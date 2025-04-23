using Distributed

# Add processes if not already added (one per available core)
if nprocs() == 1
    addprocs(Sys.CPU_THREADS - 1)  # Leave one core for system processes
    @info "Added $(nprocs() - 1) worker processes"
end

# Load packages on all processes
@everywhere using Pkg
@everywhere Pkg.activate(".")  # Ensure the same environment is used on all workers
@everywhere begin
    # Import necessary packages
    using Base.Threads
    using Korg
    using CSV
    using DataFrames
end

# Define your function on all processes
@everywhere begin
    # load the linelist
    linelist = Korg.get_GES_linelist()

    # set up sparse matrix that applies the LSF and resamples to the RVS wavelength
    wvln_lim = (8460.0, 8700.0)
    synthesis_wvln = wvln_lim[1]:0.01:wvln_lim[2]
    rvs_wvln = wvln_lim[1]:0.2:wvln_lim[2]
    LSF = Korg.compute_LSF_matrix(synthesis_wvln, rvs_wvln, 11_500)

    function synthesize_spectrum(output_path, params)
        # Make a dictionary of keys and values for keys not in the stellar parameters:
        stellar_params = ["Teff", "logg", "m_H", "a_H", "vmic", "vsini"]
        elems = Dict(key => value for (key, value) in params if !(key in stellar_params) && key != "index")

        println("Synthesizing spectrum for parameters: $params")
        println("Elements: $elems")

        # construct vector of A(X)-format abundances
        A_X = format_A_X(params["m_H"], params["a_H"], elems)
        vmic = get(params, "vmic", 0.0)  # Default to 0 km/s
        vsini = get(params, "vsini", 0.0)  # Default to 0 km/s

        # create model atmosphere by interpolating from one of a few grid
        # (RGB will be normal SDSS MARCS grid)
        atm = interpolate_marcs(params["Teff"], params["logg"], A_X)

        # this returns a "solution" object with lots of data
        spectrum = synthesize(
            atm, linelist, A_X, synthesis_wvln,
            vmic=vmic, use_MHD_for_hydrogen_lines=false,
            molecular_cross_sections=[],
        )

        # continuum normalized flux:
        rectified_flux = spectrum.flux ./ spectrum.cntm
        if vsini > 0.0
            rectified_flux = apply_rotation(rectified_flux, spectrum.wavelengths, vsini)
        end
        final_spectrum = LSF * rectified_flux

        # Save the spectrum to a CSV file (column: wavelength, flux, continuum)
        specidx = lpad(string(Int(params["index"])), 5, "0")
        output_file = joinpath(output_path, "spectrum_$(specidx).csv")
        CSV.write(output_file, DataFrame(wavelength=rvs_wvln, relative_flux=final_spectrum))

        return final_spectrum
    end
end

function read_stellar_parameters(csv_path)
    # Read the CSV file containing the grid of parameters
    df = CSV.read(csv_path, DataFrame)

    # Verify that the minimum required columns exist
    required_cols = ["Teff", "logg", "m_H", "a_H"]
    missing_cols = filter(col -> !(col in names(df)), required_cols)

    if !isempty(missing_cols)
        error("Missing required columns in CSV: $(join(missing_cols, ", "))")
    end

    # Optional columns to check for
    optional_cols = ["vmic", "vsini"]

    # Create vector of dictionaries
    data = Vector{Dict{String,Float64}}(undef, nrow(df))

    # Identify element abundance columns (any column that's not required or optional)
    all_standard_cols = [required_cols; optional_cols]
    element_cols = filter(col -> !(col in all_standard_cols), names(df))

    # Convert each row to a dictionary
    for (i, row) in enumerate(eachrow(df))
        # Initialize with required parameters
        row_params = Dict{String,Any}(
            "index" => i,
            "Teff" => row.Teff,
            "logg" => row.logg,
            "m_H" => row.m_H,
            "a_H" => row.a_H,
        )

        # Add optional parameters if they exist
        for col in optional_cols
            if col in names(df)
                row_params[col] = row[col]
            end
        end

        # Add element abundances
        for element in element_cols
            row_params[element] = row[element]
        end

        data[i] = row_params
    end

    return data
end

# Read command line arguments
if length(ARGS) < 1
    println("ERROR: Usage: julia make_korg_grid.jl <path_to_grid_csv_file>")
    exit(1)
end

csv_path = ARGS[1]
println("Reading stellar parameters from: $csv_path")

# Extract directory and base filename without extension
csv_basename = basename(csv_path)
filename_no_ext = replace(csv_basename, r".csv$" => "")

# Create output directory path
output_path = joinpath(dirname(csv_path), filename_no_ext)

# Create the directory if it doesn't exist
if !isdir(output_path)
    mkpath(output_path)
    println("Created output directory: $output_path")
else
    println("Using existing output directory: $output_path")
end

param_grid = read_stellar_parameters(csv_path)
println("Read $(length(param_grid)) rows of stellar parameters.")

results = pmap(params -> synthesize_spectrum(output_path, params), param_grid)

# If you need to clean up workers when done
# rmprocs(workers())

