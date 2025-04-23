# spectral-grid-inference
A playground for fitting spectra with grids of model spectra, and testing ideas.

1. Install the dependencies:

```
uv sync --all-extras --dev
```


2. Make the spectroscopic grid values by executing `scripts/Make-grid-values.ipynb`.

3. Run `scripts/make_korg_grid.jl` with the outputted files:

```
julia scripts/make_korg_grid.jl data/grid_params.csv
julia scripts/make_korg_grid.jl data/test_params.csv
```

4. Run / play with `scripts/Fit-grid-interpolate.ipynb`
