---
name: analyze-results
description: Analyze and plot mixtest_1d output (time series, profiles, sweep comparisons). Use when asked to plot or compare ROMS results.
---

# Analyze results

Load data with `utils.utils.open_roms_dataset` (config → NetCDF → xgcm metrics dx, dy, dz, dA, dV). Reuse helpers in `utils/utils.py` rather than reimplementing.

```bash
# Single run
python analysis/prep_timeseries.py --resolved_config runs/<name>/resolved_config.yaml --variable AKt
python analysis/prep_profiles.py   --resolved_config runs/<name>/resolved_config.yaml --variable AKt

# Compare two sweeps
python analysis/wrapper.py --sweep1 sweeps/a/manifest.yaml --sweep2 sweeps/b/manifest.yaml \
    --variables AKt temp --save-dir figures
```

- `plot_sweep.py` plots one variable across all runs in a sweep.
- Background/theory notes are in `notes/`.
