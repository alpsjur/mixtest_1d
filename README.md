# mixtest_1d

**Idealized 1D ROMS application for testing turbulence closure / mixing schemes.**

The project provides a fully scripted workflow for running the [ROMS](https://www.myroms.org/) ocean model in a single-column (1D) setup: building input files, running simulations, sweeping over parameter combinations, and analysing/plotting results.

---

## Overview

The physical setup is a single-column ocean (uniform depth, periodic horizontal boundaries) forced by a body force in the x-direction. The goal is to isolate and test turbulence closures, particularly the GLS (Generic Length Scale) scheme, as well as a novel **structure drag / mixing** parameterization that represents the effect of subgrid-scale structures (e.g. wind turbine foundations) on flow and turbulence.

### Key concepts

| Concept | Description |
|---|---|
| **GLS closure** | General turbulence closure in ROMS, configurable as k-ε, k-ω, GEN, etc. |
| **Structure area density** (`str_a`) | Frontal area per unit volume of in-water structures (m⁻¹) |
| **Structure drag** | Quadratic drag exerted by structures on the flow |
| **Structure production** | Turbulence kinetic energy produced by wakes behind structures |

---

## Repository layout

```
mixtest_1d/
├── configs/
│   ├── baseline.yaml          # Default parameter set — all runs start from here
│   └── variants/              # Per-experiment overrides (merged on top of baseline)
├── templates/
│   ├── mixtest_1d.in.j2       # Jinja2 template for the ROMS input file
│   ├── k-e_sweep.yaml         # Example sweep definition (k-epsilon variations)
│   └── gen_sweep.yaml         # Example sweep definition (GEN closure variations)
├── tools/
│   ├── prep_experiment.py     # Build input files for a single run
│   ├── run_experiment.py      # Execute a single ROMS run
│   ├── make_grd.py            # Generate the ROMS grid NetCDF
│   ├── make_ini.py            # Generate the initial conditions NetCDF
│   ├── prep_sweep.py          # Prepare a cartesian parameter sweep
│   └── run_sweep.py           # Execute all runs in a prepared sweep
├── analysis/
│   ├── prep_timeseries.py     # Volume-averaged time series of a variable
│   ├── prep_profiles.py       # Horizontally averaged vertical profile
│   ├── plot_sweep.py          # Plot a variable across all runs in a sweep
│   └── wrapper.py             # Top-level script: compare two sweeps side by side
├── tests/
│   ├── run_tests.py           # Test runner (auto-discovers test_*.py)
│   ├── test_UV_BODYFORCE.py   # Validates body-force driven acceleration
│   ├── test_STRUCTURE_DRAG.py # Validates quadratic structure drag
│   └── test_STRUCTURE_PRODUCTION.py  # Validates TKE production by structures
├── utils/
│   └── utils.py               # Shared utilities (YAML I/O, ROMS metrics, dataset loader)
├── roms/                      # ROMS executable and supporting files (do not modify)
└── environment.yml            # Conda environment specification
```

---

## Setup

### 1. Create the conda environment

```bash
conda env create -f environment.yml
conda activate mixtest_1d   # or whatever name is in the yml
```

### 2. ROMS executable

The compiled ROMS executable is expected at `roms/romsS`.  
If you need to (re-)compile ROMS, see `roms/build_roms.sh`.

---

## Workflow

### Single experiment

**Prepare** (creates `runs/<name>/` with grid, IC, and ROMS `.in` file):
```bash
python tools/prep_experiment.py configs/baseline.yaml [configs/variants/my_variant.yaml]
```

**Run:**
```bash
python tools/run_experiment.py runs/<name>/resolved_config.yaml
```

The run produces:
- `runs/<name>/output/mixtest_1d_his.nc` — ROMS history file
- `runs/<name>/logs/simulation.log` — ROMS stdout/stderr
- `runs/<name>/logs/status.yaml` — machine-readable run status

### Parameter sweep

1. Copy and edit a sweep template (e.g. `templates/k-e_sweep.yaml`) into `sweeps/`.
2. **Prepare** all runs:
   ```bash
   python tools/prep_sweep.py sweeps/my_sweep.yaml
   ```
   This creates one run directory per parameter combination and writes `sweeps/my_sweep/manifest.yaml`.
3. **Run** all prepared simulations:
   ```bash
   python tools/run_sweep.py sweeps/my_sweep/manifest.yaml
   ```
   Runs that are already marked `done` are skipped automatically.

### Running tests

```bash
# Use existing model output (runs model first if output is missing):
python tests/run_tests.py

# Force a fresh model run before each test:
python tests/run_tests.py --run-model

# Run a specific test:
python tests/run_tests.py --tests test_UV_BODYFORCE
```

### Analysis

**Time series or profile for a single run:**
```bash
python analysis/prep_timeseries.py --resolved_config runs/baseline/resolved_config.yaml --variable AKt
```

**Compare two sweeps side by side:**
```bash
python analysis/wrapper.py \
    --sweep1 sweeps/k-e_variations/manifest.yaml \
    --sweep2 sweeps/gen_variations/manifest.yaml \
    --variables AKt temp \
    --save-dir figures
```

---

## Configuration

All configuration lives in YAML files. Runs are built by **deep-merging** `baseline.yaml` with an optional variant file — the variant only needs to contain keys that differ from the baseline.

### Key config sections

| Section | Purpose |
|---|---|
| `run.name` | Name of the run directory under `runs/` |
| `grid` | Domain size (Lm, Mm, N levels), depth H0, grid spacing DX/DY |
| `vertical` | S-coordinate stretching parameters (THETA_S, THETA_B, HC) |
| `time_stepping` | NTIMES, DT (seconds), NHIS (output interval) |
| `GLS` | Turbulence closure coefficients (see ROMS manual) |
| `phys` | Coriolis F0, bottom drag RDRG2 |
| `structure` | `str_a` (area density m⁻¹), `CD` (drag coefficient), `c4` (production coefficient), `depth_zero_below` |
| `bodyforce` | Horizontal body force amplitude and ramp timing |
| `initial` | Initial temperature/salinity profile parameters |
| `files` | NetCDF file names for grid, IC, and history output |

---

## How key components fit together

```
baseline.yaml ─┐
variant.yaml  ─┴─► prep_experiment.py ──► make_grd.py  → grid NetCDF
                                      ──► make_ini.py  → IC NetCDF
                                      ──► mixtest_1d.in.j2 → ROMS input file
                                      ──► resolved_config.yaml

resolved_config.yaml ──► run_experiment.py ──► roms/romsS → history NetCDF

history NetCDF + resolved_config.yaml ──► open_roms_dataset()
                                       ──► prep_timeseries / prep_profiles
                                       ──► plot_sweep / wrapper.py
```

`utils/utils.py` is the shared foundation used by all scripts:
- `load_yaml` / `save_yaml` / `ensure_dir` — basic file I/O helpers
- `compute_stretching` / `compute_depths` — ROMS S-coordinate geometry
- `prep_ds` — attaches xgcm grid metrics (dx, dy, dz, dA, dV) to a dataset
- `open_roms_dataset` — one-call shortcut: load config → open NetCDF → prep_ds

---

## Tests

Each test verifies one physical mechanism against an analytical solution:

| Test | Physics | Pass criterion |
|---|---|---|
| `test_UV_BODYFORCE` | Uniform body force → linear acceleration u(t) = F·t | All grid points match u_analytical with rtol=1e-5 |
| `test_STRUCTURE_DRAG` | Body force balanced by structure drag → u(t) = √(F/α)·tanh(t√(Fα)) | All grid points match u_analytical with rtol=1e-4 |
| `test_STRUCTURE_PRODUCTION` | Structure TKE production balances dissipation in steady state | Domain-integrated ε matches Pd within 2% |
