---
name: run-experiment
description: Prepare and run a single 1D ROMS experiment or a parameter sweep. Use when asked to run, configure, or sweep mixtest_1d simulations.
---

# Run experiments

Config = `configs/baseline.yaml` deep-merged with an optional variant from `configs/variants/` (variants list only differing keys).

## Single run
```bash
python tools/prep_experiment.py configs/baseline.yaml [configs/variants/<variant>.yaml]
python tools/run_experiment.py runs/<name>/resolved_config.yaml
```
Outputs: `runs/<name>/output/mixtest_1d_his.nc`, `logs/simulation.log`, `logs/status.yaml`.

## Sweep
1. Copy a template (`templates/k-e_sweep.yaml`, `gen_sweep.yaml`) to `sweeps/<name>.yaml`; set `parameters` (dotted keys, e.g. `structure.c4: [0.6, 1.4]`), `run_name_template`, `output_dir`.
2. `python tools/prep_sweep.py sweeps/<name>.yaml`
3. `python tools/run_sweep.py sweeps/<name>/manifest.yaml` (skips runs marked `done`).

## Notes
- Never hand-edit generated `.in` files; change config or `templates/mixtest_1d.in.j2`.
- ROMS executable lives in `roms-related/` — don't modify it.
- Activate the conda env from `environment.yml` first.
