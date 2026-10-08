#!/usr/bin/env python3
"""
tools/prep_shear_geometry_sweep.py

Experiment: does a surface-intensified (sheared) background current change
the spar-vs-semisubmersible mixing comparison, compared to the uniform
current used in notes/floating_structure_depth_distribution.md and
notes/current_strength_sensitivity.md?

Motivation: the Norwegian Coastal Current is surface-intensified (strongest
in the upper ~10-40 m, decaying with depth -- freshwater-plume/Ekman-like
structure). A semisubmersible's drag (structure.str_a) is concentrated in
a shallow draft (~20 m) where such a profile is strongest, while a spar's
drag is spread over a much deeper draft (~90 m) including water depths
where the profile has already decayed close to its background value. Since
bodyforce.mode="profile" + structure.str_a are both depth-dependent, this
could make the two geometries respond differently to shear even though
they were built to be "matched" (same depth-integrated total_a=1.5) under
a uniform current.

This only prepares 2 *new* runs (one sheared run per geometry). The
matching uniform-current baselines already exist from the current-strength
sweep (bodyforce.U0=0.15 m/s, depth-mean-equivalent to the shear profile
used here -- see profiles/ncc_shear_surface_intensified.txt) and are
reused directly:
    runs/current_spar_u0p15
    runs/current_semisub_u0p15

Usage:
    python tools/prep_shear_geometry_sweep.py
"""
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

import yaml

from utils.utils import ensure_dir
from prep_experiment import prepare_run_from_resolved, deep_merge
from prep_sweep import write_manifest_yaml, write_manifest_csv, set_by_dotted_key

BASE_CONFIG = "configs/baseline.yaml"
APPLY_CONFIG = "configs/variants/floating_structure_mixing.yaml"

PROFILE_FILE = "profiles/ncc_shear_surface_intensified.txt"

# Same matched-total_a=1.5 geometries as the geometry/current-strength sweeps.
GEOMETRIES = {
    "spar": {"draft": 90.0, "str_a": 1.5 / 90.0},
    "semisub": {"draft": 20.0, "str_a": 1.5 / 20.0},
}

# Reused uniform-current baselines (already run, depth-mean = 0.15 m/s,
# matching the shear profile's thickness-weighted column mean):
UNIFORM_BASELINE_RUNS = {
    "spar": "current_spar_u0p15",
    "semisub": "current_semisub_u0p15",
}


def build_runs():
    runs = []
    for struct_name, geom in GEOMETRIES.items():
        name = f"shear_{struct_name}"
        overrides = {
            "structure.str_a": geom["str_a"],
            "structure.depth_zero_below": geom["draft"],
            "bodyforce.mode": "profile",
            "bodyforce.profile_file": PROFILE_FILE,
        }
        params = {
            "structure_type": struct_name,
            "draft": geom["draft"],
            "str_a": geom["str_a"],
            "forcing": "shear_profile",
            "uniform_baseline_run": UNIFORM_BASELINE_RUNS[struct_name],
        }
        runs.append((name, overrides, params))
    return runs


def main():
    sweep_out_dir = "sweeps/shear_geometry_sweep"
    ensure_dir(sweep_out_dir)

    with open(BASE_CONFIG, "r") as f:
        base = yaml.safe_load(f) or {}
    with open(APPLY_CONFIG, "r") as f:
        applied = yaml.safe_load(f) or {}
    base = deep_merge(base, applied)

    manifest_rows = []
    for index, (run_name, overrides, params) in enumerate(build_runs(), start=1):
        over = {}
        for k, v in overrides.items():
            set_by_dotted_key(over, k, v)
        cfg = deep_merge(base, over)
        cfg.setdefault("run", {})
        cfg["run"]["name"] = run_name

        result = prepare_run_from_resolved(cfg)

        manifest_rows.append({
            "hash_exact": result["hash_exact"],
            "status": "prepared",
            "run_name": run_name,
            "resolved_config": result["resolved_config"],
            "params": params,
        })
        print(f"Prepared run {index}: {run_name}  (params={params})")

    write_manifest_yaml(os.path.join(sweep_out_dir, "manifest.yaml"),
                         "shear vs uniform geometry comparison sweep", manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)


if __name__ == "__main__":
    main()
