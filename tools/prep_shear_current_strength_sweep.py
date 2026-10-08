#!/usr/bin/env python3
"""
tools/prep_shear_current_strength_sweep.py

Shear x current-strength sweep, see notes/shear_residence_time_tradeoff.md.

Combines two earlier ideas: (1) the residence-time vs. mixing-timescale
tradeoff (notes/residence_time_tradeoff.md), which used the *uniform*-
current current-strength sweep, and (2) the sheared, surface-intensified
NCC-like current profile introduced for the shear-geometry experiment
(notes/shear_geometry_interaction.md). This sweep asks: does the
residence-time tradeoff conclusion (faster currents still deliver *more*
mixing by farm exit, despite shorter residence time) still hold once the
current is realistically sheared rather than depth-uniform?

Reuses the same shape as profiles/ncc_shear_surface_intensified.txt
(depth-mean 0.15 m/s), uniformly rescaled to 4 other depth-mean targets
(0.05, 0.10, 0.20, 0.30 m/s -- matching the uniform-current sweep's
targets exactly for a fair comparison): see
profiles/ncc_shear_scale_u{0p05,0p10,0p20,0p30}.txt. The 0.15 m/s point
reuses the already-completed runs/shear_spar and runs/shear_semisub
(not re-run here). Same matched-total_a=1.5 geometries as all prior
spar/semisub comparisons in this project.

Usage:
    python tools/prep_shear_current_strength_sweep.py
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

# tag -> (profile_file, nominal depth-mean u target [m/s])
PROFILES = {
    "0p05": ("profiles/ncc_shear_scale_u0p05.txt", 0.05),
    "0p10": ("profiles/ncc_shear_scale_u0p10.txt", 0.10),
    "0p20": ("profiles/ncc_shear_scale_u0p20.txt", 0.20),
    "0p30": ("profiles/ncc_shear_scale_u0p30.txt", 0.30),
}

GEOMETRIES = {
    "spar": {"draft": 90.0, "str_a": 1.5 / 90.0},
    "semisub": {"draft": 20.0, "str_a": 1.5 / 20.0},
}


def build_runs():
    runs = []
    for struct_name, geom in GEOMETRIES.items():
        for tag, (profile_file, u_target) in PROFILES.items():
            name = f"shearcur_{struct_name}_u{tag}"
            overrides = {
                "structure.str_a": geom["str_a"],
                "structure.depth_zero_below": geom["draft"],
                "bodyforce.mode": "profile",
                "bodyforce.profile_file": profile_file,
            }
            params = {
                "structure_type": struct_name,
                "draft": geom["draft"],
                "str_a": geom["str_a"],
                "forcing": "shear_profile",
                "u_target": u_target,
            }
            runs.append((name, overrides, params))
    return runs


def main():
    sweep_out_dir = "sweeps/shear_current_strength_sweep"
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
                         "shear x current-strength sweep", manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)


if __name__ == "__main__":
    main()
