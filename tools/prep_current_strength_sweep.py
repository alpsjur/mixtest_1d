#!/usr/bin/env python3
"""
tools/prep_current_strength_sweep.py

Experiment D: background current (forcing) strength sensitivity, see
notes/current_strength_sensitivity.md.

Sweeps bodyforce.U0 -- the actual UV_BODYFORCE target velocity (see
tools/make_grd.py:build_bfrc, mode "uniform": bfrc_u = BFRC_CD * U0) -- for
the spar-like (draft=90 m) and semisub-like (draft=20 m) matched-total_a
(total_a=1.5) geometries already used in
notes/floating_structure_depth_distribution.md, spanning a representative
Norwegian Coastal Current speed range.

NOTE on an earlier bug in this script: bodyforce.BFRC_U/BFRC_V are NOT fed
into the simulated body force for mode="uniform" runs -- they are only
used by utils.analytic_mixing_timescale()'s *analytic* (pre-run) tau_mix
estimate, and only for full-depth structures. An earlier version of this
script swept BFRC_U instead of U0, which (since U0 was left at its
baseline.yaml default of 0.1 m/s throughout) produced 8 runs with
*identical* actual forcing and hence identical diagnosed speeds --
caught only by diagnosing the real model "u" field in
analysis/current_strength_sensitivity.py rather than trusting the nominal
labels. See notes/current_strength_sensitivity.md for the full writeup.

The u_target=0.1 m/s point reuses the existing geom_matched_draft90 /
geom_matched_draft20 runs from the floating_structure_geometry sweep
(bodyforce.U0=0.1, the baseline.yaml default) -- not re-run here.

Usage:
    python tools/prep_current_strength_sweep.py
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

U_TARGETS = [0.05, 0.15, 0.2, 0.3]  # 0.1 already covered by the geometry sweep

GEOMETRIES = {
    "spar": {"draft": 90.0, "str_a": 1.5 / 90.0},
    "semisub": {"draft": 20.0, "str_a": 1.5 / 20.0},
}


def build_runs():
    runs = []
    for struct_name, geom in GEOMETRIES.items():
        for u_target in U_TARGETS:
            name = f"current_{struct_name}_u{u_target:.2f}".replace(".", "p")
            overrides = {
                "structure.str_a": geom["str_a"],
                "structure.depth_zero_below": geom["draft"],
                "bodyforce.U0": u_target,
                "bodyforce.V0": 0.0,
                "initial.ubar0": u_target,
            }
            params = {"structure_type": struct_name, "draft": geom["draft"],
                      "str_a": geom["str_a"], "u_target": u_target}
            runs.append((name, overrides, params))
    return runs


def main():
    sweep_out_dir = "sweeps/current_strength_sensitivity"
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
                         "current strength sensitivity sweep", manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)


if __name__ == "__main__":
    main()
