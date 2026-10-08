#!/usr/bin/env python3
"""
tools/prep_draft_shear_sweep.py

"Goldilocks draft" sweep: under the surface-intensified NCC-like shear
profile (profiles/ncc_shear_surface_intensified.txt, see
notes/shear_geometry_interaction.md), sweep structure.depth_zero_below
("draft") continuously, keeping structure.str_a (drag per unit depth)
*fixed* -- i.e. a deeper draft means more total drag, as it would for a
real structure of some fixed cross-sectional drag density simply
extended to different depths. (This is a deliberately different
convention from the "same total depth-integrated drag budget"
total_a-matching used for the spar-vs-semisub comparisons elsewhere in
this project -- here we ask how a shear profile interacts with draft
depth alone, independent of any drag-budget normalization.) Question: is
there a draft depth that is *worst* (fastest mixing / highest impact) or
*best* (slowest mixing / least impact) for this current shear profile,
rather than mixing speed simply increasing or decreasing monotonically
with draft?

Reuses nothing from prior sweeps -- draft=90 is re-run here with the new,
realistically-derived str_a (does not exactly match the older
runs/shear_spar, which used an arbitrary total_a-matched str_a).

Usage:
    python tools/prep_draft_shear_sweep.py
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

# Fixed structure.str_a (drag per unit depth), *not* total_a-matched --
# i.e. a deeper draft means more total drag, as it would for a real
# structure of fixed cross-sectional drag density extended to different
# depths (this is a deliberately different convention from the "same
# total depth-integrated drag budget" total_a-matching used for the
# spar-vs-semisub comparisons elsewhere in this project).
#
# Value derived from *real* structure geometry rather than chosen
# arbitrarily, following the exact methodology already established in
# notes/floating_structure_depth_distribution.md Sec 2B (Hywind Tampen
# spar: 14.7 m effective hull diameter, 500x500 m turbine spacing):
#   a_real = diameter / spacing^2 = 14.7 / 500^2 = 5.88e-5 m^-1
# rescaled by the same factor already used throughout this project to
# bring real (tiny) frontal-area densities into a tractable-mixing-
# timescale regime (total_a=1.5 for the monopile/full-depth reference,
# vs. its real total_a = 8/500^2 * 150 = 4.8e-3):
#   scale = 1.5 / 4.8e-3 = 312.5  ->  str_a = 5.88e-5 * 312.5 = 0.018375
# i.e. this *is* (a rescaled version of) a real spar's frontal-area
# density, not an arbitrary round number -- draft is then varied on top
# of that fixed, physically-anchored drag density.
STR_A = 0.018375

DRAFTS = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]


def build_runs():
    runs = []
    for draft in DRAFTS:
        name = f"draft_shear_d{draft:03d}"
        overrides = {
            "structure.str_a": STR_A,
            "structure.depth_zero_below": float(draft),
            "bodyforce.mode": "profile",
            "bodyforce.profile_file": PROFILE_FILE,
        }
        params = {"draft": float(draft), "str_a": STR_A,
                  "total_a": STR_A * draft, "forcing": "shear_profile"}
        runs.append((name, overrides, params))
    return runs


def main():
    sweep_out_dir = "sweeps/draft_shear_sweep"
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
                         "Goldilocks draft sweep under shear profile", manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)


if __name__ == "__main__":
    main()
