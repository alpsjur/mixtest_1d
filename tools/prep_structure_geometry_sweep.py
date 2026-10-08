#!/usr/bin/env python3
"""
tools/prep_structure_geometry_sweep.py

Prepare the floating-structure depth-distribution sweeps (spar vs
semisubmersible vs full-depth "monopile-like" reference), see
notes/floating_structure_depth_distribution.md.

Unlike tools/prep_sweep.py's cartesian product, the runs here are defined
as an explicit list of (name, overrides) pairs -- draft and str_a are not
independent axes (str_a is derived from draft to hold the total
depth-integrated frontal area, structure.str_a * structure.depth_zero_below,
fixed or at a prescribed ratio), so a cartesian sweep is the wrong tool.

Usage:
    python tools/prep_structure_geometry_sweep.py
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


def build_runs():
    """
    Returns a list of (run_name, overrides_dict, params_dict_for_manifest)
    tuples. overrides use dotted keys (same convention as prep_sweep.py's
    sweep YAML "parameters" keys).

    total_a = structure.str_a * structure.depth_zero_below -- the
    depth-integrated frontal area density, m^-1 * m = dimensionless
    (equivalent to a depth-uniform str_a over a 1 m column, times draft).

    total_a = 1.5 is the "matched" reference, equal to the existing
    validated mixing_timescale.yaml full-depth case (str_a=0.01 x H0=150m).

    "Realistic" total_a values use the *ratio* of real frontal areas
    between structure types at a common farm spacing (500x500 m cell,
    the README worked example), derived from public geometry:
      - monopile (bottom-fixed reference): d=8 m, full column (H0=150 m)
        -> a=8/250000=3.2e-5, total_a_mono=3.2e-5*150=4.8e-3
      - spar (Hywind Tampen: draft=90 m, diameter=14.7 m, uniform cylinder
        approximation): a=14.7/250000=5.88e-5, total_a_spar=5.88e-5*90=5.29e-3
      - semisubmersible (UMaine VolturnUS-S for IEA 15MW: draft=20 m,
        3 offset columns 12.5 m diameter + 1 central 10 m diameter;
        effective frontal width averaged over current heading approximated
        as ~2 column-diameters=25 m, a simplification -- see note):
        a=25/250000=1.0e-4, total_a_semisub=1.0e-4*20=2.0e-3
    These are rescaled by (total_a_type / total_a_mono) * 1.5 to keep the
    same absolute scale as the "matched" sweep (for comparable run
    durations), while preserving the *realistic* ratio between types.
    """
    runs = []

    # --- A: draft sweep at fixed total_a=1.5 (pycnocline fixed at zt=40 m) ---
    total_a_matched = 1.5
    drafts = [15, 20, 30, 40, 50, 65, 90, 120, 150]
    for draft in drafts:
        str_a = total_a_matched / draft
        name = f"geom_matched_draft{draft}"
        overrides = {
            "structure.str_a": str_a,
            "structure.depth_zero_below": float(draft),
        }
        params = {"draft": draft, "str_a": str_a, "total_a": total_a_matched,
                  "temp_zt": 40.0, "temp_dT": 10.0}
        runs.append((name, overrides, params))

    # --- B: realistic total_a ratios (spar vs semisub vs the matched
    #         full-depth/monopile reference already in set A at draft=150) ---
    total_a_mono = 4.8e-3
    total_a_spar_real = 5.29e-3 / total_a_mono * total_a_matched
    total_a_semisub_real = 2.0e-3 / total_a_mono * total_a_matched

    for name_suffix, draft, total_a in [
        ("spar", 90, total_a_spar_real),
        ("semisub", 20, total_a_semisub_real),
    ]:
        str_a = total_a / draft
        name = f"geom_realistic_{name_suffix}"
        overrides = {
            "structure.str_a": str_a,
            "structure.depth_zero_below": float(draft),
        }
        params = {"draft": draft, "str_a": str_a, "total_a": total_a,
                  "temp_zt": 40.0, "temp_dT": 10.0}
        runs.append((name, overrides, params))

    # --- C: seasonality x structure type (spar draft=90, semisub draft=20,
    #        both at total_a=1.5 "matched", already have zt=40/dT=10 from A) ---
    seasons = [
        ("summer", 15.0, 8.0),
        ("winter", 60.0, 3.0),
    ]
    for struct_name, draft in [("spar", 90), ("semisub", 20)]:
        str_a = total_a_matched / draft
        for season_name, zt, dT in seasons:
            name = f"season_{season_name}_{struct_name}"
            overrides = {
                "structure.str_a": str_a,
                "structure.depth_zero_below": float(draft),
                "initial.temp_zt": zt,
                "initial.temp_dT": dT,
            }
            params = {"draft": draft, "str_a": str_a, "total_a": total_a_matched,
                      "temp_zt": zt, "temp_dT": dT, "season": season_name,
                      "structure_type": struct_name}
            runs.append((name, overrides, params))

    return runs


def main():
    sweep_out_dir = "sweeps/floating_structure_geometry"
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
                         "floating structure geometry/seasonality sweep", manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)


if __name__ == "__main__":
    main()
