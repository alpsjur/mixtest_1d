#!/usr/bin/env python3
"""
tools/prep_pycnocline_zt_sweep.py

Prepare the pycnocline-length-scale test sweep (see
templates/pycnocline_zt_sweep.yaml and notes/mixing_timescale_analysis.md
section 4.5/6): varies initial.temp_zt independently of grid.H0 via
explicit (H0, temp_zt) pairs, crossed with structure.c4, to test whether
L = sqrt(z_t*(H0-z_t)) collapses the dimensionless mixing time over a wider
range than the original zt=40-fixed, H0={75,150} test.

This mirrors tools/prep_mixing_timescale_sweep.py's NTIMES-derivation
logic (each run gets NTIMES ~= ntimes_margin * tau_mix_theory / DT), but
combos are built from sweep.pairs x cartesian(sweep.parameters) instead of
a single cartesian product over all parameters -- a full cartesian product
over H0 and temp_zt would include invalid combinations where the
thermocline sits too close to a boundary.

Usage:
    python tools/prep_pycnocline_zt_sweep.py [templates/pycnocline_zt_sweep.yaml]
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
import shutil
import filecmp

from utils.utils import ensure_dir, analytic_mixing_timescale
from prep_experiment import prepare_run_from_resolved, deep_merge
from prep_sweep import (
    cartesian_dict,
    format_run_name,
    write_manifest_yaml,
    write_manifest_csv,
)


def derive_ntimes(cfg: dict, margin: float) -> dict:
    """Same as prep_mixing_timescale_sweep.derive_ntimes."""
    result = analytic_mixing_timescale(cfg)
    result = {k: float(v) for k, v in result.items()}

    DT = float(cfg["time_stepping"]["DT"])
    NHIS = int(cfg["time_stepping"]["NHIS"])

    ntimes_raw = margin * result["tau_mix"] / DT
    n_hist_records = max(1, int(np.ceil(ntimes_raw / NHIS)))
    ntimes = n_hist_records * NHIS

    cfg["time_stepping"]["NTIMES"] = int(ntimes)
    result["NTIMES"] = int(ntimes)
    result["duration_days"] = ntimes * DT / 86400.0
    return result


def main():
    import numpy as np  # noqa: F401
    globals()["np"] = np

    sweep_path = sys.argv[1] if len(sys.argv) > 1 else "templates/pycnocline_zt_sweep.yaml"
    with open(sweep_path, "r") as f:
        sweep = yaml.safe_load(f)

    sdef = sweep["sweep"]
    sweep_id = sdef["id"]
    base_config_path = sdef["base_config"]
    apply_config_path = sdef.get("apply_config")
    pairs = sdef["pairs"]
    parameters = sdef.get("parameters", {})
    fixed = sdef.get("fixed", {})
    run_name_template = sdef["run_name_template"]
    sweep_out_dir = sdef["output_dir"]
    ntimes_margin = float(sdef.get("ntimes_margin", 2.5))

    ensure_dir(sweep_out_dir)

    sweep_yaml_src = os.path.abspath(sweep_path)
    sweep_yaml_dst = os.path.join(sweep_out_dir, "sweep.yaml")
    if os.path.exists(sweep_yaml_dst) and filecmp.cmp(sweep_yaml_src, sweep_yaml_dst, shallow=False):
        print(f"File already exists and is identical. No need to copy: {sweep_yaml_dst}")
    else:
        shutil.copy2(sweep_yaml_src, sweep_yaml_dst)
        print(f"Copied sweep definition to: {sweep_yaml_dst}")

    with open(base_config_path, "r") as f:
        base = yaml.safe_load(f) or {}
    if apply_config_path:
        with open(apply_config_path, "r") as f:
            applied = yaml.safe_load(f) or {}
        base = deep_merge(base, applied)

    manifest_rows = []
    index = 1
    for pair in pairs:
        for param_combo in cartesian_dict(parameters) if parameters else [{}]:
            combo = {
                "grid.H0": pair["H0"],
                "initial.temp_zt": pair["temp_zt"],
                **{k: v for k, v in fixed.items()},
                **param_combo,
            }

            over = {}
            for k, v in combo.items():
                parts = k.split(".")
                d = over
                for p in parts[:-1]:
                    d = d.setdefault(p, {})
                d[parts[-1]] = v
            cfg = deep_merge(base, over)

            run_name = format_run_name(run_name_template, cfg, combo, index)
            cfg.setdefault("run", {})
            cfg["run"]["name"] = run_name

            tau_info = derive_ntimes(cfg, ntimes_margin)

            result = prepare_run_from_resolved(cfg)

            manifest_rows.append({
                "hash_exact": result["hash_exact"],
                "status": "prepared",
                "run_name": run_name,
                "resolved_config": result["resolved_config"],
                "params": combo,
                "tau_mix_theory_days": tau_info["tau_mix"] / 86400.0,
                "u_inf": tau_info["u_inf"],
                "Pstr_theory": tau_info["Pstr"],
                "NTIMES": tau_info["NTIMES"],
                "duration_days": tau_info["duration_days"],
            })

            print(
                f"Prepared run {index}: {run_name}  "
                f"(tau_mix_theory={tau_info['tau_mix']/86400.0:.2f} d, "
                f"NTIMES={tau_info['NTIMES']}, duration={tau_info['duration_days']:.2f} d)"
            )
            index += 1

    write_manifest_yaml(os.path.join(sweep_out_dir, "manifest.yaml"), sweep_id, manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)


if __name__ == "__main__":
    main()
