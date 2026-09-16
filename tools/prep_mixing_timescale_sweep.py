#!/usr/bin/env python3
"""
tools/prep_mixing_timescale_sweep.py

Prepare the mixing-timescale sensitivity sweep (see
templates/mixing_timescale_sweep.yaml and the "Mixing time scales" notes,
following Carpenter et al. 2016).

This is a thin wrapper around tools/prep_sweep.py's cartesian-sweep
machinery, with one addition: before each run is prepared, NTIMES is
derived from analytic_mixing_timescale() so that every parameter
combination gets a run length proportional to its own predicted tau_mix,
instead of a single NTIMES sized for the slowest (or wastefully long for
the fastest) case.

Usage:
    python tools/prep_mixing_timescale_sweep.py [templates/mixing_timescale_sweep.yaml]
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
    """
    Compute tau_mix_theory for a resolved cfg and set
    cfg['time_stepping']['NTIMES'] = round_up(margin * tau_mix_theory / DT)
    to the nearest multiple of NHIS (so the last output lands on the last
    step). Returns the analytic_mixing_timescale() dict for bookkeeping.
    """
    result = analytic_mixing_timescale(cfg)
    # Cast all values to plain Python floats -- np.float64 is not
    # yaml.safe_dump-serializable and would break manifest writing.
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
    import numpy as np  # noqa: F401  (imported here so derive_ntimes can use np)
    globals()["np"] = np

    sweep_path = sys.argv[1] if len(sys.argv) > 1 else "templates/mixing_timescale_sweep.yaml"
    with open(sweep_path, "r") as f:
        sweep = yaml.safe_load(f)

    sdef = sweep["sweep"]
    sweep_id = sdef["id"]
    base_config_path = sdef["base_config"]
    apply_config_path = sdef.get("apply_config")
    parameters = sdef["parameters"]
    fixed = sdef.get("fixed", {})
    run_name_template = sdef["run_name_template"]
    sweep_out_dir = sdef["output_dir"]
    ntimes_margin = float(sdef.get("ntimes_margin", 2.5))

    ensure_dir(sweep_out_dir)

    # Copy the sweep YAML into the sweep output directory (overwrite if changed)
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
    for param_combo in cartesian_dict(parameters):
        combo = {**fixed, **param_combo}
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

        # Derive NTIMES from the analytic mixing time scale before writing
        # any run files -- avoids a fixed NTIMES that is far too short/long
        # for individual combinations.
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
