#!/usr/bin/env python3
"""
tools/prep_floating_depth_sweep.py

Prepare the floating-structure structured-zone-thickness sensitivity sweep
(see templates/floating_depth_sweep.yaml and
notes/floating_structure_sensitivity_analysis.md): varies
structure.depth_zero_below (expressed as a fraction of grid.H0, via a
synthetic "structure.depth_frac" sweep parameter) crossed with
structure.c4 and (grid.H0, initial.temp_zt) pairs (reusing the
pairs x parameters mechanism from prep_pycnocline_zt_sweep.py), with
bfrc_cb == CD fixed throughout (no shear between the structured and
floating zones -- required for analytic_mixing_timescale() to apply).

Two sizing modes are supported via sweep.ntimes_basis:
  - "tau_x" (default): NTIMES ~= margin * tau_x_theory / DT, where tau_x is
    the time to reach a fixed fractional reduction of phi (default
    x_frac=0.10, i.e. 10% mixing) -- used for the main sweep, since
    full-column mixing (tau_mix) may never be reached for a thin
    structured zone.
  - "tau_mix": NTIMES ~= margin * tau_mix_theory / DT (the original,
    full-completion time scale) -- used for a small, separate "long-run"
    sweep (see templates/floating_depth_longrun_sweep.yaml) sized to give
    phi(t) enough time to plateau, for the phi_inf/mixed-fraction
    diagnostic (see utils.detect_phi_plateau).

Usage:
    python tools/prep_floating_depth_sweep.py [templates/floating_depth_sweep.yaml]
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


def derive_ntimes(cfg: dict, margin: float, basis: str, x_frac: float) -> dict:
    """
    Compute tau_mix_theory/tau_x_theory for a resolved cfg and set
    cfg['time_stepping']['NTIMES'] = round_up(margin * tau / DT) to the
    nearest multiple of NHIS, where tau is tau_x_theory (basis="tau_x") or
    tau_mix_theory (basis="tau_mix"). Returns the analytic_mixing_timescale()
    dict for bookkeeping.
    """
    result = analytic_mixing_timescale(cfg, x_frac=x_frac)
    # Cast all values to plain Python floats -- np.float64 is not
    # yaml.safe_dump-serializable and would break manifest writing.
    result = {k: float(v) for k, v in result.items()}

    DT = float(cfg["time_stepping"]["DT"])
    NHIS = int(cfg["time_stepping"]["NHIS"])

    if basis == "tau_x":
        tau = result["tau_x"]
    elif basis == "tau_mix":
        tau = result["tau_mix"]
    else:
        raise ValueError(f"Unknown ntimes_basis {basis!r}: use 'tau_x' or 'tau_mix'.")

    ntimes_raw = margin * tau / DT
    n_hist_records = max(1, int(np.ceil(ntimes_raw / NHIS)))
    ntimes = n_hist_records * NHIS

    cfg["time_stepping"]["NTIMES"] = int(ntimes)
    result["NTIMES"] = int(ntimes)
    result["duration_days"] = ntimes * DT / 86400.0
    result["ntimes_basis"] = basis
    return result


def main():
    import numpy as np  # noqa: F401
    globals()["np"] = np

    sweep_path = sys.argv[1] if len(sys.argv) > 1 else "templates/floating_depth_sweep.yaml"
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
    ntimes_basis = sdef.get("ntimes_basis", "tau_x")
    x_frac = float(sdef.get("x_frac", 0.10))

    cb = fixed.get("structure.cb")
    CD = fixed.get("structure.CD")
    if cb is None or CD is None or float(cb) != float(CD):
        raise ValueError(
            "floating_depth_sweep requires sweep.fixed to set both "
            "structure.CD and structure.cb to the same value (no shear "
            "between the structured and floating zones), matching the "
            "analytic_mixing_timescale() assumption -- see "
            "notes/floating_structure_sensitivity_analysis.md."
        )

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
            param_combo = dict(param_combo)
            depth_frac = float(param_combo.pop("structure.depth_frac"))
            H0 = float(pair["H0"])
            depth_zero_below = depth_frac * H0

            combo = {
                "grid.H0": H0,
                "initial.temp_zt": pair["temp_zt"],
                "structure.depth_zero_below": depth_zero_below,
                **{k: v for k, v in fixed.items()},
                **param_combo,
            }
            # Kept only for manifest bookkeeping/plotting (not a real ROMS
            # config path -- do not pass depth_frac itself into `over`).
            combo_for_manifest = {**combo, "structure.depth_frac": depth_frac}

            over = {}
            for k, v in combo.items():
                parts = k.split(".")
                d = over
                for p in parts[:-1]:
                    d = d.setdefault(p, {})
                d[parts[-1]] = v
            cfg = deep_merge(base, over)

            run_name = format_run_name(run_name_template, cfg, combo_for_manifest, index)
            cfg.setdefault("run", {})
            cfg["run"]["name"] = run_name

            tau_info = derive_ntimes(cfg, ntimes_margin, ntimes_basis, x_frac)

            result = prepare_run_from_resolved(cfg)

            manifest_rows.append({
                "hash_exact": result["hash_exact"],
                "status": "prepared",
                "run_name": run_name,
                "resolved_config": result["resolved_config"],
                "params": combo_for_manifest,
                "d_struct": tau_info["d_struct"],
                "tau_mix_theory_days": tau_info["tau_mix"] / 86400.0,
                "tau_x_theory_days": tau_info["tau_x"] / 86400.0,
                "u_inf": tau_info["u_inf"],
                "Pstr_theory": tau_info["Pstr"],
                "NTIMES": tau_info["NTIMES"],
                "duration_days": tau_info["duration_days"],
                "ntimes_basis": tau_info["ntimes_basis"],
            })

            print(
                f"Prepared run {index}: {run_name}  "
                f"(depth_frac={depth_frac:.2f}, tau_x_theory={tau_info['tau_x']/86400.0:.3f} d, "
                f"tau_mix_theory={tau_info['tau_mix']/86400.0:.2f} d, "
                f"NTIMES={tau_info['NTIMES']}, duration={tau_info['duration_days']:.2f} d, "
                f"basis={ntimes_basis})"
            )
            index += 1

    write_manifest_yaml(os.path.join(sweep_out_dir, "manifest.yaml"), sweep_id, manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)


if __name__ == "__main__":
    main()
