#!/usr/bin/env python3
"""
tools/prep_wind_structure_synergy.py

Experiment E: structure-induced vs. wind-induced mixing additivity/synergy
check, see notes/surface_mixing_synergy.md.

Idealised 2x2 design (structure drag on/off) x (surface wind stress
on/off), all at the same fixed background current (U0=0.1 m/s, matched
spar-like geometry: draft=90 m, str_a=1.5/90 -- same convention as
notes/floating_structure_depth_distribution.md) and the same fixed run
duration. Since the GLS k-eps/k-omega TKE equations are nonlinear in their
production terms, structure-driven production P_d and wind-driven surface
shear production P do not have to add linearly in their effect on the
stratification decay phi(t) -- this experiment checks whether they do, in
this idealised, non-rotating, tideless setting.

Wind stress is a single "moderate" magnitude, from the standard bulk
formula tau = rho_air * Cd10 * W10^2 with W10=10 m/s, Cd10=1.3e-3,
rho_air=1.225 kg/m3 -> tau=0.159 Pa; kinematic stress (ROMS convention,
tau/rho0) = 0.159/1025 = 1.55e-4 m2/s2. Applied in the same along-current
(+x) direction as the background current, constant in time (no need for a
diurnal/storm time structure here -- see notes for why this is a
reasonable idealisation).

Runs (4 total):
    synergy_struct_off_wind_off   (control: no vertical shear, no mixing expected)
    synergy_struct_on_wind_off    (structure-induced mixing only -- numerically
                                   identical to runs/geom_matched_draft90,
                                   but re-run here for a consistent frc file/
                                   binary and an exact run-to-run baseline)
    synergy_struct_off_wind_on    (wind-induced mixing only)
    synergy_struct_on_wind_on     (both; the interesting synergy/additivity case)

Usage:
    python tools/prep_wind_structure_synergy.py
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

DRAFT = 90.0
STR_A_ON = 1.5 / DRAFT
STR_A_OFF = 0.0

# Moderate 10 m/s wind, standard bulk formula (see module docstring).
RHO_AIR = 1.225
CD10 = 1.3e-3
W10 = 10.0
RHO0 = 1025.0
TAU_PA = RHO_AIR * CD10 * W10 ** 2
SUSTR_TAU0_ON = TAU_PA / RHO0
SUSTR_TAU0_OFF = 0.0


def build_runs():
    runs = []
    for struct_on in [False, True]:
        for wind_on in [False, True]:
            name = f"synergy_struct_{'on' if struct_on else 'off'}_wind_{'on' if wind_on else 'off'}"
            overrides = {
                "structure.str_a": STR_A_ON if struct_on else STR_A_OFF,
                "structure.depth_zero_below": DRAFT,
                "forcing.sustr_tau0": SUSTR_TAU0_ON if wind_on else SUSTR_TAU0_OFF,
                "forcing.svstr_tau0": 0.0,
            }
            params = {"struct_on": struct_on, "wind_on": wind_on,
                      "draft": DRAFT, "sustr_tau0": overrides["forcing.sustr_tau0"]}
            runs.append((name, overrides, params))
    # Extra diagnostic: a 10x ("storm-force", ~35 m/s equivalent) wind
    # stress with the structure off, to check whether *any* achievable
    # wind stress produces a detectable bulk mixing signal in this setup
    # (see notes/surface_mixing_synergy.md Section 3-4).
    runs.append((
        "synergy_windcheck_strong",
        {
            "structure.str_a": STR_A_OFF,
            "structure.depth_zero_below": DRAFT,
            "forcing.sustr_tau0": 10.0 * SUSTR_TAU0_ON,
            "forcing.svstr_tau0": 0.0,
        },
        {"struct_on": False, "wind_on": "storm", "draft": DRAFT,
         "sustr_tau0": 10.0 * SUSTR_TAU0_ON},
    ))
    return runs


def main():
    sweep_out_dir = "sweeps/wind_structure_synergy"
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
                         "wind / structure mixing synergy check", manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)
    print(f"sustr_tau0 (wind on) = {SUSTR_TAU0_ON:.4e} m2/s2 (tau={TAU_PA:.4f} Pa, W10={W10} m/s)")


if __name__ == "__main__":
    main()
