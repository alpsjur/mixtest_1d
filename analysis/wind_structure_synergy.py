#!/usr/bin/env python3
# analysis/wind_structure_synergy.py
"""
Analysis/plotting for the structure-induced vs. wind-induced mixing
additivity/synergy check -- see notes/surface_mixing_synergy.md.

Usage:
    python analysis/wind_structure_synergy.py
"""
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from utils.utils import open_roms_dataset, compute_phi, compute_time_vector

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

MAIN_RUNS = [
    ("synergy_struct_off_wind_off", "structure off, wind off (control)", "k", "-"),
    ("synergy_struct_off_wind_on", "structure off, wind on (moderate, 10 m/s)", "tab:blue", "-"),
    ("synergy_struct_on_wind_off", "structure on, wind off", "tab:orange", "-"),
    ("synergy_struct_on_wind_on", "structure on, wind on (moderate)", "tab:red", "--"),
]

# Extra diagnostic run: same as struct-off/wind-off but with a storm-force
# (10x moderate) wind stress, to check whether *any* achievable wind stress
# in this setup produces a detectable bulk mixing signal.
WINDCHECK_RUN = "synergy_windcheck_strong"


def load_phi_star(run_name: str):
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    phi, _ = compute_phi(ds, grid, params)
    phi_star = (phi / phi.isel(ocean_time=0)).values
    t_days = compute_time_vector(params)
    ds.close()
    return t_days, phi_star


def load_tke_profile(run_name: str, frac_time: float = 1.0):
    """Depth profile of TKE at a given fraction through the run (default: final time)."""
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    n_t = ds.sizes["ocean_time"]
    i = min(int(n_t * frac_time), n_t - 1)
    ds_t = ds.isel(ocean_time=i)
    tke = ds_t["tke"].mean(dim=["eta_rho", "xi_rho"])
    z_w = ds_t["z_w"].mean(dim=["eta_rho", "xi_rho"]) if "z_w" in ds_t else None
    z = z_w.values if z_w is not None else np.arange(tke.sizes["s_w"])
    ds.close()
    return z, tke.values


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    ax = axes[0]
    for run_name, label, color, ls in MAIN_RUNS:
        t_days, phi_star = load_phi_star(run_name)
        ax.plot(t_days, phi_star, color=color, ls=ls, label=label)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("time (days)")
    ax.set_ylabel(r"$\phi^*(t) = \phi(t)/\phi(0)$")
    ax.set_title("Stratification decay: structure vs. wind forcing")
    ax.legend(fontsize=7)

    ax2 = axes[1]
    for run_name, label, color in [
        ("synergy_struct_off_wind_off", "wind off", "k"),
        ("synergy_struct_off_wind_on", "wind on (moderate, 10 m/s)", "tab:blue"),
        (WINDCHECK_RUN, "wind on (storm, ~35 m/s equiv., 10x stress)", "tab:red"),
    ]:
        z, tke = load_tke_profile(run_name)
        ax2.plot(tke, z, color=color, label=label)
    ax2.set_xscale("log")
    ax2.set_xlabel("TKE (m2/s2)")
    ax2.set_ylabel("z (m, 0=surface)")
    ax2.set_ylim(-30, 1)
    ax2.set_title("Near-surface TKE profile at end of run\n(structure off; wind-driven TKE stays confined near surface)")
    ax2.legend(fontsize=7)

    fig.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "wind_structure_synergy.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
