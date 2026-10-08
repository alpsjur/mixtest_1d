#!/usr/bin/env python3
# analysis/gamma_hovmoller.py
"""
Depth-time heatmap of the diagnosed local mixing coefficient

    Gamma(z,t) = -B(z,t) / epsilon(z,t)

for the spar vs. semisub shear-geometry runs (runs/shear_spar,
runs/shear_semisub -- see notes/shear_geometry_interaction.md), companion
to analysis/dissipation_hovmoller.py (which showed B and epsilon
separately) and analysis/gamma_observational_check.py (which only
produced a single bulk/epsilon-weighted Gamma number per run). This shows
the full depth-time *field*, directly answering "where/when is Gamma
large" rather than collapsing it to one number.

Gamma is only meaningful where epsilon is not vanishingly small (the ratio
is numerically unstable/meaningless in the quiescent deep interior where
both B and epsilon are ~0) -- following
analysis/gamma_observational_check.py's approach, cells with epsilon below
a small percentile-based floor are masked out (shown as grey) rather than
plotted as noisy, meaningless Gamma values.

Usage:
    python analysis/gamma_hovmoller.py
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

from utils.utils import open_roms_dataset, compute_epsilon, compute_time_vector, G, RHO0

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

RUNS = [
    ("shear_spar", 90.0, "spar (draft=90 m)"),
    ("shear_semisub", 20.0, "semisub (draft=20 m)"),
]

# Schultze et al. (2020) field/LES bulk mixing-efficiency estimate for
# monopile wakes -- same observational benchmark as
# analysis/gamma_observational_check.py, shown here as context for the
# colour scale (Gamma~0.1 is "typical"/observed, not an upper bound).
SCHULTZE_GAMMA_RANGE = (0.08, 0.14)

N_DAYS = 6.0
EPS_FLOOR_PERCENTILE = 20.0  # mask cells with epsilon below this percentile


def load_fields(run_name, draft):
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    t_days = compute_time_vector(params)
    n = np.searchsorted(t_days, N_DAYS) + 1

    drho_dz = grid.derivative(ds["rho"], "Z")
    B = (G / RHO0) * ds["AKt"] * drho_dz
    eps = compute_epsilon(ds, params)
    Gamma = -B / eps

    Gamma = Gamma.isel(eta_rho=0, xi_rho=0).values[:n]
    eps_v = eps.isel(eta_rho=0, xi_rho=0).values[:n]
    z_w = ds["z_w"].isel(eta_rho=0, xi_rho=0).values
    ds.close()
    return t_days[:n], Gamma, eps_v, z_w, draft


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    data = [load_fields(name, draft) for name, draft, _ in RUNS]

    eps_floor = np.percentile(
        np.concatenate([d[2].ravel() for d in data]), EPS_FLOOR_PERCENTILE
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)

    pc = None
    for ax, (t_days, Gamma, eps_v, z_w, draft), (name, _, label) in zip(axes, data, RUNS):
        masked = np.ma.masked_where(eps_v < eps_floor, Gamma)
        pc = ax.pcolormesh(t_days, z_w, masked.T, shading="auto",
                            cmap="RdBu_r", vmin=-0.3, vmax=0.3)
        ax.axhline(-draft, color="black", ls="--", lw=1.2, label="structure draft")
        ax.set_title(label)
        ax.set_xlabel("time (days)")
        ax.set_facecolor("lightgray")
        ax.legend(fontsize=8, loc="lower right")

    axes[0].set_ylabel("depth (m)")
    cbar = fig.colorbar(pc, ax=axes, shrink=0.85, label=r"$\Gamma = -B/\epsilon$")
    cbar.ax.axhspan(*SCHULTZE_GAMMA_RANGE, color="green", alpha=0.3)
    cbar.ax.text(3.3, np.mean(SCHULTZE_GAMMA_RANGE), "Schultze\net al. (2020)",
                 fontsize=7, va="center", color="darkgreen", transform=cbar.ax.transData,
                 clip_on=False)

    fig.suptitle(r"Local mixing coefficient $\Gamma(z,t)$"
                 "\n(grey = epsilon too small for a meaningful ratio)")
    fig.subplots_adjust(top=0.82, right=0.86, wspace=0.08)
    out_path = os.path.join(FIGURES_DIR, "gamma_hovmoller.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
