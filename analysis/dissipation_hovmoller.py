#!/usr/bin/env python3
# analysis/dissipation_hovmoller.py
"""
Depth-time ("Hovmöller") heatmaps of density and TKE dissipation rate
(epsilon), for the spar vs. semisub shear-geometry runs
(runs/shear_spar, runs/shear_semisub -- see
notes/shear_geometry_interaction.md).

Unlike the single eroding-density-line view in
analysis/mixing_race_animation.py, this shows *where in the water column*
and *when* dissipation is strongest, and how that spatial pattern relates
to where/when the density field actually erodes -- i.e. a direct visual
check of the "engine" behind the mixing seen in the line plots.

epsilon is derived via utils.utils.compute_epsilon (GLS-generic, see that
docstring) and plotted on a log colour scale since it varies over several
orders of magnitude between the structure's wake and the quiescent
interior.

Usage:
    python analysis/dissipation_hovmoller.py
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
from matplotlib.colors import LogNorm

from utils.utils import open_roms_dataset, compute_epsilon, compute_time_vector

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

# (run_name, draft, label)
RUNS = [
    ("shear_spar", 90.0, "spar (draft=90 m)"),
    ("shear_semisub", 20.0, "semisub (draft=20 m)"),
]

N_DAYS = 6.0


def load_fields(run_name, draft):
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    t_days = compute_time_vector(params)
    n = np.searchsorted(t_days, N_DAYS) + 1

    rho = ds["rho"].isel(eta_rho=0, xi_rho=0).values[:n]          # (time, s_rho)
    z_rho = ds["z_rho"].isel(eta_rho=0, xi_rho=0).values          # (s_rho,)
    eps = compute_epsilon(ds, params).isel(eta_rho=0, xi_rho=0).values[:n]  # (time, s_w)
    z_w = ds["z_w"].isel(eta_rho=0, xi_rho=0).values              # (s_w,)
    ds.close()
    return t_days[:n], rho, z_rho, eps, z_w, draft


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    data = [load_fields(name, draft) for name, draft, _ in RUNS]
    eps_all = np.concatenate([d[3].ravel() for d in data])
    eps_all = eps_all[eps_all > 0]
    vmin, vmax = np.percentile(eps_all, [1, 99.5])

    rho_all = np.concatenate([d[1].ravel() for d in data])
    rho_vmin, rho_vmax = rho_all.min(), rho_all.max()

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)

    for row, ((t_days, rho, z_rho, eps, z_w, draft), (_, _, label)) in enumerate(zip(data, RUNS)):
        ax_rho = axes[row, 0]
        ax_eps = axes[row, 1]

        pc0 = ax_rho.pcolormesh(t_days, z_rho, rho.T, shading="auto",
                                 cmap="viridis", vmin=rho_vmin, vmax=rho_vmax)
        ax_rho.axhline(-draft, color="white", ls="--", lw=1.2)
        ax_rho.set_ylabel(f"{label}\ndepth (m)")
        ax_rho.set_title("density anomaly" if row == 0 else "")

        pc1 = ax_eps.pcolormesh(t_days, z_w, eps.T, shading="auto",
                                 cmap="inferno", norm=LogNorm(vmin=vmin, vmax=vmax))
        ax_eps.axhline(-draft, color="cyan", ls="--", lw=1.2,
                        label="structure draft" if row == 0 else None)
        ax_eps.set_title(r"dissipation rate $\epsilon$" if row == 0 else "")
        if row == 0:
            ax_eps.legend(fontsize=8, loc="lower right")

    for ax in axes[-1, :]:
        ax.set_xlabel("time (days)")

    cbar0 = fig.colorbar(pc0, ax=axes[:, 0], shrink=0.85, label=r"$\rho$ (kg m$^{-3}$)")
    cbar1 = fig.colorbar(pc1, ax=axes[:, 1], shrink=0.85, label=r"$\epsilon$ (m$^2$ s$^{-3}$)")

    fig.suptitle("Where and when does mixing happen?\ndensity erosion vs. dissipation rate, spar vs. semisub (sheared current)")
    out_path = os.path.join(FIGURES_DIR, "dissipation_hovmoller.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved figure to {out_path}")


if __name__ == "__main__":
    main()
