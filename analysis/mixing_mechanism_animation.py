#!/usr/bin/env python3
# analysis/mixing_mechanism_animation.py
"""
Extended "mixing race" animation: alongside the density profile and the
phi*(t) race track (see analysis/mixing_race_animation.py), this also
animates the two ingredients of the Rennau et al. (2012) mixing
coefficient

    Gamma(z,t) = -B(z,t) / epsilon(z,t)

directly: the dissipation rate epsilon(z,t) (utils.utils.compute_epsilon)
and the buoyancy destruction of TKE B(z,t) = (g/rho0) * AKt * d(rho)/dz,
for the spar vs. semisub shear-geometry runs. Companion to
analysis/dissipation_hovmoller.py's static Hovmöller view -- this shows
the same mechanism "in motion".

Usage:
    python analysis/mixing_mechanism_animation.py
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
import matplotlib.animation as animation

from utils.utils import (
    open_roms_dataset, compute_phi, compute_time_vector, compute_epsilon, G, RHO0,
)

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

RUNS = {
    "spar (draft=90 m)": ("shear_spar", 90.0, "tab:blue"),
    "semisub (draft=20 m)": ("shear_semisub", 20.0, "tab:orange"),
}

N_DAYS = 6.0
FRAME_STRIDE = 1
GIF_FPS = 5


def load_run(run_name, draft):
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    t_days = compute_time_vector(params)

    rho = ds["rho"].isel(eta_rho=0, xi_rho=0).values          # (time, s_rho)
    z_rho = ds["z_rho"].isel(eta_rho=0, xi_rho=0).values       # (s_rho,)

    eps = compute_epsilon(ds, params).isel(eta_rho=0, xi_rho=0).values  # (time, s_w)
    drho_dz = grid.derivative(ds["rho"], "Z")
    B = ((G / RHO0) * ds["AKt"] * drho_dz).isel(eta_rho=0, xi_rho=0).values  # (time, s_w)
    z_w = ds["z_w"].isel(eta_rho=0, xi_rho=0).values           # (s_w,)

    phi, _ = compute_phi(ds, grid, params)
    phi_star = (phi / phi.isel(ocean_time=0)).values
    ds.close()

    n = np.searchsorted(t_days, N_DAYS)
    idx = np.arange(0, n, FRAME_STRIDE)
    return t_days[idx], rho[idx], z_rho, eps[idx], B[idx], z_w, phi_star[idx], draft


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    data = {label: load_run(name, draft) for label, (name, draft, _) in RUNS.items()}
    n_frames = min(len(d[0]) for d in data.values())

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    ax_rho, ax_eps = axes[0, 0], axes[0, 1]
    ax_B, ax_race = axes[1, 0], axes[1, 1]

    # --- density profile ---
    rho_lines = {}
    for label, (_, _, _, _, _, _, _, draft) in data.items():
        color = RUNS[label][2]
        (line,) = ax_rho.plot([], [], color=color, lw=2, label=label)
        rho_lines[label] = line
        ax_rho.axhspan(-draft, 0, color=color, alpha=0.08)
    all_rho = np.concatenate([d[1] for d in data.values()])
    ax_rho.set_xlim(all_rho.min() - 0.2, all_rho.max() + 0.2)
    ax_rho.set_ylim(-150, 0)
    ax_rho.set_xlabel(r"$\rho$ (kg m$^{-3}$)")
    ax_rho.set_ylabel("depth (m)")
    ax_rho.set_title("density profile")
    ax_rho.legend(fontsize=8, loc="lower right")
    ax_rho.grid(alpha=0.3)
    time_text = ax_rho.text(0.03, 0.03, "", transform=ax_rho.transAxes, fontsize=11)

    # --- dissipation rate profile (log x-axis) ---
    eps_lines = {}
    for label in data:
        color = RUNS[label][2]
        (line,) = ax_eps.plot([], [], color=color, lw=2, label=label)
        eps_lines[label] = line
    all_eps = np.concatenate([d[3].ravel() for d in data.values()])
    all_eps = all_eps[all_eps > 0]
    ax_eps.set_xscale("log")
    ax_eps.set_xlim(np.percentile(all_eps, 1), np.percentile(all_eps, 99.9))
    ax_eps.set_ylim(-150, 0)
    ax_eps.set_xlabel(r"$\epsilon$ (m$^2$ s$^{-3}$)")
    ax_eps.set_title("dissipation rate")
    ax_eps.grid(alpha=0.3)

    # --- buoyancy destruction of TKE profile ---
    B_lines = {}
    for label in data:
        color = RUNS[label][2]
        (line,) = ax_B.plot([], [], color=color, lw=2, label=label)
        B_lines[label] = line
    # Exclude the very first (transient, initial-adjustment) frame when
    # setting axis limits -- it is a one-timestep spike that would otherwise
    # swamp the meaningful later-time range.
    all_B = np.concatenate([d[4][2:].ravel() for d in data.values()])
    lo, hi = np.percentile(all_B, [0.5, 99.5])
    pad = 0.15 * (hi - lo)
    ax_B.set_xlim(lo - pad, hi + pad)
    ax_B.axvline(0, color="gray", lw=0.8)
    ax_B.set_ylim(-150, 0)
    ax_B.set_xlabel(r"$B = (g/\rho_0) A_{Kt} \partial_z\rho$ (m$^2$ s$^{-3}$)")
    ax_B.set_ylabel("depth (m)")
    ax_B.set_title("buoyancy destruction of TKE\n(negative = stratification eating turbulence)")
    ax_B.grid(alpha=0.3)

    # --- race track ---
    race_lines = {}
    race_markers = {}
    for label, (t_days, _, _, _, _, _, phi_star, _) in data.items():
        color = RUNS[label][2]
        ax_race.plot(t_days, phi_star, color=color, lw=1, alpha=0.3)
        (line,) = ax_race.plot([], [], color=color, lw=2.5, label=label)
        (marker,) = ax_race.plot([], [], "o", color=color, markersize=10)
        race_lines[label] = line
        race_markers[label] = marker
    ax_race.set_xlim(0, N_DAYS)
    ax_race.set_ylim(0, 1.02)
    ax_race.set_xlabel("time (days)")
    ax_race.set_ylabel(r"$\phi^*(t)$")
    ax_race.set_title("race to zero")
    ax_race.legend(fontsize=8)
    ax_race.grid(alpha=0.3)

    fig.suptitle(r"Anatomy of mixing: $\Gamma = -B/\epsilon$ in motion")
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    def update(frame):
        artists = []
        for label, (t_days, rho, z_rho, eps, B, z_w, phi_star, draft) in data.items():
            rho_lines[label].set_data(rho[frame], z_rho)
            eps_lines[label].set_data(np.clip(eps[frame], 1e-12, None), z_w)
            B_lines[label].set_data(B[frame], z_w)
            race_lines[label].set_data(t_days[: frame + 1], phi_star[: frame + 1])
            race_markers[label].set_data([t_days[frame]], [phi_star[frame]])
            artists += [rho_lines[label], eps_lines[label], B_lines[label],
                        race_lines[label], race_markers[label]]
        t0 = list(data.values())[0][0][frame]
        time_text.set_text(f"t = {t0:.2f} days")
        artists.append(time_text)
        return artists

    anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000 / GIF_FPS, blit=False)
    out_path = os.path.join(FIGURES_DIR, "mixing_mechanism.gif")
    anim.save(out_path, writer="pillow", fps=GIF_FPS)
    print(f"Saved animation to {out_path}")


if __name__ == "__main__":
    main()
