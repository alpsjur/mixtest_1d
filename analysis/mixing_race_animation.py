#!/usr/bin/env python3
# analysis/mixing_race_animation.py
"""
"Mixing race" animation -- a fun visualization companion to
notes/shear_geometry_interaction.md.

Animates the horizontally-averaged density profile eroding over time,
side by side for the spar (runs/shear_spar) and semisubmersible
(runs/shear_semisub) geometries under the same surface-intensified
shear current, together with a "race track" showing phi*(t) for both
with a moving marker -- so you can watch the semisubmersible catch up
(partially) once a realistic shear profile is used, instead of a uniform
current (see notes/shear_geometry_interaction.md for the quantitative
result: semisub t(phi*=0.5) drops from 8.73 to 2.55 days).

Usage:
    python analysis/mixing_race_animation.py
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

from utils.utils import open_roms_dataset, compute_phi, compute_time_vector

RUNS_DIR = os.path.join(ROOT_DIR, "runs")
FIGURES_DIR = os.path.join(ROOT_DIR, "figures")

RUNS = {
    "spar (draft=90 m)": ("shear_spar", 90.0, "tab:blue"),
    "semisub (draft=20 m)": ("shear_semisub", 20.0, "tab:orange"),
}

# Only animate the first N_DAYS days (covers both mixing timescales with
# margin -- semisub's t(phi*=0.5)=2.55 d is the slower of the two) and
# subsample frames for a reasonably small/fast GIF.
N_DAYS = 6.0
FRAME_STRIDE = 1
GIF_FPS = 5


def load_run(run_name, draft):
    cfg_path = os.path.join(RUNS_DIR, run_name, "resolved_config.yaml")
    ds, grid, params = open_roms_dataset(cfg_path)
    t_days = compute_time_vector(params)
    rho = ds["rho"].isel(eta_rho=0, xi_rho=0).values  # (time, s_rho)
    z = ds["z_rho"].isel(eta_rho=0, xi_rho=0).values  # (s_rho,) -- static column
    phi, _ = compute_phi(ds, grid, params)
    phi_star = (phi / phi.isel(ocean_time=0)).values
    ds.close()
    n = np.searchsorted(t_days, N_DAYS)
    idx = np.arange(0, n, FRAME_STRIDE)
    return t_days[idx], rho[idx], z, phi_star[idx], draft


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    data = {label: load_run(run_name, draft) for label, (run_name, draft, _) in RUNS.items()}
    n_frames = min(len(d[0]) for d in data.values())

    fig, (ax_prof, ax_race) = plt.subplots(1, 2, figsize=(10, 5))

    # --- left: density profile panel ---
    prof_lines = {}
    for label, (_, _, color) in RUNS.items():
        (line,) = ax_prof.plot([], [], color=color, lw=2, label=label)
        prof_lines[label] = line
        _, _, draft = RUNS[label][1], None, RUNS[label][1]
        ax_prof.axhspan(-draft, 0, color=color, alpha=0.08)
    all_rho = np.concatenate([d[1] for d in data.values()])
    ax_prof.set_xlim(all_rho.min() - 0.2, all_rho.max() + 0.2)
    ax_prof.set_ylim(-150, 0)
    ax_prof.set_xlabel(r"$\rho$ (density anomaly, kg m$^{-3}$)")
    ax_prof.set_ylabel("depth (m)")
    ax_prof.set_title("Density profile")
    ax_prof.legend(fontsize=8, loc="lower right")
    ax_prof.grid(alpha=0.3)
    time_text = ax_prof.text(0.03, 0.03, "", transform=ax_prof.transAxes, fontsize=11)

    # --- right: phi*(t) "race track" panel ---
    race_lines = {}
    race_markers = {}
    for label, (t_days, _, _, phi_star, _) in data.items():
        color = RUNS[label][2]
        ax_race.plot(t_days, phi_star, color=color, lw=1, alpha=0.3)
        (line,) = ax_race.plot([], [], color=color, lw=2.5, label=label)
        (marker,) = ax_race.plot([], [], "o", color=color, markersize=10)
        race_lines[label] = line
        race_markers[label] = marker
    ax_race.set_xlim(0, N_DAYS)
    ax_race.set_ylim(0, 1.02)
    ax_race.set_xlabel("time (days)")
    ax_race.set_ylabel(r"$\phi^*(t)$ (stratification remaining)")
    ax_race.set_title("Race to zero: who mixes out first?")
    ax_race.legend(fontsize=8)
    ax_race.grid(alpha=0.3)

    def update(frame):
        artists = []
        for label, (t_days, rho, z, phi_star, draft) in data.items():
            prof_lines[label].set_data(rho[frame], z)
            race_lines[label].set_data(t_days[: frame + 1], phi_star[: frame + 1])
            race_markers[label].set_data([t_days[frame]], [phi_star[frame]])
            artists += [prof_lines[label], race_lines[label], race_markers[label]]
        t0 = list(data.values())[0][0][frame]
        time_text.set_text(f"t = {t0:.2f} days")
        artists.append(time_text)
        return artists

    anim = animation.FuncAnimation(fig, update, frames=n_frames, interval=1000 / GIF_FPS, blit=False)
    out_path = os.path.join(FIGURES_DIR, "mixing_race.gif")
    anim.save(out_path, writer="pillow", fps=GIF_FPS)
    print(f"Saved animation to {out_path}")


if __name__ == "__main__":
    main()
