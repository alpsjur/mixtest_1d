import numpy as np
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import open_roms_dataset, load_velocity_profile


def analytical_u(U0, bfrc_cd, t):
    """
    Exact solution of du/dt = bfrc_cd*(U0 - u), starting from rest
    (u(0)=0), with no ramp (BFRC_TSTR=BFRC_TEND=0). U0 may be an array
    (one target per depth/level); t is broadcast against it.
    """
    return U0 * (1.0 - np.exp(-bfrc_cd * t))


def run_test(make_plots: bool = True) -> bool:
    experimentpath = "runs/test_UV_BODYFORCE_PROFILE"
    ds, grid, params = open_roms_dataset(f"{experimentpath}/resolved_config.yaml")

    NTIMES = params["time_stepping"]["NTIMES"]
    DT = params["time_stepping"]["DT"]
    NHIS = params["time_stepping"]["NHIS"]
    bfrc_cd = params["bodyforce"]["BFRC_CD"]

    dt = DT * NHIS
    T = NTIMES * DT
    t = np.arange(0, T + dt / 2, dt)

    # Target steady-state velocity profile, U0(z), from the same text file
    # used to build the body force / initial condition, interpolated onto
    # each s_rho level's (domain-mean, since the grid is horizontally
    # uniform here) depth.
    profile_file = os.path.join(ROOT_DIR, params["bodyforce"]["profile_file"])
    depth, prof_u, _ = load_velocity_profile(profile_file)

    z_r = ds.z_rho.mean(dim=["eta_rho", "xi_rho"]).values  # (s_rho,), negative down
    dist = -z_r
    U0_k = np.interp(dist, depth, prof_u)  # (s_rho,), target per level

    u_simulation = ds.u.values  # (time, s_rho, eta_rho, xi_u)

    # Tolerance is looser here than in test_UV_BODYFORCE.py (uniform case):
    # with a depth-varying target, vertical (GLS) mixing continuously
    # smooths the prescribed shear, so the achieved profile approaches
    # -- but does not exactly equal -- the per-level target, especially
    # near the surface where this profile's curvature/shear is largest.
    for k in range(ds.s_rho.size):
        u_analytical_k = analytical_u(U0_k[k], bfrc_cd, t)
        for i in range(ds.xi_u.size):
            for j in range(ds.eta_rho.size):
                assert np.isclose(u_analytical_k, u_simulation[:, k, j, i], rtol=2e-2, atol=1e-4).all(), (
                    f"Mismatch at xi={i}, eta={j}, s_rho={k} (target U0={U0_k[k]:.4f} m/s)"
                )
    print(
        "Simulated depth-varying u matches the prescribed velocity profile "
        "for all grid points with relative tolerance 2e-2."
    )

    if make_plots:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        th = t / (60 * 60)

        u_simulation_mean = u_simulation.mean(axis=(2, 3))  # (time, s_rho)
        cmap = plt.get_cmap("viridis")
        for k in range(ds.s_rho.size):
            color = cmap(k / max(ds.s_rho.size - 1, 1))
            ax.plot(th, analytical_u(U0_k[k], bfrc_cd, t), color=color, linewidth=2)
            ax.plot(th, u_simulation_mean[:, k], color=color, linestyle='--', linewidth=1)

        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('u (m/s)')
        ax.set_title(
            'UV_BODYFORCE velocity-shear profile test\n'
            '(solid: analytical per level, dashed: simulated, colored by depth)'
        )
        ax.grid()
        fig.tight_layout()
        fig.savefig("figures/test_UV_BODYFORCE_PROFILE.png")
        plt.close(fig)

    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--no-plot", action="store_true", help="Disable plotting (useful in CI)."
    )
    args = parser.parse_args()

    try:
        ok = run_test(make_plots=not args.no_plot)
        sys.exit(0 if ok else 1)
    except AssertionError as e:
        print(str(e))
        sys.exit(1)
