import numpy as np
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import open_roms_dataset


def analytical_u(F, str_a, Cd, t):
    alpha = 0.5 * Cd * str_a
    return np.sqrt(F / alpha) * np.tanh(t * np.sqrt(F * alpha))


def run_test(make_plots: bool = True) -> bool:
    """
    Control/continuity test: with bfrc_cb == CD, the Gb term (floating
    zone) and the Gd term (structured zone) apply the exact same drag law
    everywhere (same coefficient, same str_a value via the str_a_ref_omn
    reference), so the whole column -- despite str_a being zero below
    depth_zero_below -- must reproduce the single-layer analytical
    solution from test_STRUCTURE_DRAG.py at every level, with no shear.
    """
    experimentpath = "runs/test_UV_BODYFORCE_FLOATING_CB_EQ_CD"
    ds, grid, params = open_roms_dataset(f"{experimentpath}/resolved_config.yaml")

    NTIMES = params["time_stepping"]["NTIMES"]
    DT = params["time_stepping"]["DT"]
    NHIS = params["time_stepping"]["NHIS"]
    F = params["bodyforce"]["BFRC_U"] * 1e-7

    dt = DT * NHIS
    T = NTIMES * DT
    t = np.arange(0, T + dt / 2, dt)

    str_a = params["structure"]["str_a"]
    CD = params["structure"]["CD"]
    cb = params["structure"]["cb"]
    assert cb == CD, "This control test requires bfrc_cb == CD (no shear)."

    u_analytical = analytical_u(F, str_a, CD, t)
    u_simulation = ds.u.values

    for i in range(ds.xi_u.size):
        for j in range(ds.eta_rho.size):
            for k in range(ds.s_rho.size):
                assert np.isclose(
                    u_analytical, u_simulation[:, k, j, i], rtol=1e-4
                ).all(), f"Mismatch at xi={i}, eta={j}, s_rho={k}"
    print(
        "Simulated u equals the single-layer analytical solution at every "
        "grid point (relative tolerance 1e-4), confirming Gb reduces to Gd "
        "when bfrc_cb == CD."
    )

    if make_plots:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        th = t / (60 * 60)

        ax.plot(th, u_analytical, label='Analytical solution', color='blue', linewidth=2)

        u_simulation_mean = np.mean(u_simulation, axis=(1, 2, 3))
        ax.plot(
            th,
            u_simulation_mean,
            label='Numerical solution,\ndomain mean',
            color='orange',
            linestyle='--',
            linewidth=2
        )

        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('u (m/s)')
        ax.set_title(
            'Gb == Gd control case (bfrc_cb == CD):\n'
            'analytical vs. numerical solutions'
        )
        ax.legend()
        ax.grid()

        figpath = os.path.join(ROOT_DIR, "figures", "test_UV_BODYFORCE_FLOATING_CB_EQ_CD.png")
        os.makedirs(os.path.dirname(figpath), exist_ok=True)
        fig.tight_layout()
        fig.savefig(figpath, dpi=150)
        plt.close(fig)

    return True


if __name__ == "__main__":
    ok = run_test(make_plots=True)
    sys.exit(0 if ok else 1)
