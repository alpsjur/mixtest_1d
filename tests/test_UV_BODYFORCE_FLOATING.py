import numpy as np
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import open_roms_dataset, compute_z_r


def analytical_u(F, str_a, drag_coef, t):
    """Steady-drag-balanced body force solution (same ODE as
    test_STRUCTURE_DRAG.py's analytical_u), parametrized by a generic drag
    coefficient so it can be reused for both the structured zone (CD) and
    the floating/Gb zone (cb)."""
    alpha = 0.5 * drag_coef * str_a
    return np.sqrt(F / alpha) * np.tanh(t * np.sqrt(F * alpha))


def zone_levels(params):
    """
    Return (structured_levels, floating_levels): s_rho index arrays
    (0-based, matching the ROMS history file's s_rho dimension) for the
    structured (str_a>0) and floating (str_a==0, Gb-balanced) zones,
    mirroring tools/make_grd.py::build_str_a exactly.
    """
    H0 = params["grid"]["H0"]
    N = params["grid"]["N"]
    HC = params["vertical"]["HC"]
    THETA_S = params["vertical"]["THETA_S"]
    THETA_B = params["vertical"]["THETA_B"]
    depth_zero_below = params["structure"]["depth_zero_below"]

    h = np.array([[H0]])
    z_r = compute_z_r(h, HC, THETA_S, THETA_B, N)[:, 0, 0]
    dist_from_surface = -z_r
    structured_levels = np.where(dist_from_surface <= depth_zero_below)[0]
    floating_levels = np.where(dist_from_surface > depth_zero_below)[0]
    return structured_levels, floating_levels


def run_test(make_plots: bool = True) -> bool:
    """
    Validates the UV_BODYFORCE Gb balancing-drag term (floating structures)
    in a two-layer (structured/floating) water column:

    - The (majority) floating zone -- str_a==0, balanced by the new Gb term
      using bfrc_cb and a reference str_a taken from the nearest nonzero
      level above -- is compared tightly against the same closed-form
      tanh solution used in test_STRUCTURE_DRAG.py, since it dominates the
      column's volume and is only weakly perturbed by vertical mixing with
      the thin structured zone above it.

    - The (minority) structured zone -- str_a>0, balanced by the usual
      STRUCTURE_MIXING drag (str_cd) -- is checked only qualitatively:
      it must be finite/bounded, and (since CD > cb here) its mean
      velocity must be strictly lower than the floating zone's. An exact
      analytic match is not expected there because vertical turbulent
      mixing couples it to the much larger, faster floating zone below,
      which the simple single-layer power-balance ODE does not capture.
    """
    experimentpath = "runs/test_UV_BODYFORCE_FLOATING"
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
    assert cb < CD, (
        "This test assumes cb < CD (weaker drag in the floating zone) so "
        "that the floating zone reaches a higher, distinguishable steady "
        "velocity than the structured zone."
    )

    structured_levels, floating_levels = zone_levels(params)
    assert structured_levels.size > 0, "Expected at least one structured level."
    assert floating_levels.size > 0, "Expected at least one floating level."

    # Vertical turbulent mixing smooths the sharp str_a step over several
    # levels near the structured/floating interface, so only the deep
    # interior of the floating zone (well below the interface, away from
    # its influence) is expected to follow the single-layer analytical
    # solution tightly. Levels near the interface (both zones) are checked
    # only qualitatively below.
    interior_floating_levels = floating_levels[: floating_levels.size // 2]

    u_analytical_floating = analytical_u(F, str_a, cb, t)
    u_simulation = ds.u.values

    # --- Tight check: deep interior of the floating zone (majority,
    #     far from the mixing-influenced interface) -----------------------
    for k in interior_floating_levels:
        for i in range(ds.xi_u.size):
            for j in range(ds.eta_rho.size):
                assert np.isclose(
                    u_analytical_floating, u_simulation[:, k, j, i], rtol=3e-2
                ).all(), f"Floating-zone mismatch at xi={i}, eta={j}, s_rho={k}"
    print(
        "Floating-zone (Gb-balanced) u matches the single-layer analytical "
        "solution within 3% at all grid points."
    )

    # --- Qualitative checks: structured zone (minority, mixing-perturbed) -
    assert np.isfinite(u_simulation[:, structured_levels, :, :]).all(), (
        "Structured-zone velocity is not finite -- Gb/str_cd balance failed "
        "to bound the flow."
    )

    u_inf_naive_structured = np.sqrt(2.0 * F / (CD * str_a))
    u_inf_naive_floating = np.sqrt(2.0 * F / (cb * str_a))
    u_mean_structured = np.mean(u_simulation[-1, structured_levels, :, :])
    u_mean_floating = np.mean(u_simulation[-1, floating_levels, :, :])

    assert u_mean_structured < u_mean_floating, (
        "Expected the structured zone (CD) to reach a lower steady velocity "
        "than the floating zone (cb < CD), but found the opposite."
    )
    # Sanity bound: structured-zone mean should stay within the naive
    # decoupled range [CD-only estimate, cb-only estimate], since vertical
    # mixing can only pull it towards the (faster) floating zone below, not
    # push it beyond that.
    assert u_inf_naive_structured * 0.9 <= u_mean_structured <= u_inf_naive_floating * 1.05, (
        f"Structured-zone mean velocity {u_mean_structured:.4f} m/s is "
        f"outside the expected range "
        f"[{0.9*u_inf_naive_structured:.4f}, {1.05*u_inf_naive_floating:.4f}] m/s."
    )
    print(
        "Structured-zone u is finite, bounded, and consistently lower than "
        "the floating zone's, as expected from CD > cb."
    )

    if make_plots:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        th = t / (60 * 60)

        ax.plot(th, u_analytical_floating, label='Analytical, floating zone (Gb)',
                color='green', linewidth=2)

        u_sim_structured_mean = np.mean(u_simulation[:, structured_levels, :, :], axis=(1, 2, 3))
        u_sim_floating_mean = np.mean(u_simulation[:, floating_levels, :, :], axis=(1, 2, 3))
        ax.plot(th, u_sim_structured_mean, label='Numerical, structured zone,\ndomain mean',
                color='orange', linestyle='--', linewidth=2)
        ax.plot(th, u_sim_floating_mean, label='Numerical, floating zone,\ndomain mean',
                color='red', linestyle='--', linewidth=2)

        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('u (m/s)')
        ax.set_title(
            'UV_BODYFORCE floating-structure Gb term:\n'
            'analytical (floating zone) vs. numerical solutions'
        )
        ax.legend()
        ax.grid()

        figpath = os.path.join(ROOT_DIR, "figures", "test_UV_BODYFORCE_FLOATING.png")
        os.makedirs(os.path.dirname(figpath), exist_ok=True)
        fig.tight_layout()
        fig.savefig(figpath, dpi=150)
        plt.close(fig)

    return True


if __name__ == "__main__":
    ok = run_test(make_plots=True)
    sys.exit(0 if ok else 1)
