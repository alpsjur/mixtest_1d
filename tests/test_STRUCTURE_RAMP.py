import numpy as np
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import open_roms_dataset


def ramp_factor(t, T):
    return np.where(t >= T, 1.0, 0.5 * (1.0 - np.cos(np.pi * t / T)))


def integrate_u(F, alpha, T_ramp, t, nsub=200):
    """Integrate du/dt = F - fac(t)*alpha*u^2 with RK4 on the output times t."""
    def rhs(tt, u):
        return F - float(ramp_factor(tt, T_ramp)) * alpha * u * u

    u = np.zeros_like(t)
    for n in range(1, t.size):
        h = (t[n] - t[n - 1]) / nsub
        tt, un = t[n - 1], u[n - 1]
        for _ in range(nsub):
            k1 = rhs(tt, un)
            k2 = rhs(tt + h / 2, un + h / 2 * k1)
            k3 = rhs(tt + h / 2, un + h / 2 * k2)
            k4 = rhs(tt + h, un + h * k3)
            un += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            tt += h
        u[n] = un
    return u


def run_test(make_plots: bool = True) -> bool:
    experimentpath = "runs/test_STRUCTURE_RAMP"
    ds, grid, params = open_roms_dataset(f"{experimentpath}/resolved_config.yaml")

    NTIMES = params["time_stepping"]["NTIMES"]
    DT = params["time_stepping"]["DT"]
    NHIS = params["time_stepping"]["NHIS"]
    F = params["bodyforce"]["F_U"]
    str_a = params["structure"]["str_a"]
    Cd = params["structure"]["CD"]
    ramp = params["structure"]["ramp"]
    assert ramp > 0, "This test requires structure.ramp > 0"

    alpha = 0.5 * Cd * str_a
    T_ramp = ramp * DT
    dt = DT * NHIS
    t = np.arange(0, NTIMES * DT + dt / 2, dt)

    u_ramp = integrate_u(F, alpha, T_ramp, t)
    u_noramp = np.sqrt(F / alpha) * np.tanh(t * np.sqrt(F * alpha))
    u_sim = ds.u.values

    # 1. Simulation follows the ramped ODE (drag follows the cosine factor)
    for i in range(ds.xi_u.size):
        for j in range(ds.eta_rho.size):
            for k in range(ds.s_rho.size):
                assert np.isclose(
                    u_ramp, u_sim[:, k, j, i], rtol=1e-3, atol=1e-6
                ).all(), f"Mismatch with ramped solution at xi={i}, eta={j}, s_rho={k}"
    print("Simulated u matches the ramped analytical/ODE solution (rtol 1e-3).")

    # 2. Ramp slows startup: less drag early means u is larger than without ramp,
    #    and the difference vanishes once the ramp is complete and flow equilibrates.
    u_mean = u_sim.mean(axis=(1, 2, 3))
    after = t > T_ramp
    assert (u_mean[after] >= u_noramp[after] - 1e-9).all(), \
        "Ramped run should not have less velocity than the no-ramp run"
    gap = u_mean - u_noramp
    assert gap[-1] < gap.max(), "Ramp influence should decay after the ramp ends"
    print("Ramp behaviour consistent with expectation.")

    if make_plots:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        th = t / 3600
        ax.plot(th, u_noramp, label="Analytical, no ramp", color="gray")
        ax.plot(th, u_ramp, label="ODE with cosine ramp", color="blue", linewidth=2)
        ax.plot(th, u_mean, label="Numerical, domain mean", color="orange",
                linestyle="--", linewidth=2)
        ax.axvline(T_ramp / 3600, color="k", linestyle=":", label="End of ramp")
        ax.set_xlabel("Time (hours)")
        ax.set_ylabel("u (m/s)")
        ax.set_title("Structure drag with cosine ramp-up")
        ax.legend()
        ax.grid()
        fig.tight_layout()
        os.makedirs("figures", exist_ok=True)
        fig.savefig("figures/test_STRUCTURE_RAMP.png")
        plt.close(fig)

    return True


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-plot", action="store_true",
                        help="Disable plotting (useful in CI).")
    args = parser.parse_args()
    try:
        sys.exit(0 if run_test(make_plots=not args.no_plot) else 1)
    except AssertionError as e:
        print(str(e))
        sys.exit(1)
