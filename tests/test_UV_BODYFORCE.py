import xarray as xr
import numpy as np
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import prep_ds, load_yaml


def analytical_u(F, t):
    return F * t


def run_test(make_plots: bool = True) -> bool:
    experimentpath = "runs/test_UV_BODYFORCE"
    params = load_yaml(f"{experimentpath}/resolved_config.yaml")

    # Analytical solution for u
    NTIMES = params["time_stepping"]["NTIMES"]
    DT = params["time_stepping"]["DT"]
    NHIS = params["time_stepping"]["NHIS"]
    F = params["bodyforce"]["BFRC_U"] * 1e-7

    dt = DT * NHIS
    T = NTIMES * DT
    t = np.arange(0, T + dt / 2, dt)
    u_analytical = analytical_u(F, t)

    # Read the dataset and prepare it
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_dataset(
        params["io"]["output_dir"] + "/" + params["files"]["his"],
        decode_times=time_coder
    )
    ds, grid = prep_ds(ds, params)
    u_simulation = ds.u.values

    # Check that the analytical solution matches the numerical solution
    # for every xi, eta and s_rho point
    for i in range(ds.xi_u.size):
        for j in range(ds.eta_rho.size):
            for k in range(ds.s_rho.size):
                assert np.isclose(u_analytical, u_simulation[:, k, j, i], rtol=1e-5).all(), (
                    f"Mismatch at xi={i}, eta={j}, s_rho={k}"
                )
    print("Simulated u is equal to analytical u for all grid points with relative tolerance 1e-5.")

    if make_plots:
        import matplotlib.pyplot as plt

        # Also plot the analytical and numerical solutions for visual inspection
        fig, ax = plt.subplots(figsize=(6, 4))
        th = t / (60 * 60)

        # Plot the analytical solution
        ax.plot(th, u_analytical, label='Analytical solution', color='blue', linewidth=2)

        # Plot the numerical solution with spread
        u_simulation_mean = np.mean(u_simulation, axis=(1, 2, 3))
        ax.plot(
            th,
            u_simulation_mean,
            label='Numerical solution,\ndomain mean',
            color='orange',
            linestyle='--',
            linewidth=2
        )
        u_simulation_max = np.max(u_simulation, axis=(1, 2, 3))
        u_simulation_min = np.min(u_simulation, axis=(1, 2, 3))
        ax.fill_between(
            th,
            u_simulation_min,
            u_simulation_max,
            color='orange',
            alpha=0.2,
            label='Numerical solution,\nspread'
        )

        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('u (m/s)')
        ax.set_title('Comparison of analytical and numerical\nsolutions for u given a constant body force')
        ax.legend()
        ax.grid()

        # add analytical expression for u in the plot
        ax.text(
            0.5, 0.1, r'$u(t) = F t$',
            transform=ax.transAxes, fontsize=12, verticalalignment='bottom',
            horizontalalignment='center'
        )
        fig.tight_layout()
        fig.savefig("figures/test_UV_BODYFORCE.png")
        plt.close(fig)

    return True


if __name__ == "__main__":
    import argparse
    import sys

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