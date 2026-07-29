"""
The idea if this test is to see if the production term
in STRUCTURE_MIXING balance dissipation in steady state,
when no other production term is present (i.e. no shear, 
no stratification).
"""

import xarray as xr
import numpy as np
import os
import sys

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import prep_ds, load_yaml


def calc_Pd(str_a, Cd, u):
    return 0.5 * Cd * str_a * np.abs(u ** 3)


def run_test(make_plots: bool = True) -> bool:
    experimentpath = "runs/test_STRUCTURE_PRODUCTION"
    params = load_yaml(f"{experimentpath}/resolved_config.yaml")

    str_a = params["structure"]["str_a"]
    Cd = params["structure"]["CD"]

    # Read the dataset and prepare it
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_dataset(
        params["io"]["output_dir"] + "/" + params["files"]["his"],
        decode_times=time_coder
    )
    ds, grid = prep_ds(ds, params)

    # Shift variables to cell centers for easier comparison
    u = grid.interp(ds.u, "X")
    gls = grid.interp(ds.gls, "Z")
    Pd = calc_Pd(str_a, Cd, u)

    GLS = grid.integrate(gls, ("X", "Y", "Z"))
    PD = grid.integrate(Pd, ("X", "Y", "Z"))

    # Check that the domain integrated production and dissipation terms are equal
    # in steady state up to 2 percent tolerance
    diff = np.abs(GLS[-1] - PD[-1])
    accuracy = diff / np.abs(PD[-1])
    assert accuracy < 2e-2, (
        f"Mismatch between domain integrated production and dissipation terms: "
        f"{accuracy * 100:.2f}%"
    )
    print(
        f"Domain integrated production and dissipation terms are equal in steady state "
        f"with {accuracy * 100:.2f}% tolerance."
    )

    if make_plots:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(PD, label=r'Analytical $P_d$', color='blue')
        ax.plot(GLS, label=r'Numerical $\epsilon$', color='orange', linestyle='--')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Domain integrated values')
        ax.legend()
        ax.set_title('Comparing domain integrated structure production\nand dissipation terms for shear free, unstratified flow')
        ax.grid()
        fig.tight_layout()
        fig.savefig("figures/test_STRUCTURE_PRODUCTION.png")
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