#!/usr/bin/env python3
# prep_hovmoller.py
import os
import sys

# Make project utils importable
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import open_roms_dataset, compute_time_vector


def prep_hovmoller(resolved_config_path: str, var: str):
    """
    Computes the area average of a variable over the horizontal 
    domain for each time step, returning a 2D array of shape (time, depth),
    and the corresponding time and depth vectors.

    Parameters
    ----------
    resolved_config_path : str
        Path to a resolved_config.yaml produced by prep_experiment.
    var : str
        Variable name to read from the history file.

    Returns
    -------
    days : ndarray
        Time axis in days (one entry per output record).
    depth : ndarray
            Depth axis corresponding to the second dimension of `mean_da`.
    mean_da : DataArray
        Area average of `var` over the horizontal domain for each time step.
    """
    ds, grid, params = open_roms_dataset(resolved_config_path)

    if var not in ds:
        ds.close()
        raise KeyError(f"Variable '{var}' not found in dataset.")

    mean_da = grid.average(ds[var], axis=("X", "Y")).squeeze()

    # Make sure var is on rho points 
    #mean_da = grid.interp(mean_da, "Z", to="center")

    days = compute_time_vector(params)

    if "s_rho" in mean_da.dims:
        depth = grid.average(ds.z_rho, axis=("X", "Y")).squeeze().values  # depth vector
    else:
        depth = grid.average(ds.z_w, axis=("X", "Y")).squeeze().values

    return days, depth, mean_da

def main():
    import argparse
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser(description="Plot a Hovmoller diagram from a resolved config.")
    parser.add_argument("--resolved_config", type=str, default="runs/baseline/resolved_config.yaml", help="Path to resolved_config.yaml")
    parser.add_argument("--variable", type=str, default="rho", help="Variable name to plot")
    parser.add_argument("--output", type=str, default=None, help="Output plot file (optional)")
    args = parser.parse_args()

    days, depth, mean_da = prep_hovmoller(args.resolved_config, args.variable)

    print(f"days shape: {days.shape}, depth shape: {depth.shape}, mean_da shape: {mean_da.shape}")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    c = ax.pcolormesh(days, depth, mean_da.T, shading='auto')
    fig.colorbar(c, ax=ax, label=args.variable)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(f"Hovmoller diagram of area-averaged {args.variable}")

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()


if __name__ == "__main__":
    main()