#!/usr/bin/env python3
# prep_timeseries.py
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Make project utils importable
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import open_roms_dataset, compute_time_vector


def prep_timeseries(resolved_config_path: str, var: str):
    """
    Compute the volume-averaged time series of a variable.

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
    mean_da : DataArray
        Volume average of `var` over time.
    """
    ds, grid, params = open_roms_dataset(resolved_config_path)

    if var not in ds:
        ds.close()
        raise KeyError(f"Variable '{var}' not found in dataset.")

    mean_da = grid.average(ds[var], axis=("X", "Y", "Z")).squeeze()
    days = compute_time_vector(params)

    return days, mean_da

def main():
    parser = argparse.ArgumentParser(description="Plot a single timeseries from a resolved config.")
    parser.add_argument("--resolved_config", type=str, default="runs/baseline/resolved_config.yaml", help="Path to resolved_config.yaml")
    parser.add_argument("--variable", type=str, default="AKt", help="Variable name to plot")
    parser.add_argument("--output", type=str, default=None, help="Output plot file (optional)")
    args = parser.parse_args()

    days, mean_da = prep_timeseries(args.resolved_config, args.variable)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(days, mean_da.values, label=args.variable)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel(f"{args.variable} (volume-avg)")
    ax.set_title(f"Time series of domain averaged {args.variable}")
    ax.grid()

    if args.output:
        plt.savefig(args.output)
        print(f"Plot saved to {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()