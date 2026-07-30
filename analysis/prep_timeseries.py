#!/usr/bin/env python3
# prep_timeseries.py
import os
import sys
import argparse
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

# Make project utils importable (adjust paths if needed)
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import prep_ds, load_yaml


def compute_time_vector(params: dict) -> np.ndarray:
    DT = params["time_stepping"]["DT"]
    NHIS = params["time_stepping"]["NHIS"]
    NTIMES = params["time_stepping"]["NTIMES"]
    N = NTIMES // NHIS
    dt = DT * NHIS  # seconds per record
    t_seconds = np.arange(N+1, dtype=np.float64) * dt

    return t_seconds / (60 * 60 * 24)  # days



def prep_timeseries(resolved_config_path: str, var: str):

    # Load configuration
    params = load_yaml(resolved_config_path)
    out_dir = params["io"]["output_dir"]
    his_file = params["files"]["his"]
    ds_path = os.path.join(out_dir, his_file)

    # Open dataset; use_cftime to avoid warnings for non-proleptic Gregorian calendars
    time_coder = xr.coders.CFDatetimeCoder(use_cftime=True)
    ds = xr.open_dataset(ds_path, decode_times=time_coder)

    # Prepare dataset and grid (project-specific)
    ds, grid = prep_ds(ds, params)

    if var not in ds:
        ds.close()
        raise KeyError(f"Variable '{var}' not found in dataset.")

    da = ds[var]

    # Spatial average
    avg_dims = ("X", "Y", "Z")  # Default average dimensions
    mean_da = grid.average(da, axis=avg_dims).squeeze()

    # Time axis
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