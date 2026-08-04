#!/usr/bin/env python3
"""
analysis/wrapper.py

Plot time series and vertical profiles for all runs in one or more sweep manifests.

For each variable and each analysis type (timeseries / profile), one figure is
produced that overlays all runs from both sweeps. The first sweep is drawn with
solid lines, the second with dashed lines, making it easy to visually compare
two parameter sweep families side by side.

Default sweeps and variables are defined in main() — edit them there to match
your current experiment.

Usage:
    python analysis/wrapper.py
    python analysis/wrapper.py --variables AKt temp --no-show --save-dir figures
"""

import os
import sys
import argparse
import matplotlib.pyplot as plt

THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from analysis.prep_timeseries import prep_timeseries
from analysis.prep_profiles import prep_profiles
from analysis.plot_sweep import plot_sweep


# ---------------------------------------------------------------------------
# Analysis type registry
# Each entry defines: the prep function, its extra kwargs, and axis labels.
# Add new analysis types here without changing the main loop.
# ---------------------------------------------------------------------------
ANALYSIS_TYPES = {
    "timeseries": {
        "prep_fn": prep_timeseries,
        "prep_kwargs": {},
        "xlabel": "Time (days)",
        "ylabel": "Variable (volume-avg)",
    },
    "profile": {
        "prep_fn": prep_profiles,
        "prep_kwargs": {"timestep": -1, "subtract_ini": True},
        "xlabel": "Variable (area-avg, last timestep, difference from initial)",
        "ylabel": "Depth (m)",
    },
}


def plot_all(sweep1: str, sweep2: str, variables: list, save_dir: str | None, show: bool):
    """
    Produce one figure per (variable, analysis_type) combination.

    Parameters
    ----------
    sweep1 : str
        Path to the first sweep's manifest.yaml (plotted with solid lines).
    sweep2 : str
        Path to the second sweep's manifest.yaml (plotted with dashed lines).
    variables : list of str
        Variable names to plot (e.g. ["AKt", "temp"]).
    save_dir : str or None
        Directory to save figures to. Figures are not saved if None.
    show : bool
        If True, call plt.show() at the end.
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    for variable in variables:
        for analysis, info in ANALYSIS_TYPES.items():
            print(f"Plotting {analysis} of {variable}...")
            fig, ax = plt.subplots(figsize=(8, 4.5))

            plot_sweep(sweep1, variable=variable, ax=ax,
                       prep_fn=info["prep_fn"], prep_kwargs=info["prep_kwargs"],
                       linestyle="-", linewidth=2)
            plot_sweep(sweep2, variable=variable, ax=ax,
                       prep_fn=info["prep_fn"], prep_kwargs=info["prep_kwargs"],
                       linestyle="--", linewidth=1)

            ax.set_xlabel(info["xlabel"])
            ax.set_ylabel(info["ylabel"])
            ax.set_title(f"{analysis.capitalize()} of {variable}")
            ax.grid(True, alpha=0.3)
            plt.tight_layout()

            if save_dir:
                fname = os.path.join(save_dir, f"{analysis}_{variable}.png")
                fig.savefig(fname)
                print(f"  Saved: {fname}")

    if show:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot sweep comparisons for selected variables.")
    parser.add_argument("--sweep1", default="sweeps/k-e_variations/manifest.yaml",
                        help="Path to first sweep manifest (solid lines)")
    parser.add_argument("--sweep2", default="sweeps/gen_variations/manifest.yaml",
                        help="Path to second sweep manifest (dashed lines)")
    parser.add_argument("--variables", nargs="+", default=["AKt", "temp"],
                        help="Variable names to plot (default: AKt temp)")
    parser.add_argument("--save-dir", default=None,
                        help="Save figures to this directory instead of (or in addition to) showing")
    parser.add_argument("--no-show", action="store_true",
                        help="Do not call plt.show() (useful for batch/CI runs)")
    args = parser.parse_args()

    plot_all(
        sweep1=args.sweep1,
        sweep2=args.sweep2,
        variables=args.variables,
        save_dir=args.save_dir,
        show=not args.no_show,
    )


if __name__ == "__main__":
    main()