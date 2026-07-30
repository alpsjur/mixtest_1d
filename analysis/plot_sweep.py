#!/usr/bin/env python3
# plot_sweep.py
import os
import sys
import matplotlib.pyplot as plt
from typing import Callable, Dict, Any, Tuple, Optional

# Ensure local and project imports work (adjust paths if needed)
THIS_DIR = os.path.dirname(__file__)
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from utils.utils import load_yaml
from analysis.prep_timeseries import prep_timeseries  # default prep function


def plot_sweep(
    manifest_path: str,
    variable: str = "AKt",
    ax: Optional[plt.Axes] = None,
    prep_fn: Callable[..., Tuple] = prep_timeseries,
    prep_kwargs: Optional[Dict[str, Any]] = None,
    label_runs: bool = True,
    linestyle: str = "-",
    linewidth: float = 2.0,
):
    """
    Plot data for all runs in a sweep manifest (YAML) onto the provided axis,
    using a user-supplied prep function.

    Parameters
    - manifest_path: path to manifest.yaml with a 'runs' list; each run needs 'resolved_config'
                     and optionally 'run_name'.
    - variable: variable name to pass to the prep function (default: 'AKt')
    - ax: existing matplotlib Axes; if None, a new one is created and returned
    - prep_fn: callable that prepares data for plotting. Must accept
               (resolved_config_path: str, variable: str, **kwargs) and return (x, y).
               Example: prep_timeseries -> (time_days, series_values)
                        prep_profile   -> (depth, profile_values)
    - prep_kwargs: optional dict of extra keyword args passed to prep_fn
    - label_runs: include a legend entry per run (default: True)
    - linestyle, linewidth: simple style controls applied to all runs in this sweep

    Returns
    - ax: matplotlib Axes with plotted lines
    """
    if prep_kwargs is None:
        prep_kwargs = {}

    manifest = load_yaml(manifest_path)
    runs = manifest.get("runs", [])
    runs = [r for r in runs if str(r.get("status", "done")).lower() == "done"]

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 4.5))

    for i, r in enumerate(runs):
        resolved_config = r["resolved_config"]
        label = r.get("run_name", f"run_{i+1}") if label_runs else None

        # Expect (x, y) from the prep function
        x, y = prep_fn(resolved_config, variable, **prep_kwargs)

        ax.plot(
            x,
            y,
            label=label,
            linestyle=linestyle,
            linewidth=linewidth,
        )

    if label_runs:
        ax.legend(loc="best", fontsize=9, frameon=True)

    return ax


# Example usage
if __name__ == "__main__":
    # Default: use timeseries prep
    fig, ax = plt.subplots(figsize=(8, 4.5))
    plot_sweep(
        "sweeps/k-e_variations/manifest.yaml",
        variable="AKt",
        ax=ax,
        prep_fn=prep_timeseries,      # timeseries: returns (time_days, series_values)
        linestyle="-",
        linewidth=2,
    )
    plot_sweep(
        "sweeps/gen_variations/manifest.yaml",
        variable="AKt",
        ax=ax,
        prep_fn=prep_timeseries,      # same prep, different sweep
        linestyle="--",
        linewidth=1,
    )
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("AKt (volume-avg)")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()