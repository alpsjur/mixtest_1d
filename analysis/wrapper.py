#!/usr/bin/env python3
# wrapper.py
import os
import sys
import matplotlib.pyplot as plt

THIS_DIR = os.path.dirname(__file__)
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from analysis.prep_timeseries import prep_timeseries 
from analysis.prep_profiles import prep_profiles   
from analysis.plot_sweep import plot_sweep


prep_dict = {
    "timeseries": {"prep_fn": prep_timeseries,
                   "prep_kwargs": {},
                   "xlabel": "Time (days)",
                   "ylabel": "Variable (volume-avg)",
                   },
    "profile": {"prep_fn": prep_profiles,
                 "prep_kwargs": {"timestep": -1, "subtract_ini": True},
                   "xlabel": "Variable (area-avg, last timestep, difference from initial)",
                   "ylabel": "Depth (m)",
                   },
}

for variable in ["AKt", "temp"]:  
    for analysis, prep_info in prep_dict.items():  
        print(f"Running {analysis} analysis on {variable}...")
        fig, ax = plt.subplots(figsize=(8, 4.5))
        plot_sweep(
            "sweeps/k-e_variations/manifest.yaml",
            variable=variable,
            ax=ax,
            prep_fn=prep_info["prep_fn"],     
            prep_kwargs=prep_info["prep_kwargs"],
            linestyle="-",
            linewidth=2,
        )
        plot_sweep(
            "sweeps/gen_variations/manifest.yaml",
            variable=variable,
            ax=ax,
            prep_fn=prep_info["prep_fn"],      
            prep_kwargs=prep_info["prep_kwargs"],
            linestyle="--",
            linewidth=1,
        )
        ax.set_xlabel(prep_info["xlabel"])
        ax.set_ylabel(prep_info["ylabel"])
        ax.set_title(f"{analysis.capitalize()} of {variable}")
        ax.grid(True, alpha=0.3)
        plt.tight_layout()


plt.show()