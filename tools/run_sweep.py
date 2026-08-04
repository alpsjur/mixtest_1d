#!/usr/bin/env python3
"""
tools/run_sweep.py

Run all ROMS simulations listed in a sweep manifest (manifest.yaml).

Iterates the runs in order, calls run_single_resolved for each, and updates
the manifest YAML + CSV after every run (so a crash mid-sweep doesn't lose progress).

Usage:
    python tools/run_sweep.py sweeps/<sweep_id>/manifest.yaml
"""

import os
import sys
from typing import List, Dict

# Ensure project root and tools/ are importable when called from tools/
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)
if THIS_DIR not in sys.path:
    sys.path.insert(0, THIS_DIR)

from utils.utils import load_yaml, save_yaml
from run_experiment import run_single_resolved
from prep_sweep import write_manifest_csv


def run_from_manifest(manifest_yaml_path: str) -> None:
    """
    Iterate runs in a sweep manifest and execute each simulation.
    Updates manifest YAML and CSV after each run for robustness.
    Per-run status is written by run_single_resolved to runs/<name>/logs/status.yaml.
    """
    manifest = load_yaml(manifest_yaml_path)
    runs = manifest.get("runs", [])
    sweep_dir = os.path.dirname(os.path.abspath(manifest_yaml_path))
    csv_path = os.path.join(sweep_dir, "manifest.csv")
    #         print(f"[{i}/{total}] Skipping (missing resolved_config): {run_name}")
    sweep_id = os.path.basename(sweep_dir)
    print(f"Running sweep from {sweep_id}")
    total = len(runs)
    for i, row in enumerate(runs, start=1):
        run_name = row.get("run_name", "<unnamed>")
        resolved_cfg_path = row.get("resolved_config")

        if not resolved_cfg_path or not os.path.isfile(resolved_cfg_path):
            print(f"[{i}/{total}] Skipping (missing resolved_config): {run_name}")
            row["status"] = "missing_config"
            save_yaml(manifest_yaml_path, manifest)
            write_manifest_csv(csv_path, runs)
            continue

        # Skip if already successfully done
        if row.get("status") == "done" and row.get("returncode", 1) == 0:
            print(f"[{i}/{total}] Already done: {run_name}")
            continue

        print(f"[{i}/{total}] Running: {run_name}")
        result = run_single_resolved(resolved_cfg_path)

        # Update manifest row
        row["status"] = result["state"]

        # Persist after each run (safer on interruption)
        save_yaml(manifest_yaml_path, manifest)
        write_manifest_csv(csv_path, runs)

    # Final write to ensure CSV/YAML are synced
    save_yaml(manifest_yaml_path, manifest)
    write_manifest_csv(csv_path, runs)
    print("All simulations processed from manifest.")


def main():
    if len(sys.argv) != 2:
        print("Usage: python tools/run_sweep.py sweeps/<sweep_id>/manifest.yaml", file=sys.stderr)
        sys.exit(1)
    run_from_manifest(sys.argv[1])


if __name__ == "__main__":
    main()