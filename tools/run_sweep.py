#!/usr/bin/env python3
"""
tools/run_sweep_sims.py

Run all ROMS simulations listed in a sweep manifest (manifest.yaml).
- Uses tools/run_single.py to execute each simulation.
- Updates per-run status (runs/<name>/logs/status.yaml).
- Updates the sweep manifest (YAML + CSV) with status, return code, log path, and status file.

Usage:
  python tools/run_sweep_sims.py sweeps/<sweep_id>/manifest.yaml
"""

import os
import sys
import yaml
from typing import List, Dict

# Ensure project root is importable when called from tools/
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from tools.run_experiment import run_single_resolved  # reuse the fixed-exec single-runner


def load_yaml(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def save_yaml(path: str, data: dict) -> None:
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def write_manifest_csv(csv_path: str, rows: List[Dict]) -> None:
    """
    Write a CSV manifest with updated status/return codes/logs.
    Columns align with prep_sweep output.
    """
    cols = ["hash_exact", "status", "run_name", "resolved_config", "params"]
    with open(csv_path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            params_str = yaml.safe_dump(r.get("params", {}), sort_keys=True).strip().replace("\n", "; ")
            line = [
                r.get("hash_exact", ""),
                r.get("status", ""),
                r.get("run_name", ""),
                r.get("resolved_config", ""),
                f"\"{params_str}\"",
            ]
            f.write(",".join(line) + "\n")


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
        print("Usage: python tools/run_sweep_sims.py sweeps/<sweep_id>/manifest.yaml", file=sys.stderr)
        sys.exit(1)
    run_from_manifest(sys.argv[1])


if __name__ == "__main__":
    main()