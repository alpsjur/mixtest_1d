#!/usr/bin/env python3
"""
tools/run_sweep.py

Run a cartesian sweep over parameters, creating one run per combination.

- Loads a sweep YAML (single positional arg)
- For each combination:
  - resolved_cfg = deep_merge(base, overrides_for_combo)
  - set resolved_cfg['run']['name'] from run_name_template
  - delegate per-run setup to setup_experiment.prepare_run_from_resolved(resolved_cfg)
  - record manifest entries (YAML + CSV) under sweep.output_dir

Usage:
  python tools/prep_sweep.py sweeps/my_sweep.yaml
"""

import os
import re
import sys
import yaml
import itertools
from copy import deepcopy
import shutil
import filecmp

# Ensure project root is importable when running from tools/
THIS_DIR = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(THIS_DIR, ".."))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)

from prep_experiment import prepare_run_from_resolved, deep_merge


def set_by_dotted_key(d: dict, dotted: str, value):
    """Set a nested dictionary value by dotted path, creating dicts as needed."""
    parts = dotted.split(".")
    cur = d
    for p in parts[:-1]:
        if p not in cur or not isinstance(cur[p], dict):
            cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value


def get_by_dotted_key(d: dict, dotted: str):
    """Get a nested dictionary value by dotted path."""
    cur = d
    for p in dotted.split("."):
        cur = cur[p]
    return cur


def cartesian_dict(param_dict: dict):
    """Yield a dict for each cartesian combination of parameter lists."""
    keys = list(param_dict.keys())
    for combo in itertools.product(*(param_dict[k] for k in keys)):
        yield dict(zip(keys, combo))


def format_run_name(template: str, cfg: dict, params: dict, index: int) -> str:
    """Fill {dotted.path} placeholders from cfg (or params), and {index}."""
    def repl(m):
        key = m.group(1).strip()
        if key == "index":
            return str(index)
        try:
            val = get_by_dotted_key(cfg, key)
        except Exception:
            val = params.get(key, "")
        s = str(val)
        return s.replace("/", "-").replace(" ", "")
    return re.sub(r"\{([^{}]+)\}", repl, template)


def ensure_dir(path: str) -> None:
    """Create directory if missing."""
    os.makedirs(path, exist_ok=True)


def write_manifest_yaml(path: str, sweep_id: str, rows: list) -> None:
    """Write YAML manifest for the sweep."""
    with open(path, "w") as f:
        yaml.safe_dump({"sweep_id": sweep_id, "runs": rows}, f, sort_keys=False)


def write_manifest_csv(path: str, rows: list) -> None:
    """Write CSV manifest for quick analysis."""
    cols = ["hash_exact", "status", "run_name", "resolved_config", "params"]
    with open(path, "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            params_str = yaml.safe_dump(r["params"], sort_keys=True).strip().replace("\n", "; ")
            line = [
                r["hash_exact"],
                r["status"],
                r["run_name"],
                r["resolved_config"],
                f"\"{params_str}\"",
            ]
            f.write(",".join(line) + "\n")


def main():
    """Load the sweep definition, iterate combinations, delegate per-run prep, write manifests."""
    sweep_path = sys.argv[1] if len(sys.argv) > 1 else "sweeps/sweep.yaml"
    with open(sweep_path, "r") as f:
        sweep = yaml.safe_load(f)

    sdef = sweep["sweep"]
    sweep_id = sdef["id"]
    base_config_path = sdef["base_config"]
    apply_config_path = sdef.get("apply_config")
    params = sdef["parameters"]
    run_name_template = sdef["run_name_template"]
    sweep_out_dir = sdef["output_dir"]

    ensure_dir(sweep_out_dir)

    # Copy the sweep YAML into the sweep output directory (overwrite if exists)
    sweep_yaml_src = os.path.abspath(sweep_path)
    sweep_yaml_dst = os.path.join(sweep_out_dir, "sweep.yaml")

    # Check if the destination file exists and is identical to the source
    if os.path.exists(sweep_yaml_dst) and filecmp.cmp(sweep_yaml_src, sweep_yaml_dst, shallow=False):
        print(f"File already exists and is identical. No need to copy: {sweep_yaml_dst}")
    else:
        shutil.copy2(sweep_yaml_src, sweep_yaml_dst)
        print(f"Copied sweep definition to: {sweep_yaml_dst}")

    with open(base_config_path, "r") as f:
        base = yaml.safe_load(f) or {}

    if apply_config_path:
        with open(apply_config_path, "r") as f:
            applied = yaml.safe_load(f) or {}
        base = deep_merge(base, applied)

    manifest_rows = []
    index = 1
    for combo in cartesian_dict(params):
        # Build override for this combination and resolve config
        over = {}
        for k, v in combo.items():
            set_by_dotted_key(over, k, v)
        cfg = deep_merge(base, over)

        # Construct run name and set into cfg
        run_name = format_run_name(run_name_template, cfg, combo, index)
        cfg.setdefault("run", {})
        cfg["run"]["name"] = run_name

        # Delegate per-run preparation (hashing, dirs, grid/ini/bry, .in)
        result = prepare_run_from_resolved(cfg)

        # Record manifest entry (use hash from setup_experiment)
        manifest_rows.append({
            "hash_exact": result["hash_exact"],
            "status": "prepared",
            "run_name": run_name,
            "resolved_config": result["resolved_config"],
            "params": combo,
        })

        print(f"Prepared run {index}: {run_name} -> {result['run_dir']}")
        index += 1

    # Write manifests
    write_manifest_yaml(os.path.join(sweep_out_dir, "manifest.yaml"), sweep_id, manifest_rows)
    write_manifest_csv(os.path.join(sweep_out_dir, "manifest.csv"), manifest_rows)




if __name__ == "__main__":
    main()