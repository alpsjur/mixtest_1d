#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

try:
    import yaml  # PyYAML
except ImportError:
    yaml = None

# Paths
TESTS_DIR = Path(__file__).resolve().parent
ROOT = TESTS_DIR.parent  # repository root

TOOLS_DIR = ROOT / "tools"
CONFIGS_DIR = ROOT / "configs"
VARIANTS_DIR = CONFIGS_DIR / "variants"
RUNS_DIR = ROOT / "runs"

PREP_SCRIPT = TOOLS_DIR / "prep_experiment.py"
RUN_SCRIPT = TOOLS_DIR / "run_experiment.py"
BASELINE_YAML = CONFIGS_DIR / "baseline.yaml"


def die(msg: str):
    print(msg, file=sys.stderr)
    sys.exit(2)


def check_layout():
    missing = []
    for p in (TOOLS_DIR, CONFIGS_DIR, VARIANTS_DIR, RUNS_DIR, PREP_SCRIPT, RUN_SCRIPT, BASELINE_YAML):
        if not p.exists():
            missing.append(str(p.relative_to(ROOT)))
    if missing:
        die("Missing required paths under repo root:\n  - " + "\n  - ".join(missing))


def discover_tests():
    """
    Return a sorted list of test names (filenames without .py)
    for files matching tests/test_*.py, excluding this runner.
    """
    names = []
    for p in TESTS_DIR.glob("test_*.py"):
        if p.name == "run_tests.py":
            continue
        names.append(p.stem)
    return sorted(names)


def run_cmd(cmd, cwd: Path = ROOT, env=None):
    proc = subprocess.run(
        cmd,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return proc.returncode, proc.stdout


def read_yaml(path: Path):
    if yaml is None:
        raise RuntimeError("PyYAML not installed. Install with: pip install pyyaml")
    with open(path, "r") as f:
        return yaml.safe_load(f)


def to_abs_under_root(path_str: str) -> Path:
    """
    If a path from config is relative, make it absolute under ROOT.
    If it's already absolute, return as-is.
    """
    p = Path(path_str)
    return p if p.is_absolute() else (ROOT / p)


def resolved_config_path(test_name: str) -> Path:
    return RUNS_DIR / test_name / "resolved_config.yaml"


def his_path_from_resolved_config(resolved_cfg_file: Path) -> Path:
    cfg = read_yaml(resolved_cfg_file)
    out_dir = to_abs_under_root(cfg["io"]["output_dir"])
    his_name = cfg["files"]["his"]
    return out_dir / his_name


def variant_yaml(test_name: str) -> Path:
    return VARIANTS_DIR / f"{test_name}.yaml"


def ensure_outputs(test_name: str, force_model: bool = False, verbose: bool = False):
    """
    Ensure the model outputs exist for the given test.
    - If force_model is True: always run prep + run.
    - Else: only run prep + run if history file is missing (or resolved config missing).
    Returns (ran: bool, ok: bool, log: str)
    """
    logs = []
    resolved_cfg = resolved_config_path(test_name)

    need_run = force_model or (not resolved_cfg.exists())
    his_exists = False
    if resolved_cfg.exists():
        try:
            his_path = his_path_from_resolved_config(resolved_cfg)
            his_exists = his_path.exists()
            if not force_model and not his_exists:
                need_run = True
        except Exception as e:
            logs.append(f"[{test_name}] Could not inspect his file from resolved_config.yaml: {e}")
            need_run = True

    if not need_run and his_exists:
        # Outputs present — nothing to do
        return (False, True, "\n".join(logs))

    # Prepare experiment
    prep_cmd = [
        sys.executable,
        str(PREP_SCRIPT),
        str(BASELINE_YAML),
        str(variant_yaml(test_name)),
    ]
    if verbose:
        logs.append(f"[{test_name}] Prep: {' '.join(prep_cmd)}")
    rc, out = run_cmd(prep_cmd)
    logs.append(out)
    if rc != 0:
        return (True, False, "\n".join(logs))

    # Run experiment
    rcfg = resolved_config_path(test_name)
    run_cmdline = [
        sys.executable,
        str(RUN_SCRIPT),
        str(rcfg),
    ]
    if verbose:
        logs.append(f"[{test_name}] Run: {' '.join(run_cmdline)}")
    rc, out = run_cmd(run_cmdline)
    logs.append(out)
    if rc != 0:
        return (True, False, "\n".join(logs))

    # Validate history file
    try:
        his_path = his_path_from_resolved_config(rcfg)
        if not his_path.exists():
            logs.append(f"[{test_name}] Model run finished but history file not found: {his_path}")
            return (True, False, "\n".join(logs))
    except Exception as e:
        logs.append(f"[{test_name}] Error validating model output: {e}")
        return (True, False, "\n".join(logs))

    return (True, True, "\n".join(logs))


def run_test_script(test_name: str):
    """
    Run a single test script from tests/ as a subprocess with CWD=repo root.
    Returns (rc, elapsed_seconds, output)
    """
    script_path = TESTS_DIR / f"{test_name}.py"
    if not script_path.exists():
        return (127, 0.0, f"[{test_name}] Test script not found at: {script_path}")

    env = os.environ.copy()
    # Ensure the repo root is importable for 'utils', etc.
    env["PYTHONPATH"] = str(ROOT) + os.pathsep + env.get("PYTHONPATH", "")

    start = time.time()
    rc, out = run_cmd([sys.executable, str(script_path)], env=env)
    end = time.time()
    return rc, end - start, out


def main():
    check_layout()

    parser = argparse.ArgumentParser(
        description=(
            "Run ROMS tests from tests/: "
            "default analyzes existing outputs (auto-run model if missing). "
            "Use --run-model to always run model before analysis."
        )
    )
    parser.add_argument(
        "--run-model",
        action="store_true",
        help="Always prepare and run experiments before analyzing outputs.",
    )
    parser.add_argument(
        "--tests",
        nargs="*",
        help="Subset of tests to run (filenames without .py). Default: auto-discovered tests/test_*.py.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show prep/run command lines and logs.",
    )
    args = parser.parse_args()

    all_tests = discover_tests()
    if not all_tests:
        die("No tests found in tests/ (expected files like tests/test_*.py).")

    if args.tests:
        selected = [t for t in args.tests if (TESTS_DIR / f"{t}.py").exists()]
        missing = set(args.tests) - set(selected)
        if missing:
            print("Warning: these tests were not found in tests/: " + ", ".join(sorted(missing)))
    else:
        selected = all_tests

    any_fail = False
    results = []

    print(f"Project root: {ROOT}")
    print("Selected tests:", ", ".join(selected))

    for test in selected:
        print(f"\n=== {test} ===")

        # Ensure outputs (force or lazy)
        ran, ok, model_log = ensure_outputs(test, force_model=args.run_model, verbose=args.verbose)
        if args.verbose and model_log:
            print(model_log)

        if not ok:
            print(f"[{test}] Model step FAILED.")
            results.append((test, "MODEL FAILED", 0.0))
            any_fail = True
            continue
        else:
            if ran:
                print(f"[{test}] Model step completed.")
            else:
                print(f"[{test}] Using existing outputs.")

        # Run the test analysis script
        rc, elapsed, out = run_test_script(test)
        print(out.strip())
        status = "PASS" if rc == 0 else f"FAIL (rc={rc})"
        print(f"[{test}] {status} in {elapsed:.1f}s")
        if rc != 0:
            any_fail = True
        results.append((test, "PASS" if rc == 0 else "FAIL", elapsed))

    # Summary
    print("\n=== Summary ===")
    for test, status, elapsed in results:
        suffix = "" if status == "MODEL FAILED" else f" ({elapsed:.1f}s)"
        print(f"- {test}: {status}{suffix}")

    sys.exit(1 if any_fail else 0)


if __name__ == "__main__":
    main()