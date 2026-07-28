#!/usr/bin/env python3
import os
import sys
import yaml
import datetime
import subprocess
import shlex

# Fixed ROMS executable path
ROMS_EXEC = "roms/romsS"


def load_yaml(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def save_yaml(path: str, data: dict) -> None:
    with open(path, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False)


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def write_run_status(logs_dir: str, status: dict) -> str:
    """Write status to runs/<run_name>/logs/status.yaml and return its path."""
    status_path = os.path.join(logs_dir, "status.yaml")
    save_yaml(status_path, status)
    return status_path


def run_single_resolved(resolved_cfg_path: str) -> dict:
    """
    Run a ROMS simulation from a resolved config.
    Returns: dict with run_dir, log, status_file, returncode, started_at, finished_at, state.
    """
    cfg = load_yaml(resolved_cfg_path)

    run_name = cfg["run"]["name"]
    run_dir = os.path.join("runs", run_name)
    input_dir = cfg["io"]["input_dir"]
    logs_dir = os.path.join(run_dir, "logs")
    ensure_dir(logs_dir)

    in_file = os.path.join(input_dir, "mixtest_1d.in")
    log_file = os.path.join(logs_dir, "simulation.log")

    # Resolve ROMS exec to an absolute path BEFORE changing cwd
    exec_path = os.path.abspath(ROMS_EXEC)

    # Preflight checks
    if not os.path.isfile(exec_path):
        raise FileNotFoundError(f"ROMS executable not found: {exec_path}")
    if not os.path.isfile(in_file):
        raise FileNotFoundError(f"Input file not found: {in_file}")

    # Build command string for status
    cmd = [exec_path]
    cmd_str = f"{shlex.join(cmd)} < {in_file} > {log_file} (cwd={run_dir})"

    started_at = datetime.datetime.now().isoformat(timespec="seconds")
    status = {
        "state": "running",
        "started_at": started_at,
        "cmd": cmd_str,
    }
    status_file = write_run_status(logs_dir, status)

    try:
        with open(in_file, "rb") as fin, open(log_file, "wb") as flog:
            proc = subprocess.run(
                cmd,
                stdin=fin,
                stdout=flog,
                stderr=subprocess.STDOUT,
                cwd=run_dir,
                check=False,
            )
        returncode = proc.returncode
    except FileNotFoundError as e:
        returncode = -1
        with open(log_file, "ab") as flog:
            flog.write(f"\nERROR: {e}\n".encode("utf-8"))

    finished_at = datetime.datetime.now().isoformat(timespec="seconds")
    state = "done" if returncode == 0 else "failed"
    status.update({"state": state, "finished_at": finished_at, "returncode": returncode})
    write_run_status(logs_dir, status)

    return {
        "run_dir": run_dir,
        "log": log_file,
        "status_file": status_file,
        "returncode": returncode,
        "started_at": started_at,
        "finished_at": finished_at,
        "state": state,
    }


def main():
    if len(sys.argv) != 2:
        print("Usage: python tools/run_experiment.py runs/<run_name>/resolved_config.yaml", file=sys.stderr)
        sys.exit(1)
    resolved_cfg_path = sys.argv[1]
    res = run_single_resolved(resolved_cfg_path)
    print(f"Run finished: state={res['state']} returncode={res['returncode']}")
    print(f"Log: {res['log']}")
    print(f"Status: {res['status_file']}")


if __name__ == "__main__":
    main()