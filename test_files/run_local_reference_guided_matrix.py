#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def _run(cmd, cwd):
    return subprocess.run(cmd, cwd=str(cwd), check=False)


def _check_local_read_paths(job):
    missing = []
    pe = job.get("paired_end_lib") or {}
    se = job.get("single_end_lib") or {}
    for p in [pe.get("read1"), pe.get("read2"), se.get("read")]:
        if not p:
            continue
        pp = Path(os.path.expanduser(str(p)))
        if not pp.exists():
            missing.append(str(pp))
    return missing


def main():
    ap = argparse.ArgumentParser(description="Run local reference-guided matrix jobs.")
    ap.add_argument(
        "--matrix",
        default="test_files/local_reference_guided_matrix.json",
        help="Path to matrix JSON file.",
    )
    ap.add_argument(
        "--output-root",
        default="runs/local_matrix",
        help="Directory for per-case outputs.",
    )
    ap.add_argument(
        "--execute",
        action="store_true",
        help="Actually execute jobs. Without this flag, run preflight only.",
    )
    ap.add_argument(
        "--case-index",
        type=int,
        default=0,
        help="1-based matrix case index to run; 0 runs all.",
    )
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[1]
    matrix_path = (repo / args.matrix).resolve()
    output_root = (repo / args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    matrix = json.loads(matrix_path.read_text())
    if not isinstance(matrix, list):
        raise ValueError("Matrix JSON must be a list.")

    print(f"Loaded {len(matrix)} cases from {matrix_path}")
    any_fail = False

    for i, case in enumerate(matrix, start=1):
        if args.case_index and i != args.case_index:
            continue
        name = case["name"]
        job = case["job"]
        case_dir = output_root / f"{i:02d}_{name}"
        case_dir.mkdir(parents=True, exist_ok=True)
        job_file = case_dir / "job.json"
        job_file.write_text(json.dumps(job, indent=2))

        missing = _check_local_read_paths(job)
        if missing:
            any_fail = True
            print(f"[{i:02d}] {name}: missing local read files")
            for p in missing:
                print(f"  - {p}")
            continue

        print(f"[{i:02d}] {name}: preflight ok")
        if not args.execute:
            continue

        cmd = [
            sys.executable,
            str(repo / "scripts" / "run_viral_assembly.py"),
            "-j",
            str(job_file),
            "-o",
            str(case_dir),
        ]
        print(f"  running: {' '.join(cmd)}")
        rc = _run(cmd, repo)
        if rc.returncode != 0:
            any_fail = True
            print(f"  -> FAILED (exit {rc.returncode})")
        else:
            print("  -> OK")

    if any_fail:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

