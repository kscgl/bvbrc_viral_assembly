#!/usr/bin/env python3

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Tuple, List, Set, Optional


def _run(cmd: List[str], desc: str, cwd: Path) -> None:
    print(f"[Running] {desc}\n$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=str(cwd), check=True)
    print(f"[Done] {desc}\n")


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _normalize_fasta_bytes(b: bytes) -> bytes:
    """
    Normalize FASTA for comparison:
    - keep headers as-is (minus trailing whitespace)
    - remove all whitespace from sequence lines
    - ensure newline at end of each line
    """
    out: List[bytes] = []
    for raw in b.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(b">"):
            out.append(line + b"\n")
        else:
            out.append(b"".join(line.split()) + b"\n")
    return b"".join(out)


def _find_first(patterns: List[str], root: Path) -> Path:
    for pat in patterns:
        hits = sorted(root.glob(pat))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"Could not find any of {patterns} under {root}")


def _compare_fastas(old_fa: Path, new_fa: Path) -> Tuple[bool, Dict[str, str]]:
    old_b = _normalize_fasta_bytes(old_fa.read_bytes())
    new_b = _normalize_fasta_bytes(new_fa.read_bytes())
    same = old_b == new_b
    return same, {
        "old": str(old_fa),
        "new": str(new_fa),
        "old_sha256": hashlib.sha256(old_b).hexdigest(),
        "new_sha256": hashlib.sha256(new_b).hexdigest(),
    }


def _list_files(root: Path, patterns: List[str]) -> Dict[str, List[str]]:
    """
    Return sorted relative paths of files under root matching each glob pattern.
    """
    result: Dict[str, List[str]] = {}
    for pat in patterns:
        paths = sorted(str(p.relative_to(root)) for p in root.glob(pat))
        result[pat] = paths
    return result


def _compare_name_sets(
    old_root: Path,
    new_root: Path,
    patterns: List[str],
) -> Dict[str, Dict[str, List[str]]]:
    """
    For each pattern, report which relative paths are only in old, only in new, and in both.
    """
    summary: Dict[str, Dict[str, List[str]]] = {}
    for pat in patterns:
        old_paths: Set[str] = {str(p.relative_to(old_root)) for p in old_root.glob(pat)}
        new_paths: Set[str] = {str(p.relative_to(new_root)) for p in new_root.glob(pat)}
        summary[pat] = {
            "only_in_old": sorted(old_paths - new_paths),
            "only_in_new": sorted(new_paths - old_paths),
            "in_both": sorted(old_paths & new_paths),
        }
    return summary


def _brutal_compare(
    old_root: Path,
    new_root: Path,
    patterns: List[str],
) -> Dict[str, Dict]:
    """
    For each pattern:
      - compare the set of relative file paths
      - for files present in both, compare SHA-256 (raw bytes)
    """
    report: Dict[str, Dict] = {}
    for pat in patterns:
        old_paths: Set[str] = {str(p.relative_to(old_root)) for p in old_root.glob(pat)}
        new_paths: Set[str] = {str(p.relative_to(new_root)) for p in new_root.glob(pat)}
        only_old = sorted(old_paths - new_paths)
        only_new = sorted(new_paths - old_paths)
        both = sorted(old_paths & new_paths)

        mismatched: List[Dict[str, str]] = []
        matched: List[str] = []
        for rel in both:
            o = old_root / rel
            n = new_root / rel
            o_sha = _sha256_file(o)
            n_sha = _sha256_file(n)
            if o_sha != n_sha:
                mismatched.append(
                    {"path": rel, "old_sha256": o_sha, "new_sha256": n_sha}
                )
            else:
                matched.append(rel)

        report[pat] = {
            "only_in_old": only_old,
            "only_in_new": only_new,
            "matched_sha256": matched,
            "mismatched_sha256": mismatched,
        }
    return report


def _resolve_local_read_path(p: Optional[str], repo: Path, old_out_rel: str, new_out_rel: str) -> Optional[str]:
    """
    Resolve a read path from job JSON for local execution.

    Accepts:
    - absolute paths
    - paths relative to repo root
    - legacy prefixes like old_results/..., new_results/...
      (we remap these to the user-provided --old-out / --new-out locations)
    """
    if not p:
        return None

    # If it looks like a BV-BRC workspace reference, we can't fetch it locally here.
    if p.startswith("ws:"):
        raise ValueError(
            f"Job read path looks like a workspace path ({p}). "
            "For local runs, set paired_end_lib.read1/read2 (or single_end_lib.read) to local FASTQ paths."
        )

    cand = Path(p)
    if cand.is_absolute() and cand.exists():
        return str(cand)

    cand = (repo / p).resolve()
    if cand.exists():
        return str(cand)

    old_out_base = Path(old_out_rel)
    if not old_out_base.is_absolute():
        old_out_base = (repo / old_out_base).resolve()
    new_out_base = Path(new_out_rel)
    if not new_out_base.is_absolute():
        new_out_base = (repo / new_out_base).resolve()

    # Remap old_results/... and new_results/... prefixes to the provided output folders
    parts = Path(p).parts
    if parts:
        if parts[0] == "old_results":
            remap = (old_out_base / Path(*parts[1:])).resolve()
            if remap.exists():
                return str(remap)
        if parts[0] == "new_results":
            remap = (new_out_base / Path(*parts[1:])).resolve()
            if remap.exists():
                return str(remap)

    # Common case: JSON points at old_results/<file>, but caller used runs/old_results.
    # Try basename under old_out.
    remap = (old_out_base / Path(p).name).resolve()
    if remap.exists():
        return str(remap)

    raise FileNotFoundError(f"Could not resolve local read path: {p}")


def _run_new_locally(job: dict, repo: Path, new_out: Path, old_out_rel: str, new_out_rel: str) -> None:
    """
    Run the JSON-driven pipeline locally by calling scripts.reference_guided_assembly.run_reference_guided
    directly (avoids BV-BRC workspace tooling like p3-cp / p3-sra).
    """
    # Import from repo
    sys.path.insert(0, str(repo))
    from scripts.reference_guided_assembly import run_reference_guided  # type: ignore

    read1 = None
    read2 = None
    read_single = None

    pe = job.get("paired_end_lib") or {}
    se = job.get("single_end_lib") or {}

    if isinstance(pe, dict) and (pe.get("read1") or pe.get("read2")):
        read1 = _resolve_local_read_path(pe.get("read1"), repo, old_out_rel, new_out_rel)
        read2 = _resolve_local_read_path(pe.get("read2"), repo, old_out_rel, new_out_rel)
    elif isinstance(se, dict) and se.get("read"):
        read_single = _resolve_local_read_path(se.get("read"), repo, old_out_rel, new_out_rel)
    else:
        # Local-compare convenience:
        # If the job only specifies sra_id (backend mode), reuse FASTQs produced by the
        # original pipeline run under --old-out.
        srr = str(job.get("sra_id") or job.get("srr_id") or "")
        if srr:
            old_out_abs = Path(old_out_rel)
            if not old_out_abs.is_absolute():
                old_out_abs = (repo / old_out_abs).resolve()
            else:
                old_out_abs = old_out_abs.resolve()
            r1 = old_out_abs / f"{srr}_1.fastq"
            r2 = old_out_abs / f"{srr}_2.fastq"
            rs = old_out_abs / f"{srr}.fastq"
            if r1.exists() and r2.exists():
                read1, read2 = str(r1), str(r2)
            elif rs.exists():
                read_single = str(rs)

    print("[Running] JSON-driven reference-guided (local mode; no workspace fetch)")
    summary = run_reference_guided(
        job_data=job,
        output_dir=str(new_out),
        read1=read1,
        read2=read2,
        read_single=read_single,
    )
    print("[Done] JSON-driven reference-guided (local mode; no workspace fetch)\n")
    # Keep a breadcrumb for debugging
    try:
        (new_out / "run_reference_guided.summary.json").write_text(json.dumps(summary, indent=2))
    except Exception:
        pass


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run original + JSON-driven reference-guided pipeline and compare outputs."
    )
    ap.add_argument("--csv", dest="csv_path", help="CSV/TSV mapping file for original pipeline", required=True)
    ap.add_argument("--email", help="Entrez email for original pipeline (if accessions used)", default="")
    ap.add_argument("--job-json", dest="job_json", help="Job JSON for scripts/run_viral_assembly.py", required=True)
    ap.add_argument("--workdir", help="Repo root (default: script's parent)", default="")
    ap.add_argument("--old-out", dest="old_out", default="runs/old_results", help="Output dir for original pipeline")
    ap.add_argument("--new-out", dest="new_out", default="runs/new_results", help="Output dir for JSON-driven pipeline")
    ap.add_argument("--skip-run", action="store_true", help="Skip running pipelines; just compare outputs")
    args = ap.parse_args()

    repo = Path(args.workdir).resolve() if args.workdir else Path(__file__).resolve().parents[1]

    # Allow absolute output paths (useful on WSL to avoid /mnt/c I/O issues).
    old_out = Path(args.old_out)
    if not old_out.is_absolute():
        old_out = (repo / old_out).resolve()
    else:
        old_out = old_out.resolve()

    new_out = Path(args.new_out)
    if not new_out.is_absolute():
        new_out = (repo / new_out).resolve()
    else:
        new_out = new_out.resolve()

    csv_path = Path(args.csv_path).resolve()
    job_json = Path(args.job_json).resolve()

    if not args.skip_run:
        old_out.mkdir(parents=True, exist_ok=True)
        new_out.mkdir(parents=True, exist_ok=True)

        _run(
            [sys.executable, str(repo / "test_files" / "original_reference_guided_assembly.py"), "--map", str(csv_path), "--email", args.email, "--output", str(old_out)],
            "Original reference-guided (CSV)",
            cwd=repo,
        )
        # Run "new" pipeline locally (avoid BV-BRC workspace tooling)
        job = json.loads(job_json.read_text())

        # Speed up local parity: if the old run already produced a real
        # "<SRR>/<SRR>.sra" under --old-out, copy it into the new output so
        # run_reference_guided uses the same p3-sra staging path as production when sra_id is set.
        srr = str(job.get("sra_id") or job.get("srr_id") or "")
        if srr:
            old_sra = old_out / srr / f"{srr}.sra"
            new_sra = new_out / srr / f"{srr}.sra"
            if old_sra.exists() and not new_sra.exists():
                (new_out / srr).mkdir(parents=True, exist_ok=True)
                import shutil
                shutil.copy2(str(old_sra), str(new_sra))

        _run_new_locally(job, repo, new_out, str(old_out), str(new_out))

    # Load job to infer sample name and expected outputs
    job = json.loads(job_json.read_text())
    sample = str(job.get("output_file") or job.get("sra_id") or job.get("srr_id") or "sample")

    # Old consensus (original script writes either <sample>.consensus.multi.fasta or <sample>.consensus.<seg>.fasta etc.)
    old_cons = _find_first(
        [f"{sample}.consensus.multi.fasta", f"{sample}.consensus.*.fasta", f"{sample}.consensus.fasta"],
        old_out,
    )

    new_cons = _find_first(
        [f"{sample}.consensus.multi.fasta", f"{sample}.consensus.fasta"],
        new_out,
    )

    same, details = _compare_fastas(old_cons, new_cons)

    # Quick artifact checks
    new_vcf = new_out / f"{sample}.vcf.gz"
    old_vcf = old_out / f"{sample}.vcf.gz"
    artifacts = {
        "old_vcf_exists": str(old_vcf.exists()),
        "new_vcf_exists": str(new_vcf.exists()),
        "new_quast_html_exists": str((new_out / "quast_reference_guided" / "report.html").exists()),
        "new_quast_txt_exists": str((new_out / "quast_reference_guided" / "report.txt").exists()),
    }

    # File name / count comparisons for key artifact types
    patterns = [
        f"{sample}.consensus*.fasta",
        f"{sample}.vcf*",
        f"{sample}*.trim.fastq",
        f"{sample}.sam",
        f"{sample}.sorted.bam*",
        "quast_reference_guided/*",
    ]
    name_diffs = _compare_name_sets(old_out, new_out, patterns)
    brutal = _brutal_compare(old_out, new_out, patterns)

    print("=== Consensus comparison ===")
    print(json.dumps(details, indent=2))
    print(f"FASTA_match: {same}")
    print("=== Artifacts ===")
    print(json.dumps(artifacts, indent=2))
    print("=== File name comparison (old vs new) ===")
    print(json.dumps(name_diffs, indent=2))
    print("=== Brutal compare (SHA-256 for every common file) ===")
    print(json.dumps(brutal, indent=2))

    all_ok = same
    for pat, info in brutal.items():
        if info["only_in_old"] or info["only_in_new"] or info["mismatched_sha256"]:
            all_ok = False
            break

    return 0 if all_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())

