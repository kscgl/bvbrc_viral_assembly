#!/usr/bin/env python3
"""
Viral Assembly Pipeline (CSV-driven; supports segmented **and** non-segmented genomes)

Synopsis:
  Command-line pipeline that builds consensus genomes from short reads using a CSV mapping file.
  Supports local FASTA(s) or GenBank accessions, paired- or single-end reads (local or SRA), and
  segmented or non-segmented references.

Input:
  - Mapping file: CSV/TSV specified with `--map` (header required).
  - Required CSV columns:
      `SRA` or `sra_id`          : sample identifier or SRA accession (used as sample name)
      `ref_fasta` or `reference` : single reference FASTA path or accession
      `ref_fastas` or `references` : semicolon-separated list of FASTA paths/accessions (for segmented genomes)
  - Optional CSV columns:
      `segment_names`            : semicolon-separated contig names to use (when supplying a multi-FASTA)
      `region`                   : region/contig to restrict variant calling (passed to bcftools)
      Read file columns (case-insensitive):
        Paired-end: `fastq1` / `fastq2` (many name variants accepted)
        Single-end: `fastq` (or `fastq_single`)
      If read paths are provided, the pipeline copies them into the sample output directory.
      If no local reads provided, the pipeline attempts to download from SRA (`prefetch` + `fasterq-dump`).

Behavior / assumptions:
  - Reference tokens may be local FASTA paths or GenBank accessions. Accessions are resolved via Entrez.
  - If a single multi-record FASTA is provided and `segment_names` are given, headers may be rewritten.
  - Paired vs single-end is auto-detected by presence of both `fastq1` and `fastq2` (or a single `fastq`).
  - Tools required on PATH: `bwa`, `samtools`, `bcftools`, `fastp`, `prefetch`, `fasterq-dump`, `plot-vcfstats`.

Outputs (per-sample, written to `--output` directory):
  - Copied/downloaded FASTQ(s): <SAMPLE>_1.fastq, <SAMPLE>_2.fastq or <SAMPLE>.fastq
  - Trimmed FASTQ(s): <SAMPLE>.trim.fastq or <SAMPLE>_1.trim.fastq, <SAMPLE>_2.trim.fastq
  - Alignment: <SAMPLE>.sam, <SAMPLE>.sorted.bam, index
  - Variant calls: <SAMPLE>.vcf.gz (+ index)
  - Consensus sequences: <SAMPLE>.consensus.multi.fasta and per-segment FASTA(s)
  - Low-coverage mask BED (if `--depth-cutoff` > 0): <SAMPLE>.lowcov_dp<cutoff>.bed
  - Optional stats and plots: <SAMPLE>.vcf.stats, <SAMPLE>_plots/
  - Summary table: results/variant_counts.tsv (contains sample, variant counts, file paths)

Invocation:
  python viral_assembly_segmented_v4.py --map `samples.csv` --email `you@org.org` [--output `results`] [--workers N]
  Common options: `--depth-cutoff`, `--stats-plots`, `--align-threads`, `--fastp-threads`, `--version`

Notes:
  - Use exact contig/segment names matching the reference when supplying `segment_names` or when using `-r` with `bcftools consensus`.
  - Paths in the CSV may be absolute or relative; tilde (`~`) is expanded.
"""

import os
import sys
import csv
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import pandas as pd
import re

from Bio import Entrez, SeqIO

# ----------------------- utilities -----------------------

def run_command(cmd: str, desc: str, *, capture_stdout: bool = False):
    print(f"[Running] {desc}\n$ {cmd}")
    try:
        if capture_stdout:
            out = subprocess.check_output(cmd, shell=True)
            print(f"[Done] {desc}\n")
            return out
        else:
            subprocess.check_call(cmd, shell=True)
            print(f"[Done] {desc}\n")
            return None
    except subprocess.CalledProcessError as e:
        print(f"[Failed] {desc}\nExit code: {e.returncode}")
        sys.exit(1)

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

# ----------------------- FASTA helpers -----------------------

def list_contigs_from_fasta(fa_path: Path) -> List[str]:
    contigs = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.append(line[1:].strip().split()[0])
    return contigs

def build_multifasta_from_paths(fastas: List[Path], out_fa: Path, segment_names: Optional[List[str]] = None):
    """Concatenate multiple FASTA files into out_fa, optionally overriding headers with segment_names."""
    with open(out_fa, "w") as out:
        seg_idx = 0
        for fa in fastas:
            with open(fa) as fh:
                for line in fh:
                    if line.startswith(">"):
                        if segment_names is not None:
                            if seg_idx >= len(segment_names):
                                raise ValueError("segment_names shorter than number of records across input FASTAs.")
                            out.write(f">{segment_names[seg_idx]}\n")
                            seg_idx += 1
                        else:
                            out.write(line)
                    else:
                        out.write(line)
        if segment_names is not None and seg_idx != len(segment_names):
            raise ValueError("segment_names length does not match total number of records across input FASTAs.")

def rewrite_headers_for_multifasta(src_multi: Path, out_fa: Path, new_names: List[str]):
    """Rewrite headers of a single multi-record FASTA to the provided names (1:1)."""
    # Count records first
    headers = list_contigs_from_fasta(src_multi)
    if len(headers) != len(new_names):
        raise ValueError(
            f"segment_names count ({len(new_names)}) does not match number of records in {src_multi} ({len(headers)})."
        )
    with open(src_multi) as inp, open(out_fa, "w") as out:
        idx = 0
        for line in inp:
            if line.startswith(">"):
                out.write(f">{new_names[idx]}\n")
                idx += 1
            else:
                out.write(line)

def bwa_index_present(ref_fa: Path) -> bool:
    suffixes = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    return all((ref_fa.with_suffix(ref_fa.suffix + s)).exists() for s in suffixes)

def write_lowcov_bed_from_depth(depth_bytes: bytes, cutoff: int, bed_path: Path):
    """Convert `samtools depth -a` output (bytes) to merged BED regions where depth < cutoff."""
    with open(bed_path, "w") as bed:
        last_ctg = None
        start = None
        last_pos = None
        for line in depth_bytes.decode().splitlines():
            if not line:
                continue
            ctg, pos_s, dp_s = line.split()[:3]
            pos = int(pos_s)
            dp = int(dp_s)
            if dp < cutoff:
                if ctg == last_ctg and last_pos is not None and pos == last_pos + 1:
                    last_pos = pos
                else:
                    if last_ctg is not None and start is not None:
                        bed.write(f"{last_ctg}\t{start-1}\t{last_pos}\n")
                    last_ctg = ctg
                    start = pos
                    last_pos = pos
            else:
                if last_ctg is not None and start is not None:
                    bed.write(f"{last_ctg}\t{start-1}\t{last_pos}\n")
                last_ctg = None
                start = None
                last_pos = None
        if last_ctg is not None and start is not None:
            bed.write(f"{last_ctg}\t{start-1}\t{last_pos}\n")

# ----------------------- Entrez (GenBank) helpers -----------------------

def is_probable_accession(token: str) -> bool:
    """Heuristic: if token doesn't look like a local path and contains no path separators, treat as accession."""
    p = Path(os.path.expanduser(token))
    if p.exists():
        return False
    # Basic guard: accessions usually alnum + underscores/dots, no slashes
    return ("/" not in token) and ("\\" not in token)

def fetch_genbank_to_fasta(accession: str, email: str, out_fa: Path) -> Path:
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    # Some accessions may contain multiple records (rare). Use SeqIO.parse then write all out.
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    if not records:
        raise ValueError(f"No FASTA sequence returned for accession {accession}")
    with open(out_fa, "w") as fh:
        SeqIO.write(records, fh, "fasta")
    return out_fa

def resolve_refs_to_fastas(
    ref_tokens: List[str],
    cache_dir: Path,
    email: str
) -> List[Path]:
    """Resolve a list of tokens that are either local paths or GenBank accessions into local FASTA paths."""
    ensure_dir(cache_dir)
    result: List[Path] = []
    for tok in ref_tokens:
        tok = tok.strip()
        p = Path(os.path.expanduser(tok)).resolve()
        if p.exists():
            # local file
            result.append(p)
        else:
            if not is_probable_accession(tok):
                raise FileNotFoundError(f"Reference path not found: {tok}")
            # fetch accession
            out_fa = cache_dir / f"{tok}.fasta"
            if not out_fa.exists() or out_fa.stat().st_size == 0:
                fetch_genbank_to_fasta(tok, email, out_fa)
            result.append(out_fa)
    return result


# ----------------------- CSV reader (raw refs + optional reads) -----------------------

def read_rows(csv_path: Path) -> List[Tuple[str, List[str], Optional[List[str]], Optional[str], Optional[str], Optional[str], Optional[str]]]:
    """
    Read mapping file.
    Returns list of tuples:
      (sra_id, [ref_tokens], [segment_names or None], region or None, fastq1 or None, fastq2 or None, fastq_single or None)
    CSV may optionally include read file columns (local paths): fastq/fastq1/fastq2 (many name variants accepted).
    """
    rows = []
    with open(csv_path, newline="") as fh:
        head = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(head, delimiters=",\t")
        except csv.Error:
            dialect = csv.get_dialect("excel")
        reader = csv.DictReader(fh, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError("CSV/TSV appears to have no header row.")
        fm = {k.lower(): k for k in reader.fieldnames}

        sra_k = fm.get("sra") or fm.get("sra_id") or fm.get("sraaccession")

        ref_single_k = fm.get("ref_fasta") or fm.get("reference_fasta") or fm.get("reference")
        refs_multi_k = fm.get("ref_fastas") or fm.get("references") or fm.get("fasta_list")

        names_k = fm.get("segment_names")
        region_k = fm.get("region")

        # read-file candidate keys (lowercased)
        fastq1_k = fm.get("fastq1") or fm.get("fastq_1") or fm.get("reads1") or fm.get("reads_1") or fm.get("read1")
        fastq2_k = fm.get("fastq2") or fm.get("fastq_2") or fm.get("reads2") or fm.get("reads_2") or fm.get("read2")
        fastq_k = fm.get("fastq") or fm.get("fastq_single") or fm.get("reads") or fm.get("read")

        if not sra_k or not (ref_single_k or refs_multi_k):
            raise ValueError(
                "CSV must include SRA (sra/sra_id) and either a single reference column "
                "(ref_fasta/reference_fasta/reference) or a multi-reference column "
                "(ref_fastas/references/fasta_list; semicolon-separated). Entries may be paths or accessions."
            )

        for r in reader:
            sra = r[sra_k].strip()

            # Parse refs as raw tokens (path or accession)
            if refs_multi_k and r.get(refs_multi_k):
                tokens = [x.strip() for x in r[refs_multi_k].split(";") if x.strip()]
            elif ref_single_k and r.get(ref_single_k):
                val = r[ref_single_k].strip()
                tokens = [val] if val else []
            else:
                raise ValueError(f"Row for {sra}: no reference provided.")

            names = None
            if names_k and r.get(names_k):
                names = [x.strip() for x in r[names_k].split(";") if x.strip()]

            region = r[region_k].strip() if region_k and r.get(region_k) else None

            # Optional read-file columns (local paths)
            fastq1 = r[fastq1_k].strip() if fastq1_k and r.get(fastq1_k) and r[fastq1_k].strip() else None
            fastq2 = r[fastq2_k].strip() if fastq2_k and r.get(fastq2_k) and r[fastq2_k].strip() else None
            fastq_single = r[fastq_k].strip() if fastq_k and r.get(fastq_k) and r[fastq_k].strip() else None

            rows.append((sra, tokens, names, region, fastq1, fastq2, fastq_single))
    return rows

# ----------------------- core processing -----------------------

def process_sample_segmented(
    sra_id: str,
    ref_tokens: List[str],
    segment_names: Optional[List[str]],
    output_dir: Path,
    region: Optional[str],
    depth_cutoff: int,
    stats_plots: bool,
    align_threads: int,
    fastp_threads: int,
    email: str,
    read1_path: Optional[str] = None,
    read2_path: Optional[str] = None,
    read_single_path: Optional[str] = None,
) -> dict:
    """Run the full pipeline for one sample. Returns a dict with summary info."""
    outdir = ensure_dir(output_dir)
    sample_prefix = outdir / sra_id
    ref_cache = ensure_dir(outdir / "ref_cache")

    # 0) Resolve references (paths or accessions) -> local FASTA paths
    ref_fastas: List[Path] = resolve_refs_to_fastas(ref_tokens, ref_cache, email)

    # 1) Prepare multi-FASTA reference
    multi_ref = outdir / f"{sra_id}.reference.multi.fasta"

    if len(ref_fastas) == 1:
        src = ref_fastas[0]
        if segment_names:
            rewrite_headers_for_multifasta(src, multi_ref, segment_names)
        else:
            if src != multi_ref:
                shutil.copyfile(src, multi_ref)
    else:
        build_multifasta_from_paths(ref_fastas, multi_ref, segment_names)

    contigs = segment_names if segment_names else list_contigs_from_fasta(multi_ref)
    if not contigs:
        raise ValueError(f"Row {sra_id}: no sequences detected in reference FASTA(s).")

    is_segmented = len(contigs) > 1
    single_multi_input = (len(ref_fastas) == 1 and is_segmented)

    # 2) FASTQs: either provided local reads or download via SRA tools
    fq1 = outdir / f"{sra_id}_1.fastq"
    fq2 = outdir / f"{sra_id}_2.fastq"
    fqs = outdir / f"{sra_id}.fastq"

    # If user provided read paths in CSV, copy them into outdir and skip prefetch
    if read1_path or read2_path or read_single_path:
        if read1_path and read2_path:
            s1 = Path(os.path.expanduser(read1_path)).resolve()
            s2 = Path(os.path.expanduser(read2_path)).resolve()
            if not s1.exists() or not s2.exists():
                raise FileNotFoundError(f"Provided read files for {sra_id} not found: {s1}, {s2}")
            if not fq1.exists():
                shutil.copyfile(s1, fq1)
            if not fq2.exists():
                shutil.copyfile(s2, fq2)
        elif read_single_path:
            s = Path(os.path.expanduser(read_single_path)).resolve()
            if not s.exists():
                raise FileNotFoundError(f"Provided read file for {sra_id} not found: {s}")
            if not fqs.exists():
                shutil.copyfile(s, fqs)
        elif read1_path and not read2_path:
            # single file provided under fastq1 column -> treat as single-end
            s = Path(os.path.expanduser(read1_path)).resolve()
            if not s.exists():
                raise FileNotFoundError(f"Provided read file for {sra_id} not found: {s}")
            if not fqs.exists():
                shutil.copyfile(s, fqs)
        # else: nothing to do
    else:
        # No local reads provided: try to download from SRA if fastq files not already present
        if not fq1.exists() and not fqs.exists():
            run_command(f"prefetch {sra_id}", f"Prefetch {sra_id}")
            run_command(f"fasterq-dump {sra_id} -O {outdir}", f"FASTQ download {sra_id}")

    # 3) Trimming (fastp)
    trim1 = outdir / f"{sra_id}_1.trim.fastq"
    trim2 = outdir / f"{sra_id}_2.trim.fastq"
    fastp_thr_opt = f"--thread {fastp_threads}" if fastp_threads > 0 else ""
    if fq1.exists() and fq2.exists():
        run_command(
            f"fastp -i {fq1} -I {fq2} -o {trim1} -O {trim2} {fastp_thr_opt} "
            f"-h {sample_prefix}_fastp.html -j {sample_prefix}_fastp.json",
            f"Trim paired-end {sra_id}"
        )
        pe = True
    else:
        trim1 = outdir / f"{sra_id}.trim.fastq"
        run_command(
            f"fastp -i {fqs} -o {trim1} {fastp_thr_opt} -h {sample_prefix}_fastp.html -j {sample_prefix}_fastp.json",
            f"Trim single-end {sra_id}"
        )
        pe = False

    # 4) Index & align
    if not bwa_index_present(multi_ref):
        run_command(f"bwa index {multi_ref}", f"BWA index {multi_ref.name}")
    sam = outdir / f"{sra_id}.sam"
    bwa_thr_opt = f"-t {align_threads}" if align_threads > 0 else ""
    if pe:
        run_command(f"bwa mem {bwa_thr_opt} {multi_ref} {trim1} {trim2} > {sam}", f"Align {sra_id} (PE)")
    else:
        run_command(f"bwa mem {bwa_thr_opt} {multi_ref} {trim1} > {sam}", f"Align {sra_id} (SE)")

    bam_sorted = outdir / f"{sra_id}.sorted.bam"
    run_command(f"samtools view -Sb {sam} | samtools sort -o {bam_sorted}", f"Sort BAM {sra_id}")
    run_command(f"samtools index {bam_sorted}", f"Index BAM {sra_id}")

    # 5) Variant calling
    vcf_gz = outdir / f"{sra_id}.vcf.gz"
    region_opt = f"-r {region} " if region else ""
    run_command(
        f"bcftools mpileup -Ou -f {multi_ref} {region_opt}{bam_sorted} | bcftools call -mv -Oz -o {vcf_gz}",
        f"Variant calling {sra_id}"
    )
    run_command(f"bcftools index {vcf_gz}", f"Index VCF {sra_id}")

    # 6) Low-coverage mask BED (from depth)
    mask_bed = outdir / f"{sra_id}.lowcov_dp{depth_cutoff}.bed"
    depth_bytes = run_command(f"samtools depth -a {bam_sorted}", f"Depth {sra_id}", capture_stdout=True)
    if depth_cutoff > 0:
        write_lowcov_bed_from_depth(depth_bytes, depth_cutoff, mask_bed)
        mask_opt = f"-m {mask_bed}"
    else:
        mask_opt = ""

    combined_cons = outdir / f"{sra_id}.consensus.multi.fasta"
    per_seg_files: List[Path] = []

    # 7) Consensus generation (unchanged logic)
    if not is_segmented:
        only = contigs[0]
        seg_cons = outdir / f"{sra_id}.consensus.{re.sub(r'[^A-Za-z0-9._-]', '_', only)}.fasta"
        run_command(
            f"bcftools consensus -f {multi_ref} {mask_opt} {vcf_gz} > {seg_cons}",
            f"Consensus {sra_id} (single-contig: {only})"
        )
        with open(combined_cons, "w") as all_out, open(seg_cons) as fh:
            all_out.write(fh.read())
        per_seg_files.append(seg_cons)
    else:
        if single_multi_input:
            run_command(
                f"bcftools consensus -f {multi_ref} {mask_opt} {vcf_gz} > {combined_cons}",
                f"Consensus {sra_id} (multi-record single FASTA -> combined only)"
            )
        else:
            def _bcftools_supports_r() -> bool:
                try:
                    out = run_command(
                        "bcftools consensus -h 2>&1 | grep -E '\\-r|--regions' || true",
                        "Check bcftools consensus -r support",
                        capture_stdout=True
                    )
                    return bool(out and out.decode().strip())
                except Exception:
                    return False

            supports_r = _bcftools_supports_r()
            run_command(f"samtools faidx {multi_ref}", f"samtools faidx {multi_ref.name}")

            with open(combined_cons, "w") as all_out:
                for seg in contigs:
                    safe_seg = re.sub(r"[^A-Za-z0-9._-]", "_", seg)
                    seg_cons = outdir / f"{sra_id}.consensus.{safe_seg}.fasta"
                    if supports_r:
                        run_command(
                            f"bcftools consensus -f {multi_ref} -r {seg} {mask_opt} {vcf_gz} > {seg_cons}",
                            f"Consensus {sra_id}:{seg}"
                        )
                    else:
                        run_command(
                            f"samtools faidx {multi_ref} {seg} | "
                            f"bcftools consensus {mask_opt} {vcf_gz} > {seg_cons}",
                            f"Consensus {sra_id}:{seg} (streamed)"
                        )
                    per_seg_files.append(seg_cons)
                    with open(seg_cons) as fh:
                        all_out.write(fh.read())

    # 8) Optional stats & plots
    stats_file = None
    plots_dir = None
    if stats_plots:
        stats_file = outdir / f"{sra_id}.vcf.stats"
        plots_dir = outdir / f"{sra_id}_plots"
        run_command(f"bcftools stats {vcf_gz} > {stats_file}", f"VCF stats {sra_id}")
        run_command(f"plot-vcfstats {stats_file} -p {plots_dir}", f"VCF plots {sra_id}")

    # 9) Variant count
    count = int((run_command(f"bcftools view -H {vcf_gz} | wc -l", "Count variants", capture_stdout=True) or b"0").decode().strip() or 0)

    return {
        "Sample": sra_id,
        "Num_Variants": count,
        "VCF": str(vcf_gz),
        "Consensus_Multi": str(combined_cons),
        "Consensus_Segments": "" if (is_segmented and single_multi_input) else ";".join(map(str, per_seg_files)),
        "Stats_File": str(stats_file) if stats_file else "",
        "Plots_Dir": str(plots_dir) if plots_dir else "",
        "Mask_BED": str(mask_bed) if depth_cutoff > 0 else "",
    }

# ----------------------- main -----------------------

def make_parser() -> argparse.ArgumentParser:
    epilog = (
        "Examples:\n"
        "  # Local single multi-FASTA -> only combined multi\n"
        "  python viral_assembly_segmented_v3.py --map samples.csv --email you@org.org --output results\n\n"
        "  # Multiple accessions -> per-seg files + combined multi\n"
        "  python viral_assembly_segmented_v3.py --map samples_multi_acc.csv --email you@org.org --output results --workers 6\n"
    )
    parser = argparse.ArgumentParser(
        prog="viral_assembly_segmented_v3.py",
        description=(
            "Viral Assembly Pipeline: CSV-driven, supports local FASTA paths and GenBank accessions.\n"
            "Builds per-segment and combined consensuses; optimized behavior for single multi-FASTA inputs."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=epilog,
    )
    parser.add_argument("--version", action="version", version="viral-assembly-pipeline 1.4.0")
    req = parser.add_argument_group("Required")
    req.add_argument("--map", required=True, help="CSV/TSV mapping file (see docstring)")
    req.add_argument("--email", required=True, help="Email for NCBI Entrez (used when fetching accessions)")

    opt = parser.add_argument_group("Options")
    opt.add_argument("--output", default="results", help="Output directory (default: results)")
    opt.add_argument("--workers", type=int, default=cpu_count(), help="Parallel workers (default: CPU count)")
    opt.add_argument("--depth-cutoff", type=int, default=10, help="Mask bases with depth < cutoff (0 = disable masking)")
    opt.add_argument("--stats-plots", action="store_true", help="Generate bcftools stats and plot-vcfstats outputs")
    opt.add_argument("--align-threads", type=int, default=0, help="Threads for bwa mem (0 = tool default)")
    opt.add_argument("--fastp-threads", type=int, default=0, help="Threads for fastp (0 = tool default)")
    return parser

def main():
    parser = make_parser()
    args = parser.parse_args()

    outdir = ensure_dir(Path(args.output))

    rows = read_rows(Path(args.map))

    summaries = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [
            ex.submit(
                process_sample_segmented,
                sra_id, ref_tokens, names, outdir, region,
                args.depth_cutoff, args.stats_plots, args.align_threads, args.fastp_threads, args.email,
                read1, read2, read_single
            )
            for (sra_id, ref_tokens, names, region, read1, read2, read_single) in rows
        ]
        for f in as_completed(futs):
            summaries.append(f.result())

    summary_file = outdir / "variant_counts.tsv"
    pd.DataFrame(summaries).to_csv(summary_file, sep="\t", index=False)
    print("\n✅ Pipeline complete")
    print(f"- Variant summary: {summary_file}")
    print(f"- Output dir: {outdir}")

if __name__ == "__main__":
    main()
