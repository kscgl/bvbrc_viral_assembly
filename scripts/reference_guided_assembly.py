#!/usr/bin/env python3

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional, List, Dict

from Bio import Entrez, SeqIO


def run_command(cmd: str, desc: str, capture_stdout: bool = False) -> Optional[bytes]:
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
        raise


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def list_contigs_from_fasta(fa_path: Path) -> List[str]:
    contigs: List[str] = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.append(line[1:].strip().split()[0])
    return contigs


def _safe_segment_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def build_multifasta_from_paths(
    fastas: List[Path],
    out_fa: Path,
    segment_names: Optional[List[str]] = None,
) -> None:
    """Concatenate FASTA files into out_fa, optionally overriding headers with segment_names (record-level)."""
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


def rewrite_headers_for_multifasta(src_multi: Path, out_fa: Path, new_names: List[str]) -> None:
    """Rewrite headers of a single multi-record FASTA to the provided names (1:1)."""
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


def write_lowcov_bed_from_depth(depth_bytes: bytes, cutoff: int, bed_path: Path) -> None:
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

def is_probable_accession(token: str) -> bool:
    """Heuristic: if token doesn't look like a local path and contains no path separators, treat as accession."""
    p = Path(os.path.expanduser(token))
    if p.exists():
        return False
    return ("/" not in token) and ("\\" not in token)


def fetch_genbank_to_fasta(accession: str, email: str, out_fa: Path) -> Path:
    if not email:
        raise ValueError("Entrez email is required when using GenBank references for reference-guided assembly.")
    Entrez.email = email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    if not records:
        raise ValueError(f"No FASTA sequence returned for accession {accession}")
    with open(out_fa, "w") as fh:
        SeqIO.write(records, fh, "fasta")
    return out_fa


def resolve_reference_fasta(token: str, cache_dir: Path, email: Optional[str]) -> Path:
    """
    Resolve a single reference token that is either a local path or a GenBank accession into a local FASTA path.
    """
    ensure_dir(cache_dir)
    tok = token.strip()
    p = Path(os.path.expanduser(tok)).resolve()
    if p.exists():
        return p

    if not is_probable_accession(tok):
        raise FileNotFoundError(f"Reference path not found: {tok}")

    out_fa = cache_dir / f"{tok}.fasta"
    if not out_fa.exists() or out_fa.stat().st_size == 0:
        fetch_genbank_to_fasta(tok, email, out_fa)  # type: ignore[arg-type]
    return out_fa


def resolve_refs_to_fastas(ref_tokens: List[str], cache_dir: Path, email: Optional[str]) -> List[Path]:
    """Resolve a list of tokens (paths or accessions) into local FASTA paths."""
    ensure_dir(cache_dir)
    return [resolve_reference_fasta(tok, cache_dir, email) for tok in ref_tokens]


def _parse_string_list(v) -> Optional[List[str]]:
    """
    Accept either:
      - list[str]
      - semicolon-separated string
    and return a cleaned list, or None.
    """
    if v is None:
        return None
    if isinstance(v, list):
        items = [str(x).strip() for x in v if str(x).strip()]
        return items or None
    if isinstance(v, str):
        items = [x.strip() for x in v.split(";") if x.strip()]
        return items or None
    return None


def _bcftools_consensus_supports_regions() -> bool:
    """
    Detect whether `bcftools consensus` supports -r/--regions.
    We avoid shell grep for portability and tool restrictions.
    """
    try:
        p = subprocess.run(
            ["bcftools", "consensus", "-h"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        out = p.stdout or ""
        # Older bcftools builds don't support -r/--regions at all.
        # To avoid false positives (e.g. "-r" inside other text), require
        # an option line that starts with the short flag.
        for line in out.splitlines():
            stripped = line.lstrip()
            if stripped.startswith("-r ") or stripped.startswith("-r,") or stripped.startswith("--regions"):
                return True
        return False
    except Exception:
        return False


def run_reference_guided(
    job_data: Dict,
    output_dir: str,
    read1: Optional[str] = None,
    read2: Optional[str] = None,
    read_single: Optional[str] = None,
) -> Dict:
    """
    Run a single-sample reference-guided assembly using bwa/bcftools.

    Inputs:
      - job_data: job JSON dict, expected to include:
          strategy: \"reference_guided\"
          reference_type: \"genbank\" or \"fasta\"
          reference_assembly / genbank_accession / fasta_file
          email (optional; required for genbank)
      - read1/read2/read_single: local FASTQ paths prepared by the caller
    Outputs (in output_dir):
      - <sample>.consensus.fasta
      - <sample>.vcf.gz (+ index)
      - QUAST output directory with report.txt and report.html
    Returns:
      - summary dict with keys: consensus_fasta, vcf, quast_txt, quast_html
    """
    outdir = ensure_dir(Path(output_dir))

    sample_name = job_data.get("output_file") or job_data.get("srr_id") or "sample"
    sample_name = str(sample_name)

    reference_type = (job_data.get("reference_type") or "").lower()
    if reference_type not in {"genbank", "fasta"}:
        raise ValueError("reference_type must be 'genbank' or 'fasta' for reference-guided strategy.")

    email = job_data.get("email") or os.environ.get("ENTREZ_EMAIL")

    # Resolve reference FASTA(s)
    ref_cache = ensure_dir(outdir / "ref_cache")

    segment_names = _parse_string_list(job_data.get("segment_names"))
    ref_tokens = _parse_string_list(job_data.get("reference_fastas")) or _parse_string_list(job_data.get("references"))

    ref_token = None
    if reference_type == "genbank":
        if not ref_tokens:
            ref_token = job_data.get("genbank_accession") or job_data.get("reference_assembly")
            if not ref_token:
                raise ValueError("genbank_accession or reference_assembly must be provided for GenBank references.")
    elif reference_type == "fasta":
        if not ref_tokens:
            ref_token = job_data.get("fasta_file") or job_data.get("reference_assembly")
            if not ref_token:
                raise ValueError("fasta_file or reference_assembly must be provided for FASTA references.")

    # Build a canonical reference multi-FASTA in the output directory.
    # - If multiple refs are given (paths or accessions), concatenate them.
    # - If a single multi-record FASTA is given and segment_names provided, rewrite headers.
    # - Otherwise, copy the single reference.
    multi_ref = outdir / f"{sample_name}.reference.multi.fasta"

    if ref_tokens:
        ref_fastas = resolve_refs_to_fastas(ref_tokens, ref_cache, email)
        build_multifasta_from_paths(ref_fastas, multi_ref, segment_names=segment_names)
        contigs = segment_names if segment_names else list_contigs_from_fasta(multi_ref)
        single_multi_input = False
    else:
        ref_fasta = resolve_reference_fasta(str(ref_token), ref_cache, email)
        if segment_names:
            rewrite_headers_for_multifasta(ref_fasta, multi_ref, segment_names)
        else:
            if ref_fasta != multi_ref:
                shutil.copyfile(ref_fasta, multi_ref)
        contigs = segment_names if segment_names else list_contigs_from_fasta(multi_ref)
        single_multi_input = (len(contigs) > 1)

    if not contigs:
        raise ValueError(f"No sequences detected in reference FASTA {multi_ref}")

    is_paired = bool(read1 and read2)

    # Trimming with fastp
    fastp_threads = job_data.get("fastp_threads", 0)
    fastp_thr_opt = f"--thread {fastp_threads}" if fastp_threads and fastp_threads > 0 else ""

    trim1 = outdir / f"{sample_name}_1.trim.fastq"
    trim2 = outdir / f"{sample_name}_2.trim.fastq"
    trim_single = outdir / f"{sample_name}.trim.fastq"

    if is_paired:
        if not (read1 and read2):
            raise ValueError("Both read1 and read2 must be provided for paired-end reference-guided assembly.")
        run_command(
            f"fastp -i {read1} -I {read2} -o {trim1} -O {trim2} {fastp_thr_opt} "
            f"-h {outdir / (sample_name + '_fastp.html')} -j {outdir / (sample_name + '_fastp.json')}",
            f"Trim paired-end {sample_name}",
        )
    else:
        if not read_single and not read1:
            raise ValueError("A single-end read file must be provided for reference-guided assembly.")
        src = read_single or read1  # type: ignore[assignment]
        run_command(
            f"fastp -i {src} -o {trim_single} {fastp_thr_opt} "
            f"-h {outdir / (sample_name + '_fastp.html')} -j {outdir / (sample_name + '_fastp.json')}",
            f"Trim single-end {sample_name}",
        )

    # Align with bwa
    if not bwa_index_present(multi_ref):
        run_command(f"bwa index {multi_ref}", f"BWA index {multi_ref.name}")

    sam = outdir / f"{sample_name}.sam"
    align_threads = job_data.get("align_threads", 0)
    bwa_thr_opt = f"-t {align_threads}" if align_threads and align_threads > 0 else ""

    if is_paired:
        run_command(
            f"bwa mem {bwa_thr_opt} {multi_ref} {trim1} {trim2} > {sam}",
            f"Align {sample_name} (PE)",
        )
    else:
        run_command(
            f"bwa mem {bwa_thr_opt} {multi_ref} {trim_single} > {sam}",
            f"Align {sample_name} (SE)",
        )

    bam_sorted = outdir / f"{sample_name}.sorted.bam"
    run_command(
        f"samtools view -Sb {sam} | samtools sort -o {bam_sorted}",
        f"Sort BAM {sample_name}",
    )
    run_command(f"samtools index {bam_sorted}", f"Index BAM {sample_name}")

    # Variant calling with bcftools
    vcf_gz = outdir / f"{sample_name}.vcf.gz"
    region = job_data.get("region")
    region_opt = f"-r {region} " if region else ""
    run_command(
        f"bcftools mpileup -Ou -f {multi_ref} {region_opt}{bam_sorted} | "
        f"bcftools call -mv -Oz -o {vcf_gz}",
        f"Variant calling {sample_name}",
    )
    run_command(f"bcftools index {vcf_gz}", f"Index VCF {sample_name}")

    # Low-coverage mask BED (depth-based masking similar to original pipeline)
    depth_cutoff = int(job_data.get("depth_cutoff", 10))
    mask_bed = None
    if depth_cutoff > 0:
        mask_bed = outdir / f"{sample_name}.lowcov_dp{depth_cutoff}.bed"
        depth_bytes = run_command(
            f"samtools depth -a {bam_sorted}",
            f"Depth {sample_name}",
            capture_stdout=True,
        ) or b""
        write_lowcov_bed_from_depth(depth_bytes, depth_cutoff, mask_bed)
        mask_opt = f"-m {mask_bed}"
    else:
        mask_opt = ""

    # Consensus
    is_segmented = len(contigs) > 1
    emit_per_segment = bool(job_data.get("per_segment_consensus")) if single_multi_input else True

    consensus_segments: List[str] = []
    if not is_segmented:
        consensus_primary = outdir / f"{sample_name}.consensus.fasta"
        run_command(
            f"bcftools consensus -f {multi_ref} {mask_opt} {vcf_gz} > {consensus_primary}",
            f"Consensus {sample_name}",
        )
        # For backward-compatibility with the original pipeline, also emit a
        # contig-named consensus when there is exactly one contig and no
        # explicit segment_names override, e.g.
        #   <sample>.consensus.NC_001474.2.fasta
        if len(contigs) == 1:
            only = contigs[0]
            legacy_name = outdir / f"{sample_name}.consensus.{only}.fasta"
            if legacy_name != consensus_primary:
                shutil.copyfile(consensus_primary, legacy_name)
        # Also emit a multi-consensus file to mirror the original pipeline.
        consensus_multi = outdir / f"{sample_name}.consensus.multi.fasta"
        if consensus_multi != consensus_primary:
            shutil.copyfile(consensus_primary, consensus_multi)
    else:
        consensus_primary = outdir / f"{sample_name}.consensus.multi.fasta"
        if single_multi_input and not emit_per_segment:
            run_command(
                f"bcftools consensus -f {multi_ref} {mask_opt} {vcf_gz} > {consensus_primary}",
                f"Consensus {sample_name} (combined)",
            )
        else:
            supports_r = _bcftools_consensus_supports_regions()
            run_command(f"samtools faidx {multi_ref}", f"samtools faidx {multi_ref.name}")
            with open(consensus_primary, "w") as all_out:
                for seg in contigs:
                    safe_seg = _safe_segment_name(seg)
                    seg_cons = outdir / f"{sample_name}.consensus.{safe_seg}.fasta"
                    if supports_r:
                        try:
                            run_command(
                                f"bcftools consensus -f {multi_ref} -r {seg} {mask_opt} {vcf_gz} > {seg_cons}",
                                f"Consensus {sample_name}:{seg}",
                            )
                        except subprocess.CalledProcessError:
                            # Some bcftools builds print -r in help but don't accept it.
                            # Fall back to the streamed approach for this and remaining segments.
                            supports_r = False
                            run_command(
                                f"samtools faidx {multi_ref} {seg} | bcftools consensus {mask_opt} {vcf_gz} > {seg_cons}",
                                f"Consensus {sample_name}:{seg} (streamed)",
                            )
                    else:
                        run_command(
                            f"samtools faidx {multi_ref} {seg} | bcftools consensus {mask_opt} {vcf_gz} > {seg_cons}",
                            f"Consensus {sample_name}:{seg} (streamed)",
                        )
                    consensus_segments.append(str(seg_cons))
                    with open(seg_cons) as fh:
                        all_out.write(fh.read())

    # QUAST on consensus
    quast_dir = outdir / "quast_reference_guided"
    ensure_dir(quast_dir)
    quast_cmd = (
        f"quast.py -o {quast_dir} -t 4 --min-contig 200 {consensus_primary}"
    )
    try:
        run_command(quast_cmd, f"QUAST reference-guided {sample_name}")
        quast_txt = quast_dir / "report.txt"
        quast_html_rel = "quast_reference_guided/report.html"
    except Exception:
        quast_txt = None
        quast_html_rel = ""

    return {
        "sample": sample_name,
        "consensus_fasta": str(consensus_primary),
        "consensus_multi_fasta": str(consensus_primary) if is_segmented else "",
        "consensus_segments": consensus_segments,
        "vcf": str(vcf_gz),
        "quast_txt": str(quast_txt) if quast_txt and quast_txt.exists() else "",
        "quast_html": quast_html_rel,
    }

