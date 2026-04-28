#!/usr/bin/env python3

import os
import re
import sys
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


def _reference_label_from_token(token: str) -> str:
    """
    Derive a stable file-label from a reference token (path or accession).
    For paths, use the basename without extension. For accessions/other
    tokens, use the token itself.
    """
    t = str(token).strip()
    p = Path(os.path.expanduser(t))
    # Accessions like NC_045512.2 contain dots but are not file paths.
    # Only apply suffix/stem logic when this actually looks like a path.
    if p.exists() or ("/" in t) or ("\\" in t):
        return _safe_segment_name(p.stem)
    return _safe_segment_name(t)

def _rewrite_consensus_headers(cons_fa: Path, sample_name: str) -> None:
    """
    Rewrite FASTA headers in a consensus file so that each header starts with
    '<sample>.consensus.<orig_header_token>'.
    """
    if not cons_fa.exists():
        return
    lines: List[str] = []
    with open(cons_fa) as inp:
        for line in inp:
            if line.startswith(">"):
                orig = line[1:].strip().split()[0]
                # Avoid double-prefixing if already in the desired form.
                if orig.startswith(f"{sample_name}.consensus."):
                    new_header = orig
                else:
                    new_header = f"{sample_name}.consensus.{orig}"
                lines.append(f">{new_header}\n")
            else:
                lines.append(line)
    with open(cons_fa, "w") as out:
        out.writelines(lines)


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


def is_bvbrc_genome_id(token: str) -> bool:
    """Return True if token looks like a BV-BRC genome ID (digits.digits, e.g. '11060.9352')."""
    return bool(re.match(r"^\d+\.\d+$", token.strip()))


def fetch_bvbrc_genome_to_fasta(genome_id: str, out_fa: Path) -> Path:
    """
    Download the nucleotide FASTA for a BV-BRC genome ID (e.g. '11060.9352')
    using the ``p3-genome-fasta`` CLI tool (available in the BV-BRC runtime).
    """
    genome_id = genome_id.strip()

    try:
        out = subprocess.check_output(
            ["p3-genome-fasta", genome_id],
            stderr=subprocess.DEVNULL,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        raise ValueError(
            f"p3-genome-fasta failed for genome_id {genome_id!r}: {e}. "
            f"Verify that the genome exists and is publicly accessible."
        ) from e

    if not out or not out.strip().startswith(b">"):
        raise ValueError(
            f"p3-genome-fasta returned no FASTA content for genome_id {genome_id!r}. "
            f"Verify that the genome exists and is publicly accessible."
        )

    with open(out_fa, "wb") as fh:
        fh.write(out)
    print(f"[BV-BRC] Downloaded genome {genome_id} via p3-genome-fasta -> {out_fa}")
    return out_fa


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
    Resolve a single reference token into a local FASTA path.

    Resolution order:
    1. Local file path  — returned as-is if the file exists.
    2. BV-BRC genome ID (digits.digits, e.g. '11060.9352') — fetched via
       ``fetch_bvbrc_genome_to_fasta`` and cached under *cache_dir*.
    3. GenBank accession — fetched via Entrez and cached under *cache_dir*.
    """
    ensure_dir(cache_dir)
    tok = token.strip()
    p = Path(os.path.expanduser(tok)).resolve()
    if p.exists():
        return p

    if is_bvbrc_genome_id(tok):
        out_fa = cache_dir / f"bvbrc_{tok}.fasta"
        if not out_fa.exists() or out_fa.stat().st_size == 0:
            fetch_bvbrc_genome_to_fasta(tok, out_fa)
        return out_fa

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


def run_reference_guided_resolved(
    *,
    output_dir: str,
    sample_name: str,
    reference_type: str,
    reference_tokens: List[str],
    email: Optional[str],
    read1: Optional[str] = None,
    read2: Optional[str] = None,
    read_single: Optional[str] = None,
    segment_names: Optional[List[str]] = None,
    per_segment_consensus: Optional[bool] = None,
    align_threads: int = 0,
    fastp_threads: int = 0,
    region: Optional[str] = None,
    depth_cutoff: int = 10,
) -> Dict:
    """
    Preferred entrypoint: resolved inputs from run_viral_assembly (local paths, accessions).
    """
    if reference_type not in {"genome", "genbank", "fasta"}:
        raise ValueError("reference_type must be 'genome', 'genbank', or 'fasta'.")
    if not reference_tokens:
        raise ValueError("reference_tokens must contain at least one value.")

    return run_reference_guided_core(
        output_dir=output_dir,
        sample_name=sample_name,
        reference_type=reference_type,
        reference_tokens=reference_tokens,
        email=email,
        read1=read1,
        read2=read2,
        read_single=read_single,
        segment_names=segment_names,
        per_segment_consensus=per_segment_consensus,
        align_threads=align_threads,
        fastp_threads=fastp_threads,
        region=region,
        depth_cutoff=depth_cutoff,
    )


def run_reference_guided(
    job_data: Dict,
    output_dir: str,
    read1: Optional[str] = None,
    read2: Optional[str] = None,
    read_single: Optional[str] = None,
) -> Dict:
    """
    Entrypoint that accepts raw job JSON (requires reference_genbank_accession for genbank
    or reference_fasta_file for fasta — see run_viral_assembly._resolve_reference_inputs).
    Stages SRA reads via run_viral_assembly.ensure_sra_fastqs when sra_id is set and no reads were passed.
    Prefer run_reference_guided_resolved from the run script for new code.
    """
    _scripts_dir = os.path.dirname(os.path.abspath(__file__))
    if _scripts_dir not in sys.path:
        sys.path.insert(0, _scripts_dir)

    import run_viral_assembly as rva

    sra_run = job_data.get("sra_id") or job_data.get("srr_id")
    if sra_run and not (read1 or read2 or read_single):
        read1, read2, read_single = rva.ensure_sra_fastqs(
            output_dir, str(sra_run), job_data
        )

    ref_type, ref_tokens = rva._resolve_reference_inputs(job_data)
    sample_name = str(job_data.get("output_file") or job_data.get("sra_id") or job_data.get("srr_id") or "sample")
    return run_reference_guided_resolved(
        output_dir=output_dir,
        sample_name=sample_name,
        reference_type=ref_type,
        reference_tokens=ref_tokens,
        email=job_data.get("email") or os.environ.get("ENTREZ_EMAIL"),
        read1=read1,
        read2=read2,
        read_single=read_single,
        segment_names=_parse_string_list(job_data.get("segment_names")),
        per_segment_consensus=job_data.get("per_segment_consensus"),
        align_threads=int(job_data.get("align_threads", 0) or 0),
        fastp_threads=int(job_data.get("fastp_threads", 0) or 0),
        region=job_data.get("region"),
        depth_cutoff=int(job_data.get("depth_cutoff", 10) or 10),
    )


def run_reference_guided_core(
    *,
    output_dir: str,
    sample_name: str,
    reference_type: str,
    reference_tokens: List[str],
    email: Optional[str],
    read1: Optional[str] = None,
    read2: Optional[str] = None,
    read_single: Optional[str] = None,
    segment_names: Optional[List[str]] = None,
    per_segment_consensus: Optional[bool] = None,
    align_threads: int = 0,
    fastp_threads: int = 0,
    region: Optional[str] = None,
    depth_cutoff: int = 10,
) -> Dict:
    """
    Reference-guided pipeline: trim, align, call variants, consensus, QUAST.
    Expects local FASTQ paths and reference tokens (GenBank accessions or local FASTA paths).
    """
    outdir = ensure_dir(Path(output_dir))

    reference_type = (reference_type or "").lower()
    if reference_type not in {"genome", "genbank", "fasta"}:
        raise ValueError("reference_type must be 'genome', 'genbank', or 'fasta' for reference-guided strategy.")

    email = email or os.environ.get("ENTREZ_EMAIL")

    # Use a dedicated output-files folder for non-primary artifacts.
    output_files_dir = ensure_dir(outdir / "output_files")
    # Resolve/cache references under output_files.
    ref_cache = ensure_dir(output_files_dir / "ref_cache")

    if len(reference_tokens) > 1:
        ref_tokens = reference_tokens
        ref_token = None
    else:
        ref_tokens = None
        ref_token = reference_tokens[0]

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
        primary_ref_label = _reference_label_from_token(ref_tokens[0]) if len(ref_tokens) == 1 else "multi_ref"
    else:
        ref_fasta = resolve_reference_fasta(str(ref_token), ref_cache, email)
        if segment_names:
            rewrite_headers_for_multifasta(ref_fasta, multi_ref, segment_names)
        else:
            if ref_fasta != multi_ref:
                shutil.copyfile(ref_fasta, multi_ref)
        contigs = segment_names if segment_names else list_contigs_from_fasta(multi_ref)
        single_multi_input = (len(contigs) > 1)
        primary_ref_label = _reference_label_from_token(str(ref_token))

    if not contigs:
        raise ValueError(f"No sequences detected in reference FASTA {multi_ref}")

    is_paired = bool(read1 and read2)

    # Trimming with fastp
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
    region_opt = f"-r {region} " if region else ""
    run_command(
        f"bcftools mpileup -Ou -f {multi_ref} {region_opt}{bam_sorted} | "
        f"bcftools call -mv -Oz -o {vcf_gz}",
        f"Variant calling {sample_name}",
    )
    run_command(f"bcftools index {vcf_gz}", f"Index VCF {sample_name}")

    # Low-coverage mask BED (depth-based masking similar to original pipeline)
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
    emit_per_segment = bool(per_segment_consensus) if single_multi_input else True

    consensus_segments: List[str] = []
    legacy_consensus_files: List[Path] = []
    root_consensus_alias: Optional[Path] = None
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
                legacy_consensus_files.append(legacy_name)
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
        # Move report.html up one level so it sits next to consensus FASTAs.
        quast_html_src = quast_dir / "report.html"
        quast_html_top = outdir / "report.html"
        if quast_html_src.exists():
            shutil.move(str(quast_html_src), str(quast_html_top))
            quast_html_rel = "report.html"
        else:
            quast_html_rel = ""
    except Exception:
        quast_txt = None
        quast_html_rel = ""

    # Normalize consensus headers to start with "<sample>.consensus.<name>"
    _rewrite_consensus_headers(consensus_primary, sample_name)
    for legacy_path in legacy_consensus_files:
        _rewrite_consensus_headers(legacy_path, sample_name)
    for seg_path in consensus_segments:
        _rewrite_consensus_headers(Path(seg_path), sample_name)

    # For FASTA-based jobs, expose a root-level consensus named using the
    # reference FASTA label (e.g. <sample>.consensus.influenza_ref.fasta).
    if reference_type == "fasta":
        root_consensus_alias = outdir / f"{sample_name}.consensus.{primary_ref_label}.fasta"
        if root_consensus_alias != consensus_primary:
            shutil.copyfile(consensus_primary, root_consensus_alias)
            _rewrite_consensus_headers(root_consensus_alias, sample_name)

    # Single-accession GenBank/genome input: keep only contig-named consensus at root.
    single_genbank_single_ref = (
        reference_type in {"genbank", "genome"}
        and len(contigs) == 1
        and bool(ref_token)
        and not ref_tokens
    )

    # After all files are generated, move non-kept top-level files into
    # "output_files", leaving selected consensus FASTAs (and report.html) at
    # the top level.
    output_subdir = output_files_dir
    # Build a set of consensus filenames to keep at the top level.
    keep_names = set()
    if reference_type == "fasta":
        if root_consensus_alias is not None:
            keep_names.add(root_consensus_alias.name)
        keep_names.update(Path(p).name for p in consensus_segments)
    elif single_genbank_single_ref:  # genbank or genome, single contig
        keep_names.update(Path(p).name for p in legacy_consensus_files)
    else:
        keep_names.add(Path(consensus_primary).name)
        keep_names.update(Path(p).name for p in legacy_consensus_files)
        keep_names.update(Path(p).name for p in consensus_segments)
    keep_names.add("report.html")
    keep_names.add("AssemblyReport.html")
    for item in outdir.iterdir():
        if item.is_dir():
            # Keep subdirectories where they are.
            continue
        if item.name in keep_names:
            continue
        # Move everything else into output/.
        shutil.move(str(item), str(output_subdir / item.name))

    returned_consensus = str(consensus_primary)
    if reference_type == "fasta" and root_consensus_alias is not None:
        returned_consensus = str(root_consensus_alias)
    elif single_genbank_single_ref and legacy_consensus_files:  # genbank or genome, single contig
        returned_consensus = str(legacy_consensus_files[0])

    return {
        "sample": sample_name,
        "consensus_fasta": returned_consensus,
        "consensus_multi_fasta": str(consensus_primary) if is_segmented else "",
        "consensus_segments": consensus_segments,
        "vcf": str(vcf_gz),
        "quast_txt": str(quast_txt) if quast_txt and quast_txt.exists() else "",
        "quast_html": quast_html_rel,
    }

