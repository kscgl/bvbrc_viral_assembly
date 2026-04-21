#!/usr/bin/env python

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from time import strftime, localtime

from reference_guided_assembly import run_reference_guided_resolved

#
# Determine paths
#
def first_existing_path(*paths):
  for p in paths:
    if p and os.path.exists(p):
      return p
  return paths[-1]  # fallback

top = os.getenv("KB_TOP")
# Allow local/dev runs without KB_TOP by falling back to repo root.
if not top:
  top = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

ASSEMBLY_REPORT_TEMPLATE = first_existing_path(
    os.path.join(top, "lib", "viral_assembly_report_template.html"),
    os.path.join(top, "modules", "bvbrc_viral_assembly", "lib", "viral_assembly_report_template.html"),
    "/home/ac.mkuscuog/git/dev_container/modules/bvbrc_viral_assembly/lib/viral_assembly_report_template.html",
)

FLU_AD_INIT = first_existing_path(
    os.path.join(top, "lib", "flu_ad_init.sh"),
    os.path.join(top, "modules", "bvbrc_viral_assembly", "lib", "flu_ad_init.sh"),
    "/home/ac.mkuscuog/git/dev_container/modules/bvbrc_viral_assembly/lib/flu_ad_init.sh",
)

DEFAULT_STRATEGY = "IRMA"
DEFAULT_IRMA_MODULE = "FLU"
FLU_FALLBACK_MODULE = "FLU-lowQC"
SRA_FETCH_RETRIES = 5
SRA_RETRY_DELAY_SEC = 10


def fetch_file_from_ws(ws_path, local_path):
  """Fetch a file from workspace to local path using p3-cp."""
  # Local fallback: for dev/testing, allow ws_path to be a real local file.
  if ws_path:
    ws_str = str(ws_path)
    ws_expanded = os.path.expanduser(ws_str)
    local_parent = os.path.dirname(local_path)
    if local_parent:
      os.makedirs(local_parent, exist_ok=True)

    # Absolute local path?
    if os.path.isabs(ws_expanded) and os.path.exists(ws_expanded):
      shutil.copyfile(ws_expanded, local_path)
      print(f"Copied local file {ws_expanded} -> {local_path}")
      return True

    # Relative local path? Treat it as relative to repo root (KB_TOP / script parent).
    if not os.path.isabs(ws_expanded) and os.path.exists(os.path.join(top, ws_expanded)):
      src = os.path.join(top, ws_expanded)
      shutil.copyfile(src, local_path)
      print(f"Copied local relative file {src} -> {local_path}")
      return True

  cmd = ["p3-cp", f"ws:{ws_path}", local_path]
  try:
    print(f"Fetching file from {ws_path} to {local_path}")
    subprocess.run(cmd, check=True)
    print(f"File successfully copied to {local_path}")
  except subprocess.CalledProcessError as e:
    print(f"Error fetching file from {ws_path}: {e}")
    return False
  return True


def fetch_fastqs_from_sra(sra_id, temp_dir="/tmp", output_dir="sra_fastqs"):
  """Download FASTQs with ``p3-sra``. Returns (read1_or_single, read2); single-end is (path, None)."""
  os.makedirs(output_dir, exist_ok=True)
  cmd = ["p3-sra", "--id", sra_id, "--out", output_dir]

  for attempt in range(1, SRA_FETCH_RETRIES + 1):
    print(f"Attempt {attempt}: Fetching FASTQs for SRA ID {sra_id} using p3-sra")
    try:
      subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
      print(f"Error fetching FASTQs for SRA ID {sra_id}: {e}")
      return None, None

    r1 = os.path.join(output_dir, f"{sra_id}_1.fastq")
    r2 = os.path.join(output_dir, f"{sra_id}_2.fastq")
    r_single = os.path.join(output_dir, f"{sra_id}.fastq")

    if os.path.exists(r1) or os.path.exists(r2):
      return (
        r1 if os.path.exists(r1) else None,
        r2 if os.path.exists(r2) else None,
      )

    if os.path.exists(r_single):
      return r_single, None

    print(f"FASTQ files not found after attempt {attempt}. Retrying in {SRA_RETRY_DELAY_SEC} seconds...")
    time.sleep(SRA_RETRY_DELAY_SEC)

  print(f"No valid FASTQ files found for SRA ID {sra_id}.")
  return None, None


def ensure_sra_fastqs(output_dir, sra_run_accession, job_data=None):
  """
  Ensure FASTQs for *sra_run_accession* (e.g. SRRnnn) exist under *output_dir* (reuse on disk, else ``p3-sra``).

  job_data optional keys ``download_sra_from_prefetch`` / ``download_fastqs_from_sra`` (default True):
  if either is False, download is skipped and files must already exist.

  Returns (read1, read2, read_single).
  """
  if job_data is None:
    job_data = {}
  os.makedirs(output_dir, exist_ok=True)
  sid = str(sra_run_accession)

  dl_p = job_data.get("download_sra_from_prefetch")
  dl_p = True if dl_p is None else bool(dl_p)
  dl_f = job_data.get("download_fastqs_from_sra")
  dl_f = True if dl_f is None else bool(dl_f)
  allow_fetch = dl_p and dl_f

  fq1 = os.path.join(output_dir, f"{sid}_1.fastq")
  fq2 = os.path.join(output_dir, f"{sid}_2.fastq")
  fqs = os.path.join(output_dir, f"{sid}.fastq")

  have_pe = os.path.exists(fq1) and os.path.exists(fq2)
  have_se = os.path.exists(fqs)
  have_fastqs = have_pe or have_se

  if not have_fastqs:
    if not allow_fetch:
      raise RuntimeError(
        f"No FASTQs for {sid} under {output_dir} and SRA download is disabled "
        "(download_sra_from_prefetch / download_fastqs_from_sra). "
        "Provide reads or enable download."
      )
    r1, r2 = fetch_fastqs_from_sra(sid, output_dir=output_dir)
    if r2:
      return r1, r2, None
    if r1:
      return None, None, r1
    raise RuntimeError(f"p3-sra did not produce FASTQs for {sid} under {output_dir}.")

  if have_pe:
    return fq1, fq2, None
  return None, None, fqs


def _parse_string_list(v):
  if v is None:
    return None
  if isinstance(v, list):
    vals = [str(x).strip() for x in v if str(x).strip()]
    return vals or None
  if isinstance(v, str):
    vals = [x.strip() for x in v.split(";") if x.strip()]
    return vals or None
  return None


def _job_sra_run_accession(job_data):
  """SRA run accession (e.g. SRR123). Canonical JSON key ``sra_id``; legacy ``srr_id`` accepted."""
  return job_data.get("sra_id") or job_data.get("srr_id") or None


def _resolve_reference_inputs(job_data):
  """
  Read reference tokens from job JSON (runner-only). No legacy key aliases.

  - reference_type genbank: require reference_genbank_accession (string, list, or ';'-separated accessions).
  - reference_type fasta: require reference_fasta_file (path string, list, or ';'-separated paths).
  """
  reference_type = str(job_data.get("reference_type") or "").lower()
  if reference_type not in {"genbank", "fasta"}:
    raise ValueError("reference_type must be 'genbank' or 'fasta' for reference-guided strategy.")

  if reference_type == "genbank":
    tokens = _parse_string_list(job_data.get("reference_genbank_accession"))
    if not tokens:
      raise ValueError(
        "reference_genbank_accession is required for GenBank references "
        "(e.g. 'NC_045512.2' or 'NC_045512.2;MN908947.3')."
      )
    return reference_type, tokens

  tokens = _parse_string_list(job_data.get("reference_fasta_file"))
  if not tokens:
    raise ValueError(
      "reference_fasta_file is required for FASTA references (workspace path(s) or local path(s))."
    )
  return reference_type, tokens


def concatenate_fasta_files(fasta_dir, output_fasta):
    """
    Concatenate all FASTA files in the output directory into one file.

    Args:
        fasta_dir (str): Directory containing output FASTA files.
        output_fasta (str): Path to the concatenated output FASTA file.
    """
    fasta_files = [
        os.path.join(fasta_dir, f)
        for f in os.listdir(fasta_dir)
        if f.endswith(".fasta") or f.endswith(".fa")
    ]

    if not fasta_files:
        print("No FASTA files found in the output directory.")
        return False

    print(f"Concatenating {len(fasta_files)} FASTA files into {output_fasta}...")
    with open(output_fasta, "w") as outfile:
        for fasta in fasta_files:
            with open(fasta, "r") as infile:
                outfile.write(infile.read())
    print("Concatenation complete.")
    return True


def run_quast(quast_output_dir, fasta_file, threads=12, min_contig=200):
    """
    Run QUAST on the concatenated FASTA file from IRMA output.

    Args:
        quast_output_dir (str): Directory to save QUAST output.
        fasta_file (str): Path to the FASTA file to analyse.
        threads (int): Number of threads to use for QUAST. Defaults to 12.
        min_contig (int): Minimum contig size for QUAST analysis. Defaults to 200.
    """
    os.makedirs(quast_output_dir, exist_ok=True)

    quast_cmd = [
        "quast.py",
        "-o", quast_output_dir,
        "-t", str(threads),
        "--min-contig", str(min_contig),
        fasta_file
    ]

    try:
        print(f"Running QUAST with command: {' '.join(quast_cmd)}")
        subprocess.run(quast_cmd, check=True)
        print("QUAST analysis completed successfully!")
    except subprocess.CalledProcessError as e:
        print(f"Error running QUAST: {e}")


def has_assembly_output(output_dir):
  """Check if IRMA produced any non-empty FASTA files in the output directory."""
  if not os.path.isdir(output_dir):
    return False
  for f in os.listdir(output_dir):
    if (f.endswith(".fasta") or f.endswith(".fa")) and os.path.getsize(os.path.join(output_dir, f)) > 0:
      return True
  return False


def run_irma(mode, input_file1, input_file2=None, output_dir="output"):
  irma_cmd = ["IRMA", mode, input_file1]
  if input_file2:
    irma_cmd.append(input_file2)
  irma_cmd.append(output_dir)
  if mode == "FLU_AD":
    irma_cmd.extend(["--external-config", FLU_AD_INIT])

  try:
    print(f"Running IRMA with command: {' '.join(irma_cmd)}")
    result = subprocess.run(irma_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    irma_output = result.stdout
    print(irma_output)
    if has_assembly_output(output_dir):
      print("IRMA run completed successfully with assembly output!")
      return True, irma_output
    else:
      print("IRMA run completed but no assembly output was generated.")
      return False, irma_output
  except subprocess.CalledProcessError as e:
    print(f"Error running IRMA: {e}")
    return False, ""


# Tools that print version info to stderr (or with no args) rather than --version stdout.
_VERSION_STDERR_TOOLS = {"bwa"}

def get_software_version(software):
  try:
    if software in _VERSION_STDERR_TOOLS:
      # e.g. bwa: called with no args, version line appears in stderr
      result = subprocess.run(
        [software],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
      )
      output = result.stderr
    else:
      result = subprocess.run(
        [software, "--version"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
      )
      output = result.stdout or result.stderr

    # Decode tolerantly — some tools (e.g. samtools) include non-UTF-8 chars like ©
    text = output.decode("utf-8", errors="replace")

    # bwa embeds version as "Version: X.X.X" inside the usage block; strip the label
    if software in _VERSION_STDERR_TOOLS:
      for line in text.splitlines():
        if "version:" in line.lower():
          # Extract just the version string after the colon, e.g. "0.7.17-r1188"
          version = line.split(":", 1)[-1].strip()
          print(f"{software} version: {version}")
          return version
      return "Unknown"

    version = text.splitlines()[0].strip()
    print(f"{software} version: {version}")
    return version
  except Exception as e:
    print(f"Error fetching {software} version: {e}")
    return "Unknown"


def move_fasta_files(source_folder, destination_folder):
  for filename in os.listdir(source_folder):
    if filename.endswith(".fasta"):
      source_path = os.path.join(source_folder, filename)
      destination_path = os.path.join(destination_folder, filename)
      shutil.move(source_path, destination_path)
      print(f"Moved: {filename} -> {destination_folder}")


def generate_html_report(details):
    with open(ASSEMBLY_REPORT_TEMPLATE, "r") as template_file:
      template = template_file.read()

    # Set dynamic values
    report_date = strftime("%a, %d %b %Y %H:%M:%S", localtime())

    # Handle Quast section logic
    quast_html = details.get("quast_html", "")
    quast_txt_path = details.get("quast_txt", "")
    quast_section = ""
    if quast_html:
      quast_section = f"""
      <section>
          <h2>Quast Report</h2>
          <a href="{quast_html}">Quast HTML Report</a>
      """
      if quast_txt_path:
        with open(quast_txt_path, "r") as f:
          quast_txt_content = f.read()
        quast_section += f"""
        <pre class="preformatted">{quast_txt_content}</pre>
        """
      quast_section += "</section>"

    # Handle Tools table
    tools = details.get("tools", {})
    tools_table = "\n".join(
        f"<tr><td>{tool}</td><td>{version}</td></tr>"
        for tool, version in tools.items()
    )

    # Handle Errors
    errors_section = ""
    if "errors" in details and details["errors"]:
      errors_section = """
      <section style='color: #bb0505;'>
        <h2>Error</h2>
        <ul>
          <li>Assembly may have failed due to poor data quality or not having enough reads. We recommend running the FASTQC tool to assess read coverage and quality before retrying the assembly. You can access the FASTQC tool at <a href="https://www.bv-brc.org/app/FastqUtil" target="_blank">BV-BRC Fastq Utilities</a>, and the user guide is available <a href="https://www.bv-brc.org/docs/quick_references/services/fastq_utilities_service.html" target="_blank">here</a>.</li>
        </ul>
      </section>
      """

    # Handle Notes
    notes_section = ""
    notes = details.get("notes", [])
    if notes:
      notes_items = "\n".join(f"          <li>{note}</li>" for note in notes)
      notes_section = f"""
      <section style='color: #1a5276;'>
        <h2>Note</h2>
        <ul>
{notes_items}
        </ul>
      </section>
      """

    # Replace placeholders
    html_content = template.replace("{{ report_date }}", report_date)
    html_content = html_content.replace("{{ notes_section }}", notes_section)
    html_content = html_content.replace("{{ quast_section }}", quast_section)
    html_content = html_content.replace("{{ tools_table }}", tools_table)
    html_content = html_content.replace("{{ errors_section }}", errors_section)

    return html_content


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Viral Assembly Script")
  parser.add_argument("-j", "--jfile", help="JSON file for job configuration", required=True)
  parser.add_argument("-o", "--output", help="Output directory. Defaults to current directory", required=False, default=".")

  args = parser.parse_args()

  # Load job data
  job_data = None
  try:
    with open(args.jfile, "r") as j:
      job_data = json.load(j)
  except Exception as e:
    print(f"Error loading JSON job file:\n {e}")
    traceback.print_exc(file=sys.stderr)
    sys.exit(-1)

  if not job_data:
    print("Job data is empty. Exiting.")
    sys.exit(-1)

  print("Loaded job data:", job_data)

  paired_end_lib = job_data.get("paired_end_lib", {})
  single_end_lib = job_data.get("single_end_lib", {})
  sra_id = _job_sra_run_accession(job_data)
  strategy = job_data.get("strategy", DEFAULT_STRATEGY)
  strategy_lower = str(strategy).lower()
  if strategy == "auto":
    strategy = DEFAULT_STRATEGY
  module = job_data.get("module", DEFAULT_IRMA_MODULE)

  # Input validation:
  # - IRMA: exactly one of (paired_end_lib, single_end_lib, sra_id)
  # - reference_guided: allow sra_id together with paired_end_lib/single_end_lib
  #   so the backend can stage reads while still using provided FASTQs.
  if strategy_lower == "reference_guided":
    if paired_end_lib and single_end_lib:
      print("Error: Please provide only one of paired_end_lib or single_end_lib (not both).")
      sys.exit(-1)
    if not (paired_end_lib or single_end_lib or sra_id):
      print("Error: Please provide paired_end_lib, single_end_lib, or sra_id.")
      sys.exit(-1)
  else:
    inputs_provided = sum(bool(x) for x in [paired_end_lib, single_end_lib, sra_id])
    if inputs_provided != 1:
      print("Error: Please provide exactly one of Paired End Library, Single End Library, or SRA run accession (sra_id).")
      sys.exit(-1)

  # Setup output directory
  output_dir = args.output
  output_dir = os.path.abspath(output_dir)
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  os.chdir(output_dir)

  # Initialise report details (shared by both strategies)
  report_details = {
    "errors": [],
    "notes": []
  }

  if strategy_lower == "reference_guided":
    try:
      reference_type, reference_tokens = _resolve_reference_inputs(job_data)
    except Exception as e:
      print(f"Error: invalid reference input for reference-guided strategy: {e}")
      sys.exit(-1)

    # Stage workspace / local reference FASTA paths to scratch (p3-cp or dev copy).
    if reference_type == "fasta":
      ref_dir = os.path.join(output_dir, "reference_inputs")
      os.makedirs(ref_dir, exist_ok=True)
      staged_tokens = []
      for i, tok in enumerate(reference_tokens):
        ext = os.path.splitext(str(tok))[1] or ".fasta"
        local_ref = os.path.join(ref_dir, f"ref_{i}{ext}")
        if not fetch_file_from_ws(tok, local_ref):
          print(f"Error: Failed to stage reference FASTA from {tok!r}.")
          sys.exit(-1)
        staged_tokens.append(local_ref)
      reference_tokens = staged_tokens

    # Prepare local read paths for reference-guided assembly
    local_read1 = None
    local_read2 = None
    local_single = None

    if paired_end_lib:
      read1 = paired_end_lib.get("read1")
      read2 = paired_end_lib.get("read2")
      if not read1 or not read2:
        print("Error: Missing reads for paired-end library.")
        sys.exit(-1)
      local_read1 = os.path.join(output_dir, "reads_1.fastq")
      local_read2 = os.path.join(output_dir, "reads_2.fastq")
      if not (fetch_file_from_ws(read1, local_read1) and fetch_file_from_ws(read2, local_read2)):
        print("Error: Failed to fetch paired-end reads.")
        sys.exit(-1)

    elif single_end_lib:
      read = single_end_lib.get("read")
      if not read:
        print("Error: Missing read for single-end library.")
        sys.exit(-1)
      local_single = os.path.join(output_dir, "reads.fastq")
      if not fetch_file_from_ws(read, local_single):
        print("Error: Failed to fetch single-end read.")
        sys.exit(-1)

    elif sra_id:
      try:
        local_read1, local_read2, local_single = ensure_sra_fastqs(
          output_dir, str(sra_id), job_data
        )
      except Exception as e:
        print(f"Error: SRA staging failed: {e}")
        traceback.print_exc(file=sys.stderr)
        sys.exit(-1)

    # Run reference-guided pipeline
    try:
      summary = run_reference_guided_resolved(
        output_dir=output_dir,
        sample_name=str(job_data.get("output_file") or sra_id or "sample"),
        reference_type=reference_type,
        reference_tokens=reference_tokens,
        email=(job_data.get("email") or os.environ.get("ENTREZ_EMAIL")),
        read1=local_read1,
        read2=local_read2,
        read_single=local_single,
        segment_names=_parse_string_list(job_data.get("segment_names")),
        per_segment_consensus=job_data.get("per_segment_consensus"),
        align_threads=int(job_data.get("align_threads", 0) or 0),
        fastp_threads=int(job_data.get("fastp_threads", 0) or 0),
        region=job_data.get("region"),
        depth_cutoff=int(job_data.get("depth_cutoff", 10) or 10),
      )
      if summary.get("quast_txt"):
        report_details["quast_txt"] = summary["quast_txt"]
      if summary.get("quast_html"):
        report_details["quast_html"] = summary["quast_html"]
    except Exception as e:
      print(f"Reference-guided assembly failed: {e}")
      traceback.print_exc(file=sys.stderr)
      report_details["errors"].append("Reference-guided assembly failed.")

    # Retrieve software versions for reference-guided tools
    report_details["tools"] = {
      "bwa": get_software_version("bwa"),
      "samtools": get_software_version("samtools"),
      "bcftools": get_software_version("bcftools"),
      "fastp": get_software_version("fastp"),
      "quast": get_software_version("quast.py"),
    }

  else:
    # IRMA-based assembly
    assembly_output_dir = os.path.join(output_dir, strategy)
    lowqc_output_dir = os.path.join(output_dir, f"{strategy}_lowQC")

    if paired_end_lib:
      read1 = paired_end_lib.get("read1")
      read2 = paired_end_lib.get("read2")
      if not read1 or not read2:
        print("Error: Missing reads for paired-end library.")
        sys.exit(-1)
      local_read1 = os.path.join(output_dir, "read1.fastq")
      local_read2 = os.path.join(output_dir, "read2.fastq")
      if not (fetch_file_from_ws(read1, local_read1) and fetch_file_from_ws(read2, local_read2)):
        print("Error: Failed to fetch paired-end reads.")
        sys.exit(-1)
      success, irma_output = run_irma(module, local_read1, local_read2, output_dir=assembly_output_dir)
      if not success and module == "FLU":
        if "found no QC'd data" in irma_output:
          report_details["notes"].append(
            f"FLU module found no QC'd data. Assembly was re-run with the {FLU_FALLBACK_MODULE} module.")
        print(f"FLU module produced no assembly output. Retrying with {FLU_FALLBACK_MODULE} module...")
        run_irma(FLU_FALLBACK_MODULE, local_read1, local_read2, output_dir=lowqc_output_dir)
        assembly_output_dir = lowqc_output_dir

    elif single_end_lib:
      read = single_end_lib.get("read")
      if not read:
        print("Error: Missing read for single-end library.")
        sys.exit(-1)
      local_read = os.path.join(output_dir, "read.fastq")
      if not fetch_file_from_ws(read, local_read):
        print("Error: Failed to fetch single-end read.")
        sys.exit(-1)
      success, irma_output = run_irma(module, local_read, output_dir=assembly_output_dir)
      if not success and module == "FLU":
        if "found no QC'd data" in irma_output:
          report_details["notes"].append(
            f"FLU module found no QC'd data. Assembly was re-run with the {FLU_FALLBACK_MODULE} module.")
        print(f"FLU module produced no assembly output. Retrying with {FLU_FALLBACK_MODULE} module...")
        run_irma(FLU_FALLBACK_MODULE, local_read, output_dir=lowqc_output_dir)
        assembly_output_dir = lowqc_output_dir

    elif sra_id:
      r1, r2 = fetch_fastqs_from_sra(sra_id, output_dir=output_dir)
      if not r1:
        print("Error: Failed to fetch FASTQs for SRA run accession.")
        sys.exit(-1)
      success, irma_output = run_irma(module, r1, r2, output_dir=assembly_output_dir)
      if not success and module == "FLU":
        if "found no QC'd data" in irma_output:
          report_details["notes"].append(
            f"FLU module found no QC'd data. Assembly was re-run with the {FLU_FALLBACK_MODULE} module.")
        print(f"FLU module produced no assembly output. Retrying with {FLU_FALLBACK_MODULE} module...")
        run_irma(FLU_FALLBACK_MODULE, r1, r2, output_dir=lowqc_output_dir)
        assembly_output_dir = lowqc_output_dir

    fasta_file = os.path.join(output_dir, f"{job_data.get('output_file', 'sample')}_all.fasta")
    is_concatenated = concatenate_fasta_files(assembly_output_dir, fasta_file)
    if is_concatenated:
      if os.path.getsize(fasta_file) == 0:
        print("FASTA file(s) generated by IRMA is empty.")
        report_details["errors"].append("FASTA file(s) generated by IRMA is empty.")
      else:
        quast_dir = os.path.join(output_dir, "quast")
        try:
          run_quast(quast_dir, fasta_file)
          report_details["quast_txt"] = os.path.join(quast_dir, "report.txt")
          report_details["quast_html"] = "quast/report.html"
        except Exception:
          report_details["errors"].append("QUAST failed to run.")
    else:
      report_details["errors"].append("FASTA file was not generated by IRMA.")

    # Move individual segment files to the main output folder
    move_fasta_files(assembly_output_dir, output_dir)

    # Retrieve software versions
    report_details["tools"] = {
      "IRMA": get_software_version("IRMA"),
      "quast": get_software_version("quast.py")
    }

  # Generate HTML report (common to both strategies)
  html_report = generate_html_report(report_details)
  with open(os.path.join(output_dir, "AssemblyReport.html"), "w") as output_file:
    output_file.write(html_report)
