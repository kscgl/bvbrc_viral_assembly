# BV-BRC Module

## Overview

This repository holds a BV-BRC module.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

## Viral assembly entrypoint (`scripts/run_viral_assembly.py`)

`scripts/run_viral_assembly.py` is the service entrypoint that takes a **single job JSON** and runs one of:

- **Reference-guided assembly** (strategy `reference_guided`) via `scripts/reference_guided_assembly.py`
- **IRMA-based assembly** (other strategies; existing behavior)

### How it runs reference-guided assembly

For reference-guided jobs, the script:

- **Obtains input reads** from exactly one source (staging uses `p3-cp` for workspace paths in production; local paths are supported for dev):
  - **Paired-end reads**: `paired_end_lib.read1` + `paired_end_lib.read2`
  - **Single-end reads**: `single_end_lib.read`
  - **SRA**: `sra_id` (run accession, e.g. `SRR28752452`) — **`p3-sra`** runs in **`scripts/run_viral_assembly.py`** by default; assembly receives local FASTQs only.
- **Resolves reference FASTA** based on `reference_type` and runs the science pipeline in `scripts/reference_guided_assembly.py`:
  - `reference_type: genome` / `auto` — fetches the genome sequence from BV-BRC using **`p3-genome-fasta <genome_id>`** (e.g. `11060.9352`); result is cached as `bvbrc_<id>.fasta` under `output_files/ref_cache/`.
  - `reference_type: fasta` — stages the workspace FASTA into `reference_inputs/` via `p3-cp`.
  - `reference_type: genbank` — fetches from NCBI Entrez (requires `email` or `ENTREZ_EMAIL`); result is cached as `<accession>.fasta` under `output_files/ref_cache/`.
  - `fastp` trimming
  - `bwa mem` alignment
  - `samtools` sort/index + depth
  - `bcftools mpileup/call` variant calling
  - `bcftools consensus` consensus FASTA (with low-coverage masking)
  - `quast.py` report
- **Generates** `AssemblyReport.html` using the template under `lib/viral_assembly_report_template.html`.

### Environment / paths

- **`KB_TOP`**: In the BV-BRC runtime, `KB_TOP` points at the module root. For local/dev runs, the script now falls back to the repo root automatically if `KB_TOP` is unset.
- **Workspace paths**: In BV-BRC production, read paths and `reference_fasta_file` are workspace object paths and are copied with `p3-cp` (`fetch_file_from_ws` in `run_viral_assembly.py`).
  - **Dev convenience**: if the path is an existing absolute or repo-relative file, it is copied without `p3-cp`.
  - For parity testing without the full service, use `run_reference_guided(...)` with the same job keys as production, or `test_files/compare_reference_guided.py`.

### Local setup and run (quickstart)

Required tools in your env:

- `python3`, `fastp`, `bwa`, `samtools`, `bcftools`, `quast.py`
- for SRA jobs in BV-BRC: `p3-sra` on `PATH`
- for `genome`/`auto` reference jobs in BV-BRC: `p3-genome-fasta` on `PATH`
- Python package: `biopython`

Recommended local invocation:

```bash
/home/<user>/miniconda3/bin/conda run -n bvbrc python3 scripts/run_viral_assembly.py -j <job.json> -o <out_dir>
```

Example:

```bash
/home/<user>/miniconda3/bin/conda run -n bvbrc python3 scripts/run_viral_assembly.py \
  -j test_files/job_reference_guided_genbank_nc001474.json \
  -o runs/local_one_job
```

## UI → service JSON contract (reference-guided)

The form definition in `app_specs/ViralAssembly.json` matches the table below. Older reference keys (`fasta_file`, `genbank_accession`, `reference_assembly`, `reference_input`) are **not** read by the runner. For SRA-only reads, use **`sra_id`**; legacy **`srr_id`** is still accepted.

### Required inputs (every reference-guided job)

| Input | Value |
|-------|--------|
| `strategy` | `"reference_guided"` |
| `reference_type` | `"genome"`, `"auto"`, `"genbank"`, or `"fasta"` |
| Reference sequence | If `reference_type` is **`genome`** or **`auto`**: set **`reference_genome_id`** (BV-BRC genome ID, e.g. `"11060.9352"`). If **`genbank`**: set **`reference_genbank_accession`** (one accession, `"ACC1;ACC2"` for multiple, or a JSON list). If **`fasta`**: set **`reference_fasta_file`** (one workspace/local path, `"path1;path2"` or a list for multiple files). |
| Reads (pick **one**) | **`paired_end_lib`**: `read1` and `read2`, **or** **`single_end_lib`**: `read`, **or** **`sra_id`** (run accession such as `SRR28752452`). |
| Output basename | **`output_file`** (recommended). If missing, some jobs fall back to `sra_id`. |

**`genome` / `auto`**: the genome sequence is fetched from BV-BRC via `p3-genome-fasta` — no `email` needed.  
**`genbank`-only requirement:** set **`email`** (or environment variable `ENTREZ_EMAIL`) so NCBI sequences can be fetched for `reference_genbank_accession`.

### Optional inputs (tuning and segmented references)

| Input | Purpose |
|-------|--------|
| `segment_names` | For multi-record or multi-token references: list or `"a;b;..."` string; must match the number of reference contigs after resolution. Renames/relabels FASTA headers before alignment when used with a single multi-record FASTA. |
| `per_segment_consensus` | For a **single** multi-record FASTA: if `true`, also write one consensus FASTA per segment (default is combined-only for that case). |
| `align_threads`, `fastp_threads` | Thread counts for `bwa` and `fastp`. |
| `region` | Limit `bcftools mpileup` to one contig/region. |
| `depth_cutoff` | Depth below which consensus is masked with `N` (default `10`; `0` disables). |
| `download_sra_from_prefetch`, `download_fastqs_from_sra` | If either is explicitly **`false`**, skip **`p3-sra`** and require FASTQs already under the output directory (defaults: both **`true`** when `sra_id` is set). |

**Segmented / multi-reference behavior (concepts):**

- **One multi-record FASTA** (`reference_fasta_file` → one file, many `>` records): optional `segment_names` + `per_segment_consensus` control header labels and per-segment outputs.
- **Several separate references**: put multiple accessions in `reference_genbank_accession` or multiple paths in `reference_fasta_file` (list or semicolon-separated); they are concatenated into one alignment reference.

### JSON examples

**1) BV-BRC genome selection + paired-end reads** (`reference_type: genome` — user selects a specific genome):

```json
{
  "strategy": "reference_guided",
  "reference_type": "genome",
  "reference_genome_id": "11060.9694",
  "output_file": "viral_ref_guided_fasta_genome",
  "paired_end_lib": {
    "read1": "ws:/path/to/R1.fastq",
    "read2": "ws:/path/to/R2.fastq"
  }
}
```

**2) Auto-select reference via Minhash + paired-end reads** (`reference_type: auto` — UI finds the closest BV-BRC genome automatically):

```json
{
  "strategy": "reference_guided",
  "reference_type": "auto",
  "reference_genome_id": "11060.9352",
  "output_file": "viral_ref_guided_fasta_auto",
  "paired_end_lib": {
    "read1": "ws:/path/to/R1.fastq",
    "read2": "ws:/path/to/R2.fastq"
  }
}
```

**3) Minimal FASTA + paired-end reads** (only required inputs; production paths are usually workspace objects):

```json
{
  "strategy": "reference_guided",
  "reference_type": "fasta",
  "reference_fasta_file": "ws:/path/to/reference.fasta",
  "output_file": "my_sample",
  "paired_end_lib": {
    "read1": "ws:/path/to/R1.fastq",
    "read2": "ws:/path/to/R2.fastq"
  }
}
```

**4) Minimal GenBank + SRA** (`email` required if `ENTREZ_EMAIL` is not set):

```json
{
  "strategy": "reference_guided",
  "reference_type": "genbank",
  "reference_genbank_accession": "NC_001474.2",
  "output_file": "SRR27422853",
  "sra_id": "SRR27422853",
  "email": "you@org.org"
}
```

**5) FASTA + paired-end + optional segmented influenza labels** (`segment_names` and `per_segment_consensus` are **not** required; they only apply to this multi-segment FASTA use case):

```json
{
  "strategy": "reference_guided",
  "reference_type": "fasta",
  "reference_fasta_file": "test_files/Influenza_Genome.fasta",
  "output_file": "SRR28752452",
  "paired_end_lib": { "read1": "ws:/path/read1.fastq", "read2": "ws:/path/read2.fastq" },
  "segment_names": ["A_HA","A_NA","A_PB1","A_MP","A_NP","A_NS","A_PA","A_PB2"],
  "per_segment_consensus": true
}
```

**6) Multiple GenBank accessions + single-end reads** (optional `segment_names` must match the number of reference contigs—here, two accessions → two names):

```json
{
  "strategy": "reference_guided",
  "reference_type": "genbank",
  "reference_genbank_accession": "NC_045512.2;MN908947.3",
  "segment_names": "SARS2_ref1;SARS2_ref2",
  "output_file": "sampleX",
  "single_end_lib": { "read": "ws:/path/reads.fastq" },
  "email": "you@org.org"
}
```

### Expected outputs

Top-level output folder contains final human-consumable artifacts; intermediate/tool files are under `output_files/`.

Top-level (typical):

- `AssemblyReport.html`
- `report.html` (QUAST moved up from quast subdir)
- consensus FASTA(s):
  - `<sample>.consensus.<accession_or_segment>.fasta`
  - `<sample>.consensus.multi.fasta`
  - FASTA single-reference alias (when applicable): `<sample>.consensus.<ref_label>.fasta`

`output_files/` (typical):

- staged reads (`reads_1.fastq`, `reads_2.fastq` or `reads.fastq`)
- trimmed reads
- reference index files (`.amb/.ann/.bwt/.pac/.sa/.fai`)
- alignment/variant files (`.sam`, `.sorted.bam`, `.bai`, `.vcf.gz`, `.csi`)
- low-coverage BED (`*.lowcov_dp<cutoff>.bed`)
- `ref_cache/` with fetched reference FASTAs: `bvbrc_<genome_id>.fasta` for `genome`/`auto` jobs, `<accession>.fasta` for `genbank` jobs

SRA note:

- For `sra_id` jobs, **`p3-sra`** writes `{SRR}_1.fastq` / `{SRR}_2.fastq` (or `{SRR}.fastq`) under the job output directory (same naming as before; the **value** is still an SRR-style run accession).

## Local parity testing and comparison (old vs new)

The repo includes a local comparison runner that:

- runs the **original CSV-driven pipeline** (`test_files/original_reference_guided_assembly.py`)
- runs the **new JSON-driven pipeline** locally (direct `run_reference_guided(...)`, no workspace tooling)
- then compares consensus FASTA plus a “brutal” SHA-256 comparison across multiple artifact patterns.

### Output folders

All local runs should go under `runs/` (ignored by git). Example:

- `runs/old_results/` (old pipeline)
- `runs/new_results/` (new pipeline)

### Run the comparison

From the repo root:

```bash
python test_files/compare_reference_guided.py \
  --csv test_files/sample.csv \
  --email you@org.org \
  --job-json test_files/job_reference_guided_fasta_segmented.json \
  --old-out runs/old_results \
  --new-out runs/new_results
```

If you already ran both and only want to compare:

```bash
python test_files/compare_reference_guided.py \
  --csv test_files/sample.csv \
  --email you@org.org \
  --job-json test_files/job_reference_guided_fasta_segmented.json \
  --old-out runs/old_results \
  --new-out runs/new_results \
  --skip-run
```

For more notes and examples, see `test_files/reference_guided_parity_test.md`.
