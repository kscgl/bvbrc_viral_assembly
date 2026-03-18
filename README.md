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

- **Obtains input reads** from exactly one source:
  - **Paired-end reads**: `paired_end_lib.read1` + `paired_end_lib.read2`
  - **Single-end reads**: `single_end_lib.read`
  - **SRA**: `srr_id` (downloads reads via `p3-sra`)
- **Resolves the reference** and runs the pipeline in `scripts/reference_guided_assembly.py`:
  - `fastp` trimming
  - `bwa mem` alignment
  - `samtools` sort/index + depth
  - `bcftools mpileup/call` variant calling
  - `bcftools consensus` consensus FASTA (with low-coverage masking)
  - `quast.py` report
- **Generates** `AssemblyReport.html` using the template under `lib/viral_assembly_report_template.html`.

### Environment / paths

- **`KB_TOP`**: In the BV-BRC runtime, `KB_TOP` points at the module root. For local/dev runs, the script now falls back to the repo root automatically if `KB_TOP` is unset.
- **Workspace paths (`ws:`)**: In BV-BRC production, `paired_end_lib.read1/read2` and `single_end_lib.read` are typically workspace paths and are fetched with `p3-cp`.
  - Local paths (like `old_results/SRR..._1.fastq`) are **not** valid for the `p3-cp ws:` fetch step.
  - For local testing of reference-guided assembly, call `run_reference_guided(...)` directly (see below), or use the local-compare runner in `test_files/compare_reference_guided.py`.

## UI → service JSON contract (reference-guided)

At minimum, the UI must provide:

- `strategy`: `"reference_guided"`
- `reference_type`: `"fasta"` or `"genbank"`
- A reference selector:
  - FASTA: `fasta_file` or `reference_assembly`
  - GenBank: `genbank_accession` or `reference_assembly`
- Exactly one read source:
  - `paired_end_lib` **OR** `single_end_lib` **OR** `srr_id`
- A sample/output name:
  - `output_file` (preferred) or `srr_id`

### Supported optional keys (reference-guided)

- `email`: required for GenBank fetching if `ENTREZ_EMAIL` isn’t set.
- `align_threads`, `fastp_threads`: tool thread knobs.
- `region`: restricts `bcftools mpileup` to a contig/region.
- `depth_cutoff`: integer depth cutoff for masking low-coverage bases with `N` (default `10`, `0` disables).

### Segmented references / multi-contig consensus

The reference-guided implementation supports segmented references in two ways:

- **Single multi-record FASTA**: `fasta_file` points at a FASTA with multiple records.
  - If `segment_names` is provided, the FASTA headers are rewritten to those names.
  - By default this emits a combined multi-consensus FASTA; per-segment emission can be controlled via `per_segment_consensus`.
- **Multiple references**: provide `reference_fastas` (list) or `references` (semicolon-separated string) where each token is either a local FASTA path or a GenBank accession. These are concatenated into one multi-FASTA reference.

Optional segmented keys:

- `segment_names`: list or semicolon-separated string. Must match the number of FASTA records after reference resolution.
- `per_segment_consensus`: boolean. If the input is a single multi-record FASTA, default behavior matches the original pipeline (combined-only). Set `true` to also emit per-segment FASTAs.

### JSON examples

**FASTA (paired-end, segmented multi-record FASTA):**

```json
{
  "strategy": "reference_guided",
  "reference_type": "fasta",
  "fasta_file": "test_files/Influenza_Genome.fasta",
  "reference_assembly": "test_files/Influenza_Genome.fasta",
  "output_file": "SRR28752452",
  "paired_end_lib": { "read1": "ws:/path/read1.fastq", "read2": "ws:/path/read2.fastq" },
  "segment_names": ["A_HA","A_NA","A_PB1","A_MP","A_NP","A_NS","A_PA","A_PB2"],
  "per_segment_consensus": true
}
```

**GenBank (SRA reads):**

```json
{
  "strategy": "reference_guided",
  "reference_type": "genbank",
  "genbank_accession": "NC_001474.2",
  "reference_assembly": "NC_001474.2",
  "output_file": "SRR27422853",
  "srr_id": "SRR27422853",
  "email": "you@org.org"
}
```

**Multiple references (segmented tokens):**

```json
{
  "strategy": "reference_guided",
  "reference_type": "fasta",
  "references": "ref/seg1.fasta;ref/seg2.fasta;ref/seg3.fasta",
  "segment_names": "seg1;seg2;seg3",
  "output_file": "sampleX",
  "single_end_lib": { "read": "ws:/path/reads.fastq" }
}
```

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
