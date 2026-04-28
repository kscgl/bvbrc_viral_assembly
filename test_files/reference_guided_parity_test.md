## Reference-guided parity test (original vs JSON-driven)

This doc shows how to run both:
- `test_files/original_reference_guided_assembly.py` (CSV-driven original)
- `scripts/run_viral_assembly.py` + `scripts/reference_guided_assembly.py` (backend JSON-driven)

and compare their consensus FASTA outputs.

### 1) Prepare a JSON job (segmented FASTA example)

Create `test_files/job_reference_guided_fasta_segmented.json` with:

```json
{
  "strategy": "reference_guided",
  "reference_type": "fasta",
  "fasta_file": "test_files/Influenza_Genome.fasta",
  "reference_assembly": "test_files/Influenza_Genome.fasta",
  "output_file": "SRR28752452",
  "srr_id": "SRR28752452",
  "segment_names": ["A_HA", "A_NA", "A_PB1", "A_MP", "A_NP", "A_NS", "A_PA", "A_PB2"]
}
```

Notes:
- `segment_names` is optional; if provided it must match the number of FASTA records.
- If your backend provides multiple reference tokens, you can pass them as `reference_fastas` (list) or `references` (semicolon-separated string). Each token can be a local path or a GenBank accession.

### 2) Run the compare script

From the repo root:

```bash
python test_files/compare_reference_guided.py ^
  --csv test_files/sample.csv ^
  --email you@org.org ^
  --job-json test_files/job_reference_guided_fasta_segmented.json ^
  --old-out runs/old_results ^
  --new-out runs/new_results
```

If you already ran both pipelines and only want to compare:

```bash
python test_files/compare_reference_guided.py ^
  --csv test_files/sample.csv ^
  --email you@org.org ^
  --job-json test_files/job_reference_guided_fasta_segmented.json ^
  --old-out runs/old_results ^
  --new-out runs/new_results ^
  --skip-run
```

### 3) Expected outputs (new JSON-driven)

In `runs/new_results/` you should see:
- `SRR28752452.vcf.gz` (+ index)
- If reference has 1 contig: `SRR28752452.consensus.fasta`
- If reference has >1 contig: `SRR28752452.consensus.multi.fasta`
- If reference has >1 contig and per-segment emission is enabled: `SRR28752452.consensus.<segment>.fasta`
- `quast_reference_guided/report.txt` and `quast_reference_guided/report.html` (when QUAST succeeds)

### 4) JSON knobs added for segmented support

These are optional and do not break existing backend JSON:
- `segment_names`: list or semicolon-separated string
- `reference_fastas`: list of reference tokens (paths or accessions)
- `references`: semicolon-separated string of reference tokens
- `per_segment_consensus`: boolean
  - If the input is a single multi-record FASTA, default behavior matches the original (combined-only).
  - Set `per_segment_consensus: true` to force per-segment FASTAs for a single multi-record FASTA.

