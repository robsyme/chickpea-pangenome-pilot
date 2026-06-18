# 2026-06-10 — Stage 0: gate-zero pipeline

## Goal of this phase

Before committing to reassembly, confirm the deposited stLFR data is usable as
linked-reads — i.e. the per-molecule barcodes are intact in read2. Build the result
as a Nextflow pipeline (the executable lab notebook) run on Seqera Platform.

## What we did

- Identified the data: NCBI BioProject PRJNA1225167, 15 stLFR (MGI DNBSEQ-T7) WGS
  runs, one per cultivar. Pilot cultivar: **CBA Captain = SRR32381426**.
- Built a typed (DSL2 static-types) pipeline:
  - `FETCH_READS` — resolve ENA `fastq_ftp`, stream-subsample the first N read pairs
    (no full download). Hardened later with output read-count validation +
    `curl --fail --retry` + a `fastq-dl` ENA→SRA fallback (using the prebuilt
    `fastq-dl` container, not a conda re-solve, which had pulled a broken `rich`).
  - `READ_STRUCTURE` — `seqkit stats` on R1/R2.
  - `BARCODE_CHECK` — `bin/stlfr_barcode_check.py`: detects the three 10 bp barcode
    offsets in read2 **empirically** (whitelist hit-rate scan), then reports valid
    fraction, distinct-barcode count, reads-per-barcode histogram. Emits a MultiQC
    custom-content row.
  - `MULTIQC` — aggregates seqkit stats + the barcode table into one report.
- Infrastructure: org `robsyme-research` / workspace `chickpea-pangenome-pilot`,
  AWS Batch + Fusion + Wave CE `aws-batch-default`, results published to
  `s3://scidev-playground-us-east-1/robsyme/chickpea-pilot/results`. Launchpad
  pipeline `chickpea-gate-zero` (v2 parser toggle on; `dumpHashes` + `outputDir` in
  config).

## Key findings

- Read structure: R1 = 100 bp, R2 = 142 bp (100 bp insert + 42 bp barcode block).
- **Barcode offsets are 100 / 116 / 132** in read2 (layout 10 + 6-spacer + 10 +
  6-spacer + 10). This differs from BGI `split_barcode.pl`'s default (third barcode
  at 126), so the check must detect offsets rather than assume them — a hardcoded
  offset gave a false 6.8% before this was caught.
- **Result: gate PASS.** Full 2M-read subsample: **90.3% of reads carry a valid
  barcode** (1,805,844 / 2,000,000), **1,657,967 distinct barcodes** with a realistic
  reads-per-barcode tail. The data is genuinely usable as linked-reads.

## Conventions established

- Compute via Seqera Platform `tw` CLI; `TOWER_WORKSPACE_ID=66005590816034`.
- Typed syntax requires `NXF_SYNTAX_PARSER=v2` (env var, can't go in config; also a
  per-pipeline UI toggle on the Launchpad).
- `outputDir` must be set (S3) or the `output{}` results are lost to the ephemeral
  head-job dir.

## Outcome / next

Stage 0 green end-to-end. Next: stage 1 baseline assembly (`stlfr2supernova` →
Supernova), which needs full-coverage reads (subsample ~220× down to ~56×).
