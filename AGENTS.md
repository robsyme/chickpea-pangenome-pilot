# Chickpea pan-genome reassembly pilot

## What this project is

A pilot to test whether we can produce better, more honest genome assemblies for
Australian chickpea (*Cicer arietinum*) cultivars than those published in:

> Garg et al. (2025) "An Australian chickpea pan-genome provides insights into
> genome organization." *Plant Biotechnology Journal* 23:3967–3983.
> https://doi.org/10.1111/pbi.70192

The paper builds a 15-cultivar pan-genome and structural-variant catalogue from
**stLFR linked-read data** (MGI single-tube long-fragment read, co-barcoded short
reads on a DNBSEQ-T7). This is barcoded short-read data, not true long reads —
there is no PacBio HiFi or Oxford Nanopore in the study.

## The critique we want to test

The published "chromosome-scale" assemblies (scaffold N50 ~58 Mb) get their
contiguity almost entirely from **RagTag reference-guided scaffolding onto CDC
Frontier v2.0**. The actual de novo contiguity is the contig N50: 17–34 kb
(avg ~23 kb).

Because this is a pan-genome / SV study, that creates a circularity problem:
contigs are forced onto the CDC Frontier coordinate system, then SVs are called
against CDC Frontier. Reference-guided scaffolding systematically flattens the
large rearrangements (inversions, translocations, large indels) that diverge from
the reference — exactly the SVs the paper claims to discover. Fragmented contigs
(~23 kb N50) also split genes, which can show up as false "dispensable"/absent
gene calls in the core/dispensable analysis.

## What we can and cannot deliver

Can deliver (defensible):
- Reference-free / barcode-scaffolded assemblies that do not inherit CDC Frontier's
  structure — scientifically better for SV discovery even at similar contig N50.
- Transparent QC the original lacks: k-mer QV and completeness from the reads
  (Merqury-style), barcode-based misassembly detection, BUSCO on the de novo
  (not reference-guided) assembly.
- Quantified demonstration of the artefacts above (RagTag chimeric joins,
  suppressed SVs, fragmentation-driven false gene absences).

Cannot deliver without new data:
- HiFi/ONT-grade contig contiguity. Chickpea is ~45% repeats (mostly LTR
  retrotransposons). Barcodes scaffold across repeats but do not resolve repeat
  sequence longer than the ~100 bp reads. The short+linked-read repeat ceiling is
  hard; reaching contiguous, reference-free contigs needs new long-read sequencing.

So the honest deliverable is "methodologically better and more honest assemblies,"
not "more contiguous genomes." That distinction must stay explicit in any write-up.

## Approach: Nextflow as an executable lab notebook

The work is recorded as a Nextflow pipeline so every step is reproducible and the
repo itself is the record of what was done. Runs are submitted to Seqera Platform.

Pilot scope is a **single cultivar first** — CBA Captain (SRR32381426) — to
de-risk before scaling to all 15. Planned stages:

1. **Gate zero — barcode integrity.** Pull SRR32381426 and confirm the stLFR
   barcode is intact in read2. Everything downstream depends on this. The ~242 bp
   combined paired-read length suggests read2 still carries the barcode block, but
   it must be verified on real data before committing.
2. **Baseline.** Reproduce the authors' route (`stlfr2supernova` → Supernova
   v2.1.1) to get a comparable baseline. stLFR coverage here is ~220×, which
   over-saturates linked-read assemblers (Supernova targets ~56×), so subsample.
3. **Reference-free assembly + non-reference scaffolding.** Barcode-aware
   assembly/scaffolding (e.g. ARCS/ARKS+LINKS, SLR-superscaffolder, cloudSPAdes)
   that does not use CDC Frontier for ordering.
4. **QC + artefact quantification.** BUSCO, k-mer QV, misassembly detection;
   compare against the published assemblies and demonstrate the artefacts.

Decision point after the pilot: if clear artefacts and recovered SVs are
demonstrated, scale to all 15 cultivars; if the pilot mostly reproduces the
published numbers, stop having spent one genome's worth of effort.

## Data

- NCBI BioProject **PRJNA1225167** — 15 stLFR WGS runs (DNBSEQ-T7), one per
  cultivar, ~88–109 Gbp each. Submitter: Murdoch University.
- Pilot cultivar: **CBA Captain**, run **SRR32381426**, BioSample SAMN46879420.
- Published assemblies + annotations: Figshare
  https://doi.org/10.6084/m9.figshare.28632494.v1
- Reference used by the paper: CDC Frontier v2.0.

## Infrastructure

Seqera Platform (Seqera Cloud, https://api.cloud.seqera.io), interacted with via
the `tw` CLI (v0.25.0). `TOWER_ACCESS_TOKEN` is set in the environment.

| Item | Value |
|------|-------|
| Organization | `robsyme-research` (ID 151253561315404) |
| Workspace | `chickpea-pangenome-pilot` (ID 66005590816034), PRIVATE |
| Credentials | `aws-scidev-playground` (AWS) |
| Region | `us-east-1` |
| Work dir | `s3://scidev-playground-us-east-1/robsyme/chickpea-pilot/work` |
| Default CE | `aws-batch-default` (ID `5x6WMlOJHt9uwO0SzXcNBm`), AWS Batch Forge — primary |

The default compute environment is AWS Batch (Forge) on Spot (capacity-optimized),
instance families c6id/m6id/r6id (all with local NVMe), with Fusion v2 + Wave,
NVMe fast instance storage, and Fusion Snapshots (spot-interruption recovery).
EBS auto-scale is off (Fusion + NVMe replace it). Max 1000 CPUs. The r6id family
covers the memory-hungry assembly steps. An earlier `aws-cloud` CE was used only to
test credentials with `nextflow-io/hello` and was then deleted.

Recipe to recreate the default CE:

```bash
export TOWER_WORKSPACE_ID=66005590816034
tw compute-envs add aws-batch forge \
  --name aws-batch-default \
  --credentials aws-scidev-playground \
  --region us-east-1 \
  --work-dir s3://scidev-playground-us-east-1/robsyme/chickpea-pilot/work \
  --allow-buckets s3://scidev-playground-us-east-1 \
  --max-cpus 1000 \
  --provisioning-model SPOT \
  --alloc-strategy SPOT_CAPACITY_OPTIMIZED \
  --instance-types c6id,m6id,r6id \
  --fusion-v2 --wave --fast-storage --snapshots --no-ebs-auto-scale \
  --wait AVAILABLE
tw compute-envs primary set --name aws-batch-default
```

Launching a pipeline to the CE:

```bash
export TOWER_WORKSPACE_ID=66005590816034
tw launch <pipeline-url-or-name> --compute-env aws-cloud-us-east-1 --wait SUBMITTED
```

## Conventions

- All compute runs go through Seqera Platform on the workspace above; the `tw` CLI
  is the interface.
- Set `TOWER_WORKSPACE_ID=66005590816034` before `tw` commands that operate on the
  workspace.
- The Nextflow pipeline in this repo is the canonical record of the analysis — keep
  it runnable and keep parameter/version choices in code, not just prose.

## Pipeline prerequisites and conventions

- **Nextflow ≥ 26.04** (developed against 26.04.2). Required for the syntax below.
- **Strict v2 parser is mandatory for the typed syntax.** Set the environment variable
  `NXF_SYNTAX_PARSER=v2` at launch. This is read before any config is parsed, so it
  **cannot** live in `nextflow.config` — set it as a compute-environment env var on
  Seqera Platform (or prefix the command locally). Without it, typed processes fail to
  compile.
- **Enable the type system per script.** Put `nextflow.enable.types = true` at the top
  of every script (`main.nf` and each module) that uses typed processes/workflows or
  records. It is a feature-flag declaration, not a config setting. The `params {}` and
  `output {}` blocks work without the flag, but typed processes/records need it.
- **Validated typed syntax (Nextflow 26.04, preview — may change):**
  - Record: `record SraRun { accession: String }`; construct with
    `new SraRun(accession: 'SRR...')` (named args, not positional).
  - Typed process inputs: `name: Type` (e.g. `run: SraRun`, `reads: Path`,
    `n_reads: Integer`). A value channel input broadcasts across a queue channel.
  - Typed process outputs: `name: Type = <fn>` where `<fn>` is `file('...')`,
    `files('...')`, `stdout()`, `env('VAR')`, `record(...)`, or `tuple(...)`.
    Single unnamed output may omit the name.
  - Workflow: `main:` then `publish:` with `name = channel` (assignment form; `>>` is
    not valid here).
  - Output block: `output { name: Channel<T> { path '...' } }`. A `Value<T>` may be
    written as `T`.
  - Records can be shared across module files via `include { SraRun } from './types.nf'`.
- **Workflow output syntax (new).** Use the top-level `output { ... }` block plus the
  workflow `publish:` section. Do not use `publishDir`.
- **Multi-accession by design.** The pipeline takes a list of run accessions and
  parallelizes per-accession (one accession per channel element); a single-accession
  run is just the N=1 case. First pilot run uses only `SRR32381426`, but the structure
  must not assume one input.

## Gate-zero result (CBA Captain / SRR32381426)

Validated locally on a 50,000-read subsample pulled from ENA (2026-06-10):

- Read structure: R1 = 100 bp, R2 = 142 bp (insert 100 + 42 bp barcode block).
- stLFR barcodes are **intact in read2**, sitting at offsets **100 / 116 / 132**
  (layout 10 + 6-spacer + 10 + 6-spacer + 10). Note the 6 bp spacer before the third
  barcode — this differs from the `split_barcode.pl` default (`n4=0`, barcode at 126),
  so the check **detects the offsets empirically** rather than assuming them.
- **87.6%** of reads carry a valid 3-barcode combination (1-mismatch tolerance);
  43,593 distinct barcodes in the sample. Gate **PASS** (threshold 0.80).
- Conclusion: the data survived the SRA/ENA round-trip with barcodes usable as
  linked-reads. Gate zero is green for the pilot to proceed to the baseline stage.

The pipeline implements this check generically (the `barcode_check` module +
its module-scoped binary at
`modules/barcode_check/resources/usr/bin/stlfr_barcode_check.py`, unit-tested in
`tests/`). Module binaries require `nextflow.enable.moduleBinaries = true` in the
pipeline script and Wave when the work directory is on S3 (both in place).

**Formal Platform run (full 2M-read subsample):** run `2sSrLdl93wx6nA`
(`gate-zero-SRR32381426-fullsubsample`) on `aws-batch-default`, 2026-06-11,
**SUCCEEDED** — gate **PASS**. Published results:

- `valid_fraction` = **0.903** (1,805,844 / 2,000,000 reads carry a valid barcode).
- **1,657,967 distinct barcodes**, with a realistic reads-per-barcode tail (1.52M
  singletons, 120k doubletons, down to a max of 14) — the genuine co-barcoding
  signature, not random matching.
- seqkit: R1 100 bp, R2 142 bp; Q20 97.1% / 95.7%; GC 33.2% / 39.6% (R1/R2).
- Barcode offsets 100/116/132 detected on the full sample, consistent with the
  50k local run. BARCODE_CHECK ran in ~12 s on a c6id.large spot instance,
  peak RSS ~574 MB.

## Reference: local Nextflow source

Authoritative Nextflow source (all version worktrees) is checked out at
`~/dev/bears/nextflow` — match the runtime version dir (e.g. `v26.04.3`). The typed
syntax docs are under `docs/process-typed.md`, `docs/workflow-typed.md`, and
`docs/reference/syntax.md`; the grammar is `modules/nf-lang/src/main/antlr/ScriptParser.g4`.
Validate syntax empirically with `NXF_SYNTAX_PARSER=v2 nextflow run <file> -stub`.
