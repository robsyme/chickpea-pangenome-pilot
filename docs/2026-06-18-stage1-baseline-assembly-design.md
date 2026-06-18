# 2026-06-18 — Stage 1 design: baseline Supernova assembly

A design (not yet an attempt log) for the baseline assembly stage. Stage 0
(gate zero) is green; this stage reproduces the authors' route
(`stlfr2supernova` → Supernova v2.1.1) on our own reads to get a same-input
baseline to compare later reference-free assemblies against.

## Goal

Produce a defensible, same-input Supernova baseline for each pilot accession, so
that when Stage 3 (reference-free assembly) lands we can compare on identical
input rather than against the paper's published numbers. The design must allow
**multiple assemblies per accession** so alternative assemblers can be added and
compared later.

## Decisions

- **Same-input Supernova run**, not just anchoring to the paper's reported
  contig N50. Reviewers will expect the standard linked-read assembler run on the
  same reads.
- **Two coverage variants per accession**, modelled as one parameterised process:
  - `220×` — faithful to the paper. The paper ran "default parameters"; Supernova's
    default `--maxreads` (1.2 B) works out to ~243× for a 740 Mb genome, i.e.
    essentially their full ~220× data with no real subsampling. Setting `220`
    explicitly reproduces that without depending on a version-specific default.
  - `56×` — Supernova's recommended raw coverage. The over-saturation at 220× may
    itself degrade the assembly; running both tests that directly and is a
    publishable finding on its own.
- **Extend the existing `main.nf` / `chickpea-gate-zero` Launchpad pipeline**, not
  a separate pipeline. We build iteratively and want to reuse the existing task
  cache. None of the four gate-zero processes are modified, so their hashes (and
  caches) survive; assembly is added as a new downstream subgraph.
- **Bake Supernova into a single private image** (Supernova in a low layer,
  `stlfr2supernova` on top). A private, single-org registry is the same legal
  posture as runtime-fetch-from-our-own-S3 — the binary never leaves infra we
  control — so this is not redistribution to a third party. The 10x EULA accepted
  at download is the actual arbiter; standard commercial terms permit internal
  copies. (Not legal advice.)

## Architecture

Assembly hangs off the same per-accession `runs` channel as gate zero. Nothing in
the gate-zero path changes.

```
runs ─┬─ (gate-zero subgraph, unchanged) ── MULTIQC
      └─ FETCH_FULL_READS ── STLFR_CONVERT ── (combine with coverage_cutoffs) ──
             SUPERNOVA ── ASSEMBLY_STATS
```

Shape that matters: **convert once, assemble N times.** `STLFR_CONVERT`
(split_barcode → SOAPfilter → merge_barcodes) is coverage-independent and
expensive, so it runs once per accession. Its output is crossed with the coverage
list via `combine`, and a single `SUPERNOVA` process fans out over the result:

```groovy
cutoffs    = channel.fromList(params.coverage_cutoffs)   // [220, 56]
full_reads = FETCH_FULL_READS(runs)
converted  = STLFR_CONVERT(full_reads)
assemblies = SUPERNOVA(converted.combine(cutoffs))       // one task per (accession × cutoff)
stats      = ASSEMBLY_STATS(assemblies)
```

`FETCH_FULL_READS` is a **new** process, not a modification of the gate-zero
`FETCH_READS` (which only pulls the 2M-read subsample). Keeping it separate
preserves the gate-zero hash. It reuses the proven ENA → `fastq-dl` fallback
logic to pull the full run (~88–109 Gbp).

`mkoutput --style=pseudohap` runs as a tail step inside the `SUPERNOVA` process
(one scaffold per record — what the paper fed to RagTag), so the large Supernova
output dir on scratch is never re-staged.

Assembly runs for every accession in `params.accessions`; the gate-zero report
stays informational rather than auto-gating (both pilot accessions already pass).

## Data model

New record in `types.nf`:

```groovy
record Assembly {
    accession: String
    assembler: String     // 'supernova' (open for cloudSPAdes etc. later)
    coverage:  Integer    // 220, 56, ...
    fasta:     Path       // pseudohap scaffolds, gzipped
}
```

`assembler` + `coverage` travel inside the record through `ASSEMBLY_STATS` and
into both the published path and the MultiQC table label. The coverage is
therefore threaded end to end by construction, so the two variants of one
accession can never collide or be mislabelled.

## Coverage control & parameters

```groovy
params {
    genome_size: Integer = 740_000_000          // Cicer arietinum, ~740 Mb (confirm vs CDC Frontier v2.0)
    coverage_cutoffs: List<Integer> = [220, 56]
    supernova_localcores: Integer = 16
    supernova_localmem: Integer = 250            // GB
}
```

`SUPERNOVA` computes `maxreads = genome_size × coverage / read_len`, with
`read_len` measured from the converted FASTQ (the stLFR insert is ~100 bp, but
read rather than assumed). It logs target coverage and the realised `--maxreads`,
and we cross-check against the coverage Supernova reports in its own log.

Two subtleties to capture in code comments:

- The `stlfr2supernova` merge step collapses every 8 stLFR barcodes into one of
  10x's 4.79 M barcodes and silently discards overflow. So the converted pool is
  already capped before `--maxreads` applies. For `220×`, `--maxreads` likely
  exceeds the available reads (uses all — faithful); for `56×` it subsamples down.
  This asymmetry is expected and matches how the paper ran.
- `genome_size` is a param, not a constant, so the same pipeline serves the other
  13 cultivars and is portable.

## Outputs

`ASSEMBLY_STATS` runs `seqkit stats -a` (already used in gate zero) and emits two
files, mirroring `BARCODE_CHECK`: a machine-readable `stats.tsv` and a `_mqc`
custom-content table (assembler, coverage, scaffold N50, total length, scaffold
count, %N).

Published layout, keyed by accession then strategy:

```
results/assemblies/<accession>/<assembler>-<coverage>x/
    assembly.scaffolds.fa.gz
    supernova_summary.csv      # Supernova's own report, for provenance
    stats.tsv
```

Stage 1 stops at contiguity stats. BUSCO, Merqury k-mer QV, and misassembly
detection are Stage 4 — adding them now would bloat the per-task hash before the
expensive assembly is cached.

**MultiQC merge (last step):** the `_mqc` files are produced now but wired into
the existing `MULTIQC` collect only once the stats fields are settled. Adding them
changes MULTIQC's hash (it re-runs), but it is cheap.

## Resources & Fusion

Per-process resources via labels (infra config lives on the CE, not the repo):

- `SUPERNOVA`: 16 cores / 250 GB → r6id.8xlarge-class (32 vCPU, 256 GB, NVMe). The
  memory-bound step.
- `STLFR_CONVERT`: ~16 cores / 64 GB, CPU-bound over the full ~90 Gbp.
- `FETCH_FULL_READS`: network-bound, modest CPU/RAM, tens of GB of output.

No `disk` directive — it is not honoured on AWS Batch. Fusion chunks and streams
inputs and outputs against S3, so the large converted FASTQs and assembly outputs
are not bounded by local disk. The r6id family is still the right choice: its
local NVMe backs Fusion's cache and gives fast scratch for Supernova's
memory-mapped, random-access temp I/O. Whether Supernova's heavy internal scratch
performs acceptably streamed through Fusion is the one thing to watch on the smoke
run (assemblers can be sensitive to non-local temp).

**Spot reclamation is the real risk.** A ~1-day Supernova task on Spot can be
reclaimed near the end and, with `maxRetries = 2`, restart from zero. Mitigations:
the CE's Fusion Snapshots (spot-interruption recovery), or point the `SUPERNOVA`
step at on-demand capacity while the cheap steps stay on Spot. This is a CE/queue
decision at launch, not baked into the repo.

## Validation ladder

1. `-stub` all new processes first to prove the typed DAG compiles cheaply.
2. Smoke run on one accession at a tiny coverage to exercise convert → Supernova →
   mkoutput → stats end to end before any full run.
3. First real run: CBA Captain, `56×` variant only (smaller/faster than `220×`) as
   a cost/time probe, then fan out to both accessions × both cutoffs.
4. Relaunch with `--commit-id $(git rev-parse HEAD)` throughout so gate zero,
   fetch, and convert stay cached while iterating on the assembly tail.

## One-time setup prerequisites (human, before first real run)

- Download `supernova-2.1.1.tar.gz` from the 10x portal (registration + EULA).
- Build and push the private image (Supernova baked low, `stlfr2supernova` on
  top); pin by digest in the module/config.
- Confirm the genome-size param against CDC Frontier v2.0 assembly length.
- Verify the CE can place an r6id.8xlarge and decide Spot vs on-demand for the
  assembly step.

## Open items / risks

- Supernova internal scratch performance under Fusion (confirm on smoke run).
- Effective coverage into Supernova is governed by both the barcode cap and
  `--maxreads`; the realised number must be read from Supernova's log, not assumed.
- Spot reclamation on the ~1-day assembly (mitigations above).
- `read_len` after conversion is the one empirical value feeding the maxreads
  arithmetic — verify on the first run.
