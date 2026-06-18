# 2026-06-18 — Stage 1 built: baseline Supernova assembly (stub-validated)

## Goal of this phase

Extend the gate-zero pipeline with the Stage 1 baseline: reproduce the authors'
route (`stlfr2supernova` → Supernova v2.1.1) on our own reads to get a same-input
assembly to compare later reference-free work against. Design recorded in
`2026-06-18-stage1-baseline-assembly-design.md`; the task-by-task plan is in
`plans/2026-06-18-stage1-baseline-assembly.md`.

## What we did

Built the assembly subgraph as an extension of the existing `main.nf`, leaving the
four gate-zero processes untouched so their task cache survives on relaunch.

- `FETCH_FULL_READS` (`modules/fetch_full_reads.nf`) — full-run download via
  `fastq-dl` (auto ENA→SRA), separate from the subsampling `FETCH_READS` so the
  gate-zero hash is unchanged. Guards against missing paired output and verifies the
  downloads are valid gzip before the expensive stage.
- `STLFR_CONVERT` (`modules/stlfr_convert.nf`) — runs the stLFR→Chromium conversion
  once per accession (coverage-independent), emitting a `ConvertedReads` record.
- `SUPERNOVA` (`modules/supernova.nf`) — one parameterised process fanned out over
  `(accession × coverage)` by `converted.combine(cutoffs)`. It sizes `--maxreads`
  from `genome_size × coverage / read_len` (read length measured from the converted
  reads, not assumed), runs `supernova run` + `mkoutput --style=pseudohap`, and emits
  an `Assembly` record tagged with `assembler` + `coverage`.
- `ASSEMBLY_STATS` (`modules/assembly_stats/main.nf` + the module-scoped
  `assembly_stats.py`) — `seqkit stats -a` per assembly, emitting a stats TSV
  (`AssemblyStats` record) and a MultiQC custom-content table. The summariser is
  unit-tested (`tests/test_assembly_stats.py`).
- MultiQC now also aggregates the per-assembly contiguity tables (the collect gained
  `.mix(asm_stats.mqc)`), so the single report covers gate zero plus the assemblies.

Data model: `assembler` and `coverage` travel inside the typed records through to the
published path (`assemblies/<accession>/<assembler>-<coverage>x/`), so the two
variants of one accession never collide. Adding a third assembler or coverage later is
a list edit, not a graph change.

Coverage variants run as `[220, 56]`: 220× reproduces the paper's "default
parameters" run (which used essentially all ~220× of the data), 56× is Supernova's
recommended raw coverage — running both tests whether the published assembly was
degraded by over-saturation.

## Validation

The whole DAG runs green offline in stub mode (`completed=19`: gate zero + 2 full
fetches + 2 conversions + 4 Supernova + 4 stats), with the 12 assembly artefacts
published under the keyed layout and no collisions. Unit tests pass (19: 11
barcode-check + 8 assembly-stats). Built task-by-task with stub validation at each
step; the `combine`→`record(...)`→typed-input fan-out compiled without needing an
explicit record-type annotation.

Local stub validation needs `-profile stub_local`, which caps the `SUPERNOVA`
process's production request (16 CPU / 250 GB) so it fits a dev machine. Production
runs use no profile and keep the real request.

## Open items (before the first real run)

These are deliberately deferred — the proprietary tools are not yet in hand, and stub
mode never executes them:

- **Build the private image.** Download `supernova-2.1.1.tar.gz` from the 10x portal
  (EULA), bake it + `stlfr2supernova` into one private image, and point
  `supernova_image` / `stlfr2supernova_image` at the pinned digests (override on the
  Launchpad).
- **Finalise two `TODO(smoke)` script bodies** against the real tools: the
  conversion-only invocation in `STLFR_CONVERT`, and the `mkoutput` output
  filename/gzip handling in `SUPERNOVA` (the current `mv ... 2>/dev/null || mv ...`
  fallback collapses to one deterministic path once the real output name is known).
- **Confirm on the smoke run:** the read length feeding the `--maxreads` arithmetic;
  the realised coverage in Supernova's own log vs our computed `--maxreads`; whether
  Supernova's internal scratch performs acceptably streamed through Fusion; and
  Spot-reclamation behaviour on the ~1-day task (consider on-demand for `SUPERNOVA`).
- **Confirm `genome_size`** against the CDC Frontier v2.0 assembly length.

## Outcome / next

Stage 1 is code-complete and stub-green on branch `stage1-baseline-assembly`. Next:
build the private image, then a smoke run on one accession at 56× (cheaper than 220×)
to finalise the `TODO(smoke)` invocations and confirm resource behaviour, before
fanning out to both accessions × both cutoffs.
