# 2026-06-18 — Stage 1 rebuilt: open ABySS contigs baseline (stub-validated)

## Goal of this phase

Replace the (removed) Supernova baseline with a fully-open, stLFR-native contigs
assembly per accession, comparable to the paper's published Supernova contig N50
(17–34 kb) and with no proprietary dependency. Design:
`2026-06-18-stage1b-abyss-contigs-design.md`; task plan:
`plans/2026-06-18-stage1b-abyss-contigs.md`.

## Why the pivot

Supernova turned out to be a delisted 10x product: the download portal 404s and the
CDN needs a signed URL only 10x could mint, so there is no self-service way to obtain
`supernova-2.1.1.tar.gz`. Rather than chase an entitled copy of a proprietary,
discontinued binary, the baseline moved to the Birol-lab (BCGSC) open stack, which is
GPL-3, on Bioconda, stLFR-native, actively maintained, and — usefully — is itself the
reference-free route (it never touches CDC Frontier coordinates).

## What we did

Tore the Supernova subgraph back to a clean gate-zero + full-fetch base, then rebuilt
additively, stub-validating at each step (`completed` 9 → 11 → 17 → 23). Reused
`FETCH_FULL_READS`, `ASSEMBLY_STATS` (with a small edit), the MultiQC merge, and the
`combine`→record fan-out idiom; the gate-zero processes were never touched.

- Dropped `--coverage` from the `assembly_stats.py` summariser (TDD): assemblies are
  now identified by an `assembler` string, not a numeric coverage.
- `LRTK_CONVERT` (`modules/lrtk_convert.nf`, `conda bioconda::lrtk`) — `lrtk FQCONVER
  -IT stLFR` lifts the stLFR barcode out of read2 into `BX:Z:` tags, emitting insert
  reads. Runs once per accession; reused by the scaffolding stage.
- `ABYSS` (`modules/abyss.nf`, `conda bioconda::abyss`) — one parameterised process
  fanned out over a k-sweep (`converted.combine(abyss_kmers)`, default `[64, 80, 96]`),
  `abyss-pe` Bloom-filter mode, emitting contigs tagged `assembler='abyss-k<k>'`.
- `ASSEMBLY_STATS` rewritten for the new record shapes (no coverage); publishes
  per-assembly contiguity stats and feeds its `_mqc` table into the MultiQC report.

Records redefined in `types.nf`: `ConvertedReads {accession, r1, r2}`, `AssemblyJob
{accession, r1, r2, k}`, `Assembly {accession, assembler, fasta}`, `AssemblyStats
{accession, assembler, file}`. The `assembler` string carries the variant to the
published path `assemblies/<accession>/abyss-k<k>/`, so the six variants never collide.

The headline win: every process uses a Bioconda `conda` directive (lrtk, abyss,
seqkit) resolved via Wave — no private image, no EULA, no manual binary. The whole
blocker that stopped the Supernova path is gone.

## Validation

Full DAG runs green offline in stub mode: `completed=23` (7 gate-zero + 2 full-fetch +
2 LRTK + 6 ABySS [2 accessions × 3 k] + 6 stats), with the 12 contig+stats artefacts
published under the keyed layout. Unit tests pass (20: 11 barcode-check + 9
assembly-stats). Local stub runs use `-profile stub_local` to cap ABYSS's 16-CPU
request to fit a dev machine; production runs use no profile.

## Open items (smaller now — no private image needed)

- Confirm the exact `lrtk FQCONVER` flags + whitelist on a read subset (LRTK is the
  least-proven tool; this is the one real `TODO(smoke)`).
- Confirm the `abyss-pe` contigs output filename and tune `abyss_bloom` for the realised
  k-mer count; confirm the ABySS memory envelope at 220×.
- Anchor: compare ABySS contig N50 per k to the paper's published Supernova contig N50.

## Outcome / next

Stage 1 baseline is code-complete and stub-green on branch `stage1b-abyss-contigs`,
with no proprietary dependency. Next: a smoke run on one accession to finalise the
`lrtk FQCONVER` invocation and pick a sensible k, then build the reference-free
scaffolding stage (`TIGMINT → ARKS → LINKS → stLFR_GapCloser`) on top of these contigs.
