# 2026-06-18 — Stage 1 redesign: open-stack ABySS contigs baseline

A design for the baseline assembly, replacing the Supernova route. Supernova
turned out to be a delisted proprietary 10x product with no self-service download
(the portal 404s; the CDN needs a signed URL only 10x could mint), so we pivot to
a fully open, Bioconda-packaged, stLFR-native stack. This supersedes
`2026-06-18-stage1-baseline-assembly-design.md`; the Supernova pipeline built under
that design is removed (git history preserves it).

## Goal

A defensible, fully-open de novo **contigs** baseline for each pilot accession,
comparable to the paper's published Supernova contig N50 (17–34 kb), with no
proprietary dependency. This is the first slice of the larger open pipeline; the
reference-free scaffolding chain is sketched below and built next.

## Why the open stack

The Birol lab (BCGSC) tools are open (GPL-3), on Bioconda with auto-built
biocontainers (so they drop into our Wave/Fusion/AWS Batch setup via a `conda`
directive, exactly like the gate-zero modules), stLFR-native, and actively
maintained. They also *are* the reference-free route the project actually wants:
they never touch CDC Frontier coordinates. cloudSPAdes was rejected (blew past 2 TB
RAM on stLFR in a published benchmark), Aquila is reference-assisted, Athena is
metagenome-only, Tell-Link is proprietary/TELL-Seq-only.

## Build scope

Contigs first. This effort builds `LRTK_CONVERT → ABYSS → ASSEMBLY_STATS` and gets
it green + smoke-validated. Scaffolding is a separate later build. Each tool is an
integration risk; this de-risks and yields an early, defensible number.

## Architecture

The pivot keeps the whole gate-zero subgraph, `FETCH_FULL_READS`, `ASSEMBLY_STATS`,
and the MultiQC merge. It replaces the two Supernova-specific processes.

```
runs ─┬─ (gate-zero subgraph, unchanged) ── MULTIQC
      └─ FETCH_FULL_READS ── LRTK_CONVERT ── (combine with abyss_kmers) ── ABYSS ── ASSEMBLY_STATS
```

- `LRTK_CONVERT` replaces `STLFR_CONVERT`: `lrtk FQCONVER -IT stLFR` emits
  `BX:Z:`-tagged insert FASTQs (barcode lifted out of R2). Reused by scaffolding.
- `ABYSS` replaces `SUPERNOVA`: one parameterised process fanned out over a k-sweep
  (`converted.combine(abyss_kmers)`), `abyss-pe` Bloom-filter mode, emits contigs.

**Removed:** `modules/supernova.nf`, `modules/stlfr_convert.nf`, and the params
`supernova_image`, `stlfr2supernova_image`, `genome_size`, `coverage_cutoffs`. The
`stub_local` profile is repurposed to cap `ABYSS` (not `SUPERNOVA`) for local runs.

**Reused unchanged:** `FETCH_FULL_READS`, the MultiQC merge, and the
`combine`→record fan-out idiom. `ASSEMBLY_STATS` is reused with a small edit (below).

The win: every process uses a Bioconda `conda` directive (LRTK, ABySS, seqkit). No
private image, no EULA, no manual binary — the entire Supernova blocker is gone.

## Data model

The Supernova-shaped records are redefined in `types.nf` (they are unused after the
pivot):

```groovy
record ConvertedReads {       // LRTK output: BX:Z:-tagged insert reads
    accession: String
    r1: Path
    r2: Path
}

record AssemblyJob {          // one per (accession × k) after combine
    accession: String
    r1: Path
    r2: Path
    k: Integer
}

record Assembly {
    accession: String
    assembler: String         // 'abyss-k64' — encodes the full variant
    fasta: Path
}

record AssemblyStats {
    accession: String
    assembler: String
    file: Path
}
```

Changes from the Supernova version: `ConvertedReads` carries `r1`/`r2` (the two
FASTQs ABySS consumes) instead of a `dir`; the variant identity collapses into a
single `assembler` string (`abyss-k64`), dropping the `coverage` field and the
Supernova-only `summary`. The k value lives structured in `AssemblyJob` (ABySS needs
it on the command line) and is baked into the `assembler` string for labelling and
the published path `assemblies/<accession>/abyss-k<k>/`. The `assembler` string keeps
multi-assembly-per-accession extensible: scaffolding later yields `abyss-k64-arks`,
etc., and every process and path keys on it so variants never collide.

## Processes

**`LRTK_CONVERT`** (`conda 'bioconda::lrtk=2.0'`): input `reads: ReadPair`, output
`ConvertedReads`.

```
lrtk FQCONVER -IT stLFR -I1 ${reads.r1} -I2 ${reads.r2} \
    -O1 ${acc}.bx_R1.fq.gz -O2 ${acc}.bx_R2.fq.gz   # TODO(smoke): confirm exact flags + whitelist
```

LRTK is the least-proven tool in the stack; the exact `FQCONVER` flags and whether a
whitelist path is needed are confirmed on a read subset at smoke time. The I/O
contract (ReadPair in, two BX-tagged FASTQs out) is final.

**`ABYSS`** (`conda 'bioconda::abyss=2.3.10'`): input `job: AssemblyJob`, output
`Assembly`. Bloom-filter mode for memory efficiency on a 740 Mb genome.

```
abyss-pe name=${acc}_k${job.k} k=${job.k} j=${task.cpus} \
    B=${params.abyss_bloom} in='${job.r1} ${job.r2}' \
    ${acc}_k${job.k}-contigs.fa
gzip -c ${acc}_k${job.k}-contigs.fa > ${acc}.abyss-k${job.k}.contigs.fa.gz
```

emits `record(accession: job.accession, assembler: "abyss-k${job.k}",
fasta: file("${acc}.abyss-k${job.k}.contigs.fa.gz"))`.

New params: `abyss_kmers: List<Integer> = [64, 80, 96]` (drives the fan-out; all
< the 100 bp read length — k cannot exceed read length) and `abyss_bloom: String =
'20G'` (tune at smoke). Resources `cpus 16`, `memory 64.GB`.

Fan-out wiring in `main.nf`:
```groovy
kmers     = channel.fromList(params.abyss_kmers)
converted = LRTK_CONVERT(full)
jobs      = converted.combine(kmers).map { conv, k ->
    record(accession: conv.accession, r1: conv.r1, r2: conv.r2, k: k)
}
contigs   = ABYSS(jobs)
asm_stats = ASSEMBLY_STATS(contigs)
```

## ASSEMBLY_STATS change

Reused as-is except the summariser `assembly_stats.py` drops its `--coverage`
argument/column; `--assembler` now carries the full variant string and the MultiQC
Sample becomes `<accession>-<assembler>`. This is a TDD edit: update
`tests/test_assembly_stats.py` first, then the binary.

## Target chain (recorded, built next)

The contigs baseline is the trunk for the scaffolding stage, which reuses the LRTK
BX-tagged reads and the ABySS contigs:

```
LRTK_CONVERT → ABYSS ─┬─ ASSEMBLY_STATS                              (baseline contig N50)
                      └─ TIGMINT → ARKS → LINKS → stLFR_GapCloser → ASSEMBLY_STATS  (reference-free scaffolds)
```

Scaffolding adds `Assembly` variants per accession (`abyss-k64-arks`, …), making
contigs vs scaffolds directly comparable under the same model — that comparison is
the headline result (contiguity gained from linked reads, reference-free).
`stLFR_GapCloser` is the only non-Bioconda tool and gets a DIY container then.

## Testing & validation

- `assembly_stats.py` `--coverage` removal: TDD (tests first).
- Stub validation: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -profile stub_local`,
  expected `completed=23` (7 gate-zero + 2 full-fetch + 2 LRTK + 6 ABySS
  [2 accessions × 3 k] + 6 stats).
- Anchoring: compare ABySS contig N50 to the paper's published Supernova contig N50
  (17–34 kb). No same-input Supernova run.

## Open items for the smoke run

- Confirm the exact `LRTK FQCONVER` flags + whitelist on a read subset (the one real
  unknown).
- Tune `abyss_bloom` for the realised k-mer count; confirm the chosen k values behave.
- Confirm the ABySS resource envelope at 220× and the contigs output filename.
