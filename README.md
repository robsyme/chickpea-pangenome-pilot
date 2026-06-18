# chickpea-pangenome-pilot

## Overarching target

Test whether we can produce **better, more honest de novo assemblies** of Australian
chickpea (*Cicer arietinum*) cultivars than those published in Garg et al. (2025),
[*An Australian chickpea pan-genome*](https://onlinelibrary.wiley.com/doi/10.1111/pbi.70192),
working only from the **stLFR linked-read data** the paper deposited
(NCBI BioProject [PRJNA1225167](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1225167)).

The paper's "chromosome-scale" assemblies get their contiguity from **RagTag
reference-guided scaffolding onto CDC Frontier**; the real de novo contig N50 is only
17–34 kb. For a pan-genome / structural-variant study that is partly circular —
contigs are forced onto the reference, then SVs are called against the same
reference, which suppresses the very rearrangements being catalogued. Our aim is to
quantify those artefacts and deliver **reference-free** assemblies that don't inherit
CDC Frontier's structure.

**Honest ceiling (kept explicit):** stLFR is barcoded *short* reads. Barcodes help
scaffold across repeats but cannot resolve repeat sequence longer than the ~100 bp
reads, and chickpea is ~45% repeat. So the deliverable is *methodologically better
and more honest* assemblies (reference-free, transparently QC'd), **not** HiFi/ONT-grade
contiguity — that would need new long-read sequencing. This target may be adapted as
results come in; changes are recorded in the progress journal below.

## Strategy (stages)

0. **Gate zero — barcode integrity.** Confirm the stLFR barcodes survived deposition
   before committing to reassembly. ✅ done (see journal).
1. **Baseline.** Reproduce the authors' route (`stlfr2supernova` → Supernova) for a
   comparable baseline. *(next)*
2. **Reference-free assembly + barcode scaffolding** (e.g. ARCS/ARKS, SLR-superscaffolder,
   cloudSPAdes) that does not use CDC Frontier for ordering.
3. **QC + artefact quantification** (BUSCO, k-mer QV, misassembly detection; suppressed
   SVs; fragmentation-driven false gene absences).

Pilot on a **single cultivar** (CBA Captain, `SRR32381426`) first; scale to all 15 only
if the pilot demonstrates clear artefacts and recovered SVs.

## What's in this repo

- **`main.nf` + `modules/`** — the Nextflow pipeline (the executable lab notebook).
  Currently implements stage 0 (gate zero): subsample reads from ENA, `seqkit` read
  structure, stLFR barcode-integrity check, aggregated into a MultiQC report.
- **`AGENTS.md`** — working context: infrastructure (Seqera Platform / AWS), validated
  typed-syntax notes, and hard-won gotchas. Read this before editing the pipeline.
- **`docs/`** — dated progress journal; one file per attempt/session. Start here to
  catch up across sessions.
- **`tests/`** — unit tests for the barcode-check logic.

## Running (current stage: gate zero)

```bash
export NXF_SYNTAX_PARSER=v2
nextflow run . --accessions SRR32381426 --check_reads 2000000 --pass_fraction 0.80
```

Requires Nextflow ≥ 26.04.4 and the strict v2 parser (`NXF_SYNTAX_PARSER=v2`). On
Seqera Platform, launch the `chickpea-gate-zero` Launchpad pipeline. Per-accession
results publish under `results/<accession>/`, with an aggregate `multiqc_report.html`.
Tests: `python3 -m unittest discover -s tests -v`.

## Status

Gate zero is green end-to-end on Seqera Platform. CBA Captain barcodes are intact
(90.3% valid on 2M reads, offsets 100/116/132) — the data is usable as linked-reads.
See `docs/` for the full history.
