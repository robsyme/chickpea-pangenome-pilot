# Progress journal

Dated records of what we attempted, what happened, and what we decided. Newest
last. The overarching target lives in the top-level `README.md`; working context and
gotchas live in `AGENTS.md`.

**Reference (not dated — update in place):**
- `environment.md` — Seqera Platform workspace, compute environment, storage paths,
  Launchpad pipeline, and version pinning. Start here to know *where* things run.

**Attempt log:**

- `2026-06-10-gate-zero-pipeline.md` — Stage 0 (gate zero) built and green on
  Platform; stLFR barcodes confirmed intact for CBA Captain.
- `2026-06-17-nextflow-record-bug-and-module-binaries.md` — A Nextflow typed-record
  staging bug surfaced, was filed and fixed upstream (26.04.4); module binaries
  re-adopted once the fix shipped.
- `2026-06-18-stage1-baseline-assembly-design.md` — Design for Stage 1: same-input
  Supernova baseline (`stlfr2supernova` → Supernova v2.1.1), two coverage variants
  (220× faithful + 56× best-practice) per accession, extending the gate-zero
  pipeline. Task-by-task plan in `plans/2026-06-18-stage1-baseline-assembly.md`.
- `2026-06-18-stage1-baseline-assembly-build.md` — Stage 1 built and stub-validated
  end to end (the four assembly processes, the coverage fan-out, MultiQC merge).
  Awaiting the private Supernova image for the first real run.
