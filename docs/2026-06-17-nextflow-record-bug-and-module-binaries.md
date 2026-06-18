# 2026-06-17 — Nextflow typed-record staging bug; module binaries re-adopted

## Context

Part of this pilot is stress-testing Nextflow's new typed/static-types syntax. We
moved the barcode-check script to a **module-scoped binary**
(`modules/barcode_check/resources/usr/bin/`, `nextflow.enable.moduleBinaries`).
Runs then failed on Fusion with "No such file or directory".

## Investigation

- Symptom: a downstream task's input file path was written into `.command.sh` as a
  raw object-store path `/scidev-…/…fastq.gz` — **missing the `/fusion/s3` mount
  prefix** — while the file existed in S3. Direct `Path` inputs were fine.
- Reframed correctly (credit: user) — `.command.sh` is written by the head job
  *before* Fusion starts, so this is **head-job path rendering**, not a runtime/CE
  issue. (A detour considering Fusion Snapshots / CE config was dropped.)
- `getClass()` probe at script-render time: a record-field `Path` rendered as
  `nextflow.cloud.aws.nio.S3Path` (raw, unstaged) while a direct `Path` rendered as
  `nextflow.processor.TaskPath` (staged). Locally both are `TaskPath` → **the bug is
  S3-specific** and does not reproduce on a local filesystem.
- Source trace: `TaskInputResolver.normalizePath` keeps a raw `Path` when
  `holders.containsKey(value)` misses; a raw `S3Path` then bypasses
  `FusionHelper.toContainerMount` (which builds `/fusion/s3/...`).

## Root cause (per upstream)

We initially attributed the trigger to `nextflow.enable.moduleBinaries`. The Nextflow
maintainers **retitled** the issue to the true cause: *record `Path` fields are not
staged when the **record type is included from another module**.* Our repro had both
(record included from `types.nf` **and** the flag), so the flag was a confound. The
diagnosis (record `Path` fields unstaged → raw `S3Path`) was right; the trigger
attribution was not.

## Filed and fixed

- Issue: https://github.com/nextflow-io/nextflow/issues/7225
- Fix PR: #7226 "Fix staging of Path fields of included record types" (merged to
  master `72ab74963`).
- Released in **Nextflow 26.04.4** (backport `6bc3510ae`; confirmed ancestor of the
  `v26.04.4` tag and listed in its release notes).
- Minimal reproducer repo: https://github.com/robsyme/nf-record-fusion-repro

## Action taken

- While the bug was live, the pilot used a top-level `bin/` script (plain process) —
  which works because it sidesteps the included-record staging path.
- With the fix in 26.04.4, **re-adopted the module-scoped binary** for BARCODE_CHECK
  and bumped `manifest.nextflowVersion` to `>=26.04.4`. Verified end-to-end on the
  Fusion CE pinned to `NXF_VER=26.04.4` (the repro's record-field consumer now
  stages correctly).

## Lessons

- A documented typed-syntax feature (record `Path` fields staged like inputs) had a
  real Fusion staging defect; stress-testing surfaced it.
- `tw runs relaunch --pull-latest` keeps the **original run's commit** — pin
  `--commit-id $(git rev-parse HEAD)` to resume onto new code.
- Editing a Launchpad pipeline's config via the `tw` CLI clobbers the UI-only "v2
  parser" toggle; manage that pipeline's config in the UI.
