# Stage 1 Baseline Assembly Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a same-input Supernova baseline assembly stage to the existing
gate-zero Nextflow pipeline, producing two coverage variants (220× faithful, 56×
best-practice) per accession, without disturbing the gate-zero task cache.

**Architecture:** New downstream subgraph hanging off the same per-accession
`runs` channel: `FETCH_FULL_READS → STLFR_CONVERT → SUPERNOVA(jobs) →
ASSEMBLY_STATS`. The convert step runs once per accession; a `combine` with a
coverage list fans one `SUPERNOVA` process out over `(accession × coverage)`.
Coverage and assembler identity travel inside typed records so the variants never
collide. Design doc: `docs/2026-06-18-stage1-baseline-assembly-design.md`.

**Tech Stack:** Nextflow ≥ 26.04.4 typed DSL2 (`NXF_SYNTAX_PARSER=v2`), Python 3
stdlib (unittest), seqkit, a private image with `stlfr2supernova` + Supernova 2.1.1.

---

## Conventions and the two validation loops

Every module begins with `nextflow.enable.types = true`. Records are built with
`record(...)`, never `new` (the `new` form breaks task caching — see AGENTS.md).
Each process has a `stub:` block so the whole DAG can run offline.

This plan uses two test loops:

1. **Python unit tests** (for the one piece of real logic, the stats summariser):
   ```
   python3 -m unittest tests.test_assembly_stats -v
   ```
2. **Nextflow stub run** (for all pipeline wiring — validates typed syntax,
   channel topology, record flow, and the output block, offline, in seconds):
   ```
   NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub
   ```
   A green run ends with `[SUCCESS] completed=<N> failed=0`. Cached tasks from a
   prior stub run show as `cached`; that is fine. Delete `/tmp/s1_stub` between
   runs if you want a clean output tree.

The TDD rhythm for a wiring task is: edit `main.nf` to reference the new module
first (stub run fails — module/process not found), then create the module with its
stub (stub run passes), then commit.

Do **not** modify `modules/fetch_reads.nf`, `modules/read_structure.nf`,
`modules/barcode_check/`, or `modules/multiqc.nf` until Task 7. Leaving the four
gate-zero processes' script/inputs/container untouched keeps their task hashes
stable so the production cache (especially the expensive fetches) survives.

---

### Task 0: Confirm the baseline is green

**Step 1: Run the existing stub pipeline**

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS] completed=7 failed=0` with `read_structure`,
`barcode_check`, and `multiqc_report` outputs listed.

**Step 2: Run the existing unit tests**

Run: `python3 -m unittest discover -s tests -v`
Expected: `Ran 11 tests`, `OK`.

If either fails, stop and report — that is a pre-existing problem, not part of
this stage.

---

### Task 1: Assembly stats summariser (Python, TDD)

The only non-trivial logic in the stage: parse `seqkit stats -a -T` output and
render both a machine-readable TSV and a MultiQC custom-content table. Mirrors
`modules/barcode_check/resources/usr/bin/stlfr_barcode_check.py` and its tests.

**Files:**
- Create: `tests/test_assembly_stats.py`
- Create: `modules/assembly_stats/resources/usr/bin/assembly_stats.py`

**Step 1: Write the failing tests**

Create `tests/test_assembly_stats.py`:

```python
"""Tests for modules/assembly_stats/resources/usr/bin/assembly_stats.py
(stdlib unittest, no pytest needed).

Run: python3 -m unittest discover -s tests -v
"""
import importlib.util
import pathlib
import unittest

_BIN = (pathlib.Path(__file__).resolve().parent.parent
        / "modules" / "assembly_stats" / "resources" / "usr" / "bin"
        / "assembly_stats.py")
_spec = importlib.util.spec_from_file_location("asm", _BIN)
asm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(asm)

# A representative `seqkit stats -a -T` row for a FASTA (quality cols are '0').
SEQKIT = (
    "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\t"
    "Q1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)\n"
    "asm.fa.gz\tFASTA\tDNA\t1200\t700000000\t500\t583333\t5000000\t"
    "0\t0\t0\t14000000\t23000\t9000\t0\t0\t0\t33.5\n"
)


class TestParseSeqkit(unittest.TestCase):
    def test_extracts_fields_by_name(self):
        m = asm.parse_seqkit(SEQKIT)
        self.assertEqual(m["num_seqs"], 1200)
        self.assertEqual(m["sum_len"], 700000000)
        self.assertEqual(m["n50"], 23000)
        self.assertEqual(m["sum_gap"], 14000000)

    def test_missing_column_defaults_zero(self):
        text = "num_seqs\tsum_len\tN50\n10\t1000\t100\n"
        self.assertEqual(asm.parse_seqkit(text)["sum_gap"], 0)

    def test_no_data_row_raises(self):
        with self.assertRaises(ValueError):
            asm.parse_seqkit("num_seqs\tsum_len\n")


class TestPctN(unittest.TestCase):
    def test_basic(self):
        self.assertEqual(asm.pct_n(14000000, 700000000), 2.0)

    def test_zero_length(self):
        self.assertEqual(asm.pct_n(0, 0), 0.0)


class TestRenderMqc(unittest.TestCase):
    def setUp(self):
        self.m = asm.parse_seqkit(SEQKIT)
        self.out = asm.render_mqc("SRR32381426", "supernova", 56, self.m)

    def test_has_custom_content_config_header(self):
        self.assertIn("# plot_type: 'table'", self.out)
        self.assertIn("# id: 'assembly_stats'", self.out)

    def test_row_keyed_by_strategy_sample(self):
        self.assertIn("SRR32381426-supernova-56x", self.out)
        self.assertIn("23000", self.out)


class TestRenderStats(unittest.TestCase):
    def test_tsv_has_identity_and_metrics(self):
        m = asm.parse_seqkit(SEQKIT)
        out = asm.render_stats("SRR32381426", "supernova", 220, m)
        lines = out.strip().splitlines()
        self.assertEqual(lines[0].split("\t")[0], "accession")
        row = lines[1].split("\t")
        self.assertEqual(row[0], "SRR32381426")
        self.assertEqual(row[2], "220")
        self.assertEqual(row[5], "23000")  # n50 column


if __name__ == "__main__":
    unittest.main()
```

**Step 2: Run tests to verify they fail**

Run: `python3 -m unittest tests.test_assembly_stats -v`
Expected: FAIL — `FileNotFoundError` / import error (the script does not exist yet).

**Step 3: Write the implementation**

Create `modules/assembly_stats/resources/usr/bin/assembly_stats.py`:

```python
#!/usr/bin/env python3
"""Summarise an assembly FASTA from `seqkit stats -a -T` output.

Emits a machine-readable stats.tsv (augmented with assembly identity) and a
MultiQC custom-content table row. Importable for unit tests (no import-time work).
"""
import argparse


def parse_seqkit(text):
    """Parse `seqkit stats -a -T` text into the fields we report.

    Indexed by column NAME, not position, so seqkit version differences in
    column order do not break parsing. Returns ints (0 if a column is absent).
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        raise ValueError("seqkit stats output has no data row")
    cells = dict(zip(lines[0].split("\t"), lines[1].split("\t")))

    def as_int(key):
        raw = cells.get(key, "0").replace(",", "")
        try:
            return int(float(raw))
        except ValueError:
            return 0

    return {
        "num_seqs": as_int("num_seqs"),
        "sum_len": as_int("sum_len"),
        "n50": as_int("N50"),
        "sum_gap": as_int("sum_gap"),
    }


def pct_n(sum_gap, sum_len):
    """Percent of the assembly that is N (gap) bases."""
    if sum_len <= 0:
        return 0.0
    return round(100.0 * sum_gap / sum_len, 3)


def render_stats(accession, assembler, coverage, m):
    """Tidy machine-readable one-row TSV with identity + metrics."""
    header = ["accession", "assembler", "coverage", "num_scaffolds",
              "total_length", "n50", "n_bases", "pct_n"]
    row = [accession, assembler, str(coverage),
           str(m["num_seqs"]), str(m["sum_len"]), str(m["n50"]),
           str(m["sum_gap"]), str(pct_n(m["sum_gap"], m["sum_len"]))]
    return "\t".join(header) + "\n" + "\t".join(row) + "\n"


def render_mqc(accession, assembler, coverage, m):
    """MultiQC custom-content table row for one strategy."""
    sample = f"{accession}-{assembler}-{coverage}x"
    lines = [
        "# id: 'assembly_stats'",
        "# section_name: 'Assembly stats'",
        "# plot_type: 'table'",
        "\t".join(["Sample", "Assembler", "Coverage", "Scaffolds",
                   "Total length (bp)", "N50 (bp)", "N (%)"]),
        "\t".join([sample, assembler, f"{coverage}x",
                   str(m["num_seqs"]), str(m["sum_len"]),
                   str(m["n50"]), str(pct_n(m["sum_gap"], m["sum_len"]))]),
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seqkit", required=True, help="seqkit stats -a -T TSV")
    p.add_argument("--accession", required=True)
    p.add_argument("--assembler", required=True)
    p.add_argument("--coverage", required=True, type=int)
    p.add_argument("--out", required=True)
    p.add_argument("--mqc", required=True)
    args = p.parse_args(argv)

    with open(args.seqkit) as fh:
        m = parse_seqkit(fh.read())
    with open(args.out, "w") as fh:
        fh.write(render_stats(args.accession, args.assembler, args.coverage, m))
    with open(args.mqc, "w") as fh:
        fh.write(render_mqc(args.accession, args.assembler, args.coverage, m))


if __name__ == "__main__":
    main()
```

Then make it executable: `chmod +x modules/assembly_stats/resources/usr/bin/assembly_stats.py`

**Step 4: Run tests to verify they pass**

Run: `python3 -m unittest tests.test_assembly_stats -v`
Expected: `Ran 8 tests`, `OK`.

**Step 5: Commit**

```bash
git add tests/test_assembly_stats.py modules/assembly_stats/resources/usr/bin/assembly_stats.py
git commit -m "feat: assembly stats summariser + unit tests"
```

---

### Task 2: New record types

**Files:**
- Modify: `types.nf`

**Step 1: Add the records**

Append to `types.nf`:

```groovy
record ConvertedReads {
    accession: String
    dir: Path
}

record AssemblyJob {
    accession: String
    dir: Path
    coverage: Integer
}

record Assembly {
    accession: String
    assembler: String
    coverage: Integer
    fasta: Path
    summary: Path
}

record AssemblyStats {
    accession: String
    assembler: String
    coverage: Integer
    file: Path
}
```

**Step 2: Validate the file still parses**

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS] completed=7` (records are defined but unused — nothing should
break). A parse error in the new records would fail compilation here.

**Step 3: Commit**

```bash
git add types.nf
git commit -m "feat: add Stage 1 record types"
```

---

### Task 3: FETCH_FULL_READS

Downloads the FULL run (no subsample) for assembly. Deliberately separate from
the gate-zero `FETCH_READS` so that process's hash stays stable.

**Files:**
- Create: `modules/fetch_full_reads.nf`
- Modify: `main.nf`

**Step 1: Wire it into `main.nf` first (red)**

Add include near the other includes (`main.nf` top):

```groovy
include { FETCH_FULL_READS } from './modules/fetch_full_reads.nf'
```

In the workflow `main:` block, after the existing `report = MULTIQC(...)` line, add:

```groovy
    // --- Stage 1: baseline Supernova assembly --------------------------------
    full = FETCH_FULL_READS(runs)
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: FAIL — unable to resolve `./modules/fetch_full_reads.nf`.

**Step 2: Create the module (green)**

Create `modules/fetch_full_reads.nf`:

```groovy
nextflow.enable.types = true

include { SraRun ; ReadPair } from '../types.nf'

// Download the FULL run for an accession (no subsample) — needed for assembly.
// Separate from FETCH_READS (which subsamples for gate zero) so the gate-zero
// task hash stays stable and its cache survives. fastq-dl auto-falls back ENA->SRA.
process FETCH_FULL_READS {
    tag "${run.accession}"
    container 'community.wave.seqera.io/library/fastq-dl:3.0.1--fa446f61dfc85bc3'

    input:
    run: SraRun

    output:
    reads: ReadPair = record(
        accession: run.accession,
        r1: file("${run.accession}_R1.fastq.gz"),
        r2: file("${run.accession}_R2.fastq.gz")
    )

    script:
    def acc = run.accession
    """
    set -euo pipefail
    rm -rf dl && mkdir dl
    fastq-dl --accession ${acc} --outdir dl >&2
    mv "\$(ls dl/*_1.fastq.gz | head -1)" ${acc}_R1.fastq.gz
    mv "\$(ls dl/*_2.fastq.gz | head -1)" ${acc}_R2.fastq.gz
    rm -rf dl
    """

    stub:
    def acc = run.accession
    """
    printf '@r/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}_R1.fastq.gz
    printf '@r/2\\nACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIII\\n' | gzip -c > ${acc}_R2.fastq.gz
    """
}
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS]`, now including `FETCH_FULL_READS (SRR32381426)` and
`(SRR32381425)` (completed count rises by 2 to 9).

**Step 3: Commit**

```bash
git add modules/fetch_full_reads.nf main.nf
git commit -m "feat: FETCH_FULL_READS full-run download for assembly"
```

---

### Task 4: STLFR_CONVERT

Runs `stlfr2supernova` conversion (steps 0–2 only: split_barcode → SOAPfilter →
merge_barcodes) to produce Chromium-style FASTQs. Runs once per accession.

**Files:**
- Create: `modules/stlfr_convert.nf`
- Modify: `main.nf`

**Step 1: Add the image param and wiring (red)**

In the `params {}` block of `main.nf`, add:

```groovy
    // Private image carrying stlfr2supernova (+ SOAPfilter). Override on the
    // Launchpad once built/pushed; placeholder keeps the repo portable.
    stlfr2supernova_image: String = 'PLACEHOLDER/stlfr2supernova:0.1'
```

Add the include and, after the `full = ...` line:

```groovy
include { STLFR_CONVERT } from './modules/stlfr_convert.nf'
```
```groovy
    converted = STLFR_CONVERT(full)
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: FAIL — cannot resolve `./modules/stlfr_convert.nf`.

**Step 2: Create the module (green)**

Create `modules/stlfr_convert.nf`:

```groovy
nextflow.enable.types = true

include { ReadPair ; ConvertedReads } from '../types.nf'

// Convert MGI stLFR co-barcoded reads into Chromium-style FASTQs that Supernova
// accepts (stlfr2supernova steps 0-2). Coverage-independent and expensive, so it
// runs ONCE per accession; both coverage variants assemble from this one set.
//
// NOTE: the exact in-container CLI for running ONLY the conversion (stopping
// before Supernova) must be confirmed against the stlfr2supernova_pipeline scripts
// when the image is built. The wiring and I/O contract below are final; the script
// body is the best-effort invocation to finalise at smoke-test time.
process STLFR_CONVERT {
    tag "${reads.accession}"
    container params.stlfr2supernova_image

    input:
    reads: ReadPair

    output:
    converted: ConvertedReads = record(
        accession: reads.accession,
        dir: file("${reads.accession}_chromium")
    )

    script:
    def acc = reads.accession
    """
    set -euo pipefail
    mkdir -p ${acc}_chromium
    # TODO(smoke): confirm exact stlfr2supernova conversion-only invocation.
    stlfr2supernova_convert \\
        --r1 ${reads.r1} \\
        --r2 ${reads.r2} \\
        --outdir ${acc}_chromium
    """

    stub:
    def acc = reads.accession
    """
    mkdir -p ${acc}_chromium
    printf '@r 1:N:0:0\\nACGT\\n+\\nIIII\\n' | gzip -c \\
        > ${acc}_chromium/read-RA_si-AAAA_lane-001-chunk-001.fastq.gz
    """
}
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS]` with `STLFR_CONVERT` tasks (completed rises by 2 to 11).

**Step 3: Commit**

```bash
git add modules/stlfr_convert.nf main.nf
git commit -m "feat: STLFR_CONVERT stLFR->Chromium conversion (convert once)"
```

---

### Task 5: SUPERNOVA (the combine fan-out)

The riskiest typed-syntax step: `combine` a coverage list with the converted
reads, `map` into an `AssemblyJob` record, and feed one parameterised process.

**Files:**
- Create: `modules/supernova.nf`
- Modify: `main.nf`

**Step 1: Add params + fan-out wiring (red)**

In `params {}` of `main.nf` add:

```groovy
    // Stage 1 coverage variants and genome size for --maxreads sizing.
    genome_size: Integer = 740_000_000          // Cicer arietinum ~740 Mb
    coverage_cutoffs: List<Integer> = [220, 56] // faithful + best-practice
    // Private image with Supernova 2.1.1 baked in. Override on the Launchpad.
    supernova_image: String = 'PLACEHOLDER/supernova:2.1.1'
```

Add the include and, after `converted = ...`:

```groovy
include { SUPERNOVA } from './modules/supernova.nf'
```
```groovy
    cutoffs    = channel.fromList(params.coverage_cutoffs)
    jobs       = converted.combine(cutoffs).map { conv, cov ->
        record(accession: conv.accession, dir: conv.dir, coverage: cov)
    }
    assemblies = SUPERNOVA(jobs, params.genome_size)
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: FAIL — cannot resolve `./modules/supernova.nf`.

> **If, after Step 2, the stub fails on the `map { conv, cov -> record(...) }`
> line** (record type not inferred in workflow context), the fallback is to add an
> explicit annotation: `jobs = converted.combine(cutoffs).map { conv, cov ->
> AssemblyJob job = record(...); job }` or import + construct via the typed form
> the parser accepts. Validate empirically per AGENTS.md. Do NOT switch to
> `new AssemblyJob(...)` — that breaks caching.

**Step 2: Create the module (green)**

Create `modules/supernova.nf`:

```groovy
nextflow.enable.types = true

include { AssemblyJob ; Assembly } from '../types.nf'

// Assemble one (accession x coverage) job with Supernova. Sizes --maxreads from
// the target coverage, genome size, and measured read length, runs `supernova run`,
// and extracts a pseudohap scaffold FASTA with `mkoutput`.
//
// NOTE: Supernova is a proprietary 10x binary baked into a private image
// (params.supernova_image). The run/mkoutput flags are standard but confirm the
// mkoutput output filename + gzip behaviour on the smoke run; the I/O contract and
// maxreads arithmetic below are final.
process SUPERNOVA {
    tag "${job.accession}-${job.coverage}x"
    container params.supernova_image
    cpus 16
    memory 250.GB

    input:
    job: AssemblyJob
    genome_size: Integer

    output:
    assembly: Assembly = record(
        accession: job.accession,
        assembler: 'supernova',
        coverage: job.coverage,
        fasta: file("${job.accession}.${job.coverage}x.scaffolds.fa.gz"),
        summary: file("${job.accession}.${job.coverage}x.summary.csv")
    )

    script:
    def acc = job.accession
    def cov = job.coverage
    def id  = "${acc}_${cov}x"
    """
    set -euo pipefail

    # Size --maxreads to hit the requested raw coverage:
    #   maxreads = genome_size * coverage / read_len
    # read_len is measured from the converted reads (not assumed).
    one=\$(ls ${job.dir}/*.fastq.gz | head -1)
    read_len=\$(gzip -dc "\$one" | awk 'NR==2{print length(\$0); exit}')
    [ "\$read_len" -gt 0 ] || { echo "ERROR: could not measure read length" >&2; exit 1; }
    maxreads=\$(awk -v g=${genome_size} -v c=${cov} -v l=\$read_len 'BEGIN{printf "%d", (g*c)/l}')
    echo "[SUPERNOVA] ${acc} target=${cov}x read_len=\$read_len maxreads=\$maxreads" >&2

    localmem=\$(( ${task.memory.toGiga()} - 6 ))

    supernova run \\
        --id=${id} \\
        --fastqs=${job.dir} \\
        --maxreads=\$maxreads \\
        --localcores=${task.cpus} \\
        --localmem=\$localmem

    # TODO(smoke): confirm mkoutput writes <outprefix>.fasta.gz (then drop the gzip).
    supernova mkoutput \\
        --style=pseudohap \\
        --asmdir=${id}/outs/assembly \\
        --outprefix=${acc}.${cov}x.scaffolds
    [ -f ${acc}.${cov}x.scaffolds.fasta.gz ] || gzip -f ${acc}.${cov}x.scaffolds.fasta
    mv ${acc}.${cov}x.scaffolds.fasta.gz ${acc}.${cov}x.scaffolds.fa.gz 2>/dev/null || \\
        mv ${acc}.${cov}x.scaffolds.fasta ${acc}.${cov}x.scaffolds.fa.gz

    cp ${id}/outs/summary.csv ${acc}.${cov}x.summary.csv
    """

    stub:
    def acc = job.accession
    def cov = job.coverage
    """
    printf '>scaffold_1\\nACGTACGTACGT\\n' | gzip -c > ${acc}.${cov}x.scaffolds.fa.gz
    printf 'metric,value\\nest_genome_size,740000000\\n' > ${acc}.${cov}x.summary.csv
    """
}
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS]` with **four** `SUPERNOVA` tasks
(`SRR32381426-220x`, `SRR32381426-56x`, `SRR32381425-220x`, `SRR32381425-56x`) —
this confirms the `accession × coverage` fan-out (completed rises by 4 to 15).

**Step 3: Commit**

```bash
git add modules/supernova.nf main.nf
git commit -m "feat: SUPERNOVA assembly with coverage fan-out via combine"
```

---

### Task 6: ASSEMBLY_STATS

Wraps the Task 1 summariser as a module-scoped binary.

**Files:**
- Create: `modules/assembly_stats/main.nf`
- Modify: `main.nf`

**Step 1: Wire it in (red)**

Add include and, after `assemblies = ...`:

```groovy
include { ASSEMBLY_STATS } from './modules/assembly_stats/main.nf'
```
```groovy
    asm_stats = ASSEMBLY_STATS(assemblies)
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: FAIL — cannot resolve `./modules/assembly_stats/main.nf`.

**Step 2: Create the module (green)**

Create `modules/assembly_stats/main.nf`:

```groovy
nextflow.enable.types = true

include { Assembly ; AssemblyStats } from '../../types.nf'

// Contiguity summary for one assembly: scaffold count, total length, N50, %N.
// Emits a machine-readable stats.tsv (AssemblyStats record) plus a MultiQC
// custom-content row (mirrors BARCODE_CHECK). Summariser is a module-scoped
// binary (resources/usr/bin, on PATH). Requires Nextflow >= 26.04.4 for record
// Path-field staging under Fusion with module binaries (#7225/#7226).
process ASSEMBLY_STATS {
    tag "${assembly.accession}-${assembly.coverage}x"
    conda 'bioconda::seqkit=2.8.2 conda-forge::python=3.12'

    input:
    assembly: Assembly

    output:
    res: AssemblyStats = record(
        accession: assembly.accession,
        assembler: assembly.assembler,
        coverage: assembly.coverage,
        file: file("${assembly.accession}.${assembly.coverage}x.stats.tsv")
    )
    mqc: Path = file("${assembly.accession}.${assembly.coverage}x.assembly_mqc.tsv")

    script:
    def acc = assembly.accession
    def cov = assembly.coverage
    """
    seqkit stats -a -T ${assembly.fasta} > seqkit.tsv
    assembly_stats.py \\
        --seqkit seqkit.tsv \\
        --accession ${acc} \\
        --assembler ${assembly.assembler} \\
        --coverage ${cov} \\
        --out ${acc}.${cov}x.stats.tsv \\
        --mqc ${acc}.${cov}x.assembly_mqc.tsv
    """

    stub:
    def acc = assembly.accession
    def cov = assembly.coverage
    """
    printf 'accession\\tassembler\\tcoverage\\tn50\\n${acc}\\t${assembly.assembler}\\t${cov}\\t23000\\n' \\
        > ${acc}.${cov}x.stats.tsv
    printf "# id: 'assembly_stats'\\n# plot_type: 'table'\\nSample\\tN50 (bp)\\n${acc}-${assembly.assembler}-${cov}x\\t23000\\n" \\
        > ${acc}.${cov}x.assembly_mqc.tsv
    """
}
```

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS]` with four `ASSEMBLY_STATS` tasks (completed rises by 4 to 19).

**Step 3: Commit**

```bash
git add modules/assembly_stats/main.nf main.nf
git commit -m "feat: ASSEMBLY_STATS contiguity summary per assembly"
```

---

### Task 7: Publish assemblies and wire stats into MultiQC

Adds the `publish:`/`output {}` entries and folds the assembly `_mqc` rows into the
existing MultiQC report. **This is the only task that touches a gate-zero process
(MULTIQC), so it lands last** — it changes MULTIQC's hash (it re-runs once), which
is cheap, while every expensive new task is already cached by this point.

**Files:**
- Modify: `main.nf`

**Step 1: Mix assembly stats into the MultiQC collect**

Change the existing `qc_files` assignment from:

```groovy
    qc_files = structure
        .map { r -> r.file }
        .mix(barcode.mqc)
        .collect()
```

to:

```groovy
    qc_files = structure
        .map { r -> r.file }
        .mix(barcode.mqc)
        .mix(asm_stats.mqc)
        .collect()
```

**Step 2: Add publish + output entries**

In the `publish:` section add:

```groovy
    assembly       = assemblies
    assembly_stats = asm_stats.res
```

In the `output {}` block add (alongside the existing entries), and add
`Assembly ; AssemblyStats` to the `from './types.nf'` include:

```groovy
    assembly: Channel<Assembly> {
        path { a -> "assemblies/${a.accession}/${a.assembler}-${a.coverage}x" }
    }
    assembly_stats: Channel<AssemblyStats> {
        path { s -> "assemblies/${s.accession}/${s.assembler}-${s.coverage}x" }
    }
```

**Step 3: Validate**

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS]`. The `Outputs:` block now lists `assembly` and
`assembly_stats`, each with four entries published under
`assemblies/<accession>/supernova-<cov>x/`.

**Step 4: Commit**

```bash
git add main.nf
git commit -m "feat: publish assemblies + merge assembly stats into MultiQC"
```

---

### Task 8: Full-DAG verification and unit-test sweep

**Step 1: Clean stub run of the whole pipeline**

Run: `rm -rf /tmp/s1_stub && NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1_stub`
Expected: `[SUCCESS] completed=19 failed=0` (7 gate-zero + 2 fetch-full +
2 convert + 4 supernova + 4 assembly-stats), with assemblies published for both
accessions at both coverages.

**Step 2: Confirm the published layout**

Run: `find /tmp/s1_stub/assemblies -type f | sort`
Expected: under each `SRR32381426/` and `SRR32381425/`, a `supernova-220x/` and
`supernova-56x/` directory, each containing the scaffold FASTA, summary CSV, and
stats TSV.

**Step 3: Run all unit tests**

Run: `python3 -m unittest discover -s tests -v`
Expected: `Ran 19 tests` (11 barcode-check + 8 assembly-stats), `OK`.

**Step 4: Commit (if anything changed; otherwise skip)**

No code change expected — this task is verification only.

---

### Task 9: Update docs and journal

**Files:**
- Modify: `docs/README.md` (add this plan and a forthcoming attempt-log line later)
- Modify: `AGENTS.md` (Planned stages → mark Stage 1 as in progress; note the new params and the manual image prerequisite)

**Step 1: Note the human prerequisites for the first real run**

These are NOT code and block only the real run, not the stub-validated build:
- Download `supernova-2.1.1.tar.gz` from the 10x portal (EULA), build + push the
  private image, set `params.supernova_image` / `params.stlfr2supernova_image`
  (or override them on the Launchpad) to the pinned digests.
- Finalise the two `TODO(smoke)` CLI invocations (`stlfr2supernova` conversion-only;
  `mkoutput` output filename) against the real tools.
- Confirm `genome_size` against CDC Frontier v2.0.

**Step 2: Commit**

```bash
git add docs/README.md AGENTS.md
git commit -m "docs: Stage 1 build notes and prerequisites"
```

---

## First real run (after the manual prerequisites, separate from this build)

Relaunch the existing Launchpad pipeline pinned to this branch's commit so gate
zero stays cached and only the new tasks run. Probe cost on the cheaper variant
first by temporarily setting `coverage_cutoffs = [56]` for one accession, then
restore `[220, 56]` for both:

```bash
export TOWER_WORKSPACE_ID=66005590816034
tw runs relaunch -i <gate-zero-runId> --revision stage1-baseline-assembly \
  --commit-id $(git rev-parse HEAD)
```

Watch on the smoke run: realised coverage in Supernova's own log vs our computed
`--maxreads`; whether Supernova's internal scratch performs acceptably under
Fusion; and Spot-reclamation behaviour on the long task (consider on-demand for
`SUPERNOVA`).

---

## Notes for the executor

- Run every command from the worktree root
  (`.worktrees/stage1-baseline-assembly`).
- The stub run writes a `work/` dir in the worktree — gitignored, leave it.
- Keep changes cache-friendly: do not touch the four gate-zero process files until
  Task 7, and only the MULTIQC collect there.
- The typed syntax is a preview feature; the `[WARN] Static typing is a preview`
  line is expected and harmless.
