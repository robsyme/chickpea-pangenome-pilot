# Open ABySS Contigs Baseline Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the (removed) Supernova baseline with a fully-open contigs baseline:
`FETCH_FULL_READS → LRTK_CONVERT → ABYSS (k-sweep) → ASSEMBLY_STATS`, all Bioconda.

**Architecture:** Tear the Supernova assembly subgraph back to a clean gate-zero +
full-fetch base, then rebuild additively with two new Bioconda processes (LRTK barcode
extraction, ABySS contig assembly fanned out over a k-sweep), reusing `FETCH_FULL_READS`,
`ASSEMBLY_STATS`, the MultiQC merge, and the `combine`→record fan-out idiom. Design:
`docs/2026-06-18-stage1b-abyss-contigs-design.md`.

**Tech Stack:** Nextflow ≥26.04.4 typed DSL2 (`NXF_SYNTAX_PARSER=v2`), Bioconda
(`lrtk`, `abyss`, `seqkit`) via `conda` directives + Wave, Python 3 stdlib (unittest).

---

## Validation loops & conventions

Two loops, as in the prior build:
1. **Python unit tests:** `python3 -m unittest discover -s tests -v`
2. **Nextflow stub run:** `NXF_SYNTAX_PARSER=v2 nextflow run . -stub [-profile stub_local] -output-dir /tmp/<tag>`
   Green ends with `[SUCCESS] completed=<N> failed=0`. The `-profile stub_local` flag is
   only needed once `ABYSS` exists (Task 4+), to cap its 16-CPU request for a dev machine.

Every module starts with `nextflow.enable.types = true`; records are built with
`record(...)`, never `new`. Do NOT modify the gate-zero processes (`fetch_reads.nf`,
`read_structure.nf`, `barcode_check/`, `multiqc.nf`). Expected harmless warnings:
`Static typing is a preview feature` and `Access to undefined parameter` (modules reading
`params.x`).

Run all commands from the worktree root
(`.worktrees/stage1b-abyss-contigs`).

---

### Task 0: Confirm the baseline

Run: `python3 -m unittest discover -s tests` → expect `Ran 19 tests`, `OK`.
Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -profile stub_local -output-dir /tmp/s1b_base`
→ expect `[SUCCESS] completed=19 failed=0` (the merged Supernova pipeline; this is the
starting point we pivot away from). If either fails, stop and report.

---

### Task 1: Drop `--coverage` from the stats summariser (TDD)

The open stack has no per-assembly coverage variant; identity is the `assembler` string
(`abyss-k64`). Update the summariser and its tests.

**Files:**
- Modify: `tests/test_assembly_stats.py`
- Modify: `modules/assembly_stats/resources/usr/bin/assembly_stats.py`

**Step 1: Replace the test file** with this exact content:

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
        self.out = asm.render_mqc("SRR32381426", "abyss-k64", self.m)

    def test_has_custom_content_config_header(self):
        self.assertIn("# plot_type: 'table'", self.out)
        self.assertIn("# id: 'assembly_stats'", self.out)

    def test_row_keyed_by_assembler_sample(self):
        self.assertIn("SRR32381426-abyss-k64", self.out)
        self.assertIn("23000", self.out)

    def test_no_coverage_column(self):
        self.assertNotIn("Coverage", self.out)


class TestRenderStats(unittest.TestCase):
    def test_tsv_has_identity_and_metrics(self):
        m = asm.parse_seqkit(SEQKIT)
        out = asm.render_stats("SRR32381426", "abyss-k64", m)
        lines = out.strip().splitlines()
        header = lines[0].split("\t")
        row = lines[1].split("\t")
        self.assertEqual(header[0], "accession")
        self.assertNotIn("coverage", header)
        self.assertEqual(row[0], "SRR32381426")
        self.assertEqual(row[1], "abyss-k64")
        self.assertEqual(row[header.index("n50")], "23000")


if __name__ == "__main__":
    unittest.main()
```

**Step 2: Run tests to verify they fail**

Run: `python3 -m unittest tests.test_assembly_stats -v`
Expected: FAIL (the binary still requires/uses `coverage`; `render_*` signatures mismatch).

**Step 3: Replace `modules/assembly_stats/resources/usr/bin/assembly_stats.py`** with:

```python
#!/usr/bin/env python3
"""Summarise an assembly FASTA from `seqkit stats -a -T` output.

Emits a machine-readable stats.tsv (augmented with assembly identity) and a
MultiQC custom-content table row. Importable for unit tests (no import-time work).
"""
import argparse


def parse_seqkit(text):
    """Parse `seqkit stats -a -T` text into the fields we report (by column name)."""
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


def render_stats(accession, assembler, m):
    """Tidy machine-readable one-row TSV with identity + metrics."""
    header = ["accession", "assembler", "num_scaffolds",
              "total_length", "n50", "n_bases", "pct_n"]
    row = [accession, assembler,
           str(m["num_seqs"]), str(m["sum_len"]), str(m["n50"]),
           str(m["sum_gap"]), str(pct_n(m["sum_gap"], m["sum_len"]))]
    return "\t".join(header) + "\n" + "\t".join(row) + "\n"


def render_mqc(accession, assembler, m):
    """MultiQC custom-content table row for one assembly variant."""
    sample = f"{accession}-{assembler}"
    lines = [
        "# id: 'assembly_stats'",
        "# section_name: 'Assembly stats'",
        "# plot_type: 'table'",
        "\t".join(["Sample", "Assembler", "Scaffolds",
                   "Total length (bp)", "N50 (bp)", "N (%)"]),
        "\t".join([sample, assembler,
                   str(m["num_seqs"]), str(m["sum_len"]),
                   str(m["n50"]), str(pct_n(m["sum_gap"], m["sum_len"]))]),
    ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seqkit", required=True)
    p.add_argument("--accession", required=True)
    p.add_argument("--assembler", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--mqc", required=True)
    args = p.parse_args(argv)

    with open(args.seqkit) as fh:
        m = parse_seqkit(fh.read())
    with open(args.out, "w") as fh:
        fh.write(render_stats(args.accession, args.assembler, m))
    with open(args.mqc, "w") as fh:
        fh.write(render_mqc(args.accession, args.assembler, m))


if __name__ == "__main__":
    main()
```

**Step 4: Run tests to verify they pass**

Run: `python3 -m unittest tests.test_assembly_stats -v` → expect `Ran 9 tests`, `OK`.
(Full suite `python3 -m unittest discover -s tests` → `Ran 20 tests`, `OK`.)

The `assembly_stats/main.nf` module still passes `--coverage` at this point; that is
only exercised in a real run (never in stub) and is reconciled in Task 5. Stub stays green.

**Step 5: Commit**
```bash
git add tests/test_assembly_stats.py modules/assembly_stats/resources/usr/bin/assembly_stats.py
git commit -m "refactor: drop --coverage from assembly stats summariser"
```

---

### Task 2: Tear down the Supernova assembly subgraph

Reduce the pipeline to a clean gate-zero + full-fetch base. Keep `FETCH_FULL_READS`
(reused). Remove the Supernova processes, params, records, and the assembly publish/output.

**Files:**
- Delete: `modules/supernova.nf`, `modules/stlfr_convert.nf`
- Replace: `main.nf`, `types.nf`
- Modify: `nextflow.config`

**Step 1: Delete the two Supernova modules**
```bash
git rm modules/supernova.nf modules/stlfr_convert.nf
```

**Step 2: Replace `types.nf`** with (gate-zero records only — assembly records are
re-added per task as needed):
```groovy
// Shared record types.
// Included by main.nf and the process modules.

record SraRun {
    accession: String
}

record ReadPair {
    accession: String
    r1: Path
    r2: Path
}

record Result {
    accession: String
    file: Path
}
```

**Step 3: Replace `main.nf`** with:
```groovy
nextflow.enable.types = true
nextflow.enable.moduleBinaries = true

include { FETCH_READS }      from './modules/fetch_reads.nf'
include { FETCH_FULL_READS } from './modules/fetch_full_reads.nf'
include { READ_STRUCTURE }   from './modules/read_structure.nf'
include { BARCODE_CHECK }    from './modules/barcode_check/main.nf'
include { MULTIQC }          from './modules/multiqc.nf'
include { SraRun ; Result } from './types.nf'

// Gate-zero pipeline plus the open-stack assembly baseline (added in later tasks).
params {
    // Run accessions (one element per accession; parallelized). Two on purpose so
    // aggregation steps stay exercised. CBA Captain + Neelam.
    accessions: List<String> = ['SRR32381426', 'SRR32381425']

    // Read pairs to subsample per accession for the gate-zero check.
    check_reads: Integer = 2000000

    // Minimum fraction of reads carrying a valid barcode for the gate to pass.
    pass_fraction: Float = 0.80
}

workflow {
    main:
    runs = channel.fromList(params.accessions)
        .map { acc -> record(accession: acc) }

    reads     = FETCH_READS(runs, params.check_reads)
    structure = READ_STRUCTURE(reads)
    barcode   = BARCODE_CHECK(
        reads,
        file("${projectDir}/assets/stlfr_barcode_whitelist.txt"),
        params.check_reads,
        params.pass_fraction as Float
    )

    // --- Stage 1 (open-stack baseline): full reads fetched here; assembly added next.
    full = FETCH_FULL_READS(runs)

    qc_files = structure
        .map { r -> r.file }
        .mix(barcode.mqc)
        .collect()
    report = MULTIQC(qc_files, channel.value('stLFR gate-zero'))

    publish:
    read_structure = structure
    barcode_check  = barcode.res
    multiqc_report = report.report
}

output {
    read_structure: Channel<Result> {
        path { r -> r.accession }
    }
    barcode_check: Channel<Result> {
        path { r -> r.accession }
    }
    multiqc_report: Path {
        path '.'
    }
}
```

**Step 4: Edit `nextflow.config`** — remove the `stub_local` profile entirely (it capped
`SUPERNOVA`, which no longer exists; re-added for `ABYSS` in Task 4). Leave `local_conda`,
the `manifest`, and the top-level `process { errorStrategy/maxRetries }` untouched.

**Step 5: Validate** (no `-profile` needed — nothing heavy left)

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1b_t2`
Expected: `[SUCCESS] completed=9 failed=0` (7 gate-zero + 2 FETCH_FULL_READS), with only
`read_structure`/`barcode_check`/`multiqc_report` outputs. `modules/assembly_stats/` stays
on disk but is unreferenced (not compiled).

**Step 6: Commit**
```bash
git add -A
git commit -m "refactor: tear down Supernova assembly subgraph to clean base"
```

---

### Task 3: LRTK_CONVERT (stLFR barcode extraction)

**Files:**
- Modify: `types.nf` (add `ConvertedReads`)
- Create: `modules/lrtk_convert.nf`
- Modify: `main.nf` (include + one wiring line)

**Step 1: Add the record** — append to `types.nf`:
```groovy
record ConvertedReads {       // LRTK output: BX:Z:-tagged insert reads
    accession: String
    r1: Path
    r2: Path
}
```

**Step 2 (red): wire into `main.nf`** — add include with the others:
```groovy
include { LRTK_CONVERT } from './modules/lrtk_convert.nf'
```
and after `full = FETCH_FULL_READS(runs)`:
```groovy
    converted = LRTK_CONVERT(full)
```
Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1b_t3`
Expected: FAIL — cannot resolve `./modules/lrtk_convert.nf`.

**Step 3 (green): create `modules/lrtk_convert.nf`**:
```groovy
nextflow.enable.types = true

include { ReadPair ; ConvertedReads } from '../types.nf'

// Extract stLFR barcodes from read2 into BX:Z: tags, emitting barcoded insert reads
// that ABySS (and later Tigmint/ARKS) consume. Runs once per accession; reused by the
// scaffolding stage. Open Bioconda tool — no private image.
//
// NOTE: the exact `lrtk FQCONVER` flags + whether a whitelist arg is needed are
// confirmed on a read subset at smoke time (LRTK is the least-proven tool here). The
// I/O contract (ReadPair in, two BX-tagged FASTQs out) is final.
process LRTK_CONVERT {
    tag "${reads.accession}"
    conda 'bioconda::lrtk=2.0'

    input:
    reads: ReadPair

    output:
    converted: ConvertedReads = record(
        accession: reads.accession,
        r1: file("${reads.accession}.bx_R1.fq.gz"),
        r2: file("${reads.accession}.bx_R2.fq.gz")
    )

    script:
    def acc = reads.accession
    """
    set -euo pipefail
    # TODO(smoke): confirm exact FQCONVER flags + whitelist arg against the tool.
    lrtk FQCONVER -IT stLFR \\
        -I1 ${reads.r1} -I2 ${reads.r2} \\
        -O1 ${acc}.bx_R1.fq.gz -O2 ${acc}.bx_R2.fq.gz
    """

    stub:
    def acc = reads.accession
    """
    printf '@r BX:Z:AAAA-1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}.bx_R1.fq.gz
    printf '@r BX:Z:AAAA-1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}.bx_R2.fq.gz
    """
}
```
Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -output-dir /tmp/s1b_t3`
Expected: `[SUCCESS] completed=11 failed=0` (9 + 2 LRTK_CONVERT). `converted` is unused
for now (consumed in Task 4) — fine.

**Step 4: Commit**
```bash
git add types.nf modules/lrtk_convert.nf main.nf
git commit -m "feat: LRTK_CONVERT stLFR barcode extraction"
```

---

### Task 4: ABYSS contig assembly with k-sweep fan-out

**Files:**
- Modify: `types.nf` (add `AssemblyJob`, `Assembly`)
- Create: `modules/abyss.nf`
- Modify: `main.nf` (params + include + fan-out wiring)
- Modify: `nextflow.config` (re-add `stub_local` capping `ABYSS`)

**Step 1: Add the records** — append to `types.nf`:
```groovy
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
```

**Step 2: Add params to `main.nf`** (after `pass_fraction`):
```groovy
    // Stage 1 baseline: ABySS k-sweep (k must be < the 100 bp read length) and the
    // Bloom-filter size (tune at smoke for the realised k-mer count).
    abyss_kmers: List<Integer> = [64, 80, 96]
    abyss_bloom: String = '20G'
```

**Step 3: Re-add the `stub_local` profile to `nextflow.config`** (alongside `local_conda`):
```groovy
    // Cap the assembly step's resource request so the whole pipeline can be validated
    // locally in stub mode. ABYSS requests 16 CPU / 64 GB for production; this withName
    // selector overrides that only when the profile is active. Production runs (no
    // profile) keep the real request. Resource directives are not part of the task hash.
    stub_local {
        process {
            withName: ABYSS {
                cpus   = 2
                memory = 2.GB
            }
        }
    }
```

**Step 4 (red): wire the fan-out into `main.nf`** — add include:
```groovy
include { ABYSS } from './modules/abyss.nf'
```
and after `converted = LRTK_CONVERT(full)`:
```groovy
    kmers   = channel.fromList(params.abyss_kmers)
    jobs    = converted.combine(kmers).map { conv, k ->
        record(accession: conv.accession, r1: conv.r1, r2: conv.r2, k: k)
    }
    contigs = ABYSS(jobs)
```
Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -profile stub_local -output-dir /tmp/s1b_t4`
Expected: FAIL — cannot resolve `./modules/abyss.nf`.

> If the stub fails on the `map { conv, k -> record(...) }` line with a record-inference
> error, the simple form worked for the analogous Supernova fan-out, so first re-verify
> the typed-syntax exactly; only if it genuinely won't infer, annotate
> `AssemblyJob j = record(...); j` (and add `AssemblyJob` to the `from './types.nf'`
> include). Never use `new AssemblyJob(...)`. If neither works, STOP and report the error.

**Step 5 (green): create `modules/abyss.nf`**:
```groovy
nextflow.enable.types = true

include { AssemblyJob ; Assembly } from '../types.nf'

// De novo contig assembly with ABySS (Bloom-filter mode for memory efficiency on a
// ~740 Mb genome). One job per (accession x k); the k-sweep is the fan-out. Plain
// paired-end de Bruijn assembly — barcodes are not used until the scaffolding stage.
//
// NOTE: confirm the contigs output filename and tune `abyss_bloom` on the smoke run.
process ABYSS {
    tag "${job.accession}-k${job.k}"
    conda 'bioconda::abyss=2.3.10'
    cpus 16
    memory 64.GB

    input:
    job: AssemblyJob

    output:
    assembly: Assembly = record(
        accession: job.accession,
        assembler: "abyss-k${job.k}",
        fasta: file("${job.accession}.abyss-k${job.k}.contigs.fa.gz")
    )

    script:
    def acc = job.accession
    def k   = job.k
    """
    set -euo pipefail
    abyss-pe name=${acc}_k${k} k=${k} j=${task.cpus} \\
        B=${params.abyss_bloom} in='${job.r1} ${job.r2}' \\
        ${acc}_k${k}-contigs.fa
    gzip -c ${acc}_k${k}-contigs.fa > ${acc}.abyss-k${k}.contigs.fa.gz
    """

    stub:
    def acc = job.accession
    def k   = job.k
    """
    printf '>contig_1\\nACGTACGTACGT\\n' | gzip -c > ${acc}.abyss-k${k}.contigs.fa.gz
    """
}
```
Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -profile stub_local -output-dir /tmp/s1b_t4`
Expected: `[SUCCESS] completed=17 failed=0` with **six** ABYSS tasks
(`SRR32381426-k64/k80/k96`, `SRR32381425-k64/k80/k96`). `contigs` is unused until Task 5.

**Step 6: Commit**
```bash
git add types.nf modules/abyss.nf main.nf nextflow.config
git commit -m "feat: ABYSS contig assembly with k-sweep fan-out"
```

---

### Task 5: Rewire ASSEMBLY_STATS and publish

**Files:**
- Modify: `types.nf` (add `AssemblyStats`)
- Replace: `modules/assembly_stats/main.nf`
- Modify: `main.nf` (include, wiring, MultiQC mix, publish, output)

**Step 1: Add the record** — append to `types.nf`:
```groovy
record AssemblyStats {
    accession: String
    assembler: String
    file: Path
}
```

**Step 2: Replace `modules/assembly_stats/main.nf`** with (new record shapes, no coverage):
```groovy
nextflow.enable.types = true

include { Assembly ; AssemblyStats } from '../../types.nf'

// Contiguity summary for one assembly: scaffold count, total length, N50, %N. Emits a
// machine-readable stats.tsv (AssemblyStats record) plus a MultiQC custom-content row.
// Summariser is a module-scoped binary (resources/usr/bin, on PATH).
process ASSEMBLY_STATS {
    tag "${assembly.accession}-${assembly.assembler}"
    conda 'bioconda::seqkit=2.8.2 conda-forge::python=3.12'

    input:
    assembly: Assembly

    output:
    res: AssemblyStats = record(
        accession: assembly.accession,
        assembler: assembly.assembler,
        file: file("${assembly.accession}.${assembly.assembler}.stats.tsv")
    )
    mqc: Path = file("${assembly.accession}.${assembly.assembler}.assembly_mqc.tsv")

    script:
    def acc = assembly.accession
    def asm = assembly.assembler
    """
    seqkit stats -a -T ${assembly.fasta} > seqkit.tsv
    assembly_stats.py \\
        --seqkit seqkit.tsv \\
        --accession ${acc} \\
        --assembler ${asm} \\
        --out ${acc}.${asm}.stats.tsv \\
        --mqc ${acc}.${asm}.assembly_mqc.tsv
    """

    stub:
    def acc = assembly.accession
    def asm = assembly.assembler
    """
    printf 'accession\\tassembler\\tn50\\n${acc}\\t${asm}\\t23000\\n' > ${acc}.${asm}.stats.tsv
    printf "# id: 'assembly_stats'\\n# plot_type: 'table'\\nSample\\tN50 (bp)\\n${acc}-${asm}\\t23000\\n" > ${acc}.${asm}.assembly_mqc.tsv
    """
}
```

**Step 3: Wire into `main.nf`** — add include:
```groovy
include { ASSEMBLY_STATS } from './modules/assembly_stats/main.nf'
```
add `Assembly ; AssemblyStats` to the types include (so the `output {}` block can name them):
```groovy
include { SraRun ; Result ; Assembly ; AssemblyStats } from './types.nf'
```
after `contigs = ABYSS(jobs)` add:
```groovy
    asm_stats = ASSEMBLY_STATS(contigs)
```
change the `qc_files` collect to also mix the assembly tables:
```groovy
    qc_files = structure
        .map { r -> r.file }
        .mix(barcode.mqc)
        .mix(asm_stats.mqc)
        .collect()
```
add to `publish:`:
```groovy
    assembly       = contigs
    assembly_stats = asm_stats.res
```
add to `output {}`:
```groovy
    assembly: Channel<Assembly> {
        path { a -> "assemblies/${a.accession}/${a.assembler}" }
    }
    assembly_stats: Channel<AssemblyStats> {
        path { s -> "assemblies/${s.accession}/${s.assembler}" }
    }
```

**Step 4: Validate**

Run: `NXF_SYNTAX_PARSER=v2 nextflow run . -stub -profile stub_local -output-dir /tmp/s1b_t5`
Expected: `[SUCCESS] completed=23 failed=0` (17 + 6 ASSEMBLY_STATS). `Outputs:` lists
`assembly` and `assembly_stats`. Then:
`find /tmp/s1b_t5/assemblies -type l | sort` → 12 entries: each of SRR32381426/SRR32381425
× abyss-k64/k80/k96 has `.contigs.fa.gz` and `.stats.tsv` under
`assemblies/<accession>/abyss-k<k>/`.

**Step 5: Commit**
```bash
git add types.nf modules/assembly_stats/main.nf main.nf
git commit -m "feat: rewire ASSEMBLY_STATS for ABySS contigs + publish + MultiQC merge"
```

---

### Task 6: Full-DAG verification + unit-test sweep (verification only)

**Step 1:** `rm -rf /tmp/s1b_final && NXF_SYNTAX_PARSER=v2 nextflow run . -stub -profile stub_local -output-dir /tmp/s1b_final`
→ expect `[SUCCESS] completed=23 failed=0`.

**Step 2:** `find /tmp/s1b_final/assemblies -type l | sort` → 12 artifacts under
`assemblies/<accession>/abyss-k<k>/` (2 accessions × 3 k × {contigs.fa.gz, stats.tsv}).

**Step 3:** `python3 -m unittest discover -s tests -v` → `Ran 20 tests`, `OK` (11
barcode-check + 9 assembly-stats).

No commit (verification only).

---

### Task 7: Update docs and journal

**Files:**
- Create: `docs/2026-06-18-stage1b-abyss-contigs-build.md` (build journal entry)
- Modify: `docs/README.md` (register the build entry + this plan)
- Modify: `AGENTS.md` (replace the Supernova Stage-1 working-context section with the
  open-stack one: new params `abyss_kmers`/`abyss_bloom`, the `conda`-only no-image win,
  the `stub_local` cap on `ABYSS`, the LRTK `TODO(smoke)`, anchoring contig N50 to the
  paper's published numbers, and the scaffolding chain as the next stage; update the
  Planned-stages baseline note to point at the open stack)

Record the first-real-run prerequisites (much smaller now): confirm `LRTK FQCONVER` flags
on a subset, tune `abyss_bloom`, confirm ABySS resource envelope and contigs filename. No
private image needed.

**Commit:**
```bash
git add docs/ AGENTS.md
git commit -m "docs: record open ABySS contigs baseline build + refresh AGENTS.md"
```

---

## First real run (after the smoke checks, separate from this build)

Relaunch the Launchpad pipeline pinned to this branch's commit so gate zero + the full
fetch stay cached and only the new tasks run:
```bash
export TOWER_WORKSPACE_ID=66005590816034
tw runs relaunch -i <gate-zero-runId> --revision stage1b-abyss-contigs \
  --commit-id $(git rev-parse HEAD)
```
Watch on the smoke run: the `LRTK FQCONVER` output on a subset (read count, BX tags
present), ABySS contig N50 per k vs the paper's published Supernova contig N50 (17–34 kb),
and the ABySS memory envelope at 220× (tune `abyss_bloom`).
