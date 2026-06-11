# chickpea-pangenome-pilot

A Nextflow pipeline used as the executable record of a pilot investigating whether
the stLFR linked-read data behind Garg et al. (2025),
[*An Australian chickpea pan-genome*](https://onlinelibrary.wiley.com/doi/10.1111/pbi.70192),
can be reassembled de novo to address artefacts in the published reference-guided
assemblies.

This first stage is **gate zero**: confirm the stLFR barcodes are intact in read2 of
the deposited data (NCBI BioProject PRJNA1225167) before committing to reassembly.
For each run accession it subsamples reads from ENA, summarises read structure with
`seqkit`, and verifies barcode integrity against the BGI stLFR whitelist.

The barcode check (`modules/barcode_check/resources/usr/bin/stlfr_barcode_check.py`,
a module-scoped binary) detects where the three 10 bp
barcodes sit in read2 empirically rather than assuming fixed offsets, because the
deposited layout (barcodes at 100 / 116 / 132) differs from BGI's `split_barcode.pl`
default. On the pilot cultivar (CBA Captain, SRR32381426) it found 87.6% of reads
carry a valid barcode — the data is usable as linked-reads.

## Requirements

- Nextflow ≥ 26.04
- `NXF_SYNTAX_PARSER=v2` (the typed syntax needs the strict v2 parser; this is an
  environment variable, not a config setting)

## Run

```bash
export NXF_SYNTAX_PARSER=v2
nextflow run . \
    --accessions SRR32381426 \
    --check_reads 2000000 \
    --pass_fraction 0.80
```

`--accessions` takes a list; each accession is processed in parallel. Results are
published per accession under `results/<accession>/`.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `accessions` | `['SRR32381426']` | SRA/ENA run accessions to check |
| `check_reads` | `2000000` | Read pairs to subsample per accession |
| `pass_fraction` | `0.80` | Minimum valid-barcode fraction for the gate to pass |

## Tests

```bash
python3 -m unittest discover -s tests -v
```

See [`AGENTS.md`](AGENTS.md) for the full project context, infrastructure, and the
validated typed-syntax notes.
