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
