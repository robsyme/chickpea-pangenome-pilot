#!/usr/bin/env python3
"""stLFR barcode-integrity check.

Verifies that read2 of an stLFR library carries the three 10 bp barcodes the BGI
chemistry encodes after the genomic insert, using the canonical whitelist and the
`split_barcode.pl` layout: insert + b1(10) + spacer(6) + b2(10) + b3(10), matched
with 1-mismatch fault tolerance.

The insert length (and thus barcode offsets) is detected empirically by scanning
candidate lengths, because deposited read2 lengths do not always match the textbook
136/146 bp layouts.

Pure stdlib. Usable as a library (see tests) or as a CLI (see main()).
"""
import argparse
import gzip
import io
import json
import sys
from collections import Counter

_BASES = "ACGT"


def _one_mismatch_variants(seq):
    """Yield every single-substitution variant of seq (excludes seq itself)."""
    for i, base in enumerate(seq):
        for alt in _BASES:
            if alt != base:
                yield seq[:i] + alt + seq[i + 1:]


def build_fault_tolerant_map(whitelist_seqs):
    """Map barcode sequences (exact + unambiguous 1-mismatch) to 1-based indices.

    Exact sequences always win. A 1-mismatch variant is kept only if it resolves to
    exactly one barcode and is not itself an exact whitelist sequence.
    """
    exact = {}
    for idx, seq in enumerate(whitelist_seqs, start=1):
        exact[seq] = idx

    variant_owners = {}  # variant -> set of barcode indices
    for idx, seq in enumerate(whitelist_seqs, start=1):
        for variant in _one_mismatch_variants(seq):
            if variant in exact:
                continue
            variant_owners.setdefault(variant, set()).add(idx)

    ft_map = dict(exact)
    for variant, owners in variant_owners.items():
        if len(owners) == 1:
            ft_map[variant] = next(iter(owners))
    return ft_map


def match_barcode(seq, ft_map):
    """Return the 1-based barcode index for seq, or None if unmatched/ambiguous."""
    return ft_map.get(seq)


def extract_triplet(r2_seq, offsets, bc_len=10):
    """Return the 10 bp barcode substrings of read2 at the given start offsets.

    Returns a tuple of None if read2 is too short to contain any requested barcode.
    The number of offsets determines the number of barcodes returned.
    """
    out = []
    for off in offsets:
        if len(r2_seq) < off + bc_len:
            return tuple(None for _ in offsets)
        out.append(r2_seq[off:off + bc_len])
    return tuple(out)


def _valid_triplet_index(r2_seq, ft_map, offsets):
    parts = extract_triplet(r2_seq, offsets)
    if parts[0] is None:
        return None
    idx = tuple(match_barcode(p, ft_map) for p in parts)
    if all(idx):
        return idx
    return None


def detect_barcode_offsets(r2_seqs, ft_map, n_barcodes=3, bc_len=10, min_sep=10):
    """Find where the barcodes sit in read2 by scanning every window position.

    Computes the whitelist match rate for each 10 bp window, then greedily selects
    the `n_barcodes` highest-rate, non-overlapping positions. This makes the check
    robust to chemistry/spacer differences instead of assuming fixed offsets.
    """
    seqs = list(r2_seqs)
    if not seqs:
        return []
    total = len(seqs)
    read_len = min(len(s) for s in seqs)
    rates = []
    for start in range(0, read_len - bc_len + 1):
        hits = sum(1 for s in seqs if match_barcode(s[start:start + bc_len], ft_map) is not None)
        rates.append((start, hits / total))
    rates.sort(key=lambda x: -x[1])
    chosen = []
    for start, _rate in rates:
        if all(abs(start - c) >= min_sep for c in chosen):
            chosen.append(start)
        if len(chosen) == n_barcodes:
            break
    return sorted(chosen)


def analyze(r2_seqs, ft_map, offsets):
    """Summarise barcode validity for a sample of read2 sequences at given offsets."""
    n_total = 0
    n_valid = 0
    barcode_counts = Counter()
    for s in r2_seqs:
        n_total += 1
        bc = _valid_triplet_index(s, ft_map, offsets)
        if bc is not None:
            n_valid += 1
            barcode_counts[bc] += 1
    reads_per_barcode = Counter(barcode_counts.values())
    return {
        "barcode_offsets": list(offsets),
        "n_total": n_total,
        "n_valid": n_valid,
        "valid_fraction": (n_valid / n_total) if n_total else 0.0,
        "distinct_barcodes": len(barcode_counts),
        "reads_per_barcode_histogram": dict(sorted(reads_per_barcode.items())),
    }


# --- CLI -------------------------------------------------------------------

def load_whitelist(path):
    seqs = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                seqs.append(line.split()[0].upper())
    return seqs


def _open_maybe_gzip(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path)


def iter_fastq_seqs(path, limit=None):
    n = 0
    with _open_maybe_gzip(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # plus
            fh.readline()  # qual
            yield seq
            n += 1
            if limit and n >= limit:
                break


def main(argv=None):
    p = argparse.ArgumentParser(description="stLFR barcode-integrity check on read2")
    p.add_argument("--r1", required=True, help="read1 fastq(.gz)")
    p.add_argument("--r2", required=True, help="read2 fastq(.gz) carrying barcodes")
    p.add_argument("--whitelist", required=True, help="stLFR barcode whitelist")
    p.add_argument("--reads", type=int, default=2_000_000, help="reads to sample")
    p.add_argument("--detect-sample", type=int, default=20_000,
                   help="reads used to locate barcode offsets")
    # A healthy stLFR library recovers a barcode (within 1 mismatch) for ~85-90%
    # of reads; the rest carry >1 error in a barcode. 0.80 passes healthy data and
    # still catches stripped/mislaid barcodes (which score near 0% or <10%).
    p.add_argument("--pass-fraction", type=float, default=0.80)
    p.add_argument("--out", default="barcode_check.json")
    args = p.parse_args(argv)

    wl = load_whitelist(args.whitelist)
    ft_map = build_fault_tolerant_map(wl)

    r1_lens = [len(s) for s in iter_fastq_seqs(args.r1, limit=min(args.reads, 100_000))]
    r2_seqs = list(iter_fastq_seqs(args.r2, limit=args.reads))
    r2_lens = [len(s) for s in r2_seqs]

    offsets = detect_barcode_offsets(r2_seqs[:args.detect_sample], ft_map)
    res = analyze(r2_seqs, ft_map, offsets)

    res["whitelist_size"] = len(wl)
    res["r1_len_mode"] = Counter(r1_lens).most_common(1)[0][0] if r1_lens else None
    res["r2_len_mode"] = Counter(r2_lens).most_common(1)[0][0] if r2_lens else None
    res["pass_fraction_threshold"] = args.pass_fraction
    res["gate_pass"] = bool(res["valid_fraction"] >= args.pass_fraction)

    with open(args.out, "w") as fh:
        json.dump(res, fh, indent=2)

    verdict = "PASS" if res["gate_pass"] else "FAIL"
    print(f"[stLFR barcode check] {verdict}: "
          f"{res['valid_fraction']:.1%} of {res['n_total']} reads carry a valid barcode; "
          f"{res['distinct_barcodes']} distinct barcodes; "
          f"barcode offsets={offsets}, R1={res['r1_len_mode']}, R2={res['r2_len_mode']}",
          file=sys.stderr)
    return 0 if res["gate_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
