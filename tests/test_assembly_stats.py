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
