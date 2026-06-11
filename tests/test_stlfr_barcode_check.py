"""Tests for bin/stlfr_barcode_check.py (stdlib unittest, no pytest needed).

Run: python3 -m unittest discover -s tests -v
"""
import importlib.util
import pathlib
import unittest

_BIN = pathlib.Path(__file__).resolve().parent.parent / "bin" / "stlfr_barcode_check.py"
_spec = importlib.util.spec_from_file_location("sbc", _BIN)
sbc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sbc)

# Four real whitelist entries (1-based indices 1..4 in the vendored list).
WL = ["TAACAGCCAA", "CTAAGAGTCC", "TTACTGCCTT", "CGCTGAATTC"]

# Real stLFR layout observed in SRR32381426 read2: barcodes at offsets 100/116/132
# (10 + 6-bp spacer + 10 + 6-bp spacer + 10), read length 142.
REAL_OFFSETS = [100, 116, 132]
REAL_LEN = 142


def _make_r2(offsets, barcodes, length, fill="N"):
    s = list(fill * length)
    for off, bc in zip(offsets, barcodes):
        s[off:off + 10] = list(bc)
    return "".join(s)


class TestMatchBarcode(unittest.TestCase):
    def setUp(self):
        self.ftm = sbc.build_fault_tolerant_map(WL)

    def test_exact_match_returns_one_based_index(self):
        self.assertEqual(sbc.match_barcode("TAACAGCCAA", self.ftm), 1)
        self.assertEqual(sbc.match_barcode("CTAAGAGTCC", self.ftm), 2)

    def test_one_mismatch_matches(self):
        self.assertEqual(sbc.match_barcode("AAACAGCCAA", self.ftm), 1)

    def test_two_mismatches_returns_none(self):
        self.assertIsNone(sbc.match_barcode("AAGCAGCCAA", self.ftm))

    def test_unknown_returns_none(self):
        self.assertIsNone(sbc.match_barcode("GGGGGGGGGG", self.ftm))


class TestExtractTriplet(unittest.TestCase):
    def test_extracts_at_given_offsets(self):
        r2 = _make_r2(REAL_OFFSETS, ["AAAAAAAAAA", "CCCCCCCCCC", "TTTTTTTTTT"], REAL_LEN)
        b1, b2, b3 = sbc.extract_triplet(r2, REAL_OFFSETS)
        self.assertEqual(b1, "AAAAAAAAAA")
        self.assertEqual(b2, "CCCCCCCCCC")
        self.assertEqual(b3, "TTTTTTTTTT")

    def test_short_read_returns_none_triplet(self):
        b1, _, _ = sbc.extract_triplet("ACGT" * 5, REAL_OFFSETS)
        self.assertIsNone(b1)


class TestDetectBarcodeOffsets(unittest.TestCase):
    def setUp(self):
        self.ftm = sbc.build_fault_tolerant_map(WL)

    def test_finds_three_real_offsets(self):
        # 200 reads with whitelist barcodes at 100/116/132; "A"-filled elsewhere.
        reads = [_make_r2(REAL_OFFSETS, [WL[0], WL[1], WL[2]], REAL_LEN, fill="A")
                 for _ in range(200)]
        offsets = sbc.detect_barcode_offsets(reads, self.ftm)
        self.assertEqual(offsets, REAL_OFFSETS)


class TestAnalyze(unittest.TestCase):
    def setUp(self):
        self.ftm = sbc.build_fault_tolerant_map(WL)

    def test_reports_valid_fraction_and_distinct_barcodes(self):
        reads = [
            _make_r2(REAL_OFFSETS, [WL[0], WL[1], WL[2]], REAL_LEN),  # barcode A
            _make_r2(REAL_OFFSETS, [WL[0], WL[1], WL[2]], REAL_LEN),  # barcode A again
            _make_r2(REAL_OFFSETS, [WL[0], WL[1], WL[3]], REAL_LEN),  # barcode B
            _make_r2(REAL_OFFSETS, ["GGGGGGGGGG", WL[1], WL[2]], REAL_LEN),  # invalid b1
        ]
        res = sbc.analyze(reads, self.ftm, REAL_OFFSETS)
        self.assertEqual(res["n_total"], 4)
        self.assertEqual(res["n_valid"], 3)
        self.assertAlmostEqual(res["valid_fraction"], 0.75, places=6)
        self.assertEqual(res["distinct_barcodes"], 2)


class TestRenderMqc(unittest.TestCase):
    def setUp(self):
        self.res = {
            "barcode_offsets": [100, 116, 132],
            "n_total": 2000000,
            "n_valid": 1805844,
            "valid_fraction": 0.902922,
            "distinct_barcodes": 1657967,
            "r2_len_mode": 142,
            "gate_pass": True,
        }

    def test_has_custom_content_config_header(self):
        out = sbc.render_mqc("SRR32381426", self.res)
        self.assertIn("# id: 'stlfr_barcode'", out)
        self.assertIn("plot_type: 'table'", out)
        self.assertIn("section_name:", out)

    def test_data_row_keyed_by_sample_with_values(self):
        out = sbc.render_mqc("SRR32381426", self.res)
        # Locate the single data row (last non-empty, non-comment line).
        rows = [ln for ln in out.splitlines() if ln and not ln.startswith("#")]
        header, data = rows[0], rows[1]
        self.assertEqual(header.split("\t")[0], "Sample")
        cells = data.split("\t")
        self.assertEqual(cells[0], "SRR32381426")
        # valid percentage rendered to one decimal
        self.assertIn("90.3", data)
        self.assertIn("1657967", data)
        self.assertIn("100;116;132", data)
        self.assertIn("PASS", data)

    def test_gate_fail_renders_fail(self):
        res = dict(self.res, gate_pass=False)
        out = sbc.render_mqc("S2", res)
        self.assertIn("FAIL", out)


if __name__ == "__main__":
    unittest.main()
