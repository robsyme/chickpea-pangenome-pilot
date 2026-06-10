# Vendored assets provenance

## stlfr_barcode_whitelist.txt

The canonical stLFR barcode whitelist: 1536 unique 10 bp barcode sequences, each
with a 1-based index, tab-separated.

- Source: `BGI-Qingdao/stLFR_barcode_split`, file `barcode_list.txt`
  https://raw.githubusercontent.com/BGI-Qingdao/stLFR_barcode_split/master/barcode_list.txt
- Retrieved: 2026-06-10
- Unmodified copy (1536 lines).

This is the same whitelist used by BGI's `split_barcode.pl`, which reads three
10 bp barcodes from read2 with the segment layout `(n1,n2,n3,n4,n5) = (10,6,10,0,10)`
(barcode, 6 bp spacer, barcode, barcode) following the genomic insert, matched with
1-mismatch fault tolerance.
