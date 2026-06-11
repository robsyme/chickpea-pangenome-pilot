nextflow.enable.types = true

include { ReadPair ; Result } from '../../types.nf'

// Gate zero: verify the stLFR barcodes are intact in read2. Detects the barcode
// offsets empirically, then reports the fraction of reads carrying a valid
// 3-barcode combination, the distinct-barcode count, and a reads-per-barcode
// histogram. Emits both the JSON result and a MultiQC custom-content table row.
// The check script is a module-scoped binary (resources/usr/bin, on PATH).
process BARCODE_CHECK {
    tag "${reads.accession}"
    conda 'conda-forge::python=3.12'

    input:
    reads: ReadPair
    whitelist: Path
    n_reads: Integer
    pass_fraction: Float

    output:
    res: Result = record(
        accession: reads.accession,
        file: file("${reads.accession}.barcode_check.json")
    )
    mqc: Path = file("${reads.accession}.barcode_mqc.tsv")

    script:
    """
    stlfr_barcode_check.py \\
        --r1 ${reads.r1} \\
        --r2 ${reads.r2} \\
        --whitelist ${whitelist} \\
        --reads ${n_reads} \\
        --pass-fraction ${pass_fraction} \\
        --sample ${reads.accession} \\
        --out ${reads.accession}.barcode_check.json \\
        --mqc ${reads.accession}.barcode_mqc.tsv
    """

    stub:
    """
    echo '{"gate_pass": true, "valid_fraction": 1.0, "stub": true}' > ${reads.accession}.barcode_check.json
    printf "# id: 'stlfr_barcode'\\n# plot_type: 'table'\\nSample\\tGate\\n${reads.accession}\\tPASS\\n" > ${reads.accession}.barcode_mqc.tsv
    """
}
