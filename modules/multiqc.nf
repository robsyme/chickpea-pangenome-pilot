nextflow.enable.types = true

// Aggregate per-accession QC into a single report: seqkit read-structure stats
// (native MultiQC module) plus the stLFR barcode-integrity custom-content table.
process MULTIQC {
    conda 'bioconda::multiqc=1.25'

    input:
    qc_files: Bag<Path>
    title: String

    output:
    report: Path = file('multiqc_report.html')
    data: Path = file('multiqc_data')

    script:
    """
    multiqc --force --title '${title}' --filename multiqc_report.html .
    """

    stub:
    """
    mkdir -p multiqc_data
    echo '<html><body>stub</body></html>' > multiqc_report.html
    """
}
