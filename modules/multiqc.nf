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
    // MultiQC names the data dir after --filename: multiqc_report -> multiqc_report_data
    data: Path = file('multiqc_report_data')

    script:
    """
    multiqc --force --title '${title}' --filename multiqc_report.html .
    """

    stub:
    """
    mkdir -p multiqc_report_data
    echo '<html><body>stub</body></html>' > multiqc_report.html
    """
}
