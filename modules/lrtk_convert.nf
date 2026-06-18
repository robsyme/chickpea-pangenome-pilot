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
