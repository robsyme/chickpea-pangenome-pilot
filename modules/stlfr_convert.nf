nextflow.enable.types = true

include { ReadPair ; ConvertedReads } from '../types.nf'

// Convert MGI stLFR co-barcoded reads into Chromium-style FASTQs that Supernova
// accepts (stlfr2supernova steps 0-2). Coverage-independent and expensive, so it
// runs ONCE per accession; both coverage variants assemble from this one set.
//
// NOTE: the exact in-container CLI for running ONLY the conversion (stopping
// before Supernova) must be confirmed against the stlfr2supernova_pipeline scripts
// when the image is built. The wiring and I/O contract below are final; the script
// body is the best-effort invocation to finalise at smoke-test time.
process STLFR_CONVERT {
    tag "${reads.accession}"
    container params.stlfr2supernova_image

    input:
    reads: ReadPair

    output:
    converted: ConvertedReads = record(
        accession: reads.accession,
        dir: file("${reads.accession}_chromium")
    )

    script:
    def acc = reads.accession
    """
    set -euo pipefail
    mkdir -p ${acc}_chromium
    # TODO(smoke): confirm exact stlfr2supernova conversion-only invocation.
    stlfr2supernova_convert \\
        --r1 ${reads.r1} \\
        --r2 ${reads.r2} \\
        --outdir ${acc}_chromium
    """

    stub:
    def acc = reads.accession
    """
    mkdir -p ${acc}_chromium
    printf '@r 1:N:0:0\\nACGT\\n+\\nIIII\\n' | gzip -c \\
        > ${acc}_chromium/read-RA_si-AAAA_lane-001-chunk-001.fastq.gz
    """
}
