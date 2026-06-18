nextflow.enable.types = true

include { SraRun ; ReadPair } from '../types.nf'

// Download the FULL run for an accession (no subsample) — needed for assembly.
// Separate from FETCH_READS (which subsamples for gate zero) so the gate-zero
// task hash stays stable and its cache survives. fastq-dl auto-falls back ENA->SRA.
process FETCH_FULL_READS {
    tag "${run.accession}"
    container 'community.wave.seqera.io/library/fastq-dl:3.0.1--fa446f61dfc85bc3'

    input:
    run: SraRun

    output:
    reads: ReadPair = record(
        accession: run.accession,
        r1: file("${run.accession}_R1.fastq.gz"),
        r2: file("${run.accession}_R2.fastq.gz")
    )

    script:
    def acc = run.accession
    """
    set -euo pipefail
    rm -rf dl && mkdir dl
    fastq-dl --accession ${acc} --outdir dl >&2
    mv "\$(ls dl/*_1.fastq.gz | head -1)" ${acc}_R1.fastq.gz
    mv "\$(ls dl/*_2.fastq.gz | head -1)" ${acc}_R2.fastq.gz
    rm -rf dl
    """

    stub:
    def acc = run.accession
    """
    printf '@r/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}_R1.fastq.gz
    printf '@r/2\\nACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIII\\n' | gzip -c > ${acc}_R2.fastq.gz
    """
}
