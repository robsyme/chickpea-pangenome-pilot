nextflow.enable.types = true

include { ReadPair ; Result } from '../types.nf'

// Read-length and base-composition summary for R1 and R2. Confirms the expected
// stLFR asymmetry (R1 ~100 bp insert, R2 ~142 bp = insert + barcode block).
process READ_STRUCTURE {
    tag "${reads.accession}"
    conda 'bioconda::seqkit=2.8.2'

    input:
    reads: ReadPair

    output:
    res: Result = record(
        accession: reads.accession,
        file: file("${reads.accession}.seqkit_stats.tsv")
    )

    script:
    """
    seqkit stats -a -T ${reads.r1} ${reads.r2} > ${reads.accession}.seqkit_stats.tsv
    """

    stub:
    """
    printf 'file\\tnum_seqs\\tavg_len\\n${reads.r1}\\t1\\t100\\n${reads.r2}\\t1\\t142\\n' > ${reads.accession}.seqkit_stats.tsv
    """
}
