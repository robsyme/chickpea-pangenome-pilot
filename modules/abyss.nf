nextflow.enable.types = true

include { AssemblyJob ; Assembly } from '../types.nf'

// De novo contig assembly with ABySS (Bloom-filter mode for memory efficiency on a
// ~740 Mb genome). One job per (accession x k); the k-sweep is the fan-out. Plain
// paired-end de Bruijn assembly — barcodes are not used until the scaffolding stage.
//
// NOTE: confirm the contigs output filename and tune `abyss_bloom` on the smoke run.
process ABYSS {
    tag "${job.accession}-k${job.k}"
    conda 'bioconda::abyss=2.3.10'
    cpus 16
    memory 64.GB

    input:
    job: AssemblyJob

    output:
    assembly: Assembly = record(
        accession: job.accession,
        assembler: "abyss-k${job.k}",
        fasta: file("${job.accession}.abyss-k${job.k}.contigs.fa.gz")
    )

    script:
    def acc = job.accession
    def k   = job.k
    """
    set -euo pipefail
    abyss-pe name=${acc}_k${k} k=${k} j=${task.cpus} \\
        B=${params.abyss_bloom} in='${job.r1} ${job.r2}' \\
        ${acc}_k${k}-contigs.fa
    gzip -c ${acc}_k${k}-contigs.fa > ${acc}.abyss-k${k}.contigs.fa.gz
    """

    stub:
    def acc = job.accession
    def k   = job.k
    """
    printf '>contig_1\\nACGTACGTACGT\\n' | gzip -c > ${acc}.abyss-k${k}.contigs.fa.gz
    """
}
