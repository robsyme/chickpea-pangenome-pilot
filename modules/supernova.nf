nextflow.enable.types = true

include { AssemblyJob ; Assembly } from '../types.nf'

// Assemble one (accession x coverage) job with Supernova. Sizes --maxreads from
// the target coverage, genome size, and measured read length, runs `supernova run`,
// and extracts a pseudohap scaffold FASTA with `mkoutput`.
//
// NOTE: Supernova is a proprietary 10x binary baked into a private image
// (params.supernova_image). The run/mkoutput flags are standard but confirm the
// mkoutput output filename + gzip behaviour on the smoke run; the I/O contract and
// maxreads arithmetic below are final.
process SUPERNOVA {
    tag "${job.accession}-${job.coverage}x"
    container params.supernova_image
    cpus 16
    memory 250.GB

    input:
    job: AssemblyJob
    genome_size: Integer

    output:
    assembly: Assembly = record(
        accession: job.accession,
        assembler: 'supernova',
        coverage: job.coverage,
        fasta: file("${job.accession}.${job.coverage}x.scaffolds.fa.gz"),
        summary: file("${job.accession}.${job.coverage}x.summary.csv")
    )

    script:
    def acc = job.accession
    def cov = job.coverage
    def id  = "${acc}_${cov}x"
    """
    set -euo pipefail

    # Size --maxreads to hit the requested raw coverage:
    #   maxreads = genome_size * coverage / read_len
    # read_len is measured from the converted reads (not assumed).
    one=\$(ls ${job.dir}/*.fastq.gz | head -1)
    read_len=\$(gzip -dc "\$one" | awk 'NR==2{print length(\$0); exit}')
    [ "\$read_len" -gt 0 ] || { echo "ERROR: could not measure read length" >&2; exit 1; }
    maxreads=\$(awk -v g=${genome_size} -v c=${cov} -v l=\$read_len 'BEGIN{printf "%d", (g*c)/l}')
    echo "[SUPERNOVA] ${acc} target=${cov}x read_len=\$read_len maxreads=\$maxreads" >&2

    # Leave ~6 GB headroom for the JVM/OS below the instance's memory request.
    localmem=\$(( ${task.memory.toGiga()} - 6 ))
    # Floor it so a reduced-memory invocation never passes Supernova a negative value.
    [ "\$localmem" -lt 1 ] && localmem=1

    supernova run \\
        --id=${id} \\
        --fastqs=${job.dir} \\
        --maxreads=\$maxreads \\
        --localcores=${task.cpus} \\
        --localmem=\$localmem

    # TODO(smoke): confirm mkoutput writes <outprefix>.fasta.gz (then drop the gzip).
    supernova mkoutput \\
        --style=pseudohap \\
        --asmdir=${id}/outs/assembly \\
        --outprefix=${acc}.${cov}x.scaffolds
    [ -f ${acc}.${cov}x.scaffolds.fasta.gz ] || gzip -f ${acc}.${cov}x.scaffolds.fasta
    mv ${acc}.${cov}x.scaffolds.fasta.gz ${acc}.${cov}x.scaffolds.fa.gz 2>/dev/null || \\
        mv ${acc}.${cov}x.scaffolds.fasta ${acc}.${cov}x.scaffolds.fa.gz

    cp ${id}/outs/summary.csv ${acc}.${cov}x.summary.csv
    """

    stub:
    def acc = job.accession
    def cov = job.coverage
    """
    printf '>scaffold_1\\nACGTACGTACGT\\n' | gzip -c > ${acc}.${cov}x.scaffolds.fa.gz
    printf 'metric,value\\nest_genome_size,740000000\\n' > ${acc}.${cov}x.summary.csv
    """
}
