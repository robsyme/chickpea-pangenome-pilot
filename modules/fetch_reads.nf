nextflow.enable.types = true

include { SraRun ; ReadPair } from '../types.nf'

// Resolve the ENA fastq_ftp paths for a run accession and stream the first
// `n_reads` read pairs of R1 and R2. We subsample rather than download the full
// ~94 GB run: a few million pairs is enough to verify barcode integrity.
process FETCH_READS {
    tag "${run.accession}"
    conda 'conda-forge::curl conda-forge::gzip conda-forge::coreutils conda-forge::gawk'

    input:
    run: SraRun
    n_reads: Integer

    output:
    reads: ReadPair = record(
        accession: run.accession,
        r1: file("${run.accession}.subset_R1.fastq.gz"),
        r2: file("${run.accession}.subset_R2.fastq.gz")
    )

    script:
    def acc = run.accession
    def n_lines = n_reads * 4
    """
    set -e
    api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp&format=tsv"
    ftp=\$(curl -sL "\$api" | awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="fastq_ftp") c=i} NR==2{print \$c}')
    if [ -z "\$ftp" ]; then echo "ERROR: no fastq_ftp returned for ${acc}" >&2; exit 1; fi
    r1url=\$(echo "\$ftp" | cut -d';' -f1)
    r2url=\$(echo "\$ftp" | cut -d';' -f2)
    curl -sL "https://\${r1url}" | gzip -dc | head -n ${n_lines} | gzip -c > ${acc}.subset_R1.fastq.gz
    curl -sL "https://\${r2url}" | gzip -dc | head -n ${n_lines} | gzip -c > ${acc}.subset_R2.fastq.gz
    """

    stub:
    def acc = run.accession
    """
    printf '@r/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}.subset_R1.fastq.gz
    printf '@r/2\\nACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIII\\n' | gzip -c > ${acc}.subset_R2.fastq.gz
    """
}
