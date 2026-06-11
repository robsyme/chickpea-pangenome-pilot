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
    set -euo pipefail

    api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp&format=tsv"
    ftp=\$(curl -fsSL --retry 3 --retry-delay 5 --retry-all-errors "\$api" | awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="fastq_ftp") c=i} NR==2{print \$c}')
    if [ -z "\$ftp" ]; then echo "ERROR: no fastq_ftp returned for ${acc}" >&2; exit 1; fi
    r1url=\$(echo "\$ftp" | cut -d';' -f1)
    r2url=\$(echo "\$ftp" | cut -d';' -f2)

    # Stream, decompress, take the first N lines, recompress. head closing the pipe
    # sends SIGPIPE upstream (expected), so the pipe exit status is unreliable; we
    # validate the read count of the OUTPUT instead. curl --fail makes HTTP errors
    # fail loudly (not silent HTML) and --retry rides out transient ENA hiccups. A
    # short or empty output is a hard failure, which the process errorStrategy retries.
    get_subset() {
        url="\$1"; out="\$2"
        curl -fsSL --retry 3 --retry-delay 5 --retry-all-errors "https://\${url}" \\
            | gzip -dc | head -n ${n_lines} | gzip -c > "\$out" || true
        got=\$(gzip -dc "\$out" 2>/dev/null | wc -l) || got=0
        got=\$(echo "\$got" | tr -d ' ')
        if [ "\$got" -ne ${n_lines} ]; then
            echo "ERROR: \$out truncated -- got \$got of ${n_lines} expected lines (download failed)" >&2
            exit 1
        fi
    }

    get_subset "\$r1url" ${acc}.subset_R1.fastq.gz
    get_subset "\$r2url" ${acc}.subset_R2.fastq.gz
    """

    stub:
    def acc = run.accession
    """
    printf '@r/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}.subset_R1.fastq.gz
    printf '@r/2\\nACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIII\\n' | gzip -c > ${acc}.subset_R2.fastq.gz
    """
}
