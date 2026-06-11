nextflow.enable.types = true

include { SraRun ; ReadPair } from '../types.nf'

// Subsample the first `n_reads` read pairs for an accession. Primary path streams
// directly from ENA and takes the first reads (cheap, no full download). If that
// fails (e.g. ENA serving a broken object), fall back to fastq-dl, which downloads
// the full run with an automatic ENA->SRA fallback, then subsample its output.
process FETCH_READS {
    tag "${run.accession}"
    // Prebuilt fastq-dl image (verified to carry curl, gzip, awk, coreutils, python and
    // a working fastq-dl). Using the image directly avoids a Wave conda re-solve, which
    // had pulled an incompatible `rich` and crashed fastq-dl in the fallback.
    container 'community.wave.seqera.io/library/fastq-dl:3.0.1--fa446f61dfc85bc3'

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
    set -uo pipefail
    R1=${acc}.subset_R1.fastq.gz
    R2=${acc}.subset_R2.fastq.gz

    # A gz is valid only if it decompresses to exactly the requested line count.
    # (head closing the pipe SIGPIPEs upstream, so pipe exit status is unreliable;
    # we trust the output line count instead.)
    count_ok() {
        got=\$(gzip -dc "\$1" 2>/dev/null | wc -l) || got=0
        got=\$(echo "\$got" | tr -d ' ')
        [ "\$got" -eq ${n_lines} ]
    }

    # Primary: stream-subsample directly from ENA, no full download. curl --fail
    # makes HTTP errors fail loudly (not silent HTML); --retry rides out hiccups.
    stream_subset() {
        url="\$1"; out="\$2"
        [ -n "\$url" ] || return 1
        curl -fsSL --retry 3 --retry-delay 5 --retry-all-errors "https://\${url}" \\
            | gzip -dc | head -n ${n_lines} | gzip -c > "\$out" || true
        count_ok "\$out"
    }

    # Fallback: fastq-dl downloads the full run (auto ENA->SRA), then we subsample.
    fastqdl_subset() {
        echo "[FETCH] ENA streaming failed; falling back to fastq-dl (ENA->SRA)" >&2
        rm -rf dl && mkdir dl
        fastq-dl --accession ${acc} --outdir dl >&2 || return 1
        r1f=\$(ls dl/*_1.fastq.gz 2>/dev/null | head -1)
        r2f=\$(ls dl/*_2.fastq.gz 2>/dev/null | head -1)
        [ -n "\$r1f" ] && [ -n "\$r2f" ] || { echo "[FETCH] fastq-dl produced no paired files" >&2; return 1; }
        gzip -dc "\$r1f" | head -n ${n_lines} | gzip -c > "\$R1" || true
        gzip -dc "\$r2f" | head -n ${n_lines} | gzip -c > "\$R2" || true
        rm -rf dl
        count_ok "\$R1" && count_ok "\$R2"
    }

    api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${acc}&result=read_run&fields=fastq_ftp&format=tsv"
    ftp=\$(curl -fsSL --retry 3 --retry-delay 5 --retry-all-errors "\$api" | awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="fastq_ftp") c=i} NR==2{print \$c}') || true
    r1url=\$(echo "\$ftp" | cut -d';' -f1)
    r2url=\$(echo "\$ftp" | cut -d';' -f2)

    if stream_subset "\$r1url" "\$R1" && stream_subset "\$r2url" "\$R2"; then
        echo "[FETCH] ENA stream-subsample OK" >&2
    elif fastqdl_subset; then
        echo "[FETCH] fastq-dl fallback OK" >&2
    else
        echo "ERROR: both ENA streaming and fastq-dl fallback failed for ${acc}" >&2
        exit 1
    fi
    """

    stub:
    def acc = run.accession
    """
    printf '@r/1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}.subset_R1.fastq.gz
    printf '@r/2\\nACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIII\\n' | gzip -c > ${acc}.subset_R2.fastq.gz
    """
}
