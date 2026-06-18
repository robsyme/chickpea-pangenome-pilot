nextflow.enable.types = true

include { ReadPair ; ConvertedReads } from '../types.nf'

// Extract stLFR barcodes from read2 into BX:Z: tags, emitting barcoded insert reads
// that ABySS (and later Tigmint/ARKS) consume. Runs once per accession; reused by the
// scaffolding stage. Open Bioconda tool — no private image.
//
// FQCONVER (stLFR) interface verified against lrtk 2.0 (quay biocontainer) on this
// data's barcode layout (offsets 100/116/132): all reads recover a barcode (0 unmapped),
// a `BX:Z:<bc1bc2bc3>-N` tag is written into both read names, and R2 is trimmed to its
// 100 bp insert. Three non-obvious requirements, all handled below:
//   1. `-BW` must be a bwa-indexed FASTA of the barcodes (not the `barcode<TAB>index`
//      whitelist TSV); FQCONVER runs `bwa aln` against it. bwa ships in the lrtk package.
//   2. Inputs must be UNCOMPRESSED fastq (gz makes extract_bc choke on the magic byte).
//   3. Outputs are UNCOMPRESSED (the `.gz` in -O names is ignored), so we gzip after.
process LRTK_CONVERT {
    tag "${reads.accession}"
    conda 'bioconda::lrtk=2.0'
    cpus 8

    input:
    reads: ReadPair
    whitelist: Path

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
    # bwa-indexed FASTA whitelist (>index / barcode) from the repo's barcode TSV.
    awk '{print ">"\$2"\\n"\$1}' ${whitelist} > wl.fa
    bwa index wl.fa

    # lrtk needs uncompressed input.
    gzip -dc ${reads.r1} > r1.fastq
    gzip -dc ${reads.r2} > r2.fastq

    lrtk FQCONVER -IT stLFR -BW wl.fa \\
        -I1 r1.fastq -I2 r2.fastq \\
        -O1 ${acc}.bx_R1.fastq -O2 ${acc}.bx_R2.fastq \\
        -T ${task.cpus}

    # lrtk writes uncompressed; recompress to the published handoff.
    gzip -c ${acc}.bx_R1.fastq > ${acc}.bx_R1.fq.gz
    gzip -c ${acc}.bx_R2.fastq > ${acc}.bx_R2.fq.gz
    """

    stub:
    def acc = reads.accession
    """
    printf '@r BX:Z:AAAA-1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}.bx_R1.fq.gz
    printf '@r BX:Z:AAAA-1\\nACGTACGT\\n+\\nIIIIIIII\\n' | gzip -c > ${acc}.bx_R2.fq.gz
    """
}
