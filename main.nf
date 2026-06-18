nextflow.enable.types = true
nextflow.enable.moduleBinaries = true

include { FETCH_READS }    from './modules/fetch_reads.nf'
include { READ_STRUCTURE } from './modules/read_structure.nf'
include { BARCODE_CHECK }  from './modules/barcode_check/main.nf'
include { MULTIQC }        from './modules/multiqc.nf'
include { SraRun ; Result } from './types.nf'

// Gate-zero pipeline: for each SRA/ENA run accession, subsample reads, summarise
// read structure, and verify stLFR barcode integrity in read2. Designed to take
// many accessions and parallelize per accession; the pilot passes only one.
params {
    // Run accessions to check (one element per accession; parallelized).
    accessions: List<String> = ['SRR32381426']

    // Read pairs to subsample per accession for the check.
    check_reads: Integer = 2000000

    // Minimum fraction of reads that must carry a valid barcode for the gate to
    // pass. A healthy stLFR library recovers ~85-90%; broken/stripped barcodes
    // score near 0-10%.
    pass_fraction: Float = 0.80
}

workflow {
    main:
    runs = channel.fromList(params.accessions)
        .map { acc -> record(accession: acc) }

    reads     = FETCH_READS(runs, params.check_reads)
    structure = READ_STRUCTURE(reads)
    barcode   = BARCODE_CHECK(
        reads, 
        file("${projectDir}/assets/stlfr_barcode_whitelist.txt"), 
        params.check_reads,
        params.pass_fraction as Float
    )

    // Gather all per-accession QC into one MultiQC report: seqkit stats (native
    // module) plus the barcode-integrity custom-content tables.
    qc_files = structure
        .map { r -> r.file }
        .mix(barcode.mqc)
        .collect()
    report = MULTIQC(qc_files, channel.value('stLFR gate-zero'))

    publish:
    read_structure = structure
    barcode_check  = barcode.res
    multiqc_report = report.report
}

output {
    read_structure: Channel<Result> {
        path { r -> r.accession }
    }
    barcode_check: Channel<Result> {
        path { r -> r.accession }
    }
    multiqc_report: Path {
        path '.'
    }
}
