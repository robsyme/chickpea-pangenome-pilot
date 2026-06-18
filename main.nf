nextflow.enable.types = true
nextflow.enable.moduleBinaries = true

include { FETCH_READS }    from './modules/fetch_reads.nf'
include { FETCH_FULL_READS } from './modules/fetch_full_reads.nf'
include { STLFR_CONVERT } from './modules/stlfr_convert.nf'
include { READ_STRUCTURE } from './modules/read_structure.nf'
include { BARCODE_CHECK }  from './modules/barcode_check/main.nf'
include { MULTIQC }        from './modules/multiqc.nf'
include { SraRun ; Result } from './types.nf'

// Gate-zero pipeline: for each SRA/ENA run accession, subsample reads, summarise
// read structure, and verify stLFR barcode integrity in read2. Designed to take
// many accessions and parallelize per accession.
params {
    // Run accessions to check (one element per accession; parallelized).
    // The pilot uses TWO accessions on purpose so combination/aggregation steps
    // (MULTIQC now; cross-sample merges later) are always exercised and the pipeline
    // can't silently fail to scale beyond one sample. CBA Captain + Neelam.
    accessions: List<String> = ['SRR32381426', 'SRR32381425']

    // Read pairs to subsample per accession for the check.
    check_reads: Integer = 2000000

    // Minimum fraction of reads that must carry a valid barcode for the gate to
    // pass. A healthy stLFR library recovers ~85-90%; broken/stripped barcodes
    // score near 0-10%.
    pass_fraction: Float = 0.80

    // Private image carrying stlfr2supernova (+ SOAPfilter). Override on the
    // Launchpad once built/pushed; placeholder keeps the repo portable.
    stlfr2supernova_image: String = 'PLACEHOLDER/stlfr2supernova:0.1'
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

    // --- Stage 1: baseline Supernova assembly --------------------------------
    full = FETCH_FULL_READS(runs)
    converted = STLFR_CONVERT(full)

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
