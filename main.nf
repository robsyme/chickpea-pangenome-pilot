nextflow.enable.types = true

include { FETCH_READS }    from './modules/fetch_reads.nf'
include { READ_STRUCTURE } from './modules/read_structure.nf'
include { BARCODE_CHECK }  from './modules/barcode_check.nf'
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
        .map { acc -> new SraRun(accession: acc) }

    n_reads   = channel.value(params.check_reads)
    whitelist = channel.value(file("${projectDir}/assets/stlfr_barcode_whitelist.txt"))

    pass_frac = channel.value(params.pass_fraction)

    reads     = FETCH_READS(runs, n_reads)
    structure = READ_STRUCTURE(reads)
    barcode   = BARCODE_CHECK(reads, whitelist, n_reads, pass_frac)

    publish:
    read_structure = structure
    barcode_check  = barcode
}

output {
    read_structure: Channel<Result> {
        path { r -> r.accession }
    }
    barcode_check: Channel<Result> {
        path { r -> r.accession }
    }
}
