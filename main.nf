nextflow.enable.types = true
nextflow.enable.moduleBinaries = true

include { FETCH_READS }      from './modules/fetch_reads.nf'
include { FETCH_FULL_READS } from './modules/fetch_full_reads.nf'
include { READ_STRUCTURE }   from './modules/read_structure.nf'
include { BARCODE_CHECK }    from './modules/barcode_check/main.nf'
include { LRTK_CONVERT }     from './modules/lrtk_convert.nf'
include { MULTIQC }          from './modules/multiqc.nf'
include { SraRun ; Result } from './types.nf'

// Gate-zero pipeline plus the open-stack assembly baseline (added in later tasks).
params {
    // Run accessions (one element per accession; parallelized). Two on purpose so
    // aggregation steps stay exercised. CBA Captain + Neelam.
    accessions: List<String> = ['SRR32381426', 'SRR32381425']

    // Read pairs to subsample per accession for the gate-zero check.
    check_reads: Integer = 2000000

    // Minimum fraction of reads carrying a valid barcode for the gate to pass.
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

    // --- Stage 1 (open-stack baseline): full reads fetched here; assembly added next.
    full = FETCH_FULL_READS(runs)
    converted = LRTK_CONVERT(full)

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
