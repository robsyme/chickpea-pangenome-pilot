nextflow.enable.types = true

include { Assembly ; AssemblyStats } from '../../types.nf'

// Contiguity summary for one assembly: scaffold count, total length, N50, %N.
// Emits a machine-readable stats.tsv (AssemblyStats record) plus a MultiQC
// custom-content row (mirrors BARCODE_CHECK). Summariser is a module-scoped
// binary (resources/usr/bin, on PATH). Requires Nextflow >= 26.04.4 for record
// Path-field staging under Fusion with module binaries (#7225/#7226).
process ASSEMBLY_STATS {
    tag "${assembly.accession}-${assembly.coverage}x"
    conda 'bioconda::seqkit=2.8.2 conda-forge::python=3.12'

    input:
    assembly: Assembly

    output:
    res: AssemblyStats = record(
        accession: assembly.accession,
        assembler: assembly.assembler,
        coverage: assembly.coverage,
        file: file("${assembly.accession}.${assembly.coverage}x.stats.tsv")
    )
    mqc: Path = file("${assembly.accession}.${assembly.coverage}x.assembly_mqc.tsv")

    script:
    def acc = assembly.accession
    def cov = assembly.coverage
    """
    seqkit stats -a -T ${assembly.fasta} > seqkit.tsv
    assembly_stats.py \\
        --seqkit seqkit.tsv \\
        --accession ${acc} \\
        --assembler ${assembly.assembler} \\
        --coverage ${cov} \\
        --out ${acc}.${cov}x.stats.tsv \\
        --mqc ${acc}.${cov}x.assembly_mqc.tsv
    """

    stub:
    def acc = assembly.accession
    def cov = assembly.coverage
    """
    printf 'accession\\tassembler\\tcoverage\\tn50\\n${acc}\\t${assembly.assembler}\\t${cov}\\t23000\\n' \\
        > ${acc}.${cov}x.stats.tsv
    printf "# id: 'assembly_stats'\\n# plot_type: 'table'\\nSample\\tN50 (bp)\\n${acc}-${assembly.assembler}-${cov}x\\t23000\\n" \\
        > ${acc}.${cov}x.assembly_mqc.tsv
    """
}
