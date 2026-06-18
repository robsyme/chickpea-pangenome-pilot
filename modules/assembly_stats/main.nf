nextflow.enable.types = true

include { Assembly ; AssemblyStats } from '../../types.nf'

// Contiguity summary for one assembly: scaffold count, total length, N50, %N. Emits a
// machine-readable stats.tsv (AssemblyStats record) plus a MultiQC custom-content row.
// Summariser is a module-scoped binary (resources/usr/bin, on PATH).
process ASSEMBLY_STATS {
    tag "${assembly.accession}-${assembly.assembler}"
    conda 'bioconda::seqkit=2.8.2 conda-forge::python=3.12'

    input:
    assembly: Assembly

    output:
    res: AssemblyStats = record(
        accession: assembly.accession,
        assembler: assembly.assembler,
        file: file("${assembly.accession}.${assembly.assembler}.stats.tsv")
    )
    mqc: Path = file("${assembly.accession}.${assembly.assembler}.assembly_mqc.tsv")

    script:
    def acc = assembly.accession
    def asm = assembly.assembler
    """
    seqkit stats -a -T ${assembly.fasta} > seqkit.tsv
    assembly_stats.py \\
        --seqkit seqkit.tsv \\
        --accession ${acc} \\
        --assembler ${asm} \\
        --out ${acc}.${asm}.stats.tsv \\
        --mqc ${acc}.${asm}.assembly_mqc.tsv
    """

    stub:
    def acc = assembly.accession
    def asm = assembly.assembler
    """
    printf 'accession\\tassembler\\tn50\\n${acc}\\t${asm}\\t23000\\n' > ${acc}.${asm}.stats.tsv
    printf "# id: 'assembly_stats'\\n# plot_type: 'table'\\nSample\\tN50 (bp)\\n${acc}-${asm}\\t23000\\n" > ${acc}.${asm}.assembly_mqc.tsv
    """
}
