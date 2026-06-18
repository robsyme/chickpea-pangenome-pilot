// Shared record types for the gate-zero pipeline.
// Included by main.nf and the process modules.

record SraRun {
    accession: String
}

record ReadPair {
    accession: String
    r1: Path
    r2: Path
}

record Result {
    accession: String
    file: Path
}

record ConvertedReads {
    accession: String
    dir: Path
}

record AssemblyJob {
    accession: String
    dir: Path
    coverage: Integer
}

record Assembly {
    accession: String
    assembler: String
    coverage: Integer
    fasta: Path
    summary: Path
}

record AssemblyStats {
    accession: String
    assembler: String
    coverage: Integer
    file: Path
}
