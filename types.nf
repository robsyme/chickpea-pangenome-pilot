// Shared record types.
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

record ConvertedReads {       // LRTK output: BX:Z:-tagged insert reads
    accession: String
    r1: Path
    r2: Path
}

record AssemblyJob {          // one per (accession × k) after combine
    accession: String
    r1: Path
    r2: Path
    k: Integer
}

record Assembly {
    accession: String
    assembler: String         // 'abyss-k64' — encodes the full variant
    fasta: Path
}

record AssemblyStats {
    accession: String
    assembler: String
    file: Path
}
