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
