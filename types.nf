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
