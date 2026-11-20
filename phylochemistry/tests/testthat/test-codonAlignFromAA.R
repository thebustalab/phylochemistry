## Ported from phylochemistry/test_codonAlignFromAA.R (website tree).
## codonAlignFromAA() is the pure-R replacement for orthologr::codon_aln() / PAL2NAL inside
## alignSequences(mode = "codon_align"): it back-translates a gapped amino-acid alignment onto codons.

test_that("exact-translation case with a gap matches PAL2NAL", {
    aa <- c(NUC1 = "MACD", NUC2 = "MA-D")
    nt <- c(NUC1 = "ATGGCTTGTGAT", NUC2 = "ATGGCTGAT")
    res <- codonAlignFromAA(aa, nt)
    expect_identical(res[["NUC1"]], "ATGGCTTGTGAT")
    expect_identical(res[["NUC2"]], "ATGGCT---GAT")
})

test_that("peptide/nucleotide length mismatch fails loudly", {
    aa <- c(NUC1 = "MACD", NUC2 = "MA-D")
    expect_error(codonAlignFromAA(aa, c(NUC1 = "ATGGCTTGTGAT", NUC2 = "ATGGCTGA")),
                 "inconsistency between the peptide and nucleotide")
})

test_that("a peptide with no nucleotide sequence fails loudly", {
    aa <- c(NUC1 = "MACD", NUC2 = "MA-D")
    expect_error(codonAlignFromAA(aa, c(NUC1 = "ATGGCTTGTGAT")), "no nucleotide sequence named")
})
