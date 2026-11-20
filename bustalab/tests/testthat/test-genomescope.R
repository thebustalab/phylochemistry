## A synthetic diploid k-mer histogram drawn from GenomeScope's own p = 2 model
## (k = 21, 1% heterozygosity, coverage 30), plus a ramp of low-coverage error k-mers.

write_synthetic_histogram <- function(path) {
    x <- 1:400
    model <- bustalab:::predict2_1(r1 = 0.01, k = 21, d = 0.05, kmercov = 30, bias = 0.5, x = x)
    y <- round(5e7 * model + 2e7 * exp(-x / 1.5))
    writeLines(paste(x, y), path)
}

test_that("genomeScope fits a synthetic diploid histogram and restores options(warn)", {
    skip_if_not_installed("minpack.lm")
    dir <- tempfile("gs_")
    dir.create(dir)
    on.exit(unlink(dir, recursive = TRUE), add = TRUE)
    hist <- file.path(dir, "kmers.histo")
    out <- file.path(dir, "out")
    write_synthetic_histogram(hist)

    old <- options(warn = 1)
    on.exit(options(old), add = TRUE)
    capture.output(genomeScope(hist, out, ploidy = 2, kmer_length = 21))
    expect_equal(getOption("warn"), 1)

    summary_file <- file.path(out, "summary.txt")
    expect_true(file.exists(summary_file))
    s <- readLines(summary_file)
    expect_true(any(grepl("Heterozygous (ab)", s, fixed = TRUE)))
    expect_true(any(grepl("p = 2", s, fixed = TRUE)))
})

test_that("genomeScope restores options(warn) even when it errors", {
    old <- options(warn = 1)
    on.exit(options(old), add = TRUE)
    expect_error(suppressWarnings(capture.output(
        genomeScope(tempfile(), file.path(tempdir(), "gs_missing"), 2, 21))))
    expect_equal(getOption("warn"), 1)
})
