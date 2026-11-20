## Ported from phylochemistry/test_runMatrixAnalysis_distance_methods.R (website tree).
## Every distance goes through runMatrixAnalysis()'s internal .rma_distance(); each metric is
## checked against its definition (vegan where vegan has it), and the tree branches use the helper.

d <- tibble::tibble(
    sample = c("S1", "S2", "S3", "S4"),
    a = c(1, 2, 0, 5),
    b = c(2, 4, 3, 1),
    c = c(3, 6, 1, 2),
    e = c(0, 0, 4, 0)
)
vals <- c("a", "b", "c", "e")
x <- as.matrix(d[, vals]); rownames(x) <- d$sample
rma <- function(method, scale = FALSE, ...) {
    out <- runMatrixAnalysis(
        data = d, analysis = "dist", distance_method = method, scale_variance = scale,
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = "sample", ...
    )
    if (inherits(out, "dist")) as.matrix(out)[d$sample, d$sample] else out
}
same <- function(a, b) isTRUE(all.equal(unname(as.matrix(a)), unname(as.matrix(b)), tolerance = 1e-10))

test_that("each metric matches its definition", {
    expect_true(same(rma("euclidean", scale = TRUE), as.matrix(dist(scale(x)))))
    expect_true(same(rma("manhattan"), as.matrix(dist(x, "manhattan"))))
    expect_true(same(rma("correlation"), 1 - cor(t(x))))
    expect_lt(abs(rma("correlation")["S1", "S2"]), 1e-12)
    expect_lt(abs(rma("cosine")["S1", "S2"]), 1e-12)
    expect_lt(abs(rma("cosine")["S1", "S3"] - (1 - sum(x[1, ] * x[3, ]) / sqrt(sum(x[1, ]^2) * sum(x[3, ]^2)))), 1e-12)
    expect_true(same(rma("bray_curtis"), as.matrix(vegan::vegdist(x, "bray"))))
    expect_lt(abs(rma("bray_curtis")["S1", "S3"] - 8 / 14), 1e-12)
    expect_true(same(rma("jaccard"), as.matrix(vegan::vegdist(x, "jaccard", binary = TRUE))))
    expect_equal(rma("jaccard")["S1", "S4"], 0)
})

test_that("aitchison is euclidean on clr and announces its zero replacement", {
    xa <- x; xa[xa == 0] <- min(x[x > 0]) / 2
    clr <- log(xa) - rowMeans(log(xa))
    expect_true(same(suppressMessages(rma("aitchison")), as.matrix(dist(clr))))
    ait_abc <- as.matrix(runMatrixAnalysis(
        data = d, analysis = "dist", distance_method = "aitchison", scale_variance = FALSE,
        columns_w_values_for_single_analyte = c("b", "c"), columns_w_sample_ID_info = "sample"
    ))
    expect_lt(abs(ait_abc["S1", "S2"]), 1e-12)
    expect_message(rma("aitchison"), "zeros were replaced")
})

test_that("raw-value metrics refuse scaled data; unknown names are refused", {
    for (m in c("bray_curtis", "jaccard", "aitchison")) {
        expect_error(suppressMessages(rma(m, scale = TRUE)), "scale_variance = FALSE", fixed = TRUE)
    }
    expect_error(rma("maximum"), "is not available")
})

test_that("long output keeps distinct zero-distance pairs and drops only self-pairs", {
    lj <- rma("jaccard", output_format = "long")
    expect_equal(nrow(lj), 12)
    expect_true(any(lj$sample_unique_ID_sample_1 == "S1" & lj$sample_unique_ID_sample_2 == "S4" & lj$distance == 0))
})

test_that("the tree branches use the same helper", {
    nj <- runMatrixAnalysis(
        data = d, analysis = "hclust_phylo", tree_method = "neighbor_joining",
        distance_method = "bray_curtis", scale_variance = FALSE,
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = "sample"
    )
    expect_true(same(ape::cophenetic.phylo(nj)[d$sample, d$sample],
                     ape::cophenetic.phylo(ape::nj(vegan::vegdist(x, "bray")))[d$sample, d$sample]))
    lk <- suppressMessages(runMatrixAnalysis(
        data = d, analysis = "hclust_phylo", distance_method = "bray_curtis", scale_variance = FALSE,
        agglomeration_method = "average",
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = "sample"
    ))
    ref <- ape::as.phylo(hclust(vegan::vegdist(x, "bray"), method = "average"))
    expect_true(same(ape::cophenetic.phylo(lk)[d$sample, d$sample], ape::cophenetic.phylo(ref)[d$sample, d$sample]))
})

test_that("gower takes mixed columns", {
    dm <- tibble::tibble(item = c("A", "B", "C"), kind = c("x", "x", "y"), flag = c(TRUE, FALSE, TRUE), value = c(1, 2, 5))
    gm <- as.matrix(runMatrixAnalysis(
        data = dm, analysis = "dist", distance_method = "gower",
        columns_w_values_for_single_analyte = c("kind", "flag", "value"), columns_w_sample_ID_info = "item"
    ))
    ref_df <- data.frame(kind = factor(dm$kind), flag = factor(dm$flag), value = dm$value, row.names = dm$item)
    expect_true(same(gm[dm$item, dm$item], as.matrix(cluster::daisy(ref_df, metric = "gower"))))
    expect_lt(abs(gm["A", "B"] - 1.25 / 3), 1e-12)
})
