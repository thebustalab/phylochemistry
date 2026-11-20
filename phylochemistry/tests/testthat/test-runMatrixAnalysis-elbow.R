## Ported from phylochemistry/test_runMatrixAnalysis_elbow.R (website tree).
## analysis = "kmeans", parameters = "elbow" returns the k / within_cluster_variance table that book
## ch.10 plots; the flat-clustering escape room (waterfalls, station2) is graded at k = 3 on it.

d <- tibble::tibble(
    specimen = paste0("s", 1:12),
    x = c(0, 1, 0, 1,   20, 21, 20, 21,   0,  1,  0,  1),
    y = c(0, 0, 1, 1,    0,  0,  1,  1,  20, 20, 21, 21)
)
vals <- c("x", "y")
elbow_for <- function(parameters) {
    suppressMessages(runMatrixAnalysis(
        data = d, analysis = "kmeans", parameters = parameters,
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("specimen")
    ))
}

test_that("shape: k and within_cluster_variance", {
    e <- elbow_for(c("elbow", 6))
    expect_true(is.data.frame(e))
    expect_identical(colnames(e), c("k", "within_cluster_variance"))
    expect_equal(nrow(e), 6L)
    expect_identical(as.integer(e$k), 1:6)
    expect_true(is.numeric(e$within_cluster_variance))
})

test_that("it is an elbow and bends at k = 3", {
    e <- elbow_for(c("elbow", 6))
    drops <- -diff(e$within_cluster_variance)
    expect_true(all(drops >= -1e-8))
    expect_true(drops[2] > 10 * max(drops[3:length(drops)]))
})

test_that("default range, capped at nrow - 1, case-insensitive", {
    expect_equal(nrow(elbow_for("elbow")), 10L)
    expect_equal(nrow(elbow_for(c("elbow", 50))), 11L)
    expect_identical(elbow_for("Elbow")$k, elbow_for("elbow")$k)
})

test_that("the guards stop", {
    expect_error(suppressMessages(runMatrixAnalysis(
        data = d, analysis = "pca", parameters = "elbow",
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("specimen")
    )), "only makes sense with analysis", fixed = TRUE)
    expect_error(elbow_for(c("elbow", "banana")))
    expect_error(elbow_for(c("elbow", 1)))
})

test_that("the ordinary k-means path is untouched", {
    km <- suppressMessages(runMatrixAnalysis(
        data = d, analysis = "kmeans", parameters = c(3),
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("specimen")
    ))
    expect_true(is.data.frame(km))
    expect_equal(nrow(km), 12L)
    expect_true("cluster" %in% colnames(km))
    expect_equal(length(unique(km$cluster)), 3L)
    expect_true(all(table(km$cluster) == 4L))
})

test_that("k-means and its elbow honour scale_variance", {
    elbow_sv <- function(sv) suppressMessages(runMatrixAnalysis(
        data = d, analysis = "kmeans", parameters = c("elbow", 3),
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("specimen"),
        scale_variance = sv
    ))$within_cluster_variance[1]
    raw_ss <- sum(sapply(vals, function(v) sum((d[[v]] - mean(d[[v]]))^2)))
    expect_equal(elbow_sv(TRUE), (nrow(d) - 1) * length(vals))
    expect_equal(elbow_sv(FALSE), raw_ss)
    expect_equal(elbow_for(c("elbow", 3))$within_cluster_variance[1], elbow_sv(TRUE))

    d_big <- dplyr::mutate(d, big = rep(c(1000, 6000, 11000), times = 4) + rep(c(0, 1), each = 6))
    km_raw <- suppressMessages(runMatrixAnalysis(
        data = d_big, analysis = "kmeans", parameters = c(3),
        columns_w_values_for_single_analyte = c("x", "y", "big"),
        columns_w_sample_ID_info = c("specimen"), scale_variance = FALSE
    ))
    by_big <- rep(1:3, times = 4)
    expect_equal(nrow(unique(cbind(km_raw$cluster[match(d_big$specimen, km_raw$sample_unique_ID)], by_big))), 3L)
})
