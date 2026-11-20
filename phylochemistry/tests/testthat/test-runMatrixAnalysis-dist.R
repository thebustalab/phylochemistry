## Ported from phylochemistry/test_runMatrixAnalysis_dist.R (website tree).
## A call with exactly ONE columns_w_sample_ID_info used to die in the long-format dist branch
## (that column is renamed to sample_unique_ID during pre-processing).

d <- tibble::tibble(
    lake = c("Alpha", "Beta", "Gamma", "Delta"),
    park = c("P1", "P1", "P2", "P2"),
    a    = c(1, 2, 3, 4),
    b    = c(4, 3, 2, 1),
    c    = c(0, 1, 0, 1)
)
vals <- c("a", "b", "c")

test_that("one ID column works and returns the long pair table", {
    one <- runMatrixAnalysis(
        data = d, analysis = "dist", output_format = "long", scale_variance = TRUE,
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("lake")
    )
    expect_true(is.data.frame(one))
    expect_identical(colnames(one)[1:3], c("sample_unique_ID_sample_1", "sample_unique_ID_sample_2", "distance"))
    expect_equal(nrow(one), 12)
    expect_true(all(sort(unique(one$sample_unique_ID_sample_1)) == sort(d$lake)))
})

test_that("two ID columns carry their columns through and agree on the numbers", {
    one <- runMatrixAnalysis(
        data = d, analysis = "dist", output_format = "long", scale_variance = TRUE,
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("lake")
    )
    two <- runMatrixAnalysis(
        data = d, analysis = "dist", output_format = "long", scale_variance = TRUE,
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("lake", "park")
    )
    expect_identical(colnames(two)[1:3], c("sample_unique_ID_sample_1", "sample_unique_ID_sample_2", "distance"))
    expect_true(all(c("lake_sample_1", "park_sample_1", "lake_sample_2", "park_sample_2") %in% colnames(two)))
    expect_equal(nrow(two), 12)
    key1 <- paste(one$sample_unique_ID_sample_1, one$sample_unique_ID_sample_2)
    key2 <- paste(two$lake_sample_1, two$lake_sample_2)
    expect_equal(one$distance[order(key1)], two$distance[order(key2)])
})

test_that("wide output is still a dist object", {
    wide <- runMatrixAnalysis(
        data = d, analysis = "dist", scale_variance = TRUE,
        columns_w_values_for_single_analyte = vals, columns_w_sample_ID_info = c("lake")
    )
    expect_s3_class(wide, "dist")
})
