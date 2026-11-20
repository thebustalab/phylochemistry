## runMatrixAnalysis() gets bootstrap support from Bclust(). shipunov (which provides it) is archived
## on CRAN, so the package carries its own computation for when shipunov is missing. The two must
## agree exactly, random draws included.

fun <- function(.x) stats::hclust(stats::dist(.x), method = "ward.D2")
m <- as.matrix(scale(datasets::USArrests[1:12, ]))

test_that("the fallback reproduces shipunov::Bclust exactly", {
    skip_if_not_installed("shipunov")
    set.seed(42); ref <- shipunov::Bclust(m, FUN = fun, iter = 200, monitor = FALSE)
    set.seed(42); own <- phylochemistry:::.bclust_bootstrap(m, FUN = fun, iter = 200, monitor = FALSE,
                                                            bootstrap = TRUE, relative = FALSE, hclist = NULL)
    expect_identical(own$values, ref$values)
    expect_identical(own$hclust$merge, ref$hclust$merge)
    expect_equal(own$consensus, ref$consensus)
    expect_identical(class(own), class(ref))
})

test_that("the fallback gives sensible support values on its own", {
    set.seed(1)
    own <- phylochemistry:::.bclust_bootstrap(m, FUN = fun, iter = 50, monitor = FALSE,
                                              bootstrap = TRUE, relative = FALSE, hclist = NULL)
    expect_length(own$values, nrow(m) - 1)
    expect_true(all(own$values >= 0 & own$values <= 1))
    expect_equal(own$values[nrow(m) - 1], 1)  # the root cluster holds every sample in every replicate
})

test_that("hclust analyses run (with or without shipunov)", {
    wide <- tidyr::pivot_wider(dplyr::select(alaska_lake_data, -element_type),
                               names_from = "element", values_from = "mg_per_L")
    set.seed(7)
    out <- suppressMessages(runMatrixAnalysis(
        data = wide, analysis = "hclust", scale_variance = TRUE, agglomeration_method = "ward.D2",
        columns_w_values_for_single_analyte = colnames(wide)[3:15], columns_w_sample_ID_info = c("lake", "park")
    ))
    expect_true(is.data.frame(out))
    expect_true(all(c("label", "isTip") %in% colnames(out)))
})
