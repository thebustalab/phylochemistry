## Smoke test: real calls copied from the Integrated Bioanalytical Chemistry chapters, run on the
## datasets that ship with the package, with only the package's Imports.

`%>%` <- magrittr::`%>%`
alaska_lake_data_wide <- tidyr::pivot_wider(dplyr::select(alaska_lake_data, -element_type),
                                            names_from = "element", values_from = "mg_per_L")

test_that("ch7 networks: long dist -> edge list -> buildNetwork()", {
    alaska_dist_long <- runMatrixAnalysis(
        data = alaska_lake_data_wide,
        analysis = c("dist"),
        output_format = "long",
        scale_variance = TRUE,
        columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15],
        columns_w_sample_ID_info = c("lake", "park")
    )
    alaska_edges <- alaska_dist_long %>%
        dplyr::filter(distance > 0) %>%
        dplyr::mutate(similarity = 1 / (1 + distance) * 100) %>%
        dplyr::select(lake_sample_1, lake_sample_2, similarity)
    net <- suppressMessages(buildNetwork(
        edgelist = alaska_edges,
        node_attributes = dplyr::select(alaska_lake_data_wide, lake, park)
    ))
    expect_true(all(c("nodes", "edges") %in% names(net)))
    expect_gt(nrow(net$nodes), 0)
    air_net <- suppressMessages(buildNetwork(edgelist = passenger_flows, node_attributes = us_airports, directed = TRUE))
    expect_gt(nrow(air_net$edges), 0)
})

test_that("ch8 hierarchical clustering + geom_cut()", {
    set.seed(1)
    AK_lakes_clustered <- suppressMessages(runMatrixAnalysis(
        data = alaska_lake_data_wide,
        analysis = c("hclust"),
        scale_variance = TRUE,
        agglomeration_method = "ward.D2",
        columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15],
        columns_w_sample_ID_info = c("lake", "park")
    ))
    p <- ggtree::ggtree(AK_lakes_clustered) + geom_cut(height = 5)
    expect_s3_class(p, "ggplot")
    expect_true(is.numeric(ggplot2::layer_data(p, length(p$layers))$xintercept))
    expect_error(geom_cut(), "needs either")
})

test_that("ch9 PCA", {
    pca <- runMatrixAnalysis(
        data = alaska_lake_data_wide,
        analysis = c("pca"),
        columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15],
        columns_w_sample_ID_info = c("lake", "park")
    )
    expect_true(all(c("Dim.1", "Dim.2") %in% colnames(pca)) || all(c("Dim_1", "Dim_2") %in% colnames(pca)))
    expect_equal(nrow(pca), nrow(alaska_lake_data_wide))
})

test_that("ch11 comparing means: tests and pGroups()", {
    K_data <- dplyr::filter(hawaii_aquifers, analyte == "K")
    K_data_1_6 <- dplyr::filter(K_data, aquifer_code %in% c("aquifer_1", "aquifer_6"))
    expect_s3_class(K_data_1_6 %>% tTest(abundance ~ aquifer_code), "data.frame")
    expect_s3_class(K_data_1_6 %>% wilcoxTest(abundance ~ aquifer_code), "data.frame")
    expect_s3_class(K_data %>% dplyr::group_by(aquifer_code) %>% shapiroTest(abundance), "data.frame")
    expect_s3_class(K_data %>% leveneTest(abundance ~ aquifer_code), "data.frame")
    gh <- K_data %>% gamesHowellTest(abundance ~ aquifer_code) %>% pGroups()
    expect_true(all(c("treatment", "group") %in% names(gh)))
    dn <- K_data %>% dunnTest(abundance ~ aquifer_code) %>% pGroups()
    expect_true(all(c("treatment", "group") %in% names(dn)))
    expect_s3_class(K_data %>% kruskalTest(abundance ~ aquifer_code), "data.frame")
})

test_that("ch12 models: buildModel() and predictWithModel()", {
    md <- metabolomics_data
    md[] <- lapply(md, function(x) if (is.numeric(x)) replace(x, is.na(x), stats::median(x, na.rm = TRUE)) else x)
    basic_regression_model <- buildModel(
        data = md,
        model_type = "linear_regression",
        predictors = "iso-Leucine",
        response = "Valine"
    )
    expect_identical(names(basic_regression_model), c("model_type", "response", "predictors", "model", "metrics"))
    predictions <- predictWithModel(data = md, model = basic_regression_model)
    expect_true("predicted_Valine" %in% names(predictions))
    unknown_patients <- dplyr::mutate(metabolomics_unknown, patient = dplyr::row_number())
    expect_equal(nrow(predictWithModel(data = unknown_patients, model = basic_regression_model)), nrow(unknown_patients))
    expect_s3_class(as.data.frame(anovaTest(buildModel(data = md, formula = Valine ~ patient_status))), "data.frame")
})

test_that("mode(), dropNA(), readCSV()/readMonolist() round trip, writeMonolist()", {
    expect_equal(as.numeric(mode(c(1, 2, 2, 3))), 2)
    expect_identical(dropNA(c(1, NA, 3)), c(1, 3))
    f <- withr::local_tempfile(fileext = ".csv")
    writeMonolist(solvents, f)
    back <- readMonolist(f)
    expect_equal(nrow(back), nrow(solvents))
    expect_identical(colnames(back), colnames(solvents))
})

test_that("Suggests-only functions say what to install", {
    skip_if(requireNamespace("Biostrings", quietly = TRUE), "Biostrings is installed here")
    expect_error(readFasta("x.fasta"), "Biostrings")
})
