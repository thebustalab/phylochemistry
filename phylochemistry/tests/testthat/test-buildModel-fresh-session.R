## buildModel() in a FRESH R process with only library(phylochemistry) attached, as a student runs it.
## The old phylochemistry.R attached parsnip, workflows and tune, so a bare call to an unimported
## function (fit() here, once) worked there and in any test session that happened to have them
## loaded, yet failed for students with "could not find function". Running the calls in a clean
## subprocess is what catches that class of bug. Calls are the book's (12c, 12d), shrunk to run fast.

test_that("buildModel() random forests and logistic regression run with only phylochemistry attached", {
    skip_if_not_installed("callr")
    skip_if_not_installed("ranger")
    ## The subprocess loads the INSTALLED package; under devtools::load_all() that copy may be stale.
    installed <- find.package("phylochemistry", lib.loc = .libPaths(), quiet = TRUE)
    loaded <- getNamespaceInfo(asNamespace("phylochemistry"), "path")
    skip_if(length(installed) == 0 || normalizePath(installed[1]) != normalizePath(loaded),
            "needs the package under test to be the installed copy (R CMD check, or test after install)")

    res <- callr::r(function() {
        suppressPackageStartupMessages(library(phylochemistry))
        attached <- search()
        loaded_at_library <- loadedNamespaces()
        run <- function(expr) tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
        set.seed(123)
        out <- list(attached = attached, loaded_at_library = loaded_at_library)
        ## the chapter imputes missing values first; here, keep the complete columns
        md <- metabolomics_data[, colSums(is.na(metabolomics_data)) == 0]
        metabolites <- setdiff(colnames(md)[-(1:2)], "Valine")[1:10]
        out$rf_regression <- run({
            m <- suppressWarnings(buildModel(
                data = md,
                model_type = "random_forest_regression",
                predictors = metabolites,
                response = "Valine",
                fold_cross_validation = 2,
                optimization_parameters = list(n_vars_tried_at_split = c(1, 6), n_trees = 10, min_leaf_size = 1)
            ))
            p <- predictWithModel(data = md, model = m)
            stopifnot(inherits(m$model, "workflow"), "predicted_Valine" %in% names(p))
        })
        out$rf_classification <- run({
            m <- suppressWarnings(buildModel(
                data = md[c(1:34, 41:93), ],
                model_type = "random_forest_classification",
                predictors = metabolites,
                response = "patient_status",
                fold_cross_validation = 2,
                optimization_parameters = list(n_vars_tried_at_split = c(2, 5), n_trees = 10, min_leaf_size = 1)
            ))
            p <- predictWithModel(data = md[35:40, ], model = m)
            stopifnot(inherits(m$model, "workflow"), "predicted_class" %in% names(p))
        })
        out$logistic_regression <- run({
            red <- wine_quality[wine_quality$type == "red", ]
            red$excellent <- ifelse(red$quality_score >= 7, "yes", "no")
            m <- buildModel(
                data = red,
                model_type = "logistic_regression",
                predictors = c("alcohol", "volatile_acidity", "sulphates"),
                response = "excellent"
            )
            p <- predictWithModel(data = red, model = m)
            stopifnot(inherits(m$model, "glm"), "probability_yes" %in% names(p))
        })
        out
    }, libpath = .libPaths())

    ## the point of the subprocess: nothing but phylochemistry (and base R) on the search path
    expect_false(any(c("package:parsnip", "package:workflows", "package:tune", "package:generics") %in% res$attached))
    ## the book-only Imports (R/zzz_book_imports.R) must stay unloaded: loading ggtern, ggside or ggpp
    ## re-registers ggplot2's print/+ methods for the whole session
    expect_false(any(c("ggtern", "ggside", "ggpp") %in% res$loaded_at_library))
    expect_null(res$rf_regression)
    expect_null(res$rf_classification)
    expect_null(res$logistic_regression)
})
