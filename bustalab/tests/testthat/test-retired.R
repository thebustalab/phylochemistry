test_that("retired GC-MS apps stop and point at library(bustalab) + analyzeGCMSdata7()", {
    for (v in 2:6) {
        f <- get(paste0("analyzeGCMSdata", v), envir = asNamespace("bustalab"))
        err <- tryCatch(f("/some/folder"), error = function(e) e)
        expect_s3_class(err, "error")
        msg <- conditionMessage(err)
        expect_match(msg, paste0("analyzeGCMSdata", v, "() has been retired"), fixed = TRUE)
        expect_match(msg, "library(bustalab); analyzeGCMSdata7(", fixed = TRUE)
        expect_false(grepl("source(", msg, fixed = TRUE))
    }
})

test_that("the retired stubs are exported", {
    expect_true(all(paste0("analyzeGCMSdata", 2:7) %in% getNamespaceExports("bustalab")))
})
