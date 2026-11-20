## analyzeGCMSdata7() ingests the folder and, with jupyter = FALSE, returns the
## shinyApp unserved (the same route Utilities/gcms_headless uses). The input is a
## synthetic single-detector (ANDI-Chrom) CDF written here, so no data file ships.

write_synthetic_chrom_cdf <- function(path, n = 3000, interval = 0.1) {
    rt <- (seq_len(n) - 1) * interval
    signal <- 50 + 4000 * exp(-((rt - 80) / 1.5)^2) + 2500 * exp(-((rt - 160) / 2)^2)
    pt <- ncdf4::ncdim_def("point_number", "", seq_len(n), create_dimvar = FALSE)
    v_ord <- ncdf4::ncvar_def("ordinate_values", "", pt, prec = "float")
    v_delay <- ncdf4::ncvar_def("actual_delay_time", "seconds", list(), prec = "float")
    v_int <- ncdf4::ncvar_def("actual_sampling_interval", "seconds", list(), prec = "float")
    nc <- ncdf4::nc_create(path, list(v_ord, v_delay, v_int))
    ncdf4::ncvar_put(nc, v_ord, signal)
    ncdf4::ncvar_put(nc, v_delay, 0)
    ncdf4::ncvar_put(nc, v_int, interval)
    ncdf4::ncatt_put(nc, 0, "sample_name", "synthetic")
    ncdf4::nc_close(nc)
    invisible(path)
}

test_that("analyzeGCMSdata7 ingests a folder and returns a shiny app without launching it", {
    skip_if_not_installed("ncdf4")
    dir <- tempfile("gcms_")
    dir.create(dir)
    old_wd <- getwd()
    on.exit({ setwd(old_wd); unlink(dir, recursive = TRUE) }, add = TRUE)
    write_synthetic_chrom_cdf(file.path(dir, "synthetic.CDF"))

    expect_equal(bustalab:::gcms_detect_cdf_type(file.path(dir, "synthetic.CDF")), "chrom")

    app <- NULL
    capture.output(suppressMessages(
        app <- analyzeGCMSdata7(dir, n_cores = 1, ms_classifier_path = NULL, jupyter = FALSE)
    ))
    expect_s3_class(app, "shiny.appobj")
    ## documented side effects: the app works inside the data folder
    expect_equal(normalizePath(getwd()), normalizePath(dir))
    expect_true(file.exists(file.path(dir, "synthetic.CDF.csv")))
    expect_true(file.exists(file.path(dir, "chromatograms.csv")))
    expect_true(file.exists(file.path(dir, "samples_monolist.csv")))
})

test_that("analyzeGCMSdata7 refuses a folder with no CDF files", {
    dir <- tempfile("gcms_empty_")
    dir.create(dir)
    old_wd <- getwd()
    on.exit({ setwd(old_wd); unlink(dir, recursive = TRUE) }, add = TRUE)
    expect_error(analyzeGCMSdata7(dir, jupyter = FALSE), "does not contain any .CDF files")
})

test_that("the port pool mirror stays internal and unchanged", {
    pool <- bustalab:::.shiny_port_pool
    expect_equal(pool$port, c(10101:10110, 10123, 10456, 10789))
    expect_false("pick_shiny_port" %in% getNamespaceExports("bustalab"))
})
