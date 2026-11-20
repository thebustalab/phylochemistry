test_that("getLabDataset lists the lab datasets without touching the network", {
    expect_setequal(getLabDataset(), c("busta_spectral_library", "plant_phylogeny", "plant_species"))
    expect_error(getLabDataset("nope"), "Unknown lab dataset")
})

test_that("getLabDataset reads a cached copy with the old column types", {
    cache <- tempfile("labcache_")
    dir.create(cache)
    on.exit(unlink(cache, recursive = TRUE), add = TRUE)
    ## a two-row stand-in, placed where the download would land
    writeLines(c("Compound_common_name,mz_50", "007,1.5", "C24 alkane,2"),
               file.path(cache, "busta_spectral_library_v1.csv"))
    lib <- suppressMessages(getLabDataset("busta_spectral_library", cache_dir = cache))
    expect_type(lib$Compound_common_name, "character")
    expect_identical(lib$Compound_common_name[1], "007")
    expect_type(lib$mz_50, "double")
})
