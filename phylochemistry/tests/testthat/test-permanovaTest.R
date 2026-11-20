## Ported from phylochemistry/test_permanovaTest.R (website tree). permanovaTest() pinned against
## vegan::adonis2, anovaTest(), long-vs-wide input, the ANOVA chapter's two worked examples, and
## its input guards. Break check: PERMANOVA_TEST_BREAK=1 nudges one plant and must fail.

ft <- utils::read.csv(test_path("fixtures", "fertilizer_trial.csv"))
if (Sys.getenv("PERMANOVA_TEST_BREAK") == "1") ft$limonene[1] <- ft$limonene[1] + 5
`%>%` <- magrittr::`%>%`

test_that("two-way PERMANOVA matches the chapter and adonis2", {
    p2 <- ft %>% permanovaTest(cbind(limonene, pinene) ~ species * fertilizer)
    expect_identical(p2$term, c("species", "fertilizer", "species:fertilizer", "noise", "total"))
    expect_equal(p2$share, c(0.747, 0.089, 0.040, 0.124, 1))
    expect_equal(p2$F[1:3], c(72.15, 17.12, 3.85))
    ref <- vegan::adonis2(ft[c("limonene", "pinene")] ~ species * fertilizer, data = ft,
                          method = "euclidean", by = "terms", permutations = 0)
    expect_equal(p2$sum_of_squares, signif(ref$SumOfSqs, 4))
})

test_that("one response with euclidean distance equals anovaTest()", {
    a <- ft %>% anovaTest(height_cm ~ species * fertilizer)
    p1 <- ft %>% permanovaTest(height_cm ~ species * fertilizer)
    cols <- c("term", "sum_of_squares", "free_values", "variance", "F", "share")
    expect_equal(a[cols], p1[cols])
})

test_that("long (shuffled, grouped) input equals wide input", {
    p2 <- ft %>% permanovaTest(cbind(limonene, pinene) ~ species * fertilizer)
    long <- rbind(data.frame(plant = ft$plant, species = ft$species, fertilizer = ft$fertilizer, compound = "limonene", amount = ft$limonene),
                  data.frame(plant = ft$plant, species = ft$species, fertilizer = ft$fertilizer, compound = "pinene", amount = ft$pinene))
    set.seed(3); long <- long[sample(nrow(long)), ]
    pl <- long %>% dplyr::group_by(species) %>%
        permanovaTest(amount ~ species * fertilizer, profile = "compound", sample = "plant")
    expect_equal(pl[c("sum_of_squares", "F", "share")], p2[c("sum_of_squares", "F", "share")])
    expect_error(long[-1, ] %>% permanovaTest(amount ~ species, profile = "compound", sample = "plant"), "full profile")
})

test_that("hawaii example: bray 47%, euclidean 34%, 106 wells", {
    h <- as.data.frame(hawaii_aquifers)
    h <- h[!h$analyte %in% c("dissolved_solids", "HCO3"), ]
    hb <- h %>% permanovaTest(abundance ~ aquifer_code, profile = "analyte", sample = c("well_name", "longitude", "latitude"), distance = "bray")
    he <- h %>% permanovaTest(abundance ~ aquifer_code, profile = "analyte", sample = c("well_name", "longitude", "latitude"))
    expect_equal(hb$share[1], 0.473)
    expect_equal(he$share[1], 0.343)
    expect_equal(hb$free_values[3], 105)
    expect_error(h %>% permanovaTest(abundance ~ aquifer_code, profile = "analyte", sample = "well_name"), "more than one value")
})

test_that("input guards", {
    expect_error(ft %>% permanovaTest(log(limonene) ~ species), "column names")
    f2 <- ft; f2$species[1] <- NA
    expect_error(f2 %>% permanovaTest(cbind(limonene, pinene) ~ species), "missing values")
})
