## Ported from phylochemistry/test_model_first_stats.R (website tree): buildModel(formula =),
## predictWithModel(), tukeyTest()/anovaTest()/tTest()/gamesHowellTest()/pairwiseTTest() on data and
## on fitted models, and pGroups(). Same checks, now against the package rather than parsed source.

check <- function(cond, msg) expect_true(isTRUE(cond), label = msg)
err_of <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))

test_that("model-first statistics verbs behave as the course expects", {
    ## --- data: three groups with different slopes, one group name with a hyphen ------------------------
    set.seed(1)
    d <- expand.grid(x = c(20, 30, 40, 50), rep = 1:4, grp = c("A", "B-2", "C"), stringsAsFactors = FALSE)
    slope <- c(A = 0.10, `B-2` = 0.10, C = 0.30)
    d$y <- 1 + slope[d$grp] * d$x + rnorm(nrow(d), sd = 0.4)
    d$rep <- NULL
    
    ## ===================================================================================================
    ## Phase 1: buildModel(formula =)
    ## ===================================================================================================
    cat("Phase 1: buildModel(formula =)\n")
    
    ## variable-list path: its parts (response, predictors; no old output_variable part) and its numbers
    m_old <- buildModel(data = d[, c("x", "y")], model_type = "linear_regression",
                        predictors = "x", response = "y")
    check(identical(names(m_old), c("model_type", "response", "predictors", "model", "metrics")),
          "variable-list result has exactly its six parts")
    check(m_old$response == "y" && identical(m_old$predictors, "x") && is.null(m_old$output_variable),
          "response and predictors are stored; no output_variable part (old objects still predict via the fallback)")
    check(isTRUE(all.equal(m_old$metrics$value[m_old$metrics$variable == "x"],
                           unname(coef(lm(y ~ x, d))["x"]))),
          "variable-list slope equals lm()'s, at full precision")
    check(all(m_old$metrics$std_err[m_old$metrics$type == "coefficient"] ==
              unname(summary(lm(y ~ x, d))$coefficients[, 2])),
          "standard errors are not rounded")
    
    ## formula path with an interaction: per-group slopes equal separate per-group fits
    m <- buildModel(data = d, formula = y ~ x * grp)
    check(m$model_type == "linear_regression", "formula alone defaults to linear regression")
    check(m$response == "y" && identical(m$predictors, c("x", "grp")), "response and predictors read from the formula")
    check(is.factor(m$data$grp) && identical(levels(m$data$grp), c("A", "B-2", "C")),
          "text group column stored as a factor")
    check(inherits(m$formula, "formula"), "formula stored on the result")
    cf <- coef(m$model)
    fit_slope <- c(A = cf[["x"]], `B-2` = cf[["x"]] + cf[["x:grpB-2"]], C = cf[["x"]] + cf[["x:grpC"]])
    sep_slope <- sapply(c("A", "B-2", "C"), function(g) coef(lm(y ~ x, d[d$grp == g, ]))[["x"]])
    check(isTRUE(all.equal(fit_slope, sep_slope)), "interaction model slopes = separate per-group fits")
    check(any(grepl("x:grp", m$metrics$variable)), "metrics list the interaction coefficients")
    
    ## a character formula works too
    m_chr <- buildModel(data = d, formula = "y ~ x * grp")
    check(isTRUE(all.equal(coef(m_chr$model), coef(m$model))), "formula given as text gives the same model")
    
    ## columns outside the formula are ignored, including their NAs
    d_na <- d; d_na$junk <- NA
    check(is.na(err_of(buildModel(data = d_na, formula = y ~ x * grp))),
          "NA in a column the formula does not use is ignored")
    
    ## predictWithModel on a formula model
    p <- predictWithModel(data = d, model = m)
    check("predicted_y" %in% names(p) && isTRUE(all.equal(unname(p$predicted_y), unname(fitted(m$model)))),
          "predictWithModel works on a formula model")
    
    ## new data with the group as text predicts the same; an unseen group gets a plain message
    nd <- d; nd$grp <- as.character(nd$grp)
    check(isTRUE(all.equal(predictWithModel(data = nd, model = m)$predicted_y, p$predicted_y)),
          "new data with a text group column predicts fine")
    nd$grp[1] <- "Z"
    check(grepl("group the model was not built on", err_of(predictWithModel(data = nd, model = m))),
          "an unseen group in new data: plain error")
    
    ## non-syntactic names, backticked
    dn <- d; names(dn)[names(dn) == "x"] <- "iso-Leucine"
    mn <- buildModel(data = dn, formula = y ~ `iso-Leucine` * grp)
    check(isTRUE(all.equal(unname(coef(mn$model)), unname(coef(m$model)))), "backticked non-syntactic input works")
    check(isTRUE(all.equal(tukeyTest(mn, slopes_of = "iso-Leucine")$p.adj, tukeyTest(m, slopes_of = "x")$p.adj)),
          "tukeyTest slopes_of a non-syntactic input")
    
    ## plain errors
    check(grepl("not a calculation", err_of(buildModel(data = d, formula = log(y) ~ x))),
          "a calculated response is refused (it would mislabel predictions)")
    check(grepl("dot", err_of(buildModel(data = d, formula = y ~ .))), "a dot formula is refused")
    check(grepl("not both", err_of(buildModel(data = d, formula = y ~ x, predictors = "x"))),
          "formula plus predictors is refused")
    check(grepl("not both", err_of(buildModel(data = d, formula = y ~ x, input_variables = "x"))),
          "formula plus the old input_variables is refused too")
    check(grepl("linear_regression", err_of(buildModel(data = d, model_type = "random_forest_regression", formula = y ~ x))),
          "formula with a forest is refused")
    check(grepl("not present", err_of(buildModel(data = d, formula = y ~ nope))),
          "a missing formula column is named")
    
    ## predictors / response replaced input_variables / output_variable (2026-11-07); the old names are silent aliases
    cat("Aliases: input_variables / output_variable\n")
    m_new_lin <- buildModel(data = d, model_type = "linear_regression", predictors = c("x", "grp"), response = "y")
    m_alias_lin <- buildModel(data = d, model_type = "linear_regression", input_variables = c("x", "grp"), output_variable = "y")
    check(identical(coef(m_new_lin$model), coef(m_alias_lin$model)) && identical(m_new_lin$metrics, m_alias_lin$metrics),
          "old names give an identical linear fit")
    check(m_alias_lin$response == "y" && identical(m_alias_lin$predictors, c("x", "grp")),
          "old names fill $response and $predictors")
    check(identical(buildModel(d, "linear_regression", c("x", "grp"), "y")$metrics, m_new_lin$metrics),
          "positional call (data, model_type, predictors, response) works")
    d_lg <- d; d_lg$high <- ifelse(d_lg$y > median(d_lg$y), "yes", "no")
    m_new_lg <- suppressWarnings(buildModel(data = d_lg, model_type = "logistic_regression", predictors = "x", response = "high"))
    m_alias_lg <- suppressWarnings(buildModel(data = d_lg, model_type = "logistic_regression", input_variables = "x", output_variable = "high"))
    check(identical(coef(m_new_lg$model), coef(m_alias_lg$model)), "old names give an identical logistic fit")
    check(identical(predictWithModel(data = d_lg, model = m_new_lg), predictWithModel(data = d_lg, model = m_alias_lg)),
          "old-name logistic model predicts identically")
    check(grepl("only once", err_of(buildModel(data = d, model_type = "linear_regression", predictors = "x", input_variables = "x", response = "y"))),
          "predictors and input_variables both given: plain error")
    check(grepl("only once", err_of(buildModel(data = d, model_type = "linear_regression", predictors = "x", response = "y", output_variable = "y"))),
          "response and output_variable both given: plain error")
    old_shape <- list(model_type = m_new_lin$model_type, output_variable = m_new_lin$response,
                      model = m_new_lin$model, metrics = m_new_lin$metrics)   # what buildModel returned before 2026-11-07
    check(identical(names(old_shape), c("model_type", "output_variable", "model", "metrics")), "old-shaped object built (only $output_variable)")
    check(identical(predictWithModel(data = d, model = old_shape), predictWithModel(data = d, model = m_new_lin)) &&
          "predicted_y" %in% names(predictWithModel(data = d, model = old_shape)),
          "an old-shaped model object (no $response) still predicts, into predicted_y")
    check(identical(buildModel2(d, "linear_regression", c("x", "grp"), "y")$metrics, m_new_lin$metrics),
          "buildModel2 alias still works")
    
    ## ===================================================================================================
    ## Phase 2: tukeyTest(model) and pGroups()
    ## ===================================================================================================
    cat("Phase 2: tukeyTest(model), pGroups()\n")
    
    ## one-way means: the model path agrees with the data path (rstatix rounds p.adj to 3 sig figs)
    d_h <- d[d$grp != "B-2", ]                       # rstatix itself garbles hyphenated names, so leave B-2 out here
    d_h <- d_h[-c(1, 2, 3), ]                         # unbalanced on purpose
    m1 <- buildModel(data = d_h, formula = y ~ grp)
    tk_model <- tukeyTest(m1)
    tk_data <- tukeyTest(d_h, y ~ grp)
    check(identical(names(tk_model), names(tk_data)), "model path has tukey_hsd's columns")
    check(identical(tk_model$group1, tk_data$group1) && identical(tk_model$group2, tk_data$group2),
          "same pairs in the same order as the data path")
    check(isTRUE(all.equal(tk_model$estimate, tk_data$estimate, tolerance = 1e-3)),
          "same estimates and sign (group2 - group1)")
    check(all(abs(signif(tk_model$p.adj, 3) - tk_data$p.adj) < 5e-4), "same Tukey p.adj")
    check(identical(pGroups(tk_model)$group, pGroups(tk_data)$group), "same letters as the data path")
    
    ## slopes: equals emtrends + pairwise Tukey done by hand
    tk_s <- tukeyTest(m, slopes_of = "x")
    ref <- as.data.frame(summary(pairs(emmeans::emtrends(lm(y ~ x * grp, d), ~ grp, var = "x"), adjust = "tukey")))
    check(isTRUE(all.equal(tk_s$p.adj, ref$p.value)), "slope p.adj = emtrends + Tukey by hand")
    check(identical(tk_s$group1, c("A", "A", "B-2")) && identical(tk_s$group2, c("B-2", "C", "C")),
          "hyphenated group names come back intact (read from weights, not labels)")
    est <- attr(tk_s, "estimates")
    check(isTRUE(all.equal(setNames(est$estimate, est$group), sep_slope)), "per-group slopes attached")
    g <- pGroups(tk_s)
    shares <- mapply(function(a, b) length(intersect(strsplit(g[a, "group"], "")[[1]],
                                                      strsplit(g[b, "group"], "")[[1]])) > 0,
                     tk_s$group1, tk_s$group2)
    check(identical(unname(shares), tk_s$p.adj > 0.05),
          "slope letters agree with the p.adj of every pair (hyphenated name included)")
    check(isTRUE(all.equal(setNames(g$estimate, g$treatment), sep_slope[g$treatment])),
          "pGroups carries each group's slope")
    
    ## pGroups hyphen fix, on a table built WITHOUT rstatix: A ~ B-1, both differ from C
    hy <- data.frame(group1 = c("A", "A", "B-1"), group2 = c("B-1", "C", "C"), p.adj = c(0.9, 0.001, 0.001))
    gh <- pGroups(hy)
    check(gh["A", "group"] == gh["B-1", "group"] && gh["C", "group"] != gh["A", "group"],
          "pGroups gives correct letters with a hyphen in a group name")
    check(identical(names(gh), c("treatment", "group", "spaced_group")), "no estimate column without a model")
    sp <- pGroups(data.frame(group1 = c("A", "A", "B"), group2 = c("B", "C", "C"), p.adj = c(0.001, 0.001, 0.9)))
    check(all(nchar(sp$spaced_group) == 2) && !any(grepl(" ", sp$spaced_group, fixed = TRUE)) &&
          sp["C", "spaced_group"] == "\u00a0b",
          "spaced_group is padded with non-breaking spaces, so HTML tables keep the alignment")
    
    ## a bare lm works too, and plain errors
    check(isTRUE(all.equal(tukeyTest(lm(y ~ x * grp, d), slopes_of = "x")$p.adj, tk_s$p.adj)), "bare lm accepted")
    check(grepl("no factor", err_of(tukeyTest(buildModel(data = d, formula = y ~ x)))), "no factor: plain error")
    check(grepl("slopes_of", err_of(tukeyTest(m, slopes_of = "grp"))), "slopes_of must be numeric")
    check(grepl("for a model", err_of(tukeyTest(d, y ~ grp, slopes_of = "x"))), "slopes_of on a data table: plain error")
    check(grepl("already contains", err_of(tukeyTest(m, y ~ grp))), "formula given to a model: plain error")
    
    ## data path, hyphenated group names (rstatix::tukey_hsd splits its labels on "-" and garbled them)
    hy <- data.frame(g = rep(c("A", "B-2", "C-x-1"), each = 5), y = c(1:5, 3:7, 10:14) + c(0.1, -0.2, 0.3, 0, -0.1))
    tk_hy <- tukeyTest(hy, y ~ g)
    base_hy <- TukeyHSD(aov(y ~ g, hy))$g
    check(identical(paste0(tk_hy$group2, "-", tk_hy$group1), rownames(base_hy)), "data path: hyphenated names come back intact, same orientation as TukeyHSD")
    check(isTRUE(all.equal(signif(tk_hy$p.adj, 3), signif(unname(base_hy[, "p adj"]), 3))), "data path: hyphenated p.adj match TukeyHSD at 3 significant figures (older rstatix rounds, newer does not)")
    check(isTRUE(all.equal(tk_hy$estimate, unname(base_hy[, "diff"]))), "data path: hyphenated estimates match TukeyHSD")
    check(identical(pGroups(tk_hy)$group, pGroups(tukeyTest(buildModel(data = hy, formula = y ~ g)))$group), "data path letters = model path letters with hyphens")
    hy_g <- dplyr::group_by(transform(hy, site = rep(c("s1", "s2"), length.out = 15)), site)
    tk_hy_g <- tukeyTest(hy_g, y ~ g)
    check(all(c(tk_hy_g$group1, tk_hy_g$group2) %in% c("A", "B-2", "C-x-1")), "data path, grouped data: hyphenated names intact")
    plain <- data.frame(g = rep(c("a", "b", "c"), each = 5), y = c(1:5, 3:7, 10:14))
    check(identical(tukeyTest(plain, y ~ g), rstatix::tukey_hsd(plain, y ~ g)), "data path without hyphens: output identical to rstatix")
    hy2 <- transform(hy, h = rep(c("p-1", "q"), length.out = 15))
    tk_int <- tukeyTest(hy2, y ~ g * h)
    check(all(grepl("^(A|B-2|C-x-1|p-1|q)(:(p-1|q))?$", c(tk_int$group1, tk_int$group2))), "data path, two factors with hyphens: interaction labels intact")
    
    ## gamesHowellTest (2026-11-03): data path = rstatix; model path = data path plus group means; hyphen-safe
    gh_data <- gamesHowellTest(hy, y ~ g)
    check(identical(gh_data, rstatix::games_howell_test(hy, y ~ g)), "gamesHowellTest data path = rstatix::games_howell_test")
    check(all(c(gh_data$group1, gh_data$group2) %in% c("A", "B-2", "C-x-1")), "gamesHowellTest: hyphenated names intact")
    gh_model <- gamesHowellTest(buildModel(data = hy, formula = y ~ g))
    check(isTRUE(all.equal(gh_model$p.adj, gh_data$p.adj)) && identical(gh_model$group1, gh_data$group1), "gamesHowellTest model path = data path")
    check(isTRUE(all.equal(sort(pGroups(gh_model)$estimate), sort(as.numeric(tapply(hy$y, hy$g, mean))))), "gamesHowellTest model path carries group means for pGroups")
    check(grepl("one predictor, a factor", err_of(gamesHowellTest(m))), "gamesHowellTest on a slope model: plain error pointing to tukeyTest(slopes_of)")
    ## pairwiseTTest is unpooled (Welch per pair) by default since 2026-11-03
    pw2 <- data.frame(g = rep(c("a", "b", "c"), c(4, 4, 4)), y = c(1, 2, 3, 4, 2, 3, 9, 11, 5, 5.5, 6, 6.2))
    pt_course <- pairwiseTTest(pw2, y ~ g)
    pt_welch <- rstatix::pairwise_t_test(pw2, y ~ g, pool.sd = FALSE, p.adjust.method = "holm")
    check(isTRUE(all.equal(pt_course$p, pt_welch$p)), "pairwiseTTest is unpooled (Welch) by default")
    
    ## ===================================================================================================
    ## Phase 3: anovaTest(model)
    ## ===================================================================================================
    cat("Phase 3: anovaTest(model)\n")
    d2 <- d; d2$h <- rep(c("p", "q"), length.out = nrow(d2)); d2 <- d2[-c(1, 5, 9), ]   # unbalanced two-way
    a_data <- anovaTest(d2, y ~ grp * h)
    a_model <- as.data.frame(anovaTest(buildModel(data = d2, formula = y ~ grp * h)))
    check(isTRUE(all.equal(a_data, a_model, check.attributes = FALSE)), "two-way unbalanced: model path = data path")
    a_anc_data <- as.data.frame(anovaTest(d, y ~ x * grp))
    a_slope <- as.data.frame(anovaTest(m))
    check(isTRUE(all.equal(a_anc_data, a_slope, check.attributes = FALSE)), "numeric x group: model path = data path")
    msg <- tryCatch({ tukeyTest(buildModel(data = d2, formula = y ~ grp * h), compare = "grp"); "" }, message = function(e) conditionMessage(e))
    check(grepl("part of an interaction", msg), "means of a factor inside an interaction: the student is told")
    src_all <- utils::packageDescription("phylochemistry")$Imports
    check(grepl("emmeans", src_all), "emmeans is in phylochemistry's Imports")
    check("x:grp" %in% a_slope$term, "slope model ANOVA has the interaction row")
    ## The course table matches the ANOVA instrument (2026-11-05): term rows, then noise, then total; shares add to 1.
    inst <- data.frame(
        species = rep(c("A", "B", "C"), each = 10),
        fert = rep(rep(c("control", "fertilised"), each = 5), 3),
        y = c(30,34,33,36,32, 36,39,35,41,37, 44,47,42,46,45, 55,58,52,56,57, 38,41,36,40,39, 40,43,39,42,41))
    ai <- as.data.frame(anovaTest(inst, y ~ species * fert))
    check(identical(names(ai), c("term", "sum_of_squares", "free_values", "variance", "F", "p", "share")), "anovaTest columns are the instrument's")
    check(identical(ai$term, c("species", "fert", "species:fert", "noise", "total")), "rows: terms, noise, total")
    check(isTRUE(all.equal(ai$sum_of_squares, c(1164, 258.1, 98.47, 104, 1625))), "sums of squares match the instrument (4 significant figures)")
    check(identical(ai$free_values, c(2L, 1L, 2L, 24L, 29L)), "free values match the instrument (29 = 30 - 1 spoken for)")
    check(abs(sum(ai$share[ai$term != "total"]) - 1) < 0.0025, "shares of the term and noise rows add to 1 (to rounding)")
    check(abs(ai$F[1] - ai$variance[1] / ai$variance[4]) < 0.05, "F is a row's variance over the noise variance")
    check(all(is.na(ai[ai$term == "total", c("variance", "F", "p")])), "total row carries no variance, F or p")
    ## One factor: share is exactly what rstatix called ges
    one <- as.data.frame(anovaTest(inst, y ~ species))
    ges1 <- rstatix::anova_test(data = inst, y ~ species)$ges
    check(isTRUE(all.equal(one$share[1], ges1)), "one factor: share equals ges")
    check(is.data.frame(ai) && !inherits(ai, "tbl_df"), "a plain data.frame, so no column is hidden when printed")
    check(grepl("one formula", err_of(anovaTest(inst))), "data with no formula: plain error")
    check(grepl("on its own", err_of(anovaTest(m, y ~ x))), "formula given to a model: plain error")
    
    ## ===================================================================================================
    ## tTest(): data path unchanged; model + slopes_of = compares slopes unpooled (Welch), Holm-adjusted
    ## ===================================================================================================
    cat("tTest: data path and slopes_of\n")
    check(identical(tTest(d, y ~ grp), rstatix::t_test(data = d, y ~ grp)), "data + formula: identical to rstatix::t_test")
    check(identical(tTest(d, y ~ grp, var.equal = TRUE, p.adjust.method = "bonferroni"),
                    rstatix::t_test(data = d, y ~ grp, var.equal = TRUE, p.adjust.method = "bonferroni")),
          "data path passes var.equal and p.adjust.method through")
    check(identical(tTest(data = d, formula = y ~ grp, ref.group = "A"), rstatix::t_test(data = d, formula = y ~ grp, ref.group = "A")),
          "data path with named data/formula/ref.group")
    ## unequal scatter: C is five times noisier than A and B-2
    set.seed(2)
    dv <- d
    dv$y <- 1 + slope[dv$grp] * dv$x + rnorm(nrow(dv), sd = ifelse(dv$grp == "C", 2, 0.4))
    mv <- buildModel(data = dv, formula = y ~ x * grp)
    tt <- tTest(mv, slopes_of = "x")
    sep <- lapply(split(dv, dv$grp), function(s) summary(lm(y ~ x, data = s))$coefficients["x", ])
    check(identical(tt$group1, c("A", "A", "B-2")) && identical(tt$group2, c("B-2", "C", "C")), "every pair, in level order")
    check(isTRUE(all.equal(attr(tt, "estimates")$estimate, unname(sapply(sep, `[[`, "Estimate")))), "slopes = separate per-group fits")
    check(isTRUE(all.equal(attr(tt, "estimates")$std_err, unname(sapply(sep, `[[`, "Std. Error")))),
          "each slope's SE from its own group alone (no pooled variance)")
    se_pooled_A <- as.data.frame(summary(emmeans::emtrends(mv$model, specs = "grp", var = "x")))$SE[1]
    check(attr(tt, "estimates")$std_err[1] < se_pooled_A / 2, "quiet group's SE is not inflated by the noisy group (unlike tukeyTest)")
    seA <- sep$A[["Std. Error"]]; seC <- sep$C[["Std. Error"]]
    tAC <- (sep$C[["Estimate"]] - sep$A[["Estimate"]]) / sqrt(seA^2 + seC^2)
    dfAC <- (seA^2 + seC^2)^2 / (seA^4 / 14 + seC^4 / 14)
    check(isTRUE(all.equal(tt$statistic[2], tAC)) && isTRUE(all.equal(tt$df[2], dfAC)) &&
          isTRUE(all.equal(tt$p[2], 2 * pt(-abs(tAC), dfAC))), "A vs C: Welch t, Welch-Satterthwaite df, two-sided p")
    check(isTRUE(all.equal(tt$p.adj, p.adjust(tt$p, "holm"))), "p.adj is Holm by default")
    check(isTRUE(all.equal(tTest(mv, slopes_of = "x", p.adjust.method = "bonferroni")$p.adj, p.adjust(tt$p, "bonferroni"))),
          "p.adjust.method = is honoured on the model path")
    check(isTRUE(all.equal(tt$estimate, tt$estimate)) && isTRUE(all.equal(tt$estimate[2], sep$C[["Estimate"]] - sep$A[["Estimate"]])),
          "estimate = group2 slope - group1 slope (tukeyTest's orientation)")
    pg <- pGroups(tt)
    check(identical(sort(pg$treatment), c("A", "B-2", "C")) && !is.null(pg$estimate) && !is.null(pg$std_err),
          "pGroups reads it: letters plus each group's slope and SE")
    check(pg["A", "group"] != pg["C", "group"] && pg["A", "group"] == pg["B-2", "group"], "letters: C differs, A and B-2 share")
    check(isTRUE(all.equal(tTest(mv$model, slopes_of = "x")$p, tt$p)), "a bare lm works too")
    ## guards
    check(grepl("slopes_of =", err_of(tTest(mv))), "model without slopes_of: plain error")
    check(grepl("for a model", err_of(tTest(dv, y ~ grp, slopes_of = "x"))), "slopes_of with a data table: plain error")
    check(grepl("numeric predictor", err_of(tTest(mv, slopes_of = "grp"))), "slopes_of naming a factor: plain error")
    check(grepl("already contains", err_of(tTest(mv, y ~ grp, slopes_of = "x"))), "formula given to a model: plain error")
    dz <- dv; dz$z <- rnorm(nrow(dz))
    check(grepl("tukeyTest", err_of(tTest(buildModel(data = dz, formula = y ~ x * grp + z), slopes_of = "x"))),
          "an extra predictor is refused, pointing at tukeyTest")
    
})
