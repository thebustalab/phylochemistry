## stats — moved verbatim from phylochemistry.R (package conversion, commit 1)


    ## Set up better shapiroTest warning messages

        #' Shapiro-Wilk normality test by group
        #'
        #' @param data Not yet documented.
        #' @param var Not yet documented.
        #' @export
        shapiroTest <- function(data, var) {
            # Capture the variable using tidy evaluation
            var <- rlang::enquo(var)
            
            # Rename the specified column to "values"
            data <- data %>% dplyr::rename(values = !!var)
            
            # Check if any group has fewer than 3 observations
            group_sizes <- dplyr::summarize(data, size = dplyr::n())
            if (any(group_sizes$size < 3)) {
                stop("One of the groups has fewer than 3 members. A Shapiro\u2013Wilk test cannot be run on such a group. 
                    Please filter your data or choose new groups.")
            }
            
            modified <- FALSE
            
            # For each group with exactly 3 rows, check if there are duplicates in the 'values' column;
            # if yes, add a small amount of random noise to each point
            data <- data %>% dplyr::group_modify(~ {
                if (nrow(.x) == 3 && dplyr::n_distinct(.x$values) < 3) {
                    .x <- .x %>% dplyr::mutate(
                        values = values + stats::rnorm(
                            n(),
                            mean = values / 100000,
                            sd   = values / 200000
                        )
                    )
                    modified <<- TRUE
                }
                .x
            })
            
            if (modified) {
                message(
                    "Note: Found a group with exactly 3 data points where at least two were identical.\n",
                    "A tiny bit of noise was added to these data points, because the Shapiro\u2013Wilk test\n",
                    "p-value can be artificially low in such edge cases due to test mechanics.\n",
                    "Consider carefully whether you want to interpret the results of the test, even with noise added."
                )
            }
            
            # Run the Shapiro–Wilk test on the modified dataset; note that the test will use the 'values' column
            rstatix::shapiro_test(data = data, var = values)
        }


    ## Rename some stats functions so they use camelCase

        #' Levene's test for equal variances
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @export
        leveneTest <- function( data, ... ) {rstatix::levene_test( data = data, ... )}

        # tTest() takes data + formula (rstatix::t_test, unchanged), or a fitted linear model from buildModel()
        # with slopes_of =: then it compares that numeric predictor's slope between every pair of groups with
        # a Welch-style t-test. Each group's slope and standard error come from a line fitted to that group
        # alone, so no group's scatter is pooled with another's (tukeyTest(slopes_of =) pools them). p.adj is
        # Holm-adjusted by default; estimate = group2 slope - group1 slope, so pGroups() reads it like tukeyTest.
        #' t-tests on data or on a fitted model
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @param compare Not yet documented.
        #' @param slopes_of Not yet documented.
        #' @param p.adjust.method Not yet documented.
        #' @export
        tTest <- function( data, ..., compare = NULL, slopes_of = NULL, p.adjust.method = "holm" ) {
            is_model <- inherits(data, "lm") ||
                (is.list(data) && !is.data.frame(data) && !is.null(data$model_type) && !is.null(data$model))
            if (!is_model) {
                if (!is.null(compare) || !is.null(slopes_of)) {
                    stop("tTest(): compare = and slopes_of = are for a model from buildModel(). With a data table, give a formula, e.g. tTest(abundance ~ aquifer_code).", call. = FALSE)
                }
                return(rstatix::t_test( data = data, ..., p.adjust.method = p.adjust.method ))
            }
            if (length(list(...)) > 0) {
                stop("tTest(): a model already contains its formula. Use slopes_of = (which numeric predictor) and, if the model has more than one factor, compare = (which groups).", call. = FALSE)
            }
            if (is.null(slopes_of)) {
                stop("tTest(): with a model, name the numeric predictor whose slopes to compare, e.g. tTest(model, slopes_of = \"temperature\"). To compare group means, give the data and a formula instead.", call. = FALSE)
            }
            if (inherits(data, "lm")) {
                fit <- data
            } else {
                if (data$model_type != "linear_regression") {
                    stop("tTest(): only linear regression models can be t-tested.", call. = FALSE)
                }
                fit <- data$model
            }
            mf <- stats::model.frame(fit)
            inputs <- all.vars(stats::formula(fit)[[3]])
            groups <- inputs[!vapply(mf[inputs], is.numeric, logical(1))]
            if (is.null(compare)) {
                if (length(groups) == 0) {
                    stop("tTest(): this model has no factor (a text predictor) to compare. Build it with one, e.g. formula = rate ~ temperature * chemistry.", call. = FALSE)
                }
                if (length(groups) > 1) {
                    stop("tTest(): this model has more than one factor (", paste(groups, collapse = ", "), "). Choose one with compare =.", call. = FALSE)
                }
                compare <- groups
            }
            if (!compare %in% groups) {
                stop("tTest(): compare = should name a factor (a text predictor) of the model: ", paste(groups, collapse = ", "), ".", call. = FALSE)
            }
            numeric_inputs <- setdiff(inputs, groups)
            if (!slopes_of %in% numeric_inputs) {
                stop("tTest(): slopes_of = should name a numeric predictor of the model: ", paste(numeric_inputs, collapse = ", "), ".", call. = FALSE)
            }
            others <- setdiff(inputs, c(compare, slopes_of))
            if (length(others) > 0) {
                stop("tTest(): slopes_of = fits one line per group, so the model can hold only ", slopes_of, " and ", compare, ", not also ", paste(others, collapse = ", "), ". Use tukeyTest(slopes_of =) for this model.", call. = FALSE)
            }
            y <- stats::model.response(mf)
            x <- mf[[slopes_of]]
            g <- mf[[compare]]
            group_levels <- if (is.factor(g)) levels(droplevels(g)) else sort(unique(as.character(g)))
            if (length(group_levels) < 2) {
                stop("tTest(): ", compare, " has only one group, so there are no slopes to compare.", call. = FALSE)
            }
            per_group <- lapply(group_levels, function(lv) {
                keep <- as.character(g) == lv
                gx <- x[keep]
                gy <- y[keep]
                if (length(unique(gx)) < 2 || length(gx) < 3) {
                    stop("tTest(): group ", lv, " needs at least three points at two or more values of ", slopes_of, " to give a slope and its standard error.", call. = FALSE)
                }
                s <- summary(stats::lm(gy ~ gx))$coefficients
                c(slope = s["gx", "Estimate"], se = s["gx", "Std. Error"], df = length(gx) - 2, n = length(gx))
            })
            per_group <- as.data.frame(do.call(rbind, per_group))
            pairs <- combn(length(group_levels), 2)
            a <- per_group[pairs[1, ], ]
            b <- per_group[pairs[2, ], ]
            se_diff <- sqrt(a$se^2 + b$se^2)
            welch_df <- (a$se^2 + b$se^2)^2 / (a$se^4 / a$df + b$se^4 / b$df)
            statistic <- (b$slope - a$slope) / se_diff
            p <- 2 * stats::pt(-abs(statistic), welch_df)
            p_adj <- stats::p.adjust(p, method = p.adjust.method)
            out <- tibble::tibble(
                term = compare,
                group1 = group_levels[pairs[1, ]],
                group2 = group_levels[pairs[2, ]],
                n1 = a$n,
                n2 = b$n,
                estimate = b$slope - a$slope,
                statistic = statistic,
                df = welch_df,
                p = p,
                p.adj = p_adj,
                p.adj.signif = as.character(stats::symnum(p_adj,
                    cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"), na = ""))
            )
            attr(out, "estimates") <- data.frame(
                group = group_levels,
                estimate = per_group$slope,
                std_err = per_group$se
            )
            out
        }

        #' Wilcoxon rank-sum test
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @export
        wilcoxTest <- function( data, ... ) {rstatix::wilcox_test( data = data, ... )}

        # anovaTest() takes data + formula, or a fitted linear model from buildModel(), and returns the course's
        # ANOVA table, matching the ANOVA instrument column for column: one row per term, then noise, then total.
        # Sequential sums of squares (stats::anova), so the rows add up to the total and share sums to 1 (to rounding).
        # Columns: term, sum_of_squares, free_values (degrees of freedom), variance (sum_of_squares / free_values),
        # F (variance / noise variance), p, share (sum_of_squares / total sum of squares).
        #' Analysis of variance table
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @export
        anovaTest <- function( data, ... ) {
            if (inherits(data, "lm") ||
                (is.list(data) && !is.data.frame(data) && !is.null(data$model_type) && !is.null(data$model))) {
                if (length(list(...)) > 0) {
                    stop("anovaTest(): a model already contains its formula, so give it on its own: model %>% anovaTest().", call. = FALSE)
                }
                if (!inherits(data, "lm") && data$model_type != "linear_regression") {
                    stop("anovaTest(): only linear regression models have an ANOVA table.", call. = FALSE)
                }
                fit <- if (inherits(data, "lm")) data else data$model
            } else {
                dots <- list(...)
                fml <- if (!is.null(dots$formula)) dots$formula else if (length(dots) > 0 && inherits(dots[[1]], "formula")) dots[[1]] else NULL
                if (is.null(fml) || length(dots) > 1) {
                    stop("anovaTest(): give a data table and one formula, e.g. anovaTest(abundance ~ aquifer_code), or a model from buildModel().", call. = FALSE)
                }
                fit <- stats::lm(fml, data = dplyr::ungroup(data))
            }
            tab <- stats::anova(fit)
            ss <- tab[["Sum Sq"]]
            df <- tab[["Df"]]
            terms <- trimws(rownames(tab))
            terms[length(terms)] <- "noise"
            # A plain data.frame, rounded for reading, so a narrow console never hides a column (a tibble
            # truncated the term names and pushed `share` off the edge).
            data.frame(
                term = c(terms, "total"),
                sum_of_squares = signif(c(ss, sum(ss)), 4),
                free_values = c(df, sum(df)),
                variance = signif(c(ss / df, NA), 4),
                F = round(c(tab[["F value"]], NA), 2),
                p = signif(c(tab[["Pr(>F)"]], NA), 3),
                share = round(c(ss / sum(ss), 1), 3)
            )
        }

        # permanovaTest() is the ANOVA table for a whole profile (PERMANOVA, vegan::adonis2): the same rows and
        # columns as anovaTest(), but the spread is measured as distances between samples' profiles, and p comes
        # from shuffling the group labels (so it varies a little from run to run and cannot go below
        # 1 / (permutations + 1)). Two ways in:
        #   wide data, one row per sample:  permanovaTest(cbind(limonene, pinene) ~ species * fertilizer)
        #   long data, one row per value:   permanovaTest(abundance ~ aquifer_code, profile = "analyte", sample = "well_name")
        # distance = any vegan::vegdist method ("euclidean" by default; "bray" for abundance profiles). With one
        # response and "euclidean" it gives anovaTest()'s sums of squares and F exactly. Rows: one per term,
        # added in order (balanced designs: order does not matter), then noise, then total.
        # Everything is inline, not in helpers, because the WebR shim generator carries only the ...Test functions.
        #' Permutational multivariate analysis of variance
        #'
        #' @param data Not yet documented.
        #' @param formula Not yet documented.
        #' @param distance Not yet documented.
        #' @param permutations Not yet documented.
        #' @param profile Not yet documented.
        #' @param sample Not yet documented.
        #' @export
        permanovaTest <- function( data, formula, distance = "euclidean", permutations = 999, profile = NULL, sample = NULL ) {
            if (missing(formula) || !inherits(formula, "formula") || length(formula) != 3) {
                stop("permanovaTest(): give a data table and a formula, e.g. permanovaTest(cbind(limonene, pinene) ~ species * fertilizer).", call. = FALSE)
            }
            if (!requireNamespace("vegan", quietly = TRUE)) {
                stop("permanovaTest(): needs the vegan package: install.packages(\"vegan\").", call. = FALSE)
            }
            data <- as.data.frame(dplyr::ungroup(data))
            group_vars <- all.vars(formula[[3]])
            value_vars <- all.vars(formula[[2]])
            # The left side must be plain column names (one, or several inside cbind()): a transformation such
            # as log(x) would otherwise be dropped without a word, since only the names are read.
            lhs <- formula[[2]]
            lhs_parts <- if (is.call(lhs) && identical(lhs[[1]], as.name("cbind"))) as.list(lhs)[-1] else list(lhs)
            if (!all(vapply(lhs_parts, is.name, logical(1)))) {
                stop("permanovaTest(): the left side must be column names, e.g. cbind(limonene, pinene); transform a column with mutate() first.", call. = FALSE)
            }
            missing_cols <- setdiff(c(group_vars, value_vars, profile, sample), names(data))
            if (length(missing_cols) > 0) {
                stop("permanovaTest(): no column called ", paste(missing_cols, collapse = ", "), " in the data.", call. = FALSE)
            }
            if (!is.null(profile) || !is.null(sample)) {
                if (is.null(profile) || is.null(sample)) {
                    stop("permanovaTest(): long data needs both profile = (the column naming each measurement, e.g. \"analyte\") and sample = (the column or columns naming each sample).", call. = FALSE)
                }
                if (length(value_vars) != 1) {
                    stop("permanovaTest(): with long data the left side is the one value column, e.g. abundance ~ aquifer_code.", call. = FALSE)
                }
                # sample = may name several columns that together identify a sample (two wells can share a name)
                sample_id <- do.call(paste, c(lapply(data[sample], as.character), sep = " / "))
                dup <- anyDuplicated(data.frame(sample_id, data[[profile]]))
                if (dup > 0) {
                    stop("permanovaTest(): sample ", sample_id[dup], " has more than one value for ", profile, " ", data[[profile]][dup],
                         ". If the rows are exact copies, remove them first with distinct(); if they are different samples, add the column that tells them apart to sample =.", call. = FALSE)
                }
                design <- unique(data.frame(sample_id, data[group_vars], stringsAsFactors = FALSE))
                if (anyDuplicated(design$sample_id) > 0) {
                    stop("permanovaTest(): sample ", design$sample_id[anyDuplicated(design$sample_id)], " appears in more than one group; each sample must belong to one.", call. = FALSE)
                }
                wide <- tapply(data[[value_vars]], list(sample_id, as.character(data[[profile]])), function(v) v[1])
                if (anyNA(wide)) {
                    stop("permanovaTest(): some samples are missing some ", profile, " values; every sample needs the full profile.", call. = FALSE)
                }
                design <- design[match(rownames(wide), design$sample_id), group_vars, drop = FALSE]
                profiles <- wide
            } else {
                profiles <- as.matrix(data[value_vars])
                design <- data[group_vars]
                if (anyNA(profiles)) {
                    stop("permanovaTest(): the profile columns contain missing values; filter or fill them first.", call. = FALSE)
                }
            }
            if (anyNA(design)) {
                stop("permanovaTest(): ", paste(group_vars[vapply(design, anyNA, logical(1))], collapse = ", "),
                     " has missing values; filter those samples out first.", call. = FALSE)
            }
            for (v in group_vars) if (is.character(design[[v]])) design[[v]] <- factor(design[[v]])
            dists <- vegan::vegdist(profiles, method = distance)
            fml <- stats::as.formula(call("~", quote(dists), formula[[3]]))
            environment(fml) <- environment()
            tab <- vegan::adonis2(fml, data = design, permutations = permutations, by = "terms")
            ss <- tab[["SumOfSqs"]]
            df <- tab[["Df"]]
            n <- length(ss)
            terms <- trimws(rownames(tab))
            terms[n - 1] <- "noise"
            terms[n] <- "total"
            data.frame(
                term = terms,
                sum_of_squares = signif(ss, 4),
                free_values = df,
                variance = signif(c(ss[-n] / df[-n], NA), 4),
                F = round(tab[["F"]], 2),
                p = signif(tab[["Pr(>F)"]], 3),
                share = round(tab[["R2"]], 3)
            )
        }

        # tukeyTest() takes either data + formula (rstatix, tukey_hsd's first argument is x, not data) or a
        # fitted linear model from buildModel(): then it compares the model's group means, or with
        # slopes_of = the slope of that numeric predictor in each group, Tukey-adjusted via emmeans. The output
        # has tukey_hsd's columns (estimate = group2 - group1) so pGroups() reads both. Everything is inline,
        # not in helpers, because the WebR shim generator carries only the ...Test functions themselves.
        #' Tukey's honest significant difference test
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @param compare Not yet documented.
        #' @param slopes_of Not yet documented.
        #' @export
        tukeyTest <- function( data, ..., compare = NULL, slopes_of = NULL ) {
            is_model <- inherits(data, "lm") ||
                (is.list(data) && !is.data.frame(data) && !is.null(data$model_type) && !is.null(data$model))
            if (!is_model) {
                if (!is.null(compare) || !is.null(slopes_of)) {
                    stop("tukeyTest(): compare = and slopes_of = are for a model from buildModel(). With a data table, give a formula, e.g. tukeyTest(abundance ~ aquifer_code).", call. = FALSE)
                }
                # rstatix::tukey_hsd splits its "group2-group1" labels on "-", so a group name containing a
                # hyphen comes back garbled (B-2 -> "2" / "B"). Swap such levels for zero-padded tokens (same
                # sort order, so every comparison keeps its orientation), test, then put the names back.
                dots <- list(...)
                fml <- if (!is.null(dots$formula)) dots$formula else if (length(dots) > 0 && inherits(dots[[1]], "formula")) dots[[1]] else NULL
                hyphen_key <- character(0)
                if (!is.null(fml) && length(fml) == 3) {
                    rhs_vars <- intersect(all.vars(fml[[3]]), names(data))
                    for (i in seq_along(rhs_vars)) {
                        v <- rhs_vars[i]
                        x <- data[[v]]
                        if (is.numeric(x) || !any(grepl("-", as.character(x)), na.rm = TRUE)) next
                        lev <- if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
                        tok <- sprintf("v%dl%05d", i, seq_along(lev))
                        data[[v]] <- factor(tok[match(as.character(x), lev)], levels = tok)
                        hyphen_key <- c(hyphen_key, stats::setNames(lev, tok))
                    }
                }
                res <- rstatix::tukey_hsd( data, ... )
                if (length(hyphen_key) > 0) {
                    untoken <- function(s) vapply(strsplit(as.character(s), ":", fixed = TRUE), function(p) {
                        hit <- p %in% names(hyphen_key)
                        p[hit] <- hyphen_key[p[hit]]
                        paste(p, collapse = ":")
                    }, character(1))
                    res$group1 <- untoken(res$group1)
                    res$group2 <- untoken(res$group2)
                }
                return(res)
            }
            if (length(list(...)) > 0) {
                stop("tukeyTest(): a model already contains its formula. Use compare = (which groups) and, for slopes, slopes_of = (which numeric predictor).", call. = FALSE)
            }
            if (inherits(data, "lm")) {
                fit <- data
                model_data <- fit$model
            } else {
                if (data$model_type != "linear_regression") {
                    stop("tukeyTest(): only linear regression models can be Tukey-tested.", call. = FALSE)
                }
                fit <- data$model
                model_data <- if (!is.null(data$data)) data$data else fit$model
            }
            inputs <- all.vars(stats::formula(fit)[[3]])
            groups <- inputs[!vapply(model_data[inputs], is.numeric, logical(1))]
            if (is.null(compare)) {
                if (length(groups) == 0) {
                    stop("tukeyTest(): this model has no factor (a text predictor) to compare. Build it with one, e.g. formula = rate ~ temperature * chemistry.", call. = FALSE)
                }
                if (length(groups) > 1) {
                    stop("tukeyTest(): this model has more than one factor (", paste(groups, collapse = ", "), "). Choose one with compare =.", call. = FALSE)
                }
                compare <- groups
            }
            if (!compare %in% groups) {
                stop("tukeyTest(): compare = should name a factor (a text predictor) of the model: ", paste(groups, collapse = ", "), ".", call. = FALSE)
            }
            if (is.null(slopes_of)) {
                in_interaction <- grepl(":", attr(stats::terms(fit), "term.labels")) &
                    vapply(strsplit(attr(stats::terms(fit), "term.labels"), ":"), function(t) compare %in% t, logical(1))
                if (any(in_interaction)) {
                    message("tukeyTest(): ", compare, " is part of an interaction in this model, so these are its means averaged over the other predictors. If the interaction matters, compare within each level instead (filter, then test), or use slopes_of = for a numeric predictor.")
                }
                em <- emmeans::emmeans(fit, specs = compare, data = model_data)
            } else {
                numeric_inputs <- setdiff(inputs, groups)
                if (!slopes_of %in% numeric_inputs) {
                    stop("tukeyTest(): slopes_of = should name a numeric predictor of the model: ", paste(numeric_inputs, collapse = ", "), ".", call. = FALSE)
                }
                em <- emmeans::emtrends(fit, specs = compare, var = slopes_of, data = model_data)
            }
            em_table <- as.data.frame(summary(em))
            prs <- emmeans::contrast(em, method = "pairwise", adjust = "tukey")
            prs_table <- as.data.frame(summary(prs, infer = c(TRUE, TRUE)))
            # which two groups each row compares, read from the contrast weights, never from its label
            weights <- as.matrix(stats::coef(prs)[, -seq_len(ncol(stats::coef(prs)) - nrow(prs_table)), drop = FALSE])
            group_levels <- as.character(em_table[[compare]])
            lo_col <- grep("lower", names(prs_table), value = TRUE)[1]
            hi_col <- grep("upper", names(prs_table), value = TRUE)[1]
            out <- tibble::tibble(
                term = compare,
                group1 = group_levels[apply(weights, 2, function(w) which(w > 0))],
                group2 = group_levels[apply(weights, 2, function(w) which(w < 0))],
                null.value = 0,
                estimate = -prs_table$estimate,
                conf.low = -prs_table[[hi_col]],
                conf.high = -prs_table[[lo_col]],
                p.adj = prs_table$p.value,
                p.adj.signif = as.character(stats::symnum(prs_table$p.value,
                    cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                    symbols = c("****", "***", "**", "*", "ns"), na = ""))
            )
            attr(out, "estimates") <- data.frame(
                group = group_levels,
                estimate = em_table[[if (is.null(slopes_of)) "emmean" else paste0(slopes_of, ".trend")]],
                std_err = em_table$SE
            )
            out
        }

        # gamesHowellTest() compares every pair of group means using each group's OWN spread (Welch-style),
        # so it needs no equal-variance check first: the course's default pairwise test since 2026-11-03
        # (Lucas; with 3-5 per group a Levene test almost never detects unequal spreads, and Tukey then runs
        # at about double its stated false-positive rate). Takes data + formula, or a one-factor model from
        # buildModel(); because Games-Howell uses each group's own spread rather than the model's pooled
        # residual, the model path re-reads the fitted data. Output has rstatix's columns, so pGroups()
        # reads it; the model path also carries each group's mean for pGroups()' estimate column.
        # Inline, not in helpers, because the WebR shim generator carries only the ...Test functions.
        #' Games-Howell post hoc test
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @param compare Not yet documented.
        #' @export
        gamesHowellTest <- function( data, ..., compare = NULL ) {
            is_model <- inherits(data, "lm") ||
                (is.list(data) && !is.data.frame(data) && !is.null(data$model_type) && !is.null(data$model))
            if (!is_model) {
                if (!is.null(compare)) {
                    stop("gamesHowellTest(): compare = is for a model from buildModel(). With a data table, give a formula, e.g. gamesHowellTest(abundance ~ aquifer_code).", call. = FALSE)
                }
                return(rstatix::games_howell_test( data = data, ... ))
            }
            if (length(list(...)) > 0) {
                stop("gamesHowellTest(): a model already contains its formula. Use compare = to say which groups.", call. = FALSE)
            }
            if (inherits(data, "lm")) {
                fit <- data
                model_data <- fit$model
            } else {
                if (data$model_type != "linear_regression") {
                    stop("gamesHowellTest(): only linear regression models can be tested.", call. = FALSE)
                }
                fit <- data$model
                model_data <- if (!is.null(data$data)) data$data else fit$model
            }
            response <- all.vars(stats::formula(fit)[[2]])
            inputs <- all.vars(stats::formula(fit)[[3]])
            groups <- inputs[!vapply(model_data[inputs], is.numeric, logical(1))]
            if (length(inputs) != 1 || length(groups) != 1) {
                stop("gamesHowellTest(): this works on a model with one predictor, a factor, e.g. formula = abundance ~ aquifer_code. With more predictors, compare the groups within one level of the others (filter, then test), or for slopes use tukeyTest(slopes_of = ).", call. = FALSE)
            }
            if (!is.null(compare) && !identical(compare, groups)) {
                stop("gamesHowellTest(): compare = should name the model's factor: ", groups, ".", call. = FALSE)
            }
            d <- as.data.frame(model_data)[, c(response, groups)]
            names(d) <- c("response", "group")
            d$group <- as.character(d$group)
            out <- rstatix::games_howell_test(data = d, formula = response ~ group)
            out$.y. <- response
            means <- stats::aggregate(response ~ group, data = d, FUN = mean)
            ses <- stats::aggregate(response ~ group, data = d, FUN = function(v) stats::sd(v) / sqrt(length(v)))
            attr(out, "estimates") <- data.frame(group = means$group, estimate = means$response, std_err = ses$response)
            out
        }

        # kruskalTest() and dunnTest() refuse a backwards formula. rstatix's kruskal_test and dunn_test
        # do NOT error on `garden ~ value`: they rank the text column and return a confident, meaningless
        # p-value, so a student never learns the formula was reversed (the t and Wilcoxon tests do error,
        # and the console explains that error). The check is inlined in each wrapper, not a shared
        # helper, because the WebR shim generator carries only the ...Test functions themselves.
        #' Kruskal-Wallis test
        #'
        #' @param data Not yet documented.
        #' @param formula Not yet documented.
        #' @param ... Not yet documented.
        #' @export
        kruskalTest <- function( data, formula, ... ) {
            lhs <- if (length(formula) == 3) all.vars(formula[[2]]) else character(0)
            rhs <- if (length(formula) == 3) all.vars(formula[[3]]) else character(0)
            if (length(lhs) == 1 && lhs %in% rhs) {
                stop("kruskalTest(): the formula has `", lhs, "` on both sides of the ~. It reads measurement ~ groups: the numeric column you measured on the left, the column naming the groups on the right.", call. = FALSE)
            }
            if (length(lhs) == 1 && lhs %in% names(data) && !is.numeric(data[[lhs]])) {
                stop("kruskalTest(): `", lhs, "`, on the left of the ~, is not numeric. The formula reads measurement ~ groups: the numeric column you measured on the left, the column naming the groups on the right.", call. = FALSE)
            }
            rstatix::kruskal_test( data = ungroup(data), formula = formula, ... )
        }

        #' Dunn's post hoc test
        #'
        #' @param data Not yet documented.
        #' @param formula Not yet documented.
        #' @param ... Not yet documented.
        #' @export
        dunnTest <- function( data, formula, ... ) {
            lhs <- if (length(formula) == 3) all.vars(formula[[2]]) else character(0)
            rhs <- if (length(formula) == 3) all.vars(formula[[3]]) else character(0)
            if (length(lhs) == 1 && lhs %in% rhs) {
                stop("dunnTest(): the formula has `", lhs, "` on both sides of the ~. It reads measurement ~ groups: the numeric column you measured on the left, the column naming the groups on the right.", call. = FALSE)
            }
            if (length(lhs) == 1 && lhs %in% names(data) && !is.numeric(data[[lhs]])) {
                stop("dunnTest(): `", lhs, "`, on the left of the ~, is not numeric. The formula reads measurement ~ groups: the numeric column you measured on the left, the column naming the groups on the right.", call. = FALSE)
            }
            rstatix::dunn_test( data = ungroup(data), formula = formula, ... )
        }

        # pairwiseWilcoxTest() makes the same cross-group correction as pairwiseTTest() below (see
        # that comment): on grouped data every p-value is adjusted together, not group by group.
        #' Pairwise Wilcoxon tests
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @param p.adjust.method Not yet documented.
        #' @export
        pairwiseWilcoxTest <- function( data, ..., p.adjust.method = "holm" ) {
            if (!dplyr::is_grouped_df(data)) {
                return(rstatix::pairwise_wilcox_test(data = data, ..., p.adjust.method = p.adjust.method))
            }
            res <- rstatix::pairwise_wilcox_test(data = data, ..., p.adjust.method = "none")
            keep <- attributes(res)
            cols <- unclass(res)
            cols[["p.adj"]] <- stats::p.adjust(cols[["p"]], method = p.adjust.method)
            cols[["p.adj.signif"]] <- as.character(stats::symnum(cols[["p.adj"]],
                cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                symbols = c("****", "***", "**", "*", "ns"), na = ""))
            keep$names <- names(cols)
            attributes(cols) <- keep
            if (!is.null(keep$args)) attr(cols, "args")$p.adjust.method <- p.adjust.method
            cols
        }

        # pairwiseTTest() deliberately differs from rstatix::pairwise_t_test() on GROUPED data. rstatix
        # adjusts p-values within each group only, so the chapter's recipe, group_by(facet) then
        # pairwiseTTest(y ~ x) with one comparison per group, got no multiple-comparison correction at
        # all (p.adj == p). Here each group is tested unadjusted, then every p-value across all groups
        # and comparisons is adjusted together with p.adjust.method (default "holm", as in rstatix;
        # "none" is honoured), and p.adj.signif is recomputed with rstatix's cut-points. Columns, column
        # order, class and grouping are kept as rstatix returns them. Ungrouped data go straight to
        # rstatix, unchanged.
        # pool.sd = FALSE (2026-11-03, Lucas): each pair is a Welch t-test on its own two spreads, so, like
        # tTest() and gamesHowellTest(), it needs no equal-variance check first. rstatix's default pools.
        #' Pairwise t-tests
        #'
        #' @param data Not yet documented.
        #' @param ... Not yet documented.
        #' @param p.adjust.method Not yet documented.
        #' @param pool.sd Not yet documented.
        #' @export
        pairwiseTTest <- function( data, ..., p.adjust.method = "holm", pool.sd = FALSE ) {
            if (!dplyr::is_grouped_df(data)) {
                return(rstatix::pairwise_t_test(data = data, ..., pool.sd = pool.sd, p.adjust.method = p.adjust.method))
            }
            res <- rstatix::pairwise_t_test(data = data, ..., pool.sd = pool.sd, p.adjust.method = "none")
            keep <- attributes(res)
            cols <- unclass(res)
            cols[["p.adj"]] <- stats::p.adjust(cols[["p"]], method = p.adjust.method)
            cols[["p.adj.signif"]] <- as.character(stats::symnum(cols[["p.adj"]],
                cutpoints = c(0, 1e-04, 0.001, 0.01, 0.05, 1),
                symbols = c("****", "***", "**", "*", "ns"), na = ""))
            keep$names <- names(cols)
            attributes(cols) <- keep
            if (!is.null(keep$args)) attr(cols, "args")$p.adjust.method <- p.adjust.method
            cols
        }



        #### pGroups


            #' Determine which samples are identical based on p values
            #'
            #' @param data Not yet documented.
            #' @export
            pGroups <- function(data) {

                ## Letters from a pairwise table (tukeyTest, dunnTest, ...): a shared letter = no
                ## significant difference. multcompLetters() splits pair names on "-", so group names
                ## are swapped for hyphen-free tokens and swapped back (a name like "iso-Leucine"
                ## used to get the wrong letters).
                g1 <- as.character(data$group1)
                g2 <- as.character(data$group2)
                lev <- unique(c(g1, g2))
                tok <- paste0("g", seq_along(lev))
                p <- data$p.adj
                names(p) <- paste(tok[match(g1, lev)], tok[match(g2, lev)], sep = "-")

                letters_out <- multcompView::multcompLetters(p)
                spaced_group <- if (is.null(letters_out$monospacedLetters)) letters_out$Letters else letters_out$monospacedLetters
                ## pad with non-breaking spaces: HTML tables (Jupyter, knitted pages) collapse runs of plain
                ## spaces, which made "  cde" display as "cde" and look misaligned (Lucas, 2026-11-02)
                spaced_group <- gsub(" ", "\u00a0", spaced_group, fixed = TRUE)

                treatment <- lev[match(names(letters_out$Letters), tok)]
                output <- data.frame(
                    treatment = treatment,
                    group = unname(letters_out$Letters),
                    spaced_group = unname(spaced_group),
                    row.names = treatment
                )

                ## tukeyTest() on a model attaches each group's mean or slope: carry it along for plotting
                estimates <- attr(data, "estimates")
                if (!is.null(estimates)) {
                    output$estimate <- estimates$estimate[match(output$treatment, estimates$group)]
                    output$std_err <- estimates$std_err[match(output$treatment, estimates$group)]
                }
                return(output)
            }


        #### fitGaussians


            #' Fit a Gaussian curve to a set of x y data
            #'
            #' This function may work better if you feed it a smoothed signal, or add smoothing to the derivatives
            #' @param x Not yet documented.
            #' @param y Not yet documented.
            #' @param peak_detection_threshold Not yet documented.
            #' @param peak_means_stdevs_heights Not yet documented.
            #' @examples
            #' \dontrun{
            #' fitGaussians()
            #' }
            #' @export
            fitGaussians <- function(
                                x, 
                                y, 
                                peak_detection_threshold = 10,
                                peak_means_stdevs_heights = NULL
                            ) {

                ## Find peaks, if peak number is not provided.

                    if (is.null(peak_means_stdevs_heights)) {

                        ## Find "change points" maxima and minima in first derivative = peak half heights
                            deriv1_y <- as.vector(stats::predict(pspline::sm.spline(x, y), x, 1))
                            deriv1_extreme_indices <- c(cumsum(rle(abs(diff(deriv1_y)) > median(diff(deriv1_y))*2)$lengths)+1) # This from https://stackoverflow.com/questions/46029090/how-to-find-changing-points-in-a-dataset
                            deriv1_extremes <- x[deriv1_extreme_indices]

                        ## Find peaks using second derivative - better at detecting shoulder peaks v. 1st derivative
                            deriv2_y <- stats::predict(pspline::sm.spline(x, deriv1_y), x, 1)

                            peak_means <- list()
                            for (i in peak_detection_threshold:(length(x)-peak_detection_threshold)) {
                                if ( # Find all points in 2nd derivative where values on either side of it are increasing for n points in a row
                                    all(
                                        all(deriv2_y[(i-(peak_detection_threshold-1)):(i-1)] > deriv2_y[i]), 
                                        all(deriv2_y[(i+1):(i+peak_detection_threshold)] > deriv2_y[i])
                                    )
                                ) {
                                    peak_means[[i]] <- x[i]
                                }
                            }
                            peak_means <- unlist(peak_means)

                        ## Find sigma for each peak based on mean and "change points"
                            peak_sigma <- vector()
                            for (peak in 1:length(peak_means)) {
                                peak_sigma[[peak]] <- 1.5*min(abs(deriv1_extremes - peak_means[peak]))
                            }
                            
                        ## Peak heights are just the height at the mean
                            peak_heights <- y[x %in% peak_means]

                        ## Report on peaks found and plot them
                            cat(paste0("Found ", length(peak_heights), " peaks at peak_detection_threshold = ", peak_detection_threshold, "\n"))

                            plot <- ggplot() +
                                geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "red") +
                                geom_line(data = data.frame(x = x, y = y), aes(x = x, y = deriv1_y), color = "black") +
                                geom_line(data = data.frame(x = x, y = y), aes(x = x, y = deriv2_y), color = "blue") +
                                geom_vline(aes(xintercept = peak_means)) +
                                geom_vline(aes(xintercept = peak_means + peak_sigma), color = "red") +
                                geom_vline(aes(xintercept = peak_means - peak_sigma), color = "red")
                            print(plot)

                    }

                ## If peak_means_stdevs_heights provided, then use those

                    if (!is.null(peak_means_stdevs_heights)) {
                        n_peaks <- length(peak_means_stdevs_heights)/3
                        peak_means <- peak_means_stdevs_heights[1:n_peaks]
                        peak_sigma <- peak_means_stdevs_heights[(n_peaks+1):(n_peaks*2)]
                        peak_heights <- peak_means_stdevs_heights[((n_peaks*2)+1):(n_peaks*3)]
                    }

                ## Fit peaks with gaussians using nls or manually and extract fit data

                    algorithm = "port"
                    fit <- NULL

                    if (length(peak_means) == 1) {
                        try(
                            fit <-  nls(    y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                sigma1 = peak_sigma[1],
                                                C1 = peak_heights[1]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 2) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2))
                                        ),
                                data = data.frame(x = x, y = y),
                                start = list(
                                            mean1 = peak_means[1],
                                            mean2 = peak_means[2],
                                            sigma1 = peak_sigma[1],
                                            sigma2 = peak_sigma[2],
                                            C1 = peak_heights[1],
                                            C2 = peak_heights[2]
                                        ),
                                algorithm = algorithm, control = list(maxiter = 500000000))
                    , silent = TRUE)
                    }
                    
                    if (length(peak_means) == 3) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 4) {
                        try(
                            fit <-  nls(    y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 5) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2)) +
                                            C5*exp(-(x-mean5)**2/(2 * sigma5**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                mean5 = peak_means[5],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                sigma5 = peak_sigma[5],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4],
                                                C5 = peak_heights[5]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 5000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 6) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2)) +
                                            C5*exp(-(x-mean5)**2/(2 * sigma5**2)) +
                                            C6*exp(-(x-mean6)**2/(2 * sigma6**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                mean5 = peak_means[5],
                                                mean6 = peak_means[6],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                sigma5 = peak_sigma[5],
                                                sigma6 = peak_sigma[6],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4],
                                                C5 = peak_heights[5],
                                                C6 = peak_heights[6]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 10000000))
                        , silent = TRUE)
                    }

                    if (length(peak_means) == 7) {
                        try(
                            fit <- nls( y ~ (
                                            C1*exp(-(x-mean1)**2/(2 * sigma1**2)) +
                                            C2*exp(-(x-mean2)**2/(2 * sigma2**2)) +
                                            C3*exp(-(x-mean3)**2/(2 * sigma3**2)) +
                                            C4*exp(-(x-mean4)**2/(2 * sigma4**2)) +
                                            C5*exp(-(x-mean5)**2/(2 * sigma5**2)) +
                                            C6*exp(-(x-mean6)**2/(2 * sigma6**2)) +
                                            C7*exp(-(x-mean7)**2/(2 * sigma7**2))
                                        ),
                                    data = data.frame(x = x, y = y),
                                    start = list(
                                                mean1 = peak_means[1],
                                                mean2 = peak_means[2],
                                                mean3 = peak_means[3],
                                                mean4 = peak_means[4],
                                                mean5 = peak_means[5],
                                                mean6 = peak_means[6],
                                                mean7 = peak_means[7],
                                                sigma1 = peak_sigma[1],
                                                sigma2 = peak_sigma[2],
                                                sigma3 = peak_sigma[3],
                                                sigma4 = peak_sigma[4],
                                                sigma5 = peak_sigma[5],
                                                sigma6 = peak_sigma[6],
                                                sigma7 = peak_sigma[7],
                                                C1 = peak_heights[1],
                                                C2 = peak_heights[2],
                                                C3 = peak_heights[3],
                                                C4 = peak_heights[4],
                                                C5 = peak_heights[5],
                                                C6 = peak_heights[6],
                                                C7 = peak_heights[7]
                                            ),
                                    algorithm = algorithm, control = list(maxiter = 10000000))
                        , silent = TRUE)
                    }

                    if ( !is.null(fit) ) {
                        
                        cat(paste0("Fit ", length(peak_means), " peak(s) with gaussians at this peak detection threshold.\n\n"))
                        dffit <- data.frame(x = x)
                        dffit$y <- predict(fit, newdata = dffit)
                        fit.sum <- summary(fit)
                        coef.fit <- fit.sum$coefficients[,1]
                        mean.fit <- coef.fit[1:length(peak_means)]
                        sigma.fit <- coef.fit[(length(peak_means)+1):(2*length(peak_means))]
                        C.fit <- coef.fit[((2*length(peak_means))+1):(3*length(peak_means))]
                        
                    } else {
                        
                        cat(paste0("Couldn't fit gaussian(s) at this peak detection threshold\n"))
                        gaussians <- NULL
                        cat(paste0("Trying to fit manually using least squares starting from initial values\n"))

                        # ## Generate gaussian curves, using least squares to optimize

                            mean.fit <- peak_means
                            sigma.fit <- peak_sigma
                            C.fit <- peak_heights

                            init_step <- 0.9
                            step_size <- 0.05
                            final_step <- 1.1

                            grid <- expand.grid(
                                x = seq(sigma.fit[1]*init_step, sigma.fit[1]*final_step, sigma.fit[1]*step_size),
                                y = seq(sigma.fit[2]*init_step, sigma.fit[2]*final_step, sigma.fit[2]*step_size),
                                z = seq(sigma.fit[3]*init_step, sigma.fit[3]*final_step, sigma.fit[3]*step_size),
                                l = seq(mean.fit[1]*init_step, mean.fit[1]*final_step, mean.fit[1]*step_size),
                                m = seq(mean.fit[2]*init_step, mean.fit[2]*final_step, mean.fit[2]*step_size),
                                n = seq(mean.fit[3]*init_step, mean.fit[3]*final_step, mean.fit[3]*step_size),
                                a = seq(C.fit[1]*init_step, C.fit[1]*final_step, C.fit[1]*step_size),
                                b = seq(C.fit[2]*init_step, C.fit[2]*final_step, C.fit[2]*step_size),
                                c = seq(C.fit[3]*init_step, C.fit[3]*final_step, C.fit[3]*step_size)
                            )

                        #     # grid <- expand.grid(
                        #     #     x = seq(sigma.fit[1]*0.999, sigma.fit[1]*1.001, sigma.fit[1]*0.001),
                        #     #     y = seq(sigma.fit[2]*0.999, sigma.fit[2]*1.001, sigma.fit[2]*0.001),
                        #     #     z = seq(sigma.fit[3]*0.999, sigma.fit[3]*1.001, sigma.fit[3]*0.001),
                        #     #     l = seq(mean.fit[1]*0.999, mean.fit[1]*1.001, mean.fit[1]*0.001),
                        #     #     m = seq(mean.fit[2]*0.999, mean.fit[2]*1.001, mean.fit[2]*0.001),
                        #     #     n = seq(mean.fit[3]*0.999, mean.fit[3]*1.001, mean.fit[3]*0.001)
                        #     # )

                        #     # # grid <- expand.grid(
                        #     # #     x = seq(sigma.fit[1]*0.9999, sigma.fit[1]*1.0001, sigma.fit[1]*0.0001),
                        #     # #     y = seq(sigma.fit[2]*0.9999, sigma.fit[2]*1.0001, sigma.fit[2]*0.0001),
                        #     # #     z = seq(sigma.fit[3]*0.9999, sigma.fit[3]*1.0001, sigma.fit[3]*0.0001),
                        #     # #     l = seq(mean.fit[1]*0.9999, mean.fit[1]*1.0001, mean.fit[1]*0.0001),
                        #     # #     m = seq(mean.fit[2]*0.9999, mean.fit[2]*1.0001, mean.fit[2]*0.0001),
                        #     # #     n = seq(mean.fit[3]*0.9999, mean.fit[3]*1.0001, mean.fit[3]*0.0001)
                        #     # # )

                            out <- list()
                            pb <- progress::progress_bar$new(total = dim(grid)[1])
                            for (i in 1:dim(grid)[1]) {
                                pb$tick()
                                gaussians <- list()
                                for (peak in 1:length(peak_means)) {
                                    peak_data <-    data.frame(
                                                        x = x,
                                                        # y = C.fit[peak] * exp(-((x)-mean.fit[peak])**2/(2 * sigma.fit[peak]**2)),
                                                        # y = C.fit[peak] * exp(-((x)-grid[i,peak])**2/(2 * sigma.fit[peak]**2)),
                                                        # y = C.fit[peak] * exp(-((x)-mean.fit[peak])**2/(2 * grid[i,peak]**2)),
                                                        y = grid[i,peak+6] * exp(-((x)-grid[i,peak+3])**2/(2 * grid[i,peak]**2)),
                                                        peak_number = peak
                                                    )
                                    peak_data$peak_area <- sum(peak_data$y)
                                    gaussians[[peak]] <- peak_data
                                }
                                gaussians <- do.call(rbind, gaussians)

                                test_set <- gaussians
                                test_set$ref <- y
                                test_set %>%
                                    group_by(x) %>%
                                    summarize(sum = abs(sum(y) - ref)) %>%
                                    ungroup() %>%
                                    select(sum) %>%
                                    summarize(total = sum(sum)) -> qual

                                out[[i]] <- qual

                            }

                            totals <- as.numeric(as.data.frame(do.call(rbind, out))[,1])

                            out <- data.frame(
                                x = seq(1,length(totals), 1),
                                y = totals
                            )

                            fit_data <- as.numeric(grid[which.min(out$y),])
                            sigma.fit <- fit_data[1:n_peaks]
                            mean.fit <- fit_data[(n_peaks+1):(n_peaks*2)]
                            C.fit <- fit_data[((n_peaks*2)+1):(n_peaks*3)]
                    }

                ## Generate gaussian curves

                    gaussians <- list()
                    for (peak in 1:length(peak_means)) {
                        peak_data <-    data.frame(
                                            x = x,
                                            y = C.fit[peak] * exp(-((x)-mean.fit[peak])**2/(2 * sigma.fit[peak]**2)),
                                            peak_number = peak
                                        )
                        peak_data$peak_area <- sum(peak_data$y)
                        gaussians[[peak]] <- peak_data
                    }
                    gaussians <- do.call(rbind, gaussians)

                ## Plot the results

                    plot <- ggplot() +
                        geom_line(data = data.frame(x = x, y = y), aes(x = x, y = y), color = "red", size = 2) +
                        geom_vline(aes(xintercept = peak_means)) +
                        stat_summary(data = gaussians, geom = "line", fun = "sum", color = "blue", aes(x = x, y = y), size = 2) +
                        geom_line(data = gaussians, aes(x = x, y = y, group = peak_number))

                    print(plot)

                ## Return gaussians and fit parameters

                    output <- list()

                    output$fit_data <- data.frame(
                        peak_height = C.fit,
                        peak_mean = mean.fit,
                        peak_sd = sigma.fit,
                        peak_number = seq(1,length(C.fit),1)
                    )
                    output$gaussians <- gaussians
                    
                    return(output)

            }


        #### mode


            # mode <- function(x, ignore_zero = FALSE) {
                
            #     x <- as.data.frame(table(x))
            #     x <- x[order(x$Freq, decreasing = TRUE),]
            #     if (ignore_zero) {
            #       x <- filter(x, x != 0)
            #     }
                
            #     ## If it can be numeric, return numeric, otherwise, return character
            #         if(any(is.na(as.numeric(as.character(x$x))))) {
            #             return(as.character(x$x[1]))
            #         } else {
            #             return(as.numeric(as.character(x$x[1])))
            #         }
            # }

            #' mode
            #'
            #' @param x Not yet documented.
            #' @param ignore_zero Not yet documented.
            #' @examples
            #' \dontrun{
            #' mode()
            #' }
            #' @export
            mode <- function(x, ignore_zero = FALSE) {
              # Create a frequency table (do not use factor levels for ordering)
              freq_table <- as.data.frame(table(x), stringsAsFactors = FALSE)
              # Rename columns for clarity
              names(freq_table) <- c("value", "count")
              
              # If requested, filter out zeros (note: zeros are treated as character "0")
              if (ignore_zero) {
                freq_table <- subset(freq_table, value != "0")
              }
              
              # Order by frequency count in descending order
              freq_table <- freq_table[order(freq_table$count, decreasing = TRUE), ]
              
              # Get the most common value (the mode)
              mode_val <- freq_table$value[1]
              
              # Try converting to numeric if possible; if conversion fails, return character
              num_val <- suppressWarnings(as.numeric(mode_val))
              if (!is.na(num_val)) {
                return(num_val)
              } else {
                return(mode_val)
              }
            }


        #### geomSignif


            #' geomSignif
            #'
            #' @param data Not yet documented.
            #' @param orientation Not yet documented.
            #' @examples
            #' \dontrun{
            #' geomSignif()
            #' }
            #' @export
            geomSignif <- function(data, orientation) {

                if (orientation == "horizontal") {
                    return(
                        list(
                            geom_segment(
                                data = data,
                                aes(x = xmin, xend = xmax, y = y_position, yend = y_position),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = xmin, xend = xmin, y = y_position, yend = y_position - tip_length_xmin),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = xmax, xend = xmax, y = y_position, yend = y_position - tip_length_xmax),
                                inherit.aes = FALSE
                            ),
                            geom_text(
                                data = data,
                                aes(
                                    x = text_horiz_offset, y = y_position+text_vert_offset,
                                    label = text, size = text_size, hjust = hjust, vjust = vjust
                                ),
                                show.legend = FALSE, inherit.aes = FALSE
                            )
                        )
                    )
                }

                if (orientation == "vertical") {
                    return(
                        list(
                            geom_segment(
                                data = data,
                                aes(x = x_position, xend = x_position, y = ymin, yend = ymax),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = x_position, xend = x_position - tip_length_ymin, y = ymin, yend = ymin),
                                inherit.aes = FALSE
                            ),
                            geom_segment(
                                data = data,
                                aes(x = x_position, xend = x_position - tip_length_ymax, y = ymax, yend = ymax),
                                inherit.aes = FALSE
                            ),
                            geom_text(
                                data = data,
                                aes(
                                    x = x_position+text_horiz_offset, y = text_vert_offset,
                                    label = text, size = text_size, hjust = hjust, vjust = vjust
                                ),
                                show.legend = FALSE, inherit.aes = FALSE
                            )
                        )
                    )
                }
                
            }


        #### grubbsFilter


            #' grubbsFilter on grouped data
            #'
            #' @param data Not yet documented.
            #' @param threshold Not yet documented.
            #' @param col Not yet documented.
            #' @param ... Not yet documented.
            #' @examples
            #' \dontrun{
            #' grubbsFilter()
            #' }
            #' @export
            grubbsFilter <- function(data, col, threshold = 0.05, ...) {
                rlang::check_installed(c("outliers"), reason = "to use grubbsFilter()")

                if (dplyr::is_grouped_df(data)) {
                  return(dplyr::do(data, grubbsFilter(., col, ...))) 
                }

                test_vector <- as.data.frame(data[,which(colnames(data) == col)])[[1]]

                if (length(test_vector) > 3) {
                    result <- outliers::grubbs.test(test_vector)

                    if (result$p.value < threshold) {
                        # print(result)
                        tail <- gsub(" .*$", "", result$alternative)
                        if (tail == "lowest") { data <- data[-which.min(test_vector),] }
                        if (tail == "highest") { data <- data[-which.max(test_vector),] }
                    }
                }

                return(data)

                # data.frame(
                #     id = c(rep("A", 5), rep("B", 5)),
                #     val = c(-10,2,3,4,5,1,2,3,4,15)
                # ) %>%
                # as_tibble() %>%
                # group_by(id) -> data

                # grubbsFilter(data = data, col = "val")

            }


        #### anovaTestMean


            #' anova on means and standard deviations
            #'
            #' @param n Not yet documented.
            #' @param object Not yet documented.
            #' @param ybar Not yet documented.
            #' @param s Not yet documented.
            #' @param ... Not yet documented.
            #' @param ylabel Not yet documented.
            #' @examples
            #' \dontrun{
            #' anovaTestMean()
            #' }
            #' @export
            anovaMean <- function(object, n, ybar, s, ..., ylabel="ylabel") {
                ## the object argument contains the levels of the factor
                df <- sum(n)-length(ybar)
                sigmahat <- (sum((n-1)*s^2) / df)^.5
                SS.group <- sum(n*(ybar - (ybar %*% n)/sum(n))^2)
                Df <- c(length(ybar)-1, sum(n)-length(ybar))
                Sum.of.Sq <- c(SS.group, df*sigmahat^2)
                Mean.Sq <- c(SS.group/(length(ybar)-1), sigmahat^2)
                F.Value <- c(Mean.Sq[1]/Mean.Sq[2], NA)
                Pr.F <- c(1-pf(F.Value[1], Df[1], Df[2]), NA)
                result <- list(
                    Df=Df, "Sum of Sq"=Sum.of.Sq, "Mean Sq"=Mean.Sq,
                    "F value"=F.Value, "Pr(F)"=Pr.F
                )
                oldClass(result) <- c("anova","data.frame")
                attr(result,"row.names") <- format(c(ylabel,"Residuals"))
                attr(result, "heading") <-
                c("Analysis of Variance Table\n",
                paste("Response: ", ylabel, "\n", sep=""),
                "Terms added sequentially (first to last)")
                result
            }


            #' ANOVA from group means, standard deviations and sizes
            #'
            #' @param category Not yet documented.
            #' @param n Not yet documented.
            #' @param mean Not yet documented.
            #' @param stdev Not yet documented.
            #' @export
            anovaTestMean <- function(category, n, mean, stdev) {

                output <- as.data.frame(
                    anovaMean(
                        object = category,
                        n = as.numeric(n),
                        ybar = as.numeric(mean),
                        s = as.numeric(stdev),
                        ylabel = "output"
                    )
                )

                colnames(output)[5] <- "p_value"
                return(output)

            }


        #### allPairsTTestMean

            #' Pairwise t-tests from group means, standard deviations and sizes
            #'
            #' @param category Not yet documented.
            #' @param mean Not yet documented.
            #' @param stdev Not yet documented.
            #' @param n Not yet documented.
            #' @param correction_method Not yet documented.
            #' @export
            allPairsTTestMean <- function(category, mean, stdev, n, correction_method = c("BH", "holm", "bonferroni", "none")) {

                df <- data.frame(
                    category = category,
                    mean = mean,
                    stdev = stdev,
                    n = n
                )

                # Initialize list to store results
                    results <- list()

                # Compute all possible combinations
                    combinations <- combn(df$category, 2)

                # Loop over all combinations
                    for (i in seq(ncol(combinations))) {
                        # Select the two groups
                        group1 <- df[df$category == combinations[1, i],]
                        group2 <- df[df$category == combinations[2, i],]

                        # Compute standard error for each group
                        group1_se <- group1$stdev / sqrt(group1$n)
                        group2_se <- group2$stdev / sqrt(group2$n)

                        # Compute the difference and standard error for the difference between means
                        diff <- group1$mean - group2$mean
                        diff_se <- sqrt(group1_se^2 + group2_se^2)

                        # Compute p-value
                        p_value <- round((2 * pt(-abs(diff / diff_se), (group1$n + group2$n - 2))), 5)

                        # Store results
                        results[[i]] <- data.frame(
                            group1 = group1$category,
                            group2 = group2$category,
                            difference_in_means = diff,
                            difference_in_means_se = round(diff_se, 4),
                            p_value = p_value
                        )
                    }

                # Check the results
                    results <- do.call(rbind, results)
                    if (correction_method[1] == "none") {
                        results$p.adj <- results$p_value
                    } else {
                        results$p.adj <- p.adjust(results$p_value, method = correction_method[1])
                    }
                    results$signif <- ifelse(results$p.adj < 0.05, "*", "")
                    return(results)
            }

