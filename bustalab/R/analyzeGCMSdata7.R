#### Shiny port pool for jupyter = TRUE
####
#### nginx on host1 (bustalab.d.umn.edu) proxies 13 named subdomains to fixed
#### localhost ports. A Shiny app started inside a JupyterHub R session can claim
#### one; nginx then exposes it as https://shinyN.bustalab.d.umn.edu over HTTPS.
####
#### Defined here as a top-level helper so this file remains self-contained
#### (per the bustalabfunctions/AGENTS.md self-contained rule). The same
#### definition exists in phylochemistry.R; if you change the pool, keep both
#### in sync.

    .shiny_port_pool <- data.frame(
        port      = c(10101:10110, 10123, 10456, 10789),
        subdomain = c(paste0("shiny", 1:10), "shiny123", "shiny456", "shiny789"),
        stringsAsFactors = FALSE
    )

    .shiny_port_is_free <- function(port) {
        tryCatch({
            s <- httpuv::startServer("127.0.0.1", port,
                                     list(call = function(req) list(status = 200, body = "")))
            s$stop()
            TRUE
        }, error = function(e) FALSE)
    }

    pick_shiny_port <- function() {
        for (i in sample(nrow(.shiny_port_pool))) {
            p <- .shiny_port_pool$port[i]
            if (.shiny_port_is_free(p)) {
                return(list(
                    port = p,
                    url  = paste0("https://", .shiny_port_pool$subdomain[i], ".bustalab.d.umn.edu")
                ))
            }
        }
        NULL
    }


#### Co-elution split (v7)
####
#### Two compounds under one TIC hump, each with a diagnostic ion (or set of
#### ions) that mostly resolves. Rather than cut the hump with a vertical drop
#### line, the TIC is shared out SCAN BY SCAN in proportion to how much of
#### each compound the diagnostic ions say is present at that scan. The two
#### areas therefore sum exactly to the TIC area of the window, and the
#### overlapping tails go to the compound they belong to.
####
#### The model, per scan t, with a(t) and b(t) the diagnostic signal each
#### compound makes on its OWN ion:
####
####     ion_a(t) = a(t)            + xB * b(t)
####     ion_b(t) = xA * a(t)       + b(t)
####     TIC(t)   = kA * a(t)       + kB * b(t)
####
#### xA, xB are cross-talk ratios taken from each compound's spectrum (A's
#### signal on B's ion relative to A's own ion, and vice versa); 0 by
#### default, i.e. each diagnostic ion belongs to one compound. kA, kB are
#### TIC per unit diagnostic signal (1 / the fraction of the compound's TIC
#### its ion carries). They are FITTED, by least squares of the TIC on a and
#### b over every scan in the window, unless supplied.
####
#### Why fit only two numbers: an earlier version estimated a full 2x2
#### response matrix from each peak's "pure" edge. With cross-talk free to
#### vary, many splits fit the data exactly (closure 1.000), and at 5:1 and
#### 60% overlap it settled on one that put +23% on the small peak. Fixing
#### the spectral ratios and fitting only the scales has a single answer.
#### test_coelution_split.R pins this.
####
#### `fit_r2`: how well kA*a + kB*b reproduces the window's TIC. Well below 1
#### means the window, the baseline or the ions are wrong. It does NOT detect
#### a third compound with no diagnostic ion eluting between the two: that
#### one's profile is absorbed into kA and kB and its area is shared between
#### them (test_coelution_split.R shows fit_r2 0.9998 with one present). kA
#### and kB are properties of the compounds, so they should agree across
#### samples run on the same method; a sample whose fitted k is an outlier is
#### the better tell.
####
#### Three compounds (one diagnostic ion set each) use the same model with a
#### third column, TIC(t) = kA*a + kB*b + kC*c, and no cross-talk (each ion
#### set is taken to belong to its own compound). A compound whose ions have
#### no signal at all in one sample's hump keeps its row there with area 0,
#### so every sample carries the same peaks and the numbers line up. So does
#### one whose per-sample fit comes out non-positive (background only on its
#### ions): it is flagged CHECK and stored with no factor, so refresh keeps it 0.

    ## Signed per-scan sum (the v6 `area` convention: a rectangle sum over
    ## scans, not scaled by the scan interval).
    .gcms_scan_sum <- function(x) sum(x, na.rm = TRUE)

    ## Unmix the diagnostic traces (a list, one per compound) into each
    ## compound's own signal. Cross-talk (see the model above) is defined for a
    ## pair only; with three compounds each ion set is taken to belong to its
    ## own compound. NULL when the ions cannot discriminate the compounds.
    .gcms_unmix <- function(ions, crosstalk = c(0, 0)) {
        if (length(ions) == 2) {
            um <- gcms_unmix_pair(ions[[1]], ions[[2]], crosstalk)
            return(if (is.null(um)) NULL else list(um$a, um$b))
        }
        if (any(crosstalk != 0)) return(NULL)
        lapply(ions, function(x) pmax(x, 0))
    }

    ## Split the TIC across two or three co-eluting compounds. The model above,
    ## with one column per compound: TIC(t) = sum_j k_j * s_j(t).
    ##
    ##   rt            scan retention times
    ##   tic_net       baseline-corrected TIC per scan (signed; v6's TIC - baseline)
    ##   ions          list of baseline-corrected diagnostic traces, one per
    ##                 compound (one ion, or the sum of several), same scans; the
    ##                 order of elution does not matter
    ##   crosstalk     c(xA, xB), two compounds only. Default 0, 0.
    ##   scale         optional k per compound to use instead of fitting them
    ##   allow_absent  TRUE: a compound with no diagnostic signal at all in the
    ##                 window gets area 0 (and its k may be NA) instead of failing
    ##                 the split, so the sample still gets one row per compound
    ##                 and the peak numbers line up across samples. At least one
    ##                 compound must be present.
    ##
    ## Returns list(areas (sum to sum(tic_net)), shares (scans x compounds, in rt
    ## order), contrib (k_j * s_j per scan), k, present, fit_r2, status ("ok" or
    ## a reason)). On failure the areas are NA, never 0.
    gcms_split_coeluting_n <- function(rt, tic_net, ions, crosstalk = c(0, 0),
                                       scale = NULL, allow_absent = FALSE) {
        m <- length(ions)
        fail <- function(why) list(areas = rep(NA_real_, m), shares = NULL, contrib = NULL,
                                   k = rep(NA_real_, m), present = rep(NA, m),
                                   fit_r2 = NA_real_, status = why)
        if (!(m %in% 2:3)) return(fail("give ions for two or three compounds"))
        n <- length(rt)
        if (n < 5 || length(tic_net) != n || any(lengths(ions) != n))
            return(fail("traces missing or of unequal length"))
        if (any(!is.finite(c(rt, tic_net, unlist(ions)))))
            return(fail("non-finite values in the traces"))
        o <- order(rt)
        rt <- rt[o]; tic_net <- tic_net[o]; ions <- lapply(ions, `[`, o)
        if (!allow_absent && any(vapply(ions, max, numeric(1)) <= 0))
            return(fail("a diagnostic ion has no signal in the window"))
        if (length(crosstalk) != 2 || any(!is.finite(crosstalk)) || any(crosstalk < 0))
            return(fail("crosstalk must be two non-negative numbers"))
        if (m == 3 && any(crosstalk != 0))
            return(fail("cross-talk is supported for two compounds only"))

        um <- .gcms_unmix(ions, crosstalk)
        if (is.null(um))
            return(fail("diagnostic ions do not discriminate the compounds"))
        S <- do.call(cbind, um)
        present <- colSums(S) > 0
        if (!any(present)) return(fail("no compound has diagnostic signal in the window"))
        if (!allow_absent && !all(present))
            return(fail("a compound has no diagnostic signal after unmixing"))
        P <- which(present)

        if (is.null(scale)) {
            ## The profiles must differ, or the fit cannot tell their scales
            ## apart. Two real peaks at 90% overlap correlate at ~0.98 and still
            ## fit exactly, so only near-proportional traces (the same compound
            ## on two ion sets) are refused.
            if (length(P) > 1) {
                cr <- suppressWarnings(stats::cor(S[, P, drop = FALSE]))
                if (any(!is.finite(cr)) || max(cr[upper.tri(cr)]) > 0.999)
                    return(fail("diagnostic ions do not discriminate the compounds"))
            }
            kp <- unname(stats::lm.fit(S[, P, drop = FALSE], tic_net)$coefficients)
            if (any(!is.finite(kp)) || any(kp <= 0))
                return(fail(sprintf(
                    "TIC fit gave a non-positive scale (%s); check the ions or supply `scale`",
                    paste(sprintf("k%s %.3g", .gcms_split_labels[P], kp), collapse = ", "))))
            k <- rep(NA_real_, m); k[P] <- kp
        } else {
            k <- as.numeric(scale)
            if (length(k) != m) return(fail("scale must give one factor per compound"))
            ## A compound stored with no factor (it was absent, or only background,
            ## when the split was made) stays absent: area 0.
            if (allow_absent) {
                present <- present & is.finite(k)
                if (!any(present)) return(fail("no compound has both diagnostic signal and a response factor"))
                P <- which(present)
            }
            if (any(!is.finite(k[P])) || any(k[P] <= 0))
                return(fail("scale must be a positive number for every compound"))
        }

        C <- S
        for (j in seq_len(m)) C[, j] <- if (present[j]) k[j] * S[, j] else 0
        tot <- rowSums(C)
        resid <- tic_net - tot
        fit_r2 <- 1 - sum(resid^2) / sum((tic_net - mean(tic_net))^2)

        ## Scans with no diagnostic signal (the window's edges) take the shares
        ## of the nearest scans that have some.
        ok <- which(tot > 0)
        if (length(ok) < 2) return(fail("no scans carry diagnostic signal"))
        shares <- matrix(0, n, m)
        for (j in P)
            shares[, j] <- stats::approx(rt[ok], C[ok, j] / tot[ok], xout = rt,
                                         method = "constant", rule = 2, f = 0.5)$y
        areas <- vapply(seq_len(m), function(j) .gcms_scan_sum(tic_net * shares[, j]), numeric(1))

        list(areas = areas, shares = shares, contrib = C, k = k, present = present,
             fit_r2 = fit_r2, status = "ok")
    }

    ## The two-compound form, kept for its drop-line comparison. Returns a list:
    ## area_a, area_b (sum to sum(tic_net)), share_a per scan, kA/kB used,
    ## fit_r2, the ion-guided drop-line answer for comparison (drop_rt,
    ## area_a_drop, area_b_drop) and `status` ("ok" or a reason). On failure
    ## the areas are NA, never 0.
    gcms_split_coeluting <- function(rt, tic_net, ion_a, ion_b,
                                     crosstalk = c(0, 0),
                                     scale = NULL) {
        s <- gcms_split_coeluting_n(rt, tic_net, list(ion_a, ion_b), crosstalk = crosstalk, scale = scale)
        if (s$status != "ok") return(list(
            area_a = NA_real_, area_b = NA_real_, share_a = NULL,
            k_a = NA_real_, k_b = NA_real_, fit_r2 = NA_real_,
            drop_rt = NA_real_, area_a_drop = NA_real_, area_b_drop = NA_real_,
            status = s$status))
        o <- order(rt); rt <- rt[o]; tic_net <- tic_net[o]; n <- length(rt)
        share_a <- s$shares[, 1]

        ## Drop-line comparison: the first scan after A's apex where B takes
        ## the majority of the TIC.
        after <- which(seq_len(n) > which.max(s$contrib[, 1]) & share_a < 0.5)
        if (length(after)) {
            j <- after[1]
            drop_rt <- rt[j]
            area_a_drop <- .gcms_scan_sum(tic_net[seq_len(j - 1)])
            area_b_drop <- .gcms_scan_sum(tic_net[j:n])
        } else {
            drop_rt <- NA_real_; area_a_drop <- NA_real_; area_b_drop <- NA_real_
        }

        list(area_a = s$areas[1], area_b = s$areas[2], share_a = share_a,
             k_a = s$k[1], k_b = s$k[2], fit_r2 = s$fit_r2,
             drop_rt = drop_rt, area_a_drop = area_a_drop, area_b_drop = area_b_drop,
             status = "ok")
    }


    ## Unmix two diagnostic traces into each compound's own signal, given the
    ## cross-talk ratios (see the model above). NULL when the ions cannot
    ## discriminate the two compounds.
    gcms_unmix_pair <- function(ion_a, ion_b, crosstalk = c(0, 0)) {
        X <- matrix(c(1, crosstalk[1], crosstalk[2], 1), 2, 2)
        if (abs(det(X)) < 0.05) return(NULL)
        ab <- solve(X, rbind(ion_a, ion_b))
        list(a = pmax(ab[1, ], 0), b = pmax(ab[2, ], 0))
    }

    ## "57, 71" / "57;71" / "57 71" -> c(57L, 71L). Empty or unparseable -> integer(0).
    gcms_parse_ions <- function(txt) {
        if (is.null(txt) || length(txt) == 0 || is.na(txt[1])) return(integer(0))
        parts <- strsplit(trimws(as.character(txt[1])), "[,;[:space:]]+")[[1]]
        parts <- parts[nzchar(parts)]
        v <- suppressWarnings(as.numeric(parts))
        if (length(v) == 0 || any(!is.finite(v)) || any(v <= 0)) return(integer(0))
        unique(as.integer(round(v)))
    }

    ## Build the traces gcms_split_coeluting_n() needs for one sample's hump.
    ##
    ##   chrom_rows  that sample's rows of chromatograms_updated (long format:
    ##               ion == 0 is the TIC, ion == "baseline" the baseline)
    ##   allions     that sample's all-ion store (rt, mz, abundance)
    ##   start, end  the hump, in the sample's NATIVE rt
    ##   ion_sets    list of integer ion vectors, one per compound
    ##
    ## tic_net is TIC minus baseline scan by scan, so its sum is exactly the
    ## v6 `area` of [start, end]. Each diagnostic trace gets its own baseline:
    ## the median of the trace over a flank as wide as the hump on each side,
    ## interpolated across, the same idea as peak_ion_area(). The traces are
    ## put onto the TIC's scans (the grids are the same acquisition, so this
    ## is a lookup in practice, not a resampling). Returns list(rt, tic_net,
    ## base, ions), `ions` one trace per ion set.
    gcms_split_traces <- function(chrom_rows, allions, start, end, ion_sets) {
        ion_col <- as.character(chrom_rows$ion)
        in_w <- chrom_rows$rt >= start & chrom_rows$rt <= end
        tic <- chrom_rows[in_w & ion_col == "0", c("rt", "abundance")]
        bl  <- chrom_rows[in_w & ion_col == "baseline", c("rt", "abundance")]
        if (nrow(tic) < 5) return(NULL)
        tic <- tic[order(tic$rt), ]
        ## Measure everything from the hump's actual first and last scans. A
        ## Shift+G hump's bounds fall between scans; refresh later sees the
        ## snapped bounds, and flanks sized from the unsnapped ones moved the
        ## areas ~0.05% on the first Shift+Q after a Split.
        start <- min(tic$rt); end <- max(tic$rt)
        base <- bl$abundance[match(tic$rt, bl$rt)]
        base[is.na(base)] <- 0
        tic_net <- tic$abundance - base

        if (is.null(allions) || nrow(allions) == 0 || all(allions$mz == 0)) return(NULL)
        w <- end - start
        lo <- start - w; hi <- end + w
        near <- allions[allions$rt >= lo & allions$rt <= hi, ]
        near <- data.frame(rt = near$rt, mz = near$mz, abundance = near$abundance)
        grid <- sort(unique(near$rt))
        if (length(grid) < 5) return(NULL)

        trace <- function(ions) {
            sub <- near[near$mz %in% ions, ]
            y <- numeric(length(grid))
            if (nrow(sub)) {
                s <- tapply(sub$abundance, sub$rt, sum)
                y[match(as.numeric(names(s)), grid)] <- as.numeric(s)
            }
            left  <- grid < start; right <- grid > end
            inside <- !left & !right
            if (any(left) && any(right)) {
                bl_y <- stats::approx(c(stats::median(grid[left]), stats::median(grid[right])),
                                      c(stats::median(y[left]),    stats::median(y[right])),
                                      xout = grid, rule = 2)$y
            } else if (any(left)) {
                bl_y <- rep(stats::median(y[left]), length(grid))
            } else if (any(right)) {
                bl_y <- rep(stats::median(y[right]), length(grid))
            } else {
                bl_y <- rep(min(y[inside]), length(grid))
            }
            net <- pmax(y - bl_y, 0)
            stats::approx(grid, net, xout = tic$rt, rule = 2)$y
        }
        list(rt = tic$rt, tic_net = tic_net, base = base, ions = lapply(ion_sets, trace))
    }

    ## The Split button names the hump by its peak number. Numbers differ between
    ## samples (one extra early peak shifts everything after it), so the number is
    ## turned into a window in ALIGNED rt - the median start and end of peak `n`
    ## across the samples - and each sample's hump is then found by time, not by
    ## number. NULL when no sample has a peak `n`.
    gcms_window_for_peak <- function(peak_table, n) {
        if (is.null(n) || length(n) != 1 || !is.finite(n) ||
            !("peak_number_within_sample" %in% names(peak_table))) return(NULL)
        r <- peak_table[!is.na(peak_table$peak_number_within_sample) &
                        peak_table$peak_number_within_sample == round(n) &
                        is.finite(peak_table$peak_start_rt_offset) &
                        is.finite(peak_table$peak_end_rt_offset), , drop = FALSE]
        if (nrow(r) == 0) return(NULL)
        c(stats::median(r$peak_start_rt_offset), stats::median(r$peak_end_rt_offset))
    }

    ## A split's rows are labelled A, B (and C for three compounds), in the
    ## order of the ion fields.
    .gcms_split_labels <- c("A", "B", "C")

    ## Columns a split writes onto its rows. Everything needed to redo the split
    ## later lives on the rows themselves (ions, response factors, an id), so
    ## Shift+Q can refresh it without the sidebar, and two different humps can
    ## be split with different ions. split_ions_c / split_k_c are blank for a
    ## two-compound split; that blank is what says how many rows it has.
    .gcms_split_cols <- list(coelution_split = NA_character_, split_id = NA_character_,
                             split_status = NA_character_, split_ions_a = NA_character_,
                             split_ions_b = NA_character_, split_ions_c = NA_character_,
                             split_k_a = NA_real_, split_k_b = NA_real_, split_k_c = NA_real_,
                             area_tic_window = NA_real_, split_fit_r2 = NA_real_, split_k_dev = NA_real_)
    .gcms_ensure_split_cols <- function(pt) {
        for (col in names(.gcms_split_cols)) {
            if (!(col %in% names(pt))) pt[[col]] <- rep(.gcms_split_cols[[col]], nrow(pt))
            else if (is.character(.gcms_split_cols[[col]])) pt[[col]] <- as.character(pt[[col]])
            else pt[[col]] <- as.numeric(pt[[col]])
        }
        pt
    }

    ## How many compounds (rows) the split that row `r` belongs to has: three
    ## when it stores a third ion set, else two.
    .gcms_split_size <- function(pt, r) {
        ic <- if ("split_ions_c" %in% names(pt)) as.character(pt$split_ions_c[r]) else NA_character_
        if (length(ic) && !is.na(ic[1]) && nzchar(trimws(ic[1]))) 3L else 2L
    }

    ## The Split ion fields -> a list of ion sets (two, or three when the third
    ## field is filled), or a character message saying what is wrong.
    gcms_split_ion_sets <- function(a, b, c = NULL) {
        sets <- list(gcms_parse_ions(a), gcms_parse_ions(b))
        if (!is.null(c) && length(c) && !is.na(c[1]) && nzchar(trimws(c[1]))) {
            ic <- gcms_parse_ions(c)
            if (length(ic) == 0) return("Split: the ions for compound 3 could not be read; give numbers, or leave it empty for two compounds.")
            sets[[3]] <- ic
        }
        if (any(lengths(sets) == 0) || anyDuplicated(unlist(sets)))
            return("Split: give ions for compounds 1 and 2 (and 3 if there is one), with no ion in two compounds.")
        sets
    }

    ## Rewrite one split's rows from a finished split. EVERY row keeps the whole
    ## hump's bounds (Lucas, 2026-10-16): the shares are symmetric in the
    ## compounds, so nothing needs to know which elutes first - an earlier
    ## version cut the rows at a drop line, which needed an elution order the
    ## data often cannot give. `area` = each compound's split area (0 for a
    ## compound with no signal in this sample). Existing rows of the split keep
    ## their own other columns (a peak_ID typed in, say), matched by label.
    ## Returns one row per compound.
    .gcms_split_rows <- function(tmpl, tr, s) {
        m <- length(s$areas)
        labs <- .gcms_split_labels[seq_len(m)]
        idx <- rep(1L, m)
        if ("coelution_split" %in% names(tmpl)) {
            hit <- match(labs, as.character(tmpl$coelution_split))
            idx[!is.na(hit)] <- hit[!is.na(hit)]
        }
        tmpl <- tmpl[idx, , drop = FALSE]
        tmpl$peak_start <- rep(min(tr$rt), m)
        tmpl$peak_end   <- rep(max(tr$rt), m)
        tmpl$peak_start_rt_offset <- tmpl$peak_start + tmpl$rt_offset
        tmpl$peak_end_rt_offset   <- tmpl$peak_end   + tmpl$rt_offset
        tmpl$area <- s$areas
        if ("area_peak_ions" %in% names(tmpl)) tmpl$area_peak_ions <- NA_real_
        if ("peak_ions" %in% names(tmpl))      tmpl$peak_ions      <- NA_character_
        ## A whole-hump spectrum belongs to no one compound, so no agreement score.
        for (cc in c("spectrum_agreement", "spectrum_signal_frac"))
            if (cc %in% names(tmpl)) tmpl[[cc]] <- NA_real_
        tmpl$coelution_split <- labs
        tmpl$split_status    <- "split"
        tmpl$area_tic_window <- sum(s$areas)
        tmpl$split_fit_r2    <- s$fit_r2
        rownames(tmpl) <- NULL
        tmpl
    }

    ## Renumber within each sample after rows were replaced.
    .gcms_renumber <- function(pt) {
        cs <- if ("coelution_split" %in% names(pt)) as.character(pt$coelution_split) else rep(NA_character_, nrow(pt))
        pt <- pt[order(as.character(pt$path_to_cdf_csv), pt$peak_start, cs, na.last = TRUE), , drop = FALSE]
        if ("peak_number_within_sample" %in% names(pt))
            pt$peak_number_within_sample <- stats::ave(seq_len(nrow(pt)),
                                                       as.character(pt$path_to_cdf_csv), FUN = seq_along)
        rownames(pt) <- NULL
        pt
    }

    ## Largest |x| over the finite entries; NA when there are none.
    .gcms_max_dev <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(abs(x)) else NA_real_ }

    ## SPLIT (the sidebar button). In each sample the hump is the ONE peak
    ## touching `window`; two or more touching is refused, never merged (merging
    ## swallowed a neighbouring peak and the gap before it into B). A split
    ## already made with this same id is re-split (new ions, or a third compound,
    ## say); a split with a different id is left alone and reported.
    ##
    ##   get_traces  function(sample, start, end, ion_sets) -> gcms_split_traces()
    ##   window      c(start, end) in ALIGNED rt
    ##   ion_sets    list of two or three integer ion vectors, one per compound
    ##   k_mode      "sample": each sample fits its own k per compound (TIC per unit
    ##               diagnostic signal); "batch": one set of k across every sample's
    ##               hump, with samples whose own fit is > k_dev_flag off excluded and
    ##               refitted. The compound RATIO is free per sample in both modes:
    ##               it comes from each sample's own ion signals. k only converts ion
    ##               signal to TIC.
    ##
    ## A sample in which a compound's ions have no signal at all still gets that
    ## compound's row, with area 0, so every sample carries the same peaks and
    ## the numbers line up. A compound with no signal in ANY sample refuses the
    ## split (its ions are wrong). Only the split's rows are touched; every other
    ## row, and every other split, is left exactly as it was. Returns
    ## list(peak_table, log).
    gcms_apply_coelution_split <- function(peak_table, get_traces, window, ion_sets,
                                           crosstalk = c(0, 0), k_mode = c("sample", "batch"),
                                           split_id = NULL, k_dev_flag = 0.25, r2_flag = 0.98) {
        k_mode <- match.arg(k_mode)
        log <- character(0)
        say <- function(...) log <<- c(log, paste0(...))
        peak_table <- .gcms_ensure_split_cols(peak_table)
        out <- function() list(peak_table = peak_table, log = log)
        if (nrow(peak_table) == 0) return(out())
        m <- length(ion_sets)
        if (!(m %in% 2:3)) { say("Split: give ions for two or three compounds."); return(out()) }
        if (m == 3 && any(crosstalk != 0)) { say("Split: cross-talk is supported for two compounds only."); return(out()) }
        labs <- .gcms_split_labels[seq_len(m)]
        ws <- min(window); we <- max(window)
        if (is.null(split_id)) split_id <- sprintf("split@%.2f", ws)
        ions_txt <- vapply(ion_sets, paste, character(1), collapse = ", ")

        smp_col <- as.character(peak_table$path_to_cdf_csv)
        ps <- peak_table$peak_start_rt_offset; pe <- peak_table$peak_end_rt_offset
        no_off <- unique(smp_col[is.na(ps) | is.na(pe)])
        for (smp in no_off) say("  skipped ", basename(smp), ": no rt_offset, so the window cannot be placed")
        touches <- !is.na(ps) & !is.na(pe) & pe >= ws & ps <= we
        lab <- peak_table$coelution_split %in% .gcms_split_labels

        groups <- list()
        for (smp in setdiff(unique(smp_col[touches]), no_off)) {
            here <- which(smp_col == smp)
            mine <- here[lab[here] & peak_table$split_id[here] %in% split_id]
            other <- here[touches[here] & lab[here] & !(peak_table$split_id[here] %in% split_id)]
            if (length(other)) {
                say("  skipped ", basename(smp), ": the window touches peaks already split as ",
                    paste(unique(peak_table$split_id[other]), collapse = ", "), "; unsplit those first")
                next
            }
            if (length(mine)) {
                if (length(mine) >= 2 && length(mine) == .gcms_split_size(peak_table, mine[1])) {
                    groups[[smp]] <- mine; next
                }
                say("  skipped ", basename(smp), if (length(mine) == 1)
                        ": only one row of its earlier split remains; delete it and integrate the hump as one peak"
                    else ": only some rows of its earlier split remain; delete them and integrate the hump as one peak")
                next
            }
            cand <- here[touches[here]]
            ## A peak that only grazes the window is not the hump: require it to
            ## overlap at least half of the shorter of itself and the window.
            if (length(cand) == 1) {
                ov <- min(pe[cand], we) - max(ps[cand], ws)
                if (ov < 0.5 * min(pe[cand] - ps[cand], we - ws)) {
                    say("  skipped ", basename(smp), ": its only peak near the window (",
                        round(ps[cand], 2), "-", round(pe[cand], 2), ") barely overlaps it")
                    next
                }
            }
            if (length(cand) > 1) {
                say("  skipped ", basename(smp), ": ", length(cand),
                    " peaks touch the window; integrate the hump as ONE peak")
                next
            }
            groups[[smp]] <- cand
        }
        if (length(groups) == 0) {
            say("Split: no sample has a single peak at ", round(ws, 2), "-", round(we, 2), "; nothing split.")
            return(out())
        }

        humps <- list()
        for (smp in names(groups)) {
            rows <- groups[[smp]]
            st <- min(peak_table$peak_start[rows]); en <- max(peak_table$peak_end[rows])
            tr <- try(get_traces(smp, st, en, ion_sets), silent = TRUE)
            if (identical(tr, "not_loaded")) { say("  skipped ", basename(smp), ": not loaded this session"); next }
            if (inherits(tr, "try-error") || is.null(tr)) {
                say("  skipped ", basename(smp), ": could not build the TIC/ion traces"); next
            }
            um <- .gcms_unmix(tr$ions, crosstalk)
            if (is.null(um)) {
                say("  skipped ", basename(smp), ": the diagnostic ions do not discriminate the compounds"); next
            }
            S <- do.call(cbind, um)
            present <- colSums(S) > 0
            if (!any(present)) {
                say("  skipped ", basename(smp), ": none of the diagnostic ions has signal in this sample's hump"); next
            }
            P <- which(present)
            r <- if (length(P) > 1) {
                cr <- suppressWarnings(stats::cor(S[, P, drop = FALSE]))
                if (any(!is.finite(cr))) NA_real_ else max(abs(cr[upper.tri(cr)]))
            } else 0
            own_p <- unname(stats::lm.fit(S[, P, drop = FALSE], tr$tic_net)$coefficients)
            own <- rep(NA_real_, m); own[P] <- own_p
            own_ok <- all(is.finite(own_p)) && all(own_p > 0)
            ## Per-sample fallback (Lucas, 2026-10-26): a compound the fit gives a
            ## non-positive factor has only background on its ions in this sample.
            ## Treat it as absent (area 0, no stored factor) and refit the rest, so
            ## the sample keeps every row and the peak numbers still line up.
            weak <- integer(0); own_fb <- if (own_ok) own else NULL
            if (!own_ok) {
                Pf <- P; kf <- own_p
                while (length(Pf) && !(all(is.finite(kf)) && all(kf > 0))) {
                    bad <- Pf[!is.finite(kf) | kf <= 0]
                    weak <- c(weak, bad); Pf <- setdiff(Pf, bad)
                    if (length(Pf)) kf <- unname(stats::lm.fit(S[, Pf, drop = FALSE], tr$tic_net)$coefficients)
                }
                if (length(Pf)) { own_fb <- rep(NA_real_, m); own_fb[Pf] <- kf }
            }
            humps[[smp]] <- list(rows = rows, tr = tr, S = S, present = present, r = r,
                                 own = if (own_ok) own else NULL, own_fb = own_fb, weak = sort(weak))
        }
        if (length(humps) == 0) { say("Split: no sample could be split."); return(out()) }
        seen <- Reduce(`|`, lapply(humps, `[[`, "present"))
        if (!all(seen)) {
            say("Split: compound ", paste(which(!seen), collapse = " and "),
                " has no signal on its ions in any sample; check its ions. Nothing split.")
            return(out())
        }

        fit_k <- function(hs) {
            S <- do.call(rbind, lapply(hs, `[[`, "S"))
            Y <- unlist(lapply(hs, function(h) h$tr$tic_net))
            if (any(colSums(S) <= 0)) return(NULL)
            cr <- suppressWarnings(stats::cor(S))
            if (any(!is.finite(cr)) || max(abs(cr[upper.tri(cr)])) > 0.999) return(NULL)
            k <- unname(stats::lm.fit(S, Y)$coefficients)
            if (any(!is.finite(k)) || any(k <= 0)) NULL else k
        }
        batch_k <- NULL
        if (k_mode == "batch") {
            batch_k <- fit_k(humps)
            if (is.null(batch_k)) {
                say("Split: the batch fit failed (ions do not discriminate, or a non-positive factor); check the ions and the hump.")
                return(out())
            }
            devs <- vapply(humps, function(h) if (is.null(h$own)) NA_real_ else .gcms_max_dev(h$own / batch_k - 1), numeric(1))
            outl <- names(humps)[is.na(devs) | devs > k_dev_flag]
            keep <- setdiff(names(humps), outl)
            if (length(outl) && length(keep)) {
                k2 <- fit_k(humps[keep])
                if (!is.null(k2)) {
                    batch_k <- k2
                    say("  excluded from the batch fit (own fit too far off): ", paste(basename(outl), collapse = ", "))
                }
            }
            say("Split: ", paste(sprintf("k%s = %.4g", labs, batch_k), collapse = ", "), " across the batch.")
        }

        ## The k are properties of the compounds, so per-sample fits should
        ## agree; a sample far from the median is the tell for an extra compound
        ## or a wrong window, even when its own fit looks fine.
        med_k <- NULL
        if (k_mode == "sample") {
            ks <- do.call(rbind, lapply(humps, function(h) if (!is.null(h$own)) h$own))
            if (!is.null(ks) && nrow(ks) >= 3) med_k <- apply(ks, 2, stats::median, na.rm = TRUE)
        }
        new_rows <- list(); drop_rows <- integer(0)
        for (smp in names(humps)) {
            h <- humps[[smp]]
            if (k_mode == "sample") {
                if (!is.finite(h$r) || h$r > 0.999 || is.null(h$own_fb)) {
                    say("  skipped ", basename(smp), ": this sample's own response factors could not be fitted (",
                        if (!is.finite(h$r) || h$r > 0.999) "the ions rise together" else "no compound's factor came out positive",
                        "); try the batch fit")
                    next
                }
                k <- h$own_fb
            } else k <- batch_k
            s <- gcms_split_coeluting_n(h$tr$rt, h$tr$tic_net, h$tr$ions, crosstalk = crosstalk,
                                        scale = k, allow_absent = TRUE)
            if (s$status != "ok") { say("  skipped ", basename(smp), ": ", s$status); next }
            tmpl <- .gcms_split_rows(peak_table[h$rows, , drop = FALSE], h$tr, s)
            tmpl$split_id <- split_id
            tmpl$split_ions_a <- ions_txt[1]; tmpl$split_ions_b <- ions_txt[2]
            tmpl$split_ions_c <- if (m == 3) ions_txt[3] else NA_character_
            tmpl$split_k_a <- k[1]; tmpl$split_k_b <- k[2]
            tmpl$split_k_c <- if (m == 3) k[3] else NA_real_
            tmpl$split_k_dev <- if (k_mode == "batch" && !is.null(h$own)) .gcms_max_dev(h$own / k - 1)
                                else if (k_mode == "sample" && !is.null(med_k)) .gcms_max_dev(k / med_k - 1) else NA_real_
            new_rows[[smp]] <- tmpl
            drop_rows <- c(drop_rows, h$rows)
            if (any(!h$present))
                say(sprintf("  %s: compound %s has no signal on its ions here, so it is kept as a peak with area 0",
                            basename(smp), paste(which(!h$present), collapse = " and ")))
            if (k_mode == "sample" && length(h$weak))
                say(sprintf("  CHECK %s: compound %s shows only background on its ions here (its fitted factor was not positive), so it is kept as a peak with area 0",
                            basename(smp), paste(h$weak, collapse = " and ")))
            if (k_mode == "sample")
                say(sprintf("  %s: %s, fit_r2 %.3f", basename(smp),
                            paste(sprintf("k%s %.4g", labs, k), collapse = ", "), s$fit_r2))
            if (k_mode == "sample" && !is.na(tmpl$split_k_dev[1]) && tmpl$split_k_dev[1] > k_dev_flag)
                say(sprintf("  CHECK %s: its factors are %.0f%% off the median across samples", basename(smp), 100 * tmpl$split_k_dev[1]))
            if ((k_mode == "batch" && (is.na(tmpl$split_k_dev[1]) || tmpl$split_k_dev[1] > k_dev_flag)) ||
                (is.finite(s$fit_r2) && s$fit_r2 < r2_flag))
                say(sprintf("  CHECK %s: fit_r2 %.3f%s", basename(smp), s$fit_r2,
                            if (k_mode == "batch") sprintf(", own-fit k %s off the batch k",
                                if (is.na(tmpl$split_k_dev[1])) "not computable and so" else sprintf("%.0f%%", 100 * tmpl$split_k_dev[1])) else ""))
        }
        if (length(drop_rows)) {
            peak_table <- peak_table[-drop_rows, , drop = FALSE]
            peak_table <- .gcms_renumber(rbind(peak_table, do.call(rbind, new_rows)))
        }
        say(sprintf("Split (%s, %d compounds): %d sample(s) split as %s.", k_mode, m, length(new_rows), split_id))
        out()
    }

    ## REFRESH (every Shift+Q). Redoes every existing split from what is stored on
    ## its rows - same ions, same response factors - over the split's current
    ## extent, after the area recompute has just given those rows plain drop-line
    ## areas. Never creates a split. Returns list(peak_table, profiles, log);
    ## `profiles` holds each split's per-scan stacked bands for the plot.
    gcms_refresh_splits <- function(peak_table, get_traces, crosstalk = c(0, 0)) {
        log <- character(0)
        say <- function(...) log <<- c(log, paste0(...))
        peak_table <- .gcms_ensure_split_cols(peak_table)
        profiles <- list()
        lab <- which(peak_table$coelution_split %in% .gcms_split_labels)
        if (length(lab) == 0) return(list(peak_table = peak_table, profiles = NULL, log = log))
        key <- paste(peak_table$path_to_cdf_csv[lab], peak_table$split_id[lab], sep = "\r")
        new_rows <- list(); drop_rows <- integer(0)
        for (g in unique(key)) {
            rows <- lab[key == g]
            rows <- rows[order(peak_table$coelution_split[rows])]
            smp <- as.character(peak_table$path_to_cdf_csv[rows[1]])
            m <- .gcms_split_size(peak_table, rows[1])
            labs <- .gcms_split_labels[seq_len(m)]
            if (length(rows) != m || !setequal(peak_table$coelution_split[rows], labs)) {
                peak_table$split_status[rows] <- if (length(rows) == 1) "not split: only one row of the pair remains"
                                                 else "not split: only some rows of the split remain"
                say("  ", basename(smp), ": only ", if (length(rows) == 1) "one row" else "some rows",
                    " of split ", peak_table$split_id[rows[1]],
                    " remain", if (length(rows) == 1) "s" else "", "; delete ",
                    if (length(rows) == 1) "it" else "them", " and integrate the hump as one peak")
                next
            }
            k <- c(peak_table$split_k_a[rows[1]], peak_table$split_k_b[rows[1]], peak_table$split_k_c[rows[1]])[seq_len(m)]
            ion_sets <- lapply(c(peak_table$split_ions_a[rows[1]], peak_table$split_ions_b[rows[1]],
                                 peak_table$split_ions_c[rows[1]])[seq_len(m)], gcms_parse_ions)
            fail <- function(why) {
                peak_table$split_status[rows] <<- paste("not split:", why)
                ## Every row spans the whole hump, so the recompute has just given
                ## EACH of them the whole hump's area. Blank them, so nothing that
                ## sums the column counts the hump more than once.
                peak_table$area[rows] <<- NA_real_
                say("  ", basename(smp), ": split ", peak_table$split_id[rows[1]], " not refreshed (", why,
                    "); its areas are blank until it is fixed or unsplit")
            }
            if (!any(is.finite(k)) || any(lengths(ion_sets) == 0)) { fail("its stored settings are missing"); next }
            tr <- try(get_traces(smp, min(peak_table$peak_start[rows]), max(peak_table$peak_end[rows]), ion_sets), silent = TRUE)
            ## Not loaded this session (samples_monolist_subset): its areas were not
            ## recomputed either, so leave the split exactly as it was.
            if (identical(tr, "not_loaded")) next
            if (inherits(tr, "try-error") || is.null(tr)) { fail("no TIC/ion traces"); next }
            s <- gcms_split_coeluting_n(tr$rt, tr$tic_net, tr$ions, crosstalk = crosstalk,
                                        scale = k, allow_absent = TRUE)
            if (s$status != "ok") { fail(s$status); next }
            tmpl <- .gcms_split_rows(peak_table[rows, , drop = FALSE], tr, s)
            new_rows[[g]] <- tmpl
            drop_rows <- c(drop_rows, rows)
            ## Stacked bands: A from the baseline up to its share, B on top of it,
            ## then C up to the TIC, so every part of the hump is coloured by its
            ## owner. A compound with no signal here is a band of zero height.
            b0 <- if (is.null(tr$base)) rep(0, length(tr$rt)) else tr$base
            cum <- b0 + tr$tic_net * t(apply(s$shares, 1, cumsum))
            lo <- cbind(b0, cum[, -m, drop = FALSE])
            profiles[[g]] <- data.frame(
                path_to_cdf_csv = smp, split_id = tmpl$split_id[1],
                rt_rt_offset = rep(tr$rt + tmpl$rt_offset[1], m),
                part = rep(labs, each = length(tr$rt)),
                ymin = as.vector(lo), top = as.vector(cum),
                stringsAsFactors = FALSE)
        }
        if (length(drop_rows)) {
            peak_table <- peak_table[-drop_rows, , drop = FALSE]
            peak_table <- .gcms_renumber(rbind(peak_table, do.call(rbind, new_rows)))
        }
        list(peak_table = peak_table,
             profiles = if (length(profiles)) do.call(rbind, profiles) else NULL, log = log)
    }

    ## UNSPLIT (the sidebar button). Every split that peak `n` belongs to, in
    ## any sample, goes back to one peak spanning all its rows, its area the sum
    ## of theirs (the whole hump's TIC area). Returns list(peak_table, log).
    gcms_unsplit_peak <- function(peak_table, n) {
        log <- character(0)
        peak_table <- .gcms_ensure_split_cols(peak_table)
        hit <- !is.na(peak_table$peak_number_within_sample) & peak_table$peak_number_within_sample == round(n) &
               peak_table$coelution_split %in% .gcms_split_labels
        ids <- unique(peak_table$split_id[hit])
        if (length(ids) == 0)
            return(list(peak_table = peak_table, log = sprintf("Unsplit: peak %s is not part of a split.", n)))
        merged <- list(); drop_rows <- integer(0)
        rows_all <- which(peak_table$split_id %in% ids & peak_table$coelution_split %in% .gcms_split_labels)
        key <- paste(peak_table$path_to_cdf_csv[rows_all], peak_table$split_id[rows_all], sep = "\r")
        for (g in unique(key)) {
            rows <- rows_all[key == g]
            one <- peak_table[rows[which.min(peak_table$peak_start[rows])], , drop = FALSE]
            one$peak_end <- max(peak_table$peak_end[rows])
            one$peak_end_rt_offset <- one$peak_end + one$rt_offset
            one$area <- sum(peak_table$area[rows])
            for (col in names(.gcms_split_cols)) one[[col]] <- .gcms_split_cols[[col]]
            merged[[g]] <- one
            drop_rows <- c(drop_rows, rows)
        }
        peak_table <- peak_table[-drop_rows, , drop = FALSE]
        peak_table <- .gcms_renumber(rbind(peak_table, do.call(rbind, merged)))
        list(peak_table = peak_table,
             log = sprintf("Unsplit: %s merged back to one peak in %d sample(s).",
                           paste(ids, collapse = ", "), length(merged)))
    }

    ## Labels for a split's zero-area peaks (a compound with no signal in that
    ## sample: its stacked band has no height). `prof` is the band table as the
    ## plot builds it (path_to_cdf_csv, split_id, part, rt_rt_offset, ymin, top,
    ## peak_number_within_sample). One row per split with any zero-area part: x
    ## at the middle of the hump, y at the lowest point of its bottom edge (the
    ## text is hung below it), label the peak number(s), e.g. "4" or "4, 5".
    gcms_split_zero_labels <- function(prof) {
        empty <- data.frame(path_to_cdf_csv = character(0), rt_rt_offset = numeric(0),
                            y = numeric(0), label = character(0), stringsAsFactors = FALSE)
        if (is.null(prof) || nrow(prof) == 0) return(empty)
        split_key <- paste(prof$path_to_cdf_csv, prof$split_id, sep = "\r")
        out <- lapply(split(prof, split_key), function(g) {
            thick <- tapply(g$top - g$ymin, g$part, function(x) max(abs(x)))
            zero <- names(thick)[thick <= 0]
            if (!length(zero)) return(NULL)
            nums <- unique(g$peak_number_within_sample[g$part %in% zero])
            nums <- sort(nums[!is.na(nums)])
            if (!length(nums)) return(NULL)
            data.frame(path_to_cdf_csv = g$path_to_cdf_csv[1],
                       rt_rt_offset = mean(range(g$rt_rt_offset)),
                       y = min(g$ymin), label = paste(nums, collapse = ", "),
                       stringsAsFactors = FALSE)
        })
        out <- do.call(rbind, out)
        if (is.null(out)) empty else out
    }

    ## Hard-wrap labels every `width` characters, mid-word if need be (sample
    ## names rarely have spaces, so a word wrap would not wrap them at all).
    gcms_wrap_hard <- function(x, width = 14) {
        x <- as.character(x)
        if (!is.finite(width) || width < 1) return(x)
        vapply(x, function(s) {
            if (is.na(s) || nchar(s) <= width) return(s)
            starts <- seq(1, nchar(s), by = width)
            paste(substring(s, starts, pmin(starts + width - 1, nchar(s))), collapse = "\n")
        }, character(1), USE.NAMES = FALSE)
    }

    ## The peak table on screen. read.csv() reads a column with no values in it
    ## (area_peak_ions on a TIC-only folder, every split column before the first
    ## split) as logical, and rhandsontable draws a logical column as tick boxes.
    ## Show such a column as blank text instead. Display only: the file is not
    ## touched.
    gcms_peak_table_display <- function(pt) {
        for (col in names(pt))
            if (is.logical(pt[[col]]) && all(is.na(pt[[col]]))) pt[[col]] <- rep("", nrow(pt))
        pt
    }

    ## Shift+Z. peak_ID is the only column edited by hand, so the save takes the
    ## file as it is on disk and copies in only the on-screen peak_ID column.
    ## Areas, ions and every split setting therefore never pass through the
    ## widget. Returns NULL when the two tables do not line up row for row.
    gcms_apply_peak_id_edits <- function(on_disk, edited) {
        if (is.null(on_disk) || is.null(edited) || nrow(on_disk) != nrow(edited) ||
            !("peak_ID" %in% names(on_disk)) || !("peak_ID" %in% names(edited))) return(NULL)
        on_disk$peak_ID <- as.character(edited$peak_ID)
        on_disk
    }

    ## The TIC view's peak overlay rows (shaded ribbons, peak-start lines, peak-number
    ## positions) for every peak on the page, in ONE pass: a non-equi join of the peak
    ## windows onto the page's rows. It replaces a loop that filtered each sample's rows
    ## once per peak (~5 s of an ~11 s redraw at 95 samples, ~0.02 s this way). Same
    ## rows, order, columns and types as that loop, including its skip message for a
    ## peak whose window has TIC rows but no matching baseline rows. `cuf` is the
    ## page's rows BEFORE the facet factor is applied. Pinned against the frozen loop
    ## by test_overlay_onepass.R.
    ## Thin a page's line data to what the screen can show (2026-10-21). For each
    ## trace (sample x ion) the window [xs, xe] is cut into `n_bins` columns (about
    ## one per pixel) and each column keeps only its lowest and highest points, so the
    ## drawn line covers the same pixels (0.02% of pixels differ, by antialiasing, on
    ## a 20-sample page) at half the draw time, while a full-range
    ## page draws a few thousand points per trace instead of every scan (a 95-sample
    ## page held ~2.5M). Rows are kept as they are (no new values) and stay in
    ## their original order. A page with no more than 4 points per column overall
    ## is returned untouched, so a zoomed-in view is drawn exactly as before. Only
    ## the LINE layers use this; the peak ribbons and every area read the full rows.
    ## test_render_thin.R.
    gcms_thin_lines <- function(d, xs, xe, n_bins = 2000L, x = "rt_rt_offset", y = "abundance",
                                by = c("path_to_cdf_csv", "ion")) {
        if (is.null(d) || nrow(d) == 0 || !is.finite(xs) || !is.finite(xe) || xe <= xs ||
            !is.finite(n_bins) || n_bins < 1) return(d)
        if (nrow(d) <= 4 * n_bins) return(d)
        xv <- d[[x]]; yv <- d[[y]]
        bin <- pmin(floor((xv - xs) / (xe - xs) * n_bins), n_bins - 1)
        ## One integer per (trace, column): trace codes from the grouping columns.
        tr <- rep(0, nrow(d))
        for (b in intersect(by, names(d))) {
            f <- as.integer(factor(d[[b]], exclude = NULL)); tr <- tr * (max(f) + 1) + f
        }
        gid <- tr * (n_bins + 1) + bin
        ok <- !is.na(gid) & !is.na(xv) & !is.na(yv)
        keep <- !ok                                   # anything unbinnable is kept as is
        w <- which(ok)
        ## Within each column: its lowest and highest point.
        o <- w[order(gid[w], yv[w], w)]
        g <- gid[o]
        keep[o[!duplicated(g)]] <- TRUE
        keep[o[!duplicated(g, fromLast = TRUE)]] <- TRUE
        d[keep, , drop = FALSE]
    }

    gcms_overlay_rows <- function(cuf, peak_table) {
        none <- list(ribbons = dplyr::bind_rows(list()), vlines = dplyr::bind_rows(list()),
                     labels = dplyr::bind_rows(list()))
        if (is.null(cuf) || nrow(cuf) == 0 || is.null(peak_table) || nrow(peak_table) == 0) return(none)
        s <- peak_table$peak_start_rt_offset; e <- peak_table$peak_end_rt_offset
        keep_pk <- which(!is.na(s) & !is.na(e))
        if (!length(keep_pk)) return(none)
        cd <- data.table::data.table(i = seq_len(nrow(cuf)), p = as.character(cuf$path_to_cdf_csv),
                                     rt = cuf$rt_rt_offset)
        cd <- cd[!is.na(rt)]
        pt <- data.table::data.table(k = keep_pk, p = as.character(peak_table$path_to_cdf_csv)[keep_pk],
                                     s = s[keep_pk], e = e[keep_pk])
        hits <- cd[pt, on = list(p, rt >= s, rt <= e), nomatch = 0L, list(k = i.k, i = x.i)]
        if (!nrow(hits)) return(none)
        data.table::setorderv(hits, c("k", "i"))
        ion <- cuf$ion[hits$i]
        hits[, tic := ion == 0 & !is.na(ion)]
        hits[, bl := ion == "baseline" & !is.na(ion)]
        per <- hits[, list(n_tic = sum(tic), n_bl = sum(bl), first_i = i[1]), by = k]
        bad <- per[n_tic == 0 | n_tic != n_bl]
        for (j in seq_len(nrow(bad))) {
            cat(paste0("Skipping peak overlay for ", peak_table[bad$k[j], ]$path_to_cdf_csv,
                       ", peak ", peak_table$peak_number_within_sample[bad$k[j]], ": ",
                       bad$n_tic[j], " TIC rows vs ", bad$n_bl[j],
                       " baseline rows. Press Shift+Q to rebuild the baseline.\n"))
        }
        good <- per[n_tic > 0 & n_tic == n_bl]
        if (!nrow(good)) return(none)
        hg <- hits[k %in% good$k]
        pnum <- peak_table$peak_number_within_sample
        rib_i <- hg[tic == TRUE]
        rows_of <- function(idx) if (data.table::is.data.table(cuf)) cuf[idx] else cuf[idx, , drop = FALSE]
        ribbons <- rows_of(rib_i$i)
        ribbons$peak_number_within_sample <- pnum[rib_i$k]
        ribbons$baseline <- cuf$abundance[hg[bl == TRUE]$i]
        vlines <- rows_of(good$first_i)
        vlines$peak_number_within_sample <- pnum[good$k]
        lab <- data.table::data.table(k = rib_i$k, rt = ribbons$rt_rt_offset, ab = ribbons$abundance,
                                      p = ribbons$path_to_cdf_csv)
        lab <- lab[, list(x = stats::median(rt), y = max(ab), path_to_cdf_csv = p[1]), by = k]
        labels <- data.frame(peak_number_within_sample = pnum[lab$k], x = lab$x, y = lab$y,
                             path_to_cdf_csv = lab$path_to_cdf_csv, stringsAsFactors = FALSE)
        if (!data.table::is.data.table(cuf)) { rownames(ribbons) <- NULL; rownames(vlines) <- NULL }
        list(ribbons = ribbons, vlines = vlines, labels = labels)
    }

    ## TRUE when the browser's copy of the peak table carries peak_ID edits that
    ## the table it was drawn from does not. Only peak_ID is editable. Tables of
    ## different lengths are not comparable (the browser is still showing an older
    ## render), so they count as no edits rather than flashing the flag.
    gcms_peak_ids_edited <- function(edited, snapshot) {
        if (is.null(edited) || is.null(snapshot) || !("peak_ID" %in% names(edited)) ||
            !("peak_ID" %in% names(snapshot)) || nrow(edited) != nrow(snapshot)) return(FALSE)
        norm <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; trimws(x) }
        any(norm(edited$peak_ID) != norm(snapshot$peak_ID))
    }

    ## Move one sample along the retention-time axis (the Align strip, 2026-10-21):
    ## add `delta` to its rt_offset and carry it through every table that stores an
    ## offset-corrected time, so the redraw is right with no Shift+Q. Raw rt and the
    ## integration windows (peak_start / peak_end) never change: the sample's
    ## chromatogram and its peaks move together, so no area changes. Any table may
    ## be NULL. Returns the updated tables plus the sample's new offset.
    gcms_shift_sample_offset <- function(sample, delta, samples_ml, peaks = NULL, chrom = NULL, profiles = NULL) {
        if (is.null(samples_ml) || !(sample %in% samples_ml$path_to_cdf_csv) || !is.finite(delta))
            return(NULL)
        r <- samples_ml$path_to_cdf_csv == sample
        old <- suppressWarnings(as.numeric(samples_ml$rt_offset[r][1])); if (!is.finite(old)) old <- 0
        new <- round(old + delta, 6)
        samples_ml$rt_offset[r] <- new
        if (!is.null(peaks) && nrow(peaks) && "path_to_cdf_csv" %in% names(peaks)) {
            k <- as.character(peaks$path_to_cdf_csv) == sample
            if (any(k)) {
                peaks$rt_offset[k] <- new
                if ("peak_start_rt_offset" %in% names(peaks)) peaks$peak_start_rt_offset[k] <- peaks$peak_start[k] + new
                if ("peak_end_rt_offset" %in% names(peaks))   peaks$peak_end_rt_offset[k]   <- peaks$peak_end[k] + new
            }
        }
        if (!is.null(chrom) && nrow(chrom)) {
            k <- as.character(chrom$path_to_cdf_csv) == sample
            if (any(k)) { chrom$rt_offset[k] <- new; chrom$rt_rt_offset[k] <- chrom$rt[k] + new }
        }
        if (!is.null(profiles) && nrow(profiles) && "rt_rt_offset" %in% names(profiles)) {
            k <- as.character(profiles$path_to_cdf_csv) == sample
            if (any(k)) profiles$rt_rt_offset[k] <- profiles$rt_rt_offset[k] + (new - old)
        }
        list(samples = samples_ml, peaks = peaks, chrom = chrom, profiles = profiles, offset = new)
    }

    ## Align every sample on one peak (the sidebar's "Align on a peak", 2026-10-22).
    ## Each sample's peak `n` is found by ITS OWN peak number (the number drawn on
    ## its panel), and its centre is the middle of its integration window in raw
    ## rt (Lucas: window centre, not apex). The target is the median of those
    ## centres, and each sample's offset is SET to target - centre: it replaces
    ## whatever offset the sample had, so pressing Align twice gives the same
    ## answer, and the per-panel arrows then nudge on top of it. The two rows of a
    ## split pair carry consecutive numbers but span the same hump, so either number
    ## aligns on the hump (and in a sample with a split, every later number is one
    ## higher than in an unsplit sample, as drawn). Samples with no peak `n`
    ## are left as they are and listed in `missing`. Returns NULL when no sample has
    ## a peak `n`.
    gcms_align_offsets_for_peak <- function(peak_table, samples_ml, n) {
        if (is.null(n) || length(n) != 1 || !is.finite(n) || is.null(peak_table) || !nrow(peak_table) ||
            !all(c("peak_number_within_sample", "peak_start", "peak_end", "path_to_cdf_csv") %in% names(peak_table)))
            return(NULL)
        r <- peak_table[!is.na(peak_table$peak_number_within_sample) &
                        peak_table$peak_number_within_sample == round(n) &
                        is.finite(peak_table$peak_start) & is.finite(peak_table$peak_end), , drop = FALSE]
        if (!nrow(r)) return(NULL)
        ctr <- tapply((r$peak_start + r$peak_end) / 2, as.character(r$path_to_cdf_csv), mean)
        ctr <- ctr[names(ctr) %in% as.character(samples_ml$path_to_cdf_csv)]
        if (!length(ctr)) return(NULL)
        target <- stats::median(ctr)
        offsets <- stats::setNames(round(target - as.numeric(ctr), 6), names(ctr))
        list(offsets = offsets, target = target,
             missing = setdiff(as.character(samples_ml$path_to_cdf_csv), names(ctr)))
    }

    ## Set several samples' offsets at once (Align on a peak, Reset all offsets):
    ## `offsets` is a named vector, sample -> new rt_offset. The same tables and the
    ## same rule as gcms_shift_sample_offset() (raw rt and integration windows never
    ## change, so no area moves), done in one vectorised pass rather than one pass
    ## over the chromatogram table per sample. Samples not named are untouched.
    gcms_set_sample_offsets <- function(offsets, samples_ml, peaks = NULL, chrom = NULL, profiles = NULL) {
        if (is.null(samples_ml) || !length(offsets) || is.null(names(offsets))) return(NULL)
        offsets <- offsets[is.finite(offsets) & names(offsets) %in% as.character(samples_ml$path_to_cdf_csv)]
        if (!length(offsets)) return(NULL)
        old_all <- suppressWarnings(as.numeric(samples_ml$rt_offset)); old_all[!is.finite(old_all)] <- 0
        old <- stats::setNames(old_all, as.character(samples_ml$path_to_cdf_csv))[names(offsets)]
        i <- match(as.character(samples_ml$path_to_cdf_csv), names(offsets)); k <- !is.na(i)
        samples_ml$rt_offset[k] <- unname(offsets[i[k]])
        if (!is.null(peaks) && nrow(peaks) && "path_to_cdf_csv" %in% names(peaks)) {
            i <- match(as.character(peaks$path_to_cdf_csv), names(offsets)); k <- !is.na(i)
            if (any(k)) {
                peaks$rt_offset[k] <- unname(offsets[i[k]])
                if ("peak_start_rt_offset" %in% names(peaks)) peaks$peak_start_rt_offset[k] <- peaks$peak_start[k] + peaks$rt_offset[k]
                if ("peak_end_rt_offset" %in% names(peaks))   peaks$peak_end_rt_offset[k]   <- peaks$peak_end[k] + peaks$rt_offset[k]
            }
        }
        if (!is.null(chrom) && nrow(chrom)) {
            i <- match(as.character(chrom$path_to_cdf_csv), names(offsets)); k <- !is.na(i)
            if (any(k)) { chrom$rt_offset[k] <- unname(offsets[i[k]]); chrom$rt_rt_offset[k] <- chrom$rt[k] + chrom$rt_offset[k] }
        }
        if (!is.null(profiles) && nrow(profiles) && "rt_rt_offset" %in% names(profiles)) {
            i <- match(as.character(profiles$path_to_cdf_csv), names(offsets)); k <- !is.na(i)
            if (any(k)) profiles$rt_rt_offset[k] <- profiles$rt_rt_offset[k] + unname(offsets[i[k]] - old[i[k]])
        }
        list(samples = samples_ml, peaks = peaks, chrom = chrom, profiles = profiles,
             changed = names(offsets)[abs(offsets - old) > 1e-9])
    }

    ## Slide the x window by HALF its width (the Earlier / Later buttons, 2026-10-21;
    ## Lucas: 50% overlap, so 500-600 -> 550-650), so a peak on the edge of one view
    ## is whole in the next. Same width throughout. Clamped to the run [lo, hi] keeping the width, so the
    ## last step lands flush with the end rather than off it. Returns NULL when the
    ## window is unusable (not set yet, or zero width).
    gcms_slide_window <- function(xs, xe, direction, lo, hi) {
        if (length(xs) != 1 || length(xe) != 1 || !is.finite(xs) || !is.finite(xe) || xe <= xs) return(NULL)
        w <- xe - xs
        ns <- xs + sign(direction) * w / 2
        if (is.finite(lo) && is.finite(hi) && hi > lo) {
            if (w >= hi - lo) return(c(lo, hi))
            ns <- min(max(ns, lo), hi - w)
        }
        c(ns, ns + w)
    }

    ## View history (the Previous view / Subsequent view buttons, 2026-10-22; Lucas:
    ## zoom in quickly, then step back to the last zoom rather than all the way out).
    ## Browser-style: a list of views plus a position. A view is everything a zoom
    ## sets (x limits, the raw-count y pair, and the keyed y-zoom). Recording a view
    ## that differs from the one at the position drops anything after it, as a new
    ## page does in a browser, and the list is capped at `cap` views.
    gcms_view_history_record <- function(hist, view, cap = 50L) {
        if (is.null(hist)) hist <- list(views = list(), pos = 0L)
        if (is.null(view) || length(view$xs) != 1 || length(view$xe) != 1 ||
            !is.finite(view$xs) || !is.finite(view$xe)) return(hist)
        if (hist$pos > 0 && identical(hist$views[[hist$pos]], view)) return(hist)
        views <- c(hist$views[seq_len(hist$pos)], list(view))
        if (length(views) > cap) views <- views[(length(views) - cap + 1):length(views)]
        list(views = views, pos = length(views))
    }
    ## One step back (-1) or forward (+1). NULL when there is nowhere to go.
    gcms_view_history_step <- function(hist, direction) {
        if (is.null(hist) || !length(hist$views)) return(NULL)
        p <- hist$pos + sign(direction)
        if (p < 1 || p > length(hist$views)) return(NULL)
        list(views = hist$views, pos = as.integer(p))
    }

    ## Undo (the Undo button, 2026-10-26; Lucas: undo only, no redo, button only,
    ## recalculations are not steps). Everything an edit changes lives in these
    ## files in the working folder, so a step is a copy of their bytes taken just
    ## before the edit, and undoing writes them back. NULL = the file did not exist.
    gcms_undo_files <- c("peaks_monolist.csv", "samples_monolist.csv", "coelution_profiles.csv")
    gcms_undo_capture <- function(files = gcms_undo_files) {
        stats::setNames(lapply(files, function(f) {
            if (!file.exists(f)) return(NULL)
            readBin(f, "raw", n = file.info(f)$size)
        }), files)
    }
    ## Push a step unless the edit changed nothing. A step carrying a `coalesce`
    ## key equal to the top step's (a run of arrow nudges on one sample) is
    ## folded into it: the top step already holds the state before the run.
    gcms_undo_push <- function(stack, step, now, cap = 20L) {
        if (identical(step$files, now)) return(stack)
        top <- if (length(stack)) stack[[length(stack)]] else NULL
        if (!is.null(step$coalesce) && !is.null(top) && identical(top$coalesce, step$coalesce)) return(stack)
        stack <- c(stack, list(step))
        if (length(stack) > cap) stack <- stack[(length(stack) - cap + 1):length(stack)]
        stack
    }
    ## Write a step's files back, each atomically (temp file, then rename). The
    ## peak list is never deleted -- Shift+Q reads it unconditionally -- so a step
    ## taken before it existed restores an empty one, as Shift+0 writes.
    gcms_undo_restore <- function(files) {
        for (f in names(files)) {
            bytes <- files[[f]]
            if (is.null(bytes) && f == "peaks_monolist.csv")
                bytes <- charToRaw("\"peak_start\",\"peak_end\",\"peak_ID\",\"path_to_cdf_csv\",\"area\"\n")
            if (is.null(bytes)) {
                if (file.exists(f)) file.remove(f)
                next
            }
            tmp <- tempfile(pattern = paste0(".", basename(f), "."), tmpdir = dirname(f))
            writeBin(bytes, tmp)
            if (!file.rename(tmp, f)) { unlink(tmp); stop("Undo: could not write ", f) }
        }
        invisible(NULL)
    }

    ## ---- Peak spectra and assignment confidence (phase 3, 2026-10-21) ----------
    ## Every peak gets a square coloured by how well its mass
    ## spectrum agrees with the other peaks given the same peak_ID in this folder,
    ## so a mis-assigned or co-eluted peak stands out without opening each one.
    ##
    ## The spectrum is the lab library's extraction (Lucas, 2026-10-21): the mean
    ## of the apex scan and `apex_half` scans either side, minus the mean of
    ## `bg_scans` scans just outside each window edge, clipped at 0, on nominal
    ## m/z (the all-ion store's floor(mz + 0.3)) over 40-1000, scaled to max 100.
    ## Without the background step the column bleed (207/281) dominates a weak
    ## peak's spectrum and every weak peak looks like every other.
    ## `signal_frac` is the share of the apex signal left after subtraction.
    ## `store` is a sample's all-ion store (rt, mz, abundance; rt-keyed or not).
    ## Returns NULL when the window spans fewer than 3 scans or holds no signal.
    gcms_spectrum_mz <- 40:1000
    gcms_peak_spectrum <- function(store, start_rt, end_rt, apex_half = 2L, bg_scans = 3L,
                                   mz_axis = gcms_spectrum_mz, rts = NULL) {
        if (is.null(store) || !nrow(store) || !is.finite(start_rt) || !is.finite(end_rt) ||
            end_rt <= start_rt) return(NULL)
        ## `rts` (the sample's sorted unique scan times) can be passed in, so a
        ## caller scoring every peak of a sample computes it once, not per peak.
        if (is.null(rts)) rts <- sort(unique(store$rt))
        i0 <- findInterval(start_rt, rts, left.open = TRUE) + 1L     # first scan >= start
        i1 <- findInterval(end_rt, rts)                              # last scan <= end
        if (i1 - i0 < 2L) return(NULL)
        lo_scan <- max(1L, i0 - bg_scans); hi_scan <- min(length(rts), i1 + bg_scans)
        mz_lo <- min(mz_axis); mz_hi <- max(mz_axis); rt_lo <- rts[lo_scan]; rt_hi <- rts[hi_scan]
        d <- data.table::as.data.table(store)
        ## An rt-keyed store (load_allions()) is sliced by binary search, as
        ## peak_ion_area() does; scanning ~0.5M rows per peak would dominate.
        sub <- if (identical(data.table::key(d)[1], "rt")) {
            j_lo <- findInterval(rt_lo, d$rt, left.open = TRUE) + 1L; j_hi <- findInterval(rt_hi, d$rt)
            if (j_hi >= j_lo) d[j_lo:j_hi] else d[0]
        } else d[rt >= rt_lo & rt <= rt_hi]
        sub <- sub[mz >= mz_lo & mz <= mz_hi]
        if (!nrow(sub)) return(NULL)
        scan_rts <- rts[lo_scan:hi_scan]
        m <- matrix(0, nrow = length(scan_rts), ncol = length(mz_axis))
        m[cbind(match(sub$rt, scan_rts), match(sub$mz, mz_axis))] <- as.double(sub$abundance)
        win <- (i0:i1) - lo_scan + 1L
        apex <- win[which.max(rowSums(m[win, , drop = FALSE]))]
        top <- max(win[1], apex - apex_half):min(win[length(win)], apex + apex_half)
        sig <- colMeans(m[top, , drop = FALSE])
        bg_rows <- c(seq_len(i0 - lo_scan), if (hi_scan > i1) (i1 - lo_scan + 2L):(hi_scan - lo_scan + 1L))
        bg <- if (length(bg_rows)) colMeans(m[bg_rows, , drop = FALSE]) else numeric(length(mz_axis))
        s <- pmax(sig - bg, 0)
        if (sum(sig) <= 0 || max(s) <= 0) return(NULL)
        list(spectrum = 100 * s / max(s), apex_rt = scan_rts[apex], signal_frac = sum(s) / sum(sig))
    }

    ## Agreement of each peak with the other peaks of its peak_ID in the folder:
    ## cosine (sqrt-weighted, which tames the big low-mass fragments) of its
    ## spectrum to the leave-one-out MEDIAN of the rest. NA when it has no
    ## spectrum, no real ID ("", NA, unknown), or fewer than 2 other peaks with
    ## that ID to compare against. `spectra` is a list of numeric vectors or
    ## NULLs, parallel to `ids`.
    gcms_spectrum_agreement <- function(spectra, ids) {
        n <- length(spectra); out <- rep(NA_real_, n)
        ids <- as.character(ids)
        has <- !vapply(spectra, is.null, logical(1)) & !is.na(ids) & nzchar(trimws(ids)) &
               !(tolower(trimws(ids)) %in% c("unknown", "na"))
        for (id in unique(ids[has])) {
            idx <- which(has & ids == id)
            if (length(idx) < 3L) next
            W <- sqrt(do.call(rbind, spectra[idx]))
            Wn <- W / pmax(sqrt(rowSums(W^2)), 1e-12)
            cons_all <- gcms_loo_col_medians(W)
            for (j in seq_along(idx)) {
                cons <- cons_all[j, ]
                out[idx[j]] <- sum(Wn[j, ] * cons) / max(sqrt(sum(cons^2)), 1e-12)
            }
        }
        out
    }

    ## Row j = the column medians of W with row j left out, identical to
    ## apply(W[-j, ], 2, stats::median). That apply() was ~2/3 of every Shift+Q at
    ## 95 samples (961 medians per peak per ID, 2026-10-21); here each column is
    ## sorted once and row j's median is read off the sorted column, skipping j's
    ## own rank. Even counts average the two middle values through rowMeans()
    ## (long-double sum, as mean() in median.default). test_confidence.R pins it
    ## against the apply() form.
    gcms_loo_col_medians <- function(W) {
        k <- nrow(W); p <- ncol(W); m <- k - 1L
        if (m < 1L) return(matrix(NA_real_, k, p))
        ## One stable radix order over (column, value) sorts every column at once.
        o <- order(col(W), W, method = "radix")
        S <- matrix(W[o], nrow = k)
        R <- matrix(0L, k, p); R[o] <- rep(seq_len(k), p)     # rank of each cell within its column
        cols <- seq_len(p)
        out <- matrix(NA_real_, k, p)
        for (j in seq_len(k)) {
            r <- R[j, ]
            at <- function(pos) { q <- pos + (pos >= r); S[cbind(q, cols)] }
            out[j, ] <- if (m %% 2L == 1L) at((m + 1L) %/% 2L) else
                rowMeans(cbind(at(m %/% 2L), at(m %/% 2L + 1L)))
        }
        out
    }

    ## Marker class for a score: thresholds c(red_below, green_from), chosen on
    ## Lucas's leaf-wax folder (2026-10-21): IDs that are right every time (the
    ## standard, C27-C31 alkanes, most acids) median 0.95-1.0; the known problem
    ## IDs (C29/C30/C32 alcohols, fernenone, isoarborinol) median 0.5-0.7.
    gcms_confidence_class <- function(score, thresholds = c(0.75, 0.90)) {
        cls <- ifelse(is.na(score), "none",
               ifelse(score >= thresholds[2], "high",
               ifelse(score >= thresholds[1], "medium", "low")))
        factor(cls, levels = c("high", "medium", "low", "none"))
    }
    gcms_confidence_colours <- c(high = "#1a9850", medium = "#f0a30a", low = "#d73027", none = "grey60")

    ## One marker per overlaid peak: the label position from gcms_overlay_rows()
    ## (window median x, apex y) joined to the peak's score by sample + peak number.
    ## A table written before phase 3 has no score column: every marker is "none".
    gcms_overview_markers <- function(labels, peak_table, thresholds = c(0.75, 0.90)) {
        if (is.null(labels) || !nrow(labels)) return(data.frame(x = numeric(0), y = numeric(0),
            confidence = gcms_confidence_class(numeric(0), thresholds), path_to_cdf_csv = character(0)))
        score <- if ("spectrum_agreement" %in% names(peak_table)) {
            key_pt <- paste(peak_table$path_to_cdf_csv, peak_table$peak_number_within_sample, sep = "\r")
            key_lb <- paste(as.character(labels$path_to_cdf_csv), labels$peak_number_within_sample, sep = "\r")
            suppressWarnings(as.numeric(peak_table$spectrum_agreement))[match(key_lb, key_pt)]
        } else rep(NA_real_, nrow(labels))
        data.frame(x = labels$x, y = labels$y, confidence = gcms_confidence_class(score, thresholds),
                   path_to_cdf_csv = labels$path_to_cdf_csv, peak_number_within_sample = labels$peak_number_within_sample,
                   score = score, stringsAsFactors = FALSE)
    }

    ## The confidence squares as drawn (Lucas, 2026-10-22). The ribbons own the
    ## plot's fill scale (peak number), so the squares
    ## cannot MAP fill; each confidence class is its own layer with a constant fill,
    ## which needs no scale and adds no legend. Returns a list to add to a ggplot.
    gcms_detail_marker_layers <- function(markers, size = 5) {
        layers <- list()
        for (cls in names(gcms_confidence_colours)) {
            m <- markers[as.character(markers$confidence) == cls, , drop = FALSE]
            if (!nrow(m)) next
            layers[[length(layers) + 1]] <- ggplot2::geom_point(
                data = m, mapping = ggplot2::aes(x = x, y = y), inherit.aes = FALSE,
                shape = 22, size = size, stroke = 0.25, colour = "grey20",
                fill = gcms_confidence_colours[[cls]], alpha = 0.85)
        }
        layers
    }

    ## The key to those squares, shown in the button row (Lucas, 2026-10-26).
    gcms_confidence_legend <- function(thresholds = c(0.75, 0.90)) {
        ## R colour names such as "grey60" are not CSS colours: always hand CSS a hex.
        css_col <- function(x) grDevices::rgb(t(grDevices::col2rgb(x)), maxColorValue = 255)
        sw <- function(cls) shiny::tags$span(style = paste0(
            "display:inline-block; width:11px; height:11px; border:1px solid #333; margin:0 3px 0 9px;",
            " vertical-align:middle; background:", css_col(gcms_confidence_colours[[cls]]), ";"))
        lo <- format(thresholds[1]); hi <- format(thresholds[2])
        shiny::tags$span(style = "margin-left:15px; font-size:12px; white-space:nowrap;",
            title = "Each peak's square: how well its mass spectrum matches the other peaks with the same peak_ID in this folder",
            "MS agreement with same-named peaks:",
            sw("high"), paste0("\u2265 ", hi),
            sw("medium"), paste0(lo, "\u2013", hi),
            sw("low"), paste0("< ", lo),
            sw("none"), "no comparison")
    }

    ## The TIC view's ion key, moved off the plot into the button row beside the
    ## agreement key (Lucas, 2026-10-28). `pal` is the named ion -> colour map the
    ## plot's colour scale is built from, so the two cannot disagree. NULL hides it.
    gcms_ion_legend <- function(pal) {
        if (is.null(pal) || !length(pal)) return(NULL)
        css_col <- function(x) grDevices::rgb(t(grDevices::col2rgb(x)), maxColorValue = 255)
        keys <- lapply(names(pal), function(ion) shiny::tagList(
            shiny::tags$span(style = paste0(
                "display:inline-block; width:16px; height:3px; margin:0 3px 0 9px;",
                " vertical-align:middle; background:", css_col(pal[[ion]]), ";")),
            ion))
        shiny::tags$span(style = "margin-left:15px; font-size:12px; white-space:nowrap;",
            "ion:", keys)
    }

    ## "Zoom to peaks" (Lucas, 2026-10-26). Accepts "12", "12-15", "12 to 15" or
    ## "15-12"; returns c(lo, hi) as integers, or NULL when it cannot read one.
    gcms_parse_peak_range <- function(txt) {
        if (is.null(txt) || length(txt) != 1 || is.na(txt)) return(NULL)
        n <- suppressWarnings(as.numeric(regmatches(txt, gregexpr("[0-9]+(\\.[0-9]+)?", txt))[[1]]))
        n <- n[is.finite(n)]
        if (length(n) < 1 || length(n) > 2) return(NULL)
        n <- round(n)
        as.integer(c(min(n), max(n)))
    }
    ## The aligned-rt window holding peaks lo..hi: from the earliest start to the
    ## latest end of each peak's median window (gcms_window_for_peak()), padded by
    ## `pad` of its width each side. Numbers no sample has are skipped; NULL if none.
    gcms_peak_range_window <- function(peak_table, lo, hi, pad = 0.1) {
        w <- lapply(seq(lo, hi), function(n) gcms_window_for_peak(peak_table, n))
        w <- do.call(rbind, w[!vapply(w, is.null, TRUE)])
        if (is.null(w) || !nrow(w)) return(NULL)
        r <- c(min(w[, 1]), max(w[, 2]))
        r + c(-1, 1) * pad * max(diff(r), 1e-6)
    }

    ## Right-click > Remove peak (Lucas, 2026-10-26). The row of `peak_table` whose
    ## bounds (peak_start..peak_end, RAW rt, inclusive) hold `rt` in `sample`, or NA.
    ## Where peaks nest, the narrowest wins: it is the one the cursor singles out.
    gcms_peak_at <- function(peak_table, sample, rt) {
        if (is.null(peak_table) || !nrow(peak_table) || is.null(sample) || length(sample) != 1 ||
            length(rt) != 1 || !is.finite(rt)) return(NA_integer_)
        hit <- which(as.character(peak_table$path_to_cdf_csv) == sample &
                     peak_table$peak_start <= rt & peak_table$peak_end >= rt)
        if (!length(hit)) return(NA_integer_)
        hit[which.min(peak_table$peak_end[hit] - peak_table$peak_start[hit])]
    }

#### Stateless helpers lifted out of analyzeGCMSdata7() (2026-10-22)
####
#### These read nothing from the app's own state, so they live at the top
#### level where tests can load them directly. gcms_ prefix as above, so
#### sourcing the file cannot clobber (or be clobbered by) a user's objects.

    ## Detect CDF flavour: "ms" for ANDI-MS (full rt x mz cube),
    ## "chrom" for ANDI-Chrom single-detector (e.g. GC-TCD, GC-FID, GC-ECD),
    ## "unknown" otherwise. Used to fork the ingest path below.
    gcms_detect_cdf_type <- function(path) {
        nc <- ncdf4::nc_open(path)
        on.exit(ncdf4::nc_close(nc))
        vars <- names(nc$var)
        if (all(c("mass_values", "intensity_values", "scan_index") %in% vars)) {
            return("ms")
        }
        if ("ordinate_values" %in% vars) {
            return("chrom")
        }
        "unknown"
    }

    ## Read an ANDI-Chrom CDF and return a data frame matching the
    ## framedDataFile shape used downstream (mz, intensity, rt). A
    ## single-detector trace has no mz axis, so mz is a sentinel 0:
    ## the TIC pipeline sums by rt to recover the one signal channel,
    ## and any extracted-ion request returns empty (no spectrum exists).
    gcms_read_chrom_cdf_as_framed <- function(path) {
        nc <- ncdf4::nc_open(path)
        on.exit(ncdf4::nc_close(nc))
        ord <- as.numeric(ncdf4::ncvar_get(nc, "ordinate_values"))
        var_names <- names(nc$var)
        delay <- if ("actual_delay_time" %in% var_names) {
            as.numeric(ncdf4::ncvar_get(nc, "actual_delay_time"))
        } else 0
        interval <- if ("actual_sampling_interval" %in% var_names) {
            as.numeric(ncdf4::ncvar_get(nc, "actual_sampling_interval"))
        } else 1
        n <- length(ord)
        rt <- delay + (seq_len(n) - 1) * interval
        data.frame(mz = 0, intensity = ord, rt = rt)
    }

    ## ---- cache stamping -------------------------------------------------
    ## Every derived file in this folder used to be keyed on its source's
    ## NAME alone, which made replacing an acquisition invisible: drop a new
    ## run in under the same filename and the app kept plotting and
    ## integrating the OLD one, because <file>.CDF.csv already existed so it
    ## was never rebuilt and chromatograms.csv still held its rows. Deleting
    ## the .CDF.csv by hand did not help either, since chromatograms.csv was
    ## the copy actually being read. Each derived file now carries a .stamp
    ## sidecar holding its source's mtime and size.
    gcms_file_stamp <- function(path) {
        fi <- file.info(path)
        if (is.na(fi$size)) return(NA_character_)
        paste0(format(as.numeric(fi$mtime), digits = 17), ":", fi$size)
    }

    gcms_stamp_path <- function(derived) paste0(derived, ".stamp")

    gcms_stamp_is_current <- function(derived, source, tag = "") {
        if (!file.exists(derived)) return(FALSE)
        want <- paste0(gcms_file_stamp(source), tag)
        if (is.na(gcms_file_stamp(source))) return(TRUE)     # source gone: keep what we have
        sp <- gcms_stamp_path(derived)
        if (!file.exists(sp)) return(FALSE)
        identical(readLines(sp, warn = FALSE)[1], want)
    }

    gcms_write_stamp <- function(derived, source, tag = "") {
        want <- gcms_file_stamp(source)
        if (!is.na(want)) writeLines(paste0(want, tag), gcms_stamp_path(derived))
        invisible(NULL)
    }

    ## Write through a temp file in the same directory, then rename. An
    ## interrupted fwrite straight to the destination left a TRUNCATED file
    ## behind, which file.exists() then accepted for ever.
    gcms_fwrite_atomic <- function(x, path) {
        tmp <- paste0(path, ".tmp", Sys.getpid())
        data.table::fwrite(x, tmp)
        if (!file.rename(tmp, path)) {
            file.copy(tmp, path, overwrite = TRUE)
            unlink(tmp)
        }
        invisible(path)
    }

    gcms_allion_store_path <- function(cdf_csv) paste0(cdf_csv, ".allions.csv")

    ## Nominal m/z for the all-ion store. NOT round(): R rounds exact halves
    ## to the EVEN integer, and this instrument reports heavy ions at x.5 -
    ## Skylar's ion 467 sits at m/z 467.5 and went into bin 468, along with
    ## 468.5 (2026-10-16), while the extracted-ion traces (a +/-0.6 window)
    ## showed it at 467. floor(mz + 0.3) bins 281.7-282.69 as 282: whole-
    ## number data unchanged, mass defects from -0.3 to +0.69 (wax-range
    ## hydrocarbons up to ~m/z 600) kept on their nominal mass.
    ## test_nominal_binning.R pins this.
    gcms_nominal_mz <- function(mz) floor(as.numeric(mz) + 0.3)

    ## y transform for the overlay, plus the per-line opacity that goes
    ## with it. "norm" divides each ion by its own maximum in the window,
    ## which is what actually makes co-elution obvious: every ion becomes
    ## a unit-height peak, so a two-scan apex shift is unmissable. It
    ## also destroys quantitation and would draw a noise-level ion as
    ## tall as the base peak -- hence opacity carrying absolute
    ## abundance in that mode, and only in that mode.
    gcms_allion_scale <- function(w, mode) {
        w <- data.table::copy(w)
        if (identical(mode, "norm")) {
            w[, peak := max(net), by = mz]
            w[, y := ifelse(peak > 0, net / peak, 0)]
            top <- max(w$peak)
            w[, opacity := 0.10 + 0.65 * (peak / top)^0.35]
            w[, peak := NULL]
        } else {
            w[, y := switch(mode,
                raw  = net,
                sqrt = sqrt(net),
                log  = log10(pmax(net, 1)),
                net)]
            w[, opacity := 0.35]
        }
        w[]
    }

    gcms_allion_y_label <- function(mode) switch(mode,
        raw  = "Abundance (counts, baseline removed)",
        sqrt = "sqrt abundance (baseline removed)",
        log  = "log10 abundance (baseline removed)",
        norm = "Abundance, each ion scaled to its own max",
        "Abundance")

    ## Read a panel_params range whichever way this ggplot2 stores it.
    gcms_pp_range <- function(pp, axis) {
        v <- pp[[paste0(axis, ".range")]]
        if (is.null(v)) v <- tryCatch(pp[[axis]]$continuous_range, error = function(e) NULL)
        if (is.null(v)) v <- tryCatch(pp[[axis]]$dimension(), error = function(e) NULL)
        if (length(v) == 2 && all(is.finite(v)) && v[2] > v[1]) v else NULL
    }

    ## Turn a fallback brush coordinate into a true fraction of the rendered
    ## IMAGE, using the brush's OWN reported domain.
    ##
    ## Shiny's degenerate coordmap does NOT report 0..1, which is what this
    ## code assumed for three weeks. The labelled probe on host1
    ## (2026-10-02, 5 facets, 1124x700 css at a 2x pixel ratio) came back:
    ##
    ##   domain  left=-0.04 right=1.04 bottom=-0.04 top=1.04
    ##   range   left=0 right=2248 bottom=1399 top=-1
    ##   coords_img xmin=230 xmax=1222 ymin=296 ymax=364
    ##   xmin=0.0704982 xmax=0.547082 ymin=0.758429 ymax=0.810886
    ##
    ## i.e. a UNIT domain carrying ggplot's own 5% expansion, mapped across
    ## the WHOLE 2248x1400 image. Inverting Shiny's transform reproduces the
    ## raw pixels to five decimals on both axes: 230/2248 = 0.10231 =
    ## (0.0704982 + 0.04)/1.08. So reading xmin straight as a fraction is a
    ## -0.04 shift plus a 1.08 stretch -- about 4% of the window, which is
    ## why every reading was "nearly right" and why a drag inside a panel
    ## could report a y fraction ABOVE the panel top (0.8109 vs the panel's
    ## 0.8016) and get clamped, losing the top of the y-zoom.
    ##
    ## y needs no flipping: Shiny has already done it (range top=-1 < bottom
    ## =1399), so the larger value is the higher point on the picture. That
    ## is the "y fractions are measured from the BOTTOM" rule, confirmed.
    ##
    ## Returns NULL -- leaving the caller's old reading untouched -- unless
    ## the domain really is a unit domain. A surviving coordmap reports the
    ## domain in DATA units (365..4664 here), and those coordinates are
    ## already data values that must never be rescaled.
    gcms_brush_unit_fraction <- function(brush, axis = c("x", "y")) {
        axis <- match.arg(axis)
        d <- brush$domain
        if (is.null(d)) return(NULL)
        if (axis == "x") {
            lo <- d$left;   hi <- d$right; v <- c(brush$xmin, brush$xmax)
        } else {
            lo <- d$bottom; hi <- d$top;   v <- c(brush$ymin, brush$ymax)
        }
        if (is.null(lo) || is.null(hi) || length(lo) != 1 || length(hi) != 1) return(NULL)
        lo <- suppressWarnings(as.numeric(lo)); hi <- suppressWarnings(as.numeric(hi))
        if (!is.finite(lo) || !is.finite(hi)) return(NULL)
        span <- hi - lo
        ## A unit domain: span ~1 (1.08 with the expansion, 1.0 without) and
        ## an origin near 0. Anything else is a real data axis.
        if (!is.finite(span) || span < 0.9 || span > 1.2 || abs(lo) > 0.2) return(NULL)
        if (length(v) != 2 || !all(is.finite(v))) return(NULL)
        (v - lo) / span
    }

    ## Which panel (1 = top) a from-the-bottom image fraction falls in.
    ## Returns NULL when the geometry is unknown; a fraction in a gutter or
    ## just off the stack snaps to the nearest panel, because a drag that
    ## overshoots the panel edge is the normal gesture, not an error.
    gcms_panel_at_bottom_fraction <- function(pf, frac_from_bottom) {
        if (is.null(pf) || is.null(pf$panels) || nrow(pf$panels) == 0) return(NULL)
        if (!is.finite(frac_from_bottom)) return(NULL)
        lo <- 1 - pf$panels$y_hi          # panels carry fractions from the TOP
        hi <- 1 - pf$panels$y_lo
        hit <- which(frac_from_bottom >= lo & frac_from_bottom <= hi)
        if (length(hit) >= 1) {
            return(hit[which.min(abs(frac_from_bottom - (lo[hit] + hi[hit]) / 2))])
        }
        near <- which.min(pmin(abs(frac_from_bottom - lo), abs(frac_from_bottom - hi)))
        if (length(near) == 1) near else NULL
    }

    ## The panels (1 = top) a tall brush covers, for Shift+A / Shift+E across
    ## several samples at once (Lucas, 2026-10-22). A panel counts only when the
    ## brush covers its vertical MIDDLE, so the usual overshoot into a
    ## neighbouring panel does not drag that sample in. `fry` is the brush's
    ## y extent as from-the-bottom fractions (gcms_brush_unit_fraction()).
    ## integer(0) when the geometry is unknown or no middle is covered.
    gcms_panels_covering <- function(pf, fry) {
        if (is.null(pf) || is.null(pf$panels) || nrow(pf$panels) == 0) return(integer(0))
        if (length(fry) != 2 || !all(is.finite(fry))) return(integer(0))
        mid <- 1 - (pf$panels$y_lo + pf$panels$y_hi) / 2     # panels carry fractions from the TOP
        which(mid >= min(fry) & mid <= max(fry))
    }

    ## Split a multi-sample brush selection into one window per sample, each
    ## from that sample's OWN points, so raw rt follows the sample's rt_offset.
    ## Returns a data frame: path_to_cdf_csv, start, end, area (TIC minus
    ## baseline, long format: ion == 0 / ion == "baseline").
    ## Drop every peak lying wholly inside its own sample's window (inclusive).
    gcms_excise_windows <- function(peak_table, win) {
        if (!nrow(peak_table) || !nrow(win)) return(peak_table)
        i <- match(as.character(peak_table$path_to_cdf_csv), win$path_to_cdf_csv)
        gone <- (!is.na(i) & peak_table$peak_start >= win$start[i] & peak_table$peak_end <= win$end[i]) %in% TRUE
        peak_table[!gone, , drop = FALSE]
    }

    gcms_brush_windows_by_sample <- function(pp) {
        if (is.null(pp) || !nrow(pp)) return(data.frame(path_to_cdf_csv = character(0), start = numeric(0),
                                                         end = numeric(0), area = numeric(0)))
        s <- as.character(pp$path_to_cdf_csv)
        do.call(rbind, lapply(unique(s), function(k) {
            g <- pp[s == k, , drop = FALSE]
            data.frame(path_to_cdf_csv = k, start = min(g$rt), end = max(g$rt),
                       area = sum(g$abundance[g$ion == 0]) - sum(g$abundance[g$ion == "baseline"]),
                       stringsAsFactors = FALSE)
        }))
    }

    ## A peak add or a spectrum extraction belongs to ONE sample. When the
    ## facet cannot be recovered the selection spans every panel on the page,
    ## and taking $path_to_cdf_csv[1] silently attributes it to the first
    ## sample in table order rather than the brushed one -- a wrong peak
    ## written against a wrong sample, with an area summed over all twenty.
    ## Refuse and say so; the diagnostics above name what was resolved.
    ## (Shift+G is exempt by design: it writes the same window to every
    ## sample on purpose.)
    gcms_single_sample <- function(pp, what) {
        if (is.null(pp) || nrow(pp) == 0) return(NULL)
        s <- unique(as.character(pp$path_to_cdf_csv))
        if (length(s) == 1) return(s)
        cat(paste0(what, " needs a brush inside ONE sample's panel, but this selection covers ",
                   length(s), " samples (", paste(utils::head(s, 3), collapse = ", "),
                   if (length(s) > 3) ", ..." else "", "). Nothing written.\n"))
        NULL
    }

    ## ONE derivation for Sample_ID, used by both the create and the
    ## append path below. They used to disagree: the first write stored
    ## the full `foo.CDF.csv`, while later appends ran
    ## gsub("\\..*$", "", ...) -- truncating at the FIRST dot, so
    ## `WT.rep2.CDF` and `WT.rep3.CDF` both collapsed to `WT`. Strip only
    ## the .CDF.csv suffix (and any directory), as the facet labels do.
    gcms_sample_id_from_path <- function(x) {
        gsub("\\.CDF\\.csv$", "", gsub(".*/", "", x), ignore.case = TRUE)
    }

#### analyzeGCMSdata7
####
#### v7 = v6 plus the co-elution split above: the "Co-elution split" sidebar
#### fields, applied on Shift+Q by gcms_apply_coelution_split(). With the ion
#### fields empty, v7 behaves exactly as v6.
####
#### --- v6 notes, inherited ---
####
#### v6 = v5 plus the all-ion chromatogram view.
####
#### The TIC hides co-eluting compounds: two compounds under one hump sum to a
#### single smooth peak, and nothing in the trace says there are two. Drawing
#### every ion as its own thin, semi-transparent line exposes them, because the
#### ions of the two compounds apex at slightly different retention times and
#### carry different peak shapes. A pure peak draws as a tight bundle of lines;
#### an impure one splays. Set all_ion_view = TRUE, or switch "Display" in the
#### sidebar.
####
#### The same ion-selection machinery also feeds integration: `area_peak_ions`
#### sums only the ions that rise across a peak's own window, so a peak is no
#### longer charged for the column bleed and tailing neighbours sitting under
#### it. It is reported ALONGSIDE the v5 TIC `area`, never instead of it, so
#### existing numbers stay comparable.

    analyzeGCMSdata7 <- function(
        CDF_directory_path = getwd(),
        zoom_and_scroll_rate = 100,
        baseline_window = 400,
        samples_per_page = 20,
        compact = FALSE,
        x_axis_start_default = NULL,
        x_axis_end_default = NULL,
        path_to_reference_library = busta_spectral_library,
        samples_monolist_subset = NULL,
        ions = 0,
        all_ion_view = FALSE,
        all_ion_threshold = 0.02,
        all_ion_scaling = "sqrt",
        all_ion_max_points = 1200,
        integrate_peak_ions = TRUE,
        score_peak_confidence = TRUE,
        confidence_thresholds = c(0.75, 0.90),
        coelution_crosstalk = c(0, 0),
        facet_label_size = 12,
        facet_label_wrap = 14,
        n_cores = NULL,
        ms_classifier_path = "/project_data/shared/mass_spectral_library/ms_classifier.py",
        jupyter = FALSE
    ) {

        ## The TIC (ion 0) is always extracted: the baseline, integration and
        ## the default display all run on it. A call such as ions = c(523, 537)
        ## used to leave every sample with no TIC ("has 0 scans").
        ions <- unique(c(0, as.numeric(ions)))

        ## Overview (compact) mode was removed on 2026-10-22: Detail is fast enough to
        ## page through a whole folder. The argument is still accepted so an old call
        ## keeps working; it is ignored.
        if (isTRUE(compact)) message("analyzeGCMSdata7(): compact / Overview mode has been removed; drawing the normal panels.")
        ## Likewise zoom_and_scroll_rate: it set the step of the Shift+F/D/V/C pan
        ## keys, removed on 2026-10-22 (the Earlier / Later buttons pan instead).
        if (!missing(zoom_and_scroll_rate)) message("analyzeGCMSdata7(): zoom_and_scroll_rate is no longer used; pan with the Earlier / Later buttons.")




        ## ---- All-ion chromatogram support (new in v6) ----------------------
        ##
        ## Three things make an every-ion overlay affordable at 20 samples per
        ## page, which measurement on a 8434-scan run says it otherwise is not
        ## (~10M line segments, ~25 s per ggplot redraw):
        ##
        ##   * m/z is binned to NOMINAL mass. The extraction path above rounds
        ##     to 0.1, which on unit-resolution quadrupole data yields roughly
        ##     five times more channels than there are real ones.
        ##   * only ions that RISE inside the visible window are drawn. Column
        ##     bleed and persistent background are flat there, so they cost ink
        ##     and CPU without carrying information.
        ##   * retention time is decimated onto a pixel grid, keeping the
        ##     MAXIMUM in each bin so apexes survive. Nothing is gained by
        ##     drawing more x-points than the panel has pixels.
        ##
        ## Together those take the same render to ~1.5M segments and ~5 s, and a
        ## zoomed window to ~1 s.

            .allion_cache <- new.env(parent = emptyenv())
            .allion_cache_order <- character(0)   # least-recently-used first
            .allion_cache_max <- 40L




            ## Bumped whenever the binning changes, so every existing store rebuilds
            ## once instead of serving bins made the old way.
            .allion_binning_tag <- ":nominal-floor+0.3"

            ## One row per (rt, nominal m/z) with non-zero signal. Built once per
            ## sample from the framed CSV and cached on disk beside it. This is
            ## deliberately NOT folded into chromatograms.csv: that file is
            ## shared by every sample, and appending the full cube to it would
            ## turn a small index into hundreds of megabytes.
            build_allion_store <- function(cdf_csv, force = FALSE) {
                out <- gcms_allion_store_path(cdf_csv)
                if (!force && gcms_stamp_is_current(out, cdf_csv, .allion_binning_tag)) return(invisible(out))
                framed <- data.table::as.data.table(data.table::fread(cdf_csv))
                framed[, mz := gcms_nominal_mz(mz)]
                framed <- framed[!is.na(mz) & !is.na(rt) & !is.na(intensity)]
                agg <- framed[, .(abundance = sum(intensity)), by = .(rt, mz)]
                ## Keep every single-detector scan. Dropping non-positive rows is the
                ## right sparsity win for MS counts (an empty m/z cell carries nothing),
                ## but TCD/FID/ECD ordinate values are SIGNED and routinely sit below
                ## zero between peaks -- pruning them punched holes in the store, so
                ## geom_line drew chords across the gaps and peak_ion_area summed only
                ## the positive part of each peak.
                agg <- agg[mz == 0 | abundance > 0]
                data.table::setorder(agg, rt, mz)
                gcms_fwrite_atomic(agg, out)
                gcms_write_stamp(out, cdf_csv, .allion_binning_tag)
                invisible(out)
            }

            ## In-memory cache: validity is checked against the store's stamp, so a
            ## rebuilt store is not shadowed by the copy already in RAM, and the
            ## cache is bounded rather than growing one m/z cube per sample opened
            ## for the life of the session.
            allion_cache_touch <- function(cdf_csv) {
                .allion_cache_order <<- c(setdiff(.allion_cache_order, cdf_csv), cdf_csv)
                while (length(.allion_cache_order) > .allion_cache_max) {
                    drop <- .allion_cache_order[1]
                    .allion_cache_order <<- .allion_cache_order[-1]
                    if (!is.null(.allion_cache[[drop]])) rm(list = drop, envir = .allion_cache)
                }
                invisible(NULL)
            }

            load_allions <- function(cdf_csv) {
                p <- gcms_allion_store_path(cdf_csv)
                fresh <- gcms_stamp_is_current(p, cdf_csv, .allion_binning_tag)
                if (!fresh) {
                    if (!is.null(.allion_cache[[cdf_csv]])) rm(list = cdf_csv, envir = .allion_cache)
                    build_allion_store(cdf_csv, force = TRUE)
                }
                if (!is.null(.allion_cache[[cdf_csv]])) {
                    allion_cache_touch(cdf_csv)
                    return(.allion_cache[[cdf_csv]])
                }
                d <- data.table::as.data.table(data.table::fread(p))
                data.table::setkey(d, rt)
                assign(cdf_csv, d, envir = .allion_cache)
                allion_cache_touch(cdf_csv)
                d
            }

            ## Ions worth drawing in [xs, xe], with a per-ion baseline removed.
            ## The baseline is that ion's 10th percentile INSIDE the window
            ## rather than a median of the flanks either side: it costs one pass
            ## instead of three, and it degrades gracefully as the window
            ## narrows. `threshold` is a fraction of the largest baseline-
            ## corrected ion in the window, so it means the same thing whether
            ## the window holds the base peak or a trace one.
            ## Returns NULL when nothing clears the threshold.
            ## `rt_shift` is this sample's rt_offset. The decimation grid is built in
            ## OFFSET-CORRECTED space and then shifted back, so every sample on the page
            ## snaps to the SAME centres. Without it the "one x grid" rule held only
            ## within a sample: the union across samples was as dense as the raw data,
            ## which made the ion map's tile width (computed from that union) several
            ## times too narrow and drew the map as thin stripes with white gaps.
            allion_window <- function(cdf_csv, xs, xe, threshold, max_points, rt_shift = 0) {
                d <- load_allions(cdf_csv)
                w <- d[rt >= xs & rt <= xe]
                if (nrow(w) == 0) return(NULL)
                w[, base := stats::quantile(abundance, 0.1), by = mz]
                w[, net := pmax(abundance - base, 0)]
                top <- w[, .(peak = max(net)), by = mz]
                if (nrow(top) == 0 || max(top$peak) <= 0) return(NULL)
                keep <- top[peak > threshold * max(top$peak), mz]
                w <- w[mz %in% keep]
                if (nrow(w) == 0) return(NULL)
                if (is.finite(max_points) && data.table::uniqueN(w$rt) > max_points) {
                    ## Grid spans the REQUESTED window in shifted space, not this
                    ## sample's own min/max, so the centres are identical for every
                    ## sample regardless of where its data starts or what its offset is.
                    br <- seq(xs + rt_shift, xe + rt_shift, length.out = as.integer(max_points) + 1L) - rt_shift
                    if (!all(is.finite(br)) || length(br) < 2L || br[length(br)] <= br[1]) {
                        br <- seq(min(w$rt), max(w$rt), length.out = as.integer(max_points) + 1L)
                    }
                    w[, bin := findInterval(rt, br, rightmost.closed = TRUE)]
                    w <- w[w[, .I[which.max(net)], by = .(mz, bin)]$V1]
                    ## Snap to the bin centre so every ion shares ONE x grid.
                    ## Without this each ion keeps its own within-bin rt, the
                    ## union of rt across ions is as dense as the undecimated
                    ## data (so nothing is saved), and the ion map has no
                    ## regular grid to tile. Apexes move by at most half a bin,
                    ## which at these bin counts is sub-pixel.
                    centres <- (br[-length(br)] + br[-1L]) / 2
                    w[, rt := centres[pmin(pmax(bin, 1L), length(centres))]]
                    w[, bin := NULL]
                }
                w[, base := NULL]
                w[]
            }



            ## Shiny infers the brush's x/y variable names from the FIRST LAYER's
            ## aes, then brushedPoints() looks those names up as columns of the
            ## data frame it is handed -- always chromatograms_updated here.
            ## That holds in the TIC view (x = rt_rt_offset, y = abundance, both
            ## real columns) and breaks in the v6 views, where the first layer
            ## maps y to the transformed `y` (all ions) or to `mz` (ion map).
            ## Neither is a column of chromatograms_updated, so every one of the
            ## ~18 brush-driven actions -- Shift+Q's zoom, Shift+A/G peak add,
            ## Shift+1 MS extraction -- would error the moment the view changed.
            ##
            ## So: delegate to brushedPoints() when the mapping is usable, and
            ## otherwise select on x and the facet alone. The brush's y range in
            ## those modes is in sqrt / normalised / m-z units and means nothing
            ## against raw counts, so ignoring it is not a degradation -- it is
            ## the only correct reading.
            ## Where the plot PANEL sits inside the rendered image, as fractions
            ## of image width.
            ##
            ## Needed because Shiny's fallback coordmap (see brushed_chromatogram
            ## below) reports the brush as a fraction of the whole IMAGE, while
            ## the data only occupies the panel -- the y-axis labels, facet strip
            ## and colour bar all sit outside it. Verified on host1: at 960 px the
            ## panel runs 0.041..0.914, at 1400 px 0.028..0.941, i.e. constant
            ## 39 px and 83 px margins. Ignoring that is a left shift plus a
            ## stretch, which is exactly the "not quite the right region" symptom.
            ##
            ## Absolute gtable widths convert straight to inches; the panel is the
            ## "null" column that soaks up whatever is left. Cached per image size
            ## because building a second grob is not free on a million-row plot.
            .panel_cache <- new.env(parent = emptyenv())


            ## view_key and n_facets are part of the CACHE KEY, not decoration: the
            ## margins depend on the plot, not only on the image size. Measured at
            ## 960x600 with this same function -- TIC-like panel x spans
            ## 0.0585..0.9240, all-ions 0.0460..0.9136, all-ions/normalised
            ## 0.0429..0.9136 (the m/z colourbar, and the y-tick label width, both
            ## move with the view). Keying on size alone meant whichever view was
            ## brushed first at a given window size poisoned every other one.
            ## `window` (the DRAWN x and y limits) is part of the cache key as
            ## well, and has to be: the cached object carries $x_range / $y_range /
            ## $panels$y_data_*, which ARE the window. Keying on size + view + facet
            ## count alone meant the first zoom was read against the right window
            ## and every later one against the FIRST window, because the key had not
            ## changed. Field log, host1 2026-10-02: after zooming to 905..2222 the
            ## next brush still printed "drawn x 150.311-4878.66" -- the expansion of
            ## the full 365..4664 view -- and so jumped to a region the user never
            ## selected. The x-tick label widths move the panel edges too, so this is
            ## not only about the ranges.
            ## `facet_sig` likewise: the cached $panels$facet is how a mapping-less
            ## brush recovers WHICH SAMPLE it was drawn on, so it is page-specific.
            ## With samples_per_page = 20 and 100 samples, every page draws 20
            ## facets at the same size in the same view -- so paging was a cache
            ## HIT and Shift+A would have written the brushed peak against page
            ## one's sample names.
            panel_fraction <- function(p, px_w, px_h, res = 72, view_key = NULL, n_facets = NA,
                                       window = NULL, facet_sig = NULL) {
                win_tag <- if (is.null(window) || !length(window)) "" else
                    paste(signif(as.numeric(window), 10), collapse = ",")
                fct_tag <- if (is.null(facet_sig) || !length(facet_sig)) "" else
                    paste(as.character(facet_sig), collapse = "\u001f")
                key <- paste0(px_w, "x", px_h, "|", if (is.null(view_key)) "" else view_key,
                              "|", n_facets, "|", win_tag, "|", fct_tag)
                if (!is.null(.panel_cache[[key]])) return(.panel_cache[[key]])
                tf <- tempfile(fileext = ".png")
                grDevices::png(tf, width = px_w, height = px_h, res = res)
                out <- tryCatch({
                    ## One build, used for both the gtable and the panel ranges.
                    built <- ggplot2::ggplot_build(p)
                    gt <- ggplot2::ggplot_gtable(built)
                    ## No grid.draw(): the unit conversions below need an open device
                    ## (font metrics), not a drawn plot. Drawing was ~45% of this
                    ## function, paid on the first brush of every new view; the
                    ## geometry is identical without it (checked 2026-10-21, test_batch3.R).
                    ws <- gt$widths
                    n <- length(ws)
                    isnull <- grid::unitType(ws) == "null"
                    abs_in <- numeric(n)
                    if (any(!isnull)) abs_in[!isnull] <- grid::convertWidth(ws[!isnull], "in", valueOnly = TRUE)
                    total_in <- px_w / res
                    nv <- numeric(n)
                    if (any(isnull)) nv[isnull] <- as.numeric(ws[isnull])
                    rem <- max(total_in - sum(abs_in), 0)
                    win <- abs_in + if (sum(nv) > 0) rem * nv / sum(nv) else 0
                    cum <- cumsum(win)
                    ip <- grepl("^panel", gt$layout$name)
                    lcol <- min(gt$layout$l[ip]); rcol <- max(gt$layout$r[ip])
                    xf <- c(if (lcol > 1) cum[lcol - 1] else 0, cum[rcol]) / total_in

                    ## Same again down the heights. Reported as fractions from the
                    ## TOP of the image, because that is the direction image
                    ## pixels run -- and the y brush has to be flipped through it.
                    hs <- gt$heights
                    nh <- length(hs)
                    hnull <- grid::unitType(hs) == "null"
                    habs <- numeric(nh)
                    if (any(!hnull)) habs[!hnull] <- grid::convertHeight(hs[!hnull], "in", valueOnly = TRUE)
                    total_h_in <- px_h / res
                    hv <- numeric(nh)
                    if (any(hnull)) hv[hnull] <- as.numeric(hs[hnull])
                    hrem <- max(total_h_in - sum(habs), 0)
                    hin <- habs + if (sum(hv) > 0) hrem * hv / sum(hv) else 0
                    hcum <- cumsum(hin)
                    trow <- min(gt$layout$t[ip]); brow <- max(gt$layout$b[ip])
                    yf <- c(if (trow > 1) hcum[trow - 1] else 0, hcum[brow]) / total_h_in

                    ## PER-PANEL extents. yf above is the bounding box of ALL
                    ## panels, which is right for x (one panel column, shared
                    ## axis) and wrong for y: with facet_grid(path ~ .) and
                    ## samples_per_page = 20 it spans all twenty stacked panels
                    ## while the brush was drawn inside exactly one of them.
                    ## Measured on a 20-facet page at 960x2200: the stack spans
                    ## 0.0025..0.9849 and each panel is 0.0468 of the image tall,
                    ## so the same in-facet gesture in panel 1 and in panel 20
                    ## mapped to completely different zooms.
                    lay <- gt$layout[ip, , drop = FALSE]
                    lay <- lay[order(lay$t, lay$l), , drop = FALSE]
                    panels <- data.frame(
                        y_lo = vapply(lay$t, function(k) (if (k > 1) hcum[k - 1] else 0) / total_h_in, numeric(1)),
                        y_hi = vapply(lay$b, function(k) hcum[k] / total_h_in, numeric(1)),
                        x_lo = vapply(lay$l, function(k) (if (k > 1) cum[k - 1] else 0) / total_in, numeric(1)),
                        x_hi = vapply(lay$r, function(k) cum[k] / total_in, numeric(1)),
                        stringsAsFactors = FALSE
                    )

                    ## Which facet each panel holds. The gtable's panels, ordered
                    ## top-to-bottom then left-to-right, line up with the build
                    ## layout ordered by ROW then COL -- that is how the gtable is
                    ## assembled. This is what lets a mapping-less brush recover
                    ## the sample it was drawn on; the fallback's own
                    ## brush$mapping$panelvar1 is always NULL, by definition.
                    bl <- built$layout$layout
                    meta <- c("PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y", "COORD", "AXIS_X", "AXIS_Y")
                    facet_cols <- setdiff(names(bl), meta)
                    bl <- bl[order(bl$ROW, bl$COL), , drop = FALSE]
                    if (length(facet_cols) >= 1 && nrow(bl) == nrow(panels)) {
                        panels$facet <- as.character(bl[[facet_cols[1]]])
                        facet_var <- facet_cols[1]
                    } else {
                        panels$facet <- NA_character_
                        facet_var <- NA_character_
                    }

                    ## Each panel's own drawn data range, in the same display order.
                    ## Every v6 view sets explicit y limits (scale limits, or
                    ## coord_cartesian for the ion map), so scales = "free_y" is
                    ## overridden and these all agree -- but reading them per panel
                    ## means the conversion stays correct if that ever changes.
                    pps <- built$layout$panel_params
                    if (nrow(bl) == nrow(panels) && !is.null(bl$PANEL)) {
                        pidx <- as.integer(bl$PANEL)
                        panels$y_data_lo <- vapply(pidx, function(k) { v <- gcms_pp_range(pps[[k]], "y"); if (is.null(v)) NA_real_ else v[1] }, numeric(1))
                        panels$y_data_hi <- vapply(pidx, function(k) { v <- gcms_pp_range(pps[[k]], "y"); if (is.null(v)) NA_real_ else v[2] }, numeric(1))
                    } else {
                        panels$y_data_lo <- NA_real_
                        panels$y_data_hi <- NA_real_
                    }

                    ## The panel does NOT span the scale limits: ggplot expands a
                    ## continuous axis by 5% each side, so limits c(200, 3200) draw
                    ## a panel running 50..3350. Read the truth off the built plot
                    ## rather than re-deriving it from the window, or every fraction
                    ## is out by 5% of the window, worst at the panel edges.
                    pp <- built$layout$panel_params[[1]]
                    list(x = xf, y = yf, panels = panels, facet_var = facet_var,
                         x_range = gcms_pp_range(pp, "x"), y_range = gcms_pp_range(pp, "y"))
                }, error = function(e) NULL)
                grDevices::dev.off(); unlink(tf)
                ok <- function(v) length(v) == 2 && all(is.finite(v)) && v[2] > v[1]
                if (!is.null(out) && ok(out$x) && ok(out$y)) {
                    assign(key, out, envir = .panel_cache)
                    return(out)
                }
                NULL
            }



            ## The sample whose panel a plot CLICK landed in (the Align strip picks it
            ## this way). Same trust order as a brush: Shiny's panelvar1 if the coordmap
            ## survived, else the panel geometry of the picture as drawn. NULL if unknown.
            ## strict = TRUE (right-click > Remove peak): a click in a gutter or margin is
            ## NULL rather than snapped to the nearest panel, so it cannot hit a neighbour.
            clicked_sample <- function(click, strict = FALSE) {
                if (is.null(click)) return(NULL)
                if (!is.null(click$panelvar1)) return(as.character(click$panelvar1))
                if (is.null(last_brush_plot) || length(last_brush_size) != 2 ||
                    !all(is.finite(last_brush_size)) || !all(last_brush_size > 0)) return(NULL)
                pfull <- panel_fraction(last_brush_plot, last_brush_size[1], last_brush_size[2],
                                        view_key = last_brush_key, n_facets = last_brush_facets,
                                        window = last_brush_window, facet_sig = last_brush_facet_sig)
                if (is.null(pfull) || is.null(pfull$panels) || !nrow(pfull$panels)) return(NULL)
                fy <- gcms_brush_unit_fraction(list(ymin = click$y, ymax = click$y, domain = click$domain), "y")
                if (is.null(fy)) return(NULL)
                idx <- gcms_panel_at_bottom_fraction(pfull, fy[1])
                if (is.null(idx)) return(NULL)
                if (isTRUE(strict) && (fy[1] < 1 - pfull$panels$y_hi[idx] || fy[1] > 1 - pfull$panels$y_lo[idx])) return(NULL)
                as.character(pfull$panels$facet[idx])
            }

            ## The ALIGNED rt (the x axis) a plot CLICK landed on (right-click > Remove
            ## peak). A surviving coordmap reports click$x in data units already; the
            ## degenerate one reports an image fraction against a unit domain, turned
            ## into rt here as brushed_chromatogram() turns a brush's x. NULL when that
            ## cannot be done reliably: a removal refuses rather than guesses.
            clicked_rt <- function(click) {
                if (is.null(click) || is.null(click$x)) return(NULL)
                x <- suppressWarnings(as.numeric(click$x))
                if (length(x) != 1 || !is.finite(x)) return(NULL)
                frx <- gcms_brush_unit_fraction(list(xmin = x, xmax = x, domain = click$domain), "x")
                if (!is.null(click$mapping$x)) return(x)
                ## No unit domain: x is in data units if Shiny sent a domain at all
                ## (a coordmap whose mapping was lost); with no domain it is unknown.
                if (is.null(frx)) return(if (is.null(click$domain)) NULL else x)
                if (is.null(last_brush_plot) || length(last_brush_size) != 2 ||
                    !all(is.finite(last_brush_size)) || !all(last_brush_size > 0)) return(NULL)
                pfull <- panel_fraction(last_brush_plot, last_brush_size[1], last_brush_size[2],
                                        view_key = last_brush_key, n_facets = last_brush_facets,
                                        window = last_brush_window, facet_sig = last_brush_facet_sig)
                if (is.null(pfull) || length(pfull$x) != 2 || length(pfull$x_range) != 2) return(NULL)
                fr <- (frx[1] - pfull$x[1]) / (pfull$x[2] - pfull$x[1])
                if (!is.finite(fr) || fr < 0 || fr > 1) return(NULL)
                pfull$x_range[1] + fr * diff(pfull$x_range)
            }

            ## Honour a stored y-zoom only in the view it was taken in; otherwise
            ## fall back to the range the data itself occupies.
            ##
            ## zero_floor pins the BOTTOM of an abundance axis at 0 (Lucas,
            ## 2026-10-02: "y minimum is always zero"). A chromatogram read against
            ## a floating baseline is unreadable -- peak heights stop being
            ## comparable between views and a y-brush that lands low in the panel
            ## used to crop the baseline off the bottom of the picture. Only the
            ## TOP of a y-zoom survives, so a y-brush now means "show me up to
            ## here", which is the only thing anyone wants from it.
            ##
            ## NOT applied to the ion map (y is m/z, where 0 is meaningless and
            ## would waste nine tenths of the panel) nor to log10 scaling (where
            ## the drawn values are logs and 0 is a real interior value).
            y_limits_for <- function(key, fallback, zero_floor = FALSE) {
                lim <- if (!is.null(y_zoom) && length(y_zoom) == 2 &&
                           !is.null(y_zoom_key) && identical(y_zoom_key, key)) {
                    y_zoom
                } else {
                    fallback
                }
                if (isTRUE(zero_floor) && length(lim) == 2 && all(is.finite(lim))) {
                    top <- max(lim)
                    if (top > 0) lim <- c(0, top)
                }
                lim
            }

            ## The brushed y range, in the data units of whatever was last drawn.
            ##
            ## When the coordmap survives, brush$ymin/ymax are already data values.
            ## When it does not they are fractions of the image measured from the
            ## TOP (the same fallback that makes x an image fraction), so they need
            ## flipping through the panel as well as rescaling: the SMALLER
            ## fraction is nearer the top of the picture and therefore the LARGER
            ## data value.
            ##
            ## Returns NULL when there is nothing dependable to convert, in which
            ## case the y axis is simply left alone.
            brush_y_range <- function(brush) {
                if (is.null(brush) || is.null(brush$ymin) || is.null(brush$ymax)) return(NULL)
                if (!is.finite(brush$ymin) || !is.finite(brush$ymax)) return(NULL)

                yv <- if (!is.null(brush$mapping)) brush$mapping$y else NULL
                if (!is.null(yv)) return(sort(c(brush$ymin, brush$ymax)))   # already data units

                why <- NULL
                if (length(last_brush_yrange) != 2 ||
                    !all(is.finite(last_brush_yrange)) || last_brush_yrange[2] <= last_brush_yrange[1]) {
                    why <- "the y range of the last render is not known"
                } else if (is.null(last_brush_plot)) {
                    why <- "no plot recorded from the last render"
                } else if (length(last_brush_size) != 2 ||
                           !all(is.finite(last_brush_size)) || !all(last_brush_size > 0)) {
                    why <- "the rendered plot size is not known"
                } else if (!is.finite(brush$ymin) || !is.finite(brush$ymax) ||
                           brush$ymin < -0.5 || brush$ymax > 1.5) {
                    why <- "the y coordinates are not usable fractions"
                }
                ## The reported values are a UNIT DOMAIN, not 0..1 -- see
                ## gcms_brush_unit_fraction(). Without this the top of a drag inside a
                ## panel reads ABOVE the panel top and clamps, so a y-zoom lost its
                ## ceiling (host1 probe 2026-10-02: 0.8109 against a panel topping
                ## out at 0.8016, which is how "it is not getting the top y" looked).
                fry <- gcms_brush_unit_fraction(brush, "y")
                yfr <- if (!is.null(fry)) fry else c(brush$ymin, brush$ymax)
                if (!is.null(why)) {
                    cat(paste0("  y-zoom skipped (", why, "); raw y ",
                               signif(brush$ymin, 4), " to ", signif(brush$ymax, 4), ".\n"))
                    return(NULL)
                }

                pf <- panel_fraction(last_brush_plot, last_brush_size[1], last_brush_size[2],
                                     view_key = last_brush_key, n_facets = last_brush_facets,
                                     window = last_brush_window, facet_sig = last_brush_facet_sig)
                if (is.null(pf)) {
                    cat("  y-zoom skipped (panel bounds unknown).\n")
                    return(NULL)
                }

                ## Dragging a box around a peak naturally overshoots the panel top
                ## and bottom, and those coordinates come back outside 0..1.
                ## Clamping is the right reading -- the user meant "everything
                ## vertically" -- where rejecting outright silently dropped every
                ## other y-zoom in the field (2026-09-14).
                ##
                ## The y fractions are measured from the BOTTOM of the image, not
                ## the top -- the opposite of the x convention and of what the
                ## pixel direction suggests. panel_fraction() reports the panel
                ## from the top, so flip it to match. Established from the field
                ## on 2026-09-14: selecting 600-400 zoomed to 200-0 and selecting
                ## 200-0 zoomed to 600-400, a clean reflection about the axis
                ## midpoint. A selection near the middle looked almost right under
                ## either reading, which is why it first appeared to "mostly" work.
                ## ...and it has to be rescaled inside the ONE panel the gesture
                ## was made in, not across the whole stack. With a single facet the
                ## two are identical, which is the case the field evidence above
                ## validated; with twenty, the stack is ~21x the panel, so an
                ## in-facet drag in panel 1 and the same drag in panel 20 used to
                ## produce completely different zooms.
                centre <- mean(yfr)
                idx <- if (!is.null(pf$panels) && nrow(pf$panels) > 1) gcms_panel_at_bottom_fraction(pf, centre) else NULL
                if (!is.null(idx)) {
                    panel_lo <- 1 - pf$panels$y_hi[idx]
                    panel_hi <- 1 - pf$panels$y_lo[idx]
                    where <- paste0("panel ", idx, "/", nrow(pf$panels),
                                    if (!is.na(pf$panels$facet[idx])) paste0(" ", pf$panels$facet[idx]) else "")
                } else {
                    panel_lo <- 1 - pf$y[2]        # panel bottom, from the image bottom
                    panel_hi <- 1 - pf$y[1]        # panel top
                    where <- "whole panel stack"
                }
                fr <- pmin(pmax((yfr - panel_lo) / (panel_hi - panel_lo), 0), 1)

                ## The panel spans its DRAWN range, which carries ggplot's 5% axis
                ## expansion -- limits c(0, 1000) draw a panel running -50..1050.
                ## Read it off the built plot; fall back to the limits recorded at
                ## render time only when the build did not report one.
                drawn <- if (!is.null(idx) && !is.na(pf$panels$y_data_lo[idx]) && !is.na(pf$panels$y_data_hi[idx])) {
                    c(pf$panels$y_data_lo[idx], pf$panels$y_data_hi[idx])
                } else if (!is.null(pf$y_range)) {
                    pf$y_range
                } else {
                    last_brush_yrange
                }
                lo <- drawn[1]; hi <- drawn[2]
                vals <- sort(lo + fr * (hi - lo))
                cat(paste0("  y-zoom: raw ", signif(brush$ymin, 4), "-", signif(brush$ymax, 4),
                           if (!is.null(fry)) paste0(" (domain-corrected to ", signif(yfr[1], 4),
                                                     "-", signif(yfr[2], 4), ")") else "",
                           " in ", where,
                           ", panel ", signif(panel_lo, 4), "-", signif(panel_hi, 4), " (from bottom)",
                           ", against drawn y ", signif(lo, 6), "-", signif(hi, 6),
                           " -> ", signif(vals[1], 6), " to ", signif(vals[2], 6), "\n"))
                if (!all(is.finite(vals)) || vals[2] <= vals[1]) return(NULL)
                vals
            }


            ## The mass-spectrum plots have the chromatogram's coordmap problem too,
            ## and a worse symptom: with an empty mapping a bare brushedPoints()
            ## ERRORS ("not able to automatically infer `xvar`") rather than
            ## misreading, which took the whole Shift+1/2/3 block down with it --
            ## including the re-plot the user was trying to get back to. Same
            ## conversion as brushed_chromatogram, one axis, no facets.
            brushed_ms <- function(df, brush, xcol = "mz") {
                empty <- if (is.null(df)) NULL else df[0, , drop = FALSE]
                if (is.null(brush) || is.null(df) || nrow(df) == 0 || !(xcol %in% names(df))) return(empty)

                xv <- if (!is.null(brush$mapping)) brush$mapping$x else NULL
                yv <- if (!is.null(brush$mapping)) brush$mapping$y else NULL
                if (!is.null(xv) && !is.null(yv) && xv %in% names(df) && yv %in% names(df)) {
                    return(tryCatch(shiny::brushedPoints(df, brush),
                                    error = function(e) { cat(paste0("MS brush could not be read (", conditionMessage(e), ").\n")); empty }))
                }
                if (is.null(brush$xmin) || is.null(brush$xmax) ||
                    !is.finite(brush$xmin) || !is.finite(brush$xmax)) return(empty)

                rng <- range(df[[xcol]], na.rm = TRUE)
                bxmin <- brush$xmin; bxmax <- brush$xmax
                frx <- gcms_brush_unit_fraction(brush, "x")
                if ((!is.null(frx) || (bxmin >= 0 && bxmax <= 1)) && rng[2] > 1) {
                    fr <- if (!is.null(frx)) frx else c(bxmin, bxmax)
                    pfm <- NULL
                    if (!is.null(last_ms_plot) && length(last_ms_size) == 2 &&
                        all(is.finite(last_ms_size)) && all(last_ms_size > 0)) {
                        pfm <- panel_fraction(last_ms_plot, last_ms_size[1], last_ms_size[2],
                                              view_key = "massSpectra_1", n_facets = 1,
                                              window = last_ms_window)
                    }
                    if (!is.null(pfm)) {
                        fr <- pmin(pmax((fr - pfm$x[1]) / (pfm$x[2] - pfm$x[1]), 0), 1)
                        lo <- if (!is.null(pfm$x_range)) pfm$x_range[1] else rng[1]
                        hi <- if (!is.null(pfm$x_range)) pfm$x_range[2] else rng[2]
                    } else {
                        cat("  (mass-spectrum panel bounds unknown - using the uncorrected image fraction, expect an offset)\n")
                        lo <- rng[1]; hi <- rng[2]
                    }
                    bxmin <- lo + fr[1] * (hi - lo)
                    bxmax <- lo + fr[2] * (hi - lo)
                    cat(paste0("MS brush read as ", signif(bxmin, 6), " to ", signif(bxmax, 6), " ", xcol, ".\n"))
                }
                if (!is.finite(bxmin) || !is.finite(bxmax) || bxmax <= bxmin) return(empty)
                df[!is.na(df[[xcol]]) & df[[xcol]] >= bxmin & df[[xcol]] <= bxmax, , drop = FALSE]
            }

            ## multi = TRUE (Shift+A / Shift+E only): a brush that covers the
            ## middles of two or more panels selects all of those samples, and the
            ## result carries attr "gcms_multi_panel" = TRUE so the caller knows the
            ## several samples are intended, not an unresolved facet.
            brushed_chromatogram <- function(df, brush, multi = FALSE) {
                empty <- df[0, , drop = FALSE]
                multi_applied <- FALSE
                if (is.null(brush)) return(empty)

                xv <- if (!is.null(brush$mapping)) brush$mapping$x else NULL
                yv <- if (!is.null(brush$mapping)) brush$mapping$y else NULL

                ## Best case: the mapping names real columns, so brushedPoints()
                ## can do its own job (log scales, discrete axes, panels).
                if (!is.null(xv) && !is.null(yv) && xv %in% names(df) && yv %in% names(df)) {
                    return(shiny::brushedPoints(df, brush))
                }

                ## Otherwise fall back to x alone. Two reasons to land here:
                ## the y mapping names a column that does not exist in
                ## chromatograms_updated (the v6 views map y to the transformed
                ## `y`, or to `mz`), or the mapping is missing altogether.
                if (is.null(brush$mapping) || is.null(xv)) {
                    cat(paste0("Brush has no plot mapping. Fields present: ",
                               paste(names(brush), collapse = ", "), "\n"))
                    ## RAW PAYLOAD DUMP. The y convention in this fallback is still
                    ## only inferred -- x fractions were read off the field, y was
                    ## assumed to follow the same rule, and the host1 log of
                    ## 2026-10-02 contradicts it: a drag in the bottom facet came
                    ## back ymin/ymax 0.013-0.030, which is below the bottom panel's
                    ## own extent (0.047) whichever end of the image you measure
                    ## from, so it clamped to a degenerate range and the facet
                    ## always resolved to the bottom panel. coords_css / coords_img
                    ## are in PIXELS with a known origin, so they settle it outright
                    ## -- print them rather than guess again (canon: never ship an
                    ## unverified coordinate convention).
                    dump_xy <- function(tag, v) {
                        if (is.null(v)) return(invisible(NULL))
                        v <- unlist(v)
                        if (!length(v)) return(invisible(NULL))
                        cat(paste0("    ", tag, ": ",
                                   paste(names(v), signif(suppressWarnings(as.numeric(v)), 6),
                                         sep = "=", collapse = " "), "\n"))
                    }
                    cat(paste0("  RAW brush: xmin=", signif(brush$xmin, 6), " xmax=", signif(brush$xmax, 6),
                               " ymin=", signif(brush$ymin, 6), " ymax=", signif(brush$ymax, 6),
                               " direction=", if (is.null(brush$direction)) "NA" else brush$direction, "\n"))
                    dump_xy("coords_css", brush$coords_css)
                    dump_xy("coords_img", brush$coords_img)
                    dump_xy("img_css_ratio", brush$img_css_ratio)
                    dump_xy("domain", brush$domain)
                    dump_xy("range", brush$range)
                    cat(paste0("    rendered css size (w,h): ",
                               paste(signif(last_brush_size, 6), collapse = ", "),
                               "; facets drawn: ", last_brush_facets, "\n"))
                }
                xv_mapped <- if (!is.null(xv) && xv %in% names(df)) xv else NULL
                if (is.null(xv) || !(xv %in% names(df))) xv <- "rt_rt_offset"
                if (!(xv %in% names(df))) return(empty)

                ## xmin/xmax are already in DATA space when Shiny had a coordinate
                ## map to convert with, so they stay usable even when the mapping
                ## itself did not survive. Sanity-check them against the data's own
                ## x range before trusting them: if the plot had no coordinate map
                ## they are pixels, and pixels would silently select the wrong
                ## slice of the run rather than failing loudly.
                if (is.null(brush$xmin) || is.null(brush$xmax) ||
                    !is.finite(brush$xmin) || !is.finite(brush$xmax)) {
                    cat("Brush carries no usable x range - ignoring it.\n")
                    return(empty)
                }
                rng <- range(df[[xv]], na.rm = TRUE)
                bxmin <- brush$xmin; bxmax <- brush$xmax
                pfull <- NULL        # panel geometry, filled in by the fallback below

                ## Field evidence (host1, 2026-09-14): when Shiny's ggplot coordmap
                ## extraction fails it still sends a brush, but with an EMPTY
                ## mapping and xmin/xmax expressed as FRACTIONS OF THE PANEL
                ## (0..1) against a degenerate 0..1 domain. Two successive drags
                ## came back 0.259-0.411 and 0.347-0.498: identical widths, the
                ## second shifted right, exactly as two equal drags at different
                ## positions would look.
                ##
                ## Those are recoverable. The panel spans the current x-axis
                ## window, so a fraction maps straight onto it. Guarded on the
                ## data's own x range sitting clear of 0..1 -- true for retention
                ## times in seconds, and the check means real data coordinates are
                ## never reinterpreted by mistake.
                ## The guard widened with the domain correction: a unit domain runs
                ## -0.04..1.04, so a drag against the left edge legitimately reports
                ## a NEGATIVE xmin and used to fail "bxmin >= 0" and be discarded.
                ## A unit DOMAIN (frx non-NULL) identifies the fallback on its own,
                ## whatever the data range. The old "rng[1] > 1" guard applied to it
                ## as well, and ANDI-Chrom traces (FID/TCD) start at rt ~0.02 s, so
                ## every FID brush was read as raw seconds 0..1 and Shift+Q zoomed
                ## to the first second of the run (ticks 0.25/0.5/0.75) whatever was
                ## brushed. The range check now guards only the domain-less legacy
                ## reading, where a 0..1 value could genuinely be a data coordinate.
                frx <- gcms_brush_unit_fraction(brush, "x")
                if (is.null(xv_mapped) &&
                    (!is.null(frx) || (rng[1] > 1 && bxmin >= 0 && bxmax <= 1))) {
                    win_lo <- if (!is.null(x_axis_start) && length(x_axis_start) == 1 && is.finite(x_axis_start)) x_axis_start else rng[1]
                    win_hi <- if (!is.null(x_axis_end)   && length(x_axis_end)   == 1 && is.finite(x_axis_end))   x_axis_end   else rng[2]
                    if (win_hi > win_lo) {
                        fr <- if (!is.null(frx)) frx else c(brush$xmin, brush$xmax)

                        ## Re-express image fractions as PANEL fractions. Skipped
                        ## (with a warning) if the last plot or its rendered size
                        ## is unknown, in which case the old uncorrected reading
                        ## is used -- offset, but better than nothing.
                        pf <- NULL
                        if (!is.null(last_brush_plot) && length(last_brush_size) == 2 &&
                            all(is.finite(last_brush_size)) && all(last_brush_size > 0)) {
                            pfull <- panel_fraction(last_brush_plot, last_brush_size[1], last_brush_size[2],
                                                    view_key = last_brush_key, n_facets = last_brush_facets,
                                                    window = last_brush_window, facet_sig = last_brush_facet_sig)
                            if (!is.null(pfull)) pf <- pfull$x
                        }
                        if (!is.null(pf)) {
                            fr <- pmin(pmax((fr - pf[1]) / (pf[2] - pf[1]), 0), 1)
                        } else {
                            cat("  (panel bounds unknown - using the uncorrected image fraction, expect an offset)\n")
                        }

                        ## The panel spans its DRAWN x range, not the axis window:
                        ## scale_x_continuous(limits = c(200, 3200)) draws a panel
                        ## running 50..3350, because ggplot expands a continuous
                        ## axis by 5% each side. Fraction 0 is rt 50, and reading it
                        ## as 200 is a 5%-of-window error, worst at the panel edges
                        ## and zero in the middle -- the same signature as the
                        ## y-flip bug, which is why it could hide behind it.
                        if (!is.null(pfull) && !is.null(pfull$x_range)) {
                            win_lo <- pfull$x_range[1]; win_hi <- pfull$x_range[2]
                        }

                        bxmin <- win_lo + fr[1] * (win_hi - win_lo)
                        bxmax <- win_lo + fr[2] * (win_hi - win_lo)
                        cat(paste0("Brush came back as image fractions (", signif(brush$xmin, 4), "-",
                                   signif(brush$xmax, 4), ")",
                                   if (!is.null(frx)) paste0(", domain-corrected to ",
                                                             signif(frx[1], 4), "-", signif(frx[2], 4)) else "",
                                   if (!is.null(pf)) paste0(", panel spans ", signif(pf[1], 4), "-", signif(pf[2], 4)) else "",
                                   ", drawn x ", signif(win_lo, 6), "-", signif(win_hi, 6),
                                   "; reading as ", signif(bxmin, 8), " to ", signif(bxmax, 8), ".\n"))
                    }
                }

                if (bxmax < rng[1] || bxmin > rng[2]) {
                    cat(paste0("Brush x range (", signif(bxmin, 6), " to ", signif(bxmax, 6),
                               ") lies outside the data (", signif(rng[1], 6), " to ", signif(rng[2], 6),
                               "). Ignoring it.\n"))
                    return(empty)
                }

                keep <- !is.na(df[[xv]]) & df[[xv]] >= bxmin & df[[xv]] <= bxmax

                ## WHICH FACET was brushed, in descending order of trust:
                ##   1. brush$mapping$panelvar1 + brush$panelvar1 - the coordmap survived;
                ##   2. brush$panelvar1 on its own - Shiny sometimes sends the panel
                ##      value at top level even when the mapping came back empty;
                ##   3. the panel the brush's y centre falls in, from the per-panel
                ##      extents measured off the built plot.
                ## Route 1 was the only one there before, and it read the variable
                ## NAME out of brush$mapping -- whose emptiness is precisely what puts
                ## us in this fallback. So it was dead code: the selection ran across
                ## every facet on the page, and Shift+A then took $path_to_cdf_csv[1],
                ## the first sample in table order rather than the brushed one.
                pv <- if (!is.null(brush$mapping)) brush$mapping$panelvar1 else NULL
                if (is.null(pv) && !is.null(pfull) && !is.na(pfull$facet_var)) pv <- pfull$facet_var
                panel_value <- if (!is.null(brush$panelvar1)) as.character(brush$panelvar1) else NULL
                panel_source <- if (!is.null(panel_value)) "brush$panelvar1" else NULL
                if (is.null(panel_value) && !is.null(pfull) && !is.null(pfull$panels) &&
                    nrow(pfull$panels) > 1 && !is.null(brush$ymin) && !is.null(brush$ymax) &&
                    is.finite(brush$ymin) && is.finite(brush$ymax) &&
                    brush$ymin > -0.5 && brush$ymax < 1.5) {
                    ## Domain-corrected, like everything else: the uncorrected
                    ## centre sat ~2% too high, which on the host1 probe put a
                    ## gesture made in panel 4 into panel 5 and attributed the
                    ## peak to wildtype7 instead of vector7.
                    fry <- gcms_brush_unit_fraction(brush, "y")
                    ycentre <- if (!is.null(fry)) mean(fry) else mean(c(brush$ymin, brush$ymax))
                    cov <- if (isTRUE(multi) && !is.null(fry)) gcms_panels_covering(pfull, fry) else integer(0)
                    cov_values <- as.character(pfull$panels$facet[cov]); cov_values <- cov_values[!is.na(cov_values)]
                    pidx <- if (length(cov_values) >= 2) NULL else gcms_panel_at_bottom_fraction(pfull, ycentre)
                    if (length(cov_values) >= 2) {
                        panel_value <- cov_values
                        panel_source <- paste0("panels ", paste(cov, collapse = ", "), " of ", nrow(pfull$panels),
                                               ", whose middles the brush covers")
                    } else if (!is.null(pidx) && !is.na(pfull$panels$facet[pidx])) {
                        panel_value <- pfull$panels$facet[pidx]
                        panel_source <- paste0("panel ", pidx, "/", nrow(pfull$panels), " by y position")
                        ## gcms_panel_at_bottom_fraction() SNAPS to the nearest panel when the
                        ## centre lands in a gutter, which is right for a drag that
                        ## overshoots -- and indistinguishable, silently, from a y
                        ## convention that is simply wrong. The host1 log of 2026-10-02
                        ## showed a 5-facet page where every brush reported y 0.01-0.03,
                        ## outside the whole stack, and so always snapped to the bottom
                        ## panel. Say so, instead of reporting a guess as a reading.
                        p_lo <- 1 - pfull$panels$y_hi[pidx]
                        p_hi <- 1 - pfull$panels$y_lo[pidx]
                        if (ycentre < p_lo || ycentre > p_hi) {
                            cat(paste0("  WARNING: the brush's y centre (", signif(ycentre, 4),
                                       ") is OUTSIDE every panel; panel ", pidx,
                                       " (", signif(p_lo, 4), "-", signif(p_hi, 4),
                                       " from the bottom) is only the nearest one. Treat the facet",
                                       " and the y-zoom as unreliable.\n"))
                        }
                    }
                }
                if (!is.null(pv) && pv %in% names(df) && !is.null(panel_value)) {
                    keep <- keep & as.character(df[[pv]]) %in% panel_value
                    ## Marked multi only HERE, where the facet filter was actually applied:
                    ## an unfiltered selection spans every panel and must still be refused.
                    multi_applied <- length(panel_value) > 1
                    cat(paste0("Brush resolved to ", if (length(panel_value) > 1) paste0(length(panel_value), " facets: ") else "facet ",
                               paste(panel_value, collapse = ", "), " (", panel_source, ").\n"))
                } else if (!is.null(pfull) && !is.null(pfull$panels) && nrow(pfull$panels) > 1) {
                    cat("Brush facet could not be resolved - the selection spans every panel on the page.\n")
                }
                out <- df[which(keep), , drop = FALSE]
                if (nrow(out) == 0) {
                    cat(paste0("Brush (", xv, " ", signif(bxmin, 6), " to ", signif(bxmax, 6),
                               if (!is.null(brush$panelvar1)) paste0(", panel ", brush$panelvar1) else "",
                               ") matched no rows.\n"))
                } else {
                    cat(paste0("Brush selected ", nrow(out), " points, ", xv, " ",
                               signif(min(out[[xv]]), 6), " to ", signif(max(out[[xv]]), 6), ".\n"))
                }
                if (isTRUE(multi_applied)) attr(out, "gcms_multi_panel") <- TRUE
                out
            }

            ## Area over the ions a peak actually generates. v5 integrates the
            ## TIC, which charges the peak for whatever co-elutes beneath it.
            ## Restricting the sum to ions that genuinely belong to the peak
            ## removes that contribution. No decimation here -- this is a number,
            ## not a picture, so every scan counts.
            ##
            ## Rewritten 2026-09-30. The first version reused allion_window(),
            ## which was built for DRAWING, and it was wrong here in two ways:
            ##
            ##   * Its per-ion baseline is the 10th percentile INSIDE the
            ##     window. For a correctly-bounded peak the peak fills its own
            ##     window, so that percentile sits well up the flank and the
            ##     peak eats its own baseline. Measured on a clean Gaussian:
            ##     bounds at +/-3 sigma recovered 94% of the true area, +/-2
            ##     sigma 64%, +/-1.5 sigma 39%. The number therefore swung ~2x
            ##     with nothing but how tightly the bounds were drawn.
            ##     => the baseline is now taken from the FLANKS, OUTSIDE
            ##        [start, end], and interpolated across the peak, which is
            ##        what the TIC path already does.
            ##
            ##   * It selected ions by in-window AMPLITUDE. A co-eluting second
            ##     compound's ions are tall inside that window, so they cleared
            ##     the cut and were summed straight back in -- i.e. the column
            ##     did not do the one thing it exists to do.
            ##     => selection is now a SHAPE test: an ion is kept only if it
            ##        apexes with the peak and tracks its profile.
            ##
            ## Returns NA (not 0) when there is nothing dependable to integrate,
            ## so a failure reads as "no number" rather than "no signal".
            peak_ion_area <- function(cdf_csv, start_rt, end_rt, threshold,
                                      flank_frac = 1.0, apex_tol_sigmas = 0.5,
                                      min_cor = 0.8) {

                na_out <- list(area = NA_real_, ions = NA_character_)
                if (!is.finite(start_rt) || !is.finite(end_rt) || end_rt <= start_rt) return(na_out)

                d <- try(load_allions(cdf_csv), silent = TRUE)
                if (inherits(d, "try-error") || is.null(d) || nrow(d) == 0) return(na_out)

                ## Single-detector data (GC-TCD/FID/ECD) carries the mz = 0
                ## sentinel and has no spectral dimension, so an ion-restricted
                ## area is meaningless. Short-circuit, as Shift+1 and Shift+4 do.
                if (all(d$mz == 0)) return(na_out)

                width  <- end_rt - start_rt
                pad_lo <- start_rt - flank_frac * width
                pad_hi <- end_rt   + flank_frac * width
                ## load_allions() keys the store on rt, so the window is one contiguous
                ## slice: find its ends by binary search instead of scanning every row
                ## of the store (~0.5M) for every peak on every Shift+Q. Same rows, same
                ## order. An unkeyed table (test stubs) takes the plain scan.
                padded <- if (identical(data.table::key(d)[1], "rt")) {
                    i_lo <- findInterval(pad_lo, d$rt, left.open = TRUE) + 1L
                    i_hi <- findInterval(pad_hi, d$rt)
                    if (i_hi >= i_lo) d[i_lo:i_hi] else d[0]
                } else d[rt >= pad_lo & rt <= pad_hi]
                if (nrow(padded) == 0) return(na_out)

                ## Abundance and rt as double, ONCE, here: median() of an INTEGER vector is
                ## integer for an odd count and double for an even one, so by-mz groups came
                ## back mixed and data.table stopped Shift+Q ("Column 1 of result for group N
                ## is type 'integer' but expecting type 'double'"). fread reads a column as
                ## integer whenever every value is whole. test_ion_area_integer.R pins this.
                ## It used to be as.double() around each median, which also stopped
                ## data.table using its fast grouped median (GForce): ~30% of the function's
                ## time. Whole numbers are exact in double, so every number is unchanged.
                if (is.integer(padded$abundance)) padded[, abundance := as.double(abundance)]
                if (is.integer(padded$rt))        padded[, rt := as.double(rt)]

                inside <- padded[rt >= start_rt & rt <= end_rt]
                if (nrow(inside) == 0 || data.table::uniqueN(inside$rt) < 3) return(na_out)

                ## Per-ion baseline from the flanks, interpolated across the peak.
                ## Median rather than mean so a neighbouring peak intruding into
                ## a flank does not drag the level up. Falls back to whichever
                ## flank exists, then to the in-window minimum when the peak sits
                ## against the start or end of the run.
                left  <- padded[rt <  start_rt]
                right <- padded[rt >  end_rt]
                lsum <- if (nrow(left))  left[,  .(bl = median(abundance), rtl = median(rt)), by = mz] else NULL
                rsum <- if (nrow(right)) right[, .(br = median(abundance), rtr = median(rt)), by = mz] else NULL

                w <- inside
                if (!is.null(lsum)) w <- merge(w, lsum, by = "mz", all.x = TRUE)
                if (!is.null(rsum)) w <- merge(w, rsum, by = "mz", all.x = TRUE)
                if (is.null(lsum)) { w[, bl := NA_real_]; w[, rtl := NA_real_] }
                if (is.null(rsum)) { w[, br := NA_real_]; w[, rtr := NA_real_] }

                fm <- w[, .(floor_mz = min(abundance)), by = mz]
                w[fm, on = "mz", floor_mz := i.floor_mz]
                w[, base := data.table::fifelse(
                        !is.na(bl) & !is.na(br) & is.finite(rtr - rtl) & (rtr - rtl) != 0,
                        bl + (br - bl) * (rt - rtl) / (rtr - rtl),
                    data.table::fifelse(!is.na(bl), bl,
                    data.table::fifelse(!is.na(br), br, floor_mz)))]
                w[, net := pmax(abundance - base, 0)]

                ## Reference = the BASE ION, not the summed profile. The sum is
                ## a blend when something co-elutes, and its apex sits between
                ## the two compounds, which is exactly the case we are trying to
                ## resolve. The single strongest ion belongs unambiguously to one
                ## compound, so it is the honest model peak.
                tops <- w[, .(top = max(net)), by = mz]
                if (nrow(tops) == 0 || max(tops$top) <= 0) return(na_out)
                base_mz <- tops$mz[which.max(tops$top)]

                ref <- w[mz == base_mz, .(rt, ref = net)]
                data.table::setorder(ref, rt)
                if (nrow(ref) < 3 || max(ref$ref) <= 0) return(na_out)
                apex_rt <- ref$rt[which.max(ref$ref)]

                ## Peak width from the base ion's full width at half maximum.
                ## The apex tolerance has to be in units of the PEAK, not of the
                ## integration window: a window twice as wide does not make two
                ## compounds twice as resolvable. Falls back to a quarter of the
                ## window when the peak is too sparsely sampled to measure.
                half  <- ref$rt[ref$ref >= max(ref$ref) / 2]
                sigma_est <- if (length(half) >= 2) (max(half) - min(half)) / 2.355 else width / 4
                if (!is.finite(sigma_est) || sigma_est <= 0) sigma_est <- width / 4

                w <- merge(w, ref, by = "rt")
                prof <- w[, .(
                        top      = max(net),
                        ion_apex = rt[which.max(net)],
                        r        = suppressWarnings(stats::cor(net, ref))
                    ), by = mz]

                keep <- prof[
                    top > threshold * max(prof$top) &                              # amplitude floor
                    abs(ion_apex - apex_rt) <= apex_tol_sigmas * sigma_est &       # co-apexes with the base ion
                    !is.na(r) & r >= min_cor,                                      # and tracks its profile
                    mz]

                ## A real peak's base ion defines `ref`, so it passes by
                ## construction; if nothing does, the window is not a peak.
                ## Keep the strongest ion rather than silently reporting 0.
                if (length(keep) == 0) keep <- prof$mz[which.max(prof$top)]

                kept <- w[mz %in% keep]
                list(area = sum(kept$net), ions = paste(sort(unique(kept$mz)), collapse = ";"))
            }

            ## Shift+Q recomputes area_peak_ions for EVERY peak, and that was ~75% of
            ## each press (714 peaks, ~40 s on host1, 2026-10-20) although a zoom
            ## changes no peak. The result depends only on the sample's all-ion store,
            ## the window and the threshold, so remember it under exactly those: the
            ## store's stamp (mtime:size) is in the key, so a rebuilt store recomputes.
            ## Pinned by test_shiftq_speed.R.
            .ion_area_memo <- new.env(hash = TRUE, parent = emptyenv())
            .ion_area_memo_max <- 100000L
            .ion_area_keys_used <- character(0)
            .ion_area_memo_dirty <- FALSE
            ## Look in the memo BEFORE loading the store: the store cache holds 40
            ## samples, so on a bigger folder every press re-read the rest from CSV
            ## just to find a memo hit (~8% of a 95-sample press, 2026-10-21). The
            ## freshness check load_allions() would make (store stamp vs source) is
            ## made once per sample per press (.store_fresh, emptied by Shift+Q);
            ## a store that is missing or stale takes the old path, which rebuilds it.
            .store_fresh <- new.env(hash = TRUE, parent = emptyenv())
            store_is_fresh <- function(cdf_csv) {
                f <- .store_fresh[[cdf_csv]]
                if (is.null(f)) {
                    p <- gcms_allion_store_path(cdf_csv)
                    f <- file.exists(p) && gcms_stamp_is_current(p, cdf_csv, .allion_binning_tag)
                    assign(cdf_csv, f, envir = .store_fresh)
                }
                f
            }
            ion_area_key <- function(cdf_csv, start_rt, end_rt, threshold)
                paste(cdf_csv, format(start_rt, digits = 17), format(end_rt, digits = 17),
                      format(threshold, digits = 17),
                      gcms_file_stamp(gcms_allion_store_path(cdf_csv)), sep = "|")
            peak_ion_area_memo <- function(cdf_csv, start_rt, end_rt, threshold) {
                if (store_is_fresh(cdf_csv)) {
                    key <- ion_area_key(cdf_csv, start_rt, end_rt, threshold)
                    hit <- .ion_area_memo[[key]]
                    if (!is.null(hit)) { .ion_area_keys_used <<- c(.ion_area_keys_used, key); return(hit) }
                }
                ok <- !inherits(try(load_allions(cdf_csv), silent = TRUE), "try-error")
                if (!ok) return(peak_ion_area(cdf_csv, start_rt, end_rt, threshold))
                key <- ion_area_key(cdf_csv, start_rt, end_rt, threshold)
                .ion_area_keys_used <<- c(.ion_area_keys_used, key)
                hit <- .ion_area_memo[[key]]
                if (!is.null(hit)) return(hit)
                res <- peak_ion_area(cdf_csv, start_rt, end_rt, threshold)
                if (length(.ion_area_memo) >= .ion_area_memo_max)
                    rm(list = ls(.ion_area_memo, all.names = TRUE), envir = .ion_area_memo)
                assign(key, res, envir = .ion_area_memo)
                .ion_area_memo_dirty <<- TRUE
                res
            }

            ## The memo above is per session, so the FIRST Shift+Q after opening a
            ## folder still recomputed every peak (~27 s on Skylar's folder). It is
            ## saved beside peaks_monolist.csv after each Shift+Q and read back on the
            ## next session's first one. Same keys (so a rebuilt store, a moved window
            ## or a new threshold still recomputes), plus the text of peak_ion_area()
            ## and of the store builder/loader: a cache written by different code is ignored whole,
            ## never trusted. Only the current peaks' entries are written, so the file
            ## does not grow with every window ever tried. Pinned by test_shiftq_speed.R.
            .ion_area_cache_file <- "ion_area_cache.rds"
            .ion_area_cache_loaded <- FALSE
            .ion_area_cache_keys <- character(0)       # what the file on disk holds
            ## The integration code AND what builds the store it reads (binning, loading).
            ion_area_code_tag <- function() paste(c(deparse(peak_ion_area), deparse(load_allions),
                                                    deparse(build_allion_store), deparse(gcms_nominal_mz),
                                                    .allion_binning_tag), collapse = "\n")
            ion_area_cache_load <- function(path = .ion_area_cache_file) {
                if (.ion_area_cache_loaded) return(invisible(0L))
                .ion_area_cache_loaded <<- TRUE
                if (!file.exists(path)) return(invisible(0L))
                cache <- try(readRDS(path), silent = TRUE)
                if (inherits(cache, "try-error") || !is.list(cache) ||
                    !identical(cache$code, ion_area_code_tag()) || !is.list(cache$entries))
                    return(invisible(0L))
                .ion_area_cache_keys <<- names(cache$entries)
                n <- 0L
                for (k in names(cache$entries)) {
                    if (is.null(.ion_area_memo[[k]])) { assign(k, cache$entries[[k]], envir = .ion_area_memo); n <- n + 1L }
                }
                invisible(n)
            }
            ion_area_cache_save <- function(keys, path = .ion_area_cache_file) {
                keys <- unique(keys)
                if (!.ion_area_memo_dirty && file.exists(path) && setequal(keys, .ion_area_cache_keys))
                    return(invisible(FALSE))
                entries <- mget(keys, envir = .ion_area_memo, ifnotfound = list(NULL))
                entries <- entries[!vapply(entries, is.null, logical(1))]
                tmp <- paste0(path, ".tmp", Sys.getpid())
                ok <- !inherits(try(saveRDS(list(code = ion_area_code_tag(), entries = entries), tmp), silent = TRUE), "try-error")
                if (ok) ok <- file.rename(tmp, path)
                if (!ok) unlink(tmp)
                if (ok) { .ion_area_memo_dirty <<- FALSE; .ion_area_cache_keys <<- names(entries) }
                invisible(ok)
            }

            ## Phase 3: each peak's background-subtracted spectrum (gcms_peak_spectrum())
            ## for the confidence squares. Same caching shape as the ion areas:
            ## an in-session memo keyed on sample, window and the all-ion store's stamp,
            ## saved to peak_spectra_cache.rds beside the peak table and ignored whole if
            ## the extraction or store code has changed. Pinned by test_confidence.R.
            .spectrum_memo <- new.env(hash = TRUE, parent = emptyenv())
            .spectrum_keys_used <- character(0)
            .spectrum_memo_dirty <- FALSE
            .spectrum_cache_file <- "peak_spectra_cache.rds"
            .spectrum_cache_loaded <- FALSE
            spectrum_code_tag <- function() paste(c(deparse(gcms_peak_spectrum), deparse(load_allions),
                                                    deparse(build_allion_store), deparse(gcms_nominal_mz),
                                                    .allion_binning_tag), collapse = "\n")
            ## `get_store()` returns list(store, rts) for the sample, or NULL; it is
            ## only called on a miss, so a press whose spectra are all cached never
            ## reads a store (same reason as store_is_fresh() above). A stale or
            ## missing store skips the memo, so the store is reloaded (and rebuilt).
            spectrum_key <- function(cdf_csv, start_rt, end_rt)
                paste(cdf_csv, format(start_rt, digits = 17), format(end_rt, digits = 17),
                      gcms_file_stamp(gcms_allion_store_path(cdf_csv)), sep = "|")
            peak_spectrum_memo <- function(cdf_csv, start_rt, end_rt, get_store) {
                key <- spectrum_key(cdf_csv, start_rt, end_rt)
                .spectrum_keys_used <<- c(.spectrum_keys_used, key)
                if (store_is_fresh(cdf_csv) && exists(key, envir = .spectrum_memo, inherits = FALSE))
                    return(.spectrum_memo[[key]])
                st <- get_store()
                if (is.null(st)) return(NULL)
                res <- gcms_peak_spectrum(st$store, start_rt, end_rt, rts = st$rts)
                assign(key, res, envir = .spectrum_memo)      # NULL is a valid, cached answer
                .spectrum_memo_dirty <<- TRUE
                res
            }
            spectrum_cache_load <- function(path = .spectrum_cache_file) {
                if (.spectrum_cache_loaded) return(invisible(0L))
                .spectrum_cache_loaded <<- TRUE
                if (!file.exists(path)) return(invisible(0L))
                cache <- try(readRDS(path), silent = TRUE)
                if (inherits(cache, "try-error") || !is.list(cache) ||
                    !identical(cache$code, spectrum_code_tag()) || !is.list(cache$entries))
                    return(invisible(0L))
                for (k in names(cache$entries)) assign(k, cache$entries[[k]], envir = .spectrum_memo)
                invisible(length(cache$entries))
            }
            spectrum_cache_save <- function(keys, path = .spectrum_cache_file) {
                keys <- intersect(unique(keys), ls(.spectrum_memo, all.names = TRUE))
                if (!.spectrum_memo_dirty && file.exists(path)) return(invisible(FALSE))
                entries <- mget(keys, envir = .spectrum_memo)
                tmp <- paste0(path, ".tmp", Sys.getpid())
                ok <- !inherits(try(saveRDS(list(code = spectrum_code_tag(), entries = entries), tmp), silent = TRUE), "try-error")
                if (ok) ok <- file.rename(tmp, path)
                if (!ok) unlink(tmp)
                if (ok) .spectrum_memo_dirty <<- FALSE
                invisible(ok)
            }

            ## Cold Shift+Q in parallel (2026-10-21). With no cache (a new folder, a
            ## new threshold, edited code) the per-peak ion areas and spectra were
            ## ~70 s at 95 samples, one peak at a time. Before the serial loop runs,
            ## the peaks it would MISS are computed here, one sample per forked
            ## worker (parallel::mclapply), with the same functions on the same
            ## inputs, and dropped into the memos under the same keys; the serial
            ## loop then finds every one. Only fresh stores are prefetched (a stale
            ## one is rebuilt by the serial path), a failed store returns nothing,
            ## and any error in the workers just leaves the serial path to do it.
            ## Results are identical to the serial run (test_shiftq_parallel.R).
            ## `n_cores`: NULL = min(8, cores - 2); 1 turns it off. host1 is shared,
            ## so it never takes the whole machine. Forking is unavailable on Windows.
            gcms_prefetch_cores <- function() {
                n <- if (is.null(n_cores)) min(8L, max(1L, parallel::detectCores() - 2L)) else as.integer(n_cores)
                if (is.na(n) || .Platform$OS.type == "windows") 1L else max(1L, n)
            }
            gcms_prefetch_misses <- function(peak_table, threshold, want_ion, want_spec, split_col) {
                cores <- gcms_prefetch_cores()
                if (cores < 2L || nrow(peak_table) == 0) return(invisible(0L))
                todo <- list()
                for (smp in unique(as.character(peak_table$path_to_cdf_csv))) {
                    if (!store_is_fresh(smp)) next
                    r <- which(peak_table$path_to_cdf_csv == smp)
                    st <- peak_table$peak_start[r]; en <- peak_table$peak_end[r]
                    ion_need <- if (want_ion) r[vapply(seq_along(r), function(i)
                        is.null(.ion_area_memo[[ion_area_key(smp, st[i], en[i], threshold)]]), logical(1))] else integer(0)
                    spec_need <- if (want_spec) r[vapply(seq_along(r), function(i)
                        (is.na(split_col[r[i]]) || !nzchar(split_col[r[i]])) &&
                        !exists(spectrum_key(smp, st[i], en[i]), envir = .spectrum_memo, inherits = FALSE), logical(1))] else integer(0)
                    if (length(ion_need) || length(spec_need))
                        todo[[smp]] <- list(ion = ion_need, spec = spec_need)
                }
                n_miss <- sum(vapply(todo, function(t) length(t$ion) + length(t$spec), numeric(1)))
                if (n_miss < 40L) return(invisible(0L))       # not worth forking for a handful
                work <- function(smp) {
                    t <- todo[[smp]]
                    d <- try(load_allions(smp), silent = TRUE)
                    if (inherits(d, "try-error") || is.null(d)) return(NULL)
                    ion <- lapply(t$ion, function(i) peak_ion_area(smp, peak_table$peak_start[i], peak_table$peak_end[i], threshold))
                    names(ion) <- vapply(t$ion, function(i) ion_area_key(smp, peak_table$peak_start[i], peak_table$peak_end[i], threshold), "")
                    spec <- list()
                    if (length(t$spec) && nrow(d) > 0 && !all(d$mz == 0)) {
                        rts <- unique(d$rt); if (is.unsorted(rts)) rts <- sort(rts)
                        for (i in t$spec) spec[spectrum_key(smp, peak_table$peak_start[i], peak_table$peak_end[i])] <-
                            list(gcms_peak_spectrum(d, peak_table$peak_start[i], peak_table$peak_end[i], rts = rts))
                    }
                    list(ion = ion, spec = spec)
                }
                t0 <- Sys.time()
                res <- try(parallel::mclapply(names(todo), work, mc.cores = min(cores, length(todo)),
                                              mc.preschedule = FALSE), silent = TRUE)
                if (inherits(res, "try-error")) return(invisible(0L))
                n_ion <- 0L; n_spec <- 0L
                for (r in res) {
                    if (!is.list(r)) next                      # a worker that errored or returned NULL
                    for (k in names(r$ion))  { assign(k, r$ion[[k]],  envir = .ion_area_memo); n_ion  <- n_ion + 1L }
                    for (k in names(r$spec)) { assign(k, r$spec[[k]], envir = .spectrum_memo); n_spec <- n_spec + 1L }
                }
                if (n_ion)  .ion_area_memo_dirty <<- TRUE
                if (n_spec) .spectrum_memo_dirty <<- TRUE
                n <- n_ion + n_spec
                cat(sprintf("Shift+Q: %d ion areas / spectra computed on %d cores in %.1f s.\n", n,
                            min(cores, length(todo)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
                invisible(n)
            }

        setwd(CDF_directory_path)


        ## Chemstation exports detector channels as SIGNAL01.CDF, SIGNAL02.CDF, etc.
        ## Rename to <sample_name>.CDF using the file's sample_name global attribute
        ## so downstream labels are meaningful instead of channel numbers.
        signal_files <- dir()[grep("^SIGNAL.*\\.cdf$", dir(), ignore.case = TRUE)]
        for (sf in signal_files) {
            sample_name <- tryCatch({
                nc <- ncdf4::nc_open(sf)
                att <- ncdf4::ncatt_get(nc, 0, "sample_name")
                ncdf4::nc_close(nc)
                if (isTRUE(att$hasatt) && is.character(att$value)) att$value else ""
            }, error = function(e) "")
            if (!nzchar(sample_name)) next
            safe_name <- gsub("[^A-Za-z0-9._-]+", "_", sample_name)
            ext <- sub("^.*(\\.[Cc][Dd][Ff])$", "\\1", sf)
            target <- paste0(safe_name, ext)
            if (target == sf) next
            if (file.exists(target)) {
                i <- 2L
                while (file.exists(paste0(safe_name, "_", i, ext))) i <- i + 1L
                target <- paste0(safe_name, "_", i, ext)
            }
            file.rename(sf, target)
            if (file.exists(paste0(sf, ".csv"))) {
                file.rename(paste0(sf, ".csv"), paste0(target, ".csv"))
            }
            cat(paste0("Renamed ", sf, " -> ", target, "\n"))
        }

        paths_to_cdfs <- dir()[grep("\\.CDF$", dir(), ignore.case = TRUE)]
        paths_to_cdf_csvs <- paste0(paths_to_cdfs, ".csv")

        ## PREPARE DATA: Check for CDF to CSV conversion, check chromatograms

            if ( length(paths_to_cdfs) == 0 ) {
                stop("The directory specified does not contain any .CDF files.")
            } else {
                if ( !file.exists("chromatograms.csv") ) {
                    chromatograms <- list()
                } else {
                    chromatograms <- readMonolist("chromatograms.csv")
                }
                chromatograms_to_add <- list()

                for (file in seq_along(paths_to_cdfs)) {

                    ## If the cdf.csv doesn't exist for this cdf -- or is stale
                    ## against it -- create it. Stale means the .CDF has a different
                    ## mtime/size from the one the .csv was built from, i.e. the
                    ## acquisition was replaced under the same filename. Before the
                    ## stamp, that case silently kept serving the previous run.

                        cdf_csv_stale <- !gcms_stamp_is_current(paths_to_cdf_csvs[file], paths_to_cdfs[file])

                        if ( cdf_csv_stale ) {

                            if ( file.exists(paths_to_cdf_csvs[file]) ) {
                                cat(paste0(paths_to_cdfs[file], " has changed since its CSV was built - rebuilding it and everything derived from it.\n"))
                            }

                            cdf_type_this <- gcms_detect_cdf_type(paths_to_cdfs[file])

                            if (cdf_type_this == "chrom") {

                                ## ANDI-Chrom single-detector file (e.g. GC-TCD).
                                ## One row per timepoint with mz=0 sentinel.
                                cat(paste("CDF (single-detector / chromatography format) to CSV conversion. Reading data file ", paths_to_cdfs[file], "\n", sep = ""))
                                    framedDataFile <- gcms_read_chrom_cdf_as_framed(paths_to_cdfs[file])

                                cat("   Writing out data file as CSV... \n")
                                    data.table::fwrite(framedDataFile, file = paste(paths_to_cdfs[file], ".csv", sep = ""), col.names = TRUE, row.names = FALSE)

                            } else {

                                cat(paste("CDF to CSV conversion. Reading data file ", paths_to_cdfs[file], "\n", sep = ""))
                                    rawDataFile <- xcms::loadRaw(xcms::xcmsSource(paths_to_cdfs[file]))

                                cat("   Framing data file ... \n")
                                    rt <- rawDataFile$rt
                                    scanindex <- rawDataFile$scanindex

                                    filteredRawDataFile <- list()
                                    for ( i in seq_len(max(0L, length(rt) - 1L)) ) {
                                        filteredRawDataFile[[i]] <- data.frame(
                                            mz = rawDataFile$mz[(scanindex[i]+1):(scanindex[i+1])],
                                            intensity = rawDataFile$intensity[(scanindex[i]+1):(scanindex[i+1])],
                                            rt = rt[i]
                                        )
                                    }
                                    framedDataFile <- do.call(rbind, filteredRawDataFile)
                                    framedDataFile$mz <- round(framedDataFile$mz, digits = 1)
                                    framedDataFile <- drop_na(framedDataFile)

                                cat("   Merging duplicate rows ...\n")
                                    if ( dim(table(duplicated(paste(framedDataFile$mz, framedDataFile$rt, sep = "_")))) > 1 ) {
                                        framedDataFile %>% group_by(mz,rt) %>% summarize(intensity = sum(intensity), .groups = "drop") -> framedDataFile
                                        framedDataFile <- as.data.frame(framedDataFile)
                                    }

                                cat("   Writing out data file as CSV... \n")
                                    data.table::fwrite(framedDataFile, file = paste(paths_to_cdfs[file], ".csv", sep = ""), col.names = TRUE, row.names = FALSE)

                            }

                            ## Stamp the new CSV, then drop everything derived from
                            ## the old one: the all-ion store, its in-memory copy,
                            ## and this sample's rows in chromatograms.csv (which is
                            ## the table actually plotted, so leaving them there is
                            ## exactly how a replaced acquisition stayed invisible).
                            gcms_write_stamp(paths_to_cdf_csvs[file], paths_to_cdfs[file])
                            unlink(c(gcms_allion_store_path(paths_to_cdf_csvs[file]),
                                     gcms_stamp_path(gcms_allion_store_path(paths_to_cdf_csvs[file]))))
                            if (!is.null(.allion_cache[[paths_to_cdf_csvs[file]]])) {
                                rm(list = paths_to_cdf_csvs[file], envir = .allion_cache)
                            }
                            if (is.data.frame(chromatograms) && nrow(chromatograms) > 0) {
                                chromatograms <- chromatograms[chromatograms$path_to_cdf_csv != paths_to_cdf_csvs[file], , drop = FALSE]
                            }
                            if (file.exists("chromatograms.csv")) {
                                stale_rows <- readMonolist("chromatograms.csv")
                                writeMonolist(
                                    monolist = stale_rows[stale_rows$path_to_cdf_csv != paths_to_cdf_csvs[file], , drop = FALSE],
                                    monolist_out_path = "chromatograms.csv"
                                )
                            }

                        }

                    ## Build the nominal-mass all-ion store for this sample.
                    ## Done here rather than lazily at first render so the cost
                    ## lands during ingest, where there is already a progress
                    ## log, instead of freezing the first redraw.

                        if ( isTRUE(all_ion_view) && !gcms_stamp_is_current(gcms_allion_store_path(paths_to_cdf_csvs[file]), paths_to_cdf_csvs[file], .allion_binning_tag) ) {
                            cat(paste("   Building all-ion store for ", paths_to_cdf_csvs[file], "\n", sep = ""))
                            build_allion_store(paths_to_cdf_csvs[file])
                        }

                    ## If any chromatograms (tic and ion) are not present for this csv, extract them

                        ## Which ions this sample already has, read from the IN-MEMORY
                        ## tables rather than by re-reading chromatograms.csv from disk
                        ## once per sample. That read (plus the matching write below) is
                        ## what made folder loading O(N^2) in bytes: at 152 samples the
                        ## whole table was read and rewritten 152 times.
                        ## Both sides are compared as character, because the on-disk
                        ## column mixes "0" with "baseline" while `ions` is numeric.
                        ions_for_this_cdf_csv <- character(0)
                        if (is.data.frame(chromatograms) && nrow(chromatograms) > 0) {
                            ions_for_this_cdf_csv <- c(ions_for_this_cdf_csv,
                                as.character(chromatograms$ion[chromatograms$path_to_cdf_csv == paths_to_cdf_csvs[file]]))
                        }
                        if (is.data.frame(chromatograms_to_add) && nrow(chromatograms_to_add) > 0) {
                            ions_for_this_cdf_csv <- c(ions_for_this_cdf_csv,
                                as.character(chromatograms_to_add$ion[chromatograms_to_add$path_to_cdf_csv == paths_to_cdf_csvs[file]]))
                        }
                        ions_for_this_cdf_csv <- unique(ions_for_this_cdf_csv)

                        if ( length(ions_for_this_cdf_csv) > 0 ) {

                            missing_ions <- as.numeric(as.character(ions[!as.character(ions) %in% ions_for_this_cdf_csv]))
                            missing_ions <- dropNA(missing_ions)

                        } else {

                            missing_ions <- ions

                        }

                        if (length(missing_ions) > 0) {

                            cat(paste("Chromatogram extraction. Reading data file ", paths_to_cdf_csvs[file], "\n", sep = ""))
                                framedDataFile <- as.data.frame(data.table::fread(paths_to_cdf_csvs[file]))

                            cat("   Extracting chromatograms...\n")

                                if (0 %in% missing_ions) {
                                    framedDataFile$row_number <- seq(1,dim(framedDataFile)[1],1)
                                    framedDataFile %>%
                                        group_by(rt) %>% summarize(
                                        abundance = sum(intensity),
                                        ion = 0,
                                        rt_first_row_in_raw = min(row_number),
                                        rt_last_row_in_raw = max(row_number)
                                    ) -> chromatogram
                                    chromatogram <- as.data.frame(chromatogram)
                                    chromatogram$rt <- as.numeric(chromatogram$rt)
                                    chromatogram$path_to_cdf_csv <- paste(paths_to_cdfs[file], ".csv", sep = "")
                                    chromatograms_to_add <- rbind(chromatograms_to_add, chromatogram)
                                }

                                ## Only the MISSING ions: re-extracting every requested ion
                                ## duplicated the rows a cached sample already had.
                                if ( length(missing_ions[missing_ions != 0]) > 0 ) {

                                    numeric_ions <- as.numeric(as.character(missing_ions[missing_ions != 0]))
                                    for ( ion in seq_along(numeric_ions) ){
                                        framedDataFile$row_number <- seq(1,dim(framedDataFile)[1],1)
                                        framedDataFile %>%
                                            group_by(rt) %>%
                                            filter(mz > (numeric_ions[ion] - 0.6)) %>%
                                            filter(mz < (numeric_ions[ion] + 0.6)) -> signal
                                            summarize(signal,
                                                abundance = sum(intensity),
                                                ion = numeric_ions[ion],
                                                rt_first_row_in_raw = if (dim(signal)[1] > 0) { min(row_number) } else { 0 },
                                                rt_last_row_in_raw = if (dim(signal)[1] > 0) { max(row_number) } else { 0 }
                                            ) -> chromatogram
                                        chromatogram <- as.data.frame(chromatogram)
                                        chromatogram$rt <- as.numeric(chromatogram$rt)

                                        if (dim(signal)[1] > 0) {
                                            chromatogram$path_to_cdf_csv <- paste(paths_to_cdfs[file], ".csv", sep = "")
                                            chromatograms_to_add <- rbind(chromatograms_to_add, chromatogram)
                                        }
                                    }
                                }
                        }

                }   ## end for each file

                ## ONE write, after the loop. This used to sit INSIDE it, so the whole
                ## table was rebuilt and rewritten once per sample -- and
                ## unconditionally, so even a fully-cached relaunch with nothing new to
                ## extract rewrote it N times for no reason. Now a cached relaunch
                ## writes nothing at all unless there are stale rows to prune.
                ##
                ## Prune on the way out: rows for CDFs no longer in the folder were
                ## dropped from the in-memory copy but written back regardless, so a
                ## deleted acquisition stayed in chromatograms.csv for ever.
                    added_any <- is.data.frame(chromatograms_to_add) && nrow(chromatograms_to_add) > 0
                    had_any   <- is.data.frame(chromatograms) && nrow(chromatograms) > 0

                    if ( added_any ) {

                        to_write <- if (had_any) rbind(chromatograms, chromatograms_to_add) else chromatograms_to_add
                        to_write <- to_write[to_write$path_to_cdf_csv %in% paths_to_cdf_csvs, , drop = FALSE]
                        writeMonolist(monolist = to_write, monolist_out_path = "chromatograms.csv")

                    } else if ( had_any ) {

                        pruned <- chromatograms[chromatograms$path_to_cdf_csv %in% paths_to_cdf_csvs, , drop = FALSE]
                        if (nrow(pruned) != nrow(chromatograms)) {
                            cat(paste0("Dropping ", nrow(chromatograms) - nrow(pruned),
                                       " chromatogram rows for CDFs no longer in this folder.\n"))
                            writeMonolist(monolist = pruned, monolist_out_path = "chromatograms.csv")
                        }

                    } else if ( !file.exists("chromatograms.csv") ) {

                        writeMonolist(chromatograms_to_add, "chromatograms.csv")

                    }

                print("done")
            }

            ## Read in chromatograms

                chromatograms <- readMonolist("chromatograms.csv")

            ## Set up new samples monolist


                ## If it doesn't exist, create it

                    if ( !file.exists("samples_monolist.csv") ) {

                        samples_monolist <- data.frame(
                            Sample_ID = gcms_sample_id_from_path(unique(chromatograms$path_to_cdf_csv)),
                            rt_offset = 0,
                            baseline_window = baseline_window,
                            path_to_cdf_csv = unique(chromatograms$path_to_cdf_csv)
                        )

                        write.table(
                            x = samples_monolist,
                            file = "samples_monolist.csv",
                            row.names = FALSE,
                            sep = ","
                        )

                # If it exists, check to see if all cdfs in this folder are in it, if not, add them

                    } else {

                        samples_monolist <- readMonolist("samples_monolist.csv")

                        missing_from_samples_monolist <- paths_to_cdf_csvs[!paths_to_cdf_csvs %in% unique(readMonolist("samples_monolist.csv")$path_to_cdf_csv)]

                        if ( length(missing_from_samples_monolist) > 0 ) {

                            samples_monolist_additions <- data.frame(
                                Sample_ID = gcms_sample_id_from_path(missing_from_samples_monolist),
                                rt_offset = 0,
                                baseline_window = baseline_window,
                                path_to_cdf_csv = missing_from_samples_monolist
                            )

                            write.table(
                                x = samples_monolist_additions,
                                file = "samples_monolist.csv",
                                row.names = FALSE,
                                col.names = FALSE,
                                sep = ",",
                                append = TRUE
                            )

                            ## Re-read: the rows above went to the FILE, and the
                            ## in-memory copy was left one acquisition behind.
                            samples_monolist <- readMonolist("samples_monolist.csv")
                        }
                    }

            ## Filter chromatograms so only the CDFs in this folder are included

                chromatograms <- chromatograms[chromatograms$path_to_cdf_csv %in% dir()[grep("\\.CDF\\.csv$", dir(), ignore.case = TRUE)],]

            ## Set up several variables and x_axis limits if not specified in function call

                peak_data <- NULL
                peak_points <- NULL

                ## y-zoom state. y_zoom holds the brushed y range and y_zoom_key
                ## names the view it was taken in, so it is only ever re-applied
                ## to the axis it actually refers to.
                y_zoom <- NULL
                y_zoom_key <- NULL
                last_brush_plot <- NULL
                last_brush_size <- c(NA_real_, NA_real_)
                last_brush_key <- NULL
                last_brush_yrange <- c(NA_real_, NA_real_)
                ## How many facets the last render actually drew. Part of the
                ## panel-geometry cache key, because panel heights move with it.
                last_brush_facets <- NA_integer_
                ## The DRAWN window of the last render (x limits then y limits).
                ## Part of the panel-geometry cache key -- see panel_fraction().
                ## Deliberately NOT folded into last_brush_key, which is the
                ## y-zoom VIEW key and must stay window-independent or every
                ## x-zoom would silently discard the y-zoom with it.
                last_brush_window <- NULL
                ## The facet NAMES the last render drew, in display order -- the
                ## rest of the panel-geometry cache key. See panel_fraction().
                last_brush_facet_sig <- NULL
                ## Same pair for the mass-spectrum plot, used by brushed_ms().
                last_ms_plot <- NULL
                last_ms_size <- c(NA_real_, NA_real_)
                last_ms_window <- NULL

                ## Everything the handlers reach for with `<<-` is declared HERE, so
                ## the assignment binds in this function's frame. Without a local
                ## binding, `<<-` walks all the way out and writes to globalenv() --
                ## and nothing clears globalenv() when the app quits. Open folder A,
                ## quit, open folder B in the same R session and the panel showed
                ## folder A's traces, while Shift+6's guard passed on the stale
                ## object and wrote folder A's peaks into folder B's
                ## peaks_monolist.csv. Every guard on these names must therefore be
                ## is.null(), NOT exists(): once declared they always exist.
                chromatograms_updated <- NULL
                .shiftq_fp <- NULL          # what the last full Shift+Q recompute saw (shiftq_fingerprint())
                x_axis_start <- NULL
                x_axis_end <- NULL
                y_axis_start <- NULL
                y_axis_end <- NULL
                MS_out_1 <- NULL
                MS_ret_start_line <- NULL
                MS_ret_end_line <- NULL
                framedDataFile_to_subtract <- NULL
                ## What the editable peak table was last rendered FROM. Shift+Z writes
                ## the browser's copy straight back, so if the file moved on since, that
                ## write silently reverts it. See render_peak_table() / Shift+Z.
                peak_table_snapshot <- NULL
                ms_wide <- NULL
                predictions <- NULL
                ## Bounded: never draw more than one page of facets at once, so the
                ## render can't blow up the graphics device on large sample sets.
                ## Panel geometry (begin) -- test_view_switch.R extracts this block.
                ## Pixel height of one sample's panel (times the sidebar's Panel height)
                ## and the size of the peak numbers.
                panel_px <- 100
                peak_label_size <- 3.88
                ## Panel geometry (end)

                ## x-axis defaults are resolved HERE, not in the Shift+Q handler, so
                ## they are always a finite length-1 number by the time any handler
                ## reads them: Shift+Q with no brush resets the view to them, and
                ## NULL or numeric(0) there would hand the renderer's dplyr::filter()
                ## a length-0 predicate.
                ## A caller-supplied value is kept as given.
                if ( is.null(x_axis_start_default) || length(x_axis_start_default) == 0 || !is.finite(x_axis_start_default) ) {
                    x_axis_start_default <- if (nrow(chromatograms) > 0) min(chromatograms$rt, na.rm = TRUE) else 0
                }
                if ( is.null(x_axis_end_default) || length(x_axis_end_default) == 0 || !is.finite(x_axis_end_default) ) {
                    x_axis_end_default <- if (nrow(chromatograms) > 0) max(chromatograms$rt, na.rm = TRUE) else 1
                }

            ## Set up new peak monolist if it doesn't exist

            if ( !file.exists("peaks_monolist.csv") ) {

                peak_data <- data.frame(
                  peak_start = 0,
                  peak_end = 0,
                  peak_ID = "unknown",
                  path_to_cdf_csv = "a",
                  area = 0
                )

                write.table(
                  x = peak_data[-1,],
                  file = "peaks_monolist.csv",
                  append = FALSE,
                  row.names = FALSE,
                  col.names = TRUE,
                  sep = ","
                )

            }

        ## SET UP USER INTERFACE

            ui <- fluidPage(

                theme = shinythemes::shinytheme("yeti"),

                ## Draggable sidebar. sidebarLayout() emits a fixed Bootstrap
                ## col-sm-4 / col-sm-8 pair, so the split cannot be changed from
                ## R. Convert that row to flexbox and put a grab handle between
                ## the two columns. On release we fire a window resize, which is
                ## what makes Shiny re-measure the plot outputs -- that also
                ## refreshes last_brush_size, which the brush maths depends on.
                ## Busy orb (top of the sidebar, 2026-10-21). Driven entirely in the
                ## browser by Shiny's own shiny:busy / shiny:idle events, because while
                ## Shift+Q runs the R side cannot send anything -- a server-driven
                ## indicator would only update once the work was over. Green "Ready",
                ## pulsing amber "Working... N s", grey "Disconnected".
                tags$head(
                    tags$style(HTML("
                        .gcms-status { display: flex; align-items: center; gap: 8px; font-weight: bold;
                                       margin: 2px 0 10px 0; font-size: 15px; }
                        .gcms-orb    { width: 16px; height: 16px; border-radius: 50%; flex: 0 0 16px;
                                       background: #2e9d4f; box-shadow: 0 0 6px #2e9d4f; }
                        .gcms-orb.gcms-on { background: #f0a30a;
                                       animation: gcmsPulse 1.1s ease-in-out infinite; }
                        .gcms-orb.gcms-off { background: #9e9e9e !important; box-shadow: none !important;
                                       animation: none !important; }
                        @keyframes gcmsPulse { 0%, 100% { transform: scale(0.8); box-shadow: 0 0 3px #f0a30a; }
                                               50% { transform: scale(1.15); box-shadow: 0 0 14px #f0a30a; } }
                        .gcms-unsaved { color: #c62828; font-weight: bold; margin: -4px 0 8px 0; }
                        .gcms-plot-row { display: flex; align-items: flex-start; }
                        .gcms-plot-cell { flex: 1 1 auto; min-width: 0; }
                        .gcms-align-col { flex: 0 0 38px; position: relative; }
                        .gcms-align-pair { position: absolute; left: 0; display: flex; gap: 1px;
                                           transform: translateY(-50%); }
                        .gcms-align-pair button { width: 18px; height: 16px; padding: 0; line-height: 14px;
                                                  border: 1px solid #bbb; border-radius: 3px; background: #f7f7f7;
                                                  font-size: 9px !important; color: #333; cursor: pointer; }
                        .gcms-align-pair button:hover { background: #2e7dd7; color: #fff; border-color: #2e7dd7; }
                        .gcms-ctx-menu { position: fixed; display: none; z-index: 10000; min-width: 120px;
                                         background: #ffffff; border: 1px solid #999999; border-radius: 4px;
                                         box-shadow: 0 2px 8px rgba(0,0,0,0.25); padding: 3px 0; font-size: 13px; }
                        .gcms-ctx-item { padding: 4px 14px; cursor: pointer; color: #222222; }
                        .gcms-ctx-item:hover { background: #2e7dd7; color: #ffffff; }
                    ")),
                    tags$script(HTML("
                        (function () {
                          // Shift+Q is a compute then a redraw: two busy spells a moment
                          // apart. Count spells less than ~1 s apart as ONE job, so the
                          // orb does not flash green in between and the time shown is
                          // the whole wait.
                          var t0 = null, timer = null, settle = null;
                          function txt(s) { $('#gcms_status_text').text(s); }
                          $(document).on('shiny:busy', function () {
                            if (settle) { clearTimeout(settle); settle = null; }
                            if (t0 === null) t0 = Date.now();
                            $('.gcms-orb').addClass('gcms-on');
                            if (!timer) {
                              txt('Working\u2026');
                              timer = setInterval(function () {
                                txt('Working\u2026 ' + Math.round((Date.now() - t0) / 1000) + ' s');
                              }, 500);
                            }
                          });
                          $(document).on('shiny:idle', function () {
                            if (settle) clearTimeout(settle);
                            settle = setTimeout(function () {
                              settle = null;
                              if (timer) { clearInterval(timer); timer = null; }
                              var secs = t0 === null ? 0 : Math.round((Date.now() - t0) / 1000);
                              t0 = null;
                              $('.gcms-orb').removeClass('gcms-on');
                              txt(secs >= 2 ? 'Ready (last took ' + secs + ' s)' : 'Ready');
                            }, 900);
                          });
                          $(document).on('shiny:disconnected', function () {
                            if (timer) { clearInterval(timer); timer = null; }
                            if (settle) { clearTimeout(settle); settle = null; }
                            $('.gcms-orb').removeClass('gcms-on').addClass('gcms-off'); txt('Disconnected');
                          });
                        })();
                    ")),
                    ## Align arrows: one pair per panel, placed at the panel's vertical centre
                    ## (rows from align_rows_for(), sent after every chromatogram render).
                    tags$script(HTML("
                        $(function () {
                          Shiny.addCustomMessageHandler('gcms_align_rows', function (m) {
                            var col = $('#gcms_align_col');
                            col.empty().css('height', (m.height || 0) + 'px');
                            (m.rows || []).forEach(function (r) {
                              var pair = $('<div class=\"gcms-align-pair\"></div>').css('top', r.y + 'px')
                                .attr('title', r.label + ': earlier / later by the Align step');
                              [-1, 1].forEach(function (dir) {
                                var b = $('<button type=\"button\"></button>').html(dir < 0 ? '&#9664;' : '&#9654;');
                                b.on('click', function () {
                                  Shiny.setInputValue('align_nudge', {sample: r.sample, dir: dir, t: Date.now()}, {priority: 'event'});
                                });
                                pair.append(b);
                              });
                              col.append(pair);
                            });
                          });
                        });
                    "))
                ),
                tags$head(
                    tags$style(HTML("
                        /* Two independent scroll boxes (Lucas, 2026-10-22): the page itself
                           never scrolls, so the wheel over the chromatograms moves only the
                           main column and the sidebar stays put (it scrolls on its own when
                           the wheel is over it). Sticky .gcms-sticky-top sticks within
                           .gcms-main, its nearest scrolling ancestor. */
                        html, body     { height: 100%; overflow: hidden; }
                        .gcms-flex-row { display: flex !important; align-items: stretch; height: 100vh; }
                        .gcms-sidebar  { flex: 0 0 auto !important; float: none !important;
                                         max-width: none !important; display: flex; flex-direction: column;
                                         overflow-y: auto; overscroll-behavior: contain; }
                        .gcms-sidebar > .well { flex: 1 0 auto; margin-bottom: 0; }
                        .gcms-main     { flex: 1 1 auto !important; float: none !important;
                                         width: auto !important; min-width: 0;
                                         overflow-y: auto; overscroll-behavior: contain; }
                        .gcms-drag     { flex: 0 0 6px; cursor: col-resize; background: #dcdcdc;
                                         border-radius: 3px; margin: 0 6px; }
                        .gcms-drag:hover { background: #9e9e9e; }
                        body.gcms-dragging { cursor: col-resize; user-select: none; }
                        /* One font size throughout: the theme's help-text size (12px). */
                        body, body * { font-size: 12px !important; }
                        .gcms-shortcuts { padding-left: 0; list-style-position: inside; margin-bottom: 0; }
                        .gcms-shortcuts li.gcms-sep { list-style: none; }
                        #key, #message_window { margin-bottom: 6px; }
                        /* Pinned spectra + buttons; the chromatograms scroll under them. */
                        .gcms-zoom-box .form-group { margin-bottom: 0; }
                        .gcms-sticky-top { position: -webkit-sticky; position: sticky; top: 0; z-index: 100;
                                           background: #fff; padding-bottom: 4px;
                                           box-shadow: 0 3px 4px -3px rgba(0,0,0,.25); }
                    ")),
                    tags$script(HTML("
                        $(function() {
                            var well = $('.well').first();
                            if (!well.length) return;
                            var side = well.parent();
                            var row  = side.parent();
                            var main = side.next();
                            row.addClass('gcms-flex-row');
                            side.addClass('gcms-sidebar').css('width', '340px');
                            main.addClass('gcms-main');
                            var handle = $('<div></div>').addClass('gcms-drag')
                                          .attr('title', 'Drag to resize the sidebar');
                            side.after(handle);
                            var dragging = false;
                            handle.on('mousedown', function(e) {
                                dragging = true;
                                $('body').addClass('gcms-dragging');
                                e.preventDefault();
                            });
                            $(document).on('mousemove', function(e) {
                                if (!dragging) return;
                                var w = e.pageX - side.offset().left;
                                var maxw = $(window).width() - 300;
                                if (w < 160) w = 160;
                                if (w > maxw) w = maxw;
                                side.css('width', w + 'px');
                            });
                            $(document).on('mouseup', function() {
                                if (!dragging) return;
                                dragging = false;
                                $('body').removeClass('gcms-dragging');
                                $(window).trigger('resize');
                            });
                        });
                    "))
                ),

                sidebarLayout(

                    ## The sidebar and the main column scroll independently (Lucas,
                    ## 2026-10-22), see .gcms-flex-row; the grey well fills the sidebar's
                    ## full height (.gcms-sidebar > .well).
                    sidebarPanel(

                        h4("Busta Lab GCMS Analysis App"),

                        tags$div(class = "gcms-status",
                            tags$span(class = "gcms-orb"),
                            tags$span(id = "gcms_status_text", "Ready")
                        ),
                        ## Red line while the peak table holds edits Shift+Z has not saved.
                        tags$div(class = "gcms-unsaved", textOutput("unsaved_flag", inline = TRUE)),
                        ## Messages now pop up as notifications (notify()). The box stays in
                        ## the page, hidden; it does not fill while hidden, so headless drivers
                        ## time off the busy orb instead.
                        tags$div(style = "display: none;", verbatimTextOutput("message_window", placeholder = TRUE)),

                        tags$hr(),
                        strong("Keyboard Shortcuts:"),
                        tags$ul(class = "gcms-shortcuts",
                            tags$li("Shift + Q => Zoom to the box; no box = full view (recalculates first if anything changed)"),
                            tags$li("Recalculate button => recalculate, keep the view"),
                            tags$li(class = "gcms-sep", "------------------"),
                            tags$li("Shift + 6 => Detect peaks"),
                            tags$li("Shift + A => Add peak (every panel the brush spans)"),
                            tags$li("Shift + G => Add global peak"),
                            tags$li("Shift + E => Excise peak (every panel the brush spans)"),
                            tags$li("Right-click a peak => Remove it, or split it (that sample only)"),
                            tags$li("Shift + R => Remove global peak"),
                            tags$li("Shift + 0 => Remove ALL peaks"),
                            tags$li("Shift + Z => Save peak table"),
                            tags$li(class = "gcms-sep", "------------------"),
                            tags$li("Shift + 1 => Extract MS from selection"),
                            tags$li("Shift + 2 => Zoom in on selected MS portion"),
                            tags$li("Shift + 3 => Subtract selection MS from current MS"),
                            tags$li("Shift + 4 => Library search"),
                            tags$li("Shift + 5 => Save current MS")
                        ),

                        tags$hr(),
                        strong("View"),
                        numericInput(
                            "samples_per_page_ui", "Samples per page",
                            value = samples_per_page, min = 1, max = 500, step = 1
                        ),
                        sliderInput(
                            "panel_height_mult", "Panel height",
                            min = 0.5, max = 4, value = 1, step = 0.25, post = "x"
                        ),
                        helpText("Panel height stretches every panel; it redraws when you let go of the slider."),
                        numericInput("align_step", "Align step (the \u25c0 \u25b6 beside each panel)",
                                     value = 1, min = 0.001, step = 0.5, width = "100%"),

                        tags$hr(),
                        strong("Align on a peak"),
                        helpText("Give a peak number and press Align: every sample that has that peak is shifted so the centres of their integration windows line up (on the median). This replaces each sample's offset; the \u25c0 \u25b6 arrows then fine-tune from there. Samples without that peak are left as they are. Reset all offsets puts every sample back to zero."),
                        numericInput("alignall_peak", "Peak number", value = NA, min = 1, step = 1),
                        fluidRow(
                            column(6, actionButton("alignall_btn", "Align")),
                            column(6, actionButton("align_reset_btn", "Reset all offsets"))
                        ),

                        tags$hr(),
                        strong("Chromatogram display"),
                        selectInput(
                            "display_mode", NULL,
                            choices = c("Total ion current" = "tic",
                                        "All ions (thin lines)" = "all_ions",
                                        "Ion map (heatmap)" = "ion_map"),
                            selected = if (isTRUE(all_ion_view)) "all_ions" else "tic"
                        ),
                        selectInput(
                            "allion_scaling", "Y scaling",
                            choices = c("Raw counts" = "raw",
                                        "Square root" = "sqrt",
                                        "Log10" = "log",
                                        "Each ion to its own max" = "norm"),
                            selected = all_ion_scaling
                        ),
                        tags$hr(),
                        strong("All-ion settings (apply on Shift+Q)"),
                        sliderInput(
                            "allion_threshold", "Ion inclusion (fraction of the window's largest ion)",
                            min = 0.001, max = 0.2, value = all_ion_threshold, step = 0.001
                        ),
                        sliderInput(
                            "allion_max_points", "Max x-points per ion (render cost)",
                            min = 200, max = 4000, value = all_ion_max_points, step = 100
                        ),

                        tags$hr(),
                        strong("Co-elution split"),
                        helpText("Integrate the hump as ONE peak and press Shift+Q. Then give its peak number and a diagnostic ion (or several) for each of the two compounds, in any order, and press Split; fill compound 3 as well to split it three ways. Pick ions that are strong in every scan of the hump: a sparse ion gives a noisy split. The hump becomes two (or three) peaks in every sample, its TIC shared scan by scan by the ions; a sample with no signal on one compound's ions still gets that peak, with area 0, so the numbers line up. Each split keeps its own settings, and Shift+Q only refreshes it. Unsplit merges the split a peak number belongs to back into one peak."),
                        numericInput("split_peak", "Peak number", value = NA, min = 1, step = 1),
                        textInput("split_ions_a", "Ions, compound 1 (e.g. 282, 281)", value = ""),
                        textInput("split_ions_b", "Ions, compound 2 (e.g. 111, 125)", value = ""),
                        textInput("split_ions_c", "Ions, compound 3 (optional)", value = ""),
                        radioButtons("split_k_mode", "Response factors",
                                     choices = c("Fit per sample" = "sample", "Fit across the batch" = "batch"),
                                     selected = "sample", inline = TRUE),
                        fluidRow(
                            column(6, actionButton("split_btn", "Split peak")),
                            column(6, actionButton("unsplit_btn", "Unsplit peak"))
                        ),

                    ),

                    mainPanel(

                        ## Typing in a text or number box must not fire the Shift shortcuts
                        ## (a capital A in an ion box added peaks). The peak table is exempt:
                        ## Shift+Z saves it straight after an edit there.
                        ## HTML() is load-bearing: a bare string is HTML-escaped, so && arrives
                        ## as &amp;&amp;, the script fails to parse and no shortcut works
                        ## (2026-10-27, test_shortcut_text_inputs.R).
                        tags$script(HTML('
                            function gcmsTypingInField(t) {
                               if (!t || !t.tagName) return false;
                               if (t.closest && t.closest(".handsontable, .HandsontableCopyPaste")) return false;
                               return /^(INPUT|TEXTAREA|SELECT)$/.test(t.tagName) || !!t.isContentEditable;
                            }
                            $(document).on("keypress", function (e) {
                               if (gcmsTypingInField(e.target)) return;
                               // priority "event": Shiny otherwise drops a value equal to the
                               // last one sent, so a second Shift+Q in a row did nothing.
                               Shiny.setInputValue("keypress", e.which, {priority: "event"});
                            });
                        ')),

                        ## Right-click menu on the chromatograms (Lucas, 2026-10-26): Remove peak,
                        ## Split peak... Server side: the Right-click menu block.
                        tags$script(HTML("
                            (function () {
                              var menu = null, at = null;
                              function hide() { if (menu) menu.style.display = 'none'; at = null; }
                              function build() {
                                menu = document.createElement('div');
                                menu.className = 'gcms-ctx-menu';
                                [['Remove peak', 'remove'], ['Split peak\u2026', 'split']].forEach(function (it) {
                                var item = document.createElement('div');
                                item.className = 'gcms-ctx-item';
                                item.textContent = it[0];
                                item.addEventListener('click', function (e) {
                                  e.stopPropagation();
                                  var p = at; hide();
                                  if (!p || !window.Shiny) return;
                                  // Re-send the right-click's spot through Shiny's own click
                                  // handler: the same jQuery event Shiny triggers for a left
                                  // click. Cleared first, because a spot outside every panel
                                  // sends nothing and would leave an older click standing.
                                  Shiny.setInputValue('chromatogram_click', null, {priority: 'event'});
                                  $('#chromatograms').trigger($.Event('mousedown2', {which: 1, pageX: p.pageX, pageY: p.pageY}));
                                  Shiny.setInputValue('chromatogram_ctx', {action: it[1], n: Date.now()}, {priority: 'event'});
                                });
                                menu.appendChild(item);
                                });
                                document.body.appendChild(menu);
                              }
                              $(document).on('contextmenu', '#chromatograms', function (e) {
                                e.preventDefault();
                                if (!menu) build();
                                at = {pageX: e.pageX, pageY: e.pageY};
                                menu.style.left = e.clientX + 'px';
                                menu.style.top = e.clientY + 'px';
                                menu.style.display = 'block';
                              });
                              $(document).on('mousedown', function (e) {
                                if (menu && !$(e.target).closest('.gcms-ctx-menu').length) hide();
                              });
                              $(document).on('keydown', function (e) { if (e.key === 'Escape') hide(); });
                              document.addEventListener('scroll', hide, true);
                              $(window).on('blur resize', hide);
                            })();
                        ")),


                        tabsetPanel(type = "tabs",

                            tabPanel("Main",

                                ## The mass-spectrum panel and the button row stay pinned
                                ## at the top while the chromatograms scroll beneath them (Lucas,
                                ## 2026-10-22). The second spectrum panel (massSpectra_2, never
                                ## rendered to) was taken off the screen on 2026-10-26 (Lucas). Rows built by hand, not with verticalLayout(),
                                ## because that wraps every child in its own row and a sticky
                                ## element can only stick within its parent: the pinned block
                                ## has to be ONE row whose parent is the whole tab.
                                tags$div(class = "row gcms-sticky-top", column(12,

                                    plotOutput(
                                        outputId = "massSpectra_1",
                                        brush = brushOpts(
                                            id = "massSpectra_1_brush",
                                            resetOnNew = TRUE
                                        ),
                                        height = 150
                                    ),

                                    ## Pagination controls: draw only samples_per_page
                                    ## samples at a time and step through pages.
                                    fluidRow(
                                        column(12,
                                            actionButton("prev_page", paste0("\u25c0 Prev ", samples_per_page)),
                                            actionButton("next_page", paste0("Next ", samples_per_page, " \u25b6")),
                                            ## Slide the retention-time window by half its width.
                                            actionButton("pan_back", "\u25c0 Earlier", style = "margin-left:15px;"),
                                            actionButton("pan_fwd", "Later \u25b6"),
                                            ## Step back and forth through past zoom windows.
                                            actionButton("view_back", "\u25c0 Previous view", style = "margin-left:15px;"),
                                            actionButton("view_fwd", "Subsequent view \u25b6"),
                                            ## Shift+Q's recompute without its zoom: the view on screen stays.
                                            actionButton("recalc_btn", "Recalculate", style = "margin-left:15px;"),
                                            ## Put the peak list and offsets back as they were before the last edit.
                                            actionButton("undo_btn", "Undo", style = "margin-left:15px;"),
                                            ## Zoom to a run of peaks by number ("12", "12-15").
                                            tags$span(class = "gcms-zoom-box", style = "margin-left:15px; display:inline-block; width:90px; vertical-align:middle;",
                                                textInput("zoom_peaks", label = NULL, placeholder = "e.g. 12-15", width = "90px")),
                                            actionButton("zoom_peaks_btn", "Zoom to peaks"),
                                            tags$span(
                                                style = "margin-left:15px; font-weight:bold;",
                                                textOutput("page_indicator", inline = TRUE)
                                            ),
                                            ## The TIC view's ion key (filled by the renderer).
                                            uiOutput("ion_legend", inline = TRUE),
                                            ## Key to the squares on each peak (the plot cannot carry it:
                                            ## its fill scale belongs to the peak ribbons).
                                            gcms_confidence_legend(confidence_thresholds)
                                        )
                                    )
                                )),
                                fluidRow(column(12,
                                    ## Align (Lucas, 2026-10-21): a pair of arrows beside EACH panel moves
                                    ## that sample along rt by the sidebar's Align step. The column is
                                    ## filled by the renderer (gcms_align_rows message, positions from
                                    ## align_rows_for()); a click sends input$align_nudge.
                                    tags$div(class = "gcms-plot-row",
                                    tags$div(id = "gcms_align_col", class = "gcms-align-col"),
                                    tags$div(class = "gcms-plot-cell",
                                    plotOutput(
                                        outputId = "chromatograms",
                                        ## resetOnNew: the box is cleared whenever the plot
                                        ## re-renders (Lucas, 2026-10-04). A brush is drawn in the
                                        ## coordinates of the picture it was drawn ON; once Shift+Q,
                                        ## a page change or a display switch redraws it, the same
                                        ## box means a different region, and the next keypress used
                                        ## to act on it anyway. The renderer does not depend on the
                                        ## brush, so this cannot loop.
                                        brush = brushOpts(
                                            id = "chromatogram_brush",
                                            resetOnNew = TRUE
                                        ),
                                      click = "chromatogram_click",
                                      dblclick = "chromatogram_double_click",
                                      ## The height is set per render by renderPlot(height =
                                      ## chromatogram_height), so it follows the mode and page.
                                      height = "auto"
                                    ))))),

                                fluidRow(column(12,
                                    rhandsontable::rHandsontableOutput("peak_table")
                                )),

                                tags$head(
                                    HTML(
                                        "
                                        <script>
                                            var socket_timeout_interval
                                            var n = 0
                                            $(document).on('shiny:connected', function(event) {
                                            socket_timeout_interval = setInterval(function(){
                                            Shiny.onInputChange('count', n++)
                                            }, 15000)
                                            });
                                            $(document).on('shiny:disconnected', function(event) {
                                            clearInterval(socket_timeout_interval)
                                            });
                                        </script>
                                        "
                                    )
                                ),

                                textOutput("keepAlive")
                            ),

                            tabPanel("Peak Detection",
                              fluidRow(
                                column(3,
                                  sliderInput("peakStartSlope", "Peak Start Slope",
                                              min = 1, max = 100, value = 10, step = 1),
                                  sliderInput("apexSlope", "Apex Slope",
                                              min = -20, max = 20, value = 0, step = 1),
                                  sliderInput("peakEndSlope", "Peak Tail Slope",
                                              min = 1, max = 100, value = 5, step = 1),
                                  sliderInput("minSignalAboveBaseline", "Min Signal Above Baseline",
                                              min = 0, max = 50, value = 2, step = 1),
                                  sliderInput("minPeakArea", "Minimum Peak Area",
                                              min = 0, max = 1e7, value = 1e5, step = 1e4)
                                ),
                                column(9,
                                  h4("Press Shift + 6 to auto-detect peaks using these settings."),
                                  p("peakStartSlope: how steep the slope must be to call a peak start."),
                                  p("Once derivative goes < 0 => apex is reached."),
                                  p("peakEndSlope: how negative the slope can stay before we say the peak is done."),
                                  p("minSignalAboveBaseline: how close to baseline to consider the tail done."),
                                  p("minPeakArea: if the integrated area is below this, ignore the peak."),
                                  plotOutput("peak_detection_plot", height = "600px")
                                )
                              )
                            ),


                            tabPanel("MS Library",

                                verticalLayout(
                                    fluidRow(
                                        column(3,
                                            actionButton("extract_ms", "Extract Mass Spectra"),
                                            br(), br(),
                                            actionButton("classify_ms", "Classify Spectra")
                                        ),
                                        column(9,
                                            DT::dataTableOutput("ms_results")
                                        )
                                    ),

                                    plotOutput(
                                        outputId = "massSpectrumLookup",
                                        height = 1200
                                    )
                                )
                            )
                        )
                    )
                )
            )

        ## SET UP SERVER

            server <- function(input, output, session) {

                ## Keyboard map. The browser sends e.which from a keypress event, so a
                ## letter code is the UPPER-case letter (Shift+letter, or Caps Lock) and a
                ## shifted digit arrives as its symbol. The handlers follow below; Shift+Q/W
                ## share one observer and the four mass-spectrum keys share another.
                ##   Shift+Q              recompute if anything changed, then zoom    (81)
                ##                        (the Recalculate button: the same recompute, view kept)
                ##   Shift+W              always recompute                            (87)
                ##   Shift+Z              save peak_ID edits from the table           (90)
                ##   Shift+A / Shift+G    add the brushed peak: brushed samples / all (65 / 71)
                ##   Shift+E / Shift+R    remove peaks inside the brush: brushed / all (69 / 82)
                ##                        ("brushed samples": every panel whose middle the
                ##                        brush covers, else the one panel it sits in)
                ##   Shift+6 (^)          automatic peak detection                    (94)
                ##   Shift+0 ())          remove ALL peaks                            (41)
                ##   Shift+1 (!)          mass spectrum of the brushed region         (33)
                ##   Shift+2 (@)          re-plot the mass spectrum                   (64)
                ##   Shift+3 (#)          subtract the brushed region from it         (35)
                ##   Shift+4 ($)          classify every peak (ms_classifier.py)      (36)
                ## Undo is a button (no key): it reverts the last peak/offset edit, see the
                ## Undo block below. Right-click on a peak opens a menu (Remove peak; Split
                ## peak, this sample only), see the Right-click menu block. Panning is the Earlier / Later buttons (the old Shift+F/D/V/C keys were
                ## removed on 2026-10-22).

                ## Store console info:
                    message_data <- reactiveVal("")

                    # Helper function to append a new line
                    logMessage <- function(msg) {
                        msg <- trimws(gsub("\n", " ", msg))
                        if (!nzchar(msg)) return(invisible(NULL))
                        ## isolate(): notify() is also called from inside renders, which must not
                        ## come to depend on the message box.
                        old <- strsplit(isolate(message_data()), "\n", fixed = TRUE)[[1]]
                        message_data(paste(head(c(msg, old[nzchar(old)]), 8), collapse = "\n"))
                    }

                    ## An action's outcome: to the console (Jupyter output) AND the
                    ## Messages box, newest at the top.
                    notify <- function(msg) {
                        cat(msg, "\n", sep = "")
                        logMessage(msg)
                        ## The Messages box is hidden, so each outcome pops up briefly.
                        try(showNotification(msg, duration = 5, type = "message"), silent = TRUE)
                    }

                    ## Split/refresh logs: everything to the console; to the box,
                    ## all but the per-sample factor lines.
                    notify_log <- function(lines) {
                        if (!length(lines)) return(invisible(NULL))
                        cat(paste0(lines, "\n"), sep = "")
                        for (l in rev(lines[!grepl("^  [^ ]+: kA ", lines)])) logMessage(l)
                        ## The Messages box is hidden, so the outcome pops up: every headline
                        ## (unindented line) and a count of skipped samples with the first
                        ## reasons. Before this (2026-10-21) a split ran or refused silently.
                        heads <- lines[!grepl("^\\s", lines) & nzchar(trimws(lines))]
                        skips <- trimws(lines[grepl("^  skipped ", lines)])
                        for (h in heads) try(showNotification(h, duration = 10, type = "message"), silent = TRUE)
                        if (length(skips)) {
                            msg <- paste0(length(skips), " sample", if (length(skips) > 1) "s" else "", " skipped: ",
                                          paste(utils::head(skips, 3), collapse = "; "),
                                          if (length(skips) > 3) "; ... (all listed in the R console)" else "")
                            try(showNotification(msg, duration = 15, type = "warning"), silent = TRUE)
                        }
                    }

                    # Render them in the message window
                    output$message_window <- renderText({
                        message_data()
                    })

                ## Set up MS values
                    ms_data_reactive <- reactiveVal(NULL)
                    ms_results_reactive <- reactiveVal(NULL)

                ## Pagination state. current_page = which page of samples is shown;
                ## redraw_trigger is bumped by Shift+Q so the (reactive) chromatogram
                ## plot recomputes only when the data actually changes, not on every brush.
                    current_page   <- reactiveVal(1)
                    redraw_trigger <- reactiveVal(0)

                ## The ion key in the button row: the renderer sets the named ion -> colour
                ## map it drew the TIC view with (NULL in the all-ion views).
                    ion_legend_pal <- reactiveVal(NULL)
                    output$ion_legend <- renderUI(gcms_ion_legend(ion_legend_pal()))

                ## Don't let it time out
                    output$keepAlive <- renderText({
                        req(input$count)
                        paste("\nstayin' alive ", input$count)
                    })


                ## Set up null values for peak detection
                    peakDetectionPlotData <- reactiveVal(NULL)

                        output$peak_detection_plot <- renderPlot({
                        # If there's no data yet, do nothing
                        df <- peakDetectionPlotData()
                        validate(
                          need(!is.null(df), "No random peak subset to display yet. Run Shift+6 to detect peaks.")
                        )

                        # Plot them in facet_wrap:
                        ggplot(df, aes(x = rt, y = abundance, color = ion)) +
                          geom_line() +
                          geom_line(aes(y = baseline), color = "gray50") +
                          geom_ribbon(aes(ymin = baseline, ymax = abundance, fill = peak_label, group = peak_label),
                                      alpha = 0.2, color = NA) +
                          facet_wrap(~peak_label, scales = "free") +
                          theme_bw() +
                          labs(title = "Randomly Sampled Peak Regions",
                               x = "Scan (RT)", y = "Abundance (counts)") +
                          guides(fill = "none")
                    })

                ## Unsaved-edits flag: compares the browser's peak_IDs with what the table
                ## was drawn from (gcms_peak_ids_edited()). Cleared by a re-render or a save.
                    table_dirty <- reactiveVal(FALSE)
                    observeEvent(input$peak_table, {
                        ed <- try(rhandsontable::hot_to_r(input$peak_table), silent = TRUE)
                        if (inherits(ed, "try-error")) return()
                        table_dirty(gcms_peak_ids_edited(ed, peak_table_snapshot))
                    }, ignoreNULL = TRUE)
                    output$unsaved_flag <- renderText({
                        if (table_dirty()) "\u25cf Unsaved peak-table edits: Shift+Z to save" else ""
                    })

                ## ONE way to render the editable peak table. Shift+G used to render it
                ## as DT::renderDataTable while every other handler used rhandsontable --
                ## and `input$peak_table` from a DT output carries nothing hot_to_r() can
                ## read, so Shift+Z silently did nothing for the rest of the session.
                ## Each render also records what it rendered FROM, so Shift+Z can tell
                ## an edit from a revert.
                    render_peak_table <- function() {
                        current <- if (file.exists("peaks_monolist.csv")) read.csv("peaks_monolist.csv") else NULL
                        peak_table_snapshot <<- current
                        table_dirty(FALSE)      # a fresh table holds no edits
                        ## Only peak_ID is editable; everything else is computed.
                        output$peak_table <- rhandsontable::renderRHandsontable({
                            if (is.null(current)) return(NULL)
                            hot <- rhandsontable::rhandsontable(gcms_peak_table_display(current))
                            hot <- rhandsontable::hot_cols(hot, readOnly = TRUE)
                            if ("peak_ID" %in% names(current)) hot <- rhandsontable::hot_col(hot, "peak_ID", readOnly = FALSE)
                            hot
                        })
                        invisible(NULL)
                    }

                ## Save manual changes to table on Shift+Z (90)

                    observeEvent(input$keypress, {
                        if (input$keypress == 90 ) {
                            undo_begin("save peak-ID edits (Shift+Z)"); on.exit(undo_end(), add = TRUE)
                            ## Write out any modifications to peak table (i.e. sample IDs)
                                hot = isolate(input$peak_table)
                                if (is.null(hot)) {
                                    notify("Nothing to save yet. Press Shift+Q first.")
                                    return()
                                }
                                edited <- rhandsontable::hot_to_r(hot)

                                ## The browser holds a SNAPSHOT from the last render. If the file
                                ## has moved on since (a Shift+A, a Shift+6, another Shift+Q), then
                                ## writing that snapshot back reverts everything added in between.
                                ## Refuse rather than guess: a cell-level merge would need a stable
                                ## row key, and this table has none that the user cannot also edit.
                                on_disk <- if (file.exists("peaks_monolist.csv")) read.csv("peaks_monolist.csv") else NULL
                                moved_on <- !is.null(peak_table_snapshot) && !is.null(on_disk) &&
                                            !identical(dim(peak_table_snapshot), dim(on_disk))
                                if (moved_on) {
                                    notify(paste0("Not saved: the peak list changed since the table was drawn (",
                                               nrow(peak_table_snapshot), " rows then, ", nrow(on_disk),
                                               " now). Press Shift+Q, then redo your edits."))
                                    return()
                                }

                                saved <- gcms_apply_peak_id_edits(on_disk, edited)
                                if (is.null(saved)) {
                                    notify("Not saved: the table is out of date. Press Shift+Q, then redo your edits.")
                                    return()
                                }
                                writeMonolist(saved, "peaks_monolist.csv")
                                peak_table_snapshot <<- saved
                                table_dirty(FALSE)
                                notify("Peak list saved.")
                        }
                    })

                ## Delete peaks (Shift+0)
                    observeEvent(input$keypress, {
                      # SHIFT+0 => `)` => ASCII code 41
                      if (input$keypress == 41) {
                          undo_begin("remove all peaks (Shift+0)"); on.exit(undo_end(), add = TRUE)
                        cat("Removing ALL peaks...\n")
                        if (file.exists("peaks_monolist.csv")) {
                          # Overwrite with empty table
                          empty_peaks <- data.frame(
                            peak_start = numeric(),
                            peak_end   = numeric(),
                            peak_ID    = character(),
                            path_to_cdf_csv = character(),
                            area       = numeric(),
                            stringsAsFactors = FALSE
                          )
                          write.table(
                            empty_peaks,
                            "peaks_monolist.csv",
                            row.names = FALSE,
                            col.names = TRUE,
                            sep = ","
                          )
                          notify("All peaks removed.")

                          # ALSO CLEAR THE facet-plot reactiveVal
                          peakDetectionPlotData(NULL)

                        } else {
                          notify("No peaks to remove.")
                        }
                      }
                    })

                ## Peak detection
                    observeEvent(input$keypress, {
                      if (input$keypress == 94) {  # '^' = Shift+6
                          undo_begin("peak detection (Shift+6)"); on.exit(undo_end(), add = TRUE)
                        peakDetectionPlotData(NULL)
                        cat("Starting custom peak detection...\n")

                        # 1) Read thresholds from your sliders
                        peakStartSlope         <- input$peakStartSlope
                        apexSlope              <- input$apexSlope    # might be zero, used for crossing detection
                        peakTailSlope          <- input$peakEndSlope # how negative it can be before flattening
                        minSignalAboveBaseline <- input$minSignalAboveBaseline
                        minPeakArea            <- input$minPeakArea  # default = 100,000

                        # 2) Load any existing peaks (so we can append to them)
                        if (file.exists("peaks_monolist.csv")) {
                          peak_data_existing <- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                        } else {
                          peak_data_existing <- data.frame(
                            peak_start       = numeric(),
                            peak_end         = numeric(),
                            peak_ID          = character(),
                            path_to_cdf_csv  = character(),
                            area             = numeric(),
                            stringsAsFactors = FALSE
                          )
                        }

                        peak_data_new <- data.frame(
                          peak_start       = numeric(),
                          peak_end         = numeric(),
                          peak_ID          = character(),
                          path_to_cdf_csv  = character(),
                          area             = numeric(),
                          stringsAsFactors = FALSE
                        )

                        # 3) Ensure chromatograms_updated exists
                        if (is.null(chromatograms_updated)) {
                          notify("No chromatograms loaded yet. Press Shift+Q first.")
                          return()
                        }

                        # 4) If you have sample filters
                        samples_monolist <- read.csv("samples_monolist.csv", stringsAsFactors = FALSE)
                        if (length(samples_monolist_subset) > 0) {
                          samples_monolist <- samples_monolist[samples_monolist_subset, ]
                        }

                        # 5) Loop through each sample’s TIC
                        for (this_file in unique(samples_monolist$path_to_cdf_csv)) {

                          # Pull data for this sample
                          this_chrom <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == this_file)
                          # TIC only
                          tic_data <- dplyr::filter(this_chrom, ion == 0)
                          if (nrow(tic_data) < 3) next

                          # Sort by RT
                          tic_data <- tic_data[order(tic_data$rt), ]

                          # Subtract baseline
                          baseline_data <- dplyr::filter(this_chrom, ion == "baseline")
                          baseline_data <- baseline_data[order(baseline_data$rt), ]

                          merged_data <- merge(
                            tic_data,
                            baseline_data[, c("rt", "abundance")],
                            by = "rt",
                            suffixes = c("", ".bl")
                          )
                          merged_data$net_abundance <- merged_data$abundance - merged_data$abundance.bl

                          # Derivative of net_abundance
                          d1 <- c(0, diff(merged_data$net_abundance))

                          # 6) State machine
                          outsidePeak <- 0
                          ascending   <- 1
                          descending  <- 2

                          currentState <- outsidePeak
                          peak_start_idx <- NA

                          for (i in seq_along(d1)) {
                            if (currentState == outsidePeak) {
                              # Start if derivative > peakStartSlope
                              if (d1[i] > peakStartSlope) {
                                currentState <- ascending
                                peak_start_idx <- i
                              }

                            } else if (currentState == ascending) {
                              # If derivative crosses below apexSlope => apex reached => begin descending
                              if (d1[i] < apexSlope) {
                                currentState <- descending
                              }

                            } else if (currentState == descending) {
                              # If slope is now flattening out or net_abundance is near baseline => end peak
                              # i.e. derivative > -peakTailSlope OR net_abundance < minSignalAboveBaseline
                              if (d1[i] > -peakTailSlope || merged_data$net_abundance[i] < minSignalAboveBaseline) {

                                # Record the peak
                                rt_start <- merged_data$rt[peak_start_idx]
                                rt_end   <- merged_data$rt[i]
                                area_val <- sum(merged_data$net_abundance[peak_start_idx:i])
                                if (area_val < 0) area_val <- 0

                                # Apply the minPeakArea filter
                                if (area_val >= minPeakArea) {
                                  peak_data_new <- rbind(
                                    peak_data_new,
                                    data.frame(
                                      peak_start      = rt_start,
                                      peak_end        = rt_end,
                                      peak_ID         = "auto_detected",
                                      path_to_cdf_csv = this_file,
                                      area            = area_val,
                                      stringsAsFactors = FALSE
                                    )
                                  )
                                }

                                # Reset
                                currentState <- outsidePeak
                                peak_start_idx <- NA
                              }
                            }
                          }

                          # If we exit while still in ascending or descending, end at the last data point
                          if (currentState != outsidePeak && !is.na(peak_start_idx)) {
                            last_idx <- nrow(merged_data)
                            rt_start <- merged_data$rt[peak_start_idx]
                            rt_end   <- merged_data$rt[last_idx]
                            area_val <- sum(merged_data$net_abundance[peak_start_idx:last_idx])
                            if (area_val < 0) area_val <- 0

                            if (area_val >= minPeakArea) {
                              peak_data_new <- rbind(
                                peak_data_new,
                                data.frame(
                                  peak_start      = rt_start,
                                  peak_end        = rt_end,
                                  peak_ID         = "auto_detected_unfinished",
                                  path_to_cdf_csv = this_file,
                                  area            = area_val,
                                  stringsAsFactors = FALSE
                                )
                              )
                            }
                          }

                        } # end for each file

                        # 7) Append the new peaks to existing, write out
                            if (nrow(peak_data_new) > 0) {
                              # bind_rows, not rbind: once a Shift+Q has run, the on-disk table carries extra
                              # columns (peak_number_within_sample, rt_offset, area_peak_ions, ...) that
                              # peak_data_new does not, and rbind errors on the mismatch. The new rows
                              # get NA in those columns; the next Shift+Q recomputes them.
                              peak_data_combined <- dplyr::bind_rows(peak_data_existing, peak_data_new)
                              write.table(
                                peak_data_combined,
                                file = "peaks_monolist.csv",
                                row.names = FALSE,
                                col.names = TRUE,
                                sep = ","
                              )
                              notify(paste0("Peaks detected: ", nrow(peak_data_new), "."))
                            } else {
                              notify("No new peaks found.")
                              # No new peaks => might not want to do the random plot:
                              return()
                            }

                            ##
                            ## NEW: Build random subset of ~30 peaks for facet plotting
                            ##
                            # 1) Reload final table from disk (or just reuse `peak_data_combined`)
                            peak_table_final <- read.csv("peaks_monolist.csv")

                            # 2) If there are more than 30 peaks, sample 30
                            peak_table_sampled <- peak_table_final %>%
                              dplyr::arrange(desc(area)) %>%
                              dplyr::slice_head(n=30)

                            # 3) For each sampled peak, gather the chromatogram data
                            all_peaks_data <- list()

                            # Make sure we have 'chromatograms_updated' loaded
                            if (is.null(chromatograms_updated)) {
                              notify("No chromatograms loaded yet. Press Shift+Q first.")
                              return()
                            }

                            for (i in seq_len(nrow(peak_table_sampled))) {

                              rowi <- peak_table_sampled[i, ]

                              # Subset data from chromatograms_updated
                              # Only the relevant file, plus RT slice
                              df_signal <- dplyr::filter(
                                chromatograms_updated,
                                path_to_cdf_csv == rowi$path_to_cdf_csv,
                                rt >= rowi$peak_start,
                                rt <= rowi$peak_end
                              )

                              # If no data, skip
                              if (nrow(df_signal) == 0) next

                              # Optionally keep only ion==0 and ion=="baseline" if you want
                              # so the plot isn't cluttered. Or keep all if you prefer.
                              # Let's keep both TIC (ion=0) and baseline (ion="baseline"):
                              df_signal <- df_signal[df_signal$ion %in% c(0, "baseline"), ]

                              # Add a facet label, e.g. "Peak #1"
                              df_signal$peak_label <- paste0("Peak #", i)

                              all_peaks_data[[i]] <- df_signal
                            }

                            # 4) Combine into one big data frame
                            if (length(all_peaks_data) > 0) {
                              all_peaks_data <- dplyr::bind_rows(all_peaks_data)

                              # To help the plotting, let's define a separate 'baseline' column
                              # For baseline rows, abundance is already the baseline,
                              # For TIC rows, let's store the baseline in a new column.
                              # We'll match them by RT if we want. Or we can do a simpler approach:

                              # For convenience, let's do a small join to assign a 'baseline' column
                              # to the entire data, so we can do a ribbon from baseline->abundance.

                              df_baselines <- dplyr::filter(all_peaks_data, ion=="baseline") %>%
                                dplyr::select(rt, path_to_cdf_csv, peak_label, baseline=abundance)

                              # For TIC rows
                              df_signal_only <- dplyr::filter(all_peaks_data, ion==0)

                              # left_join by path_to_cdf_csv, rt, peak_label
                              df_signal_joined <- dplyr::left_join(
                                df_signal_only,
                                df_baselines,
                                by = c("rt", "path_to_cdf_csv", "peak_label")
                              )

                              # For baseline rows themselves, we can keep them separate or combine them,
                              # but the simplest might be: "baseline" is the same as "abundance" for baseline.
                              # Let’s store them in the same df, though you can do separate geoms if you prefer.
                              df_baselines$baseline <- df_baselines$baseline # unchanged
                              df_baselines$ion <- "baseline"
                              # We'll rename abundance to something else to not confuse the ribbon
                              df_baselines$abundance <- df_baselines$baseline

                              # Combine them
                              df_plot <- dplyr::bind_rows(df_signal_joined, df_baselines)

                              # 5) Store in reactiveVal => triggers the plot
                              peakDetectionPlotData(df_plot)
                            } else {
                              # If no peaks or no data slices, set to NULL
                              peakDetectionPlotData(NULL)
                            }

                      }


                    })

                ## Update chromatogram on Shift+Q (81)

                    ## Everything Shift+Q's recompute reads that can change mid-session: the two
                    ## monolists (every edit path - Shift+A/G/6/Z, Split, Align - writes one of
                    ## them), the live ion threshold and the extracted chromatograms. When none
                    ## has changed since the last recompute, Shift+Q only zooms: most presses are
                    ## zooms, and the recompute was minutes at 95 samples (2026-10-21). Content
                    ## hashes, not mtimes, so two writes in the same second cannot collide.
                    ## Shift+W always recomputes. Pinned by test_shiftq_skip.R.
                    shiftq_fingerprint <- function() {
                        md5 <- function(f) if (file.exists(f)) unname(tools::md5sum(f)) else ""
                        thr <- if (is.null(input$allion_threshold)) all_ion_threshold else input$allion_threshold
                        c(samples = md5("samples_monolist.csv"), peaks = md5("peaks_monolist.csv"),
                          threshold = format(thr, digits = 17),
                          chromatograms = paste(nrow(chromatograms), ncol(chromatograms)))
                    }

                    ## One body for Shift+Q, Shift+W and the Recalculate button (Lucas,
                    ## 2026-10-22). force = TRUE (Shift+W) recomputes even when the
                    ## fingerprint matches. keep_view = TRUE (Recalculate) recomputes when
                    ## something changed and redraws the view already on screen: it never
                    ## reads the brush, so Shift+Q is the zoom and Recalculate the refresh.
                    shiftq_run <- function(force = FALSE, keep_view = FALSE) {

                          shiftq_recompute <- force || is.null(chromatograms_updated) ||
                              !identical(shiftq_fingerprint(), .shiftq_fp)
                          if (!shiftq_recompute) {
                              cat(if (keep_view) "Recalculate: nothing changed since the last recompute (Shift+W forces one).\n" else
                                  "Shift+Q: nothing changed since the last recompute, so zoom only (Shift+W recomputes everything).\n")
                          }

                          if (shiftq_recompute) {
                            rm(list = ls(.store_fresh, all.names = TRUE), envir = .store_fresh)

                            ## Read in samples monolist and put chromatograms into chromatograms_updated

                                samples_monolist <- read.csv("samples_monolist.csv")
                                if ( length(samples_monolist_subset) > 0 ) {
                                    samples_monolist <- samples_monolist[samples_monolist_subset,]
                                }
                                chromatograms_updated <- dplyr::filter(chromatograms, path_to_cdf_csv %in% samples_monolist$path_to_cdf_csv)

                            ## Calculate baseline for each sample

                                baselined_chromatograms <- list()

                                ## Each sample's TIC rows taken from ONE split of the table, not one
                                ## dplyr::filter() of the whole table (~1.8M rows at 95 samples) per
                                ## sample. Same rows in the same order (test_shiftq_speed.R, 2026-10-21).
                                baseline_samples <- unique(chromatograms_updated$path_to_cdf_csv)
                                tic_rows_all <- chromatograms_updated[which(chromatograms_updated$ion == 0), , drop = FALSE]
                                tic_by_sample <- split(tic_rows_all, factor(tic_rows_all$path_to_cdf_csv, levels = unique(tic_rows_all$path_to_cdf_csv)))
                                for ( chrom in seq_along(baseline_samples) ) {

                                    tic <- tic_by_sample[[as.character(baseline_samples[chrom])]]
                                    if (is.null(tic)) tic <- tic_rows_all[0, , drop = FALSE]
                                    rownames(tic) <- NULL
                                    chromatogram <- tic

                                    prelim_baseline_window <- samples_monolist$baseline_window[match(chromatogram$path_to_cdf_csv[1], samples_monolist$path_to_cdf_csv)]

                                    ## A sample with fewer scans than baseline_window (a short TCD run,
                                    ## a narrow SIM method) gives floor(...) == 0; `1:0` then negative-indexes
                                    ## to numeric(0), min() returns Inf and data.frame() throws
                                    ## "arguments imply differing number of rows: 0, 1" with no hint that
                                    ## baseline_window is the knob. Clamp to at least one window, iterate with
                                    ## seq_len(), and clamp the slice so the last window cannot over-run.
                                    n_scans_in_tic <- length(tic$rt)
                                    n_prelim_baseline_windows <- max(1L, floor(n_scans_in_tic/prelim_baseline_window))
                                    if ( n_scans_in_tic < prelim_baseline_window ) {
                                        cat(paste0("Sample ", chromatogram$path_to_cdf_csv[1], " has ", n_scans_in_tic,
                                                   " scans, fewer than baseline_window (", prelim_baseline_window,
                                                   "); using one window over the whole trace. Lower baseline_window in samples_monolist.csv for a finer baseline.\n"))
                                    }
                                    prelim_baseline <- list()
                                    for ( i in seq_len(n_prelim_baseline_windows) ) {
                                        window_lo <- (prelim_baseline_window*(i-1))+1
                                        window_hi <- min(prelim_baseline_window*i, n_scans_in_tic)
                                        abundances_in_window <- tic$abundance[window_lo:window_hi]
                                        prelim_baseline[[i]] <- data.frame(
                                            rt = tic$rt[(which.min(abundances_in_window)+((i-1)*prelim_baseline_window))],
                                            min = min(abundances_in_window)
                                        )
                                    }
                                    prelim_baseline <- do.call(rbind, prelim_baseline)

                                    ## bound baseline to first and last point in chromatogram
                                    ## (index by position, not by rt value — original code accidentally
                                    ##  used the rt value as an integer index, which only happened to
                                    ##  land in-range when rt values were large; on TCD data where rt
                                    ##  starts near 0, tic$abundance[0] is length-0 and breaks data.frame())
                                    prelim_baseline <- rbind(
                                        data.frame(rt = min(tic$rt), min = tic$abundance[which.min(tic$rt)]),
                                        prelim_baseline,
                                        data.frame(rt = max(tic$rt), min = tic$abundance[which.max(tic$rt)])
                                    )

                                    tic$in_prelim_baseline <- FALSE
                                    tic$in_prelim_baseline[tic$rt %in% prelim_baseline$rt] <- TRUE

                                    y = prelim_baseline$min
                                    x = prelim_baseline$rt

                                    baseline2 <- data.frame(
                                        rt = tic$rt,
                                        y = approx(x, y, xout = tic$rt)$y
                                    )
                                    baseline2 <- baseline2[!is.na(baseline2$y),]
                                    tic <- tic[tic$rt %in% baseline2$rt,]
                                    tic$baseline <- baseline2$y

                                    baselined_chromatograms[[chrom]] <- data.frame(
                                        rt = tic$rt,
                                        abundance = tic$baseline,
                                        ion = "baseline",
                                        path_to_cdf_csv = tic$path_to_cdf_csv,
                                        rt_first_row_in_raw = tic$rt_first_row_in_raw,
                                        rt_last_row_in_raw = tic$rt_last_row_in_raw
                                    )
                                }

                                baselined_chromatograms <- do.call(rbind, baselined_chromatograms)
                                ## rbind onto chromatograms_updated (the SUBSET), not the unfiltered
                                ## `chromatograms` — rebuilding from the latter re-admitted every excluded
                                ## sample with no ion == "baseline" rows and rt_offset = NA, which killed
                                ## the render as soon as one of them landed on the visible page with a peak.
                                chromatograms_updated <- rbind(chromatograms_updated, baselined_chromatograms)

                            ## Add rt offset information for all chromatograms

                                chromatograms_updated$rt_offset <- samples_monolist$rt_offset[match(chromatograms_updated$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                chromatograms_updated$rt_rt_offset <- chromatograms_updated$rt + chromatograms_updated$rt_offset
                                chromatograms_updated <<- chromatograms_updated
                          }

                            ## Subset x_axis according to selection in chromatogram.
                            ## Recalculate skips this whole block once a view exists; before
                            ## the first one it falls through as a press with no brush.

                          view_ready <- length(x_axis_start) == 1 && length(x_axis_end) == 1 &&
                              is.finite(x_axis_start) && is.finite(x_axis_end)
                          brush <- if (keep_view) NULL else input$chromatogram_brush
                          if (keep_view && view_ready) {
                              cat(paste0("Recalculate: view kept, x-axis ", signif(x_axis_start, 8), " to ", signif(x_axis_end, 8), "\n"))
                          } else {

                                ## If null from initial start up, assign extreme values

                                    if ( is.null(x_axis_start_default) ) {
                                        x_axis_start_default <<- min(chromatograms$rt)
                                        cat(paste("x_axis_start_default is ", x_axis_start_default, "\n"))
                                    }

                                    if ( is.null(x_axis_end_default) ) {
                                        x_axis_end_default <<- max(chromatograms$rt)
                                        cat(paste("x_axis_end_default is ", x_axis_end_default, "\n"))
                                    }

                                ## If brush is null, assign default values to start and end

                                    if ( is.null(brush) ) {
                                        x_axis_start <<- x_axis_start_default
                                        x_axis_end <<- x_axis_end_default
                                        y_axis_start <<- 0
                                        y_axis_end <<- max(chromatograms$abundance)
                                        ## No brush means "reset", and that has to include y --
                                        ## otherwise a stale y window silently clips the full view.
                                        y_zoom <<- NULL
                                        y_zoom_key <<- NULL
                                    }

                                ## If brush is not null, assign brush values to start and end

                                    if ( !is.null(brush) ) {
                                        peak_points <<- isolate(brushed_chromatogram(chromatograms_updated, brush))

                                        ## An empty selection used to propagate straight through:
                                        ## min()/max() of nothing return Inf/-Inf, the axis limits
                                        ## become Inf..-Inf, the next render filters every row away,
                                        ## and ggplot dies with "Faceting variables must have at
                                        ## least one value". Because the limits are globals the app
                                        ## then stays wedged until it is restarted -- there is no
                                        ## brush you can draw to recover. Keep a usable view instead
                                        ## and say what happened.
                                        ## rt_rt_offset, not rt: the renderer filters the page on
                                        ## rt_rt_offset, so setting the limits from native rt is wrong
                                        ## by exactly rt_offset as soon as retention-time alignment is
                                        ## used -- and if the offset exceeds the brush width the window
                                        ## contains no points at all, so a brush drawn on a visible peak
                                        ## yields the "no points in range" placeholder.
                                        brush_x_col <- if ("rt_rt_offset" %in% names(peak_points)) "rt_rt_offset" else "rt"
                                        brush_ok <- nrow(peak_points) > 0 &&
                                            is.finite(min(peak_points[[brush_x_col]])) && is.finite(max(peak_points[[brush_x_col]]))

                                        if (brush_ok) {
                                            x_axis_start <<- min(peak_points[[brush_x_col]])
                                            x_axis_end <<- max(peak_points[[brush_x_col]])

                                            ## y follows the brush box itself, not the x-slice's
                                            ## full extent. Stored with a key naming the view it
                                            ## was taken in (display mode + y scaling), because
                                            ## the same number means different things in raw
                                            ## counts, square root, log and per-ion-normalised
                                            ## space -- applying it to the wrong one would clip
                                            ## the trace to nothing.
                                            ysel <- brush_y_range(input$chromatogram_brush)
                                            if (!is.null(ysel) && ysel[2] > ysel[1]) {
                                                y_zoom <<- ysel
                                                y_zoom_key <<- last_brush_key
                                                ## y_axis_start/end are the UNKEYED pair, and the renderer
                                                ## hands them to y_limits_for() as its FALLBACK. Writing a
                                                ## sqrt / log / per-ion-normalised / m-z range into them
                                                ## leaked that zoom into every other view: the key check
                                                ## correctly rejected the foreign y_zoom and then fell
                                                ## straight back to the same numbers, so a brush taken in
                                                ## all-ions/sqrt flattened the TIC through oob = squish,
                                                ## and one taken on the ion map gave the TIC a y-axis in
                                                ## m/z. They now only ever hold raw counts.
                                                if (is.null(y_zoom_key) || grepl("^tic/", y_zoom_key)) {
                                                    ## Bottom pinned at 0 -- see y_limits_for()'s
                                                    ## zero_floor. The y-brush contributes its TOP only.
                                                    y_axis_start <<- 0
                                                    y_axis_end   <<- ysel[2]
                                                } else {
                                                    cat(paste0("  y-zoom stored for ", y_zoom_key,
                                                               " only; the raw-count axis is left as it was.\n"))
                                                }
                                            } else {
                                                y_zoom <<- NULL
                                                y_zoom_key <<- NULL
                                                y_axis_start <<- 0
                                                y_axis_end <<- max(peak_points$abundance)
                                            }
                                        } else {
                                            notify("Selection held no data. View unchanged.")
                                            cat("  Re-brush on a drawn trace and press Shift+Q again, or press Shift+Q with no brush to reset.\n")
                                            if (is.null(x_axis_start) || length(x_axis_start) == 0 || !is.finite(x_axis_start)) x_axis_start <<- x_axis_start_default
                                            if (is.null(x_axis_end)   || length(x_axis_end)   == 0 || !is.finite(x_axis_end))   x_axis_end   <<- x_axis_end_default
                                            if (is.null(y_axis_start) || length(y_axis_start) == 0 || !is.finite(y_axis_start)) y_axis_start <<- 0
                                            if (is.null(y_axis_end)   || length(y_axis_end)   == 0 || !is.finite(y_axis_end))   y_axis_end   <<- max(chromatograms$abundance)
                                        }
                                    }

                                ## One line that distinguishes every way the zoom can fail: no
                                ## brush registered at all, a brush that produced no selection,
                                ## or limits that were set correctly and then ignored downstream.

                                    cat(paste0(
                                        if (keep_view) "Recalculate: no view yet, full run" else
                                            paste0("Shift+Q: brush ", if (is.null(brush)) "ABSENT" else "present"),
                                        "; x-axis now ", signif(x_axis_start, 8), " to ", signif(x_axis_end, 8),
                                        "; y-axis ",
                                        if (is.null(y_zoom)) "auto" else paste0(signif(y_zoom[1], 6), " to ", signif(y_zoom[2], 6)),
                                        "\n"))
                          }

                            ## No plot is assembled here. Shift+Q computes and writes; the
                            ## paginated output$chromatograms below draws. v5 moved the drawing
                            ## out but left v4's assembly in place, where it went on building a
                            ## chromatogram_plot that was immediately discarded -- two gsub passes
                            ## and a named-vector build over ONE element per row (~1-2M rows at
                            ## 152 samples), plus a full copy of the table attached to the dead
                            ## plot, on every single Shift+Q. Deleted 2026-10-01; nothing
                            ## downstream read chromatograms_updated_filtered or the plot.

                            ## Add peaks, if any
                          if (shiftq_recompute) {

                                peak_table <- read.csv("peaks_monolist.csv")

                                if (dim(peak_table)[1] > 0) {

                                    ## Filter out duplicate peaks and NA peaks

                                        ## Dedupe on the (start, end) PAIR. It used to run two
                                        ## independent passes -- one on peak_start, one on peak_end --
                                        ## so two genuinely distinct peaks that happened to share an
                                        ## end bound had the second silently dropped and committed to
                                        ## disk. Sorted by peak_start first, so "first wins" means the
                                        ## earlier peak rather than whichever row was appended first.
                                        peak_table <- peak_table[!is.na(peak_table$peak_start),]
                                        peak_table <- peak_table[order(peak_table$path_to_cdf_csv, peak_table$peak_start, peak_table$peak_end),]
                                        peak_table <- peak_table %>% group_by(path_to_cdf_csv) %>%
                                            ## + coelution_split: the two rows of a co-elution split share the
                                            ## hump's bounds by design and must not dedupe into one.
                                            mutate(duplicated = duplicated(paste(peak_start, peak_end,
                                                if ("coelution_split" %in% names(peak_table)) coelution_split else "", sep = "_")))
                                        peak_table <- as.data.frame(peak_table)
                                        peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
                                        peak_table <- peak_table[,!colnames(peak_table) == "duplicated"]

                                    ## Update with peak_number_within_sample

                                        peak_table_updated <- list()
                                        samples_with_peaks <- unique(as.character(peak_table$path_to_cdf_csv))
                                        for (sample_number in seq_along(samples_with_peaks)) {
                                          peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_with_peaks[sample_number],]
                                          peaks_in_this_sample <- peaks_in_this_sample[order(peaks_in_this_sample$peak_start,
                                              ## Full-length fallback: order() rejects a length-1 key (test_shiftq_unsplit_table.R).
                                              if ("coelution_split" %in% names(peaks_in_this_sample)) as.character(peaks_in_this_sample$coelution_split) else rep(NA_character_, nrow(peaks_in_this_sample)),
                                              na.last = TRUE),]
                                          peaks_in_this_sample$peak_number_within_sample <- seq(1,length(peaks_in_this_sample$path_to_cdf_csv),1)
                                          peak_table_updated[[sample_number]] <- peaks_in_this_sample
                                        }
                                        peak_table_updated <- do.call(rbind, peak_table_updated)
                                        peak_table <- peak_table_updated

                                    ## Modify peaks with RT offset

                                        # for (sample_number in 1:length(unique(samples_monolist$path_to_cdf_csv))) {

                                        #   peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_monolist$path_to_cdf_csv[sample_number],]

                                        #   rt_offsets <- samples_monolist$rt_offset[match(peaks_in_this_sample$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                        #   peak_start_rt_offsets <- peak_table$peak_start + peak_table$rt_offset
                                        #   peak_end_rt_offsets <- peak_table$peak_end + peak_table$rt_offset

                                        #   peak_table$rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- rt_offsets
                                        #   peak_table$peak_start_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_start_rt_offsets
                                        #   peak_table$peak_end_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_end_rt_offsets

                                        # }

                                        peak_table$rt_offset <- samples_monolist$rt_offset[match(peak_table$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                        peak_table$peak_start_rt_offset <- peak_table$peak_start + peak_table$rt_offset
                                        peak_table$peak_end_rt_offset <- peak_table$peak_end + peak_table$rt_offset
                                        peak_table$path_to_cdf_csv <- as.character(peak_table$path_to_cdf_csv)

                                    ## Update all peak areas in case baseline was adjusted

                                        if (isTRUE(integrate_peak_ions)) {
                                            if (is.null(peak_table$area_peak_ions)) peak_table$area_peak_ions <- NA_real_
                                            if (is.null(peak_table$peak_ions))      peak_table$peak_ions      <- NA_character_
                                        }

                                        ## Iterate over the UNIQUE sample list, and index into that
                                        ## same vector. The bound used to be length(unique(...)) while
                                        ## the index went into the raw column, so one duplicated row in
                                        ## samples_monolist.csv made the loop visit an early sample
                                        ## twice and never reach the last one -- whose areas then kept
                                        ## their stale values with no error anywhere.
                                        samples_for_area <- unique(as.character(samples_monolist$path_to_cdf_csv))
                                        ## Split the table by sample ONCE. Each peak used to subset the
                                        ## whole chromatograms_updated (TIC + baseline + every extracted
                                        ## ion, all samples) to its own sample: 714 full scans per
                                        ## Shift+Q, ~6.5 s on Skylar's folder vs 0.35 s split once.
                                        ## Same rows in the same order (test_shiftq_speed.R).
                                        chrom_by_sample <- split(chromatograms_updated,
                                                                 as.character(chromatograms_updated$path_to_cdf_csv))
                                        if (isTRUE(integrate_peak_ions)) {
                                            ion_area_cache_load()
                                            .ion_area_keys_used <<- character(0)
                                        }
                                        ## Phase 3: one spectrum per peak, parallel to peak_table's rows.
                                        ## Split rows are skipped (a whole-hump spectrum is neither compound).
                                        peak_spectra <- vector("list", nrow(peak_table))
                                        split_col <- if (is.null(peak_table$coelution_split)) rep(NA_character_, nrow(peak_table)) else
                                            as.character(peak_table$coelution_split)
                                        if (isTRUE(score_peak_confidence)) {
                                            spectrum_cache_load()
                                            .spectrum_keys_used <<- character(0)
                                        }
                                        gcms_prefetch_misses(peak_table,
                                            if (is.null(input$allion_threshold)) all_ion_threshold else input$allion_threshold,
                                            want_ion = isTRUE(integrate_peak_ions), want_spec = isTRUE(score_peak_confidence),
                                            split_col = split_col)
                                        for (sample_number in seq_along(samples_for_area)) {

                                          this_sample <- samples_for_area[sample_number]
                                          peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == this_sample,]

                                          ## Area = TIC sum over the peak window, minus the
                                          ## interpolated baseline over the same window.
                                          ##
                                          ## Fixed 2026-09-30. The baseline lives as ROWS tagged
                                          ## ion == "baseline" carrying their value in `abundance`
                                          ## (long format). It is NOT a `baseline` COLUMN -- that
                                          ## is the WIDE layout of modules/gcms.R, from which this
                                          ## block was copied verbatim. v4 changed the table to
                                          ## long form and adapted only the first sum; the second
                                          ## kept asking for $baseline, which returns NULL, and
                                          ## sum(NULL) is 0. So from v4 until now the subtraction
                                          ## was a silent no-op and every `area` ever written was
                                          ## a raw, un-baseline-corrected TIC sum. The old form
                                          ## also had no ion == 0 restriction on the baseline sum,
                                          ## so it would have summed across the TIC, the baseline
                                          ## and every extracted-ion row.
                                          ##
                                          ## seq_len, not 1:length(): a sample with no peaks yet
                                          ## gave 1:0 == c(1, 0), and the zero index handed
                                          ## dplyr::filter() a length-0 predicate, which errors.
                                          ## That is the normal state of a part-annotated folder,
                                          ## so Shift+Q died on the first un-annotated sample.
                                          areas <- numeric(nrow(peaks_in_this_sample))
                                          ## The sample's columns pulled out once; each peak's window is
                                          ## then a which() over them -- the same rows in the same order
                                          ## as the dplyr::filter() it replaces (which() drops NA like
                                          ## filter does), so the sums are bit-identical, without one
                                          ## filter() call per peak (2026-10-21).
                                          rows_for_sample <- chrom_by_sample[[this_sample]]
                                          if (is.null(rows_for_sample)) rows_for_sample <- chromatograms_updated[0, ]
                                          rt_s <- rows_for_sample$rt; ab_s <- rows_for_sample$abundance; ion_s <- rows_for_sample$ion
                                          for (peak in seq_len(nrow(peaks_in_this_sample))) {
                                            w_in <- which(rt_s >= peaks_in_this_sample$peak_start[peak] &
                                                          rt_s <= peaks_in_this_sample$peak_end[peak])
                                            ab_w <- ab_s[w_in]; ion_w <- ion_s[w_in]
                                            areas[peak] <-
                                                sum(ab_w[ion_w == 0]) -
                                                sum(ab_w[ion_w == "baseline"])
                                          }

                                          peak_table$area[peak_table$path_to_cdf_csv == this_sample] <- areas

                                          ## v6: the same area restricted to the ions the peak
                                          ## actually generates. The TIC area above charges a peak
                                          ## for everything under it -- column bleed, a tailing
                                          ## neighbour, a co-eluting second compound. Summing only
                                          ## the ions that rise across the peak window drops that.
                                          ## Reported alongside `area`, not instead of it, so
                                          ## numbers from v5 stay comparable.
                                          if (isTRUE(integrate_peak_ions) && nrow(peaks_in_this_sample) > 0) {
                                            ion_areas <- numeric(nrow(peaks_in_this_sample))
                                            ion_lists <- character(nrow(peaks_in_this_sample))
                                            for (peak in seq_len(nrow(peaks_in_this_sample))) {
                                              res <- peak_ion_area_memo(
                                                as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),
                                                peaks_in_this_sample$peak_start[peak],
                                                peaks_in_this_sample$peak_end[peak],
                                                ## The LIVE control, not the startup argument. These
                                                ## two diverge the moment the slider is touched, so
                                                ## the ion set you tune on screen was not the one
                                                ## that produced area_peak_ions.
                                                if (is.null(input$allion_threshold)) all_ion_threshold else input$allion_threshold
                                              )
                                              ion_areas[peak] <- res$area
                                              ion_lists[peak] <- res$ions
                                            }
                                            rows <- peak_table$path_to_cdf_csv == this_sample
                                            peak_table$area_peak_ions[rows] <- ion_areas
                                            peak_table$peak_ions[rows]      <- ion_lists
                                          }

                                          if (isTRUE(score_peak_confidence) && nrow(peaks_in_this_sample) > 0) {
                                            ## Loaded on the first memo miss only (see peak_spectrum_memo()).
                                            ## A store that fails, is empty or is single-detector gives NULL
                                            ## spectra; remembered so it is tried once per sample.
                                            st_cache <- new.env(parent = emptyenv())
                                            get_store <- function() {
                                              if (is.null(st_cache$v)) {
                                                store <- try(load_allions(this_sample), silent = TRUE)
                                                st_cache$v <- if (!inherits(store, "try-error") && !is.null(store) && nrow(store) > 0 &&
                                                                 !all(store$mz == 0)) {
                                                  store_rts <- unique(store$rt)      # rt-keyed, so already sorted
                                                  if (is.unsorted(store_rts)) store_rts <- sort(store_rts)
                                                  list(store = store, rts = store_rts)
                                                } else FALSE
                                              }
                                              if (isFALSE(st_cache$v)) NULL else st_cache$v
                                            }
                                            for (r in which(peak_table$path_to_cdf_csv == this_sample)) {
                                              if (!is.na(split_col[r]) && nzchar(split_col[r])) next
                                              sp <- peak_spectrum_memo(this_sample, peak_table$peak_start[r],
                                                                       peak_table$peak_end[r], get_store)
                                              if (!is.null(sp)) peak_spectra[r] <- list(sp)
                                            }
                                          }

                                          # peaks_in_this_sample$area <- areas
                                          # peak_table_updated[[sample_number]] <- peaks_in_this_sample
                                        }
                                        # peak_table_updated <- do.call(rbind, peak_table_updated)
                                        # peak_table <- peak_table_updated
                                        if (isTRUE(integrate_peak_ions)) ion_area_cache_save(.ion_area_keys_used)
                                        if (isTRUE(score_peak_confidence)) {
                                            got <- !vapply(peak_spectra, is.null, logical(1))
                                            agree <- gcms_spectrum_agreement(
                                                lapply(peak_spectra, function(x) if (is.null(x)) NULL else x$spectrum),
                                                peak_table$peak_ID)
                                            peak_table$spectrum_agreement <- round(agree, 3)
                                            peak_table$spectrum_signal_frac <- ifelse(got,
                                                round(vapply(peak_spectra, function(x) if (is.null(x)) NA_real_ else x$signal_frac, numeric(1)), 3),
                                                NA_real_)
                                            spectrum_cache_save(.spectrum_keys_used)
                                        }

                                    ## v7: refresh existing co-elution splits. Runs AFTER the area
                                    ## recompute above, which has just given every split pair plain
                                    ## drop-line areas; each pair is redone from the settings stored
                                    ## on its rows. Shift+Q never CREATES a split - that is the
                                    ## Split button.
                                        refreshed <- gcms_refresh_splits(peak_table, coelution_traces,
                                                                         crosstalk = coelution_crosstalk)
                                        peak_table <- refreshed$peak_table
                                        notify_log(refreshed$log)
                                        write_coelution_profiles(refreshed$profiles)

                                    ## Write out peaks now with assigned peak_number_within_sample and RT offset, update the peak_table in ui

                                        write.table(peak_table, file = "peaks_monolist.csv", col.names = TRUE, sep = ",", row.names = FALSE)

                                        render_peak_table()

                                }
                                ## After the write above, so the table this recompute wrote counts as seen.
                                .shiftq_fp <<- shiftq_fingerprint()
                          }

                            ## The x-window globals must be usable before the paginated renderer
                            ## filters on them. Hoisted out of the peak block 2026-10-01: it used
                            ## to sit inside `if (dim(peak_table)[1] > 0)`, so on a folder with no
                            ## peaks yet the safety net did not run at all.
                                if (length(x_axis_start) == 0) {x_axis_start <<- min(chromatograms$rt)}
                                if (length(x_axis_end) == 0) {x_axis_end <<- max(chromatograms$rt)}

                                ## Drawing moved out of Shift+Q: the chromatogram is now drawn by
                                ## the reactive output$chromatograms defined just below this handler,
                                ## which renders only the current PAGE of samples (see pagination
                                ## block). Here we only signal that fresh data is ready.
                                redraw_trigger(isolate(redraw_trigger()) + 1)

                                notify(if (keep_view) "Recalculated (view kept)." else "Chromatogram updated.")
                    }
                    observeEvent(input$keypress, {
                        if( input$keypress == 81 || input$keypress == 87 ) { # Shift+Q; Shift+W forces the full recompute
                            shiftq_run(force = input$keypress == 87)
                        }
                    })
                    observeEvent(input$recalc_btn, shiftq_run(keep_view = TRUE))

                ## ------------------------------------------------------------------
                ## Paginated chromatogram renderer + page controls
                ##
                ## A single faceted plot (one panel per sample) is what crashed v4 at
                ## 100-200 samples: render cost scales with the number of panels, and a
                ## ~20,000px image exhausts the graphics device (plot greys out, app
                ## goes unresponsive). Here we draw only samples_per_page samples at a
                ## time. Shift+Q recomputes the data and bumps redraw_trigger; Prev/Next
                ## just change current_page. Paging does NOT re-run the (expensive)
                ## baseline / area computation, so moving between views stays snappy.
                ## ------------------------------------------------------------------

                    ## Page size -- a live sidebar control; samples_per_page is only its
                    ## starting value. Held in a reactiveVal rather than read straight off
                    ## the box, so an emptied or out-of-range entry is ignored.
                        page_size_val <- reactiveVal(samples_per_page)
                        page_size <- function() page_size_val()
                        clamp_page_size <- function(v) {
                            v <- suppressWarnings(as.integer(v))
                            if (length(v) != 1 || is.na(v) || v < 1) return(NULL)
                            min(v, 500L)
                        }
                        ## Keep the first sample of the current page on screen when the
                        ## page size changes (page_bounds() clamps the result).
                        repage <- function(old_size, new_size) {
                            first <- (current_page() - 1) * old_size + 1
                            current_page(floor((first - 1) / new_size) + 1)
                        }
                        observeEvent(input$samples_per_page_ui, {
                            v <- clamp_page_size(input$samples_per_page_ui)
                            if (is.null(v)) return()
                            if (v == page_size_val()) return()
                            repage(page_size_val(), v)
                            page_size_val(v)
                        })
                        observe({
                            n <- page_size()
                            updateActionButton(session, "prev_page", label = paste0("\u25c0 Prev ", n))
                            updateActionButton(session, "next_page", label = paste0("Next ", n, " \u25b6"))
                        })
                        ## Live plot height: 200 px of axes and margins plus one panel per
                        ## sample actually on this page, so a short last page is not stretched.
                        ## Panel height multiplier (sidebar slider). Debounced: one redraw at
                        ## 95 samples takes seconds, and a slider sends values while dragged.
                        panel_height_mult <- debounce(reactive({
                            v <- suppressWarnings(as.numeric(input$panel_height_mult))
                            if (length(v) != 1 || !is.finite(v) || v <= 0) 1 else min(v, 4)
                        }), 700)
                        chromatogram_height <- function() {
                            redraw_trigger()
                            b <- page_bounds()
                            n_here <- if (!is.null(b) && b$n > 0) b$hi - b$lo + 1 else
                                min(page_size(), max(1, length(unique(chromatograms$path_to_cdf_csv))))
                            200 + round(panel_px * panel_height_mult()) * n_here
                        }

                    ## Which samples exist, and the current page's slice
                        page_bounds <- function() {
                            if (is.null(chromatograms_updated)) return(NULL)
                            all_samples <- sort(unique(as.character(chromatograms_updated$path_to_cdf_csv)))
                            n <- length(all_samples)
                            ## No samples at all: return an EMPTY page rather than
                            ## all_samples[1:0], which gives NA_character_ -- a phantom sample
                            ## that the renderer then tried to facet on, under an indicator
                            ## reading "samples 1-0 of 0".
                            if (n == 0) {
                                return(list(all_samples = character(0), n = 0L, n_pages = 1L,
                                            pg = 1L, lo = 0L, hi = 0L, page_samples = character(0)))
                            }
                            n_pages <- max(1, ceiling(n / page_size()))
                            pg <- min(max(1, current_page()), n_pages)
                            lo <- (pg - 1) * page_size() + 1
                            hi <- min(pg * page_size(), n)
                            list(all_samples = all_samples, n = n, n_pages = n_pages,
                                 pg = pg, lo = lo, hi = hi,
                                 page_samples = all_samples[seq(lo, hi)])
                        }

                    ## v7 co-elution split: the Split / Unsplit buttons. Both work on the
                    ## table and chromatograms the last Shift+Q produced, then refresh every
                    ## split (which also rebuilds the stacked-band file) and redraw.
                        coelution_traces <- function(smp, st, en, ion_sets) {
                            if (!(smp %in% chromatograms_updated$path_to_cdf_csv)) return("not_loaded")
                            gcms_split_traces(
                                chromatograms_updated[chromatograms_updated$path_to_cdf_csv == smp, ],
                                load_allions(smp), st, en, ion_sets)
                        }
                        write_coelution_profiles <- function(profiles) {
                            if (!is.null(profiles)) {
                                utils::write.csv(profiles, "coelution_profiles.csv", row.names = FALSE)
                            } else if (file.exists("coelution_profiles.csv")) {
                                file.remove("coelution_profiles.csv")
                            }
                        }
                        coelution_finish <- function(pt) {
                            refreshed <- gcms_refresh_splits(pt, coelution_traces, crosstalk = coelution_crosstalk)
                            notify_log(refreshed$log)
                            write_coelution_profiles(refreshed$profiles)
                            write.table(refreshed$peak_table, file = "peaks_monolist.csv", col.names = TRUE, sep = ",", row.names = FALSE)
                            render_peak_table()
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                        }
                        coelution_table <- function() {
                            if (is.null(chromatograms_updated) || !file.exists("peaks_monolist.csv")) {
                                notify("Split: press Shift+Q first.")
                                return(NULL)
                            }
                            pt <- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                            if (nrow(pt) == 0 || !all(c("peak_start_rt_offset", "peak_number_within_sample") %in% names(pt))) {
                                notify("Split: no numbered peaks yet. Integrate the hump and press Shift+Q first.")
                                return(NULL)
                            }
                            pt
                        }
                        observeEvent(input$split_btn, {
                            undo_begin(paste0("split peak ", input$split_peak)); on.exit(undo_end(), add = TRUE)
                            pt <- coelution_table(); if (is.null(pt)) return()
                            ion_sets <- gcms_split_ion_sets(input$split_ions_a, input$split_ions_b, input$split_ions_c)
                            n <- input$split_peak
                            if (is.character(ion_sets)) { notify(ion_sets); return() }
                            wnd <- gcms_window_for_peak(pt, n)
                            if (is.null(wnd)) { notify(paste0("Split: no sample has a peak number ", n, ".")); return() }
                            cat(sprintf("Split: peak %d is at %.2f-%.2f (aligned rt).\n", as.integer(round(n)), wnd[1], wnd[2]))
                            res <- gcms_apply_coelution_split(pt, coelution_traces, wnd, ion_sets,
                                                              crosstalk = coelution_crosstalk,
                                                              k_mode = if (is.null(input$split_k_mode)) "sample" else input$split_k_mode)
                            notify_log(res$log)
                            coelution_finish(res$peak_table)
                        })
                        observeEvent(input$unsplit_btn, {
                            undo_begin(paste0("unsplit peak ", input$split_peak)); on.exit(undo_end(), add = TRUE)
                            pt <- coelution_table(); if (is.null(pt)) return()
                            res <- gcms_unsplit_peak(pt, input$split_peak)
                            notify_log(res$log)
                            coelution_finish(res$peak_table)
                        })

                    ## Page navigation buttons
                        observeEvent(input$prev_page, {
                            current_page(max(1, current_page() - 1))
                        })
                        observeEvent(input$next_page, {
                            b <- page_bounds()
                            if (is.null(b)) { notify("No chromatograms loaded yet. Press Shift+Q first."); return() }
                            current_page(min(b$n_pages, current_page() + 1))
                        })

                    ## Earlier / Later: the same-width window half a width before or after.
                    ## Moves the view only (x limits); a y-zoom is kept, as with a brush zoom.
                        slide_view <- function(direction) {
                            if (is.null(chromatograms_updated)) { notify("No chromatograms loaded yet. Press Shift+Q first."); return() }
                            rng <- suppressWarnings(range(chromatograms_updated$rt_rt_offset, na.rm = TRUE))
                            nw <- gcms_slide_window(x_axis_start, x_axis_end, direction, rng[1], rng[2])
                            if (is.null(nw)) { notify("Zoom in first (brush + Shift+Q), then slide."); return() }
                            if (isTRUE(all.equal(nw, c(x_axis_start, x_axis_end)))) {
                                notify(if (isTRUE(all.equal(nw, as.numeric(rng)))) "Showing the whole run: zoom in first, then slide."
                                       else if (direction > 0) "Already at the end of the run." else "Already at the start of the run.")
                                return()
                            }
                            x_axis_start <<- nw[1]; x_axis_end <<- nw[2]
                            cat(sprintf("Slide: x-axis now %s to %s\n", signif(nw[1], 8), signif(nw[2], 8)))
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                        }
                        observeEvent(input$pan_back, slide_view(-1))
                        observeEvent(input$pan_fwd,  slide_view(+1))

                    ## View history (begin) -- test_view_history.R extracts this block.
                    ## Previous view / Subsequent view. Every view change (a Shift+Q zoom
                    ## or reset, Earlier / Later) ends in a redraw, so the history is
                    ## recorded there: on each redraw, a view that differs from the one at
                    ## the current position becomes a new entry. Stepping sets the position
                    ## FIRST and then restores that view, so its own redraw finds them equal
                    ## and records nothing. Plain environment, not `<<-` globals.
                        view_hist <- new.env(parent = emptyenv())
                        view_hist$h <- NULL
                        current_view <- function() list(xs = x_axis_start, xe = x_axis_end,
                                                        ys = y_axis_start, ye = y_axis_end,
                                                        yz = y_zoom, yk = y_zoom_key)
                        ## ignoreInit: the view globals live in analyzeGCMSdata7()'s frame, which
                        ## every browser session of one launch shares, so at startup they can hold
                        ## an EARLIER session's zoom. Recording it would offer a view this session
                        ## never showed (seen on host1 after a page reload).
                        observeEvent(redraw_trigger(), {
                            view_hist$h <- gcms_view_history_record(view_hist$h, current_view())
                        }, ignoreInit = TRUE)
                        step_view <- function(direction) {
                            nh <- gcms_view_history_step(view_hist$h, direction)
                            if (is.null(nh)) {
                                notify(if (direction < 0) "No earlier view to go back to." else "No later view to go forward to.")
                                return()
                            }
                            v <- nh$views[[nh$pos]]
                            view_hist$h <- nh
                            x_axis_start <<- v$xs; x_axis_end <<- v$xe
                            y_axis_start <<- v$ys; y_axis_end <<- v$ye
                            y_zoom <<- v$yz; y_zoom_key <<- v$yk
                            cat(sprintf("View %d of %d: x-axis %s to %s\n", nh$pos, length(nh$views),
                                        signif(v$xs, 8), signif(v$xe, 8)))
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                        }
                        observeEvent(input$view_back, step_view(-1))
                        observeEvent(input$view_fwd,  step_view(+1))
                    ## View history (end)

                    ## Zoom to peaks (begin) -- test_zoom_peaks.R extracts this block.
                    ## Zoom to peaks by number (the box and button in the button row).
                    ## Sets the same globals a Shift+Q brush zoom does; y is fitted to
                    ## the tallest TIC in the window, any y-zoom cleared.
                        observeEvent(input$zoom_peaks_btn, {
                            if (is.null(chromatograms_updated) || !file.exists("peaks_monolist.csv")) {
                                notify("No chromatograms loaded yet. Press Shift+Q first."); return()
                            }
                            rng <- gcms_parse_peak_range(input$zoom_peaks)
                            if (is.null(rng)) { notify("Zoom to peaks: type a peak number, or a range like 12-15."); return() }
                            pt <- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                            if (!all(c("peak_number_within_sample", "peak_start_rt_offset", "peak_end_rt_offset") %in% names(pt))) {
                                notify("Zoom to peaks: no numbered peaks yet. Press Shift+Q first."); return()
                            }
                            win <- gcms_peak_range_window(pt, rng[1], rng[2])
                            if (is.null(win)) {
                                notify(paste0("Zoom to peaks: no sample has peak ",
                                              if (rng[1] == rng[2]) rng[1] else paste0(rng[1], "-", rng[2]), "."))
                                return()
                            }
                            tic <- chromatograms_updated[!is.na(chromatograms_updated$ion) & chromatograms_updated$ion == 0 &
                                                         chromatograms_updated$rt_rt_offset >= win[1] &
                                                         chromatograms_updated$rt_rt_offset <= win[2], , drop = FALSE]
                            top <- suppressWarnings(max(tic$abundance, na.rm = TRUE))
                            x_axis_start <<- win[1]; x_axis_end <<- win[2]
                            y_zoom <<- NULL; y_zoom_key <<- NULL
                            y_axis_start <<- 0
                            y_axis_end <<- if (is.finite(top) && top > 0) top * 1.05 else max(chromatograms$abundance)
                            cat(sprintf("Zoom to peaks %d-%d: x-axis %s to %s\n", rng[1], rng[2], signif(win[1], 8), signif(win[2], 8)))
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                        })
                    ## Zoom to peaks (end)

                    ## Align: move one sample along rt. Same tables Shift+Q would rebuild,
                    ## updated in place (gcms_shift_sample_offset()), so it redraws at once.
                    ## Driven by the arrow pair beside each panel (input$align_nudge, sent by
                    ## the gcms_align_rows script), by the sidebar's Align step.
                        nudge_sample <- function(smp, delta, absolute = FALSE) {
                            ## A run of arrow presses on one sample is one undo step.
                            undo_begin(paste0("align ", basename(as.character(smp))), coalesce = paste0("nudge:", smp))
                            on.exit(undo_end(), add = TRUE)
                            fp_was_current <- !is.null(.shiftq_fp) && identical(shiftq_fingerprint(), .shiftq_fp)
                            if (is.null(chromatograms_updated)) { notify("No chromatograms loaded yet. Press Shift+Q first."); return() }
                            if (is.null(smp) || !nzchar(smp)) { notify("No sample to align."); return() }
                            if (isTRUE(table_dirty())) { notify("Save the peak table first (Shift+Z): aligning redraws it."); return() }
                            sm <- read.csv("samples_monolist.csv", stringsAsFactors = FALSE)
                            if (absolute) {
                                cur <- suppressWarnings(as.numeric(sm$rt_offset[match(smp, sm$path_to_cdf_csv)]))
                                delta <- -(if (is.finite(cur)) cur else 0)
                            }
                            pk <- if (file.exists("peaks_monolist.csv")) read.csv("peaks_monolist.csv", stringsAsFactors = FALSE) else NULL
                            pr <- if (file.exists("coelution_profiles.csv")) try(read.csv("coelution_profiles.csv", stringsAsFactors = FALSE), silent = TRUE) else NULL
                            if (inherits(pr, "try-error")) pr <- NULL
                            res <- gcms_shift_sample_offset(smp, delta, sm, pk, chromatograms_updated, pr)
                            if (is.null(res)) { notify("That sample is not in samples_monolist.csv."); return() }
                            write.table(res$samples, file = "samples_monolist.csv", row.names = FALSE, sep = ",")
                            if (!is.null(res$peaks)) write.table(res$peaks, file = "peaks_monolist.csv", col.names = TRUE, sep = ",", row.names = FALSE)
                            if (!is.null(res$profiles)) utils::write.csv(res$profiles, "coelution_profiles.csv", row.names = FALSE)
                            chromatograms_updated <<- res$chrom
                            ## The plot redraws now; the peak table (thousands of rows, ~2 s
                            ## for the browser to rebuild) once the arrows have been left alone
                            ## for a moment -- see align_table_due below. Shift+Z is safe
                            ## meanwhile: it copies only peak_ID onto the file as it is on disk,
                            ## and a nudge never changes peak_ID or the row order.
                            ## Counted from when this redraw has been SENT (onFlushed), not from
                            ## the press: the redraw itself takes ~2.5 s, so a timer started at
                            ## the press ran out before the plot even appeared.
                            session$onFlushed(function() align_table_due(isolate(align_table_due()) + 1), once = TRUE)
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                            ## Align updates every table Shift+Q would rebuild, in place, so a
                            ## recompute that was current stays current: the next Shift+Q is a zoom.
                            if (fp_was_current) .shiftq_fp <<- shiftq_fingerprint()
                            cat(sprintf("Align: %s rt_offset now %g\n", smp, res$offset))
                        }
                        ## Debounced: a run of nudges refreshes the table once, 1.5 s after the
                        ## last. Skipped if the table now holds unsaved peak_ID edits (a
                        ## re-render would discard them); the next Shift+Q refreshes it then.
                        align_table_due <- reactiveVal(0)
                        align_table_quiet <- debounce(reactive(align_table_due()), 1500)
                        observeEvent(align_table_quiet(), {
                            if (align_table_quiet() > 0 && !isTRUE(isolate(table_dirty()))) render_peak_table()
                        }, ignoreInit = TRUE)
                        ## Align on a peak / Reset all offsets: many samples in one pass
                        ## (gcms_set_sample_offsets()), then the same redraw, debounced table
                        ## refresh and fingerprint rule as a single nudge.
                        set_offsets <- function(make_offsets, what) {
                            undo_begin(tolower(what)); on.exit(undo_end(), add = TRUE)
                            fp_was_current <- !is.null(.shiftq_fp) && identical(shiftq_fingerprint(), .shiftq_fp)
                            if (is.null(chromatograms_updated)) { notify("No chromatograms loaded yet. Press Shift+Q first."); return() }
                            if (isTRUE(table_dirty())) { notify("Save the peak table first (Shift+Z): aligning redraws it."); return() }
                            sm <- read.csv("samples_monolist.csv", stringsAsFactors = FALSE)
                            pk <- if (file.exists("peaks_monolist.csv")) read.csv("peaks_monolist.csv", stringsAsFactors = FALSE) else NULL
                            plan <- make_offsets(sm, pk)
                            if (is.null(plan)) return()
                            pr <- if (file.exists("coelution_profiles.csv")) try(read.csv("coelution_profiles.csv", stringsAsFactors = FALSE), silent = TRUE) else NULL
                            if (inherits(pr, "try-error")) pr <- NULL
                            res <- gcms_set_sample_offsets(plan$offsets, sm, pk, chromatograms_updated, pr)
                            if (is.null(res)) { notify(paste0(what, ": nothing to change.")); return() }
                            write.table(res$samples, file = "samples_monolist.csv", row.names = FALSE, sep = ",")
                            if (!is.null(res$peaks)) write.table(res$peaks, file = "peaks_monolist.csv", col.names = TRUE, sep = ",", row.names = FALSE)
                            if (!is.null(res$profiles)) utils::write.csv(res$profiles, "coelution_profiles.csv", row.names = FALSE)
                            chromatograms_updated <<- res$chrom
                            session$onFlushed(function() align_table_due(isolate(align_table_due()) + 1), once = TRUE)
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                            if (fp_was_current) .shiftq_fp <<- shiftq_fingerprint()
                            msg <- paste0(what, ": ", length(res$changed), " sample", if (length(res$changed) == 1) "" else "s", " moved.")
                            if (length(plan$missing)) {
                                nm <- basename(plan$missing)
                                msg <- paste0(msg, " ", length(nm), " without that peak, left as they were: ",
                                              paste(utils::head(nm, 5), collapse = ", "), if (length(nm) > 5) ", ..." else "", ".")
                            }
                            notify(msg)
                        }
                        observeEvent(input$alignall_btn, {
                            n <- input$alignall_peak
                            set_offsets(function(sm, pk) {
                                if (is.null(pk) || !("peak_number_within_sample" %in% names(pk))) {
                                    notify("Align: no numbered peaks yet. Press Shift+Q first."); return(NULL)
                                }
                                plan <- gcms_align_offsets_for_peak(pk, sm, n)
                                if (is.null(plan)) { notify(paste0("Align: no sample has a peak number ", n, ".")); return(NULL) }
                                cat(sprintf("Align: peak %d, window centres onto %.4f (raw rt median).\n", as.integer(round(n)), plan$target))
                                plan
                            }, paste0("Aligned on peak ", n))
                        })
                        observeEvent(input$align_reset_btn, {
                            set_offsets(function(sm, pk) list(offsets = stats::setNames(rep(0, nrow(sm)), as.character(sm$path_to_cdf_csv))),
                                        "Reset all offsets")
                        })
                        observeEvent(input$align_nudge, {
                            a <- input$align_nudge
                            smp <- if (is.list(a)) a$sample else NULL
                            dir <- suppressWarnings(as.numeric(if (is.list(a)) a$dir else NA))
                            if (is.null(smp) || !is.finite(dir)) return()
                            st <- suppressWarnings(as.numeric(input$align_step)); if (!length(st) || !is.finite(st) || st <= 0) st <- 1
                            nudge_sample(as.character(smp), sign(dir) * st)
                        })

                    ## Undo (begin) -- test_undo.R extracts this block.
                    ## The Undo button. Each edit handler calls undo_begin() before it
                    ## touches anything and undo_end() on exit (on.exit, so early returns
                    ## and errors are covered); an edit that changed no file leaves no step.
                    ## Recalculations (Shift+Q/W, Recalculate) are not steps. Per browser
                    ## session, last 20 steps, gone when the app closes.
                        undo_hist <- new.env(parent = emptyenv())
                        undo_hist$stack <- list()
                        undo_hist$pending <- NULL
                        undo_label_refresh <- function() {
                            st <- undo_hist$stack
                            updateActionButton(session, "undo_btn",
                                label = if (length(st)) paste0("Undo: ", st[[length(st)]]$label) else "Undo")
                        }
                        undo_begin <- function(label, coalesce = NULL) {
                            fp <- shiftq_fingerprint()
                            undo_hist$pending <- list(label = label, coalesce = coalesce, files = gcms_undo_capture(),
                                                      fp = fp, fp_current = !is.null(.shiftq_fp) && identical(fp, .shiftq_fp))
                        }
                        undo_end <- function() {
                            step <- undo_hist$pending
                            undo_hist$pending <- NULL
                            if (is.null(step)) return(invisible(NULL))
                            undo_hist$stack <- gcms_undo_push(undo_hist$stack, step, gcms_undo_capture())
                            undo_label_refresh()
                        }
                        observeEvent(input$undo_btn, {
                            st <- undo_hist$stack
                            if (!length(st)) { notify("Nothing to undo."); return() }
                            if (isTRUE(table_dirty())) { notify("Save the peak table first (Shift+Z): undoing redraws it."); return() }
                            step <- st[[length(st)]]
                            ## The loaded chromatograms carry each sample's rt_offset; put back
                            ## the offsets of the step being restored, as Align would -- worked
                            ## out now, applied only once the files are back.
                            new_chrom <- NULL
                            old_sm <- step$files[["samples_monolist.csv"]]
                            if (!is.null(chromatograms_updated) && !is.null(old_sm) && file.exists("samples_monolist.csv")) {
                                then <- try(utils::read.csv(text = rawToChar(old_sm), stringsAsFactors = FALSE), silent = TRUE)
                                cur  <- try(utils::read.csv("samples_monolist.csv", stringsAsFactors = FALSE), silent = TRUE)
                                if (!inherits(then, "try-error") && !inherits(cur, "try-error") &&
                                    all(c("path_to_cdf_csv", "rt_offset") %in% names(then))) {
                                    off <- suppressWarnings(as.numeric(then$rt_offset)); off[!is.finite(off)] <- 0
                                    res <- gcms_set_sample_offsets(stats::setNames(off, as.character(then$path_to_cdf_csv)),
                                                                   cur, NULL, chromatograms_updated, NULL)
                                    if (!is.null(res)) new_chrom <- res$chrom
                                }
                            }
                            ok <- try(gcms_undo_restore(step$files), silent = TRUE)
                            if (inherits(ok, "try-error")) { notify(paste0("Undo failed: ", conditionMessage(attr(ok, "condition")))); return() }
                            if (!is.null(new_chrom)) chromatograms_updated <<- new_chrom
                            undo_hist$stack <- st[-length(st)]
                            ## Restored to a state the last recompute had already seen: the
                            ## next Shift+Q is a zoom again, as after an Align.
                            if (isTRUE(step$fp_current) && identical(shiftq_fingerprint(), step$fp)) .shiftq_fp <<- step$fp
                            render_peak_table()
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                            undo_label_refresh()
                            notify(paste0("Undid: ", step$label, "."))
                        })
                    ## Undo (end)

                    ## Right-click menu (begin) -- test_context_menu.R extracts this block.
                    ## Right-click on the chromatograms opens a small menu (script beside the
                    ## keypress one). Picking an item re-sends the right-click's spot as
                    ## input$chromatogram_click (NULL when it fell outside every panel), then
                    ## sends input$chromatogram_ctx with the item's action:
                    ##   remove  the peak under that spot, in that sample, goes;
                    ##   split   a dialog splits that ONE peak in that ONE sample by ions
                    ##           (the sidebar's Split does every sample), or re-splits /
                    ##           unsplits it if it is already a split pair.
                    ## A split pair is one hump on screen, so both its rows act together.
                        ctx_target <- new.env(parent = emptyenv())   # the peak the Split dialog is for
                        ctx_chr <- function(x) if (length(x) != 1 || is.na(x)) "" else as.character(x)
                        ## The peak row under the right-click: list(pt, i, smp, what), or NULL
                        ## having said why.
                        ctx_hit <- function() {
                            if (is.null(chromatograms_updated) || !file.exists("peaks_monolist.csv")) {
                                notify("No chromatograms loaded yet. Press Shift+Q first."); return(NULL)
                            }
                            click <- input$chromatogram_click
                            smp <- clicked_sample(click, strict = TRUE)
                            rt  <- clicked_rt(click)
                            if (is.null(smp) || is.null(rt)) {
                                notify("Right-click inside a sample's panel, on a peak."); return(NULL)
                            }
                            ## The x axis is aligned rt; the peak list holds raw rt.
                            off <- suppressWarnings(as.numeric(chromatograms_updated$rt_offset[
                                       match(smp, as.character(chromatograms_updated$path_to_cdf_csv))]))
                            if (length(off) != 1 || !is.finite(off)) {
                                notify("That sample's alignment offset is unknown. Press Shift+Q first."); return(NULL)
                            }
                            pt <- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                            i <- gcms_peak_at(pt, smp, rt - off)
                            if (is.na(i)) { notify("No peak there."); return(NULL) }
                            num <- if ("peak_number_within_sample" %in% names(pt))
                                       suppressWarnings(as.numeric(pt$peak_number_within_sample[i])) else NA
                            list(pt = pt, i = i, smp = smp,
                                 what = if (is.finite(num)) paste0("peak ", num) else "the peak")
                        }
                        ## Row i, or every row of the split it belongs to.
                        ctx_rows <- function(pt, i) {
                            if (!all(c("coelution_split", "split_id") %in% names(pt)) ||
                                !(ctx_chr(pt$coelution_split[i]) %in% .gcms_split_labels) || !nzchar(ctx_chr(pt$split_id[i]))) return(i)
                            rows <- which(as.character(pt$path_to_cdf_csv) == as.character(pt$path_to_cdf_csv[i]) &
                                          as.character(pt$split_id) %in% ctx_chr(pt$split_id[i]) &
                                          pt$coelution_split %in% .gcms_split_labels)
                            if (length(rows) >= 2 && length(rows) == .gcms_split_size(pt, rows[1])) rows else i
                        }
                        ## The dialog's peak, found again in the table as it is NOW (by sample and
                        ## bounds, or by split id): rows, or NULL if it has gone or changed.
                        ctx_find <- function(pt) {
                            t <- ctx_target$t
                            if (is.null(t)) return(NULL)
                            smp <- as.character(pt$path_to_cdf_csv)
                            if (nzchar(t$split_id)) {
                                if (!all(c("coelution_split", "split_id") %in% names(pt))) return(NULL)
                                rows <- which(smp == t$smp & as.character(pt$split_id) %in% t$split_id &
                                              pt$coelution_split %in% .gcms_split_labels)
                                return(if (length(rows) >= 2 && length(rows) == .gcms_split_size(pt, rows[1])) rows else NULL)
                            }
                            paired <- if ("coelution_split" %in% names(pt)) pt$coelution_split %in% .gcms_split_labels else FALSE
                            rows <- which(smp == t$smp & !paired & abs(pt$peak_start - t$start) < 1e-9 & abs(pt$peak_end - t$end) < 1e-9)
                            if (length(rows) >= 1) rows[1] else NULL
                        }
                        ctx_remove <- function() {
                            undo_begin("remove peak (right-click)"); on.exit(undo_end(), add = TRUE)
                            h <- ctx_hit(); if (is.null(h)) return()
                            rows <- ctx_rows(h$pt, h$i)
                            write.table(
                              x = h$pt[-rows, , drop = FALSE],
                              file = "peaks_monolist.csv",
                              append = FALSE,
                              row.names = FALSE,
                              col.names = TRUE,
                              sep = ","
                            )
                            notify(paste0("Removed ", h$what, if (length(rows) == 2) " (both parts of its split)" else if (length(rows) > 2) " (all three parts of its split)" else "",
                                          " from ", basename(h$smp), "."))
                            ## Unsaved table edits stay on screen; Shift+Z then refuses safely
                            ## (the row count moved) rather than write the peak back.
                            if (!isTRUE(table_dirty())) render_peak_table()
                            redraw_trigger(isolate(redraw_trigger()) + 1)
                        }
                        ctx_split_dialog <- function() {
                            h <- ctx_hit(); if (is.null(h)) return()
                            if (!all(c("peak_start_rt_offset", "peak_number_within_sample") %in% names(h$pt))) {
                                notify("Split: no numbered peaks yet. Press Shift+Q first."); return()
                            }
                            rows <- ctx_rows(h$pt, h$i)
                            pair <- length(rows) >= 2
                            r1 <- h$pt[rows[1], , drop = FALSE]
                            ctx_target$t <- list(smp = h$smp, start = r1$peak_start, end = max(h$pt$peak_end[rows]),
                                                 split_id = if (pair) ctx_chr(r1$split_id) else "")
                            ia <- if (pair) ctx_chr(r1$split_ions_a) else ctx_chr(isolate(input$split_ions_a))
                            ib <- if (pair) ctx_chr(r1$split_ions_b) else ctx_chr(isolate(input$split_ions_b))
                            ic <- if (pair) ctx_chr(r1$split_ions_c) else ctx_chr(isolate(input$split_ions_c))
                            showModal(modalDialog(
                                title = paste0(if (pair) "Re-split " else "Split ", h$what, " in ", basename(h$smp)),
                                helpText(if (pair) "Already split. Change the ions and press Split to redo it, or Unsplit to make it one peak again. This sample only."
                                         else "Split this peak into two compounds by their diagnostic ions, in this sample only (the sidebar's Split does every sample). Pick ions that are strong in every scan of the hump."),
                                textInput("ctx_split_ions_a", "Ions, compound 1 (e.g. 282, 281)", value = ia),
                                textInput("ctx_split_ions_b", "Ions, compound 2 (e.g. 111, 125)", value = ib),
                                textInput("ctx_split_ions_c", "Ions, compound 3 (optional)", value = ic),
                                footer = tagList(
                                    if (pair) actionButton("ctx_unsplit_go", "Unsplit"),
                                    modalButton("Cancel"),
                                    actionButton("ctx_split_go", "Split")),
                                easyClose = TRUE))
                        }
                        observeEvent(input$chromatogram_ctx, {
                            a <- input$chromatogram_ctx
                            act <- if (is.list(a) && is.character(a$action) && length(a$action) == 1) a$action else ""
                            if (act == "remove") ctx_remove() else if (act == "split") ctx_split_dialog()
                        })
                        observeEvent(input$ctx_split_go, {
                            undo_begin("split one peak (right-click)"); on.exit(undo_end(), add = TRUE)
                            if (isTRUE(table_dirty())) { notify("Save the peak table first (Shift+Z): splitting redraws it."); return() }
                            pt <- coelution_table(); if (is.null(pt)) return()
                            rows <- ctx_find(pt)
                            if (is.null(rows)) { removeModal(); notify("Split: that peak has changed since you right-clicked it; right-click it again."); return() }
                            ion_sets <- gcms_split_ion_sets(input$ctx_split_ions_a, input$ctx_split_ions_b, input$ctx_split_ions_c)
                            if (is.character(ion_sets)) { notify(ion_sets); return() }
                            removeModal()
                            t <- ctx_target$t
                            sub <- pt[rows, , drop = FALSE]
                            wnd <- c(min(sub$peak_start_rt_offset), max(sub$peak_end_rt_offset))
                            id <- if (nzchar(t$split_id)) t$split_id else sprintf("split@%.4f-%.4f %s", wnd[1], wnd[2], basename(t$smp))
                            res <- gcms_apply_coelution_split(sub, coelution_traces, wnd, ion_sets,
                                                              crosstalk = coelution_crosstalk, k_mode = "sample", split_id = id)
                            notify_log(res$log)
                            np <- res$peak_table
                            ions_txt <- vapply(ion_sets, paste, character(1), collapse = ", ")
                            done <- sum(np$split_id %in% id & np$coelution_split %in% .gcms_split_labels &
                                        np$split_ions_a %in% ions_txt[1] & np$split_ions_b %in% ions_txt[2] &
                                        (if (length(ion_sets) == 3) np$split_ions_c %in% ions_txt[3] else is.na(np$split_ions_c))) == length(ion_sets)
                            if (!done) return()
                            coelution_finish(.gcms_renumber(rbind(.gcms_ensure_split_cols(pt)[-rows, , drop = FALSE], np)))
                        })
                        observeEvent(input$ctx_unsplit_go, {
                            undo_begin("unsplit one peak (right-click)"); on.exit(undo_end(), add = TRUE)
                            removeModal()
                            if (isTRUE(table_dirty())) { notify("Save the peak table first (Shift+Z): unsplitting redraws it."); return() }
                            pt <- coelution_table(); if (is.null(pt)) return()
                            rows <- ctx_find(pt)
                            if (length(rows) < 2) { notify("Unsplit: that split has changed since you right-clicked it; right-click it again."); return() }
                            sub <- pt[rows, , drop = FALSE]
                            res <- gcms_unsplit_peak(sub, sub$peak_number_within_sample[1])
                            notify_log(res$log)
                            coelution_finish(.gcms_renumber(rbind(.gcms_ensure_split_cols(pt)[-rows, , drop = FALSE], res$peak_table)))
                        })
                    ## Right-click menu (end)

                    ## Page indicator text (e.g. "Page 2 / 8   (samples 21-40 of 152)")
                        output$page_indicator <- renderText({
                            redraw_trigger(); current_page()
                            b <- page_bounds()
                            if (is.null(b)) return("")
                            if (b$n == 0) return("No samples loaded")
                            paste0("Page ", b$pg, " / ", b$n_pages,
                                   "   (samples ", b$lo, "–", b$hi, " of ", b$n, ")")
                        })

                    ## The chromatogram plot itself - reactive on redraw_trigger (Shift+Q)
                    ## and current_page (Prev/Next). Axis limits come from the globals that
                    ## Shift+Q set from the brush/defaults, so brushing does not itself force
                    ## a redraw (matches v4: press Q to refresh the view).
                        ## Where each panel sits, for the Align arrows beside the plot. Panels
                        ## are equal height, so the layout is three numbers: the margin
                        ## above the first panel, the gap between panels and the band below the
                        ## last (x axis + title), all fixed in pixels by the theme. They are
                        ## measured ONCE from a 2-panel skeleton of the real plot (same
                        ## facets, x axis, theme), through panel_fraction(); a full-size skeleton
                        ## agreed with the real plot to 0 px but cost ~3 s at 95 panels.
                        ## test_align.R pins the arithmetic against a full skeleton.
                        .align_layout <- new.env(parent = emptyenv())
                        align_layout_for <- function() {
                            k <- "detail"
                            if (!is.null(.align_layout[[k]])) return(.align_layout[[k]])
                            H0 <- 600
                            sk <- ggplot(data.frame(rt_rt_offset = c(0, 1), abundance = c(0, 1),
                                                    path_to_cdf_csv = factor(c("a", "b"))),
                                         aes(x = rt_rt_offset, y = abundance)) + geom_blank() +
                                scale_x_continuous(name = "Retention (Scan number)") +
                                scale_y_continuous(name = "Abundance (counts)") +
                                facet_grid(path_to_cdf_csv ~ ., scales = "free_y", drop = FALSE) + theme_classic()
                            g <- try(panel_fraction(sk, 1000, H0, view_key = paste0("align-skeleton-", k), n_facets = 2), silent = TRUE)
                            if (inherits(g, "try-error") || is.null(g$panels) || nrow(g$panels) != 2) return(NULL)
                            pn <- g$panels[order(g$panels$y_lo), ]
                            out <- c(top = pn$y_lo[1] * H0, gap = (pn$y_lo[2] - pn$y_hi[1]) * H0,
                                     bottom = (1 - pn$y_hi[2]) * H0)
                            assign(k, out, envir = .align_layout)
                            out
                        }
                        align_rows_for <- function(samples, height) {
                            n <- length(samples); L <- align_layout_for()
                            if (!n || is.null(L) || !is.finite(height)) return(list())
                            h <- (height - L[["top"]] - L[["bottom"]] - (n - 1) * L[["gap"]]) / n
                            if (!is.finite(h) || h <= 0) return(list())
                            lapply(seq_len(n), function(i) list(
                                sample = samples[i], label = unname(gcms_sample_id_from_path(samples[i])),
                                y = round(L[["top"]] + (i - 1) * (h + L[["gap"]]) + h / 2, 1)))
                        }

                        output$chromatograms <- renderPlot({

                            redraw_trigger()
                            b <- page_bounds()
                            plot_px_h <- chromatogram_height()
                            ## One line column per drawn pixel (gcms_thin_lines()); 2000 until
                            ## the browser has reported the plot's width.
                            thin_bins <- suppressWarnings(as.numeric(session$clientData$output_chromatograms_width))
                            if (!length(thin_bins) || !is.finite(thin_bins) || thin_bins < 100) thin_bins <- 2000

                            ## Returning NULL here draws a plot with NO coordinate map, and the
                            ## browser then holds a degenerate 0..1 map for this output until the
                            ## next successful render -- one route to the fraction-coordinate
                            ## brushes seen on host1. Hand back a real ggplot with a real x
                            ## mapping instead, so the output always has a usable coordmap.
                            if (is.null(b)) {
                                return(
                                    ggplot(data.frame(rt_rt_offset = range(chromatograms$rt), abundance = c(0, 0)),
                                           aes(x = rt_rt_offset, y = abundance)) +
                                        geom_blank() +
                                        annotate("text", x = mean(range(chromatograms$rt)), y = 0, size = 5,
                                                 label = "Press Shift+Q to load chromatograms") +
                                        theme_classic()
                                )
                            }
                            page_samples <- b$page_samples

                            ## Axis limits set by Shift+Q (globals); guard against unset
                                xs <- x_axis_start; xe <- x_axis_end
                                if (length(xs) == 0 || is.null(xs)) xs <- min(chromatograms$rt)
                                if (length(xe) == 0 || is.null(xe)) xe <- max(chromatograms$rt)
                                ys <- if (!is.null(y_axis_start) && length(y_axis_start) > 0) y_axis_start else 0
                                ye <- if (!is.null(y_axis_end)   && length(y_axis_end)   > 0) y_axis_end   else max(chromatograms$abundance)

                            ## Display settings (v6). The controls may not have rendered on the
                            ## first pass, so each falls back to the function argument.
                                mode_display <- if (is.null(input$display_mode)) {
                                    if (isTRUE(all_ion_view)) "all_ions" else "tic"
                                } else input$display_mode
                                mode_scale <- if (is.null(input$allion_scaling))    all_ion_scaling    else input$allion_scaling
                                ## The two EXPENSIVE controls are read through isolate(), so this
                                ## renderer does not take a reactive dependency on them. Dragging
                                ## either slider used to queue one multi-million-segment redraw per
                                ## intermediate value, each recomputing the ion window for the whole
                                ## page. They now apply on the next Shift+Q, which is what the help
                                ## text has always said. Display mode and y scaling stay live: both
                                ## are cheap and users expect them to be instant.
                                thr        <- isolate(if (is.null(input$allion_threshold))  all_ion_threshold  else input$allion_threshold)
                                npx        <- isolate(if (is.null(input$allion_max_points)) all_ion_max_points else input$allion_max_points)

                            ## Base data for this page only.
                            ## INCLUSIVE bounds, matching integration. Integration uses >= / <=
                            ## while the overlay and this view filter used > / <, so a peak sitting
                            ## flush against the zoom edge was integrated but vanished from the
                            ## plot, and a ribbon excluded the two boundary scans its own area
                            ## included.
                                cuf <- dplyr::filter(
                                    chromatograms_updated,
                                    path_to_cdf_csv %in% page_samples,
                                    rt_rt_offset >= xs & rt_rt_offset <= xe
                                )

                                ## EVERY sample on this page gets a panel, drawn or not. The facet
                                ## column is a factor over page_samples and the facets are built
                                ## with drop = FALSE, so a sample with nothing over the ion
                                ## threshold shows as an empty panel rather than vanishing. It used
                                ## to disappear entirely: a page returned 14 panels while the
                                ## indicator still read "samples 21-40 of 152", which reads as
                                ## "those six samples failed".
                                ## The peak overlay below is built from these rows in one pass
                                ## (gcms_overlay_rows()). Every overlaid peak lies inside [xs, xe] on
                                ## a page sample, so its rows are all here. Kept BEFORE the facet
                                ## factor below, so the overlay carries the same character column the
                                ## old per-peak filter returned. Pinned by test_overlay_onepass.R.
                                cuf_rows <- cuf
                                as_page_facet <- function(x) factor(as.character(x), levels = page_samples)
                                cuf$path_to_cdf_csv <- as_page_facet(cuf$path_to_cdf_csv)

                                ## Per SAMPLE, not per ROW. labeller() only ever looks up one
                                ## entry per facet, but this built a named vector as long as the
                                ## table -- two gsub passes and a names<- over ~1-2M elements on
                                ## every redraw, including every page flip.
                                facet_labels <- gcms_sample_id_from_path(page_samples)
                                facet_labels <- gcms_wrap_hard(facet_labels, facet_label_wrap)
                                names(facet_labels) <- page_samples

                            ## Gather the all-ion data for the page, one sample at a time.
                            ## Each sample carries its own rt offset and its own baseline, so
                            ## the window is queried in that sample's native rt and the offset
                            ## is added back afterwards. A sample that yields nothing (no ion
                            ## clears the threshold in this window) is simply skipped.
                                allion <- NULL
                                if (mode_display %in% c("all_ions", "ion_map")) {
                                    npx_used <- if (mode_display == "ion_map") min(npx, 600) else npx
                                    ## Read the offsets FRESH. Shift+Q re-reads
                                    ## samples_monolist into a handler-LOCAL variable, so this
                                    ## renderer used to resolve the name to the startup copy:
                                    ## editing an rt_offset moved the TIC and the peak marks
                                    ## but not the ion lines, which in all-ions mode manufactures
                                    ## exactly the apex shift the view exists to detect. One
                                    ## small CSV per render is nothing beside the draw.
                                    offsets_now <- if (file.exists("samples_monolist.csv")) {
                                        try(read.csv("samples_monolist.csv", stringsAsFactors = FALSE), silent = TRUE)
                                    } else NULL
                                    if (inherits(offsets_now, "try-error") || is.null(offsets_now) ||
                                        !all(c("rt_offset", "path_to_cdf_csv") %in% names(offsets_now))) {
                                        offsets_now <- samples_monolist
                                    }
                                    pieces <- list()
                                    for (s in page_samples) {
                                        off <- offsets_now$rt_offset[match(s, offsets_now$path_to_cdf_csv)]
                                        if (length(off) == 0 || is.na(off)) off <- 0
                                        w <- try(allion_window(s, xs - off, xe - off, thr, npx_used, rt_shift = off), silent = TRUE)
                                        if (inherits(w, "try-error") || is.null(w) || nrow(w) == 0) next
                                        w <- gcms_allion_scale(w, mode_scale)
                                        ## Round the reconstructed grid. allion_window() snaps every
                                        ## sample to one shared set of centres in offset-corrected
                                        ## space, but it hands them back shifted into native rt, and
                                        ## (c - off) + off is not bit-identical to c. Those ~1e-15
                                        ## differences make unique() see one grid per sample again,
                                        ## which is what made the ion map's tiles too narrow. Six
                                        ## decimals is orders of magnitude finer than any real
                                        ## retention-time resolution.
                                        w$rt_rt_offset  <- round(w$rt + off, 6)
                                        w$path_to_cdf_csv <- s
                                        pieces[[length(pieces) + 1L]] <- w
                                    }
                                    if (length(pieces) > 0) allion <- as.data.frame(data.table::rbindlist(pieces))
                                    if (!is.null(allion)) allion$path_to_cdf_csv <- as_page_facet(allion$path_to_cdf_csv)

                                    ## Single-detector data (GC-TCD/FID/ECD) carries the mz = 0
                                    ## sentinel, and the all-ion views colour and position by m/z
                                    ## on a LOG scale: log10(0) is -Inf, so every line came back
                                    ## grey with an empty legend and the ion map's y axis collapsed
                                    ## to c(0, 0). Fall back to the TIC, as Shift+1 and Shift+4
                                    ## already do for the same reason.
                                    if (!is.null(allion) && all(allion$mz == 0, na.rm = TRUE)) {
                                        notify("All-ion view needs mass-spec data. Showing the TIC.")
                                        mode_display <- "tic"
                                        allion <- NULL
                                    }
                                }

                                ## Nothing to draw in the requested mode -> fall back to the TIC
                                ## rather than showing an empty panel.
                                if (mode_display %in% c("all_ions", "ion_map") && is.null(allion)) mode_display <- "tic"

                                ## Nothing in range at all -- e.g. axis limits left over from a
                                ## failed brush. facet_grid() errors on an empty frame ("Faceting
                                ## variables must have at least one value"), which in a renderPlot
                                ## leaves a broken panel and no way back, so say what happened.
                                if (nrow(cuf) == 0 && is.null(allion)) {
                                    ## Carry a REAL x mapping over the current window. The old
                                    ## placeholder was annotate(x = 0, y = 0) + theme_void(): a
                                    ## degenerate 0..0 domain with no panel, which is exactly the
                                    ## "never hand back a plot with no coordinate map" trap -- and
                                    ## it appeared precisely when the user was trying to brush
                                    ## their way out of a bad zoom, leaving the next brush with
                                    ## nothing sane to convert against.
                                    mid <- if (is.finite(xs) && is.finite(xe)) (xs + xe) / 2 else 0
                                    return(
                                        ggplot(data.frame(rt_rt_offset = c(xs, xe), abundance = c(0, 1)),
                                               aes(x = rt_rt_offset, y = abundance)) +
                                            geom_blank() +
                                            annotate("text", x = mid, y = 0.5, size = 5, lineheight = 1.2,
                                                     label = paste0(
                                                         "No chromatogram points between x = ",
                                                         signif(xs, 6), " and ", signif(xe, 6), ".\n",
                                                         "Press Shift+Q with no brush to reset the view.")) +
                                            scale_x_continuous(limits = c(xs, xe), name = "Retention (Scan number)", breaks = scales::breaks_pretty(n = 20), minor_breaks = function(lims) scales::breaks_pretty(n = 80)(lims)) +
                                            scale_y_continuous(name = "Abundance (counts)") +
                                            theme_classic()
                                    )
                                }

                            ## A y-zoom only means something in the view it was taken in --
                            ## raw counts, square root, log, per-ion-normalised and m/z are five
                            ## different axes. Tag it with the view, and ignore it elsewhere.
                            ## drawn_yrange is what the panel spans when no zoom applies, and is
                            ## also what a later brush's y fractions get measured against.
                                view_key <- paste(mode_display, mode_scale, sep = "/")
                                ## Pin the bottom of the y axis at 0 wherever y is an
                                ## abundance. The ion map's y is m/z, and log10 scaling
                                ## draws logs -- 0 is an interior value there, not a floor.
                                zero_floor_y <- if (mode_display == "ion_map") {
                                    FALSE
                                } else if (mode_display == "all_ions") {
                                    !identical(mode_scale, "log")
                                } else {
                                    TRUE
                                }
                                drawn_yrange <- if (mode_display == "ion_map" && !is.null(allion)) {
                                    range(allion$mz, na.rm = TRUE)
                                } else if (mode_display == "all_ions" && !is.null(allion)) {
                                    range(allion$y, na.rm = TRUE)
                                } else {
                                    c(ys, ye)
                                }

                            ## Base plot, per display mode
                                if (mode_display == "all_ions") {

                                    chromatogram_plot <- ggplot() +
                                        geom_line(
                                            data = allion,
                                            mapping = aes(x = rt_rt_offset, y = y, colour = mz, alpha = opacity,
                                                          group = interaction(path_to_cdf_csv, mz)),
                                            linewidth = 0.18
                                        ) +
                                        scale_alpha_identity() +
                                        ## log colour: nearly every ion in a GC-MS run sits
                                        ## between 40 and 150, so a linear m/z ramp renders the
                                        ## whole plot one shade of purple and wastes the top
                                        ## two-thirds of the scale on a handful of heavy ions.
                                        scale_colour_viridis_c(name = "m/z", trans = "log10") +
                                        scale_x_continuous(limits = c(xs, xe), name = "Retention (Scan number)", breaks = scales::breaks_pretty(n = 20), minor_breaks = function(lims) scales::breaks_pretty(n = 80)(lims)) +
                                        scale_y_continuous(name = gcms_allion_y_label(mode_scale),
                                                           limits = y_limits_for(view_key, drawn_yrange,
                                                                                 zero_floor = zero_floor_y),
                                                           oob = scales::squish) +
                                        facet_grid(path_to_cdf_csv~., scales = "free_y", drop = FALSE, labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                        theme_classic() +
                                        guides(fill = "none")

                                } else if (mode_display == "ion_map") {

                                    ## geom_tile rather than geom_raster: decimation can leave
                                    ## the rt grid with gaps, which geom_raster rejects.
                                    ## Tile width from the ACTUAL spacing of the drawn grid.
                                    ## Dividing the window by the COUNT of unique x values was
                                    ## several times too narrow whenever the union across samples
                                    ## was denser than any one sample's grid -- the map rendered as
                                    ## thin stripes with white gaps, hidden on a single-sample page.
                                    ## It is also wrong whenever decimation does not fire at all.
                                    ## Samples now share one grid (see allion_window's rt_shift),
                                    ## so the median gap is the true tile width.
                                    tile_x <- sort(unique(allion$rt_rt_offset))
                                    tile_w <- if (length(tile_x) > 1) stats::median(diff(tile_x)) else (xe - xs)
                                    if (!is.finite(tile_w) || tile_w <= 0) tile_w <- (xe - xs) / max(1, length(tile_x))
                                    chromatogram_plot <- ggplot() +
                                        geom_tile(
                                            data = allion,
                                            mapping = aes(x = rt_rt_offset, y = mz, fill = y),
                                            width = tile_w, height = 1
                                        ) +
                                        scale_fill_viridis_c(name = gcms_allion_y_label(mode_scale)) +
                                        ## coord_cartesian, not scale limits: a tile straddling
                                        ## the window edge should be clipped, not dropped.
                                        ## Both limits go through coord_cartesian, not the scales:
                                        ## tiles have width and height, so a tile straddling either
                                        ## edge should be clipped, not dropped.
                                        coord_cartesian(xlim = c(xs, xe),
                                                        ylim = y_limits_for(view_key, drawn_yrange,
                                                                            zero_floor = zero_floor_y)) +
                                        scale_x_continuous(name = "Retention (Scan number)", breaks = scales::breaks_pretty(n = 20), minor_breaks = function(lims) scales::breaks_pretty(n = 80)(lims)) +
                                        scale_y_continuous(name = "m/z") +
                                        facet_grid(path_to_cdf_csv~., scales = "free_y", drop = FALSE, labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                        theme_classic()

                                } else {

                                    ## Named, so the key in the button row (gcms_ion_legend())
                                    ## shows exactly these colours. Same order as the unnamed
                                    ## scale gave: sorted as character.
                                    ion_lv  <- sort(unique(as.character(cuf$ion[as.character(cuf$ion) != "baseline"])))
                                    ion_pal <- stats::setNames(rep_len(discrete_palette, length(ion_lv)), ion_lv)

                                    chromatogram_plot <- ggplot() +
                                        geom_line(
                                            data = gcms_thin_lines(filter(cuf, ion == "baseline"), xs, xe, thin_bins),
                                            mapping = aes(x = rt_rt_offset, y = abundance), color = "grey"
                                        ) +
                                        scale_x_continuous(limits = c(xs, xe), name = "Retention (Scan number)", breaks = scales::breaks_pretty(n = 20), minor_breaks = function(lims) scales::breaks_pretty(n = 80)(lims)) +
                                        scale_y_continuous(limits = y_limits_for(view_key, c(ys, ye),
                                                                                 zero_floor = zero_floor_y),
                                                           name = "Abundance (counts)", oob = scales::squish) +
                                        facet_grid(path_to_cdf_csv~., scales = "free_y", drop = FALSE, labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                        theme_classic() +
                                        guides(fill = "none") +
                                        scale_fill_continuous(type = "viridis") +
                                        scale_color_manual(values = ion_pal, guide = "none")

                                }
                                ion_legend_pal(if (mode_display == "tic") ion_pal else NULL)

                            ## Peak overlay (read the table Shift+Q wrote), page + x-axis filtered.
                            ## The shaded ribbons are drawn in raw counts, so they only make sense
                            ## against the TIC's raw y-axis; in the other modes the peak bounds and
                            ## numbers are still drawn, but the shading is dropped.
                                if (file.exists("peaks_monolist.csv")) {
                                    peak_table <- read.csv("peaks_monolist.csv")
                                    needed <- c("peak_start_rt_offset", "peak_end_rt_offset", "peak_number_within_sample")
                                    if (nrow(peak_table) > 0 && all(needed %in% colnames(peak_table))) {
                                        peak_table <- dplyr::filter(
                                            peak_table,
                                            path_to_cdf_csv %in% page_samples,
                                            peak_start_rt_offset >= xs & peak_end_rt_offset <= xe
                                        )
                                        ## v7: a split pair is drawn as two overlapping curves (below),
                                        ## each compound's own share of the TIC, not as two drop-line
                                        ## ribbons that would hide how the hump was shared out.
                                        split_rows <- if (mode_display == "tic" && "split_status" %in% colnames(peak_table))
                                            which(peak_table$split_status %in% "split") else integer(0)
                                        split_peaks <- peak_table[split_rows, , drop = FALSE]
                                        if (length(split_rows)) peak_table <- peak_table[-split_rows, , drop = FALSE]
                                            all_ribbons <- list(); all_vlines <- list(); all_labels <- list()
                                            if (nrow(peak_table) > 0) {
                                                ## One pass over the page's rows for every peak (gcms_overlay_rows()).
                                                overlay <- gcms_overlay_rows(cuf_rows, peak_table)
                                                all_ribbons <- overlay$ribbons
                                                all_vlines  <- overlay$vlines
                                                all_labels  <- overlay$labels
                                            ## Same facet factor as the base layer, or these land
                                            ## in panels of their own.
                                            if (nrow(all_ribbons) > 0) all_ribbons$path_to_cdf_csv <- as_page_facet(all_ribbons$path_to_cdf_csv)
                                            if (nrow(all_vlines)  > 0) all_vlines$path_to_cdf_csv  <- as_page_facet(all_vlines$path_to_cdf_csv)
                                            if (nrow(all_labels)  > 0) all_labels$path_to_cdf_csv  <- as_page_facet(all_labels$path_to_cdf_csv)
                                            if (nrow(all_ribbons) > 0) {
                                                if (mode_display == "tic") {
                                                    chromatogram_plot <- chromatogram_plot +
                                                        geom_vline(data = all_vlines, mapping = aes(xintercept = rt_rt_offset), alpha = 0.3) +
                                                        geom_ribbon(
                                                            data = all_ribbons,
                                                            mapping = aes(x = rt_rt_offset, ymax = abundance, ymin = baseline,
                                                                          fill = peak_number_within_sample, group = peak_number_within_sample),
                                                            alpha = 0.8
                                                        ) +
                                                        ## A square at each apex coloured by how well the peak's
                                                        ## spectrum agrees with the folder's other peaks of the same
                                                        ## peak_ID, drawn under the peak number so it stays readable.
                                                        gcms_detail_marker_layers(gcms_overview_markers(all_labels, peak_table, confidence_thresholds)) +
                                                        geom_text(data = all_labels, mapping = aes(label = peak_number_within_sample, x = x, y = y), color = "black", size = peak_label_size)
                                                } else {
                                                    chromatogram_plot <- chromatogram_plot +
                                                        geom_vline(data = all_vlines, mapping = aes(xintercept = rt_rt_offset), alpha = 0.3)
                                                }
                                            }
                                        }
                                        if (nrow(split_peaks) > 0 && file.exists("coelution_profiles.csv")) {
                                            prof <- try(utils::read.csv("coelution_profiles.csv", stringsAsFactors = FALSE), silent = TRUE)
                                            ## Only bands whose row is still in the table: matching on the
                                            ## sample alone drew BOTH bands after one row of a pair had been
                                            ## removed (the "residual content", 2026-10-16).
                                            if (!inherits(prof, "try-error") && nrow(prof) > 0 && "split_id" %in% names(prof)) {
                                                key_pk <- paste(split_peaks$path_to_cdf_csv, split_peaks$split_id, split_peaks$coelution_split, sep = "\r")
                                                key_pr <- paste(prof$path_to_cdf_csv, prof$split_id, prof$part, sep = "\r")
                                                prof$peak_number_within_sample <- split_peaks$peak_number_within_sample[match(key_pr, key_pk)]
                                                prof <- prof[!is.na(prof$peak_number_within_sample) &
                                                             prof$rt_rt_offset >= xs & prof$rt_rt_offset <= xe, , drop = FALSE]
                                                ## The topmost band of each split, which needs no dividing line.
                                                prof$top_part <- stats::ave(prof$part, paste(prof$path_to_cdf_csv, prof$split_id, sep = "\r"), FUN = max)
                                            } else prof <- NULL
                                            if (!is.null(prof) && nrow(prof) > 0) {
                                                prof$grp <- paste(prof$path_to_cdf_csv, prof$split_id, prof$part, sep = "\r")
                                                ## Label each band where it is thickest, halfway up it.
                                                ## A zero-area band (no signal for that compound in this sample) has
                                                ## nowhere to sit, so it gets no label here; zero_labels below put its
                                                ## number in grey just under the hump instead.
                                                split_labels <- do.call(rbind, lapply(split(prof, prof$grp), function(g) {
                                                    if (max(g$top - g$ymin) <= 0) return(NULL)
                                                    i <- which.max(g$top - g$ymin)
                                                    data.frame(path_to_cdf_csv = g$path_to_cdf_csv[i], rt_rt_offset = g$rt_rt_offset[i],
                                                               y = (g$top[i] + g$ymin[i]) / 2,
                                                               peak_number_within_sample = g$peak_number_within_sample[i])
                                                }))
                                                if (is.null(split_labels))
                                                    split_labels <- data.frame(path_to_cdf_csv = character(0), rt_rt_offset = numeric(0),
                                                                               y = numeric(0), peak_number_within_sample = numeric(0))
                                                ## Zero-area peaks of a split (Lucas, 2026-10-26): their numbers, in grey,
                                                ## centred just below the bottom of the hump, so the user sees the peak
                                                ## was made. One label per split ("4" or "4, 5").
                                                zero_labels <- gcms_split_zero_labels(prof)
                                                zero_labels$path_to_cdf_csv <- as_page_facet(zero_labels$path_to_cdf_csv)
                                                prof$path_to_cdf_csv <- as_page_facet(prof$path_to_cdf_csv)
                                                split_labels$path_to_cdf_csv <- as_page_facet(split_labels$path_to_cdf_csv)
                                                if (mode_display == "tic") {
                                                    chromatogram_plot <- chromatogram_plot +
                                                        geom_ribbon(data = prof,
                                                                    mapping = aes(x = rt_rt_offset, ymin = ymin, ymax = top,
                                                                                  fill = peak_number_within_sample, group = grp),
                                                                    alpha = 0.6) +
                                                        geom_line(data = prof[prof$part != prof$top_part, , drop = FALSE],
                                                                  mapping = aes(x = rt_rt_offset, y = top, group = grp),
                                                                  colour = "black", linewidth = 0.3) +
                                                        geom_text(data = split_labels,
                                                                  mapping = aes(label = peak_number_within_sample, x = rt_rt_offset, y = y),
                                                                  colour = "black", size = peak_label_size) +
                                                        geom_text(data = zero_labels,
                                                                  mapping = aes(label = label, x = rt_rt_offset, y = y),
                                                                  colour = "#999999", size = peak_label_size, vjust = 1.3)
                                                }
                                            }
                                        }
                                    }
                                }

                            ## TIC trace on top (page samples only). In the all-ion modes the
                            ## overlay IS the chromatogram, and a raw-counts TIC would sit far
                            ## off the transformed y-scale.
                                if (mode_display == "tic") {
                                    chromatogram_plot <- chromatogram_plot +
                                        geom_line(
                                            data = gcms_thin_lines(filter(cuf, ion != "baseline"), xs, xe, thin_bins),
                                            mapping = aes(x = rt_rt_offset, y = abundance, color = ion),
                                            alpha = 0.8
                                        )
                                }

                            ## Any theme change goes on BEFORE the plot is kept for the brush
                            ## maths, so panel_fraction() measures the panels as drawn.
                            ## Sample names (the strips on the right): horizontal, larger, and
                            ## wrapped by facet_labels above (Lucas, 2026-10-26).
                                chromatogram_plot <- chromatogram_plot +
                                    theme(strip.text.y = element_text(angle = 0, hjust = 0, size = facet_label_size))

                            ## Kept so a mapping-less brush can be re-expressed against the
                            ## panel's real bounds. clientData carries the size the browser
                            ## actually rendered at, which is what the fractions refer to.
                                last_brush_plot   <<- chromatogram_plot
                                ## The height is the one this render is drawn at (renderPlot's
                                ## height function), not clientData's, which still holds the
                                ## previous picture's height when a mode or page-size change
                                ## resizes the plot.
                                last_brush_size   <<- suppressWarnings(as.numeric(c(
                                    session$clientData$output_chromatograms_width,
                                    plot_px_h)))
                                last_brush_key    <<- view_key
                                last_brush_facet_sig <<- if (mode_display %in% c("all_ions", "ion_map") && !is.null(allion)) {
                                    sort(unique(as.character(allion$path_to_cdf_csv)))
                                } else {
                                    sort(unique(as.character(cuf$path_to_cdf_csv)))
                                }
                                last_brush_facets <<- length(last_brush_facet_sig)
                                last_brush_yrange <<- y_limits_for(view_key, drawn_yrange,
                                                                   zero_floor = zero_floor_y)
                                last_brush_window <<- suppressWarnings(as.numeric(
                                    c(xs, xe, last_brush_yrange)))
                                ## One Align arrow pair per panel, top to bottom in page order
                                ## (the facets are a factor over page_samples, drop = FALSE).
                                session$sendCustomMessage("gcms_align_rows",
                                    list(height = plot_px_h, rows = align_rows_for(page_samples, plot_px_h)))

                            chromatogram_plot
                        }, height = chromatogram_height)

                ## Transfer chromatogram_brush info to selected_peak table

                    output$selected_peak <- DT::renderDataTable(DT::datatable({

                        if ( !is.null(input$chromatogram_brush )) {
                            peak_points <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                            peak_data <-  data.frame(
                                peak_start = min(peak_points$rt),
                                peak_end = max(peak_points$rt),
                                peak_ID = "unknown",
                                path_to_cdf_csv = peak_points$path_to_cdf_csv[1],
                                area = sum(peak_points$abundance[peak_points$ion == 0])
                            )
                            peak_data
                        } else {
                            NULL
                        }

                    }))

                ## Single-peak remove => SHIFT+E => ASCII code 69
                      observeEvent(input$keypress, {
                        if( input$keypress == 69 ) { # Shift+E, "excise"
                            undo_begin("remove peak, brushed samples (Shift+E)"); on.exit(undo_end(), add = TRUE)

                          cat("Excising single peak...\n")
                          if ( !is.null(input$chromatogram_brush )) {

                            peak_points <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush, multi = TRUE)
                            if (nrow(peak_points) == 0) {
                                notify("Nothing excised: the selection held no data.")
                                return()
                            }
                            ## Several samples only when the brush covers their panels' middles
                            ## (as Shift+A); each sample is excised over its OWN raw-rt window.
                            if (!isTRUE(attr(peak_points, "gcms_multi_panel"))) {
                                path_to_cdf_csv <- gcms_single_sample(peak_points, "Shift+E")
                                if (is.null(path_to_cdf_csv)) return()
                            }
                            win <- gcms_brush_windows_by_sample(peak_points)

                            peak_table <- read.csv("peaks_monolist.csv")
                            if (nrow(peak_table) == 0) {
                                notify("No peaks to excise.")
                                return()
                            }

                            # Remove any peak that fully resides within its sample's window. INCLUSIVE
                            # bounds: with `>` and `<` a peak could not be excised by the very brush
                            # that created it, because Shift+A set its bounds to exactly this
                            # selection's min and max.
                            n_before <- nrow(peak_table)
                            peak_table <- gcms_excise_windows(peak_table, win)
                            n_gone <- n_before - nrow(peak_table)

                            write.table(
                              x = peak_table,
                              file = "peaks_monolist.csv",
                              append = FALSE,
                              row.names = FALSE,
                              col.names = TRUE,
                              sep = ","
                            )
                            notify(if (nrow(win) > 1) paste0("Excised ", n_gone, " peak", if (n_gone == 1) "" else "s",
                                                              " across ", nrow(win), " samples.") else "Peak excised.")
                          } else {
                            notify("Nothing excised: draw a selection first.")
                          }
                        }
                      })

                ## Global peak remove => SHIFT+R => ASCII code 82
                      observeEvent(input$keypress, {
                        if (input$keypress == 82) { # Shift+R, "remove" globally
                            undo_begin("remove peak, all samples (Shift+R)"); on.exit(undo_end(), add = TRUE)

                          cat("Removing selected peaks GLOBALLY...\n")
                          if ( !is.null(input$chromatogram_brush )) {

                            peak_points <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                            if (nrow(peak_points) == 0) {
                                notify("Nothing removed: the selection held no data.")
                                return()
                            }

                            # read peak table
                            peak_table <- read.csv("peaks_monolist.csv")
                            if (nrow(peak_table) == 0) {
                                notify("No peaks to remove.")
                                return()
                            }

                            # A GLOBAL removal has to compare like with like. The brush is in the
                            # brushed sample's NATIVE rt, while every other sample's peak bounds are
                            # in its own native rt -- so comparing the two across samples is wrong by
                            # each sample's rt_offset. Use the offset-corrected pair on both sides
                            # when the table carries it (it does after any Shift+Q), and fall back to
                            # native rt only on a table that has never been through one.
                            aligned <- all(c("peak_start_rt_offset", "peak_end_rt_offset") %in% names(peak_table)) &&
                                       "rt_rt_offset" %in% names(peak_points) &&
                                       !any(is.na(peak_table$peak_start_rt_offset))
                            if (aligned) {
                                selection_start <- min(peak_points$rt_rt_offset)
                                selection_end   <- max(peak_points$rt_rt_offset)
                                starts <- peak_table$peak_start_rt_offset
                                ends   <- peak_table$peak_end_rt_offset
                            } else {
                                selection_start <- min(peak_points$rt)
                                selection_end   <- max(peak_points$rt)
                                starts <- peak_table$peak_start
                                ends   <- peak_table$peak_end
                                cat("  (no rt_offset columns yet - comparing native rt; press Shift+Q first if samples are aligned)\n")
                            }

                            # remove ANY peak in ANY file that fully resides in [start, end];
                            # inclusive, for the same reason as Shift+E
                            peak_table <- peak_table[!
                              apply(cbind(
                                starts >= selection_start,
                                ends <= selection_end
                              ), 1, all)
                            ,]

                            write.table(
                              x = peak_table,
                              file = "peaks_monolist.csv",
                              append = FALSE,
                              row.names = FALSE,
                              col.names = TRUE,
                              sep = ","
                            )
                            notify("Peak removed from all samples.")
                          } else {
                            notify("Nothing removed: draw a selection first.")
                          }
                        }
                      })

                ## Append single peak with Shift+A (65)

                    observeEvent(input$keypress, {

                        # Do nothing if no selection
                            if(is.null(input$chromatogram_brush)) {
                                return()
                            }

                        # If selection and "a" is pressed, add the selection to the peak table
                            if( input$keypress == 65 ) {
                                undo_begin("add peak, brushed samples (Shift+A)"); on.exit(undo_end(), add = TRUE)

                                ## Evaluate the brush ONCE. It was called four times here,
                                ## each re-running the whole coordinate conversion and printing
                                ## its own diagnostic line.
                                pp <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush, multi = TRUE)

                                ## An empty selection must not be written. brushed_chromatogram()
                                ## returns a 0-row frame on several paths (brush outside the data,
                                ## unusable x range, nothing matched), and min/max of nothing give
                                ## Inf/-Inf -- which used to be appended to peaks_monolist.csv as
                                ## a permanent poison row that the Shift+Q dedupe does not drop.
                                if (nrow(pp) == 0) {
                                    notify("No peak added: the selection held no data.")
                                    return()
                                }
                                ## A brush covering the middles of several panels adds one peak
                                ## to EACH of those samples (Lucas, 2026-10-22), each from its own
                                ## points. Anything else must still resolve to one sample.
                                if (!isTRUE(attr(pp, "gcms_multi_panel"))) {
                                    pp_sample <- gcms_single_sample(pp, "Shift+A")
                                    if (is.null(pp_sample)) return()
                                }
                                win <- gcms_brush_windows_by_sample(pp)

                                ## area = TIC minus baseline, both read out of the LONG format
                                ## (rows tagged ion == 0 / ion == "baseline", value in `abundance`).
                                ## Fixed 2026-09-30: this asked for $tic and $baseline, which are
                                ## columns of modules/gcms.R's WIDE table, not of this one. Both
                                ## were NULL, sum(NULL) is 0, so every hand-added peak was written
                                ## with area = 0.
                                write.table(
                                    x = data.frame(
                                            peak_start = win$start,
                                            peak_end = win$end,
                                            peak_ID = "unknown",
                                            path_to_cdf_csv = win$path_to_cdf_csv,
                                            area = win$area
                                        ),
                                    file = "peaks_monolist.csv",
                                    append = TRUE,
                                    row.names = FALSE,
                                    col.names = FALSE,
                                    sep = ","
                                )

                                render_peak_table()
                                notify(if (nrow(win) > 1) paste0("Peak added to ", nrow(win), " samples.") else "Peak added.")
                            }
                    })

                ## Global append peak with Shift+G (71)

                    observeEvent(input$keypress, {

                        # Do nothing if no selection
                            if(is.null(input$chromatogram_brush)) {
                                return()
                            }

                        # If selection and "G" is pressed, add the selection to the peak table
                            if( input$keypress == 71 ) {
                                undo_begin("add peak, all samples (Shift+G)"); on.exit(undo_end(), add = TRUE)

                                ## Evaluate the brush once (see Shift+A) and refuse an empty one.
                                pp <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                                if (nrow(pp) == 0) {
                                    notify("No peak added: the selection held no data.")
                                    return()
                                }

                                x_peaks <-  data.frame(
                                                peak_start = min(pp$rt_rt_offset),
                                                peak_end = max(pp$rt_rt_offset),
                                                peak_ID = "unknown",
                                                path_to_cdf_csv = unique(chromatograms_updated$path_to_cdf_csv),
                                                area = NA_real_
                                            )

                                ## Shift+G writes the SAME window to every sample, so the area has
                                ## to be integrated per sample rather than shared. The old code
                                ## computed one scalar from the brushed sample and gave it to all
                                ## of them -- and that scalar was 0 anyway, because $tic/$baseline
                                ## are columns of modules/gcms.R's wide table, not of this one.
                                ## Fixed 2026-09-30.
                                for (row_i in seq_len(nrow(x_peaks))) {
                                    win <- dplyr::filter(
                                        chromatograms_updated,
                                        path_to_cdf_csv == x_peaks$path_to_cdf_csv[row_i],
                                        rt_rt_offset >= x_peaks$peak_start[row_i],
                                        rt_rt_offset <= x_peaks$peak_end[row_i]
                                    )
                                    x_peaks$area[row_i] <-
                                        sum(win$abundance[win$ion == 0]) -
                                        sum(win$abundance[win$ion == "baseline"])
                                }

                                x_peaks$peak_start <- x_peaks$peak_start - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]
                                x_peaks$peak_end <- x_peaks$peak_end - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]

                                write.table(
                                    x = x_peaks,
                                    file = "peaks_monolist.csv",
                                    append = TRUE,
                                    row.names = FALSE,
                                    col.names = FALSE,
                                    sep = ","
                                )

                                render_peak_table()
                                notify("Peak added to all samples.")
                            }
                    })

                ## Mass-spectrum keys: Shift+1 extract, Shift+2 re-plot, Shift+3 subtract, Shift+4 classify peaks

                    observeEvent(input$keypress, {

                        ## If "shift+1", MS from chromatogram brush -> MS_out_1

                            if( input$keypress == 33 ) {

                                ## Evaluate the brush once (it was converted three times here),
                                ## refuse an empty one, and refuse a selection that spans more than
                                ## one sample rather than extracting the spectrum from whichever
                                ## sample happens to sort first.
                                ms_pp <- brushed_chromatogram(chromatograms_updated, input$chromatogram_brush)
                                if (is.null(ms_pp) || nrow(ms_pp) == 0) {
                                    notify("No spectrum extracted: the selection held no data.")
                                    return()
                                }
                                sample_name_MS <- gcms_single_sample(ms_pp, "Shift+1")
                                if (is.null(sample_name_MS)) return()
                                ret_start_MS <- min(ms_pp$rt)
                                ret_end_MS <- max(ms_pp$rt)

                                chromatogram_updated_MS <- filter(chromatograms_updated, path_to_cdf_csv == sample_name_MS)

                                MS_ret_start_line <<- chromatogram_updated_MS$rt_first_row_in_raw[which.min(abs(
                                    chromatogram_updated_MS$rt - ret_start_MS
                                ))] + 1

                                MS_ret_end_line <<- chromatogram_updated_MS$rt_last_row_in_raw[which.min(abs(
                                    chromatogram_updated_MS$rt - ret_end_MS
                                ))] + 1

                                if (.Platform$OS.type == "unix") {

                                    ## shQuote every path: CDF_directory_path is caller-supplied and
                                    ## a space in it used to split the command, so `head`/`sed` wrote
                                    ## temp_MS.csv somewhere else entirely and the read below picked
                                    ## up the PREVIOUS peak's spectrum with no error.
                                    ms_src <- shQuote(file.path(CDF_directory_path, sample_name_MS))
                                    ms_tmp <- shQuote(file.path(CDF_directory_path, "temp_MS.csv"))
                                    system(paste0("head -1 ", ms_src, " > ", ms_tmp))
                                    system(paste0("sed -n ", MS_ret_start_line, ",", MS_ret_end_line, "p ", shQuote(sample_name_MS), " >> ", ms_tmp))
                                    framedDataFile <- readMonolist(paste0(CDF_directory_path, "/temp_MS.csv"))

                                }

                                if (.Platform$OS.type == "windows") {

                                    # system(paste0("head -1"," ",CDF_directory_path,"/",sample_name_MS," > ",CDF_directory_path,"/temp_MS.csv"))
                                    # system(paste0("sed -n ",MS_ret_start_line,",",MS_ret_end_line,"p ",sample_name_MS," >> ",CDF_directory_path,"/temp_MS.csv"))
                                    framedDataFile <-   isolate(
                                        as.data.frame(
                                            data.table::fread(sample_name_MS)
                                        )
                                    )

                                }

                                framedDataFile <- isolate(dplyr::filter(framedDataFile, rt > ret_start_MS, rt < ret_end_MS))
                                framedDataFile$mz <- round(framedDataFile$mz, 1)
                                framedDataFile <- framedDataFile %>% group_by(mz) %>% summarize(intensity = sum(intensity))
                                MS_out_1 <<- as.data.frame(framedDataFile)
                            }

                        ## If "shift+3" subtract chromatogram brush from MS_out_1 -> MS_out_1

                            if ( input$keypress == 35 ) {

                                if (is.null(MS_out_1)) {
                                    notify("No mass spectrum extracted yet. Use Shift+1 first.")
                                    return()
                                } else {

                                    ## Same again: one conversion, and the subtraction has to come
                                    ## from a single sample's panel.
                                    sub_pp <- isolate(brushed_chromatogram(chromatograms_updated, input$chromatogram_brush))
                                    if (is.null(sub_pp) || nrow(sub_pp) == 0) {
                                        notify("Nothing subtracted: the selection held no data.")
                                        return()
                                    }
                                    sub_sample <- gcms_single_sample(sub_pp, "Shift+3")
                                    if (is.null(sub_sample)) return()
                                    framedDataFile_to_subtract <- isolate(as.data.frame(
                                                        data.table::fread(sub_sample)
                                    ))
                                    framedDataFile_to_subtract <- isolate(dplyr::filter(
                                                            framedDataFile_to_subtract,
                                                            rt > min(sub_pp$rt),
                                                            rt < max(sub_pp$rt)
                                                        ))
                                    framedDataFile_to_subtract$mz <- round(framedDataFile_to_subtract$mz, 1)
                                    framedDataFile_to_subtract <- framedDataFile_to_subtract %>% group_by(mz) %>% summarize(intensity = sum(intensity))
                                    framedDataFile_to_subtract <<- as.data.frame(framedDataFile_to_subtract)
                                    joined <- left_join(MS_out_1, framedDataFile_to_subtract, by = "mz")
                                    joined$intensity.y[is.na(joined$intensity.y)] <- 0
                                    MS_out_1$intensity <- joined$intensity.x - joined$intensity.y
                                    MS_out_1$intensity[MS_out_1$intensity < 0] <- 0
                                    MS_out_1 <<- MS_out_1
                                }
                            }

                        ## If "shift+1", or "shift+2", or "shift+3, make plot based on brush if any

                            if ( input$keypress == 33 | input$keypress == 64 | input$keypress == 35 ) {

                                if (is.null(MS_out_1)) {
                                    notify("No mass spectrum extracted yet. Use Shift+1 first.")
                                    return()
                                } else if (nrow(MS_out_1) == 0 || (nrow(MS_out_1) == 1 && all(MS_out_1$mz == 0))) {
                                    # Single-detector data (e.g. GC-TCD): no mass spectrum exists.
                                    output$massSpectra_1 <- renderPlot({
                                        ggplot() +
                                            annotate("text", x = 0.5, y = 0.5,
                                                     label = "No mass spectrum available\n(non-MS detector - e.g. TCD, FID, ECD)",
                                                     size = 6, colour = "grey40") +
                                            theme_void() +
                                            xlim(0, 1) + ylim(0, 1)
                                    })
                                } else {

                                    # Normalize to max abu 100
                                        MS_out_1$intensity <- MS_out_1$intensity*100/max(MS_out_1$intensity)

                                    # Set ranges on mass spec (allow zooming in with selection), Sometimes brush returns Inf or -Inf if the range gets too small, so if that happens, just go up to the main view

                                        if (isolate(is.null(input$massSpectra_1_brush))) {
                                            MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
                                        } else {
                                            ## One conversion, not two, and through brushed_ms() so an
                                            ## empty coordmap cannot error the whole block away.
                                            ms1_sel <- isolate(brushed_ms(MS_out_1, input$massSpectra_1_brush))
                                            if (is.null(ms1_sel) || nrow(ms1_sel) == 0) {
                                                notify("Selection held no ions. Showing the full spectrum.")
                                                MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
                                            } else {
                                                MS1_low_x_limit <- min(ms1_sel$mz)
                                                MS1_high_x_limit <- max(ms1_sel$mz)
                                                in_win <- dplyr::filter(MS_out_1, mz >= MS1_low_x_limit & mz <= MS1_high_x_limit)$intensity
                                                MS1_high_y_limit <- if (length(in_win) > 0) max(in_win) + 8 else 110
                                            }
                                        }
                                        if (MS1_low_x_limit %in% c(Inf, -Inf) | MS1_high_x_limit %in% c(Inf, -Inf) | MS1_high_y_limit %in% c(Inf, -Inf)) {
                                            MS1_low_x_limit <- 0; MS1_high_x_limit <- 1200; MS1_high_y_limit <- 110
                                        }

                                    # Make plot

                                        output$massSpectra_1 <- renderPlot({

                                            ms1_plot <- ggplot() +
                                                geom_bar(
                                                    data = MS_out_1,
                                                    mapping = aes(x = mz, y = intensity),
                                                    stat = "identity", width = 0.1,
                                                    color = "black", fill = "grey"
                                                ) +
                                                theme_classic() +
                                                scale_x_continuous(expand = c(0,0)) +
                                                scale_y_continuous(expand = c(0,0)) +
                                                coord_cartesian(xlim = c(MS1_low_x_limit, MS1_high_x_limit), ylim = c(0, MS1_high_y_limit)) +
                                                geom_text(
                                                    data = # If it's less than 10 bars, label all. Otherwise, label 10 biggest ones
                                                        if (MS1_high_x_limit - MS1_low_x_limit >= 1) {
                                                            dplyr::filter(MS_out_1, mz > MS1_low_x_limit & mz < MS1_high_x_limit)[
                                                                order(
                                                                    dplyr::filter(MS_out_1, mz > MS1_low_x_limit & mz < MS1_high_x_limit)$intensity,
                                                                    decreasing = TRUE
                                                                )[1:10]
                                                            ,]
                                                        } else {
                                                            MS_out_1
                                                        },
                                                    mapping = aes(x = mz, y = intensity + 5, label = mz)
                                                )

                                            ## Kept so a mapping-less brush on this plot can be
                                            ## re-expressed against the panel's real bounds, exactly
                                            ## as the chromatogram's is.
                                            last_ms_plot <<- ms1_plot
                                            last_ms_size <<- suppressWarnings(as.numeric(c(
                                                session$clientData$output_massSpectra_1_width,
                                                session$clientData$output_massSpectra_1_height)))
                                            ## The m/z window changes from peak to peak, and the
                                            ## cached panel geometry carries it -- so it belongs in
                                            ## the cache key, or every spectrum after the first is
                                            ## read against the first one's axis.
                                            last_ms_window <<- suppressWarnings(as.numeric(
                                                c(MS1_low_x_limit, MS1_high_x_limit)))
                                            ms1_plot
                                        })
                                }
                            }

                        ## Extract and classify mass spectra on Shift+4

                            if (input$keypress == 36) {
                                undo_begin("classify peaks (Shift+4)"); on.exit(undo_end(), add = TRUE)
                                # Read in the peak table (ensure your file is in the correct format)
                                    peak_table <<- read.csv("peaks_monolist.csv", stringsAsFactors = FALSE)
                                    if(nrow(peak_table) == 0){
                                        notify("No peaks to classify.")
                                        return()
                                    }

                                ### somethign to do with the duplicate names
                                ### also all predicted peaks ID should be prefixed with "prediction"
                                ### perhaps it should only output predictions on things labelled "unknown"?

                                # Loop over each detected peak.
                                #
                                # The framed CSV is read once per RUN of peaks from the same sample,
                                # not once per peak. It used to be fread() inside this loop, so a
                                # sample with 40 peaks read its whole m/z cube 40 times -- the single
                                # dominant cost of Shift+4. Shift+Q now writes peaks_monolist.csv
                                # sorted by (sample, peak_start), so in practice that is one read per
                                # sample. Deliberately NOT re-sorted here: `predictions` is still
                                # matched to peak_table POSITIONALLY (finding A4, deferred), so
                                # changing the visiting order would change which peak gets which
                                # label. One cube is held at a time either way.
                                    ms_list <- list()
                                    loaded_file <- NULL; df <- NULL
                                    for(i in seq_len(nrow(peak_table))){
                                        # message(paste0("\n*getting spectrum for peak ", i))
                                        peak <- peak_table[i, ]
                                        # Get the corresponding .CDF.csv file
                                        cdf_file <- as.character(peak$path_to_cdf_csv)
                                        if(file.exists(cdf_file)){
                                            if (is.null(loaded_file) || !identical(loaded_file, cdf_file)) {
                                                df <- data.table::fread(cdf_file)
                                                loaded_file <- cdf_file
                                            }
                                            # Filter rows that fall within the peak’s retention time window
                                            df_subset <- df[df$rt >= peak$peak_start & df$rt <= peak$peak_end, ]
                                            if(nrow(df_subset) > 0){
                                                # Round m/z values and average intensities by m/z
                                                df_agg <- df_subset %>%
                                                    dplyr::mutate(mz = round(mz, 0)) %>%
                                                    dplyr::group_by(mz) %>%
                                                    dplyr::summarize(intensity = sum(intensity, na.rm = TRUE))
                                                # Tag the result with a peak identifier
                                                df_agg$peak_id <- peak$peak_ID
                                                df_agg$peak_unique_id <- i
                                                ms_list[[length(ms_list) + 1]] <- df_agg
                                            }
                                        }
                                    }

                                    if(length(ms_list) == 0){
                                        notify("Classification: no spectra could be extracted from the peaks.")
                                        return()
                                    }

                                # Combine all peaks into one data frame
                                    ms_data <<- dplyr::bind_rows(ms_list)

                                # Guard: single-detector data (e.g. GC-TCD) has no
                                # mass spectrum to classify, so stop here cleanly.
                                    if (all(ms_data$mz == 0)) {
                                        notify("Classification needs mass-spec data. Skipped.")
                                        return()
                                    }

                                    ms_data_reactive(ms_data)
                                    message("Mass spectra extracted successfully.")

                                # When the user clicks "Classify Spectra"
                                    ms_data <- ms_data_reactive()
                                    if(is.null(ms_data)){
                                        message("Please extract mass spectra first.")
                                        return()
                                    }

                                # Convert the long-format mass spectral data into a wide format
                                # so that each row corresponds to one peak and columns are m/z values.
                                ms_data %>%
                                    mutate(mz = round(mz)) %>%
                                    group_by(peak_unique_id, mz) %>%
                                    summarize(intensity = sum(intensity)) %>%
                                    bind_rows(data.frame(peak_unique_id = 0, mz = 0:1000, intensity = 0)) %>%
                                    pivot_wider(names_from = mz, values_from = intensity, values_fill = 0) -> ms_wide_raw
                                ms_wide_pre <- ms_wide_raw[ms_wide_raw$peak_unique_id != 0, ]
                                ms_wide_norm <- as_tibble(t(apply(ms_wide_pre[,-1], 1, function(x) (x / max(x)) * 100))) # normalize to 100
                                ms_wide <<- ms_wide_norm

                                message("ms_wide made and saved for inspection and plotting.")

                                # Classify peaks with the Python classifier (XGBoost, or sklearn
                                # HistGradientBoosting fallback). The trained model pkl must sit beside
                                # ms_classifier.py, which resolves it relative to its own __file__.
                                if (is.null(ms_classifier_path) || !file.exists(ms_classifier_path)) {
                                    message(paste0("Classification skipped - ms_classifier.py not found at: ",
                                                   ifelse(is.null(ms_classifier_path), "NULL", ms_classifier_path)))
                                    return()
                                }

                                # Export ms_wide (with peak_unique_id) for the Python classifier
                                ms_wide_export <- cbind(peak_unique_id = ms_wide_pre$peak_unique_id, ms_wide_norm)
                                write.csv(ms_wide_export, "ms_wide_for_prediction.csv", row.names = FALSE)

                                exit_code <- system(paste0("python3 ", shQuote(ms_classifier_path),
                                                           " predict ms_wide_for_prediction.csv predictions.csv"))
                                if (exit_code != 0 || !file.exists("predictions.csv")) {
                                    notify("Classification failed: the Python classifier did not run.")
                                    return()
                                }

                                # Read predictions back (columns: predicted_compound, confidence)
                                pred_df <- read.csv("predictions.csv", stringsAsFactors = FALSE)
                                predictions <<- pred_df$predicted_compound

                                message("predictions made:")

                                peak_table$peak_ID <- paste0("(predicted) ", predictions,
                                                             " [", round(pred_df$confidence * 100, 0), "%]")
                                writeMonolist(peak_table, "peaks_monolist.csv")

                                notify("Peaks classified.")
                            }

                    })

            }

        ## Call the app

            if ( jupyter == TRUE) {

                slot <- pick_shiny_port()
                if (is.null(slot)) {
                    cat("All available ports are in use. Please try again later.")
                } else {
                    cat(paste0("Connect at: ", slot$url))
                    runApp(shinyApp(ui = ui, server = server), host = "127.0.0.1", port = slot$port)
                }

            }

            if (jupyter == FALSE) {
                shinyApp(ui = ui, server = server)
            }


    }

