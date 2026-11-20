##############################################
## PHYLOCHEMISTRY -- compatibility loader   ##
##############################################
##
## phylochemistry is now an R package. New code:
##
##     install.packages("phylochemistry", repos = c("https://thebustalab.r-universe.dev",
##                                                  BiocManager::repositories()))
##     library(tidyverse); library(phylochemistry)
##
## This file keeps the old one-liner working:
##
##     source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")
##
## It installs the package if needed, attaches it with the packages the old script attached,
## and puts the course datasets in your workspace as before. Set `bustalab <- TRUE` first for
## the lab package, `phylochemistry_scope <- "full"` for the sequence / phylogenetics chapters.
##
## Options (for testing or offline mirrors):
##   phylochemistry.repos   repositories to install from (default: r-universe, Bioconductor, CRAN)
## Set by the thin wrappers at the old module URLs, not meant for users:
##   phylochemistry.compat.mode      "legacy" (default), "packages" or "datasets"
##   phylochemistry.compat.bustalab  TRUE to load the lab package without a global `bustalab`

local({

    ###### Install/load preflight (verbatim from the pre-package phylochemistry.R)

    ## Install/load preflight — stop early with a clear, actionable message
    ## instead of a cryptic error partway through installing ~280 packages.
    ## Targets the two things students trip on most: an old / pre-existing R,
    ## and a package library that R has no permission to write to.

        .phylochemistry_min_R <- "4.1.0"

        .phylochemistry_preflight <- function(min_R = .phylochemistry_min_R) {

            ## 1. R version
            if (getRversion() < min_R) {
                stop(paste0(
                    "\n\n[phylochemistry] Your R is version ", getRversion(),
                    ", but phylochemistry needs R ", min_R, " or newer.\n",
                    "  -> Install the latest R from https://cloud.r-project.org\n",
                    "  -> Then FULLY quit and reopen RStudio so it uses the new R, and source phylochemistry again.\n",
                    "     (you can check your version any time with:  R.version.string )\n"
                ), call. = FALSE)
            }

            ## 2. A writable package library must exist on the path
            libs <- .libPaths()
            if (!any(file.access(libs, mode = 2) == 0)) {
                userlib <- path.expand(Sys.getenv("R_LIBS_USER"))
                fixed <- FALSE
                if (nzchar(userlib)) {
                    if (!dir.exists(userlib))
                        try(dir.create(userlib, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
                    if (dir.exists(userlib) && file.access(userlib, mode = 2) == 0) {
                        .libPaths(c(userlib, libs))
                        message("[phylochemistry] No writable library was on the path; ",
                                "using your personal one:\n  ", userlib)
                        fixed <- TRUE
                    }
                }
                if (!fixed) {
                    stop(paste0(
                        "\n\n[phylochemistry] R has nowhere it can install packages -- every folder\n",
                        "on your library path is read-only:\n  ",
                        paste(libs, collapse = "\n  "), "\n",
                        "Create a personal library by running these two lines, then source phylochemistry again:\n",
                        "  dir.create(Sys.getenv('R_LIBS_USER'), recursive = TRUE)\n",
                        "  .libPaths(Sys.getenv('R_LIBS_USER'))\n"
                    ), call. = FALSE)
                }
            }

            ## 3. Shared machines: a writable library whose PACKAGES belong to
            ##    someone else. On the lab JupyterHub, /usr/local/lib/R/site-library
            ##    is world-writable but its package folders are owned by root and by
            ##    other users, so install.packages() can add a NEW package yet cannot
            ##    UPGRADE an existing one. That surfaces as the baffling
            ##    "installation of package 'ggplot2' had non-zero exit status" while
            ##    the folder looks perfectly writable. Fix: put the user's own library
            ##    FIRST, so upgrades land in $HOME and shadow the shared copy.
            first <- .libPaths()[1]
            if (file.access(first, mode = 2) == 0) {
                foreign <- tryCatch({
                    pkgs <- list.dirs(first, recursive = FALSE)
                    length(pkgs) > 0 && any(file.access(pkgs, mode = 2) != 0)
                }, error = function(e) FALSE)
                userlib <- path.expand(Sys.getenv("R_LIBS_USER"))
                if (foreign && nzchar(userlib) &&
                    !identical(path.expand(first), userlib)) {
                    if (!dir.exists(userlib))
                        try(dir.create(userlib, recursive = TRUE, showWarnings = FALSE), silent = TRUE)
                    if (dir.exists(userlib) && file.access(userlib, mode = 2) == 0) {
                        .libPaths(c(userlib, .libPaths()))
                        message("[phylochemistry] ", first, " holds packages owned by other users, ",
                                "which R cannot upgrade.\n  Installing to your personal library instead:\n  ", userlib)
                    }
                }
            }

            invisible(TRUE)
        }

        ## Lightweight, NON-fatal connectivity probe; only consulted when
        ## something actually needs downloading. A blocking proxy can fail this
        ## yet still allow install.packages(), so it only warns.
        .phylochemistry_online <- function(test_url = "https://cloud.r-project.org") {
            old <- options(timeout = 10); on.exit(options(old), add = TRUE)
            suppressWarnings(tryCatch({
                con <- url(test_url, open = "rb", method = "libcurl")
                on.exit(close(con), add = TRUE)
                readBin(con, "raw", 1L)
                TRUE
            }, error = function(e) FALSE))
        }

        .phylochemistry_preflight()

    ###### Settings

        ## Raise these when a script starts relying on a newer package release.
        min_version <- c(phylochemistry = "1.0.0", bustalab = "0.1.0")
        universe    <- "https://thebustalab.r-universe.dev"

        g <- globalenv()
        compat_mode <- getOption("phylochemistry.compat.mode", "legacy")
        lab <- (exists("bustalab", envir = g) && isTRUE(get("bustalab", envir = g))) ||
            isTRUE(getOption("phylochemistry.compat.bustalab"))
        scope <- if (exists("phylochemistry_scope", envir = g)) {
            get("phylochemistry_scope", envir = g)
        } else if (lab) "full" else "core"
        full <- identical(scope, "full") || lab
        skip_packages <- exists("packages", envir = g)

        ## The packages the old script attached, in the old order (lab CRAN set, core,
        ## extended, then the Bioconductor sets). Old code relies on them being attached
        ## (cowplot::plot_grid, ggrepel, ape::read.tree, ... with no library() call).
        lab_cran   <- c("minpack.lm", "XML", "DT", "rhandsontable")
        core_cran  <- c("tidyverse", "dplyr", "stringr", "data.table", "jsonlite", "progress",
                        "ggrepel", "cowplot", "ggforce", "ggside", "ggdist", "ggpmisc", "viridis", "RColorBrewer",
                        "ggnewscale", "patchwork", "gridExtra", "grid", "rstatix", "multcompView",
                        "performance", "Hmisc", "FactoMineR", "ape", "ips", "fpc", "Rtsne", "umap",
                        "shipunov", "Rfast", "recipes", "rsample", "workflows", "tune", "yardstick",
                        "parsnip", "network", "ggnetwork", "shiny", "shinyjs", "shinythemes", "httpuv",
                        "later", "parallel", "BiocManager")
        ## installed, never attached (ggtern hijacks ggplot2 printing, plyr masks dplyr)
        core_install <- c("ggtern", "maps", "mapproj", "geosphere", "eulerr", "plyr", "see",
                          "ranger", "igraph", "emmeans", "vegan")
        ext_cran   <- c("httr", "rentrez", "RCurl", "base64enc", "phangorn", "seqinr", "picante",
                        "phytools", "bio3d", "png", "exifr", "rhdf5", "googlesheets4", "googledrive",
                        "pracma", "lubridate")
        ## msa: alignSequences(); a Suggests of the package that the old file only got via bustalab
        ext_install <- c("msa")
        lab_bioc   <- c("xcms", "DESeq2", "msa", "rtracklayer", "GenomicRanges", "GenomicFeatures",
                        "Rsamtools", "shinythemes", "FNN")
        core_bioc  <- c("ggtree")
        ext_bioc   <- c("ggtreeExtra", "Biostrings", "GenomeInfoDb")

    ###### Helpers

        is_installed <- function(pkgs) {
            pkgs %in% basename(find.package(pkgs, quiet = TRUE))
        }
        outdated <- function(pkg) {
            v <- tryCatch(utils::packageVersion(pkg), error = function(e) NULL)
            is.null(v) || v < min_version[[pkg]]
        }
        repos <- function() {
            r <- getOption("phylochemistry.repos")
            if (!is.null(r)) return(r)
            cran <- getOption("repos")[["CRAN"]]
            if (is.null(cran) || identical(cran, "@CRAN@")) cran <- "https://cloud.r-project.org"
            if (!is_installed("BiocManager")) {
                message("[phylochemistry] Installing BiocManager (finds the right Bioconductor release for your R)...")
                utils::install.packages("BiocManager", repos = cran)
            }
            bioc <- tryCatch({
                op <- options(repos = c(CRAN = cran)); on.exit(options(op), add = TRUE)
                b <- suppressMessages(BiocManager::repositories())
                b[names(b) != "CRAN"]
            }, error = function(e) c(BioCsoft = "https://bioc.r-universe.dev"))
            c(thebustalab = universe, bioc, CRAN = cran)
        }
        install <- function(pkgs, what) {
            pkgs <- unique(pkgs)
            if (length(pkgs) == 0) return(invisible())
            if (is.null(getOption("phylochemistry.repos")) && !.phylochemistry_online()) {
                warning(
                    "[phylochemistry] Can't reach CRAN (https://cloud.r-project.org). ",
                    "If the install below fails, check your internet connection or campus VPN/proxy and try again.",
                    call. = FALSE, immediate. = TRUE
                )
            }
            message("[phylochemistry] Installing ", what, ": ", paste(pkgs, collapse = ", "),
                    "\n  (one time only; this can take several minutes)")
            op <- options(install.packages.compile.from.source = "never")
            on.exit(options(op), add = TRUE)
            ## An unreachable repository errors outright; report it and let the checks below
            ## give the student the plain-English message.
            tryCatch(utils::install.packages(pkgs, repos = repos(), dependencies = NA),
                     error = function(e) warning("[phylochemistry] Installing failed: ", conditionMessage(e),
                                                 call. = FALSE, immediate. = TRUE))
        }
        attach_quietly <- function(pkgs) {
            failed <- character()
            for (p in pkgs) {
                ok <- tryCatch({
                    suppressWarnings(suppressPackageStartupMessages(
                        library(p, character.only = TRUE)))
                    TRUE
                }, error = function(e) FALSE)
                if (!ok) failed <- c(failed, p)
            }
            failed
        }

    ###### 1. Make sure the packages are installed

        message("[phylochemistry] This loader now installs the phylochemistry package; ",
                "new code can use library(phylochemistry).")

        ours <- c("phylochemistry", if (lab) "bustalab")
        stale <- ours[vapply(ours, outdated, logical(1))]
        loaded_stale <- stale[vapply(stale, isNamespaceLoaded, logical(1))]

        legacy <- identical(compat_mode, "legacy") && !skip_packages
        if (legacy) {
            attach_cran <- unique(c(if (lab) lab_cran, core_cran, if (full) ext_cran))
            attach_bioc <- unique(c(if (lab) lab_bioc, core_bioc, if (full) ext_bioc))
            attach_set  <- unique(c(attach_cran, attach_bioc))
            extras <- unique(c(attach_set, core_install, if (full) ext_install))
        } else {
            extras <- character()
        }
        missing_extras <- extras[!is_installed(extras)]

        install(stale, if (lab) "the phylochemistry / bustalab packages" else "the phylochemistry package")
        install(missing_extras, "packages the old phylochemistry script loaded")

        for (p in ours) {
            if (outdated(p)) {
                stop("\n\n[phylochemistry] The '", p, "' package (version ", min_version[[p]],
                     " or newer) could not be installed. Scroll up for the first error, then try:\n",
                     "  install.packages(\"", p, "\", repos = c(\"", universe, "\", BiocManager::repositories()))\n",
                     call. = FALSE)
            }
        }
        if (length(loaded_stale) > 0) {
            warning("[phylochemistry] ", paste(loaded_stale, collapse = ", "),
                    " was updated while an older copy was loaded. Restart R (Session > Restart R) ",
                    "and source this file again.", call. = FALSE, immediate. = TRUE)
        }

    ###### 2. Attach

        if (skip_packages) {
            message("Object 'packages' exists, not loading phylochemistry packages...")
        }

        if (!identical(compat_mode, "datasets")) {
            if (legacy) {
                message("Loading packages...")
                failed <- attach_quietly(attach_set)
                if (length(failed) > 0) {
                    warning("[phylochemistry] Could not load: ", paste(failed, collapse = ", "),
                            ". Functions that need them will fail.", call. = FALSE, immediate. = TRUE)
                }
            }
            suppressPackageStartupMessages(library("phylochemistry", character.only = TRUE))
            if (lab) suppressPackageStartupMessages(library("bustalab", character.only = TRUE))
        }

        ## The old script made dplyr's verbs win over the masks from Hmisc (summarize),
        ## AnnotationDbi (select), IRanges (slice) etc. by assigning them in the workspace.
        ## Here they sit in a small attached environment on top of the search path instead,
        ## so they survive rm(list = ls()) and leave the workspace clean.
        if (identical(compat_mode, "legacy")) {
            while ("phylochemistry:compat" %in% search()) detach("phylochemistry:compat", character.only = TRUE)
            verbs <- new.env()
            for (v in c("mutate", "summarize", "filter", "select", "slice")) {
                assign(v, getExportedValue("dplyr", v), envir = verbs)
            }
            base::attach(verbs, pos = 2L, name = "phylochemistry:compat", warn.conflicts = FALSE)
        }

    ###### 3. Datasets in the workspace, as modules/datasets.R used to do

        if (!identical(compat_mode, "packages")) {
            shipped <- utils::data(package = "phylochemistry")$results[, "Item"]
            utils::data(list = shipped, package = "phylochemistry", envir = g)
            for (obj in c("cont_1", "cont_2", "discrete_palette")) {
                assign(obj, getExportedValue("phylochemistry", obj), envir = g)
            }
            ## Big or rarely used sets: downloaded (and cached) only when first touched.
            lazy <- c("ny_trees", "lake_superior_shoreline", "mushrooms",
                      if (is_installed("Biostrings")) "OSC_sequences")
            for (nm in lazy) {
                eval(substitute(delayedAssign(NM, phylochemistry::getDataset(NM),
                                              eval.env = baseenv(), assign.env = g),
                                list(NM = nm)))
            }
            if (lab) {
                for (nm in c("busta_spectral_library", "plant_phylogeny", "plant_species")) {
                    eval(substitute(delayedAssign(NM, bustalab::getLabDataset(NM),
                                                  eval.env = baseenv(), assign.env = g),
                                    list(NM = nm)))
                }
            }
        }

        options(readr.show_progress = FALSE, dplyr.summarise.inform = FALSE)

        invisible(NULL)
})
