## Every course dataset, by the name the book uses. The small ones ship in data/ and are
## lazy-loaded (type "data"); the large or rarely used ones are downloaded once from the
## website's sample_data/ folder and cached (type "download").
.dataset_catalogue <- data.frame(
    name = c("algae_data", "alaska_lake_data", "solvents", "periodic_table", "fadb_sample",
             "periodic_table_subset", "ny_trees", "metabolomics_data", "wine_grape_data",
             "hawaii_aquifers", "beer_components", "wood_smoke", "unknown_smoke", "hops_components",
             "tequila_chemistry", "chemical_blooms", "metabolomics_unknown", "metabolomics_hospital",
             "metabolomics_hospital_unknown", "wine_quality", "lake_superior_shoreline",
             "rice_proteins", "mushrooms", "unknown_mushroom", "passenger_flows", "us_airports",
             "meteorite_fragments", "OSC_sequences"),
    file = c("algae_data.csv", "alaska_lake_data.csv", "solvents.csv", "per_table.csv", "fadb_sample.csv",
             "per_table_small.csv", "ny_trees.csv", "metabolomics_data.csv", "wine_grape_data.csv",
             "hawaii_aquifers.csv", "beer_components.csv", "wood_smoke_data.csv", "unknown_smoke.csv",
             "hops_components.csv", "tequila_chemistry.csv", "chemical_blooms.csv",
             "metabolomics_unknown.csv", "metabolomics_hospital.csv",
             "metabolomics_hospital_unknown.csv", "wine_quality.csv", "lake_superior_shoreline.csv",
             "rice_proteins.csv", "mushrooms.csv", "unknown_mushroom.csv", "passenger_flows.csv",
             "us_airports.csv", "meteorite_fragments.csv", "OSCs.fasta"),
    type = c("data", "data", "data", "data", "data", "data", "download", "data", "data", "data",
             "data", "data", "data", "data", "data", "data", "data", "data", "data", "data",
             "download", "data", "download", "data", "data", "data", "data", "download"),
    stringsAsFactors = FALSE
)

.dataset_base_url <- "https://thebustalab.github.io/phylochemistry/sample_data/"

#' Get a course dataset by name
#'
#' The small datasets used in the book ship with the package and are available directly by
#' name (for example `alaska_lake_data`). The large or rarely used ones (`ny_trees`,
#' `lake_superior_shoreline`, `mushrooms`, `OSC_sequences`) are downloaded from the course
#' website the first time you ask for them and kept in a cache folder, so later calls work
#' offline. `getDataset()` also returns the shipped datasets, so one call works for every name.
#'
#' @param name The dataset's name, as a string (for example `"ny_trees"`). Leave empty to list
#'   every available dataset.
#' @param refresh If `TRUE`, download the file again even if it is already cached.
#' @return A tibble (or, for `OSC_sequences`, a Biostrings `AAStringSet`). With no `name`, a
#'   data frame listing the datasets: `name`, `file`, `type` (`"data"` = shipped with the
#'   package, `"download"` = fetched and cached on first use) and `cached`.
#' @examples
#' getDataset()
#' \dontrun{
#' ny_trees <- getDataset("ny_trees")
#' }
#' @export
getDataset <- function(name = NULL, refresh = FALSE) {
    cat_ <- .dataset_catalogue
    if (is.null(name)) {
        cat_$cached <- cat_$type == "download" & file.exists(file.path(.dataset_cache_dir(), cat_$file))
        return(cat_)
    }
    if (!is.character(name) || length(name) != 1 || !name %in% cat_$name) {
        stop("Unknown dataset ", deparse(name), ". Run getDataset() to see the available names.", call. = FALSE)
    }
    row <- cat_[cat_$name == name, ]
    if (row$type == "data") {
        env <- new.env()
        utils::data(list = name, package = "phylochemistry", envir = env)
        return(env[[name]])
    }
    if (identical(name, "OSC_sequences")) .check_bioc("Biostrings", "getDataset")
    path <- .dataset_fetch(row$file, refresh = refresh)
    if (grepl("\\.fasta$", row$file)) return(Biostrings::readAAStringSet(path))
    out <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
    out <- tibble::as_tibble(as.data.frame(out))
    attr(out, "spec") <- NULL
    out
}

.dataset_cache_dir <- function() tools::R_user_dir("phylochemistry", which = "cache")

.dataset_fetch <- function(file, refresh = FALSE) {
    dir <- .dataset_cache_dir()
    path <- file.path(dir, file)
    if (file.exists(path) && !refresh) return(path)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(fileext = paste0(".", tools::file_ext(file)))
    on.exit(unlink(tmp), add = TRUE)
    message("Downloading ", file, " (first use only; cached in ", dir, ")")
    old <- options(timeout = max(600, getOption("timeout")))
    on.exit(options(old), add = TRUE)
    status <- tryCatch(utils::download.file(paste0(.dataset_base_url, file), tmp, mode = "wb", quiet = TRUE),
                       error = function(e) e)
    if (inherits(status, "error") || !identical(status, 0L) || !file.exists(tmp) || file.size(tmp) == 0) {
        stop("Could not download ", file, " from ", .dataset_base_url,
             ". Check your internet connection and try again.", call. = FALSE)
    }
    if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not write ", path, call. = FALSE)
    path
}
