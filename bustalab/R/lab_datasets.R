## Lab datasets, formerly loaded into the global environment by modules/datasets.R when
## `bustalab <- TRUE`. They are not bundled (busta_spectral_library is ~40 MB): each is
## downloaded from the same public URL on first use and cached on disk. The readers and
## column types are the ones datasets.R used.

.lab_datasets <- list(
    busta_spectral_library = list(
        url  = "https://thebustalab.github.io/phylochemistry/sample_data/busta_spectral_library_v1.csv",
        read = function(path) readr::read_csv(path, col_types = c(Compound_common_name = "c"))
    ),
    plant_phylogeny = list(
        url  = "https://thebustalab.github.io/data/plant_phylogeny.newick",
        read = function(path) ape::read.tree(path)
    ),
    plant_species = list(
        url  = "https://thebustalab.github.io/data/plant_species.csv",
        read = function(path) readr::read_csv(path, show_col_types = FALSE)
    )
)

#' Lab reference datasets (downloaded once, then cached)
#'
#' Fetches one of the lab's reference datasets from the lab website on first use,
#' caches the file under `tools::R_user_dir("bustalab", "cache")` and reads it with
#' the same reader and column types the old `bustalab <- TRUE` loader used.
#'
#' * `"busta_spectral_library"`: the lab mass-spectral library (CSV, about 40 MB).
#' * `"plant_phylogeny"`: a plant phylogeny (`ape::read.tree()` of a Newick file).
#' * `"plant_species"`: the plant species table (CSV).
#'
#' @param name Dataset name; omit it to list the available names.
#' @param refresh If `TRUE`, download again even when a cached copy exists.
#' @param cache_dir Where cached files are kept.
#' @return The dataset (a tibble, or a `phylo` tree for `"plant_phylogeny"`); with
#'   `name` missing, a character vector of dataset names.
#' @examples
#' getLabDataset()
#' \dontrun{
#' busta_spectral_library <- getLabDataset("busta_spectral_library")
#' }
#' @export
getLabDataset <- function(name, refresh = FALSE,
                          cache_dir = tools::R_user_dir("bustalab", "cache")) {
    if (missing(name)) return(names(.lab_datasets))
    if (!is.character(name) || length(name) != 1 || !name %in% names(.lab_datasets)) {
        stop("Unknown lab dataset. Choose one of: ",
             paste(names(.lab_datasets), collapse = ", "), call. = FALSE)
    }
    spec <- .lab_datasets[[name]]
    path <- file.path(cache_dir, basename(spec$url))
    if (refresh || !file.exists(path)) {
        dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
        tmp <- tempfile(tmpdir = cache_dir, fileext = ".part")
        on.exit(unlink(tmp), add = TRUE)
        old_timeout <- options(timeout = max(600, getOption("timeout")))
        on.exit(options(old_timeout), add = TRUE)
        utils::download.file(spec$url, tmp, mode = "wb", quiet = TRUE)
        if (!file.rename(tmp, path)) stop("Could not write the cached copy to ", path, call. = FALSE)
    }
    spec$read(path)
}
