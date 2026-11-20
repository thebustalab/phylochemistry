#' phylochemistry: chemical, statistical and sequence analyses
#'
#' Teaching and research toolkit for the Integrated Bioanalytical Chemistry course book and the
#' Busta lab. See the package index for the analysis verbs and datasets.
#'
#' @keywords internal
#' @importFrom ape cophenetic.phylo drop.tip pic
#' @importFrom base64enc dataURI
#' @importFrom data.table fread
#' @importFrom dplyr all_of arrange as_tibble bind_cols everything filter
#' @importFrom dplyr full_join group_by left_join mutate n pull right_join select
#' @importFrom dplyr slice summarize sym tibble ungroup
#' @importFrom FNN get.knnx
#' @importFrom ggforce geom_link
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggplot2 aes aes_string arrow coord_cartesian coord_equal coord_fixed coord_polar CoordPolar
#' @importFrom ggplot2 element_blank element_line element_text facet_grid facet_wrap fortify geom_col geom_histogram
#' @importFrom ggplot2 geom_hline geom_line geom_point geom_rect geom_segment geom_text geom_tile geom_vline
#' @importFrom ggplot2 ggplot ggproto ggsave guide_colorbar guides is.ggplot labs layer
#' @importFrom ggplot2 scale_color_manual scale_fill_gradient scale_fill_gradient2 scale_fill_manual scale_size_continuous scale_x_continuous scale_y_continuous scale_y_discrete
#' @importFrom ggplot2 stat stat_summary theme theme_bw theme_classic theme_void unit
#' @importFrom ggtree geom_inset ggtree
#' @importFrom graphics image
#' @importFrom grDevices col2rgb dev.off png rgb
#' @importFrom httr content content_type_json GET POST timeout
#' @importFrom magrittr %>% set_colnames
#' @importFrom parsnip fit rand_forest set_args set_engine set_mode
#' @importFrom phangorn parsimony phyDat
#' @importFrom phytools fastAnc
#' @importFrom readr cols read_csv write_csv
#' @importFrom recipes recipe
#' @importFrom rlang .data warn
#' @importFrom rsample vfold_cv
#' @importFrom Rtsne Rtsne
#' @importFrom shiny actionButton br brushedPoints brushOpts column conditionalPanel div eventReactive
#' @importFrom shiny fluidPage fluidRow h2 h3 h4 helpText hr HTML
#' @importFrom shiny icon isolate mainPanel navbarPage need observe observeEvent plotOutput
#' @importFrom shiny radioButtons reactive reactiveVal reactiveValues renderPlot renderPrint renderText renderUI
#' @importFrom shiny req runApp selectInput selectizeInput shinyApp sidebarLayout sidebarPanel sliderInput
#' @importFrom shiny splitLayout stopApp tabPanel tabsetPanel tags textAreaInput textInput textOutput
#' @importFrom shiny titlePanel uiOutput updateRadioButtons validate verbatimTextOutput verticalLayout
#' @importFrom shinyjs hide useShinyjs
#' @importFrom stats approx as.formula binomial coefficients cor glm hclust lm
#' @importFrom stats median nls p.adjust pf predict pt sd setNames
#' @importFrom stringr str_count str_split str_trim
#' @importFrom tidyr drop_na pivot_longer pivot_wider unite
#' @importFrom tune collect_metrics finalize_workflow select_best tune tune_grid
#' @importFrom umap umap
#' @importFrom utils askYesNo combn file.edit head read.csv read.delim read.table setTxtProgressBar
#' @importFrom utils tail txtProgressBar unzip URLencode write.csv write.table
#' @importFrom viridis scale_color_viridis scale_fill_viridis
#' @importFrom workflows add_model add_recipe workflow
"_PACKAGE"

## Column names used through non-standard evaluation (dplyr/ggplot2), and the few names some
## interactive functions deliberately read from or write to the user's session.
utils::globalVariables(c(
    ".", "abundance", "alignment_directory_path", "amino_acid", "analyte", "Analyte",
    "atom_number", "B", "bond_direction", "bond_end_x", "bond_end_y", "bond_start_x",
    "bond_start_y", "bond_type", "build_tree_substitutions", "cc", "cluster", "cluster_k",
    "cluster_threshold", "color_n", "color_swatch", "conservation_percent", "corner_points_x",
    "corner_points_y", "dat", "Dim_1", "Dim_2", "distance", "doi", "element_or_group", "end",
    "Expect", "fortified_tree", "G", "gap_percent", "Gaps_Percent", "Gaps_Ratio", "Genus_species",
    "goi_region_name", "GPT_response", "GPT_temperature", "Greater_Than_Line", "Hit_Name", "hjust",
    "Identities_Line", "Identities_Percent", "Identities_Ratio", "image_as_df", "image_calibrated",
    "image_points", "imRed", "index", "intensity", "is_genus_in_tree", "is_species_in_tree",
    "isTip", "knn", "label", "metric", "mode_reference_value_B", "mode_reference_value_G",
    "mode_reference_value_R", "molecule_component", "mutate", "mz", "n_clusters", "name", "node",
    "peak_number", "percent_variance_explained", "photo_metadata", "picker", "pickerOutput",
    "plots", "position", "principal_component", "profile_data", "R", "Raw_Score", "ref",
    "reference_pixel_output", "relative_abundance", "renderPicker", "residue", "rt", "Sample",
    "sample_1", "sample_2", "sample_pixel_data", "sample_unique_id", "sample_unique_ID",
    "sample_unique_ID_sample_1", "sample_unique_ID_sample_2", "Score", "Score_Line", "sd", "shell",
    "start", "state", "status", "sub_data", "swatch_data", "text", "text_horiz_offset",
    "text_size", "text_vert_offset", "tip_length_xmax", "tip_length_xmin", "tip_length_ymax",
    "tip_length_ymin", "trait", "value", "Value", "values", "verbose", "vjust", "x", "x_position",
    "x_round", "xint", "xmax", "xmin", "y", "y_position", "yb", "yint", "yintercept", "ymax",
    "ymin"
))

## phylochemistry.R set these two options when it was sourced; keep them for package users, but only
## when the user has not set them.
.onLoad <- function(libname, pkgname) {
    defaults <- list(readr.show_progress = FALSE, dplyr.summarise.inform = FALSE)
    unset <- !(names(defaults) %in% names(options()))
    if (any(unset)) options(defaults[unset])
    invisible()
}

## Bioconductor packages cannot be installed by rlang::check_installed(); say how instead.
.check_bioc <- function(pkgs, fn) {
    missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop(fn, "() needs the Bioconductor package", if (length(missing) > 1) "s", " ",
             paste(missing, collapse = ", "), ". Install with:\n  ",
             "install.packages(\"BiocManager\"); BiocManager::install(c(",
             paste0("\"", missing, "\"", collapse = ", "), "))", call. = FALSE)
    }
    invisible(TRUE)
}

## ranger is parsnip's random-forest engine: buildModel() reaches it through parsnip, never by name.
.phylochemistry_engines <- function() ranger::ranger

