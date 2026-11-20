## Packages the Integrated Bioanalytical Chemistry book uses, declared in Imports so that installing
## phylochemistry installs them: students install only tidyverse and phylochemistry. The chapters
## attach them with library() or call them with pkg:: where they need them; the package itself never
## attaches them.
##
## They are referenced below, never called, so R CMD check counts them as used ("All declared Imports
## should be used"). Deliberately NOT @importFrom: an import loads the namespace at
## library(phylochemistry), and loading alone changes ggplot2 for the whole session: ggtern replaces
## print.ggplot, plot.ggplot and grid.draw.ggplot, ggside replaces `+.gg`, ggpp replaces two grid
## methods. A `::` reference here loads nothing until a chapter asks for the package.
##
##   cowplot (plot_grid), ggrepel, ggforce, viridis, RColorBrewer, ggpp (geom_table, geom_plot),
##   ggpmisc, ggside, ggdist, Hmisc (impute), performance + see (check_model), ggtern, maps +
##   mapproj (map_data, coord_map), geosphere (distHaversine), eulerr (vennAnalysis)
.book_packages <- function() {
    list(
        geosphere::distHaversine,
        ggdist::stat_halfeye,
        ggpmisc::stat_poly_eq,
        ggpp::geom_table,
        ggside::geom_xsidedensity,
        ggtern::ggtern,
        Hmisc::impute,
        mapproj::mapproject,
        maps::map,
        performance::check_model,
        RColorBrewer::brewer.pal,
        see::theme_lucid
    )
}
