# pca_visualizer.R -- compatibility wrapper.
#
# pcaVisualizer() is now part of the phylochemistry package:
#     library(phylochemistry); pcaVisualizer(data, ...)
#
# Sourcing this file installs the package if needed and attaches it.

local({
    op <- options(phylochemistry.compat.mode = "packages")
    on.exit(options(op), add = TRUE)
    source(getOption("phylochemistry.loader",
                     "https://thebustalab.github.io/phylochemistry/phylochemistry.R"))
})
