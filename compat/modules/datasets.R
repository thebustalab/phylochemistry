# datasets.R -- compatibility wrapper.
#
# The course datasets now ship with the phylochemistry package:
#     library(phylochemistry)       # alaska_lake_data etc. are then available by name
#     getDataset()                  # lists every dataset; getDataset("ny_trees") fetches a big one
#
# Sourcing this file still puts them in your workspace, as it used to. It installs the
# package first if needed (via the main phylochemistry.R loader).

local({
    op <- options(phylochemistry.compat.mode = "datasets")
    on.exit(options(op), add = TRUE)
    source(getOption("phylochemistry.loader",
                     "https://thebustalab.github.io/phylochemistry/phylochemistry.R"))
})
