#### genomescope.R -- compatibility wrapper.
####
#### genomeScope() is now part of the bustalab package:
####     library(bustalab)
####
#### Sourcing this file installs bustalab if needed and attaches it.

local({
    op <- options(phylochemistry.compat.mode = "packages", phylochemistry.compat.bustalab = TRUE)
    on.exit(options(op), add = TRUE)
    source(getOption("phylochemistry.loader",
                     "https://thebustalab.github.io/phylochemistry/phylochemistry.R"))
})
