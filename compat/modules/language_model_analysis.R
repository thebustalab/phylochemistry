# language_model_analysis.R -- compatibility wrapper.
#
# searchNCBI(), searchPubMed(), embedText(), generateText() and embedAminoAcids() are now
# part of the phylochemistry package:
#     library(phylochemistry)
# The lab's dgx client (dgxHealth(), dgxGenerate(), ...) is in the lab package:
#     library(bustalab)
#
# Sourcing this file installs and attaches phylochemistry, plus bustalab when you have set
# `bustalab <- TRUE` (as the old file did through modules/dgx.R).

local({
    op <- options(phylochemistry.compat.mode = "packages")
    on.exit(options(op), add = TRUE)
    source(getOption("phylochemistry.loader",
                     "https://thebustalab.github.io/phylochemistry/phylochemistry.R"))
})
