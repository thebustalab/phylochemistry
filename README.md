# phylochemistry

R tools for chemical, statistical and sequence analyses, used in the Busta lab and in the
*Integrated Bioanalytical Chemistry* course book.

This repository holds two R packages:

- **`phylochemistry/`**: the student and general-purpose toolkit (data visualization helpers,
  statistics verbs, matrix analysis, modelling, sequence and phylogenetic tools, course datasets).
- **`bustalab/`**: Busta lab extras (GC-MS analysis app, genome size estimation, language-model
  and DGX clients, lab datasets). Loading it also loads `phylochemistry`.

## Install

```r
install.packages("phylochemistry", repos = c("https://thebustalab.r-universe.dev", "https://cloud.r-project.org"))
```

Then, in each session:

```r
library(tidyverse)
library(phylochemistry)
```

Lab members install `bustalab` the same way and use `library(bustalab)`.

The 2021 version of this repository is preserved on the `legacy-2021` branch.
