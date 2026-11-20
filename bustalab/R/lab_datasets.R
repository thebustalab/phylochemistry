## Moved verbatim from modules/datasets.R (phylochemistry website tree).

## Busta lab specific datasets

    if (exists("bustalab")) {
        
        if (bustalab == TRUE) {

            busta_spectral_library <- read_csv("https://thebustalab.github.io/phylochemistry/sample_data/busta_spectral_library_v1.csv", col_types = c(Compound_common_name = "c"))
            plant_phylogeny <- read.tree("https://thebustalab.github.io/data/plant_phylogeny.newick")
            plant_species <- read_csv("https://thebustalab.github.io/data/plant_species.csv", show_col_types = FALSE)
        
        }
    }


