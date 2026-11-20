## Every shipped dataset loads under its book name, as a plain tibble with the expected columns.
## Expected columns are the header of the matching sample_data/ CSV on the website.

expected_columns <- list(
    alaska_lake_data = c(c("lake", "park", "water_temp", "pH", "element", "mg_per_L", "element_type")),
    algae_data = c(c("replicate", "algae_strain", "harvesting_regime", "chemical_species", "abundance")),
    beer_components = c(c("ingredient", "replicate", "analyte", "analyte_class", "abundance", "analyte_odor")),
    chemical_blooms = c(c("label", "Alkanes", "Sec_Alcohols", "Others", "Fatty_acids", "Alcohols", "Triterpenoids", "Ketones", "Other_compounds", "Aldehydes")),
    fadb_sample = c(c("fatty_acid", "Agonandra_brasiliensis", "Agonandra_silvatica", "Agonandra_excelsa", "Heisteria_silvianii", "Malania_oleifera", "Ximenia_americana", "Ongokea_gore", "Comandra_pallida", "Buckleya_distichophylla", "Nuytsia_floribunda")),
    hawaii_aquifers = c(c("aquifer_code", "well_name", "longitude", "latitude", "analyte", "abundance")),
    hops_components = c(c("hop_variety", "hop_origin", "hop_brewing_usage", "hop_aroma", "total_oil", "alpha_acids", "beta_acids", "humulene", "myrcene", "humulone", "caryophyllene", "farnesene")),
    metabolomics_data = c(c("patient_number", "patient_status", "Pyruvate", "MethylSuccinate", "Phosphoglyceric Acid", "Glucose 1-phosphate", "Malonic Acid", "Fumaric Acid", "Alpha-Ketoglutaric Acid", "Sarcosine", "Cadaverine", "5-Aminovaleric Acid", "Pipecolate", "Tyramine", "1-Methylhistidine", "isoValeric Acid", "2-Hydroxyisovaleric Acid", "3-Methyl-2-Oxovaleric Acid", "2-Aminoisobutyric acid", "Creatinine", "Betaine", "Taurine", "1-Methylhistamine", "Pyroglutamic Acid", "Creatine", "4-Hydroxyproline", "Homocysteine", "Cysteamine", "Carnitine", "Glycerate", "N-AcetylGlycine", "Guanidinoacetate", "Benzoic acid", "Citraconic Acid", "Glutaric Acid", "Oxaloacetate", "Anthranilate", "2-Hydroxyglutarate", "2-Aminoadipate", "Homogentisate", "Hydroxyphenylpyruvic acid", "Kynurenate", "Xanthurenate", "Pantothenic acid", "L-Kynurenine", "3-Hydroxykynurenine", "Cytidine", "Uridine", "1-Methyladenosine", "Hypoxanthine", "Oxypurinol", "Xanthine", "Orotate", "Allantoin", "Urate", "Xanthosine", "AMP", "IMP", "ADP", "Adenylosuccinate", "Choline", "Dimethylglycine", "Niacinamide", "Methylmalonate", "Ascorbic Acid", "4-Pyridoxic acid", "Biotin", "Pyridoxal 5'-phosphate", "Trimethylamine", "Trimethylamine N-oxide", "Shikimic Acid", "Hippuric Acid", "4-Hydroxybutyrate", "Acetylcholine", "Deoxycarnitine", "13-HODE", "18-Hydroxyarachidonic acid", "Glycochenodeoxycholate", "Acetylcarnitine", "Azelaic Acid", "Myristic Acid", "Margaric Acid", "Linoleic Acid", "Linolenic Acid", "Mannose", "Sorbitol", "Glyceraldehyde", "Erythrose", "myoinositol", "Fructose", "D-Glucoronic acid", "N-Acetylneuraminic acid", "Lactose", "Sucrose", "Ornithine", "Citrulline", "Propionate", "Oxalic Acid", "Acetoacetate", "Adipic Acid", "Ribose-5-P", "Geranyl-PP", "Glycine", "Alanine", "Serine", "Proline", "Valine", "Threonine", "Leucine", "iso-Leucine", "Asparagine", "Aspartic Acid", "Glutamine", "Lysine", "Glutamic acid", "Methionine", "Histidine", "Phenylalanine", "Arginine", "Tyrosine", "Tryptophan", "Cystine", "lactate", "3-Hydroxybenzoic acid", "Succinate", "Glucose")),
    metabolomics_hospital = c(c("patient_number", "patient_status", "metabolite", "concentration")),
    metabolomics_hospital_unknown = c(c("patient_number", "metabolite", "concentration")),
    metabolomics_unknown = c(c("Pyruvate", "MethylSuccinate", "Phosphoglyceric Acid", "Glucose 1-phosphate", "Malonic Acid", "Fumaric Acid", "Alpha-Ketoglutaric Acid", "Sarcosine", "Cadaverine", "5-Aminovaleric Acid", "Pipecolate", "Tyramine", "1-Methylhistidine", "isoValeric Acid", "2-Hydroxyisovaleric Acid", "3-Methyl-2-Oxovaleric Acid", "2-Aminoisobutyric acid", "Creatinine", "Betaine", "Taurine", "1-Methylhistamine", "Pyroglutamic Acid", "Creatine", "4-Hydroxyproline", "Homocysteine", "Cysteamine", "Carnitine", "Glycerate", "N-AcetylGlycine", "Guanidinoacetate", "Benzoic acid", "Citraconic Acid", "Glutaric Acid", "Oxaloacetate", "Anthranilate", "2-Hydroxyglutarate", "2-Aminoadipate", "Homogentisate", "Hydroxyphenylpyruvic acid", "Kynurenate", "Xanthurenate", "Pantothenic acid", "L-Kynurenine", "3-Hydroxykynurenine", "Cytidine", "Uridine", "1-Methyladenosine", "Hypoxanthine", "Oxypurinol", "Xanthine", "Orotate", "Allantoin", "Urate", "Xanthosine", "AMP", "IMP", "ADP", "Adenylosuccinate", "Choline", "Dimethylglycine", "Niacinamide", "Methylmalonate", "Ascorbic Acid", "4-Pyridoxic acid", "Biotin", "Pyridoxal 5'-phosphate", "Trimethylamine", "Trimethylamine N-oxide", "Shikimic Acid", "Hippuric Acid", "4-Hydroxybutyrate", "Acetylcholine", "Deoxycarnitine", "13-HODE", "18-Hydroxyarachidonic acid", "Glycochenodeoxycholate", "Acetylcarnitine", "Azelaic Acid", "Myristic Acid", "Margaric Acid", "Linoleic Acid", "Linolenic Acid", "Mannose", "Sorbitol", "Glyceraldehyde", "Erythrose", "myoinositol", "Fructose", "D-Glucoronic acid", "N-Acetylneuraminic acid", "Lactose", "Sucrose", "Ornithine", "Citrulline", "Propionate", "Oxalic Acid", "Acetoacetate", "Adipic Acid", "Ribose-5-P", "Geranyl-PP", "Glycine", "Alanine", "Serine", "Proline", "Valine", "Threonine", "Leucine", "iso-Leucine", "Asparagine", "Aspartic Acid", "Glutamine", "Lysine", "Glutamic acid", "Methionine", "Histidine", "Phenylalanine", "Arginine", "Tyrosine", "Tryptophan", "Cystine", "lactate", "3-Hydroxybenzoic acid", "Succinate", "Glucose")),
    meteorite_fragments = c(c("fragment", "find_site", "iron_percent", "density_g_per_cm3", "nickel_percent", "silica_percent")),
    passenger_flows = c(c("origin", "destination", "passengers")),
    periodic_table = c(c("atomic_number", "element_name", "atomic_symbol", "group_number", "period", "atomic_mass_rounded", "melting_point_C", "boiling_point_C", "state_at_RT", "density_g_per_mL", "electronegativity_pauling", "first_ionization_poten_eV", "second_ionization_poten_eV", "third_ionization_poten_eV", "electron_affinity_eV", "atomic_radius_ang", "ionic_radius_ang", "covalent_radius_ang", "atomic_volume_cm3_per_mol", "crystal_structure", "electrical_conductivity_mho_per_cm", "specific_heat_J_per_g_K", "heat_of_fusion_kJ_per_mol", "heat_of_vaporization_kJ_per_mol", "thermal_conductivity_W_per_m_K", "mg_per_kg_in_Earths_crust", "mg_per_L_in_seawater", "percent_human_body_mass", "Notes", "color", "structure", "hardness_mohs", "hydride", "oxide", "chloride", "polarizability_A_cubed", "heat_atomization_kJ_per_mol", "year_discovered", "rel_abund_solar_system_log", "abundance_earths_crust_log", "cost_pure_dollar_per_100g")),
    periodic_table_subset = c(c("atomic_number", "atomic_symbol", "group_number", "period", "atomic_mass_rounded", "melting_point_C", "boiling_point_C", "state_at_RT", "density_g_per_mL", "electronegativity_pauling", "first_ionization_poten_eV", "second_ionization_poten_eV", "third_ionization_poten_eV", "electron_affinity_eV", "atomic_radius_ang", "ionic_radius_ang", "covalent_radius_ang", "atomic_volume_cm3_per_mol", "electrical_conductivity_mho_per_cm", "specific_heat_J_per_g_K", "heat_of_fusion_kJ_per_mol", "heat_of_vaporization_kJ_per_mol", "thermal_conductivity_W_per_m_K", "polarizability_A_cubed", "heat_atomization_kJ_per_mol")),
    rice_proteins = c(c("ID_number", "Protein_name", "Protein_function")),
    solvents = c(c("solvent", "formula", "boiling_point", "melting_point", "density", "miscible_with_water", "solubility_in_water", "relative_polarity", "vapor_pressure", "CAS_number", "formula_weight", "refractive_index", "specific_gravity", "category")),
    tequila_chemistry = c(c("class", "compound", "bottle", "replicate", "value", "average")),
    unknown_mushroom = c(c("cap_diameter", "cap_shape", "cap_color", "does_bruise_or_bleed", "gill_color", "stem_height", "stem_width", "stem_color", "has_ring", "habitat", "season")),
    unknown_smoke = c(c("wood_type", "species", "abundance_mg_g", "compound_class", "compound")),
    us_airports = c(c("airport", "name", "city", "state", "region", "latitude", "longitude")),
    wine_grape_data = c(c("cultivar", "treatment", "chromatography", "metabolite", "log_abundance")),
    wine_quality = c(c("type", "quality_score", "quality_category", "fixed_acidity", "volatile_acidity", "citric_acid", "residual_sugar", "chlorides", "free_sulfur_dioxide", "total_sulfur_dioxide", "density", "pH", "sulphates", "alcohol")),
    wood_smoke = c(c("wood_type", "species", "abundance_mg_g", "compound_class", "compound"))
)

test_that("each shipped dataset loads with its columns and no readr spec", {
    for (nm in names(expected_columns)) {
        x <- getExportedValue("phylochemistry", nm)
        expect_s3_class(x, "tbl_df")
        expect_false(inherits(x, "spec_tbl_df"), label = nm)
        expect_null(attr(x, "spec"), label = nm)
        expect_identical(names(x), expected_columns[[nm]], label = nm)
        expect_gt(nrow(x), 0)
    }
})

test_that("datasets are reachable by bare name after library(phylochemistry)", {
    expect_true(is.data.frame(alaska_lake_data))
    expect_identical(getDataset("alaska_lake_data"), alaska_lake_data)
})

test_that("getDataset() lists every dataset and refuses unknown names", {
    listing <- getDataset()
    expect_true(all(c("name", "file", "type", "cached") %in% names(listing)))
    expect_true(all(names(expected_columns) %in% listing$name))
    expect_true(all(c("ny_trees", "lake_superior_shoreline", "mushrooms", "OSC_sequences") %in%
                    listing$name[listing$type == "download"]))
    expect_error(getDataset("not_a_dataset"), "Unknown dataset")
})

test_that("getDataset() downloads and caches a large dataset", {
    skip_on_cran()
    skip_if_offline("thebustalab.github.io")
    withr::local_envvar(R_USER_CACHE_DIR = withr::local_tempdir())
    x <- getDataset("lake_superior_shoreline")
    expect_s3_class(x, "tbl_df")
    expect_gt(nrow(x), 0)
    expect_true(getDataset()$cached[getDataset()$name == "lake_superior_shoreline"])
})
