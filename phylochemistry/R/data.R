## Datasets shipped with the package (lazy-loaded). Generated from sample_data/ on the website;
## the CSVs stay at their public URLs, which the @source lines point at.

#' algae_data
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`algae_data.csv`).
#'
#' @format A tibble with 180 rows and 5 columns:
#' \describe{
#'   \item{replicate}{numeric}
#'   \item{algae_strain}{character}
#'   \item{harvesting_regime}{character}
#'   \item{chemical_species}{character}
#'   \item{abundance}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/algae_data.csv>
"algae_data"

#' alaska_lake_data
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`alaska_lake_data.csv`).
#'
#' @format A tibble with 220 rows and 7 columns:
#' \describe{
#'   \item{lake}{character}
#'   \item{park}{character}
#'   \item{water_temp}{numeric}
#'   \item{pH}{numeric}
#'   \item{element}{character}
#'   \item{mg_per_L}{numeric}
#'   \item{element_type}{character}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/alaska_lake_data.csv>
"alaska_lake_data"

#' solvents
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`solvents.csv`).
#'
#' @format A tibble with 32 rows and 14 columns:
#' \describe{
#'   \item{solvent}{character}
#'   \item{formula}{character}
#'   \item{boiling_point}{numeric}
#'   \item{melting_point}{numeric}
#'   \item{density}{numeric}
#'   \item{miscible_with_water}{logical}
#'   \item{solubility_in_water}{numeric}
#'   \item{relative_polarity}{numeric}
#'   \item{vapor_pressure}{numeric}
#'   \item{CAS_number}{character}
#'   \item{formula_weight}{numeric}
#'   \item{refractive_index}{numeric}
#'   \item{specific_gravity}{numeric}
#'   \item{category}{character}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/solvents.csv>
"solvents"

#' periodic_table
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`per_table.csv`).
#'
#' @format A tibble with 118 rows and 41 columns:
#' \describe{
#'   \item{atomic_number}{numeric}
#'   \item{element_name}{character}
#'   \item{atomic_symbol}{character}
#'   \item{group_number}{character}
#'   \item{period}{numeric}
#'   \item{atomic_mass_rounded}{numeric}
#'   \item{melting_point_C}{numeric}
#'   \item{boiling_point_C}{numeric}
#'   \item{state_at_RT}{character}
#'   \item{density_g_per_mL}{numeric}
#'   \item{electronegativity_pauling}{numeric}
#'   \item{first_ionization_poten_eV}{numeric}
#'   \item{second_ionization_poten_eV}{numeric}
#'   \item{third_ionization_poten_eV}{numeric}
#'   \item{electron_affinity_eV}{numeric}
#'   \item{atomic_radius_ang}{numeric}
#'   \item{ionic_radius_ang}{numeric}
#'   \item{covalent_radius_ang}{numeric}
#'   \item{atomic_volume_cm3_per_mol}{numeric}
#'   \item{crystal_structure}{character}
#'   \item{electrical_conductivity_mho_per_cm}{numeric}
#'   \item{specific_heat_J_per_g_K}{numeric}
#'   \item{heat_of_fusion_kJ_per_mol}{numeric}
#'   \item{heat_of_vaporization_kJ_per_mol}{numeric}
#'   \item{thermal_conductivity_W_per_m_K}{numeric}
#'   \item{mg_per_kg_in_Earths_crust}{numeric}
#'   \item{mg_per_L_in_seawater}{numeric}
#'   \item{percent_human_body_mass}{numeric}
#'   \item{Notes}{character}
#'   \item{color}{character}
#'   \item{structure}{character}
#'   \item{hardness_mohs}{numeric}
#'   \item{hydride}{character}
#'   \item{oxide}{character}
#'   \item{chloride}{character}
#'   \item{polarizability_A_cubed}{numeric}
#'   \item{heat_atomization_kJ_per_mol}{numeric}
#'   \item{year_discovered}{numeric}
#'   \item{rel_abund_solar_system_log}{numeric}
#'   \item{abundance_earths_crust_log}{numeric}
#'   \item{cost_pure_dollar_per_100g}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/per_table.csv>
"periodic_table"

#' fadb_sample
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`fadb_sample.csv`).
#'
#' @format A tibble with 7 rows and 11 columns:
#' \describe{
#'   \item{fatty_acid}{character}
#'   \item{Agonandra_brasiliensis}{numeric}
#'   \item{Agonandra_silvatica}{numeric}
#'   \item{Agonandra_excelsa}{numeric}
#'   \item{Heisteria_silvianii}{numeric}
#'   \item{Malania_oleifera}{numeric}
#'   \item{Ximenia_americana}{numeric}
#'   \item{Ongokea_gore}{numeric}
#'   \item{Comandra_pallida}{numeric}
#'   \item{Buckleya_distichophylla}{numeric}
#'   \item{Nuytsia_floribunda}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/fadb_sample.csv>
"fadb_sample"

#' periodic_table_subset
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`per_table_small.csv`).
#'
#' @format A tibble with 56 rows and 25 columns:
#' \describe{
#'   \item{atomic_number}{numeric}
#'   \item{atomic_symbol}{character}
#'   \item{group_number}{numeric}
#'   \item{period}{numeric}
#'   \item{atomic_mass_rounded}{numeric}
#'   \item{melting_point_C}{numeric}
#'   \item{boiling_point_C}{numeric}
#'   \item{state_at_RT}{character}
#'   \item{density_g_per_mL}{numeric}
#'   \item{electronegativity_pauling}{numeric}
#'   \item{first_ionization_poten_eV}{numeric}
#'   \item{second_ionization_poten_eV}{numeric}
#'   \item{third_ionization_poten_eV}{numeric}
#'   \item{electron_affinity_eV}{numeric}
#'   \item{atomic_radius_ang}{numeric}
#'   \item{ionic_radius_ang}{numeric}
#'   \item{covalent_radius_ang}{numeric}
#'   \item{atomic_volume_cm3_per_mol}{numeric}
#'   \item{electrical_conductivity_mho_per_cm}{numeric}
#'   \item{specific_heat_J_per_g_K}{numeric}
#'   \item{heat_of_fusion_kJ_per_mol}{numeric}
#'   \item{heat_of_vaporization_kJ_per_mol}{numeric}
#'   \item{thermal_conductivity_W_per_m_K}{numeric}
#'   \item{polarizability_A_cubed}{numeric}
#'   \item{heat_atomization_kJ_per_mol}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/per_table_small.csv>
"periodic_table_subset"

#' metabolomics_data
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`metabolomics_data.csv`).
#'
#' @format A tibble with 93 rows and 126 columns:
#' \describe{
#'   \item{patient_number}{numeric}
#'   \item{patient_status}{character}
#'   \item{Pyruvate}{numeric}
#'   \item{MethylSuccinate}{numeric}
#'   \item{Phosphoglyceric Acid}{numeric}
#'   \item{Glucose 1-phosphate}{numeric}
#'   \item{Malonic Acid}{numeric}
#'   \item{Fumaric Acid}{numeric}
#'   \item{Alpha-Ketoglutaric Acid}{numeric}
#'   \item{Sarcosine}{numeric}
#'   \item{Cadaverine}{numeric}
#'   \item{5-Aminovaleric Acid}{numeric}
#'   \item{Pipecolate}{numeric}
#'   \item{Tyramine}{numeric}
#'   \item{1-Methylhistidine}{numeric}
#'   \item{isoValeric Acid}{numeric}
#'   \item{2-Hydroxyisovaleric Acid}{numeric}
#'   \item{3-Methyl-2-Oxovaleric Acid}{numeric}
#'   \item{2-Aminoisobutyric acid}{numeric}
#'   \item{Creatinine}{numeric}
#'   \item{Betaine}{numeric}
#'   \item{Taurine}{numeric}
#'   \item{1-Methylhistamine}{numeric}
#'   \item{Pyroglutamic Acid}{numeric}
#'   \item{Creatine}{numeric}
#'   \item{4-Hydroxyproline}{numeric}
#'   \item{Homocysteine}{numeric}
#'   \item{Cysteamine}{numeric}
#'   \item{Carnitine}{numeric}
#'   \item{Glycerate}{numeric}
#'   \item{N-AcetylGlycine}{numeric}
#'   \item{Guanidinoacetate}{numeric}
#'   \item{Benzoic acid}{numeric}
#'   \item{Citraconic Acid}{numeric}
#'   \item{Glutaric Acid}{numeric}
#'   \item{Oxaloacetate}{numeric}
#'   \item{Anthranilate}{numeric}
#'   \item{2-Hydroxyglutarate}{numeric}
#'   \item{2-Aminoadipate}{numeric}
#'   \item{Homogentisate}{numeric}
#'   \item{Hydroxyphenylpyruvic acid}{numeric}
#'   \item{Kynurenate}{numeric}
#'   \item{Xanthurenate}{numeric}
#'   \item{Pantothenic acid}{numeric}
#'   \item{L-Kynurenine}{numeric}
#'   \item{3-Hydroxykynurenine}{numeric}
#'   \item{Cytidine}{numeric}
#'   \item{Uridine}{numeric}
#'   \item{1-Methyladenosine}{numeric}
#'   \item{Hypoxanthine}{numeric}
#'   \item{Oxypurinol}{numeric}
#'   \item{Xanthine}{numeric}
#'   \item{Orotate}{numeric}
#'   \item{Allantoin}{numeric}
#'   \item{Urate}{numeric}
#'   \item{Xanthosine}{numeric}
#'   \item{AMP}{numeric}
#'   \item{IMP}{numeric}
#'   \item{ADP}{numeric}
#'   \item{Adenylosuccinate}{numeric}
#'   \item{Choline}{numeric}
#'   \item{Dimethylglycine}{numeric}
#'   \item{Niacinamide}{numeric}
#'   \item{Methylmalonate}{numeric}
#'   \item{Ascorbic Acid}{numeric}
#'   \item{4-Pyridoxic acid}{numeric}
#'   \item{Biotin}{numeric}
#'   \item{Pyridoxal 5'-phosphate}{numeric}
#'   \item{Trimethylamine}{numeric}
#'   \item{Trimethylamine N-oxide}{numeric}
#'   \item{Shikimic Acid}{numeric}
#'   \item{Hippuric Acid}{numeric}
#'   \item{4-Hydroxybutyrate}{numeric}
#'   \item{Acetylcholine}{numeric}
#'   \item{Deoxycarnitine}{numeric}
#'   \item{13-HODE}{numeric}
#'   \item{18-Hydroxyarachidonic acid}{numeric}
#'   \item{Glycochenodeoxycholate}{numeric}
#'   \item{Acetylcarnitine}{numeric}
#'   \item{Azelaic Acid}{numeric}
#'   \item{Myristic Acid}{numeric}
#'   \item{Margaric Acid}{numeric}
#'   \item{Linoleic Acid}{numeric}
#'   \item{Linolenic Acid}{numeric}
#'   \item{Mannose}{numeric}
#'   \item{Sorbitol}{numeric}
#'   \item{Glyceraldehyde}{numeric}
#'   \item{Erythrose}{numeric}
#'   \item{myoinositol}{numeric}
#'   \item{Fructose}{numeric}
#'   \item{D-Glucoronic acid}{numeric}
#'   \item{N-Acetylneuraminic acid}{numeric}
#'   \item{Lactose}{numeric}
#'   \item{Sucrose}{numeric}
#'   \item{Ornithine}{numeric}
#'   \item{Citrulline}{numeric}
#'   \item{Propionate}{numeric}
#'   \item{Oxalic Acid}{numeric}
#'   \item{Acetoacetate}{numeric}
#'   \item{Adipic Acid}{numeric}
#'   \item{Ribose-5-P}{numeric}
#'   \item{Geranyl-PP}{numeric}
#'   \item{Glycine}{numeric}
#'   \item{Alanine}{numeric}
#'   \item{Serine}{numeric}
#'   \item{Proline}{numeric}
#'   \item{Valine}{numeric}
#'   \item{Threonine}{numeric}
#'   \item{Leucine}{numeric}
#'   \item{iso-Leucine}{numeric}
#'   \item{Asparagine}{numeric}
#'   \item{Aspartic Acid}{numeric}
#'   \item{Glutamine}{numeric}
#'   \item{Lysine}{numeric}
#'   \item{Glutamic acid}{numeric}
#'   \item{Methionine}{numeric}
#'   \item{Histidine}{numeric}
#'   \item{Phenylalanine}{numeric}
#'   \item{Arginine}{numeric}
#'   \item{Tyrosine}{numeric}
#'   \item{Tryptophan}{numeric}
#'   \item{Cystine}{numeric}
#'   \item{lactate}{numeric}
#'   \item{3-Hydroxybenzoic acid}{numeric}
#'   \item{Succinate}{numeric}
#'   \item{Glucose}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/metabolomics_data.csv>
"metabolomics_data"

#' wine_grape_data
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`wine_grape_data.csv`).
#'
#' @format A tibble with 1070 rows and 5 columns:
#' \describe{
#'   \item{cultivar}{character}
#'   \item{treatment}{character}
#'   \item{chromatography}{character}
#'   \item{metabolite}{character}
#'   \item{log_abundance}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/wine_grape_data.csv>
"wine_grape_data"

#' hawaii_aquifers
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`hawaii_aquifers.csv`).
#'
#' @format A tibble with 954 rows and 6 columns:
#' \describe{
#'   \item{aquifer_code}{character}
#'   \item{well_name}{character}
#'   \item{longitude}{numeric}
#'   \item{latitude}{numeric}
#'   \item{analyte}{character}
#'   \item{abundance}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/hawaii_aquifers.csv>
"hawaii_aquifers"

#' beer_components
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`beer_components.csv`).
#'
#' @format A tibble with 3625 rows and 6 columns:
#' \describe{
#'   \item{ingredient}{character}
#'   \item{replicate}{numeric}
#'   \item{analyte}{character}
#'   \item{analyte_class}{character}
#'   \item{abundance}{numeric}
#'   \item{analyte_odor}{character}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/beer_components.csv>
"beer_components"

#' wood_smoke
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`wood_smoke_data.csv`).
#'
#' @format A tibble with 1524 rows and 5 columns:
#' \describe{
#'   \item{wood_type}{character}
#'   \item{species}{character}
#'   \item{abundance_mg_g}{numeric}
#'   \item{compound_class}{character}
#'   \item{compound}{character}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/wood_smoke_data.csv>
"wood_smoke"

#' unknown_smoke
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`unknown_smoke.csv`).
#'
#' @format A tibble with 254 rows and 5 columns:
#' \describe{
#'   \item{wood_type}{character}
#'   \item{species}{character}
#'   \item{abundance_mg_g}{numeric}
#'   \item{compound_class}{character}
#'   \item{compound}{character}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/unknown_smoke.csv>
"unknown_smoke"

#' hops_components
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`hops_components.csv`).
#'
#' @format A tibble with 38 rows and 12 columns:
#' \describe{
#'   \item{hop_variety}{character}
#'   \item{hop_origin}{character}
#'   \item{hop_brewing_usage}{character}
#'   \item{hop_aroma}{character}
#'   \item{total_oil}{numeric}
#'   \item{alpha_acids}{numeric}
#'   \item{beta_acids}{numeric}
#'   \item{humulene}{numeric}
#'   \item{myrcene}{numeric}
#'   \item{humulone}{numeric}
#'   \item{caryophyllene}{numeric}
#'   \item{farnesene}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/hops_components.csv>
"hops_components"

#' tequila_chemistry
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`tequila_chemistry.csv`).
#'
#' @format A tibble with 4720 rows and 6 columns:
#' \describe{
#'   \item{class}{character}
#'   \item{compound}{character}
#'   \item{bottle}{character}
#'   \item{replicate}{numeric}
#'   \item{value}{numeric}
#'   \item{average}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/tequila_chemistry.csv>
"tequila_chemistry"

#' chemical_blooms
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`chemical_blooms.csv`).
#'
#' @format A tibble with 78 rows and 10 columns:
#' \describe{
#'   \item{label}{character}
#'   \item{Alkanes}{numeric}
#'   \item{Sec_Alcohols}{numeric}
#'   \item{Others}{numeric}
#'   \item{Fatty_acids}{numeric}
#'   \item{Alcohols}{numeric}
#'   \item{Triterpenoids}{numeric}
#'   \item{Ketones}{numeric}
#'   \item{Other_compounds}{numeric}
#'   \item{Aldehydes}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/chemical_blooms.csv>
"chemical_blooms"

#' metabolomics_unknown
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`metabolomics_unknown.csv`).
#'
#' @format A tibble with 7 rows and 124 columns:
#' \describe{
#'   \item{Pyruvate}{numeric}
#'   \item{MethylSuccinate}{numeric}
#'   \item{Phosphoglyceric Acid}{numeric}
#'   \item{Glucose 1-phosphate}{numeric}
#'   \item{Malonic Acid}{numeric}
#'   \item{Fumaric Acid}{numeric}
#'   \item{Alpha-Ketoglutaric Acid}{numeric}
#'   \item{Sarcosine}{numeric}
#'   \item{Cadaverine}{numeric}
#'   \item{5-Aminovaleric Acid}{numeric}
#'   \item{Pipecolate}{numeric}
#'   \item{Tyramine}{numeric}
#'   \item{1-Methylhistidine}{numeric}
#'   \item{isoValeric Acid}{numeric}
#'   \item{2-Hydroxyisovaleric Acid}{numeric}
#'   \item{3-Methyl-2-Oxovaleric Acid}{numeric}
#'   \item{2-Aminoisobutyric acid}{numeric}
#'   \item{Creatinine}{numeric}
#'   \item{Betaine}{numeric}
#'   \item{Taurine}{numeric}
#'   \item{1-Methylhistamine}{numeric}
#'   \item{Pyroglutamic Acid}{numeric}
#'   \item{Creatine}{numeric}
#'   \item{4-Hydroxyproline}{numeric}
#'   \item{Homocysteine}{numeric}
#'   \item{Cysteamine}{numeric}
#'   \item{Carnitine}{numeric}
#'   \item{Glycerate}{numeric}
#'   \item{N-AcetylGlycine}{numeric}
#'   \item{Guanidinoacetate}{numeric}
#'   \item{Benzoic acid}{numeric}
#'   \item{Citraconic Acid}{numeric}
#'   \item{Glutaric Acid}{numeric}
#'   \item{Oxaloacetate}{numeric}
#'   \item{Anthranilate}{numeric}
#'   \item{2-Hydroxyglutarate}{numeric}
#'   \item{2-Aminoadipate}{numeric}
#'   \item{Homogentisate}{numeric}
#'   \item{Hydroxyphenylpyruvic acid}{numeric}
#'   \item{Kynurenate}{numeric}
#'   \item{Xanthurenate}{numeric}
#'   \item{Pantothenic acid}{numeric}
#'   \item{L-Kynurenine}{numeric}
#'   \item{3-Hydroxykynurenine}{numeric}
#'   \item{Cytidine}{numeric}
#'   \item{Uridine}{numeric}
#'   \item{1-Methyladenosine}{numeric}
#'   \item{Hypoxanthine}{numeric}
#'   \item{Oxypurinol}{numeric}
#'   \item{Xanthine}{numeric}
#'   \item{Orotate}{numeric}
#'   \item{Allantoin}{numeric}
#'   \item{Urate}{numeric}
#'   \item{Xanthosine}{numeric}
#'   \item{AMP}{numeric}
#'   \item{IMP}{numeric}
#'   \item{ADP}{numeric}
#'   \item{Adenylosuccinate}{numeric}
#'   \item{Choline}{numeric}
#'   \item{Dimethylglycine}{numeric}
#'   \item{Niacinamide}{numeric}
#'   \item{Methylmalonate}{numeric}
#'   \item{Ascorbic Acid}{numeric}
#'   \item{4-Pyridoxic acid}{numeric}
#'   \item{Biotin}{numeric}
#'   \item{Pyridoxal 5'-phosphate}{numeric}
#'   \item{Trimethylamine}{numeric}
#'   \item{Trimethylamine N-oxide}{numeric}
#'   \item{Shikimic Acid}{numeric}
#'   \item{Hippuric Acid}{numeric}
#'   \item{4-Hydroxybutyrate}{numeric}
#'   \item{Acetylcholine}{numeric}
#'   \item{Deoxycarnitine}{numeric}
#'   \item{13-HODE}{numeric}
#'   \item{18-Hydroxyarachidonic acid}{numeric}
#'   \item{Glycochenodeoxycholate}{numeric}
#'   \item{Acetylcarnitine}{numeric}
#'   \item{Azelaic Acid}{numeric}
#'   \item{Myristic Acid}{numeric}
#'   \item{Margaric Acid}{numeric}
#'   \item{Linoleic Acid}{numeric}
#'   \item{Linolenic Acid}{numeric}
#'   \item{Mannose}{numeric}
#'   \item{Sorbitol}{numeric}
#'   \item{Glyceraldehyde}{numeric}
#'   \item{Erythrose}{numeric}
#'   \item{myoinositol}{numeric}
#'   \item{Fructose}{numeric}
#'   \item{D-Glucoronic acid}{numeric}
#'   \item{N-Acetylneuraminic acid}{numeric}
#'   \item{Lactose}{numeric}
#'   \item{Sucrose}{numeric}
#'   \item{Ornithine}{numeric}
#'   \item{Citrulline}{numeric}
#'   \item{Propionate}{numeric}
#'   \item{Oxalic Acid}{numeric}
#'   \item{Acetoacetate}{numeric}
#'   \item{Adipic Acid}{numeric}
#'   \item{Ribose-5-P}{numeric}
#'   \item{Geranyl-PP}{numeric}
#'   \item{Glycine}{numeric}
#'   \item{Alanine}{numeric}
#'   \item{Serine}{numeric}
#'   \item{Proline}{numeric}
#'   \item{Valine}{numeric}
#'   \item{Threonine}{numeric}
#'   \item{Leucine}{numeric}
#'   \item{iso-Leucine}{numeric}
#'   \item{Asparagine}{numeric}
#'   \item{Aspartic Acid}{numeric}
#'   \item{Glutamine}{numeric}
#'   \item{Lysine}{numeric}
#'   \item{Glutamic acid}{numeric}
#'   \item{Methionine}{numeric}
#'   \item{Histidine}{numeric}
#'   \item{Phenylalanine}{numeric}
#'   \item{Arginine}{numeric}
#'   \item{Tyrosine}{numeric}
#'   \item{Tryptophan}{numeric}
#'   \item{Cystine}{numeric}
#'   \item{lactate}{numeric}
#'   \item{3-Hydroxybenzoic acid}{numeric}
#'   \item{Succinate}{numeric}
#'   \item{Glucose}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/metabolomics_unknown.csv>
"metabolomics_unknown"

#' metabolomics_hospital
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`metabolomics_hospital.csv`).
#'
#' @format A tibble with 200 rows and 4 columns:
#' \describe{
#'   \item{patient_number}{numeric}
#'   \item{patient_status}{character}
#'   \item{metabolite}{character}
#'   \item{concentration}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/metabolomics_hospital.csv>
"metabolomics_hospital"

#' metabolomics_hospital_unknown
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`metabolomics_hospital_unknown.csv`).
#'
#' @format A tibble with 10 rows and 3 columns:
#' \describe{
#'   \item{patient_number}{character}
#'   \item{metabolite}{character}
#'   \item{concentration}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/metabolomics_hospital_unknown.csv>
"metabolomics_hospital_unknown"

#' wine_quality
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`wine_quality.csv`).
#'
#' @format A tibble with 6497 rows and 14 columns:
#' \describe{
#'   \item{type}{character}
#'   \item{quality_score}{numeric}
#'   \item{quality_category}{character}
#'   \item{fixed_acidity}{numeric}
#'   \item{volatile_acidity}{numeric}
#'   \item{citric_acid}{numeric}
#'   \item{residual_sugar}{numeric}
#'   \item{chlorides}{numeric}
#'   \item{free_sulfur_dioxide}{numeric}
#'   \item{total_sulfur_dioxide}{numeric}
#'   \item{density}{numeric}
#'   \item{pH}{numeric}
#'   \item{sulphates}{numeric}
#'   \item{alcohol}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/wine_quality.csv>
"wine_quality"

#' rice_proteins
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`rice_proteins.csv`).
#'
#' @format A tibble with 57 rows and 3 columns:
#' \describe{
#'   \item{ID_number}{character}
#'   \item{Protein_name}{character}
#'   \item{Protein_function}{character}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/rice_proteins.csv>
"rice_proteins"

#' unknown_mushroom
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`unknown_mushroom.csv`).
#'
#' @format A tibble with 1 rows and 11 columns:
#' \describe{
#'   \item{cap_diameter}{numeric}
#'   \item{cap_shape}{character}
#'   \item{cap_color}{character}
#'   \item{does_bruise_or_bleed}{logical}
#'   \item{gill_color}{character}
#'   \item{stem_height}{numeric}
#'   \item{stem_width}{numeric}
#'   \item{stem_color}{character}
#'   \item{has_ring}{logical}
#'   \item{habitat}{character}
#'   \item{season}{character}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/unknown_mushroom.csv>
"unknown_mushroom"

#' passenger_flows
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`passenger_flows.csv`).
#'
#' @format A tibble with 192 rows and 3 columns:
#' \describe{
#'   \item{origin}{character}
#'   \item{destination}{character}
#'   \item{passengers}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/passenger_flows.csv>
"passenger_flows"

#' us_airports
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`us_airports.csv`).
#'
#' @format A tibble with 29 rows and 7 columns:
#' \describe{
#'   \item{airport}{character}
#'   \item{name}{character}
#'   \item{city}{character}
#'   \item{state}{character}
#'   \item{region}{character}
#'   \item{latitude}{numeric}
#'   \item{longitude}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/us_airports.csv>
"us_airports"

#' meteorite_fragments
#'
#' Teaching dataset from the Integrated Bioanalytical Chemistry course (`meteorite_fragments.csv`).
#'
#' @format A tibble with 37 rows and 6 columns:
#' \describe{
#'   \item{fragment}{character}
#'   \item{find_site}{character}
#'   \item{iron_percent}{numeric}
#'   \item{density_g_per_cm3}{numeric}
#'   \item{nickel_percent}{numeric}
#'   \item{silica_percent}{numeric}
#' }
#' @source <https://thebustalab.github.io/phylochemistry/sample_data/meteorite_fragments.csv>
"meteorite_fragments"

