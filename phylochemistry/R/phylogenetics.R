## phylogenetics — moved verbatim from phylochemistry.R (package conversion, commit 1)



        #### readTree


            #' Reads a phylogenetic tree
            #'
            #' @param tree_in_path The path to the phylogenetic tree
            #' @export
            readTree <- function( tree_in_path ) {

                ## Read and return the tree
                    return( ape::read.tree(file = tree_in_path) )
            }


        #### writeTree


            #' Writes a phylogenetic tree
            #'
            #' @param tree_out_path The path to where the tree should be written
            #' @param tree Not yet documented.
            #' @export
            writeTree <- function( tree, tree_out_path ) {
                ape::write.tree( phy = tree, file = tree_out_path )
            }


    ##### Tree and taxa manipulation

        #### descendantsFortified


            #' Provides a descendants list for a fortified tree
            #'
            #' @param fortified_tree The fortified tree that should be analyzed
            #' @export
            descendantsFortified <- function(fortified_tree){
        
                input <- fortified_tree
                
                if(any(duplicated(filter(input, isTip == TRUE)$label))){ stop("There are duplicate tip labels") }
                output <- list()
                entry <- 1
                root_node_number <- min(input$parent)

                for (node in 1:length(dropNA(input$node))) {
                    parent <- 0
                    if(input$isTip[node] == TRUE) {label = input$label[node]} else {label = node}
                    while(node != root_node_number) {
                        parent <- input[input$node == node,]$parent
                        output[[entry]] <- data.frame(label = label, parent = parent)
                        node <- parent
                        entry <- entry + 1
                    }
                }

                output <- do.call(rbind, output)
                return(output)
               
            }


        #### buildTree


            #' Construct various types of phylogenetic trees from alignments or other trees
            #'
            #' This function allows the construction of phylogenetic trees using different types of input data.
            #'
            #' @param scaffold_type The type of information that should be used to build the tree. Options include "amin_alignment", "nucl_alignment", "newick", or "url".
            #' @param scaffold_in_path The path to the information that should be used to build the tree. If using a URL, a recommended source is:
            #'        https://sftp.kew.org/pub/treeoflife/current_release/tree/species/
            #' @param members The tips of the tree that should be included. Default is to include everything.
            #' @param ml TRUE/FALSE whether to use maximum likelihood when constructing the tree.
            #' @param model_test TRUE/FALSE whether to test various maximum likelihood models while constructing the tree.
            #' @param bootstrap TRUE/FALSE whether to calculate bootstrap values for tree nodes.
            #' @param ancestral_states TRUE/FALSE whether to calculate ancestral states at nodes in the tree. Requires specifying a root via the 'root' parameter.
            #' @param root The tree tip to use as the root of the tree.
            #' @export
            buildTree <-    function(
                                scaffold_type = c("amin_alignment", "nucl_alignment", "newick", "url"),
                                scaffold_in_path,
                                members = NULL,
                                ml = FALSE, 
                                model_test = FALSE,
                                bootstrap = FALSE,
                                ancestral_states = FALSE,
                                root = NULL
                            ) {

                if ( scaffold_type == "nucl_alignment" ) {

                    ## Create distrance tree
                        cat(paste("Making tree with ", scaffold_in_path," ...\n", sep = ""))
                            nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path), format = "fasta", type = "DNA")
                            
                        cat("Creating neighbor-joining tree...\n")
                            dm <- phangorn::dist.ml(nucl_seqs_aligned, "F81")
                            NJ_tree <- phangorn::NJ(dm)
                            output <- NJ_tree

                    ## Make ml tree
                        if ( ml == TRUE ) {
                            if (any(duplicated(nucl_seqs_aligned))) { stop("There are identical sequences in your alignment, maximum liklihood will not work properly.") }
                            ## Test all available nucl models, use the best one to optimize for ml
                                if ( model_test == TRUE ) {
                                    cat(paste("Testing 24 maximum liklihood models... \n"))
                                    mt <- phangorn::modelTest(nucl_seqs_aligned, tree = NJ_tree, multicore = TRUE)
                                    best_nucl_model <- gsub("\\+.*$", "", mt$Model[which.max(mt$logLik)])
                                    cat(paste("Tested 24 models, using best model:", as.character(gsub("\\+.*$","",best_nucl_model)), "\n", sep = " "))
                                } else {
                                    best_nucl_model <- "GTR"
                                }
                                cat("Creating maximum liklihood tree...\n")
                                ML_tree_start <- phangorn::pml(NJ_tree, nucl_seqs_aligned, k = 4)
                                ML_tree_optimized <- phangorn::optim.pml(ML_tree_start, rearrangement = "stochastic", optInv = TRUE, optGamma = TRUE, model = as.character(best_nucl_model))
                                output <- ML_tree_optimized$tree
                        }

                    ## Run bootstrap analysis
                        if ( bootstrap == TRUE ) {
                            if ( ml == FALSE ) {
                                stop("To calculate bootstrap values, please also run maximum likelihood estimation (ml = TRUE).\n")
                            }
                            bootstraps <- phangorn::bootstrap.pml(ML_tree_optimized, bs = 100, optNni = TRUE, multicore = FALSE)
                            ML_tree_optimized$tree$node.label <- phangorn::plotBS(ML_tree_optimized$tree, bootstraps)$node.label
                            output <- ML_tree_optimized$tree
                        }

                    ## Root the tree and run ancestral states
                        if ( ancestral_states == TRUE ) {
                            if ( ml == FALSE ) {
                                cat("To enable ancestral state reconstruction, please also run maximum likelihood estimation (ml = TRUE).\n")
                            }
                            if ( ml == TRUE ) {
                                # Root the tree, then calculate ancestral_states
                                    ML_tree_optimized$tree <- ape::root(ML_tree_optimized$tree, as.character(root))
                                    output <- list()
                                    output$tree <- ML_tree_optimized$tree
                                    output$ancestral_states <- phangorn::ancestral.pml(ML_tree_optimized)
                            }
                        }

                    ## Root the tree if no ancestral_states were requested
                        if ( ancestral_states == FALSE ) {
                            if ( length(root) > 0 ) {
                                output <- ape::root(output, as.character(root))
                            }
                        }

                    ## Return the tree
                        cat("Pro tip: most tree read/write functions reset node numbers.\nFortify your tree and save it as a csv file to preserve node numbering.\nDo not save your tree as a newick or nexus file.\n")
                        return( output )
                }

                if ( scaffold_type == "amin_alignment" ) {

                    ## Make distance tree
                        cat(paste("Making tree with ", scaffold_in_path," ...\n", sep = ""))
                            amin_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path), format = "fasta", type = "AA")
                                
                        # cat("Creating neighbor-joining tree...\n")
                        dm = phangorn::dist.ml(amin_seqs_aligned, model = "JTT")
                        NJ_tree = phangorn::NJ(dm)
                        output <- NJ_tree

                    ## Make ml tree
                        if ( ml == TRUE ) {
                            if (any(duplicated(amin_seqs_aligned))) { stop("There are identical sequences in your alignment, maximum liklihood will not work properly.") }
                            if ( model_test == TRUE ) {
                                ## Test all available amino acid models and extract the best one
                                    cat(paste("Testing 24 maximum liklihood models... \n"))
                                    mt <- phangorn::modelTest(amin_seqs_aligned, tree = NJ_tree, model = "all", multicore = TRUE)
                                    best_amin_model <- gsub("\\+.*$", "", mt$Model[which.max(mt$logLik)])
                                    cat(paste("Tested 24 models, using best model:", as.character(gsub("\\+.*$", "", best_amin_model)), "\n", sep = " "))
                            } else {
                                best_amin_model <- "JTT"
                            }
                            cat("Creating maximum liklihood tree...\n")
                            ML_tree_start <- phangorn::pml(NJ_tree, amin_seqs_aligned, model = as.character(best_amin_model), k = 4, inv = .2)
                            ML_tree_optimized <- phangorn::optim.pml(ML_tree_start, rearrangement = "stochastic", optInv = TRUE, optGamma = TRUE)
                            output <- ML_tree_optimized$tree
                        }

                    ## Run bootstrap analysis
                        if ( bootstrap == TRUE ) {
                            if ( ml == FALSE ) {
                                stop("To calculate bootstrap values, please also run maximum likelihood estimation (ml = TRUE)")
                            }
                            bootstraps <- phangorn::bootstrap.pml(ML_tree_optimized, bs = 100, optNni = TRUE, multicore = FALSE)
                            ML_tree_optimized$tree$node.label <- phangorn::plotBS(ML_tree_optimized$tree, bootstraps)$node.label
                        }

                    ## Root the tree and run ancestral states
                        if ( ancestral_states == TRUE ) {
                            if ( ml == FALSE ) {
                                cat("To enable ancestral state reconstruction, please also run maximum likelihood estimation (ml = TRUE).\n")
                            }
                            if ( ml == TRUE ) {
                                if ( length(root) > 0 ) {
                                    ML_tree_optimized$tree <- ape::root(ML_tree_optimized$tree, as.character(root))
                                }
                                output <- list()
                                output$tree <- ML_tree_optimized$tree
                                ancestral_states <- phangorn::ancestral.pml(ML_tree_optimized, return = "prob")

                                ancestral_states_output <- data.frame()
                                for( i in 1:length(ancestral_states) ) {
                                    temp <- ancestral_states[[i]]
                                    colnames(temp) <- attr(ancestral_states, "levels")
                                    rownames(temp) <- NULL
                                    ancestral_states_output <- rbind(
                                        ancestral_states_output,
                                        cbind(
                                            data.frame(
                                                name = names(ancestral_states)[i],
                                                position = seq(1,dim(temp)[1],1)
                                            ),
                                            temp
                                        )
                                    )
                                }
                                
                                output$ancestral_states <- ancestral_states_output

                            }
                        }

                    ## Return tree
                        cat("Pro tip: most tree read/write functions reset node numbers.\nFortify your tree and save it as a csv file to preserve node numbering.\nDo not save your tree as a newick or nexus file.\n")
                        return ( output )
                }

                if ( scaffold_type == "newick" ) {

                    ## Read in the newick scaffold
                        if (length(grep("google", scaffold_in_path)) > 0) {
                            temp_tree <- tempfile(fileext = ".newick")
                            suppressMessages(googledrive::drive_download(
                                file = googledrive::as_id(scaffold_in_path),
                                path = temp_tree,
                                overwrite = TRUE
                            ))
                            newick <- readTree(temp_tree)
                        } else {
                            newick <- readTree( tree_in_path = scaffold_in_path )    
                        }

                    ## Are the Genus_species in your members in the newick? Are the genera in your members in the newick?
                        compatibility <- data.frame( Genus_species = unique(members), is_species_in_tree = NA, is_genus_in_tree = NA )
                        compatibility$is_species_in_tree <- compatibility$Genus_species %in% newick$tip.label
                        compatibility$is_genus_in_tree <- gsub("_.*$", "", compatibility$Genus_species) %in% gsub("_.*$", "", as.character(newick$tip.label))

                    substitutions <- list()
                    removals <- list()

                    if ( all(compatibility$is_species_in_tree) == FALSE ) {
                        ## For Genus_species in members whose genus is missing from the tree (orphans), remove them
                            orphans <- as.character(dplyr::filter(compatibility, is_species_in_tree == FALSE & is_genus_in_tree == FALSE)$Genus_species)
                            members <- members[!(members %in% orphans)]
                            if ( length(orphans) > 0 ) {
                                removals <- c(removals, orphans)
                            }

                        ## Check compatibility again
                            compatibility <- data.frame(Genus_species = unique(members), is_species_in_tree = NA, is_genus_in_tree = NA)
                            compatibility$is_species_in_tree <- compatibility$Genus_species %in% newick$tip.label
                            compatibility$is_genus_in_tree <- gsub("_.*$", "", compatibility$Genus_species) %in% gsub("_.*$", "", as.character(newick$tip.label))

                        ## For unique(members$Genus_species) in members not in the tree but whose genus in the tree (adoptees), substitute
                            adoptees <- as.character(dplyr::filter(compatibility, is_species_in_tree == FALSE & is_genus_in_tree == TRUE)$Genus_species)

                            for ( i in 1:length(adoptees) ) {
                                potential_foster_species <- newick$tip.label[gsub("_.*$", "", newick$tip.label) %in% gsub("_.*$", "", adoptees[i])] # all species in tree of the adoptees genus
                                available_foster_species <- potential_foster_species[!potential_foster_species %in% unique(members)] # potential_foster_species not already in the quantities
                                if ( length(available_foster_species) == 0) {
                                    members <- members[!members %in% adoptees[i]]
                                    removals <- c(removals, adoptees[i])
                                } else {
                                    substitutions <- c(substitutions, paste(available_foster_species[1], "substituted with", adoptees[i]))
                                    newick$tip.label[newick$tip.label == as.character(available_foster_species[1])] <- as.character(adoptees[i])
                                }
                            }
                    }

                    if (length(substitutions) > 0 || length(removals) > 0) {
                        message("IMPORTANT: Some species substitutions or removals were made as part of buildTree. Run build_tree_substitutions() to see them all.")
                    }

                    build_tree_substitutions <- function() {
                        if (length(removals) > 0) {
                            message("The following species were removed: ", paste(removals, collapse = ", "))
                        }
                        if (length(substitutions) > 0) {
                            message("The following substitutions were made: ", paste(substitutions, collapse = "; "))
                        }
                    }
                    assign("build_tree_substitutions", build_tree_substitutions, envir = globalenv())

                    ## Drop tree tips not in desired members
                        newick <- ape::drop.tip(newick, newick$tip.label[!newick$tip.label %in% unique(members)])

                    ## Sort members according to the tree
                        ordered_tip_labels <- subset(ggtree::fortify(newick), isTip)$label[order(subset(ggtree::fortify(newick), isTip)$y, decreasing = TRUE)]
                        members <- factor(members, levels = rev(ordered_tip_labels))

                    ## Return tree
                        message("Pro tip: most tree read/write functions reset node numbers. Fortify your tree and save it as a csv file to preserve node numbering.")
                        return ( newick )
                }

                if ( scaffold_type == "url" ) {

                    ## Read in the tree from the URL
                        newick <- ape::read.tree(scaffold_in_path)

                    ## Are the Genus_species in your members in the newick? Are the genera in your members in the newick?
                        compatibility <- data.frame( Genus_species = unique(members), is_species_in_tree = NA, is_genus_in_tree = NA )
                        compatibility$is_species_in_tree <- compatibility$Genus_species %in% newick$tip.label
                        compatibility$is_genus_in_tree <- gsub("_.*$", "", compatibility$Genus_species) %in% gsub("_.*$", "", as.character(newick$tip.label))

                    if ( all(compatibility$is_species_in_tree) == FALSE ) {
                        ## For Genus_species in members whose genus is missing from the tree (orphans), remove them
                            orphans <- as.character(dplyr::filter(compatibility, is_species_in_tree == FALSE & is_genus_in_tree == FALSE)$Genus_species)
                            members <- members[!(members %in% orphans)]
                            if ( length(orphans) > 0 ) {
                                message("The following species belong to a genus not found in the newick scaffold and were removed: ")
                                for ( orphan in 1:length(orphans) ) {
                                    # message("\n")
                                    message(orphans[orphan])
                                }
                            }

                        ## Check compatibility again
                            compatibility <- data.frame(Genus_species = unique(members), is_species_in_tree = NA, is_genus_in_tree = NA)
                            compatibility$is_species_in_tree <- compatibility$Genus_species %in% newick$tip.label
                            compatibility$is_genus_in_tree <- gsub("_.*$", "", compatibility$Genus_species) %in% gsub("_.*$", "", as.character(newick$tip.label))

                        ## For unique(members$Genus_species) in members not in the tree but whose genus in the tree (adoptees), substitute
                            adoptees <- as.character(dplyr::filter(compatibility, is_species_in_tree == FALSE & is_genus_in_tree == TRUE)$Genus_species)

                            for ( i in 1:length(adoptees) ) {
                                potential_foster_species <- newick$tip.label[gsub("_.*$", "", newick$tip.label) %in% gsub("_.*$", "", adoptees[i])] # all species in tree of the adoptees genus
                                available_foster_species <- potential_foster_species[!potential_foster_species %in% unique(members)] # potential_foster_species not already in the quantities
                                if ( length(available_foster_species) == 0) {
                                    members <- members[!members %in% adoptees[i]]
                                    message(paste("There aren't enough fosters to include the following species in the tree so it was removed:", adoptees[i], "\n", sep = " "))
                                } else {
                                    message(paste("Scaffold newick tip", available_foster_species[1], "substituted with", adoptees[i], sep = " "))
                                    newick$tip.label[newick$tip.label == as.character(available_foster_species[1])] <- as.character(adoptees[i])
                                }
                            }
                    }

                    ## Drop tree tips not in desired members
                        newick <- ape::drop.tip(newick, newick$tip.label[!newick$tip.label %in% unique(members)])

                    ## Sort members according to the tree
                        ordered_tip_labels <- subset(ggtree::fortify(newick), isTip)$label[order(subset(ggtree::fortify(newick), isTip)$y, decreasing = TRUE)]
                        members <- factor(members, levels = rev(ordered_tip_labels))

                    ## Return tree
                        message("Pro tip: most tree read/write functions reset node numbers.\nFortify your tree and save it as a csv file to preserve node numbering.\nDo not save your tree as a newick or nexus file.\n")
                        return ( newick )
                }

            }


        #### pruneTree


            #' Prune a tree so it only contains user-specified tips
            #'
            #' @param tree_in Either the path to the tree to manipulate, an in-memory phylo object, or a fortified tree
            #' @param tips_to_keep The tips to keep
            #' @export
            pruneTree <- function(tree_in, tips_to_keep) {

                ## Check if tree_in is a path, a phylo object, or a fortified tree
                if (is.character(tree_in)) {
                    tree <- readTree(tree_in)
                } else if (inherits(tree_in, "phylo")) {
                    tree <- tree_in
                } else if (is.data.frame(tree_in) && "node" %in% colnames(tree_in)) {
                    fortified_tree <- tree_in
                    tree <- ape::as.phylo(fortified_tree)
                } else {
                    stop("tree_in must be either a file path, a phylo object, or a fortified tree")
                }

                ## Drop all tips but those specified by user
                tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% tips_to_keep])

                ## If the input was a fortified tree, convert back to fortified and join with original
                if (exists("fortified_tree", inherits = FALSE)) {
                    pruned_fortified_tree <- ggtree::fortify(tree)
                    pruned_fortified_tree <- dplyr::right_join(fortified_tree, pruned_fortified_tree, by = "node")
                    return(pruned_fortified_tree)
                }

                return(tree)
            }


        #### collapseTree


            #' Create a high-level cladogram from a low-level tree
            #'
            #' This function takes a detailed phylogenetic tree and collapses it into a higher-level cladogram based on specified associations.
            #' It is useful for summarizing complex trees by grouping tips into broader categories.
            #'
            #' @param tree A phylogenetic tree to manipulate. It can be an object of class 'phylo' or a fortified tree.
            #' @param associations A data frame or matrix with two columns. The first column contains tree tip names, and the second column contains the
            #'        corresponding higher-level group names to which these tips should be summarized.
            #' @return A phylogenetic tree of class 'phylo' with tips renamed and collapsed according to the provided associations.
            #' @examples
            #' \dontrun{
            #' # Example usage:
            #' # tree <- readTree("path/to/tree")
            #' # associations <- data.frame(tip_name = c("tip1", "tip2"), group_name = c("group1", "group1"))
            #' # collapsed_tree <- collapseTree(tree, associations)
            #' # plot(collapsed_tree)
            #'
            #' collapseTree()
            #' }
            #' @export
            collapseTree <- function(tree, associations) {

                ## Convert associations into a data frame
                associations <- as.data.frame(associations)

                ## Check if the tree is fortified
                if (is.data.frame(tree) && "node" %in% colnames(tree)) {
                    message("Warning: Data columns are only being preserved for select members of the original data frame.")
                    fortified_tree <- tree
                    tree <- ape::as.phylo(fortified_tree)
                }

                ## Drop all species names not in the association list
                tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% associations[,1]])

                ## Subset the association list to one member per family
                tip_labels <- associations[associations[,1] %in% tree$tip.label,]

                ## Collapse tree, rename tips
                if (length(duplicated(tip_labels[,2])) > 0) {
                    collapsed_tree <- ape::drop.tip(tree, tip_labels[,1][duplicated(tip_labels[,2])])
                } else {
                    collapsed_tree <- tree
                }
                
                collapsed_tree$tip.label <- associations[,2][match(collapsed_tree$tip.label, associations[,1])]

                ## Set branch lengths to one
                collapsed_tree$edge.length <- rep(1, length(collapsed_tree$edge.length))
                message("Branch lengths have been set to one.")

                ## If the input was a fortified tree, convert back to fortified and join with original
                if (exists("fortified_tree", inherits = FALSE)) {
                    pruned_fortified_tree <- ggtree::fortify(collapsed_tree)
                    pruned_fortified_tree <- dplyr::right_join(fortified_tree, pruned_fortified_tree, by = "node")
                    return(pruned_fortified_tree)
                }

                return(collapsed_tree)
            }


    ##### Phylogenetic statistical testing

        #### runPhylogeneticAnalyses

            #' Run phylogenetic comparative analyses on traits
            #'
            #' @param traits Not yet documented.
            #' @param column_w_names_of_tiplabels Not yet documented.
            #' @param column_w_names_of_traits Not yet documented.
            #' @param column_w_values_for_traits Not yet documented.
            #' @param tree Not yet documented.
            #' @param replicates Not yet documented.
            #' @param cost Not yet documented.
            #' @export
            runPhylogeneticAnalyses <- function(
                traits,
                column_w_names_of_tiplabels,
                column_w_names_of_traits,
                column_w_values_for_traits,
                tree, ## tree should be in phylo format
                replicates = 999,
                cost = NULL
            ) {

                ## Preprocessing and TRAITS WIDE - Make tree and traits compatible

                    ## Make data frame and first column characters
                        traits <- as.data.frame(traits)
                        traits[,colnames(traits) == column_w_names_of_tiplabels] <- as.character(traits[,colnames(traits) == column_w_names_of_tiplabels])
                        data_tip_labels <- unique(traits[,colnames(traits) == column_w_names_of_tiplabels])

                    ## Drop tips and rows as necessary
                        species_dropped_from_traits <- data_tip_labels[!data_tip_labels %in% tree$tip.label]
                        if (length(species_dropped_from_traits) > 0) {
                          traits <- traits[traits[,colnames(traits) == column_w_names_of_tiplabels] %in% tree$tip.label,]
                          message(paste("Species dropped from traits by runPhylogeneticAnalyses:", paste(species_dropped_from_traits, collapse = ", ")))
                        }

                        species_dropped_from_tree <- tree$tip.label[!tree$tip.label %in% data_tip_labels]
                        if (length(species_dropped_from_tree) > 0) {
                          tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% data_tip_labels])
                          message(paste("Species dropped from tree by runPhylogeneticAnalyses:", paste(species_dropped_from_tree, collapse = ", ")))
                        }

                    ## Reorder traits to match the order of the input
                        traits[,colnames(traits) == column_w_names_of_tiplabels] <- factor(traits[,colnames(traits) == column_w_names_of_tiplabels], levels = tree$tip.label)

                    ## Fortify tree and add traits to it
                        treefort <- fortify(tree)
                        traits <- traits[,colnames(traits) %in% c(column_w_names_of_tiplabels, column_w_names_of_traits, column_w_values_for_traits)]
                        treefort_1 <- right_join(
                            treefort,
                            traits,
                            by = c("label" = colnames(traits)[colnames(traits) == column_w_names_of_tiplabels])
                        )
                        colnames(treefort_1)[colnames(treefort_1) == column_w_names_of_traits] <- "trait"
                        colnames(treefort_1)[colnames(treefort_1) == column_w_values_for_traits] <- "value"

                    ## Pivot traits wider and deal with any missing values
                        traits %>%
                            select(
                                column_w_names_of_tiplabels,
                                column_w_names_of_traits,
                                column_w_values_for_traits
                            ) -> traits_wide
                        colnames(traits_wide) <- c("sample_unique_ID", "trait", "value")
                        traits_wide <- pivot_wider(traits_wide, names_from = "trait", values_from = "value")

                        traits_wide[is.na(traits_wide)] <- 0

                ## ANCESTRAL STATES - Loop over traits and get ancestral states
                    anc_traits <- list()
                    for (i in 2:dim(traits_wide)[2]) { # i=2
                        trait <- as.numeric(unlist(traits_wide[,i]))
                        names(trait) <- as.character(unlist(traits_wide[,1]))

                        fastAnc_output <- fastAnc( tree = tree, x = trait, vars = FALSE, CI = TRUE )
                        anc_traits[[i]] <- data.frame(
                            node = as.numeric(names(fastAnc_output[[1]])),
                            trait = colnames(traits_wide)[i],
                            value = fastAnc_output[[1]]
                        )
                    }
                    anc_traits <- do.call(rbind, anc_traits)
                    treefort_2 <- right_join(treefort, anc_traits, by = "node")

                    ## Put it all together (including analyte metadata!) and return
                        treefort_3 <- rbind(
                            treefort_1[,colnames(treefort_1) %in% colnames(treefort_2)],
                            treefort_2
                        )

                ## PHYLOGENETIC SIGNAL

                    results <- list()
                    for (i in 2:dim(traits_wide)[2]) { #i=2

                        trait <- as.data.frame(traits_wide[,i])[,1]
                        names(trait) <- as.data.frame(traits_wide[,1])[,1]

                        if (class(trait) %in% c("factor", "character")) {

                            ### For discrete traits

                                ## Get the states in which the trait may exist (levels)
                                ## Check that the variable is indeed categorical by comparing the length of the levels vector to the length of the trait vector
                                    levels <- attributes(factor(trait))$levels
                                    if (length(levels) == length(trait)) { warn("Are you sure this variable is categorical?") }

                                ## Make the transition matrix
                                    if (is.null(cost)) {
                                        # By default, all transitions have a cost of 1, except transitions from a state to itself which have a cost of 0
                                        cost1 <- 1-diag(length(levels))
                                    } else {
                                        # If a custom cost matrix is provided, it should have dimensions that match the number of levels
                                    if (length(levels) != dim(cost)[1])
                                        stop("Dimensions of the character state transition matrix do not agree with the number of levels")
                                        cost1 <- t(cost)
                                    }
                                    dimnames(cost1) <- list(levels, levels)

                                ## Make the trait numeric
                                    trait_as_numeric <- as.numeric(as.factor(trait))
                                    names(trait_as_numeric) <- names(trait)

                                ## Make the phyDat object and get the parsimony score for the tree with the associated observations
                                    # obs <- t(data.frame(trait))
                                    obs <- phyDat( trait, type = "USER", levels = attributes(factor(trait))$levels )
                                    # This parsimony score represents the minimum total cost of all state changes that must have occurred on the tree given the observed trait data
                                    OBS <- parsimony( tree, obs, method = "sankoff", cost = cost1 )

                                ## Make "replicates" number of random tree-trait associations and check their parsimony score
                                    null_model <- matrix(NA, replicates, 1)
                                    for (i in 1:replicates){
                                        ## Randomize the traits and get the parsimony score for the random traits on that tree
                                            null <- sample(as.numeric(trait_as_numeric))
                                            attributes(null)$names <- attributes(trait_as_numeric)$names
                                            # null <- t(data.frame(null))
                                            null <- phyDat( null,type = "USER",levels = attributes(factor(null))$levels )
                                            null_model[i,] <- parsimony( tree, null, method = "sankoff", cost = cost1 )
                                    }

                                ## Assess observed parsimony score in the context of the random ones
                                    p_value <- sum(OBS >= null_model)/(replicates + 1)

                                ## Summarize output and report it

                                    results[[i-1]] <- data.frame(
                                        trait = colnames(traits_wide)[i],
                                        trait_type = "categorical",
                                        # n_species = length(tree$tip.label),
                                        # number_of_levels = length(attributes(factor(trait))$levels), 
                                        # evolutionary_transitions_observed = OBS,
                                        # median_evolutionary_transitions_in_randomization = median(null_model),
                                        # minimum_evolutionary_transitions_in_randomization = min(null_model),
                                        # evolutionary_transitions_in_randomization = max(null_model),
                                        phylogenetic_signal_k_value = NA,
                                        phylogenetic_signal_k_p_value = p_value,
                                        phylogenetic_signal_lambda_value = NA,
                                        phylogenetic_signal_lambda_p_value = NA
                                    )
                        }

                        if (class(trait) %in% c("numeric", "integer")) {

                            # Determine phylogenetic signal for continuous traits
                                results[[i-1]] <- data.frame(
                                    trait = colnames(traits_wide)[i],
                                    trait_type = "continuous",
                                    # n_species = length(tree$tip.label),
                                    # number_of_levels = NA,
                                    # evolutionary_transitions_observed = NA,
                                    # median_evolutionary_transitions_in_randomization = NA,
                                    # minimum_evolutionary_transitions_in_randomization = NA,
                                    # evolutionary_transitions_in_randomization = NA,
                                    phylogenetic_signal_k_value = if (sum(trait) == 0) { NA } else {
                                        as.numeric(picante::phylosignal(trait, tree)[1])
                                    },
                                    phylogenetic_signal_k_p_value = if (sum(trait) == 0) { NA } else {
                                        as.numeric(picante::phylosignal(trait, tree)[4])
                                    },
                                    phylogenetic_signal_lambda_value = if (sum(trait) == 0) { NA } else {
                                        round(phytools::phylosig(tree, trait, method = "lambda", test = TRUE)$lambda, 4)
                                    },
                                    phylogenetic_signal_lambda_p_value = if (sum(trait) == 0) { NA } else {
                                        round(phytools::phylosig(tree, trait, method = "lambda", test = TRUE)$P, 4)
                                    }
                                )
                        }
                    }

                    ## Return results
                        results <- do.call(rbind, results)
                        treefort_3 <- left_join(treefort_3, results)

                ## INDEPENDENT CONTRASTS
                    results <- list()
                    for (i in 2:(dim(traits_wide)[2])) { # i=2
                        results[[i-1]] <- data.frame(
                            pic = unname(pic(
                                x = as.numeric(as.data.frame(traits_wide)[,i]),
                                phy = tree
                            )),
                            node = filter(fortify(tree), isTip == FALSE)$node,
                            trait = colnames(traits_wide)[i]
                        )
                    }
                    results <- do.call(rbind, results)
                    return(left_join(treefort_3, results, by = c("node", "trait")))
            }


        #### harmonizeTreeTraits


            #' Harmonizes species in a set of traits (wide or long) and a tree
            #'
            #' @param traits A data frame in which the first column contains the names of the species
            #' @param tree A phylo object containing the tree
            #' @param column_w_names_of_tiplabels Not yet documented.
            #' @export
            harmonizeTreeTraits <- function( traits, column_w_names_of_tiplabels, tree ) {

                ## Make data frame and first column characters
                    traits <- as.data.frame(traits)
                    traits[,colnames(traits) == column_w_names_of_tiplabels] <- as.character(traits[,colnames(traits) == column_w_names_of_tiplabels])
                    data_tip_labels <- unique(traits[,colnames(traits) == column_w_names_of_tiplabels])

                ## Drop tips and rows as necessary
                    species_dropped_from_traits <- data_tip_labels[!data_tip_labels %in% tree$tip.label]
                    if (length(species_dropped_from_traits) > 0) {
                      traits <- traits[traits[,colnames(traits) == column_w_names_of_tiplabels] %in% tree$tip.label,]
                      message(paste("Species dropped from traits:", paste(species_dropped_from_traits, collapse = ", ")))
                    }

                    species_dropped_from_tree <- tree$tip.label[!tree$tip.label %in% data_tip_labels]
                    if (length(species_dropped_from_tree) > 0) {
                      tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% data_tip_labels])
                      message(paste("Species dropped from tree:", paste(species_dropped_from_tree, collapse = ", ")))
                    }

                ## Reorder traits to match the order of the input
                    traits[,colnames(traits) == column_w_names_of_tiplabels] <- factor(traits[,colnames(traits) == column_w_names_of_tiplabels], levels = tree$tip.label)

                return(list(tree = tree, traits = traits))
            }


        #### *phylogeneticSignal


                #' Compute phylogenetic signal for continuous and discrete traits
                #'
                #' @param tree A phylogenetic tree (class phylo) with tips that exactly match the names of trait
                #' @param replicates Number of random replications to run
                #' @param cost Optional. A specialized transition matrix
                #' @param traits Not yet documented.
                #' @param column_w_names_of_tiplabels Not yet documented.
                #' @param column_w_names_of_traits Not yet documented.
                #' @param column_w_values_for_traits Not yet documented.
                #' @export
                phylogeneticSignal <- function(
                    traits,
                    column_w_names_of_tiplabels,
                    column_w_names_of_traits,
                    column_w_values_for_traits,
                    tree,
                    replicates = 999,
                    cost = NULL
                ) {

                    ## Make tree and traits compatible, make the traits square
                        match <- harmonizeTreeTraits( 
                            traits = traits,
                            tree = tree,
                            column_w_names_of_tiplabels = column_w_names_of_tiplabels
                        )
                        traits <- match$traits
                        tree <- match$tree

                        results <- list()
                        trait_list <- data.frame(unique(traits[,colnames(traits) == column_w_names_of_traits]))[,1]

                        traits %>%
                            pivot_wider(names_from = column_w_names_of_traits, values_from = column_w_values_for_traits, values_fill = 0) -> traits
                        ## Loop over traits

                            for (i in 1:length(trait_list)) { #i=1

                                trait <- as.data.frame(traits[,colnames(traits) == trait_list[i]])[,1]
                                names(trait) <- as.data.frame(traits[,colnames(traits) == column_w_names_of_tiplabels])[,1]

                                if (class(trait) %in% c("factor", "character")) {

                                    ### For discrete traits

                                        ## Get the states in which the trait may exist (levels)

                                            levels <- attributes(factor(trait))$levels

                                        ## Check that the variable is indeed categorical by comparing the length of the levels vector to the length of the trait vector

                                            if (length(levels) == length(trait)) { warn("Are you sure this variable is categorical?") }

                                        ## Make the transition matrix

                                            if (is.null(cost)) {
                                                # By default, all transitions have a cost of 1, except transitions from a state to itself which have a cost of 0
                                                cost1 <- 1-diag(length(levels))

                                            } else {
                                                # If a custom cost matrix is provided, it should have dimensions that match the number of levels
                                            if (length(levels) != dim(cost)[1])

                                                stop("Dimensions of the character state transition matrix do not agree with the number of levels")

                                                cost1 <- t(cost)
                                            }
                                            dimnames(cost1) <- list(levels, levels)

                                        ## Make the trait numeric

                                            trait_as_numeric <- as.numeric(as.factor(trait))
                                            names(trait_as_numeric) <- names(trait)

                                        ## Make the phyDat object and get the parsimony score for the tree with the associated observations

                                            # obs <- t(data.frame(trait))
                                            obs <- phyDat( trait, type = "USER", levels = attributes(factor(trait))$levels )
                                            # This parsimony score represents the minimum total cost of all state changes that must have occurred on the tree given the observed trait data
                                            OBS <- parsimony( tree, obs, method = "sankoff", cost = cost1 )

                                        ## Make "replicates" number of random tree-trait associations and check their parsimony score

                                            null_model <- matrix(NA, replicates, 1)
                                            for (i in 1:replicates){

                                                ## Randomize the traits and get the parsimony score for the random traits on that tree
                                                    null <- sample(as.numeric(trait_as_numeric))
                                                    attributes(null)$names <- attributes(trait_as_numeric)$names
                                                    # null <- t(data.frame(null))
                                                    null <- phyDat( null,type = "USER",levels = attributes(factor(null))$levels )
                                                    null_model[i,] <- parsimony( tree, null, method = "sankoff", cost = cost1 )

                                            }

                                        ## Assess observed parsimony score in the context of the random ones

                                            p_value <- sum(OBS >= null_model)/(replicates + 1)
                                            p_value

                                        ## Summarize output and report it

                                            results[[i]] <- data.frame(
                                                trait = trait_list[i],
                                                trait_type = "categorical",
                                                n_species = length(tree$tip.label),
                                                number_of_levels = length(attributes(factor(trait))$levels), 
                                                evolutionary_transitions_observed = OBS,
                                                median_evolutionary_transitions_in_randomization = median(null_model),
                                                minimum_evolutionary_transitions_in_randomization = min(null_model),
                                                evolutionary_transitions_in_randomization = max(null_model),
                                                phylogenetic_signal_k_value = NA,
                                                phylogenetic_signal_k_p_value = p_value,
                                                phylogenetic_signal_lambda_value = NA,
                                                phylogenetic_signal_lambda_p_value = NA
                                            )
                                }

                                if (class(trait) %in% c("numeric", "integer")) {

                                    # Determine phylogenetic signal for continuous traits
                                        results[[i]] <- data.frame(
                                            trait = trait_list[i],
                                            trait_type = "continuous",
                                            n_species = length(tree$tip.label),
                                            number_of_levels = NA,
                                            evolutionary_transitions_observed = NA,
                                            median_evolutionary_transitions_in_randomization = NA,
                                            minimum_evolutionary_transitions_in_randomization = NA,
                                            evolutionary_transitions_in_randomization = NA,
                                            phylogenetic_signal_k_value = if (sum(trait) == 0) { NA } else {
                                                as.numeric(picante::phylosignal(trait, tree)[1])
                                            },
                                            phylogenetic_signal_k_p_value = if (sum(trait) == 0) { NA } else {
                                                as.numeric(picante::phylosignal(trait, tree)[4])
                                            },
                                            phylogenetic_signal_lambda_value = round(phytools::phylosig(tree, trait, method = "lambda", test = TRUE)$lambda, 4),
                                            phylogenetic_signal_lambda_p_value = round(phytools::phylosig(tree, trait, method = "lambda", test = TRUE)$P, 4)
                                        )
                                }
                            }

                    ## Return results
                        results <- do.call(rbind, results)
                        results
                        return(results)
                }


        #### independentContrasts

            #' Phylogenetically independent contrasts for traits
            #'
            #' @param traits Not yet documented.
            #' @param column_w_names_of_tiplabels Not yet documented.
            #' @param column_w_names_of_traits Not yet documented.
            #' @param column_w_values_for_traits Not yet documented.
            #' @param tree Not yet documented.
            #' @export
            independentContrasts <- function(
                traits,
                column_w_names_of_tiplabels,
                column_w_names_of_traits,
                column_w_values_for_traits,
                tree
            ) {

                ## Make tree and traits compatible
                    match <- harmonizeTreeTraits(
                        traits = traits,
                        column_w_names_of_tiplabels = column_w_names_of_tiplabels,
                        tree = tree
                    )

                ## Pivot wider and deal with missing values
                    match$traits %>%
                        select(
                            column_w_names_of_tiplabels,
                            column_w_names_of_traits,
                            column_w_values_for_traits
                        ) -> traits_wide
                    colnames(traits_wide) <- c("sample_unique_ID", "trait", "value")
                    traits_wide <- pivot_wider(traits_wide, names_from = "trait", values_from = "value")

                    traits_wide[is.na(traits_wide)] <- 0

                ## Calculate independent contrasts for every variable and return
                    results <- list()
                    for (i in 2:(dim(traits_wide)[2])) { # i=2

                        results[[i-1]] <- data.frame(
                            pic = unname(pic(
                                x = as.numeric(as.data.frame(traits_wide)[,i]),
                                phy = match$tree
                            )),
                            node = filter(fortify(match$tree), isTip == FALSE)$node,
                            trait = colnames(traits_wide)[i]
                        )

                    }
                    results <- do.call(rbind, results)
                    return(results)
            }


        #### phyloTraitDistance

            #' Compare phylogenetic and trait distances
            #'
            #' @param traits Not yet documented.
            #' @param column_w_names_of_tiplabels Not yet documented.
            #' @param column_w_values_for_traits Not yet documented.
            #' @param tree Not yet documented.
            #' @param mode Not yet documented.
            #' @export
            phyloTraitDistance <- function(
                traits,
                column_w_names_of_tiplabels,
                column_w_values_for_traits,
                tree,
                mode = c("presence_absence", "quantitative")
            ) {

                distances <- as.data.frame(cophenetic.phylo(tree))
                distances <- cbind(rownames(distances), distances)
                distances <- pivot_longer(distances, cols = 2:dim(distances)[2])
                colnames(distances) <- c("species_1", "species_2", "phylogenetic_distance")
                distances <- distances[distances$phylogenetic_distance != 0,]

                distances$trait_similarity <- 0
                for (i in 1:dim(distances)[1]) { #i=1

                    species_1_traits <- traits[which(traits[,colnames(traits) == column_w_names_of_tiplabels] == distances$species_1[i]),]
                    species_2_traits <- traits[which(traits[,colnames(traits) == column_w_names_of_tiplabels] == distances$species_2[i]),]
                    
                    ## convert to presence absence
                    if (mode[1] == "presence_absence") {
                        s1pa <- species_1_traits[,which(colnames(species_1_traits) == column_w_values_for_traits)] > 0
                        s2pa <- species_2_traits[,which(colnames(species_2_traits) == column_w_values_for_traits)] > 0
                        overlap <- s1pa == s2pa
                        distances$trait_similarity[i] <- round(sum(overlap == TRUE) / length(overlap), 4)          
                    }
                    
                    if (mode[1] == "quantitative") {
                        sum_of_differences <- sum(abs(
                            species_1_traits[,which(colnames(species_1_traits) == column_w_values_for_traits)] - 
                            species_2_traits[,which(colnames(species_2_traits) == column_w_values_for_traits)]
                        ))
                        sum_of_all_trait_values <- sum(abs(
                            species_1_traits[,which(colnames(species_1_traits) == column_w_values_for_traits)] + 
                            species_2_traits[,which(colnames(species_2_traits) == column_w_values_for_traits)]
                        ))
                        distances$trait_similarity[i] <- 1 - (sum_of_differences/sum_of_all_trait_values)
                    }
                }

                return(distances)

            }


        #### *ancestralTraits

            #' Reconstruct ancestral trait values on a tree
            #'
            #' @param traits Not yet documented.
            #' @param column_w_names_of_tiplabels Not yet documented.
            #' @param column_w_names_of_traits Not yet documented.
            #' @param column_w_values_for_traits Not yet documented.
            #' @param tree Not yet documented.
            #' @export
            ancestralTraits <- function(
                traits,
                column_w_names_of_tiplabels,
                column_w_names_of_traits,
                column_w_values_for_traits,
                tree
            ) {
                ## Make tree and traits compatible
                    match <- harmonizeTreeTraits(
                        traits = traits,
                        column_w_names_of_tiplabels = column_w_names_of_tiplabels,
                        tree = tree
                    )
                  
                ## Fortify tree and add traits to it
                    treefort <- fortify(match$tree)
                    treefort_1 <- right_join(
                        treefort,
                        match$traits,
                        by = c("label" = colnames(match$traits)[colnames(match$traits) == column_w_names_of_tiplabels])
                    )
                    colnames(treefort_1)[colnames(treefort_1) == column_w_names_of_traits] <- "trait"
                    colnames(treefort_1)[colnames(treefort_1) == column_w_values_for_traits] <- "value"
    
                ## Pivot traits wider and deal with any missing values
                    match$traits %>%
                        select(
                            column_w_names_of_tiplabels,
                            column_w_names_of_traits,
                            column_w_values_for_traits
                        ) -> traits_wide
                    colnames(traits_wide) <- c("sample_unique_ID", "trait", "value")
                    traits_wide <- pivot_wider(traits_wide, names_from = "trait", values_from = "value")
    
                    traits_wide[is.na(traits_wide)] <- 0
                        
                ## Loop over traits and get ancestral states
                    anc_traits <- list()
                    for (i in 2:dim(traits_wide)[2]) { # i=2
                        trait <- as.numeric(unlist(traits_wide[,i]))
                        names(trait) <- as.character(unlist(traits_wide[,1]))

                        fastAnc_output <- fastAnc( tree = tree, x = trait, vars = FALSE, CI = TRUE )
                        anc_traits[[i]] <- data.frame(
                            node = as.numeric(names(fastAnc_output[[1]])),
                            trait = colnames(traits_wide)[i],
                            value = fastAnc_output[[1]]
                        )
                    }
                    anc_traits <- do.call(rbind, anc_traits)
                    treefort_2 <- right_join(treefort, anc_traits, by = "node")

                ## Put it all together (including analyte metadata!) and return
                    treefort_3 <- rbind(
                        treefort_1[,colnames(treefort_1) %in% colnames(treefort_2)],
                        treefort_2
                    )
                    
                    return(treefort_3)
            }


        #### geom_ancestral_pie


            #' Plots pie charts for data at any node using columns `trait` and `value`
            #' @param data Not yet documented.
            #' @param pie_size Not yet documented.
            #' @param pie_alpha Not yet documented.
            #' @param pie_fill_colors Not yet documented.
            #' @param pie_border_size Not yet documented.
            #' @param ... Not yet documented.
            #' @export
            geom_ancestral_pie <- function(
                data,
                pie_size = 0.1,
                pie_alpha = 1,
                pie_fill_colors = discrete_palette,
                pie_border_size = 0.1,
                ...
            ) {

                ## Data wrangling 
                    data <- as.data.frame(data)
                    if (!"node" %in% colnames(data)) { stop("data should have a column 'node'...") }
                    data %>% select(node, trait, value) -> ldf
                    ldf_split <- split(data, data$node)

                ## Make pies and apply colors   
                    pie_list <- list()
                    for (i in 1:length(ldf_split)) { #i=10
                    pie_list[[i]] <- ggplot() +
                        geom_col(
                            data = ldf_split[[i]],
                            aes(y = value, x = 1, fill = trait),
                            color = "black", linewidth = pie_border_size, alpha = pie_alpha
                        ) +
                        scale_fill_manual(
                            values = pie_fill_colors,
                            guide = "none",
                            na.value = "white"
                        ) +
                        coord_polar(theta = "y") + theme_void()
                    }
                    names(pie_list) <- names(ldf_split)

                    return(geom_inset(
                        pie_list, width = pie_size, height = pie_size
                    ))
            }

