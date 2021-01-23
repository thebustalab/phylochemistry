#######################
# PHYLOCHEMISTRY V1.0 #
#######################

##### Read and write functions

    #### readMonolist

        #' Reads a monolist
        #'
        #' @param monolist_in_path The path to the monolist (in .csv format) to be read
        #' @examples
        #' @export
        #' readMonolist

        readMonolist <- function( monolist_in_path ) {
            monolist <- as.data.frame(data.table::fread(file = monolist_in_path))
            return( monolist )
        }

    #### writeMonolist

        #' Writes a monolist
        #'
        #' @param monolist_out_path The path to where the monolist should be written
        #' @examples
        #' @export
        #' writeMonolist

        writeMonolist <- function( monolist, monolist_out_path ) {
            write.table(monolist, file = monolist_out_path, sep = ",", row.names = FALSE, col.names = TRUE)
        }

    #### readPolylist

        #' Reads a polylist
        #'
        #' Reads a spreadsheet in polylist format (wide format, multiple column and row headers), and converts it into tidy format
        #' @param polylist_in_path The path to the polylist (in .csv format) to be read
        #' @param centroid The characters in the cell that defines the boundaries of the multiple row and column headers
        #' @param table_value_unit The name to be given to the column that will contain the values in the polylist table
        #' @examples
        #' @export
        #' readPolylist

        readPolylist <- function(   
                            polylist_in_path,
                            centroid = "~~~",
                            table_value_unit = "abundance"
                        ) {

            # Check for centroid
                if ( length(centroid) != 1 ) {
                    stop("Please provide a centroid")
                }

            # Import polylist
                polylist <- as.data.frame(data.table::fread(polylist_in_path, header = FALSE))

            # Identify location of centroid
                center_column <- unlist(apply(polylist, 1, function(x) grep(centroid, x)))
                center_row <- grep(centroid, polylist[,center_column])

            # Use centroid to extract vertical_monolist
                vertical_monolist <- polylist[(center_row+1):dim(polylist)[1], 1:(center_column-1)]
                colnames(vertical_monolist) <- as.character(unlist(polylist[center_row, 1:(center_column-1)]))
                vertical_monolist$URI_URI_URI <- apply(vertical_monolist, 1, function(x) paste(x, collapse = ""))

            # Use centroid to extract horizontal_monolist
                horizontal_monolist <- as.data.frame(t(polylist[-c(center_row), (center_column+1):dim(polylist)[2]]))
                rownames(horizontal_monolist) <- NULL
                colnames(horizontal_monolist) <-    c(
                                                    as.character(unlist(polylist[c(1:(center_row-1), (center_row+1):dim(polylist)[1]), center_column]))[1:(center_row-1)],
                                                    as.character(vertical_monolist$URI_URI_URI)
                                                )
                horizontal_monolist <- tidyr::gather(horizontal_monolist, URI_URI_URI, table_value_unit, (center_row):dim(horizontal_monolist)[2])
                colnames(horizontal_monolist)[colnames(horizontal_monolist) == "table_value_unit"] <- table_value_unit

            # Bind the two monolists, drop the URI column
                polylist <- cbind(horizontal_monolist,vertical_monolist[match(horizontal_monolist$URI_URI_URI, vertical_monolist[,colnames(vertical_monolist) == "URI_URI_URI"]),])
                polylist <- polylist[,colnames(polylist) != "URI_URI_URI"]

            # Add Genus_species column if not present
                # if ( any(colnames(polylist) == "Genus_species") == FALSE) {
                #     print("Adding Genus_species column")
                #     polylist$Genus_species <- paste(polylist$Genus, polylist$species, sep="_")
                # }

            # Reset row numbers, return the polylist
                rownames(polylist) <- NULL
                return(polylist)
        }

    #### writePolylist

        #' Writes a polylist
        #'
        #' @param polylist_out_path The path to where the polylist should be written
        #' @examples
        #' @export
        #' writePolylist

        writePolylist <- function( polylist, polylist_out_path ) {
            write.table( polylist, file = polylist_out_path, sep = ",", row.names = FALSE, col.names = FALSE)
        }

    #### readTree

        #' Reads a phylogenetic tree
        #'
        #' @param tree_in_path The path to the phylogenetic tree
        #' @importFrom ape read.tree
        #' @examples
        #' @export
        #' readTree

        readTree <- function( tree_in_path ) {

            ## Read and return the tree
                return( ape::read.tree(file = tree_in_path) )
        }

    #### writeTree

        #' Writes a phylogenetic tree
        #'
        #' @param tree_out_path The path to where the tree should be written
        #' @importFrom ape write.tree
        #' @examples
        #' @export
        #' writeTree

        writeTree <- function( tree, tree_out_path ) {
            ape::write.tree( phy = tree, file = tree_out_path )
        }

    #### readFasta

    	#' Returns contents of one fasta file as a StringSet object.
        #'
        #' @param fasta_in_path The path to the fasta to read
        #' @param fasta_type The type of sequence data contained in the fasta file. Options: "nonspecific", "DNA", "RNA", or "AA"
        #' @importFrom Biostrings readBStringSet readDNAStringSet readRNAStringSet readAAStringSet
        #' @examples
        #' @export
        #' readFasta

        readFasta <- function( fasta_in_path, fasta_type = "nonspecific" ) {
        	
        	if ( fasta_type == "nonspecific" ) {
        		return(Biostrings::readBStringSet( filepath = fasta_in_path ))
        	}

        	if ( fasta_type == "DNA" ) {
        		return(Biostrings::readDNAStringSet( filepath = fasta_in_path ))
        	}

        	if ( fasta_type == "RNA" ) {
        		return(Biostrings::readRNAStringSet( filepath = fasta_in_path ))
        	}

        	if ( fasta_type == "AA" ) {
        		return(Biostrings::readAAStringSet( filepath = fasta_in_path ))
        	}
        }

    #### readManyFasta

    	#' Returns the contents of many fasta files as a StringSet object.
        #'
        #' @param fasta_in_paths The paths to the fastas to read
        #' @param fasta_type The type of sequence data contained in the fasta file. Options: "nonspecific", "DNA", "RNA", or "AA"
        #' @importFrom Biostrings readBStringSet readDNAStringSet readRNAStringSet readAAStringSet
        #' @examples
        #' @export
        #' readManyFasta

        readManyFasta <- function( fasta_in_paths, fasta_type = "nonspecific" ) {

    		if ( fasta_type == "nonspecific" ) {
        		temp <- BStringSet()
        	}

        	if ( fasta_type == "DNA" ) {
        		temp <- DNAStringSet()
        	}

        	if ( fasta_type == "RNA" ) {
        		temp <- RNAStringSet()
        	}

        	if ( fasta_type == "AA" ) {
        		temp <- AAStringSet()
        	}

        	for ( fasta in 1:length(fasta_in_paths) ) {

				if ( fasta_type == "nonspecific" ) {
	        		temp <- c(temp, Biostrings::readBStringSet( filepath = fasta_in_paths[fasta] ))
	        	}

	        	if ( fasta_type == "DNA" ) {
	        		temp <- c(temp, Biostrings::readDNAStringSet( filepath = fasta_in_paths[fasta] ))
	        	}

	        	if ( fasta_type == "RNA" ) {
	        		temp <- c(temp, Biostrings::readRNAStringSet( filepath = fasta_in_paths[fasta] ))
	        	}

	        	if ( fasta_type == "AA" ) {
	        		temp <- c(temp, Biostrings::readAAStringSet( filepath = fasta_in_paths[fasta] ))
	        	}

        	}

        	return( temp )

        }

    #### writeSingleFastas

    	#' Writes one or more sequences as individual fasta file(s), each with one sequence
        #'
        #' @param XStringSet The sequence(s) to write out
        #' @param fasta_out_directory_path The path to the directory where the file(s) should be written
        #' @importFrom Biostrings writeXStringSet
        #' @examples
        #' @export
        #' writeSingleFastas

        writeSingleFastas <- function( XStringSet, fasta_out_directory_path ) {
        	for (sequence in 1:length(XStringSet)) {
        		Biostrings::writeXStringSet( x = XStringSet[sequence], filepath = paste0(fasta_out_directory_path, XStringSet@ranges@NAMES[sequence], ".fa") )
        	}
        }

    #### writeFasta

    	#' Writes one or more sequences as a single fasta file
        #'
        #' @param XStringSet The sequences to write out
        #' @param fasta_out_path The path to where the file should be written
        #' @importFrom Biostrings writeXStringSet
        #' @examples
        #' @export
        #' writeFasta

        writeFasta <- function( XStringSet, fasta_out_path ) {
        	Biostrings::writeXStringSet( x = XStringSet, filepath = fasta_out_path )
        }

    #### writeSupplementalTable

        #' Writes a human-readable supplemental table of quantitative data
        #'
        #' @param supplementalTable The data to write out
        #' @param round_level The number of decimal places summary statistics should have
        #' @param file_out_path The path to the location where the file should be written
        #' @examples
        #' @export
        #' writeSupplementalTable

        writeSupplementalTable <- function(supplementalTable, round_level, file_out_path, replicates) {
    
            # Make new object from input
                supp_table <- supplementalTable
                supp_table$abundance <- round(supp_table$abundance, round_level)
            
            # Get number of compound levels
                number_of_compound_levels <- length(grep("compound", colnames(supp_table)))

            # Assign sample_unique_id
                supp_table$sample_unique_id <- paste(supp_table$sample, supp_table$replicate, sep = "..")

            # Check that sample levels and compound levels make unique IDs
                if (number_of_compound_levels == 2) {
                    
                    supp_table$compound_unique_id <- paste(supp_table$compound_level_1, supp_table$compound_level_2, sep = "..")
                }

                if ( sd(table(supp_table$compound_unique_id)) != 0 ) {
                    stop("Compound levels do not produce unique sample IDs")
                }
            
            # Spread a subset of the table "supp_table_min"
                supp_table_min <- supp_table[colnames(supp_table) %in% c("sample_unique_id", "abundance", "compound_unique_id")]
                supp_table_min <- tidyr::spread(supp_table_min, sample_unique_id, abundance)

            # Assign supp_table_min the compound levels, then remove compound_unique_id
                if (number_of_compound_levels == 2) {
                    supp_table_min$compound_level_1 <- gsub("\\.\\..*", "", supp_table_min$compound_unique_id)
                    supp_table_min$compound_level_2 <- gsub(".*\\.\\.", "", supp_table_min$compound_unique_id)
                }
                supp_table_min <- supp_table_min[,colnames(supp_table_min) != "compound_unique_id"]

            # Order the table according to the compound levels and then the sample names
                supp_table_min <- supp_table_min[order(supp_table_min$compound_level_1, supp_table_min$compound_level_2),]
                supp_table_min <- supp_table_min[,match(c("compound_level_1", "compound_level_2", unique(supp_table$sample_unique_id)), colnames(supp_table_min))]

            # Total of each compound class
                compound_level_1_totals <- list()
                for (i in 1:length(unique(supp_table_min$compound_level_1))) {
                    data <- supp_table_min[supp_table_min$compound_level_1 == unique(supp_table_min$compound_level_1)[i],]
                    data <- data[,colnames(data) %in% unique(supp_table$sample_unique_id)]
                    compound_level_1_totals[[i]] <- data.frame(
                        compound_level_1 = paste0("Total ", as.character(unique(supp_table_min$compound_level_1)[i])),
                        compound_level_2 = " ",
                        t(colSums(data))
                    )
                }
                compound_level_1_totals <- do.call(rbind, compound_level_1_totals)
            
            # Total abundnace for each sample
                total_abundance <- cbind(
                    data.frame(
                        compound_level_1 = "Total abundance",
                        compound_level_2 = " "
                    ),
                    t(colSums(supp_table_min[,colnames(supp_table_min) %in% unique(supp_table$sample_unique_id)]))
                )
                
            # Space between values and summary stats (total unk, total wax)
                space <- t(data.frame(rep(" ", dim(supp_table_min)[2])))
                colnames(space) <- colnames(supp_table_min)

            # Bind it all together, Make everything numeric, fix rownames
                supp_table_min <- rbind(
                    supp_table_min[supp_table_min$compound_level_1 != "UNK",],
                    space,
                    compound_level_1_totals,
                    space,
                    total_abundance
                )
                supp_table_min[,3:dim(supp_table_min)[2]] <- t(apply(supp_table_min[,3:dim(supp_table_min)[2]], 1, function(x) as.numeric(x)))
                rownames(supp_table_min) <- NULL

            # Percents

                for (i in 1:length(unique(supp_table$sample_unique_id))) {

                    data <- supp_table_min[,colnames(supp_table_min) == unique(supp_table$sample_unique_id)[i]]
                    supp_table_min$percent <- round(100*data/data[length(data)], round_level)

                    colnames(supp_table_min)[colnames(supp_table_min) == "percent"] <- paste0("percent_", unique(supp_table$sample_unique_id)[i])
                
                }

            # Make averages and stdevs

                colnames(supp_table_min)[colnames(supp_table_min) %in% unique(supp_table$sample_unique_id)] <- paste0(
                    "abundance_",
                    colnames(supp_table_min)[colnames(supp_table_min) %in% unique(supp_table$sample_unique_id)]
                )
                
                for (i in 1:length(unique(supp_table$sample)) ) {

                    supp_table_min$average_1 <- round(
                        apply(supp_table_min[,grep(as.character(paste0("abundance_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, mean), round_level
                    )
                    supp_table_min$error_1 <- round(
                        apply(supp_table_min[,grep(as.character(paste0("abundance_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, sd), round_level
                    )
                    supp_table_min$error_1 <- round(supp_table_min$error_1/sqrt(replicates), round_level)

                    colnames(supp_table_min)[colnames(supp_table_min) == "average_1"] <- paste0("abundance_", unique(supp_table$sample)[i], "_average")
                    colnames(supp_table_min)[colnames(supp_table_min) == "error_1"] <- paste0("abundance_", unique(supp_table$sample)[i], "_stderror")

                    supp_table_min$average_1 <- round(
                        apply(supp_table_min[,grep(as.character(paste0("percent_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, mean), round_level
                    )
                    supp_table_min$error_1 <- round(
                        apply(supp_table_min[,grep(as.character(paste0("percent_", unique(supp_table$sample)[i])), colnames(supp_table_min))], 1, sd), round_level
                    )
                    supp_table_min$error_1 <- round(supp_table_min$error_1/sqrt(replicates), round_level)

                    colnames(supp_table_min)[colnames(supp_table_min) == "average_1"] <- paste0("percent_", unique(supp_table$sample)[i], "_average")
                    colnames(supp_table_min)[colnames(supp_table_min) == "error_1"] <- paste0("percent_", unique(supp_table$sample)[i], "_stderror")

                }

            # Order the columns
                column_order <- vector()
                for ( i in 1:length(unique(supp_table$sample)) ) {
                    column_order <- c(column_order, colnames(supp_table_min)[grep(paste0("abundance_", unique(supp_table$sample)[i]), colnames(supp_table_min))])
                    column_order <- c(column_order, colnames(supp_table_min)[grep(paste0("percent_", unique(supp_table$sample)[i]), colnames(supp_table_min))])
                }

                column_order <- match(c("compound_level_1", "compound_level_2", column_order), colnames(supp_table_min))

                supp_table_min <- supp_table_min[,column_order]

                # tibble::add_column(supp_table_min, .after = c(2,6,9))
                # grep("stdev", colnames)

            # Finalization of table
                # supp_table <- cbind(supp_table, t(space), supp_table_percent)
                supp_table_min[supp_table_min == 0] <- "n.d."
                # gsub(" ", "_", colnames(supp_table))
                write.csv(supp_table_min, file = file_out_path, row.names = FALSE)

        }

##### Polylist construction and manipulation

    #### buildPolylist

        #' Build a human-readable matrix of samples, analytes, and measurements
        #'
        #' @param samples_monolist_in_path The monolist of samples to be used in building the polylist. Should contain a column "sample_unique_ID".
        #' @param analytes_monolist_in_path The monolist of analytes to be used in building the polylist. Should contain a column "analyte_unique_ID".
        #' @param measurements_monolist_in_path The monolist of measurements to be used in building the polylist. Should contain columns "sample_unique_ID", "analyte_unique_ID", and a column of values corresponding to measurements
        #' @param centroid The characters to use in the centroid
        #' @param polylist_out_path The path to which the polylist should be written
        #' @keywords lists
        #' @examples
        #' @export
        #' buildPolylist

        buildPolylist <-    function(  
                                samples_monolist_in_path,
                                analytes_monolist_in_path,
                                measurements_monolist_in_path = NULL,
                                centroid = "~~~",
                                polylist_out_path
                            ) {

            # Read in the monolists
                samples_monolist <- readMonolist(samples_monolist_in_path)
                analytes_monolist <- readMonolist(analytes_monolist_in_path)
                centroid_location <- c((dim(analytes_monolist)[2]+1), (dim(samples_monolist)[2]+1))

            # Build polylist
                samples_monolist <- as.data.frame(t(samples_monolist))
                samples_monolist <- cbind(rownames(samples_monolist), samples_monolist)

                for (i in 1:(dim(analytes_monolist)[2])) {
                    samples_monolist <- cbind("NA", samples_monolist)
                }
                colnames(samples_monolist) <- as.character(seq(1,dim(samples_monolist)[2]))
                samples_monolist <- apply(samples_monolist, 2, FUN = as.character)

                for (i in 1:(dim(samples_monolist)[2]-dim(analytes_monolist)[2])) {
                    analytes_monolist <- cbind(analytes_monolist, "NA")
                }
                analytes_monolist <- apply(analytes_monolist, 2, FUN = as.character)
                analytes_monolist <- rbind(colnames(analytes_monolist), analytes_monolist)
                colnames(analytes_monolist) <- as.character(seq(1,dim(analytes_monolist)[2]))
                
                polylist <- rbind(samples_monolist,analytes_monolist)

            # Fill in with measurements if specified
                if ( !is.null(measurements_monolist_in_path) ) {
                    measurements_monolist <- readMonolist(monolist = measurements_monolist_in_path)
                    colnames(measurements_monolist)[!colnames(measurements_monolist) %in% c("sample_unique_ID", "analyte_unique_ID")] <- "value"

                    analyte_unique_ID_col <- grep("analyte_unique_ID", polylist[centroid_location[2],])
                    sample_unique_ID_row <- grep("sample_unique_ID", polylist[,centroid_location[1]])

                    pb <- progress::progress_bar$new(total = dim(measurements_monolist)[1])
                    for (i in 1:dim(measurements_monolist)[1]) {
                        data_point <- measurements_monolist[i,]
                        polylist[
                            match(data_point$analyte_unique_ID, polylist[,analyte_unique_ID_col]),
                            match(data_point$sample_unique_ID, polylist[sample_unique_ID_row,])
                        ] <- as.character(data_point$value)
                        pb$tick()
                    }
                }
                
            # Clean up and write out polylist
                polylist[polylist %in% c("NA", "\"NA\"")] <- ""
                polylist[centroid_location[2], centroid_location[1]] <- centroid
                writePolylist(polylist = polylist, polylist_out_path = polylist_out_path)
        }

##### Tree and taxa manipulation

    #### buildTree

        #' Construct various types of phylogenetic trees from alignments or other trees
        #'
        #' @param scaffold_type The type of information that should be used to build the tree. One of "amin_alignment", "nucl_alignment", or "newick"
        #' @param scaffold_in_path The path to the information that should be used to build the tree.
        #' @param members The tips of the tree that should be included. Default is to include everything.
        #' @param gblocks TRUE/FALSE whether to use gblocks on the alignment
        #' @param gblocks_path Path to the gblocks module, perhaps something like "/Users/lucasbusta/Documents/Science/_Lab_Notebook/Bioinformatics/programs/Gblocks_0.91b/Gblocks"
        #' @param ml TRUE/FALSE whether to use maximum likelihood when constructing the tree.
        #' @param model_test TRUE/FALSE whether to test various maximum likelihood models while constructing the tree
        #' @param bootstrap TRUE/FALSE whether to calculate bootstrap values for tree nodes.
        #' @param rois TRUE/FALSE whether to use only certain portions of an alignment for constructing the tree
        #' @param rois_data The position in the alignment to use when constructing the tree. Default = NULL, i.e. all positions.
        #' @param ancestral_states TRUE/FALSE whether to calculate ancestral states at nodes in the tree. Requires specifying a root via the 'root' parameter
        #' @param root The tree tip to use as the root of the tree
        #' @examples
        #' @export
        #' buildTree

        buildTree <-    function(
                            scaffold_type = c("amin_alignment", "nucl_alignment", "newick"),
                            scaffold_in_path,
                            members = NULL,
                            gblocks = FALSE, 
                            gblocks_path = NULL,
                            ml = FALSE, 
                            model_test = FALSE,
                            bootstrap = FALSE,
                            rois = FALSE, 
                            rois_data = NULL,
                            ancestral_states = FALSE,
                            root = NULL
                        ) {

            if ( scaffold_type == "nucl_alignment" ) {

                ## Use ROIS if specified
                    # if ( rois == FALSE ) {
                    #     cat("Skipping ROIs...\n")
                    # } 
                    # if ( rois == TRUE ) {
                    #     rois_data$region_start <- as.numeric(as.character(rois_data$region_start))
                    #     rois_data$region_end <- as.numeric(as.character(rois_data$region_end))
                    #     nucl_seqs_aligned <- Biostrings::readDNAStringSet(file = paste(scaffold_in_path), format = "fasta")
                    #     nucl_seqs_aligned_roi <- subseq(nucl_seqs_aligned, 1, 0)
                    #         for ( i in 1:dim(rois_data)[1] ) {
                    #             nucl_seqs_aligned_roi <- Biostrings::xscat(nucl_seqs_aligned_roi, subseq(nucl_seqs_aligned, (rois_data[,2][i]*3)-2, (rois_data[,3][i]*3)))
                    #         }
                    #     nucl_seqs_aligned_roi@ranges@NAMES <- nucl_seqs_aligned@ranges@NAMES
                    #     Biostrings::writeXStringSet(nucl_seqs_aligned_roi, filepath = paste(phylochemistry_directory, "/", project_name, "/alignments/", scaffold_in_path, "_roi", sep = ""))
                    #     nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(phylochemistry_directory, "/", project_name, "/alignments/", scaffold_in_path, "_roi", sep = ""), format="fasta", type="DNA")
                    #     if ( gblocks == FALSE ) {
                    #         cat(paste("Making tree with ", scaffold_in_path, "_roi...", sep = ""))
                    #     }
                    # }

                ## Use gblocks on the alignment if specified
                    if ( gblocks == FALSE ) {
                        cat("Skipping gBlocks...\n")                    
                    }
                    if ( gblocks == TRUE ) {
                        if ( rois == FALSE ) {
                            nucl_seqs_aligned <- ape::read.dna(file = paste(scaffold_in_path), format = "fasta") # Needs to be DNAbin format
                            nucl_seqs_aligned_blocked <- ips::gblocks(nucl_seqs_aligned, b5 = "a", exec = gblocks_path)
                            ape::write.dna(nucl_seqs_aligned_blocked, paste(scaffold_in_path, "_blocked", sep = ""), format = "fasta")
                            nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path, "_blocked", sep = ""), format = "fasta", type = "DNA")
                            cat(paste("Making tree with ", scaffold_in_path, "_blocked...\n", sep = ""))
                        } 
                        if ( rois == TRUE ) {
                            # nucl_seqs_aligned <- seqinr::read.alignment(file = paste(scaffold_in_path, "_roi", sep = ""), format = "fasta")
                            # nucl_seqs_aligned_blocked <- ips::gblocks(nucl_seqs_aligned, b5 = "n", exec=gblocks_path)
                            # ape::write.dna(nucl_seqs_aligned_blocked, paste(phylochemistry_directory, "/", project_name, "/alignments/", scaffold_in_path, "_roi_blocked", sep = ""), format = "fasta")
                            # nucl_seqs_aligned <- phangorn::read.phyDat(file=paste(phylochemistry_directory, "/", project_name, "/alignments/", scaffold_in_path, "_roi_blocked", sep = ""), format="fasta", type="DNA")
                            # cat(paste("Making tree with ", scaffold_in_path, "_roi_blocked...", sep = ""))
                        }
                    }

                ## If neither gBlocks nor ROIs, read in alignment as phyDat for tree creation
                    if ( rois == FALSE ) {
                        if ( gblocks == FALSE ) {
                            nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path), format = "fasta", type = "DNA")
                            cat(paste("Making tree with ", scaffold_in_path," ...\n", sep = ""))
                        }
                    }

                ## Create the tree
                    output <- list()

                    ## Create distrance tree
                        cat("Creating neighbor-joining tree...\n")
                        dm <- phangorn::dist.ml(nucl_seqs_aligned, "F81")
                        NJ_tree <- phangorn::NJ(dm)
                        output <- NJ_tree

                    ## Make ml tree
                        if ( ml == TRUE ) {
                            ## Test all available nucl models, use the best one to optimize for ml
                                if ( model_test == TRUE ) {
                                    cat(paste("Testing 24 maximum liklihood models... \n"))
                                    mt <- phangorn::modelTest(nucl_seqs_aligned, tree = NJ_tree, multicore = TRUE)
                                    best_nucl_model <- gsub("\\+.*$", "", mt$Model[which.max(mt$logLik)])
                                    cat(paste("Tested 24 models, using best model:", as.character(gsub("\\+.*$","",best_nucl_model)), "\n", sep = " "))
                                } else {
                                    best_nucl_model <- "GTR"
                                }
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
                    return( output )
            }

            if ( scaffold_type == "amin_alignment" ) {

                ## Read in data
                    if (rois == FALSE) {
                            cat("Skipping roi...\n")
                            amin_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path), format = "fasta", type = "AA")
                            cat(paste("Making tree with ", scaffold_in_path, "...\n", sep = ""))
                        } else {
                            rois_data$region_start <- as.numeric(as.character(rois_data$region_start))
                            rois_data$region_end <- as.numeric(as.character(rois_data$region_end))
                            amin_seqs_aligned <- Biostrings::readAAStringSet(filepath = scaffold_in_path, format = "fasta")
                            amin_seqs_aligned_roi <- subseq(amin_seqs_aligned, 1, 0)
                            for (i in 1:dim(rois_data)[1]) {
                                amin_seqs_aligned_roi <- Biostrings::xscat(amin_seqs_aligned_roi, subseq(amin_seqs_aligned, rois_data[,2][i], rois_data[,3][i]))
                            }
                        amin_seqs_aligned_roi@ranges@NAMES <- amin_seqs_aligned@ranges@NAMES
                        Biostrings::writeXStringSet(amin_seqs_aligned_roi, filepath = paste(scaffold_in_path, "_roi", sep = ""))
                        amin_seqs_aligned <- phangorn::read.phyDat(file = paste(scaffold_in_path, "_roi", sep = ""), format = "fasta", type = "AA")
                        cat(paste("Making tree with ", scaffold_in_path, "_roi...", sep = ""))
                    }

                ## Make distance tree
                    dm = phangorn::dist.ml(amin_seqs_aligned, model = "JTT")
                    amin_tree = phangorn::NJ(dm)

                ## Make ml tree
                    if ( ml == TRUE ) {
                        if ( model_test == TRUE ) {
                            ## Test all available amino acid models and extract the best one
                                # mt <- modelTest(amin_seqs_aligned, model = "all", multicore = TRUE)
                                # fitStart = eval(get(mt$Model[which.min(mt$BIC)], env), env)
                                # LG+G+I
                        } else {
                            model <- "LG"
                        }
                        
                        ## Optimize for maximum liklihood and write out
                            fitStart = phangorn::pml(amin_tree, amin_seqs_aligned, model = model, k = 4, inv = .2)
                            amin_tree = phangorn::optim.pml(fitStart, rearrangement = "stochastic", optInv = TRUE, optGamma = TRUE)$tree
                    }

                ## Run bootstrap analysis
                    if ( bootstrap == TRUE ) {
                        if ( ml == FALSE ) {
                            stop("To calculate bootstrap values, please also run maximum likelihood estimation (ml = TRUE)")
                        }
                        bootstraps <- phangorn::bootstrap.pml(fitStart, bs = 100, optNni = TRUE, multicore = FALSE)
                        amin_tree$node.label <- phangorn::plotBS(amin_tree, bootstraps)$node.label
                    }

                ## Return tree
                    return ( amin_tree )
            }

            if ( scaffold_type == "newick" ) {

                ## Read in the newick scaffold
                    newick <- readTree( tree_in_path = scaffold_in_path )

                ## Are the Genus_species in your members in the newick? Are the genera in your members in the newick?
                    compatibility <- data.frame( Genus_species = unique(members), is_species_in_tree = NA, is_genus_in_tree = NA )
                    compatibility$is_species_in_tree <- compatibility$Genus_species %in% newick$tip.label
                    compatibility$is_genus_in_tree <- gsub("_.*$", "", compatibility$Genus_species) %in% gsub("_.*$", "", as.character(newick$tip.label))

                if ( all(compatibility$is_species_in_tree) == FALSE ) {
                    ## For Genus_species in members whose genus is missing from the tree (orphans), remove them
                        orphans <- as.character(dplyr::filter(compatibility, is_species_in_tree == FALSE & is_genus_in_tree == FALSE)$Genus_species)
                        members <- members[!(members %in% orphans)]
                        if ( length(orphans) > 0 ) {
                            cat("The following species belong to a genus not found in the newick scaffold and were removed: ")
                            for ( orphan in 1:length(orphans) ) {
                                cat("\n")
                                cat(orphans[orphan])
                            }
                            cat("\n")
                            cat("\n")
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
                                cat(paste("There aren't enough fosters to include the following species in the tree so it was removed:", adoptees[i], "\n", sep = " "))
                            } else {
                                cat(paste("Scaffold newick tip", available_foster_species[1], "substituted with", adoptees[i], "\n", sep = " "))
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
                    return ( newick )
            }

            cat("Pro tip: most tree read/write functions reset node numbers.\n 
                Fortify your tree and save it as a csv file to preserve node numbering.\n
                Do not save your tree as a newick or nexus file."
            )
        }

    #### pruneTree

        #' Prune a tree so it only contains user-specified tips
        #'
        #' @param tree The tree to manipulate
        #' @param tips_to_keep The tips to keep
        #' @importFrom ape drop.tip
        #' @examples
        #' @export
        #' pruneTree

        pruneTree <- function( tree_in_path, tips_to_keep ) {

            ## Read tree, drop all tips but those specified by user and return
                tree <- readTree(tree_in_path)
                tree <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% tips_to_keep])
                return( tree )
        }

##### Other functions

    #### normalize

        #' Normalizes a vector of numbers to a range of zero to one.
        #'
        #' @param x The vector to normalize
        #' @param old_min The minimum of the old range
        #' @param old_max The maximum of the old range
        #' @param new_min The minimum of the new range, defaults to 0
        #' @param new_max The maximum of the new range, defaults to 1
        #' @examples
        #' @export
        #' normalize

        normalize <- function( x, old_min = NULL, old_max = NULL, new_min = 0, new_max = 1 ) {

            if ( length(old_min) == 0 & length(old_max) == 0 ) {
                
                (x - min(x)) * (new_max - new_min) / (max(x) - min(x)) + new_min                
            
            } else {

                (x - old_min) * (new_max - new_min) / (old_max - old_min) + new_min

            }

        }

    #### OsDirectoryPathCorrect

        #' OS-aware correction of directory paths
        #'
        #' @param directory_in_path The path to correct
        #' @examples
        #' @export
        #' OsDirectoryPathCorrect

        OsDirectoryPathCorrect <- function( directory_path ) {

            OS <- .Platform$OS.type

            if (OS == "unix"){

                if ( 
                    substr(
                        directory_path, 
                        nchar(directory_path), 
                        nchar(directory_path)
                    ) == "/" ) {
                } else {
                    directory_path <- paste(directory_path, "/", sep = "")
                }
                directory_path_corrected <- directory_path

            } else if (OS == "windows"){

                if ( 
                    substr(
                        directory_path, 
                        nchar(directory_path), 
                        nchar(directory_path)
                    ) == "/" ) {
                } else {
                    directory_path <- paste(directory_path, "\\", sep = "")
                }
                directory_path_corrected <- directory_path

            } else {

                warning("ERROR: OS could not be identified")

            }

            return(directory_path_corrected)

        }

##### Sequence data handling

    #### extractORFs

        #' Extract open reading frames from multifasta files
        #'
        #' @param file The multifasta file to analyze
        #' @param write_out_ORFs TRUE/FALSE whether to write out a new fasta that contains the ORFs
        #' @param overwrite TRUE/FALSE Optionally overwrite the input file with the ORF file
        #' @examples
        #' @export
        #' extractORFs

        extractORFs <- function(file, write_out_ORFs = FALSE, overwrite = FALSE) {

            fasta <- seqinr::read.fasta(file)
            ORFs <- list()
            ORFs_to_write_out <- Biostrings::DNAStringSet()

            if ( write_out_ORFs == TRUE ) {
                if (file.exists(paste(file, "_orfs", sep = ""))) { file.remove(paste(file, "_orfs", sep = "")) }
            }

            for (j in 1:length(fasta)) {
                orf_coordinates <- data.frame()
                mRNA <- unlist(fasta[j])

                ## Search for ORFs in the forward sequence
                    data <- seqinr::translate(mRNA, frame = 0, sens = "F")
                        S <- grep("\\*", data) # find all stop codons
                        M <- grep("M", data) # find all start codons
                        if (length(S)>0 & length(M)>0) { # if there are stop codons, remove all start codons after the last stop codon
                            M <- M[!M > max(S)]
                        }
                        if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                            for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                orf_coordinates <-  rbind(orf_coordinates,  data.frame(
                                                                                start_codon = (M[i]*3-2), 
                                                                                stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                direction = "forward"
                                                                            )
                                                    )
                            }
                        }
                    data <- seqinr::translate(mRNA, frame = 1, sens = "F")
                        S <- grep("\\*", data) # find all stop codons
                        M <- grep("M", data) # find all start codons
                        if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                            M <- M[!M > max(S)] # remove all start codons after the last stop codon
                        }
                        if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                            for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                start_codon = (M[i]*3-2+1), 
                                                                                stop_codon = ((M[i]*3-2+1)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                direction = "forward"
                                                                            )
                                                    )
                            }
                        }
                    data <- seqinr::translate(mRNA, frame = 2, sens = "F")
                        S <- grep("\\*", data) # find all stop codons
                        M <- grep("M", data) # find all start codons
                        if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                            M <- M[!M > max(S)] # remove all start codons after the last stop codon
                        }
                        if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                            for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                start_codon = (M[i]*3-2+2),
                                                                                stop_codon = ((M[i]*3-2+2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                direction = "forward"
                                                                            )
                                                    )
                            }
                        }

                ####### Searching for ORFs in the reverse direction leads to issues during codon alignment...
                ## Search for ORFs in the reverse direction
                    # mRNA <- rev(mRNA)
                    # data <- seqinr::translate(mRNA, frame = 0, sens = "F")
                    #     S <- grep("\\*", data) # find all stop codons
                    #     M <- grep("M", data) # find all start codons
                    #     if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                    #         M <- M[!M > max(S)] # remove all start codons after the last stop codon
                    #     }
                    #     if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                    #         for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                    #             F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                    #             orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                    #                                                             start_codon = (M[i]*3-2), 
                    #                                                             stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                    #                                                             direction = "reverse"
                    #                                                         )
                    #                                 )
                    #         }
                    #     }
                    # data <- seqinr::translate(mRNA, frame = 1, sens = "F")
                    #     S <- grep("\\*", data) # find all stop codons
                    #     M <- grep("M", data) # find all start codons
                    #     if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                    #         M <- M[!M > max(S)] # remove all start codons after the last stop codon
                    #     }
                    #     if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                    #         for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                    #             F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                    #             orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                    #                                                             start_codon = (M[i]*3-2), 
                    #                                                             stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                    #                                                             direction = "reverse"
                    #                                                         )
                    #                                 )
                    #         }
                    #     }
                    # data <- seqinr::translate(mRNA, frame = 2, sens = "F")
                    #     S <- grep("\\*", data) # find all stop codons
                    #     M <- grep("M", data) # find all start codons
                    #     if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                    #         M <- M[!M > max(S)] # remove all start codons after the last stop codon
                    #     }
                    #     if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                    #         for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                    #             F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                    #             orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                    #                                                             start_codon = (M[i]*3-2), 
                    #                                                             stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                    #                                                             direction = "reverse"
                    #                                                         )
                    #                                 )
                    #         }
                    #     }

                ## Process discovered ORFs

                    if ( dim(orf_coordinates)[1] == 0 ) {
                        ORFs[[j]] <-    data.frame(
                                            start_codon = 0, 
                                            stop_codon = 0,
                                            direction = "forward",
                                            orf_length = 0
                                        )
                    } else {
                        rownames(orf_coordinates) <- NULL # Reset rownames because of previous step
                        orf_coordinates$orf_length <- (orf_coordinates$stop_codon - orf_coordinates$start_codon + 1)
                        orf_coordinates <- orf_coordinates[order(orf_coordinates$orf_length, decreasing = TRUE),]
                        ORFs[[j]] <- orf_coordinates
                    }

                    if ( write_out_ORFs == TRUE) {
                        if ( dim(orf_coordinates)[1] == 0 ) {
                            ORFs_to_write_out[j] <- c("a","t","g","t","a","a") ## Weird, but necessary so errors are not thrown during codon alignment
                            ORFs_to_write_out@ranges@NAMES[j] <- attr(fasta[j], "name")
                        } else {
                            if ( orf_coordinates$direction[1] == "forward" ) {
                                ORFs_to_write_out <- c(ORFs_to_write_out, Biostrings::DNAStringSet(paste(mRNA[orf_coordinates$start_codon[1]:orf_coordinates$stop_codon[1]], collapse = ""))) # extract longest orf from the mRNA
                                ORFs_to_write_out@ranges@NAMES[j] <- attr(fasta[j], "name")
                            }
                            # if ( orf_coordinates$direction[1] == "reverse" ) {
                            #     orf <- rev(mRNA)[orf_coordinates$start_codon[1]:orf_coordinates$stop_codon[1]] # extract longest orf from the mRNA
                            #     seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "a") # append ORF to the list
                            # }
                        }
                    }
            }

            # seqinr::write.fasta(fake_orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "a") # append fake_ORF to the list
            if ( write_out_ORFs == TRUE) {
                if (overwrite == TRUE) {
                    Biostrings::writeXStringSet(ORFs_to_write_out, file = file)
                    # seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "w") # overwrite transcript file
                } else {
                    Biostrings::writeXStringSet(ORFs_to_write_out, file = paste(file, "_orfs", sep = ""))
                    # seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "a") # append ORF to the list        
                }
            }

            ORFs <- do.call(rbind, ORFs)
            rownames(ORFs) <- NULL
            return(ORFs)
        }

    #### blastTranscriptomes

        #' BLAST search local transcriptomes
        #'
        #' Search locally stored transcriptomes for query sequences.
        #' @param transcriptomes A list of paths to the transcriptomes that should be searched, named by taxonomic identifier (e.g. species names)
        #' @param initial_query_in_path Path to a fasta file containing the query
        #' @param sequences_of_interest_directory_path Path to a directory where blast hits should be written out as fasta files
        #' @param blast_module_directory_path Path to directory containing the BLAST+ module (perhaps something like "/usr/local/ncbi/blast/bin/")
        #' @param blast_type One of "blastn" or "dc-megablast". "blastn" is a traditional BLASTN requiring an exact match of 11. "dc-megablast" is a discontiguous megablast used to find more distant (e.g., interspecies) sequences.
        #' @param remove Names of columns to remove from the output monolist
        #' @param monolist_out_path Path to where the output monolist should be written
        #' @examples
        #' @export
        #' blastTranscriptomes

        blastTranscriptomes <- function(
                                transcriptomes, 
                                initial_query_in_path, 
                                iterative_blast = FALSE,
                                iterative_blast_length_cutoff = 700,
                                sequences_of_interest_directory_path,
                                blast_module_directory_path,
                                blast_type = c("blastn", "dc-megablast"), 
                                remove = NULL,
                                monolist_out_path
                            ) {

            ### Check paths
                # sequences_of_interest_directory_path <- OsDirectoryPathCorrect(sequences_of_interest_directory_path)
                # blast_module_directory_path <- OsDirectoryPathCorrect(blast_module_directory_path)

            ### Make sure transcriptomes object is character and build BLAST databases
                names <- names(transcriptomes)
                transcriptomes <- as.character(transcriptomes)
                names(transcriptomes) <- names
                
                for (transcriptome in 1:length(transcriptomes)) {
                    cat(paste0("\nSpecies: ", names(transcriptomes)[transcriptome]))
                    system( 
                        paste(
                            blast_module_directory_path, 
                            "makeblastdb -in ", 
                            transcriptomes[transcriptome], 
                            " -dbtype nucl", 
                            sep = "" 
                        ) 
                    )
                }

            ### Read in query, set initial hit number for iterative BLAST
                query_seqs <- Biostrings::readDNAStringSet(filepath = initial_query_in_path, format = "fasta")

                if ( iterative_blast ) {
                    change_in_hit_number <- 1
                    initial_hit_number <- length(query_seqs)
                    iteration_number <- 0
                    cat("\nStarting iterative BLAST process\n\n")
                } else {
                    change_in_hit_number <- 1
                }

            ### Start BLAST process

                while ( change_in_hit_number > 0 ) {

                    ## Iteration number
                        if ( iterative_blast ) {
                            iteration_number <- iteration_number + 1 
                            cat(paste0("Iteration number: ", iteration_number, "\n"))
                        }

                    ## Loop over each member of the query and use it to blast each transcriptome
                        monolist <- data.frame()
                        for ( query_seq in 1:length(query_seqs) ) {

                            ## Write individual files for each member of the query
                                Biostrings::writeXStringSet(
                                    query_seqs[query_seq], 
                                    filepath = paste(initial_query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""), 
                                    append = FALSE
                                )

                            ## Loop over the transcriptomes, run the blast on each, add hits to monolist
                                for (transcriptome in 1:length(transcriptomes)) {

                                    ## Run BLAST
                                        system(
                                            paste(
                                                blast_module_directory_path,
                                                "blastn -task ",
                                                blast_type,
                                                " -db ", 
                                                transcriptomes[transcriptome],
                                                " -query ",
                                                paste(initial_query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                " -out ",
                                                paste(transcriptomes[transcriptome], ".out", sep = ""),
                                                # " -outfmt '6 sacc'",
                                                " -outfmt '6 sallacc'",
                                                sep = ""
                                            )
                                        )

                                    ## Extract BLAST hits from transcriptome, add them to the monolist, write them to individual files

                                        # Read in whole transcriptome
                                            temp_seqs <- Biostrings::readDNAStringSet( 
                                                filepath = as.character(transcriptomes[transcriptome]), 
                                                format = "fasta"
                                            )
                                    
                                        # Attempt to read in hits list, subset transcriptome according to that list, add hits to monolist
                                            cat(paste0("Query ", query_seq, ": "))
                                            if ( inherits( try( read.table(file = paste(transcriptomes[transcriptome], ".out", sep = "")), silent = TRUE), "try-error") == TRUE ) { # Skips over empty files
                                                cat(paste("No BLAST hits found for ", query_seqs[query_seq]@ranges@NAMES, " in ", names(transcriptomes)[transcriptome], "\n", sep = ""))
                                            } else {
                                                temp_hits <- read.table(file = paste(transcriptomes[transcriptome], ".out", sep = ""))
                                                cat(paste(
                                                    "Found ", 
                                                    dim(temp_hits)[1], 
                                                    " BLAST hits found for ", 
                                                    query_seqs[query_seq]@ranges@NAMES, 
                                                    " in ", 
                                                    names(transcriptomes)[transcriptome], 
                                                    "\n", 
                                                    sep = "")
                                                )

                                                accession <- gsub(" .*", "", temp_seqs@ranges@NAMES)
                                                temp_seqs <- temp_seqs[accession %in% temp_hits[,1]]
                                                accession <- gsub(" .*", "", temp_seqs@ranges@NAMES)

                                                if ( all(temp_hits[,1] %in% accession) != TRUE ) {
                                                    warning("Couldn't find some BLAST hits within the transcriptome!")
                                                }

                                                # If there are hits, add them to the monolist
                                                    if (length(temp_seqs) > 0) {

                                                        # Write out blast hits
                                                            longest_ORFs <- vector()
                                                            for (temp_seq in 1:length(temp_seqs)) {
                                                                Biostrings::writeXStringSet(
                                                                    temp_seqs[temp_seq], 
                                                                    filepath = paste(sequences_of_interest_directory_path, accession[temp_seq], ".fa", sep = ""), 
                                                                    append = FALSE
                                                                )
                                                                longest_ORFs <- c(longest_ORFs, extractORFs(paste(sequences_of_interest_directory_path, accession[temp_seq], ".fa", sep = ""))$orf_length[1])
                                                            }

                                                        # Add hit information to the monolist
                                                            monolist <- rbind(monolist, data.frame(
                                                                                                accession = gsub(" .*", "", temp_seqs@ranges@NAMES),
                                                                                                Genus = as.character(gsub("_.*$", "", names(transcriptomes)[transcriptome])),
                                                                                                species = gsub(".*_", "", names(transcriptomes)[transcriptome]),
                                                                                                annotation = temp_seqs@ranges@NAMES,
                                                                                                length = temp_seqs@ranges@width,
                                                                                                subset_all = TRUE,
                                                                                                query = query_seqs@ranges@NAMES[query_seq],
                                                                                                longestORF = longest_ORFs
                                                                                            )
                                                                        )
                                                    }
                                            }
                                }
                            ## Remove individual query files
                                if (file.exists(paste(query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))) {file.remove(paste(query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))}
                        }

                    ## Write out the monolist

                        ## Optionally remove certain columns from the monolist
                            if ( length(remove) > 0  ) {
                                monolist <- monolist[,!colnames(monolist) %in% remove]
                            }

                        monolist <- unique(monolist)
                        rownames(monolist) <- NULL
                        if (file.exists(monolist_out_path)) {file.remove(monolist_out_path)}
                        writeMonolist(monolist = monolist, monolist_out_path = monolist_out_path)

                    if ( iterative_blast ) {
                        
                        ## Should iterations continue?
                            monolist <- readMonolist(monolist_out_path)
                            final_hit_number <- dim(monolist[monolist$longestORF > iterative_blast_length_cutoff,])[1]
                            change_in_hit_number <- final_hit_number - initial_hit_number
                            cat(paste0(query_seq, ": "))
                            cat(paste0(change_in_hit_number, " new hits above cutoff length!\n\n"))
                            initial_hit_number <- final_hit_number

                    } else {
                        change_in_hit_number <- 0
                    }

                }
           
                cat("\nDone!\n\n")
        }

    #### alignSequences

        #' Align sequences in a seqlist
        #'
        #' @param monolist Monolist of the sequences to be aligned. First column should be "accession"
        #' @param subset TRUE/FALSE column in monolist that specifies which sequences should be included in the alignment
        #' @param alignment_directory_path Path to where the alignment should be written
        #' @param sequences_of_interest_directory_path Path to a directory where blast hits should be written out as fasta files
        #' @param input_sequence_type One of "nucl" or "amin"
        #' @param mode One of "nucl_align", "amin_align", or "codon_align"
        #' @param base_fragment
        #' TROUBLESHOOTING: 
        #'      "ERROR: inconsistency between the following pep and nuc seqs" - usually means there are duplicate accessions numbers in the input monolist
        #' @import msa
        #' @examples
        #' @export
        #' alignSequences

        alignSequences <-   function(
                                monolist, 
                                subset, 
                                alignment_directory_path, 
                                sequences_of_interest_directory_path, 
                                input_sequence_type = c("nucl", "amin"), 
                                mode = c("nucl_align", "amin_align", "codon_align", "fragment_align"),
                                base_fragment = NULL
                            ){  

            ## Check directory_paths
                sequences_of_interest_directory_path <- OsDirectoryPathCorrect(sequences_of_interest_directory_path)
                alignment_directory_path <- OsDirectoryPathCorrect(alignment_directory_path)

            ## Get appropriate subset of the monolist
                monolist_subset <- monolist[monolist[,colnames(monolist) == as.character(subset)],]

            ## If starting with nucleotide sequences
                
                if ( input_sequence_type == "nucl" ) {

                    ## Remove existing version of this alignment and it's files
                        if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs.fa", sep = ""))) {
                            file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs.fa", sep = ""))}
                        if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_aligned.fa", sep = ""))) {
                            file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_aligned.fa", sep = ""))}

                    ## Import the subset's sequences, then write out *_nucl_seqs.fa
                        nucl_seqs_set <- Biostrings::DNAStringSet()
                        for (i in 1:dim(monolist_subset)[1]) {
                            nucl_seqs_set <- c(nucl_seqs_set, Biostrings::readDNAStringSet(paste(sequences_of_interest_directory_path, as.character(monolist_subset$accession[i]), ".fa", sep = "")))
                        }
                        Biostrings::writeXStringSet(nucl_seqs_set, filepath = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""))

                    ## Run plain nucleotide alignment, if requested
                        if ( mode == "nucl_align") {
                            nucl_seqs_set_aligned <- msa::msa(nucl_seqs_set, order = c("input"))
                            msa <- msa::msaConvert(nucl_seqs_set_aligned, type = "seqinr::alignment")
                            n <- dim(as.data.frame(msa$seq))[1]
                            for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                            seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_nucl_seqs_aligned.fa", sep = ""))
                        }

                    ## Run codon alignment, if requested
                        if ( mode == "codon_align") {

                            ## Remove existing codon alignment
                                if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_codon_aligned.fa", sep = ""))) {
                                file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_codon_aligned.fa", sep = ""))}
                            
                            ## Extarct ORFs
                                extractORFs(file = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), write_out_ORFs = TRUE)
                            
                            ## Translate the nucleotide sequences and write out *_amin_seqs.fa
                                nucl_seqs_set <- Biostrings::readDNAStringSet(paste(alignment_directory_path, subset, "_nucl_seqs.fa_orfs", sep = ""))
                                amin_seqs_set <- Biostrings::translate(nucl_seqs_set, if.fuzzy.codon = "solve")
                                if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
                                Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

                            ## Run amino acid alignment and write out *_amin_seqs.fa
                                amin_seqs_set_aligned <- msa::msa(amin_seqs_set, order = c("input"))
                                msa <- msa::msaConvert(amin_seqs_set_aligned, type = "seqinr::alignment")
                                n <- dim(as.data.frame(msa$seq))[1]
                                for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                                if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))) {
                                file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))}
                                seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_amin_seqs_aligned.fa", sep = ""))

                            ## Run codon alignment and write to nucl_seqs_aligned.fa
                                nucl_seqs_codon_aligned <-  orthologr::codon_aln(
                                                                file_aln = paste(alignment_directory_path, subset, "_amin_seqs_aligned.fa", sep = ""), 
                                                                file_nuc = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), 
                                                                get_aln = TRUE
                                                            )
                                msa <- nucl_seqs_codon_aligned
                                n <- dim(as.data.frame(msa$seq))[1]
                                for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                                seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_nucl_seqs_codon_aligned.fa", sep = ""))
                        }

                    ## Run fragment alignment, if requested
                        if ( mode == "fragment_align") {

                            ## Remove existing fragment alignment
                                if (file.exists(paste(alignment_directory_path, as.character(subset), "_fragment_seqs_aligned.fa", sep = ""))) {
                                file.remove(paste(alignment_directory_path, as.character(subset), "_fragment_seqs_aligned.fa", sep = ""))}

                            ## Extarct ORFs
                                extractORFs(file = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), write_out_ORFs = TRUE)
                            
                            ## Translate the nucleotide sequences and write out *_amin_seqs.fa
                                nucl_seqs_set <- Biostrings::readDNAStringSet(paste(alignment_directory_path, subset, "_nucl_seqs.fa_orfs", sep = ""))
                                amin_seqs_set <- Biostrings::translate(nucl_seqs_set, if.fuzzy.codon = "solve")
                                if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
                                Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

                            ## Define base fragment
                                base_fragment_seq <- amin_seqs_set[amin_seqs_set@ranges@NAMES == base_fragment]
                                fragments_seq_set <- amin_seqs_set[amin_seqs_set@ranges@NAMES != base_fragment]
                                current_base_fragment_seq <- base_fragment_seq
                            
                            ## Align each fragment on its own with the base fragment
                                aligned_fragments <- AAStringSet()
                                for ( fragment in 1:length(fragments_seq_set) ) {
                                    fragment_pair_seq_set <- c(base_fragment_seq, fragments_seq_set[fragment])
                                    fragment_pair_seq_set_aligned <- msa::msa(fragment_pair_seq_set, order = c("input"))
                                    fragment_pair_seq_set_aligned <- AAStringSet(fragment_pair_seq_set_aligned)
                                    current_base_fragment_seq <- fragment_pair_seq_set_aligned[fragment_pair_seq_set_aligned@ranges@NAMES == base_fragment]
                                    aligned_fragments <- c(aligned_fragments, fragment_pair_seq_set_aligned[fragment_pair_seq_set_aligned@ranges@NAMES != base_fragment])
                                }

                            ## Perform final alignment and write it out
                                fragment_seq_set <- c(current_base_fragment_seq, aligned_fragments)
                                fragment_seq_set <- AAStringSet(gsub("-", "Z", fragment_seq_set))
                                fragment_seq_set_aligned <- msa::msa(fragment_seq_set, order = c("input"))
                                fragment_seq_set_aligned <- AAStringSet(fragment_seq_set_aligned)
                                writeXStringSet(fragment_seq_set_aligned, filepath = paste(alignment_directory_path, subset, "_fragment_seqs_aligned.fa", sep = ""))
                        }
                }

                if ( input_sequence_type == "amin" ) {

                    if ( mode == "amin_align" ) {

                        ## Remove existing version of this alignment and it's files
                            if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
                                file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
                            if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))) {
                                file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs_aligned.fa", sep = ""))}

                        ## Import the subset's sequences, then write out *_amin_seqs.fa
                            amin_seqs_set <- Biostrings::DNAStringSet()
                            for (i in 1:dim(monolist_subset)[1]) {
                                amin_seqs_set <- c(amin_seqs_set, Biostrings::readDNAStringSet(paste(sequences_of_interest_directory_path, as.character(monolist_subset$accession[i]), ".fa", sep = "")))
                            }
                            Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

                        ## Run the alignment
                            amin_seqs_set_aligned <- msa::msa(amin_seqs_set, order = c("input"))
                            msa <- msa::msaConvert(amin_seqs_set_aligned, type = "seqinr::alignment")
                            n <- dim(as.data.frame(msa$seq))[1]
                            for (i in 1:n){msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))}
                            seqinr::write.fasta(as.list(msa$seq),as.list(msa$nam), file.out = paste(alignment_directory_path, subset, "_amin_seqs_aligned.fa", sep = ""))
                    }
                }
        }

    #### readGFFs

        #' readGFFs
        #'
        #' Read multiple genome feature files into a tidy dataframe
        #' @param input_frame A dataframe with columns: Genus_species, GFF_in_path, region_reach, concatenation_spacing
        #' @param query A .fa file containing the query
        #' @param phylochem The phylochem object to which the data should be bound
        #' @examples
        #' @export
        #' readGFFs

            readGFFs <- function(
                            input_frame, # Dataframe with columns
                            subset_mode = c("none", "goi_only", "foi_only", "goi_region", "foi_region", "name_check"), 
                            goi = NULL, # Character vector of gene names. Labelling by species is not necessary
                            foi = NULL,
                            region_reach = 200000,
                            concatenate_by_species = TRUE,
                            concatenation_spacing = 10000,
                            omit = NULL
                        ) {

                framed_GFFs <- list()

                for ( i in 1:dim(input_frame)[1] ) {

                    # Load GFFs
                        cat(paste("Loading GFF for ", input_frame$Genus_species[i], "...\n", sep = ""))
                        gff_as_GRange <- rtracklayer::import.gff(as.character(input_frame$GFF_in_path[i]))
                        chr <- unique(gff_as_GRange@seqnames@values)
                        cat(paste("Number of chromosomes and/or scaffolds: ", length(chr), "\n", sep = ""))

                    # Make the GFF into dataframes "chrom_scaff", and "ranges"
                        library(GenomicRanges)
                        chrom_scaff <- as.data.frame(GenomicRanges::seqinfo(gff_as_GRange))
                        chrom_scaff$chrom_scaff_name <- rownames(chrom_scaff)
                        chrom_scaff$chrom_scaff_number <- seq(1,dim(chrom_scaff)[1],1)
                        rownames(chrom_scaff) <- NULL

                        ranges <- as.data.frame(gff_as_GRange)
                        ranges <- ranges[,colnames(ranges) != "Parent"]
                        
                    # Custom gene name editing on species-by-species basis
                        if ( input_frame$Genus_species[i] == "Zea_mays" ) {
                            # ranges$Name <- gsub("_.*$", "", ranges$Name)
                            ranges$Name <- ranges$locus_tag
                        }

                        if ( input_frame$Genus_species[i] == "Helianthus_annuus" ) {
                            ranges$Name <- ranges$ID
                        }

                        if ( input_frame$Genus_species[i] == "Solanum_lycopersicum" ) {
                            ranges$Name <- gsub("0\\.1", "0", ranges$Name)
                            ranges$Name <- gsub("0\\.2", "0", ranges$Name)
                            ranges$Name <- gsub("0\\.3", "0", ranges$Name)
                        }

                    # No subsetting
                        if ( subset_mode == "none" ) {
                            subsetted_ranges <- ranges
                        }

                    # Name check mode
                        if ( subset_mode == "name_check" ) {
                            framed_GFFs[[i]] <- unique(data.frame(
                                                    Genus_species = input_frame$Genus_species[i],
                                                    names = ranges$Name
                                                ))
                        }

                    # Subset "chrom_scaff" and "ranges"
                        if ( subset_mode == "goi_region" ) {

                            # Set up the genes_of_interest_frame for this species
                                goi_sp <- goi[goi %in% ranges$Name]
                                goi_sp_frame <- data.frame(goi_sp = goi_sp, chrom_scaff_of_interest = ranges$seqnames[match(goi_sp, ranges$Name)], to_subset = TRUE)

                            # Remove chrom_scaffs and ranges not on chrom_scaffs that contain goi
                                chrom_scaff <- chrom_scaff[chrom_scaff$chrom_scaff_name %in% goi_sp_frame$chrom_scaff_of_interest,]
                                ranges <- ranges[ranges$seqnames %in% goi_sp_frame$chrom_scaff_of_interest,]
                                
                            # Remove cols to omit
                                if ( length(omit) > 0 ) {
                                    ranges <- ranges[, !colnames(ranges) %in% omit]    
                                }
                                
                            # Remove ranges not within region_reach of a goi, if goi are within region_reach of eachother, expand the size of that region iteratively
                                subsetted_ranges <- data.frame()
                                goi_region_number <- 1
                                while ( any(goi_sp_frame$to_subset) == TRUE ) {
                                    j <- which(goi_sp_frame$to_subset == TRUE)[1]
                                    
                                    # Define the start of the goi_start_sites and the goi inside it
                                        ranges_on_this_chrom_scaff <- ranges[ranges$seqnames == goi_sp_frame$chrom_scaff_of_interest[j],]
                                        goi_start_sites <- ranges_on_this_chrom_scaff[ranges_on_this_chrom_scaff$Name %in% goi_sp_frame$goi_sp[j],]$start
                                        known_goi_in_this_region <- vector()
                                        known_goi_in_this_region <- c(known_goi_in_this_region, as.character(goi_sp_frame$goi_sp[j]))
                                        region <- data.frame(start = min(goi_start_sites) - input_frame$region_reach[i], end = max(goi_start_sites) + input_frame$region_reach[i])
                                        region$start[region$start < 0] <- 0

                                    # Expand the range if it includes another goi
                                        finished <- FALSE
                                        while ( finished == FALSE) {
                                            goi_region_bounds <- vector()
                                            for (k in 1:dim(region)[1]) {
                                                goi_region_bounds <- c(goi_region_bounds,seq(region[k,1], region[k,2], 1))
                                            }
                                            ranges_in_this_goi_region <- ranges_on_this_chrom_scaff[ranges_on_this_chrom_scaff$start %in% goi_region_bounds,]
                                            goi_in_current_goi_region <- goi_sp_frame$goi_sp[goi_sp_frame$goi_sp %in% ranges_in_this_goi_region$Name]
                                            if ( all(goi_in_current_goi_region %in% known_goi_in_this_region) ) {
                                                finished <- TRUE
                                            } else {
                                                temp <- ranges_on_this_chrom_scaff$Name %in% goi_in_current_goi_region[!goi_in_current_goi_region %in% known_goi_in_this_region]
                                                addl_gene <- ranges_on_this_chrom_scaff[temp,]
                                                addl_gene <- addl_gene[1,]
                                                goi_start_sites <- c(goi_start_sites, addl_gene$start)
                                                region <- data.frame(start = min(goi_start_sites) - input_frame$region_reach[i], end = max(goi_start_sites) + input_frame$region_reach[i])
                                                region$start[region$start < 0] <- 0

                                                # Mark that the addl_gene is now known in this region
                                                    known_goi_in_this_region <- c(known_goi_in_this_region, as.character(addl_gene$Name))
                                            }
                                        }
                                    goi_sp_frame$to_subset[goi_sp_frame$goi_sp %in% goi_in_current_goi_region] <- FALSE
                                    ranges_in_this_goi_region$goi_region_name <- paste(input_frame$Genus_species[i], goi_region_number, sep="_")
                                    ranges_in_this_goi_region$goi_region_number <- goi_region_number
                                    goi_region_number <- goi_region_number + 1
                                    subsetted_ranges <- rbind(subsetted_ranges, ranges_in_this_goi_region)
                                }

                            # Reorient each goi_region so it begins at zero
                                subsetted_ranges <- subsetted_ranges[subsetted_ranges$type != "region",] # remove "ranges" that are whole chrom_scaffs
                                subsetted_ranges <- plyr::ddply(subsetted_ranges, .(goi_region_name), mutate, end = end - min(start), start = start - min(start))
                                subsetted_ranges <- plyr::ddply(subsetted_ranges, .(goi_region_name), mutate, goi_region_length = max(end))

                            # Concatenate range_scaffs of this species
                                if ( concatenate_by_species == TRUE ) {
                                    if ( length(unique(subsetted_ranges$goi_region_number)) > 1 ) {
                                        end_previous_goi_region <- max(subsetted_ranges[subsetted_ranges$goi_region_number == 1,]$end)
                                        for (this_goi_region_number in 2:length(unique(subsetted_ranges$goi_region_number))) {
                                            amount_to_advance_this_goi_region <- end_previous_goi_region + input_frame$concatenation_spacing[i]
                                            subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$start <- subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$start + amount_to_advance_this_goi_region
                                            subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$end <- subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$end + amount_to_advance_this_goi_region
                                            end_previous_goi_region <- max(subsetted_ranges[subsetted_ranges$goi_region_number == this_goi_region_number,]$end)
                                        }
                                    }
                                }

                            # Bind to output
                                subsetted_ranges$Genus_species <- input_frame$Genus_species[i]
                                subsetted_ranges$Genus <- gsub("_.*$", "", input_frame$Genus_species[i])
                                subsetted_ranges$species <- gsub(".*_", "", input_frame$Genus_species[i])
                                subsetted_ranges$y_position <- input_frame$y_position[i]
                                framed_GFFs[[i]] <- subsetted_ranges
                        }
                }

                framed_GFFs <- do.call(plyr::rbind.fill, framed_GFFs)

                if ( subset_mode != "name_check" ) {

                    center <- (max(framed_GFFs$end) - min(framed_GFFs$start))/2
                    framed_GFFs <- plyr::ddply(framed_GFFs, .(Genus_species), mutate, center_start = (start + (center-(max(end)-min(start))/2)), center_end = (end + (center-(max(end)-min(start))/2)) )

                }

                return(framed_GFFs)
            }

##### Chemical data handling

    #### convertCDFstoCSVs

        #' Convert mass spectral datafiles (CDF) into a csv file
        #'
        #' @param paths_to_cdfs Paths to CDF files
        #' @param min_mz Smallest m/z value to process (acts like a filter)
        #' @param max_mz Largest m/z value to process (acts like a filter)
        #' @param min_rt (optional) Smallest retention time to process (acts like a filter)
        #' @param max_rt (optional) Largest retention time to process (acts like a filter)
        #' @param force Convert all CDF files, even if they've already been converted previously
        #' @examples
        #' @export
        #' convertCDFstoCSVs

            convertCDFstoCSVs <- 	function(
						            	paths_to_cdfs, 
						            	min_mz = 50, 
						            	max_mz = 800, 
						            	min_rt = NULL, 
						            	max_rt = NULL, 
						            	force = FALSE
						            ) {

                if ( force == FALSE ) {
                    paths_to_cdfs <- paths_to_cdfs[!file.exists(paste0(paths_to_cdfs, ".csv"))]
                }

                if (length(paths_to_cdfs) > 0) {

	                for (file in 1:length(paths_to_cdfs)) {

	                    cat(paste("Reading data file ", paths_to_cdfs[file], "\n", sep = ""))
	                        rawDataFile <- xcms::loadRaw(xcms::xcmsSource(paths_to_cdfs[file]))

	                    cat("   Framing data file... \n")
	                        rt <- rawDataFile$rt
	                        scanindex <- rawDataFile$scanindex

	                        filteredRawDataFile <- list()
	                        for ( i in 1:(length(rt)-1) ) {
	                            filteredRawDataFile[[i]] <- data.frame(
	                                mz = rawDataFile$mz[(scanindex[i]+1):(scanindex[i+1])],
	                                intensity = rawDataFile$intensity[(scanindex[i]+1):(scanindex[i+1])],
	                                rt = rt[i]
	                            )
	                        }
	                        framedDataFile <- do.call(rbind, filteredRawDataFile)
	                        framedDataFile$mz <- round(framedDataFile$mz, digits = 1)

	                    cat("   Merging duplicate rows ...\n")
	                        library(plyr)
	                        if ( dim(table(duplicated(paste(framedDataFile$mz, framedDataFile$rt, sep = "_")))) > 1 ) {
	                            framedDataFile <- plyr::ddply(framedDataFile, .(mz, rt), summarize, intensity = sum(intensity))
	                        }

	                    cat("   Filling in blank m/z values ...\n")
	                        spreadDataFile <- tidyr::spread(framedDataFile, rt, intensity, fill = 0)
	                        
	                        numbers <- round(seq(min_mz, max_mz, 0.1), digits = 1)
	                        add <- numbers[!numbers %in% spreadDataFile$mz]
	                        addition <- spreadDataFile[1:length(add),]
	                        
	                        if ( length(add) > 0 ) {
	                            addition$mz <- add     
	                        }
	                        
	                        addition <- tidyr::gather(addition, rt, intensity, 2:dim(addition)[2])
	                        addition$intensity <- 0
	                        addition$rt <- as.numeric(addition$rt)
	                        
	                        framedDataFile <- tidyr::gather(spreadDataFile, rt, intensity, 2:dim(spreadDataFile)[2])
	                        framedDataFile <- rbind(framedDataFile, addition)
	                        framedDataFile$rt <- as.numeric(framedDataFile$rt)
	                        framedDataFile <- framedDataFile[sort.list(framedDataFile$mz),]
	                        framedDataFile <- framedDataFile[sort.list(framedDataFile$rt),]
	                        rownames(framedDataFile) <- NULL

	                    if ( length(min_rt) > 0 ) {
	                        cat("   Filtering by minimum retention time thresholds ...\n")
	                            framedDataFile <- dplyr::filter(framedDataFile, rt > min_rt)
	                    }

	                    if ( length(max_rt) > 0 ) {
	                        cat("   Filtering by maximum retention time thresholds ...\n")
	                            framedDataFile <- dplyr::filter(framedDataFile, rt < max_rt)    
	                    }

	                    cat("   Writing out data file as CSV... \n")
	                        data.table::fwrite(framedDataFile, file = paste(paths_to_cdfs[file], ".csv", sep = ""), col.names = TRUE, row.names = FALSE)
	                }

	            }
            }

    #### extractChromatogramsFromCSVs

        #' Extract total ion chromatograms from csv files containing mass spectral data
        #'
        #' @param paths_to_cdf_csvs Paths to the cdf.csvs
        #' @examples
        #' @export
        #' extractChromatogramsFromCSVs

            extractChromatogramsFromCSVs <- function( paths_to_cdf_csvs, force = FALSE ) {

                chromatograms <- list()
                for ( file in 1:length(paths_to_cdf_csvs) ) {

                    cat(paste("Reading data file ", paths_to_cdf_csvs[file], "\n", sep = ""))
                        framedDataFile <- as.data.frame(data.table::fread(paths_to_cdf_csvs[file]))

                    cat("   Extracting the total ion chromatogram...\n")
                        library(plyr)
                        chromatogram <- plyr::ddply(framedDataFile, .(rt), summarize, tic = sum(intensity))
                        chromatogram$rt <- as.numeric(chromatogram$rt)
                        chromatogram$path_to_cdf_csv <- paths_to_cdf_csvs[file]

                    cat("   Appending chromatogram to list...\n")
                        chromatograms[[file]] <- chromatogram
                }

                chromatograms <- do.call(rbind, chromatograms)

                return( chromatograms )
            }

    #### integrationApp

        #' A Shiny app to integrate GC-FID and GC-MS data
        #'
        #' @param chromatograms A data frame containing columns: "rt", "tic", and "path_to_cdf_csv", which contain retention time, total ion chromatogram intensities, and paths to CDF.csv files generated by the convertCDFstoCSVs function.
        #' @param x_axis_start A numeric value for the lower x-axis bounds on the plot generated by the app. Defaults to full length.
        #' @param x_axis_end A numeric value for the upper x-axis bounds on the plot generated by the app. Defaults to full length.
        #' @param samples_monolist_path A path to a .csv file containing metadata for the samples you wish to analyze. Requied columns are: "rt_offset", "baseline_window", and "path_to_cdf_csv", which are for aligning chromatograms, adjusting baseline determination, and defining the path to the CDF.csv files for each sample, respectively.
        #' @param samples_monolist_subset Optional, a numeric vector (for example, "c(1:10)"), defining a subset of samples to be loaded.
        #' @param peaks_monolist_path A path to a .csv file containing metadata for all peaks in the sample set. Required columns are: peak_start", "peak_end", "path_to_cdf_csv", "area", "peak_number_within_sample", "rt_offset", "peak_start_rt_offset", "peak_end_rt_offset". This file is automatically generated by the app.
        #' @param zoom_and_scroll_rate Defines intervals of zooming and scrolling movement while running the app
        #' @examples
        #' @export
        #' integrationApp

        integrationApp <-   function(
                chromatograms,
                x_axis_start = NULL,
                x_axis_end = NULL,
                samples_monolist_path,
                create_new_samples_monolist = FALSE,
                samples_monolist_subset = NULL,
                peaks_monolist_path,
                create_new_peak_monolist = FALSE,
                zoom_and_scroll_rate
            ) {

                library(shiny)
                library(ggplot2)
                library(plyr)
                library(dplyr)
                library(DT)
                library(RColorBrewer)

                ## Messages for setup

	            	if( length(samples_monolist_path) == 0 ) {
	            		stop("You must specify a samples_monolist_path")
	            	}
	            	if( length(peaks_monolist_path) == 0 ) {
	            		stop("You must specify a peaks_monolist_path")
	            	}

                ## Set up new samples monolist

                    if ( create_new_samples_monolist == TRUE ) {
                        
                        samples_monolist <- data.frame(
                            Sample_ID = gsub("\\..*$", "", gsub(".*/", "", unique(chromatograms$path_to_cdf_csv))),
                            extraction_time = NA,
                            plant_age = 0,
                            relative_area = 0,
                            rt_offset = 0,
                            baseline_window = 400,
                            path_to_cdf_csv = unique(chromatograms$path_to_cdf_csv)
                        )

                        write.table(
                            x = samples_monolist,
                            file = samples_monolist_path,
                            row.names = FALSE,
                            sep = ","
                        )

                    }

                ## Set up several variables, plot_height, and x_axis limits if not specified in function call
                    
                    peak_data <- NULL
                    peak_points <- NULL
                    if ( length(samples_monolist_subset) > 0 ) {
                        plot_height <- 1000 + 100*length(samples_monolist_subset)    
                    } else {
                        plot_height <- 1000 + 100*dim(samples_monolist)[1]
                    }

                    if (length(x_axis_start) == 0) {
                        x_axis_start <<- min(chromatograms$rt)
                    }

                    if (length(x_axis_end) == 0) {
                        x_axis_end <<- max(chromatograms$rt)
                    }

                ## Set up new peak monolist
                    
                    if ( create_new_peak_monolist == TRUE ) {
                        
                        peak_data <- data.frame(
                          peak_start = 0,
                          peak_end = 0,
                          path_to_cdf_csv = "a",
                          area = 0
                        )

                        write.table(
                          x = peak_data[-1,],
                          file = peaks_monolist_path,
                          append = FALSE,
                          row.names = FALSE,
                          col.names = TRUE,
                          sep = ","
                        )

                    }

                ## Set up user interface

                    ui <- fluidPage(

                        tags$script('
                        $(document).on("keypress", function (e) {
                           Shiny.onInputChange("keypress", e.which);
                        });
                        '), 

                        verticalLayout(

                            plotOutput(
                                "massSpectra_z",
                                height = 150
                            ),

                            plotOutput(
                                "massSpectra_x",
                                height = 150
                            ),

                            plotOutput(
                              "chromatograms",
                              brush = "chromatogram_brush", 
                              click = "chromatogram_click", 
                              dblclick = "chromatogram_double_click",
                              height = plot_height
                            ),

                            verbatimTextOutput("key", placeholder = TRUE),

                            # DT::dataTableOutput("selected_peak"), ## Useful for troubleshooting

                            DT::dataTableOutput("peak_table")
                        )
                    )

                ## Set up the server

                    server <- function(input, output, session) {

                        ## Check keystoke value
                            output$key <- renderPrint({
                                input$keypress
                            })

                        ## Keys to move chromatogram view - zoom in and out, move L and R
                            observeEvent(input$keypress, {
                                if( input$keypress == 102 ) { x_axis_start <<- x_axis_start + zoom_and_scroll_rate } # Forward on "f"
                                if( input$keypress == 102 ) { x_axis_end <<- x_axis_end + zoom_and_scroll_rate } # Forward on "f"
                                if( input$keypress == 100 ) { x_axis_start <<- x_axis_start - zoom_and_scroll_rate } # Backward on "d"
                                if( input$keypress == 100 ) { x_axis_end <<- x_axis_end - zoom_and_scroll_rate } # Backward on "d"
                                if( input$keypress == 118 ) { x_axis_start <<- x_axis_start - zoom_and_scroll_rate } # Wider on "v"
                                if( input$keypress == 118 ) { x_axis_end <<- x_axis_end + zoom_and_scroll_rate } # Wider on "v"
                                if( input$keypress == 99 ) { x_axis_start <<- x_axis_start + zoom_and_scroll_rate } # Closer on "c"
                                if( input$keypress == 99 ) { x_axis_end <<- x_axis_end - zoom_and_scroll_rate } # Closer on "c"
                            })

                        ## Key to update chromatogram
                            
                            observeEvent(input$keypress, {      
                                
                                if( input$keypress == 113 ) { # Update on "q"
                                    
                                    output$chromatograms <- renderPlot({

                                        ## Read in samples monolist and put chromatograms into chromatograms_updated
                                            samples_monolist <- read.csv(samples_monolist_path)
                                            if ( length(samples_monolist_subset) > 0 ) {
                                                samples_monolist <- samples_monolist[samples_monolist_subset,]    
                                            }
                                            chromatograms_updated <- dplyr::filter(chromatograms, path_to_cdf_csv %in% samples_monolist$path_to_cdf_csv)

                                        ## Calculate baseline for each sample

                                            baselined_chromatograms <- list()

                                            for ( chrom in 1:length(unique(chromatograms_updated$path_to_cdf_csv)) ) {
                                      
                                                chromatogram <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == unique(chromatograms_updated$path_to_cdf_csv)[chrom])

                                                prelim_baseline_window <- samples_monolist$baseline_window[match(chromatogram$path_to_cdf_csv[1], samples_monolist$path_to_cdf_csv)]

                                                n_prelim_baseline_windows <- floor(length(chromatogram$rt)/prelim_baseline_window)
                                                prelim_baseline <- list()
                                                for ( i in 1:n_prelim_baseline_windows ) {
                                                  min <- min(chromatogram$tic[((prelim_baseline_window*(i-1))+1):(prelim_baseline_window*i)])
                                                  prelim_baseline[[i]] <-     data.frame(
                                                                          rt = chromatogram$rt[chromatogram$tic == min],
                                                                          min = min
                                                                      )
                                                }
                                                prelim_baseline <- do.call(rbind, prelim_baseline)
                                                chromatogram$in_prelim_baseline <- FALSE
                                                chromatogram$in_prelim_baseline[chromatogram$rt %in% prelim_baseline$rt] <- TRUE

                                                y = prelim_baseline$min
                                                x = prelim_baseline$rt

                                                baseline2 <- data.frame(
                                                              rt = chromatogram$rt,
                                                              y = approx(x, y, xout = chromatogram$rt)$y
                                                          )
                                                baseline2 <- baseline2[!is.na(baseline2$y),]
                                                chromatogram <- chromatogram[chromatogram$rt %in% baseline2$rt,]
                                                chromatogram$baseline <- baseline2$y

                                                baselined_chromatograms[[chrom]] <- data.frame(
                                                    rt = chromatogram$rt,
                                                    tic = chromatogram$tic,
                                                    path_to_cdf_csv = chromatogram$path_to_cdf_csv,
                                                    in_prelim_baseline = chromatogram$in_prelim_baseline,
                                                    baseline = chromatogram$baseline
                                                )

                                            }

                                            chromatograms_updated <- do.call(rbind, baselined_chromatograms)

                                        ## Add rt offset information for all chromatograms

                                            chromatograms_updated$rt_offset <- samples_monolist$rt_offset[match(chromatograms_updated$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                            chromatograms_updated$rt_rt_offset <- chromatograms_updated$rt + chromatograms_updated$rt_offset
                                            chromatograms_updated <<- chromatograms_updated

                                        ## Plot with peaks, if any
                                            
                                            peak_table <- read.csv(peaks_monolist_path)
                                    
                                            if (dim(peak_table)[1] > 0) {

                                                ## Filter out duplicate peaks and NA peaks
                                                    
                                                    peak_table <- plyr::ddply(peak_table, .(path_to_cdf_csv), mutate, duplicated = duplicated(peak_start))
                                                    peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
                                                    peak_table <- plyr::ddply(peak_table, .(path_to_cdf_csv), mutate, duplicated = duplicated(peak_end))
                                                    peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
                                                    peak_table <- peak_table[,!colnames(peak_table) == "duplicated"]
                                                    peak_table <- peak_table[!is.na(peak_table$peak_start),]

                                                ## Update with peak_number_within_sample
                                                    
                                                    peak_table_updated <- list()
                                                    for (sample_number in 1:length(unique(peak_table$path_to_cdf_csv))) {
                                                      peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == unique(peak_table$path_to_cdf_csv)[sample_number],]
                                                      peaks_in_this_sample <- peaks_in_this_sample[order(peaks_in_this_sample$peak_start),]
                                                      peaks_in_this_sample$peak_number_within_sample <- seq(1,length(peaks_in_this_sample$path_to_cdf_csv),1)
                                                      peak_table_updated[[sample_number]] <- peaks_in_this_sample
                                                    }
                                                    peak_table_updated <- do.call(rbind, peak_table_updated)
                                                    peak_table <- peak_table_updated

                                                ## Modify peaks with RT offset
                                                    
                                                    # for (sample_number in 1:length(unique(samples_monolist$path_to_cdf_csv))) {
                                                      
                                                    #   peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_monolist$path_to_cdf_csv[sample_number],]
                                                      
                                                    #   rt_offsets <- samples_monolist$rt_offset[match(peaks_in_this_sample$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                                    #   peak_start_rt_offsets <- peak_table$peak_start + peak_table$rt_offset
                                                    #   peak_end_rt_offsets <- peak_table$peak_end + peak_table$rt_offset
                                                      
                                                    #   peak_table$rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- rt_offsets
                                                    #   peak_table$peak_start_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_start_rt_offsets
                                                    #   peak_table$peak_end_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_end_rt_offsets

                                                    # }

                                                    peak_table$rt_offset <- samples_monolist$rt_offset[match(peak_table$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
                                                    peak_table$peak_start_rt_offset <- peak_table$peak_start + peak_table$rt_offset
                                                    peak_table$peak_end_rt_offset <- peak_table$peak_end + peak_table$rt_offset
                                                    peak_table$path_to_cdf_csv <- as.character(peak_table$path_to_cdf_csv)
                                        
                                                ## Update all peak areas in case baseline was adjusted
                                                    # peak_table_updated <- list()
                                                    for (sample_number in 1:length(unique(samples_monolist$path_to_cdf_csv))) {
                                                      
                                                      peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_monolist$path_to_cdf_csv[sample_number],]
                                                      
                                                      areas <- vector()
                                                      for (peak in 1:length(peaks_in_this_sample$peak_number_within_sample)) {
                                                        areas <- append(areas, 
                                                          sum(dplyr::filter(
                                                            chromatograms_updated[chromatograms_updated$path_to_cdf_csv == as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),], 
                                                            rt >= peaks_in_this_sample$peak_start[peak] & rt <= peaks_in_this_sample$peak_end[peak])$tic
                                                          ) - 
                                                          sum(dplyr::filter(
                                                            chromatograms_updated[chromatograms_updated$path_to_cdf_csv == as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),], 
                                                            rt >= peaks_in_this_sample$peak_start[peak] & rt <= peaks_in_this_sample$peak_end[peak])$baseline
                                                          )
                                                        )
                                                      }

                                                      peak_table$area[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- areas

                                                      # peaks_in_this_sample$area <- areas
                                                      # peak_table_updated[[sample_number]] <- peaks_in_this_sample
                                                    }
                                                    # peak_table_updated <- do.call(rbind, peak_table_updated)
                                                    # peak_table <- peak_table_updated
                                      
                                                ## Write out peaks now with assigned peak_number_within_sample and RT offset
                                        
                                                    write.table(peak_table, file = peaks_monolist_path, col.names = TRUE, sep = ",", row.names = FALSE)


                                                ## Create the plot with subsetted data to make it faster

                                                    ## Subset the chromatograms and peaks

                                                        chromatograms_updated <- dplyr::filter(chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end)
                                                        peak_table <- dplyr::filter(peak_table, peak_start_rt_offset > x_axis_start & peak_end_rt_offset < x_axis_end)

                                                    ## Make chromatogram plot object

                                                        # Make their labels easy to read
                                                            facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated$path_to_cdf_csv))
                                                            names(facet_labels) <- chromatograms_updated$path_to_cdf_csv

                                                            p <-    ggplot() + 
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = tic)) +
                                                                scale_x_continuous(limits = c(x_axis_start, x_axis_end)) +
                                                                facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # facet_grid(path_to_cdf_csv~., labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # scale_y_continuous(limits = c(0, max(dplyr::filter(chromatograms_updated, rt > x_axis_start & rt < x_axis_end)$tic))) +
                                                                theme_classic() +
                                                                scale_fill_continuous(type = "viridis") +
                                                                theme(
                                                                  legend.position = 'none',
                                                                  legend.title = element_blank()
                                                                )

                                                    ## Add peaks

                                                        for (peak in 1:dim(peak_table)[1]) {
                                                            
                                                            signal_for_this_peak <- dplyr::filter(
                                                                chromatograms_updated[chromatograms_updated$path_to_cdf_csv == peak_table[peak,]$path_to_cdf_csv,], 
                                                                rt_rt_offset > peak_table[peak,]$peak_start_rt_offset, 
                                                                rt_rt_offset < peak_table[peak,]$peak_end_rt_offset
                                                            )

                                                            if (dim(signal_for_this_peak)[1] > 0) {

                                                                signal_for_this_peak$peak_number_within_sample <- peak_table$peak_number_within_sample[peak]
                                                                
                                                                p <- p +  geom_vline(data = signal_for_this_peak[1,], mapping = aes(xintercept = rt_rt_offset), alpha = 0.3) +
                                                                          geom_ribbon(data = signal_for_this_peak, mapping = aes(x = rt_rt_offset, ymax = tic, ymin = baseline, fill = peak_number_within_sample)) +
                                                                          geom_text(data = signal_for_this_peak, mapping = aes(label = peak_number_within_sample, x = median(rt_rt_offset), y = max(tic)))
                                                            }
                                                        }

                                            } else {

                                                ## Make chromatogram plot object

                                                    chromatograms_updated <- dplyr::filter(chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end)

                                                        # Make their labels easy to read
                                                            facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated$path_to_cdf_csv))
                                                            names(facet_labels) <- chromatograms_updated$path_to_cdf_csv

                                                        p <-  ggplot() + 
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
                                                                geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = tic)) +
                                                                scale_x_continuous(limits = c(x_axis_start, x_axis_end)) +
                                                                facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # facet_grid(path_to_cdf_csv~., labeller = labeller(path_to_cdf_csv = facet_labels)) +
                                                                # scale_y_continuous(limits = c(0, max(dplyr::filter(chromatograms_updated, rt > x_axis_start & rt < x_axis_end)$tic))) +
                                                                theme_classic() +
                                                                scale_fill_continuous(type = "viridis") +
                                                                theme(
                                                                  legend.position = 'none',
                                                                  legend.title = element_blank()
                                                                )
                                            }

                                        ## Draw the plot
                                            
                                            p
                                    })
                                }
                            })

                        ## Transfer chromatogram_brush info to selected_peak table
                            output$selected_peak <- DT::renderDataTable(DT::datatable({

                                if ( !is.null(input$chromatogram_brush )) {
                                    peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
                                    peak_data <-  data.frame(
                                        peak_start = min(peak_points$rt),
                                        peak_end = max(peak_points$rt),
                                        path_to_cdf_csv = peak_points$path_to_cdf_csv[1],
                                        area = sum(peak_points$tic)
                                    )
                                    peak_data
                                } else {
                                    NULL
                                }

                            }))

                        ## Remove selected peaks
                            observeEvent(input$keypress, {      
                                if( input$keypress == 114 ) { # Update on "r"

                                    if ( !is.null(input$chromatogram_brush )) {

                                        peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
                                        selection_start = min(peak_points$rt)
                                        selection_end = max(peak_points$rt)
                                        path_to_cdf_csv = peak_points$path_to_cdf_csv[1]
                                        peak_table <- read.csv(peaks_monolist_path)

                                        peak_table <- peak_table[
                                            !apply(cbind(
                                                peak_table$peak_start > selection_start,
                                                peak_table$peak_end < selection_end,
                                                peak_table$path_to_cdf_csv == as.character(peak_points$path_to_cdf_csv[1])
                                            ), 1, all)
                                        ,]

                                        write.table(
                                            x = peak_table,
                                            file = peaks_monolist_path,
                                            append = FALSE,
                                            row.names = FALSE,
                                            col.names = TRUE,
                                            sep = ","
                                        )
                                    }
                                }
                            })

                        ## Append single peak with "a" keystroke 
                            observeEvent(input$keypress, {

                                # Do nothing if no selection
                                    if(is.null(input$chromatogram_brush)) {
                                        return()
                                    }

                                # If selection and "a" is pressed, add the selection to the peak table
                                    if( input$keypress == 97 ) {
                                    
                                        write.table(
                                            x = data.frame(
                                                    peak_start = min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                    peak_end = max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                    path_to_cdf_csv = brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1],
                                                    area = sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$tic) - sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$baseline)
                                                ),
                                            file = peaks_monolist_path,
                                            append = TRUE,
                                            row.names = FALSE,
                                            col.names = FALSE,
                                            sep = ","
                                        )

                                        output$peak_table <- DT::renderDataTable(DT::datatable({
                                            peak_table <- read.csv(peaks_monolist_path)
                                            peak_table
                                        }))
                                    }
                            })

                        ## Global append peak with "g" keystroke 
                            observeEvent(input$keypress, {

                                # Do nothing if no selection
                                    if(is.null(input$chromatogram_brush)) {
                                        return()
                                    }

                                # If selection and "g" is pressed, add the selection to the peak table
                                    if( input$keypress == 103 ) {
                                    
                                        x_peaks <-  data.frame(
                                                        peak_start = min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt_rt_offset),
                                                        peak_end = max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt_rt_offset),
                                                        path_to_cdf_csv = unique(chromatograms_updated$path_to_cdf_csv),
                                                        area = sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$tic) - sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$baseline)
                                                    )

                                        x_peaks$peak_start <- x_peaks$peak_start - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]
                                        x_peaks$peak_end <- x_peaks$peak_end - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]

                                        write.table(
                                            x = x_peaks,
                                            file = peaks_monolist_path,
                                            append = TRUE,
                                            row.names = FALSE,
                                            col.names = FALSE,
                                            sep = ","
                                        )

                                        output$peak_table <- DT::renderDataTable(DT::datatable({
                                            peak_table <- read.csv(peaks_monolist_path)
                                            peak_table
                                        }))
                                    }
                            })

                        ## Show mass spectrum on "1" keypress
                            observeEvent(input$keypress, {

                                # Do nothing if no selection
                                    if(is.null(input$chromatogram_brush)) {
                                        return()
                                    }

                                # If selection and "1" is pressed, extract and print mass spectra
                                    if( input$keypress == 49 ) {

                                        output$massSpectra_z <- renderPlot({
                                            framedDataFile <- isolate(as.data.frame(
                                                                data.table::fread(as.character(
                                                                    brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                                                                ))
                                            ))
                                            framedDataFile <- isolate(dplyr::filter(
                                                                    framedDataFile, 
                                                                    rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                                    rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                                                ))
                                            framedDataFile <- plyr::ddply(framedDataFile, .(mz), summarize, intensity = sum(intensity))
                                            framedDataFile_1 <<- framedDataFile
                                            framedDataFile$intensity <- framedDataFile$intensity*100/max(framedDataFile$intensity)
                                            ggplot() + 
                                                geom_bar(data = framedDataFile, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
                                                theme_classic() +
                                                scale_x_continuous(expand = c(0,0)) +
                                                scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
                                                geom_text(data = dplyr::filter(framedDataFile, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
                                        })
                                    }
                            })

                        ## Show subtracted mass spectrum for "1" on "3" keypress
                            observeEvent(input$keypress, {

                                # Do nothing if no selection
                                    if(is.null(input$chromatogram_brush)) {
                                        return()
                                    }

                                # If selection and "51" is pressed, extract and print mass spectra
                                    if( input$keypress == 51 ) {

                                        output$massSpectra_z <- renderPlot({
                                            framedDataFile_to_subtract <- isolate(as.data.frame(
                                                                data.table::fread(as.character(
                                                                    brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                                                                ))
                                            ))
                                            framedDataFile_to_subtract <- isolate(dplyr::filter(
                                                                    framedDataFile_to_subtract, 
                                                                    rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                                    rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                                                ))
                                            framedDataFile_to_subtract <- plyr::ddply(framedDataFile_to_subtract, .(mz), summarize, intensity = sum(intensity))
                                            framedDataFile_1$intensity <- framedDataFile_1$intensity - framedDataFile_to_subtract$intensity
                                            framedDataFile_1 <<- framedDataFile_1
                                            framedDataFile_1$intensity <- framedDataFile_1$intensity*100/max(framedDataFile_1$intensity)
                                            ggplot() + 
                                                geom_bar(data = framedDataFile_1, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
                                                theme_classic() +
                                                scale_x_continuous(expand = c(0,0)) +
                                                scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
                                                geom_text(data = dplyr::filter(framedDataFile_1, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
                                        })
                                    }
                            })

                        ## Save mass spectrum in "1" on "5" keypress
                            observeEvent(input$keypress, {

                                # Do nothing if no selection
                                    if(is.null(input$chromatogram_brush)) {
                                        return()
                                    }

                                # If selection and "51" is pressed, extract and print mass spectra
                                    if( input$keypress == 53 ) {
                                        framedDataFile_1$intensity[framedDataFile_1$intensity < 0] <- 0
                                        write.csv(framedDataFile_1, "integration_app_spectrum_1.csv", row.names = FALSE)
                                    }
                            })

                        ## Show mass spectrum on "2" keypress
                            observeEvent(input$keypress, {

                                # Do nothing if no selection
                                    if(is.null(input$chromatogram_brush)) {
                                        return()
                                    }

                                # If selection and "2" is pressed, extract and print mass spectra
                                    if( input$keypress == 50 ) {

                                        output$massSpectra_x <- renderPlot({
                                            framedDataFile <- isolate(as.data.frame(
                                                                data.table::fread(as.character(
                                                                    brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                                                                ))
                                        ))
                                        framedDataFile <- isolate(dplyr::filter(
                                                                    framedDataFile, 
                                                                    rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                                    rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                                                ))
                                        framedDataFile <- plyr::ddply(framedDataFile, .(mz), summarize, intensity = sum(intensity))
                                        framedDataFile_2 <<- framedDataFile
                                        framedDataFile$intensity <- framedDataFile$intensity*100/max(framedDataFile$intensity)
                                        ggplot() + 
                                            geom_bar(data = framedDataFile, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
                                            theme_classic() +
                                            scale_x_continuous(expand = c(0,0)) +
                                            scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
                                            geom_text(data = dplyr::filter(framedDataFile, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
                                        })
                                    }
                            })

                            ## Show subtracted mass spectrum for "2" on "4" keypress
                                observeEvent(input$keypress, {

                                    # Do nothing if no selection
                                        if(is.null(input$chromatogram_brush)) {
                                            return()
                                        }

                                    # If selection and "51" is pressed, extract and print mass spectra
                                        if( input$keypress == 52 ) {

                                            output$massSpectra_x <- renderPlot({
                                                framedDataFile_to_subtract <- isolate(as.data.frame(
                                                                    data.table::fread(as.character(
                                                                        brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
                                                                    ))
                                                ))
                                                framedDataFile_to_subtract <- isolate(dplyr::filter(
                                                                        framedDataFile_to_subtract, 
                                                                        rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
                                                                        rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
                                                                    ))
                                                framedDataFile_to_subtract <- plyr::ddply(framedDataFile_to_subtract, .(mz), summarize, intensity = sum(intensity))
                                                framedDataFile_2$intensity <- framedDataFile_2$intensity - framedDataFile_to_subtract$intensity
                                                framedDataFile_2 <<- framedDataFile_2
                                                framedDataFile_2$intensity <- framedDataFile_2$intensity*100/max(framedDataFile_2$intensity)
                                                ggplot() + 
                                                    geom_bar(data = framedDataFile_2, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
                                                    theme_classic() +
                                                    scale_x_continuous(expand = c(0,0)) +
                                                    scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
                                                    geom_text(data = dplyr::filter(framedDataFile_2, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
                                            })
                                        }
                                })

                            ## Save mass spectrum in "2" on "6" keypress
                                observeEvent(input$keypress, {

                                    # Do nothing if no selection
                                        if(is.null(input$chromatogram_brush)) {
                                            return()
                                        }

                                    # If selection and "51" is pressed, extract and print mass spectra
                                        if( input$keypress == 54 ) {
                                            framedDataFile_2$intensity[framedDataFile_2$intensity < 0] <- 0
                                            write.csv(framedDataFile_2, "integration_app_spectrum_2.csv")
                                        }
                                })
                    }

                ## Call the app
                    
                    shinyApp(ui = ui, server = server)

            }

    #### mergePeakLists

	    #' Merge multiple peaklists
	    #'
	    #' @param analysis_directory_path Path to the directory above all the peak lists
	    #' @examples
	    #' @export
	    #' mergePeakLists

	    mergePeakLists <- function(analysis_directory_path) {

	        analysis_directory_path <- "/Users/lucasbusta/Desktop/test"

	        paths_to_peak_monolists <- paste(
	            analysis_directory_path,
	            dir(analysis_directory_path),
	            "peaks_monolist.csv",
	            sep = "/"
	        )

	        peak_list <- list()
	        for (i in 1:length(paths_to_peak_monolists)) {
	            peak_list[[i]] <- readMonolist(paths_to_peak_monolists[i])
	        }

	        do.call(rbind, peak_list)

	    }

	#### integrationAppLite

	    #' A Shiny app to integrate GC-FID and GC-MS data
	    #'
	    #' @param chromatograms A data frame containing columns: "rt", "tic", and "path_to_cdf_csv", which contain retention time, total ion chromatogram intensities, and paths to CDF.csv files generated by the convertCDFstoCSVs function.
	    #' @param x_axis_start A numeric value for the lower x-axis bounds on the plot generated by the app. Defaults to full length.
	    #' @param x_axis_end A numeric value for the upper x-axis bounds on the plot generated by the app. Defaults to full length.
	    #' @param samples_monolist_path A path to a .csv file containing metadata for the samples you wish to analyze. Requied columns are: "rt_offset", "baseline_window", and "path_to_cdf_csv", which are for aligning chromatograms, adjusting baseline determination, and defining the path to the CDF.csv files for each sample, respectively.
	    #' @param samples_monolist_subset Optional, a numeric vector (for example, "c(1:10)"), defining a subset of samples to be loaded.
	    #' @param peaks_monolist_path A path to a .csv file containing metadata for all peaks in the sample set. Required columns are: peak_start", "peak_end", "path_to_cdf_csv", "area", "peak_number_within_sample", "rt_offset", "peak_start_rt_offset", "peak_end_rt_offset". This file is automatically generated by the app.
	    #' @param zoom_and_scroll_rate Defines intervals of zooming and scrolling movement while running the app
	    #' @examples
	    #' @export
	    #' integrationAppLite

	    integrationAppLite <-   function(
	            CDF_directory_path,
	            create_new_samples_monolist = FALSE,
	            create_new_peak_monolist = FALSE,
	            zoom_and_scroll_rate = 100,
	            baseline_window = 400
	        ) {

	            library(shiny)
	            library(tidyverse)
	            library(DT)
	            library(RColorBrewer)

	            ## Set up monolists

	            	samples_monolist_subset = NULL
	            	samples_monolist_path <- paste0(CDF_directory_path, "/samples_monolist.csv")
	        		peaks_monolist_path <- paste0(CDF_directory_path, "/peaks_monolist.csv")
	            	
	            ## Covert CDF to csv, if necessary

	                paths_to_cdfs <- paste(
	                    CDF_directory_path,
	                    dir(CDF_directory_path)[grep("*.CDF$", dir(CDF_directory_path))],
	                    sep = "/"
	                )

	                convertCDFstoCSVs(paths_to_cdfs)

	            ## Recreate chromatograms.csv, if necessary

	                paths_to_cdf_csvs <- paste0(paths_to_cdfs, ".csv")

	                # If chromatograms doesn't exist, make it 
	                    if (!file.exists(paste(CDF_directory_path, "chromatograms.csv", sep = "/"))) {
	                        chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
	                        writeMonolist(chromatograms, paste(CDF_directory_path, "chromatograms.csv", sep = "/"))
	                    } else {
	                        # If it exists, check to see if all cdfs in this folder are in it, if not, recreate
	                        chromatograms <- readMonolist(paste(CDF_directory_path, "chromatograms.csv", sep = "/"))
	                        if (!all(paths_to_cdf_csvs %in% unique(chromatograms$path_to_cdf_csv))) {
	                            chromatograms <- extractChromatogramsFromCSVs(paths_to_cdf_csvs)
	                            writeMonolist(chromatograms, paste(CDF_directory_path, "chromatograms.csv", sep = "/"))
	                        }
	                    }

	            ## Set up new samples monolist

	                if ( create_new_samples_monolist == TRUE ) {
	                    
	                    samples_monolist <- data.frame(
	                        Sample_ID = gsub("\\..*$", "", gsub(".*/", "", unique(chromatograms$path_to_cdf_csv))),
	                        rt_offset = 0,
	                        baseline_window = baseline_window,
	                        path_to_cdf_csv = unique(chromatograms$path_to_cdf_csv)
	                    )

	                    write.table(
	                        x = samples_monolist,
	                        file = samples_monolist_path,
	                        row.names = FALSE,
	                        sep = ","
	                    )

	                }

	            ## Set up several variables, plot_height, and x_axis limits if not specified in function call
	                
	                peak_data <- NULL
	                peak_points <- NULL
	                if ( length(samples_monolist_subset) > 0 ) {
	                    plot_height <- 200 + 100*length(samples_monolist_subset)    
	                } else {
	                    plot_height <- 200 + 100*dim(samples_monolist)[1]
	                }

	                x_axis_start <<- min(chromatograms$rt)
	                x_axis_end <<- max(chromatograms$rt)
	                
	            ## Set up new peak monolist
	                
	                if ( create_new_peak_monolist == TRUE ) {
	                    
	                    peak_data <- data.frame(
	                      peak_start = 0,
	                      peak_end = 0,
	                      peak_ID = "unknown",
	                      path_to_cdf_csv = "a",
	                      area = 0
	                    )

	                    write.table(
	                      x = peak_data[-1,],
	                      file = peaks_monolist_path,
	                      append = FALSE,
	                      row.names = FALSE,
	                      col.names = TRUE,
	                      sep = ","
	                    )

	                }

	            ## Set up user interface

	                ui <- fluidPage(

	                    tags$script('
	                    $(document).on("keypress", function (e) {
	                       Shiny.onInputChange("keypress", e.which);
	                    });
	                    '), 

	                    verticalLayout(

	                        plotOutput(
	                            "massSpectra_z",
	                            height = 150
	                        ),

	                        plotOutput(
	                            "massSpectra_x",
	                            height = 150
	                        ),

	                        plotOutput(
	                          "chromatograms",
	                          brush = "chromatogram_brush", 
	                          click = "chromatogram_click", 
	                          dblclick = "chromatogram_double_click",
	                          height = plot_height
	                        ),

	                        verbatimTextOutput("key", placeholder = TRUE),

	                        # DT::dataTableOutput("selected_peak"), ## Useful for troubleshooting

	                        rhandsontable::rHandsontableOutput("peak_table")
	                    )
	                )

	            ## Set up the server

	                server <- function(input, output, session) {

	                    ## Check keystoke value
	                        output$key <- renderPrint({
	                            input$keypress
	                        })

	                    ## Keys to move chromatogram view - zoom in and out, move L and R
	                        observeEvent(input$keypress, {
	                            if( input$keypress == 70 ) { x_axis_start <<- x_axis_start + zoom_and_scroll_rate } # Forward on "F"
	                            if( input$keypress == 70 ) { x_axis_end <<- x_axis_end + zoom_and_scroll_rate } # Forward on "F"
	                            if( input$keypress == 68 ) { x_axis_start <<- x_axis_start - zoom_and_scroll_rate } # Backward on "D"
	                            if( input$keypress == 68 ) { x_axis_end <<- x_axis_end - zoom_and_scroll_rate } # Backward on "D"
	                            if( input$keypress == 86 ) { x_axis_start <<- x_axis_start - zoom_and_scroll_rate } # Wider on "V"
	                            if( input$keypress == 86 ) { x_axis_end <<- x_axis_end + zoom_and_scroll_rate } # Wider on "V"
	                            if( input$keypress == 67 ) { x_axis_start <<- x_axis_start + zoom_and_scroll_rate } # Closer on "C"
	                            if( input$keypress == 67 ) { x_axis_end <<- x_axis_end - zoom_and_scroll_rate } # Closer on "C"
	                        })

	                    ## Save manual changes to table on "Z" (90) keystroke

	                        observeEvent(input$keypress, {

	                            if (input$keypress == 90 ) {

	                                ## Write out any modifications to peak table (i.e. sample IDs)
	                                    
	                                    hot = isolate(input$peak_table)
	                                    if (!is.null(hot)) {
	                                        writeMonolist(hot_to_r(input$peak_table), peaks_monolist_path)
	                                        print(peaks_monolist_path)
	                                    }

	                            }

	                        })

	                    ## Update chromatogram on "Q" (81) keystroke
	                        
	                        observeEvent(input$keypress, {      
	                            
	                            if( input$keypress == 81 ) { # Update on "Q"
	                                
	                                output$chromatograms <- renderPlot({

	                                    ## Read in samples monolist and put chromatograms into chromatograms_updated
	                                        
	                                        samples_monolist <- read.csv(samples_monolist_path)
	                                        if ( length(samples_monolist_subset) > 0 ) {
	                                            samples_monolist <- samples_monolist[samples_monolist_subset,]    
	                                        }
	                                        chromatograms_updated <- dplyr::filter(chromatograms, path_to_cdf_csv %in% samples_monolist$path_to_cdf_csv)

	                                    ## Calculate baseline for each sample

	                                        baselined_chromatograms <- list()

	                                        for ( chrom in 1:length(unique(chromatograms_updated$path_to_cdf_csv)) ) {
	                                  
	                                            chromatogram <- dplyr::filter(chromatograms_updated, path_to_cdf_csv == unique(chromatograms_updated$path_to_cdf_csv)[chrom])

	                                            prelim_baseline_window <- samples_monolist$baseline_window[match(chromatogram$path_to_cdf_csv[1], samples_monolist$path_to_cdf_csv)]

	                                            n_prelim_baseline_windows <- floor(length(chromatogram$rt)/prelim_baseline_window)
	                                            prelim_baseline <- list()
	                                            for ( i in 1:n_prelim_baseline_windows ) {
	                                              min <- min(chromatogram$tic[((prelim_baseline_window*(i-1))+1):(prelim_baseline_window*i)])
	                                              prelim_baseline[[i]] <-     data.frame(
	                                                                      rt = chromatogram$rt[chromatogram$tic == min],
	                                                                      min = min
	                                                                  )
	                                            }
	                                            prelim_baseline <- do.call(rbind, prelim_baseline)
	                                            chromatogram$in_prelim_baseline <- FALSE
	                                            chromatogram$in_prelim_baseline[chromatogram$rt %in% prelim_baseline$rt] <- TRUE

	                                            y = prelim_baseline$min
	                                            x = prelim_baseline$rt

	                                            baseline2 <- data.frame(
	                                                          rt = chromatogram$rt,
	                                                          y = approx(x, y, xout = chromatogram$rt)$y
	                                                      )
	                                            baseline2 <- baseline2[!is.na(baseline2$y),]
	                                            chromatogram <- chromatogram[chromatogram$rt %in% baseline2$rt,]
	                                            chromatogram$baseline <- baseline2$y

	                                            baselined_chromatograms[[chrom]] <- data.frame(
	                                                rt = chromatogram$rt,
	                                                tic = chromatogram$tic,
	                                                path_to_cdf_csv = chromatogram$path_to_cdf_csv,
	                                                in_prelim_baseline = chromatogram$in_prelim_baseline,
	                                                baseline = chromatogram$baseline
	                                            )

	                                        }

	                                        chromatograms_updated <- do.call(rbind, baselined_chromatograms)

	                                    ## Add rt offset information for all chromatograms

	                                        chromatograms_updated$rt_offset <- samples_monolist$rt_offset[match(chromatograms_updated$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
	                                        chromatograms_updated$rt_rt_offset <- chromatograms_updated$rt + chromatograms_updated$rt_offset
	                                        chromatograms_updated <<- chromatograms_updated

	                                    ## Plot with peaks, if any
	                                        
	                                        peak_table <- read.csv(peaks_monolist_path)
	                                
	                                        if (dim(peak_table)[1] > 0) {

	                                            ## Filter out duplicate peaks and NA peaks
	                                                
	                                                peak_table <- plyr::ddply(peak_table, .(path_to_cdf_csv), mutate, duplicated = duplicated(peak_start))
	                                                peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
	                                                peak_table <- plyr::ddply(peak_table, .(path_to_cdf_csv), mutate, duplicated = duplicated(peak_end))
	                                                peak_table <- dplyr::filter(peak_table, duplicated == FALSE)
	                                                peak_table <- peak_table[,!colnames(peak_table) == "duplicated"]
	                                                peak_table <- peak_table[!is.na(peak_table$peak_start),]

	                                            ## Update with peak_number_within_sample
	                                                
	                                                peak_table_updated <- list()
	                                                for (sample_number in 1:length(unique(peak_table$path_to_cdf_csv))) {
	                                                  peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == unique(peak_table$path_to_cdf_csv)[sample_number],]
	                                                  peaks_in_this_sample <- peaks_in_this_sample[order(peaks_in_this_sample$peak_start),]
	                                                  peaks_in_this_sample$peak_number_within_sample <- seq(1,length(peaks_in_this_sample$path_to_cdf_csv),1)
	                                                  peak_table_updated[[sample_number]] <- peaks_in_this_sample
	                                                }
	                                                peak_table_updated <- do.call(rbind, peak_table_updated)
	                                                peak_table <- peak_table_updated

	                                            ## Modify peaks with RT offset
	                                                
	                                                # for (sample_number in 1:length(unique(samples_monolist$path_to_cdf_csv))) {
	                                                  
	                                                #   peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_monolist$path_to_cdf_csv[sample_number],]
	                                                  
	                                                #   rt_offsets <- samples_monolist$rt_offset[match(peaks_in_this_sample$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
	                                                #   peak_start_rt_offsets <- peak_table$peak_start + peak_table$rt_offset
	                                                #   peak_end_rt_offsets <- peak_table$peak_end + peak_table$rt_offset
	                                                  
	                                                #   peak_table$rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- rt_offsets
	                                                #   peak_table$peak_start_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_start_rt_offsets
	                                                #   peak_table$peak_end_rt_offset[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- peak_end_rt_offsets

	                                                # }

	                                                peak_table$rt_offset <- samples_monolist$rt_offset[match(peak_table$path_to_cdf_csv, samples_monolist$path_to_cdf_csv)]
	                                                peak_table$peak_start_rt_offset <- peak_table$peak_start + peak_table$rt_offset
	                                                peak_table$peak_end_rt_offset <- peak_table$peak_end + peak_table$rt_offset
	                                                peak_table$path_to_cdf_csv <- as.character(peak_table$path_to_cdf_csv)
	                                    
	                                            ## Update all peak areas in case baseline was adjusted
	                                                
	                                                for (sample_number in 1:length(unique(samples_monolist$path_to_cdf_csv))) {
	                                                  
	                                                  peaks_in_this_sample <- peak_table[peak_table$path_to_cdf_csv == samples_monolist$path_to_cdf_csv[sample_number],]
	                                                  
	                                                  areas <- vector()
	                                                  for (peak in 1:length(peaks_in_this_sample$peak_number_within_sample)) {
	                                                    areas <- append(areas, 
	                                                      sum(dplyr::filter(
	                                                        chromatograms_updated[chromatograms_updated$path_to_cdf_csv == as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),], 
	                                                        rt >= peaks_in_this_sample$peak_start[peak] & rt <= peaks_in_this_sample$peak_end[peak])$tic
	                                                      ) - 
	                                                      sum(dplyr::filter(
	                                                        chromatograms_updated[chromatograms_updated$path_to_cdf_csv == as.character(peaks_in_this_sample$path_to_cdf_csv[peak]),], 
	                                                        rt >= peaks_in_this_sample$peak_start[peak] & rt <= peaks_in_this_sample$peak_end[peak])$baseline
	                                                      )
	                                                    )
	                                                  }

	                                                  peak_table$area[peak_table$path_to_cdf_csv == as.character(samples_monolist$path_to_cdf_csv[sample_number])] <- areas

	                                                  # peaks_in_this_sample$area <- areas
	                                                  # peak_table_updated[[sample_number]] <- peaks_in_this_sample
	                                                }
	                                                # peak_table_updated <- do.call(rbind, peak_table_updated)
	                                                # peak_table <- peak_table_updated
	                                  
	                                            ## Write out peaks now with assigned peak_number_within_sample and RT offset, update the peak_table in ui
	                                    
	                                                write.table(peak_table, file = peaks_monolist_path, col.names = TRUE, sep = ",", row.names = FALSE)

	                                                output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
				                                        peak_table2 <- read.csv(peaks_monolist_path)
				                                        peak_table2
				                                    }))

	                                            ## Create the plot with subsetted data to make it faster

	                                                ## Subset the chromatograms and peaks

	                                                    chromatograms_updated <- dplyr::filter(chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end)
	                                                    peak_table <- dplyr::filter(peak_table, peak_start_rt_offset > x_axis_start & peak_end_rt_offset < x_axis_end)

	                                                ## Make chromatogram plot object

	                                                    # Make their labels easy to read
	                                                        facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated$path_to_cdf_csv))
	                                                        names(facet_labels) <- chromatograms_updated$path_to_cdf_csv

	                                                        p <-    ggplot() + 
	                                                            geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
	                                                            geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = tic)) +
	                                                            scale_x_continuous(limits = c(x_axis_start, x_axis_end)) +
	                                                            facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
	                                                            # facet_grid(path_to_cdf_csv~., labeller = labeller(path_to_cdf_csv = facet_labels)) +
	                                                            # scale_y_continuous(limits = c(0, max(dplyr::filter(chromatograms_updated, rt > x_axis_start & rt < x_axis_end)$tic))) +
	                                                            theme_classic() +
	                                                            scale_fill_continuous(type = "viridis") +
	                                                            theme(
	                                                              legend.position = 'none',
	                                                              legend.title = element_blank()
	                                                            )

	                                                ## Add peaks

	                                                    for (peak in 1:dim(peak_table)[1]) {
	                                                        
	                                                        signal_for_this_peak <- dplyr::filter(
	                                                            chromatograms_updated[chromatograms_updated$path_to_cdf_csv == peak_table[peak,]$path_to_cdf_csv,], 
	                                                            rt_rt_offset > peak_table[peak,]$peak_start_rt_offset, 
	                                                            rt_rt_offset < peak_table[peak,]$peak_end_rt_offset
	                                                        )

	                                                        if (dim(signal_for_this_peak)[1] > 0) {

	                                                            signal_for_this_peak$peak_number_within_sample <- peak_table$peak_number_within_sample[peak]
	                                                            
	                                                            p <- p +  geom_vline(data = signal_for_this_peak[1,], mapping = aes(xintercept = rt_rt_offset), alpha = 0.3) +
	                                                                      geom_ribbon(data = signal_for_this_peak, mapping = aes(x = rt_rt_offset, ymax = tic, ymin = baseline, fill = peak_number_within_sample)) +
	                                                                      geom_text(data = signal_for_this_peak, mapping = aes(label = peak_number_within_sample, x = median(rt_rt_offset), y = max(tic)))
	                                                        }
	                                                    }

	                                        } else {

	                                            ## Make chromatogram plot object

	                                                chromatograms_updated <- dplyr::filter(chromatograms_updated, rt_rt_offset > x_axis_start & rt_rt_offset < x_axis_end)
	                                        	
	                                                    # Make their labels easy to read
	                                                        facet_labels <- gsub(".CDF.csv", "", gsub(".*/", "", chromatograms_updated$path_to_cdf_csv))
	                                                        names(facet_labels) <- chromatograms_updated$path_to_cdf_csv

	                                                    p <-  ggplot() + 
	                                                            geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = baseline), color = "grey") +
	                                                            geom_line(data = chromatograms_updated, mapping = aes(x = rt_rt_offset, y = tic)) +
	                                                            scale_x_continuous(limits = c(x_axis_start, x_axis_end)) +
	                                                            facet_grid(path_to_cdf_csv~., scales = "free_y", labeller = labeller(path_to_cdf_csv = facet_labels)) +
	                                                            # facet_grid(path_to_cdf_csv~., labeller = labeller(path_to_cdf_csv = facet_labels)) +
	                                                            # scale_y_continuous(limits = c(0, max(dplyr::filter(chromatograms_updated, rt > x_axis_start & rt < x_axis_end)$tic))) +
	                                                            theme_classic() +
	                                                            scale_fill_continuous(type = "viridis") +
	                                                            theme(
	                                                              legend.position = 'none',
	                                                              legend.title = element_blank()
	                                                            )
	                                        }

	                                    ## Draw the plot
	                                        
	                                        p
	                                })
	                            }
	                        })

	                    ## Transfer chromatogram_brush info to selected_peak table
	                        
	                        output$selected_peak <- DT::renderDataTable(DT::datatable({

	                            if ( !is.null(input$chromatogram_brush )) {
	                                peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
	                                peak_data <-  data.frame(
	                                    peak_start = min(peak_points$rt),
	                                    peak_end = max(peak_points$rt),
	                                    peak_ID = "unknown",
	                                    path_to_cdf_csv = peak_points$path_to_cdf_csv[1],
	                                    area = sum(peak_points$tic)
	                                )
	                                peak_data
	                            } else {
	                                NULL
	                            }

	                        }))

	                    ## Remove selected peaks with "R" (82) keystroke
	                        
	                        observeEvent(input$keypress, {      
	                            if( input$keypress == 82 ) { # Update on "R"

	                                if ( !is.null(input$chromatogram_brush )) {

	                                    peak_points <- brushedPoints(chromatograms_updated, input$chromatogram_brush)
	                                    selection_start = min(peak_points$rt)
	                                    selection_end = max(peak_points$rt)
	                                    path_to_cdf_csv = peak_points$path_to_cdf_csv[1]
	                                    peak_table <- read.csv(peaks_monolist_path)

	                                    peak_table <- peak_table[
	                                        !apply(cbind(
	                                            peak_table$peak_start > selection_start,
	                                            peak_table$peak_end < selection_end,
	                                            peak_table$path_to_cdf_csv == as.character(peak_points$path_to_cdf_csv[1])
	                                        ), 1, all)
	                                    ,]

	                                    write.table(
	                                        x = peak_table,
	                                        file = peaks_monolist_path,
	                                        append = FALSE,
	                                        row.names = FALSE,
	                                        col.names = TRUE,
	                                        sep = ","
	                                    )
	                                }
	                            }
	                        })

	                    ## Append single peak with "A" (65) keystroke 
	                        
	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                            # If selection and "a" is pressed, add the selection to the peak table
	                                if( input$keypress == 65 ) {
	                                
	                                    write.table(
	                                        x = data.frame(
	                                                peak_start = min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
	                                                peak_end = max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
	                                                peak_ID = "unknown",
	                                                path_to_cdf_csv = brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1],
	                                                area = sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$tic) - sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$baseline)
	                                            ),
	                                        file = peaks_monolist_path,
	                                        append = TRUE,
	                                        row.names = FALSE,
	                                        col.names = FALSE,
	                                        sep = ","
	                                    )

	                                    output$peak_table <- rhandsontable::renderRHandsontable(rhandsontable::rhandsontable({
	                                        peak_table2 <- read.csv(peaks_monolist_path)
	                                        peak_table2
	                                    }))
	                                }
	                        })

	                    ## Global append peak with "G" (71) keystroke 
	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                            # If selection and "G" is pressed, add the selection to the peak table
	                                if( input$keypress == 71 ) {
	                                
	                                    x_peaks <-  data.frame(
	                                                    peak_start = min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt_rt_offset),
	                                                    peak_end = max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt_rt_offset),
	                                                    peak_ID = "unknown",
	                                                    path_to_cdf_csv = unique(chromatograms_updated$path_to_cdf_csv),
	                                                    area = sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$tic) - sum(brushedPoints(chromatograms_updated, input$chromatogram_brush)$baseline)
	                                                )

	                                    x_peaks$peak_start <- x_peaks$peak_start - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]
	                                    x_peaks$peak_end <- x_peaks$peak_end - chromatograms_updated$rt_offset[match(x_peaks$path_to_cdf_csv, chromatograms_updated$path_to_cdf_csv)]

	                                    write.table(
	                                        x = x_peaks,
	                                        file = peaks_monolist_path,
	                                        append = TRUE,
	                                        row.names = FALSE,
	                                        col.names = FALSE,
	                                        sep = ","
	                                    )

	                                    output$peak_table <- DT::renderDataTable(DT::datatable({
	                                        peak_table <- read.csv(peaks_monolist_path)
	                                        peak_table
	                                    }))
	                                }
	                        })

	                    ## Extract selected MS to MS panel one on "shift+1" (33) keypress
	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                            # If selection and "shift+1" is pressed, extract and print mass spectra
	                                if( input$keypress == 33 ) {

	                                    output$massSpectra_z <- renderPlot({
	                                        framedDataFile <- isolate(as.data.frame(
	                                                            data.table::fread(as.character(
	                                                                brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
	                                                            ))
	                                        ))
	                                        framedDataFile <- isolate(dplyr::filter(
	                                                                framedDataFile, 
	                                                                rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
	                                                                rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
	                                                            ))
	                                        framedDataFile <- plyr::ddply(framedDataFile, .(mz), summarize, intensity = sum(intensity))
	                                        framedDataFile_1 <<- framedDataFile
	                                        framedDataFile$intensity <- framedDataFile$intensity*100/max(framedDataFile$intensity)
	                                        ggplot() + 
	                                            geom_bar(data = framedDataFile, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
	                                            theme_classic() +
	                                            scale_x_continuous(expand = c(0,0)) +
	                                            scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
	                                            geom_text(data = dplyr::filter(framedDataFile, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
	                                    })
	                                }
	                        })

	                    ## Subtract selected MS as background from MS panel one on "shift+3" (35) keypress
	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                            # If selection and "35" is pressed, extract and print mass spectra
	                                if( input$keypress == 35 ) {

	                                    output$massSpectra_z <- renderPlot({
	                                        framedDataFile_to_subtract <- isolate(as.data.frame(
	                                                            data.table::fread(as.character(
	                                                                brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
	                                                            ))
	                                        ))
	                                        framedDataFile_to_subtract <- isolate(dplyr::filter(
	                                                                framedDataFile_to_subtract, 
	                                                                rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
	                                                                rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
	                                                            ))
	                                        framedDataFile_to_subtract <- plyr::ddply(framedDataFile_to_subtract, .(mz), summarize, intensity = sum(intensity))
	                                        framedDataFile_1$intensity <- framedDataFile_1$intensity - framedDataFile_to_subtract$intensity
	                                        framedDataFile_1 <<- framedDataFile_1
	                                        framedDataFile_1$intensity <- framedDataFile_1$intensity*100/max(framedDataFile_1$intensity)
	                                        ggplot() + 
	                                            geom_bar(data = framedDataFile_1, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
	                                            theme_classic() +
	                                            scale_x_continuous(expand = c(0,0)) +
	                                            scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
	                                            geom_text(data = dplyr::filter(framedDataFile_1, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
	                                    })
	                                }
	                        })

	                    ## Save MS in MS panel one on "shift+5" (37) keypress
	                        
	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                            # If selection and "37" is pressed, extract and print mass spectra
	                                if( input$keypress == 37 ) {
	                                    framedDataFile_1$intensity[framedDataFile_1$intensity < 0] <- 0
	                                    write.csv(framedDataFile_1, "integration_app_spectrum_1.csv", row.names = FALSE)
	                                }
	                        })

	                    ## Extract selected MS to MS panel two on "shift+2" (64) keypress
	                        
	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                            # If selection and "shift+22" is pressed, extract and print mass spectra
	                                if( input$keypress == 64 ) {

	                                    output$massSpectra_x <- renderPlot({
	                                        framedDataFile <- isolate(as.data.frame(
	                                                            data.table::fread(as.character(
	                                                                brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
	                                                            ))
	                                    ))
	                                    framedDataFile <- isolate(dplyr::filter(
	                                                                framedDataFile, 
	                                                                rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
	                                                                rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
	                                                            ))
	                                    framedDataFile <- plyr::ddply(framedDataFile, .(mz), summarize, intensity = sum(intensity))
	                                    framedDataFile_2 <<- framedDataFile
	                                    framedDataFile$intensity <- framedDataFile$intensity*100/max(framedDataFile$intensity)
	                                    ggplot() + 
	                                        geom_bar(data = framedDataFile, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
	                                        theme_classic() +
	                                        scale_x_continuous(expand = c(0,0)) +
	                                        scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
	                                        geom_text(data = dplyr::filter(framedDataFile, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
	                                    })
	                                }
	                        })

	                    ## Subtract selected MS as background from MS panel one on "shift+4" (36) keypress
	                        
	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                                # If selection and "shift+4" is pressed, extract and print mass spectra
	                                    if( input$keypress == 36 ) {

	                                        output$massSpectra_x <- renderPlot({
	                                            framedDataFile_to_subtract <- isolate(as.data.frame(
	                                                                data.table::fread(as.character(
	                                                                    brushedPoints(chromatograms_updated, input$chromatogram_brush)$path_to_cdf_csv[1]
	                                                                ))
	                                            ))
	                                            framedDataFile_to_subtract <- isolate(dplyr::filter(
	                                                                    framedDataFile_to_subtract, 
	                                                                    rt > min(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt),
	                                                                    rt < max(brushedPoints(chromatograms_updated, input$chromatogram_brush)$rt)
	                                                                ))
	                                            framedDataFile_to_subtract <- plyr::ddply(framedDataFile_to_subtract, .(mz), summarize, intensity = sum(intensity))
	                                            framedDataFile_2$intensity <- framedDataFile_2$intensity - framedDataFile_to_subtract$intensity
	                                            framedDataFile_2 <<- framedDataFile_2
	                                            framedDataFile_2$intensity <- framedDataFile_2$intensity*100/max(framedDataFile_2$intensity)
	                                            ggplot() + 
	                                                geom_bar(data = framedDataFile_2, mapping = aes(x = mz, y = intensity), stat = "identity", width = 1) +
	                                                theme_classic() +
	                                                scale_x_continuous(expand = c(0,0)) +
	                                                scale_y_continuous(expand = c(0,0), limits = c(0,110)) +
	                                                geom_text(data = dplyr::filter(framedDataFile_2, intensity > 10), mapping = aes(x = mz, y = intensity+5, label = mz))
	                                        })
	                                    }
	                            })

	                    ## Save MS in MS panel one on "shift+6" (94) keypress

	                        observeEvent(input$keypress, {

	                            # Do nothing if no selection
	                                if(is.null(input$chromatogram_brush)) {
	                                    return()
	                                }

	                            # If selection and "shift+6" is pressed, extract and print mass spectra
	                                if( input$keypress == 96 ) {
	                                    framedDataFile_2$intensity[framedDataFile_2$intensity < 0] <- 0
	                                    write.csv(framedDataFile_2, "integration_app_spectrum_2.csv")
	                                }
	                        })
	                }

	            ## Call the app
	                
	                shinyApp(ui = ui, server = server)

	        }

    #### readChromatograms

        #' Import chromatograms stored as ChemStation exported .csv files
        #'
        #' Allows the user to import chromatograms stored as .csv files. Files should have originated from an Agilent GC system running ChemStation
        #' Files should be in the format of character~var1-var2.csv
        #' "character" is the name of the chromatogram, var1 and var2 are variables you want associated with that chromatogram.
        #' @param dir A directory containing (ONLY) the chromatogram(s) that are to be plotted
        #' @param normalize_level Level at which to normalize chromatogram
        #' @param normalize_range y range to consider during normalization. Useful for normalizing but excluding solvent peak(s), for example
        #' @export
        #' @examples
        #' readChromatograms()

        readChromatograms <- function( dir, normalize_level = c("none", "var1", "var2"), normalize_range = NULL ) {

            ## Get paths to the .csv files in the directory "dir"
                list <- paste(dir, dir(dir), sep = "/")

            ## Read in time scale = retention times
                temp <- read.csv(list[1])
                s <- dim(temp)[1]
                ret <- as.numeric(as.character(temp[3:s,1]))

            ## Read in and process data into a gathered dataframe
                data <- data.frame(V1 = as.numeric(1), value = as.numeric(1), character = as.character("1"), var1 = as.character("1"), var2 = as.character("1"))

                    for (i in 1:length(list)){
                        temp <- read.csv(list[i])
                        s <- dim(temp)[1]
                        temp <- temp[3:s,1:2]
                        value <- as.numeric(as.character(temp[,2]))
                        value <- c(rep(0,as.numeric(gsub("-.*$", "", gsub(".*~", "", list[i])))),value)
                        
                        # character is whatever comes before the "~" in the filename
                            character <- as.character(rep(gsub("~.*$", "", gsub(".*/", "", list[i])),length(value)))
                        
                        # var1 is whatever is in between the "~" and the "-" in the filename
                            var1 <- as.character(rep(gsub("-.*$","", gsub(".*~", "", list[i])),length(value)))

                        # var2 is whatever comes after the "-" in the filename
                            var2 <- as.character(rep(substr(list[i],regexpr("-",list[i])[1]+1,nchar(list[i])-4),length(value)))
                        
                        data <- rbind(data,cbind(ret[1:length(value)],value, character, var1, var2))
                    }

                    data <- data[-1,]
                    colnames(data) <- c("ret", "value", "character", "var1", "var2")
                    rownames(data) <- NULL
                    data$ret <- as.numeric(as.character(data$ret))
                    data$value <- as.numeric(as.character(data$value))

            ## Remove entries with var2 = NA
                data <- data[!is.na(data$ret),]

            ## Normalize by var1
                if (normalize_level == "var1") {
                    
                    if (length(normalize_range) > 0) {
                        xmin <- normalize_range[1]
                        xmax <- normalize_range[2]
                        data <- data[data$ret >= xmin & data$ret <= xmax,]
                    }

                    data$value_norm <- data$value
                    for (i in 1:length(unique(data$var1))) {
                        data$value_norm[data$var1==unique(data$var1)[i]] <- phylochemistry::normalize(
                            data$value[data$var1==unique(data$var1)[i]],
                            old_min = min(as.numeric(data$value[data$var1==unique(data$var1)[i]])),
                            old_max = max(as.numeric(data$value[data$var1==unique(data$var1)[i]])),
                            new_min = 0, new_max = 100)
                    }
                }

            ## Normalize by var2
                if (normalize_level == "var2") {
                    
                    if (length(normalize_range) > 0) {
                        xmin <- normalize_range[1]
                        xmax <- normalize_range[2]
                        data <- data[data$ret >= xmin & data$ret <= xmax,]
                    }

                    data$value_norm <- data$value
                    for (i in 1:length(unique(data$var2))) {
                        data$value_norm[data$var2==unique(data$var2)[i]] <- phylochemistry::normalize(
                            data$value[data$var2==unique(data$var2)[i]],
                            old_min = min(as.numeric(data$value[data$var2==unique(data$var2)[i]])),
                            old_max = max(as.numeric(data$value[data$var2==unique(data$var2)[i]])),
                            new_min = 0, new_max = 100)
                    }
                }

            return(data)
        }

    #### readSpectra

        #' Import mass spectral data
        #'
        #' Used to import one or more mass spectra and collect the data in a dataframe.
        #' The names of the files containing the mass spectra must have the following format: "var2~var3-var4.csv", with var1 being the file's position in the directory
        #' @param dir The directory containing the mass spectral (.csv) files of interest
        #' @export
        #' @examples
        #' readSpectra()

        readSpectra <- function( dir ) {

            ## Set working directory and prepare the data frame
                setwd(dir)
                data <- data.frame(
                    mz = as.numeric(1), abu = as.numeric(1), var1 = as.character("1"),
                    var2 = as.character("1"), var3 = as.character("1"), var4 = as.character("1")
                )

            ## Read in all data
                for (i in 1:length(dir())) {
                    temp <- read.csv(dir()[i])
                    temp <- temp[3:dim(temp)[1],1:2]
                    mz <- as.numeric(as.character(temp[,1]))
                    abu <- 100*as.numeric(as.character(temp[,2]))/(max(as.numeric(as.character(temp[,2]))))
                    var1 <- as.character(rep(as.character(i))) # Var1 is the position in the directory
                    var2 <- as.character(rep(gsub("~.*$", "", dir()[i])))
                    var3 <- as.character(rep(gsub("-.*$", "", gsub(".*~", "", dir()[i]))))
                    var4 <- as.character(rep(gsub("\\..*$", "", gsub(".*-", "", dir()[i]))))
                    data <- rbind(data,cbind(mz,abu,var1,var2,var3,var4))
                }

            ## Clean up the data frame and correct variable types
                data <- data[-1,]
                data$mz <- as.numeric(data$mz)
                data$abu <- as.numeric(data$abu)
                data$var1 <- factor(data$var1, levels=unique(data$var1))
                data$var2 <- factor(data$var2, levels=unique(data$var2))
                data$var3 <- factor(data$var3, levels=unique(data$var3))
                data$var4 <- factor(data$var4, levels=unique(data$var4))

            ## Return data
                return( data )
        }

##### Plotting
    
    #### drawAlignment

        #' Generate a multiple alignment graphic using ggplot
        #'
        #' Allows the user to plot a multiple alignent using ggplot's grammar of graphics
        #' @param infile The alignment to use (fasta file)
        #' @param seqlist A dataframe with metadata for the alignment
        #' @param alignment_labels The name of the column to label the entries in the alignment with
        #' @param wrap TRUE/FALSE whether to use a column of plots to show the alignment
        #' @param wrap_length Length at which to cut the wrap
        #' @param roi TRUE/FALSE whether to only show regions of interest in the alignment
        #' @param roi_data Dataframe defining the regions to show
        #' @param consensus TRUE/FALSE whether to include a consensus line in the alignment
        #' @param consensus_height Height of the consensus line plot
        #' @param funct_assoc TRUE/FALSE whether to show functional associations in the alignment
        #' @param funct_assoc_data Columns in the metadata dataframe to use in searching for functional association (e.g. "c(sterol~other, sterol~cyclo)")
        #' @param funct_assoc_height Height of the functional association lines
        #' @param highlights TRUE/FALSE whether to highlight certain sites in the alignment
        #' @param highlights_data Dataframe specifying which site to highlight
        #' @param hlines TRUE/FALSE whether to divide the alignment with horizontal lines
        #' @param hlines_data Dataframe specifying where to divide alignment
        #' @param tick_spacing Spacing between x-axis ticks
        #' @param ticks_text_size Size of x-axis ticks text
        #' @param order Order in which to display the members of the alignment
        #' @param color_pal Color palette to use
        #' @importFrom phangorn read.aa
        #' @import tidyr
        #' @import ggplot2
        #' @import ggtree  
        #' @export
        #' @examples
        #' drawAlignment()

        drawAlignment <- function(   
                            infile,
                            seqlist,
                            alignment_labels,

                            wrap = FALSE,
                            wrap_length = NULL,

                            roi = FALSE,
                            roi_data = NULL,

                            consensus = FALSE,
                            consensus_height = 5,

                            funct_assoc = FALSE,
                            funct_assoc_data = NULL,
                            funct_assoc_height = 5,

                            highlights = FALSE,
                            highlights_data = NULL,

                            hlines = FALSE,
                            hlines_data = NULL,

                            tick_spacing = 5,
                            ticks_text_size = 30,

                            order = NULL,
                            color_pal = NULL
                        ) {

            if (wrap == FALSE) {
                if (roi == FALSE) {
                    print("Please specify either wrap OR roi")
                    stop()
                }
            }

            if (wrap == TRUE) {
                if (roi == TRUE) {
                    print("Please specify either wrap OR roi")
                    stop()
                }
            }

            # Read alignment, create basic plottable
                AA_phydat <- as.data.frame(as.character(phangorn::read.aa(file = infile, format = "fasta")))
                AA_phydat <- as.data.frame(t(AA_phydat))
                colnames(AA_phydat) <- as.character(seq(1,dim(AA_phydat)[2],1))
                AA_alignment_df <-  cbind(  
                                        data.frame(protein = paste(
                                            seqlist[,which(colnames(seqlist) == alignment_labels)][match(rownames(as.matrix(AA_phydat)), seqlist$accession)]
                                        )),
                                        data.frame(funct = seqlist[,colnames(seqlist)=="Function"][match(rownames(as.matrix(AA_phydat)), seqlist[,1])]),
                                        data.frame(y = rep(0,dim(AA_phydat)[1])),
                                        AA_phydat
                                    )
                AA_alignment_df_plottable <- tidyr::gather(AA_alignment_df, position, residue, 4:dim(AA_alignment_df)[2])
                AA_alignment_df_plottable$position <- as.numeric(as.character(AA_alignment_df_plottable$position))
                AA_alignment_df_plottable <- AA_alignment_df_plottable[order(AA_alignment_df_plottable$protein, AA_alignment_df_plottable$position),]
                head(AA_alignment_df_plottable)

            # Make wrapping list
                if (wrap == TRUE) {
                    wrap_break_points <- list()
                    for (i in 1:(dim(AA_phydat)[2]%/%wrap_length)) {
                        wrap_break_points <- c(wrap_break_points, 1+(wrap_length*(i)))
                    }
                    head(wrap_break_points)
                }

                if (roi == TRUE) {
                    wrap_break_points <- dim(AA_phydat)[2]
                    head(wrap_break_points)
                    wrap_length = dim(AA_phydat)[2]
                }

            # ROI
                if (roi == TRUE) {
                    #Define the ROI and Modify the alignment plottable
                        roi_ranges <- rep("NA",dim(AA_phydat)[2])
                            range_names <- roi_data[,1]
                            for (i in 1:dim(roi_data)[1]) {
                                roi_ranges[as.numeric(as.character(roi_data[i,2])):as.numeric(as.character(roi_data[i,3]))] <- as.character(roi_data[i,1])
                            }
                        AA_alignment_df_plottable$roi_ranges <- rep(roi_ranges, length(unique(AA_alignment_df_plottable$protein)))
                        AA_alignment_df_plottable <- subset(AA_alignment_df_plottable, roi_ranges %in% as.character(range_names))
                }

            # ORDER the sequences to optionally match the alignment with a tree or something
                if (!is.null(order)) {
                    AA_alignment_df_plottable$protein <- factor(AA_alignment_df_plottable$protein, levels=order)
                }
                
            # Initiate the plot(s)
                plot_list <- list()
                for (i in 1:length(wrap_break_points)) {
                    plot_list[[i]] <-   ggplot2::ggplot(data = AA_alignment_df_plottable[AA_alignment_df_plottable$position %in% 
                                                        as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                        as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                    aes(x = position, y = y)
                                        )   
                }

            # HIGHLIGHTS
                if (highlights == TRUE) {
                    if (roi == TRUE) {
                        highlights_data <- cbind(highlights_data, data.frame(roi_ranges=roi_ranges[highlights_data$xint]))
                    }
                    for (i in 1:length(wrap_break_points)) {
                        plot_list[[i]] <- plot_list[[i]] +  geom_vline(
                                                                data=highlights_data[highlights_data$xint %in% 
                                                                    as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                                    as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                                aes(xintercept=xint, color=as.character(xint)),
                                                                size=3
                                                            )
                    }
                }
                
            # CONSENSUS
                if (consensus == TRUE) {
                    # Calculate consensus scores
                        #Calculate max possible score:
                            AA_phydat_char_distrib <- t(as.data.frame(lapply(apply(AA_phydat,2,table), sd)))
                            AA_phydat_char_distrib[is.na(AA_phydat_char_distrib)] <- 0
                            max_consensus_score <- ceiling(max(AA_phydat_char_distrib))

                        # Calculate consensus scores:
                            AA_phydat_char_distrib <- t(as.data.frame(lapply(apply(AA_phydat,2,table), sd)))
                            AA_phydat_char_distrib[is.na(AA_phydat_char_distrib)] <- max_consensus_score
                            colnames(AA_phydat_char_distrib) <- "consensus_score"
                            AA_phydat_char_distrib <- normalize(AA_phydat_char_distrib)*consensus_height
                            
                    # Bind consensus score
                        consensus_protein <- data.frame(
                            protein = "consensus",
                            funct="none",
                            y=AA_phydat_char_distrib[,1],
                            position=seq(1,dim(AA_phydat)[2]),
                            residue="+"
                        )
                        head(consensus_protein)

                    # Truncate the consensus_protein to fit roi
                        if (roi == TRUE) {
                            consensus_protein <- cbind(
                                consensus_protein, 
                                data.frame(roi_ranges=  rep(    roi_ranges,
                                                                length(unique(consensus_protein$protein))
                                                        )
                                )
                            )
                            consensus_protein <- subset(consensus_protein, roi_ranges %in% range_names)
                        }

                    # Add consensus_protein to the plot(s) in plot_list
                        for (i in 1:length(wrap_break_points)) {
                            plot_list[[i]] <-  plot_list[[i]] + layer(  data=consensus_protein[consensus_protein$position %in% 
                                                                                as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                                                as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                                            aes(x=position, y=y), geom="line", stat="identity", position="identity"
                                                                        )
                        }
                }

            # FUNCTIONAL ASSOCIATION
                if (funct_assoc == TRUE) {
                    for (j in 1:length(funct_assoc_data)) {
                        # Subset alignment matrices for the two variables
                            alignment_matrix <- toupper(t(as.data.frame(as.character(phangorn::read.aa(file = infile, format = "fasta")))))
                            
                            spec1 <- list()
                            spec1[[1]] <- seqlist[,colnames(seqlist) == funct_assoc_data[j]] == as.character(gsub("~.*$", "", funct_assoc_data[j]))
                            alignment_matrix_dim1 <- alignment_matrix[Reduce("&", spec1),]
                            dim(alignment_matrix_dim1)

                            spec2 <- list()
                            spec2[[1]] <- seqlist[,colnames(seqlist)==funct_assoc_data[j]] == as.character(gsub(".*~", "", funct_assoc_data[j]))
                            alignment_matrix_dim2 <- alignment_matrix[Reduce("&", spec2),]
                            dim(alignment_matrix_dim2)

                        #Set up fixed substitution table
                            fixed_sub_table <- data.frame(funct = c(as.character(gsub("~.*$", "", funct_assoc_data[j])), as.character(gsub(".*~", "", funct_assoc_data[j]))), code=c(1,2))

                        #Get frequency stats from the first alignment subset
                            alignment_matrix_dim1_char_distrib <- lapply(apply(alignment_matrix_dim1,2,table), sort, decreasing=TRUE)
                            alignment_matrix_dim2_char_distrib <- lapply(apply(alignment_matrix_dim2,2,table), sort, decreasing=TRUE)

                        #Set up empty assoc matrix
                            association_list <- list()

                            for (i in 1:dim(alignment_matrix)[2]){

                                # Set up dim1 residue substitution table
                                    base_sub_table_dim1 <- data.frame(letter=LETTERS, number=rep(2,26))
                                    base_sub_table_dim1 <- rbind(base_sub_table_dim1, data.frame(letter="-", number=2))
                                    sub_position_dim1 <- data.frame(consensus=names(alignment_matrix_dim1_char_distrib[[i]][1]), number=1)
                                    positional_sub_table_dim1 <- base_sub_table_dim1
                                    positional_sub_table_dim1$number <- sub_position_dim1$number[match(positional_sub_table_dim1$letter, sub_position_dim1$consensus)]
                                    positional_sub_table_dim1[is.na(positional_sub_table_dim1$number),]$number <- 2
                                    positional_sub_table_dim1

                                # Substitute in the first alignment subset
                                    positional_data_frame_dim1 <- data.frame(
                                        enzyme=names(alignment_matrix_dim1[,i]),
                                        funct=as.character(gsub("~.*$", "", funct_assoc_data[j])),
                                        x_num=as.character(gsub("~.*$", "", funct_assoc_data[j])),
                                        position=i,
                                        residue=alignment_matrix_dim1[,i],
                                        y_num=alignment_matrix_dim1[,i]
                                    )

                                    positional_data_frame_dim1$x_num <- fixed_sub_table$code[match(positional_data_frame_dim1$x_num, fixed_sub_table$funct)]
                                    positional_data_frame_dim1$y_num <- positional_sub_table_dim1$number[match(positional_data_frame_dim1$y_num, positional_sub_table_dim1$letter)]
                                    positional_data_frame_dim1

                                # Set up dim2 residue substitution table
                                    base_sub_table_dim2 <- data.frame(letter=LETTERS, number=rep(1,26))
                                    base_sub_table_dim2 <- rbind(base_sub_table_dim2, data.frame(letter="-", number=1))
                                    sub_position_dim2 <- data.frame(consensus=names(alignment_matrix_dim2_char_distrib[[i]][1]), number=2)
                                    if (as.character(sub_position_dim2$consensus) == as.character(sub_position_dim1$consensus)) {sub_position_dim2$number <- 1}
                                    positional_sub_table_dim2 <- base_sub_table_dim2
                                    positional_sub_table_dim2$number <- sub_position_dim2$number[match(positional_sub_table_dim2$letter, sub_position_dim2$consensus)]
                                    positional_sub_table_dim2[is.na(positional_sub_table_dim2$number),]$number <- 1
                                    positional_sub_table_dim2

                                #Substitute in the second alignment subset
                                    positional_data_frame_dim2 <- data.frame(
                                        enzyme=names(alignment_matrix_dim2[,i]),
                                        funct=as.character(gsub(".*~", "", funct_assoc_data[j])),
                                        x_num=as.character(gsub(".*~", "", funct_assoc_data[j])),
                                        position=i,
                                        residue=alignment_matrix_dim2[,i],
                                        y_num=alignment_matrix_dim2[,i]
                                    )

                                    positional_data_frame_dim2$x_num <- fixed_sub_table$code[match(positional_data_frame_dim2$x_num, fixed_sub_table$funct)]
                                    positional_data_frame_dim2$y_num <- positional_sub_table_dim2$number[match(positional_data_frame_dim2$y_num, positional_sub_table_dim2$letter)]

                                #Add the positional_data_frame to the list
                                    association_list[[i]] <- rbind(positional_data_frame_dim1, positional_data_frame_dim2)
                            }

                            # USE THIS LINE TO MANUALLY INSPECT SITE ASSOCIATAIONS
                                # association_list[[273]]

                            cor_list <- list()
                            for (i in 1:length(association_list)) {
                                cor_list[i] <- cor(cbind(association_list[[i]]$x_num,association_list[[i]]$y_num))[1,2]
                            }

                            cor_list[is.na(cor_list)] <- 0
                            cor_list[cor_list < 0] <- 0
                            correlation_protein <-  data.frame(
                                                        protein = as.character(funct_assoc_data[j]), 
                                                        funct = "NA", 
                                                        y = t(as.data.frame(cor_list))*funct_assoc_height, 
                                                        position = seq(1,length(cor_list),1), 
                                                        residue = "NA"
                                                    )
                            head(correlation_protein)

                            # Truncate the correlation_protein plottable so it only contains things in the roi
                                if (roi == TRUE) {
                                    correlation_protein <- cbind(
                                        correlation_protein, 
                                        data.frame(roi_ranges=rep(
                                            roi_ranges,
                                            length(unique(correlation_protein$protein)))
                                        )
                                    )
                                    correlation_protein <- subset(correlation_protein, roi_ranges %in% range_names)
                                }

                            # Assign correlation protein to global environment
                                assign("correlation_protein", correlation_protein, envir = .GlobalEnv)

                        # Add the correlation_protein to the plot(s) in plot_list
                            for (i in 1:length(wrap_break_points)) {
                                plot_list[[i]] <-  plot_list[[i]] + layer(data=correlation_protein[correlation_protein$position %in% 
                                                                            as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                                                            as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                                                        aes(x=position, y=y), geom="line", stat="identity", position="identity"
                                                                    )
                            }
                    }
                }

            # HLINES
                if (hlines == TRUE) {
                    for (i in 1:length(wrap_break_points)) {
                        plot_list[[i]] <- plot_list[[i]] +  geom_hline(
                                                                data = hlines_data,
                                                                aes(yintercept = yint),
                                                                size = 1,
                                                                linetype = 2
                                                            )
                    }
                }

            # Build the rest of the alignment
                #Specify facets based on presence/absence of ROIs
                    if (roi == TRUE) { 
                        for (i in 1:length(wrap_break_points)) {
                            plot_list[[i]] <- plot_list[[i]] + facet_grid(protein~roi_ranges, scales="free", space="free", switch="y") 
                        }
                    }
                    if (wrap == TRUE) { 
                        for (i in 1:length(wrap_break_points)) {
                            plot_list[[i]] <- plot_list[[i]] + facet_grid(protein~., scales="free", space="free", switch="y") 
                        }
                    }
                    
                # Make rest of plot
                    for (i in 1:length(wrap_break_points)) {
                        plot_list[[i]] <- plot_list[[i]] + theme_classic()
                        plot_list[[i]] <- plot_list[[i]] +  scale_x_continuous(
                                                                expand = c(0.005,0.005),
                                                                name = "",
                                                                # breaks = as.numeric(generateTicks(seq(0,5000,tick_spacing))$all_breaks),
                                                                # labels = as.character(generateTicks(seq(0,5000,tick_spacing))$all_labels)
                                                            )
                        plot_list[[i]] <- plot_list[[i]] + scale_y_continuous(name = "")
                        plot_list[[i]] <- plot_list[[i]] + scale_color_manual(values = color_pal)
                        plot_list[[i]] <- plot_list[[i]] + theme(
                                    panel.spacing.y = unit(0.01, "lines"),
                                    panel.spacing.x = unit(0.2, "cm"),
                                    panel.border = element_blank(),
                                    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                                    axis.text.y = element_blank(),
                                    strip.text.y = element_text(angle = 180),
                                    text = element_text(size = ticks_text_size),
                                    legend.position = "none",
                                    axis.title = element_text(color = "#737373", face = "bold", size = 40),
                                    axis.ticks.length = unit(0.2, "cm"),
                                    axis.ticks = element_line(color = "#737373", size = 1, lineend = 6),
                                    axis.text = element_text(color = "#737373", face = "bold", size = ticks_text_size),
                                    axis.line = element_line(color = "#737373", size = 1),
                                    strip.text = element_text(hjust = 1, size = ticks_text_size),
                                    strip.background = element_blank(),
                                    strip.placement = "outside"
                                )
                    }

            # Return the Alignment
                if (wrap == TRUE) { ## IF USING WRAP, RETURNS A LIST THAT MUST BE PLOTTED USING do.call(gridExtra::grid.arrange,c(alignment, ncol=1))
                    plots <- list()
                    for (i in 1:length(wrap_break_points)) {
                        plots[[i]] <- ggplot2::ggplot_gtable(
                            ggplot2::ggplot_build(plot_list[[i]] + layer(
                                data = AA_alignment_df_plottable[AA_alignment_df_plottable$position %in% 
                                    as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                                    as.numeric(as.character(unlist(wrap_break_points)[i])),],
                                geom = "text", 
                                mapping = aes(label = residue),
                                stat = "identity",
                                position = "identity")
                            )
                        )
                    }
                    return(plots)
                }

                if (roi == TRUE) { ## IF USING ROI, RETURNS A GGPLOT THAT CAN BE PLOTTED USING NORMAL GGPLOT PLOTTING TECHNIQUES
                    plot_list[[i]] <- plot_list[[i]] + layer(
                        data = AA_alignment_df_plottable[AA_alignment_df_plottable$position %in% 
                            as.numeric(as.character(unlist(wrap_break_points)[i]-wrap_length)):
                            as.numeric(as.character(unlist(wrap_break_points)[i])),],
                        geom = "text", 
                        mapping = aes(label = residue),
                        stat = "identity",
                        position = "identity")
                    return(plot_list[[1]])
                }           
        }

    #### generateTicks

        #' Create major and minor axes ticks and labels
        #'
        #' @param major_ticks Values at which major ticks should be generated
        #' @param minor_freq Number of minor ticks between each major tick 
        #' @examples
        #' @export
        #' generateTicks

        generateTicks <-    function(major_ticks, minor_freq = 4, major_tick_size, minor_tick_size) {
                                
                                major_labels <- vector()
                                
                                for (tick in 1:(length(major_ticks)-1)) {
                                    major_labels <- c(major_labels, major_ticks[tick])
                                    major_labels <- c(major_labels, rep("", (minor_freq)))
                                }
                                
                                all_ticks <- vector()
                                
                                for (tick in 1:(length(major_ticks)-1)) {
                                    all_ticks <- c(all_ticks, major_ticks[tick])
                                    minor_tick_values <- vector()
                                    for (minor_tick in 1:minor_freq) {
                                        minor_tick_values <- c(minor_tick_values, major_ticks[tick] + minor_tick*(((major_ticks[tick+1] - major_ticks[tick])/(minor_freq+1))))
                                    }
                                    all_ticks <- c(all_ticks, minor_tick_values)
                                }
                                
                                major_labels <- c(major_labels, major_ticks[length(major_ticks)])
                                all_ticks <- c(all_ticks, major_ticks[length(major_ticks)])

                                return <- data.frame(all_ticks = as.numeric(all_ticks), major_labels = as.character(major_labels))
                                return$tick_size <- minor_tick_size
                                return$tick_size[return$major_labels != ""] <- major_tick_size

                                return(return)
                            }

##### Phylogenetic statistical testing

    #### phylogeneticSignal

        #' Convert mass spectral datafiles (CDF) into a csv file
        #'
        #' @param trait A vector of the variable, named according to which tree tip it comes from
        #' @param tree A phylogenetic tree with tips that exactly match the names of trait
        #' @param replicates Number of random replications to run
        #' @param cost Optional. A specialized transition matrix
        #' @examples
        #' @export
        #' phylogeneticSignal

            phylogeneticSignal <- function( trait, tree, replicates = 999, cost = NULL ) {

                ### For discrete traits

                    ## Get the states in which the trait may exist (levels)

                        levels <- attributes(factor(trait))$levels

                    ## Chech that the variable is indeed categorical
                                    
                        if (length(levels) == length(trait)) {
                        
                            warn("Are you sure this variable is categorical?")

                        }

                    ## Make the transition matrix
                        
                        if (is.null(cost)) {
                        
                            cost1 <- 1-diag(length(levels))
                        
                        } else {
                        
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
                        
                        output <- data.frame(
                            number_of_levels = length(attributes(factor(trait))$levels), 
                            evolutionary_transitions_observed = OBS,
                            median_evolutionary_transitions_in_randomization = median(null_model),
                            minimum_evolutionary_transitions_in_randomization = min(null_model),
                            evolutionary_transitions_in_randomization = max(null_model),
                            p_value = p_value
                        )
                        
                    return(output)

                }

