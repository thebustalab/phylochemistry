## sequence — moved verbatim from phylochemistry.R (package conversion, commit 1)

        # tukeyHSD <- function( data, ... ) {rstatix::tukey_hsd( data = data, ... )}

    ## Set up lookups

        ## Phred_ascii_33 lookup
            #' Phred+33 quality score lookup table
            #'
            #' @format A data frame mapping Phred+33 characters to quality scores.
            #' @export
            phred33_lookup <- data.frame(rbind(
                c("!","0"),c("\"","1"),c("#","2"),c("$","3"),c("%","4"),c("&","5"),c("'","6"),c("(","7"),c(")","8"),c("*","9"),c("+","10"),c(",","11"),c("-","12"),c(".","13"),c("/","14"),c("0","15"),c("1","16"),c("2","17"),c("3","18"),c("4","19"),c("5","20"),c("6","21"),c("7","22"),c("8","23"),c("9","24"),c(":","25"),c(";","26"),c("<","27"),c("=","28"),c(">","29"),c("?","30"),c("@","31"),c("A","32"),c("B","33"),c("C","34"),c("D","35"),c("E","36"),c("F","37"),c("G","38"),c("H","39"),c("I","40"),c("J","41"),c("K","42"),c("L","43"),c("M","44"),c("N","45"),c("O","46"),c("P","47"),c("Q","48"),c("R","49"),c("S","50"),c("T","51"),c("U","52"),c("V","53"),c("W","54"),c("X","55"),c("Y","56"),c("Z","57"),c("[","58"),c("\\","59"),c("]","60"),c("^","61"),c("_","62"),c("`","63"),c("a","64"),c("b","65"),c("c","66"),c("d","67"),c("e","68"),c("f","69"),c("g","70"),c("h","71"),c("i","72"),c("j","73"),c("k","74"),c("l","75"),c("m","76"),c("n","77"),c("o","78"),c("p","79"),c("q","80"),c("r","81"),c("s","82"),c("t","83"),c("u","84"),c("v","85"),c("w","86"),c("x","87"),c("y","88"),c("z","89"),c("{","90"),c("|","91"),c("}","92"),c("~","93")
            ))


        #### readAlignment


            #' Reads an alignment into a tidy dataframe
            #'
            #' @param alignment_in_path The path to the alignment
            #' @param type Not yet documented.
            #' @export
            readAlignment <- function( alignment_in_path, type = c("DNA", "AA")) {

                ## DNA v AA

                    if (length(type) > 1) {
                        stop("Please specify a single type of alignment: DNA or AA")
                    }

                ## Read and return the alignment
                    
                    if ( type == "DNA" ) {
                        alignment <- phangorn::read.phyDat(
                            file = alignment_in_path,
                            format = "fasta", type = "DNA"
                        )
                    }
                    if ( type == "AA" ) {
                        alignment <- phangorn::read.phyDat(
                            file = alignment_in_path,
                            format = "fasta", type = "AA"
                        )   
                    }
                    alignment <- t(as.matrix(as.data.frame(alignment)))
                    alignment <- as_tibble(data.frame(cbind(rownames(alignment), alignment)))
                    alignment <- pivot_longer(alignment, cols = c(2:dim(alignment)[2]), names_to = "position", values_to = "state")
                    alignment$position <- as.numeric(as.character(gsub("X", "", alignment$position)))
                    colnames(alignment)[1] <- "name"
                    
                    return(alignment)
            }


        #### writeAlignment

            #' Write Alignment DataFrame to FASTA File
            #'
            #' This function writes a long-format alignment DataFrame to a FASTA file. The DataFrame must contain
            #' specific columns: 'name', 'position', and 'state'. The function can automatically detect the type
            #' of sequence data (DNA or AA) if not specified.
            #'
            #' @param alignment A DataFrame containing the alignment data with columns 'name', 'position', and 'state'.
            #' @param alignment_out_path A string specifying the file path where the FASTA file will be saved.
            #' @param type A string indicating the type of sequence data: "DNA", "AA", or "auto" for automatic detection.
            #'             Default is "auto".
            #' @return None. The function writes the alignment to a specified FASTA file.
            #' @examples
            #' \dontrun{
            #' \dontrun{
            #' # Write a DNA alignment to a FASTA file
            #' writeAlignment(alignment_df, "output.fasta", type = "DNA")
            #'
            #' # Automatically detect the type and write to a FASTA file
            #' writeAlignment(alignment_df, "output.fasta")
            #' }
            #' }
            #' @export
            writeAlignment <- function(alignment, alignment_out_path, type = "auto") {
                .check_bioc(c("Biostrings"), "writeAlignment")
                
                # Input validation
                if (!all(c("name", "position", "state") %in% colnames(alignment))) {
                    stop("Alignment DataFrame must contain columns: 'name', 'position', and 'state'")
                }

                if (length(type) > 1) {
                    stop("Please specify a single type of alignment: DNA, AA, or auto")
                }

                # Auto-detect type if set to "auto"
                if (type == "auto") {
                    dna_letters <- c("A", "T", "C", "G", "N", "-")
                    aa_letters <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*", "-")
                    
                    unique_states <- unique(toupper(alignment$state))
                    if (all(unique_states %in% dna_letters)) {
                        type <- "DNA"
                    } else if (all(unique_states %in% aa_letters)) {
                        type <- "AA"
                    } else {
                        stop("Unable to auto-detect type. Please specify 'DNA' or 'AA'.")
                    }
                }

                # Convert long format to wide format
                alignment_wide <- pivot_wider(
                    alignment,
                    names_from = position,
                    values_from = state
                )

                # Combine sequence positions into single strings
                sequences <- apply(alignment_wide[,-1], 1, paste, collapse="")
                names(sequences) <- alignment_wide$name

                # Convert to appropriate StringSet type
                if (type == "DNA") {
                    sequences <- Biostrings::DNAStringSet(sequences)
                } else if (type == "AA") {
                    sequences <- Biostrings::AAStringSet(sequences)
                }

                # Write to file
                writeFasta(sequences, alignment_out_path, fasta_type = type)
            }


        #### downloadFasta


            #' Returns contents of one fasta file as a StringSet object.
            #'
            #' @param share_link The path to the fasta to read. Can be local or a google drive share link.
            #' @export
            downloadFasta <- function(share_link) {
                rlang::check_installed(c("googledrive"), reason = "to use downloadFasta()")
                cat("Downloading... this may take a while...")
                suppressMessages(googledrive::drive_download(
                    file = googledrive::as_id(share_link),
                    path = "TEMP___fasta",
                    overwrite = TRUE
                ))
            }


        #### readFasta


            #' Returns contents of one fasta file as a StringSet object.
            #'
            #' @param fasta_in_path The path to the fasta to read. Can be local or a google drive share link.
            #' @param fasta_type The type of sequence data contained in the fasta file. Options: "nonspecific", "DNA", "RNA", or "AA"
            #' @export
            readFasta <- function( fasta_in_path, fasta_type = "nonspecific" ) {
                .check_bioc(c("Biostrings"), "readFasta")

                ## If its on google drive, download it and make a temp file
                    if (length(grep("http", fasta_in_path)) > 0) {
                        
                        downloadFasta(share_link = fasta_in_path)

                        temp_fasta <- Biostrings::readBStringSet("TEMP___fasta")
                        
                        file.remove("TEMP___fasta")

                        return(temp_fasta)
                    }

                ## Read the file in
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


        #### readManyFastas


            #' Returns the contents of many fasta files as a StringSet object.
            #'
            #' @param fasta_in_paths The paths to the fastas to read
            #' @param fasta_type The type of sequence data contained in the fasta file. Options: "nonspecific", "DNA", "RNA", or "AA"
            #' @export
            readManyFastas <- function( fasta_in_paths, fasta_type = "nonspecific" ) {
                .check_bioc(c("Biostrings"), "readManyFastas")

        		if ( fasta_type == "nonspecific" ) {
            		temp <- Biostrings::BStringSet()
            	}

            	if ( fasta_type == "DNA" ) {
            		temp <- Biostrings::DNAStringSet()
            	}

            	if ( fasta_type == "RNA" ) {
            		temp <- Biostrings::RNAStringSet()
            	}

            	if ( fasta_type == "AA" ) {
            		temp <- Biostrings::AAStringSet()
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


        #### writeManyFastas


            #' Writes one or more sequences as individual fasta file(s), each with one sequence
            #'
            #' @param XStringSet The sequence(s) to write out
            #' @param fasta_out_directory_path The path to the directory where the file(s) should be written
            #' @export
            writeManyFastas <- function( XStringSet, fasta_out_directory_path ) {
                .check_bioc(c("Biostrings"), "writeManyFastas")
            	for (sequence in 1:length(XStringSet)) {
            		Biostrings::writeXStringSet( x = XStringSet[sequence], filepath = paste0(fasta_out_directory_path, XStringSet@ranges@NAMES[sequence], ".fa") )
            	}
            }


        #### writeFasta


            #' Writes one or more sequences as a single fasta file
            #'
            #' @param XStringSet The sequences to write out
            #' @param fasta_out_path The path to where the file should be written
            #' @param append Not yet documented.
            #' @param fasta_type Not yet documented.
            #' @export
            writeFasta <- function( XStringSet, fasta_out_path, append = FALSE, fasta_type = c("DNA", "AA")) {
                .check_bioc(c("Biostrings"), "writeFasta")

                if (fasta_type[1] == "DNA") {
            	   Biostrings::writeXStringSet( x = Biostrings::DNAStringSet(XStringSet), filepath = fasta_out_path, append = append )
                }

                if (fasta_type[1] == "AA") {
                    Biostrings::writeXStringSet( x = XStringSet, filepath = fasta_out_path, append = append )
                }

            }


        #### gblocks


            #' Subset an alignment
            #'
            #' @param alignment_in_path The path to the information that should be used to build the tree.
            #' @param max_gap_percent The maximum percentage of gaps allowed in a column of the alignment. Default is 0.5.
            #' @param min_conservation_percent The minimum percentage of conservation allowed in a column of the alignment. Default is 0.3.
            #' @param sequence_type The type of sequence in the alignment, either "DNA" or "AA"
            #' @export
            gblocks <- function(
                alignment_in_path = NULL,
                max_gap_percent = 0.5,
                min_conservation_percent = 0.3,
                sequence_type = c("DNA", "AA")
            ) {
                .check_bioc(c("Biostrings"), "gblocks")

                ## Read in the alignment and make it a matrix
                    if (sequence_type == "DNA") {
                        nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(alignment_in_path), format = "fasta", type = "DNA")
                        nucl_seqs_aligned_matrix <- t(as.matrix(as.data.frame(nucl_seqs_aligned)))
                        cat(paste0("Alignment is ", dim(nucl_seqs_aligned_matrix)[2], " positions long.\n"))
                    }
                    if (sequence_type == "AA") {
                        nucl_seqs_aligned <- phangorn::read.phyDat(file = paste(alignment_in_path), format = "fasta", type = "AA")
                        nucl_seqs_aligned_matrix <- t(as.matrix(as.data.frame(nucl_seqs_aligned)))
                        cat(paste0("Alignment is ", dim(nucl_seqs_aligned_matrix)[2], " positions long.\n"))
                    }

                ## Analyze composition of each position
                    position_scores <- list()
                    for (i in 1:dim(nucl_seqs_aligned_matrix)[2]) {
                        position_scores[[i]] <- data.frame(
                            position = i,
                            gap_percent = sum(nucl_seqs_aligned_matrix[,i] == "-", na.rm = TRUE) / dim(nucl_seqs_aligned_matrix)[1],
                            conservation_percent = suppressWarnings(sum(nucl_seqs_aligned_matrix[,i] == mode(nucl_seqs_aligned_matrix[,i]), na.rm = TRUE) / dim(nucl_seqs_aligned_matrix)[1])
                        )
                    }
                    position_scores <- do.call(rbind, position_scores)

                    print(
                        ggplot() +
                            geom_line(data = position_scores, aes(x = position, y = -gap_percent+1.5-max_gap_percent), color = "gold", alpha = 0.8) +
                            geom_line(data = position_scores, aes(x = position, y = conservation_percent+0.5-min_conservation_percent), color = "maroon", alpha = 0.4) +
                            geom_hline(yintercept = 0.5) +
                            coord_cartesian(xlim = c(
                                which(movingAverage(position_scores$conservation_percent) < 0.9)[1],
                                (dim(position_scores)[1] - which(rev(movingAverage(position_scores$conservation_percent) < 0.9))[1])
                            )) +
                            geom_text(aes(x = dim(position_scores)[1]*0.45, y = 1, label = "Positions with gold and red above the line are kept.")) +
                            theme_bw()
                    )

                ## Ask user if these settings are okay

                    if ( !if (interactive()) askYesNo("Do these gblocks parameters look okay?") ) {
                        stop()
                    }

                ## Find all positions with less than threshold conservation and remove them
                    conservation_reject <- position_scores$position[position_scores$conservation_percent < min_conservation_percent]
                    cat(paste0("Removing ", length(conservation_reject), " positions that are too conserved.\n"))
                    gap_reject <- position_scores$position[position_scores$gap_percent > max_gap_percent]
                    cat(paste0("Removing ", length(gap_reject), " positions that contain too many gaps.\n"))
                    
                    reject_positions <- unique(c(conservation_reject, gap_reject))
                    reject_positions[order(reject_positions)]

                    nucl_seqs_aligned_matrix <- nucl_seqs_aligned_matrix[,-c(reject_positions)]
                    cat(paste0("Remaining positions: ", dim(nucl_seqs_aligned_matrix)[2]), "\n")

                ## Write out as fasta
                    if (sequence_type == "DNA") {
                        blocked_alignment <- Biostrings::DNAStringSet()
                        for (i in 1:dim(nucl_seqs_aligned_matrix)[1]) {
                            blocked_alignment <- c(blocked_alignment, Biostrings::DNAStringSet(paste0(nucl_seqs_aligned_matrix[i,], collapse = "")))
                        }
                        blocked_alignment@ranges@NAMES <- rownames(nucl_seqs_aligned_matrix)
                        writeFasta(blocked_alignment, paste0(alignment_in_path, "_blocked"))
                    }
                    if (sequence_type == "AA") {
                        blocked_alignment <- Biostrings::AAStringSet()
                        for (i in 1:dim(nucl_seqs_aligned_matrix)[1]) {
                            blocked_alignment <- c(blocked_alignment, Biostrings::AAStringSet(paste0(nucl_seqs_aligned_matrix[i,], collapse = "")))
                        }
                        blocked_alignment@ranges@NAMES <- rownames(nucl_seqs_aligned_matrix)
                        writeFasta(blocked_alignment, paste0(alignment_in_path, "_blocked"))
                    }

                    blocked_alignment <- t(as.matrix(as.data.frame(blocked_alignment)))
                    position_scores <- list()
                    for (i in 1:dim(blocked_alignment)[2]) {
                        position_scores[[i]] <- data.frame(
                            position = i,
                            gap_percent = sum(blocked_alignment[,i] == "-", na.rm = TRUE) / dim(blocked_alignment)[1],
                            conservation_percent = suppressWarnings(sum(blocked_alignment[,i] == mode(blocked_alignment[,i]), na.rm = TRUE) / dim(nucl_seqs_aligned_matrix)[1])
                        )
                    }
                    position_scores <- do.call(rbind, position_scores)

            }


    ##### Sequence data handling

        #### annotateAminoAcidSequence    

            #' Annotate an amino acid sequence with InterProScan
            #'
            #' @param sequence Not yet documented.
            #' @param email Not yet documented.
            #' @param applications Not yet documented.
            #' @export
            annotateAminoAcidSequence <- function(sequence, email, applications = "PfamA") {
              base <- "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
              
              # submit (EBI requires an email as the abuse contact)
              job <- content(POST(paste0(base, "/run"),
                                  body = list(email = email, stype = "p",
                                              appl = applications, goterms = "true",
                                              sequence = sequence),
                                  encode = "form"), "text", encoding = "UTF-8")
              message("job = ", job)
              
              # poll
              repeat {
                Sys.sleep(5)
                st <- content(GET(paste0(base, "/status/", job)), "text", encoding = "UTF-8")
                if (st == "FINISHED") break
                if (st %in% c("ERROR", "FAILURE", "NOT_FOUND")) stop("InterProScan ", st, ": ", job)
              }
              
              # fetch the TSV result and parse to a tibble
              tsv <- content(GET(paste0(base, "/result/", job, "/tsv")), "text", encoding = "UTF-8")
              cols <- c("query","md5","length","analysis","sig_acc","sig_desc","start","stop",
                        "score","status","date","interpro_acc","interpro_desc","go","pathways")
              as_tibble(read.delim(text = tsv, header = FALSE, sep = "\t",
                                   col.names = cols, fill = TRUE, quote = ""))
            }


        #### blastNCBI

            #' Run a BLAST search at NCBI
            #'
            #' @param query_sequence Not yet documented.
            #' @param program Not yet documented.
            #' @param database Not yet documented.
            #' @export
            blastNCBI <- function(query_sequence, program, database) {

                # Reference: https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html#urlapi
                # Adjust program parameter if necessary
                if (program == "megablast") { 
                    program <- "blastn&MEGABLAST=on" 
                } else if (program == "rpsblast") { 
                    program <- "blastp&SERVICE=rpsblast" 
                } else if (program == "blastp") {
                    program <- "blastp"  # Enable protein->protein BLAST
                }

                # Allow friendly aliases for the clustered nr database (ClusteredNR).
                # Searches run against one representative sequence per 90/90 cluster, and
                # each hit title is the cluster's well-annotated representative protein.
                if (database %in% c("clusteredNR", "clustered_nr", "clusterednr")) {
                    database <- "nr_cluster_seq"
                }
                
                # Build the request and post the initial request to start the search
                url <- "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
                args <- list(
                    CMD = "Put", PROGRAM = program, DATABASE = database,
                    QUERY = query_sequence
                )
                response <- POST(url, body = args, encode = "form")
                
                # Parse the RID (request ID) and RTOE (how long we have to wait)
                content_text <- content(response, "text", encoding = "UTF-8")
                rid <- sub(".*RID = (\\S+).*", "\\1", content_text)
                rtoe <- as.numeric(sub(".*RTOE = (\\d+).*", "\\1", content_text))
                Sys.sleep(rtoe)
                
                # Poll for results
                repeat {
                    Sys.sleep(5)
                    poll_response <- GET(url, query = list( CMD = "Get", FORMAT_OBJECT = "SearchInfo", RID = rid ))
                    poll_content <- content(poll_response, "text", encoding = "UTF-8")
                    if (grepl("\\s+Status=WAITING", poll_content)) { next }
                    if (grepl("\\s+Status=FAILED", poll_content)) { stop(sprintf("Search %s failed; please report to blast-help@ncbi.nlm.nih.gov.", rid)) }
                    if (grepl("\\s+Status=UNKNOWN", poll_content)) { stop(sprintf("Search %s expired.", rid)) }
                    if (grepl("\\s+Status=READY", poll_content)) { if (grepl("\\s+ThereAreHits=yes", poll_content)) { break } else { stop("No hits found.") } }
                    stop("Unexpected status encountered.") # Unexpected status
                }
                
                # Retrieve and display results
                result_response <- GET(url, query = list( CMD = "Get", FORMAT_TYPE = "Text", RID = rid ))
                result_content <- content(result_response, "text", encoding = "UTF-8")
                split_lines <- unlist(strsplit(result_content, "\n"))
                
                # Process human readable results file
                identities_indices <- grep("Identities", split_lines)
                results <- data.frame(Score_Line = character(), Identities_Line = character(), Greater_Than_Line = character(), stringsAsFactors = FALSE)
                
                for (index in identities_indices) {
                    
                    # Get the score line
                    identity_line <- split_lines[index]
                    
                    # Look upwards for the first line that starts with ">"
                    greater_than_line <- NA
                    for (i in seq(index - 1, 1, by = -1)) {
                        if (startsWith(split_lines[i], ">")) {
                            greater_than_line <- split_lines[i]
                            break
                        }
                    }
                    
                    # Look upwards for the first line that starts with " Score"
                    score_line <- NA
                    for (i in seq(index - 1, 1, by = -1)) {
                        if (startsWith(split_lines[i], " Score")) {
                            score_line <- split_lines[i]
                            break
                        }
                    }
                    
                    # Append to results data frame
                    results <- rbind(results, data.frame(Score_Line = score_line, Identities_Line = identity_line, Greater_Than_Line = greater_than_line, stringsAsFactors = FALSE))
                }
                
                # Parse results
                parsed_results <- results %>%
                    mutate(
                        Score = as.numeric(sub(".*Score = ([0-9]+\\.?[0-9]*) bits.*", "\\1", Score_Line)),
                        Raw_Score = as.numeric(sub(".*Score = [0-9]+\\.?[0-9]* bits \\(([0-9]+)\\).*", "\\1", Score_Line)),
                        Expect = as.numeric(sub(".*Expect = ([^,]+).*", "\\1", Score_Line)),
                        Identities_Ratio = sub(".*Identities = ([0-9]+/[0-9]+).*", "\\1", Identities_Line),
                        Identities_Percent = as.numeric(sub(".*Identities = [0-9]+/[0-9]+ \\(([^%]+)%\\).*", "\\1", Identities_Line)),
                        Gaps_Ratio = sub(".*Gaps = ([0-9]+/[0-9]+).*", "\\1", Identities_Line),
                        Gaps_Percent = as.numeric(sub(".*Gaps = [0-9]+/[0-9]+ \\(([^%]+)%\\).*", "\\1", Identities_Line)),
                        Hit_Name = Greater_Than_Line
                    ) %>%
                    select(Hit_Name, Score, Raw_Score, Expect, Identities_Ratio, Identities_Percent, Gaps_Ratio, Gaps_Percent)
                parsed_results$Hit_Name <- gsub(">", "", parsed_results$Hit_Name)
                
                # Display the resulting parsed data frame
                return(as_tibble(parsed_results))
            }


        #### buildDomainLibrary


            #' Build a Hidden Markov Model (HMM) database from a set of sequence alignments
            #'
            #' This function takes a list of sequence alignment files and constructs a HMM database.
            #' It uses the HMMER suite to build individual HMMs from each alignment and then combines them
            #' into a single HMM database file. The resulting database is then indexed for efficient searching.
            #'
            #' @param alignment_in_paths A character vector containing the file paths to the sequence alignment files.
            #'        Each file should be in a format compatible with the HMMER 'hmmbuild' tool.
            #' @param domain_library_out_path A character string specifying the file path where the resulting HMM database
            #'        will be saved. This file will contain all the HMMs built from the input alignments.
            #' @return This function does not return a value. It creates and indexes an HMM database file at the specified path.
            #' @details The function first removes any existing files at the output path to prevent appending to old data.
            #'          It then iterates over each alignment file, builds a temporary HMM file using 'hmmbuild', and appends
            #'          it to the output database file. After processing all alignments, the temporary file is deleted, and
            #'          the 'hmmpress' command is used to index the database for faster access during searches.
            #' @examples
            #' \dontrun{
            #' \dontrun{
            #' buildDomainLibrary(
            #'     alignment_in_paths = c("alignment1.sto", "alignment2.sto"),
            #'     domain_library_out_path = "output_hmm_database.hmm"
            #' )
            #' }
            #' }
            #' @export
            buildDomainLibrary <- function(
                alignment_in_paths,
                domain_library_out_path
            ) {

                # Remove any existing files at the output path to avoid appending to old data
                system(paste0("rm ", domain_library_out_path, "*"))

                # Iterate over each alignment file to build and append HMMs
                for ( i in 1:length(alignment_in_paths) ) {
                    # Build a temporary HMM from the current alignment
                    system(paste0(
                        "hmmbuild temp.hmm ", alignment_in_paths[i]
                    ), intern = TRUE)
                    # Append the temporary HMM to the output database
                    system(paste0("cat temp.hmm >> ", domain_library_out_path))
                }
                
                # Remove the temporary HMM file
                system("rm temp.hmm")
                
                # Index the HMM database for efficient searching
                system(paste0(
                    "hmmpress ", domain_library_out_path
                ))
            }


        #### predictDomains

            #' Predict PFAM domains from sequences in a FASTA file using HMMER
            #'
            #' This function utilizes the HMMER software suite to scan a given multifasta file for PFAM domains
            #' using a specified HMM domain library. The results are parsed and returned as a data frame.
            #'
            #' @param fasta_in_path A character string specifying the path to the multifasta file containing the sequences to be analyzed.
            #'        The file should be in a format compatible with the HMMER 'hmmscan' tool.
            #' @param domain_library_in_path A character string specifying the path to the HMM file that contains a library of domains.
            #'        By default, it uses the Pfam-A.hmm library located at "/project_data/shared/general_lab_resources/hmmer/Pfam-A.hmm".
            #' @return A data frame containing the parsed results of the domain predictions. The data frame includes columns such as
            #'         domain name, domain accession, domain length, query name, query accession, sequence e-value, and more.
            #' @details The function first removes any existing temporary output file to ensure clean results. It then runs the 'hmmscan'
            #'          command to predict domains in the input sequences. The output is parsed into a structured data frame using
            #'          predefined column types. The function handles the parsing of domain descriptions separately and appends them
            #'          to the resulting data frame.
            #' @examples
            #' \dontrun{
            #' \dontrun{
            #' predictDomains(
            #'     fasta_in_path = "example_sequences.fasta",
            #'     domain_library_in_path = "custom_domain_library.hmm"
            #' )
            #' }
            #' }
            #' @export
            predictDomains <- function(
                fasta_in_path,
                domain_library_in_path = "/project_data/shared/general_lab_resources/hmmer/Pfam-A.hmm"
            ) {
                # Define the path for the temporary output file
                temp_file_path <- paste0(getwd(), "/predict_domains.txt")
                # Remove the temporary file if it already exists
                if (file.exists(temp_file_path)) {file.remove(temp_file_path)}

                ## Run HMM3SCAN to predict domains
                system(paste0(
                    "hmmscan --domtblout ", temp_file_path,
                    " ", domain_library_in_path, " ",
                    fasta_in_path
                ))

                ## Parse the output using predefined column types
                col_types <- readr::cols(
                    domain_name         = readr::col_character(),
                    domain_accession    = readr::col_character(),
                    domain_len          = readr::col_integer(),
                    query_name          = readr::col_character(),
                    query_accession     = readr::col_character(),
                    qlen                = readr::col_integer(),
                    sequence_evalue     = readr::col_double(),
                    sequence_score      = readr::col_double(),
                    sequence_bias       = readr::col_double(),
                    domain_N            = readr::col_integer(),
                    domain_of           = readr::col_integer(),
                    domain_cevalue      = readr::col_double(),
                    domain_ievalue      = readr::col_double(),
                    domain_score        = readr::col_double(),
                    domain_bias         = readr::col_double(),
                    hmm_from            = readr::col_integer(),
                    hmm_to              = readr::col_integer(),
                    ali_from            = readr::col_integer(),
                    ali_to              = readr::col_integer(),
                    env_from            = readr::col_integer(),
                    env_to              = readr::col_integer(),
                    acc                 = readr::col_double()
                )

                N <- length(col_types$cols)

                # Read the lines from the temporary output file
                lines <- readr::read_lines(temp_file_path, lazy=FALSE, progress=FALSE)

                # Parse the table from the output lines
                table <- sub(
                    pattern = sprintf("(%s).*", paste0(rep('\\S+', N), collapse=" +")),
                    replacement = '\\1',
                    x=lines,
                    perl = TRUE
                ) %>%
                gsub(pattern="  *", replacement="\t") %>%
                paste0(collapse="\n") %>%
                readr::read_tsv(
                    col_names=names(col_types$cols),
                    comment='#',
                    na='-',
                    col_types = col_types,
                    lazy=FALSE,
                    progress=FALSE
                )

                # Extract and append domain descriptions
                descriptions <- lines[!grepl("^#", lines, perl=TRUE)] %>%
                sub(
                    pattern = sprintf("%s *(.*)", paste0(rep('\\S+', N), collapse=" +")),
                    replacement = '\\1',
                    perl = TRUE
                )

                table$description <- descriptions[!grepl(" *#", descriptions, perl=TRUE)]

                return(table)
            }


        #### hmmAlign

            #' Align sequences to a Hidden Markov Model (HMM)
            #'
            #' This function aligns sequences from a FASTA file to a specified domain in an HMM file.
            #' It extracts the relevant positions based on the RF (Reference) sequence and returns
            #' the aligned sequences as an AAStringSet object.
            #'
            #' @param fasta_in_path The path to the input FASTA file containing sequences to be aligned.
            #' @param hmm_in_path The path to the HMM file containing the domain model.
            #' @param domain_name The name of the domain to fetch from the HMM file for alignment.
            #' @return An AAStringSet object containing the aligned sequences.
            #' @examples
            #' \dontrun{
            #' # Example usage:
            #' # aligned_sequences <- hmmAlign("sequences.fasta", "domains.hmm", "domain_name")
            #' # print(aligned_sequences)
            #' }
            #' @export
            hmmAlign <- function(fasta_in_path, hmm_in_path, domain_name) {
                .check_bioc(c("Biostrings"), "hmmAlign")
    
                # Fetch HMM and perform alignment
                    hmm_file <- tempfile()
                    alignment_file <- tempfile()
                    system(sprintf("hmmfetch %s %s > %s", hmm_in_path, domain_name, hmm_file))
                    system(sprintf("hmmalign %s %s > %s", hmm_file, fasta_in_path, alignment_file))

                # Read alignment file, Extract RF sequence and determine positions to keep
                    lines <- readLines(alignment_file)
                    rf_sequence <- paste0(sub("^#=GC RF\\s+", "", grep("^#=GC RF", lines, value = TRUE)), collapse = "")
                    positions_to_keep <- which(strsplit(rf_sequence, "")[[1]] == "x")
                
                # Extract and concatenate sequences per header
                    sequence_lines <- lines[!grepl("^#", lines) & grepl("^[A-Za-z0-9_]", lines)]
                    parts <- strsplit(sequence_lines, "\\s+")
                    headers <- sapply(parts, `[[`, 1)
                    sequences <- sapply(parts, function(x) x[length(x)])
                    sequences_per_header <- tapply(sequences, headers, paste0, collapse = "")
                
                # Crop sequences based on RF and assemble into AAStringSet
                    cropped_sequences <- sapply(sequences_per_header, function(seq) {
                        paste0(strsplit(seq, "")[[1]][positions_to_keep], collapse = "")
                    })
                    sequences <- Biostrings::AAStringSet(cropped_sequences)
                    names(sequences) <- names(cropped_sequences)
                    return(sequences)
            }


        #### extractORFs


            #' Extract open reading frames (ORFs) from multifasta files
            #'
            #' This function identifies and extracts open reading frames (ORFs) from sequences in a multifasta file.
            #' It can optionally write the ORFs to a new fasta file and overwrite the input file with the ORF file.
            #' The function can also search for alternative ORFs in the complementary, reverse, and reverse-complementary sequences.
            #' Note: Searching for ORFs in the reverse direction leads to issues during codon alignment...
            #'
            #' @param file_in_path The path to the multifasta file to analyze.
            #' @param write_out_ORFs Logical. If TRUE, writes out a new fasta file containing the ORFs.
            #' @param overwrite Logical. If TRUE, optionally overwrites the input file with the ORF file.
            #' @param alternative_ORFs Logical. If TRUE, searches for alternative ORFs in the complementary, reverse, and reverse-complementary sequences.
            #' @return A data frame containing the accession, start codon position, stop codon position, direction, and length of each ORF.
            #' @examples
            #' \dontrun{
            #' # Example usage:
            #' # Extract ORFs from a multifasta file and write them to a new file
            #' # orfs <- extractORFs("sequences.fasta", write_out_ORFs = TRUE, overwrite = FALSE, alternative_ORFs = TRUE)
            #' # print(orfs)
            #' }
            #' @export
            extractORFs <- function(
                file_in_path,
                write_out_ORFs = FALSE,
                overwrite = FALSE,
                alternative_ORFs = FALSE
            ) {
                .check_bioc(c("Biostrings"), "extractORFs")
                rlang::check_installed(c("seqinr"), reason = "to use extractORFs()")

                fasta <- seqinr::read.fasta(file_in_path)
                ORFs <- list()
                ORFs_to_write_out <- Biostrings::DNAStringSet()

                if ( write_out_ORFs == TRUE ) {
                    if (file.exists(paste(file_in_path, "_orfs", sep = ""))) { file.remove(paste(file_in_path, "_orfs", sep = "")) }
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

                    if (alternative_ORFs) {

                        # Search for ORFs in the compliment mRNA
                            mRNA_comp <- chartr("ATGC","TACG", mRNA)
                            mRNA_comp <- chartr("atcg","tacg", mRNA_comp)
                            data <- seqinr::translate(mRNA_comp, frame = 0, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_comp, frame = 1, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_comp, frame = 2, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "compliment"
                                                                                    )
                                                            )
                                    }
                                }

                        # Search for ORFs in the reverse mRNA
                            mRNA_rev <- rev(mRNA)
                            data <- seqinr::translate(mRNA_rev, frame = 0, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev, frame = 1, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev, frame = 2, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse"
                                                                                    )
                                                            )
                                    }
                                }

                        # Search for ORFs in the reverse compliment mRNA
                            mRNA_rev_comp <- chartr("ATGC","TACG", mRNA)
                            mRNA_rev_comp <- chartr("atcg","tacg", mRNA_rev_comp)
                            mRNA_rev_comp <- rev(mRNA_rev_comp)
                            data <- seqinr::translate(mRNA_rev_comp, frame = 0, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse_compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev_comp, frame = 1, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse_compliment"
                                                                                    )
                                                            )
                                    }
                                }
                            data <- seqinr::translate(mRNA_rev_comp, frame = 2, sens = "F")
                                S <- grep("\\*", data) # find all stop codons
                                M <- grep("M", data) # find all start codons
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    M <- M[!M > max(S)] # remove all start codons after the last stop codon
                                }
                                if (length(S)>0 & length(M)>0) { # if there are stop codons, extract orfs
                                    for (i in 1:length(M)) { # extract all cds in this frame and put them in a list
                                        F1_amin <- c(M[i],cbind(S, S > M[i])[,1][cbind(S, S > M[i])[,2]==1][1])
                                        orf_coordinates <-  rbind(orf_coordinates,   data.frame(
                                                                                        start_codon = (M[i]*3-2), 
                                                                                        stop_codon = ((M[i]*3-2)+((F1_amin[2] - F1_amin[1])+1)*3-1),
                                                                                        direction = "reverse_compliment"
                                                                                    )
                                                            )
                                    }
                                }

                    }

                    ## Process discovered ORFs

                        if ( dim(orf_coordinates)[1] == 0 ) {
                            ORFs[[j]] <-    data.frame(
                                                accession = names(fasta[j]),
                                                start_codon = 0, 
                                                stop_codon = 0,
                                                direction = "forward",
                                                orf_length = 0
                                            )
                        } else {
                            rownames(orf_coordinates) <- NULL # Reset rownames because of previous step
                            orf_coordinates$orf_length <- (orf_coordinates$stop_codon - orf_coordinates$start_codon + 1)
                            orf_coordinates <- orf_coordinates[order(orf_coordinates$orf_length, decreasing = TRUE),]
                            ORFs[[j]] <- cbind(data.frame(accession = names(fasta[j])), orf_coordinates)
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
                        Biostrings::writeXStringSet(ORFs_to_write_out, file = file_in_path)
                        # seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "w") # overwrite transcript file
                    } else {
                        Biostrings::writeXStringSet(ORFs_to_write_out, file = paste(file_in_path, "_orfs", sep = ""))
                        # seqinr::write.fasta(orf, names = attr(fasta[j], "name"), file = paste(file, "_orfs", sep = ""), open = "a") # append ORF to the list        
                    }
                }

                ORFs <- do.call(rbind, ORFs)
                rownames(ORFs) <- NULL
                return(ORFs)
            }


        #### polyBlast


            #' BLAST search local sequence collections
            #'
            #' Search locally stored sequence collections for query sequences. Download Blast+ here: https://www.ncbi.nlm.nih.gov/books/NBK279671/
            #' @param named_subjects_list A list of paths to the transcriptomes that should be searched, named by taxonomic identifier (e.g. species names)
            #' @param query_in_path Path to a fasta file containing the query or queries
            #' @param sequences_of_interest_directory_path Path to a directory where blast hits should be written out as fasta files
            #' @param blast_module_directory_path Path to directory containing the BLAST+ module (perhaps something like "/usr/local/ncbi/blast/bin/")
            #' @param blast_mode One of "blastn" or "dc-megablast". "blastn" is a traditional BLASTN requiring an exact match of 11. "dc-megablast" is a discontiguous megablast used to find more distant (e.g., interspecies) sequences.
            #' @param e_value_cutoff e-value cutoff to apply to results. Defaults to 1.
            #' @param word_size Word size for the BLAST search. Defaults to 4.
            #' @param queries_in_output Should the queries be included in the fasta output?
            #' @param monolist_out_path Path to where the output monolist should be written
            #' @return A data frame (monolist) containing the BLAST search results with columns for accession, Genus, species, annotation, length, longestORF, length_aligned_with_query, percent_identity, e_value, bitscore, subset_all, word_size, query, query_length, and query_longestORF.
            #' @details This function performs a BLAST search on a set of query sequences against a list of subject transcriptomes. It supports both nucleotide and protein BLAST searches, and can handle both Unix and Windows operating systems. The function first translates the transcriptomes if necessary, builds BLAST databases, and then performs the BLAST search. The results are processed to remove duplicates and are written to a monolist, which is also saved to a specified output path.
            #' @examples
            #' \dontrun{
            #' named_subjects_list <- list("species1" = "path/to/species1.fasta", "species2" = "path/to/species2.fasta")
            #' query_in_path <- "path/to/query.fasta"
            #' sequences_of_interest_directory_path <- "path/to/output/directory"
            #' blast_module_directory_path <- "/usr/local/ncbi/blast/bin/"
            #' blast_mode <- "blastn"
            #' e_value_cutoff <- 1
            #' word_size <- 4
            #' queries_in_output <- TRUE
            #' monolist_out_path <- "path/to/monolist_output"
            #' polyBlast(named_subjects_list, query_in_path, sequences_of_interest_directory_path, blast_module_directory_path, blast_mode, e_value_cutoff, word_size, queries_in_output, monolist_out_path)
            #' }
            #' @export
            polyBlast <- function(
                                    named_subjects_list,
                                    query_in_path,
                                    sequences_of_interest_directory_path,
                                    blast_module_directory_path,
                                    blast_mode = c("nnblastn", "ntblastp", "pnblastp"), 
                                    e_value_cutoff = 1,
                                    word_size = 4,
                                    queries_in_output = TRUE,
                                    monolist_out_path
                                ) {
                .check_bioc(c("Biostrings"), "polyBlast")
                rlang::check_installed(c("seqinr"), reason = "to use polyBlast()")

                ### Make sure transcriptomes object is character and has names
                    transcriptomes <- named_subjects_list
                    if( length(names(transcriptomes)) != length(transcriptomes) ) {
                        stop("Please provide a set of named paths to the argument `transcriptomes`")
                    }
                    
                    names <- names(transcriptomes)
                    transcriptomes <- as.character(transcriptomes)
                    names(transcriptomes) <- names

                ### If *tblast*, then translate all transcriptomes and rename the translations with "_trans" suffix

                    if ( substr(blast_mode, 2, 2) == "t" ) {
                        cat("Translating transcriptomes...\n\n")
                        for ( transcriptome in 1:length(transcriptomes) ) {
                            writeFasta(
                                XStringSet = Biostrings::translate(Biostrings::readDNAStringSet(transcriptomes[transcriptome]), if.fuzzy.codon = "solve"),
                                fasta_out_path = paste(transcriptomes[transcriptome], "trans", sep = "_"),
                                fasta_type = "AA"
                            )
                            transcriptomes[transcriptome] <- paste(transcriptomes[transcriptome], "trans", sep = "_")
                        }
                    } else if ( substr(blast_mode, 2, 2) == "n" ) {

                    } else {
                        stop("Please specify t or n to indicate whether the subjects should be translated.")
                    }

                ### Build blast database(s)

                    if ( substr(blast_mode, 1, 1) %in% c("n") ) {
                        dbtype = "nucl"
                    } else if ( substr(blast_mode, 1, 1) %in% c("p") ) {
                        dbtype = "prot"
                    } else {
                        stop("Please specify a subject type of n or p")
                    }

                    if ( .Platform$OS == "unix" ) {

                        for (transcriptome in 1:length(transcriptomes)) {
                            cat(paste0("\nSpecies: ", names(transcriptomes)[transcriptome]))
                            system(
                                paste(
                                    blast_module_directory_path,
                                    "makeblastdb -in ",
                                    transcriptomes[transcriptome],
                                    " -dbtype ",
                                    dbtype,
                                    sep = ""
                                ) 
                            )
                        }

                    } else if ( .Platform$OS == "windows" ) {

                        for (transcriptome in 1:length(transcriptomes)) {
                            cat(paste0("\nSpecies: ", names(transcriptomes)[transcriptome]))
                            shell(
                                paste(
                                    blast_module_directory_path,
                                    "makeblastdb -in ",
                                    transcriptomes[transcriptome],
                                    " -dbtype ",
                                    dbtype,
                                    sep = ""
                                )
                            )
                        }
                    }
                    
                ### Start BLAST process

                    query_seqs <- Biostrings::readBStringSet(filepath = query_in_path, format = "fasta")
                    if ( substr(blast_mode, 8, 8) == "p") { blast_type <- "blastp" }
                    if ( substr(blast_mode, 8, 8) == "n") { blast_type <- "blastn" }

                        ## Loop over each member of the query and use it to blast each transcriptome
                            
                            monolist <- data.frame()
                            cat("\n\n")
                            
                            for ( query_seq in 1:length(query_seqs) ) { # query_seq = 1

                                ## Write individual files for each member of the query
                                    
                                    Biostrings::writeXStringSet(
                                        query_seqs[query_seq], 
                                        filepath = paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""), 
                                        append = FALSE
                                    )

                                ## Loop over the transcriptomes, run the blast on each, add hits to monolist
                                    
                                    for (transcriptome in 1:length(transcriptomes)) { # transcriptome = 1

                                        ## Run BLASTs on unix system

                                            if ( .Platform$OS == "unix") {

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt '6 sallacc'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_length", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt '6 length'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_pident", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt '6 pident'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_evalue", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt '6 evalue'",
                                                        sep = ""
                                                    )
                                                )

                                                system(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt '6 bitscore'",
                                                        sep = ""
                                                    )
                                                )
                                            }

                                        ## Run BLASTs on windows system

                                            if ( .Platform$OS == "windows") {

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt \"6 sallacc\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_length", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt \"6 length\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_pident", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt \"6 pident\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_evalue", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt \"6 evalue\"",
                                                        sep = ""
                                                    )
                                                )

                                                shell(
                                                    paste(
                                                        blast_module_directory_path,
                                                        blast_type,
                                                        " -task ",
                                                        blast_type,
                                                        " -db ", 
                                                        transcriptomes[transcriptome],
                                                        " -query ",
                                                        paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""),
                                                        " -out ",
                                                        paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""),
                                                        " -evalue ",
                                                        e_value_cutoff,
                                                        " -word_size ",
                                                        word_size,
                                                        " -outfmt \"6 bitscore\"",
                                                        sep = ""
                                                    )
                                                )
                                            }

                                        ## Extract BLAST hits from transcriptome, add them to the monolist, write them to individual files

                                            # Read in whole transcriptome
                                                if ( substr(blast_mode, 8, 8) == "p" ) {
                                                    temp_seqs <- Biostrings::readAAStringSet( 
                                                        filepath = as.character(transcriptomes)[transcriptome], 
                                                        format = "fasta"
                                                    )
                                                } else if ( substr(blast_mode, 8, 8) == "n" ) {
                                                    temp_seqs <- Biostrings::readDNAStringSet( 
                                                        filepath = as.character(gsub("_trans$", "", transcriptomes)[transcriptome]), 
                                                        format = "fasta"
                                                    )
                                                }
                                        
                                            # Attempt to read in hits list, subset transcriptome according to that list, add hits to monolist
                                                cat(paste0("Query ", query_seq, ": "))
                                                if ( inherits( try( read.table(file = paste(transcriptomes[transcriptome], ".out", sep = "")), silent = TRUE), "try-error") == TRUE ) { # Skips over empty files
                                                    cat(paste("No BLAST hits found for ", query_seqs[query_seq]@ranges@NAMES, " in ", names(transcriptomes)[transcriptome], "\n", sep = ""))
                                                } else {

                                                    ## Remove duplicate hits from all output files
                                                        temp_hits <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out", sep = ""))[,1])
                                                        temp_hits_length <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_length", sep = ""))[,1])
                                                        temp_hits_pident <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_pident", sep = ""))[,1])
                                                        temp_hits_evalue <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_evalue", sep = ""))[,1])
                                                        temp_hits_bitscore <- as.character(read.table(file = paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""))[,1])
                                                        
                                                        duplicate_indeces <- duplicated(temp_hits)

                                                        temp_hits <- temp_hits[!duplicate_indeces]
                                                        temp_hits_length <- temp_hits_length[!duplicate_indeces]
                                                        temp_hits_pident <- temp_hits_pident[!duplicate_indeces]
                                                        temp_hits_evalue <- temp_hits_evalue[!duplicate_indeces]
                                                        temp_hits_bitscore <- temp_hits_bitscore[!duplicate_indeces]
                                                        
                                                        writeMonolist(temp_hits, paste(transcriptomes[transcriptome], ".out", sep = ""))
                                                        writeMonolist(temp_hits_length, paste(transcriptomes[transcriptome], ".out_length", sep = ""))
                                                        writeMonolist(temp_hits_pident, paste(transcriptomes[transcriptome], ".out_pident", sep = ""))
                                                        writeMonolist(temp_hits_evalue, paste(transcriptomes[transcriptome], ".out_evalue", sep = ""))
                                                        writeMonolist(temp_hits_bitscore, paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""))

                                                    temp_hits <- readMonolist(paste(transcriptomes[transcriptome], ".out", sep = ""))
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
                                                    temp_seqs@ranges@NAMES <- gsub(" .*", "", temp_seqs@ranges@NAMES)

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
                                                                    longestORF = longest_ORFs,
                                                                    length_aligned_with_query = readMonolist(paste(transcriptomes[transcriptome], ".out_length", sep = ""))[,1],
                                                                    percent_identity = readMonolist(paste(transcriptomes[transcriptome], ".out_pident", sep = ""))[,1],
                                                                    e_value = readMonolist(paste(transcriptomes[transcriptome], ".out_evalue", sep = ""))[,1],
                                                                    bitscore = readMonolist(paste(transcriptomes[transcriptome], ".out_bitscore", sep = ""))[,1],
                                                                    subset_all = TRUE,
                                                                    word_size = word_size,
                                                                    query = query_seqs@ranges@NAMES[query_seq],
                                                                    query_length = query_seqs@ranges@width[query_seq],
                                                                    query_longestORF = extractORFs(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))$orf_length[1]
                                                                ))
                                                        }
                                                }
                                    }

                                ## Optionally add the queries to the output

                                    if (queries_in_output == TRUE) {

                                        if (blast_type == "blastp") {
                                                longestORF <- 3*query_seqs@ranges@width[query_seq] 
                                            } else {
                                                longestORF <- extractORFs(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))$orf_length[1]
                                            }
                                        
                                        monolist <- rbind(monolist, data.frame(
                                            accession = query_seqs@ranges@NAMES[query_seq],
                                            Genus = "Query",
                                            species = "query",
                                            annotation = "query",
                                            length = query_seqs@ranges@width[query_seq],
                                            longestORF = longestORF,
                                            length_aligned_with_query = query_seqs@ranges@width[query_seq],
                                            percent_identity = 100,
                                            e_value = 0,
                                            bitscore = 1000,
                                            subset_all = TRUE,
                                            word_size = word_size,
                                            query = query_seqs@ranges@NAMES[query_seq],
                                            query_length = query_seqs@ranges@width[query_seq],
                                            query_longestORF = extractORFs(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))$orf_length[1]
                                        ))

                                        Biostrings::writeXStringSet(
                                            query_seqs[query_seq], 
                                            filepath = paste(sequences_of_interest_directory_path, query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""), 
                                            append = FALSE
                                        )
                                    }

                                ## Remove individual query files
                                    if (file.exists(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = "")) ) {
                                        file.remove(paste(query_in_path, "_", query_seqs@ranges@NAMES[query_seq], ".fa", sep = ""))
                                    }
                            }

                        ## Write out the monolist

                            monolist <- unique(monolist)
                            rownames(monolist) <- NULL
                            if (file.exists(monolist_out_path)) {file.remove(monolist_out_path)}
                            writeMonolist(monolist = monolist, monolist_out_path = monolist_out_path)
               
                    cat("\nDone!\n\n")
            }


        #### codonAlignFromAA

            #' Back-translate an amino-acid alignment onto its codons
            #'
            #' Given an aligned set of amino-acid sequences and the *ungapped* coding
            #' nucleotide sequences they were translated from, produce the corresponding
            #' codon (nucleotide) alignment: each non-gap residue becomes its three
            #' nucleotides and each gap ("-") becomes "---". This is the in-frame core of
            #' PAL2NAL and replaces the former orthologr::codon_aln() dependency. It assumes
            #' each peptide is an exact, in-frame translation of its nucleotide sequence
            #' (equal codon counts, no frameshift) -- which holds inside alignSequences(),
            #' where the peptides are translated from the same ORFs. If that correspondence
            #' is broken it stops with a clear message rather than misaligning silently.
            #'
            #' @param aa_alignment A named character vector (or object coercible with
            #'   as.character, e.g. an AAStringSet) of aligned, gapped amino-acid sequences.
            #' @param nucl_seqs A named DNAStringSet (or named character vector) of the
            #'   ungapped coding nucleotides, one per aligned sequence, matched by name.
            #' @return A named character vector of codon-aligned nucleotide sequences.
            #' @export
            codonAlignFromAA <- function(aa_alignment, nucl_seqs) {

                ## as.character() drops names on a base vector but keeps them on an
                ## XStringSet, so restore them from the original where needed.
                aa <- as.character(aa_alignment)
                if (is.null(names(aa))) names(aa) <- names(aa_alignment)
                nt_all <- as.character(nucl_seqs)
                if (is.null(names(nt_all))) names(nt_all) <- names(nucl_seqs)

                if (is.null(names(aa)) || is.null(names(nt_all))) {
                    stop("codonAlignFromAA: both aa_alignment and nucl_seqs must be named.", call. = FALSE)
                }

                out <- setNames(character(length(aa)), names(aa))

                for (nm in names(aa)) {

                    if (!nm %in% names(nt_all)) {
                        stop(paste0(
                            "codonAlignFromAA: no nucleotide sequence named '", nm,
                            "' to match the aligned peptide. Peptide and nucleotide names must correspond."
                        ), call. = FALSE)
                    }

                    residues <- strsplit(aa[[nm]], "")[[1]]
                    nt <- nt_all[[nm]]
                    n_res <- sum(residues != "-")

                    ## Guard: the peptide must be an exact, in-frame translation of the nucleotides.
                    if (nchar(nt) != n_res * 3L) {
                        stop(paste0(
                            "codonAlignFromAA: inconsistency between the peptide and nucleotide sequences for '", nm, "'.\n",
                            "  The aligned peptide has ", n_res, " residues (excluding gaps), needing ",
                            n_res * 3L, " nucleotides,\n",
                            "  but the coding sequence has ", nchar(nt), " nucleotides.\n",
                            "  Codon alignment expects each peptide to be an exact, in-frame translation of its\n",
                            "  nucleotide sequence (no frameshifts, no extra leading/trailing bases). Common causes:\n",
                            "  duplicate accessions in the monolist, an off-by-one ORF, or a stop codon present on\n",
                            "  one side only."
                        ), call. = FALSE)
                    }

                    codons <- character(length(residues))
                    pos <- 1L
                    for (j in seq_along(residues)) {
                        if (residues[j] == "-") {
                            codons[j] <- "---"
                        } else {
                            codons[j] <- substr(nt, pos, pos + 2L)
                            pos <- pos + 3L
                        }
                    }
                    out[nm] <- paste0(codons, collapse = "")
                }

                out
            }


        #### alignSequences

            #' Align sequences in a seqlist
            #'
            #' This function aligns sequences from a given monolist based on the specified mode and input sequence type.
            #' It supports nucleotide, amino acid, codon, and fragment alignments. The aligned sequences are written to the specified output path.
            #' ERROR: inconsistency between the following pep and nuc seqs" - usually means there are duplicate accessions numbers in the input monolist
            #'
            #' @param monolist A data frame containing the sequences to be aligned. The first column should be "accession".
            #' @param subset A logical vector indicating which sequences in the monolist should be included in the alignment.
            #' @param alignment_out_path A character string specifying the path to the file where the alignment should be written.
            #' @param sequences_of_interest_directory_path A character string specifying the path to a directory where blast hits should be written out as fasta files.
            #' @param input_sequence_type A character string specifying the type of input sequences. Options are "nucl" for nucleotide sequences or "amin" for amino acid sequences.
            #' @param mode A character string specifying the alignment mode. Options are "nucl_align" for nucleotide alignment, "amin_align" for amino acid alignment, "codon_align" for codon alignment, or "fragment_align" for fragment alignment.
            #' @param base_fragment A character string specifying the base fragment for fragment alignment. Default is NULL.
            #' @return This function does not return a value. It writes the aligned sequences to the specified output path.
            #' @examples
            #' \dontrun{
            #' # Example usage:
            #' # Align nucleotide sequences
            #' alignSequences(
            #'     monolist = my_monolist,
            #'     subset = my_subset,
            #'     alignment_out_path = "output_alignment.fa",
            #'     sequences_of_interest_directory_path = "sequences/",
            #'     input_sequence_type = "nucl",
            #'     mode = "nucl_align"
            #' )
            #' }
            #' @export
            alignSequences <-   function(
                        monolist, 
                        subset, 
                        alignment_out_path,
                        sequences_of_interest_directory_path,
                        input_sequence_type = c("nucl", "amin"), 
                        mode = c("nucl_align", "amin_align", "codon_align", "fragment_align"),
                        base_fragment = NULL
                    ){  
                .check_bioc(c("Biostrings", "msa"), "alignSequences")
                rlang::check_installed(c("seqinr"), reason = "to use alignSequences()")

                ## Check directory_paths, Get appropriate subset of the monolist
                    sequences_of_interest_directory_path <- OsDirectoryPathCorrect(sequences_of_interest_directory_path)
                    monolist_subset <- monolist[unlist(monolist[,as.character(colnames(monolist)) == as.character(subset)]),]

                ## If starting with nucleotide sequences

                    if ( input_sequence_type == "nucl" ) {

                        ## Remove existing version of this alignment and its files
                            if (file.exists(alignment_out_path)) {file.remove(alignment_out_path)}

                        ## Import the subset's sequences, then write out *_nucl_seqs.fa
                            nucl_seqs_set <- Biostrings::DNAStringSet()
                            for (i in 1:dim(monolist_subset)[1]) {
                                nucl_seqs_set <- c(nucl_seqs_set, Biostrings::readDNAStringSet(paste(sequences_of_interest_directory_path, as.character(monolist_subset$accession[i]), ".fa", sep = "")))
                            }

                        ## Run plain nucleotide alignment, if requested
                            if ( mode == "nucl_align") {
                                nucl_seqs_set_aligned <- msa::msa(nucl_seqs_set, order = c("input"))
                                msa <- msa::msaConvert(nucl_seqs_set_aligned, type = "seqinr::alignment")
                                for (i in 1:dim(as.data.frame(msa$seq))[1]) {
                                    msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))
                                }
                                seqinr::write.fasta(
                                    as.list(msa$seq),
                                    as.list(msa$nam),
                                    file.out = alignment_out_path
                                )
                            }

                        ## Run codon alignment, if requested
                            if ( mode == "codon_align") {

                                ## Extract ORFs, read in ORFs, translate them, then align them

                                    Biostrings::writeXStringSet(nucl_seqs_set, filepath = paste0(alignment_out_path, "_unaligned"))
                                    extractORFs(file_in_path = paste0(alignment_out_path, "_unaligned"), write_out_ORFs = TRUE)
                                    if (file.exists(paste0(alignment_out_path, "_unaligned"))) {file.remove(paste0(alignment_out_path, "_unaligned"))}
                                
                                    ## The ORF nucleotides are what the peptides are translated from, so they
                                    ## are the correct sequences to codon-align against (matched by name).
                                    orf_seqs_set <- Biostrings::readDNAStringSet(paste0(alignment_out_path, "_unaligned_orfs"))
                                    if (file.exists(paste0(alignment_out_path, "_unaligned_orfs"))) {file.remove(paste0(alignment_out_path, "_unaligned_orfs"))}

                                    amin_seqs_set <- Biostrings::translate(orf_seqs_set, if.fuzzy.codon = "solve")

                                    amin_seqs_set_aligned <- msa::msa(amin_seqs_set, order = c("input"))
                                    msa <- msa::msaConvert(amin_seqs_set_aligned, type = "seqinr::alignment")
                                    for (i in 1:dim(as.data.frame(msa$seq))[1]) {
                                        msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))
                                    }
                                    seqinr::write.fasta(
                                        as.list(msa$seq),
                                        as.list(msa$nam),
                                        file.out = paste0(alignment_out_path, "_amin_seqs_aligned.fa")
                                    )

                                ## Back-translate the peptide alignment onto codons (pure-R, replaces the
                                ## former orthologr::codon_aln() / PAL2NAL dependency), then write out.
                                    aa_aligned <- setNames(as.character(msa$seq), msa$nam)
                                    nucl_seqs_codon_aligned <- codonAlignFromAA(aa_aligned, orf_seqs_set)
                                    seqinr::write.fasta(
                                        as.list(unname(nucl_seqs_codon_aligned)),
                                        as.list(names(nucl_seqs_codon_aligned)),
                                        file.out = alignment_out_path
                                    )
                            }

                        ## Run fragment alignment, if requested
                            if ( mode == "fragment_align") {

                                ## Remove existing fragment alignment
            #                         if (file.exists(paste(alignment_directory_path, as.character(subset), "_fragment_seqs_aligned.fa", sep = ""))) {
            #                         file.remove(paste(alignment_directory_path, as.character(subset), "_fragment_seqs_aligned.fa", sep = ""))}

            #                     ## Extract ORFs
            #                         extractORFs(file_in_path = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), write_out_ORFs = TRUE)

            #                     ## Translate the nucleotide sequences and write out *_amin_seqs.fa
            #                         nucl_seqs_set <- Biostrings::readDNAStringSet(paste(alignment_directory_path, subset, "_nucl_seqs.fa_orfs", sep = ""))
            #                         amin_seqs_set <- Biostrings::translate(nucl_seqs_set, if.fuzzy.codon = "solve")
            #                         if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))) {
            #                             file.remove(paste(alignment_directory_path, as.character(subset), "_", "amin_seqs.fa", sep = ""))}
            #                         Biostrings::writeXStringSet(amin_seqs_set, filepath = paste(alignment_directory_path, subset, "_amin_seqs.fa", sep = ""))

            #                     ## Define base fragment
            #                         base_fragment_seq <- amin_seqs_set[amin_seqs_set@ranges@NAMES == base_fragment]
            #                         fragments_seq_set <- amin_seqs_set[amin_seqs_set@ranges@NAMES != base_fragment]
            #                         current_base_fragment_seq <- base_fragment_seq

            #                     ## Align each fragment on its own with the base fragment
            #                         aligned_fragments <- AAStringSet()
            #                         for ( fragment in 1:length(fragments_seq_set) ) {
            #                             fragment_pair_seq_set <- c(base_fragment_seq, fragments_seq_set[fragment])
            #                             fragment_pair_seq_set_aligned <- msa::msa(fragment_pair_seq_set, order = c("input"))
            #                             fragment_pair_seq_set_aligned <- AAStringSet(fragment_pair_seq_set_aligned)
            #                             current_base_fragment_seq <- fragment_pair_seq_set_aligned[fragment_pair_seq_set_aligned@ranges@NAMES == base_fragment]
            #                             aligned_fragments <- c(aligned_fragments, fragment_pair_seq_set_aligned[fragment_pair_seq_set_aligned@ranges@NAMES != base_fragment])
            #                         }

            #                     ## Perform final alignment and write it out
            #                         fragment_seq_set <- c(current_base_fragment_seq, aligned_fragments)
            #                         fragment_seq_set <- AAStringSet(gsub("-", "Z", fragment_seq_set))
            #                         fragment_seq_set_aligned <- msa::msa(fragment_seq_set, order = c("input"))
            #                         fragment_seq_set_aligned <- AAStringSet(fragment_seq_set_aligned)
            #                         writeXStringSet(fragment_seq_set_aligned, filepath = paste(alignment_directory_path, subset, "_fragment_seqs_aligned.fa", sep = ""))
                            }

                        ## Run amin alignment, if requested
                            if ( mode == "amin_align") {

                                ## Remove existing amin alignment
                                    if (file.exists(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_amin_aligned.fa", sep = ""))) {
                                    file.remove(paste(alignment_directory_path, as.character(subset), "_", "nucl_seqs_amin_aligned.fa", sep = ""))}

                                ## Extract ORFs
                                    extractORFs(file_in_path = paste(alignment_directory_path, subset, "_nucl_seqs.fa", sep = ""), write_out_ORFs = TRUE)

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
                            }
                    }

                ## If starting with amino acid sequences

                    if ( input_sequence_type == "amin" ) {

                        if ( mode == "amin_align" ) {

                            ## Remove existing version of this alignment and its files
                                # if (file.exists(paste0(alignment_out_path, "_unaligned"))) {
                                #     file.remove(paste0(alignment_out_path, "_unaligned")) }    
                                if (file.exists(alignment_out_path)) {
                                    file.remove(alignment_out_path)}

                            ## Import the subset's sequences, then write out *_amin_seqs.fa
                                amin_seqs_set <- Biostrings::AAStringSet()
                                for (i in 1:dim(monolist_subset)[1]) {
                                    amin_seqs_set <- c(
                                        amin_seqs_set,
                                        Biostrings::readAAStringSet(
                                            paste0(
                                                sequences_of_interest_directory_path,
                                                as.character(monolist_subset$accession[i]),
                                                ".fa"
                                            )
                                        )
                                    )
                                }
                                # Biostrings::writeXStringSet(amin_seqs_set, filepath = paste0(alignment_out_path, "_unaligned"))

                            ## Run the alignment
                                amin_seqs_set_aligned <- msa::msa(amin_seqs_set, order = c("input"))
                                msa <- msa::msaConvert(amin_seqs_set_aligned, type = "seqinr::alignment")
                                for (i in 1:dim(as.data.frame(msa$seq))[1]) {
                                    msa$seq[i] <- substr(msa$seq[i],0,nchar(msa$seq[1]))
                                }
                                seqinr::write.fasta(
                                    as.list(msa$seq),
                                    as.list(msa$nam),
                                    file.out = alignment_out_path
                                )
                        }
                    }
            }


        #### analyzeAlignment
        

            #' Analyze a multiple sequence alignment and subset it for phylogeny building
            #'
            #' This function reads a multiple sequence alignment, analyzes the composition of each position, and subsets the alignment based on user-defined thresholds for gap percentage and conservation percentage. It also provides a user interface for visualizing the alignment and the resulting phylogenetic trees.
            #'
            #' @param alignment_in_path The path to the alignment file to analyze, in fasta format.
            #' @param type The type of sequences in the alignment, either "DNA" or "AA".
            #' @param jupyter Logical. If TRUE, the function will run in a Jupyter notebook environment.
            #' @return A Shiny app for visualizing and subsetting the alignment.
            #' @examples
            #' \dontrun{
            #' # Example usage:
            #' # analyzeAlignment("path/to/alignment.fasta", type = "DNA", jupyter = FALSE)
            #' }
            #' @export
            analyzeAlignment <- function(
                alignment_in_path,
                type = c("DNA", "AA"),
                jupyter = FALSE
            ) {
                .check_bioc(c("Biostrings"), "analyzeAlignment")

                ## Prelim stuff

                ## Read in and process alignment
                alignment <- readAlignment(
                    alignment_in_path = alignment_in_path,
                    type = type
                )

                if (type == "DNA") { scaffold_type <- "nucl_alignment" }
                if (type == "AA") { scaffold_type <- "amin_alignment" }

                tree <- fortify(buildTree(
                    scaffold_type = scaffold_type,
                    scaffold_in_path = alignment_in_path
                ))

                alignment$name <- factor(alignment$name, levels = arrange(filter(tree, isTip == TRUE), y)$label)

                ## Analyze composition of each position
                alignment %>%
                  group_by(position) %>%
                  summarize(
                    gap_percent = sum(str_count(state, "-"))/n()*100,
                    conservation_percent = suppressWarnings(sum(state == mode(state))/n()*100)
                  ) %>%
                  pivot_longer(cols = c(2,3), names_to = "metric", values_to = "value") -> alignment_stats

                ## User interface

                ui <- fluidPage(

                    sidebarLayout(

                        sidebarPanel(

                            sliderInput(
                                inputId = "x_start",
                                label = "x_start",
                                min = 0,
                                max = max(alignment$position),
                                value = 0
                            ),

                            sliderInput(
                                inputId = "x_end",
                                label = "x_end",
                                min = 0,
                                max = max(alignment$position),
                                value = max(alignment$position)
                            ),

                            sliderInput(
                                inputId = "gap_threshold",
                                label = "Maximum Gap Percent",
                                min = 0,
                                max = 100,
                                value = 100
                            ),

                            sliderInput(
                                inputId = "conservation_threshold",
                                label = "Minimum Conservation Percent",
                                min = 0,
                                max = 100,
                                value = 0
                            ),

                            "Here is a NJ representation of the original alignment:",
                            plotOutput(
                                outputId = "tree",
                                height = "300px"
                            ),

                            "Here is a NJ representation of the alignment AFTER trimming:",
                            plotOutput(
                                outputId = "filtered_tree",
                                height = "300px"
                            ),

                            tags$head(
                                HTML(
                                    "
                                    <script>
                                        var socket_timeout_interval;
                                        var n = 0;
                                        $(document).on('shiny:connected', function(event) {
                                          socket_timeout_interval = setInterval(function(){
                                            Shiny.onInputChange('count', n++)
                                          }, 15000)
                                        });
                                        $(document).on('shiny:disconnected', function(event) {
                                          clearInterval(socket_timeout_interval)
                                        });
                                    </script>
                                    "
                                )
                            ),

                            textOutput("keepAlive"),
                            br(),
                            actionButton("return_btn", "Return Trimmed Alignment & Close App"),

                            width = 3
                        ),

                        mainPanel(

                            plotOutput(
                                outputId = "plot",
                                height = "1000px"
                            )
                        )
                    )
                )

                ## Server

                server <- function(input, output) {

                    ## Reactive to compute trimmed alignment from filtering criteria
                    trimmedAlignment <- reactive({

                        # Determine positions to keep based on filtering criteria
                        pos <- alignment_stats %>%
                            pivot_wider(names_from = "metric", values_from = "value") %>%
                            filter(
                                gap_percent < input$gap_threshold &
                                conservation_percent > input$conservation_threshold
                            ) %>%
                            pull(position)
                        validate(
                            need(length(pos) > 0, "No positions meet the filtering criteria. Please adjust your gap or conservation thresholds.")
                        )

                        # Create a wide format alignment with only the positions to keep
                        alignment_wide <- alignment %>%
                            filter(position %in% pos) %>%
                            pivot_wider(names_from = position, values_from = state)
                        validate(
                            need(ncol(alignment_wide) >= 2, "Trimmed alignment does not contain any sequence positions. Check your filtering criteria.")
                        )

                        # Build the trimmed alignment based on type
                        if (type == "DNA") {
                            trimmed_alignment <- Biostrings::DNAStringSet()
                            for (i in 1:nrow(alignment_wide)) {
                                seq_chars <- alignment_wide[i, 2:ncol(alignment_wide)]
                                trimmed_seq <- paste0(as.character(seq_chars), collapse = "")
                                trimmed_alignment <- c(trimmed_alignment, Biostrings::DNAStringSet(trimmed_seq))
                            }
                            # Convert factor to character before assigning names
                            trimmed_alignment@ranges@NAMES <- as.character(alignment_wide[[1]])
                            writeFasta(trimmed_alignment, paste0(alignment_in_path, "_trimmed"), fasta_type = "DNA")
                        }

                        if (type == "AA") {
                            trimmed_alignment <- Biostrings::AAStringSet()
                            for (i in 1:nrow(alignment_wide)) {
                                seq_chars <- alignment_wide[i, 2:ncol(alignment_wide)]
                                trimmed_seq <- paste0(as.character(seq_chars), collapse = "")
                                trimmed_alignment <- c(trimmed_alignment, Biostrings::AAStringSet(trimmed_seq))
                            }
                            # Convert factor to character before assigning names
                            trimmed_alignment@ranges@NAMES <- as.character(alignment_wide[[1]])
                            writeFasta(trimmed_alignment, paste0(alignment_in_path, "_trimmed"), fasta_type = "AA")
                        }

                        trimmed_alignment
                    })

                    ## Don't let it time out
                    output$keepAlive <- renderText({
                        req(input$count)
                        paste("\nstayin' alive ", input$count)
                    })

                    ## Original tree
                    output$tree <- renderPlot({
                        ggtree(tree)
                    })

                    ## Plot alignment
                    output$plot <- renderPlot({

                        ## Main plots
                        stats_plot <- ggplot(alignment_stats, aes(x = position, y = value, color = metric)) +
                            geom_line(size = 1) +
                            geom_hline(data = data.frame(
                                metric = "conservation_percent",
                                yintercept = input$conservation_threshold
                            ), aes(yintercept = yintercept)) +
                            geom_hline(data = data.frame(
                                metric = "gap_percent",
                                yintercept = input$gap_threshold
                            ), aes(yintercept = yintercept)) +
                            scale_color_viridis(discrete = TRUE, guide = "none") +
                            theme_classic() +
                            scale_y_continuous(name = "Percent") +
                            scale_x_continuous(expand = c(0,0)) +
                            coord_cartesian(xlim = c(input$x_start, input$x_end)) +
                            facet_grid(metric ~ .) +
                            theme(
                                axis.text.y = element_blank()
                            )

                        alignment_plot <- ggplot(alignment, aes(x = position, y = name)) +
                            geom_tile(aes(fill = state)) +
                            scale_fill_viridis(discrete = TRUE, guide = "none") +
                            scale_y_discrete(name = "") +
                            scale_x_continuous(expand = c(0,0)) +
                            theme_classic() +
                            coord_cartesian(xlim = c(input$x_start, input$x_end)) +
                            theme(
                                axis.text.y = element_blank()
                            )

                        alignment_stats %>%
                            pivot_wider(names_from = "metric", values_from = "value") %>%
                            filter(
                                gap_percent < input$gap_threshold &
                                conservation_percent > input$conservation_threshold
                            ) %>% 
                            select(position) %>% as.data.frame() -> positions_to_keep
                        positions_to_keep <- positions_to_keep[,1]

                        filtered_alignment_plot <- alignment %>%
                            filter(position %in% positions_to_keep) %>%
                            ggplot(aes(x = position, y = name)) +
                                geom_tile(aes(fill = state)) +
                                scale_fill_viridis(discrete = TRUE, guide = "none") +
                                scale_y_discrete(name = "") +
                                scale_x_continuous(expand = c(0,0)) +
                                theme_classic() +
                                coord_cartesian(xlim = c(input$x_start, input$x_end)) +
                                theme(
                                    axis.text.y = element_blank()
                                )

                        stats_plot + alignment_plot + filtered_alignment_plot + patchwork::plot_layout(ncol = 1)
                    })

                    ## Process and draw tree from filtered alignment
                    output$filtered_tree <- renderPlot({

                        # Trigger the reactive trimmed alignment (which writes the file)
                        trimmedAlignment()

                        trimmed_tree <- fortify(buildTree(
                            scaffold_type = scaffold_type,
                            scaffold_in_path = paste0(alignment_in_path, "_trimmed")
                        ))

                        ggtree(trimmed_tree)
                    })

                    ## Observer to handle the return button; when clicked, compute the trimmed alignment and close the app.
                    observeEvent(input$return_btn, {
                        trimmed <- trimmedAlignment()
                        stopApp(trimmed)
                    })
                }

                ## Call the app

                if (jupyter == TRUE) {

                    slot <- pick_shiny_port()
                    if (is.null(slot)) {
                        cat("All available ports are in use. Please try again later.")
                    } else {
                        cat(paste0("Connect at: ", slot$url))
                        return(runApp(shinyApp(ui = ui, server = server), host = "127.0.0.1", port = slot$port))
                    }

                }

                if (jupyter == FALSE) {
                    shinyApp(ui = ui, server = server)
                }
            }


        #### readGFFs


                #' readGFFs
                #'
                #' Read multiple genome feature files into a tidy dataframe
                #'
                #' This function reads multiple genome feature files (GFFs) and processes them into a tidy dataframe. It allows for various modes of subsetting and concatenation of GFFs based on user-defined parameters.
                #'
                #' @param input_frame A dataframe with columns: Genus_species, GFF_in_path, region_reach (20000), concatenation_spacing (10000)
                #' @param subset_mode A character vector of length 1. One of: "none", "goi_only", "goi_region", "name_check"
                #' @param goi A character vector of gene names. Labelling by species is not necessary
                #' @param concatenate_by_species A logical. If TRUE, concatenate all GFFs for a species into a single GFF
                #' @param chromosomes_only A logical. If TRUE, only keep chromosome features
                #' @param omit A character vector of columns to omit from the output dataframe
                #' @return A tidy dataframe containing the processed GFF data
                #' @examples
                #' \dontrun{
                #' # Example usage:
                #' # Read GFFs without any subsetting
                #' gff_data <- readGFFs(input_frame, subset_mode = "none")
                #'
                #' # Read GFFs and subset by genes of interest (goi)
                #' gff_data <- readGFFs(input_frame, subset_mode = "goi_only", goi = c("gene1", "gene2"))
                #'
                #' # Read GFFs and subset by regions around genes of interest (goi)
                #' gff_data <- readGFFs(input_frame, subset_mode = "goi_region", goi = c("gene1", "gene2"), region_reach = 20000)
                #'
                #' # Read GFFs and check for unique gene names
                #' gff_data <- readGFFs(input_frame, subset_mode = "name_check")
                #' }
                #' @export
                readGFFs <- function(
                                input_frame, # Dataframe with columns
                                subset_mode = c("none", "goi_only", "goi_region", "name_check"), 
                                goi = NULL, # Character vector of gene names. Labelling by species is not necessary
                                concatenate_by_species = TRUE,
                                chromosomes_only = FALSE,
                                omit = NULL
                            ) {
                    .check_bioc(c("rtracklayer", "GenomicRanges"), "readGFFs")
                    rlang::check_installed(c("plyr"), reason = "to use readGFFs()")

                    framed_GFFs <- list()

                    for ( i in 1:dim(input_frame)[1] ) {

                        # Load GFFs
                            cat(paste("Loading GFF for ", input_frame$Genus_species[i], "...\n", sep = ""))
                            gff_as_GRange <- suppressWarnings(rtracklayer::import.gff(as.character(input_frame$GFF_in_path[i]), version = "3"))
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

                        # Remove cols to omit
                            if ( length(omit) > 0 ) {
                                ranges <- ranges[, !colnames(ranges) %in% omit]    
                            }
                            
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
                                framed_GFFs <- ranges
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
                                    library(plyr)
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
                                    framed_GFFs <- do.call(rbind, framed_GFFs)
                            }
                    }

                    ## Concatenate range_scaffs for single species
                        if ( concatenate_by_species == TRUE & dim(input_frame)[1] == 1) {

                            ## Filter for only main chromosomes if specified, then number seqnames
                                if (chromosomes_only == TRUE) {
                                    chromosome_seqnames <- as.character(framed_GFFs$seqnames[framed_GFFs$chromosome %in% c(dropNA(as.numeric(unique(framed_GFFs$chromosome))), "X", "Y")])
                                    message("NA introduced by coersion here is OKAY!")
                                    framed_GFFs <- framed_GFFs[framed_GFFs$seqnames %in% chromosome_seqnames,]
                                    framed_GFFs$seqname_number <- 0
                                    framed_GFFs$seqname_number <- match(framed_GFFs$seqnames, chromosome_seqnames)
                                } else {
                                    chromosome_seqnames <- unique(framed_GFFs$seqnames)
                                    framed_GFFs$seqname_number <- 0
                                    framed_GFFs$seqname_number <- match(framed_GFFs$seqnames, chromosome_seqnames)
                                }
                            if ( length(unique(framed_GFFs$seqnames)) > 1 ) {
                                end_previous_seqname_number_region <- max(framed_GFFs[framed_GFFs$seqname_number == 1,]$end)
                                for (this_seqname_number in 2:length(unique(framed_GFFs$seqname_number))) {
                                    amount_to_advance_this_seqname_number_region <- end_previous_seqname_number_region + input_frame$concatenation_spacing[i]
                                    framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$start <- framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$start + amount_to_advance_this_seqname_number_region
                                    framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$end <- framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$end + amount_to_advance_this_seqname_number_region
                                    end_previous_seqname_number_region <- max(framed_GFFs[framed_GFFs$seqname_number == this_seqname_number,]$end)
                                }
                            }
                        }


                    ## Center the scaffolds, merging the GFFs if there is more than one
                        if (dim(input_frame)[1] > 1) {
                            framed_GFFs <- do.call(plyr::rbind.fill, framed_GFFs)
                            if ( subset_mode != "name_check" ) {
                                center <- (max(framed_GFFs$end) - min(framed_GFFs$start))/2
                                library(plyr)
                                framed_GFFs <- plyr::ddply(framed_GFFs, .(Genus_species), mutate, center_start = (start + (center-(max(end)-min(start))/2)), center_end = (end + (center-(max(end)-min(start))/2)) )
                            }
                        } else {
                            if ( subset_mode != "name_check" ) {
                                center <- (max(framed_GFFs$end) - min(framed_GFFs$start))/2
                                framed_GFFs$center_start = framed_GFFs$start + (center-(max(framed_GFFs$end)-min(framed_GFFs$start))/2)
                                framed_GFFs$center_end = framed_GFFs$end + (center-(max(framed_GFFs$end)-min(framed_GFFs$start))/2)
                            }
                        }

                    return(framed_GFFs)
                }


        #### fastxQC


            #' Run QC on fastx file(s)
            #'
            #' @param paths_to_fastxs A list of the paths to the fastx file(s) you want to analyze
            #' @param max_n_seqs The maximum number of sequences PER FILE to analyze. Default is no limit.
            #' @param type Not yet documented.
            #' @param mode Not yet documented.
            #' @param progress Not yet documented.
            #' @export
            fastxQC <-  function(
                            paths_to_fastxs,
                            type = c("fastq", "fasta"),
                            mode = c("fast", "slow"),
                            max_n_seqs = NULL,
                            progress = FALSE
                        ) {

                output <- list()
                for ( path_to_fastx in 1:length(paths_to_fastxs)) {

                    cat(paste("Analyzing file ", paths_to_fastxs[path_to_fastx], "\n", sep = ""))

                    if (length(grep("http", paths_to_fastxs[path_to_fastx])) > 0) {
                        downloadFasta(paths_to_fastxs[path_to_fastx])
                        paths_to_fastxs[path_to_fastx] <- "TEMP___fasta"
                    }

                    ## Determine max number of reads

                        if (length(max_n_seqs) > 0) {
                            n_reads <- max_n_seqs
                        } else {
                            if (type == "fasta") {
                                n_reads <- (dim(data.table::fread(paths_to_fastxs[path_to_fastx], sep = NULL, header = FALSE))[1]/2)
                            }
                            if (type == "fastq") {
                                n_reads <- (dim(data.table::fread(paths_to_fastxs[path_to_fastx], sep = NULL, header = FALSE))[1]/4)
                            }
                        }

                    ## Get length distribution. If fastq, convert to fasta first. If slow, do quality scores

                        if (type == "fasta") {

                            output[[length(output)+1]] <- cbind(
                                data.frame(
                                    file = rep(paths_to_fastxs[path_to_fastx], n_reads)
                                ),
                                read.table(Rsamtools::indexFa(
                                    paths_to_fastxs[path_to_fastx]
                                ))[1:n_reads,1:2]
                            )

                        }

                        if (type == "fastq" & mode == "fast") {
                            
                            cat(paste("Converting to fasta and getting sequence lengths...", "\n", sep = ""))

                            fastqToFasta(
                                file_in_path = paths_to_fastxs[path_to_fastx],
                                file_out_path = paste0(paths_to_fastxs[path_to_fastx], ".fasta")
                            )

                            output[[length(output)+1]] <- cbind(
                                data.frame(
                                    file = rep(paths_to_fastxs[path_to_fastx], n_reads)
                                ),
                                read.table(Rsamtools::indexFa(
                                    paste0(paths_to_fastxs[path_to_fastx], ".fasta")
                                ))[1:n_reads,1:2]
                            )

                        }

                        if (type == "fastq" & mode == "slow") {

                            cat(paste("Using Phred_ASCII_33 score conversions...", "\n", sep = ""))
                            cat(paste("Analyzing quality scores...", "\n", sep = ""))

                            if( progress ) { pb <- progress::progress_bar$new(total = n_reads) }
                            for (i in 1:n_reads) {

                                data <- data.table::fread(
                                    paths_to_fastxs[path_to_fastx], skip = (4*(i-1)),
                                    nrows = 4, sep = NULL, header = FALSE
                                )

                                # if (!length(strsplit(as.character(unlist(data[2,])), split = "")[[1]]) == length(strsplit(as.character(unlist(data[4,])), split = "")[[1]])) {
                                #     stop(paste0("Problem with quality score encoding on line ", i))
                                # }

                                read_score <- mean(as.numeric(phred33_lookup$X2[match(
                                                    strsplit(as.character(unlist(data[4,])), split = "")[[1]],
                                                    phred33_lookup$X1
                                                )]))

                                if (is.na(read_score)) {

                                    bad_characters <- strsplit(as.character(unlist(data[4,])), split = "")[[1]][which(is.na(match(
                                        strsplit(as.character(unlist(data[4,])), split = "")[[1]],
                                        phred33_lookup$X1
                                    )))]

                                    stop(paste0(
                                        "Problem with quality score encoding on line ", i, ".", "\n",
                                        "Unknown characters: ", bad_characters))

                                }
                                
                                output[[length(output)+1]] <- data.frame(
                                    file = paths_to_fastxs[path_to_fastx],
                                    name = data[1],
                                    length = nchar(data[2,]),
                                    score = read_score,
                                    accuracy = (1-10^(-read_score/10))*100
                                )

                                if( progress ) { pb$tick() }
                            
                            }
                        }
                }

                output <- do.call(rbind, output)
                rownames(output) <- NULL
                colnames(output) <- c("file", "name", "length", "score", "accuracy")[1:dim(output)[2]]

                return(as_tibble(output))
            }


        #### reverseTranslate

            #' Reverse translate a protein sequence to its corresponding nucleotide sequence
            #'
            #' This function takes a protein sequence from a fasta file and reverse translates it into a nucleotide sequence using a specified codon table. The resulting nucleotide sequence is then written to an output fasta file.
            #'
            #' @param fasta_in_path A character string specifying the path to the input fasta file containing the protein sequence(s).
            #' @param fasta_out_path A character string specifying the path to the output fasta file where the nucleotide sequence(s) will be written.
            #' @param organism A character string specifying the organism whose codon table should be used for reverse translation. Options are "Arabidopsis_thaliana", "Zea_mays", "Nicotiana_tabacum", or "Saccharomyces_cerevisiae". The default is to use the codon table for "Arabidopsis_thaliana".
            #' @return This function does not return a value. It writes the reverse-translated nucleotide sequences to the specified output fasta file.
            #' @details The function reads a codon table from a predefined URL and filters it based on the specified organism. It then constructs a substitution table to map amino acids to their most frequently used codons in the organism. Each protein sequence in the input file is reverse translated, ensuring a terminal stop codon is present, and the resulting nucleotide sequences are saved to the output file.
            #' @examples
            #' \dontrun{
            #' # Reverse translate a protein sequence for Arabidopsis thaliana
            #' reverseTranslate("https://thebustalab.github.io/R_For_Chemists/sample_data/peptide.fasta", "peptide.fasta", "Arabidopsis_thaliana")
            #' }
            #' @export
            reverseTranslate <- function(
                fasta_in_path,
                fasta_out_path,
                organism = c("Arabidopsis_thaliana", "Zea_mays", "Nicotiana_tabacum", "Saccharomyces_cerevisiae")
            ) {
                .check_bioc(c("Biostrings"), "reverseTranslate")

                # Read the codon table for the specified organism
                codon_tables <- tryCatch(
                    suppressWarnings(readMonolist("https://thebustalab.github.io/phylochemistry/sample_data/codon_tables.csv")),
                    error = function(e) stop("reverseTranslate() needs the codon usage table sample_data/codon_tables.csv, ",
                                             "which is not available on the course website, so this function cannot run.", call. = FALSE)
                )
                codon_tables$Genus_species <- paste(codon_tables$Genus, codon_tables$species, sep = "_")
                codon_tables <- dplyr::filter(codon_tables, Genus_species == organism)
                    
                # Create a substitution table mapping amino acids to codons
                substitution_table <- list()
                for (i in 1:length(unique(codon_tables$amino_acid))) {
                    this_codon <- dplyr::filter(codon_tables, amino_acid == unique(codon_tables$amino_acid)[i])
                    substitution_table[[i]] <- this_codon[which.max(this_codon$frequency),]
                }
                substitution_table <- do.call(rbind, substitution_table)
                rownames(substitution_table) <- NULL

                # Read the amino acid sequences from the input fasta file
                amin_seqs <- Biostrings::readAAStringSet(fasta_in_path)
                output <- Biostrings::DNAStringSet()
                for (i in 1:length(amin_seqs)) {
                    this_amin_seq <- strsplit(as.character(amin_seqs[i]), split = NULL)[[1]]

                    # Add a terminal stop codon if not present
                    if( this_amin_seq[length(this_amin_seq)] != "*" ) {
                        this_amin_seq <- c(this_amin_seq, "*")
                    }
                    
                    # Determine codons for each amino acid and concatenate them
                    output[[i]] <- paste(substitution_table$codon[match(this_amin_seq, substitution_table$amino_acid)], collapse = "")

                }
                # Set the names for the output sequences and write them to the output file
                output@ranges@NAMES <- amin_seqs@ranges@NAMES
                Biostrings::writeXStringSet(output, fasta_out_path)
            }


        #### fastqToFasta


            #' Convert fastq to fasta
            #'
            #' @param file_in_path The path to the fastq file
            #' @param file_out_path The path to the fasta file to write
            #' @examples
            #' \dontrun{
            #' fastqToFasta("https://thebustalab.github.io/R_For_Chemists/sample_data/reads.fastq", "reads.fasta")
            #' }
            #' @export
            fastqToFasta <- function(file_in_path, file_out_path) {
                .check_bioc(c("Biostrings"), "fastqToFasta")
                data <- readLines(file_in_path)
                fasta <- Biostrings::DNAStringSet(data[seq(1,length(data),4)+1])
                names(fasta) <- data[seq(1,length(data),4)]
                Biostrings::writeXStringSet(fasta, file_out_path)
            }


        #### canuHistogramToKmerTable


            #' Convert a kmer report file from Canu into a table that can be accepted by genomescope
            #' Use http://qb.cshl.edu/genomescope/genomescope2.0/
            #'
            #' @param file_in_path The path to the kmer report file
            #' @param file_out_path The path to the table file to write
            #' @examples
            #' \dontrun{
            #' canuHistogramToKmerTable("https://thebustalab.github.io/R_For_Chemists/sample_data/kmer_report.txt", "kmer_report.txt")
            #' }
            #' @export
            canuHistogramToKmerTable <- function(file_in_path, file_out_path) {

                dat <- readLines(file_in_path)
                dat <- dat[11:length(dat)]
                dat <- gsub("1\\.0.*$", "", gsub("0\\..*$", "", dat))
                dat <- strsplit(paste(dat, collapse = " "), " ")[[1]]
                dat <- dat[dat != ""]

                dat <- data.frame(
                    x = as.numeric(dat[seq(1,length(dat), 2)]),
                    y = as.numeric(dat[seq(2,length(dat), 2)])
                )

                head(dat)
                tail(dat)

                write.table(
                    dat, file = file_out_path,
                    sep = " ", row.names = FALSE, col.names = FALSE
                )

            }


        #### drawAlignment


            #' Generate a multiple alignment graphic using ggplot
            #'
            #' Allows the user to plot a multiple alignent using ggplot's grammar of graphics
            #' @param infile The alignment to use (fasta file)
            #' @param monolist A dataframe with metadata for the alignment
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
            #' @examples
            #' \dontrun{
            #' drawAlignment()
            #' }
            #' @export
            drawAlignment <- function(   
                                infile,
                                monolist,
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
                    AA_phydat <- as.data.frame(as.character(as.matrix(ape::read.FASTA(file = infile, type = "AA"))))
                    AA_phydat <- as.data.frame(t(AA_phydat))
                    colnames(AA_phydat) <- as.character(seq(1,dim(AA_phydat)[2],1))
                    AA_alignment_df <-  cbind(  
                                            data.frame(protein = paste(
                                                monolist[,which(colnames(monolist) == alignment_labels)][match(rownames(as.matrix(AA_phydat)), monolist$accession)]
                                            )),
                                            # data.frame(funct = monolist[,colnames(monolist)=="Function"][match(rownames(as.matrix(AA_phydat)), monolist[,1])]),
                                            data.frame(y = rep(0,dim(AA_phydat)[1])),
                                            AA_phydat
                                        )
                    # AA_alignment_df_plottable <- tidyr::gather(AA_alignment_df
                    first_number_col <- min(which(is.na(suppressWarnings(as.numeric(colnames(AA_alignment_df)))) == FALSE))
                    AA_alignment_df_plottable <- as.data.frame(tidyr::pivot_longer(AA_alignment_df, cols = first_number_col:dim(AA_alignment_df)[2], names_to = "position", values_to = "residue"))
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
                                alignment_matrix <- toupper(t(as.data.frame(as.character(as.matrix(ape::read.FASTA(file = infile, type = "AA"))))))
                                
                                spec1 <- list()
                                spec1[[1]] <- monolist[,colnames(monolist) == funct_assoc_data[j]] == as.character(gsub("~.*$", "", funct_assoc_data[j]))
                                alignment_matrix_dim1 <- alignment_matrix[Reduce("&", spec1),]
                                dim(alignment_matrix_dim1)

                                spec2 <- list()
                                spec2[[1]] <- monolist[,colnames(monolist)==funct_assoc_data[j]] == as.character(gsub(".*~", "", funct_assoc_data[j]))
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
                                                                    name = ""
                                                                    # breaks = as.numeric(generateTicks(seq(0,5000,tick_spacing))$all_breaks),
                                                                    # labels = as.character(generateTicks(seq(0,5000,tick_spacing))$all_labels)
                                                                )
                            plot_list[[i]] <- plot_list[[i]] + scale_y_continuous(name = "", breaks = 0)
                            plot_list[[i]] <- plot_list[[i]] + scale_color_manual(values = color_pal)
                            plot_list[[i]] <- plot_list[[i]] + theme(
                                        panel.spacing.y = unit(0.01, "lines"),
                                        panel.spacing.x = unit(0.2, "cm"),
                                        panel.border = element_blank(),
                                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                                        axis.text.y = element_blank(),
                                        text = element_text(size = ticks_text_size),
                                        legend.position = "none",
                                        axis.title = element_text(color = "#737373", face = "bold", size = 40),
                                        axis.ticks.length = unit(0.2, "cm"),
                                        axis.ticks = element_line(color = "#737373", size = 1, lineend = 6),
                                        axis.text = element_text(color = "#737373", face = "bold", size = ticks_text_size),
                                        axis.line = element_line(color = "#737373", size = 1),
                                        strip.text = element_text(hjust = 1, size = ticks_text_size),
                                        strip.text.y.left = element_text(angle = 0),
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


        #### annotatePDB

            #### annotatePDB


                    #' Import chemcial data
                    #'
                    #' Allows the user to plot a multiple alignent using ggplot's grammar of graphics
                    #' @param pdb_in_path Not yet documented.
                    #' @param alignment_in_path Not yet documented.
                    #' @param pdb_name_in_alignment Not yet documented.
                    #' @param highlights_data Not yet documented.
                    #' @param out_file Not yet documented.
                    #' @examples
                    #' \dontrun{
                    #' annotatePDB()
                    #' }
                    #' @export
                    annotatePDB <- function(pdb_in_path, alignment_in_path, pdb_name_in_alignment, highlights_data, out_file) {
                        rlang::check_installed(c("bio3d"), reason = "to use annotatePDB()")

                        ## Read in PDB file
                            pdb <- bio3d::read.pdb(pdb_in_path)
                            pdb$atom
                            pdb_sequence <- bio3d::aa321(pdb$atom$resid[bio3d::atom.select(pdb, "calpha")$atom])
                            pdb_frame <- data.frame(pdb_sequence=pdb_sequence, pdb_position=seq(1,length(pdb_sequence),1))
                            pdb_frame$pdb_position <- as.numeric(pdb_frame$pdb_position)
                            head(pdb_frame)

                        ## Read in alignment
                            AA_phydat <- as.data.frame(as.character(as.matrix(ape::read.FASTA(file = alignment_in_path, type = "AA"))))
                            alignment_sequence <- as.data.frame(AA_phydat)[,grep(pdb_name_in_alignment, colnames(as.data.frame(AA_phydat)))]
                            alignment_frame <- data.frame(alignment_sequence = alignment_sequence, alignment_position = seq(1,length(alignment_sequence),1))
                            head(alignment_frame)

                        ## Add spaces at the front of the PDB frame until its first AA lines up with alignment

                            head(pdb_frame)
                            head(alignment_frame)

                            space <- data.frame(pdb_sequence = "-", pdb_position = "NA")
                            space$pdb_position <- as.numeric(as.character(space$pdb_position))

                            for (i in 1:length(pdb_sequence)) {
                                if ( as.character(alignment_frame[i,1]) == as.character(pdb_frame[i,1]) ) {
                                } else {
                                    pdb_frame <- rbind(pdb_frame[0:(i-1),],space,pdb_frame[i:dim(pdb_frame)[1],])
                                }
                            }

                            head(pdb_frame, 30)
                            head(alignment_frame, 30)

                        ## Add spaces to the end of the PDB frame until it lines up with alignment

                            tail(pdb_frame)
                            tail(alignment_frame)

                            final_space <- data.frame(pdb_sequence="-", pdb_position="NA")
                            final_space$pdb_position <- as.numeric(as.character(final_space$pdb_position))        
                            
                            if (((dim(alignment_frame)[1]-dim(pdb_frame)[1])-1) > 0) {
                                for (i in 1:((dim(alignment_frame)[1]-dim(pdb_frame)[1])-1)) {
                                    final_space <- rbind(final_space, space)
                                }
                                pdb_frame <- rbind(pdb_frame, final_space)
                            } else {
                                pdb_frame <- rbind(pdb_frame, final_space)
                            }            

                            tail(pdb_frame, 20)
                            tail(alignment_frame, 20)
                            
                        ## Bind them together - they should have the same dimensions at this point
                            
                            dim(pdb_frame)
                            dim(alignment_frame)
                            translational_frame <- cbind(alignment_frame, pdb_frame)

                        ## Write it out

                            if (file.exists(out_file)) { file.remove(out_file) }

                            write("## set background", file=out_file, append=TRUE)
                            write("set bg_rgb, white", file=out_file, append=TRUE)
                            write("", file=out_file, append=TRUE)
                            write("## Set base color", file=out_file, append=TRUE)
                            write("set_color myblue = [55,126,184]", file=out_file, append=TRUE)
                            write("select base, resi 1-800", file=out_file, append=TRUE)
                            write("color myblue, base", file=out_file, append=TRUE)
                            write("", file=out_file, append=TRUE)
                            write("## General image improvement", file=out_file, append=TRUE)
                            write("set spec_reflect, 0.2", file=out_file, append=TRUE)

                            for (i in 1:dim(highlights_data)[1]) {
                                color <- paste("[", col2rgb(highlights_data$color[i])[1], ",", col2rgb(highlights_data$color[i])[2], ",", col2rgb(highlights_data$color[i])[3], "]", sep = "") 
                                write(paste("set_color ", paste("COLOR", i, sep = ""), " = ", color, sep = ""), file=out_file, append=TRUE)
                                write(paste("select ", paste("RES", i, sep = ""), paste(", resi", translational_frame$pdb_position[match(highlights_data$xint[i], translational_frame$alignment_position)]), sep = ""), file=out_file, append=TRUE)
                                write(paste("color ", paste("COLOR", i, sep = ""), paste(", RES", i, sep = ""), sep = ""), file=out_file, append=TRUE)
                                write(paste("show sticks", paste(", RES", i, sep = ""), sep = ""), file=out_file, append=TRUE)
                                write(paste(""), file=out_file, append=TRUE)
                            }

                            return(translational_frame)
                    }

