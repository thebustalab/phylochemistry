## io — moved verbatim from phylochemistry.R (package conversion, commit 1)


    ##### Read and write functions

        #### readCSV


            #' Interactive selection of a CSV file to read
            #'
            #' This function allows the user to interactively select a CSV file to read.
            #' The file can be delimited by either a comma or a space.
            #'
            #' @param sep A character vector indicating the delimiter to use when reading the file.
            #'            The default is a comma (","), but a space (" ") can also be used.
            #' @return A data frame containing the contents of the selected CSV file.
            #' @examples
            #' \dontrun{
            #' \dontrun{
            #' # To read a CSV file with comma delimiter
            #' data <- readCSV()
            #'
            #' # To read a CSV file with space delimiter
            #' data <- readCSV(sep = " ")
            #' }
            #' }
            #' @export
            readCSV <- function(sep = c(",", " ")) {
                if (sep[1] == ",") {
                    return(readr::read_csv(file.choose()))
                }  
                if (sep[1] == " ") {
                    return(utils::read.delim(file.choose(), delim = " "))
                }
            }


        #### writeCSV


            #' Write a monolist to a CSV file
            #'
            #' This function allows the user to interactively select a location to save a monolist as a CSV file.
            #'
            #' @param monolist A data frame representing the monolist to be written to a CSV file.
            #' @return None. The function writes the monolist to a file and does not return a value.
            #' @examples
            #' \dontrun{
            #' \dontrun{
            #' # To write a monolist to a CSV file
            #' writeCSV(monolist)
            #' }
            #' }
            #' @export
            writeCSV <- function(monolist) {
                writeMonolist(monolist = monolist, monolist_out_path = file.choose(new = TRUE))
            }


        #### readMonolist


            #' Read a monolist from a CSV file or URL
            #'
            #' This function reads a monolist from a specified path or URL. The monolist should be in CSV format.
            #'
            #' @param monolist_in_path A string specifying the path or URL to the monolist CSV file.
            #' @return A data frame containing the contents of the monolist.
            #' @examples
            #' \dontrun{
            #' \dontrun{
            #' # To read a monolist from a local file
            #' monolist <- readMonolist("path/to/monolist.csv")
            #'
            #' # To read a monolist from a URL
            #' monolist <- readMonolist("http://example.com/monolist.csv")
            #' }
            #' }
            #' @export
            readMonolist <- function( monolist_in_path ) {
                if (length(grep("http", monolist_in_path)) > 0) {
                    monolist <- read_csv(monolist_in_path, show_col_types = FALSE)        
                } else {
                    monolist <- as.data.frame(data.table::fread(file = monolist_in_path))
                }
                return( monolist )
            }


        #### writeMonolist


            #' Writes a monolist
            #'
            #' @param monolist The monolist to write out
            #' @param monolist_out_path The path to where the monolist should be written
            #' @param ... Not yet documented.
            #' @export
            writeMonolist <- function( monolist, monolist_out_path, ... ) {
                write.table(monolist, file = monolist_out_path, sep = ",", row.names = FALSE, ...)
                # write_csv(monolist, file = monolist_out_path)
            }


        #### readPolylist


            #' Reads a polylist
            #'
            #' Reads a spreadsheet in polylist format (wide format, multiple column and row headers), and converts it into tidy format
            #' @param polylist_in_path The path to the polylist (in .csv format) to be read
            #' @param centroid The characters in the cell that defines the boundaries of the multiple row and column headers
            #' @param table_value_unit The name to be given to the column that will contain the values in the polylist table
            #' @export
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
                    center_row <- as.numeric(grep(centroid, polylist[,center_column]))

                # Use centroid to extract vertical_monolist
                    vertical_monolist <- as.data.frame(polylist[(center_row+1):dim(polylist)[1], 1:(center_column-1)])
                    colnames(vertical_monolist) <- as.character(unlist(polylist[center_row, 1:(center_column-1)]))
                    vertical_monolist$URI_URI_URI <- apply(vertical_monolist, 1, function(x) paste(x, collapse = ""))

                # Use centroid to extract horizontal_monolist
                    horizontal_monolist <- as.data.frame(t(polylist[-c(center_row), (center_column+1):dim(polylist)[2]]))
                    rownames(horizontal_monolist) <- NULL
                    colnames(horizontal_monolist) <-    c(
                                                        as.character(unlist(polylist[c(1:(center_row-1), (center_row+1):dim(polylist)[1]), center_column]))[1:(center_row-1)],
                                                        as.character(vertical_monolist$URI_URI_URI)
                                                    )
                    horizontal_monolist <- pivot_longer(horizontal_monolist, cols = center_row:ncol(horizontal_monolist), names_to = "URI_URI_URI", values_to = "table_value_unit")

                    colnames(horizontal_monolist)[colnames(horizontal_monolist) == "table_value_unit"] <- table_value_unit

                # Bind the two monolists, drop the URI column
                    polylist <- cbind(horizontal_monolist,vertical_monolist[match(horizontal_monolist$URI_URI_URI, vertical_monolist[,colnames(vertical_monolist) == "URI_URI_URI"]),])
                    polylist <- polylist[,colnames(polylist) != "URI_URI_URI"]

                # Convert to numeric
                    for ( i in 1:dim(polylist)[2] ) {
                        if ( numbersOnly(polylist[,i]) ) {
                            polylist[,i] <- as.numeric(polylist[,i])
                        }
                    }

                # Reset row numbers, return the polylist
                    rownames(polylist) <- NULL
                    return(polylist)
            }


        #### writePolylist


            #' Writes a polylist
            #'
            #' This function allows the user to interactively select a location to save a polylist as a CSV file.
            #'
            #' @param polylist The polylist to write out
            #' @param polylist_out_path The path to where the polylist should be written
            #' @export
            writePolylist <- function( polylist, polylist_out_path ) {
                write.table( polylist, file = polylist_out_path, sep = ",", row.names = FALSE, col.names = FALSE)
            }


        #### writeSupplementalTable


            #' Writes a human-readable supplemental table of quantitative data
            #'
            #' @param supplementalTable The data to write out
            #' @param round_level The number of decimal places summary statistics should have
            #' @param file_out_path The path to the location where the file should be written
            #' @param replicates Not yet documented.
            #' @export
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


        #### renameFile

            #' Rename (move) a file, creating the destination folder if needed
            #'
            #' @param from Not yet documented.
            #' @param to Not yet documented.
            #' @export
            renameFile <- function(from, to) {
                todir <- dirname(to)
                if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
                file.rename(from = from,  to = to)
            }


        #### openRGD


                #' Open an RScript stored on Google Drive
                #'
                #' @param drive_share_link Not yet documented.
                #' @export
                openRGD <- function(drive_share_link) {
                    rlang::check_installed(c("googledrive"), reason = "to use openRGD()")

                    file_name <- googledrive::drive_get(
                        id = googledrive::as_id(drive_share_link)
                    )$name

                    googledrive::drive_rename(
                        file = googledrive::as_id(drive_share_link),
                        name = paste0("IN_USE___", file_name)
                    )

                    googledrive::drive_download(
                        file = googledrive::as_id(drive_share_link),
                        path = paste0("IN_USE___", file_name),
                        overwrite = TRUE
                    )

                    file.edit(paste0("IN_USE___", file_name))

                }


        #### closeRGD


            #' Save and close the current RScript to Google Drive
            #'
            #' @param drive_share_link Not yet documented.
            #' @export
            closeRGD <- function(drive_share_link){
                rlang::check_installed(c("rstudioapi", "googledrive"), reason = "to use closeRGD()")

                rstudioapi::documentSave()

                file_name <- googledrive::drive_get(
                    id = googledrive::as_id(drive_share_link)
                )$name

                googledrive::drive_update(
                    file = googledrive::as_id(drive_share_link),
                    media = paste0(file_name)
                )

                googledrive::drive_rename(
                    file = googledrive::as_id(drive_share_link),
                    name = paste0(gsub("IN_USE___", "", file_name))
                )

                rstudioapi::documentClose()

                file.remove(paste0(file_name))

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
            #' @export
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
                        measurements_monolist <- readMonolist(monolist_in_path = measurements_monolist_in_path)
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


    ##### Other functions

        #### OsDirectoryPathCorrect


            #' OS-aware correction of directory paths
            #'
            #' @param directory_path The path to correct
            #' @export
            OsDirectoryPathCorrect <- function( directory_path ) {

                ## Detect OS

                    OS <- .Platform$OS.type

                ## Add terminal slash if missing

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
                            ) == "\\\\" ) {
                        } else {
                            directory_path <- paste(directory_path, "\\\\", sep = "")
                        }
                        directory_path_corrected <- directory_path

                    } else {

                        warning("ERROR: OS could not be identified")

                    }

                ## Return
                    
                    return(directory_path_corrected)

            }


        #### OsPathCorrect


            #' OS-aware correction of paths
            #'
            #' @param path The path to correct
            #' @export
            OsPathCorrect <- function( path ) {

                ## Detect OS

                    OS <- .Platform$OS.type

                ## Replace single windows slashes with doubles

                    if (OS == "unix") {

                        path_corrected <- path

                    } else if (OS == "windows"){

                        path_corrected <- gsub("\\", "\\\\", path)

                    } else {

                        warning("ERROR: OS could not be identified")

                    }

                ## Return
                    
                    return(path_corrected)

            }


        #### dropNA


            #' drops NAs from a vector
            #'
            #' @param x The vector to drop NAs from
            #' @export
            dropNA <- function( x ) { return(x[!is.na(x)]) }


        #### numbersOnly

            #' Test whether every element can be read as a number
            #'
            #' @param x Not yet documented.
            #' @export
            numbersOnly <- function(x) { all(suppressWarnings(!is.na(as.numeric(as.character(x))))) }


        #### generateQRcodes


            #' Create major and minor axes ticks and labels
            #'
            #' @param share_link Share link to Google Sheet with QR code metadata
            #' @export
            generateQRcodes <- function( share_link ) {
                rlang::check_installed(c("qrencoder", "cowplot", "ggpubr", "png", "googlesheets4"), reason = "to use generateQRcodes()")

                ## Read in the Google docs spreadsheet

                    data <- googlesheets4::read_sheet(share_link)

                ## Create QR codes for the plants

                    plot_list <- list()
                    for (i in 1:dim(data)[1]) {
                        
                        metadata <- unite(data, col = "metadata")$metadata

                        png(paste0(metadata, ".png"))
                            image(
                                1 - qrencoder::qrencode(metadata),
                                asp = 1, xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1),
                                col = c("black", "white")
                            )
                        dev.off()

                        plot_list[[i]] <- cowplot::plot_grid(
                            ggpubr::text_grob(metadata, size = 5),
                            grid::rasterGrob(png::readPNG(paste0(metadata, ".png"))),
                            ncol = 1, rel_heights = c(1,30)
                        )

                        file.remove(paste0(metadata, ".png"))

                    }

                ## return
                    return(plot_list)
                    message("Consider using 'do.call(\"grid.arrange\", c(plot_list))' to plot all the plots")
        }

