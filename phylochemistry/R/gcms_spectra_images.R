## gcms_spectra_images — moved verbatim from phylochemistry.R (package conversion, commit 1)


    ##### Chemical data handling

        #### convertCDFstoCSVs

            #' Convert mass spectral data files (CDF) into a CSV file
            #'
            #' This function processes mass spectral data stored in CDF files and converts them into CSV format. It allows filtering based on m/z values and retention times, and can optionally force conversion even if CSV files already exist.
            #'
            #' @param paths_to_cdfs A character vector containing the paths to the CDF files that need to be converted.
            #' @param min_mz A numeric value specifying the smallest m/z (mass-to-charge ratio) value to process. This acts as a filter to exclude data below this threshold.
            #' @param max_mz A numeric value specifying the largest m/z value to process. This acts as a filter to exclude data above this threshold.
            #' @param min_rt (optional) A numeric value specifying the smallest retention time to process. This acts as a filter to exclude data below this threshold. Default is NULL, meaning no minimum retention time filter is applied.
            #' @param max_rt (optional) A numeric value specifying the largest retention time to process. This acts as a filter to exclude data above this threshold. Default is NULL, meaning no maximum retention time filter is applied.
            #' @param force A logical value indicating whether to convert all CDF files even if they have already been converted previously. Default is FALSE, meaning only files without existing CSV counterparts will be converted.
            #' @examples
            #' # Convert CDF files to CSV with default m/z range and no retention time filtering
            #' convertCDFstoCSVs(paths_to_cdfs = c("file1.cdf", "file2.cdf"))
            #'
            #' # Convert CDF files to CSV with specified m/z and retention time ranges, forcing conversion
            #' convertCDFstoCSVs(paths_to_cdfs = c("file1.cdf", "file2.cdf"), min_mz = 100, max_mz = 700, min_rt = 5, max_rt = 20, force = TRUE)
            #'
            #' @export
            convertCDFstoCSVs <- function(
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

                        cat(paste("Reading data file ", paths_to_cdfs[file], sep = ""))
                            rawDataFile <- xcms::loadRaw(xcms::xcmsSource(paths_to_cdfs[file]))

                        cat("   Framing data file ... \n")
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
                            if ( dim(table(duplicated(paste(framedDataFile$mz, framedDataFile$rt, sep = "_")))) > 1 ) {
                                framedDataFile %>% group_by(mz,rt) %>% summarize(intensity = sum(intensity), .groups = "drop") -> framedDataFile
                                framedDataFile <- as.data.frame(framedDataFile)
                            }

                        cat("   Filling in blank m/z values ...\n")
                            spreadDataFile <- tidyr::spread(framedDataFile, rt, intensity, fill = 0)
                            
                            numbers <- round(seq(min_mz, max_mz, 0.1), digits = 1)
                            add <- numbers[!numbers %in% spreadDataFile$mz]
                            addition <- spreadDataFile[1:length(add),]
                            
                            if ( length(add) > 0 ) {
                                addition$mz <- add     
                            }
                            
                            addition <- tidyr::pivot_longer(addition, cols = 2:dim(addition)[2], names_to = "rt", values_to = "intensity")
                            addition$intensity <- 0
                            addition$rt <- as.numeric(addition$rt)
                            
                            framedDataFile <- tidyr::pivot_longer(spreadDataFile, cols = 2:dim(spreadDataFile)[2], names_to = "rt", values_to = "intensity")

                        cat("   still working ...\n")
                            framedDataFile <- rbind(framedDataFile, addition)

                        cat("   still working ... patience please\n")
                            framedDataFile$rt <- as.numeric(framedDataFile$rt)
                            framedDataFile <- framedDataFile[sort.list(framedDataFile$mz),]

                        cat("   still working ... almost done\n")
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

                        cat("   Writing out data file as CSV... \n\n")
                            data.table::fwrite(framedDataFile, file = paste(paths_to_cdfs[file], ".csv", sep = ""), col.names = TRUE, row.names = FALSE)
                    }
                }
            }


        #### extractChromatogramsFromCSVs

            #' Extract total ion chromatograms from csv files containing mass spectral data
            #'
            #' @param paths_to_cdf_csvs Paths to the cdf.csvs
            #' @param path_to_existing_chromatograms_csv Path to an existing chromatograms.csv file. The function will append new chromatograms to this file that are not already in it.
            #' @examples
            #' @export
            #' extractChromatogramsFromCSVs

                extractChromatogramsFromCSVs <- function( paths_to_cdf_csvs, path_to_existing_chromatograms_csv = NULL ) {

                    ## If the chromatograms.csv already exists, figure out which chromatograms are already in there and don't analyze those
                        if ( is.null(path_to_existing_chromatograms_csv) == FALSE ) {
                            paths_to_cdf_csvs <- paths_to_cdf_csvs[!paths_to_cdf_csvs %in% unique(readMonolist(path_to_existing_chromatograms_csv)$path_to_cdf_csv)]
                        }

                    ## Extract the chromatograms from the full dataset
                        chromatograms <- list()
                        for ( file in 1:length(paths_to_cdf_csvs) ) {

                            cat(paste("Reading data file ", paths_to_cdf_csvs[file], "\n", sep = ""))
                                framedDataFile <- as.data.frame(data.table::fread(paths_to_cdf_csvs[file]))

                            cat("   Extracting the total ion chromatogram...\n")
                                # library(plyr)
                                # chromatogram <- plyr::ddply(framedDataFile, .(rt), summarize, tic = sum(intensity))
                                framedDataFile %>% group_by(rt) %>% summarize(tic = sum(intensity)) -> chromatogram
                                chromatogram$rt <- as.numeric(chromatogram$rt)
                                chromatogram$path_to_cdf_csv <- paths_to_cdf_csvs[file]

                            cat("   Appending chromatogram to list...\n")
                                chromatograms[[file]] <- chromatogram
                        }

                        chromatograms <- do.call(rbind, chromatograms)

                    ## If the chromatograms file already exists, append to it, otherwise, return the chromatograms 
                        if ( is.null(path_to_existing_chromatograms_csv) == FALSE ) {
                            writeMonolist(
                                monolist = rbind(
                                    readMonolist(path_to_existing_chromatograms_csv),
                                    chromatograms
                                ),
                                monolist_out_path = path_to_existing_chromatograms_csv
                            )
                        } else {
                            return( chromatograms )        
                        }                
                }


        #### mergePeakLists

            #' Merge multiple peaklists into a single data frame
            #'
            #' This function reads multiple peak list files from a specified directory, processes them, and merges them into a single data frame. Each peak list is expected to be stored in a CSV file named "peaks_monolist.csv" within subdirectories of the specified directory.
            #'
            #' @param analysis_directory_path A character string specifying the path to the directory containing subdirectories with peak list files. Each subdirectory should contain a "peaks_monolist.csv" file.
            #' @return A data frame containing the merged peak lists. The data frame includes the first nine columns of each peak list, with the 'area' column converted to numeric.
            #' @examples
            #' # Example usage:
            #' # Merge peak lists from a specified directory
            #' merged_peaks <- mergePeakLists("/path/to/analysis_directory")
            #'
            #' @export
            mergePeakLists <- function(analysis_directory_path) {

                # Construct paths to each "peaks_monolist.csv" file within the specified directory
                paths_to_peak_monolists <- paste(
                    analysis_directory_path,
                    dir(analysis_directory_path),
                    "peaks_monolist.csv",
                    sep = "/"
                )

                peak_list <- list()
                # Iterate over each path to read and process the peak list
                for (i in 1:length(paths_to_peak_monolists)) {
                    
                    if (file.exists(paths_to_peak_monolists[i])) {

                        # Read the peak list and ensure the 'area' column is numeric
                        data <- readMonolist(paths_to_peak_monolists[i])
                        data$area <- as.numeric(data$area)

                        # If the data frame is not empty, store the first nine columns in the list
                        if (dim(data)[1] > 0) {
                            peak_list[[i]] <- data[,1:9]
                        }
                    
                    }
                
                }

                # Combine all the peak lists into a single data frame
                do.call(rbind, peak_list)

            }


        #### readChromatograms

            #' Import chromatograms stored as ChemStation exported .csv files
            #'
            #' Allows the user to import chromatograms stored as .csv files. Files should have originated from an Agilent GC system running ChemStation
            #' Files should be in the format of character~var1-var2.csv
            #' "character" is the name of the chromatogram, var1 and var2 are NUMERIC variables you want associated with that chromatogram.
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
                    if (normalize_level[1] == "var1") {
                        
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
                    if (normalize_level[1] == "var2") {
                        
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
            #' @param dir The directory containing the mass spectral (.csv) files of interest and nothing else
            #' @export
            #' @examples
            #' readSpectra()

            readSpectra <- function( dir, source = c("phylochemistry", "agilent") ) {

                ## Set working directory and prepare the data frame
                    setwd(dir)
                    data <- data.frame(
                        mz = as.numeric(1), abu = as.numeric(1), var1 = as.character("1"),
                        var2 = as.character("1"), var3 = as.character("1"), var4 = as.character("1")
                    )

                ## Read in all data
                    for (i in 1:length(dir())) {
                        temp <- readMonolist(dir()[i])
                        if (source[1] == "phylochemistry" ) {
                            mz <- as.numeric(as.character(temp$mz))
                            abu <- 100*as.numeric(as.character(temp$abu))/(max(as.numeric(as.character(temp$abu))))
                        }
                        if (source[1] == "agilent") {
                            temp <- temp[3:dim(temp)[1],1:2]
                            mz <- as.numeric(as.character(temp[,1]))
                            abu <- 100*as.numeric(as.character(temp[,2]))/(max(as.numeric(as.character(temp[,2]))))    
                        }
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


        #### analyzeMassSpectralImages

            #' Use image analysis to obtain csv of mass spectrum from an image.
            #'
            #' Read a png with left-most pixel low mass, right-most pixel high mass
            #' Filename: reference~compoundName_lowmz-rightmz
            #' @param image_directory_in_path The directory containing the images to be read.
            #' @export
            #' @examples
            #' analyzeMassSpectralImages()

            analyzeMassSpectralImages <- function(image_directory_in_path) {

                images_to_analyze <- dir(image_directory_in_path)[grep(".png", dir(image_directory_in_path))]
                
                ## Skip comparison images, if necessary
                    if (length(grep("comparison.png", images_to_analyze)) > 0) {
                        images_to_analyze <- images_to_analyze[-c(grep("comparison.png", images_to_analyze))]
                    }

                merged_spectral_output <- list()
                for (image_to_analyze in 1:length(images_to_analyze)) {

                    ## Load and rotate image, filter for just one color channel

                        cat(paste0("\n", "Analyzing ", images_to_analyze[image_to_analyze], "\n", "\n"))

                        left_mz <- as.numeric(gsub("-.*$", "", gsub(".*_", "", images_to_analyze[image_to_analyze])))
                        right_mz <- as.numeric(gsub("\\..*$", "", gsub(".*-", "", images_to_analyze[image_to_analyze])))
                        reference <- gsub("~.*$", "", images_to_analyze[image_to_analyze])

                        image <- imager::load.image(paste0(image_directory_in_path, "/", images_to_analyze[image_to_analyze]))
                        image <- as_tibble(as.data.frame(as.cimg(image)))
                        image <- image[image$cc == 1,]
                        image$y <- -as.numeric(image$y)
                        image$y <- image$y + -min(image$y)
                        
                        # ggplot(image, aes(x = x, y = y, fill = value)) + geom_tile()

                    ## Set to monochrome, normalize coordinates, sort output
                        
                        filtered_image <- filter(image, value < 0.8)
                        filtered_image$value <- 1
                        filtered_image$x <- normalize(
                        filtered_image$x, 
                        old_min = min(filtered_image$x),
                        old_max = max(filtered_image$x),
                        new_min = left_mz,
                        new_max = right_mz
                    )
                    filtered_image <- filtered_image[order(filtered_image$x),]

                    # ggplot(filtered_image, aes(x = x, y = y)) + geom_point(size = 0.1) + coord_cartesian(xlim = c(330,350))

                ## Evaluate continuity and create processed image, do some clean up with median and negative numbers

                    processed_image_y <- list()
                    processed_image_x <- list()
                    pb <- progress::progress_bar$new(total = length(unique(filtered_image$x)))
                    for( i in 1:length(unique(filtered_image$x)) ) {
                        
                        data <- filtered_image[filtered_image$x == unique(filtered_image$x)[i],]
                        
                        # # One way to find gaps
                        # if (any(!diff(data$y) == -1)) {
                        #   start <- data$y[sum(head(rle(diff(data$y))$lengths, -1))+1]
                        # } else {
                        #   start <- max(data$y)
                        # }

                        # Another way.. both are slow
                        start <- min(data$y)
                        while((start + 1) %in% data$y) {
                            start <- start + 1
                        }

                        processed_image_y[[i]] <- start
                        processed_image_x[[i]] <- unique(filtered_image$x)[i]
                        pb$tick()
                    }

                    processed_image <- data.frame(
                        x = do.call(rbind, processed_image_x),
                        y = do.call(rbind, processed_image_y)
                    )

                    processed_image$y <- processed_image$y - median(processed_image$y)
                    processed_image$y[processed_image$y < 0] <- 0

                    # ggplot(processed_image, aes(x = x, y = y)) + geom_col()

                ## Bin x coordinates to nominal resolution, set negatives to zero

                    processed_image$x_round <- round(processed_image$x, 0)
                    processed_image <- group_by(processed_image, x_round)
                    processed_image <- summarize(processed_image, y = sum(y))
                    processed_image$y <- processed_image$y / max(processed_image$y) * 100
                    colnames(processed_image)[colnames(processed_image) == "x_round"] <- "x"
                    

                    # ggplot(processed_image, aes(x = x, y = y)) + geom_col()

                ## Deal with duplicate peaks by resolving them in the way that creates the highest number of odd peaks

                    badness <- 0
                    for( i in 1:(length(processed_image$x)-1) ) {
                        if (processed_image$y[i] > 0 ) {
                            if( processed_image$y[i]/processed_image$y[i+1] > 0.9 & processed_image$y[i]/processed_image$y[i+1] < 1.1 ) {
                                if (processed_image$y[i] > 25) {
                                    badness <- badness + 1
                                }
                            }
                        }
                    }
                    # warning(paste0("Badness: ", badness, "\n"))

                ## Return

                    colnames(processed_image) <- c("mz", "relative_abundance")
                    write_csv(processed_image, paste0(image_directory_in_path, "/", gsub(".png", "", images_to_analyze[image_to_analyze]), ".csv"))

                ## Make plot of image versus processed data

                    plot_output <- ggplot() +
                        geom_col(data = processed_image, aes(x = mz, y = -relative_abundance), width = 0.8, color = "black") +
                        geom_text(data = 
                            processed_image[apply(cbind(
                                processed_image$relative_abundance > 2,
                                processed_image$mz > max(processed_image$mz)*0.999,
                                processed_image$mz < min(processed_image$mz)*1.001
                            ), 1, any),],
                            aes(x = mz, y = -relative_abundance - 5, label = mz)
                        ) +
                        geom_point(data = filtered_image, aes(x = x, y = y/max(y)*100), size = 0.1) +
                        theme_bw()
                    print(plot_output)
                    ggsave(filename = paste0(image_directory_in_path, "/", gsub(".png", "", images_to_analyze[image_to_analyze]), "_comparison.png"), plot = plot_output)

                ## Add to merged_spectral_output

                    merged_spectral_output[[image_to_analyze]] <- data.frame(
                        reference = reference,
                        compound_name = gsub("_.*$", "", gsub(".*~", "", gsub(".png", "", images_to_analyze[image_to_analyze]))),
                        mz = processed_image$mz,
                        relative_abundance = processed_image$relative_abundance
                    )
                }

                    ## Write out merged_spectral_output

                        merged_spectral_output <- do.call(rbind, merged_spectral_output)
                        write_csv(merged_spectral_output, paste0(image_directory_in_path, "/_merged_spectral_output.csv"))
            }


        #### crystallinity

            #' Extract 720 and 720 wavenumber intensities from FTIR data
            #'
            #' Filenames: 123~123-1.csv
            #' To use plots: library(gridExtra); nCol <- floor(sqrt(length(plots))); do.call("grid.arrange", c(plots, ncol=nCol))
            #' @param data_directory_in_path The directory containing the data to be read and processed.
            #' @export
            #' @examples
            #' crystallinity()

            crystallinity <- function(data_directory_in_path, make_plots = TRUE) {

                data_directory_in_path <- OsDirectoryPathCorrect(data_directory_in_path)

                files <- dir(data_directory_in_path)
                out <- list()

                if (make_plots == TRUE) {
                    plots <- list()
                }
                
                for (i in 1:length(files)) {

                    cat(paste0("File ", files[i], "\n"))
                    spectrum <- read_csv(paste0(data_directory_in_path, files[[i]]), col_names = FALSE, col_types = cols())
                    colnames(spectrum) <- c("x", "y")

                    spectrum <- dplyr::filter(spectrum, x < 800 & x > 600)

                    ## Baseline the signal using windows. Doesn't work as well as window-defined line below
                        # baseline <- findBaseline(
                        #   x = spectrum$x,
                        #   y = spectrum$y,
                        #   window_width = 80
                        # ) 
                        # spectrum <- spectrum[spectrum$x %in% baseline$x,]
                        # spectrum$b <- baseline$y

                    ## Baseline the signal using a window-defined line
                        baseline_guide <- rbind(
                            dplyr::filter(spectrum, x < 760 & x > 740),
                            dplyr::filter(spectrum, x < 710 & x > 680)
                        )
                        spectrum$b <- lm(y~x, data = baseline_guide)$coefficients[2]*spectrum$x + lm(y~x, data = baseline_guide)$coefficients[1]
                    
                    ## Subtract the baseline, smooth, and plot the signal
                    ## If you dont filter, you will get extra little crappy peaks when you try to model!!
                        spectrum$y <- signal::sgolayfilt(spectrum$y, p = 1, n = 11)
                        spectrum$yb <- spectrum$y - spectrum$b

                    ## filter and just find maxima
                        peak1 <- max(dplyr::filter(spectrum, x > 728 & x < 740)$yb)
                        peak2 <- max(dplyr::filter(spectrum, x > 715 & x < 722)$yb)
                        wave1 <- dplyr::filter(spectrum, x > 728 & x < 740)[dplyr::filter(spectrum, x > 728 & x < 740)$yb == peak1,]$x
                        wave2 <- dplyr::filter(spectrum, x > 715 & x < 722)[dplyr::filter(spectrum, x > 715 & x < 722)$yb == peak2,]$x
                        ratio <- peak1/peak2
                        print(ratio)

                    ## Plot to inspect

                        if (make_plots == TRUE) {
                            plots[[i]] <- ggplot(data = spectrum) +
                                # geom_line(aes(x = x, y = y)) +
                                # geom_point(aes(x = x, y = y)) +
                                # geom_line(aes(x = x, y = b)) +
                                geom_vline(xintercept = wave1) +
                                geom_vline(xintercept = wave2) +
                                scale_x_continuous(breaks = seq(0,1000,10)) +
                                geom_line(aes(x = x, y = yb)) +
                                coord_cartesian(xlim = c(800, 600))
                        }

                    ## Filter down and fit gaussians. Just use heights. It's more reliable
                        # spectrum <- dplyr::filter(spectrum, x > 705 & x < 745)
                        # gaussians <- fitGaussians(
                        #   x = spectrum$x, 
                        #   y = spectrum$yb,
                        #   peak_detection_threshold = as.numeric(gsub(".CSV", "", gsub(".*~", "", files[i])))
                        # )

                    ## Get areas
                        # gaussians$file <- files[i]
                        # gaussians %>%
                        #   as_tibble() %>%
                        #   group_by(peak_number) %>%
                        #   summarize(
                        #       wavenumber = x[which.max(y)],
                        #       area = unique(peak_area),
                        #       file = unique(file),
                        #       height = max(y)
                        #   ) -> out[[i]]

                    cat("\n")
                    # readline(prompt="Press [enter] to continue")

                    out[[i]] <- 
                        data.frame(
                            sample = files[i],
                            ratio = ratio
                        )
                }

                plots <<- plots

                out <- do.call(rbind, out)

                out$rep <- gsub(".CSV", "", gsub(".*-", "", out$sample))
                out$name <- gsub("-.*$", "", gsub("", "", out$sample))

                out %>%
                    group_by(name) %>%
                    dplyr::summarize(mean_ratio = mean(ratio)) -> ratio_order

                out$name <- factor(
                    out$name, 
                    levels = ratio_order$name[order(ratio_order$mean_ratio, decreasing = TRUE)]
                )

                return(out)

            }


    ##### Image Analysis

        #### analyzeImages

            #' Create major and minor axes ticks and labels
            #'
            #' @param major_ticks Values at which major ticks should be generated
            #' @param minor_freq Number of minor ticks between each major tick 
            #' @examples
            #' @export
            #' analyzeImages

            analyzeImages <- function(
                                folder_URL,
                                monolist_out_path
                            ) {

                files <- googledrive::drive_ls(folder_URL)
                file <- 1
                prop <- 0.4
                corner_points_x <<- vector()
                corner_points_y <<- vector()

                ui <- fluidPage(
                    verticalLayout(
                        splitLayout(
                            verticalLayout(
                                splitLayout(
                                    sliderInput("file", "File Number:",
                                        min = 1, max = dim(files)[1],
                                        value = 1, step = 1
                                    ),
                                    actionButton("download_image", "Download Image")
                                ),
                                verbatimTextOutput("photo_unique_id", placeholder = TRUE),
                                verbatimTextOutput("photo_metadata1", placeholder = TRUE),
                                verbatimTextOutput("photo_metadata2", placeholder = TRUE),
                                splitLayout(
                                    verticalLayout(
                                        sliderInput("rotation", "Rotation:",
                                            min = 0, max = 360,
                                            value = 0, step = 90
                                        ),
                                        verbatimTextOutput('selected')
                                    ),
                                    verticalLayout(
                                        actionButton("draw_image", "Draw Image"),
                                        actionButton("colorchecker_select", "ColorChecker Selected. Calibrate."),
                                        actionButton("rgb_select", "RGB Selected. Calibrate."),
                                        actionButton("reset", "Reset.")
                                    )
                                )
                            ),
                            splitLayout(
                                verticalLayout(
                                    actionButton("write_out", "Append to Spreadsheet")
                                ),
                                plotOutput(
                                    outputId = "sample_window",
                                    height = 400
                                )
                            )
                        ),
                        splitLayout(
                            plotOutput(
                                outputId = "photograph",
                                brush = brushOpts(
                                    id = "photograph_brush"
                                ),
                                click = "photograph_click", 
                                dblclick = "photograph_double_click",
                                height = 500
                            ),
                            pickerOutput(
                                outputId = "calibrated_photograph",
                                width = '100%',
                                height = "500px"
                            )
                        )
                    )
                )

                server <- function(input, output, session) {

                    ## Download image

                        observeEvent(input$download_image, {

                            googledrive::drive_download(
                                file = files$id[input$file],
                                path = "temp.JPEG",
                                overwrite = TRUE
                            )

                            path_to_image <- "temp.JPEG"

                            imRed <<- imager::load.image(path_to_image)

                            ## Extract and print QR code data and timestamp

                                photo_metadata <- list()
                                qr_location_data <- quadrangle::qr_scan_cpp(magick::image_read(path_to_image), lighten = TRUE, darken = TRUE)
                                photo_metadata$photo_unique_id <- as.character(files$id[input$file])
                                # photo_metadata$photo_unique_id <- as.character(files$id[file])
                                photo_metadata$name <- quadrangle::qr_scan_js_from_corners(magick::image_read(path_to_image), qr_location_data$points)$data
                                if( length(exifr::read_exif(path_to_image)$SubSecCreateDate) > 0 ){ photo_metadata$time_stamp <- exifr::read_exif(path_to_image)$SubSecCreateDate } else { photo_metadata$time_stamp <- exifr::read_exif(path_to_image)$FileModifyDate }
                                output$photo_unique_id <- renderText({ photo_metadata$photo_unique_id })
                                output$photo_metadata1 <- renderText({ photo_metadata$name })
                                output$photo_metadata2 <- renderPrint({ photo_metadata$time_stamp })
                                photo_metadata$date <- gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp))
                                photo_metadata$year <- lubridate::year(gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp)))
                                photo_metadata$month <- lubridate::month(gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp)))
                                photo_metadata$day <- lubridate::day(gsub(":", "/", gsub(" .*$", "", photo_metadata$time_stamp)))
                                photo_metadata <<- photo_metadata
                        })

                    ## Draw or rotate photograph when asked

                        observeEvent(input$draw_image, {

                            imRed <<- imager::imrotate(imRed, input$rotation)
                            
                            image_as_df <- as_tibble(as.data.frame(as.cimg(imRed)))
                            image_as_df$y <- -as.numeric(image_as_df$y)
                            image_as_df$y <- image_as_df$y + -min(image_as_df$y)
                            image_as_df <- pivot_wider(image_as_df, names_from = cc, values_from = value)
                            names(image_as_df)[3:5] <- c("R", "G", "B")
                            image_as_df$hex <- rgb(image_as_df$R, image_as_df$G, image_as_df$B)
                            image_as_df <<- image_as_df

                            photograph <- ggplot() +
                                geom_tile(data = filter(image_as_df), aes(x = x, y = y), fill = image_as_df$hex) + theme_void()
                                photograph
                            output$photograph <- renderPlot({photograph})

                            corner_points_x <<- vector()
                            corner_points_y <<- vector()

                            # output$selected <- isolate(renderPrint({
                            #     corner_points_x
                            # }))

                        })

                    ## Extract sample spot

                        output$sample_window <- renderPlot({ 

                            if ( !is.null(input$calibrated_photograph_selected_points) ) {

                                sample_pixel_data <<- image_calibrated[input$calibrated_photograph_selected_points,]

                                ggplot() +
                                    geom_tile(
                                        data = sample_pixel_data,
                                        aes(x = x, y = y),
                                        fill = sample_pixel_data$hex
                                    ) +
                                    theme_void()

                                # sample_pixel_output <<- data.frame(
                                #     photo_unique_id = photo_metadata$photo_unique_id,
                                #     name = photo_metadata$name,
                                #     time_stamp = photo_metadata$time_stamp,
                                #     x = sample_pixel_data$x,
                                #     y = sample_pixel_data$y,
                                #     color = sample_pixel_data$hex
                                # )

                            } else { NULL}

                        })
                                    
                    ## Run ColorChecker correction

                        output$selected <- renderPrint({
                            corner_points_x <<- c(corner_points_x, input$photograph_click$x)
                            corner_points_y <<- c(corner_points_y, input$photograph_click$y)
                            corner_points_x
                        })

                        # observeEvent(input$reset, {
                            

                        #     output$selected <- renderPrint({
                        #         corner_points_x
                        #     })
                        # })

                        observeEvent(input$colorchecker_select, {

                            ## Pre-calculations
                                
                                mR <- raster::as.matrix(imager::R(imRed)) * 255
                                mG <- raster::as.matrix(imager::G(imRed)) * 255
                                mB <- raster::as.matrix(imager::B(imRed)) * 255
                                rR <- raster::raster(mR)
                                rG <- raster::raster(mG)
                                rB <- raster::raster(mB)
                                raster::extent(rR) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
                                raster::extent(rG) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
                                raster::extent(rB) <- c(0, dim(imRed)[2], 0, dim(imRed)[1])
                                rR <- raster::flip(t(rR), "y")
                                rG <- raster::flip(t(rG), "y")
                                rB <- raster::flip(t(rB), "y")

                            ## Get the border of the colorchecker and swatch positions

                                xy <- list()
                                xy$x <- corner_points_x
                                xy$y <- corner_points_y

                                xyDF <- as.data.frame(xy)
                                line1 <- xyDF[1:2, ]
                                line2 <- xyDF[2:3, ]
                                line3 <- xyDF[3:4, ]
                                line4 <- xyDF[c(4, 1), ]
                                xdiff <- (line1$x[2] - line1$x[1])/6
                                ydiff <- (line1$y[2] - line1$y[1])/6
                                xmin <- line1$x[1]
                                ymin <- line1$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2, (xmin + xdiff * 4 + xmin + 
                                    xdiff * 5)/2, (xmin + xdiff * 5 + xmin + xdiff * 
                                    6)/2), y = c((ymin + ymin + ydiff)/2, (ymin + 
                                    ydiff + ymin + ydiff * 2)/2, (ymin + ydiff * 
                                    2 + ymin + ydiff * 3)/2, (ymin + ydiff * 3 + 
                                    ymin + ydiff * 4)/2, (ymin + ydiff * 4 + ymin + 
                                    ydiff * 5)/2, (ymin + ydiff * 5 + ymin + ydiff * 
                                    6)/2))
                                xySubDF_1A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop, xmin + xdiff * 5 - (xdiff/2) * 
                                    prop, xmin + xdiff * 6 - (xdiff/2) * prop), y = c(ymin + 
                                    ydiff - (ydiff/2) * prop, ymin + ydiff * 2 - 
                                    (ydiff/2) * prop, ymin + ydiff * 3 - (ydiff/2) * 
                                    prop, ymin + ydiff * 4 - (ydiff/2) * prop, ymin + 
                                    ydiff * 5 - (ydiff/2) * prop, ymin + ydiff * 
                                    6 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop, xmin + xdiff * 4 + (xdiff/2) * prop, xmin + 
                                    xdiff * 5 + (xdiff/2) * prop), y = c(ymin + (ydiff/2) * 
                                    prop, ymin + ydiff + (ydiff/2) * prop, ymin + 
                                    ydiff * 2 + (ydiff/2) * prop, ymin + ydiff * 
                                    3 + (ydiff/2) * prop, ymin + ydiff * 4 + (ydiff/2) * 
                                    prop, ymin + ydiff * 5 + (ydiff/2) * prop))
                                xySubDF_1Ba <- as.data.frame(xySubBa)
                                xySubDF_1Bb <- as.data.frame(xySubBb)
                                xdiff <- (line3$x[2] - line3$x[1])/6
                                ydiff <- (line3$y[2] - line3$y[1])/6
                                xmin <- line3$x[1]
                                ymin <- line3$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2, (xmin + xdiff * 4 + xmin + 
                                    xdiff * 5)/2, (xmin + xdiff * 5 + xmin + xdiff * 
                                    6)/2), y = c((ymin + ymin + ydiff)/2, (ymin + 
                                    ydiff + ymin + ydiff * 2)/2, (ymin + ydiff * 
                                    2 + ymin + ydiff * 3)/2, (ymin + ydiff * 3 + 
                                    ymin + ydiff * 4)/2, (ymin + ydiff * 4 + ymin + 
                                    ydiff * 5)/2, (ymin + ydiff * 5 + ymin + ydiff * 
                                    6)/2))
                                xySubDF_3A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop, xmin + xdiff * 5 - (xdiff/2) * 
                                    prop, xmin + xdiff * 6 - (xdiff/2) * prop), y = c(ymin + 
                                    ydiff - (ydiff/2) * prop, ymin + ydiff * 2 - 
                                    (ydiff/2) * prop, ymin + ydiff * 3 - (ydiff/2) * 
                                    prop, ymin + ydiff * 4 - (ydiff/2) * prop, ymin + 
                                    ydiff * 5 - (ydiff/2) * prop, ymin + ydiff * 
                                    6 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop, xmin + xdiff * 4 + (xdiff/2) * prop, xmin + 
                                    xdiff * 5 + (xdiff/2) * prop), y = c(ymin + (ydiff/2) * 
                                    prop, ymin + ydiff + (ydiff/2) * prop, ymin + 
                                    ydiff * 2 + (ydiff/2) * prop, ymin + ydiff * 
                                    3 + (ydiff/2) * prop, ymin + ydiff * 4 + (ydiff/2) * 
                                    prop, ymin + ydiff * 5 + (ydiff/2) * prop))
                                xySubDF_3Ba <- as.data.frame(xySubBa)
                                xySubDF_3Bb <- as.data.frame(xySubBb)
                                xdiff <- (line2$x[2] - line2$x[1])/4
                                ydiff <- (line2$y[2] - line2$y[1])/4
                                xmin <- line2$x[1]
                                ymin <- line2$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2), y = c((ymin + ymin + ydiff)/2, 
                                    (ymin + ydiff + ymin + ydiff * 2)/2, (ymin + 
                                      ydiff * 2 + ymin + ydiff * 3)/2, (ymin + ydiff * 
                                      3 + ymin + ydiff * 4)/2))
                                xySubDF_2A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop), y = c(ymin + ydiff - (ydiff/2) * 
                                    prop, ymin + ydiff * 2 - (ydiff/2) * prop, ymin + 
                                    ydiff * 3 - (ydiff/2) * prop, ymin + ydiff * 
                                    4 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop), y = c(ymin + (ydiff/2) * prop, ymin + 
                                    ydiff + (ydiff/2) * prop, ymin + ydiff * 2 + 
                                    (ydiff/2) * prop, ymin + ydiff * 3 + (ydiff/2) * 
                                    prop))
                                xySubDF_2Ba <- as.data.frame(xySubBa)
                                xySubDF_2Bb <- as.data.frame(xySubBb)
                                xdiff <- (line4$x[2] - line4$x[1])/4
                                ydiff <- (line4$y[2] - line4$y[1])/4
                                xmin <- line4$x[1]
                                ymin <- line4$y[1]
                                xySubA <- list(x = c((xmin + xmin + xdiff)/2, (xmin + 
                                    xdiff + xmin + xdiff * 2)/2, (xmin + xdiff * 
                                    2 + xmin + xdiff * 3)/2, (xmin + xdiff * 3 + 
                                    xmin + xdiff * 4)/2), y = c((ymin + ymin + ydiff)/2, 
                                    (ymin + ydiff + ymin + ydiff * 2)/2, (ymin + 
                                      ydiff * 2 + ymin + ydiff * 3)/2, (ymin + ydiff * 
                                      3 + ymin + ydiff * 4)/2))
                                xySubDF_4A <- as.data.frame(xySubA)
                                xySubBa <- list(x = c(xmin + xdiff - (xdiff/2) * 
                                    prop, xmin + xdiff * 2 - (xdiff/2) * prop, xmin + 
                                    xdiff * 3 - (xdiff/2) * prop, xmin + xdiff * 
                                    4 - (xdiff/2) * prop), y = c(ymin + ydiff - (ydiff/2) * 
                                    prop, ymin + ydiff * 2 - (ydiff/2) * prop, ymin + 
                                    ydiff * 3 - (ydiff/2) * prop, ymin + ydiff * 
                                    4 - (ydiff/2) * prop))
                                xySubBb <- list(x = c(xmin + (xdiff/2) * prop, xmin + 
                                    xdiff + (xdiff/2) * prop, xmin + xdiff * 2 + 
                                    (xdiff/2) * prop, xmin + xdiff * 3 + (xdiff/2) * 
                                    prop), y = c(ymin + (ydiff/2) * prop, ymin + 
                                    ydiff + (ydiff/2) * prop, ymin + ydiff * 2 + 
                                    (ydiff/2) * prop, ymin + ydiff * 3 + (ydiff/2) * 
                                    prop))
                                xySubDF_4Ba <- as.data.frame(xySubBa)
                                xySubDF_4Bb <- as.data.frame(xySubBb)

                                sub_data <<- rbind(xySubDF_1A, xySubDF_1Ba, xySubDF_1Bb, xySubDF_2A, xySubDF_2Ba, xySubDF_2Bb, xySubDF_3A, xySubDF_3Ba, xySubDF_3Bb, xySubDF_4A, xySubDF_4Ba, xySubDF_4Bb)

                                labels <- list(
                                    c(1, 7, 13, 19), c(2, 8, 14, 20), 
                                    c(3, 9, 15, 21), c(4, 10, 16, 22),
                                    c(5, 11, 17, 23), c(6, 12, 18, 24)
                                )

                            ## Get swatch positions and labels
                                
                                swatch_labels <- list()
                                xyTot <- c()
                                for (e in 1:nrow(xySubDF_1A)) {
                                    xyLine1 <- xySubDF_1A[e, ]
                                    xyLine2 <- xySubDF_3A[6:1, ][e, ]
                                    xdiff <- (xyLine2$x - xyLine1$x)/4
                                    ydiff <- (xyLine2$y - xyLine1$y)/4
                                    xmin <- xyLine1$x
                                    ymin <- xyLine1$y
                                    xySub <- list(x = c((xmin + xmin + xdiff)/2, 
                                      (xmin + xdiff + xmin + xdiff * 2)/2, (xmin + 
                                        xdiff * 2 + xmin + xdiff * 3)/2, (xmin + 
                                        xdiff * 3 + xmin + xdiff * 4)/2), y = c((ymin + 
                                      ymin + ydiff)/2, (ymin + ydiff + ymin + ydiff * 
                                      2)/2, (ymin + ydiff * 2 + ymin + ydiff * 3)/2, 
                                      (ymin + ydiff * 3 + ymin + ydiff * 4)/2))
                                    xySubDF <- as.data.frame(xySub)

                                    swatch_labels[[e]] <- xySubDF
                                    
                                    xySubDFLabel <- cbind(xySubDF, label = labels[[e]])
                                    xyLine1a <- xySubDF_1Ba[e, ]
                                    xyLine1b <- xySubDF_1Bb[e, ]
                                    xyLine3a <- xySubDF_3Bb[6:1, ][e, ]
                                    xyLine3b <- xySubDF_3Ba[6:1, ][e, ]
                                    xdiffa <- (xyLine3a$x - xyLine1a$x)/4
                                    ydiffa <- (xyLine3a$y - xyLine1a$y)/4
                                    xdiffb <- (xyLine3b$x - xyLine1b$x)/4
                                    ydiffb <- (xyLine3b$y - xyLine1b$y)/4
                                    xmina <- xyLine1a$x
                                    ymina <- xyLine1a$y
                                    xminb <- xyLine1b$x
                                    yminb <- xyLine1b$y
                                    xySubAa <- list(x = c(xmina + xdiffa - (xdiffa/2) * 
                                      prop, xmina + xdiffa * 2 - (xdiffa/2) * prop, 
                                      xmina + xdiffa * 3 - (xdiffa/2) * prop, xmina + 
                                        xdiffa * 4 - (xdiffa/2) * prop), y = c(ymina + 
                                      ydiffa - (ydiffa/2) * prop, ymina + ydiffa * 
                                      2 - (ydiffa/2) * prop, ymina + ydiffa * 3 - 
                                      (ydiffa/2) * prop, ymina + ydiffa * 4 - (ydiffa/2) * 
                                      prop))
                                    xySubAb <- list(x = c(xmina + (xdiffa/2) * prop, 
                                      xmina + xdiffa + (xdiffa/2) * prop, xmina + 
                                        xdiffa * 2 + (xdiffa/2) * prop, xmina + xdiffa * 
                                        3 + (xdiffa/2) * prop), y = c(ymina + (ydiffa/2) * 
                                      prop, ymina + ydiffa + (ydiffa/2) * prop, ymina + 
                                      ydiffa * 2 + (ydiffa/2) * prop, ymina + ydiffa * 
                                      3 + (ydiffa/2) * prop))
                                    xySubBa <- list(x = c(xminb + xdiffb - (xdiffb/2) * 
                                      prop, xminb + xdiffb * 2 - (xdiffb/2) * prop, 
                                      xminb + xdiffb * 3 - (xdiffb/2) * prop, xminb + 
                                        xdiffb * 4 - (xdiffb/2) * prop), y = c(yminb + 
                                      ydiffb - (ydiffb/2) * prop, yminb + ydiffb * 
                                      2 - (ydiffb/2) * prop, yminb + ydiffb * 3 - 
                                      (ydiffb/2) * prop, yminb + ydiffb * 4 - (ydiffb/2) * 
                                      prop))
                                    xySubBb <- list(x = c(xminb + (xdiffb/2) * prop, 
                                      xminb + xdiffb + (xdiffb/2) * prop, xminb + 
                                        xdiffb * 2 + (xdiffb/2) * prop, xminb + xdiffb * 
                                        3 + (xdiffb/2) * prop), y = c(yminb + (ydiffb/2) * 
                                      prop, yminb + ydiffb + (ydiffb/2) * prop, yminb + 
                                      ydiffb * 2 + (ydiffb/2) * prop, yminb + ydiffb * 
                                      3 + (ydiffb/2) * prop))
                                    xySubAaDF <- as.data.frame(xySubAa)
                                    xySubAbDF <- as.data.frame(xySubAb)
                                    xySubBaDF <- as.data.frame(xySubBa)
                                    xySubBbDF <- as.data.frame(xySubBb)
                                    xySubRow <- cbind(xySubAaDF, xySubAbDF, xySubBaDF, 
                                      xySubBbDF, label = labels[[e]])
                                    colnames(xySubRow) <- c("x1", "y1", "x2", "y2", 
                                      "x3", "y3", "x4", "y4", "label")
                                    xyTot <- rbind(xyTot, xySubRow)
                                }

                                swatch_labels <- do.call(rbind, swatch_labels)
                                swatch_labels$label <- unlist(labels)

                            ## Get swatch colors

                                xyTot$imR <- NA
                                xyTot$imG <- NA
                                xyTot$imB <- NA
                                xyTot <- xyTot[order(xyTot$label), ]
                                polygons <- list()
                                for (e in 1:nrow(xyTot)) {
                                    
                                    polygons[[e]] <- data.frame(
                                        xmax = xyTot$x1[e],
                                        xmin = xyTot$x2[e],
                                        ymax = xyTot$y1[e],
                                        ymin = xyTot$y3[e],
                                        label = e
                                    )

                                    image_as_df %>%
                                        filter(
                                            x < xyTot$x1[e] &
                                            x > xyTot$x2[e] &
                                            y < xyTot$y1[e] &
                                            y > xyTot$y3[e]
                                        ) -> polygon_data

                                    xyTot$imR[e] <- mode(polygon_data$R) * 255
                                    xyTot$imG[e] <- mode(polygon_data$G) * 255
                                    xyTot$imB[e] <- mode(polygon_data$B) * 255

                                }
                                polygons <- do.call(rbind, polygons)


                                l1 <- c(1, 115, 82, 68)
                                l2 <- c(2, 194, 150, 130)
                                l3 <- c(3, 98, 122, 157)
                                l4 <- c(4, 87, 108, 67)
                                l5 <- c(5, 133, 128, 177)
                                l6 <- c(6, 103, 189, 170)
                                l7 <- c(7, 214, 126, 44)
                                l8 <- c(8, 80, 91, 166)
                                l9 <- c(9, 193, 90, 99)
                                l10 <- c(10, 94, 60, 108)
                                l11 <- c(11, 157, 188, 64)
                                l12 <- c(12, 224, 163, 46)
                                l13 <- c(13, 56, 61, 150)
                                l14 <- c(14, 70, 148, 73)
                                l15 <- c(15, 175, 54, 60)
                                l16 <- c(16, 231, 199, 31)
                                l17 <- c(17, 187, 86, 149)
                                l18 <- c(18, 8, 133, 161)
                                l19 <- c(19, 243, 243, 242)
                                l20 <- c(20, 200, 200, 200)
                                l21 <- c(21, 160, 160, 160)
                                l22 <- c(22, 122, 122, 121)
                                l23 <- c(23, 85, 85, 85)
                                l24 <- c(24, 52, 52, 52)
                                ColorCheckerRGB <- as.data.frame(rbind(l1, l2, l3, 
                                    l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, 
                                    l15, l16, l17, l18, l19, l20, l21, l22, l23, 
                                    l24)
                                )
                                colnames(ColorCheckerRGB) <- c("label", "sR", "sG", "sB")

                                dat <<- merge(xyTot, ColorCheckerRGB, by = "label")

                                ## Calculate regression for the colors

                                    print("Calculating polynomial regression...")
                                        sR <- dat$sR
                                        sG <- dat$sG
                                        sB <- dat$sB
                                        imR <- dat$imR
                                        imG <- dat$imG
                                        imB <- dat$imB
                                        modelR <- lm(sR ~ imR + imG + imB + imR^2 + imG^2 + 
                                            imB^2)
                                        modelG <- lm(sG ~ imR + imG + imB + imR^2 + imG^2 + 
                                            imB^2)
                                        modelB <- lm(sB ~ imR + imG + imB + imR^2 + imG^2 + 
                                            imB^2)
                                        dfIm = data.frame(imR = matrix(mR, ncol = 1), imG = matrix(mG, 
                                            ncol = 1), imB = matrix(mB, ncol = 1))

                                ## Calibrate colors

                                    print("Calibrating colors...")
                                        prR <- predict(modelR, dfIm)
                                        prG <- predict(modelG, dfIm)
                                        prB <- predict(modelB, dfIm)
                                        dfCal <- as.data.frame(cbind(prR, prG, prB))
                                    
                                    print("Rebuilding image...")
                                        Ri = matrix(dfCal$prR, nrow = dim(imRed)[1])
                                        Gi = matrix(dfCal$prG, nrow = dim(imRed)[1])
                                        Bi = matrix(dfCal$prB, nrow = dim(imRed)[1])
                                        imCal = array(dim = dim(imRed))
                                        imCal[, , , 1] = Ri
                                        imCal[, , , 2] = Gi
                                        imCal[, , , 3] = Bi
                                        imCal <- imager::as.cimg(imCal)
                                        
                                        image_calibrated <<- as_tibble(as.data.frame(as.cimg(imCal)))
                                        image_calibrated$y <- -as.numeric(image_calibrated$y)
                                        image_calibrated$y <- image_calibrated$y + -min(image_calibrated$y)
                                        image_calibrated <- pivot_wider(image_calibrated, names_from = cc, values_from = value)
                                        names(image_calibrated)[3:5] <- c("R", "G", "B")

                                        image_calibrated$R[image_calibrated$R/255 < 0] <- 0
                                        image_calibrated$G[image_calibrated$G/255 < 0] <- 0
                                        image_calibrated$B[image_calibrated$B/255 < 0] <- 0

                                        image_calibrated$R[image_calibrated$R/255 > 1] <- 1
                                        image_calibrated$G[image_calibrated$G/255 > 1] <- 1
                                        image_calibrated$B[image_calibrated$B/255 > 1] <- 1

                                        image_calibrated$hex <- rgb(image_calibrated$R/255, image_calibrated$G/255, image_calibrated$B/255)
                                        image_calibrated <<- image_calibrated

                                output$photograph <- renderPlot({
                                    
                                    ggplot() +
                                        geom_tile(data = filter(image_as_df), aes(x = x, y = y), fill = image_as_df$hex) + theme_void() +
                                        geom_point(
                                            data = sub_data,
                                            aes(x = x, y = y), shape = 21, size = 3, color = "black", fill = "maroon"
                                        ) +
                                        geom_text(
                                            data = swatch_labels,
                                            aes(x = x, y = y, label = label), color = "gold"
                                        ) +
                                        geom_rect(
                                            data = polygons,
                                            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "gold", alpha = 0
                                        )

                                })

                                output$calibrated_photograph <- renderPicker({

                                    picker(
                                        coords = as.data.frame(image_calibrated[,1:2]),
                                        colors = as.data.frame(image_calibrated[,6])[[1]],
                                        # cluster_colors = image_as_df[,6],
                                        labels = as.character(as.data.frame(image_calibrated[,1])[[1]]) 
                                        # labels = NULL
                                        # label_coords = label_coords,
                                        # polygons = polygons, 
                                        # text_props = text_props,
                                        # point_color_polygons = 'white',
                                        # grid_legend_items = grid_legend_items
                                    )
                                    
                                    # ggplot() +
                                    #     geom_tile(data = filter(image_calibrated), aes(x = x, y = y), fill = image_calibrated$hex) + theme_void()
                                        # geom_point(
                                        #     data = sub_data,
                                        #     aes(x = x, y = y), shape = 21, size = 3, color = "black", fill = "maroon"
                                        # ) +
                                        # geom_text(
                                        #     data = swatch_labels,
                                        #     aes(x = x, y = y, label = label), color = "gold"
                                        # ) +
                                        # geom_rect(
                                        #     data = polygons,
                                        #     aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = "gold"
                                        # )

                                })

                        })

                    ## Extract rgb_window

                        observeEvent(input$rgb_select, {

                            if ( isolate(!is.null(input$photograph_brush)) ) {
                                
                                image_points <<- isolate(brushedPoints(image_as_df, input$photograph_brush))
                                
                                Rmax <- cbind(
                                    data.frame(color_swatch = "red"),
                                    image_points[order(image_points$R, decreasing = TRUE),][1:50,]
                                )
                                Gmax <- cbind(
                                    data.frame(color_swatch = "green"),
                                    image_points[order(image_points$G, decreasing = TRUE),][1:50,]
                                )
                                Bmax <- cbind(
                                    data.frame(color_swatch = "blue"),
                                    image_points[order(image_points$B, decreasing = TRUE),][1:50,]
                                )
                                
                                rgb_window <- ggplot() +
                                    geom_tile(data = image_points, aes(x = x, y = y), fill = image_points$hex) +
                                    geom_rect(aes(xmin = min(Rmax$x), xmax = max(Rmax$x), ymin = min(Rmax$y), ymax = max(Rmax$y)), alpha = 0.5) +
                                    geom_rect(aes(xmin = min(Gmax$x), xmax = max(Gmax$x), ymin = min(Gmax$y), ymax = max(Gmax$y)), alpha = 0.5) +
                                    geom_rect(aes(xmin = min(Bmax$x), xmax = max(Bmax$x), ymin = min(Bmax$y), ymax = max(Bmax$y)), alpha = 0.5)

                                output$rgb_window <- renderPlot({ rgb_window })

                                swatch_data <<- rbind(
                                    select(Rmax, color_swatch, x, y, R) %>%
                                    set_colnames(c("color_swatch", "x", "y", "value")),
                                    select(Gmax, color_swatch, x, y, G) %>%
                                    set_colnames(c("color_swatch", "x", "y", "value")),
                                    select(Bmax, color_swatch, x, y, B) %>%
                                    set_colnames(c("color_swatch", "x", "y", "value"))
                                )
                                swatch_data$color_swatch <- factor(swatch_data$color_swatch, levels = c("red", "green", "blue"))

                                rgb_histogram <- ggplot(data = swatch_data) + 
                                    geom_histogram(aes(x = value, fill = color_swatch), binwidth = 0.01, color = "black") +
                                    scale_fill_manual(values = c("red", "green", "blue")) +
                                    facet_grid(color_swatch~.) +
                                    scale_x_continuous(limits = c(-0.1,1.1)) +
                                    theme_bw()

                                output$rgb_histogram <- renderPlot({ rgb_histogram })

                                # mode_reference_value_R <<- mode(swatch_data$value[swatch_data$color_swatch == "red"])
                                # mode_reference_value_G <<- mode(swatch_data$value[swatch_data$color_swatch == "green"])
                                # mode_reference_value_B <<- mode(swatch_data$value[swatch_data$color_swatch == "blue"])

                                mode_reference_value_R <<- mean(swatch_data$value[swatch_data$color_swatch == "red"])
                                mode_reference_value_G <<- mean(swatch_data$value[swatch_data$color_swatch == "green"])
                                mode_reference_value_B <<- mean(swatch_data$value[swatch_data$color_swatch == "blue"])

                                output$mode_reference_value_R <- renderText({ mode_reference_value_R })
                                output$mode_reference_value_G <- renderText({ mode_reference_value_G })
                                output$mode_reference_value_B <- renderText({ mode_reference_value_B })

                                reference_pixel_output <<- data.frame(
                                    QR = if(is.null(photo_metadata$name)) {"NA"} else {photo_metadata$name},
                                    Year = photo_metadata$year,
                                    Month = photo_metadata$month,
                                    Day = photo_metadata$day,
                                    Date = photo_metadata$date,
                                    type = "reference",
                                    x = image_points$x,
                                    y = image_points$y,
                                    R = image_points$R,
                                    G = image_points$G,
                                    B = image_points$B
                                )

                            } else { NULL }

                        })

                    ## Write out

                        ## Here be sure to use a unique_photo_ID and have the app check to see if its already in output before writing out!

                        ## Maybe a share link could be the unique photo ID

                        observeEvent(input$write_out, {

                            if (file.exists(monolist_out_path)) {

                                writeMonolist(
                                    cbind(as.data.frame(photo_metadata), sample_pixel_data),
                                    append = TRUE,
                                    monolist_out_path = monolist_out_path,
                                    col.names = FALSE
                                )
                                cat("Appended to Spreadsheet")

                            } else {

                                writeMonolist(
                                    cbind(as.data.frame(photo_metadata), sample_pixel_data),
                                    append = TRUE,
                                    monolist_out_path = monolist_out_path,
                                    col.names = TRUE
                                )
                                cat("Appended to Spreadsheet")

                            }
                                
                        })

                    ## Meta analysis tab

                            # observeEvent(input$metadata_file, {

                            #     ## Read in metadata file
                            #         metadata_file <<- input$metadata_file
                            #         print(metadata_file)
                            #         prelim_data <<- readMonolist(metadata_file$datapath)
                            #         print(head(prelim_data))
                            
                            #     ## Extract swatches
                            #         prelim_data %>%
                            #             pivot_longer(cols = 9:11, names_to = "color", values_to = "value") %>%
                            #             group_by(QR, color) %>%
                            #             arrange(desc(value), .by_group = TRUE) %>%
                            #             dplyr::slice(1:25) -> swatch_data

                            #     ## Determine global means
                            #         swatch_data$global_mean <- 0  
                            #         swatch_data$global_mean[swatch_data$color == "R"] <- mean(swatch_data$value[swatch_data$color == "R"])
                            #         swatch_data$global_mean[swatch_data$color == "G"] <- mean(swatch_data$value[swatch_data$color == "B"])
                            #         swatch_data$global_mean[swatch_data$color == "B"] <- mean(swatch_data$value[swatch_data$color == "G"])

                            #     ## Create dataframe according to which normalization should be performed
                            #         swatch_data %>%
                            #         group_by(QR, color) %>%
                            #         summarize(difference = unique(global_mean) - unique(mean(value))) %>%
                            #         # mutate(unique_id = paste(QR, "_", color)) %>%
                            #         pivot_wider(names_from = color, values_from = difference, names_prefix = "adjust_") -> normalization_data

                            #     ## Modify swatch data according to global means
                            #         swatch_data %>%
                            #         group_by(QR, color) %>%
                            #         mutate(modified_value = value + (global_mean - mean(value))) %>%
                            #         pivot_longer(cols = 10:12, names_to = "value_type", values_to = "value") %>%
                            #         ggplot() +
                            #             geom_histogram(aes(x = value, fill = value_type), bins = 30, alpha = 0.7, color = "black") +
                            #             facet_grid(QR~color) +
                            #             theme_bw() -> global_means_plot

                            #         output$metaplot <- renderPlot({ global_means_plot })

                            # })

                        ## Modify reference data according to global means
                            
                            # prelim_data %>%
                            #     mutate(
                            #         R_adjusted = R + normalization_data$adjust_R[match(QR, normalization_data$QR)],
                            #         G_adjusted = G + normalization_data$adjust_G[match(QR, normalization_data$QR)],
                            #         B_adjusted = B + normalization_data$adjust_B[match(QR, normalization_data$QR)]
                            #     ) -> prelim_data_processed

                            # prelim_data_processed$R_adjusted[prelim_data_processed$R_adjusted < 0] <- 0
                            # prelim_data_processed$G_adjusted[prelim_data_processed$G_adjusted < 0] <- 0
                            # prelim_data_processed$B_adjusted[prelim_data_processed$B_adjusted < 0] <- 0

                            # prelim_data_processed$R_adjusted[prelim_data_processed$R_adjusted > 1] <- 1
                            # prelim_data_processed$G_adjusted[prelim_data_processed$G_adjusted > 1] <- 1
                            # prelim_data_processed$B_adjusted[prelim_data_processed$B_adjusted > 1] <- 1

                            # prelim_data_processed %>%
                            # filter(type == "reference") %>%
                            # mutate(
                            # hex_before = rgb(R, G, B),
                            # hex_after = rgb(R_adjusted, G_adjusted, B_adjusted)
                            # ) %>%
                            # pivot_longer(cols = c("hex_before", "hex_after"), names_to = "state", values_to = "hex") %>%
                            # filter(QR %in% unique(prelim_data_processed$QR)) -> plot_data

                            # ggplot() +
                            # geom_tile(data = plot_data, aes(x = x, y = y), fill = plot_data$hex) +
                            # facet_grid(QR~state, space = "free", scales = "free") -> plot

                            # png(filename = "test.png", width = 10, height = 10, units = "in", res = 300)
                            # plot
                            # dev.off()

                            # ## Plot swatches before and after
                            #     ggplot() +
                            #         geom_tile(
                            #         data = filter(swatch_data, value_type == "value"),
                            #         aes(x = x, y = y, fill = hex)
                            #     )

                            # swatch_data %>%
                            #     group_by(QR, color) %>%
                            #     mutate(modified_value = value + (global_mean - mean(value))) %>%
                            #     pivot_longer(cols = 11:13, names_to = "value_type", values_to = "value") %>%
                            #     filter(value_type == "value")
                            #     ggplot() +
                            #     geom_tile(aes(x = x, y = y, fill = hex)) +
                            #     facet_grid(QR~color) +
                            #     theme_bw()

                            #     out <- list()
                            #     seq <- seq(0.6,0.8,0.001)
                            #     for (i in 1:length(seq)) {
                            #     test$black_v_white <- 0
                            #     test$black_v_white[test$greyscale > seq[i]] <- 1
                            #     out[[i]] <- data.frame(
                            #     ratio = 
                            #     as.data.frame(table(test$black_v_white))$Freq[1]/
                            #     as.data.frame(table(test$black_v_white))$Freq[2],
                            #     i = i
                            #     )
                            #     }
                            #     out <- do.call(rbind, out)

                            #     ggplot(out, aes(x = i, y = ratio)) + geom_col()



                            #     plot_a <- test %>%
                            #     ggplot(aes(x = x, y = y)) +
                            #     geom_tile(
                            #     aes(fill = greyscale),
                            #     show.legend = FALSE
                            #     ) +
                            #     scale_fill_gradient(low = "black", high = "white")

                            #     plot_b <- test %>%
                            #     ggplot(aes(x = x, y = y)) +
                            #     geom_tile(
                            #     aes(fill = black_v_white),
                            #     show.legend = FALSE
                            #     ) +
                            #     scale_fill_gradient(low = "black", high = "white")

                            #     plot_c <- ggplot(test) +
                            #     geom_histogram(aes(x = greyscale), bins = 80) +
                            #     geom_vline(xintercept = threshold)

                            #     d <- as.data.frame(table(test$black_v_white))

                            #     plot_d <- ggplot(data = d) + geom_col(aes(x = 1, y = Freq, fill = Var1))

                            #     cowplot::plot_grid(plot_a, plot_b, plot_c, plot_d)

                }

                shinyApp(ui = ui, server = server)

            }


        #### drawBaseline

            #' Draws a baseline on a signal
            #'
            #' @param data x, y dataframe
            #' @param prelim_baseline_window interval at which to draw baseline
            #' @examples
            #' @export
            #' drawBaseline

            drawBaseline <- function(data, prelim_baseline_window) {
      
                n_prelim_baseline_windows <- floor(length(data$x)/prelim_baseline_window)
                prelim_baseline <- list()
                for ( i in 1:n_prelim_baseline_windows ) {
                  abundances_in_window <- data$y[((prelim_baseline_window*(i-1))+1):(prelim_baseline_window*i)]
                  prelim_baseline[[i]] <- data.frame(
                      x = data$x[(which.min(abundances_in_window)+((i-1)*prelim_baseline_window))],
                      min = min(abundances_in_window)
                  )
                }
                prelim_baseline <- do.call(rbind, prelim_baseline)
                data$in_prelim_baseline <- FALSE
                data$in_prelim_baseline[data$x %in% prelim_baseline$x] <- TRUE

                y = prelim_baseline$min
                x = prelim_baseline$x

                baseline2 <- data.frame(
                  x = data$x,
                  y = approx(x, y, xout = data$x)$y
                )
                baseline2 <- baseline2[!is.na(baseline2$y),]
                data <- data[data$x %in% baseline2$x,]
                data$baseline <- baseline2$y

                return(data)

            }

