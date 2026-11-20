## Moved verbatim from phylochemistry.R (phylochemistry website tree).

            genomeScope <- function( input_histogram_file, output_dir, ploidy, kmer_length ) {

                ## General set up stuff

                    # library('argparse') ## Load libraries for non-linear least squares and argument parser
                    # library('genomescope') ## Load the genomescope library

                    ## Define default arguments

                        arguments <- vector()
                        arguments$input <- input_histogram_file # parser$add_argument("-i", "--input", help = "input histogram file")
                        arguments$output <- output_dir # parser$add_argument("-o", "--output", help = "output directory name")
                        arguments$ploidy <- ploidy # parser$add_argument("-p", "--ploidy", type = "integer", default = 2, help = "ploidy (1, 2, 3, 4, 5, or 6) for model to use [default 2]")
                        arguments$kmer_length <- kmer_length # parser$add_argument("-k", "--kmer_length", type = "integer", default = 21, help = "kmer length used to calculate kmer spectra [default 21]")
                        arguments$version <- FALSE
                        arguments$name_prefix <- "test" # parser$add_argument("-n", "--name_prefix", default = "", help = "optional name_prefix for output files")
                        arguments$fitted_hist <- FALSE # parser$add_argument("--fitted_hist", action="store_true", default=FALSE, help = "ADVANCED: generates a fitted histogram for kmer multiplicity 0-4 and a lookup table of probabilities")
                        arguments$initial_heterozygosities <- -1 # parser$add_argument("--initial_heterozygosities", type="character", default = -1, help = "ADVANCED: flag to set initial values for nucleotide heterozygosity rates")
                        arguments$initial_repetitiveness <- -1 # parser$add_argument("--initial_repetitiveness", type="character", default = -1, help = "ADVANCED: flag to set initial value for repetitiveness")
                        arguments$lambda <- -1 # parser$add_argument("-l", "--lambda", "--kcov", "--kmercov", type = "integer", default=-1, help = "optional initial kmercov estimate for model to use")
                        arguments$max_kmercov <- -1 # parser$add_argument("-m", "--max_kmercov", type = "integer", default=-1, help = "optional maximum kmer coverage threshold (kmers with coverage greater than max_kmercov are ignored by the model)")
                        arguments$name_prefix <- ""
                        arguments$no_unique_sequence <- FALSE # parser$add_argument("--no_unique_sequence", action="store_true", default=FALSE, help = "optional flag to turn off yellow unique sequence line in plots")
                        arguments$num_rounds <- 4 # parser$add_argument("--num_rounds", type = "integer", default = 4, help = "ADVANCED: parameter for the number of optimization rounds")
                        arguments$start_shift <- 5 # parser$add_argument("--start_shift", type = "integer", default=START_SHIFT, help = "ADVANCED: coverage shifts to exclude between fitting rounds")
                        arguments$testing <- FALSE # parser$add_argument("--testing", action="store_true", default=FALSE, help = "ADVANCED: flag to create testing.tsv file with model parameters")
                        arguments$topology <- 0 # parser$add_argument("-t", "--topology", type = "integer", default = 0, help = "ADVANCED: flag for topology for model to use")
                        arguments$trace_flag <- FALSE # parser$add_argument("--trace_flag", action="store_true", default=FALSE, help = "ADVANCED: flag to turn on printing of iteration progress of nlsLM function")
                        arguments$transform_exp <- 1 # parser$add_argument("--transform_exp", type="integer", default=1, help = "ADVANCED: parameter for the exponent when fitting a transformed (x**transform_exp*y vs. x) kmer histogram [default 1]")
                        arguments$true_params <- -1 # parser$add_argument("--true_params", type="character", default = -1, help = "ADVANCED: flag to state true simulated parameters for testing mode")
                        arguments$typical_error <- 15 # parser$add_argument("--typical_error", type = "integer", default=TYPICAL_ERROR, help = "ADVANCED: typical level of sequencing error")
                        arguments$verbose <- FALSE # parser$add_argument("--verbose", action="store_true", default=FALSE, help = "optional flag to print messages during execution")
                        arguments$version <- FALSE # parser$add_argument("-v", "--version", action="store_true", default=FALSE, help="print the version and exit")
                        version_message <- "GenomeScope 2.0\n"

                        if (is.null(arguments$input) | is.null(arguments$output)) { stop("Please provide input and output") }

                    ## Transfer arguments to variable names
                        
                        histfile <- arguments$input
                        foldername  <- arguments$output
                        p           <- arguments$ploidy
                        k           <- arguments$kmer_length
                        estKmercov  <- arguments$lambda
                        max_kmercov <- arguments$max_kmercov
                        VERBOSE     <- arguments$verbose
                        NO_UNIQUE_SEQUENCE <- arguments$no_unique_sequence
                        topology    <- arguments$topology
                        d_init      <- arguments$initial_repetitiveness
                        r_inits     <- arguments$initial_heterozygosities
                        transform_exp <- arguments$transform_exp
                        TESTING     <- arguments$testing
                        TRUE_PARAMS <- arguments$true_params
                        TRACE_FLAG <- arguments$trace_flag
                        NUM_ROUNDS <- arguments$num_rounds
                        FITTED_HIST <- arguments$fitted_hist
                        START_SHIFT <- arguments$start_shift
                  
                    ## Define some other variables
                        
                        NUM_ROUNDS = 4 ## Number of rounds before giving up
                        START_SHIFT = 5 ## Coverage steps to trim off between rounds
                        TYPICAL_ERROR = 15 ## Typical cutoff for sequencing error
                        MAX_ITERATIONS=200 ## Max rounds on NLS
                        SCORE_CLOSE = 0.20 ## Overrule if two scores are within this percent (0.05 = 5%) but larger difference in het
                        SCORE_HET_FOLD_DIFFERENCE = 10 ## Overrule heterozygosity if there is a large difference in het rate
                    
                    ## Suppress the warnings if the modeling goes crazy, those are in try/catch blocks anyways
                        options(warn = -1)

                    ## Colors for plots
                        COLOR_BGCOLOR  = "light grey"
                        COLOR_HIST     = "#56B4E9"
                        COLOR_2pPEAK   = "black"
                        COLOR_pPEAK    = "#F0E442"
                        COLOR_ERRORS   = "#D55E00"
                        COLOR_KMERPEAK = "black"
                        COLOR_RESIDUAL = "purple"
                        COLOR_COVTHRES = "red"

                    ## Given mean +/- stderr, report min and max value within 2 SE

                        min_max <- function(table) {
                          return (c(max(0,table[1] - 2*table[2]), table[1]+ 2*table[2]))
                        }

                        min_max1 <- function(table) {
                          return (c(max(0,table[1] - 2*table[2]), min(1, table[1]+ 2*table[2])))
                        }

                ## Program

                    cat(paste("GenomeScope analyzing ", histfile, " p=", p, " k=", k, " outdir=", foldername, "\n", sep=""))

                    dir.create(foldername, showWarnings=FALSE)

                    ## Initialize the status
                    progressFilename <- paste(foldername,"/", arguments$name_prefix, "progress.txt",sep="")
                    cat("starting", file=progressFilename, sep="\n")

                    kmer_prof <- read.csv(file=histfile,sep="", header=FALSE,colClasses=c("numeric","numeric"))

                    minkmerx = 1;
                    if (kmer_prof[1,1] == 0) {
                    if (VERBOSE) {cat("Histogram starts with zero, reseting minkmerx\n")}
                    minkmerx = 2;
                    }

                    kmer_prof <- kmer_prof[c(minkmerx:(length(kmer_prof[,2])-1)),] #get rid of the last position
                    kmer_prof_orig <- kmer_prof

                    ## try to find the local minimum between errors and the first (heterozygous) peak
                        kmer_trans = as.numeric(kmer_prof[,1])**transform_exp*as.numeric(kmer_prof[,2])
                        start <- tail(which(kmer_trans[1:TYPICAL_ERROR]==min(kmer_trans[1:TYPICAL_ERROR])),n=1)
                        start_max <- start + which(kmer_trans[start:length(kmer_trans)]==max(kmer_trans[start:length(kmer_trans)])) - 1

                    maxCovIndex = -1

                    ## Figure out which kmers to exclude, if any
                        if(max_kmercov == -1) {
                            maxCovIndex <- length(kmer_prof[,1])
                            max_kmercov <- kmer_prof[maxCovIndex,1]
                        } else {
                    
                    ## Figure out the index we should use for this coverage length
                            x <- kmer_prof[,1]
                            maxCovIndex <- length(x[x<=max_kmercov])
                        }

                    if (VERBOSE) {cat(paste("using max_kmercov:", max_kmercov, " with index:", maxCovIndex, "\n"))}

                    # terminate after NUM_ROUND iterations, store best result so far in container
                    round <- 0
                    best_container <- list(NULL,0)

                ## The loop

                    while (round < NUM_ROUNDS) {
                        
                        cat(paste("round", round, "trimming to", start, "trying 2p peak model... "), file=progressFilename, sep="", append=TRUE)
                        if (VERBOSE) {
                            cat(paste("round", round, "trimming to", start, "trying 2p peak model... \n"))
                        }

                        ## Reset the input trimming off low frequency error kmers
                            kmer_prof=kmer_prof_orig[1:maxCovIndex,]
                            x <- kmer_prof[start:maxCovIndex,1]
                            y <- kmer_prof[start:maxCovIndex,2]

                                ## Below is a non-subfunction implementation of estimate_Genome_peakp, to simplify code

                                    kmer_hist_orig = kmer_prof

                                    if (topology==-1) {
                                        p_to_num_topologies = c(1, 1, 1, 2, 5, 16)
                                        num_topologies = p_to_num_topologies[p]
                                        topologies = 1:num_topologies
                                    } else {
                                        num_topologies = 1
                                        topologies = c(topology)
                                    }
                                    
                                    numofKmers = sum(as.numeric(x)*as.numeric(y))
                                    
                                    if (estKmercov==-1) {
                                        #In situations with low heterozygosity, the peak with highest amplitude typically corresponds to the homozygous peak (i.e. the p-th peak).
                                        #However, with increasing heterozygosity, the highest amplitude peak may be an earlier peak.
                                        #Thus, when setting the estimated kmer coverage, we will need to iterate through these possibilities.
                                        #num_peak_indices indicates how many possibilities we need to iterate through.
                                        num_peak_indices = p
                                        y_transform = as.numeric(x)**transform_exp*as.numeric(y)
                                        estKmercov1 = x[which(y_transform==max(y_transform))][1]
                                    } else {
                                        # When the user sets the estimated kmer coverage, we only need to iterate through one possibility
                                        num_peak_indices = 1
                                        ## We set the estimated kmer coverage to be the user specified value
                                        estKmercov1 = estKmercov
                                    }
                                    
                                    estLength1 = numofKmers/estKmercov1

                                    nls00 = NULL
                                    peak_indices = 1:num_peak_indices
                                    for (i in peak_indices) {
                                        nls0 = NULL
                                        top_count = 0
                                        ## We see what happens when we set the estimated kmer coverage to be 1/i times the x-coordinate where the max peak occurs (1 <= i <= p if the user doesn't set the estimated kmer coverage, and i=1 if they do)
                                        estKmercov2 = estKmercov1 / i
                                        estLength2 = numofKmers/estKmercov2

                                        if (VERBOSE) {cat(paste("trying with kmercov: ", estKmercov2, "\n"))}

                                        for (top in topologies) {
                                            if (VERBOSE) {cat(paste("trying with topology: ", top, "\n"))}
                                            top_count = top_count + 1
                                            

                                            ## Here is a non-subfunction implementation of nls_peak, to simplify code
    
                                                # Initiate variables
                                                    
                                                    model = NULL
                                                    best_deviance = Inf
                                                    d_min = 0
                                                    if (d_init!=-1) {
                                                      d_initial = d_init
                                                    } else {
                                                      d_initial = 0.10
                                                    }
                                                    d_max = 1
                                                    r_min = 0.00001
                                                    if (top==0) {
                                                      p_to_num_r = c(0, 1, 2, 4, 6, 10)
                                                    } else {
                                                      p_to_num_r = c(0, 1, 2, 3, 4, 5)
                                                    }
                                                    num_r = p_to_num_r[p]
                                                    r_max = 1
                                                    kmercov_min = 0
                                                    kmercov_initial = estKmercov2
                                                    kmercov_max = Inf
                                                    bias_min = 0
                                                    bias_initial = 0.5
                                                    bias_max = Inf
                                                    length_min = 0
                                                    length_initial = estLength2/p
                                                    length_max = Inf

                                                # Determine what formula to use, based on p
                                                    
                                                    if (p==1) {
                                                      r_text = ""
                                                    } else {
                                                      r_text = paste(paste(lapply(1:(num_r), function(x) paste("r", as.character(x), sep="")), collapse=", "), ", ")
                                                    }
                                                    
                                                    x = x[1:min(2000,length(x))]
                                                    y = y[1:min(2000,length(y))]
                                                    y_transform = as.numeric(x)**transform_exp*as.numeric(y)
                                                    formula = as.formula(paste("y_transform ~ x**transform_exp*length*predict",p,"_",top,"(",r_text, "k, d, kmercov, bias, x)",sep=""))

                                                    if (VERBOSE) {cat("trying nlsLM algorithm (Levenberg-Marquardt)\n")}

                                                    if (r_inits!=-1) {
                                                        r_initials = unlist(lapply(strsplit(r_inits,","),as.numeric))
                                                        if (length(r_initials)!=num_r) {
                                                            stop("Incorrect number of initial rates supplied.")
                                                        }
                                                        r_initials_list = list(r_initials)
                                                    } else {
                                                        r_initials_list = list(rep(0.001, num_r), 0.001*(1:num_r), 0.001*(num_r:1), rep(0.01, num_r), 0.01*(1:num_r), 0.01*(num_r:1))
                                                    }

                                                    for (r_initials in r_initials_list) {

                                                        model1 = NULL
                                                        r_start = vector("list", num_r)
                                                        if (p > 1) {
                                                            names(r_start) = paste("r", 1:(num_r), sep="")
                                                            for (i in 1:(num_r)) {
                                                                r_start[[paste("r",i,sep="")]] = r_initials[i]
                                                            }
                                                        }

                                                        try(model1 <- nlsLM(formula = formula,
                                                            start   = c(list(d = d_initial), r_start, list(kmercov = kmercov_initial, bias = bias_initial, length = length_initial)),
                                                            lower   = c(c(d_min), rep(r_min, num_r), c(kmercov_min, bias_min, length_min)),
                                                            upper   = c(c(d_max), rep(r_max, num_r), c(kmercov_max, bias_max, length_max)),
                                                            control = list(minFactor=1e-12, maxiter=MAX_ITERATIONS, factor=0.1), trace=TRACE_FLAG), silent = TRUE
                                                        )

                                                        if (!is.null(model1)) {
                                                            current_deviance = model1$m$deviance()
                                                            #cat("Model deviance: ", current_deviance, "\n")
                                                            if (current_deviance < best_deviance) {
                                                              model = model1
                                                              best_deviance = current_deviance
                                                            }
                                                        } else {
                                                            #print("Model did not converge.")
                                                        }

                                                    }

                                                    if (!is.null(model)) {
                                                        model_sum    = summary(model)
                                                        model$p      = p
                                                        model$top = top
                                                        
                                                        if (p==1) {
                                                            model$hets = list(c(0, 0))
                                                        } else {
                                                            model$hets = lapply(1:(num_r), function(x) min_max1(model_sum$coefficients[paste('r', x, sep=""),]))
                                                        }

                                                        #model$het = c(1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 1))), 1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 2))))
                                                        model$het = c(sum(sapply(model$hets, '[[', 1)), sum(sapply(model$hets, '[[', 2)))
                                                        model$homo = 1-model$het
                                                        model$dups   = min_max(model_sum$coefficients['bias',])
                                                        model$kcov   = min_max(model_sum$coefficients['kmercov',])
                                                        model$mlen   = min_max(model_sum$coefficients['length',])
                                                        model$md     = min_max1(model_sum$coefficients['d',])
                                                        
                                                        if (p==1) {
                                                            model$ahets = list(c(0))
                                                        } else {
                                                            model$ahets = lapply(1:(num_r), function(x) model_sum$coefficients[paste('r', x, sep=""),][[1]])
                                                        }
                                                        
                                                        #model$ahet = 1-Reduce("*", 1-unlist(model$ahets))
                                                        model$ahet = Reduce("+", model$ahets)
                                                        model$ahomo = 1-model$ahet
                                                        model$adups = model_sum$coefficients['bias',][[1]]
                                                        model$akcov = model_sum$coefficients['kmercov',][[1]]
                                                        model$amlen = model_sum$coefficients['length',][[1]]
                                                        model$amd   = model_sum$coefficients['d',][[1]]
                                                    }

                                            nls1 = model ## here, model is the output of the non-subfunction version of nls_peak
                                            nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments, score_close = SCORE_CLOSE)[[1]]
                                        }
                                        if (i < num_peak_indices) { #if this is not the last evaluation
                                            nls00 = eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments, score_close = SCORE_CLOSE)[[1]]
                                        }
                                    }

                                    model_peaks <- eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments, score_close = SCORE_CLOSE)


                        if (!is.null(model_peaks[[1]])) {
                            cat(paste("converged. score: ", model_peaks[[2]]$all[[1]]), file=progressFilename, sep="\n", append=TRUE)

                            if (VERBOSE) {
                                # mdir = paste(foldername, "/round", round, sep="")
                                # dir.create(mdir, showWarnings = FALSE)
                                # report_results(kmer_prof,kmer_prof_orig, k, p, model_peaks, mdir, arguments, TRUE)
                            }
                        } else { cat(paste("unconverged"), file = progressFilename, sep="\n", append=TRUE) }

                        #check if this result is better than previous
                        if (!is.null(model_peaks[[1]])) {
                            if (is.null(best_container[[1]])) {
                                if (VERBOSE) {cat("no previous best, updating best\n")}
                                best_container = model_peaks
                            }
                            else {
                                best_container_score = best_container[[1]]$m$deviance()
                                model_peaks_score = model_peaks[[1]]$m$deviance()
                                pdiff = abs(model_peaks_score - best_container_score) / max(model_peaks_score, best_container_score)

                                if (pdiff < SCORE_CLOSE) {
                                  hetm = model_peaks[[1]]$ahet
                                  hetb = best_container[[1]]$ahet

                                  #if (hetb * SCORE_HET_FOLD_DIFFERENCE < hetm) {
                                  if (hetb + 0.01 < hetm) {
                                    if (VERBOSE) {cat("model has significantly higher heterozygosity but similar score, overruling\n")}
                                  }
                                  #else if (hetm * SCORE_HET_FOLD_DIFFERENCE < hetb) {
                                  else if (hetm + 0.01 < hetb) {
                                    if (VERBOSE) {cat("previous best has significantly higher heterozygosity and similar score, keeping\n")}
                                    best_container = model_peaks
                                  }
                                  else if (model_peaks_score < best_container_score) {
                                    if (VERBOSE) {cat("score is marginally better but het rate is not extremely different, updating\n")}
                                    best_container = model_peaks
                                  }
                                }
                                else if (model_peaks_score < best_container_score) {
                                      if (VERBOSE) {cat("score is significantly better, updating\n")}
                                      best_container = model_peaks
                                }
                            }
                        }

                        ## Ignore a larger number of kmers as errors
                            start <- start + START_SHIFT
                            round <- round + 1
                    }

                ## Report the results, note using the original full profile
                    ## This is a non-subfunction implementation of report_results
                    
                    ## arguments to report_results
                        kmer_hist = kmer_prof
                        container = best_container
                        IN_VERBOSE = FALSE
  
                    output <- list()

                    ## Some initial stuff

                        kmer_hist = kmer_prof
                        kmer_hist_orig = kmer_prof_orig
                        k = k
                        p = p
                        container = best_container
                        foldername = foldername
                        arguments = arguments
                        IN_VERBOSE = FALSE

                        x=kmer_hist_orig[[1]]
                        y_orig=kmer_hist_orig[[2]]
                        y = as.numeric(x)**transform_exp*as.numeric(y_orig)
                        kmer_hist_transform = kmer_hist_orig
                        kmer_hist_transform$V2 = as.numeric(kmer_hist_transform$V1)**transform_exp * as.numeric(kmer_hist_transform$V2)
                        model = container[[1]]

                    # automatically zoom into the relevant regions of the plot, ignore first 15 positions
                    ## I think this is just for plotting and can be turned off
                        
                        xmax = length(x)
                        start_orig = which(y_orig == min(y_orig[1:TYPICAL_ERROR]))
                        start = which(y == min(y[1:TYPICAL_ERROR]))
                        zoomx = x[start:(xmax-1)]
                        zoomy_orig = y_orig[start_orig:(xmax-1)]
                        zoomy = y[start:(xmax-1)]
            
                    ## allow for a little space above max value past the noise
                    ## I think this is just for plotting and can be turned off
                        
                        y_limit_orig = max(zoomy_orig[start_orig:length(zoomy_orig)])*1.1
                        y_limit = max(zoomy[start:length(zoomy)])*1.1

                        x_limit_orig = which(y_orig == max(y_orig[start_orig:length(zoomx)])) * 3
                        x_limit = which(y == max(y[start:length(zoomx)])) * 3

                        if (min(zoomy_orig) > zoomy_orig[1]){
                          x_limit_orig=max(which(zoomy_orig<zoomy_orig[1])[2],600)
                        }
                        if (min(zoomy) > zoomy[1]){
                          x_limit=max(which(zoomy<zoomy[1])[2],600)
                        }

                        if (!is.null(model))
                        {
                          model_sum=summary(model)
                          kcov = min_max(model_sum$coefficients['kmercov',])[1]
                          x_limit_orig = max(kcov*(2*p+1.1), x_limit_orig)
                          x_limit = max(kcov*(2*p+1.1), x_limit)
                          if (model$top==0) {
                            p_to_num_r = c(0, 1, 2, 4, 6, 10)
                          } else {
                            p_to_num_r = c(0, 1, 2, 3, 4, 5)
                          }
                        } else {
                          if (topology==0) {
                            p_to_num_r = c(0, 1, 2, 4, 6, 10)
                          } else {
                            p_to_num_r = c(0, 1, 2, 3, 4, 5)
                          }
                        }
            
                    ## Uncomment this to enforce a specific number
                    ## I think this is just for plotting and can be turned off
                    # x_limit=150
            
                    ## Define some initial values

                        het=c(-1,-1)
                        homo=c(-1,-1)
                        num_r = p_to_num_r[p]
                        if (p > 1) {
                          hets = lapply(1:(num_r), function(x) c(-1, -1))
                          ahets = lapply(1:(num_r), function(x) -1)
                        }
                        amd = -1
                        akcov = -1
                        adups = -1
                        amlen = -1
                        atotal_len = -1
                        top = -1
                        total_len=c(-1,-1)
                        repeat_len=c(-1,-1)
                        unique_len=c(-1,-1)
                        dups=c(-1,-1)
                        error_rate=c(-1,-1)
                        model_status="fail"

                        model_fit_unique      = c(0,0,0)
                        model_fit_full        = c(0,0,0)
                        model_fit_all         = c(0,0,0)
                        model_fit_allscore    = c(0,0,0)
                        model_fit_fullscore   = c(0,0,0)
                        model_fit_uniquescore = c(0,0,0)

                        plot_size=2000
                        font_size=1.2
                        resolution=300

                    ## Plot the distribution, and hopefully with the model fit (currently turned off)
                        
                        # ylabel_orig = "Frequency"
                        # if (transform_exp == 1) {
                        # ylabel_transform = "Coverage*Frequency"
                        # } else {
                        # ylabel_transform = paste("Coverage^", transform_exp, "*Frequency", sep="")
                        # }
                        # png(paste(foldername, "/", arguments$name_prefix, "linear_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_orig, ylim=c(0,y_limit_orig), xlim=c(0,x_limit_orig),
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
                        # rect(0, 0, x_limit_orig*1.1 , y_limit_orig*1.1, col=COLOR_BGCOLOR)
                        # points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
                        # #  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
                        # #    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # #  }
                        # box(col="black")

                        # png(paste(foldername, "/", arguments$name_prefix, "transformed_linear_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_transform, ylim=c(0,y_limit), xlim=c(0,x_limit),
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
                        # rect(0, 0, x_limit*1.1 , y_limit*1.1, col=COLOR_BGCOLOR)
                        # points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
                        # #  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
                        # #    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # #  }
                        # box(col="black")

                        # ## Make a second plot in log space over entire range
                        # png(paste(foldername, "/", arguments$name_prefix, "log_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_orig, log="xy",
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # rect(1e-10, 1e-10, max(kmer_hist_orig[,1])*10 , max(kmer_hist_orig[,2])*10, col=COLOR_BGCOLOR)
                        # points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
                        # if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
                        # abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # }
                        # box(col="black")

                        # png(paste(foldername, "/", arguments$name_prefix, "transformed_log_plot.png", sep=""),
                        # width=plot_size, height=plot_size, res=resolution)
                        # par(mar = c(5.1,4.1,6.1,2.1))
                        # plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
                        # xlab="Coverage", ylab=ylabel_transform, log="xy",
                        # cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
                        # rect(1e-10, 1e-10, max(kmer_hist_transform[,1])*10 , max(kmer_hist_transform[,2])*10, col=COLOR_BGCOLOR)
                        # points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
                        # if(length(kmer_hist[,1])!=length(kmer_hist_transform[,1])){
                        # abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
                        # }
                        # box(col="black")

                    ## Now the modeling
                            
                            if(!is.null(model)) {
                                
                                x=kmer_hist[[1]]
                                y=kmer_hist[[2]]
                                y_transform = as.numeric(x)**transform_exp*as.numeric(y)

                                ## The model converged!
                                    pred=predict(model, newdata=data.frame(x))

                                ## Compute the genome characteristics
                                    model_sum=summary(model)
                                    #print(model_sum)

                                ## save the model to a file
                                    capture.output(model_sum, file=paste(foldername,"/", arguments$name_prefix, "model.txt", sep=""))

                                ## Identify key values
                                    top   = model$top
                                    hets  = model$hets
                                    het   = model$het
                                    ahets = model$ahets
                                    ahet  = model$ahet
                                    homo  = model$homo
                                    ahomo = model$ahomo

                                    dups = model$dups
                                    kcov = model$kcov
                                    mlen = model$mlen
                                    md   = model$md

                                    adups = model$adups
                                    akcov = model$akcov
                                    amlen = model$amlen
                                    amd   = model$amd

                                ## Compute error rate, by counting kmers unexplained by model through first peak
                                ## truncate errors as soon as it goes to zero, dont allow it to go back up
                                    error_xcutoff = max(1, floor(kcov[1]))
                                    error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
                                    if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

                                    error_kmers = x[1:error_xcutoff_ind]**(-transform_exp)*(y_transform[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind])

                                    first_zero = -1

                                    for (i in 1:error_xcutoff_ind) {
                                      if (first_zero == -1)
                                      {
                                        if (error_kmers[i] < 1.0)
                                        {
                                          first_zero = i
                                          if (VERBOSE) {cat(paste("Truncating errors at", i, "\n"))}
                                        }
                                      }
                                      else
                                      {
                                        error_kmers[i] = 0
                                      }
                                    }

                                    if (first_zero == -1) {
                                      first_zero = error_xcutoff_ind
                                      if (VERBOSE) {cat(paste("Truncating errors at", error_xcutoff_ind, "\n"))}
                                    }

                                ## Rather than "0", set to be some very small number so log-log plot looks okay
                                    error_kmers = pmax(error_kmers, 1e-10)

                                    total_error_kmers = sum(as.numeric(error_kmers) * as.numeric(x[1:error_xcutoff_ind]))

                                    total_kmers = sum(as.numeric(x)*as.numeric(y))

                                    error_rate = 1-(1-(total_error_kmers/total_kmers))**(1/k)
                                    error_rate = c(error_rate, error_rate)

                                    total_len = (total_kmers-total_error_kmers)/(p*kcov)
                                    atotal_len = (total_kmers-total_error_kmers)/(p*akcov)

                                ## find kmers that fit the p peak model (no repeats)
                                if (p==1) {
                                  unique_hist = amlen*predict1_1_unique(k, amd, akcov, adups, x)
                                }
                                if (p==2) {
                                  unique_hist = amlen*predict2_1_unique(ahets[[1]], k, amd, akcov, adups, x)
                                }
                                if (p==3) {
                                  unique_hist = amlen*predict3_1_unique(ahets[[1]], ahets[[2]], k, amd, akcov, adups, x)
                                }
                                if (p==4) {
                                  if (top==0) {
                                    unique_hist = amlen*predict4_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)
                                  } else {
                                    unique_hist = eval(parse(text = paste("amlen*predict4_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], k, amd, akcov, adups, x)", sep="")))
                                  }
                                }
                                if (p==5) {
                                  if (top==0) {
                                    unique_hist = amlen*predict5_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], k, amd, akcov, adups, x)
                                  } else {
                                    unique_hist = eval(parse(text = paste("amlen*predict5_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)", sep="")))
                                  }
                                }
                                if (p==6) {
                                  if (top==0) {
                                    unique_hist = amlen*predict6_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], k, amd, akcov, adups, x)
                                  } else {
                                    unique_hist = eval(parse(text = paste("amlen*predict6_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], k, amd, akcov, adups, x)", sep="")))
                                  }
                                }

                                r0 = 1-ahet #aa
                                t0 = r0**k #AA
                                s0 = t0 #AA
                                s1 = 1-t0 #AB
                                alpha_1 = (1-amd)*(2*s1) + amd*(2*s0*s1 + 2*s1**2)
                                alpha_2 = (1-amd)*(s0) + amd*(s1**2)
                                alpha_3 = amd*(2*s0*s1)
                                alpha_4 = amd*(s0**2)

                                one_hist = alpha_1 * dnbinom(x, size = akcov*1 / adups, mu = akcov*1)
                                two_hist = alpha_2 * dnbinom(x, size = akcov*p / adups, mu = akcov*p)
                                thr_hist = alpha_3 * dnbinom(x, size = akcov*3 / adups, mu = akcov*3)
                                fou_hist = alpha_4 * dnbinom(x, size = akcov*2*p / adups, mu = akcov*2*p)

                                unique_hist_transform = x**transform_exp*unique_hist

                                unique_kmers = sum(as.numeric(x)*as.numeric(unique_hist))
                                repeat_kmers = max(0, total_kmers - unique_kmers - total_error_kmers)

                                repeat_len=repeat_kmers/(p*kcov)
                                if (repeat_kmers == 0) {
                                  unique_len = total_len
                                } else {
                                  unique_len=unique_kmers/(p*kcov)
                                }

                                score = container[[2]]

                                model_fit_allscore    = score$allscore
                                model_fit_fullscore   = score$fullscore
                                model_fit_uniquescore = score$uniquescore

                                model_fit_all    = score$all
                                model_fit_full   = score$full
                                model_fit_unique = score$unique

                                residual_transform = y_transform - pred
                                residual = x**(-transform_exp)*residual_transform

                                hetline_simple = paste0("heterozygosity: ", format(100*ahet, digits=3), "%")

                                if (p==1) {
                                  hetline = paste0("a:", format(100*ahomo, digits=3), "%")
                                }
                                if (p==2) {
                                  hetline = paste0("aa:", format(100*ahomo,      digits=3), "% ",
                                                   "ab:", format(100*ahets[[1]], digits=3), "%")
                                }
                                if (p==3) {
                                  hetline = paste0("aaa:", format(100*ahomo,      digits=3), "% ",
                                                   "aab:", format(100*ahets[[1]], digits=3), "% ",
                                                   "abc:", format(100*ahets[[2]], digits=3), "%")
                                }
                                if (p==4) {
                                  if (top==0) {
                                    hetline = paste0("aaaa:", format(100*ahomo,      digits=3), "% ",
                                                     "aaab:", format(100*ahets[[1]], digits=3), "% ",
                                                     "aabb:", format(100*ahets[[2]], digits=3), "% ",
                                                     "aabc:", format(100*ahets[[3]], digits=3), "% ",
                                                     "abcd:", format(100*ahets[[4]], digits=3), "%")
                                  } else {
                                    hetline = paste0("aaaa:",                       format(100*ahomo,      digits=3), "% ",
                                                     switch(top, "aaab:", "aabb:"), format(100*ahets[[1]], digits=3), "% ",
                                                     "aabc:",                       format(100*ahets[[2]], digits=3), "% ",
                                                     "abcd:",                       format(100*ahets[[3]], digits=3), "%")
                                  }
                                }
                                if (p==5) {
                                  if (top==0) {
                                    hetline = paste0("aaaaa:", format(100*ahomo,      digits=3), "% ",
                                                     "aaaab:", format(100*ahets[[1]], digits=3), "% ",
                                                     "aaabb:", format(100*ahets[[2]], digits=3), "% ",
                                                     "aaabc:", format(100*ahets[[3]], digits=3), "% ",'\n',
                                                     "aabbc:", format(100*ahets[[4]], digits=3), "% ",
                                                     "aabcd:", format(100*ahets[[5]], digits=3), "% ",
                                                     "abcde:", format(100*ahets[[6]], digits=3), "%")
                                  } else {
                                    hetline = paste0("aaaaa:",                                                      format(100*ahomo,      digits=3), "% ",
                                                     switch(top, "aaaab:", "aaaab:", "aaabb:", "aaabb:", "aaabb:"), format(100*ahets[[1]], digits=3), "% ",
                                                     switch(top, "aaabc:", "aabbc:", "aaabc:", "aabcc:", "aabcc:"), format(100*ahets[[2]], digits=3), "% ",'\n',
                                                     switch(top, "aabcd:", "aabcd:", "aabcd:", "aabcd:", "abcdd:"), format(100*ahets[[3]], digits=3), "% ",
                                                     "abcde:",                                                      format(100*ahets[[4]], digits=3), "%")
                                  }
                                }
                                if (p==6) {
                                  if (top==0) {
                                    hetline = paste0("aaaaaa:", format(100*ahomo, digits=3), "% ",
                                                     "aaaaab:", format(100*ahets[[1]], digits=3), "% ",
                                                     "aaaabb:", format(100*ahets[[2]], digits=3), "% ",
                                                     "aaabbb:", format(100*ahets[[3]], digits=3), "% ",'\n',
                                                     "aaaabc:", format(100*ahets[[4]], digits=3), "% ",
                                                     "aaabbc:", format(100*ahets[[5]], digits=3), "% ",
                                                     "aabbcc:", format(100*ahets[[6]], digits=3), "% ",
                                                     "aaabcd:", format(100*ahets[[7]], digits=3), "% ",'\n',
                                                     "aabbcd:", format(100*ahets[[8]], digits=3), "% ",
                                                     "aabcde:", format(100*ahets[[9]], digits=3), "% ",
                                                     "abcdef:", format(100*ahets[[10]], digits=3), "%")
                                  } else {
                                    hetline = paste0("aaaaaa:", format(100*ahomo, digits=3), "% ",
                                                     switch(top, "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaabbb:", "aaabbb:", "aaabbb:"), format(100*ahets[[1]], digits=3), "% ",
                                                     switch(top, "aaaabc:", "aaaabc:", "aaabbc:", "aaabbc:", "aaabbc:", "aaaabc:", "aaaabc:", "aaabcc:", "aaabcc:", "aaabcc:", "aabbcc:", "aabbcc:", "aabbcc:", "aaabbc:", "aaabbc:", "aaabbc:"), format(100*ahets[[2]], digits=3), "% ",'\n',
                                                     switch(top, "aaabcd:", "aabbcd:", "aaabcd:", "aabccd:", "aabccd:", "aaabcd:", "aabbcd:", "aaabcd:", "aabcdd:", "aabcdd:", "aabbcd:", "aabcdd:", "aabcdd:", "aaabcd:", "aabccd:", "aabccd:"), format(100*ahets[[3]], digits=3), "% ",
                                                     switch(top, "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdde:", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdde:"), format(100*ahets[[4]], digits=3), "% ",
                                                     "abcdef:", format(100*ahets[[5]], digits=3), "%")
                                  }
                                }

                                if (p >= 5) {
                                  hetline = hetline_simple
                                }

                                if (!IN_VERBOSE) {
                                  cat(paste0(hetline,"\n"))
                                }

                                # dev.set(dev.next())

                                ### THIS IS ALL NOW TURNED OFF

                            # ## finish plots (off)

                            #   ## Finish Linear Plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, x**(-transform_exp)*pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if (NO_UNIQUE_SEQUENCE) {
                            #       legend(.62 * x_limit_orig, 1.0 * y_limit_orig,
                            #       legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     } else {
                            #       legend(.62 * x_limit_orig, 1.0 * y_limit_orig,
                            #       legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     }

                            #     dev.set(dev.next())

                            #   ## Finish Linear Plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist_transform, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], (x[1:error_xcutoff_ind]**transform_exp)*error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual_transform, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if (NO_UNIQUE_SEQUENCE) {
                            #       legend(.62 * x_limit, 1.0 * y_limit,
                            #       legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     } else {
                            #       legend(.62 * x_limit, 1.0 * y_limit,
                            #       legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #       lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #       lwd=c(3,3,3,3,2),
                            #       col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #       bg="white")
                            #     }

                            #     dev.set(dev.next())

                            #   ## Finish Log plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, x**(-transform_exp)*pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if(length(kmer_hist[,1])==length(kmer_hist_orig[,1]))
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       } else {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       }
                            #     }
                            #     else
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       } else {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       }
                            #     }

                            #     dev.set(dev.next())

                            #   ## Finish Log plot
                            #     title(paste("\n\nlen:",  prettyNum(total_len[1], big.mark=","),
                            #     "bp",
                            #     " uniq:", format(100*(unique_len[1]/total_len[1]), digits=3),
                            #     "% ", "\n",
                            #     hetline, "\n",
                            #     " kcov:", format(akcov, digits=3),
                            #     " err:",   format(100*error_rate[1], digits=3),
                            #     "% ",
                            #     " dup:",  format(adups, digits=3),
                            #     " ",
                            #     " k:",   format(k, digits=3),
                            #     " p:",   format(p, digits=3),
                            #     sep=""),
                            #     cex.main=.85)

                            #     ## Mark the modes of the peaks
                            #     abline(v=akcov * (1:(2*p)), col=COLOR_KMERPEAK, lty=2)

                            #     ## Draw just the unique portion of the model
                            #     if (!NO_UNIQUE_SEQUENCE) {
                            #       lines(x, unique_hist_transform, col=COLOR_pPEAK, lty=1, lwd=3)
                            #     }
                            #     lines(x, pred, col=COLOR_2pPEAK, lwd=3)
                            #     lines(x[1:error_xcutoff_ind], (x[1:error_xcutoff_ind]**transform_exp)*error_kmers, lwd=3, col=COLOR_ERRORS)

                            #     if (VERBOSE) {
                            #       lines(x, residual_transform, col=COLOR_RESIDUAL, lwd=3)
                            #     }

                            #     ## Add legend
                            #     if(length(kmer_hist[,1])==length(kmer_hist_orig[,1]))
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       } else {
                            #         legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed"),
                            #         lwd=c(3,3,3,3,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK),
                            #         bg="white")
                            #       }
                            #     }
                            #     else
                            #     {
                            #       if (NO_UNIQUE_SEQUENCE) {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       } else {
                            #         legend("topright",
                            #         ##legend(exp(.62 * log(max(x))), 1.0 * max(y),
                            #         legend=c("observed", "full model", "unique sequence", "errors", "kmer-peaks","cov-threshold"),
                            #         lty=c("solid", "solid", "solid", "solid", "dashed", "dashed"),
                            #         lwd=c(3,3,3,3,2,3),
                            #         col=c(COLOR_HIST, COLOR_2pPEAK, COLOR_pPEAK, COLOR_ERRORS, COLOR_KMERPEAK, COLOR_COVTHRES),
                            #         bg="white")
                            #       }
                            #     }

                                ## other things (also turned off)

                            # model_status="done"

                            # if (!IN_VERBOSE) {
                            #   cat(paste("Model converged het:", format(ahet, digits=3),
                            #   " kcov:", format(akcov, digits=3),
                            #   " err:", format(error_rate[1], digits=3),
                            #   " model fit:", format(adups, digits=3),
                            #   " len:", round(total_len[1]), "\n", sep=""))
                            # }
                            # }
                            # else
                            # {
                            #   title("\nFailed to converge")
                            #   dev.set(dev.next())
                            #   title("\nFailed to converge")
                            #   cat("Failed to converge.", file=paste(foldername,"/", arguments$name_prefix, "model.txt", sep=""))
                            #   cat("Failed to converge.\n")
                            # }

                            # dev.off()
                            # dev.off()
                            # dev.off()
                            # dev.off()

                                ## Write key values to summary file

                                    summaryFile <- paste(foldername,"/", arguments$name_prefix, "summary.txt",sep="")

                                    format_column_1 = "%-30s"
                                    format_column_2 = "%-18s"
                                    format_column_3 = "%-18s"

                                    cat(paste("GenomeScope version 2.0", sep=""), file=summaryFile, sep="\n")
                                    cat(paste("input file = ", arguments$input, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste("output directory = ", arguments$output, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste("p = ", p,sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste("k = ", k,sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    if (arguments$name_prefix!="") {
                                      cat(paste("name prefix = ", substring(arguments$name_prefix,1,nchar(arguments$name_prefix)-1), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (arguments$lambda!=-1) {
                                      cat(paste("initial kmercov estimate = ", arguments$lambda, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (arguments$max_kmercov!=-1) {
                                      cat(paste("max_kmercov = ", arguments$max_kmercov, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (VERBOSE) {
                                      cat(paste("VERBOSE set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (NO_UNIQUE_SEQUENCE) {
                                      cat(paste("NO_UNIQUE_SEQUENCE set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (topology!=0) {
                                      cat(paste("topology = ", topology, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (d_init!=-1) {
                                      cat(paste("initial repetitiveness = ", d_init, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (r_inits!=-1) {
                                      cat(paste("initial heterozygosities = ", r_inits, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (transform_exp != 1) {
                                      cat(paste("TRANSFORM_EXP = ", transform_exp, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (TESTING) {
                                      cat(paste("TESTING set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (TRUE_PARAMS != -1) {
                                      cat(paste("TRUE_PARAMS = ", TRUE_PARAMS, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (TRACE_FLAG) {
                                      cat(paste("TRACE_FLAG set to TRUE", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (NUM_ROUNDS != 4) {
                                      cat(paste("NUM_ROUNDS = ", NUM_ROUNDS, sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    }

                                    property <- vector()
                                    min <- vector()
                                    max <- vector()

                                    cat(paste("\n",sprintf(format_column_1,"property"),               sprintf(format_column_2,"min"),                              sprintf(format_column_3,"max"), sep=""),                                     file=summaryFile, sep="\n", append=TRUE)
                                    if (p==1) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (a)"),            sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)

                                    }
                                    if (p==2) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aa)"),           sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (ab)"),         sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      property <- c(property, "Homozygous (aa)", "Heterozygous (ab)")
                                      min <- c(min, homo[2]*100, hets[[1]][1]*100)
                                      max <- c(max, homo[1]*100, hets[[1]][2]*100)
                                    }
                                    if (p==3) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaa)"),          sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaa)"),    sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])),  sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"aab"),                       sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"abc"),                       sprintf(format_column_2,percentage_format(hets[[2]][1])),         sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                    }
                                    if (p==4) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaaa)"),                    sprintf(format_column_2,percentage_format(homo[2])),      sprintf(format_column_3,percentage_format(homo[1])),      sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaaa)"),              sprintf(format_column_2,percentage_format(het[1])),       sprintf(format_column_3,percentage_format(het[2])),       sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1, switch(top+1, "aaab", "aaab", "aabb")), sprintf(format_column_2,percentage_format(hets[[1]][1])), sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1, switch(top+1, "aabb", "aabc", "aabc")), sprintf(format_column_2,percentage_format(hets[[2]][1])), sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1, switch(top+1, "aabc", "abcd", "abcd")), sprintf(format_column_2,percentage_format(hets[[3]][1])), sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      if (top == 0) {
                                        cat(paste(sprintf(format_column_1,"abcd"),                               sprintf(format_column_2,percentage_format(hets[[4]][1])), sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      }
                                    }
                                    if (p==5) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaaaa)"),                                                sprintf(format_column_2,percentage_format(ahomo)), sep=""),      file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaaaa)"),                                          sprintf(format_column_2,percentage_format(ahet)),  sep=""),      file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aaaab", "aaaab", "aaaab", "aaabb", "aaabb", "aaabb")), sprintf(format_column_2,percentage_format(hets[[1]][1])), sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aaabb", "aaabc", "aabbc", "aaabc", "aabcc", "aabcc")), sprintf(format_column_2,percentage_format(hets[[2]][1])), sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aaabc", "aabcd", "aabcd", "aabcd", "aabcd", "abcdd")), sprintf(format_column_2,percentage_format(hets[[3]][1])), sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1,"aabbc", "abcde", "abcde", "abcde", "abcde", "abcde")), sprintf(format_column_2,percentage_format(hets[[4]][1])), sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                      if (top == 0) {
                                        #cat(paste(sprintf(format_column_1,"aabcd"),                     sprintf(format_column_2,percentage_format(hets[[5]][1])),         sprintf(format_column_3,percentage_format(hets[[5]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"abcde"),                     sprintf(format_column_2,percentage_format(hets[[6]][1])),         sprintf(format_column_3,percentage_format(hets[[6]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      }
                                    }
                                    if (p==6) {
                                      cat(paste(sprintf(format_column_1,"Homozygous (aaaaaa)"),       sprintf(format_column_2,percentage_format(ahomo)), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      cat(paste(sprintf(format_column_1,"Heterozygous (not aaaaaa)"), sprintf(format_column_2,percentage_format(ahet)), sep=""),                 file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaaaab", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaaab:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaaabb:", "aaabbb:", "aaabbb:", "aaabbb:")),                    sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaaabb", "aaaabc:", "aaaabc:", "aaabbc:", "aaabbc:", "aaabbc:", "aaaabc:", "aaaabc:", "aaabcc:", "aaabcc:", "aaabcc:", "aabbcc:", "aabbcc:", "aabbcc:", "aaabbc:", "aaabbc:", "aaabbc:")),                    sprintf(format_column_2,percentage_format(hets[[2]][1])),         sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaabbb", "aaabcd:", "aabbcd:", "aaabcd:", "aabccd:", "aabccd:", "aaabcd:", "aabbcd:", "aaabcd:", "aabcdd:", "aabcdd:", "aabbcd:", "aabcdd:", "aabcdd:", "aaabcd:", "aabccd:", "aabccd:")),                    sprintf(format_column_2,percentage_format(hets[[3]][1])),         sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaaabc", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdde:", "aabcde:", "aabcde:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdee:", "aabcde:", "aabcde:", "abcdde:")),                    sprintf(format_column_2,percentage_format(hets[[4]][1])),         sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      #cat(paste(sprintf(format_column_1, switch(top+1, "aaabbc", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:", "abcdef:")),                    sprintf(format_column_2,percentage_format(hets[[5]][1])),         sprintf(format_column_3,percentage_format(hets[[5]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                      if (top==0) {
                                        #cat(paste(sprintf(format_column_1,"aabbcc"),                    sprintf(format_column_2,percentage_format(hets[[6]][1])),         sprintf(format_column_3,percentage_format(hets[[6]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"aaabcd"),                    sprintf(format_column_2,percentage_format(hets[[7]][1])),         sprintf(format_column_3,percentage_format(hets[[7]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"aabbcd"),                    sprintf(format_column_2,percentage_format(hets[[8]][1])),         sprintf(format_column_3,percentage_format(hets[[8]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"aabcde"),                    sprintf(format_column_2,percentage_format(hets[[9]][1])),         sprintf(format_column_3,percentage_format(hets[[9]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
                                        #cat(paste(sprintf(format_column_1,"abcdef"),                    sprintf(format_column_2,percentage_format(hets[[10]][1])),        sprintf(format_column_3,percentage_format(hets[[10]][2])), sep=""),               file=summaryFile, sep="\n", append=TRUE)
                                      }
                                    }
                                    cat(paste(sprintf(format_column_1,"Genome Haploid Length"), sprintf(format_column_2,bp_format(total_len[2])),                  sprintf(format_column_3,bp_format(total_len[1])), sep=""),                   file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Genome Repeat Length"),  sprintf(format_column_2,bp_format(repeat_len[2])),                 sprintf(format_column_3,bp_format(repeat_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Genome Unique Length"),  sprintf(format_column_2,bp_format(unique_len[2])),                 sprintf(format_column_3,bp_format(unique_len[1])), sep=""),                  file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Model Fit "),            sprintf(format_column_2,percentage_format(model_fit_allscore[1])), sprintf(format_column_3,percentage_format(model_fit_fullscore[1])), sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    cat(paste(sprintf(format_column_1,"Read Error Rate"),       sprintf(format_column_2,percentage_format(error_rate[1])),         sprintf(format_column_3,percentage_format(error_rate[2])), sep=""),          file=summaryFile, sep="\n", append=TRUE)

                                    property <- c(property, "Genome Haploid Length", "Genome Repeat Length", "Genome Unique Length", "Model Fit", "Read Error Rate")
                                    min <- c(min, total_len[2], repeat_len[2], unique_len[2], model_fit_allscore[1], error_rate[1])
                                    max <- c(max, total_len[1], repeat_len[1], unique_len[1], model_fit_fullscore[1], error_rate[2])

                                    output$summary <- data.frame(
                                      property = property,
                                      min = min,
                                      max = max
                                    )

                                    # if (VERBOSE)
                                    # {
                                    #   cat(paste("\nPercent Kmers Modeled (All Kmers) = ",  percentage_format(model_fit_allscore[1]),    " [", model_fit_allscore[2],    ", ", model_fit_allscore[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Percent Kmers Modeled (Full Model) = ",   percentage_format(model_fit_fullscore[1]),   " [", model_fit_fullscore[2],   ", ", model_fit_fullscore[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Percent Kmers Modeled (Unique Kmers) = ", percentage_format(model_fit_uniquescore[1]), " [", model_fit_uniquescore[2], ", ", model_fit_uniquescore[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)

                                    #   cat(paste("\nModel RSSE (All Kmers) = ",  model_fit_all[1],    " [", model_fit_all[2],    ", ", model_fit_all[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Model RSSE (Full Model) = ",   model_fit_full[1],   " [", model_fit_full[2],   ", ", model_fit_full[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
                                    #   cat(paste("Model RSSE (Unique Model) = ", model_fit_unique[1], " [", model_fit_unique[2], ", ", model_fit_unique[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE) 
                                    # }

                                    ## Merfin probabilities (off)

                                # colors <- c("black","red","green","purple","blue")

                                # if (p<=2 & FITTED_HIST){

                                #   if (!(p==1)){

                                #   fitted_hist=data.frame(cbind(one_hist,two_hist,thr_hist,fou_hist))

                                #   } else {

                                #   fitted_hist=data.frame(cbind(two_hist,fou_hist))

                                #   }

                                #   ## Fitted histogram  

                                #     png(paste(foldername, "/fitted_hist.png", sep=""), height = plot_size, width = plot_size, res=resolution)
                                #     layout(matrix(c(1,2), nrow=2, byrow = TRUE),heights=lcm(c(11,5.5)))
                                #     par(mar=c(0,5,1,1))
                                  
                                #     plot(kmer_hist_orig, type="n", xlab=NULL, ylab=ylabel_orig, ylim=c(0,y_limit_orig), xlim=c(0,akcov*(p*2.5)),
                                #       cex.main=font_size, cex.sub=font_size, cex.axis=font_size, cex.lab=font_size,
                                #        xaxt='n', yaxt='n')
                                #     myTicks = axTicks(4)
                                #     myTicks[1]<-0
                                #     axis(2, at=myTicks, labels=myTicks, cex.axis=font_size,
                                #        cex.lab=font_size)
                                  
                                #     peaks<-akcov * 1:6

                                #     points(kmer_hist, type="h", col=COLOR_HIST, lwd=2)
                                #     lines(error_kmers, col="black", lwd=1.5)
                                  
                                #     abline(v=peaks, col="black", lty=2, lwd=0.3)
                                  
                                #     for (i in 1:ncol(fitted_hist)) {
                                #     lines(fitted_hist[,i]*amlen, type="l", lwd=1.5, col=colors[i+1])
                                #     }

                                #     lines(rowSums(fitted_hist*amlen), col="darkgray", lwd=2, lty=2)
                                  
                                #     legend_names=c("0-copy", "1-copy", "2-copy", "3-copy", "4-copy")
                                #     legend_lty=c("solid", "solid", "solid", "solid", "solid")
                                #     legend_lwd=c(3,3,3,3,3)
                                  
                                #     legend("topright",
                                #        ##legend(exp(.65 * log(max(x))), 1.0 * max(y),
                                #        legend=c("Observed",legend_names[1:(ncol(fitted_hist)+1)],"Full model"),
                                #        lty=c("solid",legend_lty[1:(ncol(fitted_hist)+1)], "dashed"),
                                #        lwd=c(2, legend_lwd[1:(ncol(fitted_hist)+1)], 3),
                                #        col=c(COLOR_HIST,colors[1:(ncol(fitted_hist)+1)], "darkgray"),
                                #        bg="white")

                                #     ## Generate lookup_table
                                  
                                #       lookup_table <- NULL
                                #       plot_table <- NULL
                                    
                                #       fitted_hist[(akcov*(p*2.5)):length(one_hist),1:ncol(fitted_hist)] <- 0
                                    
                                #       fitted_hist <- na.zero(fitted_hist)
                                    
                                #       for (i in 1:(akcov*(p*2.5)-1)){
                                    
                                #       totalP<-sum(na.zero(error_kmers[i]), rowSums(fitted_hist[i,]*amlen))
                                   
                                #       prob<-c(na.zero(error_kmers[i]/totalP),as.numeric(fitted_hist[i,]*amlen)/totalP)
                                    
                                #       plot_table<-rbind(plot_table,prob)
                                      
                                #       max.p<-max(prob)  
                                #       readK<-which.max(prob)-1
                                    
                                #       lookup_table<-rbind(lookup_table,c(readK,max.p))
                                    
                                #       }
                                    
                                #       lookup_table<-data.frame(lookup_table)
                                #       rownames(plot_table) <- make.names(plot_table[,1], unique = TRUE)
                                    
                                #       write.table(lookup_table, paste(foldername, "/lookup_table.txt", sep=""), sep=",", row.names = FALSE, col.names = FALSE)

                                #     ## Plot lookup values
                                  
                                #       par(mar=c(5,5,0,1))
                                    
                                #       plot_table<-data.frame(plot_table)
                                #       plot(plot_table$X1,type="n", xlab="Coverage", ylab="Probability",xlim=c(0,akcov*(p*2.5)), ylim=c(0,1), cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size, yaxt='n')

                                #       abline(v=peaks, col="black", lty=2, lwd=0.3)
                                    
                                #       for (i in 1:(ncol(fitted_hist)+1)) {
                                #       lines(plot_table[,i], type="l", lwd=1.5, col=colors[i])
                                #       }
                                #       axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), cex.axis=font_size,
                                #          cex.lab=font_size)
                                  
                                #     dev.off()
                                  
                                # }

                                    ## Some sort of advanced testing (turned off)

                            # if (TESTING) {

                            #   if (TRUE_PARAMS!=-1) {
                            #     testingFile <- paste(foldername, "/", arguments$name_prefix, "SIMULATED_testing.tsv",sep="")
                            #   } else {
                            #     testingFile <- paste(foldername,"/SIMULATED_testing.tsv",sep="")
                            #   }
                            #   if (p==1) {
                            #     cat(paste(amd, akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #   }
                            #   if (p==2) {
                            #     cat(paste(amd, ahets[[1]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #   }
                            #   if (p==3) {
                            #     if (TRUE_PARAMS!=-1) {
                            #       true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #       cat(paste(amd, ahets[[1]], ahets[[2]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #     } else {
                            #       cat(paste(amd, ahets[[1]], ahets[[2]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #     }
                            #   }
                            #   if (p==4) {
                            #     if (topology==0) {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     } else {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, switch(top, ahets[[1]], 0), switch(top, 0, ahets[[1]]), ahets[[2]], ahets[[3]], akcov, adups, atotal_len, top, true_params[1], switch(true_params[5], true_params[2], 0), switch(true_params[5], 0, true_params[2]), true_params[3], true_params[4], true_params[5], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     }
                            #   }
                            #   if (p==5) {
                            #     if (topology==0) {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], true_params[7], true_params[8], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     } else {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, switch(top, ahets[[1]], ahets[[1]], 0, 0, 0), switch(top, 0, 0, ahets[[1]], ahets[[1]], ahets[[1]]), switch(top, ahets[[2]], 0, ahets[[2]], 0, 0), switch(top, 0, ahets[[2]], 0, 0, 0), switch(top, 0, 0, 0, ahets[[2]], ahets[[2]]), switch(top, ahets[[3]], ahets[[3]], ahets[[3]], ahets[[3]], 0), switch(top, 0, 0, 0, 0, ahets[[3]]), ahets[[4]], akcov, adups, atotal_len, top, true_params[1], switch(true_params[6], true_params[2], true_params[2], 0, 0, 0), switch(true_params[6], 0, 0, true_params[2], true_params[2], true_params[2]), switch(true_params[6], true_params[3], 0, true_params[3], 0, 0), switch(true_params[6], 0, true_params[3], 0, 0, 0), switch(true_params[6], 0, 0, 0, true_params[3], true_params[3]), switch(true_params[6], true_params[4], true_params[4], true_params[4], true_params[4], 0), switch(true_params[6], 0, 0, 0, 0, true_params[4]), true_params[5], true_params[6], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     }
                            #   }
                            #   if (p==6) {
                            #     if (topology==0) {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], akcov, adups, atotal_len, top, true_params[1], true_params[2], true_params[3], true_params[4], true_params[5], true_params[6], true_params[7], true_params[8], true_params[9], true_params[10], true_params[11], true_params[12], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], ahets[[7]], ahets[[8]], ahets[[9]], ahets[[10]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     } else {
                            #       if (TRUE_PARAMS!=-1) {
                            #         true_params = unlist(lapply(strsplit(TRUE_PARAMS, ","), as.numeric))
                            #         cat(paste(amd, ifelse(top %in% c(1,2,3,4,5), ahets[[1]], 0), ifelse(top %in% c(6,7,8,9,10,11,12,13), ahets[[1]], 0), ifelse(top %in% c(14,15,16), ahets[[1]], 0), ifelse(top %in% c(1,2,6,7), ahets[[2]], 0), ifelse(top %in% c(3,4,5,14,15,16), ahets[[2]], 0), ifelse(top %in% c(8,9,10), ahets[[2]], 0), ifelse(top %in% c(11,12,13), ahets[[2]], 0), ifelse(top %in% c(1,3,6,8,14), ahets[[3]], 0), ifelse(top %in% c(2,7,11), ahets[[3]], 0), ifelse(top %in% c(4,5,15,16), ahets[[3]], 0), ifelse(top %in% c(9,10,12,13), ahets[[3]], 0), ifelse(top %in% c(1,2,3,4,6,7,8,9,11,12,14,15), ahets[[4]], 0), ifelse(top %in% c(5,16), ahets[[4]], 0), ifelse(top %in% c(10,13), ahets[[4]], 0), ahets[[5]], akcov, adups, atotal_len, top, true_params[1], ifelse(true_params[7] %in% c(1,2,3,4,5), true_params[2], 0), ifelse(true_params[7] %in% c(6,7,8,9,10,11,12,13), true_params[2], 0), ifelse(true_params[7] %in% c(14,15,16), true_params[2], 0), ifelse(true_params[7] %in% c(1,2,6,7), true_params[3], 0), ifelse(true_params[7] %in% c(3,4,5,14,15,16), true_params[3], 0), ifelse(true_params[7] %in% c(8,9,10), true_params[3], 0), ifelse(true_params[7] %in% c(11,12,13), true_params[3], 0), ifelse(true_params[7] %in% c(1,3,6,8,14), true_params[4], 0), ifelse(true_params[7] %in% c(2,7,11), true_params[4], 0), ifelse(true_params[7] %in% c(4,5,15,16), true_params[4], 0), ifelse(true_params[7] %in% c(9,10,12,13), true_params[4], 0), ifelse(true_params[7] %in% c(1,2,3,4,6,7,8,9,11,12,14,15), true_params[5], 0), ifelse(true_params[7] %in% c(5,16), true_params[5], 0), ifelse(true_params[7] %in% c(10,13), true_params[5], 0), true_params[6], true_params[7], sep="\t"), file=testingFile, sep="\n", append=FALSE)
                            #       } else {
                            #         cat(paste(amd, ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], akcov, adups, atotal_len, top, sep="\t"), file=testingFile, sep="\n", append=TRUE)
                            #       }
                            #     }
                            #   }
                            # }

                            }

                    ## Make a ggplot-friendly version of the plots

                        output$plots <- data.frame(
                            kmer_hist_x = kmer_hist_orig[,1],
                            kmer_hist_y = kmer_hist_orig[,2],
                            kmer_hist_unique = unique_hist,
                            kmer_hist_predicted = x**(-transform_exp)*pred,
                            kmer_hits_error = c(error_kmers, rep(0,length(x)-length(error_kmers))),
                            kmer_hist_transformed_x = kmer_hist_transform[,1],
                            kmer_hist_transformed_y = kmer_hist_transform[,2],
                            kmer_hist_unique_transformed = unique_hist_transform,
                            kmer_hist_predicted_transformed = x**(-transform_exp)*pred,
                            kmer_hits_error = c(((x[1:error_xcutoff_ind]**transform_exp)*error_kmers), rep(0,length(x)-length(error_kmers)))
                        )

                    ## the ggplots

                        # thing1 <- ggplot() +
                        #   geom_col(data = kmer_hist_orig, aes(x = V1, y = V2), color = "black", fill = "white") +
                        #   geom_line(data = data.frame(
                        #     x = x, y = unique_hist
                        #   ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x, y = x**(-transform_exp)*pred
                        #   ), aes(x = x, y = y), color = "darkblue", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x[1:error_xcutoff_ind], y = error_kmers
                        #   ), aes(x = x, y = y), color = "red", size = 2, alpha = 0.6) +
                        #   scale_x_continuous(limits = c(0,x_limit_orig)) +
                        #   scale_y_continuous(limits = c(0,y_limit_orig)) +
                        #   theme_bw()

                        # thing2 <- ggplot() +
                        #   geom_col(data = kmer_hist_transform, aes(x = V1, y = V2), color = "black", fill = "white") +
                        #   geom_line(data = data.frame(
                        #     x = x, y = unique_hist_transform
                        #   ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x, y = pred
                        #   ), aes(x = x, y = y), color = "darkblue", size = 2, alpha = 0.6) +
                        #   geom_line(data = data.frame(
                        #     x = x[1:error_xcutoff_ind], y = (x[1:error_xcutoff_ind]**transform_exp)*error_kmers
                        #   ), aes(x = x, y = y), color = "red", size = 2, alpha = 0.6) +
                        #   scale_x_continuous(limits = c(0,x_limit)) +
                        #   scale_y_continuous(limits = c(0,y_limit)) +
                        #   theme_bw()

                        # thing3 <- ggplot() +
                        #   geom_col(data = kmer_hist_orig, aes(x = log(V1), y = log(V2)), color = "black", fill = "white") +
                        #   # geom_line(data = data.frame(
                        #   #   x = x, y = unique_hist
                        #   # ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   # scale_x_continuous(limits = c(0,x_limit_orig)) +
                        #   # scale_y_continuous(limits = c(0,y_limit_orig)) +
                        #   theme_bw()

                        # thing4 <- ggplot() +
                        #   geom_col(data = kmer_hist_transform, aes(x = log(V1), y = log(V2)), color = "black", fill = "white") +
                        #   # geom_line(data = data.frame(
                        #   #   x = x, y = unique_hist_transform
                        #   # ), aes(x = x, y = y), color = "darkgreen", size = 2, alpha = 0.6) +
                        #   # scale_x_continuous(limits = c(0,x_limit)) +
                        #   # scale_y_continuous(limits = c(0,y_limit)) +
                        #   theme_bw()

                        output
            }

