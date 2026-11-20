###### GENOMESCOPE FUNCTIONS #######

## eval_model (medium)

  #' Evaluate distinct model forms, in order to resolve ambiguity of which peak is the homozygous peak
  #'
  #' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
  #' @param nls0,nls1 The nlsLM model objects to evaluate and compare.
  #' @param p An integer corresponding to the ploidy.
  #' @param round An integer corresponding to the iteration number (0, 1, 2, 3) for the fitting process.
  #' @param foldername A character vector corresponding to the name of the output directory.
  #' @param arguments A data frame of the user-specified inputs.
  #' @return A list (nls, nlsscore) where nls is the nlsLM model object (with some additional components)
  #' and nlsscore is the score (model RSSE) corresponding to the best fit (of the p forms).
  #' @export

  eval_model<-function(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments, score_close) {
    
    VERBOSE <- FALSE
    nls0score = -1
    nls1score = -1
    transform_exp = 1
    SCORE_CLOSE <- score_close

    ## Evaluate the score the nls0
    if (!is.null(nls0)) {

        ## below is a non-subfunction implementation of score_model, there is one more, below
        ## score_model<-function(kmer_hist_orig, nls, round, foldername, transform_exp) {

          nls = nls0
          round = round+0.1
          x = kmer_hist_orig[[1]]
          y = kmer_hist_orig[[2]]
          y_transform = as.numeric(x)**transform_exp*as.numeric(y)

          pred=predict(nls, newdata=data.frame(x))
          model_sum=summary(nls)
          p=nls$p
          kcovfloor = max(1, floor(min_max(model_sum$coefficients['kmercov',])[[1]]))

          ## Compute error rate, by counting kmers unexplained by model through first peak
          ## truncate errors as soon as it goes to zero, dont allow it to go back up
            error_xcutoff = kcovfloor
            error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
            if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

          error_kmers = x[1:error_xcutoff_ind]**(-transform_exp)*(y_transform[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind])

          first_zero = -1

          for (i in 1:error_xcutoff_ind) {
            if (first_zero == -1) {
              if (error_kmers[i] < 1.0) {
                first_zero = i
              }
            }
            else {
              error_kmers[i] = 0
            }
          }

          if (first_zero == -1) {
            first_zero = error_xcutoff_ind
          }

          y_fit = y_transform

          ## The fit is residual sum of square error, excluding sequencing errors
            model_fit_all    = c(sum(as.numeric(y_fit[first_zero:length(y_fit)]                          - pred[first_zero:length(y_fit)])                           ** 2), first_zero, x[length(y_fit)])
            model_fit_full   = c(sum(as.numeric(y_fit[first_zero:(min(length(y_fit),(2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]) ** 2), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
            model_fit_unique = c(sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                      - pred[first_zero:((p+1)*kcovfloor)])                       ** 2), first_zero, ((p+1)*kcovfloor))

          ## The score is the percentage of kmers correctly modeled, excluding sequencing errors
            model_fit_allscore    = c(1-sum(abs(as.numeric(y_fit[first_zero:length(y_fit)]                           - pred[first_zero:length(y_fit)])))                           / sum(as.numeric(y_fit[first_zero:length(y_fit)])),                           first_zero, x[length(y_fit)])
            model_fit_fullscore   = c(1-sum(abs(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]))) / sum(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))])), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
            model_fit_uniquescore = c(1-sum(abs(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                       - pred[first_zero:((p+1)*kcovfloor)])))                       / sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)])),                       first_zero, ((p+1)*kcovfloor))

        nls0score = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
                           full = model_fit_full,     fullscore = model_fit_fullscore,
                           unique = model_fit_unique, uniquescore = model_fit_uniquescore)

      # if (VERBOSE) {
      #   mdir = paste(foldername, "/round", round, ".1", sep="")
      #   dir.create(mdir, showWarnings=FALSE)
      #   report_results(kmer_prof_orig,kmer_prof_orig, k, p, (list(nls0, nls0score)) , mdir, arguments, TRUE)
      # }
    }
    else {
      if (VERBOSE) {cat("nls0score failed to converge\n")}
    }

    ## Evaluate the score of nls1
    if (!is.null(nls1)) {

      ## below is a non-subfunction implementation of score_model, there is one more, above
      ## score_model<-function(kmer_hist_orig, nls, round, foldername, transform_exp) {

        nls = nls1
        round = round+0.2
        x = kmer_hist_orig[[1]]
        y = kmer_hist_orig[[2]]
        y_transform = as.numeric(x)**transform_exp*as.numeric(y)

        pred=predict(nls, newdata=data.frame(x))
        model_sum=summary(nls)
        p=nls$p
        kcovfloor = max(1, floor(min_max(model_sum$coefficients['kmercov',])[[1]]))

        ## Compute error rate, by counting kmers unexplained by model through first peak
        ## truncate errors as soon as it goes to zero, dont allow it to go back up
          error_xcutoff = kcovfloor
          error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
          if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

        error_kmers = x[1:error_xcutoff_ind]**(-transform_exp)*(y_transform[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind])

        first_zero = -1

        for (i in 1:error_xcutoff_ind) {
          if (first_zero == -1) {
            if (error_kmers[i] < 1.0) {
              first_zero = i
            }
          }
          else {
            error_kmers[i] = 0
          }
        }

        if (first_zero == -1) {
          first_zero = error_xcutoff_ind
        }

        y_fit = y_transform

        ## The fit is residual sum of square error, excluding sequencing errors
          model_fit_all    = c(sum(as.numeric(y_fit[first_zero:length(y_fit)]                          - pred[first_zero:length(y_fit)])                           ** 2), first_zero, x[length(y_fit)])
          model_fit_full   = c(sum(as.numeric(y_fit[first_zero:(min(length(y_fit),(2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]) ** 2), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
          model_fit_unique = c(sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                      - pred[first_zero:((p+1)*kcovfloor)])                       ** 2), first_zero, ((p+1)*kcovfloor))

        ## The score is the percentage of kmers correctly modeled, excluding sequencing errors
          model_fit_allscore    = c(1-sum(abs(as.numeric(y_fit[first_zero:length(y_fit)]                           - pred[first_zero:length(y_fit)])))                           / sum(as.numeric(y_fit[first_zero:length(y_fit)])),                           first_zero, x[length(y_fit)])
          model_fit_fullscore   = c(1-sum(abs(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]))) / sum(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))])), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
          model_fit_uniquescore = c(1-sum(abs(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                       - pred[first_zero:((p+1)*kcovfloor)])))                       / sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)])),                       first_zero, ((p+1)*kcovfloor))

      nls1score = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
                         full = model_fit_full,     fullscore = model_fit_fullscore,
                         unique = model_fit_unique, uniquescore = model_fit_uniquescore)

      if(VERBOSE) {cat(paste("nls1score$all:\t", nls1score$all[[1]], "\n"))}

      # if (VERBOSE) {
      #   mdir = paste(foldername, "/round", round, ".2", sep="")
      #   dir.create(mdir, showWarnings=FALSE)
      #   report_results(kmer_prof_orig, kmer_prof_orig, k, p, (list(nls1, nls1score)) , mdir, arguments, FALSE)
      # }
    }
    else {
      if (VERBOSE) {cat("nls1score failed to converge\n")}
    }

    ## Return the better of the scores
    if (!is.null(nls0)) {
      if (!is.null(nls1)) {
        pdiff = abs(nls0score$all[[1]] - nls1score$all[[1]]) / max(nls0score$all[[1]], nls1score$all[[1]])

        if (pdiff < SCORE_CLOSE) {
          het0 = nls0$ahet
          het1 = nls1$ahet

          #if (het1 * SCORE_HET_FOLD_DIFFERENCE < het0) {
          if (het1 + 0.01 < het0) {
            if (VERBOSE) {cat("returning nls0, similar score, higher het\n")}
            return (list(nls1, nls1score))
          }
          #else if (het0 * SCORE_HET_FOLD_DIFFERENCE < het1) {
          else if (het0  + 0.01 < het1) {
            if (VERBOSE) {cat("returning nls1, similar score, higher het\n")}
            return (list(nls0, nls0score))
          }
        }

        if (nls0score$all[[1]] < nls1score$all[[1]]) {
          if (VERBOSE) {cat("returning nls0, better score\n")}
          return (list(nls0, nls0score))
        }
        else {
          if (VERBOSE) {cat("returning nls1, better score\n")}
          return (list(nls1, nls1score))
        }
      }
      else {
        if (VERBOSE) {cat("returning nls0, nls1 fail\n")}
        return (list(nls0, nls0score))
      }
    }

    if (VERBOSE) {cat("returning nls1 by default\n")}
    return (list(nls1, nls1score))
  }

## predict functions (big)

  #' Produce model estimated (p=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict1_0 = function(k, d, kmercov, bias, x)
  {
    predict1_1(k, d, kmercov, bias, x)
  }

  #' Produce model estimated (p=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict1_0_unique = function(k, d, kmercov, bias, x)
  {
    predict1_1_unique(k, d, kmercov, bias, x)
  }

  #' Produce model estimated (p=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict1_1 = function(k, d, kmercov, bias, x)
  {
    r0 = 1
    if (d > 1) {return(0)}
    t0 = r0**k
    s0 = t0
    alpha_1 = (1-d)*(s0)
    alpha_2 = d*(s0)
    alpha_1 * dnbinom(x, size=kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size=kmercov*2 / bias, mu = kmercov*2)
  }

  #' Produce model estimated (p=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict1_1_unique = function(k, d, kmercov, bias, x)
  {
    r0 = 1
    if (d > 1) {return(0)}
    t0 = r0**k
    s0 = t0
    alpha_1_unique = (1-d)*(s0)
    alpha_1_unique * dnbinom(x, size=kmercov*1 / bias, mu = kmercov*1)
  }

  #AB -> AA
  #' Produce model estimated (p=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1 A numeric corresponding to the nucleotide heterozygosity ab.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict2_0 = function(r1, k, d, kmercov, bias, x)
  {
    predict2_1(r1, k, d, kmercov, bias, x)
  }

  #AB -> AA
  #' Produce model estimated (p=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1 A numeric corresponding to the nucleotide heterozygosity ab.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict2_0_unique = function(r1, k, d, kmercov, bias, x)
  {
    predict2_1_unique(r1, k, d, kmercov, bias, x)
  }

  #AB -> AA
  #' Produce model estimated (p=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1 A numeric corresponding to the nucleotide heterozygosity ab.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict2_1 = function(r1, k, d, kmercov, bias, x)
  {
    r0 = 1-r1 #aa
    if (r0 < 0 || d > 1) {return(0)}
    t0 = r0**k #AA
    s0 = t0 #AA
    s1 = 1-t0 #AB
    alpha_1 = (1-d)*(2*s1) + d*(2*s0*s1 + 2*s1**2)
    alpha_2 = (1-d)*(s0) + d*(s1**2)
    alpha_3 = d*(2*s0*s1)
    alpha_4 = d*(s0**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
  }

  #AB -> AA
  #' Produce model estimated (p=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1 A numeric corresponding to the nucleotide heterozygosity ab.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict2_1_unique = function(r1, k, d, kmercov, bias, x)
  {
    r0 = 1-r1
    if (r0 < 0 || d > 1) {return(0)}
    t0 = r0**k
    s0 = t0
    s1 = 1-t0
    alpha_1_unique = (1-d)*(2*s1)
    alpha_2_unique = (1-d)*(s0)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)
  }

  #ABC -> AAB -> AAA
  #' Produce model estimated (p=3) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1,r2 Numerics corresponding to the nucleotide heterozygosities aab and abc respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict3_0 = function(r1, r2, k, d, kmercov, bias, x)
  {
    predict3_1(r1, r2, k, d, kmercov, bias, x)
  }

  #ABC -> AAB -> AAA
  #' Produce model estimated (p=3, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1,r2 Numerics corresponding to the nucleotide heterozygosities aab and abc respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict3_0_unique = function(r1, r2, k, d, kmercov, bias, x)
  {
    predict3_1_unique(r1, r2, k, d, kmercov, bias, x)
  }

  #ABC -> AAB -> AAA
  #' Produce model estimated (p=3) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1,r2 Numerics corresponding to the nucleotide heterozygosities aab and abc respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict3_1 = function(r1, r2, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2 #aaa
    if (r0 < 0 || d > 1) {return(0)}
    t0 = r0**k #AAA
    t1 = (r0+r1)**k #AAA + AAB
    s0 = t0 #AAA
    s1 = t1-t0 #AAB
    s2 = 1-t1 #ABC
    alpha_1 = (1-d)*(s1+3*s2) + d*(2*s0*s1 + 4*s0*s2 + 2*s1**2 + 6*s1*s2 + 4*s2**2)
    alpha_2 = (1-d)*(s1) + d*(s2**2)
    alpha_3 = (1-d)*(s0) + d*(2*s1*s2)
    alpha_4 = d*(2*s0*s2 + s1**2)
    alpha_5 = d*(2*s0*s1)
    alpha_6 = d*(s0**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #ABC -> AAB -> AAA
  #' Produce model estimated (p=3, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param r1,r2 Numerics corresponding to the nucleotide heterozygosities aab and abc respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict3_1_unique = function(r1, r2, k, d, kmercov, bias, x)
  {
    r0 = 1-r1-r2
    if (r0 < 0 || d > 1) {return(0)}
    t0 = r0**k
    t1 = (r0+r1)**k
    s0 = t0
    s1 = t1-t0
    s2 = 1-t1
    alpha_1_unique = (1-d)*(s1+3*s2)
    alpha_2_unique = (1-d)*(s1)
    alpha_3_unique = (1-d)*(s0)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)
  }

  #AAAA -> (AAAB, AABB) -> AABC -> ABCD
  #' Produce model estimated (p=4, full model) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaab,raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabb, aabc, and abcd respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict4_0 = function(raaab, raabb, raabc, rabcd, k, d, kmercov, bias, x)
  {
    raaaa = 1-raaab-raabb-raabc-rabcd
    if (raaaa < 0 || d > 1) {return(0)}
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raaab+raabb+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAAAB-tAABB+tAAAA
    sABCD = 1-tAABC
    alpha_1 = (1-d)*(sAAAB + 2*sAABC + 4*sABCD) + d*(2*sAAAA*sAAAB + 2*sAAAB**2 + 2*sAAAB*sAABB + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABB*sAABC + 4*sAABC**2 + 6*sAAAA*sABCD + 8*sAAAB*sABCD + 6*sAABB*sABCD + 10*sAABC*sABCD + 6*sABCD**2)
    alpha_2 = (1-d)*(2*sAABB + sAABC) + d*(2*sAAAA*sAABB + 2*sAAAB*sAABB + 2*sAABB**2 + 2*sAABB*sAABC + 2*sAABB*sABCD + sABCD**2)
    alpha_3 = (1-d)*(sAAAB) + d*(2*sAABB*sABCD + 2*sAABC*sABCD)
    alpha_4 = (1-d)*(sAAAA) + d*(sAABB**2 + 2*sAABB*sAABC + sAABC**2 + 2*sAAAB*sABCD)
    alpha_5 = d*(2*sAAAB*sAABB + 2*sAAAB*sAABC + 2*sAAAA*sABCD)
    alpha_6 = d*(sAAAB**2 + 2*sAAAA*sAABB + 2*sAAAA*sAABC)
    alpha_7 = d*(2*sAAAA*sAAAB)
    alpha_8 = d*(sAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
  }

  #AAAA -> (AAAB, AABB) -> AABC -> ABCD
  #' Produce model estimated (p=4, full model, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaab,raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabb, aabc, and abcd respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict4_0_unique = function(raaab, raabb, raabc, rabcd, k, d, kmercov, bias, x)
  {
    raaaa = 1-raaab-raabb-raabc-rabcd
    if (raaaa < 0 || d > 1) {return(0)}
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raaab+raabb+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAAAB-tAABB+tAAAA
    sABCD = 1-tAABC
    alpha_1_unique = (1-d)*(sAAAB + 2*sAABC + 4*sABCD)
    alpha_2_unique = (1-d)*(2*sAABB + sAABC)
    alpha_3_unique = (1-d)*(sAAAB)
    alpha_4_unique = (1-d)*(sAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
  }

  #AAAA -> AAAB -> AABC -> ABCD
  #' Produce model estimated (p=4, topology=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaab,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabc, and abcd respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict4_1 = function(raaab, raabc, rabcd, k, d, kmercov, bias, x)
  {
    raaaa = 1-raaab-raabc-rabcd
    if (raaaa < 0 || d > 1) {return(0)}
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABC = (raaaa+raaab+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABC = tAABC-tAAAB
    sABCD = 1-tAABC
    alpha_1 = (1-d)*(sAAAB + 2*sAABC + 4*sABCD) + d*(2*sAAAA*sAAAB + 2*sAAAB**2 + 4*sAAAA*sAABC + 6*sAAAB*sAABC + 4*sAABC**2 + 6*sAAAA*sABCD + 8*sAAAB*sABCD + 10*sAABC*sABCD + 6*sABCD**2)
    alpha_2 = (1-d)*(sAABC) + d*(sABCD**2)
    alpha_3 = (1-d)*(sAAAB) + d*(2*sAABC*sABCD)
    alpha_4 = (1-d)*(sAAAA) + d*(sAABC**2 + 2*sAAAB*sABCD)
    alpha_5 = d*(2*sAAAB*sAABC + 2*sAAAA*sABCD)
    alpha_6 = d*(sAAAB**2 + 2*sAAAA*sAABC)
    alpha_7 = d*(2*sAAAA*sAAAB)
    alpha_8 = d*(sAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
  }

  #AAAA -> AAAB -> AABC -> ABCD
  #' Produce model estimated (p=4, topology=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaab,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aaab, aabc, and abcd respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict4_1_unique = function(raaab, raabc, rabcd, k, d, kmercov, bias, x)
  {
    raaaa = 1-raaab-raabc-rabcd
    if (raaaa < 0 || d > 1) {return(0)}
    tAAAA = raaaa**k
    tAAAB = (raaaa+raaab)**k
    tAABC = (raaaa+raaab+raabc)**k
    sAAAA = tAAAA
    sAAAB = tAAAB-tAAAA
    sAABC = tAABC-tAAAB
    sABCD = 1-tAABC
    alpha_1_unique = (1-d)*(sAAAB + 2*sAABC + 4*sABCD)
    alpha_2_unique = (1-d)*(sAABC)
    alpha_3_unique = (1-d)*(sAAAB)
    alpha_4_unique = (1-d)*(sAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
  }

  #AAAA -> AABB -> AABC -> ABCD
  #' Produce model estimated (p=4, topology=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aabb, aabc, and abcd respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict4_2 = function(raabb, raabc, rabcd, k, d, kmercov, bias, x)
  {
    raaaa = 1-raabb-raabc-rabcd
    if (raaaa < 0 || d > 1) {return(0)}
    tAAAA = raaaa**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raabb+raabc)**k
    sAAAA = tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAABB
    sABCD = 1-tAABC
    alpha_1 = (1-d)*(2*sAABC + 4*sABCD) + d*(4*sAAAA*sAABC + 4*sAABB*sAABC + 4*sAABC**2 + 6*sAAAA*sABCD + 6*sAABB*sABCD + 10*sAABC*sABCD + 6*sABCD**2)
    alpha_2 = (1-d)*(2*sAABB + sAABC) + d*(2*sAAAA*sAABB + 2*sAABB**2 + 2*sAABB*sAABC + 2*sAABB*sABCD + sABCD**2)
    alpha_3 = (1-d)*(0) + d*(2*sAABB*sABCD + 2*sAABC*sABCD)
    alpha_4 = (1-d)*(sAAAA) + d*(sAABB**2 + 2*sAABB*sAABC + sAABC**2)
    alpha_5 = d*(2*sAAAA*sABCD)
    alpha_6 = d*(2*sAAAA*sAABB + 2*sAAAA*sAABC)
    alpha_7 = d*(0)
    alpha_8 = d*(sAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)
  }

  #AAAA -> AABB -> AABC -> ABCD
  #' Produce model estimated (p=4, topology=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raabb,raabc,rabcd Numerics corresponding to the nucleotide heterozygosities aabb, aabc, and abcd respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict4_2_unique = function(raabb, raabc, rabcd, k, d, kmercov, bias, x)
  {
    raaaa = 1-raabb-raabc-rabcd
    if (raaaa < 0 || d > 1) {return(0)}
    tAAAA = raaaa**k
    tAABB = (raaaa+raabb)**k
    tAABC = (raaaa+raabb+raabc)**k
    sAAAA = tAAAA
    sAABB = tAABB-tAAAA
    sAABC = tAABC-tAABB
    sABCD = 1-tAABC
    alpha_1_unique = (1-d)*(2*sAABC + 4*sABCD)
    alpha_2_unique = (1-d)*(2*sAABB + sAABC)
    alpha_3_unique = (1-d)*(0)
    alpha_4_unique = (1-d)*(sAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)
  }

  #' Produce model estimated (p=5, full model) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaab,raaabb,raaabc,raabbc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_0 = function(raaaab, raaabb, raaabc, raabbc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaaab-raaabb-raaabc-raabbc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = (raaaaa)**k
    tAAAAB = (raaaaa+raaaab)**k
    tAAABB = (raaaaa+raaabb)**k
    tAAABC = (raaaaa+raaaab+raaabb+raaabc)**k
    tAABBC = (raaaaa+raaaab+raabbc)**k
    tAABCD = (raaaaa+raaaab+raaabb+raaabc+raabbc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAAABC = tAAABC-tAAAAB-tAAABB+tAAAAA
    sAABBC = tAABBC-tAAAAB
    sAABCD = tAABCD-tAAABC-tAABBC+tAAAAB
    sABCDE = 1-tAABCD
    alpha_1  = (1-d)*(sAAAAB + 2*sAAABC + sAABBC + 3*sAABCD + 5*sABCDE) + d*(2*sAAAAA*sAAAAB + 2*sAAAAB**2 + 2*sAAAAB*sAAABB + 4*sAAAAA*sAAABC + 6*sAAAAB*sAAABC + 4*sAAABB*sAAABC + 4*sAAABC**2 + 2*sAAAAA*sAABBC + 4*sAAAAB*sAABBC + 2*sAAABB*sAABBC + 6*sAAABC*sAABBC + 2*sAABBC**2 + 6*sAAAAA*sAABCD + 8*sAAAAB*sAABCD + 6*sAAABB*sAABCD + 10*sAAABC*sAABCD + 8*sAABBC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 10*sAAAAB*sABCDE + 8*sAAABB*sABCDE + 12*sAAABC*sABCDE + 10*sAABBC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
    alpha_2  = (1-d)*(sAAABB + 2*sAABBC + sAABCD) + d*(2*sAAAAA*sAAABB + 2*sAAAAB*sAAABB + 2*sAAABB**2 + 2*sAAABB*sAAABC + 2*sAAAAA*sAABBC + 2*sAAAAB*sAABBC + 4*sAAABB*sAABBC + 2*sAAABC*sAABBC + 2*sAABBC**2 + 2*sAAABB*sAABCD + 2*sAABBC*sAABCD + 2*sAAABB*sABCDE + 2*sAABBC*sABCDE + sABCDE**2)
    alpha_3  = (1-d)*(sAAABB + sAAABC) + d*(2*sAABBC*sABCDE + 2*sAABCD*sABCDE)
    alpha_4  = (1-d)*(sAAAAB) + d*(sAABBC**2 + 2*sAABBC*sAABCD + sAABCD**2 + 2*sAAABB*sABCDE + 2*sAAABC*sABCDE)
    alpha_5  = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABBC + 2*sAAABC*sAABBC + 2*sAAABB*sAABCD + 2*sAAABC*sAABCD + 2*sAAAAB*sABCDE)
    alpha_6  = d*(sAAABB**2 + 2*sAAABB*sAAABC + sAAABC**2 + 2*sAAAAB*sAABBC + 2*sAAAAB*sAABCD + 2*sAAAAA*sABCDE)
    alpha_7  = d*(2*sAAAAB*sAAABB + 2*sAAAAB*sAAABC + 2*sAAAAA*sAABBC + 2*sAAAAA*sAABCD)
    alpha_8  = d*(sAAAAB**2 + 2*sAAAAA*sAAABB + 2*sAAAAA*sAAABC)
    alpha_9  = d*(2*sAAAAA*sAAAAB)
    alpha_10 = d*(sAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)
  }

  #' Produce model estimated (p=5, full model, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaab,raaabb,raaabc,raabbc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_0_unique = function(raaaab, raaabb, raaabc, raabbc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaaab-raaabb-raaabc-raabbc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = (raaaaa)**k
    tAAAAB = (raaaaa+raaaab)**k
    tAAABB = (raaaaa+raaabb)**k
    tAAABC = (raaaaa+raaaab+raaabb+raaabc)**k
    tAABBC = (raaaaa+raaaab+raabbc)**k
    tAABCD = (raaaaa+raaaab+raaabb+raaabc+raabbc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAAABC = tAAABC-tAAAAB-tAAABB+tAAAAA
    sAABBC = tAABBC-tAAAAB
    sAABCD = tAABCD-tAAABC-tAABBC+tAAAAB
    sABCDE = 1-tAABCD
    alpha_1_unique  = (1-d)*(sAAAAB + 2*sAAABC + sAABBC + 3*sAABCD + 5*sABCDE)
    alpha_2_unique  = (1-d)*(sAAABB + 2*sAABBC + sAABCD)
    alpha_3_unique  = (1-d)*(sAAABB + sAAABC)
    alpha_4_unique  = (1-d)*(sAAAAB)
    alpha_5_unique  = (1-d)*(sAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
  }

  #' Produce model estimated (p=5, topology=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaab,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aaabc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_1 = function(raaaab, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaaab-raaabc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAAABC = (raaaaa+raaaab+raaabc)**k
    tAABCD = (raaaaa+raaaab+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAAABC = tAAABC-tAAAAB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
    alpha_1 = (1-d)*(sAAAAB + 2*sAAABC + 3*sAABCD + 5*sABCDE) + d*(2*sAAAAA*sAAAAB + 2*sAAAAB**2 + 4*sAAAAA*sAAABC + 6*sAAAAB*sAAABC + 4*sAAABC**2 + 6*sAAAAA*sAABCD + 8*sAAAAB*sAABCD + 10*sAAABC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 10*sAAAAB*sABCDE + 12*sAAABC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
    alpha_2 = (1-d)*(sAABCD) + d*(sABCDE**2)
    alpha_3 = (1-d)*(sAAABC) + d*(2*sAABCD*sABCDE)
    alpha_4 = (1-d)*(sAAAAB) + d*(sAABCD**2 + 2*sAAABC*sABCDE)
    alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABC*sAABCD + 2*sAAAAB*sABCDE)
    alpha_6 = d*(sAAABC**2 + 2*sAAAAB*sAABCD + 2*sAAAAA*sABCDE)
    alpha_7 = d*(2*sAAAAB*sAAABC + 2*sAAAAA*sAABCD)
    alpha_8 = d*(sAAAAB**2 + 2*sAAAAA*sAAABC)
    alpha_9 = d*(2*sAAAAA*sAAAAB)
    alpha_10= d*(sAAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
  }

  #' Produce model estimated (p=5, topology=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaab,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aaabc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_1_unique = function(raaaab, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaaab-raaabc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAAABC = (raaaaa+raaaab+raaabc)**k
    tAABCD = (raaaaa+raaaab+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAAABC = tAAABC-tAAAAB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
    alpha_1_unique = (1-d)*(sAAAAB + 2*sAAABC + 3*sAABCD + 5*sABCDE)
    alpha_2_unique = (1-d)*(sAABCD)
    alpha_3_unique = (1-d)*(sAAABC)
    alpha_4_unique = (1-d)*(sAAAAB)
    alpha_5_unique = (1-d)*(sAAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
  }

  #' Produce model estimated (p=5, topology=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaab,raabbc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aabbc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_2 = function(raaaab, raabbc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaaab-raabbc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAABBC = (raaaaa+raaaab+raabbc)**k
    tAABCD = (raaaaa+raaaab+raabbc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAABBC = tAABBC-tAAAAB
    sAABCD = tAABCD-tAABBC
    sABCDE = 1-tAABCD
    alpha_1 = (1-d)*(sAAAAB + sAABBC + 3*sAABCD + 5*sABCDE) + d*(2*sAAAAA*sAAAAB + 2*sAAAAB**2 + 2*sAAAAA*sAABBC + 4*sAAAAB*sAABBC + 2*sAABBC**2 + 6*sAAAAA*sAABCD + 8*sAAAAB*sAABCD + 8*sAABBC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 10*sAAAAB*sABCDE + 10*sAABBC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
    alpha_2 = (1-d)*(2*sAABBC + sAABCD) + d*(2*sAAAAA*sAABBC + 2*sAAAAB*sAABBC + 2*sAABBC**2 + 2*sAABBC*sAABCD + 2*sAABBC*sABCDE + sABCDE**2)
    alpha_3 = (1-d)*(0) + d*(2*sAABBC*sABCDE + 2*sAABCD*sABCDE)
    alpha_4 = (1-d)*(sAAAAB) + d*(sAABBC**2 + 2*sAABBC*sAABCD + sAABCD**2)
    alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAAAB*sABCDE)
    alpha_6 = d*(2*sAAAAB*sAABBC + 2*sAAAAB*sAABCD + 2*sAAAAA*sABCDE)
    alpha_7 = d*(2*sAAAAA*sAABBC + 2*sAAAAA*sAABCD)
    alpha_8 = d*(sAAAAB**2)
    alpha_9 = d*(2*sAAAAA*sAAAAB)
    alpha_10= d*(sAAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
  }

  #' Produce model estimated (p=5, topology=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaab,raabbc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaaab, aabbc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_2_unique = function(raaaab, raabbc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaaab-raabbc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAAAB = (raaaaa+raaaab)**k
    tAABBC = (raaaaa+raaaab+raabbc)**k
    tAABCD = (raaaaa+raaaab+raabbc+raabcd)**k
    sAAAAA = tAAAAA
    sAAAAB = tAAAAB-tAAAAA
    sAABBC = tAABBC-tAAAAB
    sAABCD = tAABCD-tAABBC
    sABCDE = 1-tAABCD
    alpha_1_unique = (1-d)*(sAAAAB + sAABBC + 3*sAABCD + 5*sABCDE)
    alpha_2_unique = (1-d)*(2*sAABBC + sAABCD)
    alpha_3_unique = (1-d)*(0)
    alpha_4_unique = (1-d)*(sAAAAB)
    alpha_5_unique = (1-d)*(sAAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
  }

  #' Produce model estimated (p=5, topology=3) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabb,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aaabc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_3 = function(raaabb, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaabb-raaabc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAAABC = (raaaaa+raaabb+raaabc)**k
    tAABCD = (raaaaa+raaabb+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAAABC = tAAABC-tAAABB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
    alpha_1 = (1-d)*(2*sAAABC + 3*sAABCD + 5*sABCDE) + d*(4*sAAAAA*sAAABC + 4*sAAABB*sAAABC + 4*sAAABC**2 + 6*sAAAAA*sAABCD + 6*sAAABB*sAABCD + 10*sAAABC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 8*sAAABB*sABCDE + 12*sAAABC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
    alpha_2 = (1-d)*(sAAABB + sAABCD) + d*(2*sAAAAA*sAAABB + 2*sAAABB**2 + 2*sAAABB*sAAABC + 2*sAAABB*sAABCD + 2*sAAABB*sABCDE + sABCDE**2)
    alpha_3 = (1-d)*(sAAABB + sAAABC) + d*(2*sAABCD*sABCDE)
    alpha_4 = (1-d)*(0) + d*(sAABCD**2 + 2*sAAABB*sABCDE + 2*sAAABC*sABCDE)
    alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABCD + 2*sAAABC*sAABCD)
    alpha_6 = d*(sAAABB**2 + 2*sAAABB*sAAABC + sAAABC**2 + 2*sAAAAA*sABCDE)
    alpha_7 = d*(2*sAAAAA*sAABCD)
    alpha_8 = d*(2*sAAAAA*sAAABB + 2*sAAAAA*sAAABC)
    alpha_9 = d*(0)
    alpha_10= d*(sAAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
  }

  #' Produce model estimated (p=5, topology=3, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabb,raaabc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aaabc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_3_unique = function(raaabb, raaabc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaabb-raaabc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAAABC = (raaaaa+raaabb+raaabc)**k
    tAABCD = (raaaaa+raaabb+raaabc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAAABC = tAAABC-tAAABB
    sAABCD = tAABCD-tAAABC
    sABCDE = 1-tAABCD
    alpha_1_unique = (1-d)*(2*sAAABC + 3*sAABCD + 5*sABCDE)
    alpha_2_unique = (1-d)*(sAAABB + sAABCD)
    alpha_3_unique = (1-d)*(sAAABB + sAAABC)
    alpha_4_unique = (1-d)*(0)
    alpha_5_unique = (1-d)*(sAAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
  }

  #' Produce model estimated (p=5, topology=4) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabb,raabcc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_4 = function(raaabb, raabcc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaabb-raabcc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tAABCD = (raaaaa+raaabb+raabcc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sAABCD = tAABCD-tAABCC
    sABCDE = 1-tAABCD
    alpha_1 = (1-d)*(sAABCC + 3*sAABCD + 5*sABCDE) + d*(2*sAAAAA*sAABCC + 2*sAAABB*sAABCC + 2*sAABCC**2 + 6*sAAAAA*sAABCD + 6*sAAABB*sAABCD + 8*sAABCC*sAABCD + 6*sAABCD**2 + 8*sAAAAA*sABCDE + 8*sAAABB*sABCDE + 10*sAABCC*sABCDE + 14*sAABCD*sABCDE + 8*sABCDE**2)
    alpha_2 = (1-d)*(sAAABB + 2*sAABCC + sAABCD) + d*(2*sAAAAA*sAAABB + 2*sAAABB**2 + 2*sAAAAA*sAABCC + 4*sAAABB*sAABCC + 2*sAABCC**2 + 2*sAAABB*sAABCD + 2*sAABCC*sAABCD + 2*sAAABB*sABCDE + 2*sAABCC*sABCDE + sABCDE**2)
    alpha_3 = (1-d)*(sAAABB) + d*(2*sAABCC*sABCDE + 2*sAABCD*sABCDE)
    alpha_4 = (1-d)*(0) + d*(sAABCC**2 + 2*sAABCC*sAABCD + sAABCD**2 + 2*sAAABB*sABCDE)
    alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABCC + 2*sAAABB*sAABCD)
    alpha_6 = d*(sAAABB**2 + 2*sAAAAA*sABCDE)
    alpha_7 = d*(2*sAAAAA*sAABCC + 2*sAAAAA*sAABCD)
    alpha_8 = d*(2*sAAAAA*sAAABB)
    alpha_9 = d*(0)
    alpha_10= d*(sAAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
  }

  #' Produce model estimated (p=5, topology=4, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabb,raabcc,raabcd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, aabcd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_4_unique = function(raaabb, raabcc, raabcd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaabb-raabcc-raabcd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tAABCD = (raaaaa+raaabb+raabcc+raabcd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sAABCD = tAABCD-tAABCC
    sABCDE = 1-tAABCD
    alpha_1_unique = (1-d)*(sAABCC + 3*sAABCD + 5*sABCDE)
    alpha_2_unique = (1-d)*(sAAABB + 2*sAABCC + sAABCD)
    alpha_3_unique = (1-d)*(sAAABB)
    alpha_4_unique = (1-d)*(0)
    alpha_5_unique = (1-d)*(sAAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
  }

  #' Produce model estimated (p=5, topology=5) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabb,raabcc,rabcdd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, abcdd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_5 = function(raaabb, raabcc, rabcdd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaabb-raabcc-rabcdd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tABCDD = (raaaaa+raaabb+raabcc+rabcdd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sABCDD = tABCDD-tAABCC
    sABCDE = 1-tABCDD
    alpha_1 = (1-d)*(sAABCC + 3*sABCDD + 5*sABCDE) + d*(2*sAAAAA*sAABCC + 2*sAAABB*sAABCC + 2*sAABCC**2 + 4*sAAAAA*sABCDD + 4*sAAABB*sABCDD + 6*sAABCC*sABCDD + 4*sABCDD**2 + 8*sAAAAA*sABCDE + 8*sAAABB*sABCDE + 10*sAABCC*sABCDE + 12*sABCDD*sABCDE + 8*sABCDE**2)
    alpha_2 = (1-d)*(sAAABB + 2*sAABCC + sABCDD) + d*(2*sAAAAA*sAAABB + 2*sAAABB**2 + 2*sAAAAA*sAABCC + 4*sAAABB*sAABCC + 2*sAABCC**2 + 2*sAAAAA*sABCDD + 4*sAAABB*sABCDD + 4*sAABCC*sABCDD + 3*sABCDD**2 + 2*sAAABB*sABCDE + 2*sAABCC*sABCDE + 4*sABCDD*sABCDE + sABCDE**2)
    alpha_3 = (1-d)*(sAAABB) + d*(2*sAABCC*sABCDD + 2*sAABCC*sABCDE)
    alpha_4 = (1-d)*(0) + d*(sAABCC**2 + 2*sAAABB*sABCDD + 2*sAAABB*sABCDE)
    alpha_5 = (1-d)*(sAAAAA) + d*(2*sAAABB*sAABCC)
    alpha_6 = d*(sAAABB**2 + 2*sAAAAA*sABCDD + 2*sAAAAA*sABCDE)
    alpha_7 = d*(2*sAAAAA*sAABCC)
    alpha_8 = d*(2*sAAAAA*sAAABB)
    alpha_9 = d*(0)
    alpha_10= d*(sAAAAA**2)
    alpha_1 * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2 * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3 * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4 * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5 * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6 * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7 * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8 * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9 * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10* dnbinom(x, size = kmercov*10/ bias, mu = kmercov*10)
  }

  #' Produce model estimated (p=5, topology=5, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabb,raabcc,rabcdd,rabcde Numerics corresponding to the nucleotide heterozygosities aaabb, aabcc, abcdd, and abcde respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict5_5_unique = function(raaabb, raabcc, rabcdd, rabcde, k, d, kmercov, bias, x)
  {
    raaaaa = 1-raaabb-raabcc-rabcdd-rabcde
    if (raaaaa < 0 || d > 1) {return(0)}
    tAAAAA = raaaaa**k
    tAAABB = (raaaaa+raaabb)**k
    tAABCC = (raaaaa+raaabb+raabcc)**k
    tABCDD = (raaaaa+raaabb+raabcc+rabcdd)**k
    sAAAAA = tAAAAA
    sAAABB = tAAABB-tAAAAA
    sAABCC = tAABCC-tAAABB
    sABCDD = tABCDD-tAABCC
    sABCDE = 1-tABCDD
    alpha_1_unique = (1-d)*(sAABCC + 3*sABCDD + 5*sABCDE)
    alpha_2_unique = (1-d)*(sAAABB + 2*sAABCC + sABCDD)
    alpha_3_unique = (1-d)*(sAAABB)
    alpha_4_unique = (1-d)*(0)
    alpha_5_unique = (1-d)*(sAAAAA)
    alpha_1_unique * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)
  }

  #' Produce model estimated (p=6, full model) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaaabb,raaabbb,raaaabc,raaabbc,raabbcc,raaabcd,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_0 = function(raaaaab, raaaabb, raaabbb, raaaabc, raaabbc, raabbcc, raaabcd, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa = 1-raaaaab-raaaabb-raaabbb-raaaabc-raaabbc-raabbcc-raaabcd-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    tAAAAAA = (raaaaaa)**k
    tAAAAAB = (raaaaaa+raaaaab)**k
    tAAAABB = (raaaaaa+raaaabb)**k
    tAAABBB = (raaaaaa+raaabbb)**k
    tAAAABC = (raaaaaa+raaaaab+raaaabb+raaaabc)**k
    tAAABBC = (raaaaaa+raaaaab+raaabbb+raaabbc)**k
    tAABBCC = (raaaaaa+raaaabb+raabbcc)**k
    tAAABCD = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raaabcd)**k
    tAABBCD = (raaaaaa+raaaaab+raaaabb+raaaabc+raabbcc+raabbcd)**k
    tAABCDE = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raabbcc+raaabcd+raabbcd+raabcde)**k
    sAAAAAA = tAAAAAA
    sAAAAAB = tAAAAAB-tAAAAAA
    sAAAABB = tAAAABB-tAAAAAA
    sAAABBB = tAAABBB-tAAAAAA
    sAAAABC = tAAAABC-tAAAABB-tAAAAAB+tAAAAAA
    sAAABBC = tAAABBC-tAAABBB-tAAAAAB+tAAAAAA
    sAABBCC = tAABBCC-tAAAABB
    sAAABCD = tAAABCD-tAAABBC-tAAAABC+tAAAAAB
    sAABBCD = tAABBCD-tAABBCC-tAAAABC+tAAAABB
    sAABCDE = tAABCDE-tAABBCD-tAAABCD+tAAAABC
    sABCDEF = 1-tAABCDE
    alpha_1 = (1-d)*(sAAAAAB + 2*sAAAABC + sAAABBC + 3*sAAABCD + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAB*sAAAABB + 4*sAAAAAA*sAAAABC + 6*sAAAAAB*sAAAABC + 4*sAAAABB*sAAAABC + 4*sAAAABC**2 + 2*sAAAAAB*sAAABBB + 4*sAAAABC*sAAABBB + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAAABB*sAAABBC + 6*sAAAABC*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 6*sAAAAAA*sAAABCD + 8*sAAAAAB*sAAABCD + 6*sAAAABB*sAAABCD + 10*sAAAABC*sAAABCD + 6*sAAABBB*sAAABCD + 8*sAAABBC*sAAABCD + 6*sAAABCD**2 + 2*sAAAAAB*sAABBCC + 4*sAAAABC*sAABBCC + 2*sAAABBC*sAABBCC + 6*sAAABCD*sAABBCC + 4*sAAAAAA*sAABBCD + 6*sAAAAAB*sAABBCD + 4*sAAAABB*sAABBCD + 8*sAAAABC*sAABBCD + 4*sAAABBB*sAABBCD + 6*sAAABBC*sAABBCD + 10*sAAABCD*sAABBCD + 4*sAABBCC*sAABBCD + 4*sAABBCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 8*sAAAABB*sAABCDE + 12*sAAAABC*sAABCDE + 8*sAAABBB*sAABCDE + 10*sAAABBC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABBCC*sAABCDE + 12*sAABBCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 10*sAAAABB*sABCDEF + 14*sAAAABC*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 16*sAAABCD*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABBCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2 = (1-d)*(sAAAABB + sAAABBC + 3*sAABBCC + 2*sAABBCD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAAAB*sAAAABB + 2*sAAAABB**2 + 2*sAAAABB*sAAAABC + 2*sAAAABB*sAAABBB + 2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 4*sAAAABB*sAAABBC + 2*sAAAABC*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAAABB*sAAABCD + 2*sAAABBC*sAAABCD + 4*sAAAAAA*sAABBCC + 4*sAAAAAB*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAAAABC*sAABBCC + 4*sAAABBB*sAABBCC + 6*sAAABBC*sAABBCC + 4*sAAABCD*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABBCD + 2*sAAAAAB*sAABBCD + 4*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAAABBB*sAABBCD + 4*sAAABBC*sAABBCD + 2*sAAABCD*sAABBCD + 6*sAABBCC*sAABBCD + 2*sAABBCD**2 + 2*sAAAABB*sAABCDE + 2*sAAABBC*sAABCDE + 4*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAABBC*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + sABCDEF**2)
    alpha_3 = (1-d)*(2*sAAABBB + sAAABBC + sAAABCD) + d*(2*sAAAAAA*sAAABBB + 2*sAAAAAB*sAAABBB + 2*sAAAABB*sAAABBB + 2*sAAAABC*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAAABCD + 2*sAAABBB*sAABBCC + 2*sAAABBB*sAABBCD + 2*sAAABBB*sAABCDE + 2*sAAABBB*sABCDEF + 2*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4 = (1-d)*(sAAAABB + sAAAABC) + d*(sAABBCC**2 + 2*sAABBCC*sAABBCD + sAABBCD**2 + 2*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + sAABCDE**2 + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF + 2*sAAABCD*sABCDEF)
    alpha_5 = (1-d)*(sAAAAAB) + d*(2*sAAABBB*sAABBCC + 2*sAAABBC*sAABBCC + 2*sAAABCD*sAABBCC + 2*sAAABBB*sAABBCD + 2*sAAABBC*sAABBCD + 2*sAAABCD*sAABBCD + 2*sAAABBB*sAABCDE + 2*sAAABBC*sAABCDE + 2*sAAABCD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAAABC*sABCDEF)
    alpha_6 = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2 + 2*sAAABBB*sAAABCD + 2*sAAABBC*sAAABCD + sAAABCD**2 + 2*sAAAABB*sAABBCC + 2*sAAAABC*sAABBCC + 2*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAAAABB*sAABCDE + 2*sAAAABC*sAABCDE + 2*sAAAAAB*sABCDEF)
    alpha_7 = d*(2*sAAAABB*sAAABBB + 2*sAAAABC*sAAABBB + 2*sAAAABB*sAAABBC + 2*sAAAABC*sAAABBC + 2*sAAAABB*sAAABCD + 2*sAAAABC*sAAABCD + 2*sAAAAAB*sAABBCC + 2*sAAAAAB*sAABBCD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
    alpha_8 = d*(sAAAABB**2 + 2*sAAAABB*sAAAABC + sAAAABC**2 + 2*sAAAAAB*sAAABBB + 2*sAAAAAB*sAAABBC + 2*sAAAAAB*sAAABCD + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCDE)
    alpha_9 = d*(2*sAAAAAB*sAAAABB + 2*sAAAAAB*sAAAABC + 2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC + 2*sAAAAAA*sAAABCD)
    alpha_10 = d*(sAAAAAB**2 + 2*sAAAAAA*sAAAABB + 2*sAAAAAA*sAAAABC)
    alpha_11 = d*(2*sAAAAAA*sAAAAAB)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, full model, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaaabb,raaabbb,raaaabc,raaabbc,raabbcc,raaabcd,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_0_unique = function(raaaaab, raaaabb, raaabbb, raaaabc, raaabbc, raabbcc, raaabcd, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa = 1-raaaaab-raaaabb-raaabbb-raaaabc-raaabbc-raabbcc-raaabcd-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    tAAAAAA = (raaaaaa)**k
    tAAAAAB = (raaaaaa+raaaaab)**k
    tAAAABB = (raaaaaa+raaaabb)**k
    tAAABBB = (raaaaaa+raaabbb)**k
    tAAAABC = (raaaaaa+raaaaab+raaaabb+raaaabc)**k
    tAAABBC = (raaaaaa+raaaaab+raaabbb+raaabbc)**k
    tAABBCC = (raaaaaa+raaaabb+raabbcc)**k
    tAAABCD = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raaabcd)**k
    tAABBCD = (raaaaaa+raaaaab+raaaabb+raaaabc+raabbcc+raabbcd)**k
    tAABCDE = (raaaaaa+raaaaab+raaaabb+raaabbb+raaaabc+raaabbc+raabbcc+raaabcd+raabbcd+raabcde)**k
    sAAAAAA = tAAAAAA
    sAAAAAB = tAAAAAB-tAAAAAA
    sAAAABB = tAAAABB-tAAAAAA
    sAAABBB = tAAABBB-tAAAAAA
    sAAAABC = tAAAABC-tAAAABB-tAAAAAB+tAAAAAA
    sAAABBC = tAAABBC-tAAABBB-tAAAAAB+tAAAAAA
    sAABBCC = tAABBCC-tAAAABB
    sAAABCD = tAAABCD-tAAABBC-tAAAABC+tAAAAAB
    sAABBCD = tAABBCD-tAABBCC-tAAAABC+tAAAABB
    sAABCDE = tAABCDE-tAABBCD-tAAABCD+tAAAABC
    sABCDEF = 1-tAABCDE
    alpha_1 = (1-d)*(sAAAAAB + 2*sAAAABC + sAAABBC + 3*sAAABCD + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2 = (1-d)*(sAAAABB + sAAABBC + 3*sAABBCC + 2*sAABBCD + sAABCDE)
    alpha_3 = (1-d)*(2*sAAABBB + sAAABBC + sAAABCD)
    alpha_4 = (1-d)*(sAAAABB + sAAAABC)
    alpha_5 = (1-d)*(sAAAAAB)
    alpha_6 = (1-d)*(sAAAAAA)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=1) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_1 = function(raaaaab, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaaabc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taaabcd  = (raaaaaa+raaaaab+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAAAAB + 2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 4*sAAAAAA*sAAAABC + 6*sAAAAAB*sAAAABC + 4*sAAAABC**2 + 6*sAAAAAA*sAAABCD + 8*sAAAAAB*sAAABCD + 10*sAAAABC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 12*sAAAABC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 14*sAAAABC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAABCDE) + d*(sABCDEF**2)
    alpha_3  = (1-d)*(sAAABCD) + d*(2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABC) + d*(sAABCDE**2 + 2*sAAABCD*sABCDEF)
    alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABCD*sAABCDE + 2*sAAAABC*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCD**2 + 2*sAAAABC*sAABCDE + 2*sAAAAAB*sABCDEF)
    alpha_7  = d*(2*sAAAABC*sAAABCD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABC**2 + 2*sAAAAAB*sAAABCD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAB*sAAAABC + 2*sAAAAAA*sAAABCD)
    alpha_10 = d*(sAAAAAB**2 + 2*sAAAAAA*sAAAABC)
    alpha_11 = d*(2*sAAAAAA*sAAAAAB)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=1, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_1_unique = function(raaaaab, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaaabc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taaabcd  = (raaaaaa+raaaaab+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAAAAB + 2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAABCDE)
    alpha_3_unique  = (1-d)*(sAAABCD)
    alpha_4_unique  = (1-d)*(sAAAABC)
    alpha_5_unique  = (1-d)*(sAAAAAB)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=2) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aabbcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_2 = function(raaaaab, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaaabc-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taabbcd  = (raaaaaa+raaaaab+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAAAAB + 2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 4*sAAAAAA*sAAAABC + 6*sAAAAAB*sAAAABC + 4*sAAAABC**2 + 4*sAAAAAA*sAABBCD + 6*sAAAAAB*sAABBCD + 8*sAAAABC*sAABBCD + 4*sAABBCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 12*sAAAABC*sAABCDE + 12*sAABBCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 14*sAAAABC*sABCDEF + 14*sAABBCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(2*sAABBCD + sAABCDE) + d*(2*sAAAAAA*sAABBCD + 2*sAAAAAB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAABBCD**2 + 2*sAABBCD*sAABCDE + 2*sAABBCD*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(0) + d*(2*sAABBCD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABC) + d*(sAABBCD**2 + 2*sAABBCD*sAABCDE + sAABCDE**2)
    alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAAABC*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABC*sAABBCD + 2*sAAAABC*sAABCDE + 2*sAAAAAB*sABCDEF)
    alpha_7  = d*(2*sAAAAAB*sAABBCD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABC**2 + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAB*sAAAABC)
    alpha_10 = d*(sAAAAAB**2 + 2*sAAAAAA*sAAAABC)
    alpha_11 = d*(2*sAAAAAA*sAAAAAB)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=2, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaaabc, aabbcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_2_unique = function(raaaaab, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaaabc-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaaabc  = (raaaaaa+raaaaab+raaaabc)**k
    taabbcd  = (raaaaaa+raaaaab+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaaab+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAAABC  = taaaabc-taaaaab
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAAAAB + 2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(2*sAABBCD + sAABCDE)
    alpha_3_unique  = (1-d)*(0)
    alpha_4_unique  = (1-d)*(sAAAABC)
    alpha_5_unique  = (1-d)*(sAAAAAB)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=3) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_3 = function(raaaaab, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaabbc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taaabcd  = (raaaaaa+raaaaab+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAAAAB + sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 6*sAAAAAA*sAAABCD + 8*sAAAAAB*sAAABCD + 8*sAAABBC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 10*sAAABBC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 12*sAAABBC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAABBC + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 2*sAAABBC*sAAABCD + 2*sAAABBC*sAABCDE + 2*sAAABBC*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(sAAABBC + sAAABCD) + d*(2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(0) + d*(sAABCDE**2 + 2*sAAABBC*sABCDEF + 2*sAAABCD*sABCDEF)
    alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABBC*sAABCDE + 2*sAAABCD*sAABCDE)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBC**2 + 2*sAAABBC*sAAABCD + sAAABCD**2 + 2*sAAAAAB*sABCDEF)
    alpha_7  = d*(2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(2*sAAAAAB*sAAABBC + 2*sAAAAAB*sAAABCD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAA*sAAABBC + 2*sAAAAAA*sAAABCD)
    alpha_10 = d*(sAAAAAB**2)
    alpha_11 = d*(2*sAAAAAA*sAAAAAB)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=3, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_3_unique = function(raaaaab, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaabbc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taaabcd  = (raaaaaa+raaaaab+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAAAAB + sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAABBC + sAABCDE)
    alpha_3_unique  = (1-d)*(sAAABBC + sAAABCD)
    alpha_4_unique  = (1-d)*(0)
    alpha_5_unique  = (1-d)*(sAAAAAB)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=4) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_4 = function(raaaaab, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaabbc-raabccd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 6*sAAAAAB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 8*sAAAAAA*sAABCDE + 10*sAAAAAB*sAABCDE + 10*sAAABBC*sAABCDE + 12*sAABCCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAAAAB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAABBC*sAABCDE + 2*sAABCCD*sAABCDE + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(sAAABBC) + d*(2*sAABCCD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAABCCD*sAABCDE + sAABCDE**2 + 2*sAAABBC*sABCDEF)
    alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABBC*sAABCCD + 2*sAAABBC*sAABCDE)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBC**2 + 2*sAAAAAB*sABCDEF)
    alpha_7  = d*(2*sAAAAAB*sAABCCD + 2*sAAAAAB*sAABCDE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(2*sAAAAAB*sAAABBC + 2*sAAAAAA*sAABCCD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAA*sAAABBC)
    alpha_10 = d*(sAAAAAB**2)
    alpha_11 = d*(2*sAAAAAA*sAAAAAB)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=4, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_4_unique = function(raaaaab, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaabbc-raabccd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaaaab+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE)
    alpha_3_unique  = (1-d)*(sAAABBC)
    alpha_4_unique  = (1-d)*(0)
    alpha_5_unique  = (1-d)*(sAAAAAB)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=5) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, abcdde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_5 = function(raaaaab, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaabbc-raabccd-rabcdde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaaaab+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
    alpha_1  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAAAAB + 2*sAAAAAB**2 + 2*sAAAAAA*sAAABBC + 4*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 6*sAAAAAB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 6*sAAAAAA*sABCDDE + 8*sAAAAAB*sABCDDE + 8*sAAABBC*sABCDDE + 10*sAABCCD*sABCDDE + 6*sABCDDE**2 + 10*sAAAAAA*sABCDEF + 12*sAAAAAB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 16*sABCDDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAAAAB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAAAAB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAAAAA*sABCDDE + 2*sAAAAAB*sABCDDE + 4*sAAABBC*sABCDDE + 4*sAABCCD*sABCDDE + 3*sABCDDE**2 + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + 4*sABCDDE*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(sAAABBC) + d*(2*sAABCCD*sABCDDE + 2*sAABCCD*sABCDEF)
    alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAAABBC*sABCDDE + 2*sAAABBC*sABCDEF)
    alpha_5  = (1-d)*(sAAAAAB) + d*(2*sAAABBC*sAABCCD)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBC**2 + 2*sAAAAAB*sABCDDE + 2*sAAAAAB*sABCDEF)
    alpha_7  = d*(2*sAAAAAB*sAABCCD + 2*sAAAAAA*sABCDDE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(2*sAAAAAB*sAAABBC + 2*sAAAAAA*sAABCCD)
    alpha_9  = d*(2*sAAAAAA*sAAABBC)
    alpha_10 = d*(sAAAAAB**2)
    alpha_11 = d*(2*sAAAAAA*sAAAAAB)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=5, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaaab,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaaab, aaabbc, aabccd, abcdde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_5_unique = function(raaaaab, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaaab-raaabbc-raabccd-rabcdde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaaab  = (raaaaaa+raaaaab)**k
    taaabbc  = (raaaaaa+raaaaab+raaabbc)**k
    taabccd  = (raaaaaa+raaaaab+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaaaab+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAAAAB  = taaaaab-taaaaaa
    sAAABBC  = taaabbc-taaaaab
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
    alpha_1_unique  = (1-d)*(sAAAAAB + sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE)
    alpha_3_unique  = (1-d)*(sAAABBC)
    alpha_4_unique  = (1-d)*(0)
    alpha_5_unique  = (1-d)*(sAAAAAB)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=6) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_6 = function(raaaabb, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaaabc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taaabcd  = (raaaaaa+raaaabb+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAAAABC + 4*sAAAABB*sAAAABC + 4*sAAAABC**2 + 6*sAAAAAA*sAAABCD + 6*sAAAABB*sAAABCD + 10*sAAAABC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 12*sAAAABC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 14*sAAAABC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAABB*sAAAABC + 2*sAAAABB*sAAABCD + 2*sAAAABB*sAABCDE + 2*sAAAABB*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(sAAABCD) + d*(2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB + sAAAABC) + d*(sAABCDE**2 + 2*sAAABCD*sABCDEF)
    alpha_5  = (1-d)*(0) + d*(2*sAAABCD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAAABC*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCD**2 + 2*sAAAABB*sAABCDE + 2*sAAAABC*sAABCDE)
    alpha_7  = d*(2*sAAAABB*sAAABCD + 2*sAAAABC*sAAABCD + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAABB*sAAAABC + sAAAABC**2 + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAA*sAAABCD)
    alpha_10 = d*(2*sAAAAAA*sAAAABB + 2*sAAAAAA*sAAAABC)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=6, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaaabc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_6_unique = function(raaaabb, raaaabc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaaabc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taaabcd  = (raaaaaa+raaaabb+raaaabc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAAABCD  = taaabcd-taaaabc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(2*sAAAABC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + sAABCDE)
    alpha_3_unique  = (1-d)*(sAAABCD)
    alpha_4_unique  = (1-d)*(sAAAABB + sAAAABC)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=7) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aabbcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_7 = function(raaaabb, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaaabc-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taabbcd  = (raaaaaa+raaaabb+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAAAABC + 4*sAAAABB*sAAAABC + 4*sAAAABC**2 + 4*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 8*sAAAABC*sAABBCD + 4*sAABBCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 12*sAAAABC*sAABCDE + 12*sAABBCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 14*sAAAABC*sABCDEF + 14*sAABBCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + 2*sAABBCD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAABB*sAAAABC + 2*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAABBCD**2 + 2*sAAAABB*sAABCDE + 2*sAABBCD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAABBCD*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(0) + d*(2*sAABBCD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB + sAAAABC) + d*(sAABBCD**2 + 2*sAABBCD*sAABCDE + sAABCDE**2)
    alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEF + 2*sAAAABC*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCD + 2*sAAAABC*sAABBCD + 2*sAAAABB*sAABCDE + 2*sAAAABC*sAABCDE)
    alpha_7  = d*(2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAABB*sAAAABC + sAAAABC**2 + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(0)
    alpha_10 = d*(2*sAAAAAA*sAAAABB + 2*sAAAAAA*sAAAABC)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=7, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaaabc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaaabc, aabbcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_7_unique = function(raaaabb, raaaabc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaaabc-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaaabc  = (raaaaaa+raaaabb+raaaabc)**k
    taabbcd  = (raaaaaa+raaaabb+raaaabc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raaaabc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAAABC  = taaaabc-taaaabb
    sAABBCD  = taabbcd-taaaabc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(2*sAAAABC + 2*sAABBCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + 2*sAABBCD + sAABCDE)
    alpha_3_unique  = (1-d)*(0)
    alpha_4_unique  = (1-d)*(sAAAABB + sAAAABC)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=8) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaabcc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_8 = function(raaaabb, raaabcc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaabcc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taaabcd  = (raaaaaa+raaaabb+raaabcc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAAABCD  = taaabcd-taaabcc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAABCC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABCC + 2*sAAAABB*sAAABCC + 2*sAAABCC**2 + 6*sAAAAAA*sAAABCD + 6*sAAAABB*sAAABCD + 8*sAAABCC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 10*sAAABCC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 12*sAAABCC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + sAAABCC + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAAAA*sAAABCC + 4*sAAAABB*sAAABCC + 2*sAAABCC**2 + 2*sAAAABB*sAAABCD + 2*sAAABCC*sAAABCD + 2*sAAAABB*sAABCDE + 2*sAAABCC*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAABCC*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(sAAABCC + sAAABCD) + d*(2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB) + d*(sAABCDE**2 + 2*sAAABCC*sABCDEF + 2*sAAABCD*sABCDEF)
    alpha_5  = (1-d)*(0) + d*(2*sAAABCC*sAABCDE + 2*sAAABCD*sAABCDE + 2*sAAAABB*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCC**2 + 2*sAAABCC*sAAABCD + sAAABCD**2 + 2*sAAAABB*sAABCDE)
    alpha_7  = d*(2*sAAAABB*sAAABCC + 2*sAAAABB*sAAABCD + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAA*sAAABCC + 2*sAAAAAA*sAAABCD)
    alpha_10 = d*(2*sAAAAAA*sAAAABB)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=8, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaabcc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_8_unique = function(raaaabb, raaabcc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaabcc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taaabcd  = (raaaaaa+raaaabb+raaabcc+raaabcd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAAABCD  = taaabcd-taaabcc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAABCC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + sAAABCC + sAABCDE)
    alpha_3_unique  = (1-d)*(sAAABCC + sAAABCD)
    alpha_4_unique  = (1-d)*(sAAAABB)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=9) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaabcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_9 = function(raaaabb, raaabcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaabcc-raabcdd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABCC + 2*sAAAABB*sAAABCC + 2*sAAABCC**2 + 4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAAABCC*sAABCDD + 4*sAABCDD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 10*sAAABCC*sAABCDE + 12*sAABCDD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 12*sAAABCC*sABCDEF + 14*sAABCDD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAAAA*sAAABCC + 4*sAAAABB*sAAABCC + 2*sAAABCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAAABCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAABB*sAABCDE + 2*sAAABCC*sAABCDE + 2*sAABCDD*sAABCDE + 2*sAAAABB*sABCDEF + 2*sAAABCC*sABCDEF + 2*sAABCDD*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(sAAABCC) + d*(2*sAABCDD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB) + d*(sAABCDD**2 + 2*sAABCDD*sAABCDE + sAABCDE**2 + 2*sAAABCC*sABCDEF)
    alpha_5  = (1-d)*(0) + d*(2*sAAABCC*sAABCDD + 2*sAAABCC*sAABCDE + 2*sAAAABB*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCC**2 + 2*sAAAABB*sAABCDD + 2*sAAAABB*sAABCDE)
    alpha_7  = d*(2*sAAAABB*sAAABCC + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABCDD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAA*sAAABCC)
    alpha_10 = d*(2*sAAAAAA*sAAAABB)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=9, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaabcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_9_unique = function(raaaabb, raaabcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaabcc-raabcdd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raaabcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sAABCDE)
    alpha_3_unique  = (1-d)*(sAAABCC)
    alpha_4_unique  = (1-d)*(sAAAABB)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=10) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaabcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, abcdee, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_10 = function(raaaabb, raaabcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaabcc-raabcdd-rabcdee-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raaabcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
    alpha_1  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sABCDEE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABCC + 2*sAAAABB*sAAABCC + 2*sAAABCC**2 + 4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAAABCC*sAABCDD + 4*sAABCDD**2 + 6*sAAAAAA*sABCDEE + 6*sAAAABB*sABCDEE + 8*sAAABCC*sABCDEE + 10*sAABCDD*sABCDEE + 6*sABCDEE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 12*sAAABCC*sABCDEF + 14*sAABCDD*sABCDEF + 16*sABCDEE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sABCDEE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 2*sAAAAAA*sAAABCC + 4*sAAAABB*sAAABCC + 2*sAAABCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAAABCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAAAA*sABCDEE + 4*sAAAABB*sABCDEE + 4*sAAABCC*sABCDEE + 4*sAABCDD*sABCDEE + 3*sABCDEE**2 + 2*sAAAABB*sABCDEF + 2*sAAABCC*sABCDEF + 2*sAABCDD*sABCDEF + 4*sABCDEE*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(sAAABCC) + d*(2*sAABCDD*sABCDEE + 2*sAABCDD*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB) + d*(sAABCDD**2 + 2*sAAABCC*sABCDEE + 2*sAAABCC*sABCDEF)
    alpha_5  = (1-d)*(0) + d*(2*sAAABCC*sAABCDD + 2*sAAAABB*sABCDEE + 2*sAAAABB*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABCC**2 + 2*sAAAABB*sAABCDD)
    alpha_7  = d*(2*sAAAABB*sAAABCC + 2*sAAAAAA*sABCDEE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABCDD)
    alpha_9  = d*(2*sAAAAAA*sAAABCC)
    alpha_10 = d*(2*sAAAAAA*sAAAABB)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=10, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raaabcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aaabcc, aabcdd, abcdee, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_10_unique = function(raaaabb, raaabcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raaabcc-raabcdd-rabcdee-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taaabcc  = (raaaaaa+raaaabb+raaabcc)**k
    taabcdd  = (raaaaaa+raaaabb+raaabcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raaabcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAAABCC  = taaabcc-taaaabb
    sAABCDD  = taabcdd-taaabcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
    alpha_1_unique  = (1-d)*(sAAABCC + 2*sAABCDD + 4*sABCDEE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + sAAABCC + 2*sAABCDD + sABCDEE)
    alpha_3_unique  = (1-d)*(sAAABCC)
    alpha_4_unique  = (1-d)*(sAAAABB)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=11) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raabbcc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabbcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_11 = function(raaaabb, raabbcc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raabbcc-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabbcd  = (raaaaaa+raaaabb+raabbcc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABBCD  = taabbcd-taabbcc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(2*sAABBCD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 4*sAABBCC*sAABBCD + 4*sAABBCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 8*sAABBCC*sAABCDE + 12*sAABBCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABBCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABBCD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 4*sAAAAAA*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABBCD + 4*sAAAABB*sAABBCD + 6*sAABBCC*sAABBCD + 2*sAABBCD**2 + 2*sAAAABB*sAABCDE + 4*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + 2*sAAAABB*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(0) + d*(2*sAABBCC*sABCDEF + 2*sAABBCD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB) + d*(sAABBCC**2 + 2*sAABBCC*sAABBCD + sAABBCD**2 + 2*sAABBCC*sAABCDE + 2*sAABBCD*sAABCDE + sAABCDE**2)
    alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCC + 2*sAAAABB*sAABBCD + 2*sAAAABB*sAABCDE)
    alpha_7  = d*(2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABBCD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(0)
    alpha_10 = d*(2*sAAAAAA*sAAAABB)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=11, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raabbcc,raabbcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabbcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_11_unique = function(raaaabb, raabbcc, raabbcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raabbcc-raabbcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabbcd  = (raaaaaa+raaaabb+raabbcc+raabbcd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabbcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABBCD  = taabbcd-taabbcc
    sAABCDE  = taabcde-taabbcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(2*sAABBCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABBCD + sAABCDE)
    alpha_3_unique  = (1-d)*(0)
    alpha_4_unique  = (1-d)*(sAAAABB)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=12) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raabbcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_12 = function(raaaabb, raabbcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raabbcc-raabcdd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(2*sAABCDD + 4*sAABCDE + 6*sABCDEF) + d*(4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAABBCC*sAABCDD + 4*sAABCDD**2 + 8*sAAAAAA*sAABCDE + 8*sAAAABB*sAABCDE + 8*sAABBCC*sAABCDE + 12*sAABCDD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABCDD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sAABCDE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 4*sAAAAAA*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAABBCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAABB*sAABCDE + 4*sAABBCC*sAABCDE + 2*sAABCDD*sAABCDE + 2*sAAAABB*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(0) + d*(2*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB) + d*(sAABBCC**2 + 2*sAABBCC*sAABCDD + sAABCDD**2 + 2*sAABBCC*sAABCDE + 2*sAABCDD*sAABCDE + sAABCDE**2)
    alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCC + 2*sAAAABB*sAABCDD + 2*sAAAABB*sAABCDE)
    alpha_7  = d*(2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABCDD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(0)
    alpha_10 = d*(2*sAAAAAA*sAAAABB)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=12, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raabbcc,raabcdd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_12_unique = function(raaaabb, raabbcc, raabcdd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raabbcc-raabcdd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    taabcde  = (raaaaaa+raaaabb+raabbcc+raabcdd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sAABCDE  = taabcde-taabcdd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(2*sAABCDD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sAABCDE)
    alpha_3_unique  = (1-d)*(0)
    alpha_4_unique  = (1-d)*(sAAAABB)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=13) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raabbcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, abcdee, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_13 = function(raaaabb, raabbcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raabbcc-raabcdd-rabcdee-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raabbcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
    alpha_1  = (1-d)*(2*sAABCDD + 4*sABCDEE + 6*sABCDEF) + d*(4*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 4*sAABBCC*sAABCDD + 4*sAABCDD**2 + 6*sAAAAAA*sABCDEE + 6*sAAAABB*sABCDEE + 6*sAABBCC*sABCDEE + 10*sAABCDD*sABCDEE + 6*sABCDEE**2 + 10*sAAAAAA*sABCDEF + 10*sAAAABB*sABCDEF + 10*sAABBCC*sABCDEF + 14*sAABCDD*sABCDEF + 16*sABCDEE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sABCDEE) + d*(2*sAAAAAA*sAAAABB + 2*sAAAABB**2 + 4*sAAAAAA*sAABBCC + 6*sAAAABB*sAABBCC + 4*sAABBCC**2 + 2*sAAAAAA*sAABCDD + 4*sAAAABB*sAABCDD + 6*sAABBCC*sAABCDD + 2*sAABCDD**2 + 2*sAAAAAA*sABCDEE + 4*sAAAABB*sABCDEE + 6*sAABBCC*sABCDEE + 4*sAABCDD*sABCDEE + 3*sABCDEE**2 + 2*sAAAABB*sABCDEF + 4*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF + 4*sABCDEE*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(0) + d*(2*sAABBCC*sABCDEE + 2*sAABCDD*sABCDEE + 2*sAABBCC*sABCDEF + 2*sAABCDD*sABCDEF)
    alpha_4  = (1-d)*(sAAAABB) + d*(sAABBCC**2 + 2*sAABBCC*sAABCDD + sAABCDD**2)
    alpha_5  = (1-d)*(0) + d*(2*sAAAABB*sABCDEE + 2*sAAAABB*sABCDEF)
    alpha_6  = (1-d)*(sAAAAAA) + d*(2*sAAAABB*sAABBCC + 2*sAAAABB*sAABCDD)
    alpha_7  = d*(2*sAAAAAA*sABCDEE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(sAAAABB**2 + 2*sAAAAAA*sAABBCC + 2*sAAAAAA*sAABCDD)
    alpha_9  = d*(0)
    alpha_10 = d*(2*sAAAAAA*sAAAABB)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=13, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaaabb,raabbcc,raabcdd,rabcdee,rabcdef Numerics corresponding to the nucleotide heterozygosities aaaabb, aabbcc, aabcdd, abcdee, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_13_unique = function(raaaabb, raabbcc, raabcdd, rabcdee, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaaabb-raabbcc-raabcdd-rabcdee-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaaabb  = (raaaaaa+raaaabb)**k
    taabbcc  = (raaaaaa+raaaabb+raabbcc)**k
    taabcdd  = (raaaaaa+raaaabb+raabbcc+raabcdd)**k
    tabcdee  = (raaaaaa+raaaabb+raabbcc+raabcdd+rabcdee)**k
    sAAAAAA  = taaaaaa
    sAAAABB  = taaaabb-taaaaaa
    sAABBCC  = taabbcc-taaaabb
    sAABCDD  = taabcdd-taabbcc
    sABCDEE  = tabcdee-taabcdd
    sABCDEF  = 1-tabcdee
    alpha_1_unique  = (1-d)*(2*sAABCDD + 4*sABCDEE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAAABB + 3*sAABBCC + 2*sAABCDD + sABCDEE)
    alpha_3_unique  = (1-d)*(0)
    alpha_4_unique  = (1-d)*(sAAAABB)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=14) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabbb,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_14 = function(raaabbb, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaabbb-raaabbc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taaabcd  = (raaaaaa+raaabbb+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 6*sAAAAAA*sAAABCD + 6*sAAABBB*sAAABCD + 8*sAAABBC*sAAABCD + 6*sAAABCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAABBB*sAABCDE + 10*sAAABBC*sAABCDE + 14*sAAABCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 16*sAAABCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAABBC + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAABBC*sAAABCD + 2*sAAABBC*sAABCDE + 2*sAAABBC*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(2*sAAABBB + sAAABBC + sAAABCD) + d*(2*sAAAAAA*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAAABCD + 2*sAAABBB*sAABCDE + 2*sAAABBB*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(0) + d*(sAABCDE**2 + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF + 2*sAAABCD*sABCDEF)
    alpha_5  = (1-d)*(0) + d*(2*sAAABBB*sAABCDE + 2*sAAABBC*sAABCDE + 2*sAAABCD*sAABCDE)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2 + 2*sAAABBB*sAAABCD + 2*sAAABBC*sAAABCD + sAAABCD**2)
    alpha_7  = d*(2*sAAAAAA*sABCDEF)
    alpha_8  = d*(2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC + 2*sAAAAAA*sAAABCD)
    alpha_10 = d*(0)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=14, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabbb,raaabbc,raaabcd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aaabcd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_14_unique = function(raaabbb, raaabbc, raaabcd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaabbb-raaabbc-raaabcd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taaabcd  = (raaaaaa+raaabbb+raaabbc+raaabcd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raaabcd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAAABCD  = taaabcd-taaabbc
    sAABCDE  = taabcde-taaabcd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAABBC + 3*sAAABCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAABBC + sAABCDE)
    alpha_3_unique  = (1-d)*(2*sAAABBB + sAAABBC + sAAABCD)
    alpha_4_unique  = (1-d)*(0)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=15) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabbb,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_15 = function(raaabbb, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaabbb-raaabbc-raabccd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
    alpha_1  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 4*sAAABBB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 8*sAAAAAA*sAABCDE + 8*sAAABBB*sAABCDE + 10*sAAABBC*sAABCDE + 12*sAABCCD*sAABCDE + 8*sAABCDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 18*sAABCDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAABBB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAABBC*sAABCDE + 2*sAABCCD*sAABCDE + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(2*sAAABBB + sAAABBC) + d*(2*sAAAAAA*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAABCCD + 2*sAAABBB*sAABCDE + 2*sAAABBB*sABCDEF + 2*sAABCCD*sABCDEF + 2*sAABCDE*sABCDEF)
    alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAABCCD*sAABCDE + sAABCDE**2 + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF)
    alpha_5  = (1-d)*(0) + d*(2*sAAABBB*sAABCCD + 2*sAAABBC*sAABCCD + 2*sAAABBB*sAABCDE + 2*sAAABBC*sAABCDE)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2)
    alpha_7  = d*(2*sAAAAAA*sABCDEF)
    alpha_8  = d*(2*sAAAAAA*sAABCCD + 2*sAAAAAA*sAABCDE)
    alpha_9  = d*(2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC)
    alpha_10 = d*(0)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=15, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabbb,raaabbc,raabccd,raabcde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, aabcde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_15_unique = function(raaabbb, raaabbc, raabccd, raabcde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaabbb-raaabbc-raabccd-raabcde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    taabcde  = (raaaaaa+raaabbb+raaabbc+raabccd+raabcde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sAABCDE  = taabcde-taabccd
    sABCDEF  = 1-taabcde
    alpha_1_unique  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sAABCDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sAABCDE)
    alpha_3_unique  = (1-d)*(2*sAAABBB + sAAABBC)
    alpha_4_unique  = (1-d)*(0)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

  #' Produce model estimated (p=6, topology=16) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabbb,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, abcdde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_16 = function(raaabbb, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaabbb-raaabbc-raabccd-rabcdde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaabbb+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
    alpha_1  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 4*sAAAAAA*sAABCCD + 4*sAAABBB*sAABCCD + 6*sAAABBC*sAABCCD + 4*sAABCCD**2 + 6*sAAAAAA*sABCDDE + 6*sAAABBB*sABCDDE + 8*sAAABBC*sABCDDE + 10*sAABCCD*sABCDDE + 6*sABCDDE**2 + 10*sAAAAAA*sABCDEF + 10*sAAABBB*sABCDEF + 12*sAAABBC*sABCDEF + 14*sAABCCD*sABCDEF + 16*sABCDDE*sABCDEF + 10*sABCDEF**2)
    alpha_2  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE) + d*(2*sAAAAAA*sAAABBC + 2*sAAABBB*sAAABBC + 2*sAAABBC**2 + 2*sAAAAAA*sAABCCD + 2*sAAABBB*sAABCCD + 4*sAAABBC*sAABCCD + 2*sAABCCD**2 + 2*sAAAAAA*sABCDDE + 2*sAAABBB*sABCDDE + 4*sAAABBC*sABCDDE + 4*sAABCCD*sABCDDE + 3*sABCDDE**2 + 2*sAAABBC*sABCDEF + 2*sAABCCD*sABCDEF + 4*sABCDDE*sABCDEF + sABCDEF**2)
    alpha_3  = (1-d)*(2*sAAABBB + sAAABBC) + d*(2*sAAAAAA*sAAABBB + 2*sAAABBB**2 + 2*sAAABBB*sAAABBC + 2*sAAABBB*sAABCCD + 2*sAAABBB*sABCDDE + 2*sAABCCD*sABCDDE + 2*sAAABBB*sABCDEF + 2*sAABCCD*sABCDEF)
    alpha_4  = (1-d)*(0) + d*(sAABCCD**2 + 2*sAAABBB*sABCDDE + 2*sAAABBC*sABCDDE + 2*sAAABBB*sABCDEF + 2*sAAABBC*sABCDEF)
    alpha_5  = (1-d)*(0) + d*(2*sAAABBB*sAABCCD + 2*sAAABBC*sAABCCD)
    alpha_6  = (1-d)*(sAAAAAA) + d*(sAAABBB**2 + 2*sAAABBB*sAAABBC + sAAABBC**2)
    alpha_7  = d*(2*sAAAAAA*sABCDDE + 2*sAAAAAA*sABCDEF)
    alpha_8  = d*(2*sAAAAAA*sAABCCD)
    alpha_9  = d*(2*sAAAAAA*sAAABBB + 2*sAAAAAA*sAAABBC)
    alpha_10 = d*(0)
    alpha_11 = d*(0)
    alpha_12 = d*(sAAAAAA**2)
    alpha_1  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)+
    alpha_7  * dnbinom(x, size = kmercov*7 / bias, mu = kmercov*7)+
    alpha_8  * dnbinom(x, size = kmercov*8 / bias, mu = kmercov*8)+
    alpha_9  * dnbinom(x, size = kmercov*9 / bias, mu = kmercov*9)+
    alpha_10 * dnbinom(x, size = kmercov*10 / bias, mu = kmercov*10)+
    alpha_11 * dnbinom(x, size = kmercov*11 / bias, mu = kmercov*11)+
    alpha_12 * dnbinom(x, size = kmercov*12 / bias, mu = kmercov*12)
  }

  #' Produce model estimated (p=6, topology=16, unique portion) y-coordinates of the kmer spectra given the kmer size, repetitiveness, average polyploid kmer coverage, bias, and x-coordinates of the kmer spectra.
  #'
  #' @param raaabbb,raaabbc,raabccd,rabcdde,rabcdef Numerics corresponding to the nucleotide heterozygosities aaabbb, aaabbc, aabccd, abcdde, and abcdef respectively.
  #' @param k An integer corresponding to the kmer length.
  #' @param d A numeric corresponding to the repetitiveness.
  #' @param kmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  #' @param bias A numeric corresponding to the overdispersion of the negative binomial distribution.
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @return A numeric vector of the model estimated y-coordinates of the kmer spectra.
  #' @export
  predict6_16_unique = function(raaabbb, raaabbc, raabccd, rabcdde, rabcdef, k, d, kmercov, bias, x)
  {
    raaaaaa  = 1-raaabbb-raaabbc-raabccd-rabcdde-rabcdef
    if (raaaaaa < 0 || d > 1) {return(0)}
    taaaaaa  = (raaaaaa)**k
    taaabbb  = (raaaaaa+raaabbb)**k
    taaabbc  = (raaaaaa+raaabbb+raaabbc)**k
    taabccd  = (raaaaaa+raaabbb+raaabbc+raabccd)**k
    tabcdde  = (raaaaaa+raaabbb+raaabbc+raabccd+rabcdde)**k
    sAAAAAA  = taaaaaa
    sAAABBB  = taaabbb-taaaaaa
    sAAABBC  = taaabbc-taaabbb
    sAABCCD  = taabccd-taaabbc
    sABCDDE  = tabcdde-taabccd
    sABCDEF  = 1-tabcdde
    alpha_1_unique  = (1-d)*(sAAABBC + 2*sAABCCD + 4*sABCDDE + 6*sABCDEF)
    alpha_2_unique  = (1-d)*(sAAABBC + 2*sAABCCD + sABCDDE)
    alpha_3_unique  = (1-d)*(2*sAAABBB + sAAABBC)
    alpha_4_unique  = (1-d)*(0)
    alpha_5_unique  = (1-d)*(0)
    alpha_6_unique  = (1-d)*(sAAAAAA)
    alpha_1_unique  * dnbinom(x, size = kmercov*1 / bias, mu = kmercov*1)+
    alpha_2_unique  * dnbinom(x, size = kmercov*2 / bias, mu = kmercov*2)+
    alpha_3_unique  * dnbinom(x, size = kmercov*3 / bias, mu = kmercov*3)+
    alpha_4_unique  * dnbinom(x, size = kmercov*4 / bias, mu = kmercov*4)+
    alpha_5_unique  * dnbinom(x, size = kmercov*5 / bias, mu = kmercov*5)+
    alpha_6_unique  * dnbinom(x, size = kmercov*6 / bias, mu = kmercov*6)
  }

## Format numbers (tiny)
  

  bp_format<-function(num) {paste(formatC(round(num),format="f",big.mark=",", digits=0), "bp",sep=" ")}

  percentage_format<-function(num) {paste(signif(num,6)*100,"%",sep="")}

  X_format<-function(num) {paste(signif(num,4),"X",sep="")}

  na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
  }

## report_results (big)

  #' Report results and make plots
  #'
  #' @param kmer_hist A data frame of the original histogram data (starting at 1 and going up to the max kmer coverage threshold).
  #' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
  #' @param k An integer corresponding to the kmer length.
  #' @param p An integer corresponding to the ploidy.
  #' @param container A list (nls, nlsscore) where nls is the nlsLM model object (with some additional components)
  #' and nlsscore is the score (model RSSE) corresponding to the best fit (of the p forms).
  #' @param foldername A character vector corresponding to the name of the output directory.
  #' @param arguments A data frame of the user-specified inputs.
  #' @param IN_VERBOSE A boolean flag to designate whether report_results is being called in a VERBOSE block.
  #' @export
  
    report_results <- function(kmer_hist,kmer_hist_orig, k, p, container, foldername, arguments, IN_VERBOSE) {
      
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

        #automatically zoom into the relevant regions of the plot, ignore first 15 positions
        xmax=length(x)
        start_orig=which(y_orig == min(y_orig[1:TYPICAL_ERROR]))
        start=which(y == min(y[1:TYPICAL_ERROR]))
        zoomx=x[start:(xmax-1)]
        zoomy_orig=y_orig[start_orig:(xmax-1)]
        zoomy=y[start:(xmax-1)]

        ## allow for a little space above max value past the noise
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
        # x_limit=150

      ## The first part of the plots
        
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

        # thing1/thing2

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

        # (thing1+thing3) / (thing2+thing4)

        ## Plot the distribution, and hopefully with the model fit
          ylabel_orig = "Frequency"
          if (transform_exp == 1) {
            ylabel_transform = "Coverage*Frequency"
          } else {
            ylabel_transform = paste("Coverage^", transform_exp, "*Frequency", sep="")
          }
          png(paste(foldername, "/", arguments$name_prefix, "linear_plot.png", sep=""),
          width=plot_size, height=plot_size, res=resolution)
          par(mar = c(5.1,4.1,6.1,2.1))
          plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
          xlab="Coverage", ylab=ylabel_orig, ylim=c(0,y_limit_orig), xlim=c(0,x_limit_orig),
          cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
          #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
          rect(0, 0, x_limit_orig*1.1 , y_limit_orig*1.1, col=COLOR_BGCOLOR)
          points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
          #  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
          #    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
          #  }
          box(col="black")

          png(paste(foldername, "/", arguments$name_prefix, "transformed_linear_plot.png", sep=""),
          width=plot_size, height=plot_size, res=resolution)
          par(mar = c(5.1,4.1,6.1,2.1))
          plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
          xlab="Coverage", ylab=ylabel_transform, ylim=c(0,y_limit), xlim=c(0,x_limit),
          cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
          #rect(0, 0, max(kmer_hist_orig[[1]])*1.1 , max(kmer_hist_orig[[2]])*1.1, col=COLOR_BGCOLOR)
          rect(0, 0, x_limit*1.1 , y_limit*1.1, col=COLOR_BGCOLOR)
          points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
          #  if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
          #    abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
          #  }
          box(col="black")

          ## Make a second plot in log space over entire range
          png(paste(foldername, "/", arguments$name_prefix, "log_plot.png", sep=""),
          width=plot_size, height=plot_size, res=resolution)
          par(mar = c(5.1,4.1,6.1,2.1))
          plot(kmer_hist_orig, type="n", main="GenomeScope Profile\n\n\n",
          xlab="Coverage", ylab=ylabel_orig, log="xy",
          cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
          rect(1e-10, 1e-10, max(kmer_hist_orig[,1])*10 , max(kmer_hist_orig[,2])*10, col=COLOR_BGCOLOR)
          points(kmer_hist_orig, type="h", col=COLOR_HIST, lwd=2)
          if(length(kmer_hist[,1])!=length(kmer_hist_orig[,1])){
            abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
          }
          box(col="black")

          png(paste(foldername, "/", arguments$name_prefix, "transformed_log_plot.png", sep=""),
          width=plot_size, height=plot_size, res=resolution)
          par(mar = c(5.1,4.1,6.1,2.1))
          plot(kmer_hist_transform, type="n", main="GenomeScope Profile\n\n\n",
          xlab="Coverage", ylab=ylabel_transform, log="xy",
          cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size)
          rect(1e-10, 1e-10, max(kmer_hist_transform[,1])*10 , max(kmer_hist_transform[,2])*10, col=COLOR_BGCOLOR)
          points(kmer_hist_transform, type="h", col=COLOR_HIST, lwd=2)
          if(length(kmer_hist[,1])!=length(kmer_hist_transform[,1])){
            abline(v=length(kmer_hist[,1]),col=COLOR_COVTHRES,lty="dashed", lwd=3)
          }
          box(col="black")

      ## Now the modeling
        if(!is.null(model))
        {
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

        for (i in 1:error_xcutoff_ind)
        {
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

        if (first_zero == -1)
        {
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
        if (p==1)
        {
          unique_hist = amlen*predict1_1_unique(k, amd, akcov, adups, x)
        }
        if (p==2)
        {
          unique_hist = amlen*predict2_1_unique(ahets[[1]], k, amd, akcov, adups, x)
        }
        if (p==3)
        {
          unique_hist = amlen*predict3_1_unique(ahets[[1]], ahets[[2]], k, amd, akcov, adups, x)
        }
        if (p==4)
        {
          if (top==0) {
            unique_hist = amlen*predict4_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)
          } else {
            unique_hist = eval(parse(text = paste("amlen*predict4_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], k, amd, akcov, adups, x)", sep="")))
          }
        }
        if (p==5)
        {
          if (top==0) {
            unique_hist = amlen*predict5_0_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], ahets[[5]], ahets[[6]], k, amd, akcov, adups, x)
          } else {
            unique_hist = eval(parse(text = paste("amlen*predict5_", top, "_unique(ahets[[1]], ahets[[2]], ahets[[3]], ahets[[4]], k, amd, akcov, adups, x)", sep="")))
          }
        }
        if (p==6)
        {
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

        # ## finish plots

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

        #   ## other things

            model_status="done"

            if (!IN_VERBOSE) {
              cat(paste("Model converged het:", format(ahet, digits=3),
              " kcov:", format(akcov, digits=3),
              " err:", format(error_rate[1], digits=3),
              " model fit:", format(adups, digits=3),
              " len:", round(total_len[1]), "\n", sep=""))
            }
            }
            else
            {
              title("\nFailed to converge")
              dev.set(dev.next())
              title("\nFailed to converge")
              cat("Failed to converge.", file=paste(foldername,"/", arguments$name_prefix, "model.txt", sep=""))
              cat("Failed to converge.\n")
            }

            dev.off()
            dev.off()
            dev.off()
            dev.off()

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
        if (p==1)
        {
          cat(paste(sprintf(format_column_1,"Homozygous (a)"),            sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)

        }
        if (p==2)
        {
          cat(paste(sprintf(format_column_1,"Homozygous (aa)"),           sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1,"Heterozygous (ab)"),         sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
          property <- c(property, "Homozygous (aa)", "Heterozygous (ab)")
          min <- c(min, homo[2]*100, hets[[1]][1]*100)
          max <- c(max, homo[1]*100, hets[[1]][2]*100)
        }
        if (p==3)
        {
          cat(paste(sprintf(format_column_1,"Homozygous (aaa)"),          sprintf(format_column_2,percentage_format(homo[2])),         sprintf(format_column_3,percentage_format(homo[1])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1,"Heterozygous (not aaa)"),    sprintf(format_column_2,percentage_format(het[1])),          sprintf(format_column_3,percentage_format(het[2])),  sep=""),                file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1,"aab"),                       sprintf(format_column_2,percentage_format(hets[[1]][1])),         sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1,"abc"),                       sprintf(format_column_2,percentage_format(hets[[2]][1])),         sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""),                file=summaryFile, sep="\n", append=TRUE)
        }
        if (p==4)
        {
          cat(paste(sprintf(format_column_1,"Homozygous (aaaa)"),                    sprintf(format_column_2,percentage_format(homo[2])),      sprintf(format_column_3,percentage_format(homo[1])),      sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1,"Heterozygous (not aaaa)"),              sprintf(format_column_2,percentage_format(het[1])),       sprintf(format_column_3,percentage_format(het[2])),       sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1, switch(top+1, "aaab", "aaab", "aabb")), sprintf(format_column_2,percentage_format(hets[[1]][1])), sprintf(format_column_3,percentage_format(hets[[1]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1, switch(top+1, "aabb", "aabc", "aabc")), sprintf(format_column_2,percentage_format(hets[[2]][1])), sprintf(format_column_3,percentage_format(hets[[2]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste(sprintf(format_column_1, switch(top+1, "aabc", "abcd", "abcd")), sprintf(format_column_2,percentage_format(hets[[3]][1])), sprintf(format_column_3,percentage_format(hets[[3]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
          if (top == 0) {
            cat(paste(sprintf(format_column_1,"abcd"),                               sprintf(format_column_2,percentage_format(hets[[4]][1])), sprintf(format_column_3,percentage_format(hets[[4]][2])), sep=""), file=summaryFile, sep="\n", append=TRUE)
          }
        }
        if (p==5)
        {
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
        if (p==6)
        {
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

        return(output)

        if (VERBOSE)
        {
          cat(paste("\nPercent Kmers Modeled (All Kmers) = ",  percentage_format(model_fit_allscore[1]),    " [", model_fit_allscore[2],    ", ", model_fit_allscore[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste("Percent Kmers Modeled (Full Model) = ",   percentage_format(model_fit_fullscore[1]),   " [", model_fit_fullscore[2],   ", ", model_fit_fullscore[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste("Percent Kmers Modeled (Unique Kmers) = ", percentage_format(model_fit_uniquescore[1]), " [", model_fit_uniquescore[2], ", ", model_fit_uniquescore[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE)

          cat(paste("\nModel RSSE (All Kmers) = ",  model_fit_all[1],    " [", model_fit_all[2],    ", ", model_fit_all[3],    "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste("Model RSSE (Full Model) = ",   model_fit_full[1],   " [", model_fit_full[2],   ", ", model_fit_full[3],   "]", sep=""), file=summaryFile, sep="\n", append=TRUE)
          cat(paste("Model RSSE (Unique Model) = ", model_fit_unique[1], " [", model_fit_unique[2], ", ", model_fit_unique[3], "]", sep=""), file=summaryFile, sep="\n", append=TRUE) 
        }

        ## Finalize the progress
        progressFilename=paste(foldername, "/", arguments$name_prefix, "progress.txt",sep="")
        cat(model_status, file=progressFilename, sep="\n", append=TRUE)

      ## Merfin probabilities

        colors <- c("black","red","green","purple","blue")

        if (p<=2 & FITTED_HIST){

          if (!(p==1)){

          fitted_hist=data.frame(cbind(one_hist,two_hist,thr_hist,fou_hist))
        
          } else {

          fitted_hist=data.frame(cbind(two_hist,fou_hist))
        
          }

          ## Fitted histogram  

            png(paste(foldername, "/fitted_hist.png", sep=""), height = plot_size, width = plot_size, res=resolution)
            layout(matrix(c(1,2), nrow=2, byrow = TRUE),heights=lcm(c(11,5.5)))
            par(mar=c(0,5,1,1))
          
            plot(kmer_hist_orig, type="n", xlab=NULL, ylab=ylabel_orig, ylim=c(0,y_limit_orig), xlim=c(0,akcov*(p*2.5)),
              cex.main=font_size, cex.sub=font_size, cex.axis=font_size, cex.lab=font_size,
               xaxt='n', yaxt='n')
            myTicks = axTicks(4)
            myTicks[1]<-0
            axis(2, at=myTicks, labels=myTicks, cex.axis=font_size,
               cex.lab=font_size)
          
            peaks<-akcov * 1:6

            points(kmer_hist, type="h", col=COLOR_HIST, lwd=2)
            lines(error_kmers, col="black", lwd=1.5)
          
            abline(v=peaks, col="black", lty=2, lwd=0.3)
          
            for (i in 1:ncol(fitted_hist)) {
            lines(fitted_hist[,i]*amlen, type="l", lwd=1.5, col=colors[i+1])
            }

            lines(rowSums(fitted_hist*amlen), col="darkgray", lwd=2, lty=2)
          
            legend_names=c("0-copy", "1-copy", "2-copy", "3-copy", "4-copy")
            legend_lty=c("solid", "solid", "solid", "solid", "solid")
            legend_lwd=c(3,3,3,3,3)
          
            legend("topright",
               ##legend(exp(.65 * log(max(x))), 1.0 * max(y),
               legend=c("Observed",legend_names[1:(ncol(fitted_hist)+1)],"Full model"),
               lty=c("solid",legend_lty[1:(ncol(fitted_hist)+1)], "dashed"),
               lwd=c(2, legend_lwd[1:(ncol(fitted_hist)+1)], 3),
               col=c(COLOR_HIST,colors[1:(ncol(fitted_hist)+1)], "darkgray"),
               bg="white")

            ## Generate lookup_table
          
              lookup_table <- NULL
              plot_table <- NULL
            
              fitted_hist[(akcov*(p*2.5)):length(one_hist),1:ncol(fitted_hist)] <- 0
            
              fitted_hist <- na.zero(fitted_hist)
            
              for (i in 1:(akcov*(p*2.5)-1)){
            
              totalP<-sum(na.zero(error_kmers[i]), rowSums(fitted_hist[i,]*amlen))
           
              prob<-c(na.zero(error_kmers[i]/totalP),as.numeric(fitted_hist[i,]*amlen)/totalP)
            
              plot_table<-rbind(plot_table,prob)
              
              max.p<-max(prob)  
              readK<-which.max(prob)-1
            
              lookup_table<-rbind(lookup_table,c(readK,max.p))
            
              }
            
              lookup_table<-data.frame(lookup_table)
              rownames(plot_table) <- make.names(plot_table[,1], unique = TRUE)
            
              write.table(lookup_table, paste(foldername, "/lookup_table.txt", sep=""), sep=",", row.names = FALSE, col.names = FALSE)

            ## Plot lookup values
          
              par(mar=c(5,5,0,1))
            
              plot_table<-data.frame(plot_table)
              plot(plot_table$X1,type="n", xlab="Coverage", ylab="Probability",xlim=c(0,akcov*(p*2.5)), ylim=c(0,1), cex.lab=font_size, cex.axis=font_size, cex.main=font_size, cex.sub=font_size, yaxt='n')

              abline(v=peaks, col="black", lty=2, lwd=0.3)
            
              for (i in 1:(ncol(fitted_hist)+1)) {
              lines(plot_table[,i], type="l", lwd=1.5, col=colors[i])
              }
              axis(2, at=c(0,0.5,1), labels=c(0,0.5,1), cex.axis=font_size,
                 cex.lab=font_size)
          
            dev.off()
          
        }

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

## score model (now implemented directly in eval_model())

  #' Score nlsLM model by number and percent of residual errors after excluding sequencing errors
  #'
  #' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
  #' @param nls An nlsLM model object.
  #' @param round An integer corresponding to the iteration number (0, 1, 2, 3) for the fitting process.
  #' @param foldername A character vector corresponding to the name of the output directory.
  #' @return A data frame where the variables "all", "full", and "unique" correspond to the model residual sum of square error (excluding sequencing error)
  #' including x-values up to the end of kmer_hist_orig, up to (2p+1)*lambda, and up to (p+1)*lambda respectively.
  #' Additionally, the variables "allscore", "fullscore", and "uniquescore" are the
  #' corresponding percentage of kmers that are correctly modeled.
  #' @export
  # score_model<-function(kmer_hist_orig, nls, round, foldername, transform_exp) {

  #     x = kmer_hist_orig[[1]]
  #     y = kmer_hist_orig[[2]]
  #     y_transform = as.numeric(x)**transform_exp*as.numeric(y)

  #     pred=predict(nls, newdata=data.frame(x))
  #     model_sum=summary(nls)
  #     p=nls$p
  #     kcovfloor = max(1, floor(min_max(model_sum$coefficients['kmercov',])[[1]]))

  #     ## Compute error rate, by counting kmers unexplained by model through first peak
  #     ## truncate errors as soon as it goes to zero, dont allow it to go back up
  #       error_xcutoff = kcovfloor
  #       error_xcutoff_ind = tail(which(x<=error_xcutoff),n=1)
  #       if (length(error_xcutoff_ind)==0) {error_xcutoff_ind=1}

  #     error_kmers = x[1:error_xcutoff_ind]**(-transform_exp)*(y_transform[1:error_xcutoff_ind] - pred[1:error_xcutoff_ind])

  #     first_zero = -1

  #     for (i in 1:error_xcutoff_ind) {
  #       if (first_zero == -1) {
  #         if (error_kmers[i] < 1.0) {
  #           first_zero = i
  #         }
  #       }
  #       else {
  #         error_kmers[i] = 0
  #       }
  #     }

  #     if (first_zero == -1) {
  #       first_zero = error_xcutoff_ind
  #     }

  #     #if (TRANSFORM) {
  #     #  y_fit = y_transform
  #     #} else {
  #     #  y_fit = y
  #     #}
  #     y_fit = y_transform

  #     ## The fit is residual sum of square error, excluding sequencing errors
  #       model_fit_all    = c(sum(as.numeric(y_fit[first_zero:length(y_fit)]                          - pred[first_zero:length(y_fit)])                           ** 2), first_zero, x[length(y_fit)])
  #       model_fit_full   = c(sum(as.numeric(y_fit[first_zero:(min(length(y_fit),(2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]) ** 2), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
  #       model_fit_unique = c(sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                      - pred[first_zero:((p+1)*kcovfloor)])                       ** 2), first_zero, ((p+1)*kcovfloor))

  #     ## The score is the percentage of kmers correctly modeled, excluding sequencing errors
  #       model_fit_allscore    = c(1-sum(abs(as.numeric(y_fit[first_zero:length(y_fit)]                           - pred[first_zero:length(y_fit)])))                           / sum(as.numeric(y_fit[first_zero:length(y_fit)])),                           first_zero, x[length(y_fit)])
  #       model_fit_fullscore   = c(1-sum(abs(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))] - pred[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))]))) / sum(as.numeric(y_fit[first_zero:(min(length(y_fit), (2*p+1)*kcovfloor))])), first_zero, (min(length(y_fit), (2*p+1)*kcovfloor)))
  #       model_fit_uniquescore = c(1-sum(abs(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)]                       - pred[first_zero:((p+1)*kcovfloor)])))                       / sum(as.numeric(y_fit[first_zero:((p+1)*kcovfloor)])),                       first_zero, ((p+1)*kcovfloor))

  #     fit = data.frame(all  = model_fit_all,      allscore  = model_fit_allscore,
  #                      full = model_fit_full,     fullscore = model_fit_fullscore,
  #                      unique = model_fit_unique, uniquescore = model_fit_uniquescore)

  #     return (fit)
  # }

## nls_peak (now implemented directly in genomeScope())

  # #' Uses nlsLM to fit 2p peak model
  # #'
  # #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  # #' @param y A numeric vector of the y-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  # #' @param k An integer corresponding to the kmer length.
  # #' @param p An integer corresponding to the ploidy.
  # #' @param top An integer corresponding to the topology.
  # #' @param estKmercov A numeric corresponding to the estimated average kmer coverage of the polyploid genome.
  # #' @param estLength A numeric corresponding to the estimated polyploid genome length.
  # #' @param max_iterations An integer corresponding to the maximum number iterations to use for nlsLM.
  # #' @return An nlsLM model object with some additional components.
  # #' @export
  # nls_peak<-function(x, y, k, p, top, estKmercov, estLength, max_iterations) {
  #   #Initiate variables
  #   model = NULL
  #   best_deviance = Inf
  #   d_min = 0
  #   if (d_init!=-1) {
  #     d_initial = d_init
  #   } else {
  #     d_initial = 0.10
  #   }
  #   d_max = 1
  #   r_min = 0.00001
  #   if (top==0) {
  #     p_to_num_r = c(0, 1, 2, 4, 6, 10)
  #   } else {
  #     p_to_num_r = c(0, 1, 2, 3, 4, 5)
  #   }
  #   num_r = p_to_num_r[p]
  #   r_max = 1
  #   kmercov_min = 0
  #   kmercov_initial = estKmercov
  #   kmercov_max = Inf
  #   bias_min = 0
  #   bias_initial = 0.5
  #   bias_max = Inf
  #   length_min = 0
  #   length_initial = estLength/p
  #   length_max = Inf

  #   #Determine what formula to use, based on p
  #   if (p==1) {
  #     r_text = ""
  #   } else {
  #     r_text = paste(paste(lapply(1:(num_r), function(x) paste("r", as.character(x), sep="")), collapse=", "), ", ")
  #   }
  #   x = x[1:min(2000,length(x))]
  #   y = y[1:min(2000,length(y))]
  #   y_transform = as.numeric(x)**transform_exp*as.numeric(y)
  #   formula = as.formula(paste("y_transform ~ x**transform_exp*length*predict",p,"_",top,"(",r_text, "k, d, kmercov, bias, x)",sep=""))

  #   if (VERBOSE) {cat("trying nlsLM algorithm (Levenberg-Marquardt)\n")}

  #   if (r_inits!=-1) {
  #     r_initials = unlist(lapply(strsplit(r_inits,","),as.numeric))
  #     if (length(r_initials)!=num_r) {
  #       stop("Incorrect number of initial rates supplied.")
  #     }
  #     r_initials_list = list(r_initials)
  #   } else {
  #     r_initials_list = list(rep(0.001, num_r), 0.001*(1:num_r), 0.001*(num_r:1), rep(0.01, num_r), 0.01*(1:num_r), 0.01*(num_r:1))
  #   }

  #   for (r_initials in r_initials_list) {

  #     model1 = NULL
  #     r_start = vector("list", num_r)
  #     if (p > 1) {
  #       names(r_start) = paste("r", 1:(num_r), sep="")
  #       for (i in 1:(num_r)) {
  #         r_start[[paste("r",i,sep="")]] = r_initials[i]
  #       }
  #     }

  #     try(model1 <- nlsLM(formula = formula,
  #                        start   = c(list(d = d_initial), r_start, list(kmercov = kmercov_initial, bias = bias_initial, length = length_initial)),
  #                        lower   = c(c(d_min), rep(r_min, num_r), c(kmercov_min, bias_min, length_min)),
  #                        upper   = c(c(d_max), rep(r_max, num_r), c(kmercov_max, bias_max, length_max)),
  #                        control = list(minFactor=1e-12, maxiter=max_iterations, factor=0.1), trace=TRACE_FLAG), silent = TRUE)

  #     if (!is.null(model1)) {
  #       current_deviance = model1$m$deviance()
  #       #cat("Model deviance: ", current_deviance, "\n")
  #       if (current_deviance < best_deviance) {
  #         model = model1
  #         best_deviance = current_deviance
  #       }
  #     } else {
  #       #print("Model did not converge.")
  #     }

  #   }

  #   if (!is.null(model))
  #   {
  #     model_sum    = summary(model)
  #     model$p      = p
  #     model$top = top
  #     if (p==1) {
  #       model$hets = list(c(0, 0))
  #     } else {
  #       model$hets = lapply(1:(num_r), function(x) min_max1(model_sum$coefficients[paste('r', x, sep=""),]))
  #     }
  #     #model$het = c(1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 1))), 1-Reduce("*", 1-unlist(lapply(model$hets, '[[', 2))))
  #     model$het = c(sum(sapply(model$hets, '[[', 1)), sum(sapply(model$hets, '[[', 2)))
  #     model$homo = 1-model$het
  #     model$dups   = min_max(model_sum$coefficients['bias',])
  #     model$kcov   = min_max(model_sum$coefficients['kmercov',])
  #     model$mlen   = min_max(model_sum$coefficients['length',])
  #     model$md     = min_max1(model_sum$coefficients['d',])
  #     if (p==1) {
  #       model$ahets = list(c(0))
  #     } else {
  #       model$ahets = lapply(1:(num_r), function(x) model_sum$coefficients[paste('r', x, sep=""),][[1]])
  #     }
  #     #model$ahet = 1-Reduce("*", 1-unlist(model$ahets))
  #     model$ahet = Reduce("+", model$ahets)
  #     model$ahomo = 1-model$ahet
  #     model$adups = model_sum$coefficients['bias',][[1]]
  #     model$akcov = model_sum$coefficients['kmercov',][[1]]
  #     model$amlen = model_sum$coefficients['length',][[1]]
  #     model$amd   = model_sum$coefficients['d',][[1]]
  #   }

  #   #print(model)
  #   #print(model$m$deviance())

  #   return(model)
  # }

## estimate_Genome_peakp (now implemented directly in genomeScope())


  #' Function to fit 2p peak model, with p forms
  #'
  #' @param kmer_hist_orig A data frame of the original histogram data (starting at 1 and with last position removed).
  #' @param x An integer vector of the x-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @param y A numeric vector of the y-coordinates of the histogram (after filtering out low coverage errors and high coverage kmers).
  #' @param k An integer corresponding to the kmer length.
  #' @param p An integer corresponding to the ploidy.
  #' @param topology An integer corresponding to the topology to use.
  #' @param estKmercov An integer corresponding to the estimated kmer coverage of the polyploid genome.
  #' Set to -1 if not specified by user.
  #' @param round An integer corresponding to the iteration number (0, 1, 2, 3) for the fitting process.
  #' @param foldername A character vector corresponding to the name of the output directory.
  #' @param arguments A data frame of the user-specified inputs.
  #' @return A list (nls, nlsscore) where nls is the nlsLM model object (with some additional components)
  #' and nlsscore is the score (model RSSE) corresponding to the best fit (of the p forms).
  #' @export
    # estimate_Genome_peakp <- function(kmer_hist_orig, x, y, k, p, topology, estKmercov, round, foldername, arguments) {
    #   if (topology==-1) {
    #     p_to_num_topologies = c(1, 1, 1, 2, 5, 16)
    #     num_topologies = p_to_num_topologies[p]
    #     topologies = 1:num_topologies
    #   }
    #   else {
    #     num_topologies = 1
    #     topologies = c(topology)
    #   }
    #   numofKmers = sum(as.numeric(x)*as.numeric(y))
    #   if (estKmercov==-1) {
    #     #In situations with low heterozygosity, the peak with highest amplitude typically corresponds to the homozygous peak (i.e. the p-th peak).
    #     #However, with increasing heterozygosity, the highest amplitude peak may be an earlier peak.
    #     #Thus, when setting the estimated kmer coverage, we will need to iterate through these possibilities.
    #     #num_peak_indices indicates how many possibilities we need to iterate through.
    #     num_peak_indices = p
    #     y_transform = as.numeric(x)**transform_exp*as.numeric(y)
    #     estKmercov1 = x[which(y_transform==max(y_transform))][1]
    #   }
    #   else {
    #     # When the user sets the estimated kmer coverage, we only need to iterate through one possibility
    #     num_peak_indices = 1
    #     ## We set the estimated kmer coverage to be the user specified value
    #     estKmercov1 = estKmercov
    #   }
    #   estLength1 = numofKmers/estKmercov1

    #   nls00 = NULL
    #   peak_indices = 1:num_peak_indices
    #   for (i in peak_indices) {
    #     nls0 = NULL
    #     top_count = 0
    #     ## We see what happens when we set the estimated kmer coverage to be 1/i times the x-coordinate where the max peak occurs (1 <= i <= p if the user doesn't set the estimated kmer coverage, and i=1 if they do)
    #     estKmercov2 = estKmercov1 / i
    #     estLength2 = numofKmers/estKmercov2

    #     if (VERBOSE) {cat(paste("trying with kmercov: ", estKmercov2, "\n"))}

    #     for (top in topologies) {
    #       if (VERBOSE) {cat(paste("trying with topology: ", top, "\n"))}
    #       top_count = top_count + 1
    #       nls1 = nls_peak(x, y, k, p, top, estKmercov2, estLength2, MAX_ITERATIONS)
    #       nls0 = eval_model(kmer_hist_orig, nls0, nls1, p, round, foldername, arguments)[[1]]
    #     }
    #     if (i < num_peak_indices) { #if this is not the last evaluation
    #       nls00 = eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments)[[1]]
    #     }
    #   }

    #   return(eval_model(kmer_hist_orig, nls00, nls0, p, round, foldername, arguments))
    # }
