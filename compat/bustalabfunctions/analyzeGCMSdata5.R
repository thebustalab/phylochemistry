#### analyzeGCMSdata5 (retired)
####
#### This version has been retired in favour of analyzeGCMSdata7(), which is now part of
#### the bustalab package. This stub stays at the old URL so anyone still sourcing it is
#### told what to use instead.

    analyzeGCMSdata5 <- function(...) {
        stop(
            "analyzeGCMSdata5() has been retired. Please use analyzeGCMSdata7() instead.\n",
            "Run it with:\n",
            "  library(bustalab); analyzeGCMSdata7(\"/path/to/your/CDF/folder\")\n",
            "(If bustalab is not installed yet, run  bustalab <- TRUE  and then\n",
            "  source(\"https://thebustalab.github.io/phylochemistry/phylochemistry.R\")  once.)",
            call. = FALSE
        )
    }

    message("Note: analyzeGCMSdata5() has been retired. Please use library(bustalab); analyzeGCMSdata7() instead.")
