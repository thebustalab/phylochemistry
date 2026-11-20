## Moved verbatim from phylochemistry.R (phylochemistry website tree).

        #### analyzeGCMSdata2-6 (retired)

            ## Versions 2-6 are retired in favour of analyzeGCMSdata7(). These stubs stop with a
            ## pointer to v7. The old code is in z_archive/phylochemistry_analyzeGCMSdata2_3.R and
            ## bustalabfunctions/z_archive/. Plain analyzeGCMSdata() is not defined here: it comes from
            ## modules/gcms.R (old embedded copy in z_archive/phylochemistry_analyzeGCMSdata1.R).

            .retiredGCMSapp <- function(version) {
                force(version)
                function(...) {
                    stop(
                        "analyzeGCMSdata", version, "() has been retired. Please use analyzeGCMSdata7() instead.\n",
                        "Load it with:\n",
                        "  source(\"https://thebustalab.github.io/phylochemistry/bustalabfunctions/analyzeGCMSdata7.R\")\n",
                        "then run, for example:\n",
                        "  analyzeGCMSdata7(\"/path/to/your/CDF/folder\")",
                        call. = FALSE
                    )
                }
            }

            analyzeGCMSdata2 <- .retiredGCMSapp(2)
            analyzeGCMSdata3 <- .retiredGCMSapp(3)
            analyzeGCMSdata4 <- .retiredGCMSapp(4)
            analyzeGCMSdata5 <- .retiredGCMSapp(5)
            analyzeGCMSdata6 <- .retiredGCMSapp(6)

