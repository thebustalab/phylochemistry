## Moved from phylochemistry.R (phylochemistry website tree); the message now
## points at the bustalab package rather than the old source() URL.

        #### analyzeGCMSdata2-6 (retired)

            ## Versions 2-6 are retired in favour of analyzeGCMSdata7(). These stubs stop with a
            ## pointer to v7. Plain analyzeGCMSdata() is not defined here: it comes from
            ## modules/gcms.R on the lab website.

            .retiredGCMSapp <- function(version) {
                force(version)
                function(...) {
                    stop(
                        "analyzeGCMSdata", version, "() has been retired. Please use analyzeGCMSdata7() instead.\n",
                        "Run it with:\n",
                        "  library(bustalab); analyzeGCMSdata7(\"/path/to/your/CDF/folder\")",
                        call. = FALSE
                    )
                }
            }

            #' Retired GC-MS apps (versions 2 to 6)
            #'
            #' These versions have been retired. Each one stops with a message pointing at
            #' [analyzeGCMSdata7()], the current app.
            #'
            #' @param ... Ignored.
            #' @return Nothing; always signals an error.
            #' @name analyzeGCMSdata-retired
            #' @examples
            #' try(analyzeGCMSdata6("/path/to/your/CDF/folder"))
            NULL

            #' @rdname analyzeGCMSdata-retired
            #' @export
            analyzeGCMSdata2 <- .retiredGCMSapp(2)
            #' @rdname analyzeGCMSdata-retired
            #' @export
            analyzeGCMSdata3 <- .retiredGCMSapp(3)
            #' @rdname analyzeGCMSdata-retired
            #' @export
            analyzeGCMSdata4 <- .retiredGCMSapp(4)
            #' @rdname analyzeGCMSdata-retired
            #' @export
            analyzeGCMSdata5 <- .retiredGCMSapp(5)
            #' @rdname analyzeGCMSdata-retired
            #' @export
            analyzeGCMSdata6 <- .retiredGCMSapp(6)
