## shiny_ports — moved verbatim from phylochemistry.R (package conversion, commit 1)

#######################
## PHYLOCHEMISTRY v1 ##
#######################

###### Shiny port pool for jupyter = TRUE Shiny launchers


    .shiny_port_pool <- data.frame(
        port      = c(10101:10110, 10123, 10456, 10789),
        subdomain = c(paste0("shiny", 1:10), "shiny123", "shiny456", "shiny789"),
        stringsAsFactors = FALSE
    )


    .shiny_port_is_free <- function(port) {
        tryCatch({
            s <- httpuv::startServer("127.0.0.1", port,
                                     list(call = function(req) list(status = 200, body = "")))
            s$stop()
            TRUE
        }, error = function(e) FALSE)
    }


    pick_shiny_port <- function() {
        for (i in sample(nrow(.shiny_port_pool))) {
            p <- .shiny_port_pool$port[i]
            if (.shiny_port_is_free(p)) {
                return(list(
                    port = p,
                    url  = paste0("https://", .shiny_port_pool$subdomain[i], ".bustalab.d.umn.edu")
                ))
            }
        }
        NULL
    }

