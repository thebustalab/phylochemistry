## Every name a package function CALLS (f(...)) must resolve in the package, its imports or base R.
## R CMD check's codetools pass misses a call whose name is also a LOCAL VARIABLE: buildModel() did
## `fit <- lm(...)` in one branch and called fit() (parsnip, never imported) in another, so check was
## clean while the random-forest paths failed for every student. Here a local counts only if it is
## assigned a function (function, reactive(), match.fun(), ...) or is an argument, never a plain value.

.call_heads <- function(e) {
    out <- character()
    walk <- function(x) {
        if (is.call(x)) {
            if (is.symbol(x[[1]])) out <<- c(out, as.character(x[[1]]))
            for (i in seq_along(x)) { el <- x[[i]]; if (!missing(el)) walk(el) }
        } else if (is.pairlist(x) || is.list(x)) {
            for (el in x) if (!missing(el)) walk(el)
        }
    }
    walk(e)
    unique(out)
}

## names bound to something callable inside a function body: `f <- function`, `f <- reactive(...)`,
## and the arguments of the function and of every nested function
.local_callables <- function(fn) {
    makers <- c("function", "reactive", "eventReactive", "reactiveVal", "reactiveValues", "match.fun", "Vectorize")
    out <- names(formals(fn))
    walk <- function(x) {
        if (is.call(x)) {
            h <- x[[1]]
            if (is.symbol(h) && as.character(h) %in% c("<-", "=", "<<-") && length(x) == 3 && is.symbol(x[[2]])) {
                rhs <- x[[3]]
                if (is.call(rhs) && is.symbol(rhs[[1]]) && as.character(rhs[[1]]) %in% makers) {
                    out <<- c(out, as.character(x[[2]]))
                }
            }
            if (is.symbol(h) && identical(as.character(h), "function")) out <<- c(out, names(x[[2]]))
            for (i in seq_along(x)) { el <- x[[i]]; if (!missing(el)) walk(el) }
        } else if (is.pairlist(x) || is.list(x)) {
            for (el in x) if (!missing(el)) walk(el)
        }
    }
    walk(body(fn))
    unique(out)
}

test_that("every bare function call in the package resolves without the user attaching anything", {
    ns <- asNamespace("phylochemistry")
    base_pkgs <- c("base", "stats", "utils", "graphics", "grDevices", "methods")
    known <- unique(c(
        unlist(lapply(base_pkgs, getNamespaceExports)), ls(baseenv(), all.names = TRUE),
        ls(ns, all.names = TRUE), ls(parent.env(ns), all.names = TRUE)
    ))
    ## Deliberate exceptions, each reached only on a path that sets it up:
    allowed <- c(
        "picker", "pickerOutput", "renderPicker", # analyzeImages(): the non-CRAN picker package
        "shell",                                  # polyBlast(): base R on Windows only
        ".",                                      # readGFFs(): library(plyr) runs first
        "embedding_model"                         # predictWithModel(): the torch module model$architecture() returns
    )
    unresolved <- character()
    for (nm in ls(ns, all.names = TRUE)) {
        f <- get(nm, envir = ns)
        if (!is.function(f) || is.primitive(f)) next
        bad <- setdiff(c(.call_heads(body(f)), .call_heads(formals(f))),
                       c(known, allowed, .local_callables(f)))
        if (length(bad)) unresolved <- c(unresolved, paste0(nm, "(): ", bad))
    }
    expect_identical(unresolved, character())
})
