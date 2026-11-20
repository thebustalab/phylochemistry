## Bootstrap support for hierarchical clustering, used by runMatrixAnalysis() (it calls Bclust()
## by that name for the "hclust" analyses).
##
## Bclust() comes from the shipunov package, which has been archived on CRAN, so it cannot be a
## hard dependency of this package. When shipunov is installed, its Bclust() is used unchanged.
## Otherwise this file computes the same bootstrap: the variables (columns) are resampled with
## replacement `iter` times, in the same order and with the same random draws, each resample is
## clustered with `FUN`, and every cluster of the original tree scores the fraction of resampled
## trees that contain exactly the same cluster. tests/testthat/test-bclust.R checks the two agree.

Bclust <- function(data, method.d = "euclidean", method.c = "ward.D",
                   FUN = function(.x) stats::hclust(stats::dist(.x, method = method.d), method = method.c),
                   iter = 1000, mc.cores = 1, monitor = TRUE, bootstrap = TRUE, relative = FALSE,
                   hclist = NULL) {
    if (requireNamespace("shipunov", quietly = TRUE)) {
        return(shipunov::Bclust(data = data, method.d = method.d, method.c = method.c, FUN = FUN,
                                iter = iter, mc.cores = mc.cores, monitor = monitor,
                                bootstrap = bootstrap, relative = relative, hclist = hclist))
    }
    .bclust_bootstrap(data = data, FUN = FUN, iter = iter, monitor = monitor,
                      bootstrap = bootstrap, relative = relative, hclist = hclist)
}

## Cluster membership of an hclust tree: one row per merge, one column per observation,
## 1 where the observation belongs to the cluster formed at that merge.
.hclust_membership <- function(tree) {
    merges <- tree$merge
    n_obs <- nrow(merges) + 1L
    member <- matrix(0L, nrow = nrow(merges), ncol = n_obs)
    for (step in seq_len(nrow(merges))) {
        for (side in 1:2) {
            node <- merges[step, side]
            if (node < 0L) {
                member[step, -node] <- 1L
            } else {
                member[step, ] <- pmax(member[step, ], member[node, ])
            }
        }
    }
    member
}

.bclust_bootstrap <- function(data, FUN, iter, monitor, bootstrap, relative, hclist) {
    if (!isTRUE(bootstrap) || isTRUE(relative) || !is.null(hclist)) {
        rlang::check_installed("shipunov", reason = "for Bclust() with bootstrap = FALSE, relative = TRUE or hclist")
    }
    fun <- match.fun(FUN)
    tree <- fun(data)
    original <- .hclust_membership(tree)
    n_obs <- ncol(original)
    key <- function(m) apply(m, 1, paste, collapse = "")
    original_keys <- key(original)
    hits <- numeric(nrow(original))
    together <- matrix(0, nrow = nrow(original), ncol = n_obs)
    for (rep in seq_len(iter)) {
        resampled <- .hclust_membership(fun(data[, sample.int(ncol(data), replace = TRUE)]))
        hits <- hits + (original_keys %in% key(resampled))
        together <- together + resampled
        if (monitor) cat(".")
    }
    if (monitor) cat("\n")
    consensus <- t(together) / iter
    row.names(consensus) <- row.names(data)
    out <- list(values = hits / iter, hclust = tree, consensus = consensus, meth = "Bootstrap", iter = iter)
    class(out) <- "Bclust"
    out
}
