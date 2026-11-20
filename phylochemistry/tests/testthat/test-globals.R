## Functions must not depend on, or silently write, objects in the user's global environment.

tree <- ape::read.tree(text = "((A:1,B:1):1,(C:1,(D:1,E:1):1):1);")

test_that("pruneTree() ignores a global fortified_tree", {
    assign("fortified_tree", data.frame(node = 1), envir = globalenv())
    on.exit(rm("fortified_tree", envir = globalenv()), add = TRUE)
    out <- pruneTree(tree, tips_to_keep = c("A", "B", "C"))
    expect_s3_class(out, "phylo")
    expect_setequal(out$tip.label, c("A", "B", "C"))
})

test_that("collapseTree() ignores a global fortified_tree", {
    assign("fortified_tree", data.frame(node = 1), envir = globalenv())
    on.exit(rm("fortified_tree", envir = globalenv()), add = TRUE)
    assoc <- data.frame(tip = c("A", "B", "C", "D", "E"), group = c("g1", "g1", "g2", "g3", "g3"))
    out <- suppressMessages(collapseTree(tree, assoc))
    expect_s3_class(out, "phylo")
    expect_setequal(out$tip.label, c("g1", "g2", "g3"))
})

test_that("pruneTree() on a fortified tree returns a fortified tree", {
    fort <- ggtree::fortify(tree)
    out <- pruneTree(fort, tips_to_keep = c("A", "B", "C"))
    expect_true(is.data.frame(out))
    expect_true("node" %in% colnames(out))
    expect_equal(nrow(out), nrow(ggtree::fortify(ape::drop.tip(tree, c("D", "E")))))
})

test_that("shapiroTest() keeps its noise flag local", {
    if (exists("modified", envir = globalenv(), inherits = FALSE)) rm("modified", envir = globalenv())
    d <- dplyr::group_by(data.frame(g = rep(c("a", "b"), each = 3), v = c(1, 1, 2, 3, 4, 6)), g)
    set.seed(1)
    expect_message(shapiroTest(d, v), "noise")
    expect_false(exists("modified", envir = globalenv(), inherits = FALSE))
})

test_that("pick_shiny_port() returns a pool port and URL, or NULL", {
    p <- pick_shiny_port()
    if (!is.null(p)) {
        expect_true(p$port %in% c(10101:10110, 10123, 10456, 10789))
        expect_match(p$url, "^https://shiny")
    } else {
        succeed()
    }
})
