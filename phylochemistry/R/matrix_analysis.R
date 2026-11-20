## matrix_analysis — moved verbatim from phylochemistry.R (package conversion, commit 1)


    ##### Networks

        #### extractModules

            #' Extract modules around bait genes from an edge list
            #'
            #' @param data Not yet documented.
            #' @param baits Not yet documented.
            #' @param module_size Not yet documented.
            #' @export
            extractModules <- function(data, baits, module_size) {
                ## All this code assumes a Gene1-Gene2-Metric format    
                    connections_to_baits <- list()
                    for( i in 1:length(baits) ) {
                        ## Find entries in the coexpression table that include the bait as Gene1 or Gene2
                        connections_to_bait <- unique(rbind(data[data[,1] == baits[i],], data[data[,2] == baits[i],]))
                        ## Order entries according to "Metric" and then take the top n (module_size) entries
                        connections_to_baits[[i]] <- connections_to_bait[order(connections_to_bait[,3], decreasing = TRUE),][1:module_size,]
                    }
                    connections_to_baits <- do.call(rbind, connections_to_baits)

                ## Append all connections between module members so you dont get a sun shape
                    data[apply(cbind(
                        data[,1] %in% c(connections_to_baits[,1], connections_to_baits[,2]),
                        data[,2] %in% c(connections_to_baits[,1], connections_to_baits[,2])
                    ), 1, all),] -> module_member_connections
                    connections <- unique(rbind(connections_to_baits, module_member_connections))
                    return(connections[!is.na(connections[,3]),])
            }


        #### buildNetwork_OG


            #' Build a network from an edgelist
            #'
            #' @param edgelist A data frame with names of origin nodes in first column, destination nodes in second column, and an optional third column of edgeweights.
            #' @param node_attributes A dataframe of attributes associated with the nodes.
            #' @param facet_variable Include in the output graphic a line plot corresponding to the degree of conservation at each site in the alignment.
            #' @examples
            #' \dontrun{
            #' buildNetwork_OG()
            #' }
            #' @export
            buildNetwork_OG <- function( edgelist, node_attributes = NULL, facet_variable = NULL ) {

                ## Make sure the edgelist is a dataframe and that its first two columns are characters
                    edgelist <- as.data.frame(edgelist)
                    edgelist[,1] <- as.character(edgelist[,1])
                    edgelist[,2] <- as.character(edgelist[,2])

                ## If edgelist has more than two columns, assume that 3rd column is edge weights and build matrix from adjacency matrix
                    if (dim(edgelist)[2] > 2) {
                        adjacency_matrix <- as.matrix(as.data.frame(tidyr::pivot_wider(edgelist, names_from = 1, values_from = 3)))
                        rownames(adjacency_matrix) <- adjacency_matrix[,1]
                        adjacency_matrix <- adjacency_matrix[,-1]
                        network_object <- network::network(adjacency_matrix, matrix.type = "adjacency")
                        cat("If you get a network building error you might have missing edges from filtering the dataset before network building. Just filter afterword instead.")
                    }

                ## Just build a network with the edgelist
                    if (dim(edgelist)[2] == 2) {
                        network_object <- network::network(edgelist, matrix.type = "edgelist", ignore.eval = FALSE)
                    }

                ## Fortify the network object into a dataframe, modify it with new colnames and additional information, crack into a list of edges and nodes
                    combined_network_frame <- ggnetwork::ggnetwork(network_object, arrow.gap = 0, layout = "kamadakawai", stringsAsFactors = TRUE)
                    # combined_network_frame <- ggnetwork::ggnetwork(network_object, arrow.gap = 0)
                    combined_network_frame$type <- "edge"
                    colnames(combined_network_frame)[colnames(combined_network_frame) == "vertex.names"] <- "start_node"
                    combined_network_frame$type[apply(cbind(combined_network_frame$x == combined_network_frame$xend, combined_network_frame$y == combined_network_frame$yend), 1, all)] <- "node"
                    combined_network_frame$end_node <- combined_network_frame$start_node[match(combined_network_frame$xend, combined_network_frame$x)]

                    network_frame <- list()
                    network_frame$edges <- combined_network_frame[combined_network_frame$type == "edge",]
                    network_frame$edges <- network_frame$edges[,colnames(network_frame$edges) %in% c("x", "y", "start_node", "xend", "yend", "end_node")]
                    network_frame$nodes <- combined_network_frame[combined_network_frame$type == "node",]
                    network_frame$nodes <- network_frame$nodes[,colnames(network_frame$nodes) %in% c("x", "y", "start_node")]
                    colnames(network_frame$nodes)[3] <- "node_name"

                    combined_network_frame$end_node <- as.character(combined_network_frame$end_node)

                ## Bind edge attributes to the network_frame
                    if (dim(edgelist)[2] > 2) {
                        edge_attributes <- list()
                        for ( edge in 1:dim(edgelist)[1] ) {
                            network_frame_entry <- which(apply(cbind(combined_network_frame$start_node == edgelist[edge,1], combined_network_frame$end_node == edgelist[edge,2]), 1, all))
                            edge_attributes[[network_frame_entry]] <- as.character(edgelist[edge, 3:dim(edgelist)[2]])
                        }
                        edge_attributes <- as.data.frame(do.call(rbind, edge_attributes))
                        colnames(edge_attributes) <- colnames(edgelist)[3:dim(edgelist)[2]]
                        network_frame$edges <- cbind(network_frame$edges, edge_attributes)
                    }

                ## Merge node_attributes frame with network_frame$nodes
                    if (!is.null(node_attributes)) {
                        network_frame$nodes <- cbind(network_frame$nodes, node_attributes[match(network_frame$nodes$node_name, node_attributes[,1]),])
                    }

                ## Create network_frame_expanded, which contains entries for each facet_variable
                    # network_frame_expanded <- data.frame()
                    # facet_variables <- as.character(unique(annotation_frame[,which(colnames(annotation_frame) == facet_variable)]))

                    # if ( length(facet_variables) > 0 ) {
                    #     for (i in 1:length(facet_variables)) {
                    #         frame_subset <- annotation_frame[annotation_frame[,which(colnames(annotation_frame) == facet_variable)] == facet_variables[i],]
                    #         frame_subset <- frame_subset[frame_subset$abundance > 0,]
                            
                    #         # Show only nodes present in that facet_variable
                    #             # nodes <- network_frame[apply(cbind(network_frame$type == "node", network_frame$start_node %in% frame_subset$compound_name), 1, all),]

                    #         # Show all nodes
                    #             nodes <- network_frame[network_frame$type == "node",]
                            
                    #         # Show only edges present in that facet variable
                    #             edges <- network_frame[apply(   cbind(
                    #                                                 network_frame$type == "edge", 
                    #                                                 network_frame$start_node %in% frame_subset$compound_name,
                    #                                                 network_frame$end_node %in% frame_subset$compound_name
                    #                                             ), 1, all),]

                    #         temp <- rbind(nodes, edges)
                    #         temp$facet_variable <- facet_variables[i]
                    #         colnames(temp)[colnames(temp) == "facet_variable"] <- facet_variable

                    #         network_frame_expanded <- rbind(network_frame_expanded, temp)
                    #     }
                    # } else {
                    #     network_frame_expanded <- network_frame
                    # }

                return(network_frame)
            }


        #### collapseReciprocalEdges


            #' Collapse duplicate and reciprocal edges in an edgelist (used by buildNetwork)
            #'
            #' Rows that are identical in every column are always dropped. Beyond that, each
            #' row's node pair is put in alphabetical order so A->B and B->A share a key. If
            #' every key's rows agree on all attribute columns (columns 3+; numbers compared
            #' with a small relative tolerance), the input is symmetric and is collapsed to one
            #' row per pair (the first occurrence, in its original direction). If any key's rows
            #' disagree, the input is treated as directed and returned whole. It is all or
            #' nothing across the table, never pair by pair.
            #'
            #' @param edgelist A data frame, node names (character) in columns 1 and 2.
            #' @param directed NULL (decide from the data, as above), TRUE (never collapse
            #'   reciprocal pairs), or FALSE (always collapse, keeping the first occurrence).
            #' @param tol Relative tolerance for comparing numeric attribute columns.
            #' @export
            collapseReciprocalEdges <- function(edgelist, directed = NULL, tol = sqrt(.Machine$double.eps)) {

                n_input <- nrow(edgelist)
                edgelist <- edgelist[!duplicated(edgelist), , drop = FALSE]
                n_exact <- n_input - nrow(edgelist)
                if (n_exact > 0) {
                    message("buildNetwork: dropped ", n_exact, " edge row(s) that were exact copies of another row.")
                }

                pair_key <- paste(
                    pmin(edgelist[[1]], edgelist[[2]]),
                    pmax(edgelist[[1]], edgelist[[2]]),
                    sep = "\r"
                )
                repeated <- duplicated(pair_key)
                if (isTRUE(directed) || !any(repeated)) {
                    rownames(edgelist) <- NULL
                    return(edgelist)
                }

                # For every row, does it match the first row sharing its pair key on all attributes?
                first_row <- match(pair_key, pair_key)
                agrees <- rep(TRUE, nrow(edgelist))
                if (ncol(edgelist) > 2) {
                    for (i in 3:ncol(edgelist)) {
                        x <- edgelist[[i]]
                        if (is.numeric(x)) {
                            ref <- x[first_row]
                            same <- (is.na(x) & is.na(ref)) |
                                (!is.na(x) & !is.na(ref) & abs(x - ref) <= tol * pmax(1, abs(ref)))
                        } else {
                            x <- as.character(x)
                            ref <- x[first_row]
                            same <- (is.na(x) & is.na(ref)) | (!is.na(x) & !is.na(ref) & x == ref)
                        }
                        agrees <- agrees & same
                    }
                }
                n_conflicting_pairs <- length(unique(pair_key[!agrees]))

                if (is.null(directed) && n_conflicting_pairs > 0) {
                    message(
                        "buildNetwork: ", n_conflicting_pairs, " node pair(s) appear more than once with ",
                        "different edge attributes (e.g. A->B and B->A with different weights), so all ",
                        nrow(edgelist), " edges were kept and the network is treated as directed. ",
                        "Set directed = FALSE to keep one edge per pair (the first occurrence) instead."
                    )
                    rownames(edgelist) <- NULL
                    return(edgelist)
                }

                collapsed <- edgelist[!repeated, , drop = FALSE]
                rownames(collapsed) <- NULL
                if (n_conflicting_pairs > 0) {
                    message(
                        "buildNetwork: directed = FALSE, so ", nrow(edgelist), " edges were collapsed to ",
                        nrow(collapsed), ", one per node pair. ", n_conflicting_pairs, " pair(s) had ",
                        "differing edge attributes; the first occurrence of each was kept."
                    )
                } else {
                    message(
                        "buildNetwork: every repeated node pair (A->B and B->A) had matching edge ",
                        "attributes, so the network is treated as undirected and ", nrow(edgelist),
                        " edges were collapsed to ", nrow(collapsed), ". ",
                        "Set directed = TRUE to keep both directions."
                    )
                }
                collapsed
            }


        #### buildNetwork


            #' Build a network from an edgelist
            #'
            #' @param edgelist A data frame with origin nodes in the first column, destination nodes in the second column, and optional edge attributes in columns 3+, which are carried through to the returned edges frame as-is. No `weight` column is invented: for the force-directed layout, a column named `edgeweight` or `weight` is used if present, otherwise column 3.
            #' @param node_attributes A dataframe of attributes associated with the nodes. First column must contain node names.
            #' @param facet_variable Currently unused; retained for API compatibility.
            #' @param directed How to handle node pairs that appear more than once. Rows identical in every column are always dropped. NULL (default) decides from the data: if every A->B / B->A pair (and any same-direction repeat) agrees on all edge attributes, as in a long-format distance or correlation matrix, the edgelist is collapsed to one edge per pair; if any pair disagrees, every edge is kept. TRUE always keeps both directions; FALSE always collapses, keeping the first occurrence of each pair. A message reports what was done. Note that when both directions are kept, their segments are drawn on top of one another, so a plot will only show one of them.
            #' @examples
            #' \dontrun{
            #' buildNetwork()
            #' }
            #' @export
            buildNetwork <- function(edgelist, node_attributes = NULL, facet_variable = NULL, directed = NULL) {

                edgelist <- as.data.frame(edgelist, stringsAsFactors = FALSE)
                if (ncol(edgelist) < 2) {
                    stop("`edgelist` must have at least two columns: start node and end node.")
                }
                if (!is.null(directed) && !(isTRUE(directed) || isFALSE(directed))) {
                    stop("`directed` must be NULL, TRUE or FALSE.")
                }

                edgelist[[1]] <- as.character(edgelist[[1]])
                edgelist[[2]] <- as.character(edgelist[[2]])

                edgelist <- collapseReciprocalEdges(edgelist, directed = directed)

                # Build from an explicit edgelist so row-level edge attributes are preserved.
                # multiple = TRUE: a same-direction repeat with differing attributes survives
                # collapseReciprocalEdges() and would otherwise make network() stop on parallel edges.
                network_object <- network::network(
                    edgelist[, 1:2, drop = FALSE],
                    matrix.type = "edgelist",
                    ignore.eval = TRUE,
                    multiple = TRUE
                )

                if (ncol(edgelist) > 2) {
                    for (i in 3:ncol(edgelist)) {
                        attr_name <- colnames(edgelist)[i]
                        if (is.null(attr_name) || is.na(attr_name) || attr_name == "") {
                            attr_name <- paste0("edge_attr_", i - 2)
                        }
                        network::set.edge.attribute(network_object, attr_name, edgelist[[i]])
                    }
                }

                # Prefer a weighted force-directed layout when usable edge weights are present.
                layout_coords <- NULL
                if (ncol(edgelist) > 2 && requireNamespace("igraph", quietly = TRUE)) {
                    weight_col <- c("edgeweight", "weight")
                    weight_col <- weight_col[weight_col %in% colnames(edgelist)]
                    if (length(weight_col) == 0) {
                        weight_col <- colnames(edgelist)[3]
                    } else {
                        weight_col <- weight_col[1]
                    }

                    layout_weights <- suppressWarnings(as.numeric(edgelist[[weight_col]]))
                    if (any(!is.na(layout_weights))) {
                        min_pos <- suppressWarnings(min(layout_weights[layout_weights > 0], na.rm = TRUE))
                        if (!is.finite(min_pos)) {
                            min_pos <- 1e-06
                        }
                        layout_weights[is.na(layout_weights) | layout_weights <= 0] <- min_pos

                        graph_for_layout <- igraph::graph_from_data_frame(
                            data.frame(
                                from = edgelist[[1]],
                                to = edgelist[[2]],
                                weight = layout_weights,
                                stringsAsFactors = FALSE
                            ),
                            directed = FALSE
                        )

                        coords <- igraph::layout_with_fr(
                            graph_for_layout,
                            weights = igraph::E(graph_for_layout)$weight,
                            niter = 2000
                        )
                        coords <- as.matrix(coords)
                        rownames(coords) <- igraph::V(graph_for_layout)$name

                        vertex_names <- network::get.vertex.attribute(network_object, "vertex.names")
                        coords <- coords[match(vertex_names, rownames(coords)), , drop = FALSE]

                        if (!any(is.na(coords))) {
                            layout_coords <- coords
                        }
                    }
                }

                combined_network_frame <- ggnetwork::ggnetwork(
                    network_object,
                    arrow.gap = 0,
                    layout = if (is.null(layout_coords)) "kamadakawai" else layout_coords,
                    stringsAsFactors = FALSE
                )

                is_node <- (combined_network_frame$x == combined_network_frame$xend) &
                    (combined_network_frame$y == combined_network_frame$yend)

                nodes <- combined_network_frame[is_node, c("x", "y", "vertex.names"), drop = FALSE]
                colnames(nodes)[3] <- "node_name"
                nodes <- nodes[!duplicated(nodes$node_name), , drop = FALSE]

                node_lookup <- nodes[, c("node_name", "x", "y"), drop = FALSE]
                colnames(node_lookup) <- c("node_name", "node_x", "node_y")

                edge_attributes <- NULL
                if (ncol(edgelist) > 2) {
                    edge_attributes <- edgelist[, 3:ncol(edgelist), drop = FALSE]
                }

                edges <- data.frame(
                    start_node = edgelist[[1]],
                    end_node = edgelist[[2]],
                    stringsAsFactors = FALSE
                )
                if (!is.null(edge_attributes)) {
                    edges <- cbind(edges, edge_attributes)
                }

                start_idx <- match(edges$start_node, node_lookup$node_name)
                end_idx <- match(edges$end_node, node_lookup$node_name)

                edges$x <- node_lookup$node_x[start_idx]
                edges$y <- node_lookup$node_y[start_idx]
                edges$xend <- node_lookup$node_x[end_idx]
                edges$yend <- node_lookup$node_y[end_idx]

                if (any(is.na(edges$x) | is.na(edges$y) | is.na(edges$xend) | is.na(edges$yend))) {
                    warning("Some edges could not be mapped to node coordinates; check edgelist node labels.")
                }

                edge_cols <- c("x", "y", "start_node", "xend", "yend", "end_node")
                extra_cols <- setdiff(colnames(edges), edge_cols)
                edges <- edges[, c(edge_cols, extra_cols), drop = FALSE]

                if (!is.null(node_attributes)) {
                    node_attributes <- as.data.frame(node_attributes, stringsAsFactors = FALSE)
                    node_key <- as.character(node_attributes[[1]])
                    matched_attrs <- node_attributes[match(nodes$node_name, node_key), , drop = FALSE]

                    # Avoid duplicating the node identifier column.
                    if (ncol(matched_attrs) > 0) {
                        matched_attrs <- matched_attrs[, -1, drop = FALSE]
                    }
                    nodes <- cbind(nodes, matched_attrs)
                }

                if (!is.null(facet_variable)) {
                    warning("`facet_variable` is currently not implemented in buildNetwork().")
                }

                return(list(edges = edges, nodes = nodes))
            }


    ##### Mathematics, Statistical Testing, Modeling, Signal Processing

        #### makeRectangular

            #' Fill in missing combinations of a long table
            #'
            #' @param dataframe Not yet documented.
            #' @param names_from Not yet documented.
            #' @param values_from Not yet documented.
            #' @param values_fill Not yet documented.
            #' @export
            makeRectangular <- function(dataframe, names_from, values_from, values_fill = c(0, NA)) {
                wide_data <- pivot_wider(dataframe, names_from = names_from, values_from = values_from, values_fill = values_fill[1])
                names <- unique(unlist(dataframe[,colnames(dataframe) == names_from]))
                out <- pivot_longer(wide_data, cols = names, values_to = values_from, names_to = names_from)
                return(out)  
            }


        #### normalize


            #' Normalizes a vector of numbers to a range of zero to one.
            #'
            #' @param x The vector to normalize
            #' @param old_min The minimum of the old range
            #' @param old_max The maximum of the old range
            #' @param new_min The minimum of the new range, defaults to 0
            #' @param new_max The maximum of the new range, defaults to 1
            #' @param na_zero Not yet documented.
            #' @export
            normalize <- function( x, old_min = NULL, old_max = NULL, new_min = 0, new_max = 1, na_zero = FALSE ) {

                if ( length(old_min) == 0 & length(old_max) == 0 ) {
                    
                    x <- (x - min(x)) * (new_max - new_min) / (max(x) - min(x)) + new_min
                
                } else {

                    x <- (x - old_min) * (new_max - new_min) / (old_max - old_min) + new_min

                }

                if (na_zero == TRUE) {
                    x[is.na(x)] <- 0
                }

                return(x)

            }


        #### movingAverage


            #' Compute a moving average
            #'
            #' movingAverage
            #' @param x Not yet documented.
            #' @param n Not yet documented.
            #' @export
            movingAverage <- function(x, n = 5){
            
                stats::filter(x, rep(1 / n, n), sides = 2)
            
            }


        #### bootstrap


            #' Runs bootstrapping on categorical data; from github sgibb/bootstrap
            #'
            #' @param x Not yet documented.
            #' @param fun Not yet documented.
            #' @param n
            #' bootstrap
            #' @param mc.cores Not yet documented.
            #' @export
            bootstrap <- function(x, fun, n = 1000L, mc.cores = getOption("mc.cores", 2L)) {
                # Match the function if given as a character string
                fun <- match.fun(fun)

                # Internal function to convert hclust object to binary matrix
                as.binary.matrix.hclust <- function(x) {
                    nr <- as.integer(nrow(x$merge))
                    m <- matrix(0L, nrow = nr, ncol = nr + 1L)
                    for (i in seq_len(nr)) {
                      left <- x$merge[i, 1L]
                      if (left < 0L) {
                        m[i, -left] <- 1L
                      } else {
                        m[i, ] <- m[left, ]
                      }
                      right <- x$merge[i, 2L]
                      if (right < 0L) {
                        m[i, -right] <- 1L
                      } else {
                        m[i, ] <- m[i, ] | m[right, ]
                      }
                    }
                    return(m)
                }

                # Internal function to perform clustering
                .clust <- function(x, fun) {
                    hc <- fun(x)
                    return(as.binary.matrix.hclust(hc))
                }

                # Internal function to calculate matches
                .calculateMatches <- function(origin, current, nc = ncol(origin)) {
                    one <- tcrossprod(origin, current)
                    zero <- tcrossprod(1 - origin, 1 - current)
                    return(rowSums(one + zero == nc))
                }

                # Internal function for resampling
                .resample <- function(x, size = ncol(x)) {
                    sel <- sample.int(ncol(x), size = size, replace = TRUE)
                    return(x[, sel])
                }

                origin <- .clust(x, fun = fun)

                v <- parallel::mclapply(seq_len(n), function(y) .calculateMatches(origin, .clust(.resample(x, size = ncol(x)), fun = fun), ncol(origin)), mc.cores = mc.cores)

                return(colSums(do.call(rbind, v)) / n)
            }


        #### runMatrixAnalyses


                ## ── RETIRED 2026-09-22 ────────────────────────────────────────────────────────────────
                ## runMatrixAnalyses() (plural) was a wide-data-only fork of runMatrixAnalysis()
                ## (singular). Two near-identical names for "the simple one" and "the full one" was a
                ## reliable source of student confusion, and the two returned DIFFERENT TYPES for the
                ## same analysis = "dist", which is worse than confusing.
                ##
                ## Everything it could do, the singular now does:
                ##   - tree_method = c("linkage_dendrogram", "neighbor_joining")  -> ported
                ##   - distance_method = "gower" via cluster::daisy               -> ported
                ##   - analysis = "dist" returning a long data frame              -> ported, now opt-in
                ##     via output_format = "long" (the default "wide" gives a base R `dist` object)
                ## and the singular additionally takes long-format input, transpose,
                ## unknown_sample_ID_info and columns_w_additional_analyte_info, none of which the
                ## plural ever had. The plural also carried a live bug: output_format = "long"
                ## referenced two variables that did not exist in it, so it always errored.
                ##
                ## It is kept as a stop() rather than deleted so that an old script gets a sentence
                ## telling it what to do, instead of "could not find function".
                #' RETIRED. Use runMatrixAnalysis() (singular) instead.
                #'
                #' Retired 2026-09-22. This was a wide-data-only fork of runMatrixAnalysis(); every
                #' capability it had has been ported into the singular, which also handles long-format
                #' input. Calling it raises an error explaining the migration.
                #'
                #' @param ... Ignored; the function always stops.
                #' @export
                runMatrixAnalyses <- function(...) {
                    stop(
                        "runMatrixAnalyses() (plural) was retired on 2026-09-22. Use runMatrixAnalysis() ",
                        "(singular) instead \u2014 it does everything the plural did and more.\n\n",
                        "  - Wide data: the arguments are the same. Rename the function and it will run.\n",
                        "  - analysis = \"dist\": the singular returns a base R `dist` object by default. ",
                        "For the long, one-row-per-pair data frame the plural returned, add ",
                        "output_format = \"long\".\n",
                        "  - Hierarchical clustering: tree_method now defaults to \"linkage_dendrogram\" ",
                        "rather than \"neighbor_joining\". Pass tree_method explicitly if you want NJ.\n",
                        call. = FALSE
                    )
                }


        #### runMatrixAnalysis

            #' The ggproto Stat behind geom_cut()
            #'
            #' Places the cut line on a dendrogram from the plot's own `x` data; used by [geom_cut()].
            #' @format A ggplot2 `ggproto` object.
            #' @export
            StatCut <- ggplot2::ggproto("StatCut", ggplot2::Stat,
                required_aes = "x",
                compute_panel = function(data, scales, height = NULL, k = NULL) {

                    if (is.null(height) && is.null(k)) {
                        stop("geom_cut() needs either the height the tree was cut at -- geom_cut(height = 5) -- or the number of groups -- geom_cut(k = 4).")
                    }
                    if (!is.null(height) && !is.null(k)) {
                        stop("geom_cut() takes height OR k, not both.")
                    }

                    tip_x <- max(data$x, na.rm = TRUE)

                    if (!is.null(height)) {
                        ## The tips sit at max(height)/2, so a merge at `height` is that far back.
                        return(data.frame(xintercept = tip_x - height / 2))
                    }

                    ## k groups means k-1 joins have been severed, so the line sits between the
                    ## (k-1)th and kth join counting in from the root. Tips share max(x); anything
                    ## to the left of that is a join.
                    joins <- sort(data$x[data$x < tip_x])
                    k <- as.integer(k)
                    if (is.na(k) || k < 1) stop("geom_cut(k = ) needs a positive number of groups.")
                    if (k == 1)               return(data.frame(xintercept = min(joins) - 0.01))
                    if (k > length(joins) + 1) stop(paste0("This tree cannot be cut into ", k, " groups."))
                    lo <- joins[k - 1]
                    hi <- if (k <= length(joins)) joins[k] else tip_x
                    data.frame(xintercept = (lo + hi) / 2)
                }
            )


            #' geom_cut
            #'
            #' Draw the cut on a dendrogram, in the right place, from the plot's own data.
            #'
            #' WHY THIS EXISTS (2026-09-23). A cut IS a vertical line across the tree, and saying so
            #' with `geom_vline(xintercept = 5)` is wrong twice over: `ape::as.phylo()` HALVES an
            #' hclust's merge heights, and `ggtree` then measures x from the ROOT, so the tips sit at
            #' `max(height)/2` and a cut at height 5 actually belongs at `max(x) - 5/2`. Drawn the naive
            #' way the line lands out past the tip labels. `geom_cut()` does that conversion, so the
            #' student writes the same height they cut at and the picture agrees with the `cluster`
            #' column.
            #'
            #' It reads only `x`, which ggtree always maps, so nothing extra has to be threaded into
            #' the layer. Tips are the rightmost nodes, so everything left of `max(x)` is a join.
            #'
            #' @param height Cut height, the same value passed to `parameters = c(height = ...)`.
            #' @param k Number of groups, the same value passed to `parameters = c(n_clusters = k)`
            #'   (or the older `c(k)` / `c(k = k)`). The line is drawn midway between the two joins
            #'   that bracket that number of groups.
            #' @param linetype Not yet documented.
            #' @param colour Not yet documented.
            #' @param ... Not yet documented.
            #' @examples
            #' \dontrun{
            #'   ggtree(out) + geom_tippoint(aes(fill = cluster)) + geom_cut(height = 5)
            #'   ggtree(out) + geom_tippoint(aes(fill = cluster)) + geom_cut(k = 4)
            #' }
            #' @export
            geom_cut <- function(height = NULL, k = NULL, linetype = "dashed",
                                 colour = "grey40", ...) {

                ## Validate HERE, not in the Stat. ggplot2 turns a Stat error into a WARNING and
                ## draws the plot anyway, so a mistyped call would quietly produce a tree with no
                ## line on it -- the one failure a student would not notice. Checking at call time
                ## makes it a real error. (The "k is larger than this tree" case still has to wait
                ## for the data, so it stays a warning.)
                if (is.null(height) && is.null(k)) {
                    stop("geom_cut() needs either the height the tree was cut at -- geom_cut(height = 5) -- or the number of groups -- geom_cut(k = 4).")
                }
                if (!is.null(height) && !is.null(k)) {
                    stop("geom_cut() takes height OR k, not both.")
                }
                if (!is.null(height) && (!is.numeric(height) || is.na(height[1]) || height[1] <= 0)) {
                    stop("geom_cut(height = ) needs a positive height.")
                }
                if (!is.null(k) && (!is.numeric(k) || is.na(k[1]) || k[1] < 1)) {
                    stop("geom_cut(k = ) needs a positive number of groups.")
                }

                ggplot2::layer(
                    stat = StatCut, geom = ggplot2::GeomVline,
                    data = NULL, mapping = NULL, position = "identity",
                    show.legend = FALSE, inherit.aes = TRUE,
                    params = list(height = height, k = k, linetype = linetype,
                                  colour = colour, na.rm = TRUE, ...)
                )
            }


                #' Runs a matrix analysis (clustering, kmeans, pca, distances).
                #'
                #' The single matrix-analysis verb for the course. It absorbed runMatrixAnalyses()
                #' (plural) on 2026-09-22; the plural is retired and now stops with a migration message.
                #'
                #' @param data The data frame or tibble to analyze. May be WIDE (one column per analyte,
                #'   named in columns_w_values_for_single_analyte) or LONG (one row per sample per
                #'   analyte, named in column_w_names_of_multiple_analytes plus
                #'   column_w_values_for_multiple_analytes).
                #' @param analysis The type of analysis to run: one of "pca", "pca_ord", "pca_dim",
                #'   "mca", "mca_ord", "mca_dim", "mds", "mds_ord", "mds_dim", "tsne", "dbscan",
                #'   "kmeans", "hclust", "hclust_phylo", "hclust_cat" or "dist".
                #' @param parameters Analysis-specific parameters. For kmeans, the number of clusters;
                #'   for dbscan, c(eps, MinPts). Required for both.
                #' @param column_w_names_of_multiple_analytes Long input: the column holding analyte names.
                #' @param column_w_values_for_multiple_analytes Long input: the column holding the values.
                #' @param columns_w_values_for_single_analyte Wide input: the columns holding the values.
                #' @param columns_w_additional_analyte_info Per-analyte metadata to carry through.
                #' @param columns_w_sample_ID_info The columns that together identify a sample.
                #' @param transpose Analyse analytes as rows rather than samples. Default FALSE.
                #' @param distance_method What "similar" means. "euclidean" (default) and "manhattan"
                #'   (absolute differences; less pulled by one extreme variable); "gower" (mixed numeric
                #'   and categorical columns, via cluster::daisy); "correlation" (1 - Pearson r between
                #'   samples: same pattern, regardless of size) and "cosine" (1 - cosine similarity: same
                #'   direction, regardless of length); "bray_curtis" (abundance data; shared absences do
                #'   not count as similarity), "jaccard" (presence/absence: any value > 0 is present) and
                #'   "aitchison" (compositional data: euclidean on centred log-ratios; zeros are replaced
                #'   with half the smallest positive value, with a message). bray_curtis, jaccard and
                #'   aitchison need raw non-negative values, so pass scale_variance = FALSE with them.
                #'   All are computed in base R (no vegan), so every one also runs in WebR.
                #' @param agglomeration_method The linkage rule for hierarchical clustering. Only used
                #'   when tree_method = "linkage_dendrogram".
                #' @param tree_method How a tree is built for "hclust"/"hclust_phylo".
                #'   "linkage_dendrogram" (the DEFAULT) is hierarchical clustering proper: it honours
                #'   agglomeration_method and returns bootstrap support. "neighbor_joining" is the
                #'   phylogenetic method, works from the distance matrix alone, IGNORES
                #'   agglomeration_method, and returns no bootstrap values.
                #' @param unknown_sample_ID_info Sample IDs to treat as unknowns for nearest-neighbour
                #'   matching against the rest (used by the GC-MS workflows).
                #' @param components_to_return How many components to return. Default 2.
                #' @param scale_variance Center each column and divide it by its standard deviation
                #'   before computing distances. Default TRUE for every analysis (including "hclust",
                #'   which defaulted to FALSE until 2026-09-28). Pass FALSE for compositional data --
                #'   where every column is a share of the same whole and is already on one scale.
                #' @param na_replacement "mean" (default), "none", "zero" or "drop".
                #' @param output_format "wide" (default) or "long". For analysis = "dist" this chooses
                #'   the RETURN TYPE: "wide" gives a base R `dist` object, "long" gives one row per pair
                #'   (sample_1, sample_2, distance) plus the sample ID columns, which is the shape an
                #'   edge list wants.
                #' @param ... Not yet documented.
                #' @export
                runMatrixAnalysis <-    function(
                                            data,
                                            analysis = c(
                                                "pca", "pca_ord", "pca_dim",
                                                "mca", "mca_ord", "mca_dim",
                                                "mds", "mds_ord", "mds_dim",
                                                "tsne", "dbscan", "kmeans",
                                                "hclust", "hclust_phylo", "hclust_cat", "dist"
                                            ),
                                            parameters = NULL,
                                            column_w_names_of_multiple_analytes = NULL,
                                            column_w_values_for_multiple_analytes = NULL,
                                            columns_w_values_for_single_analyte = NULL,
                                            columns_w_additional_analyte_info = NULL,
                                            columns_w_sample_ID_info = NULL,
                                            transpose = FALSE,
                                            distance_method = c(
                                                "euclidean", "manhattan", "gower",
                                                "correlation", "cosine",
                                                "bray_curtis", "jaccard", "aitchison"
                                            ),
                                            agglomeration_method = c(
                                                "ward.D2", "ward.D", "single", "complete",
                                                "average", # (= UPGMA)
                                                "mcquitty", # (= WPGMA)
                                                "median", # (= WPGMC)
                                                "centroid" # (= UPGMC)
                                            ),
                                            ## NOTE the order: "linkage_dendrogram" is FIRST and so is the
                                            ## default. That is deliberate and differs from the retired
                                            ## runMatrixAnalyses(), which defaulted to neighbour joining.
                                            ## Every existing caller of runMatrixAnalysis -- the escape
                                            ## rooms, the GC-MS functions, the book -- was written against
                                            ## the linkage/Bclust behaviour this function has always had,
                                            ## so defaulting to NJ would have silently changed trees that
                                            ## puzzles already grade against. Nothing live relied on the
                                            ## plural's NJ default (its hclust callers all passed
                                            ## tree_method explicitly).
                                            tree_method = c(
                                                "linkage_dendrogram", # hclust -> as.phylo; honours agglomeration_method; bootstrapped
                                                "neighbor_joining"    # ape::nj; additive, non-ultrametric; IGNORES agglomeration_method
                                            ),
                                            unknown_sample_ID_info = NULL,
                                            components_to_return = 2,
                                            scale_variance = NULL, ## default = TRUE for every analysis; pass FALSE for compositional data
                                            na_replacement = c("mean", "none", "zero", "drop"),
                                            output_format = c("wide", "long"),
                                            ...
                                        ) {

                    # Check that argument names are spelled correctly

                        passed_args <- names(c(as.list(environment()), list(...)))

                        if (!all(passed_args %in%
                            c(
                                "data",
                                "analysis",
                                "parameters",
                                "column_w_names_of_multiple_analytes",
                                "column_w_values_for_multiple_analytes",
                                "columns_w_values_for_single_analyte",
                                "columns_w_additional_analyte_info",
                                "columns_w_sample_ID_info",
                                "transpose",
                                "distance_method",
                                "agglomeration_method",
                                "tree_method",
                                "unknown_sample_ID_info",
                                "components_to_return",
                                "scale_variance",
                                "na_replacement",
                                "output_format",
                                "..."
                            ))
                        ) {stop("One of your argument names is misspelled, please double check spelling.")}

                    # distance_method must be one we can compute

                        ## Checked up front: an unknown name used to reach stats::dist(), which
                        ## either errors deep inside a branch or (for its own extra methods such
                        ## as "maximum") quietly computes something the course never teaches.
                        .rma_distance_methods <- c(
                            "euclidean", "manhattan", "gower", "correlation", "cosine",
                            "bray_curtis", "jaccard", "aitchison"
                        )
                        if (!distance_method[1] %in% .rma_distance_methods) {
                            stop(paste0(
                                "distance_method = \"", distance_method[1], "\" is not available. ",
                                "Choose one of: ", paste0("\"", .rma_distance_methods, "\"", collapse = ", "), "."
                            ))
                        }

                    # parameters = "elbow" is a k-means-only request

                        ## Caught HERE, at the top, rather than in the k-means branch, because every
                        ## other analysis returns long before that branch is reached -- so without
                        ## this check, asking a pca for an elbow would quietly hand back a pca. The
                        ## branch itself (and the reasoning for the whole feature) is further down,
                        ## under "## k-means".
                        if (
                            length(parameters) > 0 &&
                            identical(tolower(as.character(parameters[1])), "elbow") &&
                            analysis[1] != "kmeans"
                        ) {
                            stop(paste0(
                                "parameters = \"elbow\" only makes sense with analysis = \"kmeans\". ",
                                "You asked for analysis = \"", analysis[1], "\"."
                            ))
                        }

                    # Check that column names are spelled correctly

                        if( any(
                            !c(
                                column_w_names_of_multiple_analytes,
                                column_w_values_for_multiple_analytes,
                                columns_w_values_for_single_analyte,
                                columns_w_additional_analyte_info,
                                columns_w_sample_ID_info
                            ) %in% colnames(data)
                            ) == TRUE
                        ) {
                            stop("There is a mismatch in the column names delivered to the command and the column names in your data. Please double check the spelling of your column names you gave to the command.")
                        }

                        if (all(c(
                            is.null(column_w_values_for_multiple_analytes),
                            is.null(columns_w_values_for_single_analyte)
                        ))) { stop("You need to specify at least one column with values for analytes.")}

                    # Pre-process data

                        # Add analyte_unique_ID_column if necessary

                            # if( length(columns_for_analyte_unique_ID) > 1 ) {
                                #add analyte_unique_ID column if necessary
                            # }

                        # Remove columns that are not included in input column lists

                            if (length(
                                    which(!colnames(data) %in% 
                                        c(
                                            column_w_names_of_multiple_analytes,
                                            column_w_values_for_multiple_analytes,
                                            columns_w_values_for_single_analyte,
                                            columns_w_additional_analyte_info,
                                            columns_w_sample_ID_info
                                        )
                                    )
                                ) > 0 
                            ) {
                                data <- data[,-which(!colnames(data) %in% 
                                    c(
                                        column_w_names_of_multiple_analytes,
                                        column_w_values_for_multiple_analytes,
                                        columns_w_values_for_single_analyte,
                                        columns_w_additional_analyte_info,
                                        columns_w_sample_ID_info
                                    )
                                )]
                            }

                        # Check for duplicate analyte names

                            if( length(columns_w_values_for_single_analyte) > 0 ) {
                                if( any(duplicated(columns_w_values_for_single_analyte)) ) {
                                    stop("There are duplicate analyte names in columns_w_values_for_single_analyte.")
                                }
                            }

                            if( length(column_w_names_of_multiple_analytes) == 1 ) {
                                n_unique_samples <- table(duplicated(data[,colnames(data) %in% columns_w_sample_ID_info]))[1]

                                analyte_breakdown <- table(duplicated(data[,colnames(data) == column_w_names_of_multiple_analytes]))
                                n_unique_groups_of_analytes <- sum(analyte_breakdown)/analyte_breakdown[1]

                                if (n_unique_samples != n_unique_groups_of_analytes) {
                                    stop("Either:\n1. The variables you have specified in columns_w_sample_ID_info do not define groups of samples with unique sets of analytes and you probably need to add more variables to columns_w_sample_ID_info.\n\nOR\n\n2. Not all your samples have the same number of analytes in column_w_names_of_multiple_analytes/column_w_values_of_multiple_analytes. Please fix that before continuing.")
                                }
                            }

                        # Remove analyte annotation columns before pivoting

                            if( length(columns_w_additional_analyte_info) > 0 ) {
                                analyte_annotation_free_data <- data[,-match(columns_w_additional_analyte_info, colnames(data))]
                            } else {
                                analyte_annotation_free_data <- data
                            }

                        # If no pivot required, skip pivoting

                            if( length(column_w_names_of_multiple_analytes) == 0 & length(columns_w_values_for_single_analyte) >= 1 ) {
                                data_wide <- analyte_annotation_free_data
                                analyte_columns <- columns_w_values_for_single_analyte
                                data_wide <- unique(data_wide)
                            }

                        # If pivoting required, pivot_wider any long-style data

                            if( length(column_w_names_of_multiple_analytes) == 1 ) {
                                data_wide <- pivot_wider(
                                    analyte_annotation_free_data,
                                    names_from = all_of(column_w_names_of_multiple_analytes),
                                    values_from = all_of(column_w_values_for_multiple_analytes)
                                )
                                analyte_columns <- unlist(unique(analyte_annotation_free_data[,colnames(analyte_annotation_free_data) == column_w_names_of_multiple_analytes]))
                                analyte_columns <- c(columns_w_values_for_single_analyte, analyte_columns)
                            }

                        # Check to see if analyte columns are numeric and compatible with analysis

                            which_analyte_columns <- which(colnames(data_wide) %in% analyte_columns)
                            are_they_numeric <- list()
                            for( i in which_analyte_columns ) {
                              are_they_numeric <- c(are_they_numeric, is.numeric(data_wide[[i]]))
                            }

                            # Should selected analysis proceed?

                                if ( all(unlist(are_they_numeric)) ) {
                                    if ( analysis %in% c("pca", "pca_dim", "pca_ord") ) {
                                        # cat("Analytes are all numeric and compatible with the analysis selected.\n")
                                    }
                                    if ( analysis %in% c("mca", "mca_ord", "mca_dim") ) {
                                        stop("Analytes are all numeric, but the analysis selected is for categorical variables. Please choose a different analysis method.\n")
                                    }
                                }

                                if ( !all(unlist(are_they_numeric)) ) {
                                    if (analysis %in% c("mca", "mca_ord", "mca_dim")) {
                                        # cat("Analytes are all categorical and compatible with the analysis selected.\n")
                                    }
                                    if ( analysis %in% c("pca", "pca_dim", "pca_ord") ) {
                                        stop("Analytes are all categorical, but the analysis selected is for numeric variables. Please choose a different analysis method.\n")
                                    }

                                }

                        # Add sample_unique_ID_column if necessary, or just change column name of existing sample_unique_ID column

                            if( length(columns_w_sample_ID_info) > 1 ) {
                                sample_unique_IDs <- apply(
                                    data_wide[,match(columns_w_sample_ID_info, colnames(data_wide))],
                                    1, paste, collapse = "_"
                                )
                                if( any(duplicated(sample_unique_IDs)) ) {stop("columns_w_sample_ID_info specified do not lead to unique sample IDs")}
                                data_wide$sample_unique_ID <- sample_unique_IDs
                            } else {
                                colnames(data_wide)[colnames(data_wide) == columns_w_sample_ID_info] <- "sample_unique_ID"
                                if( any(duplicated(data_wide$sample_unique_ID)) ) {stop("columns_w_sample_ID_info specified do not lead to unique sample IDs")}
                            }

                            # Make sure "sample_unique_ID" is character
                                data_wide$sample_unique_ID <- as.character(data_wide$sample_unique_ID)

                        # Prepare the matrix

                            matrix <- as.data.frame(data_wide[,match(analyte_columns, colnames(data_wide))])
                            rownames(matrix) <- data_wide$sample_unique_ID

                        # Handle NAs

                            if( na_replacement[1] == "none") {
                            }

                            if( na_replacement[1] == "drop" ) {
                                message("Dropping any variables in your dataset that have NA as a value.\nVariables dropped:\n")
                                if (length(names(which(apply(is.na(matrix), 2, any)))) > 0) {
                                    message(names(which(apply(is.na(matrix), 2, any))))    
                                } else {
                                    message("none")
                                }
                                cat("\n")
                                matrix <- matrix[,!apply(is.na(matrix), 2, any)]
                            }
                            if( na_replacement[1] %in% c("zero", "mean") ) {

                                if( any(is.na(matrix)) ) {
                                    
                                    message(paste0("Replacing NAs in your data with ", na_replacement[1]), "\n")

                                        for( column in 1:dim(matrix)[2]) {
                                            
                                            if( any(is.na(matrix[,column])) ) {

                                                if( na_replacement[1] == "mean" ) {
                                                    replacement <- mean(matrix[,column], na.rm = TRUE)
                                                }
                                                if( na_replacement[1] == "zero" ) {
                                                    replacement <- 0
                                                }
                                                if( !any(na_replacement %in% c("mean", "zero")) ) {
                                                    stop("Your data contains NAs. Please specify how to deal with them using na_replacement. \n")
                                                }
                                                
                                                matrix[,column][is.na(matrix[,column])] <- as.numeric(replacement)

                                            } else {}
                                        }
                                }
                            }

                        # Transpose matrix, if requested
                            if( transpose == TRUE ) { matrix <- t(matrix) }

                        # Run unknown, if requested

                            if( length(unknown_sample_ID_info) > 0 ) {
                            
                                ## Use na_replacement != "drop"

                                    if( na_replacement[1] != "drop" ) {
                                      stop("It is highly recommended that you use na_replacement = \"drop\" when matching an unknown, since anything with an NA value will be sent to the bottom of the matching list.")
                                    }

                                ## Separate unknown and knowns matrix

                                    sample_ID_of_unknown <- unknown_sample_ID_info
                                    index_of_unknown <- which(rownames(matrix) == sample_ID_of_unknown)
                                    unknown <- matrix[index_of_unknown,]
                                    matrix_minus_unknown <- matrix[-c(index_of_unknown),]

                                ## Identify how many knowns there are to test

                                    indices_to_test <- seq(1, dim(matrix_minus_unknown)[1], 1)

                                ## Loop and select 20 nearest neighbors until there are just 100 left

                                    while( length(indices_to_test) > 100 ) { ## Has been optimized for 50, 100, and 1000 indeces bites. 100 and 1000 both take ~5s.

                                        ## Select 100 entries at random, subset matrix_minus_unknown for that, append the unknowns
                                        
                                            indices_being_tested <- sample(indices_to_test, 100, replace = FALSE)

                                            matrix_to_be_tested <- matrix_minus_unknown[indices_being_tested,]
                                            matrix_to_be_tested <- rbind(unknown, matrix_to_be_tested)
                                            dist <- Rfast::Dist(matrix_to_be_tested)
                                            results <- data.frame(
                                              distance = dist[2:dim(dist)[1], 1],
                                              index = indices_being_tested
                                            )

                                        ## Keep closest 20 indices, modify indices_to_test

                                            indices_to_keep <- results[order(results$distance),]$index[1:20]
                                            indices_to_toss <- indices_being_tested[!indices_being_tested %in% indices_to_keep]
                                            indices_to_test <- indices_to_test[!indices_to_test %in% indices_to_toss]

                                    }

                                    ## Test the last 100 neighbors and return the closest 10

                                        final_matrix_to_be_tested <- rbind(unknown, matrix_minus_unknown[indices_to_test,])
                                        dist <- Rfast::Dist(final_matrix_to_be_tested)
                                        results <- data.frame(
                                          distance = dist[2:dim(dist)[1], 1],
                                          index = indices_to_test
                                    )

                                    if( length(indices_to_test) < 10 ) {
                                        number_for_final_tree <- length(indices_to_test)
                                    } else {
                                        number_for_final_tree <- 10
                                    }

                                    indices_to_keep <- results[order(results$distance),]$index[1:number_for_final_tree]
                                    indices_to_keep

                                ## Subset the matrix

                                    matrix <- rbind(unknown, matrix_minus_unknown[indices_to_keep,])

                                    ### THINK CAREFULLY HERE! THINGS SOMEWHAT CLOSE TO THE UNKNOWN 
                                    ### MIGHT NOT CLUSTER DIRECTLY WITH IT
                                    ### DEPENDING ON HOW RELATED THEY ARE TO OTHER THINGS

                            }

                    # Run the matrix analysis selected

                        # Scale data, unless not requested

                            ## Default TRUE for EVERY analysis, hclust included (2026-09-28).
                            ## hclust used to default to FALSE, alone among the analyses. That was a
                            ## silent footgun: unscaled euclidean distance is dominated by whichever
                            ## column is reported in the largest numbers, so the Alaska lakes tree
                            ## was really a chloride tree and looked entirely reasonable. Worse, it
                            ## meant pca and hclust on the SAME table used different geometry for no
                            ## visible reason. Book ch.8 already passed scale_variance = TRUE on
                            ## every call precisely because the old default defeated the lesson.
                            ## Callers that genuinely want unscaled distances -- compositional data,
                            ## where every column is a percentage of the same whole and scaling just
                            ## inflates the near-constant minor components -- must now pass
                            ## scale_variance = FALSE explicitly.
                            if ( is.null(scale_variance) ) {
                                scale_variance <- TRUE
                            }

                            ## Gower range-normalizes every column itself, so z-scaling first changes
                            ## nothing -- and scale() cannot take the categorical columns Gower is
                            ## for. So Gower always gets the matrix as it is (2026-10-27).
                            if( scale_variance == TRUE & !analysis %in% c("mca", "mca_ord", "mca_dim") & distance_method[1] != "gower") {
                                
                                scaled_matrix <- scale(matrix)

                                if( any(is.na(scaled_matrix)) ) {
                                    message("Some analytes have zero variance and will be assigned a value of zero in the scaled matrix.")
                                    scaled_matrix[is.na(scaled_matrix)] <- 0
                                }

                            }

                            if( scale_variance == FALSE | distance_method[1] == "gower" ) {
                                scaled_matrix <- matrix
                            }

                        ## FORTIFY, WITHOUT THE LIBRARY CHATTER (2026-09-28).
                        ## Newer tidytree (the version WebR installs; desktop 0.4.6 does NOT do
                        ## this, which is why it never showed up in testing here) emits
                        ##   "Invaild edge matrix for <phylo>. A <tbl_df> is returned."
                        ## -- its own typo -- from as_tibble.phylo when its edge-matrix validity
                        ## check is unhappy. Lucas hit it mid-exercise in canyon, TWICE per call,
                        ## on top of a perfectly correct answer.
                        ## It is harmless HERE and only here: everything below reads plain
                        ## columns (label, x, y, isTip), so tbl_df vs tbl_tree makes no
                        ## difference, and the browser smoke test proves ggtree() still plots the
                        ## frame. So muffle THAT message only -- matched on its text -- rather
                        ## than wrapping the call in a blanket suppressMessages(), which would
                        ## also swallow a genuine complaint from somewhere else in the chain.
                        .rma_quiet_fortify <- function(tr) {
                            withCallingHandlers(
                                ggtree::fortify(tr),
                                message = function(m) {
                                    if (grepl("edge matrix", conditionMessage(m), fixed = TRUE))
                                        invokeRestart("muffleMessage")
                                },
                                warning = function(w) {
                                    if (grepl("edge matrix", conditionMessage(w), fixed = TRUE))
                                        invokeRestart("muffleWarning")
                                }
                            )
                        }

                        ## ONE DISTANCE FUNCTION FOR EVERY BRANCH (2026-10-27).
                        ## The dist output, neighbour joining and the bootstrapped linkage tree
                        ## all call this, so a metric means the same thing wherever it is asked
                        ## for. Base R only: vegan has no guaranteed WebR build, and these are a
                        ## few lines each. Rows are samples, columns are variables.
                        .rma_distance <- function(m, method) {
                            if (method == "gower") {
                                ## Gower exists for MIXED columns, but daisy() rejects character
                                ## (and treats logical as asymmetric binary): hand it factors.
                                m <- as.data.frame(m)
                                for (j in seq_along(m)) {
                                    if (is.character(m[[j]]) || is.logical(m[[j]])) m[[j]] <- factor(m[[j]])
                                }
                                return(cluster::daisy(m, metric = "gower"))
                            }
                            x <- as.matrix(m)
                            if (method %in% c("bray_curtis", "jaccard", "aitchison") && any(x < 0, na.rm = TRUE)) {
                                stop(paste0(
                                    "distance_method = \"", method, "\" needs raw, non-negative values, but this ",
                                    "matrix has negative ones. Scaling (scale_variance = TRUE, the default) ",
                                    "centres every column on zero, so pass scale_variance = FALSE with this metric."
                                ))
                            }
                            d <- switch(method,
                                euclidean = stats::dist(x, method = "euclidean"),
                                manhattan = stats::dist(x, method = "manhattan"),
                                correlation = stats::as.dist(1 - stats::cor(t(x))),
                                cosine = {
                                    len <- sqrt(rowSums(x^2))
                                    stats::as.dist(1 - tcrossprod(x) / outer(len, len))
                                },
                                bray_curtis = {
                                    num <- as.matrix(stats::dist(x, method = "manhattan"))
                                    den <- outer(rowSums(x), rowSums(x), "+")
                                    bc <- num / den
                                    bc[den == 0] <- 0  # two all-zero samples: identical
                                    stats::as.dist(bc)
                                },
                                jaccard = stats::dist(x > 0, method = "binary"),
                                aitchison = {
                                    if (any(x == 0, na.rm = TRUE)) {
                                        pseudo <- min(x[x > 0], na.rm = TRUE) / 2
                                        message("Aitchison distance cannot take the log of zero; zeros were replaced with ",
                                                signif(pseudo, 3), " (half the smallest positive value).")
                                        x[x == 0] <- pseudo
                                    }
                                    lx <- log(x)
                                    stats::dist(lx - rowMeans(lx), method = "euclidean")
                                }
                            )
                            if (any(is.na(d))) {
                                message("Some distances are NA. With \"correlation\" or \"cosine\" this means a sample whose ",
                                        "values are all identical (or all zero), which has no pattern to compare.")
                            }
                            d
                        }

                        ## HCLUST, HCLUST_PHYLO ##

                            if( analysis == "hclust" | analysis == "hclust_phylo") {

                                ## Two ways to turn the data into a tree. `linkage_dendrogram` is
                                ## hierarchical clustering proper -- the step-by-step merging the book
                                ## teaches -- and it is the only one that honours agglomeration_method
                                ## and that yields bootstrap support. `neighbor_joining` is the
                                ## phylogenetic method; it works from the distance matrix alone and
                                ## IGNORES agglomeration_method entirely.
                                if (tree_method[1] == "linkage_dendrogram") {

                                    ## BClust approach to bootstrapped hclust

                                        bclust <- Bclust(
                                            scaled_matrix, method.d = distance_method[1],
                                            method.c = agglomeration_method[1],
                                            ## Bclust's own default FUN calls stats::dist, which
                                            ## knows none of correlation/cosine/bray_curtis/
                                            ## jaccard/aitchison (or gower). Every bootstrap
                                            ## replicate goes through the shared helper instead.
                                            FUN = function(.x) hclust(
                                                .rma_distance(.x, distance_method[1]),
                                                method = agglomeration_method[1]
                                            ),
                                            monitor = FALSE
                                        )
                                        phylo <- ape::as.phylo(bclust$hclust)

                                } else if (tree_method[1] == "neighbor_joining") {

                                    ## The distance matrix is computed further down for the other
                                    ## analyses, so neighbour joining makes its own here.
                                    nj_dist_matrix <- .rma_distance(scaled_matrix, distance_method[1])
                                    phylo <- ape::nj(nj_dist_matrix)

                                } else {
                                    stop("tree_method must be one of \"neighbor_joining\" or \"linkage_dendrogram\".")
                                }

                                if( analysis == "hclust_phylo" ) {
                                    return(phylo)
                                    stop("Returning hclust_phylo.")
                                }
                                clustering <- .rma_quiet_fortify(phylo)
                                clustering$sample_unique_ID <- clustering$label
                                clustering$bootstrap <- NA

                                ## Add bootstrap values starting from the furthest node to the highest
                                ## node. Only the linkage route has them -- neighbour joining is not
                                ## bootstrapped here, so its bootstrap column stays NA.
                                if (tree_method[1] == "linkage_dendrogram") {
                                    bs_vals <- data.frame(
                                        xval = clustering$x[clustering$isTip != TRUE],
                                        bs_val = NA
                                    )
                                    for (i in 1:length(bclust$value)) { # i=1
                                        bs_vals$bs_val[
                                            order(bs_vals$xval, decreasing = TRUE)[i]
                                        ] <- bclust$values[i]
                                    }
                                    clustering$bootstrap[clustering$isTip != TRUE] <- bs_vals$bs_val
                                }

                                ## CUTTING THE TREE (2026-09-23). `parameters` cuts the dendrogram
                                ## into groups and returns them as a `cluster` column, deliberately the
                                ## SAME shape k-means already returns -- one column, values "cluster_1",
                                ## "cluster_2", ... -- so that chapters 8 and 10 rhyme instead of
                                ## diverging. Without it there is no course-native way to cut a tree at
                                ## all: the hclust branch returns ggtree coordinates, not an hclust
                                ## object, so stats::cutree() on the output fails with "invalid 'tree'
                                ## ('merge' component)" and the student has to drop to base R -- which
                                ## is exactly the move the book tells them not to make.
                                ##
                                ## TWO WAYS TO ASK, because they are different questions:
                                ##   parameters = c(n_clusters = 4) -- cut into 4 groups, named
                                ##                                     explicitly (preferred -- symmetric
                                ##                                     with the height form below)
                                ##   parameters = c(4)              -- the same, bare (kept working: the
                                ##                                     escape rooms' graded hints use it,
                                ##                                     and k-means already reads a bare
                                ##                                     number this way)
                                ##   parameters = c(k = 4)          -- the same, the older explicit name
                                ##                                     (also kept working)
                                ##   parameters = c(height = 5)     -- cut ACROSS the tree at height 5,
                                ##                                     and however many groups that
                                ##                                     leaves is the answer
                                ## The height form is the one the canyon scenario's flood-gate escape
                                ## re-poses on the world (water level = cut height), and it is the only
                                ## form that can answer "how many groups are there at this height" --
                                ## with a count the number is an input, so "there is no fourth group"
                                ## cannot be expressed.
                                ##
                                ## Internal nodes get NA, like `bootstrap` does for tips: only tips are
                                ## samples, and only samples belong to a group.
                                if ( length(parameters) > 0 ) {

                                    if (tree_method[1] != "linkage_dendrogram") {
                                        stop(paste0(
                                            "Cutting a tree into groups needs tree_method = ",
                                            "\"linkage_dendrogram\". Neighbour joining does not build ",
                                            "the tree by merging groups, so there are no merge heights ",
                                            "to cut at."
                                        ))
                                    }

                                    param_name <- names(parameters)[1]
                                    if (is.null(param_name)) { param_name <- "" }
                                    if (is.na(param_name)) { param_name <- "" }

                                    if (param_name %in% c("h", "height")) {

                                        cut_height <- as.numeric(parameters[1])
                                        if (is.na(cut_height) || cut_height <= 0) {
                                            stop("parameters = c(height = ...) needs a positive height to cut the tree at.")
                                        }
                                        cut_groups <- stats::cutree(bclust$hclust, h = cut_height)

                                    } else if (param_name %in% c("", "k", "n_clusters")) {

                                        n_groups <- as.integer(parameters[1])
                                        if (is.na(n_groups) || n_groups < 1) {
                                            stop("parameters must be the number of groups to cut the tree into, e.g. parameters = c(n_clusters = 3) (or c(3) / c(k = 3)), or a height, e.g. parameters = c(height = 5).")
                                        }
                                        if (n_groups > nrow(scaled_matrix)) {
                                            stop(paste0(
                                                "Cannot cut ", nrow(scaled_matrix), " samples into ",
                                                n_groups, " groups."
                                            ))
                                        }
                                        cut_groups <- stats::cutree(bclust$hclust, k = n_groups)

                                    } else {
                                        stop(paste0(
                                            "parameters was named \"", param_name, "\". For a tree it must be ",
                                            "either a number of groups -- parameters = c(n_clusters = 3) ",
                                            "(or c(3) / c(k = 3)) -- or a height to cut at -- ",
                                            "parameters = c(height = 5)."
                                        ))
                                    }

                                    clustering$cluster <- NA_character_
                                    hit <- match(clustering$label, names(cut_groups))
                                    clustering$cluster[!is.na(hit)] <-
                                        paste0("cluster_", cut_groups[hit[!is.na(hit)]])
                                }
                            }

                            if (analysis == "hclust_cat") {

                                ## "bootstrap" approach
                                    temp <- as.data.frame(lapply(scaled_matrix, function(x) if(is.character(x)) factor(x) else x))
                                    rownames(temp) <- rownames(scaled_matrix)
                                    scaled_matrix <- temp
                                    createHclustObject <- function(x)hclust(cluster::daisy(x, metric = distance_method[1]), method = agglomeration_method[1])
                                    b <- bootstrap(scaled_matrix, fun = createHclustObject, n = 100L)
                                    phylo <- ape::as.phylo(createHclustObject(scaled_matrix))

                                    clustering <- .rma_quiet_fortify(phylo)
                                    clustering$sample_unique_ID <- clustering$label
                                    clustering$bootstrap <- NA

                                ## Add bootstrap values starting from the furthest node to the highest node
                                    bs_vals <- data.frame( xval = clustering$x[clustering$isTip != TRUE], bs_val = NA )
                                    for (i in 1:length(b)) { bs_vals$bs_val[order(bs_vals$xval, decreasing = TRUE)[i]] <- b[i] }

                                clustering$bootstrap[clustering$isTip != TRUE] <- bs_vals$bs_val

                            }

                            # if( analysis == "hclust" | analysis == "hclust_phylo" ) {
                            #     scaled_matrix <- as.data.frame(lapply(scaled_matrix, function(x) if(is.character(x)) factor(x) else x))
                            #     createHclustObject <- function(x)hclust(cluster::daisy(x, metric = distance_method[1]), method = agglomeration_method[1])
                            #     b <- bootstrap(scaled_matrix, fun = createHclustObject, n = 100L)
                            #     phylo <- ape::as.phylo(createHclustObject(scaled_matrix))
                            #     if( analysis == "hclust_phylo" ) {
                            #         return(phylo)
                            #         stop("Returning hclust_phylo.")
                            #     }
                            #     clustering <- ggtree::fortify(phylo)
                            #     clustering$sample_unique_ID <- clustering$label
                            #     clustering$bootstrap <- NA

                            #     ## Add bootstrap values starting from the furthest node to the highest node
                            #         bs_vals <- data.frame(
                            #             xval = clustering$x[clustering$isTip != TRUE],
                            #             bs_val = NA
                            #         )
                            #         for (i in 1:length(b)) { # i=1
                            #             bs_vals$bs_val[
                            #                 order(bs_vals$xval, decreasing = TRUE)[i]
                            #             ] <- b[i]
                            #         }

                            #     clustering$bootstrap[clustering$isTip != TRUE] <- bs_vals$bs_val
                            # }

                        # Generate distance matrix

                            if(  !analysis %in% c("mca", "mca_ord", "mca_dim") ) {

                                ## Through the shared helper (defined above the hclust branch), so
                                ## every metric -- gower via cluster::daisy, the rest in base R --
                                ## is computed the same way here as in the trees.
                                dist_matrix <- .rma_distance(scaled_matrix, distance_method[1])

                                if( analysis == "dist") {

                                    ## Two shapes, chosen with output_format. The default "wide" gives
                                    ## the `dist` object base R expects (hclust(), cmdscale(), print()).
                                    ## "long" gives one row per PAIR -- sample_1, sample_2, distance,
                                    ## plus the sample ID columns suffixed for each side -- which is the
                                    ## shape an edge list wants, and what the networks chapter plots.
                                    if( output_format[1] == "long" ) {
                                        df <- as.data.frame(as.table(as.matrix(dist_matrix)))
                                        names(df) <- c("sample_1", "sample_2", "distance")
                                        ## Drop each sample paired with ITSELF, not every zero:
                                        ## under jaccard, correlation or cosine two DIFFERENT
                                        ## samples are routinely at distance 0, and filtering on
                                        ## the value silently deleted those pairs (2026-10-27).
                                        df <- filter(df, as.character(sample_1) != as.character(sample_2))
                                        clustering <- df

                                        ## WITH A SINGLE columns_w_sample_ID_info, THAT COLUMN NO LONGER EXISTS.
                                        ## The pre-processing above (~line 14316) RENAMES a lone ID column to
                                        ## "sample_unique_ID" rather than adding one, so asking for
                                        ## c("lake", "sample_unique_ID") matched NA and the subset blew up with an
                                        ## error naming neither the argument nor the reason. Two or more ID columns
                                        ## survive as themselves and a new sample_unique_ID is pasted together, which
                                        ## is why the 2-column call worked and the 1-column call did not (found in
                                        ## WebR, 2026-09-23). Take only the ID columns that actually survived.
                                        dist_id_cols <- unique(c(
                                            intersect(columns_w_sample_ID_info, colnames(data_wide)),
                                            "sample_unique_ID"
                                        ))
                                        dist_id_frame <- data_wide[, match(dist_id_cols, colnames(data_wide)), drop = FALSE]

                                        clustering <- right_join(
                                          dist_id_frame, clustering,
                                          by = c("sample_unique_ID" = "sample_1"), suffix = c("_sample_1", "_sample_2")
                                        )

                                        clustering <- right_join(
                                          dist_id_frame, clustering,
                                          by = c("sample_unique_ID" = "sample_2"), suffix = c("_sample_1", "_sample_2")
                                        )

                                        colnames(clustering)[colnames(clustering) == "sample_unique_ID"] <- "sample_unique_ID_sample_1"
                                        clustering <- select(clustering, sample_unique_ID_sample_1, sample_unique_ID_sample_2, distance, everything())
                                        return(clustering)
                                    }

                                    return(dist_matrix)
                                }

                            }

                            if( analysis %in% c("mca", "mca_ord", "mca_dim") ) { scaled_matrix <- matrix }

                        ## Dimensionality reduction

                            ## MDS

                                if( analysis == "mds" ) {
                                    coords <- stats::cmdscale(dist_matrix)
                                    colnames(coords) <- c("Dim_1", "Dim_2")
                                    clustering <- as_tibble(coords)
                                    clustering$sample_unique_ID <- rownames(coords)
                                }
                        
                            ## tSNE

                                if( analysis == "tsne" ) {
                                    clustering <- Rtsne(scaled_matrix, theta = 0.0, perplexity = 2)
                                    clustering <- as.data.frame(clustering$Y)
                                    clustering$sample_unique_ID <- rownames(scaled_matrix)
                                    colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                                    clustering <- select(clustering, sample_unique_ID, Dim_1, Dim_2)
                                    rownames(clustering) <- NULL
                                }

                            ## umap

                                if( analysis == "umap" ) {
                                    clustering <- as.data.frame(umap(scaled_matrix)$layout)
                                    clustering$sample_unique_ID <- rownames(clustering)
                                    colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                                    clustering <- select(clustering, sample_unique_ID, Dim_1, Dim_2)
                                    rownames(clustering) <- NULL
                                }

                            ## MCA, MCA_ORD, MCA_DIM ##

                                if( analysis == "mca" ) {
                                    message("Running Multiple Correspondence Analysis, extracting sample coordinates...\n")
                                    coords <- FactoMineR::MCA(matrix, graph = FALSE)$ind$coord[,c(1:components_to_return)]
                                    clustering <- as_tibble(coords)
                                    clustering$sample_unique_ID <- rownames(coords)
                                    colnames(clustering) <- c("Dim_1", "Dim_2", "sample_unique_ID")
                                    message("Done!\n")
                                }

                                if( analysis == "mca_ord" ) {
                                    message("Running Multiple Correspondence Analysis, extracting ordination plot...\n")
                                    coords <- FactoMineR::MCA(matrix, graph = FALSE)$var$eta2[,c(1,components_to_return)]
                                    clustering <- as_tibble(coords)
                                    clustering$analyte <- rownames(coords)
                                    colnames(clustering) <- c("Dim_1", "Dim_2", "analyte")
                                    clustering <- select(clustering, analyte, Dim_1, Dim_2)
                                    return(clustering)
                                    stop("Returning ordination plot coordinates. \nDone!")
                                }

                                if( analysis == "mca_dim" ) {
                                    message("Running Multiple Correspondence Analysis, extracting dimensional contributions...\n")
                                    coords <- FactoMineR::MCA(matrix, graph = FALSE)$eig[,2]
                                    clustering <- tibble::enframe(coords, name = NULL)
                                    clustering$principal_component <- names(coords)
                                    clustering$principal_component <- as.numeric(gsub("dim ", "", clustering$principal_component))
                                    colnames(clustering)[colnames(clustering) == "value"] <- "percent_variance_explained"
                                    clustering <- select(clustering, principal_component, percent_variance_explained)
                                    return(clustering)
                                    stop("Returning eigenvalues. \nDone!")
                                }

                            ## PCA, PCA_ORD, PCA_DIM ## 

                                if( analysis == "pca" ) {
                                    coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$ind$coord[,c(1:components_to_return)]
                                    clustering <- as_tibble(coords)
                                    ## Dim_1, Dim_2, ... to match mds/tsne/umap/mca (FactoMineR names them Dim.1)
                                    colnames(clustering) <- sub("^Dim\\.", "Dim_", colnames(clustering))
                                    clustering$sample_unique_ID <- rownames(coords)
                                }

                                if( analysis == "pca_ord" ) {
                                    coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$var$coord[,c(1:components_to_return)]
                                    clustering <- as_tibble(coords)
                                    colnames(clustering) <- sub("^Dim\\.", "Dim_", colnames(clustering))
                                    clustering$analyte <- rownames(coords)
                                    clustering <- select(clustering, analyte, c(paste0("Dim_", seq(1,components_to_return,1))))
                                    return(clustering)
                                    stop("Returning ordination plot coordinates.")
                                }

                                if( analysis == "pca_dim" ) {
                                    coords <- FactoMineR::PCA(scaled_matrix, graph = FALSE, scale.unit = FALSE)$eig[,2]
                                    clustering <- tibble::enframe(coords, name = NULL)
                                    clustering$principal_component <- names(coords)
                                    clustering$principal_component <- as.numeric(gsub("comp ", "", clustering$principal_component))
                                    colnames(clustering)[colnames(clustering) == "value"] <- "percent_variance_explained"
                                    clustering <- select(clustering, principal_component, percent_variance_explained)
                                    return(clustering)
                                    stop("Returning eigenvalues.")
                                }
                    
                        ## Clustering

                            if(  !analysis %in% c("mca", "mca_ord", "mca_dim") ) {

                                if( any(is.na(scaled_matrix)) == TRUE ) {
                                    stop("clustering cannot handle NA. Please choose an option for na_replacement.")
                                }

                            }

                            ## DBSCAN

                                if( analysis == "dbscan" ) {

                                    if ( length(parameters) > 0 ) {
                                        cluster_k <- parameters[1]
                                        cluster_threshold <- parameters[2]
                                    }

                                    if ( length(parameters) == 0 ) {
                                        findClusterParameters(dist_matrix = dist_matrix, matrix = matrix, analysis = "dbscan")
                                    }

                                    message("Using", cluster_k, "as a value for k.\n")
                                    message("Using", cluster_threshold, "as a value for threshold.\n")
                                    clustering <- as_tibble(data.frame(
                                        sample_unique_ID = colnames(as.matrix(dist_matrix)),
                                        cluster = paste0("cluster_", fpc::dbscan(dist_matrix, eps = as.numeric(cluster_threshold), MinPts = as.numeric(cluster_k), scale = FALSE, method = "dist")[[1]])
                                    ))
                                    clustering$cluster[clustering$cluster == "cluster_0"] <- NA

                                }

                            ## k-means

                                if (analysis == "kmeans") {

                                    ## parameters = "elbow" -- return the ELBOW DATA, not a clustering.
                                    ##
                                    ## WHY THIS EXISTS (2026-09-30). Choosing k by the elbow is a taught
                                    ## step -- book ch.10 explains it with a figure and then says to find
                                    ## the bend -- but the only tool that ever DREW the elbow was the
                                    ## Shiny picker inside findClusterParameters(), which needs desktop R.
                                    ## So on every WebR surface (the escape rooms, the book's cells, the
                                    ## sandbox) the course taught the idea and then withheld the
                                    ## instrument, leaving a student to hand-roll
                                    ##     sapply(1:6, function(k) kmeans(m, k, nstart = 25)$tot.withinss)
                                    ## which runMatrixAnalysis() has no way to express and the book never
                                    ## shows. This returns exactly that table, ready to pipe to ggplot:
                                    ##
                                    ##     runMatrixAnalysis(..., analysis = "kmeans", parameters = "elbow") %>%
                                    ##       ggplot(aes(x = k, y = within_cluster_variance)) + geom_line() + geom_point()
                                    ##
                                    ## SHAPE. One row per candidate k; columns `k` and
                                    ## `within_cluster_variance` (the total within-cluster sum of squares
                                    ## -- the quantity ch.10 calls the within-group variance, and the same
                                    ## quantity the retired picker plotted). It is a TERMINAL return: an
                                    ## elbow table has no samples in it, so none of the annotation/join
                                    ## machinery below applies. `pca_dim` returns early for the same reason.
                                    ##
                                    ## RANGE. 1..10 by default, capped at nrow - 1 (at k = n every point
                                    ## is its own cluster, the spread is exactly zero, and the row carries
                                    ## no information). Override the maximum with a second element:
                                    ## parameters = c("elbow", 6).
                                    ##
                                    ## MATRIX. Deliberately the SAME `scaled_matrix` the clustering call
                                    ## below uses -- otherwise the elbow a student reads would belong to
                                    ## a different geometry than the clusters they then ask for, which is
                                    ## the worst possible failure here because it is invisible. Both
                                    ## honour `scale_variance` (default TRUE, as for every analysis); if
                                    ## the clustering call is ever moved to another matrix, MOVE THIS
                                    ## WITH IT. (Until 2026-10-04 both used the UNSCALED `matrix` and
                                    ## silently ignored scale_variance, though the docs said the default
                                    ## applied to every analysis.)
                                    if (
                                        length(parameters) > 0 &&
                                        identical(tolower(as.character(parameters[1])), "elbow")
                                    ) {

                                        if ( length(parameters) > 1 ) {
                                            max_k <- suppressWarnings(as.integer(parameters[2]))
                                            if (is.na(max_k) || max_k < 2) {
                                                stop(paste0(
                                                    "parameters = c(\"elbow\", n): n is the largest number of ",
                                                    "clusters to try and must be a whole number of at least 2."
                                                ))
                                            }
                                        } else {
                                            max_k <- 10L
                                        }

                                        max_k <- min(max_k, nrow(matrix) - 1L)

                                        if (is.na(max_k) || max_k < 2) {
                                            stop(paste0(
                                                "An elbow plot needs at least three samples to be worth ",
                                                "drawing; this data has ", nrow(matrix), "."
                                            ))
                                        }

                                        ## `tibble::as_tibble`, QUALIFIED -- unlike the bare `as_tibble`
                                        ## twenty lines below. The surrounding code assumes dplyr is
                                        ## attached, which is true on the desktop toolkit and true in a
                                        ## room whose scenario `setup` attaches it. This branch is the
                                        ## one a student may reach from a bare WebR sandbox or a book
                                        ## cell with no setup, and "could not find function as_tibble"
                                        ## in a timed room is not a failure worth risking for the sake
                                        ## of matching the line below. tibble is installed on every
                                        ## surface (BASE_PACKAGES).
                                        elbow <- tibble::as_tibble(data.frame(
                                            k = seq_len(max_k),
                                            within_cluster_variance = vapply(
                                                seq_len(max_k),
                                                function(i) stats::kmeans(
                                                    x = scaled_matrix, centers = i,
                                                    nstart = 25, iter.max = 1000
                                                )$tot.withinss,
                                                numeric(1)
                                            )
                                        ))

                                        return(elbow)
                                        stop("Returning elbow data.")

                                    }

                                    if ( length(parameters) > 0 ) {
                                        n_clusters <- parameters[1]
                                    }

                                    if ( length(parameters) == 0 ) {
                                        findClusterParameters(dist_matrix = dist_matrix, matrix = matrix, analysis = "kmeans")
                                    }

                                    message("Using", n_clusters, "as a value for cluster_number.\n")
                                    clustering <- as_tibble(data.frame(
                                        sample_unique_ID = colnames(as.matrix(dist_matrix)),
                                        cluster = paste0("cluster_", stats::kmeans(x = scaled_matrix, centers = as.numeric(n_clusters), nstart = 25, iter.max = 1000)$cluster)
                                    ))

                                }

                            ## OPTICS

                                # out <- dbscan::optics(scaled_matrix, minPts = 5)
                                # out <- data.frame(
                                #     order = out$order,
                                #     reach_dist = out$reachdist,
                                #     name = rownames(scaled_matrix)
                                # )
                                # out$name <- factor(out$name, levels = rev(rownames(scaled_matrix)[out$order]))
                                # ggplot(out[2:19,]) +
                                #     geom_col(aes(x = name, y = reach_dist))                                

                    # If transpose = TRUE, then return without annotations

                        if( transpose == TRUE ) {
                            return(clustering)
                            stop("Returning transposed cluster output. Make sure all your variables have the same units!")
                        }
                        
                    # Post processing and return.

                        if( !analysis %in% c("hclust_phylo")) {

                            ## Add back annotations to the output

                                ## The empty branch is LOAD-BEARING, not dead code. With exactly one
                                ## columns_w_sample_ID_info the pre-processing above RENAMED that column to
                                ## sample_unique_ID, so it is no longer in data_wide and there is nothing to
                                ## join back -- the annotation is already there. Do not "tidy" this into an
                                ## unconditional join: that is the same mistake the dist/long branch made
                                ## (~line 14605), which threw an NA-subscript error naming neither the
                                ## argument nor the reason. Regression test: test_runMatrixAnalysis_dist.R
                                if( length(columns_w_sample_ID_info) == 1 ) {
                                } else {
                                clustering <-   right_join(
                                                    data_wide[,match(
                                                        c(columns_w_sample_ID_info, "sample_unique_ID"),
                                                        colnames(data_wide))
                                                    ], clustering, by = "sample_unique_ID"
                                                )
                                }

                                rownames_matrix <- tibble::enframe(rownames(scaled_matrix), name = NULL)
                                colnames(rownames_matrix)[1] <- "sample_unique_ID"

                                # if (analysis != "pca") { ## Don't do this for pca for some reason?? I don't understand why...
                                    clustering <- full_join(
                                        clustering,
                                        as_tibble(cbind(rownames_matrix, as_tibble(matrix))),
                                        by = "sample_unique_ID"
                                    )
                                # }
                                # clustering

                                ## Order the returned matrix so that the sample_unique_ID comes first

                                    clustering <- select(clustering, sample_unique_ID, everything())

                            # Annotate internal nodes in tree output if all its descendants share a property

                                if( analysis == "hclust" ) {
                                    for( node in dplyr::filter(clustering, isTip == FALSE)$node ) {
                                        for (sample_property in colnames(clustering)[colnames(clustering) %in% columns_w_sample_ID_info] ) {
                                            descends <- clustering[clustering$node %in% ips::descendants(phylo, node),]
                                            if (length( unlist(unique(descends[,colnames(descends) == sample_property])) ) == 1 ) {
                                                clustering[
                                                    which(clustering$node == node),
                                                    which(colnames(clustering) == sample_property)
                                                ] <- unlist(descends[,colnames(descends) == sample_property])[1]
                                            }
                                        }
                                    }
                                }

                            # Return results

                                if( output_format[1] == "long" ) {
                                    clustering <- pivot_longer(
                                        clustering,
                                        cols = c(which(colnames(clustering) == analyte_columns[1]): dim(clustering)[2]),
                                        names_to = "analyte_name", 
                                        values_to = "value"
                                    )
                                    analyte_annotation_frame <- unique(select(ungroup(data), all_of(c(column_w_names_of_multiple_analytes, columns_w_additional_analyte_info))))
                                    clustering <- left_join(clustering, analyte_annotation_frame, by = c("analyte_name" = column_w_names_of_multiple_analytes))
                                }

                                return( clustering )
                        }
                }


        #### findClusterParameters


            #' Interactively find clustering parameters from a distance matrix and source matrix
            #'
            #' @param dist_matrix Not yet documented.
            #' @param matrix Not yet documented.
            #' @param analysis Not yet documented.
            #' @export
            findClusterParameters <- function(dist_matrix, matrix, analysis = c("dbscan", "kmeans")) {

                ui <-   navbarPage(title = "Cluster Parameter Selection", id = "navbar",
                            
                            tabPanel(title = "Main",

                                sidebarLayout(

                                    sidebarPanel(

                                        conditionalPanel( "input.selection == 'kmeans'",
                                        
                                            sliderInput(
                                                inputId = "n_clusters",
                                                label = "Number of clusters",
                                                min = 0, max = 25, value = 1, step = 1
                                            )
                                        
                                        ),
                                
                                        conditionalPanel( "input.selection == 'dbscan'",

                                            sliderInput(
                                                inputId = "cluster_threshold", label = "Neighborhood size (epsilon or 'eps')",
                                                min = 0, max = 50, value = 2, step = 0.05
                                            ),

                                            sliderInput(
                                                inputId = "cluster_k", label = "Number of neighbors required to form a cluster (MinPts)",
                                                min = 1, max = 40, value = 2
                                            )


                                        ),

                                        plotOutput(
                                            outputId = "plot1",
                                            width = 200, height = 200
                                        )
                                    ),

                                    mainPanel(
                                            plotOutput(
                                                outputId = "plot2", width = "100%"
                                                # height = 800
                                            )
                                    )
                                )

                            ),

                            tabPanel(title = "Quit", value = "stop", icon = icon("circle-notch"),
                            
                                selectInput(
                                    inputId = "selection",
                                    label = "Do not touch this",
                                    choices = c("dbscan", "kmeans"),
                                    selected = analysis
                                )

                            )

                )
                            
                picked <- NULL
                server <- function(input, output, session) {

                    output$plot1 <- renderPlot({

                        ## read by runMatrixAnalysis() after the app closes
                        assign("n_clusters", input$n_clusters, envir = globalenv())
                        assign("cluster_k", input$cluster_k, envir = globalenv())
                        assign("cluster_threshold", input$cluster_threshold, envir = globalenv())
                        picked <<- list(n_clusters = input$n_clusters, cluster_k = input$cluster_k,
                                        cluster_threshold = input$cluster_threshold)

                        if (analysis == "dbscan") {
                            distances <- dbscan::kNNdist(dist_matrix, k = input$cluster_k)
                        }

                        if (analysis == "kmeans") {
                            kmeans_results <- list()
                            for( i in 1:(dim(matrix)[1]-1) ) {
                                kmeans_results[[i]] <- sum(stats::kmeans(x = matrix, centers = i, nstart = 25, iter.max = 1000)$withinss)
                            }
                            distances <- do.call(rbind, kmeans_results)
                        }

                        distances <- as.numeric(distances)[rev(order(as.numeric(distances)))]

                        ggplot(
                            data.frame(
                                index = seq(1,length(distances),1),
                                distance = distances
                            ), aes(x = index, y = distance)
                        ) + geom_point() + theme_bw() +
                        geom_hline(yintercept = input$n_clusters)

                    })

                    output$plot2 <- renderPlot ({
                        if (analysis == "dbscan") {
                            clustering <- as_tibble(data.frame(
                                sample_unique_ID = colnames(as.matrix(dist_matrix)),
                                cluster = paste0("cluster_", fpc::dbscan(
                                    dist_matrix, 
                                    eps = input$cluster_threshold, 
                                    MinPts = input$cluster_k, 
                                    scale = FALSE, method = "dist"
                                )[[1]])
                            ))
                            clustering$cluster[clustering$cluster == "cluster_0"] <- NA
                        }

                        if (analysis == "kmeans") {
                            kmeans_clusters <- stats::kmeans(
                                x = matrix, centers = input$n_clusters, 
                                nstart = 25, iter.max = 1000
                            )$cluster
                            clustering <- as_tibble(data.frame(
                                sample_unique_ID = names(kmeans_clusters), 
                                cluster = paste0("cluster_", kmeans_clusters)
                            ))
                        }

                        clustering <- cbind(
                            clustering,
                            matrix[match(rownames(matrix), clustering$sample_unique_ID),]
                        )

                        plot <- ggplot(clustering, aes_string(
                            x = colnames(matrix)[1], 
                            y = colnames(matrix)[2], 
                            fill = "cluster"
                        ))

                        if (analysis == "dbscan") {
                            plot <- plot + # draw eps neighborhoods
                                ggforce::geom_circle(
                                    aes(
                                        x0 = !!sym(colnames(matrix)[1]), 
                                        y0 = !!sym(colnames(matrix)[2]), 
                                        r = input$cluster_threshold
                                    ),
                                    color = "black", inherit.aes = FALSE, alpha = 0.1
                                )
                        }

                        plot <- plot +
                            geom_point(shape = 21, size = 5, alpha = 0.6) +
                            theme_bw() +
                            scale_fill_manual(values = discrete_palette) +
                            coord_fixed()
                        plot
                    })

                    session$onSessionEnded(function() { stopApp() })

                    observe({ if (input$navbar == "stop") stopApp() })

                }

                runApp(shinyApp(ui, server))
                invisible(picked)

            }

