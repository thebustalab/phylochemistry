## viz — moved verbatim from phylochemistry.R (package conversion, commit 1)


        #### rotate_coord

            #' rotates a set of x,y coordinates
            #'
            #' @param x The x values
            #' @param y The y values
            #' @param angle The angle to rotate (in degrees)
            #' @param center The center of rotation (default = 0,0)
            #' @examples
            #' @export
            #' rotate_coord

            rotate_coord <- function( x, y, angle, center = c(0,0) ){
  
                x <- x - center[1]
                y <- y - center[2]
                angle <- angle*pi/180

                conversionmatrix <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)), ncol=2, nrow=2)
                xy <- cbind(x,y)%*%conversionmatrix

                xy[,1] <- xy[,1]+center[1]
                xy[,2] <- xy[,2]+center[2]

                return(xy)

            }


        #### coord_radar

            coord_radar <- function (theta = "x", start = 0, direction = 1) 
                {
                    theta <- match.arg(theta, c("x", "y"))
                    r <- if (theta == "x") 
                        "y"
                    else "x"
                    ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
                        direction = sign(direction),
                        is_linear = function(coord) TRUE)
                }


    ##### Data Visualization
        
        #### -.gg

            #' Insert a geom layer below others
            #'
            #' Allows the user to plot a multiple alignent using ggplot's grammar of graphics
            #' @param plot The plot to which the layer should be added
            #' @param layer The layer to be added
            #' @export
            #' @examples
            #' -.gg()

            `-.gg` <- function(plot, layer) {
                if (missing(layer)) {
                    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
                }
                if (!is.ggplot(plot)) {
                    stop('Need a plot on the left side')
                }
                plot$layers = c(layer, plot$layers)
                plot
            }


        #### generateTicks

            #' Create major and minor axes ticks and labels
            #'
            #' @param major_ticks Values at which major ticks should be generated
            #' @param minor_freq Number of minor ticks between each major tick 
            #' @examples
            #' @export
            #' generateTicks

            generateTicks <-    function(
                                    major_ticks,
                                    minor_freq = 4,
                                    major_tick_size = 2,
                                    minor_tick_size = 1
                                ) {
                                    
                major_labels <- vector()
                
                for (tick in 1:(length(major_ticks)-1)) {
                    major_labels <- c(major_labels, major_ticks[tick])
                    major_labels <- c(major_labels, rep("", (minor_freq)))
                }
                
                all_ticks <- vector()
                
                for (tick in 1:(length(major_ticks)-1)) {
                    all_ticks <- c(all_ticks, major_ticks[tick])
                    minor_tick_values <- vector()
                    for (minor_tick in 1:minor_freq) {
                        minor_tick_values <- c(minor_tick_values, major_ticks[tick] + minor_tick*(((major_ticks[tick+1] - major_ticks[tick])/(minor_freq+1))))
                    }
                    all_ticks <- c(all_ticks, minor_tick_values)
                }
                
                major_labels <- c(major_labels, major_ticks[length(major_ticks)])
                all_ticks <- c(all_ticks, major_ticks[length(major_ticks)])

                return <- data.frame(all_ticks = as.numeric(all_ticks), major_labels = as.character(major_labels))
                return$tick_size <- minor_tick_size
                return$tick_size[return$major_labels != ""] <- major_tick_size

                return(return)
            }


        #### drawMolecules

            #' Draw chemical structures from CSV
            #'
            #' @param csv_in_path Path to the csv with the molecular coordinates.
            #' @param data Data that the functions should plot.
            #' @export
            #' @examples
            #' drawMolecules()

            drawMolecules <-    function(
                                    csv_in_path = NULL, 
                                    data = NULL, 
                                    atom_color_column = "default",
                                    bond_color_column = "default",
                                    molecule_name = TRUE,
                                    atom_numbers = TRUE
                                ) {

                ## Check input data / path

                    if (is.null(csv_in_path) & is.null(data)) {
                        stop("Please provide a data source using either the `data` argument OR the `csv_in_path` argument.")
                    }

                    if (!is.null(csv_in_path) & !is.null(data)) {
                        stop("Please provide a SINGLE data source using either the `data` argument OR the `csv_in_path` argument.")
                    }

                    if (is.null(csv_in_path) & !is.null(data)) {}

                    if (!is.null(csv_in_path) & is.null(data)) {
                        data <- read_csv(csv_in_path, col_types = cols())
                    }

                ## Clean up data
                    output <- list()
                    for (i in 1:length(unique(data$molecule_name))) {

                        temp <- filter(data, molecule_name == unique(data$molecule_name)[i])
                        temp$bond_start_x <- temp$x[match(temp$bond_start_atom, temp$atom_number)]
                        temp$bond_start_y <- temp$y[match(temp$bond_start_atom, temp$atom_number)]
                        temp$bond_end_x <- temp$x[match(temp$bond_end_atom, temp$atom_number)]
                        temp$bond_end_y <- temp$y[match(temp$bond_end_atom, temp$atom_number)]
                        temp$bond_start_element_or_group <- temp$element_or_group[match(temp$bond_start_atom, temp$atom_number)]
                        temp$bond_end_element_or_group <- temp$element_or_group[match(temp$bond_end_atom, temp$atom_number)]
                        output[[i]] <- temp
                    
                    }

                    plot_data <- do.call(rbind, output)

                ## Define color schemes

                    atom_colors_pre <- as.data.frame(rbind(
                        c("C", "black"),
                        c("CH2", "grey"),
                        c("CH2OH", "#4daf4a"), # green
                        c("CH3", "#ff7f00"), # orange
                        c("COOH", "#ffff33"), # yellow
                        c("H", "white"),
                        c("O", "#e41a1c"), # red
                        c("OH", "#377eb8"), # blue
                        c("CHO", "grey20"),
                        c("OCOCH3", "grey30"),
                        c("COOCH3", "grey40"),
                        c("CH2CH3", "grey50")
                    ))
                    atom_colors <- atom_colors_pre[,2]
                    names(atom_colors) <- atom_colors_pre[,1]

                    #984ea3 purple
                    #a65628 brown

                ## Plotting

                    plot <- ggplot()

                        ## Add bonds under default color scheme

                            if (bond_color_column == "default") {

                                ## Add achiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" &  bond_type == "single" & bond_direction == "flat"
                                        ),
                                        aes(
                                            x = bond_start_x, y = bond_start_y, xend = bond_end_x, yend = bond_end_y
                                        ), size = 2, color = "grey30"
                                    )

                                ## Add down chiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" & bond_type == "single" & bond_direction == "down"
                                        ),
                                        aes(
                                            x = bond_start_x, y = bond_start_y, xend = bond_end_x, yend = bond_end_y,
                                            size = stat(index)
                                        ), color = "grey60"
                                    )
                                
                                ## Add up chiral single bonds
                                    plot <- plot + geom_link(
                                        data = filter(plot_data, 
                                            molecule_component == "bond" & bond_type == "single" & bond_direction == "up"
                                        ),
                                        aes(
                                            x = bond_start_x, y = bond_start_y, xend = bond_end_x, yend = bond_end_y,
                                            size = stat(index)
                                        ), color = "grey0"
                                    )

                                ## Add steep non-horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & bond_start_y != bond_end_y & abs(bond_start_y- bond_end_y) != 1),
                                        aes(x = bond_start_x-0.25, y = bond_start_y-0.15, xend = bond_end_x-0.25, yend = bond_end_y-0.15),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & bond_start_y != bond_end_y & abs(bond_start_y- bond_end_y) != 1),
                                        aes(x = bond_start_x+0.25, y = bond_start_y+0.15, xend = bond_end_x+0.25, yend = bond_end_y+0.15),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Add nearly horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & abs(bond_start_y- bond_end_y) == 1),
                                        aes(x = bond_start_x, y = bond_start_y-0.3, xend = bond_end_x, yend = bond_end_y-0.3),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x & abs(bond_start_y- bond_end_y) == 1),
                                        aes(x = bond_start_x, y = bond_start_y+0.3, xend = bond_end_x, yend = bond_end_y+0.3),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Add horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_y == bond_end_y),
                                        aes(x = bond_start_x, y = bond_start_y-0.2, xend = bond_end_x, yend = bond_end_y-0.2),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_y == bond_end_y),
                                        aes(x = bond_start_x, y = bond_start_y+0.2, xend = bond_end_x, yend = bond_end_y+0.2),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Add vertical double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes(x = bond_start_x-0.25, y = bond_start_y, xend = bond_end_x-0.25, yend = bond_end_y),
                                        size = 1.5, alpha = 0.7
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes(x = bond_start_x+0.25, y = bond_start_y, xend = bond_end_x+0.25, yend = bond_end_y),
                                        size = 1.5, alpha = 0.7
                                    )

                                ## Color scheme
                                    # plot <- plot + scale_color_gradient(low = "black", high = "black")
                           
                            }

                        ## Add bonds under custom color scheme

                            if (bond_color_column != "default") {

                                if (dim(unique(filter(plot_data, molecule_component == "bond") %>% select(bond_color_column)))[1] == 1 ) {
                                    ## Add achiral single bonds
                                        plot <- plot + geom_link(
                                            data = filter(plot_data, 
                                                molecule_component == "bond" &  bond_type == "single" & bond_direction == "flat"
                                            ),
                                            aes_string(
                                                x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y"
                                            ), size = 2, color = "#440154FF"
                                        )
                                    ## Add down chiral single bonds
                                        plot <- plot + geom_link(
                                            data = filter(plot_data, 
                                                molecule_component == "bond" & bond_type == "single" & bond_direction == "down"
                                            ),
                                            aes_string(
                                                x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                                size = "stat(index)"
                                            ), color = "#440154FF"
                                        )
                                    ## Add up chiral single bonds
                                        plot <- plot + geom_link(
                                            data = filter(plot_data, 
                                                molecule_component == "bond" & bond_type == "single" & bond_direction == "up"
                                            ),
                                            aes_string(
                                                x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                                size = "stat(index)"
                                            ), color = "#440154FF"
                                        )
                                } else {
                                    ## Add achiral single bonds
                                        plot <- plot + geom_link(
                                            data = filter(plot_data, 
                                                molecule_component == "bond" &  bond_type == "single" & bond_direction == "flat"
                                            ),
                                            aes_string(
                                                x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                                color = bond_color_column
                                            ), size = 2
                                        )
                                    ## Add down chiral single bonds
                                        plot <- plot + geom_link(
                                            data = filter(plot_data, 
                                                molecule_component == "bond" & bond_type == "single" & bond_direction == "down"
                                            ),
                                            aes_string(
                                                x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                                size = "stat(index)", color = bond_color_column
                                            )
                                        )
                                    ## Add up chiral single bonds
                                        plot <- plot + geom_link(
                                            data = filter(plot_data, 
                                                molecule_component == "bond" & bond_type == "single" & bond_direction == "up"
                                            ),
                                            aes_string(
                                                x = "bond_start_x", y = "bond_start_y", xend = "bond_end_x", yend = "bond_end_y",
                                                size = "stat(index)", color = bond_color_column
                                            )
                                        )
                                }

                                ## Add nearly horizontal double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x),
                                        aes_string(x = "bond_start_x", y = "bond_start_y-0.2", xend = "bond_end_x", yend = "bond_end_y-0.2"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x != bond_end_x),
                                        aes_string(x = "bond_start_x", y = "bond_start_y+0.2", xend = "bond_end_x", yend = "bond_end_y+0.2"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    )

                                ## Add vertical double bonds
                                    plot <- plot + geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes_string(x = "bond_start_x-0.18", y = "bond_start_y", xend = "bond_end_x-0.18", yend = "bond_end_y"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    ) +
                                    geom_segment(
                                        data = filter(plot_data, molecule_component == "bond" & bond_type == "double" & bond_start_x == bond_end_x),
                                        aes_string(x = "bond_start_x+0.18", y = "bond_start_y", xend = "bond_end_x+0.18", yend = "bond_end_y"),
                                        size = 1.5, alpha = 0.7, color = "blue"
                                    )
                            }

                        
                        ## Add atom circles 
                          
                            if (atom_color_column == "default") {
                                plot <- plot + geom_point(
                                    data = filter(plot_data, molecule_component == "atom"),
                                    aes(x = x, y = y, fill = element_or_group),
                                    shape = 21, size = 4
                                ) +
                                scale_fill_manual(values = atom_colors, name = "")
                            }

                            if (atom_color_column != "default") {
                                if (dim(unique(filter(plot_data, molecule_component == "atom") %>% select(atom_color_column)))[1] == 1) {
                                    plot <- plot + geom_point(
                                        data = filter(plot_data, molecule_component == "atom"),
                                        aes_string(x = "x", y = "y"), fill = "#440154FF",
                                        shape = 21, size = 4
                                    )
                                } else {
                                    plot <- plot + geom_point(
                                        data = filter(plot_data, molecule_component == "atom"),
                                        aes_string(x = "x", y = "y", fill = atom_color_column),
                                        shape = 21, size = 4
                                    ) +
                                    scale_fill_viridis(na.value = "#440154FF")
                                }
                            }

                        ## Add atom number labels

                            if (atom_numbers) {
                          
                                plot <- plot + geom_text(
                                    data = filter(plot_data, molecule_component == "atom", atom_number %in% c(10, 50, 71)),
                                    aes(x = x, y = y, label = atom_number),
                                    size = 2, color = "black"
                                ) +

                                geom_text(
                                    data = filter(plot_data, molecule_component == "atom", !atom_number %in% c(10, 50, 71)),
                                    aes(x = x, y = y, label = atom_number),
                                    size = 2, color = "white"
                                )

                            }

                        ## Add molecule name as label
                            if (molecule_name) {

                                plot <- plot + geom_text(
                                    data = drop_na(unique(select(plot_data, molecule_name))),
                                    aes(x = 1, y = 13, label = molecule_name),
                                    size = 4, color = "black", hjust = 0
                                )

                            }

                        ## Scales and theme
                            plot <- plot + 
                            scale_size_continuous(range = c(0,4)) +
                            scale_x_continuous(breaks = seq(0,20,1)) +
                            scale_color_viridis(na.value = "#440154FF") +
                            # scale_linetype_manual(values = bond_line_types) +
                            scale_y_continuous(breaks = seq(0,20,1)) +
                            facet_wrap(.~molecule_name, ncol = 3) +
                            # facet_wrap(.~molecule_name) +
                            theme_void() +
                            guides(size = "none", alpha = "none") +
                            coord_fixed() +
                            theme(
                                strip.text = element_blank()
                            )

                    return(plot)
            }           


        #### vennAnalysis

            #' Get coordinates for quantitative Venn Diagram
            #'
            #' @param data Data frame of TRUE FALSE for whether each observation belongs to each category
            #' @examples
            #' vennAnalysis()

            vennAnalysis <- function(data) {
                # str(plot(euler(data, shape = "circle"), quantities = TRUE))
                # plot(euler(data, shape = "circle"), quantities = TRUE)$data$centers
                venn_output <- eulerr::euler(data, shape = "circle")
                ellipses <- cbind(
                    data.frame(category = rownames(venn_output$ellipses)),
                    venn_output$ellipses[,c(1:3)]
                )
                colnames(ellipses) <- c("category", "x", "y", "r")
                
                overlaps <- plot(venn_output, quantities = TRUE)$data$centers[,c(1,2)]
                overlaps <- overlaps[(dim(venn_output$ellipses)[1]+1):(dim(overlaps)[1]),c(1,2)]
                overlaps <- cbind(
                    data.frame(category = rownames(overlaps)),
                    overlaps,
                    data.frame(r = NA)
                )
                venn_circle_data <- rbind(ellipses, overlaps)
                original_values <- as.data.frame(venn_output$original.values)
                original_values$category <- rownames(original_values)
                colnames(original_values)[1] <- "value"
                venn_circle_data$count <- original_values$value[match(venn_circle_data$category, original_values$category)]

                return(venn_circle_data)
            }


        #### points3D

            #' Get coordinates for a 3D scatter plot
            #'
            #' @param data Data frame with columns x, y, z, and sample_unique_ID
            #' @param angle Angle of the z axis
            #' @param tick_round Integer to which tick labels should be rounded
            #' @param x_tick_interval Interval of ticks on x axis
            #' @param y_tick_interval Interval of ticks on y axis
            #' @param z_tick_interval Interval of ticks on z axis
            #' @examples
            #' points3D()

            points3D <- function( data, angle, tick_round, x_tick_interval, y_tick_interval, z_tick_interval ) {

                ## Initial calculations

                    output <- list()

                    ymin = plyr::round_any(min(data$y), tick_round, f = floor)
                    ymax = plyr::round_any(max(data$y), tick_round, f = ceiling)
                    ylength = ymax - ymin

                    zmin = plyr::round_any(min(data$z), tick_round, f = floor)
                    zmax = plyr::round_any(max(data$z), tick_round, f = ceiling)
                    zlength = zmax - zmin

                    xmin = plyr::round_any(min(data$x), tick_round, f = floor)
                    xmax = plyr::round_any(max(data$x), tick_round, f = ceiling)
                    xend = xmax+(zlength*cos(angle))
                    xlength = xmax - xmin

                    xintervals = seq(xmin, xmax, x_tick_interval)
                    yintervals = seq(ymin, ymax, y_tick_interval)
                    zintervals = seq(zmin, zmax, z_tick_interval)

                ## Draw grid

                    grid_x <- data.frame(
                        y = ((zintervals-zmin) * sin(angle)),
                        yend = ((zintervals-zmin) * sin(angle)),
                        x = ((zintervals-zmin) * cos(angle)) + xmin,
                        xend = (((zintervals-zmin) * cos(angle)) + xmax)
                    )

                    grid_z <- data.frame(
                        y = 0,
                        yend = zlength*sin(angle),
                        x = xintervals,
                        xend = zlength*cos(angle) + xintervals
                    )

                    output$grid <- rbind(grid_z, grid_x)

                ## Draw ticks

                    ticks_x <- data.frame(
                        y = 0,
                        yend = zlength*0.02*sin(angle),
                        x = xintervals,
                        xend = zlength*0.02*cos(angle) + xintervals
                    )

                    ticks_x2 <- data.frame(
                        y = 0,
                        yend = -zlength*0.02*sin(angle),
                        x = xintervals,
                        xend = -zlength*0.02*cos(angle) + xintervals
                    )

                    ticks_y <- data.frame( ### HERE
                        y = yintervals-ymin,
                        yend = yintervals-ymin,
                        x = -0.02*xlength + xmin,
                        xend = 0.02*xlength + xmin
                    )
                    
                    ticks_z <- data.frame(
                        y = ((zintervals-zmin) * sin(angle)),
                        yend = ((zintervals-zmin) * sin(angle)),
                        x = (((zintervals-zmin) * cos(angle)) + xmax) - 0.02*xlength,
                        xend = (((zintervals-zmin) * cos(angle)) + xmax) + 0.02*xlength
                    )

                    output$ticks <- rbind(ticks_x, ticks_x2, ticks_y, ticks_z)

                ## Draw tick labels

                    labels_x <- data.frame(
                        y = -ylength*0.04,
                        x = xintervals-0.04*xlength,
                        label = xintervals
                    )

                    labels_y <- data.frame(
                        y = yintervals-ymin,
                        x = -xlength*0.08 + xmin,
                        label = yintervals
                    )

                    labels_z <- data.frame(
                        y = ((zintervals-zmin) * sin(angle)),
                        x = (((zintervals-zmin) * cos(angle)) + xmax) + xlength*0.08,
                        label = zintervals
                    )

                    output$labels <- rbind(labels_x, labels_y, labels_z)

                ## Draw axes

                    output$axes <- data.frame(
                        x = c(xmin, xmin, xmax),
                        xend = c(xmax, xmin, xend),
                        y = c(0, 0, 0),
                        yend = c(0, ylength, zlength*sin(angle))
                    )
                
                ## Draw points and their dashed lines

                    output$point_segments <- data.frame(
                        x = data$x + (data$z-zmin)*cos(angle),
                        xend = data$x + (data$z-zmin)*cos(angle),
                        y = data$y - ymin + (data$z-zmin)*sin(angle),
                        yend = (data$z-zmin)*sin(angle)
                    )

                    output$points <- data.frame(
                        x = data$x + (data$z-zmin)*cos(angle),
                        y = data$y - ymin + (data$z-zmin)*sin(angle),
                        sample_unique_ID = data$sample_unique_ID
                    )

                return(output)
            }

