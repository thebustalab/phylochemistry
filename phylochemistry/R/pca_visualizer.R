## pcaVisualizer — moved verbatim from modules/pca_visualizer.R (library() lines dropped)

# Launch with pcaVisualizerApp(data_frame)
pcaVisualizer <- function(data,
                             columns_w_sample_ID_info = NULL,
                             columns_w_values_for_single_analyte = NULL) {
  
  stopifnot(is.data.frame(data))
  
  # sensible defaults if not supplied
  if (is.null(columns_w_sample_ID_info)) {
    # guess: keep non-numeric columns as sample-ID info
    columns_w_sample_ID_info <- names(data)[!sapply(data, is.numeric)]
  }
  if (is.null(columns_w_values_for_single_analyte)) {
    # guess: numeric columns are analytes
    columns_w_values_for_single_analyte <- names(data)[sapply(data, is.numeric)]
  }
  
  base_font_size <- 16
  ui <- fluidPage(
    tags$head(
      tags$style(HTML(sprintf(
        paste(
          "body { font-size: %dpx; }",
          ".shiny-input-container > label, .control-label { font-size: %dpx; }",
          ".radio label { font-size: %dpx; }",
          ".help-block { font-size: %dpx; }",
          "h4 { font-size: %dpx; }",
          sep = " "
        ),
        base_font_size,
        base_font_size,
        base_font_size,
        base_font_size,
        base_font_size + 4
      )))
    ),
    titlePanel("PCA Visualizer"),
    
    sidebarLayout(
      sidebarPanel(
        width = 2,
        # Color/shape aesthetics (must be among ID columns)
        selectInput(
          "color_var",
          "Color mapped to",
          choices = c("None" = "__none__", columns_w_sample_ID_info),
          selected = "__none__"
        ),
        selectInput(
          "shape_var",
          "Shape mapped to",
          choices = c("None" = "__none__", columns_w_sample_ID_info),
          selected = "__none__"
        ),
        
        tags$hr(),
        radioButtons(
          "heatmap_order",
          "Arrange Heat Map by:",
          choices = c(
            "By Dim_1" = "dim1",
            "By Dim_2" = "dim2"
          ),
          selected = "dim1"
        ),
        sliderInput(
          "loading_threshold",
          "Filter ordination plot",
          min = 0,
          max = 1,
          value = 0.8,
          step = 0.05
        )
      ),
      
      mainPanel(
        width = 10,
        fluidRow(
          column(
            width = 6,
            h4("PCA Plot"),
            helpText("also sometimes called 'scores plot'"),
            plotOutput("scores_plot", height = 520)
          ),
          column(
            width = 6,
            h4("Ordination Plot"),
            helpText("also sometimes called 'loadings plot' or 'eigenvalues'"),
            plotOutput("loadings_plot", height = 520)
          )
        ),
        br(),
        h4("Analyte Abundance Heatmap"),
        plotOutput("heatmap_plot", height = 420)
      )
    )
  )
  
  server <- function(input, output, session) {
    # Sample-ID columns are fixed to the provided defaults
    id_cols <- reactive({
      intersect(columns_w_sample_ID_info, names(data))
    })

    resolve_loading_labels <- function(ld_df) {
      if (!nrow(ld_df)) {
        return(character())
      }
      label_col <- intersect(names(ld_df), c("variable", "feature", "name", "label"))
      if (length(label_col)) {
        return(as.character(ld_df[[label_col[1]]]))
      }
      rn <- rownames(ld_df)
      if (!is.null(rn) && any(nzchar(rn))) {
        rn_clean <- rn
        empty_idx <- which(!nzchar(rn_clean))
        if (length(empty_idx)) {
          rn_clean[empty_idx] <- paste0("Var ", empty_idx)
        }
        return(rn_clean)
      }
      if (length(columns_w_values_for_single_analyte) >= nrow(ld_df)) {
        return(columns_w_values_for_single_analyte[seq_len(nrow(ld_df))])
      }
      paste0("Var ", seq_len(nrow(ld_df)))
    }

    filter_indices_by_loading <- function(ld_df, threshold) {
      if (!nrow(ld_df)) {
        return(integer())
      }
      dim_cols <- intersect(c("Dim_1", "Dim_2"), names(ld_df))
      if (length(dim_cols) < 2) {
        return(seq_len(nrow(ld_df)))
      }
      load_mat <- ld_df[, dim_cols[1:2], drop = FALSE]
      magnitudes <- sqrt(rowSums(load_mat^2, na.rm = TRUE))
      magnitudes[!is.finite(magnitudes)] <- 0
      which(magnitudes >= threshold)
    }
    
    # Reactive: run PCA (scores)
    pca_scores <- reactive({
      cols <- id_cols()
      shiny::validate(
        shiny::need(length(cols) >= 1, "Data must include at least one sample-ID column."),
        shiny::need(length(columns_w_values_for_single_analyte) >= 2, "Data must include at least two analyte columns.")
      )

      runMatrixAnalysis(
        data = data,
        analysis = "pca",
        columns_w_sample_ID_info = cols,
        columns_w_values_for_single_analyte = columns_w_values_for_single_analyte
      )
    })
    
    # Reactive: run PCA (loadings/ordination)
    pca_loadings <- reactive({
      cols <- id_cols()
      shiny::validate(shiny::need(length(cols) >= 1, "Data must include at least one sample-ID column."))

      runMatrixAnalysis(
        data = data,
        analysis = "pca_ord",
        columns_w_sample_ID_info = cols,
        columns_w_values_for_single_analyte = columns_w_values_for_single_analyte
      )
    })
    
    heatmap_data <- reactive({
      sample_cols <- id_cols()
      shiny::validate(shiny::need(length(sample_cols) >= 1, "Data must include at least one sample-ID column."))
      all_analyte_cols <- columns_w_values_for_single_analyte
      ld_all <- pca_loadings()
      labels_all <- resolve_loading_labels(ld_all)
      threshold <- if (is.null(input$loading_threshold)) 0 else input$loading_threshold
      keep_idx <- filter_indices_by_loading(ld_all, threshold)
      if (length(all_analyte_cols) >= max(c(keep_idx, 0))) {
        analyte_cols <- all_analyte_cols[keep_idx]
      } else {
        analyte_cols <- character(0)
      }
      if (!length(analyte_cols)) {
        possible <- labels_all[keep_idx]
        possible <- possible[possible %in% names(data)]
        analyte_cols <- unique(possible)
      }
      analyte_cols <- analyte_cols[analyte_cols %in% names(data)]
      if (!length(analyte_cols)) {
        return(NULL)
      }

      dim_cols <- intersect(c("Dim_1", "Dim_2"), names(ld_all))
      order_choice <- if (is.null(input$heatmap_order)) "original" else input$heatmap_order
      if (order_choice %in% c("dim1", "dim2")) {
        axis_col <- if (order_choice == "dim1") "Dim_1" else "Dim_2"
        if (axis_col %in% names(ld_all)) {
          axis_values <- ld_all[[axis_col]]
          axis_values <- suppressWarnings(as.numeric(axis_values))
          axis_values <- axis_values[keep_idx]
          axis_values[!is.finite(axis_values)] <- NA_real_
          analyte_order <- order(axis_values, decreasing = TRUE, na.last = TRUE)
          keep_idx <- keep_idx[analyte_order]
          analyte_cols <- analyte_cols[analyte_order]
        }
      }
      if (!length(analyte_cols)) {
        return(NULL)
      }
      selected <- data[, unique(c(sample_cols, analyte_cols)), drop = FALSE]
      sample_labels <- if (length(sample_cols) == 1) {
        as.character(selected[[sample_cols]])
      } else {
        apply(selected[sample_cols], 1, function(row) paste(row, collapse = " | "))
      }
      analyte_matrix <- data.matrix(selected[analyte_cols])
      scaled_matrix <- apply(analyte_matrix, 2, function(col) {
        center <- col - mean(col, na.rm = TRUE)
        max_abs <- suppressWarnings(max(abs(center), na.rm = TRUE))
        if (!is.finite(max_abs) || max_abs == 0) {
          rep(0, length(col))
        } else {
          pmax(pmin(center / max_abs, 1), -1)
        }
      })
      if (is.null(dim(scaled_matrix))) {
        scaled_matrix <- matrix(scaled_matrix, ncol = 1)
      }
      colnames(scaled_matrix) <- analyte_cols
      order_choice <- if (is.null(input$heatmap_order)) "original" else input$heatmap_order
      sample_levels <- unique(sample_labels)
      if (order_choice %in% c("dim1", "dim2")) {
        scores_df <- pca_scores()
        req(all(c("Dim_1", "Dim_2") %in% names(scores_df)))
        req(all(sample_cols %in% names(scores_df)))
        score_labels <- if (length(sample_cols) == 1) {
          as.character(scores_df[[sample_cols]])
        } else {
          apply(scores_df[sample_cols], 1, function(row) paste(row, collapse = " | "))
        }
        ordering_df <- data.frame(
          Sample = score_labels,
          Dim1 = scores_df[["Dim_1"]],
          Dim2 = scores_df[["Dim_2"]],
          stringsAsFactors = FALSE
        )
        ordering_df <- ordering_df[!is.na(ordering_df$Sample), , drop = FALSE]
        ordering_df <- ordering_df[!duplicated(ordering_df$Sample), , drop = FALSE]
        if (order_choice == "dim1") {
          ordering_df <- ordering_df[order(ordering_df$Dim1, ordering_df$Dim2, na.last = TRUE), , drop = FALSE]
        } else {
          ordering_df <- ordering_df[order(ordering_df$Dim2, ordering_df$Dim1, na.last = TRUE), , drop = FALSE]
        }
        sample_levels <- unique(c(ordering_df$Sample, sample_levels))
      }
      num_samples <- nrow(analyte_matrix)
      analyte_count <- length(analyte_cols)
      data.frame(
        Sample = factor(rep(sample_labels, each = analyte_count), levels = sample_levels),
        Analyte = factor(rep(analyte_cols, times = num_samples), levels = analyte_cols),
        Value = as.vector(t(scaled_matrix)),
        stringsAsFactors = FALSE
      )
    })
    
    # Scores plot (samples)
    output$scores_plot <- renderPlot({
      df <- pca_scores()
      req(all(c("Dim_1", "Dim_2") %in% names(df)))
      req(input$color_var, input$shape_var)

      color_var <- if (is.null(input$color_var)) "__none__" else input$color_var
      shape_var <- if (is.null(input$shape_var)) "__none__" else input$shape_var
      use_color <- !identical(color_var, "__none__")
      use_shape <- !identical(shape_var, "__none__")

      aes_args <- list(x = "Dim_1", y = "Dim_2")
      if (use_color) {
        aes_args$colour <- color_var
      }
      if (use_shape) {
        aes_args$shape <- shape_var
      }
      mapping <- do.call(aes_string, aes_args)

      point_args <- list(alpha = 0.85, size = 2.8)
      if (!use_color) {
        point_args$colour <- "black"
      }
      if (!use_shape) {
        point_args$shape <- 1
      }

      plot <- ggplot(df, mapping) +
        do.call(geom_point, point_args) +
        labs(x = "PC1 (Dim_1)", y = "PC2 (Dim_2)") +
        theme_bw(base_size = base_font_size)

      if (use_color) {
        plot <- plot +
          labs(color = color_var) +
          scale_color_manual(values = discrete_palette) ## don't mess with this line!
      }

      if (use_shape) {
        plot <- plot + labs(shape = shape_var)
      }

      plot
    })
    
    # Loadings plot (variable vectors)
    output$loadings_plot <- renderPlot({
      ld_all <- pca_loadings()
      req(all(c("Dim_1", "Dim_2") %in% names(ld_all)))
      labels_all <- resolve_loading_labels(ld_all)
      threshold <- if (is.null(input$loading_threshold)) 0 else input$loading_threshold
      keep_idx <- filter_indices_by_loading(ld_all, threshold)
      shiny::validate(shiny::need(length(keep_idx) > 0, "No analytes meet the loading threshold."))
      order_choice <- if (is.null(input$heatmap_order)) "original" else input$heatmap_order
      if (order_choice %in% c("dim1", "dim2")) {
        axis_col <- if (order_choice == "dim1") "Dim_1" else "Dim_2"
        if (axis_col %in% names(ld_all)) {
          axis_values <- ld_all[[axis_col]]
          axis_values <- suppressWarnings(as.numeric(axis_values))
          axis_values <- axis_values[keep_idx]
          axis_values[!is.finite(axis_values)] <- NA_real_
          reorder_idx <- order(axis_values, decreasing = TRUE, na.last = TRUE)
          keep_idx <- keep_idx[reorder_idx]
        }
      }
      ld <- ld_all[keep_idx, , drop = FALSE]
      lbl <- labels_all[keep_idx]
      lbl[is.na(lbl)] <- ""
      ld$`__label` <- lbl

      g <- ggplot(ld, aes(x = 0, y = 0, xend = .data[["Dim_1"]], yend = .data[["Dim_2"]])) +
        geom_segment(arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.5, color = "#636363") +
        coord_equal() +
        labs(x = "PC1 loading (Dim_1)", y = "PC2 loading (Dim_2)") +
        theme_bw(base_size = base_font_size)

      g + ggrepel::geom_text_repel(
          aes(
            x = .data[["Dim_1"]], y = .data[["Dim_2"]],
            label = paste(.data[["analyte"]], "(", round(.data[["Dim_1"]], 3), ", ", round(.data[["Dim_2"]], 3), ")")
          ),
          direction = "x",
          color = "#333333",
          size = 3,
          seed = 123,
          box.padding = 0.25,
          point.padding = 0.15,
          max.overlaps = Inf,
          min.segment.length = 0
        )
    })
    
    output$heatmap_plot <- renderPlot({
      hm <- heatmap_data()
      shiny::validate(shiny::need(!is.null(hm) && nrow(hm) > 0, "No analytes meet the loading threshold."))
      ggplot(hm) +
        geom_tile(aes(x = Analyte, y = Sample, fill = Value), color = NA) +
        # geom_text(aes(x = Analyte, y = Sample, label = round(Value, 2)), size = 3, color = "black") +
        scale_fill_gradient2(
          low = "#08519c",
          mid = "#ffffbf",
          high = "#cb181d",
          midpoint = 0,
          limits = c(-1, 1),
          na.value = "#f0f0f0",
          oob = scales::squish
        ) +
        labs(x = "Analyte", y = "Sample", fill = "Scaled Value (-1 to 1)") +
        theme_bw(base_size = base_font_size) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1)
        )
    })

  }
  
  shinyApp(ui, server)
}
