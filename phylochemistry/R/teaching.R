## teaching — moved verbatim from phylochemistry.R (package conversion, commit 1)


        #### question

            #' Create a question

            question <- function(question, distractors, correct, no, fb = "", print_question = TRUE){
                allanswers <- c(distractors, correct)[sample.int(length(distractors)+1)]
                correctanswer <- which(allanswers == correct)
                
                answercode <- paste0(sapply(1:length(allanswers), function(i){
                
                    x <- allanswers[i]
                    paste0('<div class="radio">\n  <label>\n    <input type="radio" name="question', no, '" id="opt', i,'" value="', i, '" onchange="check_answer', no, '()">\n    ', x, '\n  </label>\n</div>')
                
                }), collapse = "\n\n")

                out <- paste0(question, "\n\n", answercode,
                    '<div class="collapse" id="collapseExample', no,'">
                    <div class="card card-body" id="answerFeedback', no, '">
                    </div>
                    </div>',
                    paste0('<script type="text/javascript">
                    function check_answer', no, '()
                    {
                    var radioButtons', no, ' = document.getElementsByName("question', no, '");
                    document.getElementById("answerFeedback', no, '").innerHTML = "Try selecting an answer!!";
                    for(var i = 0; i < radioButtons', no, '.length; i++)
                    {
                    if(radioButtons', no, '[i].checked == true)
                    {
                    var feedback', no, ' = "<p style=\'color:red\'>Wrong', ifelse(fb == "", ".", paste0("; ", fb)),  '</p>";
                    if(radioButtons', no, '[i].value == "', correctanswer, '") {
                    feedback', no, ' = "<p style=\'color:green\'>Correct!</p>"
                    }
                    document.getElementById("answerFeedback', no, '").innerHTML = feedback', no, ';
                    return true;
                    }
                    }
                    }
                    </script>
                    '))
                
                if(print_question){
                    cat(out)
                } else {
                    return(out)
                }
            }


        #### questionnaire

            questionnaire <- function(x, shuffle = TRUE, print_question = TRUE){
                
                if(inherits(x, "character")) x <- read.csv(x, stringsAsFactors = FALSE, fileEncoding="UTF-8-BOM")
                
                if(!all(names(x) == c("question", "distractors", "correct", "fb"))){
                    stop("Incorrect column names")
                }
                
                if(shuffle){
                    x <- x[sample.int(nrow(x)), ]
                }
                
                out <- ""
                for(i in 1:nrow(x)){
                    out <- paste0(out,
                    "**Question ",
                    i,
                    ":**\n",
                    question(question = x$question[i],
                    distractors = eval(parse(text = x$distractors[i])),
                    correct = x$correct[i],
                    no = i,
                    fb = ifelse(is.na(x$fb[i]), "", x$fb[i]),
                    print_question = FALSE),
                    "\n\n"
                    )
                }
                if(print_question){
                    cat(out)
                } else {
                    return(out)
                }
            }


    ##### Language models and text data handling (embedText in language_model_analysis.R)

        #### buildCaption

            buildCaption <- function(jupyter = FALSE) {

                ui <- fluidPage(
                  theme = shinythemes::shinytheme("darkly"),
                  useShinyjs(),
                  sidebarLayout(
                    sidebarPanel(
                      h2("Caption"),
                      textAreaInput("caption_input", "Enter your figure caption:", "", rows = 5),
                      h2("Caption Element Checklist"),
                      
                      # TITLE SECTION
                      h3("TITLE SECTION"),
                      div(id = "elem_figure_title",
                          radioButtons("figure_title",
                                       label = "FIGURE TITLE: Provide a concise, descriptive title that summarizes the overall message or purpose of the figure, including the main topic, experimental system, or hypothesis.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      div(id = "elem_panel_id",
                          radioButtons("panel_id",
                                       label = "PANEL IDENTIFICATION: Label each panel (e.g., A, B, C) and refer to these labels consistently.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      actionButton("reset_title", "Reset Title Section"),
                      
                      br(), br(),
                      
                      # SUBPANEL SECTION
                      h3("SUBPANEL SECTION"),
                      div(id = "elem_graph_layout",
                          radioButtons("graph_layout",
                                       label = "GRAPH TYPE AND LAYOUT: State the type of plot (e.g., line plot, bar chart, scatter plot, histogram, etc.) and describe any special features such as insets, overlays, or embedded plots.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      div(id = "elem_axes_variables",
                          radioButtons("axes_variables",
                                       label = "AXES AND VARIABLE MAPPINGS: Clearly define what is on each axis (including units) and describe any additional dimensions (e.g., color coding, marker sizes, or symbols) used in the figure.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      div(id = "elem_data_representation",
                          radioButtons("data_representation",
                                       label = "DATA REPRESENTATIONS: Describe what the individual data points, bars, or error bars represent (e.g., individual measurements, means, medians, etc.), including any error metrics (such as standard error, standard deviation, or 95% confidence intervals) and details on any trend or regression lines.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      div(id = "elem_sample_size",
                          radioButtons("sample_size",
                                       label = "SAMPLE SIZE: Indicate the number of independent samples or experimental replicates underlying each element of the graph.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      div(id = "elem_stat_analysis",
                          radioButtons("stat_analysis",
                                       label = "STATISTICAL ANALYSES AND COMPARISONS: Describe any control experiments or baseline comparisons, specify the statistical tests used (e.g., t-test, ANOVA, regression analysis) and significance levels, and detail any data normalization, transformation, or curve fitting applied.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      actionButton("reset_subpanel", "Reset Subpanel Section"),
                      
                      br(), br(),
                      
                      # FINAL CHECKS SECTION
                      h3("FINAL CHECKS SECTION"),
                      div(id = "elem_caption_text",
                          radioButtons("caption_text",
                                       label = "AVOID OR DEFINE ABBREVIATIONS: Avoid using unexplained abbreviations or jargon; if abbreviations are necessary, provide definitions at first use.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      div(id = "elem_data_methodology",
                          radioButtons("data_methodology",
                                       label = "DATA METHODOLOGY: Outline the methods used to generate the data, including experimental conditions (e.g., treatment concentrations, temperature), computational parameters if applicable, and any image processing adjustments critical for interpreting the figure.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      div(id = "elem_data_source",
                          radioButtons("data_source",
                                       label = "DATA SOURCE: Indicate the origin of the data (e.g., experimental assays, clinical samples, simulations, databases, publications) and include relevant references if applicable.",
                                       choices = list("Doesn't Apply" = "na", "Complete" = "complete"),
                                       selected = character(0),
                                       inline = TRUE)
                      ),
                      actionButton("reset_final", "Reset Final Checks Section")
                    ),
                    
                    mainPanel(
                      h3("Image Input"),
                      textInput("img_path", label = "Enter image path (PDF, PNG, JPG)"),
                      actionButton("read_btn", label = "Read Image"),
                      hr(),
                      uiOutput("figureImage")
                    )
                  )
                )

                server <- function(input, output, session) {
                  
                  # Reactive image path updated on 'Read Image' button click.
                  imagePath <- eventReactive(input$read_btn, {
                    input$img_path
                  })
                  
                  # Render the image or error message if needed.
                  output$figureImage <- renderUI({
                    req(imagePath())
                    
                    if (!file.exists(imagePath())) {
                      return(tags$div(style = "color:red;", 
                                      sprintf("Error: File not found. Please check the path: %s", imagePath())))
                    }
                    
                    ext <- tolower(tools::file_ext(imagePath()))
                    
                    if (ext == "pdf") {
                      # Use iframe for PDFs.
                      tags$iframe(style = "height:600px; width:100%; border: none;", src = imagePath())
                    } else if (ext %in% c("png", "jpg", "jpeg", "gif")) {
                      mime_type <- switch(ext,
                                          "png" = "image/png",
                                          "jpg" = "image/jpeg",
                                          "jpeg" = "image/jpeg",
                                          "gif" = "image/gif",
                                          "application/octet-stream")
                      
                      img_data <- tryCatch({
                        dataURI(file = imagePath(), mime = mime_type)
                      }, error = function(e) {
                        NULL
                      })
                      
                      if (is.null(img_data)) {
                        return(tags$div(style = "color:red;", "Error: Unable to read or encode the image file."))
                      }
                      
                      tags$img(src = img_data, style = "max-width:100%; max-height:600px; height:auto;")
                    } else {
                      tags$div(style = "color:red;", 
                               "Error: Unrecognized file type. Please use PDF, PNG, JPG, JPEG, or GIF.")
                    }
                  })
                  
                  # Hide checklist items when a selection is made.
                  observeEvent(input$figure_title, {
                    if (!is.null(input$figure_title) && input$figure_title != "")
                      hide("elem_figure_title")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$caption_text, {
                    if (!is.null(input$caption_text) && input$caption_text != "")
                      hide("elem_caption_text")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$panel_id, {
                    if (!is.null(input$panel_id) && input$panel_id != "")
                      hide("elem_panel_id")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$graph_layout, {
                    if (!is.null(input$graph_layout) && input$graph_layout != "")
                      hide("elem_graph_layout")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$axes_variables, {
                    if (!is.null(input$axes_variables) && input$axes_variables != "")
                      hide("elem_axes_variables")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$data_representation, {
                    if (!is.null(input$data_representation) && input$data_representation != "")
                      hide("elem_data_representation")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$sample_size, {
                    if (!is.null(input$sample_size) && input$sample_size != "")
                      hide("elem_sample_size")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$stat_analysis, {
                    if (!is.null(input$stat_analysis) && input$stat_analysis != "")
                      hide("elem_stat_analysis")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$data_methodology, {
                    if (!is.null(input$data_methodology) && input$data_methodology != "")
                      hide("elem_data_methodology")
                  }, ignoreInit = TRUE)
                  
                  observeEvent(input$data_source, {
                    if (!is.null(input$data_source) && input$data_source != "")
                      hide("elem_data_source")
                  }, ignoreInit = TRUE)
                  
                  # Reset Title Section
                  observeEvent(input$reset_title, {
                    updateRadioButtons(session, "figure_title", selected = character(0))
                    updateRadioButtons(session, "panel_id", selected = character(0))
                    shinyjs::show("elem_figure_title")
                    shinyjs::show("elem_panel_id")
                  })
                  
                  # Reset Subpanel Section
                  observeEvent(input$reset_subpanel, {
                    updateRadioButtons(session, "graph_layout", selected = character(0))
                    updateRadioButtons(session, "axes_variables", selected = character(0))
                    updateRadioButtons(session, "data_representation", selected = character(0))
                    updateRadioButtons(session, "sample_size", selected = character(0))
                    updateRadioButtons(session, "stat_analysis", selected = character(0))
                    shinyjs::show("elem_graph_layout")
                    shinyjs::show("elem_axes_variables")
                    shinyjs::show("elem_data_representation")
                    shinyjs::show("elem_sample_size")
                    shinyjs::show("elem_stat_analysis")
                  })
                  
                  # Reset Final Checks Section
                  observeEvent(input$reset_final, {
                    updateRadioButtons(session, "caption_text", selected = character(0))
                    updateRadioButtons(session, "data_methodology", selected = character(0))
                    updateRadioButtons(session, "data_source", selected = character(0))
                    shinyjs::show("elem_caption_text")
                    shinyjs::show("elem_data_methodology")
                    shinyjs::show("elem_data_source")
                  })
                  
                }

                if ( jupyter == TRUE) {

                    slot <- pick_shiny_port()
                    if (is.null(slot)) {
                        cat("All available ports are in use. Please try again later.")
                    } else {
                        cat(paste0("Connect at: ", slot$url))
                        runApp(shinyApp(ui = ui, server = server), host = "127.0.0.1", port = slot$port)
                    }

                }

                if (jupyter == FALSE) { shinyApp(ui = ui, server = server) }
            }

