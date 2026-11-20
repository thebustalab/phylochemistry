## literature — moved verbatim from phylochemistry.R (package conversion, commit 1)


        #### parseRISFile

                #' Parse an RIS reference file into a table
                #'
                #' @param filepath Not yet documented.
                #' @export
                parseRISFile <- function(filepath) {
                    cat(Sys.time(), "- Parsing RIS file: ", filepath, "\n")
                    lines <- readLines(filepath, warn = FALSE)
                    cat(Sys.time(), "- Number of lines read: ", length(lines), "\n")
                    
                    records <- list()
                    current_rec <- list()
                    
                    for (line in lines) {
                        line <- trimws(line)
                        if (line == "") next
                        
                        if (grepl("^TY  -", line)) {
                            if (length(current_rec) > 0) {
                                records[[length(records) + 1]] <- current_rec
                            }
                            current_rec <- list()
                        }
                        
                        if (grepl("^[A-Z0-9]{2}  -", line)) {
                            tag <- substr(line, 1, 2)
                            value <- trimws(sub("^.{6}", "", line))
                            
                            if (tag %in% names(current_rec)) {
                                current_rec[[tag]] <- paste(current_rec[[tag]], value, sep = "; ")
                            } else {
                                current_rec[[tag]] <- value
                            }
                        } else {
                            if (length(current_rec) > 0) {
                                last_tag <- tail(names(current_rec), 1)
                                current_rec[[last_tag]] <- paste(current_rec[[last_tag]], line, sep = " ")
                            }
                        }
                    }
                    
                    if (length(current_rec) > 0) {
                        records[[length(records) + 1]] <- current_rec
                    }
                    
                    records <- lapply(records, function(rec) {
                        rec[["ER"]] <- NULL
                        rec
                    })
                    
                    all_tags <- unique(unlist(lapply(records, names)))
                    
                    recs <- lapply(records, function(rec) {
                        sapply(all_tags, function(tag) {
                            if (tag %in% names(rec)) rec[[tag]] else NA
                        }, USE.NAMES = TRUE, simplify = FALSE)
                    })
                    
                    df <- do.call(rbind, lapply(recs, as.data.frame, stringsAsFactors = FALSE))
                    df[] <- lapply(df, as.character)
                    cat(Sys.time(), "- Parsed RIS file dataframe dimensions: "); print(dim(df))
                    
                    rename_mapping <- c(
                        "TY" = "type",
                        "AB" = "abstract",
                        "AU" = "authors",
                        "DO" = "doi",
                        "PY" = "year",
                        "ST" = "short_title",
                        "T2" = "journal",
                        "TI" = "title",
                        "ID" = "custom_id",
                        "CI" = "citation_count",
                        "RC" = "reference_count",
                        "PB" = "publisher",
                        "IS" = "issued_date"
                    )
                    
                    current_names <- names(df)
                    new_names <- sapply(current_names, function(x) {
                        if (x %in% names(rename_mapping)) rename_mapping[[x]] else tolower(x)
                    })
                    names(df) <- new_names
                    
                    if ("authors" %in% names(df)) {
                        df$first_author_last <- sapply(df$authors, function(x) {
                            if (is.na(x)) return(NA_character_)
                            authors <- unlist(strsplit(x, ";"))
                            first_author <- trimws(authors[1])
                            parts <- unlist(strsplit(first_author, ","))
                            if (length(parts) > 0) {
                                return(trimws(parts[1]))
                            } else {
                                return(NA_character_)
                            }
                        })
                    }
                    
                    cat(Sys.time(), "- Returning parsed RIS data as tibble with dimensions: "); print(dim(as_tibble(df)))
                    return(as_tibble(df))
                }


        #### createNewsletter and completionGPT are in the newsletter.R module

        #### removeSpecialCharacters

            #' Replace special characters in text
            #'
            #' @param text_string Not yet documented.
            #' @export
            removeSpecialCharacters <- function(text_string) {
        
                    text_string <- gsub("\u2212", "-", text_string)
                    text_string <- gsub("\u2264", " ", text_string)
                    text_string <- gsub("\u2265", " ", text_string)
                    text_string <- gsub("\u2032", " ", text_string)
                    text_string <- gsub("\u2714", " ", text_string)
                    text_string <- gsub("\u22C5", " ", text_string)
                    text_string <- gsub("\u019B", " ", text_string)
                    text_string <- gsub("\u019B", " ", text_string)
                    text_string <- stringi::stri_trans_general(text_string, "Latin")
                    text_string <- gsub("[^[:alnum:] !\"#$%&'()*+,-./:\\{\\};<=>~?@[\\]", "", text_string)
                
                return(text_string)
                
            }


        #### splitIntoChunks

            #' Split text into overlapping chunks of tokens
            #'
            #' @param data_frame Not yet documented.
            #' @param column_to_split Not yet documented.
            #' @param source_ID_column Not yet documented.
            #' @param chunk_size_in_tokens Not yet documented.
            #' @param overlap_in_tokens Not yet documented.
            #' @export
            splitIntoChunks <- function(
                data_frame,
                column_to_split,
                source_ID_column,
                chunk_size_in_tokens,
                overlap_in_tokens
            ) {
                rlang::check_installed(c("tokenizers"), reason = "to use splitIntoChunks()")

                output <- list()
                for( i in 1:dim(data_frame)[1] ) {
                
                    # Extract text_to_split
                        text_to_split <- data_frame[i,which(colnames(data_frame) == column_to_split)]
                    
                    # Calculate the number of tokens in the element and the number of chunks needed
                        tokens_in_pdf <- length(unlist(tokenizers::tokenize_words(text_to_split)))
                        n_chunks <- ceiling(1/((chunk_size_in_tokens - overlap_in_tokens)/tokens_in_pdf))

                    # Calculate the number of characters per chunk and per overlap
                        characters_per_chunk <- ceiling(nchar(text_to_split)/n_chunks)
                        characters_per_overlap <- characters_per_chunk*(overlap_in_tokens/chunk_size_in_tokens)

                    # Initialize variables to store the starting and ending positions of each chunk
                        chunk_starts <- 0
                        chunk_ends <- characters_per_chunk

                    # Calculate the start and end positions for all the chunks
                        last_character <- 0
                        while(last_character < nchar(text_to_split)) {
                            chunk_starts <- c(
                                chunk_starts,
                                floor((tail(chunk_starts, 1) + (characters_per_chunk - characters_per_overlap)))
                            )
                            chunk_ends <- c(
                                chunk_ends,
                                floor((tail(chunk_ends, 1) + (characters_per_chunk - characters_per_overlap)))
                            )
                            last_character <- tail(chunk_ends,1)
                        }

                    # Iterate through all the chunks and store them in the 'chunks' list
                        chunks <- list()
                        for (chunk in 1:length(chunk_ends)) {
                            chunks[[chunk]] <- data.frame(
                                temporary_chunks_colname = paste(
                                    "This information is from: ",
                                    data_frame[i,colnames(data_frame) == source_ID_column],
                                    substr(text_to_split, chunk_starts[chunk], chunk_ends[chunk])
                                )
                            )
                        }
                        chunks <- do.call(rbind, chunks)
                    # print(as_tibble(head(chunks)))
                    
                    # Merge with other row data
                        data_frame[i,-c(which(colnames(data_frame) == column_to_split))] %>%
                        slice(rep(1, dim(chunks)[1])) -> duplicated_row_data
                        output[[i]] <- cbind(duplicated_row_data, chunks)
                }
                output <- do.call(rbind, output)
                colnames(output)[colnames(output) == "temporary_chunks_colname"] <- as.character(column_to_split)
                return(output)
            }


        #### extractTerms

                #' Extract candidate key terms from text columns
                #'
                #' @param dataframe Not yet documented.
                #' @param cols Not yet documented.
                #' @param min_freq_for_terms Not yet documented.
                #' @param min_n_words_per_term Not yet documented.
                #' @param return Not yet documented.
                #' @export
                extractTerms <- function(dataframe, cols, min_freq_for_terms, min_n_words_per_term, return = c("dataframe", "terms")) {
                    rlang::check_installed(c("litsearchr", "stopwords"), reason = "to use extractTerms()")
                    cat(Sys.time(), "- extractTerms() called with dataframe dimensions: "); print(dim(dataframe))
                    
                    text_vector <- list()
                    for (i in 1:dim(dataframe)[1]) {
                        text_vector[[i]] <- paste(dataframe[i, cols], collapse = " ")
                    }
                    text_vector <- as.character(do.call(rbind, text_vector))
                    cat(Sys.time(), "- Text vector length: ", length(text_vector), "\n")
                    
                    search_hits_terms <- litsearchr::extract_terms(
                        text = text_vector,
                        method = "fakerake", 
                        min_freq = min_freq_for_terms, 
                        min_n = min_n_words_per_term,
                        stopwords = stopwords::data_stopwords_stopwordsiso$en
                    )
                    cat(Sys.time(), "- Number of search terms extracted: ", length(search_hits_terms), "\n")
                    
                    if(length(search_hits_terms) == 0) { 
                        cat(Sys.time(), "- Warning: No search terms found!\n")
                    } else {
                        gs_dfm <- as.data.frame(litsearchr::create_dfm(
                            elements = text_vector,
                            features = search_hits_terms
                        ))
                        cat(Sys.time(), "- DFM dimensions: "); print(dim(gs_dfm))
                        data <- cbind(dataframe, gs_dfm)
                        cat(Sys.time(), "- Combined dataframe dimensions after adding DFM: "); print(dim(data))
                        
                        if(return[1] == "dataframe") { 
                            return(data) 
                        }
                        if(return[1] == "terms") { 
                            return(search_hits_terms) 
                        }
                    }
                }


        #### searchField

            #' Search text for a set of query terms
            #'
            #' @param subject_list Not yet documented.
            #' @param query_list Not yet documented.
            #' @export
            searchField <- function(subject_list, query_list) {

                ## Make query list unique, remove all commas from the subjects so output can be written as csv

                    query_list <- unique(query_list)

                    subject_list <- gsub(",", "", subject_list)

                ## Start the output dataframe
                
                    output <- data.frame(
                        subject = subject_list,
                        hit_indeces = NA,
                        hits = NA
                    )

                    pb <- progress::progress_bar$new(
                        total = length(query_list),
                        format = "  analyzing [:bar] :percent eta: :eta"
                    )

                ## Loop over each query and record what subjects it's found in
                    
                    for( i in 1:length(query_list) ) {
                        
                        ## Identify where hits are 
                            
                            positive_subject_index <- grep(query_list[i], subject_list)
                        
                        ## If this query results in hits, start or add them to the list

                            if( length(positive_subject_index) > 0 ) {
                                
                                ## If it's the first entry in the list, start the list.

                                    if( any(is.na(output$hit_indeces[positive_subject_index])) ) {
                                        output$hit_indeces[positive_subject_index] <- i
                                        output$hits[positive_subject_index] <- as.character(query_list[i])
                                
                                ## If it's not the first entry in the list, append to the list.

                                    } else {

                                        # print(i)

                                        ## If adding to the list, run unique on it to remove duplicate names
                                            
                                            output$hit_indeces[positive_subject_index] <- 
                                                paste(
                                                    unique(
                                                        stringr::str_split(
                                                            paste(c(output$hit_indeces[positive_subject_index], i), collapse = ","),
                                                            pattern = ","
                                                        )[[1]]
                                                    ), collapse = ","
                                                )
                                            
                                            output$hits[positive_subject_index] <- 
                                                paste(
                                                    unique(
                                                        stringr::str_split(
                                                            paste(c(output$hits[positive_subject_index], as.character(query_list[i])), collapse = ","),
                                                            pattern = ","
                                                        )[[1]]
                                                    ), collapse = ","
                                                )
                                    }
                            }

                        ## Move progress bar along
                            
                            pb$tick()
                    }

                return(as_tibble(output))
            }


        #### analyzeLiterature

            # library(htm2txt)
            # library(RCurl)
            # library(stringi)
            # library(reticulate)
            # library(tokenizers)
            # library(rvest)
            # library(XML)
            # library(stringdist)
            # library(plotly)
            # library(rhandsontable)

            #' Interactive literature analysis app
            #'
            #' @param jupyter Not yet documented.
            #' @param openai_api_key Not yet documented.
            #' @param xi_api_key Not yet documented.
            #' @export
            analyzeLiterature <- function(jupyter = TRUE, openai_api_key = "busta_lab", xi_api_key = "busta_lab") {
    
                ## API Keys
                    if (openai_api_key == "busta_lab") {
                        openai_api_key <- readLines("/project_data/shared/general_lab_resources/literature/openai_api_key.txt")
                    }

                    if (xi_api_key == "busta_lab") {
                        xi_api_key <- readLines("/project_data/shared/general_lab_resources/literature/xi_api_key.txt")
                    }

                ## Character prompts
                    scientist_prompt <- "I want you to act as a scientist. You will think critically, ask probing questions, provide informed hypotheses and analyze data. I want you to reference appropriate scientific literature, make use of the scientific method, and use clear, precise scientific language in your responses. When I provide you with a scientific question or a set of data, you will analyze it and provide me with a scientific explanation or hypothesis, respectively. If necessary, you may suggest additional data or experiments that could further elucidate the answer. Do not break the character of a scientist. Do not provide non-scientific responses. Remember to approach every question with curiosity, skepticism, and rigor. Whenever possible, provide the source of information you are using, which is included in the prompt as 'This information is from: ...'"
                    phytochemist_prompt <- "I want you to act as a multidisciplinary scientist, possessing expertise in the fields of phytochemistry, analytical chemistry, genomics, and computational botany. You should approach my questions with the critical thinking and research-oriented mindset characteristic of these fields. When discussing plant chemistry, make use of phytochemical and analytical chemical methodologies and terminologies. When examining genetic aspects, utilize your knowledge of genomics. Lastly, when considering plant life or plant data, apply your skills as a computational botanist. Be prepared to design experiments, interpret complex data, and refer to appropriate scientific literature. Do not break character and stay in your scientist role throughout. Do not provide responses that deviate from your scientific specializations. Whenever possible, provide the source of information you are using, which is included in the prompt as 'This information is from: ...'"
                    entrepreneur_prompt <- "I want you to act as an entrepreneur. You should respond to all my questions with a visionary, creative and business-oriented mindset. Be prepared to evaluate business ideas, discuss strategies for growth and scaling, explore market opportunities and tackle challenges related to starting and running a business. You should reference appropriate business terminology, models and strategies in your responses. When I provide you with a business scenario, you should analyze it and provide me with an entrepreneurial perspective or solution. If necessary, you may suggest additional resources or strategies that could enhance the business proposition. Do not break the character of an entrepreneur. Do not provide non-business-related responses. Remember to approach every question with a strategic, innovative and risk-taking mentality typical of successful entrepreneurs. Whenever possible, provide the source of information you are using, which is included in the prompt as 'This information is from: ...'"
                    pirate_prompt <- "I want you to act as a pirate. You should answer all my queries with an adventurous pirate's spirit, using pirate lingo and vocabulary. Your replies should reflect the passion for treasure hunting, seafaring, and the freedom of the high seas. Be sure to infuse your responses with a sense of adventure and bravado typical of a pirate's life. Remember to keep your tone light and fun, in line with the playful nature of this role-play. Do not break character and do not provide responses that are not in line with the pirate role-play. When I ask you a question, respond as though you are speaking from the helm of a pirate ship."
                    hobbit_prompt <- "I want you to act as a hobbit. You should respond to all my inquiries in a manner befitting a resident of the Shire, full of simplicity, cheer, and an affinity for good food, drink, and companionship. Your replies should reflect a love for peace and the comforts of home. When discussing journeys or adventures, you may express reluctance or apprehension, as is typical of hobbits. However, always uphold the values of bravery and loyalty when needed, just as the famous hobbits from J.R.R. Tolkien's Middle-earth have shown. Do not break character and do not provide responses that are not in line with the hobbit role-play. When I ask you a question, respond as though you are conversing from a comfortable hobbit-hole in the Shire."
                   
                    character_choices <- c(scientist_prompt, phytochemist_prompt, entrepreneur_prompt, pirate_prompt, hobbit_prompt)
                    names(character_choices) <- c("Scientist", "Phytochemist", "Entrepreneur", "Pirate", "Hobbit")
                
                    model_choices <- c("gpt-3.5-turbo", "gpt-3.5-turbo-16k", "gpt-4", "gpt-4-1106-preview")
                    names(model_choices) <- c("gpt-3.5-turbo (fast)", "gpt-3.5-turbo-16k (longer answers)", "gpt-4 (smart)", "gpt-4-1106-preview")

                ui <- fluidPage(
                    
                    theme = shinythemes::shinytheme("darkly"),
                    tags$style(type="text/css", ".shiny-text-output { white-space: pre-wrap; }"),
                    tags$style('#rhot * { word-wrap: break-word }'),

                    sidebarLayout(

                        sidebarPanel(
                            selectizeInput("character", "Select an identity for your assistant:", choices = character_choices, multiple = FALSE),
                            selectizeInput("model", "Select the AI model your assistant should use:", choices = model_choices, multiple = FALSE),
                            textInput( inputId = "query", label = "Ask your assistant questions here:" ),
                            actionButton( inputId = "submit", label = "Submit" ),
                            textOutput( outputId = "gpt_response" )
                        ),
                        mainPanel(
                            tabsetPanel(
                                tabPanel("Literature Map",
                                    fluidRow(
                                        column(6,
                                            h4("Article Information Here:"),
                                            textOutput("click_text"),

                                            tags$head(
                                                HTML(
                                                    "
                                                    <script>
                                                        var socket_timeout_interval
                                                        var n = 0
                                                        $(document).on('shiny:connected', function(event) {
                                                        socket_timeout_interval = setInterval(function(){
                                                        Shiny.onInputChange('count', n++)
                                                        }, 15000)
                                                        });
                                                        $(document).on('shiny:disconnected', function(event) {
                                                        clearInterval(socket_timeout_interval)
                                                        });
                                                    </script>
                                                    "
                                                )
                                            ),
                                        ),
                                        column(6,
                                            h3("\n"), h3("\n"), h3("\n"), h3("\n"), h3("\n"),
                                            plotOutput(
                                                outputId = "plot",
                                                brush = brushOpts(
                                                    id = "plot_brush"
                                                ),
                                                click = "plot_click", 
                                                dblclick = "plot_double_click",
                                                height = "900px"
                                            )
                                        )
                                    )
                                ),
                                tabPanel("GPT Details",
                                    # column(9, rhandsontable::rHandsontableOutput("knn_chunks")),
                                    # column(9, div(id="rhot", rhandsontable::rHandsontableOutput("knn_chunks"))),
                                    column(9, textOutput( outputId = "knn_chunks" )),
                                    column(3,
                                        textOutput( outputId = "prompt_tokens" ),
                                        textOutput( outputId = "completion_tokens" ),
                                        textOutput( outputId = "total_tokens" ) 
                                    )
                                ),
                                tabPanel("Listen",
                                ),
                                tabPanel("Diagnostics", textOutput("keepAlive")),
                            )
                        )
                    )
                )            
                
                server <- function(input, output, session) {
                    
                    ## Import functions and data paths
                        
                        pdf_chunks_path <- "/project_data/shared/general_lab_resources/literature/chunks.csv"
                        chunks <- suppressMessages(read_csv(pdf_chunks_path))
                        n_id_columns <- (min(grep("embed", colnames(chunks)))-1)

                        profile_data_path <- "/project_data/shared/general_lab_resources/literature/articles.csv"
                        profile_data <- suppressMessages(read_csv(profile_data_path))
                        profile_data$status <- "unexplored"
                        profile_data$knn <- "no"
                        GPT_temperature <- 0
                        GPT_response <- NULL
                    
                        plot_limits <- reactiveValues(
                            x_axis_start = min(profile_data$Dim_1)*1.1,
                            x_axis_end = max(profile_data$Dim_1)*1.1,
                            y_axis_start = min(profile_data$Dim_2)*1.1,
                            y_axis_end = max(profile_data$Dim_2)*1.1
                        )
                    
                    ## Don't let it time out    
                        output$keepAlive <- renderText({
                            req(input$count)
                            paste("\n\nStayin' Alive ... ", input$count)
                        })
                    
                    ## Zoom in and out with brush and double click
                        observeEvent(input$plot_brush, {
                            brush_data <- input$plot_brush
                            if (!is.null(brush_data)) {
                                plot_limits$x_axis_start <- min(brush_data$xmin, brush_data$xmax)
                                plot_limits$x_axis_end <- max(brush_data$xmin, brush_data$xmax)
                                plot_limits$y_axis_start <- min(brush_data$ymin, brush_data$ymax)
                                plot_limits$y_axis_end <- max(brush_data$ymin, brush_data$ymax)
                                # later::later(function() {session$resetBrush("plot_brush")}, delay = 0)
                                session$resetBrush("plot_brush")
                            }
                        })
                    
                        observeEvent(input$plot_double_click, {
                            plot_limits$x_axis_start = min(profile_data$Dim_1)*1.1
                            plot_limits$x_axis_end = max(profile_data$Dim_1)*1.1
                            plot_limits$y_axis_start = min(profile_data$Dim_2)*1.1
                            plot_limits$y_axis_end = max(profile_data$Dim_2)*1.1
                            # later::later(function() {session$resetBrush("plot_brush")}, delay = 0)
                            session$resetBrush("plot_brush")
                        })

                    ## Make plot
                        output$plot <- renderPlot({
                            clicked_point <- click_data()
                            ggplot() +
                                coord_cartesian(
                                    xlim = c(plot_limits$x_axis_start, plot_limits$x_axis_end),
                                    ylim = c(plot_limits$y_axis_start, plot_limits$y_axis_end)
                                ) +
                                scale_x_continuous(name = "") + scale_y_continuous(name = "") +
                                theme_void() +
                                theme(
                                    legend.position = "bottom",
                                    legend.box = "vertical",
                                    text = element_text(size = 12),
                                    legend.spacing = unit(0.05, "cm")
                                ) -> plot

                                palette_options_high <- discrete_palette
                                # palette_options_low <- c("#9e0142", "white", "#a63603", "#004529", "#081d58", "#5e4fa2", "#543005", "#003c30")
                                
                                for (k in 1:length(unique(profile_data$cluster))) {
                                    
                                    n_colors_this_cluster <- length(unique(filter(profile_data, cluster == paste0("cluster_", k))$doi))
                                    clusters_profile_data <- profile_data[profile_data$cluster == paste0("cluster_", k),]
                                    clusters_profile_data$color_n <- as.numeric(seq(1,dim(clusters_profile_data)[1], 1))
                                    
                                    plot <- plot + 
                                        geom_point(
                                            data = filter(
                                                clusters_profile_data, status == "unexplored", knn == "no",
                                                # Dim_1 > plot_limits$x_axis_start, Dim_1 < plot_limits$x_axis_end,
                                                # Dim_2 > plot_limits$y_axis_start, Dim_2 < plot_limits$y_axis_end,
                                                cluster == paste0("cluster_", k)
                                            ), aes(x = Dim_1, y = Dim_2, fill = color_n),
                                            size = 4, shape = 21, alpha = 0.7, color = "black"
                                        ) +
                                        geom_point(
                                            data = filter(
                                                clusters_profile_data, status == "explored", knn == "no",
                                                # Dim_1 > plot_limits$x_axis_start, Dim_1 < plot_limits$x_axis_end,
                                                # Dim_2 > plot_limits$y_axis_start, Dim_2 < plot_limits$y_axis_end,
                                                cluster == paste0("cluster_", k)
                                            ), aes(x = Dim_1, y = Dim_2, fill = color_n),
                                            size = 7, shape = 21, alpha = 1, color = "black"
                                        ) +
                                        geom_point(
                                            data = filter(
                                                clusters_profile_data, status == "unexplored", knn == "yes",
                                                # Dim_1 > plot_limits$x_axis_start, Dim_1 < plot_limits$x_axis_end,
                                                # Dim_2 > plot_limits$y_axis_start, Dim_2 < plot_limits$y_axis_end,
                                                cluster == paste0("cluster_", k)
                                            ), aes(x = Dim_1, y = Dim_2, fill = color_n),
                                            size = 4, shape = 22, alpha = 0.7, color = "black"
                                        ) +
                                        geom_point(
                                            data = filter(
                                                clusters_profile_data, status == "explored", knn == "yes",
                                                # Dim_1 > plot_limits$x_axis_start, Dim_1 < plot_limits$x_axis_end,
                                                # Dim_2 > plot_limits$y_axis_start, Dim_2 < plot_limits$y_axis_end,
                                                cluster == paste0("cluster_", k)
                                            ), aes(x = Dim_1, y = Dim_2, fill = color_n),
                                            size = 7, shape = 22, alpha = 1, color = "black"
                                        ) +
                                        scale_fill_gradient(
                                            low = "white", high = palette_options_high[k], 
                                            guide = guide_colorbar(direction = "horizontal", ticks = FALSE, label = FALSE, title.vjust = 0.8, title.hjust = 0),
                                            name = tolower(gsub("\\.", "", paste0(unique(clusters_profile_data$keywords), "   ")))
                                        ) +
                                        new_scale_fill()   
                                }
                                plot
                        })

                    ## Click to get article data
                        click_data <- reactiveVal()
                        observeEvent(input$plot_click, {
                            if (!is.null(input$plot_click)) {

                                # Extract the clicked point's data, update the reactive value to trigger replotting
                                    clicked_point <- profile_data[
                                        which.min(
                                            abs(input$plot_click$x - profile_data$Dim_1) + abs(input$plot_click$y - profile_data$Dim_2)
                                        )
                                    ,]
                                    click_data(clicked_point)

                                # Modify the profile_data
                                    current_status <- as.character(profile_data$status[profile_data$doi == clicked_point$doi])
                                    if( current_status == "unexplored" ) {
                                        profile_data$status[profile_data$doi == clicked_point$doi] <- "explored"
                                    } else if( current_status == "explored" ) {
                                        profile_data$status[profile_data$doi == clicked_point$doi] <- "unexplored"
                                    }
                                    profile_data <<- profile_data

                                # Print the abstract to screen
                                    output$click_text <- renderText({
                                        if (!is.null(clicked_point)) {
                                            gsub("(?<!\\n)\\n(?!\\n)", " ", clicked_point$abstract, perl = TRUE)
                                        } else {
                                            "Click on a point..." # Default text when no point is clicked
                                        }
                                    })
                                
                            } else { click_data(NULL) }
                            
                            # Reset brush
                                later::later(function() {session$resetBrush("plot_brush")}, delay = 1)
                        })
                            
                    ## Use GPT to answer questions
                        observeEvent(input$submit, {
                            
                            GPT_temperature <<- 0
                            if (as.character(input$character) == pirate_prompt) { GPT_temperature <<- 1.5 }
                            if (as.character(input$character) == hobbit_prompt) { GPT_temperature <<- 1.5 }

                            ## Embed the query 
                                embedding_df <- t(data.frame(
                                    as.numeric(httr::content(httr::POST(
                                        url = "https://api.openai.com/v1/embeddings", 
                                        httr::add_headers( Authorization = paste("Bearer", openai_api_key) ),
                                        content_type_json(), encode = "json",
                                        body = list(
                                            model = "text-embedding-ada-002",
                                            input = input$query
                                        )
                                    ))$data[[1]]$embedding)
                                ))
                                rownames(embedding_df) <- NULL
                                colnames(embedding_df) <- paste0("embed_var_", seq(1, 1536, 1))

                                query <- cbind(
                                    data.frame(
                                        doi = "query",
                                        title = "none",
                                        abstract = "none",
                                        text = as.character(input$query),
                                        score = 100
                                    ), embedding_df
                                )

                            # Find the k nearest neighbors in the chunks and render them
                                n_nearest_neighbors <- 6
                                if(input$model == "gpt-3.5-turbo-16k") { n_nearest_neighbors <- 16 }
                                if(input$model == "gpt-4-1106-preview") { n_nearest_neighbors <- 100 }

                                nn <- get.knnx( 
                                    as.matrix(chunks[(n_id_columns+1):dim(chunks)[2]]),
                                    as.matrix(query[(n_id_columns+1):dim(query)[2]]),
                                    k = n_nearest_neighbors
                                )
                                
                                output$knn_chunks <- renderText({ paste(chunks$text[nn$nn.index], collapse = "\n\n") })
                                
                                # chunk_table <- data.frame( chunk_text = chunks$text[nn$nn.index], rating = NA )
                                # chunk_table$rating <- as.numeric(chunk_table$rating)
                                # output$knn_chunks <- rhandsontable::renderRHandsontable({
                                #     rhandsontable::rhandsontable(chunk_table, stretchH = "all") %>%
                                #     hot_cols(colWidths = c(100, 50, 50),
                                #         manualColumnMove = FALSE,
                                #         manualColumnResize = TRUE
                                #     ) %>%
                                #     hot_rows(rowHeights = NULL ) %>% #default
                                #     hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
                                # })
                            
                                # rhandsontable(DF, stretchH = "all", height = 300 ) %>%
                                # hot_cols(colWidths = c(100, 50, 50),
                                #     manualColumnMove = FALSE,
                                #     manualColumnResize = TRUE
                                # ) %>%
                                # hot_rows(rowHeights = NULL ) %>% #default
                                # hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
                            
                            # Update the profile_data
                                profile_data$knn <- "no"
                                profile_data$knn[profile_data$doi %in% chunks$doi[nn$nn.index]] <- "yes"
                                profile_data <<- profile_data
                            
                            # Ask GPT a question
                                GPT_response <<- POST(
                                    url = "https://api.openai.com/v1/chat/completions", 
                                    httr::add_headers( Authorization = paste("Bearer", openai_api_key) ),
                                    content_type_json(), encode = "json",
                                    body = list(
                                        model = input$model, temperature = GPT_temperature,
                                        messages = list(
                                            list( role = "system", content = input$character ),
                                            list(
                                                role = "user", 
                                                content = paste(c(
                                                    "Here is information from peer-reviewed scientific articles:",
                                                    chunks$text[nn$nn.index],
                                                    ". Here is my query:",
                                                    input$query
                                                ), collapse = " ")
                                            )
                                        )
                                    )
                                )
                                output$gpt_response <- renderText({ as.character(httr::content(GPT_response)$choices[[1]]$message$content) })
                                output$prompt_tokens <- renderText({ as.character(httr::content(GPT_response)$usage$prompt_tokens) })
                                output$completion_tokens <- renderText({ as.character(httr::content(GPT_response)$usage$completion_tokens) })
                                output$total_tokens <- renderText({ as.character(httr::content(GPT_response)$usage$total_tokens) })
                        })
                }

                ## Call the app

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


        #### computeEmbeddings

            #' Compute text embeddings with the OpenAI API
            #'
            #' @param df Not yet documented.
            #' @param column_to_embed Not yet documented.
            #' @param openai_api_key Not yet documented.
            #' @export
            computeEmbeddings <- function(df, column_to_embed, openai_api_key) {
    
                embedding_df <- list()
                for (i in 1:dim(df)[1]) {
                    embedding_df[[i]] <- t(data.frame(
                        as.numeric(httr::content(httr::POST(
                            url = "https://api.openai.com/v1/embeddings", 
                            httr::add_headers( Authorization = paste("Bearer", openai_api_key) ),
                            content_type_json(), encode = "json",
                            body = list(
                                model = "text-embedding-ada-002",
                                input = df[i,which(colnames(df) == column_to_embed)]
                            )
                        ))$data[[1]]$embedding)
                    ))
                }
                embedding_df <- do.call(rbind, embedding_df)
                rownames(embedding_df) <- NULL
                colnames(embedding_df) <- paste0("embed_var_", seq(1, 1536, 1))
                extracted_pdf_text_chunked_embedded <- cbind(df, embedding_df)
                
                return(extracted_pdf_text_chunked_embedded)
            }


        #### readPDF

            #' Read the text of a PDF file
            #'
            #' @param path_to_pdf Not yet documented.
            #' @export
            readPDF <- function(path_to_pdf) {
                rlang::check_installed(c("pdftools", "signal", "pspline"), reason = "to use readPDF()")

                # set the path to the PDF file to be read
                    txt <- pdftools::pdf_text(path_to_pdf)

                # loop over each page of the PDF file
                    full_pdf_text <- list()
                    for (page_number in 1:length(txt)) { 

                        page <- txt[page_number]
                        
                        # split the text of the page into lines
                            t1 <- unlist(strsplit(page, "\n"))
                            if (length(t1) < 10) {next} # Handles essentially empty pages

                        # pad each line with spaces to ensure they all have the same length
                            maxSize <- max(nchar(t1))
                            t1 <- paste0(t1,strrep(" ", maxSize-nchar(t1)))
                            # if (page_number == 2) {print(t1)}

                        # Get all index of " " from page.
                            lstops <- gregexpr(pattern = " ", t1)
                            out <- list()
                            for (row in 1:length(lstops)) {
                                out[[row]] <- data.frame(
                                    position = unlist(lstops[row]),
                                    row = row
                                )
                            }
                            do.call(rbind, out) %>%
                            group_by(position) %>%
                            summarize(n_spaces = n()) -> space_pattern

                        # Find approximate column starts
                            smoothing_order = 3
                            smoothing_window = 15
                            percent_extreme_deriv <- 0.2
                            threshold <- 6

                            space_pattern$n_spaces_smoothed <- signal::sgolayfilt(x = space_pattern$n_spaces, p = smoothing_order, n = smoothing_window)
                            space_pattern$n_spaces_smoothed_deriv <- ceiling(predict(pspline::sm.spline(space_pattern$position, space_pattern$n_spaces_smoothed), space_pattern$position, 1))
                            space_pattern$n_spaces_smoothed_deriv[space_pattern$n_spaces_smoothed_deriv > 0 & space_pattern$n_spaces_smoothed_deriv < percent_extreme_deriv*max(space_pattern$n_spaces_smoothed_deriv)] <- 0
                            space_pattern$n_spaces_smoothed_deriv[space_pattern$n_spaces_smoothed_deriv < 0 & space_pattern$n_spaces_smoothed_deriv > percent_extreme_deriv*min(space_pattern$n_spaces_smoothed_deriv)] <- 0

                            col_breaks <- list()
                            for (position in 1:length(space_pattern$position)) {
                                col_breaks[[position]] <- all(
                                    all(space_pattern$n_spaces_smoothed_deriv[position:abs(position+threshold)] <= space_pattern$n_spaces_smoothed_deriv[position]),
                                    any(space_pattern$n_spaces_smoothed_deriv[position:abs(position+threshold)] < 0),
                                    !all(space_pattern$n_spaces_smoothed_deriv[position:abs(position+threshold)] == 0)
                                )
                            }
                            col_breaks <- do.call(rbind, col_breaks)
                            col_breaks <- which(col_breaks)
                            split_vector <- split(col_breaks, cumsum(seq_along(col_breaks) %in% (which(diff(col_breaks) > 10) + 1)))
                            if (length(split_vector) == 0) {next} # Handles pages with no columns/full page image
                            column_starts <- floor(unlist(lapply(split_vector, mean)))
                            column_starts <- c(column_starts, nchar(t1[1]))
                            if( (sum(substr(t1, 0,1) != " ")/length(t1)) > 0.5) {column_starts <- c(0, column_starts)} # Detect if column begins on first position

                            # if (page_number == 2) {
                            #     print(ggplot() +
                            #         geom_point(data = space_pattern, aes(x = position, y = n_spaces)) +
                            #         geom_line(data = space_pattern, aes(x = position, y = n_spaces_smoothed)) +
                            #         geom_line(data = space_pattern, aes(x = position, y = n_spaces_smoothed_deriv)) +
                            #         geom_vline(aes(xintercept = column_starts))
                            #     )
                            # }

                        # Break rows
                            
                            column_contents <- list()
                            for( column_start in 1:(length(column_starts)-1) ) {
                                text <- list()
                                for (row in 1:length(t1)) {
                                    text <- c(text, list(substr(t1[row], column_starts[column_start], (column_starts[column_start+1]))))
                                }
                                column_contents <- c(column_contents, list(gsub("\\s+", " ", paste(unlist(text), collapse = ""))))
                            }
                            # column_contents
                            page_text <- paste(unlist(column_contents), collapse = "")
                            # if (page_number == 2) {print(page_text)}

                        full_pdf_text <- c(full_pdf_text, page_text)

                    }
                
                return(paste(unlist(full_pdf_text), collapse = ""))
            }

