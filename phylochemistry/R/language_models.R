## language_models — searchNCBI, searchPubMed, embedText, generateText, embedAminoAcids moved verbatim
## from modules/language_model_analysis.R (the versions book ch13/14 run); the rest of
## literature tooling is in literature.R.

    #### searchNCBI

        #' Search NCBI for protein (or nucleotide) records and return their sequences
        #'
        #' @param search_term The term to search for in the NCBI database.
        #' @param retmax The maximum number of results to return. Default is 5.
        #' @param type `"protein"` (the default) or `"nucleotide"`.
        #' @return An `AAStringSet` (protein) or `DNAStringSet` (nucleotide), or `NULL` if nothing matched.
        #' @export
        searchNCBI <- function(search_term, retmax = 5, type = "protein") {
            .check_bioc(c("Biostrings"), "searchNCBI")
            # `type` comes from the phylochemistry.R copy of this function; the protein path is the
            # modules/language_model_analysis.R version (the one book ch14 ran), unchanged.
            if (!identical(type, "protein")) {
                db <- ifelse(type == "protein", "protein", "nucleotide")
                search_results <- rentrez::entrez_search(db = db, term = search_term, retmax = retmax)
                if (length(search_results$ids) == 0) {
                    message("No sequences found for the search term.")
                    return(NULL)
                }
                temp_fasta <- tempfile(fileext = ".fasta")
                write(rentrez::entrez_fetch(db = db, id = search_results$ids, rettype = "fasta"), file = temp_fasta)
                return(Biostrings::readDNAStringSet(temp_fasta))
            }

            # Search, message if no results
                search_results <- rentrez::entrez_search(db = "protein", term = search_term, retmax = retmax)
                if (length(search_results$ids) == 0) {
                    message("No proteins found for the search term.")
                    return(NULL)
                }

            # Write results to temp file then read in as stringset
                temp_fasta <- tempfile(fileext = ".fasta")
                write(
                    rentrez::entrez_fetch(db = "protein", id = search_results$ids, rettype = "fasta"),
                    file = temp_fasta
                )
                    return(Biostrings::readAAStringSet(temp_fasta))
        }

    #### searchPubMed

        #' Search PubMed and return article metadata
        #'
        #' @param search_terms Not yet documented.
        #' @param pubmed_api_key Not yet documented.
        #' @param sort Not yet documented.
        #' @param retmax_per_term Not yet documented.
        #' @export
        searchPubMed <- function(search_terms, pubmed_api_key, sort = c("date", "relevance"), retmax_per_term = 20) {

              pm_entries <- character()
              term_vector <- character()

              for (i in 1:length(search_terms)) { 
                  search_output <- rentrez::entrez_search(
                    db = "pubmed", term = as.character(search_terms[i]), 
                    retmax = retmax_per_term, use_history = TRUE, sort = sort[1]
                  )
                  
                  # Initialize variables for retry mechanism
                  success <- FALSE
                  attempts <- 0
                  max_attempts <- 3  # Maximum number of retry attempts
                  
                  while (!success && attempts < max_attempts) {
                      attempts <- attempts + 1
                      
                      # Attempt to fetch data
                      query_output <- try(rentrez::entrez_fetch(
                        db = "pubmed", web_history = search_output$web_history, 
                        rettype = "xml", retmax = retmax_per_term, 
                        api_key = pubmed_api_key, timeout = 60), silent = TRUE)
                      
                      # Check if the attempt was successful
                      if (inherits(query_output, "try-error")) {
                          message("Error encountered. Attempt ", attempts, " of ", max_attempts, ". Retrying in 5 seconds...")
                          Sys.sleep(5)  # Wait before retrying
                      } else {
                          # Parse and store the entries if successful
                          current_pm_entries <- XML::xmlToList(XML::xmlParse(query_output))
                          pm_entries <- c(pm_entries, current_pm_entries)
                          term_vector <- c(term_vector, rep(as.character(search_terms[i]), length(current_pm_entries)))
                          success <- TRUE
                      }
                  }
                  
                  if (!success) {
                      message("Failed to fetch data for term: ", as.character(search_terms[i]), " after ", max_attempts, " attempts.")
                  }
                  
                  Sys.sleep(4)  # Wait between different search terms to respect API rate limits
              }

              unique_indices <- !duplicated(pm_entries)
              pm_entries <- pm_entries[unique_indices]
              term_vector <- term_vector[unique_indices]

              pm_results <- list()
              for (i in 1:length(pm_entries)) { # i=1

                  if (length(pm_entries[[i]]) == 1) { next }
                  if (is.null(pm_entries[[i]]$MedlineCitation$Article$ELocationID$text)) { next }

                  options <- which(names(pm_entries[[i]]$MedlineCitation$Article) == "ELocationID")
                  for (option in options) { # option = 4
                      if (grepl("10\\.", pm_entries[[i]]$MedlineCitation$Article[[option]]$text)) {
                          doi <- pm_entries[[i]]$MedlineCitation$Article[[option]]$text
                          break
                      } else {next}
                  }

                  pm_results[[i]] <- data.frame(
                      entry_number = as.numeric(i),
                      term = term_vector[[i]],
                      date = lubridate::as_date(paste(
                          pm_entries[[i]]$MedlineCitation$DateRevised$Year,
                          pm_entries[[i]]$MedlineCitation$DateRevised$Month,
                          pm_entries[[i]]$MedlineCitation$DateRevised$Day,
                      sep = "-"
                      )),
                      journal = pm_entries[[i]]$MedlineCitation$Article$Journal$Title,
                      title = paste0(pm_entries[[i]]$MedlineCitation$Article$ArticleTitle, collapse = ""),
                      doi = doi,
                      abstract = paste0(pm_entries[[i]]$MedlineCitation$Article$Abstract$AbstractText, collapse = "")
                  )
              }
              return(as_tibble(do.call(rbind, pm_results)))
        }

    #### embedText

        #' Embed text with a language model
        #'
        #' @param df Not yet documented.
        #' @param column_name Not yet documented.
        #' @param hf_api_key Not yet documented.
        #' @param path_to_glove_file Not yet documented.
        #' @param local Not yet documented.
        #' @param model_id Not yet documented.
        #' @param batch_size Not yet documented.
        #' @param max_retries Not yet documented.
        #' @param timeout_sec Not yet documented.
        #' @param server_url Not yet documented.
        #' @param server_token Not yet documented.
        #' @export
        embedText <- function(
          df,
          column_name,
          hf_api_key = Sys.getenv("HF_TOKEN"),
          path_to_glove_file = NULL,
          local = FALSE,
          model_id = "BAAI/bge-small-en-v1.5",
          batch_size = 16,
          max_retries = 5,
          timeout_sec = 60,
          server_url = Sys.getenv("EMBED_SERVER_URL", ""),
          server_token = Sys.getenv("EMBED_SERVER_TOKEN", "")
        ) {
          # deps: imported through the package NAMESPACE (no longer attached to the session)
          verbose <- FALSE  # read by post_with_retries(); embedText() has no verbose argument

          # --- prep input ----
          df <- as_tibble(df)
          if (!column_name %in% colnames(df)) stop("`column_name` not found in df.")
          text_vector <- as.character(df[[column_name]])
          text_vector[is.na(text_vector)] <- ""  # avoid NAs to API

          if (!is.null(path_to_glove_file) && nzchar(path_to_glove_file)) {
            local <- TRUE
          }

          if (local) {
            # ---------- LOCAL (GloVe) ----------
            if (is.null(path_to_glove_file) || !nzchar(path_to_glove_file)) {
              stop("GloVe file path is required when using local embeddings.")
            }
            if (!file.exists(path_to_glove_file)) {
              stop("GloVe file not found at: ", path_to_glove_file)
            }
            # load GloVe
            glove <- fread(path_to_glove_file, data.table = FALSE, quote = "")
            rownames(glove) <- glove[, 1]
            glove <- glove[, -1, drop = FALSE]

            # simple tokenizer
            word_tokenizer <- function(x) {
              # split on non-letters/numbers, lowercase
              toks <- str_split(tolower(x), "[^a-z0-9_]+", simplify = FALSE)[[1]]
              toks[nzchar(toks)]
            }

            pb <- txtProgressBar(min = 0, max = length(text_vector), style = 3)
            embeddings_list <- vector("list", length(text_vector))
            for (i in seq_along(text_vector)) {
              tokens <- word_tokenizer(text_vector[i])
              if (length(tokens) == 0) {
                # empty -> zeros
                embeddings_list[[i]] <- rep(0, ncol(glove))
              } else {
                present <- tokens[tokens %in% rownames(glove)]
                if (length(present) == 0) {
                  embeddings_list[[i]] <- rep(0, ncol(glove))
                } else {
                  embeddings_list[[i]] <- colMeans(glove[present, , drop = FALSE], na.rm = TRUE)
                }
              }
              setTxtProgressBar(pb, i)
            }
            close(pb)

            embeddings_df <- as.data.frame(do.call(rbind, embeddings_list))
            colnames(embeddings_df) <- paste0("embedding_", seq_len(ncol(embeddings_df)))
            out <- bind_cols(df, embeddings_df)
            return(out)
          }

          # ---------- REMOTE (self-hosted embed proxy) ----------
          # When `server_url` is set (arg or EMBED_SERVER_URL), embed against a
          # self-hosted, CPU-only embed proxy instead of Hugging Face. The proxy
          # takes {"texts": [...]} and returns {"embeddings": [[...]], "dim": N}.
          # No metered API, no key spend; requires a shared class bearer token.
          if (!is.null(server_url) && nzchar(server_url)) {
            server_url <- str_trim(server_url[1])
            embed_endpoint <- paste0(sub("/+$", "", server_url), "/embed")
            # Explicit browser-like User-Agent: the proxy sits behind Cloudflare,
            # whose Bot Fight Mode 403s default library agents (e.g. Python-urllib).
            req_headers <- c(
              `Content-Type` = "application/json",
              `User-Agent` = "Mozilla/5.0 (phylochemistry embedText)"
            )
            if (nzchar(server_token)) {
              req_headers <- c(req_headers, Authorization = paste0("Bearer ", str_trim(server_token[1])))
            }

            n <- length(text_vector)
            idx <- split(seq_len(n), ceiling(seq_len(n) / batch_size))
            pb <- txtProgressBar(min = 0, max = length(idx), style = 3)
            all_rows <- vector("list", length(idx))
            out_dim <- NULL

            for (b in seq_along(idx)) {
              batch_text <- as.list(text_vector[idx[[b]]])
              payload <- jsonlite::toJSON(list(texts = batch_text), auto_unbox = TRUE)
              resp <- httr::POST(
                url = embed_endpoint,
                httr::add_headers(.headers = req_headers),
                body = payload, encode = "raw", timeout(timeout_sec)
              )
              code <- resp$status_code
              if (code < 200 || code >= 300) {
                msg <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
                stop(sprintf("Embed proxy request failed (HTTP %s): %s", code, msg))
              }
              parsed <- jsonlite::fromJSON(
                httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = FALSE
              )
              emb <- parsed$embeddings
              mat <- do.call(rbind, lapply(emb, function(v) as.numeric(unlist(v, recursive = TRUE))))
              if (is.null(out_dim)) out_dim <- ncol(mat)
              if (ncol(mat) != out_dim) stop("Inconsistent embedding dimensions returned by the proxy.")
              all_rows[[b]] <- mat
              setTxtProgressBar(pb, b)
            }
            close(pb)

            embeddings_mat <- do.call(rbind, all_rows)
            colnames(embeddings_mat) <- paste0("embedding_", seq_len(ncol(embeddings_mat)))
            embeddings_df <- as.data.frame(embeddings_mat, stringsAsFactors = FALSE)
            out <- bind_cols(df, embeddings_df)
            return(out)
          }

          # ---------- REMOTE (Hugging Face Router) ----------
          if (is.null(hf_api_key) || hf_api_key == "") {
            stop("Please provide an HF API key via `hf_api_key` or set HF_TOKEN env var.")
          }
          if (length(hf_api_key) > 1) hf_api_key <- hf_api_key[1]
          hf_api_key <- str_trim(hf_api_key)

          base_url <- sprintf(
            "https://router.huggingface.co/hf-inference/models/%s/pipeline/feature-extraction",
            URLencode(model_id, reserved = TRUE)
          )

          # helper: POST with retries/backoff
          post_with_retries <- function(payload_json) {
            delay <- 1
            for (attempt in seq_len(max_retries)) {
              resp <- httr::POST(
                url = base_url,
                httr::add_headers(
                  Authorization = paste0("Bearer ", hf_api_key),
                  `Content-Type` = "application/json"
                ),
                body = payload_json,
                encode = "raw",
                timeout(timeout_sec)
              )

              code <- resp$status_code
              if (code >= 200 && code < 300) return(resp)

              # handle 429 / 5xx with backoff
              if (code == 429 || (code >= 500 && code < 600)) {
                if (verbose) {
                  msg_retry <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
                  message(sprintf("Attempt %d retryable error (HTTP %s): %s", attempt, code, msg_retry))
                }
                retry_after_raw <- httr::headers(resp)[["retry-after"]]
                retry_after <- suppressWarnings(as.numeric(retry_after_raw))
                if (length(retry_after) == 1 && is.finite(retry_after)) {
                  Sys.sleep(retry_after)
                } else {
                  Sys.sleep(delay)
                  delay <- min(delay * 2, 30)
                }
                next
              }

              # Other errors: stop with message body
              msg <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
              stop(sprintf("HF request failed (HTTP %s): %s", code, msg))
            }
            stop("HF request failed after retries.")
          }

          # batching
          n <- length(text_vector)
          idx <- split(seq_len(n), ceiling(seq_len(n) / batch_size))

          pb <- txtProgressBar(min = 0, max = length(idx), style = 3)
          all_rows <- vector("list", length(idx))
          out_dim <- NULL

          for (b in seq_along(idx)) {
            batch_idx <- idx[[b]]
            batch_text <- as.list(text_vector[batch_idx])

            # Body: either single string or list of strings is accepted.
            # We always send a list to keep parsing simple.
            body <- list(inputs = batch_text)
            payload <- jsonlite::toJSON(body, auto_unbox = TRUE)

            resp <- post_with_retries(payload)
            txt <- httr::content(resp, as = "text", encoding = "UTF-8")
            parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)

            # Expected: list(list(numeric ...), list(numeric ...), ...)
            # Some models may return matrix-like nested lists; handle both.
            # Ensure we end up with a matrix rows = items, cols = dims.
            # If single item came back as a single vector, coerce to list-of-one.
            if (!is.list(parsed[[1]]) && is.numeric(unlist(parsed))) {
              # single vector -> wrap
              parsed <- list(parsed)
            }

            # Convert each element to numeric vector
            mat <- do.call(rbind, lapply(parsed, function(v) as.numeric(unlist(v, recursive = TRUE))))
            if (is.null(out_dim)) out_dim <- ncol(mat)
            if (ncol(mat) != out_dim) {
              stop("Inconsistent embedding dimensions returned by the model.")
            }

            all_rows[[b]] <- mat
            setTxtProgressBar(pb, b)
          }
          close(pb)

          embeddings_mat <- do.call(rbind, all_rows)
          colnames(embeddings_mat) <- paste0("embedding_", seq_len(ncol(embeddings_mat)))
          embeddings_df <- as.data.frame(embeddings_mat, stringsAsFactors = FALSE)

          # preserve original row order
          out <- bind_cols(df, embeddings_df)
          return(out)
        }

    #### generateText

        #' Generate text with a language model
        #'
        #' @param df Not yet documented.
        #' @param prompt_column Not yet documented.
        #' @param system_column Not yet documented.
        #' @param hf_api_key Not yet documented.
        #' @param model_id Not yet documented.
        #' @param provider Not yet documented.
        #' @param max_retries Not yet documented.
        #' @param timeout_sec Not yet documented.
        #' @param temperature Not yet documented.
        #' @param verbose Not yet documented.
        #' @param max_new_tokens Not yet documented.
        #' @param server_url Not yet documented.
        #' @param server_token Not yet documented.
        #' @export
        generateText <- function(
          df,
          prompt_column,
          system_column = NULL,
          hf_api_key = Sys.getenv("HF_TOKEN"),
          model_id = "Qwen/Qwen2.5-14B-Instruct",
          provider = "featherless-ai",
          max_retries = 5,
          timeout_sec = 120,
          temperature = NULL,
          verbose = FALSE,
          max_new_tokens = 512,
          server_url = Sys.getenv("GENERATE_SERVER_URL", ""),
          server_token = Sys.getenv("GENERATE_SERVER_TOKEN", "")
        ) {
          # deps: imported through the package NAMESPACE (no longer attached to the session)

          df <- as_tibble(df)
          if (!prompt_column %in% colnames(df)) stop("`prompt_column` not found in df.")
          if (!is.null(system_column) && !system_column %in% colnames(df)) {
            stop("`system_column` not found in df.")
          }

          # ---- Class generation proxy branch -------------------------------------------
          # When `server_url` is set (arg or GENERATE_SERVER_URL), generate against the lab's
          # course proxy (https://generate.lbusta.org) instead of HuggingFace. Same interface,
          # same return shape — the df plus a `generation` column — so nothing calling this
          # needs to change. Mirrors embedText's server branch exactly.
          #
          # This is the supported path for CHEM class work: HuggingFace generation is dead
          # (which is why the chapter 12 generateText chunk was disabled), and the proxy is
          # token-authenticated, output-capped and budget-limited on the lab's side.
          if (!is.null(server_url) && nzchar(server_url)) {
            # Extract prompts here: the shared extraction below happens AFTER the HF-key
            # check, which this branch deliberately skips (the proxy needs no HF key).
            prompts <- as.character(df[[prompt_column]]); prompts[is.na(prompts)] <- ""
            systems <- NULL
            if (!is.null(system_column)) {
              systems <- as.character(df[[system_column]]); systems[is.na(systems)] <- ""
            }
            server_url <- str_trim(server_url[1])
            gen_endpoint <- paste0(sub("/+$", "", server_url), "/generate")
            req_headers <- c(
              `Content-Type` = "application/json",
              # Cloudflare's Bot Fight Mode 403s default library user-agents; curl-alikes with
              # a browser UA pass. A 403 here means the edge, not your token.
              `User-Agent` = "Mozilla/5.0 (phylochemistry generateText)"
            )
            if (nzchar(server_token)) {
              req_headers <- c(req_headers,
                               Authorization = paste0("Bearer ", str_trim(server_token[1])))
            }
            generations <- character(length(prompts))
            for (i in seq_along(prompts)) {
              payload <- list(prompt = prompts[i], max_tokens = max_new_tokens)
              if (!is.null(systems) && nzchar(systems[i])) payload$system <- systems[i]
              resp <- httr::POST(
                url = gen_endpoint,
                httr::add_headers(.headers = req_headers),
                body = jsonlite::toJSON(payload, auto_unbox = TRUE),
                encode = "raw", httr::timeout(timeout_sec)
              )
              code <- httr::status_code(resp)
              if (code == 401) {
                stop("Generation server rejected the token. Set GENERATE_SERVER_TOKEN to the ",
                     "class token (ask Lucas).")
              }
              if (code == 429) {
                stop("Generation server is rate-limited or the class daily budget is used up. ",
                     "Wait a little and try again.")
              }
              if (code >= 300) {
                stop("Generation server said ", code, ": ",
                     substr(httr::content(resp, "text", encoding = "UTF-8"), 1, 200))
              }
              parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
              generations[i] <- if (is.null(parsed$completion)) NA_character_ else parsed$completion
              if (verbose) message("  generated ", i, "/", length(prompts))
            }
            return(bind_cols(df, tibble(generation = generations)))
          }
          # ---- end class proxy branch ---------------------------------------------------

          if (is.null(hf_api_key) || hf_api_key == "") {
            stop("Please provide an HF API key via `hf_api_key` or set HF_TOKEN env var.")
          }
          if (length(hf_api_key) > 1) hf_api_key <- hf_api_key[1]
          hf_api_key <- str_trim(hf_api_key)
          provider <- str_trim(provider)

          prompts <- as.character(df[[prompt_column]])
          prompts[is.na(prompts)] <- ""
          systems <- NULL
          if (!is.null(system_column)) {
            systems <- as.character(df[[system_column]])
            systems[is.na(systems)] <- ""
          }

          if (nzchar(provider)) {
            base_url <- sprintf("https://router.huggingface.co/%s/v1/chat/completions", provider)
          } else {
            base_url <- "https://router.huggingface.co/v1/chat/completions"
          }

          post_with_retries <- function(payload_json) {
            delay <- 1
            for (attempt in seq_len(max_retries)) {
              resp <- httr::POST(
                url = base_url,
                httr::add_headers(
                  Authorization = paste0("Bearer ", hf_api_key),
                  `Content-Type` = "application/json"
                ),
                body = payload_json,
                encode = "raw",
                timeout(timeout_sec)
              )

              code <- resp$status_code
              if (code >= 200 && code < 300) return(resp)

              if (code == 429 || (code >= 500 && code < 600)) {
                if (verbose) {
                  msg_retry <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
                  message(sprintf("Attempt %d retryable error (HTTP %s): %s", attempt, code, msg_retry))
                }
                retry_after_raw <- httr::headers(resp)[["retry-after"]]
                retry_after <- suppressWarnings(as.numeric(retry_after_raw))
                if (length(retry_after) == 1 && is.finite(retry_after)) {
                  Sys.sleep(retry_after)
                } else {
                  Sys.sleep(delay)
                  delay <- min(delay * 2, 30)
                }
                next
              }

              msg <- tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) "")
              if (verbose) {
                message(sprintf("Attempt %d failed (HTTP %s): %s", attempt, code, msg))
              }
              stop(sprintf("HF request failed (HTTP %s): %s", code, msg))
            }
            stop("HF request failed after retries.")
          }

          generations <- character(length(prompts))
          pb <- txtProgressBar(min = 0, max = length(prompts), style = 3)

          for (i in seq_along(prompts)) {
            user_msg <- list(role = "user", content = prompts[i])
            msgs <- list(user_msg)
            if (!is.null(systems)) {
              sys_val <- str_trim(systems[i])
              if (nzchar(sys_val)) {
                msgs <- append(list(list(role = "system", content = sys_val)), msgs)
              }
            }

            payload <- list(
              model = model_id,
              messages = msgs,
              max_tokens = max_new_tokens
            )
            if (!is.null(temperature)) {
              payload$temperature <- temperature
            }

            resp <- post_with_retries(jsonlite::toJSON(payload, auto_unbox = TRUE))
            txt <- httr::content(resp, as = "text", encoding = "UTF-8")
            parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)

            if (is.list(parsed$choices) &&
                length(parsed$choices) >= 1 &&
                is.list(parsed$choices[[1]]) &&
                !is.null(parsed$choices[[1]]$message) &&
                !is.null(parsed$choices[[1]]$message$content)) {
              generations[i] <- as.character(parsed$choices[[1]]$message$content)
            } else if (!is.null(parsed$generated_text)) {
              generations[i] <- as.character(parsed$generated_text)
            } else {
              generations[i] <- NA_character_
            }

            setTxtProgressBar(pb, i)
          }

          close(pb)

          out <- bind_cols(df, tibble(generation = generations))
          return(out)
        }

    #### embedAminoAcids

        #' Embed amino acid sequences with a protein language model
        #'
        #' @param amino_acid_stringset Not yet documented.
        #' @param dataframe Not yet documented.
        #' @param biolm_api_key Not yet documented.
        #' @param nvidia_api_key Not yet documented.
        #' @param platform Not yet documented.
        #' @param input_type Not yet documented.
        #' @param model_name Not yet documented.
        #' @param seq_column Not yet documented.
        #' @param local_script For `platform = "local"`: path to a Python script that embeds the
        #'   sequences in a FASTA file and writes one row per sequence to a CSV file. Required for that platform.
        #' @param python_path For `platform = "local"`: the Python executable to run `local_script` with.
        #' @export
        embedAminoAcids <- function(
                amino_acid_stringset = NULL,
                dataframe = NULL,
                biolm_api_key = NULL,
                nvidia_api_key = NULL,
                platform = c("nvidia", "biolm", "local"),
                input_type = c("amino_acid_stringset", "dataframe"),
                model_name = c("esm2-650m", "esm2-8m", "esm2-35m", "esm2-150m", "esm2_t30_150M_UR50D"),
                seq_column = NULL,
                local_script = NULL,
                python_path = "python3"
            ) {
              # Select first option for parameters
              platform <- platform[1]
              input_type <- input_type[1]
              model_name <- model_name[1]
              
              # Validate that the provided model_name is allowed for the chosen platform
              if (platform == "biolm") {
                allowed_models <- c("esm2-650m", "esm2-8m", "esm2-35m", "esm2-150m")
                if (!model_name %in% allowed_models) {
                  stop("For platform 'biolm', model_name must be one of: ", paste(allowed_models, collapse = ", "))
                }
              }
              if (platform == "nvidia") {
                if (model_name != "esm2-650m") {
                  stop("For platform 'nvidia', only model_name 'esm2-650m' is allowed.")
                }
              }
              if (platform == "local") {
                if (model_name != "esm2_t30_150M_UR50D") {
                  stop("For platform 'local', only model_name 'esm2_t30_150M_UR50D' is allowed.")
                }
              }
              
              # Trim any accidental whitespace/newlines in provided keys
              if (!is.null(biolm_api_key)) biolm_api_key <- trimws(biolm_api_key[1])
              if (!is.null(nvidia_api_key)) nvidia_api_key <- trimws(nvidia_api_key[1])
              
              # Validate input type
              if (input_type == "amino_acid_stringset") { 
                if (!inherits(amino_acid_stringset, "XStringSet")) {
                  stop("For input_type 'amino_acid_stringset', amino_acid_stringset must be of class 'XStringSet'.")
                }
                # Convert the XStringSet to a character vector and preserve names
                sequences <- as.character(amino_acid_stringset)
                seq_names <- names(amino_acid_stringset)
              }
              
              # BIOLM    
              if (platform == "biolm") {
                if (is.null(biolm_api_key) || biolm_api_key == "") {
                  stop("biolm_api_key is missing or empty.")
                }
                # Build the items string manually by iterating over the sequences.
                # If needed, escape any internal quotes in the sequence.
                items_str <- paste0(
                  sapply(sequences, function(seq) {
                    seq_escaped <- gsub("\"", "\\\\\"", seq)
                    sprintf("{\"sequence\": \"%s\"}", seq_escaped)
                  }),
                  collapse = ","
                )
                
                # Create the final JSON payload string.
                # The formatting (line breaks and indentations) here matches the manual example.
                payload <- sprintf("{\n    \"params\": {\n        \"include\": [\n            \"mean\",\n            \"contacts\",\n            \"logits\",\n            \"attentions\"\n        ]\n    },\n    \"items\": [\n        %s\n    ]\n}", items_str)
                
                headers <- c(
                  'Authorization' = paste('Token', biolm_api_key),
                  'Content-Type'  = 'application/json',
                  'Accept'        = 'application/json'
                )
                
                # Construct the API endpoint URL dynamically based on model_name
                biolm_url <- sprintf("https://biolm.ai/api/v3/%s/encode/", model_name)
                response_content <- RCurl::postForm(biolm_url, .opts = list(postfields = payload, httpheader = headers, followlocation = TRUE), style = "httppost")
                
                # Parse the JSON response
                result_data <- jsonlite::fromJSON(response_content)
                
                # Check for any error messages in the API response
                if (!is.null(result_data$error)) {
                  stop("API Error: ", result_data$error)
                }
                
                embeddings <- list()
                for (i in 1:length(result_data$results$embeddings)) { # i=1
                  embeddings[[i]] <- result_data$results$embeddings[[i]]$embedding[[1]]
                }
                embeddings <- as_tibble(as.data.frame(do.call(rbind, embeddings)))
                colnames(embeddings) <- paste0("embedding_", seq(1, dim(embeddings)[2], 1))
                embeddings <- as_tibble(as.data.frame(cbind(data.frame(name = seq_names), embeddings)))
                return(embeddings)
              }
              
              if (platform == "nvidia") {
                if (is.null(nvidia_api_key) || nvidia_api_key == "") {
                  stop("nvidia_api_key is missing or empty.")
                }
                # NVIDIA API accepts a list of sequences
                payload <- list(
                  sequences = as.list(as.data.frame(amino_acid_stringset)[[1]]),
                  format = "h5"
                )
                response <- httr::POST(
                  url = "https://health.api.nvidia.com/v1/biology/meta/esm2-650m",
                  httr::add_headers(
                    `Content-Type` = "application/json",
                    Authorization = paste("Bearer", nvidia_api_key)
                  ),
                  body = jsonlite::toJSON(payload, auto_unbox = TRUE),
                  encode = "json"
                )
                code <- httr::status_code(response)
                if (code >= 300) {
                  body_txt <- httr::content(response, as = "text", encoding = "UTF-8")
                  if (code == 401 || code == 403) {
                    stop("NVIDIA API request failed (auth): HTTP ", code, " - check that nvidia_api_key is valid and not expired. Response: ", body_txt)
                  }
                  stop("NVIDIA API request failed: HTTP ", code, " Response: ", body_txt)
                }

                # Save response to a temporary file; NVIDIA returns an HDF5 file, not a zip
                file_path <- tempfile(fileext = ".h5")
                resp_raw <- httr::content(response, as = "raw")
                if (length(resp_raw) == 0) {
                  stop("NVIDIA API returned an empty response.")
                }
                writeBin(resp_raw, file_path)

                # Some future versions could wrap the H5 in a zip; handle both
                read_path <- file_path
                ctype <- httr::headers(response)[["content-type"]]
                if (!is.null(ctype) && grepl("zip", ctype, ignore.case = TRUE)) {
                  temp_dir <- tempdir()
                  unzipped_files <- unzip(file_path, exdir = temp_dir)
                  if (length(unzipped_files) == 0) stop("Failed to unzip NVIDIA embeddings payload.")
                  read_path <- unzipped_files[1]
                }

                embeddings <- as_tibble(t(rhdf5::h5read(read_path, "embeddings")))
                colnames(embeddings) <- paste0("embedding_", seq(1, dim(embeddings)[2], 1))
                embeddings <- as_tibble(as.data.frame(cbind(data.frame(name = seq_names), embeddings)))
                return(embeddings)
              }
              if (platform == "local") {
                # Runs a local encoder script as: python_path local_script input.fasta output.csv model_name
                # The script writes one row of embedding values per sequence to output.csv.
                if (is.null(local_script) || !nzchar(local_script)) {
                  stop("platform = 'local' needs local_script = the path to a protein-encoding Python script ",
                       "(called as: python local_script input.fasta output.csv model_name).", call. = FALSE)
                }
                if (!file.exists(local_script)) stop("local_script not found: ", local_script, call. = FALSE)
                fasta_file <- tempfile(fileext = ".fasta")
                output_file <- tempfile(fileext = ".csv")
                on.exit(unlink(c(fasta_file, output_file)), add = TRUE)
                Biostrings::writeXStringSet(amino_acid_stringset, fasta_file)
                system2(python_path, c(shQuote(local_script), shQuote(fasta_file), shQuote(output_file), model_name))
                if (!file.exists(output_file)) stop("The local script did not produce an output file.", call. = FALSE)
                embeddings <- utils::read.csv(output_file)
                embeddings <- as_tibble(as.data.frame(cbind(data.frame(name = seq_names), embeddings)))
                return(embeddings)
              }
            } 

