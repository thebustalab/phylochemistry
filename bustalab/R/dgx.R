#################################################
## DGX ENDPOINTS (research students, JupyterHub) ##
#################################################

# Split out of modules/language_model_analysis.R on 2026-09-08.
#
# WHY THIS IS ITS OWN MODULE: language_model_analysis.R is sourced ONLY on the lab
# build (bustalab = TRUE). A research student sourcing phylochemistry.R therefore
# never got dgxGenerate/runModelGrid/dgxHealth, and dgx_example_R.ipynb failed on
# its first cell with `could not find function "dgxHealth"`. This file is sourced
# unconditionally from phylochemistry.R, and it deliberately needs only packages
# that are already in the CORE scope (httr, jsonlite, readr, dplyr, tibble) so it
# adds no install burden for a chapters-1-11 student.

message("Loading dgx module...")

        # These reach the lab's dgx through SSH tunnels that already exist on host1:
        #   9000 = ESM2 protein embeddings, 9001 = ESMC, 9002 = Gemma text generation.
        # They are NOT reachable from a laptop — only from inside JupyterHub.
        #
        # THE ONE THING TO UNDERSTAND: the dgx holds ONE model at a time. Asking for a
        # different model (or a different quantisation) evicts the old one and loads the new,
        # which takes about 3 minutes. Answering a prompt then takes under a second. So never
        # loop model-inside-prompt; use runModelGrid(), which groups by model.

        dgxGenerate <- function(
          prompt,
          model = NULL,          # e.g. "google/gemma-4-12B-it"; NULL = server default
          quant = NULL,          # 4 = nf4, 8 = int8, 16 = bf16 (full precision)
          max_tokens = 256,
          system = NULL,
          # Opt-in token probabilities from the model's RAW distribution (before temperature or
          # top_p). 1 = log-probability of each generated token; 2 or 3 = that plus the 1-2 most
          # likely other tokens at each step (3 = chosen + next two). The result gains a
          # `logprobs` data frame (one row per token: token, id, logprob, and for 2-3 a `top`
          # list-column). Needs max_tokens <= 512; for classification use max_tokens of 1-4.
          logprobs = NULL,
          server_url = Sys.getenv("DGX_GENERATE_URL", "http://127.0.0.1:9002"),
          timeout_sec = 900,     # a cold model load is ~180s; 60 would fail on every first call
          # If the dgx is too full to load the model right now (e.g. a big model while a video
          # render is running), it answers 503 "memory_admission". That means "not yet", so
          # WAIT and retry, for up to this long. 0 = fail straight away. Added 2026-09-18.
          max_wait_sec = as.numeric(Sys.getenv("DGX_MEMORY_WAIT_SEC", 4 * 3600)),
          on_wait = NULL         # internal: called before and after each memory wait
        ) {
          suppressPackageStartupMessages({ library(httr); library(jsonlite) })
          payload <- list(prompt = prompt, max_tokens = max_tokens)
          if (!is.null(model))  payload$model  <- model
          if (!is.null(quant))  payload$quant  <- as.integer(quant)
          if (!is.null(system)) payload$system <- system
          if (!is.null(logprobs) && logprobs > 0) payload$logprobs <- as.integer(logprobs)
          waited <- 0
          repeat {
            resp <- httr::POST(
              url = paste0(sub("/+$", "", server_url), "/generate"),
              httr::add_headers(`Content-Type` = "application/json",
                                `User-Agent` = "Mozilla/5.0 (phylochemistry dgxGenerate)"),
              body = jsonlite::toJSON(payload, auto_unbox = TRUE),
              encode = "raw", httr::timeout(timeout_sec)
            )
            txt <- httr::content(resp, "text", encoding = "UTF-8")
            if (httr::status_code(resp) != 503 || !grepl("memory_admission", txt, fixed = TRUE)) break
            wait <- suppressWarnings(as.numeric(httr::headers(resp)[["retry-after"]]))
            if (length(wait) != 1 || is.na(wait)) wait <- 300
            why <- tryCatch(jsonlite::fromJSON(txt)$detail, error = function(e) txt)
            why <- sub("memory_admission: ", "", why, fixed = TRUE)
            if (waited + wait > max_wait_sec)
              stop("Gave up after waiting ", round(waited / 60), " min for memory: ", why, call. = FALSE)
            if (waited == 0) {
              message("The dgx is too full to load ", model, " @ ", quant, "-bit right now:\n  ", why,
                      "\nWaiting and retrying every ", round(wait / 60), " min (up to ",
                      round(max_wait_sec / 3600), " h). Nothing is lost.")
            } else {
              message("  still waiting for memory (", round(waited / 60), " min so far)...")
            }
            if (is.function(on_wait)) on_wait()
            Sys.sleep(wait)
            waited <- waited + wait
            if (is.function(on_wait)) on_wait()
          }
          if (httr::status_code(resp) >= 300) {
            stop("dgx said ", httr::status_code(resp), ": ",
                 substr(httr::content(resp, "text", encoding = "UTF-8"), 1, 300))
          }
          jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
        }

        # ---- SHARING: the dgx has ONE model slot for the whole lab ------------------------
        # Calls on port 9002 go through the dgx broker, which queues them and keeps a loaded
        # model serving its work before swapping, so two grids take turns instead of evicting
        # each other on every prompt. runModelGrid() still refuses an unreachable server or an
        # absurdly large grid, and stops a run failing every prompt. These guards mirror
        # dgx_client.py exactly — keep the two sides in step.

        dgxMaxGridCells  <- function() as.numeric(Sys.getenv("DGX_MAX_GRID_CELLS", "5000"))

        # Seconds since the generate server last answered anyone; NA if we cannot reach it.
        .dgxServerIdle <- function(
          server_url = Sys.getenv("DGX_GENERATE_URL", "http://127.0.0.1:9002")
        ) {
          tryCatch({
            resp <- httr::GET(paste0(sub("/+$", "", server_url), "/health"), httr::timeout(10))
            if (httr::status_code(resp) != 200) return(NA_real_)
            idle <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))$idle_sec
            if (is.null(idle)) NA_real_ else as.numeric(idle)
          }, error = function(e) NA_real_)
        }

        # Can I start a grid right now? The broker queues everyone's calls, so another grid
        # running is not a reason to wait: yours takes its turn.
        dgxGridStatus <- function() {
          idle <- .dgxServerIdle()
          if (is.na(idle)) message("Server: could not reach the generate endpoint.")
          else if (idle < 300)
            message("Server: up, last answered a call ", round(idle), " s ago. Calls are queued, ",
                    "so you can start anyway — yours will take its turn.")
          else message("Server: up and quiet — idle for ", round(idle / 60), " min.")
          invisible(NULL)
        }

        # Refuse to start a grid that would monopolise the box or run against a dead server.
        .dgxPreflight <- function(total, force, verbose) {
          problems <- character(0)
          if (total > dgxMaxGridCells()) problems <- c(problems, paste0(
            "This grid is ", format(total, big.mark = ","), " generations. That is days of GPU\n",
            "    time for the whole lab. Pilot on 20-50 prompts across the full grid first,\n",
            "    decide which (model, quant) actually works, then run everything through that\n",
            "    ONE cell."))
          if (is.na(.dgxServerIdle())) problems <- c(problems, paste0(
            "Cannot reach the generate server. Run dgxHealth() first — starting a grid\n",
            "    against a dead server just fills the CSV with errors."))
          if (verbose) { message("Pre-flight check (same as dgxGridStatus()):"); dgxGridStatus() }
          if (!length(problems)) return(invisible(NULL))
          bar <- strrep("-", 90)
          body <- paste0("  * ", problems, collapse = "\n\n")
          if (force) {
            if (verbose) message("\n", bar, "\nrunModelGrid WARNINGS (overridden by force = TRUE):\n\n",
                                 body, "\n", bar)
            return(invisible(NULL))
          }
          stop("\n", bar, "\nrunModelGrid refused to start:\n\n", body,
               "\n\nIf you are sure, call runModelGrid(..., force = TRUE).\n", bar, call. = FALSE)
        }

        # Run every prompt against every (model, quantisation) and write a tidy CSV.
        # Batched BY CELL so the ~3-minute load is paid once per cell, not once per prompt.
        # RESUMABLE: appends after every prompt and skips rows that SUCCEEDED, so if the
        # kernel dies four hours in, re-running the same cell continues where it stopped.
        # Rows that FAILED are retried — a failure is not a result. A retried prompt appends a
        # second row; the later row is the one that counts.
        # IT WILL REFUSE TO START if the server is unreachable or the grid is larger than
        # DGX_MAX_GRID_CELLS (pass force = TRUE if you are sure about the size). Someone else's
        # grid running is fine: the dgx broker queues both and they take turns.
        runModelGrid <- function(
          prompts,
          models = c("google/gemma-4-E2B-it", "google/gemma-4-E4B-it",
                     "google/gemma-4-12B-it", "google/gemma-4-26B-A4B-it",
                     "google/gemma-4-31B-it"),
          quants = c(4),
          out_csv = "grid_results.csv",
          max_tokens = 16,
          system = NULL,
          labels = NULL,
          verbose = TRUE,
          force = FALSE,
          max_consecutive_errors = as.numeric(Sys.getenv("DGX_MAX_CONSECUTIVE_ERRORS", "10")),
          # Opt-in, 1-3 (see dgxGenerate). Adds two columns: `logprobs` (a JSON string, one
          # entry per generated token) and `first_token_logprob` (handy for one-token
          # classification). Changing it needs a fresh out_csv, since the columns differ.
          logprobs = NULL
        ) {
          suppressPackageStartupMessages({ library(readr); library(dplyr); library(tibble) })
          if (is.null(labels)) labels <- as.character(seq_along(prompts))
          if (length(labels) != length(prompts)) stop("labels must match prompts in length")

          total <- length(models) * length(quants) * length(prompts)
          .dgxPreflight(total, force, verbose)

          use_lp <- !is.null(logprobs) && logprobs > 0
          # A row counts as done only if it has no error. The LAST row for a key wins, so a
          # successful retry supersedes the failure that came before it.
          done_keys <- character(0); retrying <- 0
          if (file.exists(out_csv)) {
            prev <- suppressMessages(readr::read_csv(out_csv, show_col_types = FALSE))
            if (use_lp != ("logprobs" %in% names(prev)))
              stop(out_csv, " has different columns than this run (logprobs on/off changes them). ",
                   "Use a new out_csv.", call. = FALSE)
            if (nrow(prev) > 0) {
              prev <- prev %>%
                mutate(.key = paste(model, quant, as.character(prompt_id)),
                       .ok  = is.na(error) | !nzchar(trimws(as.character(error))),
                       .row = row_number()) %>%
                group_by(.key) %>% filter(.row == max(.row)) %>% ungroup()
              done_keys <- prev$.key[prev$.ok]
              retrying  <- sum(!prev$.ok)
              if (verbose) {
                message("Resuming: ", length(done_keys), " successful results already in ", out_csv)
                if (retrying > 0) message("          ", retrying, " failed rows will be RETRIED.")
              }
            }
          }

          n <- length(done_keys); consecutive_errors <- 0

          for (m in models) for (q in quants) {
            keep <- !(paste(m, q, labels) %in% done_keys)
            if (!any(keep)) { if (verbose) message("[skip] ", m, " @ ", q, "-bit"); next }
            if (verbose) message("\n[cell] ", m, " @ ", q, "-bit — ", sum(keep),
                                 " prompts\n       loading the model, ~3 minutes...")
            for (i in which(keep)) {
              failed <- FALSE
              row <- tryCatch({
                r <- dgxGenerate(prompts[i], model = m, quant = q,
                                 max_tokens = max_tokens, system = system,
                                 logprobs = if (use_lp) logprobs else NULL)
                tibble(model = m, quant = q, prompt_id = labels[i], prompt = prompts[i],
                       completion = r$completion, duration_sec = r$duration_sec,
                       error = NA_character_, timestamp = format(Sys.time())) %>%
                  { if (use_lp) {
                      lp <- r$logprobs
                      mutate(., logprobs = if (length(lp)) as.character(jsonlite::toJSON(lp, auto_unbox = TRUE)) else "[]",
                             first_token_logprob = if (length(lp)) lp$logprob[1] else NA_real_)
                    } else . }
              }, error = function(e) {    # one bad prompt must not kill a long run
                failed <<- TRUE
                tibble(model = m, quant = q, prompt_id = labels[i], prompt = prompts[i],
                       completion = NA_character_, duration_sec = NA_real_,
                       error = substr(conditionMessage(e), 1, 300), timestamp = format(Sys.time())) %>%
                  { if (use_lp) mutate(., logprobs = NA_character_, first_token_logprob = NA_real_) else . }
              })
              readr::write_csv(row, out_csv, append = file.exists(out_csv))
              n <- n + 1
              consecutive_errors <- if (failed) consecutive_errors + 1 else 0

              if (consecutive_errors >= max_consecutive_errors) {
                message("\nSTOPPED after ", consecutive_errors, " failures in a row. ",
                        "The last one was:\n    ", row$error[1],
                        "\nNothing is being generated, so there is no point continuing — the ",
                        "rest of the grid\nwould just fill ", out_csv, " with the same error.",
                        "\nCheck the server with dgxHealth(), then re-run this same cell: the ",
                        "failed rows\nwill be retried automatically.")
                stop(consecutive_errors, " consecutive failures — run stopped. Last error: ",
                     row$error[1], call. = FALSE)
              }
              if (verbose && (n %% 10 == 0)) message("       ", n, "/", total, " overall")
            }
          }
          if (verbose) message("\nDone. Results in ", out_csv)
          invisible(out_csv)
        }

        dgxHealth <- function() {
          suppressPackageStartupMessages(library(httr))
          urls <- c(
            generate = Sys.getenv("DGX_GENERATE_URL", "http://127.0.0.1:9002"),
            esm2     = Sys.getenv("DGX_ESM2_URL",     "http://127.0.0.1:9000"),
            esmc     = Sys.getenv("DGX_ESMC_URL",     "http://127.0.0.1:9001")
          )
          # Loop over NAMES, not values: `for (nm in urls)` binds nm to the URL and
          # loses the label, so every line printed the address instead of the service.
          for (nm in names(urls)) {
            ok <- tryCatch(httr::status_code(httr::GET(paste0(sub("/+$", "", urls[[nm]]), "/health"),
                                                       httr::timeout(10))) == 200,
                           error = function(e) FALSE)
            message(nm, " (", urls[[nm]], "): ",
                    if (ok) "ok" else "UNREACHABLE (are you on JupyterHub?)")
          }
          invisible(NULL)
        }

message("Done with dgx loading!")
