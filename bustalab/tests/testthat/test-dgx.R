## No dgx server is reachable in tests: every endpoint is pointed at a closed local
## port, so each request is refused at once rather than timing out.

local_dead_dgx <- function(env = parent.frame()) {
    withr::local_envvar(c(DGX_GENERATE_URL = "http://127.0.0.1:9",
                          DGX_ESM2_URL = "http://127.0.0.1:9",
                          DGX_ESMC_URL = "http://127.0.0.1:9",
                          DGX_MAX_GRID_CELLS = NA), .local_envir = env)
}

test_that("the dgx public API is exported with its frozen signatures", {
    ns <- asNamespace("bustalab")
    expect_true(all(c("dgxHealth", "dgxGenerate", "dgxGridStatus", "runModelGrid",
                      "dgxMaxGridCells") %in% getNamespaceExports("bustalab")))
    expect_false(any(c(".dgxServerIdle", ".dgxPreflight") %in% getNamespaceExports("bustalab")))
    expect_identical(names(formals(get("dgxGenerate", ns))),
                     c("prompt", "model", "quant", "max_tokens", "system", "logprobs",
                       "server_url", "timeout_sec", "max_wait_sec", "on_wait"))
    expect_identical(names(formals(get("runModelGrid", ns))),
                     c("prompts", "models", "quants", "out_csv", "max_tokens", "system",
                       "labels", "verbose", "force", "max_consecutive_errors", "logprobs"))
})

test_that("dgxMaxGridCells reads its environment variable", {
    local_dead_dgx()
    expect_equal(dgxMaxGridCells(), 5000)
    withr::local_envvar(DGX_MAX_GRID_CELLS = "12")
    expect_equal(dgxMaxGridCells(), 12)
})

test_that("dgxHealth and dgxGridStatus report an unreachable server without erroring", {
    local_dead_dgx()
    msgs <- character(0)
    res <- withCallingHandlers(dgxHealth(), message = function(m) {
        msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
    expect_null(res)
    expect_equal(sum(grepl("UNREACHABLE", msgs)), 3)
    expect_message(dgxGridStatus(), "could not reach")
})

test_that("runModelGrid refuses to start against an unreachable server", {
    local_dead_dgx()
    out <- tempfile(fileext = ".csv")
    expect_error(runModelGrid("hello", models = "m", out_csv = out, verbose = FALSE),
                 "Cannot reach the generate server")
    expect_false(file.exists(out))
})

test_that("runModelGrid refuses an oversized grid", {
    local_dead_dgx()
    withr::local_envvar(DGX_MAX_GRID_CELLS = "2")
    expect_error(runModelGrid(c("a", "b", "c"), models = "m", out_csv = tempfile(), verbose = FALSE),
                 "This grid is 3 generations")
})

test_that("dgxGenerate fails with an error when nothing is listening", {
    local_dead_dgx()
    expect_error(dgxGenerate("hello", timeout_sec = 5, max_wait_sec = 0))
})
