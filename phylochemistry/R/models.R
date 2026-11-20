## models — moved verbatim from phylochemistry.R (package conversion, commit 1)


        #### buildModel

            # Was buildModel2 until 2026-10-30 (Lucas): the older buildModel it replaced had no callers and was
            # removed. buildModel2 stays as an alias below so existing code keeps running.
            buildModel <- function(
                data,
                model_type = c(
                    "linear_regression", "random_forest_regression",
                    "random_forest_classification", "logistic_regression",
                    "contrastive_learning"
                ),
                predictors = NULL,
                response = NULL,
                fold_cross_validation = 10,
                optimization_parameters = list(
                    n_vars_tried_at_split = seq(10,20,5),
                    n_trees = seq(100,200,50),
                    min_leaf_size = 1
                ),
                formula = NULL,
                scale_variance = FALSE,
                input_variables = NULL,
                output_variable = NULL
            ) {

                ## Prepare the data

                    ## predictors / response were input_variables / output_variable until 2026-11-07
                    ## (Lucas). The old names stay accepted, silently and last in the argument list, so
                    ## saved code keeps running; they are resolved here once and never used again.
                        if (!is.null(input_variables)) {
                            if (!is.null(predictors)) stop("Please give the predictors only once.")
                            predictors <- input_variables
                        }
                        if (!is.null(output_variable)) {
                            if (!is.null(response)) stop("Please give the response only once.")
                            response <- output_variable
                        }

                    ## Correct the formula and the data, start the output
                        model_type <- match.arg(model_type)

                        ## formula = is the alternative to predictors / response (linear regression
                        ## only), e.g. formula = rate ~ temperature * chemistry. It comes after the other
                        ## arguments so positional calls are unaffected. The variables it names are worked
                        ## out here so every check below runs unchanged.
                        user_formula <- formula
                        if (!is.null(user_formula)) {
                            if (is.character(user_formula)) user_formula <- stats::as.formula(user_formula)
                            if (!inherits(user_formula, "formula") || length(user_formula) != 3) {
                                stop("formula should look like response ~ predictors, for example rate ~ temperature * chemistry.")
                            }
                            if (model_type != "linear_regression") {
                                stop("formula = works with model_type = \"linear_regression\" only. For other model types, use predictors and response.")
                            }
                            if (!is.null(predictors) || !is.null(response)) {
                                stop("Give either formula, or predictors and response, not both.")
                            }
                            if (!is.name(user_formula[[2]])) {
                                stop("The left side of the formula should be one column name, not a calculation. Make the column first, e.g. mutate(log_rate = log(rate)), then use log_rate ~ ...")
                            }
                            response <- all.vars(user_formula[[2]])
                            predictors <- all.vars(user_formula[[3]])
                            if (length(response) != 1) {
                                stop("The left side of the formula should be the one column you want to predict.")
                            }
                            if ("." %in% predictors) {
                                stop("Please name the predictor columns in the formula rather than using a dot.")
                            }
                            missing_vars <- setdiff(c(response, predictors), colnames(data))
                            if (length(missing_vars) > 0) {
                                stop(paste("The following formula variables are not present in the data:", paste(missing_vars, collapse = ", ")))
                            }
                            ## only the columns the model uses are checked for missing values and fitted
                            data <- as.data.frame(data)[, c(response, predictors), drop = FALSE]
                            ## text predictors become factors so later tests (tukeyTest) see the groups
                            for (v in predictors) {
                                if (is.character(data[[v]])) data[[v]] <- factor(data[[v]])
                            }
                        }

                        if (is.null(predictors) || length(predictors) == 0) {
                            stop("Please provide at least one predictor.")
                        }
                        if (is.null(response) || !nzchar(response)) {
                            stop("Please provide a valid response.")
                        }
                        if (!response %in% colnames(data)) {
                            stop("The specified response is not present in the data.")
                        }
                        if (!all(predictors %in% colnames(data))) {
                            missing_predictors <- predictors[!predictors %in% colnames(data)]
                            stop(
                                paste(
                                    "The following predictors are not present in the data:",
                                    paste(missing_predictors, collapse = ", ")
                                )
                            )
                        }
                        if (!is.numeric(fold_cross_validation) || fold_cross_validation < 2) {
                            stop("fold_cross_validation must be a numeric value greater than or equal to 2.")
                        }
                        if (model_type %in% c("random_forest_regression", "random_forest_classification")) {
                            required_optimization_args <- c("n_vars_tried_at_split", "n_trees", "min_leaf_size")
                            if (!all(required_optimization_args %in% names(optimization_parameters))) {
                                stop("optimization_parameters must include n_vars_tried_at_split, n_trees, and min_leaf_size.")
                            }
                            if (any(lengths(optimization_parameters[required_optimization_args]) == 0)) {
                                stop("optimization_parameters entries must contain at least one value.")
                            }
                        }
                        data <- as.data.frame(data)
                        if( any(is.na(data)) ) { stop("Yo - there are missing data in your data set, please deal with that before proceeding.") }
                        ## scale_variance = TRUE: centre each numeric predictor on its mean and divide it by its standard
                        ## deviation, as runMatrixAnalysis does, so the coefficients can be compared across predictors.
                        ## The response column is never scaled, so predictions stay in its units.
                        if (!is.logical(scale_variance) || length(scale_variance) != 1 || is.na(scale_variance)) {
                            stop("scale_variance should be TRUE or FALSE.")
                        }
                        if (isTRUE(scale_variance)) {
                            if (!is.null(user_formula)) {
                                ## a named function in the formula (log, I, factor, poly) would be applied to the
                                ## SCALED column, so refuse it; operators such as + * : ^ are fine
                                if (length(grep("^[A-Za-z._]", setdiff(all.names(user_formula[[3]]), all.vars(user_formula[[3]])))) > 0) {
                                    stop("scale_variance = TRUE needs plain column names in the formula (for example y ~ a + b). Make any calculated column with mutate() first.")
                                }
                            }
                            scaled_inputs <- setdiff(predictors, response)
                            scaled_inputs <- scaled_inputs[vapply(scaled_inputs, function(v) is.numeric(data[[v]]), logical(1))]
                            if (length(scaled_inputs) == 0) {
                                stop("scale_variance = TRUE scales numeric predictors, but none of the predictors are numbers.")
                            }
                            input_means <- vapply(scaled_inputs, function(v) mean(data[[v]]), numeric(1))
                            input_sds <- vapply(scaled_inputs, function(v) stats::sd(data[[v]]), numeric(1))
                            if (any(is.na(input_sds) | input_sds == 0)) {
                                stop(paste0("scale_variance = TRUE cannot scale ", paste(scaled_inputs[is.na(input_sds) | input_sds == 0], collapse = ", "), ": every value in it is the same."))
                            }
                            for (v in scaled_inputs) data[[v]] <- as.numeric((data[[v]] - input_means[[v]]) / input_sds[[v]])
                            predictor_scaling <- data.frame(variable = scaled_inputs, mean = unname(input_means), sd = unname(input_sds))
                        }
                        formula <- paste0("`", response, "` ~ `", paste(predictors, collapse = "` + `"), "`")
                        if (!is.null(user_formula)) formula <- user_formula
                        output <- list()
                        output$model_type <- model_type
                        output$response <- response   # predictWithModel names its column from this
                        output$predictors <- predictors
                        if (isTRUE(scale_variance)) output$predictor_scaling <- predictor_scaling
                        ## only on the formula path, so the variable-list output gains no extra parts
                        if (!is.null(user_formula)) {
                            output$formula <- user_formula
                            output$data <- data
                        }

                    ## Various checkes
                            if (dim(data)[1] < length(predictors)) {
                                stop("Heyyyyy, there are fewer observations in the data than predictors. You should use fewer predictors or make more observations.")
                            }
                            if (!is.numeric(data[[response]]) && !model_type %in% c("random_forest_classification", "logistic_regression")) {
                                stop("Uh oh! You cannot predict the value of a categorical variable with a regression model. Choose a classification model instead.")
                            }
                            if (model_type == "logistic_regression") {
                                outcome_values <- data[[response]]
                                if (!(is.character(outcome_values) || is.factor(outcome_values) || is.logical(outcome_values))) {
                                    stop(paste0("Uh oh! Logistic regression predicts a yes/no response, but ", response, " is numbers. Make it a yes/no column first, for example mutate(", response, " = ifelse(", response, " > 0, \"yes\", \"no\")), or use linear_regression to predict a number."))
                                }
                                classes <- levels(droplevels(factor(outcome_values)))
                                if (length(classes) != 2) {
                                    stop(paste0("Uh oh! Logistic regression needs a response with exactly two values (for example \"yes\" and \"no\"), but ", response, " has ", length(classes), ": ", paste(classes, collapse = ", "), "."))
                                }
                                ## The positive class (the one whose probability is predicted): "yes" for a
                                ## yes/no column (any capitalisation), TRUE for a logical column, otherwise the
                                ## second factor level (for text, the second in alphabetical order).
                                positive <- if (setequal(tolower(classes), c("yes", "no"))) classes[tolower(classes) == "yes"] else classes[2]
                                ## glm() models the probability of the SECOND level, and predictWithModel reads
                                ## the two levels back from the fitted model, so the order is fixed here.
                                data[[response]] <- factor(as.character(outcome_values), levels = c(setdiff(classes, positive), positive))
                            }
                            if (model_type == "random_forest_classification" && !is.factor(data[[response]])) {
                                data[[response]] <- as.factor(data[[response]])
                            }

                ## Modeling

                    if (model_type == "linear_regression") {

                        ## Run the fit and start the output
                            fit <- lm(
                                data = data,
                                formula = formula,
                                x = TRUE, y = TRUE, model = TRUE, qr = TRUE
                            )
                            output$model <- fit

                            input_y <- fit$y[order(as.numeric(names(fit$y)))]
                            total_sum_squares <- sum((input_y - mean(input_y, na.rm = TRUE))^2, na.rm = TRUE)
                            residual_sum_squares <- sum((summary(fit)$residuals)^2, na.rm = TRUE)

                            metrics <- rbind(
                                data.frame(
                                    variable = c("r_squared", "total_sum_squares", "residual_sum_squares"),
                                    value = c(summary(fit)$r.squared, total_sum_squares, residual_sum_squares),
                                    std_err = NA,
                                    type = "statistic",
                                    p_value = NA,
                                    p_value_adj = NA
                                ),
                                data.frame(
                                    variable = names(coefficients(fit)),
                                    value = as.numeric(coefficients(fit)),
                                    std_err = as.numeric(summary(fit)$coefficients[,2]),
                                    type = "coefficient",
                                    p_value = summary(fit)$coefficients[,4],
                                    p_value_adj = p.adjust(summary(fit)$coefficients[,4])
                                )
                            )
                            # Full precision: rounding to 4 dp collapsed small slopes (~1e-4) and their
                            # standard errors to 1e-04 or 0. Round at display time instead.
                            metrics$value <- as.numeric(metrics$value)
                            rownames(metrics) <- NULL
                            output$metrics <- metrics
                    }

                    if (model_type == "logistic_regression") {

                        ## Run the fit and start the output. Coefficients are on the log-odds scale of
                        ## the positive class; after scaling the predictors, larger absolute values still
                        ## mean a predictor matters more, as in linear regression.
                            fit <- glm(
                                data = data,
                                formula = formula,
                                family = binomial,
                                x = TRUE, y = TRUE, model = TRUE
                            )
                            output$model <- fit

                            metrics <- rbind(
                                data.frame(
                                    variable = c("null_deviance", "residual_deviance", "aic"),
                                    value = c(fit$null.deviance, fit$deviance, fit$aic),
                                    std_err = NA,
                                    type = "statistic",
                                    p_value = NA,
                                    p_value_adj = NA
                                ),
                                data.frame(
                                    variable = names(coefficients(fit)),
                                    value = as.numeric(coefficients(fit)),
                                    std_err = as.numeric(summary(fit)$coefficients[,2]),
                                    type = "coefficient",
                                    p_value = summary(fit)$coefficients[,4],
                                    p_value_adj = p.adjust(summary(fit)$coefficients[,4])
                                )
                            )
                            rownames(metrics) <- NULL
                            output$metrics <- metrics
                    }

                    if (model_type == "random_forest_regression") {

                        ## Create a workflow that has a recipe and a model
                            workflow() %>%
                                add_recipe(
                                    recipe( as.formula(formula), data = data ) # %>%
                                    # step_normalize(all_numeric()) # %>% step_impute_knn(all_inputs())
                                ) %>%
                                add_model(
                                    rand_forest() %>% # specify that the model is a random forest
                                    set_args(mtry = tune(), trees = tune(), min_n = tune()) %>% # specify that the `mtry` and `trees` parameters needs to be tuned
                                    set_engine("ranger", importance = "impurity") %>% # select the engine/package that underlies the model
                                    set_mode("regression") # choose either the continuous regression or binary classification mode
                                ) -> workflow_for_training

                        ## Tune the model on the training set
                            tune_results <- tune_grid(
                                workflow_for_training,
                                resamples = vfold_cv(data, v = fold_cross_validation), #CV object
                                grid = expand.grid(
                                    mtry = optimization_parameters$n_vars_tried_at_split,
                                    trees = optimization_parameters$n_trees,
                                    min_n = optimization_parameters$min_leaf_size
                                ), # grid of values to try
                                metrics = yardstick::metric_set(yardstick::rmse)
                            )

                        # Check model parameters if you want
                            output$metrics <- collect_metrics(tune_results)
                            names(output$metrics)[names(output$metrics) == "mtry"] <- "n_vars_tried_at_split"
                            names(output$metrics)[names(output$metrics) == "trees"] <- "n_trees"
                            names(output$metrics)[names(output$metrics) == "min_n"] <- "min_leaf_size"
                            if ("n" %in% names(output$metrics)) {
                                output$metrics$fold_cross_validation <- output$metrics$n
                                output$metrics$n <- NULL
                            }

                        ## Apply the best parameters to the workflow to creat the output model
                            output$model <- fit(
                                finalize_workflow(workflow_for_training, select_best(tune_results, metric = "rmse")),
                                data
                            )
                    }

                    if (model_type == "random_forest_classification") {

                        ## Create a workflow that has a recipe and a model
                            workflow() %>%
                                add_recipe(
                                    recipe( as.formula(formula), data = data ) # %>%
                                    # step_normalize(all_numeric()) # %>% step_impute_knn(all_inputs())
                                ) %>%
                                add_model(
                                    rand_forest() %>% # specify that the model is a random forest
                                    set_args(mtry = tune(), trees = tune(), min_n = tune()) %>% # specify that the `mtry` and `trees` parameters needs to be tuned
                                    set_engine("ranger", importance = "impurity") %>% # select the engine/package that underlies the model
                                    set_mode("classification") # choose either the continuous regression or binary classification mode
                                ) -> workflow_for_training

                        ## Tune the model on the training set
                            tune_results <- tune_grid(
                                workflow_for_training,
                                resamples = vfold_cv(data, v = fold_cross_validation), #CV object
                                grid = expand.grid(
                                    mtry = optimization_parameters$n_vars_tried_at_split,
                                    trees = optimization_parameters$n_trees,
                                    min_n = optimization_parameters$min_leaf_size
                                ), # grid of values to try
                                metrics = yardstick::metric_set(yardstick::accuracy)
                            )

                        # Check model parameters if you want
                            output$metrics <- collect_metrics(tune_results)
                            names(output$metrics)[names(output$metrics) == "mtry"] <- "n_vars_tried_at_split"
                            names(output$metrics)[names(output$metrics) == "trees"] <- "n_trees"
                            names(output$metrics)[names(output$metrics) == "min_n"] <- "min_leaf_size"
                            if ("n" %in% names(output$metrics)) {
                                output$metrics$fold_cross_validation <- output$metrics$n
                                output$metrics$n <- NULL
                            }
                            # select_best(tune_results, metric = "rmse")

                        ## Apply the best parameters to the workflow to creat the output model
                            output$model <- fit(
                                finalize_workflow(workflow_for_training, select_best(tune_results, metric = "accuracy")),
                                data
                            )
                    }

                ## Return a model trained on all the input data
                    
                    if (isTRUE(scale_variance)) attr(output$model, "predictor_scaling") <- predictor_scaling
                    return(output)
            }



        #### buildModel2 (alias, kept for back-compatibility)

            buildModel2 <- function(...) buildModel(...)


        #### predictWithModel

            predictWithModel <- function(
                data,
                model,
                model_type = NULL,
                id_cols = NULL
            ) {

                ## Lucas, 2026-10-30: the model type comes from the model itself, and the result is the
                ## input data with the prediction(s) appended, like the cluster column runMatrixAnalysis
                ## adds. `model` may be the whole buildModel() result or its $model; `model_type` is an
                ## optional override kept for old code; `id_cols` is accepted and ignored (every input
                ## column now comes back).

                ## Check to see that required data are provided
                    if (missing(data)) {
                        stop("Input data is required.")
                    }
                    if (missing(model)) {
                        stop("Please provide a model (the output of buildModel(), or its $model).")
                    }

                ## Unpack a whole buildModel() result
                    outcome <- NULL
                    if (is.list(model) && !is.null(model$model_type) && !is.null(model$model)) {
                        if (is.null(model_type)) model_type <- model$model_type
                        ## $response; models built before 2026-11-07 carry only $output_variable
                        outcome <- if (!is.null(model$response)) model$response else model$output_variable
                        model <- model$model
                    }
                    predictor_scaling <- attr(model, "predictor_scaling")
                    if (is.null(predictor_scaling)) predictor_scaling <- attr(model, "input_scaling")  # models built before the 2026-11-08 rename

                ## Infer the model type from the fitted object
                    if (is.null(model_type)) {
                        if (inherits(model, "glm")) {
                            model_type <- "logistic_regression"
                        } else if (inherits(model, "lm")) {
                            model_type <- "linear_regression"
                        } else if (inherits(model, "workflow")) {
                            mode <- model$fit$actions$model$spec$mode
                            model_type <- if (identical(mode, "classification")) "random_forest_classification" else "random_forest_regression"
                        } else if (is.list(model) && !is.null(model$architecture)) {
                            model_type <- "contrastive_learning"
                        } else {
                            stop("predictWithModel() could not tell what kind of model this is. Pass the output of buildModel().")
                        }
                    }
                    model_type <- match.arg(model_type, c(
                        "linear_regression", "random_forest_regression",
                        "random_forest_classification", "logistic_regression",
                        "contrastive_learning"
                    ))

                ## The response's name, for the column
                    if (is.null(outcome)) {
                        if (inherits(model, "lm")) outcome <- names(model$model)[1]
                        if (inherits(model, "workflow")) outcome <- names(model$pre$mold$outcomes)[1]
                    }
                    pred_col <- paste0("predicted_", if (is.null(outcome) || !nzchar(outcome)) "value" else outcome)

                ## Predictions, appended to the input data
                    output <- as_tibble(data)
                    if (!is.null(predictor_scaling)) {
                        ## built with scale_variance = TRUE: scale the new predictors with the TRAINING means and
                        ## standard deviations (never their own), so the model sees them as it was built
                        data <- as.data.frame(data)
                        for (i in seq_len(nrow(predictor_scaling))) {
                            v <- predictor_scaling$variable[i]
                            if (v %in% colnames(data)) data[[v]] <- as.numeric((data[[v]] - predictor_scaling$mean[i]) / predictor_scaling$sd[i])
                        }
                    }

                    if (model_type == "linear_regression") {
                        output[[pred_col]] <- unname(as.numeric(tryCatch(
                            stats::predict(model, data),
                            error = function(e) {
                                if (grepl("new level", conditionMessage(e))) {
                                    stop("predictWithModel(): the new data contain a group the model was not built on (",
                                         conditionMessage(e), "). A model can only predict for groups it has seen.", call. = FALSE)
                                }
                                stop(e)
                            }
                        )))
                    }

                    if (model_type == "logistic_regression") {
                        ## buildModel() fits the response as a factor whose second level is the
                        ## positive class; the cut-off is 0.5 (move it with mutate() on the probability)
                        response <- model$model[[1]]
                        classes <- if (is.factor(response)) levels(response) else levels(factor(response))
                        if (length(classes) != 2) {
                            stop("predictWithModel() could not find the two response values of this logistic model. Pass the output of buildModel().")
                        }
                        probability <- unname(as.numeric(tryCatch(
                            stats::predict(model, data, type = "response"),
                            error = function(e) {
                                if (grepl("new level", conditionMessage(e))) {
                                    stop("predictWithModel(): the new data contain a group the model was not built on (",
                                         conditionMessage(e), "). A model can only predict for groups it has seen.", call. = FALSE)
                                }
                                stop(e)
                            }
                        )))
                        output[[paste0("probability_", classes[2])]] <- probability
                        output[["predicted_class"]] <- ifelse(probability >= 0.5, classes[2], classes[1])
                    }

                    if (model_type == "random_forest_regression") {
                        output[[pred_col]] <- unname(as.numeric(unlist(stats::predict(model, new_data = data))))
                    }

                    if (model_type == "random_forest_classification") {
                        probs <- as.data.frame(stats::predict(model, new_data = data, type = "prob"))
                        classes <- sub("^\\.pred_", "", colnames(probs))
                        for (k in seq_along(classes)) {
                            output[[paste0("probability_", classes[k])]] <- probs[[k]]
                        }
                        output[["predicted_class"]] <- classes[max.col(as.matrix(probs), ties.method = "first")]
                    }

                    if (model_type == "contrastive_learning") {
                        # Convert data to the required matrix format
                        x <- data[,colnames(data) %in% model$input_variables]
                        embedding_matrix <- as.matrix(x)

                        # Load model architecture and parameters
                        embedding_model <- model$architecture()
                        embedding_model$load_state_dict(model$state_dict)

                        # Pass the data through the model to get embeddings, one column per dimension
                        embeddings <- embedding_model(torch::torch_tensor(embedding_matrix))
                        emb <- as.matrix(embeddings$detach()$cpu())  # Detach and bring to CPU if running on GPU
                        for (k in seq_len(ncol(emb))) {
                            output[[paste0("embedding_", k)]] <- emb[, k]
                        }
                    }

                ## Return
                    return(output)

            }

