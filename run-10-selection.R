run_10_selection <- function(LV = NULL, LVEP = NULL, LVHT2 = NULL, parameter_count,
                             basis_function, true_model,
                             test_data = NULL, test_parameters = NULL,
                             test_noise = NULL, test_seed = NULL,
                             initial_conditions = c(5, 3),
                             time_range = c(0, 2), n_timepoints = 1000,
                             n_observations = 22, burn_fraction = 0.25,
                             alpha = 0.05, maxit = 2000,
                             print_results = TRUE, kernel_bw = "adaptive", dim_basis = 1,
                             candidate = c(1L, 2L, 3L)) {
  start <- proc.time()[["elapsed"]]
  model_names <- c("LV", "LVEP", "LVHT2")
  if (is.numeric(candidate) && all(candidate %in% seq_along(model_names)))
    candidate <- model_names[candidate]
  if (!is.character(candidate) || length(candidate) < 2L ||
      anyDuplicated(candidate) || !all(candidate %in% model_names))
    stop("candidate must identify two or three distinct models: 1 = LV, 2 = LVEP, 3 = LVHT2.")
  candidate <- model_names[model_names %in% candidate]
  model_results <- list(LV = LV, LVEP = LVEP, LVHT2 = LVHT2)
  if (any(vapply(model_results[candidate], function(x) length(x) < 10L, logical(1))))
    stop("Supply results for all 10 datasets for each candidate model.")
  inputs <- c("prediction1", "prediction2", "residual12")

  if (is.null(test_data)) {
    times <- seq(time_range[1], time_range[2], length.out = n_timepoints)
    test_data <- generate_noisy_test_data(
      model = true_model,
      parameters = test_parameters,
      noise = test_noise,
      seed = test_seed,
      initial_conditions = initial_conditions,
      times = times
    )
  } else {
    times <- test_data$time
  }

  drm_runs <- lapply(1:10, function(run) {
    predictions <- get_predictions(run, LV, LVEP, LVHT2, times, burn_fraction, candidate)
    current_test <- if (!is.null(test_data$data)) test_data else test_data[[run]]

    lapply(inputs, function(input) {
      samples <- residual_transform(predictions, current_test, input)
      drm_start <- proc.time()[["elapsed"]]
      result <- drm_selection(samples, basis_function, parameter_count, alpha, maxit,
                              kernel_bw = kernel_bw, dim_basis = dim_basis)
      result$runtime_seconds <- proc.time()[["elapsed"]] - drm_start
      result$pairwise$dataset <- run
      result$pairwise$input <- input
      result$individual$dataset <- run
      result$individual$input <- input
      result
    })
  })

  selections <- do.call(rbind, lapply(1:10, function(run) {
    selected <- vapply(drm_runs[[run]], `[[`, character(1), "selected_model")
    data.frame(
      dataset = run,
      prediction1 = selected[1],
      prediction2 = selected[2],
      residual12 = selected[3],
      row.names = NULL
    )
  }))

  detailed <- do.call(rbind, lapply(1:10, function(run) {
    do.call(rbind, lapply(seq_along(inputs), function(i) {
      result <- drm_runs[[run]][[i]]
      pair_names <- paste(result$pairwise$model_1, result$pairwise$model_2, sep = "_vs_")
      p_values <- setNames(result$pairwise$p_value, paste0("p_", pair_names))
      wald_values <- setNames(result$individual$wald, paste0("wald_", result$individual$model))

      data.frame(
        dataset = run,
        input = inputs[i],
        basis_dimension = result$basis_dimension,
        result$diagnostics,
        runtime_seconds = result$runtime_seconds,
        as.list(c(p_values, wald_values)),
        selected_model = result$selected_model,
        row.names = NULL,
        check.names = FALSE
      )
    }))
  }))

  wald_columns <- paste0("wald_", candidate)
  detailed$minimum_wald_model <- apply(detailed[, wald_columns], 1, function(x) {
    if (all(is.na(x))) return(NA_character_)
    sub("wald_", "", names(x)[which.min(x)])
  })
  detailed$true_model_minimum_wald <- detailed$minimum_wald_model == true_model

  true_model_min_wald <- do.call(rbind, lapply(split(detailed, detailed$input), function(x) {
    data.frame(
      input = x$input[1],
      count = sum(x$true_model_minimum_wald, na.rm = TRUE),
      total = nrow(x),
      proportion = sum(x$true_model_minimum_wald, na.rm = TRUE) / nrow(x),
      row.names = NULL
    )
  }))
  rownames(true_model_min_wald) <- NULL

  benchmark <- benchmark_aic_bic(
    LV, LVEP, LVHT2, parameter_count,
    n_observations, burn_fraction, true_model, candidate = candidate
  )

  all_selections <- merge(benchmark$selections, selections, by = "dataset")
  methods <- c("AIC", "BIC", "prediction1", "prediction2", "residual12")
  accuracy <- data.frame(
    method = methods,
    correct = vapply(methods, function(x) sum(all_selections[[x]] == true_model, na.rm = TRUE), numeric(1)),
    accuracy = vapply(methods, function(x) sum(all_selections[[x]] == true_model, na.rm = TRUE) / nrow(all_selections), numeric(1))
  )

  if (print_results) {
    print(all_selections)
    print(accuracy)
  }

  list(
    test_data = test_data,
    candidate = candidate,
    runtime = list(
      total_seconds = proc.time()[["elapsed"]] - start,
      drm_total_seconds = sum(detailed$runtime_seconds),
      drm_mean_seconds = mean(detailed$runtime_seconds)
    ),
    selections = all_selections,
    accuracy = accuracy,
    detailed = list(
      results = detailed,
      true_model_min_wald = true_model_min_wald
    ),
    benchmark_scores = benchmark$scores,
    drm_fits = drm_runs
  )
}

default_basis_functions <- function() {
  signed_log <- function(y) sign(y) * log1p(abs(y))

  list(
    "q(y) = y" = function(y) c(y = y),
    "q(y) = y^2" = function(y) c(y2 = y^2),
    "q(y) = (y, y^2)" = function(y) c(y = y, y2 = y^2),
    "q(y) = y^3" = function(y) c(y3 = y^3),
    "q(y) = (y, y^2, y^3)" = function(y) c(y = y, y2 = y^2, y3 = y^3),
    "q(y) = (y, y^3)" = function(y) c(y = y, y3 = y^3),
    "q(y) = signed_log(y)" = function(y) c(logy = signed_log(y)),
    "q(y) = (y, signed_log(y))" = function(y) c(y = y, logy = signed_log(y)),
    "q(y) = (y, y^2, signed_log(y))" = function(y) c(y = y, y2 = y^2, logy = signed_log(y))
  )
}

run_multiple_basis <- function(LV = NULL, LVEP = NULL, LVHT2 = NULL, parameter_count,
                               true_model, test_parameters, test_noise,
                               test_seed, basis_functions = default_basis_functions(),
                               initial_conditions = c(5, 3),
                               time_range = c(0, 2), n_timepoints = 1000,
                               n_observations = 22, burn_fraction = 0.25,
                               alpha = 0.05, maxit = 2000,
                               kernel_bw = "adaptive", dim_basis = 1, candidate = c(1L, 2L, 3L)) {
  test_data <- generate_noisy_test_data(
    model = true_model,
    parameters = test_parameters,
    noise = test_noise,
    seed = test_seed,
    initial_conditions = initial_conditions,
    times = seq(time_range[1], time_range[2], length.out = n_timepoints)
  )

  results <- lapply(basis_functions, function(basis_function) {
    run_10_selection(
      LV = LV,
      LVEP = LVEP,
      LVHT2 = LVHT2,
      parameter_count = parameter_count,
      basis_function = basis_function,
      true_model = true_model,
      test_data = test_data,
      n_observations = n_observations,
      burn_fraction = burn_fraction,
      alpha = alpha,
      maxit = maxit,
      kernel_bw = kernel_bw,
      dim_basis = dim_basis,
      candidate = candidate,
      print_results = FALSE
    )
  })

  accuracy_by_basis <- do.call(rbind, lapply(names(results), function(basis_name) {
    accuracy <- results[[basis_name]]$accuracy
    accuracy$basis <- basis_name
    accuracy[, c("basis", "method", "correct", "accuracy")]
  }))

  drm_inputs <- c("prediction1", "prediction2", "residual12")
  scores <- accuracy_by_basis[accuracy_by_basis$method %in% drm_inputs, ]
  rownames(scores) <- NULL
  best_basis <- do.call(rbind, lapply(drm_inputs, function(input) {
    current <- scores[scores$method == input, ]
    current[current$accuracy == max(current$accuracy), ]
  }))
  rownames(best_basis) <- NULL
  best_results <- setNames(lapply(drm_inputs, function(input) {
    basis_names <- best_basis$basis[best_basis$method == input]
    results[basis_names]
  }), drm_inputs)

  benchmark_accuracy <- results[[1]]$accuracy
  benchmark_accuracy <- benchmark_accuracy[benchmark_accuracy$method %in% c("AIC", "BIC"), ]
  summary <- rbind(
    data.frame(
      method = benchmark_accuracy$method,
      accuracy = benchmark_accuracy$accuracy,
      basis = NA_character_,
      row.names = NULL
    ),
    do.call(rbind, lapply(drm_inputs, function(input) {
      current <- best_basis[best_basis$method == input, ]
      data.frame(
        method = input,
        accuracy = current$accuracy[1],
        basis = paste(current$basis, collapse = " | "),
        row.names = NULL
      )
    }))
  )

  print(summary)

  list(
    test_data = test_data,
    summary = summary,
    scores = scores,
    accuracy_by_basis = accuracy_by_basis,
    best_basis = best_basis,
    best_results = best_results,
    all_results = results
  )
}
