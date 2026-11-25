library(drmdel)

drm_selection <- function(samples, basis_function, parameter_count,
                          alpha = 0.05, maxit = 2000,
                          kernel_bw = "adaptive", dim_basis = 1) {
  models <- if (is.null(names(samples))) c("LV", "LVEP", "LVHT2") else names(samples)[-1L]
  stage <- "basis_construction"
  basis_ready <- FALSE
  convergence <- d <- NA_integer_
  optimizer_message <- error_call <- error_message <- NA_character_
  warnings <- character()
  basis_label <- if (is.function(basis_function)) "custom" else as.character(basis_function)

  result <- tryCatch(withCallingHandlers({
    if (identical(basis_function, "adaptive"))
      basis_function <- fit_adaptive_basis(samples, kernel_bw = kernel_bw, dim_basis = dim_basis)$fun
    basis_ready <- TRUE
    stage <- "drm_setup"
    fitter <- drmdel::drmdel
    environment(fitter) <- list2env(list(optim = function(...) {
      stage <<- "drm_optimization"
      opt <- stats::optim(...)
      convergence <<- opt$convergence
      if (!is.null(opt$message)) optimizer_message <<- opt$message
      if (opt$convergence != 0L)
        warning("DRM optimizer did not converge; optim code = ", opt$convergence, call. = FALSE)
      stage <<- "drm_postprocessing"
      opt
    }), parent = environment(fitter))
    all_y <- unlist(samples, use.names = FALSE)
    fit <- fitter(x = all_y, n_samples = lengths(samples), basis_func = basis_function,
                  control = list(maxit = maxit))
    d <- fit$drm_info$d
    stage <- "covariance"
    block_size <- d + 1L
    covariance <- as.matrix(meleCov(fit))
    n_total <- length(all_y)
    block_index <- function(i) ((i - 1L) * block_size + 1L):(i * block_size)
    beta_index <- function(i) block_index(i)[-1]

    stage <- "wald_individual"
    individual <- do.call(rbind, lapply(seq_along(models), function(i) {
      index <- beta_index(i)
      beta <- fit$mele[index]
      wald <- n_total * drop(t(beta) %*% solve(covariance[index, index, drop = FALSE], beta))
      data.frame(model = models[i], wald = wald,
                 p_value = pchisq(wald, d, lower.tail = FALSE))
    }))

    pairs <- combn(seq_along(models), 2, simplify = FALSE)
    stage <- "wald_pairwise"
    pairwise <- do.call(rbind, lapply(pairs, function(pair) {
      index_1 <- beta_index(pair[1])
      index_2 <- beta_index(pair[2])
      difference <- fit$mele[index_1] - fit$mele[index_2]
      contrast <- matrix(0, d, length(fit$mele))
      contrast[cbind(seq_len(d), index_1)] <- 1
      contrast[cbind(seq_len(d), index_2)] <- -1
      difference_covariance <- contrast %*% covariance %*% t(contrast)
      wald <- n_total * drop(t(difference) %*% solve(difference_covariance, difference))
      p_value <- pchisq(wald, d, lower.tail = FALSE)
      model_1 <- models[pair[1]]
      model_2 <- models[pair[2]]
      wald_1 <- individual$wald[individual$model == model_1]
      wald_2 <- individual$wald[individual$model == model_2]
      loser <- if (p_value < alpha) {
        if (wald_1 > wald_2) model_1 else model_2
      } else {
        NA_character_
      }

      data.frame(
        model_1 = model_1,
        model_2 = model_2,
        pairwise_wald = wald,
        df = d,
        p_value = p_value,
        individual_wald_1 = wald_1,
        individual_wald_2 = wald_2,
        eliminated = loser
      )
    }))

    stage <- "model_selection"
    eliminated <- unique(na.omit(pairwise$eliminated))
    remaining <- setdiff(models, eliminated)
    minimum_parameters <- min(parameter_count[remaining])
    selected <- remaining[parameter_count[remaining] == minimum_parameters]
    selected_model <- if (length(selected) == 1L) selected else paste(selected, collapse = " + ")

    list(
      fit = fit,
      basis_dimension = d,
      pairwise = pairwise,
      individual = individual,
      eliminated = eliminated,
      remaining = remaining,
      selected_model = selected_model
    )
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }, error = function(e) {
    calls <- vapply(sys.calls(), function(x) paste(deparse(x[[1]]), collapse = ""), character(1))
    if ("adaptive_bw" %in% calls) stage <<- "adaptive_bandwidth"
    else if ("adaptive_dim" %in% calls) stage <<- "adaptive_dimension"
    else if ("fit_adaptive_basis" %in% calls && "eigen" %in% calls) stage <<- "basis_eigensystem"
    error_call <<- if (is.null(conditionCall(e))) stage
    else substr(paste(deparse(conditionCall(e)), collapse = " "), 1L, 400L)
  }), error = function(e) {
    error_message <<- conditionMessage(e)
    pairs <- combn(models, 2)
    list(
      fit = NULL, basis_dimension = d,
      pairwise = data.frame(model_1 = pairs[1, ], model_2 = pairs[2, ],
        pairwise_wald = NA_real_, df = d, p_value = NA_real_,
        individual_wald_1 = NA_real_, individual_wald_2 = NA_real_, eliminated = NA_character_),
      individual = data.frame(model = models, wald = NA_real_, p_value = NA_real_),
      eliminated = character(), remaining = character(), selected_model = NA_character_
    )
  })

  valid <- is.na(error_message)
  error_type <- NA_character_
  if (!valid) {
    error_type <- if (!is.na(convergence) && convergence != 0L) "nonconvergence" else "other"
    if (grepl("infinite|missing values|non-finite|NA/NaN/Inf", error_message)) error_type <- "nonfinite_values"
    if (grepl("singular|reciprocal condition", error_message)) error_type <- "singular_matrix"
    if (grepl("dimension of the basis", error_message)) error_type <- "invalid_dimension"
  }
  result$diagnostics <- data.frame(
    basis = basis_label, valid = valid, failure_stage = if (valid) NA_character_ else stage,
    basis_ready = basis_ready, drm_convergence = convergence, optimizer_message = optimizer_message,
    error_message = error_message, error_type = error_type, error_call = error_call,
    warning_message = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_
  )
  result
}
