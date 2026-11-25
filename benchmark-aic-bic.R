benchmark_aic_bic <- function(LV, LVEP, LVHT2, parameter_count,
                              n_observations = 22, burn_fraction = 0.25,
                              true_model = "LV", candidate = c("LV", "LVEP", "LVHT2")) {
  models <- list(LV = LV, LVEP = LVEP, LVHT2 = LVHT2)

  models <- models[candidate]

  scores <- do.call(rbind, lapply(1:10, function(run) {
    do.call(rbind, lapply(names(models), function(model) {
      mcmc <- unwrap_mcmc(models[[model]][[run]])
      log_likelihood <- as.numeric(mcmc$lLRec)
      first_keep <- floor(length(log_likelihood) * burn_fraction) + 1L
      maximum_log_likelihood <- max(log_likelihood[first_keep:length(log_likelihood)], na.rm = TRUE)
      k <- parameter_count[model]

      data.frame(
        dataset = run,
        model = model,
        log_likelihood = maximum_log_likelihood,
        AIC = -2 * maximum_log_likelihood + 2 * k,
        BIC = -2 * maximum_log_likelihood + log(n_observations) * k,
        row.names = NULL
      )
    }))
  }))

  selections <- do.call(rbind, lapply(split(scores, scores$dataset), function(x) {
    data.frame(
      dataset = x$dataset[1],
      AIC = x$model[which.min(x$AIC)],
      BIC = x$model[which.min(x$BIC)]
    )
  }))

  accuracy <- data.frame(
    method = c("AIC", "BIC"),
    correct = c(sum(selections$AIC == true_model), sum(selections$BIC == true_model)),
    accuracy = c(mean(selections$AIC == true_model), mean(selections$BIC == true_model))
  )

  list(scores = scores, selections = selections, accuracy = accuracy)
}
