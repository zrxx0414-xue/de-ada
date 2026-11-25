unwrap_mcmc <- function(x) {
  if (all(c("paramsRec", "lLRec") %in% names(x))) return(x)
  if ("paramsMCMC" %in% names(x)) return(unwrap_mcmc(x$paramsMCMC))
  unwrap_mcmc(x[[1]])
}

solve_ode_model <- function(model, parameters, initial_conditions, times) {
  ode_model <- function(time, species, parameters) {
    x1 <- species[1]
    x2 <- species[2]

    derivative <- switch(
      model,
      LV = c(
        x1 * (parameters[1] - parameters[2] * x2),
        -x2 * (parameters[3] - parameters[4] * x1)
      ),
      LVEP = c(
        x1 * (parameters[1] - parameters[2] * x2 - parameters[5] * x1),
        -x2 * (parameters[3] - parameters[4] * x1)
      ),
      LVHT2 = c(
        parameters[1] * x1 - parameters[2] * x1 * x2 / (1 + parameters[5] * x1),
        -parameters[3] * x2 + parameters[4] * x1 * x2 / (1 + parameters[5] * x1)
      )
    )

    list(derivative)
  }

  solution <- deSolve::ode(
    y = initial_conditions,
    times = times,
    func = ode_model,
    parms = parameters
  )

  list(
    x1 = as.numeric(solution[, 2]),
    x2 = as.numeric(solution[, 3])
  )
}

get_prediction <- function(mcmc_result, model, times, burn_fraction = 0.25) {
  mcmc <- unwrap_mcmc(mcmc_result)
  log_likelihood <- as.numeric(mcmc$lLRec)
  first_keep <- floor(length(log_likelihood) * burn_fraction) + 1L
  keep <- first_keep:length(log_likelihood)
  ml_iteration <- keep[which.max(log_likelihood[keep])]
  n_parameters <- c(LV = 4L, LVEP = 5L, LVHT2 = 5L)[model]
  ml_values <- as.numeric(mcmc$paramsRec[ml_iteration, ])
  parameters <- ml_values[seq_len(n_parameters)]
  initial_conditions <- ml_values[n_parameters + 1:2]
  prediction <- solve_ode_model(model, parameters, initial_conditions, times)

  c(
    prediction,
    list(
      parameters = parameters,
      initial_conditions = initial_conditions,
      iteration = ml_iteration,
      log_likelihood = log_likelihood[ml_iteration]
    )
  )
}

get_predictions <- function(run, LV, LVEP, LVHT2, times,
                            burn_fraction = 0.25,
                            candidate = c("LV", "LVEP", "LVHT2")) {
  models <- list(LV = LV, LVEP = LVEP, LVHT2 = LVHT2)
  predictions <- lapply(candidate, function(model) {
    get_prediction(models[[model]][[run]], model, times, burn_fraction)
  })
  setNames(predictions, candidate)
}
