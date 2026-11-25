generate_noisy_test_data <- function(model, parameters, noise, seed,
                                     initial_conditions = c(5, 3),
                                     times = seq(0, 2, length.out = 11),
                                     save_file = NULL,
                                     object_name = "test_data") {
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

  set.seed(seed)
  solution <- deSolve::ode(
    y = initial_conditions,
    times = times,
    func = ode_model,
    parms = parameters
  )
  truth <- as.matrix(solution[, -1, drop = FALSE])
  observation <- truth + matrix(
    rnorm(length(truth), mean = 0, sd = noise),
    nrow = nrow(truth),
    ncol = ncol(truth)
  )
  observation[observation < 0] <- 0.001

  test_data <- list(
    data = observation,
    data.true = truth,
    time = times,
    params = parameters,
    noise = noise,
    seed = seed,
    model = model,
    initial.conditions = initial_conditions
  )

  if (!is.null(save_file)) {
    save_environment <- new.env(parent = emptyenv())
    assign(object_name, test_data, envir = save_environment)
    save(list = object_name, file = save_file, envir = save_environment)
  }

  test_data
}
