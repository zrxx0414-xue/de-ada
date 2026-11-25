residual_transform <- function(predictions, test_data,
                               input = c("prediction1", "prediction2", "residual12")) {
  input <- match.arg(input)
  test <- as.matrix(test_data$data)
  truth <- as.matrix(test_data$data.true)

  if (input == "prediction1") {
    return(c(list(Baseline = test[, 1]), lapply(predictions, `[[`, "x1")))
  }

  if (input == "prediction2") {
    return(c(list(Baseline = test[, 2]), lapply(predictions, `[[`, "x2")))
  }

  baseline <- c(test[, 1] - truth[, 1], test[, 2] - truth[, 2])
  residuals <- lapply(predictions, function(prediction) {
    c(test[, 1] - prediction$x1, test[, 2] - prediction$x2)
  })
  c(list(Baseline = baseline), residuals)
}
