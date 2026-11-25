# Change to 2 to run the two-dimensional adaptive basis experiment.
dim_basis <- 1
kernel_bw <- "adaptive"
candidate <- c(1L, 2L, 3L)  # 1 = LV, 2 = LVEP, 3 = LVHT2
source("adaptive_basis.R")

source("generate-noisy-test-data.R")
source("get-predictions.R")
source( "residual-transform.R")
source("drm-selection.R")
source("benchmark-aic-bic.R")
source("run-10-selection.R")

# LV true with noise sd=0.1581139
resultLV1 <- run_10_selection(
  LV = get("LV-LVtrue-0.1581139"),
  LVEP = get("LVEP-LVtrue-0.1581139"),
  LVHT2 = get("LVHT2-LVtrue-0.1581139"),
  parameter_count = c(LV = 4, LVEP = 5, LVHT2 = 5),
  basis_function = "adaptive",
  kernel_bw = kernel_bw,
  dim_basis = dim_basis,
  candidate = candidate,
  true_model = "LV",
  test_parameters = c(2,1,4,1),
  test_noise = 0.1581139,
  test_seed = 12345,
  initial_conditions = c(5, 3),
  time_range = c(0, 2),
  n_timepoints = 200
)

# LV true with noise sd=0.28374
resultLV2 <- run_10_selection(
  LV = get("LV-LVtrue-0.28374"),
  LVEP = get("LVEP-LVtrue-0.28374"),
  LVHT2 = get("LVHT2-LVtrue-0.28374"),
  parameter_count = c(LV = 4, LVEP = 5, LVHT2 = 5),
  basis_function = "adaptive",
  kernel_bw = kernel_bw,
  dim_basis = dim_basis,
  candidate = candidate,
  true_model = "LV",
  test_parameters = c(2,1,4,1),
  test_noise = 0.28374,
  test_seed = 12345,
  initial_conditions = c(5, 3),
  time_range = c(0, 2),
  n_timepoints = 1000
)

#LVHT2 true with noise sd=1.121186
resultLVHT1 <- run_10_selection(
  LV = get("LV-LVHT2true-1.121186"),
  LVEP = get("LVEP-LVHT2true-1.121186"),
  LVHT2 = get("LVHT2-LVHT2true-1.121186"),
  parameter_count = c(LV = 4, LVEP = 5, LVHT2 = 5),
  basis_function = "adaptive",
  kernel_bw = kernel_bw,
  dim_basis = dim_basis,
  candidate = candidate,
  true_model = "LVHT2",
  test_parameters = c(2.8, 3.5, 1, 2.5, 1),
  test_noise = 1.121186,
  test_seed = 12345,
  initial_conditions = c(5, 3),
  time_range = c(0, 2),
  n_timepoints = 1000
)

#LVHT2 true with noise sd=1.581139
resultLVHT2 <- run_10_selection(
  LV = get("LV-LVHT2true-1.581139"),
  LVEP = get("LVEP-LVHT2true-1.581139"),
  LVHT2 = get("LVHT2-LVHT2true-1.581139"),
  parameter_count = c(LV = 4, LVEP = 5, LVHT2 = 5),
  basis_function = "adaptive",
  kernel_bw = kernel_bw,
  dim_basis = dim_basis,
  candidate = candidate,
  true_model = "LVHT2",
  test_parameters = c(2.8, 3.5, 1, 2.5, 1),
  test_noise = 1.581139,
  test_seed = 12345,
  initial_conditions = c(5, 3),
  time_range = c(0, 2),
  n_timepoints = 500
)

