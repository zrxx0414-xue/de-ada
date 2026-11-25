# This R code contains the main function adaptivebasis, functions for quantile and density estimation, and some other functions that are called in the main function

####################################################################################################################################
# adaptivebasis: produces the adaptive basis function

# Inputs: 
# x: a scalar at which the user wants to evaluate the adaptive basis function

# data: a list that contains m+1 samples, where list[[i]] should contain a vector representing the (i+1)-th sample.

# kernel_bw: a vector that gives user-specified kernel bandwidths used in the step of estimation of eigenfunctions,
             # whose length should be the same as the number of populations considered;
             # if not supplied, the default is the proposed adaptive way of choosing such kernel bandwidths (see the paper for details).

# dim_basis: a scalar that gives user-specified dimension of the adaptive basis function;
             # if not supplied, the default dimension is the proposed adaptive way of choosing such dimension (see the paper for details).

# Output: 
# a vector containing the values of adaptive basis function at x

adaptivebasis <- function(x, data, kernel_bw = "adaptive", dim_basis = "adaptive"){
  m <- length(data) - 1 # m+1 is the total number of populations 
  sizes <- numeric(m+1) # sample sizes 
  combined_data <- NULL # the pooled sample
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)
  
  
  # aggregate the same population data at their unique values to avoid excessive computational cost 
  unique_data <- list() 
  for (i in 1:(m+1)){ 
    unique_pt <- unique(data[[i]]) 
    num_unique <- length(unique_pt)
    # the first column is the value of the unique points, the second column is the number of obs that equal to that value  
    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2) 
    for (j in 1:num_unique){ 
      unique_data[[i]][j,1] <- unique_pt[j] 
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j]) 
    } 
  } 
  
  if (kernel_bw == "adaptive") kernel_bw <- adaptive_bw(data = data, ref_model = "normal")
  if (dim_basis == "adaptive") dim_basis <- adaptive_dim(data = data, kernel_bw = kernel_bw)
  
  ##############################################################################################
  # Construct the M matrix used for obtaining eigenvalues and vectors 
  M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1)) 
  # Record the log density ratios
  ff_tmp_seq <- list() 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
  for (i in 1:length(combined_data)){ 
    x_M <- combined_data[i] 
    log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
    for (j in 1:(m+1)){ 
      ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M 
    }  
  } 
  
  # Obtain the centered log density ratio 
  f_sum_M <- numeric(length(combined_data))  
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){ 
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j] 
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]  
  }  
  f_avg_M <- f_sum_M/(m+1) 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M 
  
  # Obtain the matrix M_mat and its eigen-system
  for (i_1 in 1:(m+1)){ 
    for (i_2 in i_1:(m+1)){ 
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]  
      M_mat[i_1, i_2] <- mean(ff_tmp_mat) 
      M_mat[i_2, i_1] <- M_mat[i_1, i_2] 
    }
  }    
  eigen_sys <- eigen(M_mat) # The (m+1)-th eigenvalue of M_mat will always be 0 
  
  ##############################################################################################
  # Specify the adaptive basis function
  log_base <- log(sum(unique_data[[1]][,2]*dnorm(x = (x-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
  f_sum <- 0 
  f_naive_LDR <- numeric(m+1) 
  f_naive_LDR[1] <- 0 
  for (i in 2:(m+1)){ 
    f_naive_LDR[i] <- log(sum(unique_data[[i]][,2]*dnorm(x = (x-unique_data[[i]][,1]), mean = 0, sd = kernel_bw[i]))/sum(unique_data[[i]][,2])) - log_base 
    # centerize the log density ratios 
    f_naive_LDR[i] <- f_naive_LDR[i] - mean_log_DR[i]
    f_sum <- f_sum + f_naive_LDR[i]  
  }  
  f_avg <- f_sum/(m+1) 
  f_naive_LDR <- f_naive_LDR - f_avg 
  
  fpca_basis_fn <- numeric(dim_basis)  
  for (j in 1:dim_basis){ 
    fpca_basis_fn[j] <- (eigen_sys$values[j]^(-1/2))*sum((eigen_sys$vectors[(1:(m+1)),j])*f_naive_LDR) 
  } 
  
  return(fpca_basis_fn) 
}

# Prepare once per replicate; reuse q_fun for all evaluations.
fit_adaptivebasis <- function(data, kernel_bw = "adaptive", dim_basis = "adaptive"){
  m <- length(data) - 1 # m+1 is the total number of populations 
  sizes <- numeric(m+1) # sample sizes 
  combined_data <- NULL # the pooled sample
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)
  
  
  # aggregate the same population data at their unique values to avoid excessive computational cost 
  unique_data <- list() 
  for (i in 1:(m+1)){ 
    unique_pt <- unique(data[[i]]) 
    num_unique <- length(unique_pt)
    # the first column is the value of the unique points, the second column is the number of obs that equal to that value  
    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2) 
    for (j in 1:num_unique){ 
      unique_data[[i]][j,1] <- unique_pt[j] 
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j]) 
    } 
  } 
  
  if (identical(kernel_bw, "adaptive")) kernel_bw <- adaptive_bw(data = data, ref_model = "normal")
  if (dim_basis == "adaptive") dim_basis <- adaptive_dim(data = data, kernel_bw = kernel_bw)
  
  ##############################################################################################
  # Construct the M matrix used for obtaining eigenvalues and vectors 
  M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1)) 
  # Record the log density ratios
  ff_tmp_seq <- list() 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
  for (i in 1:length(combined_data)){ 
    x_M <- combined_data[i] 
    log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
    for (j in 1:(m+1)){ 
      ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M 
    }  
  } 
  
  # Obtain the centered log density ratio 
  f_sum_M <- numeric(length(combined_data))  
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){ 
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j] 
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]  
  }  
  f_avg_M <- f_sum_M/(m+1) 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M 
  
  # Obtain the matrix M_mat and its eigen-system
  for (i_1 in 1:(m+1)){ 
    for (i_2 in i_1:(m+1)){ 
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]  
      M_mat[i_1, i_2] <- mean(ff_tmp_mat) 
      M_mat[i_2, i_1] <- M_mat[i_1, i_2] 
    }
  }    
  eigen_sys <- eigen(M_mat) # The (m+1)-th eigenvalue of M_mat will always be 0 
  
  ##############################################################################################
  # Specify the adaptive basis function
  q_fun <- function(x) {
  log_base <- log(sum(unique_data[[1]][,2]*dnorm(x = (x-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
  f_sum <- 0 
  f_naive_LDR <- numeric(m+1) 
  f_naive_LDR[1] <- 0 
  for (i in 2:(m+1)){ 
    f_naive_LDR[i] <- log(sum(unique_data[[i]][,2]*dnorm(x = (x-unique_data[[i]][,1]), mean = 0, sd = kernel_bw[i]))/sum(unique_data[[i]][,2])) - log_base 
    # centerize the log density ratios 
    f_naive_LDR[i] <- f_naive_LDR[i] - mean_log_DR[i]
    f_sum <- f_sum + f_naive_LDR[i]  
  }  
  f_avg <- f_sum/(m+1) 
  f_naive_LDR <- f_naive_LDR - f_avg 
  
  fpca_basis_fn <- numeric(dim_basis)  
  for (j in 1:dim_basis){ 
    fpca_basis_fn[j] <- (eigen_sys$values[j]^(-1/2))*sum((eigen_sys$vectors[(1:(m+1)),j])*f_naive_LDR) 
  } 
  
  return(fpca_basis_fn)
  }
  list(d = dim_basis, fun = q_fun)
}

####################################################################################################################################



####################################################################################################################################
# quantile_est: produces quantile estimates under the density ratio model with the adaptive basis function

# Inputs: 
# data: a list that contains m+1 samples, where list[[i]] should contain a vector representing the (i+1)-th sample.

# kernel_bw: a vector that gives user-specified kernel bandwidths used in the step of estimation of eigenfunctions,
             # whose length should be the same as the number of populations considered;
             # if not supplied, the default is the proposed adaptive way of choosing such kernel bandwidths (see the paper for details).

# dim_basis: a scalar that gives user-specified dimension of the adaptive basis function;
             # if not supplied, the default dimension is the proposed adaptive way of choosing such dimension (see the paper for details).

# index: a vector of indices of the populations whose quantiles are of interest, each of which is a positive integer;
         # population G_0 corresponds to an index being 1, G_1 corresponds to an index being 2, and so on...

# prob: a vector of levels of quantiles of interest, each of which is a numeric value between 0 and 1 (exclusively);
        # it should be of the same length as index, 
        # and the i-th element of prob should correspond to the i-th element of index.

# Output: 
# a vector whose i-th element is the quantile estimate for the index[i]-th population with level prob[i]

quantile_est <- function(data, kernel_bw = "adaptive", dim_basis = "adaptive", index, prob){
  m <- length(data) - 1 # m+1 is the total number of populations 
  sizes <- numeric(m+1) # sample sizes 
  combined_data <- NULL # the pooled sample
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)
  
  
  # aggregate the same population data at their unique values to avoid excessive computational cost 
  unique_data <- list() 
  for (i in 1:(m+1)){ 
    unique_pt <- unique(data[[i]]) 
    num_unique <- length(unique_pt)
    # the first column is the value of the unique points, the second column is the number of obs that equal to that value  
    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2) 
    for (j in 1:num_unique){ 
      unique_data[[i]][j,1] <- unique_pt[j] 
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j]) 
    } 
  } 
  
  #######################################################################################
  # first obtain the adaptive basis function
  # we do not call the adaptivebasis() function directly here to save computational cost, as some values computed in this process can be reused later
  
  if (kernel_bw == "adaptive") kernel_bw <- adaptive_bw(data = data, ref_model = "normal")
  if (dim_basis == "adaptive") dim_basis <- adaptive_dim(data = data, kernel_bw = kernel_bw)
  
  # Construct the M matrix used for obtaining eigenvalues and vectors, using KDE
  M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1)) 
  # Record the log density ratios
  ff_tmp_seq <- list() 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
  for (i in 1:length(combined_data)){ 
    x_M <- combined_data[i] 
    log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
    for (j in 1:(m+1)){ 
      ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M 
    }  
  } 
  
  # Obtain the centered log density ratio 
  f_sum_M <- numeric(length(combined_data))  
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){ 
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j] 
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]  
  }  
  f_avg_M <- f_sum_M/(m+1) 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M 
  
  # Obtain the matrix M_mat and its eigen-system
  for (i_1 in 1:(m+1)){ 
    for (i_2 in i_1:(m+1)){ 
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]  
      M_mat[i_1, i_2] <- mean(ff_tmp_mat) 
      M_mat[i_2, i_1] <- M_mat[i_1, i_2] 
    }
  }    
  eigen_sys <- eigen(M_mat) # The (m+1)-th eigenvalue of M_mat will always be 0 
  
  # Evaluate and store the values of the basis function at every data points to avoid excessive amount of basis function calculations when fitting the DRM
  eigenfun_unique_data <- matrix(NA, nrow = length(unique(combined_data)), ncol = dim_basis) 
  for (ii in 1:length(unique(combined_data))){
    x <- unique(combined_data)[ii]
    x_loc <- which(combined_data == x)[1] #only use the first one to avoid duplicate issue
    log_DR <- numeric(m+1) 
    for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc] #already calculated and recorded when obtaining M_mat
    # obtain the estimated eigenfunction evaluated at the unique data with location x_loc
    fpca_basis_fn <- numeric(dim_basis)
    for (jj in 1:dim_basis) fpca_basis_fn[jj] <- (eigen_sys$values[jj]^(-1))*sum((eigen_sys$vectors[(1:(m+1)),jj])*log_DR) 
    eigenfun_unique_data[ii,] <- fpca_basis_fn
  }
  
  #######################################################################################
  # secondly, fit the DRM to the data using the adaptive basis function of different dimensions
  # create the discrete function of the basis function
  eigenfun_discrete <- function(x){ 
    loc_unique_data <- which(x == unique(combined_data))  
    return(eigenfun_unique_data[loc_unique_data,]) 
  } 
  
  # fit the DRM to the data using various basis functions
  # use a function in the package "drmdel"
  drmfit_FPCA <- drmdel(x=combined_data, n_samples=sizes, basis_func=eigenfun_discrete) 
  
  #######################################################################################
  # finally, obtain quantile estimates 
  quan_drmFPCA <- quantileDRM(k = (index-1), p = prob, drmfit = drmfit_FPCA, cov = FALSE)$est 
  
  return(quan_drmFPCA)
}
####################################################################################################################################



####################################################################################################################################
# density_est: produces density estimates under the density ratio model with the adaptive basis function

# Inputs: 
# x_eval: a scalar at which the user wants to evaluate the density estimates

# data: a list that contains m+1 samples, where list[[i]] should contain a vector representing the (i+1)-th sample.

# kernel_bw: a vector that gives user-specified kernel bandwidths used in the step of estimation of eigenfunctions,
             # whose length should be the same as the number of populations considered;
             # if not supplied, the default is the proposed adaptive way of choosing such kernel bandwidths (see the paper for details).

# dim_basis: a scalar that gives user-specified dimension of the adaptive basis function;
             # if not supplied, the default dimension is the proposed adaptive way of choosing such dimension (see the paper for details).

# index: a vector of indices of the populations whose densities are of interest, each of which is a positive integer;
         # population G_0 corresponds to an index being 1, G_1 corresponds to an index being 2, and so on...

# Outputs: 
# den_at_x: a vector whose i-th element is the density estimate for the index[i]-th population evaluated at x

# density_grid_x: a vector of 512 equally spaced points in the range of pooled sample

# density_grid_y: a matrix whose i-th row is the vector of the density estimates for the index[i]-th population evaluated at density_grid_x

density_est <- function(x_eval, data, kernel_bw = "adaptive", dim_basis = "adaptive", index){
  m <- length(data) - 1 # m+1 is the total number of populations 
  sizes <- numeric(m+1) # sample sizes 
  combined_data <- NULL # the pooled sample
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)
  
  
  # aggregate the same population data at their unique values to avoid excessive computational cost 
  unique_data <- list() 
  for (i in 1:(m+1)){ 
    unique_pt <- unique(data[[i]]) 
    num_unique <- length(unique_pt)
    # the first column is the value of the unique points, the second column is the number of obs that equal to that value  
    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2) 
    for (j in 1:num_unique){ 
      unique_data[[i]][j,1] <- unique_pt[j] 
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j]) 
    } 
  } 
  
  #######################################################################################
  # first obtain the adaptive basis function
  # we do not call the adaptivebasis() function directly here to save computational cost, as some values computed in this process can be reused later
  
  if (kernel_bw == "adaptive") kernel_bw <- adaptive_bw(data = data, ref_model = "normal")
  if (dim_basis == "adaptive") dim_basis <- adaptive_dim(data = data, kernel_bw = kernel_bw)
  
  # Construct the M matrix used for obtaining eigenvalues and vectors, using KDE
  M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1)) 
  # Record the log density ratios
  ff_tmp_seq <- list() 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
  for (i in 1:length(combined_data)){ 
    x_M <- combined_data[i] 
    log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
    for (j in 1:(m+1)){ 
      ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M 
    }  
  } 
  
  # Obtain the centered log density ratio 
  f_sum_M <- numeric(length(combined_data))  
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){ 
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j] 
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]  
  }  
  f_avg_M <- f_sum_M/(m+1) 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M 
  
  # Obtain the matrix M_mat and its eigen-system
  for (i_1 in 1:(m+1)){ 
    for (i_2 in i_1:(m+1)){ 
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]  
      M_mat[i_1, i_2] <- mean(ff_tmp_mat) 
      M_mat[i_2, i_1] <- M_mat[i_1, i_2] 
    }
  }    
  eigen_sys <- eigen(M_mat) # The (m+1)-th eigenvalue of M_mat will always be 0 
  
  # Evaluate and store the values of the basis function at every data points to avoid excessive amount of basis function calculations when fitting the DRM
  eigenfun_unique_data <- matrix(NA, nrow = length(unique(combined_data)), ncol = dim_basis) 
  for (ii in 1:length(unique(combined_data))){
    x <- unique(combined_data)[ii]
    x_loc <- which(combined_data == x)[1] #only use the first one to avoid duplicate issue
    log_DR <- numeric(m+1) 
    for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc] #already calculated and recorded when obtaining M_mat
    # obtain the estimated eigenfunction evaluated at the unique data with location x_loc
    fpca_basis_fn <- numeric(dim_basis)
    for (jj in 1:dim_basis) fpca_basis_fn[jj] <- (eigen_sys$values[jj]^(-1))*sum((eigen_sys$vectors[(1:(m+1)),jj])*log_DR) 
    eigenfun_unique_data[ii,] <- fpca_basis_fn
  }
  
  #######################################################################################
  # secondly, fit the DRM to the data using the adaptive basis function of different dimensions
  # create the discrete function of the basis function
  eigenfun_discrete <- function(x){ 
    loc_unique_data <- which(x == unique(combined_data))  
    return(eigenfun_unique_data[loc_unique_data,]) 
  } 
  
  # fit the DRM to the data using various basis functions
  # use a function in the package "drmdel"
  drmfit_FPCA <- drmdel(x=combined_data, n_samples=sizes, basis_func=eigenfun_discrete) 
  
  #######################################################################################
  # finally, obtain density estimates 
  den_at_x <- numeric(length(index))
  density_grid_x <- seq(from = min(combined_data), to = max(combined_data), length.out = 512)
  density_grid_y <- matrix(nrow = length(index), ncol = 512)
  for (kk in index){
    # obtain the quantile estimate
    den_drmFPCA <- densityDRM(k = (kk-1), drmfit = drmfit_FPCA, cut = 0) 
    # use the suggested rule of thumb kernel bandwidth
    den_drmFPCA <- densityDRM(k = (kk-1), drmfit = drmfit_FPCA, cut = 0, bw = den_drmFPCA$bw*0.9/1.06) 
    p_est_FPCA <- drmfit_FPCA$p_est[drmfit_FPCA$p_est[,1] == (kk-1), 2:3]
    
    density_grid_y[which(index == kk),] <- den_drmFPCA$y
    den_at_x[which(index == kk)] <- sum(p_est_FPCA[,2]*dnorm(x = (x_eval-p_est_FPCA[,1]), mean = 0, sd = den_drmFPCA$bw))
  }
  
  return(list(den_at_x = den_at_x, density_grid_x = density_grid_x, density_grid_y = density_grid_y))
}
####################################################################################################################################



####################################################################################################################################
# functions that are called in the main function adaptivebasis

####################################################################################################################################
# adaptive_bw: find the adaptive kernel bandwidths used in the step of estimation of eigenfunctions  

# Inputs: 
# data: a list that contains m+1 samples, where list[[i]] should contain a vector representing the (i+1)-th sample.

# ref_model: a reference model for the data generating distributions, can be either "normal" or "gamma"

# Output: 
# a vector of the adaptive kernel bandwidths that has the same length as the number of populations

adaptive_bw <- function(data, ref_model){
  m <- length(data) - 1 # m+1 is the total number of populations 
  sizes <- numeric(m+1) # sample sizes 
  combined_data <- NULL # the pooled sample
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)
  
  # rule of thumb bandwidth suggested by Silverman
  bw_RoT <- numeric(m+1)
  for (i in 1:(m+1)) bw_RoT[i] <- 0.9*(sizes[i]^(-1/5))*min(sd(data[[i]]), IQR(data[[i]])/1.34)
  
  # aggregate the same population data at their unique values to avoid excessive computational cost 
  unique_data <- list() 
  for (i in 1:(m+1)){ 
    unique_pt <- unique(data[[i]]) 
    num_unique <- length(unique_pt)
    # the first column is the value of the unique points, the second column is the number of obs that equal to that value  
    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2) 
    for (j in 1:num_unique){ 
      unique_data[[i]][j,1] <- unique_pt[j] 
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j]) 
    } 
  } 
  
  
  #######################################################################################
  # Construct the M matrix used for obtaining eigenvalues and vectors, with a reference model of Normal or Gamma
  M_mat_ref <- matrix(NA, nrow = (m+1), ncol = (m+1)) 
  # Record the log density ratios
  if (ref_model == "normal"){
    ff_tmp_seq <- list() 
    for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
    for (i in 1:length(combined_data)){ 
      x_M <- combined_data[i] 
      log_base_M <- log(dnorm(x = x_M, mean = median(data[[1]]), sd = min(sd(data[[1]]), IQR(data[[1]])/1.34))) 
      for (j in 1:(m+1)){ 
        ff_tmp_seq[[j]][i] <- log(dnorm(x = x_M, mean = median(data[[j]]), sd = min(sd(data[[j]]), IQR(data[[j]])/1.34))) - log_base_M 
      }  
    }  
  }
  
  if (ref_model == "gamma"){
    ff_tmp_seq <- list() 
    for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
    for (i in 1:length(combined_data)){ 
      x_M <- combined_data[i] 
      log_base_M <- log(dgamma(x = x_M, shape = mean(data[[1]])^2/var(data[[1]]), scale = var(data[[1]])/mean(data[[1]]))) 
      for (j in 1:(m+1)){ 
        ff_tmp_seq[[j]][i] <- log(dgamma(x = x_M, shape = mean(data[[j]])^2/var(data[[j]]), scale = var(data[[j]])/mean(data[[j]]))) - log_base_M 
      }  
    }  
  }
  
  # Obtain the centered log density ratio 
  f_sum_M <- numeric(length(combined_data))  
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){ 
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j] 
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]  
  }  
  f_avg_M <- f_sum_M/(m+1) 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M 
  
  # Obtain the matrix M_mat and its eigen-system
  for (i_1 in 1:(m+1)){ 
    for (i_2 in i_1:(m+1)){ 
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]  
      M_mat_ref[i_1, i_2] <- mean(ff_tmp_mat) 
      M_mat_ref[i_2, i_1] <- M_mat_ref[i_1, i_2] 
    }
  }    
  eigen_sys_ref <- eigen(M_mat_ref) # The (m+1)-th eigenvalue of M_mat will always be 0 
  
  # Evaluate and store the values of the basis function at every data points to avoid excessive amount of basis function calculations when fitting the DRM
  eigenfun_unique_data_ref <- matrix(NA, nrow = length(unique(combined_data)), ncol = 2) 
  for (ii in 1:length(unique(combined_data))){
    x <- unique(combined_data)[ii]
    x_loc <- which(combined_data == x)[1] #only use the first one to avoid duplicate issue
    log_DR <- numeric(m+1) 
    for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc] #already calculated and recorded when obtaining M_mat
    # obtain the estimated eigenfunction evaluated at the unique data with location x_loc
    fpca_basis_fn <- numeric(2)  
    for (jj in 1:2) fpca_basis_fn[jj] <- (eigen_sys_ref$values[jj]^(-1))*sum((eigen_sys_ref$vectors[(1:(m+1)),jj])*log_DR) 
    eigenfun_unique_data_ref[ii,] <- fpca_basis_fn
  }
  
  
  #######################################################################################
  # Find the adaptive kernel bandwidth
  obj_fun <- function(bw_factor){
    if ((min(bw_factor) < 1) | (max(bw_factor) > 5)){
      return(Inf)
    } else{
      # Construct the M matrix used for obtaining eigenvalues and vectors, using KDE
      kernel_bw <- bw_factor*bw_RoT
      M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1)) 
      # Record the log density ratios
      ff_tmp_seq <- list() 
      for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
      for (i in 1:length(combined_data)){ 
        x_M <- combined_data[i] 
        log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
        for (j in 1:(m+1)){ 
          ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M 
        }  
      } 
      
      # Obtain the centered log density ratio 
      f_sum_M <- numeric(length(combined_data))  
      mean_log_DR <- numeric(m+1)
      for (j in 1:(m+1)){ 
        mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
        ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j] 
        f_sum_M <- f_sum_M + ff_tmp_seq[[j]]  
      }  
      f_avg_M <- f_sum_M/(m+1) 
      for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M 
      
      # Obtain the matrix M_mat and its eigen-system
      for (i_1 in 1:(m+1)){ 
        for (i_2 in i_1:(m+1)){ 
          ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]  
          M_mat[i_1, i_2] <- mean(ff_tmp_mat) 
          M_mat[i_2, i_1] <- M_mat[i_1, i_2] 
        }
      }    
      eigen_sys <- eigen(M_mat) # The (m+1)-th eigenvalue of M_mat will always be 0 
      
      # Evaluate and store the values of the basis function at every data points to avoid excessive amount of basis function calculations 
      eigenfun_unique_data <- matrix(NA, nrow = length(unique(combined_data)), ncol = 2) 
      for (ii in 1:length(unique(combined_data))){
        x <- unique(combined_data)[ii]
        x_loc <- which(combined_data == x)[1] #only use the first one to avoid duplicate issue
        log_DR <- numeric(m+1) 
        for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc] #already calculated and recorded when obtaining M_mat
        # obtain the estimated eigenfunction evaluated at the unique data with location x_loc
        fpca_basis_fn <- numeric(2)
        for (jj in 1:2) fpca_basis_fn[jj] <- (eigen_sys$values[jj]^(-1))*sum((eigen_sys$vectors[(1:(m+1)),jj])*log_DR) 
        eigenfun_unique_data[ii,] <- fpca_basis_fn
      }
      
      return(sum(colMeans((eigenfun_unique_data - eigenfun_unique_data_ref)^2)))
    }
  }
  
  # Optimization step
  opt <- optimize(f = obj_fun, interval = c(1, 5))
  adaptive_bw_factor <- opt$minimum
  
  return(adaptive_bw_factor*bw_RoT)
}
####################################################################################################################################



####################################################################################################################################
# adaptive_dim: find the adaptive dimension of the adaptive basis function
library(drmdel)

# Inputs: 
# data: a list that contains m+1 samples, where list[[i]] should contain a vector representing the (i+1)-th sample.

# kernel_bw: a vector of the kernel bandwidths used in the step of estimation of eigenfunctions, whose length should be the same as the number of populations

# Output: 
# a scalar that gives the adaptive dimension of the adaptive basis function

adaptive_dim <- function(data, kernel_bw){
  m <- length(data) - 1 # m+1 is the total number of populations 
  sizes <- numeric(m+1) # sample sizes 
  combined_data <- NULL # the pooled sample
  for (i in 1:(m+1)){
    sizes[i] <- length(data[[i]])
    combined_data <- c(combined_data, data[[i]])
  }
  n <- sum(sizes)
  
  # aggregate the same population data at their unique values to avoid excessive computational cost 
  unique_data <- list() 
  for (i in 1:(m+1)){ 
    unique_pt <- unique(data[[i]]) 
    num_unique <- length(unique_pt)
    # the first column is the value of the unique points, the second column is the number of obs that equal to that value  
    unique_data[[i]] <- matrix(nrow = num_unique, ncol = 2) 
    for (j in 1:num_unique){ 
      unique_data[[i]][j,1] <- unique_pt[j] 
      unique_data[[i]][j,2] <- sum(data[[i]] == unique_pt[j]) 
    } 
  } 
  
  
  #######################################################################################
  # Construct the M matrix used for obtaining eigenvalues and vectors, using KDE
  M_mat <- matrix(NA, nrow = (m+1), ncol = (m+1)) 
  # Record the log density ratios
  ff_tmp_seq <- list() 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- numeric(length(combined_data))  
  for (i in 1:length(combined_data)){ 
    x_M <- combined_data[i] 
    log_base_M <- log(sum(unique_data[[1]][,2]*dnorm(x = (x_M-unique_data[[1]][,1]), mean = 0, sd = kernel_bw[1]))/sum(unique_data[[1]][,2])) 
    for (j in 1:(m+1)){ 
      ff_tmp_seq[[j]][i] <- log(sum(unique_data[[j]][,2]*dnorm(x = (x_M-unique_data[[j]][,1]), mean = 0, sd = kernel_bw[j]))/sum(unique_data[[j]][,2])) - log_base_M 
    }  
  } 
  
  # Obtain the centered log density ratio 
  f_sum_M <- numeric(length(combined_data))  
  mean_log_DR <- numeric(m+1)
  for (j in 1:(m+1)){ 
    mean_log_DR[j] <- mean(ff_tmp_seq[[j]])
    ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - mean_log_DR[j] 
    f_sum_M <- f_sum_M + ff_tmp_seq[[j]]  
  }  
  f_avg_M <- f_sum_M/(m+1) 
  for (j in 1:(m+1)) ff_tmp_seq[[j]] <- ff_tmp_seq[[j]] - f_avg_M 
  
  # Obtain the matrix M_mat and its eigen-system
  for (i_1 in 1:(m+1)){ 
    for (i_2 in i_1:(m+1)){ 
      ff_tmp_mat <- ff_tmp_seq[[i_1]]*ff_tmp_seq[[i_2]]  
      M_mat[i_1, i_2] <- mean(ff_tmp_mat) 
      M_mat[i_2, i_1] <- M_mat[i_1, i_2] 
    }
  }    
  eigen_sys <- eigen(M_mat) # The (m+1)-th eigenvalue of M_mat will always be 0 
  # proportions of explained variation
  eigenvalues <- eigen_sys$values
  eigenvalues_ratio <- numeric(m)
  for (j in 1:m) eigenvalues_ratio[j] <- sum(eigenvalues[1:j])/sum(eigenvalues)
  
  # Evaluate and store the values of the basis function at every data points to avoid excessive amount of basis function calculations 
  eigenfun_unique_data <- matrix(NA, nrow = length(unique(combined_data)), ncol = m) 
  for (ii in 1:length(unique(combined_data))){
    x <- unique(combined_data)[ii]
    x_loc <- which(combined_data == x)[1] #only use the first one to avoid duplicate issue
    log_DR <- numeric(m+1) 
    for (j in 1:(m+1)) log_DR[j] <- ff_tmp_seq[[j]][x_loc] #already calculated and recorded when obtaining M_mat
    # obtain the estimated eigenfunction evaluated at the unique data with location x_loc
    fpca_basis_fn <- numeric(m)
    for (jj in 1:m) fpca_basis_fn[jj] <- (eigen_sys$values[jj]^(-1))*sum((eigen_sys$vectors[(1:(m+1)),jj])*log_DR) 
    eigenfun_unique_data[ii,] <- fpca_basis_fn
  }
  
  
  #######################################################################################
  # fit the DRM to the data using the adaptive basis function of different dimensions
  # create the discrete function of the basis function
  eigenfun_discrete <- function(x, nFPC){ 
    loc_unique_data <- which(x == unique(combined_data))  
    return(eigenfun_unique_data[loc_unique_data,(1:nFPC)]) 
  } 
  
  # record the BIC values of different DRM models
  BIC_candidates <- numeric(m)
  for (j in 1:m){
    drmfit_FPCA <- drmdel(x=combined_data, n_samples=sizes, basis_func=function(x) return(eigenfun_discrete(x, nFPC=j))) 
    likelihood_FPCA <- -drmfit_FPCA$negldl - n*log(n)
    BIC_candidates[j] <- -2*likelihood_FPCA + j*m*log(n)
  }
  
  
  #######################################################################################
  # report the adaptive dimension
  
  # method 1: explained variation
  nFPC_95 <- min(m, which(eigenvalues_ratio >= 0.95)[1])
  
  # method 2: BIC
  nFPC_BIC <- which.min(BIC_candidates)
  
  # Final adaptive method: the maximum of the two methods
  nFPC_adaptive <- max(nFPC_BIC, nFPC_95)
  
  return(nFPC_adaptive)
}
####################################################################################################################################



############################################# End of file ############################################# 


