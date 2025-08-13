# =============================================================================
# CASE I and II(b): HETEROGENEOUS DATA WITH NPIs
# Generates fitting and prediction plots for heterogeneous vs homogeneous models
# =============================================================================

library(tidyverse)
library(deSolve)
library(MASS)
library(gridExtra)
library(grid)

# Source required functions
source("MLE_functions_paper.R")

# Define transformation functions if not already available
if (!exists("logit")) {
  logit <- function(p) log(p/(1-p))
}
if (!exists("expit")) {
  expit <- function(x) 1/(1+exp(-x))
}

cat("=== STREAMLINED CASE II(b): HETEROGENEOUS DATA WITH NPIs ===\n")
sdfhsdf
# =============================================================================
# STEP 1: FIXED PARAMETERS AND SETUP
# =============================================================================

# Fixed model parameters
N <- 100000
I0_fixed <- 40
E0_fixed <- I0_fixed * 2.5
n_replicates <- 200

# Model parameters matching your paper
R0_spec <- 3.0
delta_spec <- 1/5.5
rho_spec <- 0.5
gamma_spec <- 1/4
t0_spec <- 15
t1_spec <- 20
t2_spec <- 99
t3_spec <- 100
tfinal_spec <- 100
c_value1_spec <- 1
c_value2_spec <- 0.3 # NPI strength (30% reduction)
c_value3_spec <- 1

# Case II(b): Simulate heterogeneous data (CV appox  1.414)  CASE set CV_simulation=0
CV_simulation <-sqrt(1/0.5)  #=1.414  case II b. for case I b set value to be 0.
true_CV <- CV_simulation

cat("Simulation parameters:\n")
cat("- CV (heterogeneity):", round(CV_simulation, 3), "\n")
cat("- R0:", R0_spec, "\n")
cat("- NPI strength:", c_value2_spec, "\n")
cat("- I0 fixed:", I0_fixed, "\n")

# Create initial state
initial_state <- c(S = N - E0_fixed - I0_fixed, E = E0_fixed, I = I0_fixed, R = 0, C = 0)
assign("initial_state", initial_state, envir = .GlobalEnv)

# =============================================================================
# STEP 2: DATA SIMULATION  
# =============================================================================

cat("\n=== SIMULATING HETEROGENEOUS DATA WITH NPIs ===\n")

set.seed(12375)
sim_data_all <- NULL

for (j in 1:n_replicates) {
  if (j %% 50 == 0) cat("Simulating dataset", j, "of", n_replicates, "\n")
  set.seed(1000*j+n_replicates)
  # Simulate with heterogeneous susceptibility
  dsimdiscre <- simulate_cases_reduced_model(
    R0 = R0_spec,
    delta = delta_spec,
    rho = rho_spec, 
    gamma = gamma_spec,
    v = CV_simulation,  # Heterogeneous simulation
    N = N,
    E0 = E0_fixed,
    I0 = I0_fixed,
    t0 = t0_spec,
    t1 = t1_spec,
    t2 = t2_spec,
    t3 = t3_spec,
    c_value1 = c_value1_spec,
    c_value2 = c_value2_spec,
    c_value3 = c_value3_spec,
    tfinal = tfinal_spec
  ) 
  
  sim.data_1 <- dsimdiscre$sim_data
  sim.data_1$data_set <- j
  sim_data_all <- if (is.null(sim_data_all)) sim.data_1 else bind_rows(sim_data_all, sim.data_1)
}

# Set global times for fitting functions
times <- sim.data_1$time
assign("times", times, envir = .GlobalEnv)

cat("Data simulation complete:", nrow(sim_data_all), "total observations\n")

#write.csv(sim_data_all,"sim_data_all_hom", row.names = F)

# =============================================================================
# STEP 3: FIT HETEROGENEOUS MODEL (EXACT STRUCTURE MATCH)
# =============================================================================

cat("\n=== FITTING HETEROGENEOUS MODEL ===\n")

if (exists("results_het")) rm(results_het)
het_parameter_estimates <- list()        # Store NATURAL scale parameters for empirical covariance
het_trans_parameter_estimates <- list()  # Store TRANSFORMED parameters  
het_individual_covariances <- list()     # Store individual covariance matrices

for (j in 1:n_replicates) {
  if (j %% 50 == 0) cat("Het fitting dataset", j, "of", n_replicates, "\n")
  
  # Extract current dataset
  sim.data_1 <- sim_data_all %>% 
    filter(data_set == j) %>% 
    dplyr::select(time, reports)
  
  # Fit the reduced model with error handling
  z_mle <- tryCatch({
    fit4_reducedm_loglik.NPI(dat = sim.data_1)
  }, error = function(e) {
    cat("Error in fitting dataset", j, ":", e$message, "\n")
    return(NULL)
  })
  
  # Skip this dataset if fitting failed
  if (is.null(z_mle)) {
    cat("Skipping dataset", j, "due to fitting error\n")
    next
  }
  
  # Set up results dataframe for this iteration
  i_results <- as.data.frame(matrix(z_mle$parms, nrow = 1))
  colnames(i_results) <- names(z_mle$parms)
  
  # Initialize values for Hessian results
  z_se <- numeric(length(z_mle$trans_parms))
  z_cor <- c(0, 0, 0)
  z_hess <- 0
  z_pd <- 0
  z_ratio <- 0
  
  # Process Hessian matrix if available
  if (!is.null(z_mle$trans_hessian)) {
    tryCatch({
      z_hess <- 1
      z_eigen <- eigen(z_mle$trans_hessian)
      z_ratios <- z_eigen$values[1] / z_eigen$values
      z_ratio <- z_ratios[length(z_ratios)]
      
      if (all(z_eigen$values > 0)) {
        z_pd <- 1
        z_variance <- solve(z_mle$trans_hessian)
        z_d <- diag(1 / sqrt(diag(z_variance)), nrow = nrow(z_variance))
        z_correlation <- z_d %*% (z_variance %*% z_d)
        z_se <- sqrt(diag(z_variance))
        
        # Extract key correlations (R0-v, R0-t0, v-c_value2)
        z_cor <- c(
          z_correlation[2, 1],  # R0_v_trans_cor
          z_correlation[3, 1],  # R0_t0_trans_cor
          z_correlation[4, 2]   # v_c_value2_trans_cor
        )
        
        # Store parameter estimates and covariance for empirical calculation
        param_vector <- c(z_mle$parms[1], z_mle$parms[2], z_mle$parms[3], z_mle$parms[4])
        names(param_vector) <- c("R0", "v", "t0", "c_value2")
        het_parameter_estimates[[j]] <- param_vector
        
        # Store TRANSFORMED parameters 
        trans_param_vector <- z_mle$trans_parms
        names(trans_param_vector) <- c("R0_trans", "v_trans", "t0_trans", "c_value2_trans")
        het_trans_parameter_estimates[[j]] <- trans_param_vector
        
        # Store individual covariance matrix (transformed parameter space)
        het_individual_covariances[[j]] <- z_variance
        
        # Calculate confidence intervals
        par_ucl <- z_mle$trans_parms + 1.96 * sqrt(diag(z_variance))
        par_lcl <- z_mle$trans_parms - 1.96 * sqrt(diag(z_variance))
        
        C_intervals <- as.data.frame(matrix(c(
          exp(par_lcl[1]), exp(par_ucl[1]),          # R0
          exp(par_lcl[2]), exp(par_ucl[2]),          # v
          exp(par_lcl[3]), exp(par_ucl[3]),          # t0
          expit(par_lcl[4]), expit(par_ucl[4])       # c_value2
        ), nrow = 1, byrow = TRUE))
        
        colnames(C_intervals) <- c(
          "R0_lcl", "R0_ucl", 
          "v_lcl", "v_ucl", 
          "t0_lcl", "t0_ucl", 
          "c_value2_lcl", "c_value2_ucl"
        )
      }
    }, error = function(e) {
      # Keep default values if Hessian processing fails
    })
  }
  
  # Create dataframes for transformed parameters
  z1 <- as.data.frame(matrix(z_mle$trans_parms, nrow = 1))
  colnames(z1) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans", sep = "_")
  
  z_se <- as.data.frame(matrix(z_se, nrow = 1))
  colnames(z_se) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans_se", sep = "_")
  
  z_cor <- as.data.frame(matrix(z_cor, nrow = 1))
  colnames(z_cor) <- c("R0_v_trans_cor", "R0_t0_trans_cor", "v_c_value2_trans_cor")
  
  # Combine all results for this iteration
  i_results <- i_results %>% bind_cols(z1, z_se, z_cor)
  i_results <- if (exists("C_intervals")) bind_cols(i_results, C_intervals) else i_results
  
  # Add diagnostic values and dataset ID
  i_results$hess_exists <- z_hess
  i_results$hess_pd <- z_pd
  i_results$ratio_max_min_evalue <- z_ratio
  i_results$dataset_id <- j
  
  # Append to results
  results_het <- if (exists("results_het") && nrow(results_het) > 0) bind_rows(results_het, i_results) else i_results
}



#write.csv(results_het,"results_hetdata_case_II_b_het_parameters.csv", row.names = F)
# Filter valid results
valid_het <- results_het %>% 
  filter(convergence == 0 & hess_pd == 1 & !is.na(R0) & !is.na(v) & !is.na(t0) & !is.na(c_value2))

cat("Het model: Valid results:", nrow(valid_het), "out of", nrow(results_het), "\n")

saveRDS(het_parameter_estimates,"het_parameter_estimates.rds")
saveRDS(het_trans_parameter_estimates,"het_trans_parameter_estimates.rds")
saveRDS(het_individual_covariances,"het_individual_covariances.rds")

# =============================================================================
# STEP 4: FIT HOMOGENEOUS MODEL (EXACT STRUCTURE MATCH)
# =============================================================================

cat("\n=== FITTING HOMOGENEOUS MODEL ===\n")
v_spec=0
if (exists("results_hom")) rm(results_hom)
hom_parameter_estimates <- list()        # Store NATURAL scale parameters for empirical covariance
hom_trans_parameter_estimates <- list()  # Store TRANSFORMED parameters
hom_individual_covariances <- list()     # Store individual covariance matrices

for (i in 1:n_replicates) {
  if (i %% 50 == 0) print(paste("Fitting dataset", i))
  
  # Extract current dataset
  sim.data_1 <- sim_data_all %>% 
    filter(data_set == i) %>% 
    dplyr::select(time, reports)
  
  # Fit the homogeneous model with error handling
  z_mle <- tryCatch({
    fit3_hom_1epic_loglikwithNPI(dat = sim.data_1)
  }, error = function(e) {
    cat("Error in fitting dataset", i, ":", e$message, "\n")
    return(NULL)
  })
  
  # Skip this dataset if fitting failed
  if (is.null(z_mle)) {
    cat("Skipping dataset", i, "due to fitting error\n")
    next
  }
  
  # Set up results dataframe for this iteration
  i_results <- as.data.frame(matrix(z_mle$parms, nrow = 1))
  colnames(i_results) <- names(z_mle$parms)
  
  # Initialize values for Hessian results
  z_se <- numeric(length(z_mle$trans_parms))
  z_cor <- c(0, 0, 0)
  z_hess <- 0
  z_pd <- 0
  z_ratio <- 0
  
  # Process Hessian matrix if available
  if (!is.null(z_mle$trans_hessian)) {
    tryCatch({
      z_hess <- 1
      z_eigen <- eigen(z_mle$trans_hessian)
      z_ratios <- z_eigen$values[1] / z_eigen$values
      z_ratio <- z_ratios[length(z_ratios)]
      
      if (all(z_eigen$values > 0)) {
        print(i)
        
        z_pd <- 1
        # Covariance matrix calculation
        z_variance <- solve(z_mle$trans_hessian)
        # Correlations of transformed parameters
        z_d <- diag(1/sqrt(diag(z_variance)), nrow=nrow(z_variance))
        z_correlation <- z_d %*% (z_variance %*% z_d)
        z_se <- sqrt(diag(z_variance)) # Standard errors for limits
        z_cor <- c(z_correlation[2,1], z_correlation[2,3], z_correlation[1,3])  
        
        # Store parameter estimates for empirical calculation
        param_vector <- c(z_mle$parms[1], z_mle$parms[2], z_mle$parms[3])
        names(param_vector) <- c("R0", "t0", "c_value2")
        hom_parameter_estimates[[i]] <- param_vector
        
        # Store TRANSFORMED parameters
        trans_param_vector <- z_mle$trans_parms
        names(trans_param_vector) <- c("R0_trans", "t0_trans", "c_value2_trans")
        hom_trans_parameter_estimates[[i]] <- trans_param_vector
        
        # Store individual covariance matrix (transformed parameter space)
        hom_individual_covariances[[i]] <- z_variance
        
        # Calculate confidence intervals
        par_ucl <- z_mle$trans_parms + 1.96 * sqrt(diag(z_variance))
        par_lcl <- z_mle$trans_parms - 1.96 * sqrt(diag(z_variance))
        C_intervals <- as.data.frame(matrix(c(
          exp(par_lcl[1]), 
          exp(par_ucl[1]) ,
          exp(par_lcl[2]), exp(par_ucl[2]),
          expit(par_lcl[3]), expit(par_ucl[3])
        ), nrow = 1, byrow = TRUE))
        colnames(C_intervals) <- c("R0_lcl", "R0_ucl", "t0_lcl", "t0_ucl", 
                                   "c_value2_lcl", "c_value2_ucl")
      }
    }, error = function(e) {
      cat("Error processing Hessian for dataset", i, ":", e$message, "\n")
      # Default values for z_se, z_cor, etc. remain unchanged
    })
  }
  
  # Create dataframes for transformed parameters, standard errors, and correlations
  z1 <- as.data.frame(matrix(z_mle$trans_parms, nrow=1))
  colnames(z1) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans", sep="_")
  
  z_se <- as.data.frame(matrix(z_se, nrow=1))
  colnames(z_se) <- paste(names(z_mle$parms)[1:length(z_mle$trans_parms)], "trans_se", sep="_")
  
  z_cor <- as.data.frame(matrix(z_cor, nrow=1))
  colnames(z_cor) <- c("R0_t0_trans_cor", "t0_c_value2_trans_cor", "R0_c_value2_trans_cor")
  
  # Combine results
  i_results <- i_results %>% bind_cols(z1, z_se, z_cor)
  i_results <- if (exists("C_intervals")) bind_cols(i_results, C_intervals) else i_results
  
  # Add diagnostics and dataset ID
  i_results$hess_exists <- z_hess
  i_results$hess_pd <- z_pd
  i_results$ratio_max_min_evalue <- z_ratio
  i_results$dataset_id <- i  # Add dataset ID for tracking
  
  # Append to results
  results_hom <- if (exists("results_hom")) bind_rows(results_hom, i_results) else i_results
}




# SAVE THE HOM MODEL RESULTS


#write.csv(results_hom,"results_hetdata_case_II_b_hom_parameters.csv", row.names = F)

# Filter valid results  
valid_het <- results_het %>% 
  filter(convergence == 0 & hess_pd == 1 & !is.na(R0) & !is.na(v) & !is.na(t0) & !is.na(c_value2))

valid_hom <- results_hom %>% 
  filter(convergence == 0 & hess_pd == 1 & !is.na(R0) & !is.na(t0) & !is.na(c_value2))

cat("Het model: Valid results:", nrow(valid_het), "out of", nrow(results_het), "\n")
cat("Hom model: Valid results:", nrow(valid_hom), "out of", nrow(results_hom), "\n")


saveRDS(hom_parameter_estimates,"hom_parameter_estimates.rds")
saveRDS(hom_trans_parameter_estimates,"hom_trans_parameter_estimates.rds")
saveRDS(hom_individual_covariances,"hom_individual_covariances.rds")


# =============================================================================
# STEP 5: CALCULATE EMPIRICAL COVARIANCES (USING STORED ESTIMATES)
# =============================================================================

cat("\n=== CALCULATING EMPIRICAL COVARIANCES FROM STORED ESTIMATES ===\n")

# Extract valid parameter estimates (NATURAL SCALE for trajectory generation)
valid_het_params <- do.call(rbind, het_parameter_estimates[valid_het$dataset_id])
valid_hom_params <- do.call(rbind, hom_parameter_estimates[valid_hom$dataset_id])

# Extract valid TRANSFORMED parameter estimates (for empirical covariance)
valid_het_trans_params <- do.call(rbind, het_trans_parameter_estimates[valid_het$dataset_id])
valid_hom_trans_params <- do.call(rbind, hom_trans_parameter_estimates[valid_hom$dataset_id])

# Calculate empirical means and covariances IN NATURAL SPACE (for trajectory generation)
het_natural_mean <- colMeans(valid_het_params, na.rm = TRUE)
hom_natural_mean <- colMeans(valid_hom_params, na.rm = TRUE)

# Calculate empirical means and covariances IN TRANSFORMED SPACE (for sampling)
het_trans_empirical_mean <- colMeans(valid_het_trans_params, na.rm = TRUE)
het_trans_empirical_cov <- cov(valid_het_trans_params, use = "complete.obs")

hom_trans_empirical_mean <- colMeans(valid_hom_trans_params, na.rm = TRUE)
hom_trans_empirical_cov <- cov(valid_hom_trans_params, use = "complete.obs")

cat("Het natural scale mean:", round(het_natural_mean, 3), "\n")
cat("Hom natural scale mean:", round(hom_natural_mean, 3), "\n")
cat("Het transformed empirical mean:", round(het_trans_empirical_mean, 3), "\n")
cat("Hom transformed empirical mean:", round(hom_trans_empirical_mean, 3), "\n")

# =============================================================================
# STEP 6: TRAJECTORY GENERATION FUNCTION
# =============================================================================


# 3. CORRECTED GENERATE_TRAJECTORY FUNCTION
generate_trajectory <- function(R0_val, v_val, t0_val, c_val, times_vec) {
  params <- c(
    R0 = R0_val, v = v_val, rho = rho_spec, delta = delta_spec, gamma = gamma_spec,
    N = N, t0 = t0_val, t1 = t1_spec, t2 = t2_spec, t3 = t3_spec,
    c_value1 = c_value1_spec, c_value2 = c_val, c_value3 = c_value3_spec, 
    tfinal = max(times_vec)
  )
  
  # Get initial conditions - check various possible names
  if (exists("E0_1") && exists("I0_1")) {
    E0_use <- E0_1
    I0_use <- I0_1
  } else if (exists("E0_fixed") && exists("I0_fixed")) {
    E0_use <- E0_fixed
    I0_use <- I0_fixed
  } else if (exists("E0") && exists("I0")) {
    E0_use <- E0
    I0_use <- I0
  } else {
    stop("Initial conditions (E0, I0) not found in global environment")
  }
  
  initial_state <- c(S = N - E0_use - I0_use, E = E0_use, I = I0_use, R = 0, C = 0)
  
  # Check if we need to add time 0
  if (min(times_vec) > 0) {
    times_with_zero <- c(0, times_vec)
    added_zero <- TRUE
  } else {
    times_with_zero <- times_vec
    added_zero <- FALSE
  }
  
  out <- as.data.frame(ode(y = initial_state, times = times_with_zero, 
                           func = Reduced.m_intervene, parms = params))
  
  # Calculate daily incidence
  daily_incidence <- c(0, diff(out$C))
  
  # Return only for requested time points
  if (added_zero) {
    return(daily_incidence[-1])  # Remove the value for time 0
  } else {
    return(daily_incidence)
  }
}






# Time vectors
times_plot <- seq(1, tfinal_spec, by = 1)
times_extended <- seq(1, 250, by = 1)

# Generate median fitted trajectories using natural scale means
fitted_het <- generate_trajectory(as.numeric(het_natural_mean[1]), as.numeric(het_natural_mean[2]), 
                                  as.numeric(het_natural_mean[3]),as.numeric( het_natural_mean[4]), times_plot)
fitted_hom <- generate_trajectory(as.numeric(hom_natural_mean[1]), 0, 
                                  as.numeric(hom_natural_mean[2]), as.numeric(hom_natural_mean[3]), times_plot)

# Generate extended trajectories for prediction
fitted_het_extended <- generate_trajectory(as.numeric(het_natural_mean[1]), as.numeric(het_natural_mean[2]), 
                                           as.numeric(het_natural_mean[3]), as.numeric(het_natural_mean[4]), times_extended)
fitted_hom_extended <- generate_trajectory(as.numeric(hom_natural_mean[1]), 0, 
                                           as.numeric(hom_natural_mean[2]),as.numeric( hom_natural_mean[3]), times_extended)

cat("Median fitted trajectories generated\n")

# =============================================================================
# STEP 7: GENERATE CONFIDENCE BANDS (SAMPLING IN TRANSFORMED SPACE)
# =============================================================================

cat("\n=== GENERATING CONFIDENCE BANDS FROM TRANSFORMED SPACE ===\n")

n_samples <- 100  # Adjust for speed

# HETEROGENEOUS MODEL CONFIDENCE BANDS
# Sample from transformed space using empirical covariance
set.seed(456)
het_trans_samples <- mvrnorm(n_samples, mu = het_trans_empirical_mean, Sigma = het_trans_empirical_cov)

het_predictions <- matrix(NA, nrow = length(times_plot), ncol = n_samples)
het_predictions_ext <- matrix(NA, nrow = length(times_extended), ncol = n_samples)
valid_het_samples <- 0

for(i in 1:n_samples) {
  # Transform back to natural scale
  R0_val <-as.numeric(exp(het_trans_samples[i, 1]))      # exp(log(R0))
  v_val <- as.numeric(exp(het_trans_samples[i, 2]))      # exp(log(v))
  t0_val <- as.numeric(exp(het_trans_samples[i, 3]))      # exp(log(t0))
  c_val <- as.numeric(expit(het_trans_samples[i, 4]))     # expit(logit(c))
  
  # Quality control check
  if(R0_val > 0 && R0_val < 3.5 && 
     v_val >= 0 && v_val < 3 &&
     t0_val > 0 && t0_val < 20 &&
     c_val > 0 && c_val <= 1) {
    
    try({
      temp_traj <- generate_trajectory(R0_val, v_val, t0_val, c_val, times_plot)
      temp_traj_ext <- generate_trajectory(R0_val, v_val, t0_val, c_val, times_extended)
      if(all(is.finite(temp_traj)) && all(is.finite(temp_traj_ext))) {
        het_predictions[, i] <- temp_traj
        het_predictions_ext[, i] <- temp_traj_ext
        valid_het_samples <- valid_het_samples + 1
      }
    }, silent = TRUE)
  }
}

# HOMOGENEOUS MODEL CONFIDENCE BANDS
# Sample from transformed space using empirical covariance
hom_trans_samples <- mvrnorm(n_samples, mu = hom_trans_empirical_mean, Sigma = hom_trans_empirical_cov)

hom_predictions <- matrix(NA, nrow = length(times_plot), ncol = n_samples)
hom_predictions_ext <- matrix(NA, nrow = length(times_extended), ncol = n_samples)
valid_hom_samples <- 0

for(i in 1:n_samples) {
  # Transform back to natural scale
  R0_val <- as.numeric(exp(hom_trans_samples[i, 1]))      # exp(log(R0))
  t0_val <- as.numeric(exp(hom_trans_samples[i, 2]))      # exp(log(t0))
  c_val <- as.numeric(expit(hom_trans_samples[i, 3]))     # expit(logit(c))
  
  # Quality control check
  if(R0_val > 0 && R0_val < 10 &&
     t0_val > 0 && t0_val < 50 &&
     c_val > 0 && c_val <= 1) {
    
    try({
      temp_traj <- generate_trajectory(R0_val, 0, t0_val, c_val, times_plot)
      temp_traj_ext <- generate_trajectory(R0_val, 0, t0_val, c_val, times_extended)
      if(all(is.finite(temp_traj)) && all(is.finite(temp_traj_ext))) {
        hom_predictions[, i] <- temp_traj
        hom_predictions_ext[, i] <- temp_traj_ext
        valid_hom_samples <- valid_hom_samples + 1
      }
    }, silent = TRUE)
  }
}

# Calculate CI bounds
het_lower_ci <- apply(het_predictions, 1, quantile, 0.025, na.rm = TRUE)
het_upper_ci <- apply(het_predictions, 1, quantile, 0.975, na.rm = TRUE)
hom_lower_ci <- apply(hom_predictions, 1, quantile, 0.025, na.rm = TRUE)
hom_upper_ci <- apply(hom_predictions, 1, quantile, 0.975, na.rm = TRUE)

# Extended predictions
het_lower_ext <- apply(het_predictions_ext, 1, quantile, 0.025, na.rm = TRUE)
het_upper_ext <- apply(het_predictions_ext, 1, quantile, 0.975, na.rm = TRUE)
hom_lower_ext <- apply(hom_predictions_ext, 1, quantile, 0.025, na.rm = TRUE)
hom_upper_ext <- apply(hom_predictions_ext, 1, quantile, 0.975, na.rm = TRUE)

cat("Valid samples: Het =", valid_het_samples, ", Hom =", valid_hom_samples, "\n")

# =============================================================================
# STEP 8: CREATE FITTING PLOT (IMAGE 1)
# =============================================================================

cat("\n=== CREATING FITTING PLOT ===\n")

# Select representative data (first 3 datasets combined)
plot_data <- sim_data_all %>% 
  filter(data_set %in% c(1, 2, 3), time > 0) %>%
  group_by(time) %>%
  summarise(mean_reports = mean(reports, na.rm = TRUE), .groups = 'drop')

# Colors matching your images
het_color <- "#1f77b4"  # Blue
hom_color <- "#ff7f0e"  # Orange

fitting_plot <- ggplot() +
  # Confidence bands
  geom_ribbon(aes(x = times_plot, ymin = het_lower_ci, ymax = het_upper_ci), 
             fill = het_color, alpha = 0.3) +
  geom_ribbon(aes(x = times_plot, ymin = hom_lower_ci, ymax = hom_upper_ci), 
              fill = hom_color, alpha = 0.3) +
  
  # Fitted lines
  geom_line(aes(x = times_plot, y = fitted_het, color = "Heterogeneous"), size = 1.2) +
  geom_line(aes(x = times_plot, y = fitted_hom, color = "Homogeneous"), size = 1.2) +
  
  # Observed data points
  geom_point(data = plot_data, aes(x = time, y = mean_reports), size = 1, alpha = 0.8) +
  
  # True intervention line
  geom_vline(xintercept = t0_spec, linetype = "dashed", color = "darkgreen", size = 1) +
  
  # Styling
  scale_color_manual(values = c("Heterogeneous" = het_color, "Homogeneous" = hom_color)) +
  labs(
    title = "Case II b: Heterogeneous data with NPIs",
    subtitle = "Models fitted to first 100 days, forecasting next 150 days",
    x = "Time (days)", y = "Daily cases", color = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold")) +
  
  # Parameter annotations
  annotate("text", x = max(times_plot) * 0.6, y = max(plot_data$mean_reports) * 0.9,
           label = paste0("True: R0=", R0_spec, ", CV=", round(true_CV, 3), ", t0=", t0_spec, ", NPI=", c_value2_spec),
           hjust = 0, size = 3, color = "black") +
  annotate("text", x = max(times_plot) * 0.6, y = max(plot_data$mean_reports) * 0.8,
           label = paste0("Het: R0=", round(het_natural_mean[1], 2), ", CV=", round(het_natural_mean[2], 3), 
                          ", t0=", round(het_natural_mean[3], 1), ", NPI=", round(het_natural_mean[4], 3)),
           hjust = 0, size = 3, color = het_color) +
  annotate("text", x = max(times_plot) * 0.6, y = max(plot_data$mean_reports) * 0.7,
           label = paste0("Hom: R0=", round(hom_natural_mean[1], 2), ", t0=", round(hom_natural_mean[2], 1), 
                          ", NPI=", round(hom_natural_mean[3], 3)),
           hjust = 0, size = 3, color = hom_color) +
  annotate("text", x = t0_spec, y = max(plot_data$mean_reports) * 1.05,
           label = "True t0", hjust = 0.5, size = 3, color = "darkgreen")

print(fitting_plot)

# =============================================================================
# STEP 9: CREATE PREDICTION PLOT (IMAGE 2)
# =============================================================================

cat("\n=== CREATING PREDICTION PLOT ===\n")

# Prepare extended plot data with observations for fitting period
extended_plot_data <- data.frame(
  time = times_extended,
  het_fitted = fitted_het_extended,
  hom_fitted = fitted_hom_extended,
  het_lower = het_lower_ext,
  het_upper = het_upper_ext,
  hom_lower = hom_lower_ext,
  hom_upper = hom_upper_ext,
  observed = NA
)

# Add observed data for fitting period
fitting_period <- which(times_extended <= tfinal_spec)
for(i in 1:min(length(fitting_period), nrow(plot_data))) {
  day <- times_extended[fitting_period[i]]
  if(day %in% plot_data$time) {
    extended_plot_data$observed[fitting_period[i]] <- 
      plot_data$mean_reports[plot_data$time == day]
  }
}

prediction_plot <- ggplot(extended_plot_data, aes(x = time)) +
  # Confidence bands
  geom_ribbon(aes(ymin = het_lower, ymax = het_upper), fill = het_color, alpha = 0.3) +
  geom_ribbon(aes(ymin = hom_lower, ymax = hom_upper), fill = hom_color, alpha = 0.3) +
  
  # Fitted lines
  geom_line(aes(y = het_fitted, color = "Heterogeneous"), size = 1.2) +
  geom_line(aes(y = hom_fitted, color = "Homogeneous"), size = 1.2) +
  
  # Observed data points
  geom_point(aes(y = observed), size = 1.2, alpha = 0.8) +
  
  # Forecast line
  geom_vline(xintercept = tfinal_spec, linetype = "dotted", color = "purple", size = 1) +
  
  # Styling
  scale_color_manual(values = c("Heterogeneous" = het_color, "Homogeneous" = hom_color)) +
  labs(
    title = "Case II (b): Prediction trajectories for Heterogeneous data with NPIs",
    subtitle = "Models fitted to first 100 days, forecasting next 150 days",
    x = "Time (days)", y = "Daily Cases", color = "Model"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(size = 14, face = "bold")) +
  annotate("text", x = tfinal_spec, y = max(extended_plot_data$hom_upper, na.rm = TRUE) * 0.8,
           label = "Forecast begins", hjust = 0, vjust = 0, angle = 90, size = 3, color = "purple")

print(prediction_plot)

# =============================================================================
# STEP 10: SAVE PLOTS
# =============================================================================

cat("\n=== SAVING PLOTS ===\n")

# Create figures directory if it doesn't exist
dir.create("figures", showWarnings = FALSE)

# Save plots
ggsave("figures/Case_IIaii_fitting_plot.png", fitting_plot, width = 12, height = 8, dpi = 300)
ggsave("figures/Case_IIb_prediction_plot.png", prediction_plot, width = 12, height = 8, dpi = 300)

# Save combined plot
combined_plots <- grid.arrange(fitting_plot, prediction_plot, ncol = 1)
ggsave("figures/Case_IIb_combined_plots.png", combined_plots, width = 12, height = 16, dpi = 300)

print(combined_plots)
# ===============================================================================
# PLOTTING CODE -
# ===============================================================================

# Load required libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)

# ===============================================================================
# PART 1: CONVERT LISTS TO DATAFRAMES
# ===============================================================================

# Convert heterogeneous parameter estimates list to dataframe
het_params_df <- data.frame()
for (i in 1:length(het_parameter_estimates)) {
  if (!is.null(het_parameter_estimates[[i]])) {
    temp_row <- data.frame(
      dataset_id = i,
      R0 = het_parameter_estimates[[i]]["R0"],
      v = het_parameter_estimates[[i]]["v"],
      t0 = het_parameter_estimates[[i]]["t0"],
      c_value2 = het_parameter_estimates[[i]]["c_value2"]
    )
    het_params_df <- rbind(het_params_df, temp_row)
  }
}

# Convert homogeneous parameter estimates list to dataframe
hom_params_df <- data.frame()
for (i in 1:length(hom_parameter_estimates)) {
  if (!is.null(hom_parameter_estimates[[i]])) {
    temp_row <- data.frame(
      dataset_id = i,
      R0 = hom_parameter_estimates[[i]]["R0"],
      t0 = hom_parameter_estimates[[i]]["t0"],
      c_value2 = hom_parameter_estimates[[i]]["c_value2"]
    )
    hom_params_df <- rbind(hom_params_df, temp_row)
  }
}

# ===============================================================================
# PART 2: FILTER VALID RESULTS (adjust criteria as needed)
# ===============================================================================

# Filter valid results - adjust these criteria based on your needs
valid_het <- het_params_df %>%
  filter(!is.na(R0), !is.na(v), !is.na(t0), !is.na(c_value2),
         R0 > 0, v >= 0, t0 > 0, c_value2 >= 0, c_value2 <= 1)

valid_hom <- hom_params_df %>%
  filter(!is.na(R0), !is.na(t0), !is.na(c_value2),
         R0 > 0, t0 > 0, c_value2 >= 0, c_value2 <= 1)

# Print summary
cat("Valid heterogeneous fits:", nrow(valid_het), "out of", nrow(het_params_df), "\n")
cat("Valid homogeneous fits:", nrow(valid_hom), "out of", nrow(hom_params_df), "\n")


# ===============================================================================
# PART 3: CREATE DENSITY PLOTS FOR PARAMETER DISTRIBUTIONS
# ===============================================================================

# Prepare data for density plots - combine both models
# Add model type column to each dataset
valid_het_plot <- valid_het %>%
  mutate(Model = "Heterogeneous") %>%
  dplyr:: select(R0, v, t0, c_value2, Model)

# For homogeneous, add v = 0 for consistency
valid_hom_plot <- valid_hom %>%
  mutate(Model = "Homogeneous", v = 0) %>%
 dplyr:: select(R0, v, t0, c_value2, Model)

# Combine data
combined_data <- rbind(valid_het_plot, valid_hom_plot)

# Define true values (adjust these to your actual true values)
true_R0 <- R0_spec
true_v <- CV_simulation  
true_t0 <- t0_spec
true_c <- c_value2_spec

# Create individual density plots
# Colors matching your images
het_color <- "#1f77b4"  # Blue
hom_color <- "#ff7f0e"  # Orange

# R0 density plot
R0_plot <- ggplot(combined_data, aes(x = R0, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = true_R0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = true_R0, y = Inf, label = paste("True =", true_R0), 
           vjust = 2, hjust = 1.1, size = 4) +
  scale_fill_manual(values = c("Heterogeneous" = het_color, "Homogeneous" = hom_color)) +
  theme_minimal() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "R0", x = "R0", y = "Density")

# CV density plot (only heterogeneous model has varying CV)
CV_plot <- ggplot(valid_het_plot, aes(x = v)) +
  geom_density(fill = het_color, alpha = 0.7) +
  geom_vline(xintercept = true_v, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = true_v, y = Inf, label = paste("True =", true_v), 
           vjust = 2, hjust = 1.1, size = 4) +
  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "CV", x = "CV", y = "Density")

# t0 density plot
t0_plot <- ggplot(combined_data, aes(x = t0, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = true_t0, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = true_t0, y = Inf, label = paste("True =", true_t0), 
           vjust = 2, hjust = -0.1, size = 4) +
  scale_fill_manual(values = c("Heterogeneous" = het_color, "Homogeneous" = hom_color)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "t0", x = "t0 (days)", y = "Density")

# NPI (c_value2) density plot
NPI_plot <- ggplot(combined_data, aes(x = c_value2, fill = Model)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = true_c, linetype = "dashed", color = "black", size = 1) +
  annotate("text", x = true_c, y = Inf, label = paste("True =", true_c), 
           vjust = 2, hjust = 1.1, size = 4) +
  scale_fill_manual(values = c("Heterogeneous" = het_color, "Homogeneous" = hom_color)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "NPI", x = "NPI", y = "Density")

# Combine density plots
density_combined <- grid.arrange(R0_plot, CV_plot, t0_plot, NPI_plot, 
                                 ncol = 2, nrow = 2,
                                 top = "Case II(aii): Heterogeneous data with NPIs")





density_combined <- grid.arrange(R0_plot, CV_plot, t0_plot, NPI_plot, 
                                 ncol = 2, nrow = 2,
                                 top = textGrob("Case II(aii): Heterogeneous data without NPIs", 
                                                gp = gpar(fontsize = 16, fontface = "bold")))

# Save density plots
ggsave("figures/density_plots_case_IIaii.png", density_combined, width = 12, height = 10, dpi = 300)

# ===============================================================================
#  PRINT SUMMARY STATISTICS
# ===============================================================================

# Summary for heterogeneous model
cat("\n=== HETEROGENEOUS MODEL SUMMARY ===\n")
cat("R0: mean =", mean(valid_het$R0), ", sd =", sd(valid_het$R0), "\n")
cat("CV: mean =", mean(valid_het$v), ", sd =", sd(valid_het$v), "\n")
cat("t0: mean =", mean(valid_het$t0), ", sd =", sd(valid_het$t0), "\n")
cat("NPI: mean =", mean(valid_het$c_value2), ", sd =", sd(valid_het$c_value2), "\n")

# Summary for homogeneous model
cat("\n=== HOMOGENEOUS MODEL SUMMARY ===\n")
cat("R0: mean =", mean(valid_hom$R0), ", sd =", sd(valid_hom$R0), "\n")
cat("t0: mean =", mean(valid_hom$t0), ", sd =", sd(valid_hom$t0), "\n")
cat("NPI: mean =", mean(valid_hom$c_value2), ", sd =", sd(valid_hom$c_value2), "\n")



# ===============================================================================
#  part 4:CORRELATION HEATMAP CODE
# 
# ===============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)

# ===============================================================================
# STEP 1: PREPARE DATA 
# ===============================================================================


# ===============================================================================
# STEP 2: CREATE CORRELATION MATRICES
# ===============================================================================

# HETEROGENEOUS MODEL
# Calculate correlation matrix
het_cor_matrix <- valid_het %>% 
  dplyr::select(R0, v, t0, c_value2) %>%
  cor(use = "complete.obs")

# Rename columns and rows for display
colnames(het_cor_matrix) <- c("R0", "CV", "t0", "NPI")
rownames(het_cor_matrix) <- c("R0", "CV", "t0", "NPI")

# Convert to long format for ggplot
het_cor_df <- as.data.frame(het_cor_matrix)
het_cor_df$Var1 <- rownames(het_cor_df)
het_cor_df <- het_cor_df %>%
  pivot_longer(cols = -Var1, names_to = "Var2", values_to = "Correlation") %>%
  mutate(
    Var1 = factor(Var1, levels = c("R0", "CV", "t0", "NPI")),
    Var2 = factor(Var2, levels = c("R0", "CV", "t0", "NPI"))
  )

# HOMOGENEOUS MODEL
# Calculate correlation matrix
hom_cor_matrix <- valid_hom %>% 
  dplyr::select(R0, t0, c_value2) %>%
  cor(use = "complete.obs")

# Rename columns and rows for display
colnames(hom_cor_matrix) <- c("R0", "t0", "NPI")
rownames(hom_cor_matrix) <- c("R0", "t0", "NPI")

# Convert to long format for ggplot
hom_cor_df <- as.data.frame(hom_cor_matrix)
hom_cor_df$Var1 <- rownames(hom_cor_df)
hom_cor_df <- hom_cor_df %>%
  pivot_longer(cols = -Var1, names_to = "Var2", values_to = "Correlation") %>%
  mutate(
    Var1 = factor(Var1, levels = c("R0", "t0", "NPI")),
    Var2 = factor(Var2, levels = c("R0", "t0", "NPI"))
  )

# ===============================================================================
#  CREATE CORR HEATMAPS
# ===============================================================================

# HETEROGENEOUS MODEL HEATMAP
het_corr_plot <- ggplot(het_cor_df, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Correlation, 3)), 
            color = "white", size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#053061",    # Dark blue
    mid = "white",      # White
    high = "#67001f",   # Dark red
    midpoint = 0, 
    limit = c(-1, 1), 
    space = "Lab",
    name = "Correlation"
  ) +
  labs(
    #title = "Heterogeneous Model Parameter Correlations",
    subtitle = "Case II (b): Heterogeneous data with NPIs",
    x = "Parameter", 
    y = "Parameter"
  ) +
  theme_minimal() +
  theme(
    # Axis text
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, 
                               size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    # Remove grid and ticks
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    # Legend
    legend.position = "none"  # We'll add a shared legend later
  ) +
  coord_fixed()  # Ensures square tiles

# HOMOGENEOUS MODEL HEATMAP
hom_corr_plot <- ggplot(hom_cor_df, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Correlation, 3)), 
            color = "white", size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#053061",    # Dark blue
    mid = "white",      # White
    high = "#67001f",   # Dark red
    midpoint = 0, 
    limit = c(-1, 1), 
    space = "Lab",
    name = "Correlation"
  ) +
  labs(
   # title = "Homogeneous Model Parameter Correlations",
    subtitle = "Case II (b): Homogeneous data with NPIs",
    x = "Parameter", 
    y = "Parameter"
  ) +
  theme_minimal() +
  theme(
    # Axis text
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, 
                               size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    # Titles
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    # Remove grid and ticks
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    # Legend
    legend.position = "none"
  ) +
  coord_fixed()  # Ensures square tiles



# ===============================================================================
# COMBINE HEATMAPS
# ===============================================================================

# Create a plot with legend to extract it
legend_plot <- ggplot(het_cor_df, aes(x = Var2, y = Var1, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#053061", 
    mid = "white", 
    high = "#67001f",
    midpoint = 0, 
    limit = c(-1, 1), 
    space = "Lab",
    name = "Correlation"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.width = unit(2, "cm")
  )

# Extract legend
library(cowplot)
legend <- get_legend(legend_plot)

# Combine plots without individual legends
combined_corr <- grid.arrange(
  het_corr_plot,
  hom_corr_plot,
  ncol = 2,
  top = textGrob("Case II (b): Parameter Correlation Heatmaps", 
                 gp = gpar(fontsize = 16, fontface = "bold"))
)

# Final combined plot with shared legend
final_corr_plot <- grid.arrange(
  combined_corr,
  legend,
  heights = c(0.9, 0.1),
  nrow = 2
)



# Display the plots
print(het_corr_plot)
print(hom_corr_plot)
print(final_corr_plot)


