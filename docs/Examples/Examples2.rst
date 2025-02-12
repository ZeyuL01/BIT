BIT Simulation
====================

This example demonstrates how to generate simulation plots from raw data as presented in the manuscript. In manuscript, we considered three simulation cases examining :math:`\mu`, :math:`\tau^2`, and :math:`\sigma_0^2`, along with two distribution variants: Gamma and Student’s :math:`t` distributions. Below, we focus on the case of :math:`\mu`.

Generating simulation data
------------------------

To generate the raw simulation data, we run the simulation script ``SIMU_MU.R``, First, we import the required packages and configure the test values for :math:`\mu`:

.. code-block:: r

  library(BIT)
  library(truncdist)

  # Configuration parameters ---------------------------------------------------
  TR_SIZES <- c(500, 1000, 1500)
  MU_VALUES <- c(-5, -4.5, -4, -3.5, -3, -2.5)
  TAU_SQ_VALUES <- c(0.5, 0.7, 0.9, 1.1, 1.3, 1.5)
  SIGMA0_VALUES <- c(1, 1.2, 1.4, 1.6, 1.8, 2.0)
  BASE_DIR <- "./simulation/MU"
  N_VALUE <- 30000
  ITERATIONS <- 5000
  BURN_IN <- 2500


Next, we define helper functions used to generate simulation data for TRs with either single or multiple observations (datasets):

.. code-block:: r

  # Helper functions --------------------------------------------------------
  logistic_transform <- function(theta) {
    exp(theta) / (1 + exp(theta))
  }

  generate_TR_size <- function(num_TRs) {
    ceiling(rtrunc(num_TRs, spec = "lnorm",
                   meanlog = 1, sdlog = 1.75,
                   a = 0, b = 1000))
  }

  generate_multi_obs_TR <- function(n, mu, tau_sq, population) {
    TR_data <- list()
    TR_data$theta_i <- rnorm(1, mu, sqrt(tau_sq))
    TR_data$sigma_i <- rgamma(1, 0.75, 100)
    TR_data$theta_ij <- rnorm(n, TR_data$theta_i, sqrt(TR_data$sigma_i))
    TR_data$x_ij <- rbinom(n, population, logistic_transform(TR_data$theta_ij))
    TR_data$n_ij <- rep(population, n)
    TR_data
  }

  generate_single_obs_TR <- function(mu, tau_sq, sigma0_sq, population) {
    TR_data <- list()
    TR_data$theta_i <- rnorm(1, mu, sqrt(tau_sq))
    TR_data$theta_ij <- rnorm(1, TR_data$theta_i, sqrt(sigma0_sq))
    TR_data$x_ij <- rbinom(1, population, logistic_transform(TR_data$theta_ij))
    TR_data$n_ij <- population
    TR_data
  }


With these helper functions in place, we define the main function for generating simulation data:

.. code-block:: r

  # Data generation functions -----------------------------------------------
  generate_simulation_data <- function(
      num_TRs = 1000,
      population = 30000,
      mu = -4,
      tau_sq = 0.75,
      sigma0_sq = 1.5) {

    TR_sizes <- generate_TR_size(num_TRs)

    multi_obs_count <- sum(TR_sizes > 1)
    single_obs_count <- num_TRs - multi_obs_count
    ordered_sizes <- sort(TR_sizes[TR_sizes > 1])

    simulation_data <- vector("list", num_TRs)

    # Generate multi-observation TRs
    for (i in seq_len(multi_obs_count)) {
      simulation_data[[i]] <- generate_multi_obs_TR(
        ordered_sizes[i], mu, tau_sq, population
      )
    }

    # Generate single-observation TRs
    for (j in (multi_obs_count + 1):num_TRs) {
      simulation_data[[j]] <- generate_single_obs_TR(
        mu, tau_sq, sigma0_sq, population
      )
    }

    simulation_data
  }

  structure_simulation_data <- function(raw_data) {
    structured_data <- list(
      xij = unlist(lapply(raw_data, `[[`, "x_ij")),
      nij = unlist(lapply(raw_data, `[[`, "n_ij")),
      label_vec = rep(seq_along(raw_data), lengths(lapply(raw_data, `[[`, "x_ij"))),
      theta_i = unlist(lapply(raw_data, `[[`, "theta_i")),
      theta_ij = unlist(lapply(raw_data, `[[`, "theta_ij"))
    )
    structured_data
  }

The simulation workflow consists of generating the data, running the main analysis, and saving the results:

.. code-block:: r

  # Simulation workflow ----------------------------------------------------
  run_simulation <- function(mu_value, iterations, num_TRs, simulation_id) {
    data_dir <- file.path(BASE_DIR, "SIMU_DATA")
    result_dir <- file.path(BASE_DIR, "SIMU_RESULTS")
    log_dir <- file.path(BASE_DIR, "LOG")

    dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

    # Generate and save simulation data
    simulated_data <- generate_simulation_data(
      num_TRs, N_VALUE, mu_value, 0.75, 1.5
    )
    structured_data <- structure_simulation_data(simulated_data)

    data_path <- file.path(data_dir, sprintf("id_%d_data_sim_mu_%g_I_%d.rds",
                                             simulation_id, mu_value, num_TRs))
    log_path <- file.path(log_dir, sprintf("id_%d_data_sim_mu_%g_I_%d.txt",
                                           simulation_id, mu_value, num_TRs))
    saveRDS(structured_data, data_path)

    # Run main analysis
    analysis_results <- Main_Sampling(
      iterations,
      structured_data$xij,
      structured_data$nij,
      structured_data$label_vec,
      log_path
    )

    # Process and save results
    final_results <- list(
      mu = mean(analysis_results$mu0[(iterations - BURN_IN):iterations]),
      theta_i = rowMeans(analysis_results$theta_i[, (iterations - BURN_IN):iterations]),
      label_vec = structured_data$label_vec
    )

    result_path <- file.path(result_dir, sprintf("id_%d_res_sim_mu_%g_I_%d.rds",
                                                 simulation_id, mu_value, num_TRs))
    saveRDS(final_results, result_path)
  }


We run the simulation 100 times for each setting, using different total numbers of TRs (500, 1000, and 1500), and varying :math:`\mu` range from :math:`[-5,-2.5]`, while keeping :math:`\tau^2=0.75` and :math:`\sigma_0^2=1.5` as default values.

.. code-block:: r

  # Execution block ------------------------------------------------------
  for (mu_value in MU_VALUES){
    for (sim_id in seq_len(100)) {
      for (tr_size in TR_SIZES) {
        run_simulation(
          mu_value = mu_value,
          iterations = ITERATIONS,
          num_TRs = tr_size,
          simulation_id = sim_id
        )
      }
    }
  }









