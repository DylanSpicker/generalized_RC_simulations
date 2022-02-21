# _targets.R file
library(targets)
library(tarchetypes)

# Load in User Defined Functions
source("R/functions.R")
source("R/corrections/generalized_rc.R")
source("R/corrections/generalized_simex.R")
source("R/corrections/standard_simex.R")
source("R/corrections/standard_rc.R")

source("R/simulation_1/data_generation.R")
source("R/simulation_1/model_fit.R")
source("R/simulation_1/summary.R")

source("R/simulation_2/data_generation.R")
source("R/simulation_2/model_fit.R")
source("R/simulation_2/summary.R")

source("R/simulation_3/data_generation.R")
source("R/simulation_3/model_fit.R")
source("R/simulation_3/summary.R")

# Target Options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = c("MASS", "tidyverse", "rootSolve", "ggthemes", "grid", "gridExtra"), "error" = "continue")

## Global Options
total_replicates <- 10
n_batches <- 5

# Simulation Parameters
parms_1 <- tibble::tibble(
    scenario = c("n=100", "n=500", "n=1000", "n=5000", "n=10000"),
    n = c(100, 500, 1000, 5000, 10000)
)

parms_2 <- tibble::tibble(
    scenario = c("L,L", "L,M", "L,H", "M,L", "M,M", "M,H", "H,L", "H,M", "H,H"),
    n = 5000,
    ev1 = c(0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5),
    ev2 = c(0.2, 0.5, 1, 0.2, 0.5, 1, 0.2, 0.5, 1)
)

parms_3 <- tibble::tibble(
    scenario = c("U,U", "U,P", "U,N", "P,U", "P,P", "P,N", "N,U", "N,P", "N,N"),
    n = 5000,
    eta0 = c(0,0,0,0.5,0.5,0.5,-0.5,-0.5,-0.5),
    eta1 = c(0,0.5,-0.5,0,0.5,-0.5,0,0.5,-0.5)
)

experiment_1_test <- function(sim1_data) {
    sim1_data[1,]
}

list(
    # Pipeline Calls to Simulate the Datasets for each Replicated Dataset
    tar_map_rep(sim1_data, command = generate_data_simulation_1(n), batches = n_batches, reps = floor(total_replicates/n_batches), names = tidyselect::any_of("scenario"), values = parms_1),
    tar_map_rep(sim2_data, command = generate_data_simulation_2(n, ev1, ev2), batches = n_batches, reps = floor(total_replicates/n_batches), names = tidyselect::any_of("scenario"), values = parms_2),
    tar_map_rep(sim3_data, command = generate_data_simulation_3(n, eta0, eta1), batches = n_batches, reps = floor(total_replicates/n_batches), names = tidyselect::any_of("scenario"), values = parms_3),

    # Pipeline Calls to Analyze the Datasets for each Replicated Dataset
    # First Simulation Analyses
    tar_rep2(sim1_results_RC, experiment_1_RC(sim1_data), sim1_data),
    tar_rep2(sim1_results_generalized_RC, experiment_1_generalized_RC(sim1_data), sim1_data),
    tar_rep2(sim1_results_weighted_generalized_RC, experiment_1_weighted_generalized_RC(sim1_data), sim1_data),
    tar_rep2(sim1_results_SIMEX, experiment_1_SIMEX(sim1_data), sim1_data),
    tar_rep2(sim1_results_ESIMEX, experiment_1_ESIMEX(sim1_data), sim1_data),
    tar_rep2(sim1_results_generalized_SIMEX, experiment_1_generalized_SIMEX(sim1_data), sim1_data),
    tar_rep2(sim1_results_combine_first_generalized_SIMEX, experiment_1_combine_first_generalized_SIMEX(sim1_data), sim1_data),

    # Second Simulation Analyses
    tar_rep2(sim2_results_truth, experiment_2_truth(sim2_data), sim2_data),
    tar_rep2(sim2_results_RC, experiment_2_RC(sim2_data), sim2_data),
    tar_rep2(sim2_results_generalized_RC, experiment_2_generalized_RC(sim2_data), sim2_data),
    tar_rep2(sim2_results_generalized_RC_noZ, experiment_2_generalized_RC_noZ(sim2_data), sim2_data),
    tar_rep2(sim2_results_generalized_RC_wtd, experiment_2_generalized_RC_wtd(sim2_data), sim2_data),
    tar_rep2(sim2_results_generalized_RC_wtd_noZ, experiment_2_generalized_RC_wtd_noZ(sim2_data), sim2_data),
    tar_rep2(sim2_results_standard_simex, experiment_2_standard_simex(sim2_data), sim2_data),
    tar_rep2(sim2_results_empirical_simex, experiment_2_empirical_simex(sim2_data), sim2_data),
    tar_rep2(sim2_results_generalized_SIMEX, experiment_2_generalized_SIMEX(sim2_data), sim2_data),
    tar_rep2(sim2_results_generalized_SIMEX_combine_first, experiment_2_generalized_SIMEX_combine_first(sim2_data), sim2_data),

    # Third Simulation Analyses
    tar_rep2(sim3_results_all, experiment_3_all(sim3_data), sim3_data),
    tar_rep2(sim3_results_reps, experiment_3_reps(sim3_data), sim3_data),
    tar_rep2(sim3_results_all_wt, experiment_3_all_wt(sim3_data), sim3_data),
    tar_rep2(sim3_results_reps_wt, experiment_3_reps_wt(sim3_data), sim3_data),
    tar_rep2(sim3_results_generalized_sim, experiment_3_generalized_sim(sim3_data), sim3_data),
    tar_rep2(sim3_results_generalized_sim_reps, experiment_3_generalized_sim_reps(sim3_data), sim3_data),
    tar_rep2(sim3_results_generalized_sim_combine_first, experiment_3_generalized_sim_combine_first(sim3_data), sim3_data),
    tar_rep2(sim3_results_generalized_sim_combine_first_reps, experiment_3_generalized_sim_combine_first_reps(sim3_data), sim3_data),
    tar_rep2(sim3_results_simex, experiment_3_simex(sim3_data), sim3_data),
    tar_rep2(sim3_results_empirical_simex, experiment_3_empirical_simex(sim3_data), sim3_data),

    # Pipeline Calls to Summarize the Simulated Dataset Results
    # Produces Plots and Summaries for Simulation 1
    tar_target(sim1_results, combine_results_1(
        sim1_results_RC, 
        sim1_results_generalized_RC,
        sim1_results_weighted_generalized_RC,
        sim1_results_SIMEX,
        sim1_results_ESIMEX,
        sim1_results_generalized_SIMEX,
        sim1_results_combine_first_generalized_SIMEX
    )),
    tar_target(sim1_plots, generate_plot_1(sim1_results), format = "file"),
    tar_target(sim1_tables, generate_table_1(sim1_results)),

    # Produces Plots and Summaries for Simulation 2
    tar_target(sim2_results, combine_results_2(
        sim2_results_truth,
        sim2_results_RC,
        sim2_results_generalized_RC,
        sim2_results_generalized_RC_noZ,
        sim2_results_generalized_RC_wtd,
        sim2_results_generalized_RC_wtd_noZ,
        sim2_results_standard_simex,
        sim2_results_empirical_simex,
        sim2_results_generalized_SIMEX,
        sim2_results_generalized_SIMEX_combine_first
    )),
    tar_target(sim2_plots, generate_plot_2(sim2_results), format = "file"),
    tar_target(sim2_tables, generate_table_2(sim2_results)),

    # Produce Plots and Summaries for Simulation 3
    tar_target(sim3_results, combine_results_3(
        sim3_results_all,
        sim3_results_reps,
        sim3_results_all_wt,
        sim3_results_reps_wt,
        sim3_results_generalized_sim,
        sim3_results_generalized_sim_reps,
        sim3_results_generalized_sim_combine_first,
        sim3_results_generalized_sim_combine_first_reps,
        sim3_results_simex,
        sim3_results_empirical_simex
    )),
    tar_target(sim3_plots_bp, generate_plot_3_boxplot(sim3_results), format = "file"),
    tar_target(sim3_plots_pp, generate_plot_3_probplot(sim3_results), format = "file")
)