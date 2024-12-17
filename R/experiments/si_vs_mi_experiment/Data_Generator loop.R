#================================= Create DataSet =================================
library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
library(miceadds)
source("si_vs_mi_experiment/cohort_sim_utils.R")

config_name <- "config1n3000"
n_datasets <- 100
m <- 1000
for (i in 1:n_datasets) {
    set.seed(i)
    # generate a cohort DAG and synthetic data
    cohort <- simulate_ampute_cohort(
      n = 3000,
      p = 33,
      num_tiers = 5,
      in_tier_prob = 0.1,
      cross_tier_prob = 0.15,
      lB = 0.1, 
      uB = 1,
      ampute_proportion_total = 0.9,
      ampute_proportion_tier = 0.1,
      ampute_proportion_dropout = 0.1,
      n_dropout_tiers_to_be_excluded = 1,
      n_amputation_dependencies = 2
    )

    data <- cohort$data
    amputed_data <- cohort$amputed_data
    ordering <- cohort$tier_ord
    true_dag <- cohort$dag
    true_tmpdag <- cohort$tmpdag
    dataset_filename <- paste0("si_vs_mi_experiment/data/", config_name, "/dataset", i, ".csv")
    write.table(data, file = dataset_filename, row.names = FALSE, na = "", col.names = FALSE, sep = ",")

    dataset_path <- file.path(dataset_filename, fsep = .Platform$file.sep)
    imputed_data <- mice(amputed_data, m = m, method = "norm", printFlag = FALSE, remove.collinear = TRUE)

    # write mids obj
    save(imputed_data, file = paste0("si_vs_mi_experiment/imputed_datasets_cohort_", config_name, "imp", i, ".Rdata")
}
