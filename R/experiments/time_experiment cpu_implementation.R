library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
library(micd)
library(miceadds)
library(dplyr)
source("R/MeeksRules.R")
source("R/tcheckTriple.R")
source("R/tpc_cons_intern.R")
source("R/tskeleton.R")
source("R/tskeleton_parallel.R")
source("R/cuTPCMI.R")
source("R/fixedMItest.R")
source("R/cuPCMI.R")
library(ggplot2)

# Define dataset paths
dataset_path_100 <- "dataset_imputed_100/dataset_imputed_100.Rdata"
dataset_path_90 <- "dataset_imputed_90/dataset_imputed_90.Rdata"
dataset_path_80 <- "dataset_imputed_80/dataset_imputed_80.Rdata"
dataset_path_70 <- "dataset_imputed_70/dataset_imputed_70.Rdata"
dataset_path_60 <- "dataset_imputed_60/dataset_imputed_60.Rdata"
dataset_path_50 <- "dataset_imputed_50/dataset_imputed_50.Rdata"
dataset_path_25 <- "dataset_imputed_25/dataset_imputed_25.Rdata"
dataset_path_10 <- "dataset_imputed_10/dataset_imputed_10.Rdata"

datasets <- list(
  list(path = dataset_path_100, p = 100),
  list(path = dataset_path_90, p = 90),
  list(path = dataset_path_80, p = 80),
  list(path = dataset_path_70, p = 70),
  list(path = dataset_path_60, p = 60),
  list(path = dataset_path_50, p = 50),
  list(path = dataset_path_25, p = 25),
  list(path = dataset_path_10, p = 10)
)

# Experimental parameters
alphas <- c(0.01, 0.05, 0.1)
m_values <- c(1, 10, 100, 1000)
m_max <<- max(m_values)
repeats <- 10  # Number of repetitions for each run

# Initialize a data frame to store runtime results
results <- data.frame(
  p = integer(),
  m = integer(),
  alpha = numeric(),
  run_time = numeric()
)

# Function to run the experiment
run_experiment <- function(dataset, p, alpha, m) {
  # Load data
  imputed_data <- load.Rdata(dataset, "imputed_data")
  
  suffStatMI_full <- micd::getSuff(imputed_data, test="gaussMItest")
  if (m > 1){
    suffStatMI <- suffStatMI_full[-(1:(m_max - m))]
  }
  else {
    suffStatMI <- list(C = suffStatMI_full[[1]], n = suffStatMI_full[[length(suffStatMI_full)]])
  }
  for (i in 1:repeats) {
    print(format(c(p,alpha,m), scientific=F))
    start_time <- Sys.time()
    cuda_tPC_MI <- pc(suffStatMI, indepTest = gaussMItest, p = p, alpha = alpha, m.max = 12)
    end_time <- Sys.time()
    run_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    results <<- rbind(results, data.frame(p = p, m = m, alpha = alpha, run_time = run_time))
  }
}

# Conduct the experiment
for (dataset_info in datasets) {
  p <- dataset_info$p
  dataset_path <- dataset_info$path
  for (alpha in alphas) {
    for (m in m_values) {
      run_experiment(dataset_path, p, alpha, m)
    }
  }
}

# Calculate mean and standard deviation of runtimes
summary_results <- aggregate(run_time ~ p + m + alpha, data = results, FUN = function(x) c(mean = mean(x), sd = sd(x)))
summary_results <- do.call(data.frame, summary_results)
colnames(summary_results) <- c("p", "m", "alpha", "mean_runtime", "sd_runtime")

# Plot results with ggplot2
# Plot: Runtime vs m for each p and alpha
ggplot(summary_results, aes(x = log10(m), y = log10(mean_runtime), color = factor(p), shape = factor(alpha))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(p, alpha))) +
  #geom_errorbar(aes(ymin = log10(mean_runtime - sd_runtime), ymax = log10(mean_runtime + sd_runtime)), width = 0.075) +
  labs(title = "Log runtime vs log number of imputations for cuPC-MI", x = "log10(m)", y = "log10(Mean Runtime (seconds))",
       color = "p", shape = "alpha") +
  theme_minimal()

# Plot: Runtime vs alpha for each p and m
ggplot(summary_results, aes(x = alpha, y = log10(mean_runtime), shape = factor(p), color = factor(m))) +
  geom_point(size = 3) +
  scale_shape_manual(values=c(15,16,17,18,4,3,7,8)) +
  geom_line(aes(group = interaction(p, m)), linetype = "solid") +
  #geom_errorbar(aes(ymin = mean_runtime - sd_runtime, ymax = mean_runtime + sd_runtime), width = 0.01) +
  labs(title = "Log runtime vs alpha for cu_pc_MI", x = "alpha", y = "log10(Mean Runtime (seconds))",
       color = "m", shape = "p") +
  theme_minimal()

# Plot: Runtime vs p for each m and alpha
ggplot(summary_results, aes(x = p, y = log10(mean_runtime), shape = factor(m), color = factor(alpha))) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(m, alpha)), linetype = "solid") +
  #geom_errorbar(aes(ymin = mean_runtime - sd_runtime, ymax = mean_runtime + sd_runtime), width = 1) +
  labs(title = "Runtime vs p for cu_pc_MI", x = "p", y = "Mean Runtime (seconds)",
       color = "m", shape = "alpha") +
  theme_minimal()
  
# Display a summary table
summary_results
