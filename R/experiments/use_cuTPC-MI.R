# ================================= A test case =================================
library(pcalg)
library(graph)
library(tictoc)
library(mice)
library(micd)
library(dplyr)
source("R/MeeksRules.R")
source("R/tcheckTriple.R")
source("R/tpc_cons_intern.R")
source("R/tskeleton.R")
source("R/tskeleton_parallel.R")
source("R/cuTPCMI.R")
source("R/fixedMItest.R")
source("R/experiments/si_vs_mi_experiment/cohort_sim_utils.R")

set.seed(1337)

alpha <- 0.01
max_order <- 20
m <- 10

# generate cohort data
cat("Generating amputed cohort data...\n")

cohort <- simulate_ampute_cohort(
    n = 1000,
    p = 20,
    num_tiers = 3,
    in_tier_prob = 0.1,
    cross_tier_prob = 0.15,
    lB = 0.1,
    uB = 1,
    ampute_proportion_total = 0.9,
    ampute_proportion_tier = 0.1,
    ampute_proportion_dropout = 0.1,
    n_dropout_tiers_to_be_excluded = 1,
    n_amputation_dependencies = 3
)

data <- cohort$data
amputed_data <- cohort$amputed_data
tiers <- cohort$tier_ord
true_dag <- cohort$dag
true_tmpdag <- cohort$tmpdag

# impute data
cat("Generating multiple imputations...\n")

imputed_data <- mice(amputed_data, m = m, method = "norm", printFlag = FALSE, remove.collinear = TRUE)

# make suffStat
cat("Making sufficient statistics...\n")
suffStatMI <- micd::getSuff(imputed_data, test="gaussMItest")
suffStatSI <- list(C = suffStatMI[[1]], n = suffStatMI[[length(suffStatMI)]])

p <- imputed_data[[1]] %>% length()

# fit 
cat("Fitting cuTPC-MI...\n")
cat("\n")
tic()
cuda_tPC_MI <- tpc(suffStatMI, indepTest = fixedGaussMItest, p = p, alpha = alpha, m.max = max_order, skel.method = "cuda", verbose = FALSE, tiers = tiers)
print(cuda_tPC_MI)
cat("\n")
cat("Time consumed:\n")
toc()
cat("\n")

# fit with m=1
cat("Fitting cuTPC-SI...\n")
cat("\n")
tic()
cuda_tPC_SI <- tpc(suffStatSI, indepTest = fixedGaussMItest, p = p, alpha = alpha, m.max = max_order, skel.method = "cuda", verbose = FALSE, tiers = tiers)
print(cuda_tPC_SI)
cat("\n")
cat("Time consumed:\n")
toc()
cat("\n")
