# ================================= A test case =================================
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
source("R/tpc.R")
source("R/fixedMItest.R")

shdsum <- 0
hdsum <- 0
alpha <- 0.01
max_order <- 20
tiers_ <- rep(c(1,2,3,4), times = c(10,10,10,10))
i <- 1
# read data as imputed_data
dataset_path <- file.path(paste0("imputed_datasets/imp", i, "/imp", i, ".Rdata"),
                                    fsep = .Platform$file.sep)
load.Rdata(dataset_path, "imputed_data")

# make suffStat
suffStatMI <- micd::getSuff(imputed_data, test="gaussMItest")
# input params to pc
p <- imputed_data[[1]] %>% length()

# cat("Fitting with alpha =", alpha, "\n")
cat("######################################################################################\n")
cat("\n")
cat("cuda+cuda\n")
tic()
cuda_tPC_cuda <- tpc(suffStatMI, indepTest = fixedGaussMItest, p = p, alpha = alpha, m.max = max_order, skel.method = "cuda", verbose = FALSE, tiers = tiers_, tpc_cons_intern = "cuda")
print(cuda_tPC_cuda)
cat("\n")
cat("Time consumed:\n")
toc()
cat("\n")

cat("\n")
cat("cuda+standard\n")
tic()
cuda_tPC_standard <- tpc(suffStatMI, indepTest = fixedGaussMItest, p = p, alpha = alpha, m.max = max_order, skel.method = "cuda", verbose = FALSE, tiers = tiers_, tpc_cons_intern = "standard")
print(cuda_tPC_standard)
cat("\n")
cat("Time consumed:\n")
toc()
cat("\n")

cat("shd:", shd(cuda_tPC_standard, cuda_tPC_cuda), "hd:", shd(ugraph(cuda_tPC_standard@graph), ugraph(cuda_tPC_cuda@graph)), "\n")
cat("######################################################################################\n\n\n")