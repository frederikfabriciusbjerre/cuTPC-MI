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
source("R/cuTPCMI.R")
source("R/fixedMItest.R")

shdsum <- 0
hdsum <- 0
alpha <- 0.01
max_order <- 20
tiers_ <- rep(c(1,2,3,4), times = c(10,10,10,10))
# tiers_ <- NULL

for (i in 1:100){
    # read data as imputed_data
    dataset_path <- file.path(paste0("imputed_datasets/imp", i, "/imp", i, ".Rdata"),
                                     fsep = .Platform$file.sep)
    load.Rdata(dataset_path, "imputed_data")

    # make suffStat
    suffStatMI <- micd::getSuff(imputed_data, test="gaussMItest")
    suffStatSI <- list(C = suffStatMI[[1]], n = suffStatMI[[length(suffStatMI)]])
    # input params to pc
    p <- imputed_data[[1]] %>% length()

    # cat("Fitting with alpha =", alpha, "\n")
    cat("######################################################################################\n")
    cat("\n")
    cat("MI\n")
    tic()
    cuda_tPC_MI <- tpc(suffStatMI, indepTest = fixedGaussMItest, p = p, alpha = alpha, m.max = max_order, skel.method = "cuda", verbose = FALSE, tiers = tiers_)
    print(cuda_tPC@graph)
    cat("\n")
    cat("cuda_tPC ord =", cuda_tPC@max.ord, "\n")
    cat("alpha    =", alpha, "\n\n")
    cat("Time consumed:\n")
    toc()
    cat("\n")

    cat("\n")
    cat("SI\n")
    tic()
    cuda_tPC_SI <- tpc(suffStatSI, indepTest = fixedGaussMItest, p = p, alpha = alpha, m.max = max_order, skel.method = "cuda", verbose = FALSE, tiers = tiers_)
    print(cuda_tPC_SI@graph)
    cat("\n")
    cat("cuda_tPC_SI ord =", cuda_tPC_SI@max.ord, "\n")
    cat("alpha    =", alpha, "\n\n")
    cat("Time consumed:\n")
    toc()
    cat("\n")

    # mask <- near(micd_tPC@pMax, cuda_tPC@pMax, tol = 0.001)
    # lst1 <- micd_tPC@pMax[!mask]
    # lst2 <- cuda_tPC@pMax[!mask]
    # h <- 0
    # j <- 0
    # for (i in 1:length(lst1)){ # nolint
    #     if (lst2[i] != -Inf){
    #         # cat(lst1[i], lst2[i], "\n")
    #         h <- h + 1
    #     }
    #     j <- j + 1
    # }
    # cat("Number of pvals not equal             : ", j, "\n")
    # cat("Number of pvals not equal without -Inf: ", h, "\n")
    shdsum <- shdsum + shd(cuda_tPC_SI, cuda_tPC_MI)
    hdsum <- hdsum + shd(ugraph(cuda_tPC_SI@graph), ugraph(cuda_tPC_MI@graph))
    cat("shd:", shd(cuda_tPC_SI, cuda_tPC_MI), "hd:", shd(ugraph(cuda_tPC_SI@graph), ugraph(cuda_tPC_MI@graph)), "\n")
    cat("######################################################################################\n\n\n")
}
cat("shd:", shdsum, "\n")
cat("hd: ", hdsum, "\n")

