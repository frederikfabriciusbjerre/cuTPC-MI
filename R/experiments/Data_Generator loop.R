#================================= Create DataSet =================================
library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
library(miceadds)

p <- 40 #number of nodes
probability <- 0.25
n <- 1000 #number of sample
vars <- c(paste0(1:p))


# mice params
prob_miss <- 0.1
m <- 100
method <- "norm"
n_datasets <- 100
# system("rm -rf imputed_datasets")
# system("mkdir imputed_datasets")


for (i in 1:n_datasets) {
    tic()
    set.seed(i)
    # Generate a random DAG and synthetic data
    gGtrue <- randomDAG(p, prob = probability, lB = 0.1, uB = 1, V = vars)
    N1 <- runif(p, 0.5, 1.0)
    Sigma1 <- matrix(0, p, p)
    diag(Sigma1) <- N1
    eMat <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma1)
    gmG <- list(x = rmvDAG(n, gGtrue, errMat = eMat), g = gGtrue)
    dataset_filename <- paste0("data/dataset", i, ".csv")
    write.table(gmG$x, file = dataset_filename, row.names = FALSE, na = "", col.names = FALSE, sep = ",")

    dataset_path <- file.path(dataset_filename, fsep = .Platform$file.sep)
    data <- read.table(dataset_path, sep = ",")

    data_missing <- ampute(data, prop = prob_miss, mech = "MAR", bycases = TRUE)$amp
    imputed_data <- mice(data_missing, m = m, method = method, printFlag = FALSE, remove.collinear = TRUE)
imputed_datasets
    setwd("/imputed_datasets")
    # write mids obj
    write.mice.imputation(imputed_data, paste0("imp", i), dattype = "csv", mids2spss = FALSE)
    setwd("../")
}
