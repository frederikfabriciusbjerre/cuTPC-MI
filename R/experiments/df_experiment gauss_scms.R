library(mice)
library(micd)
library(pcalg)
library(graph)
library(tidyverse)
library(Metrics)

source("R/MeeksRules.R")
source("R/tcheckTriple.R")
source("R/tpc_cons_intern.R")
source("R/tskeleton.R")
source("R/tskeleton_parallel.R")
source("R/cuTPCMI.R")
source("R/fixedMItest.R")
source("cuda/cuPCMI.R")

zStatMI <- function (x, y, S, C, n)
{
  r <- pcalg::pcorOrder(x, y, S, C)
  res <- 0.5 * log_q1pm(r)
  if (is.na(res))
    0
  else res
}

log_q1pm <- function(r) log1p(2 * r / (1 - r))

df.reiter <- function(B, U, m, dfcom){
  # modified from https://github.com/mwheymans/miceafter/blob/main/R/pool_bftest.R

  t <- (m - 1)
  r <- (1+1/m)*B/U
  a <- r*t/(t-2)
  vstar <- ( (dfcom+1) / (dfcom+3) ) * dfcom
  
  c0 <- 1 / (t-4)
  c1 <- vstar - 2 * (1+a)
  c2 <- vstar - 4 * (1+a)
  
  z <- 1 / c2 +
    c0 * (a^2 * c1 / ((1+a)^2 * c2)) +
    c0 * (8*a^2 * c1 / ((1+a) * c2^2) + 4*a^2 / ((1+a) * c2)) +
    c0 * (4*a^2 / (c2 * c1) + 16*a^2 * c1 / c2^3) +
    c0 * (8*a^2 / c2^2)
  
  v <- 4 + 1/z
  
  return(v)
}

lambdas <- c()
dfs <- c()
# modified gaussMItest_df_corrected function
gaussMItest_df_corrected <- function (x, y, S, suffStat, df_correction_method = 'df_br') {
  
  # number of imputations
  M <- length(suffStat) - 1
  # sample size
  n <- suffStat[[M+1]]
  suffStat[[M+1]] <- NULL
  
  z <- sapply(suffStat, function(j) {
    zStatMI(x, y, S, C=j, n=n)
  })
  
  # 1. Average of M imputed data sets
  avgz <- mean(z)
  
  # 2. Average of completed-data variance
  W <- 1 / (n - length(S) - 3)
  
  # 3. Between variance
  B <- sum((z - avgz)^2) / (M - 1)
  
  # 4. Total variance
  TV <- W + (1 + 1 / M) * B
  
  # 5. Test statistic
  ts <- avgz / sqrt(TV)
  
  # 6. Degrees of freedom
  lambda <- (B + B/M) / TV
  lambdas <<- c(lambdas, c(lambda))
  df_old <- (M - 1) / lambda^2
  df_com <- n - (length(S) + 3)
  df_obs <- (1 - lambda) * ((df_com + 1) / (df_com + 3)) * df_com
  df_br <- df_old * df_obs / (df_old + df_obs)
  df_reiter <- df.reiter(B, W, M, df_com)
  
  # handle NaN values
  if (is.nan(df_br)){
    df_br <- ( (df_com+1) / (df_com+3) ) * df_com
  }
  
  # choose degrees of freedom based on correction method
  if (df_correction_method == 'df_old'){
    df_ <- df_old
  } else if (df_correction_method == 'df_br'){
    df_ <- df_br
  } else if (df_correction_method == 'df_reiter'){
    df_ <- df_reiter
  } else {
    stop("Invalid df_correction_method specified")
  }
  dfs <<- c(dfs, c(df_))
  # 7. pvalue
  pvalue <- 2 * stats::pt(abs(ts), df = df_, lower.tail = FALSE)
  
  return(pvalue)
}

i <- 0
alpha <- 0.05
for (epoch in 1:100){
for (p in c(10)) {
  for (prob_dag in c(0.01, 0.1, 0.25, 0.5)) {
    for (m in c(10, 100)) {
      for (n in c(25, 100, 1000)) {
        for (prop_miss in c(0.1, 0.5, 0.9)){
          cat("n =", n, ", m =", m, ", p =", p, ", prob_dag = ", prob_dag, "prop_miss = ", prop_miss, "epoch :", epoch, "\n")
          # Generate true DAG
          dag_true <- pcalg::randomDAG(p, prob = prob_dag) # true DAG
          cpdag_true <- pcalg::dag2cpdag(dag_true) # true CPDAG
          
          # simulate Gaussian data from the true DAG
          data <- pcalg::rmvDAG(n, dag_true, errDist = "normal", mix = 0.3)
          
          # ampute
          data_missing <- mice::ampute(data, prop = prop_miss, 
                                       mech = "MAR", 
                                       bycases = TRUE)$amp
          i <- i + 1
          set.seed(i)
          # multiple imputations
          tryCatch({
            imputed_data <- mice::mice(data_missing, m = m, method = "norm",
                                       visitSequence = "monotone", printFlag = FALSE)
          }, error = function(e) {
            message(paste("Error in iteration", i, ":", e$message))
          })

          suff_imputed <- micd::getSuff(imputed_data, test="gaussMItest")
          suff_complete <- list(C = cor(data), n = nrow(data))
          
          # run pc on all data
          pc_complete <- pcalg::pc(suffStat = suff_complete, indepTest = pcalg::gaussCItest, 
                                   alpha = alpha, p = p)
          adjmat_complete <- as(pc_complete@graph, "matrix")
          
          # run PC algorithm on imputed data with each df_correction_method
          for (df_correction_method in c('df_old', 'df_br', 'df_reiter')) {
            
            pc_imputed <- pcalg::pc(suffStat = suff_imputed, indepTest = function (x, y, S, suffStat) {
              gaussMItest_df_corrected(x, y, S, suffStat, df_correction_method)
            }, alpha = alpha, p = p)
            mean_lambda <- mean(lambdas)
            mean_df <- mean(dfs)
            lambdas <- c()
            dfs <- c()
            
            adjmat_imputed <- as(pc_imputed@graph, "matrix")

            # Flatten adjacency matrices
            adj_complete_vector <- as.vector((adjmat_complete + t(adjmat_complete)) > 0)
            adj_imputed_vector <- as.vector((adjmat_imputed + t(adjmat_imputed)) > 0)

            recall_ <- Metrics::recall(adj_complete_vector, adj_imputed_vector)
            precision_ <- Metrics::precision(adj_complete_vector, adj_imputed_vector)
            hamming_distance <- pcalg::shd(pc_complete@graph %>% ugraph(), pc_imputed@graph %>% ugraph())
            shd_value <- pcalg::shd(pc_complete@graph, pc_imputed@graph)
            
            # store results
            results <- rbind(results, data.frame(p = p,
                                                 prob_dag = prob_dag,
                                                 prop_miss = prop_miss,
                                                 m = m,
                                                 n = n,
                                                 epoch = epoch,
                                                 mean_lambda = mean_lambda,
                                                 df_correction_method = df_correction_method,
                                                 mean_df = mean_df,
                                                 recall = recall_,
                                                 precision = precision_,
                                                 hamming_distance = hamming_distance,
                                                 shd = shd_value,
                                                 stringsAsFactors = FALSE))
            # write.csv(results, "results.csv")
          }
        }
      }
    }
  }
  }
}
