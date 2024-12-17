library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
library(micd)
library(miceadds)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Metrics)
source("si_vs_mi_experiment/cohort_sim_utils.R")
source("R/MeeksRules.R")
source("R/tcheckTriple.R")
source("R/tpc_cons_intern.R")
source("R/tskeleton.R")
source("R/tskeleton_parallel.R")
source("R/tpc.R")
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
  df_old <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2
  df_com <- n - (length(S) + 3)
  df_obs <- (1 - lambda) * ((df_com + 1) / (df_com + 3)) * df_com
  df_br <- df_old * df_obs / (df_old + df_obs)
  df_reiter <- df.reiter(B, W, M, df_com)
  
  # Handle NaN values
  if (is.nan(df_br)){
    df_br <- ( (df_com+1) / (df_com+3) ) * df_com
  }
  
  # Choose degrees of freedom based on correction method
  if (df_correction_method == 'old'){
    df_ <- df_old
  } else if (df_correction_method == 'br'){
    df_ <- df_br
  } else if (df_correction_method == 'reiter'){
    df_ <- df_reiter
  } else {
    stop("Invalid df_correction_method specified")
  }
  dfs <<- c(dfs, c(df_))
  # 7. pvalue
  pvalue <- 2 * stats::pt(abs(ts), df = df_, lower.tail = FALSE)
  
  return(pvalue)
}

# Initialize a data frame to store results
results <- data.frame(p = integer(),
                      prob_dag = numeric(),
                      alpha = numeric(),
                      epoch = integer(),
                      mean_lambda = numeric(),
                      df_correction_method = character(),
                      mean_df = numeric(),
                      recall = numeric(),
                      precision = numeric(),
                      hamming_distance = integer(),
                      shd = integer(),
                      stringsAsFactors = FALSE)

filename <- "results_si_vs_mi_config2n3000_cpu.csv"

alpha <- 0.05
p <- 33
n <- 3000
m_max <- 1000
for (i in 1:99){  
  print(i)
  set.seed(i)
  imp_i <- paste0("imp", i)
  # load imputed data
  if (i < 1){
    dataset_imp_path<- paste0("si_vs_mi_experiment/imputed_datasets_cohort_config2n3000/", imp_i, "/", imp_i, ".Rdata")
  } else{
    dataset_imp_path<- paste0("si_vs_mi_experiment/imputed_datasets_cohort_config2n3000/", imp_i, ".Rdata")
  }
  imputed_data <- get(load(dataset_imp_path))
  
  # load complete data
  dataset_complete_path <- paste0("si_vs_mi_experiment/data/config2n3000/dataset", i, ".csv")
  complete_data <- read.csv(dataset_complete_path, header = FALSE)
  
  # get suff for all imputed datasets
  suffStatMI_full <- micd::getSuff(imputed_data, test="gaussMItest")
  
  # get suff for complete
  suff_complete <- list(C = cor(complete_data), n = nrow(complete_data))
  for (m in c(10, 100)) {
    #cat("m =", m, ", epoch :", epoch, "\n")
    # Get suffStat for imputed data
    if (m == 1000){
      print(m)
      suff_imputed <- suffStatMI_full
    }
    else if (m > 1){
      print(m)
      suff_imputed <- c(sample(head(suffStatMI_full, -1), m), c(n))
    }
    else {
      print(m)
      suff_imputed <- list(C = sample(head(suffStatMI_full, -1), 1)[[1]], 
                           n = suffStatMI_full[[length(suffStatMI_full)]])
    }
    
    # run PC-algorithm on complete data (without missing values)
    pc_complete <- pcalg::pc(suffStat = suff_complete, indepTest = pcalg::gaussCItest, alpha = alpha, p = p)
    adjmat_complete <- as(pc_complete@graph, "matrix")
    
    # run PC-algorithm on imputed data with each df_correction_method
    if (m > 1){
      for (df_correction_method in c('old')) {
        
        pc_imputed <- pcalg::pc(suffStat = suff_imputed, indepTest = function(x, y, S, suffStat) {
                                            gaussMItest_df_corrected(x, y, S, suffStat, df_correction_method)
                                        }, alpha = alpha, p = p)
        # pc_imputed <- cu_pc_MI(suffStat = suff_imputed, alpha = alpha, p = p, df_method = df_correction_method)
        mean_lambda <- mean(lambdas)
        mean_df <- mean(dfs)
        lambdas <- c()
        dfs <- c()
        
        adjmat_imputed <- as(pc_imputed@graph, "matrix")
        
        adj_complete_vector <- (adjmat_complete + t(adjmat_complete)) > 0 %>% as.numeric()
        adj_imputed_vector <- (adjmat_imputed + t(adjmat_imputed)) > 0 %>% as.numeric()
        
        recall_ <- Metrics::recall(adj_complete_vector, adj_imputed_vector)
        precision_ <- Metrics::precision(adj_complete_vector, adj_imputed_vector)
        hamming_distance <- pcalg::shd(pc_complete@graph %>% ugraph(), pc_imputed@graph %>% ugraph())
        shd_value <- pcalg::shd(pc_complete@graph, pc_imputed@graph)
      
        # store results
        results <- rbind(results, data.frame(p = p,
                                             m = m,
                                             n = n,
                                             mean_lambda = mean_lambda,
                                             df_correction_method = df_correction_method,
                                             mean_df = mean_df,
                                             recall = recall_,
                                             precision = precision_,
                                             hamming_distance = hamming_distance,
                                             shd = shd_value,
                                             stringsAsFactors = FALSE))
        write.csv(results, filename)
      }
    }
    # m = 1
    else {
      pc_imputed <- pcalg::pc(suffStat = suff_imputed, indepTest = pcalg::gaussCItest, alpha = alpha, p = p)
      mean_lambda <- NA # remember to comment on this in report
      mean_df <- NA
      lambdas <- c()
      dfs <- c()
      adjmat_imputed <- as(pc_imputed@graph, "matrix")

      adj_complete_vector <- (adjmat_complete + t(adjmat_complete)) > 0 %>% as.numeric()
      adj_imputed_vector <- (adjmat_imputed + t(adjmat_imputed)) > 0 %>% as.numeric()
      
      recall_ <- Metrics::recall(adj_complete_vector, adj_imputed_vector)
      precision_ <- Metrics::precision(adj_complete_vector, adj_imputed_vector)
      hamming_distance <- pcalg::shd(pc_complete@graph %>% ugraph(), pc_imputed@graph %>% ugraph())
      shd_value <- pcalg::shd(pc_complete@graph, pc_imputed@graph)
      
      # store results
      results <- rbind(results, data.frame(p = p,
                                           m = m,
                                           n = n,
                                           mean_lambda = mean_lambda,
                                           df_correction_method = mean_df,
                                           mean_df = mean_df,
                                           recall = recall_,
                                           precision = precision_,
                                           hamming_distance = hamming_distance,
                                           shd = shd_value,
                                           stringsAsFactors = FALSE))
      #write.csv(results, filename)
    }
  }
}
