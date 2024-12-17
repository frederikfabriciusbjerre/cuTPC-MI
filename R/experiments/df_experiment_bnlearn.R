library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(mice)
library(micd)
library(miceadds)
library(dplyr)
library(bnlearn)
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
dfs <- c()
lambdas <- c()

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
  df_old <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2
  df_com <- n - (length(S) + 3)
  df_obs <- (1 - lambda) * ((df_com + 1) / (df_com + 3)) * df_com
  df_br <- df_old * df_obs / (df_old + df_obs)
  df_reiter <- df.reiter(B, W, M, df_com)
  
  # handle NaN values
  if (is.nan(df_br) || is.infinite(df_br)){
    df_br <- ((df_com + 1) / (df_com + 3)) * df_com
  }
    if (is.nan(df_reiter) || is.infinite(df_reiter)){
    df_reiter <- ((df_com + 1) / (df_com + 3)) * df_com
  }
  
  # choose degrees of freedom based on correction method
  if (df_correction_method == 'old'){
    df_ <- df_old
  } else if (df_correction_method == 'br'){
    df_ <- df_br
  } else if (df_correction_method == 'reiter'){
    df_ <- df_reiter
  } else {
    stop("Invalid df_correction_method specified")
  }
  
  lambdas <<- c(lambdas, c(lambda))
  dfs <<- c(dfs, c(df_))
  # 7. pvalue
  pvalue <- 2 * stats::pt(abs(ts), df = df_, lower.tail = FALSE)
  
  return(pvalue)
}


# list of networks and their names
network_list <- list(
  ecoli =       readRDS("bnlearn_networks/ecoli70.rds"),
  magic_irri =  readRDS("bnlearn_networks/magic-irri.rds"),
  magic_niab =  readRDS("bnlearn_networks/magic-niab.rds"),
  arth150 =     readRDS("bnlearn_networks/arth150.rds")
)

results <- data.frame()

# Set alpha value for the PC algorithm
alpha <- 0.05

df_correction_methods <- c('old', 'br', 'reiter') 

i <- 0
for (network_name in names(network_list)) {
  network <- network_list[[network_name]]
  
  cat("Processing network:", network_name, "\n")
  
  p <- length(nodes(network))
  
  dag_graphNEL <- as.graphNEL(network)
  dag_adjmat <- as(dag_graphNEL, "matrix")
  
  for (m in c(100, 10)) {
    for (n in c(500, 100, 25)) {
        for (prop_miss in c(0.1, 0.5, 0.9)){
            for (epoch in 1:10) {
                i <- i + 1
                cat("Network:", network_name, ", n =", n, ", m =", m, ", prop_miss", prop_miss, ", epoch =", epoch, "\n")
                
                set.seed(i)  
                data_sample <- rbn(network, n = n)
                
                true_cpdag <- cpdag(network)
                
                cat("imputing... \n")
                # ampute 
                data_missing <- mice::ampute(data_sample, prop = prop_miss, 
                                            mech = "MAR", 
                                            bycases = TRUE)$amp
                
                tryCatch({
                    imputed_data <- mice::mice(data_missing, m = m, method = "norm",
                                        visitSequence = "monotone", printFlag = FALSE)
                }, error = function(e) {
                    message(paste("Error in iteration", i, ":", e$message))
                })
                
                suff_imputed <- micd::getSuff(imputed_data, test = "gaussMItest")
                
                cat("pcalg fitting... \n")
                # run PC-algorithm on complete data (without missing values)
                suff_complete <- list(C = cor(data_sample), n = nrow(data_sample))
                pc_complete <- pcalg::pc(suffStat = suff_complete, indepTest = pcalg::gaussCItest, 
                                 alpha = alpha, p = p)
                adjmat_complete <- as(pc_complete@graph, "matrix")
                
                # run PC-algorithm on imputed data with each df_correction_method
                for (df_correction_method in df_correction_methods) {
                    cat("pcalg fitting with", df_correction_method, "... \n")
                    pc_imputed <- pcalg::pc(suffStat = suff_imputed, indepTest = function(x, y, S, suffStat) {
                                            gaussMItest_df_corrected(x, y, S, suffStat, df_correction_method)
                                        }, alpha = alpha, p = p)
                    
                    mean_df <- mean(dfs, na.rm = TRUE)
                    mean_lambda <- mean(lambdas, na.rm = TRUE)
                    dfs <- c()
                    lambdas <- c()

                    adjmat_imputed <- as(pc_imputed@graph, "matrix")
                    adjmat_imputed_skeleton <- (adjmat_imputed + t(adjmat_imputed)) > 0
                    
                    adj_complete_vector <- as.vector((adjmat_complete + t(adjmat_complete)) > 0)
                    adj_imputed_vector <- as.vector((adjmat_imputed + t(adjmat_imputed)) > 0)
                    
                    recall_ <- Metrics::recall(adj_complete_vector, adj_imputed_vector)
                    precision_ <- Metrics::precision(adj_complete_vector, adj_imputed_vector)
                    shd_value <- pcalg::shd(pc_imputed@graph, pc_complete@graph)
                    hamming_distance <- pcalg::shd(ugraph(pc_imputed@graph), ugraph(pc_complete@graph))
                    
                    # store results
                    results <- rbind(results, data.frame(
                        network = network_name,
                        p = p,
                        m = m,
                        n = n,
                        epoch = epoch,
                        df_correction_method = df_correction_method,
                        mean_lambda = mean_lambda,
                        mean_df = mean_df,
                        prop_miss = prop_miss, 
                        recall = recall_,
                        precision = precision_,
                        hamming_distance = hamming_distance,
                        shd = shd_value,
                        stringsAsFactors = FALSE
                    ))
                    write.table(
                        results[nrow(results), ], 
                        file = "results_with_lambda.csv", 
                        sep = ",", 
                        row.names = FALSE, 
                        col.names = !file.exists("results_with_lambda.csv"),  # Write header only if file doesn't exist
                        append = TRUE
                    )
                    }
                }
            }
        }
    }
}
