# Load packages needed
library(RBGL)
library(pcalg)
library(bnlearn)

# Define functions (with better argument titles than earlier)

# in_par: Probability of inserting an edge between two nodes in the same tier
# cr_par: Probability of inserting a cross-tier edge between two nodes in different tiers
# tino: Vector containing number of nodes in each tier. 
# i.e. c(4,3) has 2 tiers, 4 nodes in tier 1 and 3 nodes in tier 2. 
# lB and uB. Lower and upper bound for entries in the weighted adjacency matrix. Chosen uniformly
# NoEmptyDAG: If TRUE, will NOT return a empty DAG
# Output: weighted adjacency matrix
randtierDAG <- function(in_par, cr_par, tino, lB = 0, uB = 1, NoEmptyDAG = T){
  numno <- sum(tino)
  adjmat <- matrix(0,nrow = numno,ncol = numno)
  #tier order
  tior <- c()
  for (o in 1:length(tino)){
    tior <- c(tior,rep(o,tino[o]))
  }
  
  emptyDAG <- T
  while (emptyDAG) {
    #In tier DAGs
    for (i in 1:length(tino)){
      if (tino[i] > 1){
        tempDAG <- r.gauss.pardag(p = tino[i], prob = in_par, lbe = lB, ube = uB, neg.coef = F)
        adjmat[which(tior == i),which(tior == i)] <- t(tempDAG$weight.mat())
      }
    }
    # Make across tier edges
    for (n in 1:(length(tino)-1)){
      plustier <- which(tior > n)
      tier <- which(tior == n)
      adjmat[plustier,tier] <- rbinom(n = length(plustier)*length(tier), size = 1,prob = cr_par)*runif(length(plustier)*length(tier), min = lB, max = uB)
    }
    emptyDAG <- ifelse(NoEmptyDAG,(sum(adjmat) == 0),F)
    
  }
  
  return(adjmat)
  
}

# n: number of observations
# w_adjmat: weighted adjacency matrix. 
# Output: dataframe with n rows, and columns according to w_adjmat. Each column names as an integer
rmvTDAG <- function(n, w_adjmat){
  graphnelDAG <- as(t(w_adjmat), "graphNEL")
  graphnelDAG <- subGraph(tsort(graphnelDAG),graphnelDAG)
  data <- rmvDAG(n,graphnelDAG, errDist = "normal")
  data <- data[,as.character(sort(as.integer(colnames(data))))]
  return(data)
}

# Here is an extra function that converts an adjacency matrix of a DAG to the 
# adjacency matrix of the corresponding Tiered MPDAG
# Function: Convert DAG to Tiered MPDAG
# adjmatrix: adjacency matrix (entries are 0 and 1's)
# odr: Vector containing indices of tiers for each node. SORRY: THIS IS ANOTHER TYPE THAN FOR randtierDAG()
#i.e. c(1,1,1,1,2,2,2) has 2 tiers, 4 nodes in tier 1 and 3 nodes in tier 2. 
# Output: Adjacency matrix
dag2tmpdag <- function(adjmatrix,odr){
  cpdag <- dag2cpdag(as(t(adjmatrix), "graphNEL"))
  cpadjmat <- t(as(cpdag,"matrix"))
  number_nodes <- length(odr)
  node_numbers <- 1:number_nodes
  Forbidden.edges.matrix <- matrix(0,ncol = number_nodes, nrow = number_nodes)  #list of nodes all with integer(0) entry
  for (nod in node_numbers){
    Forbidden.edges.matrix[nod,] <- odr[nod]<odr
  }
  tmadjmat <- addBgKnowledge(cpadjmat*!Forbidden.edges.matrix, checkInput = F)
  return(tmadjmat)
}

simulate_cohort <- function(n, p, num_tiers, in_tier_prob, cross_tier_prob, lB=0, uB=1){
  in_par <- in_tier_prob
  ac_par <- cross_tier_prob
  
  num_no <- p #Number of node
  num_ti <- num_tiers #Number of tiers
  if (num_no < num_ti){
    stop("Number of nodes is less than number of tiers")
  }
  # Makes a vector with integer representing tier
  tiered_ordering <- sort(c
                          (sample(c(1:num_ti), size = (num_no-num_ti), replace = T),
                            1:num_ti)
                          )
  
  # Converts to type neccessary for randtierDAG
  tino <- c()
  for (o in 1:num_ti){
    tino <- c(tino,sum(tiered_ordering == o))
  }
  
  # Simulate DAG from chosen tiers, nodes and density parameters
  weightadjmat <- randtierDAG(in_par,ac_par,tino, lB = lB, uB = uB)
  
  # Generate data from DAG
  data <- rmvTDAG(n, weightadjmat)
  
  # Convert weighted adjeceny matrix to adjacency matrix (Only 0,1 entries)
  adjmat <- ceiling(weightadjmat/uB)
  DAG <- as(t(adjmat), "graphNEL")
 
  # convert to tMPDAG
  t_adjmat <- dag2tmpdag(adjmatrix = adjmat,odr = tiered_ordering)
  T_MPDAG <- as(t(t_adjmat), "graphNEL")
  return (list(data = as.data.frame(data), tier_ord = tiered_ordering, dag = DAG, tmpdag = T_MPDAG))
}

adj_matrix2ampute_matrix <- function(adj_mat, n_amputation_dependencies) {
  # replace n_amputation_dependencies neighbors with weights, set rest to 0
  for (i in 1:nrow(adj_mat)) {
    row <- adj_mat[i, ]
    ones_indices <- which(row == 1)
    
    # if no. neighbors > n_amputation_dependencies, set extra 1s to 0
    if (length(ones_indices) > n_amputation_dependencies) {
      indices_to_zero <- sample(ones_indices, length(ones_indices) - n_amputation_dependencies)
      row[indices_to_zero] <- 0
    }
    
    # replace remaining ones with random numbers between 0 and 1
    ones_indices <- which(row == 1)
    row[ones_indices] <- runif(length(ones_indices))
    
    # normalize
    row_sum <- sum(row)
    if (row_sum > 0) {
      row <- row / row_sum
    }
    
    # update matrix
    adj_mat[i, ] <- row
  }
  
  return(adj_mat)
}

ampute_cohort <- function(data, 
                          dag, 
                          tier_order, 
                          ampute_proportion_total, 
                          ampute_proportion_tier, 
                          ampute_proportion_dropout, 
                          n_dropout_tiers_to_be_excluded = 1,
                          n_amputation_dependencies = 2) {
  n_rows <- nrow(data)
  n_cols <- ncol(data)
  data_amputed <- data
  
  # we do not ampute first or last tier in tier amputation
  tiers <- sort(unique(tier_order))
  tiers_excl_first_last <- tiers[tiers != min(tiers) & tiers != max(tiers)]
  
  ## 1. Tier Amputation ##
  
  num_rows_tier <- ceiling(ampute_proportion_tier * n_rows)
  if (num_rows_tier > 0) {
    rows_tier <- sample(1:n_rows, num_rows_tier, replace = FALSE)
    selected_tiers <- sample(tiers_excl_first_last, num_rows_tier, replace = TRUE)
    for (i in seq_along(rows_tier)) {
      row_idx <- rows_tier[i]
      tier <- selected_tiers[i]
      vars_in_tier <- which(tier_order == tier)
      
      # ampute
      data_amputed[row_idx, vars_in_tier] <- NA
    }
  } else {
    rows_tier <- integer(0)
  }
  
  ## 2. Dropout Amputation ##
  
  # exclude the first n_dropout tiers
  tiers_dropout <- tiers[-(1:n_dropout_tiers_to_be_excluded)]
  
  # exclude rows already selected in tier amputation
  available_rows <- setdiff(1:n_rows, rows_tier)
  num_rows_dropout <- ceiling(ampute_proportion_dropout * n_rows)
  if (num_rows_dropout > 0 && length(available_rows) >= num_rows_dropout) {
    rows_dropout <- sample(available_rows, num_rows_dropout, replace = FALSE)
    selected_tiers_dropout <- sample(tiers_dropout, num_rows_dropout, replace = TRUE)
    for (i in seq_along(rows_dropout)) {
      row_idx <- rows_dropout[i]
      tier <- selected_tiers_dropout[i]
      
      # identify tiers from the selected tier onwards
      tiers_to_ampute <- tiers[tiers >= tier]
      vars_to_ampute <- which(tier_order %in% tiers_to_ampute)
      
      # ampute
      data_amputed[row_idx, vars_to_ampute] <- NA
    }
  } else {
    rows_dropout <- integer(0)
  }
  
  ## 3. MAR Amputation using mice::ampute ##
  
  # get available rows / rows without amputations already
  rows_already_amputed <- unique(c(rows_tier, rows_dropout))
  available_rows_mar <- setdiff(1:n_rows, rows_already_amputed)
  
  # pull out completely observed cases
  data_not_amputed <- data_amputed[available_rows_mar, ]
  
  # calc proportion of complete data observations
  complete_case_proportion <- 1 - (ampute_proportion_tier + ampute_proportion_dropout)
  
  # since we ampute on a subset, we need to increase proportion according to 
  # proportion of missingness occurring due other types
  
  ampute_proportion_mar <- ampute_proportion_total - (ampute_proportion_tier + ampute_proportion_dropout)
  
  if (ampute_proportion_mar < 0) {
    stop("ampute_proportion_total should be larger than the sum of tier and dropout")
  }
  
  ampute_proportion_mar <- ampute_proportion_mar / complete_case_proportion

  # create the moral graph
  dag_bn <- as.bn(dag)
  moral_graph <- moral(dag_bn)
  
  # this will be used as markov blankets
  moral_adj_matrix <- amat(moral_graph)
  
  # convert moral_adj_matrix to missingness weights for MAR amputation
  weights <- adj_matrix2ampute_matrix(moral_adj_matrix, n_amputation_dependencies)
  # ampute using mice with given weight matrix
  ampute_result <- ampute(
    data = data_not_amputed,
    prop = ampute_proportion_mar, 
    weights = weights,
    mech = "MAR",
    bycases = TRUE
  )
  
  # replace in original data
  data_amputed[available_rows_mar, ] <- ampute_result$amp
  
  # mice gotta have these colnames for amputed datasets :)
  colnames(data_amputed) <- colnames(ampute_result$amp)
  
  return(data_amputed)
}

simulate_ampute_cohort <- function(n, 
                                   p, 
                                   num_tiers, 
                                   in_tier_prob, 
                                   cross_tier_prob, 
                                   lB, 
                                   uB,
                                   ampute_proportion_total,
                                   ampute_proportion_tier,
                                   ampute_proportion_dropout,
                                   n_dropout_tiers_to_be_excluded,
                                   n_amputation_dependencies){
  cohort <- simulate_cohort(
    n = n, 
    p = p,
    num_tiers = num_tiers,
    in_tier_prob = in_tier_prob,
    cross_tier_prob = cross_tier_prob, 
    lB = lB, 
    uB = lB
  )
  
  data <- cohort$data
  tier_ord <- cohort$tier_ord
  dag <- cohort$dag
  tmpdag <- cohort$tmpdag
  
  amputed_data <- ampute_cohort(
    data = data,
    dag = dag,
    tier_order = tier_ord,
    ampute_proportion_total = ampute_proportion_total,
    ampute_proportion_tier = ampute_proportion_tier,
    ampute_proportion_dropout = ampute_proportion_dropout,
    n_dropout_tiers_to_be_excluded = n_dropout_tiers_to_be_excluded,
    n_amputation_dependencies = n_amputation_dependencies
  )
  return ((list(data = data, 
                amputed_data = amputed_data, 
                tier_ord = tier_ord, 
                dag = dag, 
                tmpdag = tmpdag)))
}