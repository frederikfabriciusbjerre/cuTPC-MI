#' PC Algorithm Accounting for a Partial Node Ordering
#'
#' Like [pcalg::pc()], but takes into account a user-specified partial
#' ordering of the nodes/variables. This has two effects:
#' 1) The conditional independence between \code{x} and \code{y} given \code{S} is
#' ot tested if any variable in \code{S} lies in the future of both \code{x} and \code{y};
#' 2) edges cannot be oriented from a higher-order to a lower-order node. In addition,
#' the user may specify individual forbidden edges and context variables.
#'
#' @param suffStat A [base::list()] of sufficient statistics, containing all necessary
#' elements for the conditional independence decisions in the function [indepTest()].
#' @param indepTest A function for testing conditional independence. It is internally
#' called as \code{indepTest(x,y,S,suffStat)}, and tests conditional independence of
#' \code{x} and \code{y} given \code{S}. Here, \code{x} and \code{y} are variables, and
#' \code{S} is a (possibly empty) vector of variables (all variables are denoted by their
#' (integer) column positions in the adjacency matrix). \code{suffStat} is a list,
#' see the argument above. The return value of \code{indepTest} is the p-value of the
#' test for conditional independence.
#' @param alpha significance level (number in \emph{(0,1)} for the individual conditional
#' independence tests.
#' @param labels (optional) character vector of variable (or "node") names.
#' Typically preferred to specifying \code{p}.
#' @param p (optional) number of variables (or nodes). May be specified if \code{labels}
#' are not, in which case \code{labels} is set to \code{1:p}.
#' @param skel.method Character string specifying method; the default, "stable" provides
#' an order-independent skeleton, see [tpc::tskeleton()].
#' @param forbEdges A logical matrix of dimension p*p. If \code{[i,j]} is TRUE, then the
#' directed edge i->j is forbidden. If both \code{[i,j]} and \code{[j,i]} are TRUE, then any type of
#' edge between i and j is forbidden.
#' @param m.max Maximal size of the conditioning sets that are considered in the
#' conditional independence tests.
#' @param conservative Logical indicating if conservative PC should be used.
#' Defaults to FALSE. See [pcalg::pc()] for details.
#' @param maj.rule Logical indicating if the majority rule should be used. Defaults
#' to TRUE. See [pcalg::pc()] for details.
#' @param tiers Numeric vector specifying the tier / time point for each variable.
#' Must be of length 'p', if specified, or have the same length as 'labels', if specified.
#' A smaller number corresponds to an earlier tier / time point.
#' @param context.all Numeric or character vector. Specifies the positions or names
#' of global context variables. Global context variables have no incoming edges, i.e.
#' no parents, and are themselves parents of all non-context variables in the graph.
#' @param context.tier Numeric or character vector. Specifies the positions or
#' names of tier-specific context variables. Tier-specific context variables have no
#' incoming edges, i.e. no parents, and are themselves parents of all non-context variables
#' in the same tier.
#' @param verbose if \code{TRUE}, detailed output is provided.
#' @param numCores The numbers of CPU cores to be used.
#' @param cl.type The cluster type. Default value is \code{"PSOCK"}.
#' For High-performance clusters use \code{"MPI"}. See also \code{parallel::\link[parallel]{makeCluster}}.
#' @param clusterexport Character vector. Lists functions to be exported to nodes if numCores > 1.
#'
#' @details See \code{pcalg::\link[pcalg]{pc}} for further information on the PC algorithm.
#' The PC algorithm is named after its developers Peter Spirtes and Clark Glymour
#' (Spirtes et al., 2000).
#'
#' Specifying a tier for each variable using the \code{tier} argument has the
#' following effects:
#' 1) In the skeleton phase and v-structure learing phases,
#' conditional independence testing is restricted such that if x is in tier t(x)
#' and y is in t(y), only those variables are allowed in the conditioning set whose
#' tier is not larger than t(x).
#' 2) Following the v-structure phase, all
#' edges that were found between two tiers are directed into the direction of the
#' higher-order tier. If context variables are specified using \code{context.all}
#' and/or \code{context.tier}, the corresponding orientations are added in this step.
#'
#' @return An object of \code{\link[base]{class}} "\code{pcAlgo}"
#' (see [pcalg::pcalgo] containing an estimate of the equivalence class of
#' the underlying DAG.
#'
#' @author Original code by Markus Kalisch, Martin Maechler, and Diego Colombo.
#' Modifications by Janine Witte (Kalisch et al., 2012).
#'
#' @references   M. Kalisch, M. Maechler, D. Colombo, M.H. Maathuis and P. Buehlmann (2012).
#' Causal Inference Using Graphical Models with the R Package pcalg.
#' Journal of Statistical Software 47(11): 1--26.
#'
#' P. Spirtes, C. Glymour and R. Scheines (2000). Causation, Prediction,
#' and Search, 2nd edition. The MIT Press. https://philarchive.org/archive/SPICPA-2.
#'
#' @export
#'
#' @examples
#' # load simulated cohort data
#' data(dat_sim)
#' n <- nrow(dat_sim)
#' lab <- colnames(dat_sim)
#'
#' # estimate skeleton without taking background information into account
#' tpc.fit <- tpc(suffStat = list(C = cor(dat_sim), n = n),
#'                indepTest = gaussCItest, alpha = 0.01, labels = lab)
#' pc.fit <- pcalg::pc(suffStat = list(C = cor(dat_sim), n = n),
#'                     indepTest = gaussCItest, alpha = 0.01, labels = lab,
#'                     maj.rule = TRUE, solve.conf = TRUE)
#' identical(pc.fit@graph, tpc.fit@graph) # TRUE
#' # estimate skeleton with temporal ordering as background information
#' tiers <- rep(c(1,2,3), times=c(3,3,3))
#' tpc.fit2 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
#'                 indepTest = gaussCItest, alpha = 0.01, labels = lab, tiers = tiers)
#'
#' tpc.fit3 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
#'                 indepTest = gaussCItest, alpha = 0.01, labels = lab, tiers = tiers,
#'                 skel.method = "stable.parallel",
#'                 numCores = 2, clusterexport = c("cor", "ecdf"))
#'
#' if(requireNamespace("Rgraphviz", quietly = TRUE)){
#'  data("true_sim")
#'  oldpar <- par(mfrow = c(1,3))
#'  plot(true_sim, main = "True DAG")
#'  plot(tpc.fit, main = "PC estimate")
#'  plot(tpc.fit2, main = "tPC estimate")
#'  par(oldpar)
#'  }
#'
#'  # require that there is no edge between A1 and A1, and that any edge between A2 and B2
#'  # or A2 and C2 is directed away from A2
#'  forb <- matrix(FALSE, nrow=9, ncol=9)
#'  rownames(forb) <- colnames(forb) <- lab
#'  forb["A1","A3"] <- forb["A3","A1"] <- TRUE
#'  forb["B2","A2"] <- TRUE
#'  forb["C2","A2"] <- TRUE
#'
#'  tpc.fit3 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
#'                  indepTest = gaussCItest, alpha = 0.01,labels = lab,
#'                  forbEdges = forb, tiers = tiers)
#'
#'  if (requireNamespace("Rgraphviz", quietly = TRUE)) {
#'  # compare estimated CPDAGs
#'    data("true_sim")
#'    oldpar <- par(mfrow = c(1,2))
#'    plot(tpc.fit2, main = "old tPC estimate")
#'    plot(tpc.fit3, main = "new tPC estimate")
#'    par(oldpar)
#'  }
#'  # force edge from A1 to all other nodes measured at time 1
#'  # into the graph (note that the edge from A1 to A2 is then
#'  # forbidden)
#'  tpc.fit4 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
#'                  indepTest = gaussCItest, alpha = 0.01, labels = lab,
#'                  tiers = tiers, context.tier = "A1")
#'
#'  if (requireNamespace("Rgraphviz", quietly = TRUE)) {
#'  # compare estimated CPDAGs
#'   data("true_sim")
#'   plot(tpc.fit4, main = "alternative tPC estimate")
#'  }
#'
#'  # force edge from A1 to all other nodes into the graph
#'  tpc.fit5 <- tpc(suffStat = list(C = cor(dat_sim), n = n),
#'                  indepTest = gaussCItest, alpha = 0.01, labels = lab,
#'                  tiers = tiers, context.all = "A1")
#'
#'  if (requireNamespace("Rgraphviz", quietly = TRUE)) {
#'  # compare estimated CPDAGs
#'  data("true_sim")
#'  plot(tpc.fit5, main = "alternative tPC estimate")
#'  }
#'
tpc <- function (suffStat, indepTest, alpha, labels, p,
                 skel.method = c("stable", "stable.parallel", "cuda"),
                 forbEdges = NULL, m.max = Inf,
                 conservative = FALSE, maj.rule = TRUE,
                 tiers = NULL, context.all = NULL, context.tier = NULL,
                 verbose = FALSE,
                 numCores = NULL, cl.type = "PSOCK",
                 clusterexport = NULL,
                 tpc_cons_intern = "standard"){
  cl <- match.call()
  if (!missing(p)) {
    stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
  }
  if (missing(labels)) {
    if (missing(p))
      stop("need to specify 'labels' or 'p'")
    labels <- as.character(seq_len(p))
  } else {
    stopifnot(is.character(labels))
    if (missing(p)) {
      p <- length(labels)
    } else if (p != length(labels))
      stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
    else message("No need to specify 'p', when 'labels' is given")
  }

  if (is.null(tiers)) {
    ## if no tiers are specified, everything is tier 1
    tiers <- rep(1, p)
  } else {
    ## check if 'tiers' are correctly specified
    if (!is.numeric(tiers)) {stop("'tiers' must be a numeric vector")}
    if (length(tiers) != p) {stop("length of 'tiers' does not match 'p' or length of 'labels'")}
  }

  if (!is.null(context.all)) {
    if (is.character(context.all)) {
      if (!all(context.all %in% labels)) {stop("'context.all' includes variable names not in 'labels'")}
      context.all <- which(labels %in% context.all)
    }

    if (is.numeric(context.all)) {
      if (!all(context.all %in% (1:p))) {stop("'context.all' contains elements that are smaller than 1 or larger than 'p'")}
      if (!all(tiers[context.all]==min(tiers))) {stop("'context.all' variables must be in the first tier")}
    } else {
      stop("'context.all' must be an integer vector or character vector")
    }
  }

  if (!is.null(context.tier)) {
    if (is.character(context.tier)) {
      if (!all(context.tier %in% labels)) {stop("'context.tier' includes variable names not in 'labels'")}
      context.tier <- which(labels %in% context.tier)
    }
    if (is.numeric(context.tier)) {
      if (!all(context.tier %in% 1:p)) {stop("'context.tier' contains elements that are smaller than 1 or larger than 'p'")}
    } else {
      stop("'context.tier' must be a numeric or character vector")
    }
  }

  if ( !is.null(context.tier) && !is.null(context.all) ) {
    if (length(intersect(context.tier, context.all)) > 0) {
      stop(paste("The following variables are in both 'context.tier' and 'context.all': ",
                 paste(intersect(context.tier, context.all), collapse=",")))
    }
  }

  if (conservative && maj.rule) {
    stop("Choose either conservative PC or majority rule PC!")
  }
  if ((!conservative) && (!maj.rule)) {
    stop("Choose one of conservative PC and majority rule PC!")
  }

  if(!is.null(numCores) && is.null(cl.type)) stop("Specify cluster type.")


  ## generate fixedEdges and fixedGaps according to context.all and context.tier
  fixedEdges <- matrix(FALSE, p, p)
  fixedGaps <- matrix(FALSE, p, p)

  ## context.all
  if (!is.null(context.all)) {
    for (i in context.all) {
      fixedEdges[i, ] <- TRUE
      fixedEdges[ ,i] <- TRUE
      for (j in context.all) {
        fixedEdges[i,j] <- FALSE
        fixedEdges[j,i] <- FALSE
        fixedGaps[i,j] <- TRUE
        fixedGaps[j,i] <- TRUE
      }
    }
  }
  ## context.tier
  if (!is.null(context.tier)) {
    for (i in context.tier) {
      for (j in c(context.tier, context.all)) {
        fixedGaps[i,j] <- TRUE
        fixedGaps[j,i] <- TRUE
      }
      k <- tiers[i]
      fixedEdges[i,tiers==k] <- TRUE
      fixedEdges[tiers==k,i] <- TRUE
      fixedGaps[i,tiers!=k] <- TRUE
      fixedGaps[tiers!=k,i] <- TRUE
    }
  }

  if (!is.null(forbEdges)) {

    ## check if forbEdges contradicts context.tier or context.all
    checkMatrix <- fixedEdges
    for (i in context.all) {
      checkMatrix[ ,i] <- FALSE
    }
    for (i in context.tier) {
      k <- tiers[i]
      checkMatrix[tiers==k,i] <- FALSE
    }

    if ( sum(checkMatrix*forbEdges) > 0 ) {
      ConflictList <- which((checkMatrix*forbEdges) > 0, arr.ind=TRUE)
      ConflictList[ ,1] <- labels[ConflictList[ ,1]]
      ConflictList[ ,2] <- labels[as.numeric(ConflictList[ ,2])]
      colnames(ConflictList) <- NULL
      cat("Note: 'forbEdges' overrules 'context.tier' and 'context.all'.
          Edges between the following pairs of nodes are forbidden by
          'forbEdges' even though 'context.tier' and/or 'context.all'
          suggest they should be present:\n")
      print(ConflictList)
    }

    ## matrix of forbidden adjacencies (no type of edge allowed between a and
    ## b):
    forbAdj <- forbEdges * t(forbEdges)

    ## modify fixedEdges and fixedGaps according to forbEdges
    fixedGaps <- ( fixedGaps + forbAdj ) > 0
    fixedEdges <- ( fixedEdges - forbAdj ) > 0

    ## generate list of forbidden arrows (where the other direction is allowed)
    forbArrows <- forbEdges - forbEdges * t(forbEdges)
    forbArrowsL <- which(forbArrows > 0, arr.ind=TRUE)
    forbArrowList <- lapply(seq_len(nrow(forbArrowsL)),
                            function(i) forbArrowsL[i,])
  } else {
    forbArrowList <- list()
  }

  # skel.method <- "stable"
  skel.method <- match.arg(skel.method)

  if (skel.method == "cuda"){
    skel <- tskeleton_cuda_MI(suffStat, indepTest, alpha, labels = labels,
                              method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                              m.max = m.max, verbose = verbose, tiers = tiers)
  }
  else if (skel.method == "stable.parallel") {
    if (is.null(numCores)) {stop("Please specify 'numCores'.")}
    skel <- tskeleton_parallel(suffStat, indepTest, alpha, labels = labels,
                               method = skel.method,
                               fixedGaps = fixedGaps,
                               fixedEdges = fixedEdges,
                               tiers = tiers,
                               m.max = m.max, verbose = verbose,
                               numCores = numCores,
                               cl.type = cl.type,
                               clusterexport = clusterexport)
  } else {
    skel <- tskeleton(suffStat, indepTest, alpha, labels = labels,
                      method = skel.method, fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                      m.max = m.max, verbose = verbose, tiers = tiers)
  }

  skel@call <- cl

  if (numEdges(skel@graph) == 0) {
    return(skel)
  }

  ## step II, orientation of v-structures:

  # Prepare arguments for the .C call to tpcConsIntern
  if (tpc_cons_intern == "cuda"){

    p <- length(labels)
    M <- length(suffStat) - 1
    nrows <- suffStat[[length(suffStat)]]
    cat(p, M, nrows)
    C_list <- head(suffStat, M)
    C_array <- array(0, dim = c(p, p, M))
    print(dim(C_array))    # Should be (90, 90, M) for p = 90
    print(length(C_list)) 
    for (i in 1:M) {
      C_array[, , i] <- C_list[[i]]
    }
    C_array[is.na(C_array)] <- 0.0

    G_matrix <- as(skel@graph, "matrix")
    G_host <- as.integer(G_matrix)

    tiers_host <- as.integer(tiers)
    context_all_host <- if (!is.null(context.all)) as.integer(context.all - 1) else integer(0)
    n_context_all <- length(context_all_host)
    context_tier_host <- if (!is.null(context.tier)) as.integer(context.tier - 1) else integer(0)
    n_context_tier <- length(context_tier_host)

    # Prepare SepSet
    SepSet_flat_host <- integer(0)
    SepSet_offsets_host <- integer(p * p)
    SepSet_lengths_host <- integer(p * p)
    offset <- 0
    for (a in 1:p) {
      for (c in 1:p) {
        idx <- (a - 1) * p + (c - 1)
        sepset_ac <- skel@sepset[[a]][[c]]
        if (!is.null(sepset_ac)) {
          SepSet_offsets_host[idx + 1] <- offset
          SepSet_lengths_host[idx + 1] <- length(sepset_ac)
          SepSet_flat_host <- c(SepSet_flat_host, sepset_ac)
          offset <- offset + length(sepset_ac)
        } else {
          SepSet_offsets_host[idx + 1] <- -1
          SepSet_lengths_host[idx + 1] <- 0
        }
      }
    }
    SepSet_flat_length <- length(SepSet_flat_host)
    max_unfTriple <- 1000  # Adjust as needed
    ambiguous_triples_host <- integer(max_unfTriple * 3)
    unfTriple_count_host <- integer(1)
    version_unf0 <- 2L
    version_unf1 <- 1L
    dyn.load("cuda/tpcConsIntern.so")
    # Make the .C call
    result <- .C("tpcConsIntern",
                as.integer(p),
                as.integer(M),
                as.integer(nrows),
                G_host = as.integer(G_host),
                tiers_host = as.integer(tiers_host),
                context_all_host = as.integer(context_all_host),
                as.integer(n_context_all),
                context_tier_host = as.integer(context_tier_host),
                as.integer(n_context_tier),
                C_host = as.double(C_array),
                as.double(alpha),
                as.integer(max_unfTriple),
                ambiguous_triples_host = as.integer(ambiguous_triples_host),
                unfTriple_count_host = as.integer(unfTriple_count_host),
                as.integer(version_unf0),
                as.integer(version_unf1),
                as.logical(maj.rule),
                SepSet_flat_host = as.integer(SepSet_flat_host),
                SepSet_offsets_host = as.integer(SepSet_offsets_host),
                SepSet_lengths_host = as.integer(SepSet_lengths_host),
                as.integer(SepSet_flat_length)
    )

    # Retrieve the updated G_host and SepSet_flat_host
    G_updated <- matrix(result$G_host, nrow = p, ncol = p)
    skel@graph <- as(G_updated, "graphNEL")

    # Update the SepSet in skel
    SepSet_flat_host <- result$SepSet_flat_host
    SepSet_offsets_host <- result$SepSet_offsets_host
    SepSet_lengths_host <- result$SepSet_lengths_host
    for (a in 1:p) {
      for (c in 1:p) {
        idx <- (a - 1) * p + (c - 1)
        len <- SepSet_lengths_host[idx + 1]
        if (len > 0) {
          start <- SepSet_offsets_host[idx + 1] + 1
          end <- start + len - 1
          sepset_ac <- SepSet_flat_host[start:end] + 1
          skel@sepset[[a]][[c]] <- sepset_ac
        } else {
          skel@sepset[[a]][[c]] <- NULL
        }
      }
    }

    # Handle ambiguous triples if needed
    unfTriple_count <- result$unfTriple_count_host
    if (unfTriple_count > 0) {
      ambiguous_triples <- matrix(result$ambiguous_triples_host[1:(unfTriple_count * 3)], ncol = 3, byrow = TRUE)
      ambiguous_triples <- ambiguous_triples + 1  # Adjust indices back to R's 1-based indexing
    } else {
      ambiguous_triples <- matrix(0, nrow = 0, ncol = 3)
    }

    skelII <- list(sk = skel, unfTripl = ambiguous_triples)
  }
  else{
     skelII <- tpc.cons.intern(skel, suffStat, indepTest, alpha, version.unf = c(2, 1), maj.rule = maj.rule,
                               verbose = verbose, tiers=tiers, context.all=context.all, context.tier=context.tier,
                               forbEdges = forbEdges)
  }
  
  ## step III, orientation of edges between tiers:
  gIII <- as(skelII$sk@graph, "matrix")
  for (t in unique(tiers)) {
    gIII[tiers>t, tiers==t] <- 0
  }
  ## context variables:
  for (i in context.all) {
    gIII[ ,i] <- 0
  }
  for (i in context.tier) {
    k <- tiers[i]
    gIII[tiers==k,i] <- 0
  }
  ## edges in forbArrowList:
  for ( r in forbArrowList ) {
    gIII[r[1],r[2]] <- 0
  }

  # step IV, Meek's rules
  skelIII <- skelII$sk
  skelIII@graph <- as(gIII, "graphNEL")
  source("R/MeeksRules.R")
  MeekRules(skelIII, verbose = verbose, unfVect = skelII$unfTripl,
            solve.confl = TRUE)
}
library(Rfast)
library(mice)
tskeleton_cuda_MI <- function (suffStat, indepTest, alpha, labels, p,
                       method = c("cuda"), m.max = Inf,
                       fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                       tiers = NULL, verbose = FALSE) {
     
     cl <- match.call()
     if (!missing(p))
        stopifnot(is.numeric(p), length(p <- as.integer(p)) ==
                     1, p >= 2)
     if (missing(labels)) {
        if (missing(p))
           stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
     }   else {
        stopifnot(is.character(labels))
        if (missing(p))
           p <- length(labels)
        else if (p != length(labels))
           stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
     }
     seq_p <- seq_len(p)
     method <- match.arg(method)
     if (is.null(fixedGaps)) {
        G <- matrix(TRUE, nrow = p, ncol = p)
     } else if (!identical(dim(fixedGaps), c(p, p)))
        stop("Dimensions of the dataset and fixedGaps do not agree.")
     else if (!identical(fixedGaps, t(fixedGaps)))
        stop("fixedGaps must be symmetric")
     else G <- !fixedGaps
     diag(G) <- FALSE
     #################################################
     ## if no tiers are specified, everything is tier 0
     if (is.null(tiers)) {
        tiers <- rep(0, p)
     } else {
        ## check if 'tiers' are correctly specified
        if (!is.numeric(tiers)) {stop("'tiers' must be a numeric vector")}
        if (length(tiers) != p) {stop("length of 'tiers' does not match 'p' or length of 'labels'")}
     }
     #################################################
     if (any(is.null(fixedEdges))) {
        fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
     }
     else if (!identical(dim(fixedEdges), c(p, p)))
        stop("Dimensions of the dataset and fixedEdges do not agree.")
     else if (!identical(fixedEdges, t(fixedEdges)))
        stop("fixedEdges must be symmetric")
    
    sepset <- lapply(seq_p, function(.) vector("list", p))
    # pMax is a matrix with one p-value per edge, at the beginning all p-values
    # are -Inf
    pMax <- matrix(0, nrow = p, ncol = p)
    n <- suffStat[length(suffStat)]
    m <- length(suffStat) - 1
    C_list <- head(suffStat, m)
    C_array <- array(0, dim = c(p, p, m))


    for (i in 1:m) {
        C_array[, , i] <- C_list[[i]]
    }
    # replace NA with 0.0, this is how it is handled in pcalg package
    C_array[is.na(C_array)] <- 0.0
    C_vector <- as.vector(C_array) # is this needed?
    
    # Initialize adjacency matrix G
    G <- matrix(TRUE, nrow = p, ncol = p)
    diag(G) <- FALSE
    ord <- 0
    G <- G * 1 # Convert logical to integer

    # Determine maximum levels
    if (m.max == Inf) {
        max_level <- 32
    } else {
        max_level <- m.max
    }
    sepsetMatrix <- matrix(-1, nrow = p * p, ncol = 32)
    dyn.load("cuda/SkeletonMI.so")
    start_time <- proc.time()
    z <- .C("SkeletonMI",
        C = as.double(C_vector),
        p = as.integer(p),
        Nrows = as.integer(n),
        m = as.integer(m),
        G = as.integer(G),
        Alpha = as.double(alpha),
        l = as.integer(ord),
        max_level = as.integer(max_level),
        pmax = as.double(pMax),
        sepsetmat = as.integer(sepsetMatrix),
        tiers = as.integer(tiers)
    )
    ord <- z$l
    G <- (matrix(z$G, nrow = p, ncol = p)) > 0

    pMax <- (matrix(z$pmax, nrow = p, ncol = p))
    pMax[which(pMax == -100000)] <- -Inf

    sepsetMatrix <- t(matrix(z$sepsetmat, nrow = 32, ncol = p^2))
    #print(sepsetMatrix)
    index_of_cuted_edge <- row(sepsetMatrix)[which(sepsetMatrix != -1)]
    for (i in index_of_cuted_edge) {
        edge_idx <- i - 1  # Adjust for R's 1-based indexing
        x <- (edge_idx %% p) + 1
        y <- (edge_idx %/% p) + 1

        # find the last non -1 entry in the sepset row
        sepset_entries <- sepsetMatrix[i, ]
        #cat("x", x, "y", y, "sepset_entries", sepset_entries, "\n")
        valid_entries <- sepset_entries[sepset_entries != -1]

        # assign the separation set
        sepset[[x]][[y]] <- sepset[[y]][[x]] <- if (any(valid_entries == 0)) {
                                integer(0)
                            } else {
                                valid_entries
                            }
        # sepset[[x]][[y]] <- if (any(valid_entries == 0)) {
        #                         integer(0)
        #                     } else {
        #                         valid_entries
        #                     }
    }

   Gobject <- if (sum(G) == 0) {
      new("graphNEL", nodes = labels)
   } else {
      colnames(G) <- rownames(G) <- labels
      as(G, "graphNEL")
   }
   new("pcAlgo", graph = Gobject, call = cl, n = integer(0),
       max.ord = as.integer(ord - 1), n.edgetests = 0,
       sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1))
}
