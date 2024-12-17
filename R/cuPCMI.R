library(pcalg)
library(Rfast)
library(mice)

cu_pc_MI <- function(suffStat, indepTest, alpha, labels, p,
                     fixedGaps = NULL, fixedEdges = NULL, NAdelete = TRUE,
                     m.max = Inf, u2pd = c("relaxed", "rand", "retry"),
                     skel.method = c("stable", "original", "stable.fast"),
                     conservative = FALSE, maj.rule = FALSE,
                     solve.confl = FALSE, verbose = FALSE, df_method = c("old", "br", "reiter")) {
    ## Initial Checks
    cl <- match.call()
    if (!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
    if (missing(labels)) {
        if (missing(p)) stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
    } else { ## use labels ==> p  from it
        stopifnot(is.character(labels))
        if (missing(p)) {
            p <- length(labels)
        } else if (p != length(labels)) {
            stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
        } else {
            message("No need to specify 'p', when 'labels' is given")
        }
    }
    seq_p <- seq_len(p)

    u2pd <- match.arg(u2pd)
    skel.method <- match.arg(skel.method)
    if (u2pd != "relaxed") {
        if (conservative || maj.rule) {
            stop("Conservative PC and majority rule PC can only be run with 'u2pd = relaxed'")
        }

        if (solve.confl) {
            stop("Versions of PC using lists for the orientation rules (and possibly bi-directed edges)\n
            can only be run with 'u2pd = relaxed'")
        }
    }

    if (conservative && maj.rule) stop("Choose either conservative PC or majority rule PC!")

    ## Skeleton
    skel <- cu_skeleton_MI(suffStat, indepTest, alpha,
        labels = labels, NAdelete = NAdelete, m.max = m.max, verbose = verbose, df_method = df_method
    )
    skel@call <- cl # so that makes it into result
    ## Orient edges
    if (!conservative && !maj.rule) {
        switch (u2pd,
                "rand" = udag2pdag(skel),
                "retry" = udag2pdagSpecial(skel)$pcObj,
                "relaxed" = udag2pdagRelaxed(skel, verbose = verbose, solve.confl = solve.confl))
    } else { ## u2pd "relaxed" : conservative _or_ maj.rule
        ## version.unf defined per default
        ## Tetrad CPC works with version.unf=c(2,1)
        ## see comment on pc.cons.intern for description of version.unf
        pc. <- pc.cons.intern(skel, suffStat, indepTest, alpha,
                            version.unf = c(2,1), maj.rule = maj.rule, verbose = verbose)
        udag2pdagRelaxed(pc.$sk, verbose = verbose,
                        unfVect = pc.$unfTripl, solve.confl = solve.confl)
    }
} ## {pc}

cu_skeleton_MI <- function(suffStat, indepTest, alpha, labels, p, m.max = Inf, NAdelete = TRUE, verbose = FALSE, df_method) {
    cl <- match.call()
    if (!missing(p)) stopifnot(is.numeric(p), length(p <- as.integer(p)) == 1, p >= 2)
    if (missing(labels)) {
        if (missing(p)) stop("need to specify 'labels' or 'p'")
        labels <- as.character(seq_len(p))
    } else { ## use labels ==> p  from it
        stopifnot(is.character(labels))
        if (missing(p)) {
            p <- length(labels)
        } else if (p != length(labels)) {
            stop("'p' is not needed when 'labels' is specified, and must match length(labels)")
        } else {
            message("No need to specify 'p', when 'labels' is given")
        }
    }

    seq_p <- seq_len(p)
    pval <- NULL
    # Convert SepsetMatrix to sepset
    sepset <- lapply(seq_p, function(.) vector("list", p)) # a list of lists [p x p]
    # save maximal p value
    pMax <- matrix(0, nrow = p, ncol = p)
    
    m <- length(suffStat) - 1
    n <- suffStat[length(suffStat)]
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
    done <- TRUE
    G <- G * 1 # Convert logical to integer

    # Determine maximum levels
    if (m.max == Inf) {
        max_level <- 32
    } else {
        max_level <- m.max
    }

    sepsetMatrix <- matrix(-1, nrow = p * p, ncol = 32)
    dyn.load("cuda/SkeletonMI.so")


    if (df_method == "old"){
        df_method_int <- 0
    }
    else if (df_method == "br"){
        df_method_int <- 1
    }
    else if (df_method == "reiter"){
        df_method_int <- 2
    }
    else{
        stop("df method is not specified correctly. Should be either 'old', 'br', or 'reiter'")
    }
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
        tiers = as.integer(rep(0, p)),
        DF_method = as.integer(df_method_int)
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
        # sepset[[x]][[y]] <- sepset[[y]][[x]] <- valid_entries + 1
        sepset[[x]][[y]] <- if (any(valid_entries == 0)) {
                                integer(0)
                            } else {
                                valid_entries
                            }
    }

    # print(ord)
    ## transform matrix to graph object :
    Gobject <-
        if (sum(G) == 0) {
            new("graphNEL", nodes = labels)
        } else {
            colnames(G) <- rownames(G) <- labels
            as(G, "graphNEL")
        }
    ## final object
    new("pcAlgo",
        graph = Gobject, call = cl, n = integer(0),
        max.ord = as.integer(ord - 1), n.edgetests = 0,
        sepset = sepset, pMax = pMax, zMin = matrix(NA, 1, 1)
    )
} ## end{ skeleton }

# copied and changed from micd github
getSuffCU <- function(X) {
    if (!(mice::is.mids(X) || is.list(X))) {
        stop("Data is neither a list nor a mids object.")
    }

    if (inherits(X, "mids")) {
        X <- mice::complete(X, action = "all")
        if (any(sapply(X, is.factor))) {
            stop("Data must be all numeric.")
        }
    }

    C_list <- lapply(X, cor)
    p <- ncol(X[[1]])
    m <- length(C_list)

    C_array <- array(0, dim = c(p, p, m))

    for (i in 1:m) {
        C_array[, , i] <- C_list[[i]]
    }
    # potentially problematic, but I don't think so
    # replace NA with 0.0, this is how it is handled in pcalg package
    C_array[is.na(C_array)] <- 0.0
    return(list(C = C_array, n = nrow(X[[1]]), m = m))
}
