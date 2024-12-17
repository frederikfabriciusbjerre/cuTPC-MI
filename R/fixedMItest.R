fixedGaussMItest <- function (x, y, S, suffStat) 
{
    # number of imputations
    M <- length(suffStat) - 1
    # sample size
    n <- suffStat[[length(suffStat)]]

    z <- sapply(head(suffStat, -1), function(j) {
        zStatMI(x, y, S, C=j)
    })

    # 1. Average of M imputed data sets
    avgz <- mean(z)

    # 2. Average of completed-data variance
    W <- 1 / (n - length(S) - 3)

    # 3. Between variance
    B <- sum( ( z - avgz )^2 ) / (M-1)

    ################################
    # if single imputation / m = 1
    if (M == 1){
        return (2 * stats::pnorm(abs(avgz/sqrt(W)), lower.tail = FALSE))
    }
    ################################
    
    # 4. Total variance
    TV <- W + (1 + 1 / M) * B

    # 5. Test statistic
    ts <- avgz / sqrt(TV)

    # 6. Degrees of freedom
    df <- (M - 1) * (1 + (W / B) * (M/(M + 1)))^2

    # 7. pvalue
    pvalue <- 2 * stats::pt(abs(ts), df = df, lower.tail = FALSE)

    return(pvalue)
}

zStatMI <- function (x, y, S, C)
{
    r <- pcalg::pcorOrder(x, y, S, C)
    res <- 0.5 * log_q1pm(r)
    if (is.na(res))
        0
    else res
}

log_q1pm <- function(r) log1p(2 * r / (1 - r))
