library(survival) # for Surv() and easy estimates via survfit()
library(MASS) # for ginv()
library(dplyr) # for %>%, group_by() etc.
library(purrr) # for cross()

# function: simple_surv() determines basic survival quantities in a faster 
#           manner than survfit()
# input: time - survival time
#        cens - indicator for censoring (1=censored, 0=not censored)
# output: matrix of dimensions (n, 4) to speed up later calculations
simple_surv <- function(time, cens) {
  sur <- Surv(time, cens)
  n.all <- length(sur)
  tab <- table(sur[,1], factor(sur[,2], levels = c(0,1)))
  d <- tab[,1] # number of events in time points
  w <- tab[,2] # number of withdrawals in time points
  n <- c(n.all, n.all - cumsum(d) - cumsum(w)) # calculate risk set
  n <- n[-length(n)]
  s <- cumprod(1 - (w /  n)) # calculate Kaplan-Meier-Estimator
  matrix(c(as.numeric(row.names(tab)), s, w, n), ncol = 4)
}

# function: RMST() estimates Restricted Mean Survival Time and its variance 
#           for a sample of tupel (time, cens), given a time point tau
# input: time - survival time
#        cens - indicator for censoring (1=censored, 0=not censored)
#        tau - restriction time point
# output: rmst - estimated restricted mean survival time for tau
#         var_rmst - estimated variance of the rmst
RMST <- function(time, cens, tau) {
  n <- length(time) # number of observation in the beginning of study
  survtab <- simple_surv(time, cens) # fit convenient model for quantities
  
  t_i <- survtab[,1] <= tau # identify time points t_i <= tau
  t_sel <- survtab[,1][t_i] # select relevent time points
  S_km <- survtab[,2][t_i][-sum(t_i)] # calculate Kaplan-Meier estimator
  
  w.factor <- diff(c(0, t_sel)) # width of the area under S(t) from t_0=0
  rmst <- sum(w.factor * c(1, S_km)) # calculate area under S(t) up to t_i=tau
  
  area <- rev(cumsum(rev(diff(t_sel) * S_km)))^2 # calculate areas under S(t) from t_i to tau for
  d <- survtab[,3][t_i][-sum(t_i)] # determine number of events
  Y <- survtab[,4][t_i][-sum(t_i)] # determine number of individuals under risk
  var_rmst <- sum(((n * d) / Y^2) * area) # calculate final variance
  
  c(rmst, var_rmst)
}

# function: RMSTstatistic() calculates the test statistic
# input: X - event data, with factorial indicators
#        tau - restriction time point
#        nall - number of whole sample
#        Tmat - contrast matrix for specific null hypothesis
# output: value of test statistic
RMSTstatistic <- function(X, tau, nall, Tmat) {
  # order data in consistent way
  X <- group_by(X, A, B)
  
  # calculate sample size by group combination, try to use count()
  nk <- summarise(X, n = n(), .groups = "keep")
  nk <- nk$n
  
  # estimate RMST and its variance
  mom <- summarise(X, rmst = list(RMST(time, cens, tau)), .groups = "keep")
  mom <- matrix(unlist(mom$rmst), nrow = 6, ncol = 2, byrow = TRUE)
  
  # precalculate matrix
  # (seminar_rmst_simulateH0_seq_29836.1592505663-5.out: Error in Tmat %*% mom[, 1] : non-conformable arguments)
  Tmat_mu <- try(Tmat %*% mom[, 1, drop = FALSE], silent = TRUE)
  if(!is.matrix(Tmat_mu)) {
    save(X, file = "dasproblem.RData")
    return(X)
  }
  
  # estimate covariance matrix of RMST by group
  sigma <- diag((nall / nk) * mom[,2])
  
  # calculate teststatistic
  drop(nall * (t(Tmat_mu) %*% ginv(Tmat %*% sigma %*% t(Tmat)) %*% Tmat_mu))
}

# helper functions for constructing the contrast matrix
J <- function(dim) matrix(1, ncol = dim, nrow = dim)
P <- function(dim) diag(dim) - J(dim) / dim

# function: testRMST() tests for different for main and interaction effects 
#           in different hypotheses
# input: dat - dataframe with relevant information
#        time - column name of event time variable as character
#        cens - column name of censoring indicator as character
#        fA - column name of factor A as character
#        fB - column name of factor B as character
# output: 
testRMST <- function(dat, time, cens, fA, fB, alternative, tau, 
                     nperm = 1999, 
                     permDist = FALSE,
                     perm = FALSE) {
  
  # prepare dataframe
  X <- dat %>% dplyr::select(time = time, cens = cens, A = fA, B = fB) %>% 
    mutate(A = as.factor(A), B = as.factor(B))

  # calculate number of observations
  nall <- nrow(X)
  
  # determine number of factor levels
  a <- length(levels(X$A))
  b <- length(levels(X$B))
  
  # determine the contrast matrix
  if (alternative == "A") H <- P(a) %x% (J(b) / b)
  if (alternative == "B") H <- (J(a) / a) %x% P(b)
  if (alternative == "AB") H <- P(a) %x% P(b)
  
  # calculate projection matrix
  Tmat <- t(H) %*% ginv(H %*% t(H)) %*% H
  
  # calculate rank for asymptotic test decision
  rankT <- qr(Tmat)$rank
  
  # calculate the test statistic for the observed group affiliation
  trueW <- RMSTstatistic(X = X, nall = nall, tau = tau, Tmat = Tmat)
  
  if(perm) {
  # calculate the permutation distribution of the test statistic (approach 1)
    permW <- replicate(nperm, RMSTstatistic(X = X %>% mutate(A = sample(A), 
                       B = sample(B)), nall = nall, tau = tau, Tmat = Tmat))

  # calculate p-values (asymptotic and permutation)
    res <- list(asymp = pchisq(trueW, df = rankT, lower.tail = FALSE),
                perm = mean(trueW < permW))
    if(permDist) {
      res$perm_dist <- permW
      res$asymp_rank <- rankT
      res$trueW <- trueW
    }
  } else {
    res <- list(asymp = pchisq(trueW, df = rankT, lower.tail = FALSE))
  }
  
  return(res)
}

# function: calc_U() simulates event times to approximate the upper bound of
#           unified distribution of censoring times to obtain a certain
#           censoring rate
# input: T_cens - target censoring rate 
#        T_dist - sample of specific event time distribution of size n
# output: U - upper bound for the unified distribution for fixed censoring rate
calc_U <- function(T_cens, T_dist) {
  n <- length(T_dist)
  
  # sort data for following calculations
  dat <- sort(T_dist)
  
  # approximate first integral by empirical cumulative distribution function (ecdf)
  first_int <- (cumsum(dat) / dat) / n 
  
  # approximate second integral by 1 - ecdf
  second_int <- ((n - 1):0) / n # noticable faster than 1 - ecdf(dat)(dat)
  
  # sum up integrals
  final_int <- sort(first_int + second_int)
  
  # find bound, which leads to correct censoring rate 
  i <- findInterval(T_cens, final_int)
  
  # select the upper bound 
  U <- dat[n - i]
  
  if(i == 0) {
    # for low censoring right and therefore high U the second integral is 
    # approximately 0 and the first one can be approximated by E(T)/U and
    # therefore the distance between mean(T)/U and t_cens should be minimized
    # for large values of T \in [max(t), 100*max(t)]
    U <- optimize(function(u) (mean(dat)/u - T_cens)^6, 
                  lower = dat[n], upper = 100 * dat[n])$minimum
  }
  
  U
}

# survival function and pseudo random number generation of the piecewise exponential hazard
pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
pweR <- function(n, salt, h1, h2) {
  pre_salt <- rexp(n, rate = h1-h2)
  pre_salt <- ifelse(pre_salt <= salt, pre_salt, Inf)
  post_salt <- rexp(n, rate = h2)
  apply(cbind(pre_salt, post_salt), 1, min)
}
# calculation of the RMST for a piecewise constant hazard set up with variable
# discontinuity point and fixed tau
pwRMST <- function(x, end_tau) {
  integrate(function (t) pweS(t, salt  = x, h1 = 3, h2 = 12/35), 
            lower = 0, upper = end_tau, subdivisions = 10^8L)$value
}

# function: simDataH0() simulates data under H0
# input: nk - vector of sample sizes of the k groups
#        crate - vector of censoring rates in the k groups
#        dist - event time distribution
#        alternative - type of the hypothesis
#        U_j1 - pre calculated upper limit for uniform distributed censoring time for alternative
#        U_j2 - pre calculated upper limit for uniform distributed censoring time for null
# output: list of 2 variables:
#           dat - dataframe of 4 variables
#             time - time to event
#             cens - indicator for censoring (1=censored, 0=not censored)
#             A - factor A as character
#             B - factor B as character
#           tau - value for tau depending on the distribution
simDataH0 <- function(nk,
                      crate,
                      dist,
                      alternative, 
                      U_j1, 
                      U_j2) {
  # indicator for which to "disturb" the distribution of event time
  if(alternative == "A") alt_i <- c(1, 2, 3)
  if(alternative == "AB") alt_i <- c(1, 2)
  
  if (dist == "exp") {
    tau <- 4
    time <- unlist(sapply(nk, function(x) rexp(x, rate = 1)))
    cens_time <- unlist(sapply(nk, function(x) rexp(x, rate = crate / (1 - crate))))
  }
  
  if (dist == "exp_cross") {
    saltus <- 1 / 6
    hr1 <- 3
    hr2 <- 12 / 35
    t_rate <- 1
    tau <- optimize(function(x) (integrate(function(x) 1 - pweS(x, 1/6, 3, 12/35), lower = 0, upper = x)$value -
                          integrate(pexp, lower = 0, upper = x)$value)^30,
             lower = 0, upper = 5, tol = .Machine$double.eps)$minimum
    
    time1 <- unlist(sapply(nk[alt_i], function(x) pweR(x, salt = saltus, h1 = hr1, h2 = hr2)))
    cens_time1 <- unlist(apply(rbind(nk[alt_i], U_j1), 2 ,function(x) runif(x[1], 0, x[2])))
    
    time2 <- unlist(sapply(nk[-alt_i], function(x) rexp(x, rate = t_rate)))
    cens_time2 <- unlist(apply(rbind(nk[-alt_i], U_j2), 2 ,function(x) runif(x[1], 0, x[2])))
    
    time <- c(time1, time2)
    cens_time <- c(cens_time1, cens_time2)
  }
  
  if (dist == "weibull") {
    t_shape <- 0.4
    t_scale <- 0.8
    tau <- 1
    time <- unlist(sapply(nk, function(x) rweibull(x, shape = t_shape, scale = t_scale)))
    cens_time <- unlist(apply(rbind(nk, U_j1), 2, function(x) runif(x[1], 0, x[2])))
  }
  
  if (dist == "lnorm") {
    t_mean <- 0
    t_sd <- 1/4
    tau <- 1.5
    time <- unlist(sapply(nk, function(x) rlnorm(x, meanlog = t_mean, sdlog = t_sd)))
    cens_time <- unlist(apply(rbind(nk, U_j1), 2 ,function(x) runif(x[1], 0, x[2])))
  }
  
  grid <- as.matrix(expand.grid(A = factor(c("1", "2")), B = factor(c("1", "2", "3"))))
  grid <- grid[order(grid[, 1]), ]
  grid <- as.data.frame(matrix(rep(c(grid), c(nk, nk)), ncol = 2))
  names(grid) <- c("A", "B")
  res <- data.frame(time = apply(cbind(c(time), c(cens_time)), 1, min), 
                    cens = as.numeric(cens_time >= time))
  res <- cbind(res, grid)
  
  return(list(dat = res, tau = tau))
}

# function: simDataH1_shift() simulates data under H1 by 
# input: nk - vector of sample sizes of the k groups
#        crate - vector of censoring rates in the k groups
#        alternative - type of the hypothesis
#        saltus - point for the jump in piecewise constant hazard rate
#        U_j1 - pre calculated upper limit for uniform distributed censoring time for alternative
#        U_j2 - pre calculated upper limit for uniform distributed censoring time for null
# output: list of 2 variables:
#           dataframe of 4 variables
#             time - time to event exp(1)-distributed disturbed by piecewise exponential distributions
#             cens - indicator for censoring (1=censored, 0=not censored)
#             A - factor A as character
#             B - factor B as character
#           value for tau depending on the distribution
simDataH1_shift <- function(nk,
                            crate,
                            alternative,
                            saltus, 
                            U_j1, 
                            U_j2) {
  if (alternative == "A") alt_i <- c(1, 2, 3)
  if (alternative == "AB") alt_i <- c(1, 2)
  
  hr1 <- 3
  hr2 <- 12 / 35
  t_rate <- 1
  tau <- optimize(function(x) (integrate(function(x) 1 - pweS(x, 1/6, 3, 12/35), lower = 0, upper = x)$value -
                                 integrate(pexp, lower = 0, upper = x)$value)^30,
                  lower = 0, upper = 5, tol = .Machine$double.eps)$minimum

  time1 <- unlist(sapply(nk[alt_i], function(x) pweR(x, salt = saltus, h1 = hr1, h2 = hr2)))
  cens_time1 <- unlist(apply(rbind(nk[alt_i], U_j1), 2 ,function(x) runif(x[1], 0, x[2])))
  time2 <- unlist(sapply(nk[-alt_i], function(x) rexp(x, rate = t_rate)))
  cens_time2 <- unlist(apply(rbind(nk[-alt_i], U_j2), 2 ,function(x) runif(x[1], 0, x[2])))
    
  time <- c(time1, time2)
  cens_time <- c(cens_time1, cens_time2)
  grid <- as.matrix(expand.grid(A = factor(c("1", "2")), B = factor(c("1", "2", "3"))))
  grid <- grid[order(grid[, 1]), ]
  grid <- as.data.frame(matrix(rep(c(grid), c(nk, nk)), ncol = 2))
  names(grid) <- c("A", "B")
  res <- data.frame(time = apply(cbind(c(time), c(cens_time)), 1, min), 
                    cens = as.numeric(cens_time >= time))
  res <- cbind(res, grid)
  
  return(list(dat = res, tau = tau))
}