source("seminar_rmst_test.r")
library(parallel) # obvious
library(rlecuyer) # for meaningful pseudo RNG

nsim <- 5000
iter <- 1999
ncores <- 4
# replicate=1,1-108
# 512 Mb should be sufficient

# pre-calculate tau for discontinuity points 
p_tau <- optimize(function(x) (integrate(function(x) 1 - pweS(x, 1/6, 3, 12/35), lower = 0, upper = x)$value -
                               integrate(pexp, lower = 0, upper = x)$value)^30,
                lower = 0, upper = 5, tol = .Machine$double.eps^2)$minimum

# calculate RMST for Exp(1)
expRMST <- integrate(function(x) 1 - pexp(x), lower = 0, upper = p_tau)$value

p <- as.integer(Sys.getenv("PBS_ARRAYID"))
print(p)

p_list <- list(
  # determine alternatives
  hyp = c("A", "AB"),
  
  # factor to increase group size
  fac = c(1, 2, 4),
  
  # group sizes
  n_i = list(c(8, 12, 14, 10, 8, 8)), # unbalanced
  
  # censoring rates
  c_i = list(
    lo = c(0.08, 0.07, 0.1, 0.09, 0.06, 0.06), # low
             mid = c(0.16, 0.25, 0.21, 0.2, 0.15, 0.23)
    # , # intermediate
             # hi = c(0.41, 0.45, 0.29, 0.35, 0.4, 0.34)
             ), # high
  
  # discontinuity point for pwe hazard so the RMST increases linear
  salt = sapply(seq(0, 0.4, 0.05), function(i) 
                optimize(function(s) (i - abs(expRMST - pwRMST(s, p_tau)))^30, 
                         lower = 1/6, upper = p_tau, tol = .Machine$double.eps)$minimum) 
) %>% cross()

.lec.SetPackageSeed(c(6,13,73,4,52,1)) # rlecuyer Aequivalent zu set.seed()
nstreams <- length(p_list) # Anzahl der Zufallszahlenstreams
names <- paste("myrngstream",1:nstreams,sep="") # Irgendein Name fur den RNG Stream
.lec.CreateStream(names) # Zufallszahlenstreams erstellen
.lec.CurrentStream(names[p]) # Auswahl des p-ten Streams.

if (p_list[[p]]$hyp == "A") alt_i <- c(1, 2, 3)
if (p_list[[p]]$hyp == "AB") alt_i <- c(1, 2)
U_j1 <- unlist(sapply(p_list[[p]]$c_i[alt_i], function(x) calc_U(x, T_dist = pweR(10^6, salt = p_list[[p]]$salt, h1 = 3, h2 = 12/35))))
U_j2 <- unlist(sapply(p_list[[p]]$c_i[-alt_i], function(x) calc_U(x, T_dist = rexp(10^6, rate = 1))))
  
RNGkind("L'Ecuyer-CMRG")
res <- mclapply(1:nsim, function(x) {
  that_dat <- simDataH1_shift(nk = p_list[[p]]$fac * p_list[[p]]$n_i, 
                              crate = p_list[[p]]$c_i, 
                              saltus = p_list[[p]]$salt, 
                              alternative = p_list[[p]]$hyp,
                              U_j1 = U_j1,
                              U_j2 = U_j2);
    testRMST(time = "time", 
             cens = "cens", 
             fA = "A", 
             fB = "B", 
             dat = that_dat$dat, 
             tau = that_dat$tau, 
             alternative = p_list[[p]]$hyp, 
             nperm = iter, 
             perm = TRUE)
}, mc.cores = ncores, mc.set.seed = TRUE)

res <- matrix(unlist(res), nrow = 2)
res <- data.frame(asymp = res[1,], perm = res[2,])
res$hyp <- p_list[[p]]$hyp
res$fac <- p_list[[p]]$fac
res$n_i <- paste0(p_list[[p]]$n_i, collapse = "")
res$c_i <- paste0(p_list[[p]]$c_i, collapse = "")
res$saltus <- p_list[[p]]$salt
  
namethis <- paste0("f_resH1_shifted_", p, "_par")
assign(namethis, res)
  
save(list = namethis, file=paste0("f_resH1_shifted_", p, "_par.RData"))
gc()