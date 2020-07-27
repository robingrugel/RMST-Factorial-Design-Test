source("seminar_rmst_test.r")
library(parallel) # obvious
library(rlecuyer) # for meaningful pseudo RNG

nsim <- 5000
iter <- 1999
ncores <- 4
# replicate=1,1-144
# 512 Mb should be sufficient

p <- as.integer(Sys.getenv("PBS_ARRAYID"))
print(p)

p_list <- list(
  # determine alternatives
  hyp = c("A","AB"),
  
  # factor to increase group sizes
  fac = c(1, 2, 4),
  
  # group size designs
  n_i = list(n1 = c(10, 10, 10, 10, 10, 10), # balanced
             n2 = c(8, 12, 14, 10, 8, 8)), # unbalanced
  
  # censoring rates
  c_i = list(lo = c(0.08, 0.07, 0.1, 0.09, 0.06, 0.06), # low
             mid = c(0.16, 0.25, 0.21, 0.2, 0.15, 0.23), # intermediate
             hi = c(0.41, 0.45, 0.29, 0.35, 0.4, 0.34)), # high
  
  t_dist = c("exp", 
             "exp_cross", 
             "weibull", 
             "lnorm")
) %>% cross()

.lec.SetPackageSeed(c(6,13,73,4,52,1)) # rlecuyer Aequivalent zu set.seed()
nstreams <- length(p_list) # Anzahl der Zufallszahlenstreams
names <- paste("myrngstream",1:nstreams,sep="") # Irgendein Name fur den RNG Stream
.lec.CreateStream(names) # Zufallszahlenstreams erstellen
.lec.CurrentStream(names[p]) # Auswahl des p-ten Streams.

if (p_list[[p]]$hyp == "A") alt_i <- c(1, 2, 3)
if (p_list[[p]]$hyp == "AB") alt_i <- c(1, 2)

if(p_list[[p]]$t_dist == "exp") {
  U_j1 <- 1
  U_j2 <- 1
}
if(p_list[[p]]$t_dist == "exp_cross") {
  U_j1 <- unlist(sapply(p_list[[p]]$c_i[alt_i], function(x) calc_U(x, T_dist = pweR(10^6, salt = 1/6, h1 = 3, h2 = 12/35))))
  U_j2 <- unlist(sapply(p_list[[p]]$c_i[-alt_i], function(x) calc_U(x, T_dist = rexp(10^6, rate = 1))))
}
if(p_list[[p]]$t_dist == "weibull") {
  U_j1 <- unlist(sapply(p_list[[p]]$c_i, function(x) calc_U(x, T_dist = rweibull(10^6, shape = 0.4, scale = 0.8))))
  U_j2 <- U_j1
}
if(p_list[[p]]$t_dist == "lnorm") {
  U_j1 <- unlist(sapply(p_list[[p]]$c_i, function(x) calc_U(x, T_dist = rlnorm(10^6, meanlog = 0, sdlog = 1/4))))
  U_j2 <- U_j1
}

RNGkind("L'Ecuyer-CMRG")
res <- mclapply(1:nsim, function(x) {
  that_dat <- simDataH0(nk = p_list[[p]]$fac * p_list[[p]]$n_i, 
                        crate = p_list[[p]]$c_i, 
                        dist = p_list[[p]]$t_dist, 
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
res$dist <- p_list[[p]]$t_dist

namethis <- paste0("f_resH0_", p, "_par")
assign(namethis, res)

save(list = namethis, file=paste0("f_resH0_", p, "_par.RData"))
gc()