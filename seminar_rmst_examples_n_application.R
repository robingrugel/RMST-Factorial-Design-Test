source("seminar_rmst_test.r")
library(tables) # for latex tables
library(ggplot2) # for nice plots
library(condSURV) # for data example
library(ggpubr)
library(tidyr) # for pivot_longer
library(survminer)
library(binom)
library(xtable)
# library(modelr) # no idea

################################################################################
# motivation using real data example
################################################################################

Exp_setup1 <- function() {
  ggplot() + 
    stat_function(fun = function(x) 1 - pexp(x, rate = 1), 
                  n = 10000, aes(color = "$1$")) + 
    stat_function(fun = function(x) 1 - pexp(x, rate = 0.5), 
                  n = 10000, aes(color = "$0.5$")) + 
    stat_function(fun = function(x) 1 - pexp(x, rate = 1.5), 
                  n = 10000, aes(color = "$1.5$")) + 
    scale_colour_manual(name = "$\\lambda$", values = 1:3) +
    xlim(c(0, 5)) + xlab("$t$") + ylab("$f(t)$") + theme_bw() 
}

################################################################################
############################# Distribution of W(T) #############################
################################################################################
plotDistributions <- function(iter = c(19, 199)) {
  set.seed(1234)
  res1 <- sapply(iter, function(l) testRMST(time = "time", cens = "cens", fA = "A", fB = "B",
            dat = simDataH0(nk = rep(10, 6), crate = c(0.08, 0.07, 0.1, 0.09, 0.06, 0.06), 
            dist = "exp", alternative = "A")$dat, tau = 4, alternative = "A", nperm = l,
            permDist = TRUE, perm = TRUE)$perm_dist)
  res2 <- sapply(iter, function(l) testRMST(time = "time", cens = "cens", fA = "A", fB = "B",
              dat = simDataH0(nk = rep(1000, 6), crate = c(0.08, 0.07, 0.1, 0.09, 0.06, 0.06), 
              dist = "exp", alternative = "A")$dat, tau = 1, alternative = "A", nperm = l,
              permDist = TRUE, perm = TRUE)$perm_dist) 
  res3 <- data.frame(x = c(unlist(res1), unlist(res2)), 
                     no_perm = as.factor(rep(rep(c(19, 199), iter), 2)), 
                     n = as.factor(rep(c(10,1000), each = sum(iter))))
  levels(res3$no_perm) <- c("$n_{\\text{perm}} = 19$", "$n_{\\text{perm}} = 199$")
  levels(res3$n) <- c("$n_i = 10$", "$n_i = 1000$")
  ggplot(res3, aes(x, col = "Permutation")) + stat_ecdf() + facet_grid(no_perm ~ n) +
    stat_function(fun = pchisq, args= list(df = 1), aes(col = "$\\chi^2_{\\text{rank}(T)}$")) + 
    theme_bw() + 
    xlab("$W_n(T)$") + ylab("$F_{W_n(T)}(\\cdot)$") +
    scale_color_manual(name = "Distribution", values = c("red", "black"))
}


################################################################################
############################# H_0 simulation setup #############################
################################################################################
plotH0_setup1 <- function() {
  ggplot() + 
  stat_function(fun = function(x) 1 - pexp(x, rate = 1), 
                n = 10000, aes(color = "$\\text{Exp}(1)$")) + 
  stat_function(fun = function(x) 1 - pweibull(x, shape = 0.4, scale = 0.8), 
                n = 10000, aes(color = "$\\text{WB}(0.4, 0.8)$")) +
  stat_function(fun = function(x) 1 - plnorm(x, meanlog = 0, sdlog = 1/4), 
                n = 10000, aes(color = "$\\text{logN}(0, 0.25)$")) + 
  scale_colour_manual(name = "event time distribution", values = 1:3) +
  xlim(c(0, 3)) + xlab("$t$") + ylab("$S_i(t)$") + theme_bw() 
}



# optimize(function(ti) (integrate(function(x) pweS(x, 1/6, 3, 12/35), lower = 0, upper = ti)$value -
#                         integrate(function(a) 1 - pexp(a), lower = 0, upper = ti)$value)^30,
#          lower = 0, upper = 5, tol = .Machine$double.eps^2)$minimum
# 
# fu <- sapply(seq(0, 5, length.out=100000), function(ti) (integrate(function(x) pweS(x, 1/6, 3, 12/35), lower = 0, upper = ti)$value -
#                                                            integrate(function(a) 1 - pexp(a), lower = 0, upper = ti)$value)^30)
# plot(seq(0, 5, length.out=100000), fu, ylim = c(0, 0.004^100))



plotH0_setup2 <- function() {
  # pre-calculate tau for discontinuity points 
  p_tau <- optimize(function(x) (integrate(function(x) 1 - pweS(x, 1/6, 3, 12/35), lower = 0, upper = x)$value -
                                   integrate(pexp, lower = 0, upper = x)$value)^30,
                    lower = 0, upper = 5, tol = .Machine$double.eps)$minimum
  
  ggplot() + 
  stat_function(fun = function(x) 1 - pexp(x, rate = 1), 
                n = 10000, aes(color = "$\\text{Exp}(1)$")) + 
  stat_function(fun = function(x) pweS(x, salt = 1/6, h1 = 3, h2 = 12/35),
                n = 10000, aes(color = "$\\text{pwExp}\\left(\\frac{1}{6}, 3,\\frac{12}{35}\\right)$")) + 
  scale_colour_manual(name = "event time distribution", values = c(1,2)) +
  xlim(c(0, 3)) + ylim(c(0, 1)) + xlab("$t$") + ylab("$S_i(t)$") + 
  geom_segment(aes(x=p_tau,xend=p_tau,
                   y=0,yend=pweS(p_tau, salt = 1/6, h1 = 3, h2 = 12/35), 
                   linetype = "$\\tau \\approx 1.478929$")) + 
  scale_linetype_manual("end point for RMST", values = 2) + theme_bw() 
# +
#   geom_segment(aes(x=0,xend=1.478929,
#                    y=0,yend=0, linetype = "$\\tau=1.478929$")) + 
#   geom_segment(aes(x=0,xend=0,
#                    y=0,yend=1, 
#                    linetype = "$\\tau=1.478929$"))
}

################################################################################
############################ H_0 simulation results ############################
################################################################################
tableH0 <- function(h = "$AB$") {
  flsi <- list.files()
  flsi_i <- grep("resH0_", flsi, value = TRUE)
  lapply(flsi_i, load, .GlobalEnv)
  simresH0 <- do.call(rbind, lapply(paste0("f_resH0_", 1:144, "_par"), get, envir = globalenv()))
  redresH0 <- simresH0 %>% 
    group_by(hyp, fac, n_i, c_i, dist) %>%
    summarize(err1asymp = mean(asymp < 0.05), err1perm = mean(perm < 0.05))
    
  giveDirection <- function(x) {
    jk <- binom::binom.confint(250, 5000)
    if(x > jk[2, 5:6]$lower & x < jk[5, 5:6]$upper) {
      return("\\cellcolor{green!25}")
    }
    if(x > jk[2, 5:6]$upper) {
      return("\\cellcolor{red!25}")
    }
    if(x < jk[2, 5:6]$lower) {
      return("\\cellcolor{blue!25}")
    }
  }
  redresH0 <- redresH0 %>%mutate(Asymp = paste0(err1asymp, sapply(err1asymp, giveDirection)),
                                 Perm = paste0(err1perm, sapply(err1perm, giveDirection)))
  
  redresH0$n_i <- as.factor(redresH0$n_i)
  levels(redresH0$n_i) <- c("$n_{\\text{bal}}$", "$n_{\\text{unbal}}$")
  redresH0$c_i <- as.factor(redresH0$c_i)
  levels(redresH0$c_i) <- c("$c_{\\text{low}}$", "$c_{\\text{medium}}$", "$c_{\\text{high}}$")
  redresH0$hyp <- as.factor(redresH0$hyp)
  levels(redresH0$hyp) <- c("$A$", "$AB$")
  redresH0$dist <- as.factor(redresH0$dist)
  levels(redresH0$dist) <- c("Exp", "Exp-pwExp", "LogN", "Weibull")
  redresH0$fac <- as.factor(redresH0$fac)
  
  redresH0 <- redresH0 %>% filter(hyp == h) %>% droplevels()
  jo <-
  tabular(
    Heading("Distribution") * dist* Heading("Design") * n_i * Heading("$q$") * 
      fac ~ Heading() * c_i * Heading() *identity * (Asymp + Perm),
    data = redresH0
  )
toLatex(jo)
}
# tableH0("$A$")
# tableH0("$AB$")

################################################################################
############################# H_1 simulation setup #############################
################################################################################
plotH1_setup <- function() {
  # pre-calculate tau for discontinuity points 
  p_tau <- optimize(function(x) (integrate(function(x) 1 - pweS(x, 1/6, 3, 12/35), lower = 0, upper = x)$value -
                                   integrate(pexp, lower = 0, upper = x)$value)^30,
                    lower = 0, upper = 5, tol = .Machine$double.eps)$minimum
  
  # calculate RMST for Exp(1)
  expRMST <- integrate(function(x) 1 - pexp(x), lower = 0, upper = p_tau)$value
  
  salti <- sapply(seq(0, 0.4, 0.05), function(i) 
    optimize(function(s) (i - abs(expRMST - pwRMST(s, p_tau)))^30, 
             lower = 1/6, upper = p_tau, tol = .Machine$double.eps)$minimum) 
  
  g <- ggplot() + 
  stat_function(fun = function(x) 1 - pexp(x, rate = 1), 
                n = 10000, aes(color = "$\\text{Exp}(1)$")) +
  xlim(c(0, 3)) + ylim(c(0, 1)) + xlab("$t$") + ylab("$S_i(t)$") + 
  stat_function(fun = function(x) pweS(x, salt = salti[1], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[2], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[3], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[4], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[5], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[6], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[7], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[8], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  stat_function(fun = function(x) pweS(x, salt = salti[9], h1 = 3, h2 = 12/35), 
                n = 10000, aes(color = paste("$\\text{pwExp}\\left(s, 3,\\frac{12}{35}\\right)$"))) +
  scale_colour_manual(name = "event time distribution", values = c(1, rep(2, 8))) + theme_bw() 

  g <- g + geom_segment(aes(x=1.478929,xend=1.478929,
            y=0,yend=pweS(1.478929, salt = 1/6, h1 = 3, h2 = 12/35), 
            linetype = "$\\tau=1.478929$")) + 
    scale_linetype_manual(name = "end point for RMST", values = 2)

  g
}


################################################################################
############################ H_1 simulation results ############################
################################################################################
tableH1 <- function() {
  flsi <- list.files()
  flsi_i <- grep("resH1_", flsi, value = TRUE)
  lapply(flsi_i, load, .GlobalEnv)
  simresH1 <- do.call(rbind, lapply(paste0("f_resH1_shifted_", 1:108, "_par"), get, envir = globalenv()))
  redresH1 <- simresH1 %>% 
    group_by(hyp, fac, n_i, c_i, saltus) %>%
    summarize(Asymp = mean(asymp < 0.05), Perm = mean(perm < 0.05))
  
  redresH1$n_i <- as.factor(redresH1$n_i)
  levels(redresH1$n_i) <- c("$n_{\\text{unbal}}$")
  redresH1$c_i <- as.factor(redresH1$c_i)
  levels(redresH1$c_i) <- c("$c_{\\text{low}}$", "$c_{\\text{medium}}$")
  redresH1$hyp <- as.factor(redresH1$hyp)
  levels(redresH1$hyp) <- c("$A$", "$AB$")
  redresH1$fac <- as.factor(redresH1$fac)
  redresH1$saltus <- as.factor(redresH1$saltus)
  levels(redresH1$saltus) <- seq(0, 0.4, 0.05)
  
  jo <-
    tabular(
      Heading("$q$") * fac * Heading("$\\delta$") * saltus ~ 
        Heading() * c_i * Heading() * identity * (Asymp + Perm),
      data = redresH1
    )
  toLatex(jo)
}
# tableH1()

plotH1_res <- function() {
  flsiH1 <- list.files()
  flsi_iH1 <- grep("resH1_shifted_", flsiH1, value = TRUE)
  lapply(flsi_iH1, load, .GlobalEnv)
  simresH1_shifted <- do.call(rbind, lapply(paste0("f_resH1_shifted_", 1:108, "_par"), 
                          get, envir = globalenv()))
  
  redresH1 <- simresH1_shifted %>% 
    group_by(hyp, fac, n_i, c_i, saltus) %>%
    summarize(powperm = mean(perm < 0.05), powasymp = mean(asymp < 0.05)) %>%
    pivot_longer(cols = c(powperm, powasymp), names_to = "testtype")
  
  redresH1$fac <- as.factor(redresH1$fac)
  levels(redresH1$fac) <- paste0("$q=", c(1,2,4), "$")
  redresH1$hyp <- as.factor(redresH1$hyp)
  levels(redresH1$hyp) <- c("$A$", "$AB$")
  
  redresH1$saltus <- as.factor(redresH1$saltus)
  levels(redresH1$saltus) <- as.factor(seq(0, 0.4, 0.05))
  redresH1$saltus <- as.numeric(as.character(redresH1$saltus))
  redresH1$testtype <- as.factor(redresH1$testtype)
  levels(redresH1$testtype) <- c("Asymp", "Perm")
  
  redresH1$c_i <- as.factor(redresH1$c_i)
  levels(redresH1$c_i) <- c("$c_{\\text{low}}$", "$c_{\\text{medium}}$")
  
  ggplot(redresH1, aes(y = value, x = saltus, col = testtype, linetype = c_i)) +
    geom_line() + geom_vline(xintercept = 0) + facet_grid(hyp ~ fac) + 
    geom_hline(aes(yintercept = 0.05)) + theme_bw() + 
    scale_color_manual(name = "Test", values = 1:2) +
    scale_linetype_manual(name = "Censoring rate", values = 1:2) +
    # xlab("$\\text{RMST}_{\\text{Exp}}(\\tau) - \\text{RMST}_{\\text{pwExp}_i}(\\tau)$") + 
    xlab("$\\Delta$") + ylab("Power") + ylim(c(0, 1))
}

################################################################################
############################## Real Data Example ###############################
################################################################################
ftab <- table(colonCS$sex, colonCS$rx)
row.names(ftab) <- c("\texttt{female}", "\texttt{male}")
colnames(ftab) <- c("\texttt{Obs}", "\texttt{Lev}", "\texttt{Lev+5FU}")
# print(xtable(ftab, digits = 2), include.rownames = T,
#       sanitize.colnames.function = identity, sanitize.rownames.function = identity)

in_fit <- survfit(Surv(Stime, event) ~ sex + rx, data = colonCS)
p_fit <- ggsurvplot_facet(in_fit, facet.by = "sex", data = colonCS, 
                          ggtheme = theme_bw(), censor = FALSE,
                          panel.labs = list(sex = c("female", "male"))) + 
  xlab("$t$") + ylab("$\\widehat{S}_i(t)$") + geom_vline(xintercept = 1500)


RMST_k <- matrix(colonCS %>% group_by(sex, rx) %>% 
                 summarise(RMST(Stime, event, 1500)[1]) %>% 
                 pull(), nrow = 2, byrow = 2)
row.names(RMST_k) <- c("\\texttt{female}", "\\texttt{male}")
colnames(RMST_k) <- c("\\texttt{Obs}", "\\texttt{Lev}", "\\texttt{Lev+5FU}")
# print(xtable(RMST_k, digits = 0), include.rownames = T,
#       sanitize.colnames.function = identity, sanitize.rownames.function = identity)

cens_k <- matrix(colonCS %>% group_by(sex, rx) %>% 
                   summarise(mean(!event)) %>% 
                   pull(), nrow = 2, byrow = 2)
row.names(cens_k) <- c("\\texttt{female}", "\\texttt{male}")
colnames(cens_k) <- c("\\texttt{Obs}", "\\texttt{Lev}", "\\texttt{Lev+5FU}")
# print(xtable(cens_k, digits = 2), include.rownames = T,
#       sanitize.colnames.function = identity, sanitize.rownames.function = identity)

# resA <- testRMST(time = "Stime", cens = "event", fA = "sex", fB = "rx",
#                  dat = colonCS, tau = 1500, alternative = "A", nperm = 1999,
#                  permDist = TRUE, perm = TRUE)
# resB <- testRMST(time = "Stime", cens = "event", fA = "sex", fB = "rx",
#                  dat = colonCS, tau = 1500, alternative = "B", nperm = 1999,
#                  permDist = TRUE, perm = TRUE)
# resAB <- testRMST(time = "Stime", cens = "event", fA = "sex", fB = "rx",
#                   dat = colonCS, tau = 1500, alternative = "AB", nperm = 1999,
#                   permDist = TRUE, perm = TRUE)

plotColon <- function() {
  
  pA <- ggplot(data = data.frame(x = resA$perm_dist), 
               aes(x = x, col = "Permutation")) + stat_ecdf() + 
    stat_function(fun = pchisq, args= list(df = resA$asymp_rank), 
                  aes(col = "$\\chi^2_{\\text{rank}(\\boldsymbol{T})}$")) + 
    theme_bw() + xlim(c(0, max(resA$perm_dist, resB$perm_dist, resAB$perm_dist))) +
  geom_segment(aes(x=resA$trueW, xend=resA$trueW, y = 0, yend= 1 - resA$asymp,
                   linetype = "observed $w(\\boldsymbol{T})$"), size = 1) +
    xlab("$W_n(\\boldsymbol{T})$") + ylab("$F_{W_n(\\boldsymbol{T})}(\\cdot)$") +
    scale_color_manual(name = "Distribution", values = c("red", "black")) +
    scale_linetype_manual(name = "for $\\boldsymbol{H}$", values = 2)
  
  pB <- ggplot(data = data.frame(x = resB$perm_dist), 
               aes(x = x, col = "Permutation")) + stat_ecdf() + 
    stat_function(fun = pchisq, args= list(df = resB$asymp_rank), 
                  aes(col = "$\\chi^2_{\\text{rank}(\\boldsymbol{T})}$")) + 
    theme_bw() + xlim(c(0, max(resA$perm_dist, resB$perm_dist, resAB$perm_dist))) +
    geom_segment(aes(x=resB$trueW, xend=resB$trueW, y = 0, yend= 1 - resB$asymp,
                     linetype = "observed $w(\\boldsymbol{T})$"), size = 1) +
    xlab("$W_n(\\boldsymbol{T})$") + ylab("$F_{W_n(\\boldsymbol{T})}(\\cdot)$") +
    scale_color_manual(name = "Distribution", values = c("red", "black")) +
    scale_linetype_manual(name = "for $\\boldsymbol{H}$", values = 2)
  
  pAB <- ggplot(data = data.frame(x = resAB$perm_dist), 
                aes(x = x, col = "Permutation")) + stat_ecdf() + 
    stat_function(fun = pchisq, args= list(df = resAB$asymp_rank), 
                  aes(col = "$\\chi^2_{\\text{rank}(\\boldsymbol{T})}$")) + 
    theme_bw() + xlim(c(0, max(resA$perm_dist, resB$perm_dist, resAB$perm_dist))) +
    geom_segment(aes(x=resAB$trueW, xend=resAB$trueW, y = 0, yend= 1 - resAB$asymp,
                     linetype = "observed $w(\\boldsymbol{T})$"), size = 1) +
    xlab("$W_n(\\boldsymbol{T})$") + ylab("$F_{W_n(\\boldsymbol{T})}(\\cdot)$") +
    scale_color_manual(name = "Distribution", values = c("red", "black")) +
    scale_linetype_manual(name = "for $\\boldsymbol{H}$", values = 2)

  jo <- get_legend(pA)
  pA <- pA + theme(legend.position = "none")
  pB <- pB + theme(legend.position = "none")
  pAB <- pAB + theme(legend.position = "none")
  
  ggarrange(pA, pB, pAB, jo, ncol = 2, nrow = 2)
}

# plab <- c("$\\mathcal{H}_0(\\bm{H}_{\\texttt{sex}})$",
#           "$\\mathcal{H}_0(\\bm{H}_{\\texttt{rx}})$",
#           "$\\mathcal{H}_0(\\bm{H}_{\\texttt{sex, rx}})$")
# ptab <- matrix(paste0(
#   round(c(resA$perm, resB$perm, resAB$perm,
#           resA$asymp, resB$asymp, resAB$asymp), digits = 3), " (",
#   round(c(p.adjust(c(resA$perm, resB$perm, resAB$perm), method = "holm"),
#           p.adjust(c(resA$asymp, resB$asymp, resAB$asymp), method = "holm")), digits = 3), ")"),
# nrow = 2, byrow = TRUE)
# row.names(ptab) <- c("Perm", "Asymp")
# colnames(ptab) <- plab
# print(xtable(ptab, digits = 4), include.rownames = T,
#       sanitize.colnames.function = identity, sanitize.rownames.function = identity)
