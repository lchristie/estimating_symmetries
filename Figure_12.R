library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("ico_syms.R")
source("test_function.R")


scenario_2_perm_varying_K <- read.csv("sims/scenario_2_LCE_perms_varying_K_2.csv")[,-1]
scenario_4_perm_varying_K <- read.csv("sims/scenario_4_LCE_perms_varying_K_2.csv")[,-1]


num_sims <- 100
num_perms <- 100
alpha <- 0.05
ns <- c(10, 20, 30, 40, 50, 75, 100, 150, 200)
mu_X <- rep(0, 3)
sigma <- 0.01


par(mfrow=c(1,2))

m <- length(ns)
MSPE_baseline <- scenario_2_perm_varying_K[,(4*m+1):(5*m)]
MSPE_full <- scenario_2_perm_varying_K[,(5*m+1):(6*m)]
MSPE_smaller_with_x <- scenario_2_perm_varying_K[,(6*m+1):(7*m)]
MSPE_smaller <- scenario_2_perm_varying_K[,(7*m+1):(8*m)]
MSPE_smallest <- scenario_2_perm_varying_K[,(8*m+1):(9*m)]

plot( ns, colMeans(MSPE_baseline), type = "l", ylim = c(-sigma, 0.7), 
      ylab = "MSPE", xlab = "Sample Size", lty = 1)
lines( ns - 1, colMeans(MSPE_full), col = "blue", lty = 1)
lines( ns - 2, colMeans(MSPE_smaller_with_x), col = "red", lty = 1)
lines( ns - 2, colMeans(MSPE_smaller), col = "green", lty = 1)
lines( ns - 2, colMeans(MSPE_smallest), col = "purple", lty = 1)

MSPE_baseline.sd <- 1.96 * apply(MSPE_baseline, 2, sd) / sqrt( num_sims )
MSPE_full.sd <-  1.96 * apply(MSPE_full, 2, sd) / sqrt( num_sims )
MSPE_smaller_with_x.sd <-  1.96 * apply(MSPE_smaller_with_x, 2, sd) / sqrt( num_sims )
MSPE_smaller.sd <-  1.96 * apply(MSPE_smaller, 2, sd) / sqrt( num_sims )
MSPE_smallest.sd <-  1.96 * apply(MSPE_smallest, 2, sd) / sqrt( num_sims )

arrows(x0=ns, y0= colMeans(MSPE_baseline) - MSPE_baseline.sd - 0.01, x1=ns, y1=colMeans(MSPE_baseline) + MSPE_baseline.sd + 0.01, code=3, angle=90, length=0.05)
arrows(x0=ns - 1, y0= colMeans(MSPE_full) - MSPE_full.sd - 0.01, x1=ns - 1, y1=colMeans(MSPE_full) + MSPE_full.sd + 0.01, code=3, angle=90, length=0.05, col = "blue")
arrows(x0=ns - 2, y0= colMeans(MSPE_smaller_with_x) - MSPE_smaller_with_x.sd - 0.01, x1=ns - 2, y1=colMeans(MSPE_smaller_with_x) + MSPE_smaller_with_x.sd + 0.01, code=3, angle=90, length=0.05, col = "red")
arrows(x0=ns - 2, y0= colMeans(MSPE_smaller) - MSPE_smaller.sd - 0.01, x1=ns - 2, y1=colMeans(MSPE_smaller) + MSPE_smaller.sd + 0.01, code=3, angle=90, length=0.05, col = "green")
arrows(x0=ns - 3, y0= colMeans(MSPE_smallest) - MSPE_smallest.sd - 0.01, x1=ns - 3, y1=colMeans(MSPE_smallest) + MSPE_smallest.sd + 0.01, code=3, angle=90, length=0.05, col = "purple")


m <- length(ns)
MSPE_baseline <- scenario_4_perm_varying_K[,(4*m+1):(5*m)]
MSPE_full <- scenario_4_perm_varying_K[,(5*m+1):(6*m)]
MSPE_smaller_with_x <- scenario_4_perm_varying_K[,(6*m+1):(7*m)]
MSPE_smaller <- scenario_4_perm_varying_K[,(7*m+1):(8*m)]
MSPE_smallest <- scenario_4_perm_varying_K[,(8*m+1):(9*m)]

plot( ns, colMeans(MSPE_baseline), type = "l", ylim = c(-sigma, 0.7), 
      ylab = "MSPE", xlab = "Sample Size", lty = 1)
lines( ns - 1, colMeans(MSPE_full), col = "blue", lty = 1)
lines( ns - 2, colMeans(MSPE_smaller_with_x), col = "red", lty = 1)
lines( ns - 2, colMeans(MSPE_smaller), col = "green", lty = 1)
lines( ns - 2, colMeans(MSPE_smallest), col = "purple", lty = 1)

MSPE_baseline.sd <- 1.96 * apply(MSPE_baseline, 2, sd) / sqrt( num_sims )
MSPE_full.sd <-  1.96 * apply(MSPE_full, 2, sd) / sqrt( num_sims )
MSPE_smaller_with_x.sd <-  1.96 * apply(MSPE_smaller_with_x, 2, sd) / sqrt( num_sims )
MSPE_smaller.sd <-  1.96 * apply(MSPE_smaller, 2, sd) / sqrt( num_sims )
MSPE_smallest.sd <-  1.96 * apply(MSPE_smallest, 2, sd) / sqrt( num_sims )

arrows(x0=ns, y0= colMeans(MSPE_baseline) - MSPE_baseline.sd - 0.01, x1=ns, y1=colMeans(MSPE_baseline) + MSPE_baseline.sd + 0.01, code=3, angle=90, length=0.05)
arrows(x0=ns - 1, y0= colMeans(MSPE_full) - MSPE_full.sd - 0.01, x1=ns - 1, y1=colMeans(MSPE_full) + MSPE_full.sd + 0.01, code=3, angle=90, length=0.05, col = "blue")
arrows(x0=ns - 2, y0= colMeans(MSPE_smaller_with_x) - MSPE_smaller_with_x.sd - 0.01, x1=ns - 2, y1=colMeans(MSPE_smaller_with_x) + MSPE_smaller_with_x.sd + 0.01, code=3, angle=90, length=0.05, col = "red")
arrows(x0=ns - 2, y0= colMeans(MSPE_smaller) - MSPE_smaller.sd - 0.01, x1=ns - 2, y1=colMeans(MSPE_smaller) + MSPE_smaller.sd + 0.01, code=3, angle=90, length=0.05, col = "green")
arrows(x0=ns - 3, y0= colMeans(MSPE_smallest) - MSPE_smallest.sd - 0.01, x1=ns - 3, y1=colMeans(MSPE_smallest) + MSPE_smallest.sd + 0.01, code=3, angle=90, length=0.05, col = "purple")

