library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("ico_syms.R")
source("test_function.R")



num_sims <- 100
num_perms <- 100
alpha <- 0.05
ns <- c(10, 20, 30, 40, 50, 75, 100, 150, 200)
mu_X <- rep(0, 3)
sigma <- 0.01

scenario_1_perm <- read.csv( "sims/scenario_1_LCE_perms.csv", header = TRUE )[,-1]
scenario_2_perm <- read.csv( "sims/scenario_2_LCE_perms.csv", header = TRUE )[,-1]
scenario_3_perm <- read.csv( "sims/scenario_3_LCE_perms.csv", header = TRUE )[,-1]
scenario_4_perm <- read.csv( "sims/scenario_4_LCE_perms.csv", header = TRUE )[,-1]

#### Plots

scenario <- scenario_1_perm

m <- length(ns)
MSPE_baseline <- scenario[,(2*m+1):(3*m)]
MSPE_full <- scenario[,(3*m+1):(4*m)]
MSPE_split <- scenario[,(4*m+1):(5*m)]

plot( ns, colMeans(MSPE_baseline), type = "l", ylim = c(-sigma, 0.7), 
      ylab = "MSPE", xlab = "Sample Size",  lty = 2)
lines( ns - 1, colMeans(MSPE_full), col = "blue", lty = 2)
lines( ns - 2, colMeans(MSPE_split), col = "red", lty = 2)

points( ns, colMeans(MSPE_baseline), pch = 15)
points( ns - 1, colMeans(MSPE_full), col = "blue", pch = 16)
points( ns - 2, colMeans(MSPE_split), col = "red", pch = 17)

MSPE_baseline.sd <- 1.96 * apply(MSPE_baseline, 2, sd) / sqrt( num_sims )
MSPE_full.sd <-  1.96 * apply(MSPE_full, 2, sd) / sqrt( num_sims )
MSPE_split.sd <-  1.96 * apply(MSPE_split, 2, sd) / sqrt( num_sims )

arrows(x0=ns, y0= colMeans(MSPE_baseline) - MSPE_baseline.sd - 0.01, x1=ns, y1=colMeans(MSPE_baseline) + MSPE_baseline.sd + 0.01, code=3, angle=90, length=0.05)
arrows(x0=ns - 1, y0= colMeans(MSPE_full) - MSPE_full.sd - 0.01, x1=ns - 1, y1=colMeans(MSPE_full) + MSPE_full.sd + 0.01, code=3, angle=90, length=0.05, col = "blue")
arrows(x0=ns - 2, y0= colMeans(MSPE_split) - MSPE_split.sd - 0.01, x1=ns - 2, y1=colMeans(MSPE_split) + MSPE_split.sd + 0.01, code=3, angle=90, length=0.05, col = "red")



#### Proportion Counts


inds <- c(2, 5, 7, 9)
colMeans( scenario_1_perm[, inds] == 16)
colMeans( scenario_2_perm[, inds] == 13)
colMeans( scenario_3_perm[, inds] == 0)
colMeans( scenario_4_perm[, inds] == 0)

