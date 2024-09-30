library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("test_function.R")
source("Example_6.1_Helpers.R")


dimension <- 2
alpha <- 0.05
ns <- c(20,30,40,50,60,70,80,90,100,125,150,200,300,400,500,750)

mu_X <- rep(0, dimension)
Sigma_X = diag(rep(2, dimension))

sigma <- 0.05

a_B <- 100

num_sims <- 200

avt.file.name <- paste0( "sims/sim_results_avt_f", dimension, ".csv" )
pv.file.name <- paste0( "sims/sim_results_pv_f", dimension, ".csv" )
pv_lambda.file.name <- paste0( "sims/sim_results_pv_lambda_f", dimension, ".csv" )

#### Plots for f_d

## Read Sims
rejections_fd_avt <- as.matrix( read.csv( avt.file.name , header = TRUE) )[,-1]
names( rejections_fd_avt ) <- NULL
rownames( rejections_fd_avt ) <- NULL
colnames( rejections_fd_avt ) <- NULL
rejections_fd_pv <- as.matrix( read.csv( pv.file.name , header = TRUE ) )[,-1]
names( rejections_fd_pv ) <- NULL
rownames( rejections_fd_pv ) <- NULL
colnames( rejections_fd_pv ) <- NULL
rejections_fd_pv_lambda <- as.matrix( read.csv( pv_lambda.file.name , header = TRUE ) )[,-1]
names( rejections_fd_pv ) <- NULL
rownames( rejections_fd_pv ) <- NULL
colnames( rejections_fd_pv ) <- NULL

## Confidence Region Plots

k_hat_errs_avt <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_avt <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_pv <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_lambda <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_pv_lambda <- matrix( 0, ncol = length(ns), nrow = num_sims )

for (k in 1:length(ns)) {
  for (i in 1:num_sims) {
    ind <- (k-1)*num_sims + i
    k_hat_errs_avt[i, k] <- K_hat_sym_diff( rejections_fd_avt[,ind] > alpha )
    k_hat_errs_pv[i, k] <- K_hat_sym_diff( rejections_fd_pv[,ind] > alpha )
    k_hat_errs_pv_lambda[i, k] <- K_hat_sym_diff( rejections_fd_pv_lambda[,ind] > alpha )
    k_hat_props_contain_g_max_avt[i, k] <- ( sum( ( rejections_fd_avt[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 )
    k_hat_props_contain_g_max_pv[i, k] <- ( sum( ( rejections_fd_pv[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 ) 
    k_hat_props_contain_g_max_pv_lambda[i, k] <- ( sum( ( rejections_fd_pv_lambda[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 ) 
  }
}

plot( ns, colMeans(k_hat_errs_avt), pch = 15, type = "p", ylim = c(0,6), ylab = "Confidence Region Excess", xlab = "Sample Size (n)", log = "x")
points( ns-1, colMeans(k_hat_errs_pv), pch = 17, col = "red")
points( ns-2, colMeans(k_hat_errs_pv_lambda), pch = 19, col = "blue")
lines( ns, colMeans(k_hat_errs_avt), lty = 2)
lines( ns-1, colMeans(k_hat_errs_pv), lty = 4, col = "red")
lines( ns-2, colMeans(k_hat_errs_pv_lambda), lty = 3, col = "blue")


avt.sd <- 1.96 * apply(k_hat_errs_avt, 2, sd) / sqrt( num_sims )
pv.sd <- 1.96 * apply(k_hat_errs_pv, 2, sd) / sqrt( num_sims )
pv_lambda.sd <- 1.96 * apply(k_hat_errs_pv_lambda, 2, sd) / sqrt( num_sims )

arrows(x0=ns, y0= colMeans(k_hat_errs_avt) - avt.sd - 0.01, x1=ns, y1=colMeans(k_hat_errs_avt) + avt.sd + 0.01, code=3, angle=90, length=0.1)
arrows(x0=ns-1, y0= colMeans(k_hat_errs_pv) - pv.sd - 0.01, x1=ns-1, y1=colMeans(k_hat_errs_pv) + pv.sd + 0.01, code=3, angle=90, length=0.1, col = "red")
arrows(x0=ns-2, y0= colMeans(k_hat_errs_pv_lambda) - pv_lambda.sd - 0.01, x1=ns-2, y1=colMeans(k_hat_errs_pv_lambda) + pv_lambda.sd + 0.01, code=3, angle=90, length=0.1, col = "blue")


## Proportions for Table 5.1

colMeans( k_hat_errs_avt == 0 )[c(1,4,9,12,15)]
colMeans( k_hat_errs_pv == 0 )[c(1,4,9,12,15)]
colMeans( k_hat_errs_pv_lambda == 0 )[c(1,4,9,12,15)]

