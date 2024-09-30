library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("test_function.R")

dimension <- 20
alpha <- 0.05
ns <- c(20,50,100,200,500) #,300,400,500,750)

mu_X <- rep(0, dimension)
Sigma_X = diag(rep(2, dimension))

sigma <- 0.05

a_B <- 200

num_sims <- 100

avt.file.name <- paste0( "sims/sim_results_avt_f", dimension, "_v2.csv" )
pvt_no_proj.file.name <- paste0( "sims/sim_results_pv_f", dimension, "_no_proj_v2.csv" )
pvt_no_proj_lambda.file.name <-  paste0( "sims/sim_results_pv_f", dimension, "_no_proj_lambda_v2.csv" )
pv_with_proj.file.name <- paste0( "sims/sim_results_pv_f", dimension, "_with_proj_dim_2_v2.csv" )

## Read Sims
rejections_fd_avt <- as.matrix( read.csv( avt.file.name , header = TRUE) )[,-1]
names( rejections_fd_avt ) <- NULL
rownames( rejections_fd_avt ) <- NULL
colnames( rejections_fd_avt ) <- NULL
rejections_fd_pv_no_proj <- as.matrix( read.csv( pvt_no_proj.file.name , header = TRUE ) )[,-1]
names( rejections_fd_pv ) <- NULL
rownames( rejections_fd_pv ) <- NULL
colnames( rejections_fd_pv ) <- NULL
rejections_fd_pvt_no_proj_lambda <- as.matrix( read.csv( pvt_no_proj_lambda.file.name , header = TRUE) )[,-1]
names( rejections_fd_avt ) <- NULL
rownames( rejections_fd_avt ) <- NULL
colnames( rejections_fd_avt ) <- NULL
rejections_fd_pv_with_proj <- as.matrix( read.csv( pv_with_proj.file.name , header = TRUE ) )[,-1]
names( rejections_fd_pv ) <- NULL
rownames( rejections_fd_pv ) <- NULL
colnames( rejections_fd_pv ) <- NULL

## Confidence Region Plots

par(mfrow = c(1,1))

k_hat_errs_avt <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_avt <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_no_proj <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_pv_no_proj <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_no_proj_lambda <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_pv_no_proj_lambda <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_with_proj <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_pv_with_proj <- matrix( 0, ncol = length(ns), nrow = num_sims )



for (k in 1:length(ns)) {
  for (i in 1:num_sims) {
    ind <- (k-1)*num_sims + i
    k_hat_errs_avt[i, k] <- K_hat_sym_diff( rejections_fd_avt[,ind] > alpha )
    k_hat_errs_pv_no_proj[i, k] <- K_hat_sym_diff( rejections_fd_pv_no_proj[,ind] > alpha )
    k_hat_errs_pv_no_proj_lambda[i, k] <- K_hat_sym_diff( rejections_fd_pvt_no_proj_lambda[,ind] > alpha )
    k_hat_errs_pv_with_proj[i, k] <- K_hat_sym_diff( rejections_fd_pv_with_proj[,ind] > alpha )
    k_hat_props_contain_g_max_avt[i, k] <- ( sum( ( rejections_fd_avt[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 )
    k_hat_props_contain_g_max_pv_no_proj[i, k] <- ( sum( ( rejections_fd_pv_no_proj[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 ) 
    k_hat_props_contain_g_max_pv_no_proj_lambda[i, k] <- ( sum( ( rejections_fd_pvt_no_proj_lambda[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 ) 
    k_hat_props_contain_g_max_pv_with_proj[i, k] <- ( sum( ( rejections_fd_pv_with_proj[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 ) 
  }
}

plot( ns, colMeans(k_hat_errs_avt), pch = 15, type = "p", ylim = c(0,6), ylab = "Confidence Region Excess", xlab = "Sample Size (n)", log = "x")
points( ns, colMeans(k_hat_errs_pv_no_proj), pch = 17, col = "red")
points( ns, colMeans(k_hat_errs_pv_no_proj_lambda), pch = 18, col = "blue")
points( ns, colMeans(k_hat_errs_pv_with_proj), pch = 19, col = "green")

lines( ns, colMeans(k_hat_errs_avt), lty = 2)
lines( ns, colMeans(k_hat_errs_pv_no_proj), lty = 2, col = "red")
lines( ns, colMeans(k_hat_errs_pv_no_proj_lambda), lty = 2, col = "blue")
lines( ns, colMeans(k_hat_errs_pv_with_proj), lty = 2, col = "green")

avt.sd <- 1.96 * apply(k_hat_errs_avt, 2, sd) / sqrt( num_sims )
pv_no_proj.sd <- 1.96 * apply(k_hat_errs_pv_no_proj, 2, sd) / sqrt( num_sims )
pv_no_proj_lambda.sd <- 1.96 * apply(k_hat_errs_pv_no_proj_lambda, 2, sd) / sqrt( num_sims )
pv_with_proj.sd <- 1.96 * apply(k_hat_errs_pv_with_proj, 2, sd) / sqrt( num_sims )

arrows(x0=ns, y0= colMeans(k_hat_errs_avt) - avt.sd - 0.01, x1=ns, y1=colMeans(k_hat_errs_avt) + avt.sd + 0.01, code=3, angle=90, length=0.1)
arrows(x0=ns, y0= colMeans(k_hat_errs_pv_no_proj) - pv_no_proj.sd - 0.01, x1=ns, y1=colMeans(k_hat_errs_pv_no_proj) + pv_no_proj.sd + 0.01, code=3, angle=90, length=0.1, col = "red")
arrows(x0=ns, y0= colMeans(k_hat_errs_pv_no_proj_lambda) - pv_no_proj_lambda.sd - 0.01, x1=ns, y1=colMeans(k_hat_errs_pv_no_proj_lambda) + pv_no_proj_lambda.sd + 0.01, code=3, angle=90, length=0.1, col = "blue")
arrows(x0=ns, y0= colMeans(k_hat_errs_pv_with_proj) - pv_with_proj.sd - 0.01, x1=ns, y1=colMeans(k_hat_errs_pv_with_proj) + pv_with_proj.sd + 0.01, code=3, angle=90, length=0.1, col = "green")


