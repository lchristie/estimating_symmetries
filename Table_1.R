library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("test_function.R")
source("Example_6.1_Helpers.R")

alpha <- 0.05
ns <- c(20,30,40,50,60,70,80,90,100,125,150,200,300,400,500,750)
num_sims <- 200




avt.file.name_2 <- paste0( "sims/sim_results_avt_f", 2, ".csv" )
pv.file.name_2 <- paste0( "sims/sim_results_pv_f", 2, ".csv" )
pv_lambda.file.name_2 <- paste0( "sims/sim_results_pv_lambda_f", 2, ".csv" )

avt.file.name_4 <- paste0( "sims/sim_results_avt_f", 4, ".csv" )
pv.file.name_4 <- paste0( "sims/sim_results_pv_f", 4, ".csv" )
pv_lambda.file.name_4 <- paste0( "sims/sim_results_pv_lambda_f", 4, ".csv" )

avt.file.name_6 <- paste0( "sims/sim_results_avt_f", 6, ".csv" )
pv.file.name_6 <- paste0( "sims/sim_results_pv_f", 6, ".csv" )
pv_lambda.file.name_6 <- paste0( "sims/sim_results_pv_lambda_f", 6, ".csv" )


#### Plots for f_d

## Read Sims
rejections_fd_avt_2 <- as.matrix( read.csv( avt.file.name_2 , header = TRUE) )[,-1]
names( rejections_fd_avt_2 ) <- NULL
rownames( rejections_fd_avt_2 ) <- NULL
colnames( rejections_fd_avt_2 ) <- NULL
rejections_fd_pv_2 <- as.matrix( read.csv( pv.file.name_2 , header = TRUE ) )[,-1]
names( rejections_fd_pv_2 ) <- NULL
rownames( rejections_fd_pv_2 ) <- NULL
colnames( rejections_fd_pv_2 ) <- NULL
rejections_fd_pv_lambda_2 <- as.matrix( read.csv( pv_lambda.file.name_2 , header = TRUE ) )[,-1]
names( rejections_fd_pv_2 ) <- NULL
rownames( rejections_fd_pv_2 ) <- NULL
colnames( rejections_fd_pv_2 ) <- NULL


rejections_fd_avt_4 <- as.matrix( read.csv( avt.file.name_4 , header = TRUE) )[,-1]
names( rejections_fd_avt_4 ) <- NULL
rownames( rejections_fd_avt_4 ) <- NULL
colnames( rejections_fd_avt_4 ) <- NULL
rejections_fd_pv_4 <- as.matrix( read.csv( pv.file.name_4 , header = TRUE ) )[,-1]
names( rejections_fd_pv_4 ) <- NULL
rownames( rejections_fd_pv_4 ) <- NULL
colnames( rejections_fd_pv_4 ) <- NULL
rejections_fd_pv_lambda_4 <- as.matrix( read.csv( pv_lambda.file.name_4 , header = TRUE ) )[,-1]
names( rejections_fd_pv_lambda_4 ) <- NULL
rownames( rejections_fd_pv_lambda_4 ) <- NULL
colnames( rejections_fd_pv_lambda_4 ) <- NULL


rejections_fd_avt_6 <- as.matrix( read.csv( avt.file.name_6 , header = TRUE) )[,-1]
names( rejections_fd_avt_6 ) <- NULL
rownames( rejections_fd_avt_6 ) <- NULL
colnames( rejections_fd_avt_6 ) <- NULL
rejections_fd_pv_6 <- as.matrix( read.csv( pv.file.name_6 , header = TRUE ) )[,-1]
names( rejections_fd_pv_6 ) <- NULL
rownames( rejections_fd_pv_6 ) <- NULL
colnames( rejections_fd_pv_6 ) <- NULL
rejections_fd_pv_lambda_6 <- as.matrix( read.csv( pv_lambda.file.name_6 , header = TRUE ) )[,-1]
names( rejections_fd_pv_lambda_6 ) <- NULL
rownames( rejections_fd_pv_lambda_6 ) <- NULL
colnames( rejections_fd_pv_lambda_6 ) <- NULL


k_hat_errs_avt_2 <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_2 <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_lambda_2 <- matrix( 0, ncol = length(ns), nrow = num_sims )

for (k in 1:length(ns)) {
  for (i in 1:num_sims) {
    ind <- (k-1)*num_sims + i
    k_hat_errs_avt_2[i, k] <- K_hat_sym_diff( rejections_fd_avt_2[,ind] > alpha )
    k_hat_errs_pv_2[i, k] <- K_hat_sym_diff( rejections_fd_pv_2[,ind] > alpha )
    k_hat_errs_pv_lambda_2[i, k] <- K_hat_sym_diff( rejections_fd_pv_lambda_2[,ind] > alpha )
  }
}

k_hat_errs_avt_4 <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_4 <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_lambda_4 <- matrix( 0, ncol = length(ns), nrow = num_sims )

for (k in 1:length(ns)) {
  for (i in 1:num_sims) {
    ind <- (k-1)*num_sims + i
    k_hat_errs_avt_4[i, k] <- K_hat_sym_diff( rejections_fd_avt_4[,ind] > alpha )
    k_hat_errs_pv_4[i, k] <- K_hat_sym_diff( rejections_fd_pv_4[,ind] > alpha )
    k_hat_errs_pv_lambda_4[i, k] <- K_hat_sym_diff( rejections_fd_pv_lambda_4[,ind] > alpha )
  }
}

k_hat_errs_avt_6 <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_6 <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv_lambda_6 <- matrix( 0, ncol = length(ns), nrow = num_sims )

for (k in 1:length(ns)) {
  for (i in 1:num_sims) {
    ind <- (k-1)*num_sims + i
    k_hat_errs_avt_6[i, k] <- K_hat_sym_diff( rejections_fd_avt_6[,ind] > alpha )
    k_hat_errs_pv_6[i, k] <- K_hat_sym_diff( rejections_fd_pv_6[,ind] > alpha )
    k_hat_errs_pv_lambda_6[i, k] <- K_hat_sym_diff( rejections_fd_pv_lambda_6[,ind] > alpha )
  }
}

## Proportions for Table 5.1 From Example 6.1

Table_1 <- matrix( ncol = 8 )
colnames(Table_1) =  c( "Example", "Scenario", "Test", "20", "50", "100", "200", "500")
  
Table_1[1,] <- c( "6.1", "Dimension 2", "AVT", colMeans( k_hat_errs_avt_2 == 0 )[c(1,4,9,12,15)] )
Table_1 <- rbind( Table_1, c( "6.1", "Dimension 2", "PV", colMeans( k_hat_errs_pv_2 == 0 )[c(1,4,9,12,15)] ) )
Table_1 <- rbind( Table_1, c( "6.1", "Dimension 2", "PV, \\hat{\\lambda}", colMeans( k_hat_errs_pv_lambda_2 == 0 )[c(1,4,9,12,15)] ) )

Table_1 <- rbind( Table_1, c( "6.1", "Dimension 4", "AVT", colMeans( k_hat_errs_avt_4 == 0 )[c(1,4,9,12,15)] ) )
Table_1 <- rbind( Table_1, c( "6.1", "Dimension 4", "PV", colMeans( k_hat_errs_pv_4 == 0 )[c(1,4,9,12,15)] ) ) 
Table_1 <- rbind( Table_1, c( "6.1", "Dimension 4", "PV, \\hat{\\lambda}", colMeans( k_hat_errs_pv_lambda_4 == 0 )[c(1,4,9,12,15)] ) )

Table_1 <- rbind( Table_1, c( "6.1", "Dimension 6", "AVT", colMeans( k_hat_errs_avt_6 == 0 )[c(1,4,9,12,15)] ) )
Table_1 <- rbind( Table_1, c( "6.1", "Dimension 6", "PV", colMeans( k_hat_errs_pv_6 == 0 )[c(1,4,9,12,15)] ) ) 
Table_1 <- rbind( Table_1, c( "6.1", "Dimension 6", "PV, \\hat{\\lambda}", colMeans( k_hat_errs_pv_lambda_6 == 0 )[c(1,4,9,12,15)] ) )


#### Example 6.2 Rows for Table 2


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



#### Proportion Counts


inds <- c(2, 5, 7, 9)
Table_1 <- rbind( Table_1, c( "6.2", "Scenario 1", "PV", colMeans( scenario_1_perm[, inds] == 16), "-" ) )
Table_1 <- rbind( Table_1, c( "6.2", "Scenario 2", "PV", colMeans( scenario_2_perm[, inds] == 16), "-" ) )
Table_1 <- rbind( Table_1, c( "6.2", "Scenario 3", "PV", colMeans( scenario_3_perm[, inds] == 16), "-" ) )
Table_1 <- rbind( Table_1, c( "6.2", "Scenario 4", "PV", colMeans( scenario_4_perm[, inds] == 16), "-" ) )


Table_1

