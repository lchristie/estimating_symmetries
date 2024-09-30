library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("test_function.R")
source("Example_6.1_Helpers.R")

alpha <- 0.05
ns <- c(20,30,40,50,60,70,80,90,100,125,150,200,300,400,500,750)

mu_X <- rep(0, dimension)
Sigma_X = diag(rep(2, dimension))

sigma <- 0.05

a_B <- 100

num_sims <- 200


#### Dimension 2

dimension <- 2

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

## Power and Size Tables

rej_props_fd_avt <- matrix(0, nrow = 6, ncol = length(ns) )
rej_props_fd_pv <- matrix(0, nrow = 6, ncol = length(ns) )
rej_props_fd_pv_lambda <- matrix(0, nrow = 6, ncol = length(ns) )

for (k in 1:length(ns)) {
  start_ind <- (k-1)*num_sims + 1
  stop_ind <- (k)*num_sims
  rej_props_fd_avt[1,k] <- mean( rejections_fd_avt[1, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[2,k] <- mean( rejections_fd_avt[2, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[3,k] <- mean( rejections_fd_avt[3, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[4,k] <- mean( rejections_fd_avt[4, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[5,k] <- mean( rejections_fd_avt[5, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[6,k] <- mean( rejections_fd_avt[6, start_ind:stop_ind] < alpha )
  
  rej_props_fd_pv[1,k] <- mean( rejections_fd_pv[1, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[2,k] <- mean( rejections_fd_pv[2, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[3,k] <- mean( rejections_fd_pv[3, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[4,k] <- mean( rejections_fd_pv[4, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[5,k] <- mean( rejections_fd_pv[5, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[6,k] <- mean( rejections_fd_pv[6, start_ind:stop_ind] < alpha )
  
  rej_props_fd_pv_lambda[1,k] <- mean( rejections_fd_pv_lambda[1, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[2,k] <- mean( rejections_fd_pv_lambda[2, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[3,k] <- mean( rejections_fd_pv_lambda[3, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[4,k] <- mean( rejections_fd_pv_lambda[4, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[5,k] <- mean( rejections_fd_pv_lambda[5, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[6,k] <- mean( rejections_fd_pv_lambda[6, start_ind:stop_ind] < alpha )
}

## Used in Table 3 and 4
# rej_props_fd_avt[,c(1,4,9,12,15)]
# rej_props_fd_pv[,c(1,4,9,12,15)]
# rej_props_fd_pv_lambda[,c(1,4,9,12,15)]

Table_3 <- matrix( c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_avt[1,c(1,4,9,12,15)] ), ncol = 8 )
colnames( Table_3 ) <- c( "Dimension", "Hypothesis", "Group", "20", "50", "100", "200", "500")
Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_avt[2,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_avt[3,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_avt[4,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_avt[5,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_avt[6,c(1,4,9,12,15)] ))

Table_4 <- matrix( c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_pv[1,c(1,4,9,12,15)] ), ncol = 8 )
colnames( Table_4 ) <- c( "Dimension", "Hypothesis", "Group", "20", "50", "100", "200", "500")
Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_pv[2,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_pv[3,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_pv[4,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_pv[5,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_pv[6,c(1,4,9,12,15)] ))

Table_5 <- matrix( c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_pv_lambda[1,c(1,4,9,12,15)] ), ncol = 8)
colnames( Table_5 ) <- c( "Dimension", "Hypothesis", "Group", "20", "50", "100", "200", "500")
Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_pv_lambda[2,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_pv_lambda[3,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_pv_lambda[4,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_pv_lambda[5,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_pv_lambda[6,c(1,4,9,12,15)] ))

#### Dimension 4

dimension <- 4

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

## Power and Size Tables

rej_props_fd_avt <- matrix(0, nrow = 6, ncol = length(ns) )
rej_props_fd_pv <- matrix(0, nrow = 6, ncol = length(ns) )
rej_props_fd_pv_lambda <- matrix(0, nrow = 6, ncol = length(ns) )

for (k in 1:length(ns)) {
  start_ind <- (k-1)*num_sims + 1
  stop_ind <- (k)*num_sims
  rej_props_fd_avt[1,k] <- mean( rejections_fd_avt[1, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[2,k] <- mean( rejections_fd_avt[2, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[3,k] <- mean( rejections_fd_avt[3, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[4,k] <- mean( rejections_fd_avt[4, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[5,k] <- mean( rejections_fd_avt[5, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[6,k] <- mean( rejections_fd_avt[6, start_ind:stop_ind] < alpha )
  
  rej_props_fd_pv[1,k] <- mean( rejections_fd_pv[1, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[2,k] <- mean( rejections_fd_pv[2, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[3,k] <- mean( rejections_fd_pv[3, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[4,k] <- mean( rejections_fd_pv[4, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[5,k] <- mean( rejections_fd_pv[5, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[6,k] <- mean( rejections_fd_pv[6, start_ind:stop_ind] < alpha )
  
  rej_props_fd_pv_lambda[1,k] <- mean( rejections_fd_pv_lambda[1, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[2,k] <- mean( rejections_fd_pv_lambda[2, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[3,k] <- mean( rejections_fd_pv_lambda[3, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[4,k] <- mean( rejections_fd_pv_lambda[4, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[5,k] <- mean( rejections_fd_pv_lambda[5, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[6,k] <- mean( rejections_fd_pv_lambda[6, start_ind:stop_ind] < alpha )
}

## Used in Table 3 and 4
# rej_props_fd_avt[,c(1,4,9,12,15)]
# rej_props_fd_pv[,c(1,4,9,12,15)]
# rej_props_fd_pv_lambda[,c(1,4,9,12,15)]

Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_avt[1,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_avt[2,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_avt[3,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_avt[4,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_avt[5,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_avt[6,c(1,4,9,12,15)] ))

Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_pv[1,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_pv[2,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_pv[3,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_pv[4,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_pv[5,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_pv[6,c(1,4,9,12,15)] ))

Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_pv_lambda[1,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_pv_lambda[2,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_pv_lambda[3,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_pv_lambda[4,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_pv_lambda[5,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_pv_lambda[6,c(1,4,9,12,15)] ))



#### Dimension 6

dimension <- 6

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

## Power and Size Tables

rej_props_fd_avt <- matrix(0, nrow = 6, ncol = length(ns) )
rej_props_fd_pv <- matrix(0, nrow = 6, ncol = length(ns) )
rej_props_fd_pv_lambda <- matrix(0, nrow = 6, ncol = length(ns) )

for (k in 1:length(ns)) {
  start_ind <- (k-1)*num_sims + 1
  stop_ind <- (k)*num_sims
  rej_props_fd_avt[1,k] <- mean( rejections_fd_avt[1, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[2,k] <- mean( rejections_fd_avt[2, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[3,k] <- mean( rejections_fd_avt[3, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[4,k] <- mean( rejections_fd_avt[4, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[5,k] <- mean( rejections_fd_avt[5, start_ind:stop_ind] < alpha )
  rej_props_fd_avt[6,k] <- mean( rejections_fd_avt[6, start_ind:stop_ind] < alpha )
  
  rej_props_fd_pv[1,k] <- mean( rejections_fd_pv[1, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[2,k] <- mean( rejections_fd_pv[2, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[3,k] <- mean( rejections_fd_pv[3, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[4,k] <- mean( rejections_fd_pv[4, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[5,k] <- mean( rejections_fd_pv[5, start_ind:stop_ind] < alpha )
  rej_props_fd_pv[6,k] <- mean( rejections_fd_pv[6, start_ind:stop_ind] < alpha )
  
  rej_props_fd_pv_lambda[1,k] <- mean( rejections_fd_pv_lambda[1, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[2,k] <- mean( rejections_fd_pv_lambda[2, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[3,k] <- mean( rejections_fd_pv_lambda[3, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[4,k] <- mean( rejections_fd_pv_lambda[4, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[5,k] <- mean( rejections_fd_pv_lambda[5, start_ind:stop_ind] < alpha )
  rej_props_fd_pv_lambda[6,k] <- mean( rejections_fd_pv_lambda[6, start_ind:stop_ind] < alpha )
}

## Used in Table 3 and 4
# rej_props_fd_avt[,c(1,4,9,12,15)]
# rej_props_fd_pv[,c(1,4,9,12,15)]
# rej_props_fd_pv_lambda[,c(1,4,9,12,15)]

Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_avt[1,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_avt[2,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_avt[3,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_avt[4,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_avt[5,c(1,4,9,12,15)] ))
Table_3 <- rbind( Table_3, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_avt[6,c(1,4,9,12,15)] ))

Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_pv[1,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_pv[2,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_pv[3,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_pv[4,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_pv[5,c(1,4,9,12,15)] ))
Table_4 <- rbind( Table_4, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_pv[6,c(1,4,9,12,15)] ))

Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_h \\rangle", rej_props_fd_pv_lambda[1,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_v \\rangle", rej_props_fd_pv_lambda[2,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_0", "\\langle R_\\pi \\rangle", rej_props_fd_pv_lambda[3,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_/ \\rangle", rej_props_fd_pv_lambda[4,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_\\ \\rangle", rej_props_fd_pv_lambda[5,c(1,4,9,12,15)] ))
Table_5 <- rbind( Table_5, c( dimension, "H_1", "\\langle R_\\pi/2 \\rangle", rej_props_fd_pv_lambda[6,c(1,4,9,12,15)] ))


Table_3
Table_4
Table_5









