library(rdist)
library(MASS)
library(np)
library(ds4psy)
library(ggplot2)
# library(sphereplot)
world <- map_data("world")
library(mapproj)
library(rotasym)
library( dgof )
library(fields)

source("util_funcs.R")
source("test_function.R")
source("ico_syms.R")


set.seed(1010)
data("sunspots_births")
sunspots_births$X <-
  cbind(cos(sunspots_births$phi) * cos(sunspots_births$theta),
        cos(sunspots_births$phi) * sin(sunspots_births$theta),
        sin(sunspots_births$phi))

# head(sunspots_births$X)

p_vals <- matrix(0, nrow = 13, ncol = 16)

set.seed(1011)

num_samps <- 250
a_B <- 1000
a_lambda <- -1


for (i in 1:13) {
  m_cycle <- i + 11
  inds <- which( sunspots_births$cycle == m_cycle, arr.ind = TRUE )
  all_data <- cbind( sunspots_births$X[inds, ], sunspots_births$date[inds] )
  m_inds <- sample( dim(all_data)[1], num_samps, replace = FALSE )
  m_data <- all_data[m_inds, ]
  
  m_data <- cbind( m_data, local_average( m_data[, -4], all_data, 0.2 ) )

  p_vals[i, 1] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_1, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 2] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_2, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 3] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_3, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 4] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_4, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 5] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_5, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 6] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_6, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 7] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_7, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 8] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_8, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 9] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_9, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 10] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_10, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 11] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_11, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 12] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_12, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 13] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_13, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 14] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_14, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 15] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_15, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  p_vals[i, 16] <- PermVarTest( m_data[,1:3] , m_data[,5], sample_so3, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
  
  print( paste("Solar Cycle: ", m_cycle) )
  print( p_vals[i, ] )

}

alpha <- 0.02/16
p_vals >= alpha

#### Table 2
rowSums( p_vals >= alpha )


# 
# 
# #### Testing Zone
# m_cycle <- 14
# inds <- which( sunspots_births$cycle == m_cycle, arr.ind = TRUE )
# all_data <- cbind( sunspots_births$X[inds, ], sunspots_births$date[inds] )
# m_inds <- sample( dim(all_data)[1], num_samps, replace = FALSE )
# m_data <- all_data[m_inds, ]
# 
# m_data <- cbind( m_data, local_average( m_data[, -4], all_data, 0.2 ) )
# 
# 
# PermVarTest( m_data[,1:3] , m_data[,5], sample_ico_unif_axis_1, a_B, g_dot, plot_flag = TRUE, m_lambda = a_lambda, d_X = "Great Circle")
# 
# 

