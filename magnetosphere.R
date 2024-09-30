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

#### Setup Pars

d_X <- function ( x_i, X ) {
  return( sqrt( colSums( ( t(X) - as.numeric(x_i) )^2 ) ) )
}

norm_Y <- function ( y_i, Y ) {
  return( sqrt( colSums( (y_i - t(Y) )^2 ) ) )
}

#### Symmetries of Magnetosphere Data

m_data <- read.csv("SWARM_DATA.csv", header = TRUE)

m_data <- cbind( m_data, cos( m_data[,"LAT"] * pi / 180 ) * cos (m_data[,"LONG"] * pi / 180 ) )
m_data <- cbind( m_data, cos( m_data[,"LAT"] * pi / 180 ) * sin (m_data[,"LONG"] * pi / 180 ) )
m_data <- cbind( m_data, sin( m_data[,"LAT"] * pi / 180 )  )

colnames(m_data) <- c( colnames(m_data)[1:(dim(m_data)[2] - 3)] , "x_coord", "y_coord", "z_coord" )

for (i in 1:dim(icosahedron)[1]) {
  m_data <- cbind( m_data, acos( icosahedron[i,1] * m_data["x_coord"] + 
                                   icosahedron[i,2] * m_data["y_coord"] + 
                                   icosahedron[i,3] * m_data["z_coord"] ) ) 
  colnames(m_data)[length(colnames(m_data))] <- paste0( "proj_coord", i) 
}

set.seed(2021)
num_samples <- 1000
inds <- sample( 1:(dim(m_data)[1] * 1), num_samples, replace = FALSE)

m_data_used <- m_data[inds,]

## Symmetry Tests

a_B <- 1000
alpha <- 0.05
a_lambda <- -1
set.seed(2021)

X <- as.matrix( cbind( m_data_used["x_coord"], m_data_used["y_coord"], m_data_used["z_coord"] ) )
F_NT_data <- t(t(as.matrix( m_data_used["F_NT"])))

p_vals_F <- matrix(0, nrow = 1, ncol = 16)

p_vals_F[1] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_1, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[2] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_2, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[3] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_3, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[4] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_4, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[5] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_5, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[6] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_6, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[7] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_7, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[8] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_8, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[9] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_9, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[10] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_10, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[11] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_11, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[12] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_12, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[13] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_13, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[14] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_14, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[15] <- PermVarTest(X, F_NT_data, sample_ico_unif_axis_15, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")
p_vals_F[16] <- PermVarTest(X, F_NT_data, sample_so3, a_B, g_dot, plot_flag = FALSE, m_lambda = a_lambda, d_X = "Great Circle")


I(p_vals_F > alpha)

