library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("ico_syms.R")
source("test_function.R")

G_hat <- function( rejection_vec ) {
  ## first 9 components relate to \langle R_{2 pi / 11}^{u_i} \rangle, then SO(3), then SL(3) 
  if (sum( rejection_vec[1:9] ) == 0) {
    return(0)
  }
  if (prod( rejection_vec[1:9] ) * rejection_vec[10]  == 0) {
    accepted <- which( rejection_vec[1:9] == 1, arr.ind = TRUE)
    ico_ind <- accepted[sample( length(accepted), 1)]
    if (ico_ind <= 3) {return( ico_ind + 12) }
    else {return( ico_ind - 3 ) }
  } 
  return(16)
}


#### Simulations

## Set-up Pars

g_dot <- function( g, X ) {
  return( g %*% X )
}

d_X <- function ( x_i, X ) {
  return( sqrt( colSums( ( t(X) - as.numeric(x_i) )^2 ) ) )
}

norm_Y <- function ( y_i, Y ) {
  return( sqrt( colSums( (y_i - t(Y) )^2 ) ) )
}

a_V <- function( d ) {
  return( d )
}

f1 <- function( X ) {
  norm <- sqrt( rowSums(X^2) )
  return( sin( norm ) )
}

f2 <- function ( X ) {
  return ( sin( abs( X[,1]) ) + cos( sqrt( X[,2]^2 + X[,3]^2) ) )
}

f3 <- function ( X ) {
  return ( sin( abs( X[,1]) ) + 0.2 * X[,2] - cos( X[,3] ) )
  
}

f4 <- function ( X ) {
  rot_mat <- rotation_matrix(0.1, c(0,0,1) )
  n <- dim(X)[1]
  X_out <- matrix( 0, nrow = n, ncol = dim(X)[2] ) 
  for (i in 1:n) {
    X_out[i,] <- t( rot_mat %*% X[i,] )
  }
  return ( sin( abs( X_out[,1]) ) + cos( sqrt( X_out[,2]^2 + X_out[,3]^2) ) )
}

a_p_t <- function ( m_thres ) {
  return( 2 * exp( - (m_thres / sigma) ^2 / 4) / ( (m_thres / sigma) * sqrt(2 * pi )) )
}


#### Sims

run_sims_perm_var <- function( f, sigma,  Sigma_X, Sigma_X_test, ns, num_sims, num_perms ) {
  G_hats <- matrix(0, ncol = length(ns), nrow = num_sims)
  G_hats_split <- matrix(0, ncol = length(ns), nrow = num_sims)
  
  MSEs_Baseline <-  matrix(0, ncol = length(ns), nrow = num_sims)
  MSEs_Full <-  matrix(0, ncol = length(ns), nrow = num_sims)
  MSEs_Split <-  matrix(0, ncol = length(ns), nrow = num_sims)
  
  pb <- txtProgressBar(min = 0, max = sum(ns) * num_sims, style = 3, width = 50, char = "=") 
  counter <- 0
  
  for (k in 1:length(ns)) {
    n <- ns[k]
    m <- ns[k]
    p_vals <- matrix(0, nrow = num_sims, ncol = 11)
    p_vals_split <- matrix(0, nrow = num_sims, ncol = 11)
    rejections <- matrix(0, nrow = num_sims, ncol = 11)
    rejections_split <- matrix(0, nrow = num_sims, ncol = 11)
    
    for (i in 1:num_sims) {
      X <- mvrnorm(n, mu_X, Sigma_X)
      Y <- as.matrix( f(X) + rnorm(n, 0, sigma) )
      X.test <- mvrnorm(n, mu_X, Sigma_X_test) 
      Y.test <- as.matrix( f(X.test) + rnorm(n, 0, sigma) )

      p_vals[i, 1] <- PermVarTest(X, Y, sample_ico_unif_axis_13, num_perms, g_dot)
      p_vals[i, 2] <- PermVarTest(X, Y, sample_ico_unif_axis_14, num_perms, g_dot)
      p_vals[i, 3] <- PermVarTest(X, Y, sample_ico_unif_axis_15, num_perms, g_dot)
      p_vals[i, 4] <- PermVarTest(X, Y, sample_ico_unif_axis_1, num_perms, g_dot)
      p_vals[i, 5] <- PermVarTest(X, Y, sample_ico_unif_axis_2, num_perms, g_dot)
      p_vals[i, 6] <- PermVarTest(X, Y, sample_ico_unif_axis_3, num_perms, g_dot)
      p_vals[i, 7] <- PermVarTest(X, Y, sample_ico_unif_axis_4, num_perms, g_dot)
      p_vals[i, 8] <- PermVarTest(X, Y, sample_ico_unif_axis_5, num_perms, g_dot)
      p_vals[i, 9] <- PermVarTest(X, Y, sample_ico_unif_axis_6, num_perms, g_dot)
      p_vals[i, 10] <- PermVarTest(X, Y, sample_so3, num_perms, g_dot)
      
      rejections[i, ] <- I(p_vals[i, ] > alpha)
      G_hats[i, k] <- G_hat(rejections[i, ])
      
      p_vals_split[i, 1] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_13, num_perms, g_dot)
      p_vals_split[i, 2] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_14, num_perms, g_dot)
      p_vals_split[i, 3] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_15, num_perms, g_dot)
      p_vals_split[i, 4] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_1, num_perms, g_dot)
      p_vals_split[i, 5] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_2, num_perms, g_dot)
      p_vals_split[i, 6] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_3, num_perms, g_dot)
      p_vals_split[i, 7] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_4, num_perms, g_dot)
      p_vals_split[i, 8] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_5, num_perms, g_dot)
      p_vals_split[i, 9] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_ico_unif_axis_6, num_perms, g_dot)
      p_vals_split[i, 10] <- PermVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), sample_so3, num_perms, g_dot)

      rejections_split[i, ] <- I(p_vals_split[i, ] > alpha)
      G_hats_split[i, k] <- G_hat(rejections_split[i, ])    
      
      MSEs_Baseline[i, k] <- find_MSE_np(X, Y, X.test, Y.test)
      
      if ( G_hats[i, k] == 0 ) {
        MSEs_Full[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats[i,k] < 16 ) {
        MSEs_Full[i,k] <- find_MSE_S1u(X, Y, X.test, Y.test, G_hats[i,k] )
      } else if ( G_hats[i, k] == 16 ) {
        MSEs_Full[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else {
        MSEs_Full[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test) 
      }
      
      if ( G_hats_split[i, k] == 0 ) {
        MSEs_Split[i, k] <- find_MSE_np(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else if ( G_hats_split[i, k] < 16 ) {
        MSEs_Split[i, k] <- find_MSE_S1u(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test, G_hats_split[i, k] )
      } else if ( G_hats_split[i, k] == 16 ) {
        MSEs_Split[i, k] <- find_MSE_SO3(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } 
      
      counter <- counter + n
      setTxtProgressBar(pb, counter)
    }
    
    print("----------")
    print( n )
    print( mean( MSEs_Baseline[, k] ) )
    print( mean( MSEs_Full[, k] ) )
    print( mean( MSEs_Split[, k] ) )
    print( table( G_hats[ ,k ]))
    print( table( G_hats_split[ , k]))
    
  }
  
  return( cbind( G_hats, G_hats_split, MSEs_Baseline, MSEs_Full, MSEs_Split ) )
}

num_sims <- 100
num_perms <- 100
alpha <- 0.05
ns <- c(10, 20, 30, 40, 50, 75, 100, 150, 200)
mu_X <- rep(0, 3)
sigma <- 0.01



set.seed(1010)

scenario_1_perm <- run_sims_perm_var(f1, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_1_perm, file = "sims/scenario_1_LCE_perms.csv")
set.seed(1010)
scenario_2_perm <- run_sims_perm_var(f2, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_2_perm, file = "sims/scenario_2_LCE_perms.csv")
set.seed(1010)
scenario_3_perm <- run_sims_perm_var(f3, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_3_perm, file = "sims/scenario_3_LCE_perms.csv")
set.seed(1010)
scenario_4_perm <- run_sims_perm_var(f4, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_4_perm, file = "sims/scenario_4_LCE_perms.csv")


set.seed(1010)

scenario_2_perm_varying_K <- run_sims_perm_var_varying_K(f2, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_2_perm_varying_K, file = "sims/scenario_2_LCE_perms_varying_K_2.csv")
set.seed(1010)
scenario_4_perm_varying_K <- run_sims_perm_var_varying_K(f4, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_4_perm_varying_K, file = "sims/scenario_4_LCE_perms_varying_K_2.csv")



