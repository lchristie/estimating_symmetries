library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("test_function.R")

dimension <- 6
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

#### Group Elements
identity_mat <- diag( rep(1, dimension) )
R_h <-  diag( c( -1, rep(1, dimension - 1) ) )
R_v <-  diag( c( 1, -1, rep(1, dimension - 2) ) )
R_diag_forw <-  diag( c( 0, 0, rep(1, dimension - 2) ) )
R_diag_forw[1,2] <- 1
R_diag_forw[2,1] <- 1
R_diag_back <-  diag( c( 0, 0, rep(1, dimension - 2) ) )
R_diag_back[1,2] <- -1
R_diag_back[2,1] <- -1
R_90 <-  diag( c( 0, 0, rep(1, dimension - 2) ) )
R_90[1,2] <- -1
R_90[2,1] <- 1
R_180 <-  diag( c( -1, -1,  rep(1, dimension - 2) ) )
R_270 <-  diag( c( 0, 0, rep(1, dimension - 2) ) )
R_270[1,2] <- 1
R_270[2,1] <- -1


#### Helper Functions

sample_Rh <- function() {
  power = sample(2,1)
  if (power == 2) {
    return( identity_mat )
  } else {
    return( R_h )
  }
}
sample_Rv <- function() {
  power = sample(2,1)
  if (power == 2) {
    return( identity_mat )
  } else {
    return( R_v )
  }
}
sample_Rforwardslash <- function() {
  power = sample(2,1)
  if (power == 2) {
    return( identity_mat )
  } else {
    return( R_diag_forw )
  }
}
sample_Rbackslash <- function() {
  power = sample(2,1)
  if (power == 2) {
    return( identity_mat )
  } else {
    return( R_diag_back )
  }
}
sample_R180 <- function() {
  power = sample(2,1)
  if (power == 2) {
    return( identity_mat )
  } else {
    return( R_180 )
  }
}
sample_R90 <- function() {
  power = sample(4,1)
  if (power == 1) {
    return( R_90 )
  } 
  if (power == 2) {
    return( R_180 )
  }
  if (power == 3) {
    return( R_270 )
  }
  if (power == 4) {
    return( identity_mat )
  }
}


K_hat_sym_diff <- function ( TestRejections ) {
  ## This function tracks the error of the confidence region
  ## The true Test rejection vector is (1,1,1,-1,-1,-1)
  ## This can be sped up as a direct computation rather than conditionals, but this is fast enough for now

  mathcal_E <- 0
  if ( sum( TestRejections[1:3] ) < 3 ) {
    ## This case covers the error of $\langle R_h, R_v \rangle, $R_h$, $R_v$ and $R_\pi$
    mathcal_E <- mathcal_E + 1 + (3 - sum( TestRejections[1:3] ))
  }
  if ( TestRejections[6] * TestRejections[3] == 1) {
    ## This is the error of $\langle R_{\pi / 2} \rangle$
    mathcal_E <- mathcal_E + 1
  }
  if ( sum( TestRejections[3:5] ) == 3) {
    ## This is the error of $\langle R_/, R_\ \rangle$, $R_/$, $R_\$, and $R_\pi$
    mathcal_E <- mathcal_E + 1 + sum( TestRejections[4:5] )
  }
  if ( sum( TestRejections ) == 6 ) {
    ## This case covers the error of $D_4$
    mathcal_E <- mathcal_E + 1
  }
  return( mathcal_E )
}



G_hat_from_TR <- function ( TestRejections ) {
  ## This returns a vector indicating inclusion from the lattice base, G_hat is the join of these elements
  if ( sum(TestRejections) <= 1 ) {
    ## This case is the easiest, either returning $I$ or the singular accepted group (except R_270)
    output <- TestRejections
    output[6] <- 0
    return( output )
  }
  if ( TestRejections[3] == 0 ) {
    ## This case has only order two maximal subgroups
    possibilities = which( TestRejections[1:5] == 1, arr.ind = TRUE )
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    output[ind] <- 1
    return( output )
  } 
  if ( sum(TestRejections) == 2) {
    ## Recall that one of the accepted symmetries is G_3
    if ( TestRejections[6] == 1 ) {
      return( TestRejections )
    } else {
      possibilities = which( TestRejections == 1, arr.ind = TRUE )
      output <- rep(0, 6)
      ind <- true_sample( possibilities )
      output[ind] <- 1
      return( output )
    }
  }
  if ( (sum( TestRejections * c(1,1,1,0,0,0) ) == 3) * (sum(TestRejections) == 3) ) {
    return( TestRejections )
  }
  if ( (sum( TestRejections * c(0,0,1,1,1,0) ) == 3) * (sum(TestRejections) == 3) ) {
    return( TestRejections )
  }
  if ( (sum( TestRejections * c(1,1,1,1,1,0) ) == 5) * (sum(TestRejections) == 5) ) {
    if ( sample(2,1) == 1 ) {
      return( c(1,1,1,0,0,0) )
    } else {
      return( c(0,0,1,1,1,0) )
    }
  }
  if ( sum(TestRejections) == 6 ) {
    ## This is the case where $G_hat = D_4$
    return(TestRejections)
  }
  if ( sum(TestRejections) == 3 ) {
    ## Recall TestRejections[3] == 1, but cannot form a large maximal group
    possibilities = which( TestRejections == 1, arr.ind = TRUE )
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    output[ind] <- 1
    return( output )
  }
  if ( sum(TestRejections[1:3]) == 3) {
    possibilities = which( TestRejections[3:6] == 1, arr.ind = TRUE ) + 2
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    output[ind] <- 1
    if (output[3] == 1) {
      output[1] <- 1
      output[2] <- 1
    }
    return( output )
  }
  if ( sum(TestRejections[3:5]) == 3) {
    possibilities = which( TestRejections[c(1,2,3,6)] == 1, arr.ind = TRUE )
    output <- rep(0, 6)
    ind <- true_sample( possibilities )
    if (ind == 4) {ind <- 6}
    output[ind] <- 1
    if (output[3] == 1) {
      output[4] <- 1
      output[5] <- 1
    }
    return( output )
  }
  return( rep(-1, 6) ) ## This is an error code
}



#### Simulations

## Set-up Pars

d_X <- function ( x_i, X ) {
  # print( t(x_i) )
  # print( X )
  # return( sqrt( rowSums( (t(x_i) - X )^2 ) ) )
  # print( sweep( X, 2, t(x_i)) )
  return(sqrt( rowSums( sweep( X, 2, t(x_i))^2 )))
}

norm_Y <- function ( y_i, Y ) {
  return( sqrt( colSums( (y_i - t(Y) )^2 ) ) )
}

a_V_e <- function( d ) {
  return( exp(-1) * d )
}

f_d <- function( X ) {
  norm <- sqrt( X[,1]^2 )
  return( exp( -norm ) )
}

a_p_t <- function ( m_thres ) {
  return( 2 * exp( - (m_thres / 0.05) ^2 / 4) / ( (m_thres / 0.05) * sqrt(2 * pi )) )
}

g_dot <- function( g, X ) {
  return( g %*% X )
}

#### f_d Sims

run_sims <- function( num_sims ) {
  
  rejections_fd_avt <- matrix(0, ncol = num_sims * length(ns), nrow = 6)
  rejections_fd_pv <- matrix(0, ncol = num_sims * length(ns), nrow = 6)
  rejections_fd_pv_lambda <- matrix(0, ncol = num_sims * length(ns), nrow = 6)
  
  pb <- txtProgressBar(min = 0, max = sum(ns) * num_sims, style = 3, width = 50, char = "=") 
  counter <- 0
  
  for (k in 1:length(ns)) {
    n <- ns[k]
    m <- ns[k]
    
    for (i in 1:num_sims) {
      X <- mvrnorm(n, mu_X, Sigma_X)
      Y <- as.matrix( f_d(X) + rnorm(n, 0, sigma) )
      
      rejections_fd_avt[1, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, sample_Rh, a_p_t, 2*sigma, m, g_dot)
      rejections_fd_avt[2, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, sample_Rv, a_p_t, 2*sigma, m, g_dot)
      rejections_fd_avt[3, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, sample_R180, a_p_t, 2*sigma, m, g_dot)
      rejections_fd_avt[4, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, sample_Rforwardslash, a_p_t, 2*sigma, m, g_dot)
      rejections_fd_avt[5, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, sample_Rbackslash, a_p_t, 2*sigma, m, g_dot)
      rejections_fd_avt[6, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, sample_R90, a_p_t, 2*sigma, m, g_dot)
      
      rejections_fd_pv[1, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rh, a_B, g_dot)
      rejections_fd_pv[2, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rv, a_B, g_dot)
      rejections_fd_pv[3, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_R180, a_B, g_dot)
      rejections_fd_pv[4, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rforwardslash, a_B, g_dot)
      rejections_fd_pv[5, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rbackslash, a_B, g_dot)
      rejections_fd_pv[6, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_R90, a_B, g_dot)
      
      rejections_fd_pv_lambda[1, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rh, a_B, g_dot, m_lambda = -1)
      rejections_fd_pv_lambda[2, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rv, a_B, g_dot, m_lambda = -1)
      rejections_fd_pv_lambda[3, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_R180, a_B, g_dot, m_lambda = -1)
      rejections_fd_pv_lambda[4, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rforwardslash, a_B, g_dot, m_lambda = -1)
      rejections_fd_pv_lambda[5, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_Rbackslash, a_B, g_dot, m_lambda = -1)
      rejections_fd_pv_lambda[6, (k-1)*num_sims + i] <- PermVarTest(X, Y, sample_R90, a_B, g_dot, m_lambda = -1)
      
      counter <- counter + n
      setTxtProgressBar(pb, counter)
    }
    
    print(paste("--Acceptances for n=", n, "--"))
    print("Asymmetric Variation Test")
    print( rowMeans( rejections_fd_avt[,((k-1)*num_sims + 1):(k*num_sims)] > alpha ) )
    print("Permutation Variant")
    print( rowMeans( rejections_fd_pv[,((k-1)*num_sims + 1):(k*num_sims)] > alpha ) )
    print("Permutation Variant, non zero Lambda")
    print( rowMeans( rejections_fd_pv_lambda[,((k-1)*num_sims + 1):(k*num_sims)] > alpha ) )
  }
  write.csv( rejections_fd_avt, file = avt.file.name)
  write.csv( rejections_fd_pv, file = pv.file.name )
  write.csv( rejections_fd_pv_lambda, file = pv_lambda.file.name )
  
}

set.seed(1010)
run_sims(num_sims)



#### Timing Tests
set.seed(1010)
num_sims = 100
alpha <- 0.05

mu_X <- rep(0, dimension)
Sigma_X = diag(c(2,1,0.5,0.25))

sigma <- 0.05

a_B <- 100

n <- 100
m <- 100
X <- mvrnorm(n, mu_X, Sigma_X)
Y <- as.matrix( f_d(X) + rnorm(n, 0, sigma) )

system.time( rejections_fd_avt[1, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, 1, a_p_t, 2*sigma, m) )
## Time = 0.004 seconds

system.time( PermVarTest(X, Y, 1, a_B) )
## Time = 0.118 Seconds




