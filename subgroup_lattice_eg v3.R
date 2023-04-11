library(rdist)
library(MASS)
library(np)

gold_ratio <- ( 1 + sqrt(5) ) / 2
icosahedron <- matrix( 0 , nrow = 9, ncol = 3)
icosahedron[1, ] <- c(0, 1, gold_ratio)
icosahedron[2, ] <- c(0, 1, -gold_ratio)
icosahedron[3, ] <- c(1, gold_ratio, 0)
icosahedron[4, ] <- c(1, -gold_ratio, 0)
icosahedron[5, ] <- c(gold_ratio, 0, 1)
icosahedron[6, ] <- c(gold_ratio, 0, -1)
icosahedron <- icosahedron / sqrt( 1 + gold_ratio^2 )
icosahedron[7, ] <- c(1,0,0)
icosahedron[8, ] <- c(0,1,0)
icosahedron[9, ] <- c(0,0,1)

s1_offset <- 5

#### Helper Functions

find_second_smallest <- function (m_arr) {
  smallest_ind <- 0
  second_ind <- 0
  smallest <- Inf
  second <- Inf
  
  for (i in 1:length(m_arr)) {
    if (m_arr[i] < smallest) {
      second_ind <- smallest_ind
      second <- smallest
      smallest_ind <- i
      smallest <- m_arr[i]
    } else if (m_arr[i] < second) {
      second_ind <- i
      second <- m_arr[i]
    } 
  } 
  return(second_ind)
}

find_smallest_non_zero <- function (m_arr) {
  smallest_ind <- 0
  smallest <- Inf
  for (i in 1:length(m_arr)) {
    if (m_arr[i] < smallest) {
      if(m_arr[i] > 0) {
        smallest <- m_arr[i]
        smallest_ind <- i
      }
    }
  } 
  return(smallest_ind)
}

row_min_inds <- function( m_mat, print_flag = FALSE ) {
  output <- rep(0, dim(m_mat)[1])
  for (i in 1:dim(m_mat)[1]) {
    output[i] <- find_smallest_non_zero(m_mat[i,])
  }
  if (print_flag) {
    counts <- sort(table( output ), decreasing=T )
    print(counts)
  }
  return(output)
}

select_from_ind <- function( m_mat, m_inds ) {
  output <- rep(0, dim(m_mat)[1])
  for (i in 1:dim(m_mat)[1] ) {
    output[i] <- m_mat[i, m_inds[i]]
  }
  return(output)
}

rotation_matrix <- function( angle, axis ) {
    ## axis must be a three vector
    u <- axis / sqrt( sum( axis^2 ) )
    phi <- angle
    R <- matrix(0, nrow = 3, ncol = 3)
    R[1,1] <- cos(phi) + u[1]^2 * (1 - cos(phi))
    R[1,2] <- u[1] * u[2] * (1 - cos(phi)) - u[3] * sin(phi)
    R[1,3] <- u[1] * u[3] * (1 - cos(phi)) + u[2] * sin(phi)
    R[2,1] <- u[1] * u[2] * (1 - cos(phi)) + u[3] * sin(phi)
    R[2,2] <- cos(phi) + u[2]^2 * (1 - cos(phi))
    R[2,3] <- u[2] * u[3] * (1 - cos(phi)) - u[1] * sin(phi)
    R[3,1] <- u[3] * u[1] * (1 - cos(phi)) - u[2] * sin(phi)
    R[3,2] <- u[2] * u[3] * (1 - cos(phi)) + u[1] * sin(phi)
    R[3,3] <- cos(phi) + u[3]^2 * (1 - cos(phi))
    return(R)
}

G_hat <- function( rejection_vec ) {
  ## rejection_vec is an 11 vector from the simulations
  ## first 6 components relate to \langle R_{2 pi / 11}^{u_i} \rangle, then SO(3), then SL(3) 
  if (prod( rejection_vec[1:dim(icosahedron)[1]] ) == 0) { 
    if ( sum( rejection_vec[1:dim(icosahedron)[1]]) == 1 ) {
      return( which(rejection_vec[1:dim(icosahedron)[1]] == 1) + s1_offset )
    }
    return(1)
  }
  if ( rejection_vec[11] == 0 ) {
    return(2)
  } 
  return(3)
}

sample_G <- function( index ) {
  if (index <= dim(icosahedron)[1] ) {
    # Finite order 11 rotation around icosahedral axis
    power <- sample(1:10, size = 1)
    axis <- icosahedron[index, ]
    return( rotation_matrix( 2 * pi * power / 11, axis ) )
  }
  if (index == dim(icosahedron)[1] + 1) {
    # Samples from SO(3)
    angle <- runif( 1, min = 0, max = pi)
    axis <- rnorm(3)
    axis <- axis / sqrt( sum( axis^2 ) )
    return( rotation_matrix( angle , axis ) )
  }
  if (index == dim(icosahedron)[1] + 2 ) {
    # Samples from SL(3, R)
    A <- matrix( 0, nrow = 3, ncol = 3)
    while( det(A) == 0 ) {
      A <- matrix( rnorm(9), nrow = 3, ncol = 3)
    }
    a <- sign( det(A) ) * ( sign(det(A)) * det(A) )^(-1/3)
    return( a * A )
  }
}

find_MSE_np <- function (X, Y, X.test, Y.test, actual = FALSE) {
  m_df <- as.data.frame( cbind(X,Y) )
  colnames(m_df) <- c( "X1", "X2", "X3", "Y" )
  m_df_T <- as.data.frame( cbind(X.test,Y.test) )
  colnames(m_df_T) <- c( "X1", "X2", "X3", "Y" )
  
  bw <- npregbw(Y ~ X1 + X2 + X3, data = m_df, regtype = "lc")
  mod <- npreg(bw, data = m_df)
  fits <- predict( mod, data = m_df, newdata = m_df_T )
  
  if (actual) { return( mean( (fits - f(X.test))^2 ) ) }
  else { return( mean( (fits - Y.test)^2 ) ) }
}

find_MSE_S1 <- function(X, Y, X.test, Y.test, axis_ind, actual = FALSE) {
  axis <- icosahedron[axis_ind, ]
  m_df <- as.data.frame( cbind( rowSums( t( t(X) * axis ) ) , 
                                sqrt( rowSums(X^2) - rowSums(  t( t(X) * axis ) )^2 ), 
                                Y ) )
  colnames(m_df) <- c( "X1", "X2", "Y" )
  m_df_T <- as.data.frame( cbind( rowSums( t( t(X.test) * axis ) ) , 
                                  sqrt( rowSums(X.test^2) - rowSums( t( t(X.test) * axis ) )^2 ), 
                                  Y.test ) )
  colnames(m_df_T) <- c( "X1", "X2", "Y" )
  
  bw <- npregbw(Y ~ X1 + X2, data = m_df, regtype = "lc")
  mod <- npreg(bw, data = m_df)
  fits <- predict( mod, data = m_df, newdata = m_df_T )
  
  if (actual) { return( mean( (fits - f(X.test))^2 ) ) }
  else { return( mean( (fits - Y.test)^2 ) ) }
}

find_MSE_SO3 <- function (X, Y, X.test, Y.test, actual = FALSE) {
  m_df <- as.data.frame( cbind( sqrt( rowSums(X^2) ) , Y ) )
  colnames(m_df) <- c( "X1", "Y" )
  m_df_T <- as.data.frame( cbind( sqrt( rowSums(X.test^2) ) , Y.test) )
  colnames(m_df_T) <- c( "X1", "Y" )
  
  bw <- npregbw(Y ~ X1, data = m_df, regtype = "lc")
  mod <- npreg(bw, data = m_df)
  fits <- predict( mod, data = m_df, newdata = m_df_T )
  
  if (actual) { return( mean( (fits - f(X.test))^2 ) ) }
  else { return( mean( (fits - Y.test)^2 ) ) }
}

find_MSE_SL3 <- function (X, Y, X.test, Y.test, actual = FALSE) {
  fits <- mean( Y )
  if (actual) { return( mean( (fits - f(X.test))^2 ) ) }
  else { return( mean( (fits - Y.test)^2 ) ) }
}

#### Testing Functions

AsymVarTest <- function( m_X, m_Y, d_X, norm_Y, m_V, m_G_index, m_p_t, m_thresh, m_m, m_g_dot) {
  n <- dim(m_X)[1]
  dist_mat <- matrix(0, nrow = m_m , ncol = n)
  dif_mat <- matrix(0, nrow = m_m, ncol = n)
  
  for (i in 1:m_m) {
    g_i <- sample_G( m_G_index )
    ind <- sample(1:n, 1)
    X_i <- m_X[ind,]
    Y_i <- m_Y[ind,]
    g_cdot_X_i <- array( m_g_dot(g_i, X_i) )
    dist_mat[i,] <- d_X( g_cdot_X_i, m_X )
    dif_mat[i,] <- norm_Y( Y_i, m_Y )
  }
  inds <- row_min_inds(dist_mat)
  dists <- select_from_ind(dist_mat, inds)
  difs <- select_from_ind(dif_mat, inds)
  
  p_val <- pbinom(sum(I( m_V(dists) < difs - m_thresh) ) - 1, 
                  sum(I(difs>0)), m_p_t(m_thresh), lower.tail = FALSE )
   
  return(p_val)
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
  return( sin( - norm ) )
}

f2 <- function( X ) {
  norm <- sqrt( rowSums(X^2) )
  return( sin( - norm ) )
}

f3 <- function( X ) {
  return( sin( - abs(X[,3]) ) )
}

f4 <- function( X ) {
  return( sin( - abs(X[,1]) ) + cos( sqrt( X[,2]^2 + X[,3]^2 ))) 
}

f5 <- function( X ) {
  return (sin( - abs(X[,1] ) ) + 0.2 * X[,2] - cos(  X[,3]  ) )
}

a_p_t <- function ( m_thres ) {
  return( 2 * exp( - (m_thres / sigma) ^2 / 4) / ( (m_thres / sigma) * sqrt(2 * pi )) )
}

## Sims

run_sims <- function( f, sigma,  Sigma_X, Sigma_X_test, ns, num_sims, actual = FALSE ) {
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
      
      p_vals[i, 1] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 1, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 2] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 2, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 3] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 3, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 4] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 4, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 5] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 5, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 6] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 6, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 7] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 7, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 8] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 8, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 9] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 9, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 10] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 10, a_p_t, 2*sigma, m, g_dot)
      p_vals[i, 11] <- AsymVarTest(X, Y, d_X, norm_Y, a_V, 11, a_p_t, 2*sigma, m, g_dot)
      
      rejections[i, ] <- I(p_vals[i, ] > alpha)
      G_hats[i, k] <- G_hat(rejections[i, ])

      p_vals_split[i, 1] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 1, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 2] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 2, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 3] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 3, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 4] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 4, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 5] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 5, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 6] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 6, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 7] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 7, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 8] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 8, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 9] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 9, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 10] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 10, a_p_t, 2*sigma, m, g_dot)
      p_vals_split[i, 11] <- AsymVarTest(X[1:n/2, ], t(t(Y[1:n/2, ])), d_X, norm_Y, a_V, 11, a_p_t, 2*sigma, m, g_dot)
      
      rejections_split[i, ] <- I(p_vals_split[i, ] > alpha)
      G_hats_split[i, k] <- G_hat(rejections_split[i, ])    
      
      MSEs_Baseline[i, k] <- find_MSE_np(X, Y, X.test, Y.test)
      
      if ( G_hats[i, k] == 1 ) {
        MSEs_Full[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats[i, k] == 2 ) {
        MSEs_Full[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else if ( G_hats[i, k] == 3 ) {
        MSEs_Full[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test)
      } else {
        MSEs_Full[i, k] <- find_MSE_S1(X, Y, X.test, Y.test, G_hats[i, k] - s1_offset )
      }
      
      if ( G_hats_split[i, k] == 1 ) {
        MSEs_Split[i, k] <- find_MSE_np(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else if ( G_hats_split[i, k] == 2 ) {
        MSEs_Split[i, k] <- find_MSE_SO3(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else if ( G_hats_split[i, k] == 3 ) {
        MSEs_Split[i, k] <- find_MSE_SL3(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else {
        MSEs_Split[i, k] <- find_MSE_S1(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test, G_hats_split[i, k] - s1_offset )
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

num_sims = 100
alpha <- 0.05
ns <- c(10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000) ## Sample Sizes for Example 5.2
ns <- c( 50, 100, 250, 500, 750, 1000) ## Sample Sizes for Appendix 
mu_X <- rep(0, 3)
sigma <- 0.01

set.seed(1011)
scenario_1 <- run_sims(f2, sigma, diag(c(0.1,0.1,2)), diag(c(2,2,2)), ns, num_sims)
write.csv(scenario_2, file = "scenario_1_LCE.csv")

set.seed(38)
scenario_2 <- run_sims(f4, sigma, diag(c(0.1,0.1,2)), diag(c(2,2,2)), ns, num_sims)
write.csv(scenario_4, file = "scenario_2_LCE.csv")

set.seed(24)
scenario_3 <- run_sims(f5, sigma, diag(c(0.1,0.1,2)), diag(c(2,2,2)), ns, num_sims)
write.csv(scenario_6, file = "scenario_3_LCE.csv")

set.seed(1010)
scenario_4 <- run_sims(f1, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims)
write.csv(scenario_1, file = "scenario_4_LCE.csv")

set.seed(1012)
scenario_5 <- run_sims(f4, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims)
write.csv(scenario_3, file = "scenario_5_LCE.csv")

set.seed(15)
scenario_6 <- run_sims(f5, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims)
write.csv(scenario_5, file = "scenario_6_LCE.csv")


#### Plots

## Plots of MSPEs

scenario <- scenario_1

m <- length(ns)

plot( ns, colMeans(scenario[,(2*m+1):(3*m)]), type = "l", ylim = c(-sigma, 0.9), 
      ylab = "MSPE", xlab = "Sample Size",  lty = 2)
lines( ns, colMeans(scenario[,(3*m+1):(4*m)]), col = "blue", lty = 2)
lines( ns, colMeans(scenario[,(4*m+1):(5*m)]), col = "red", lty = 2)

points( ns, colMeans(scenario[,(2*m+1):(3*m)]), pch = 15)
points( ns, colMeans(scenario[,(3*m+1):(4*m)]), col = "blue", pch = 16)
points( ns, colMeans(scenario[,(4*m+1):(5*m)]), col = "red", pch = 17)


## Plots of \hat(G) = G_max

scenario <- scenario_5
g_hat_correct_ind <- 1
m <- length(ns)

plot( ns, colMeans( I( scenario[,1:m] == g_hat_correct_ind ) ), type = "l", ylim = c(-0.1, 1.1), 
      ylab ="", xlab = "Sample Size",  lty = 2, col = "blue")
lines( ns, colMeans( I( scenario[,(m+1):(2*m)] == g_hat_correct_ind ) ), col = "red", lty = 2)
lines( c(0,1000), c(1,1), lty = 3)

points( ns, colMeans( I( scenario[,1:m] == g_hat_correct_ind ) ), col = "blue", pch = 16)
points( ns, colMeans( I( scenario[,(m+1):(2*m)] == g_hat_correct_ind ) ), col = "red", pch = 17)

#### Extra for interest

## Finding mean scaling from sample_G(8)
m <- 100000
ratios <- rep(0, m)
for (i in 1:m) {
  a_X <- as.numeric(mvrnorm(1, mu_X, Sigma_X))
  g_X <- as.numeric(sample_G(9) %*% a_X)
  ratios[i] <- sqrt( sum( g_X^2 ) ) / sqrt( sum( a_X^2 ) )
}
hist(ratios, breaks = 30, xlim = c(-0.1,10) )
lines(c(1,1), c(0,m), lty = 1, col = "blue")
lines(c(mean(ratios),mean(ratios)), c(0,m), lty = 1, col = "red")
mean(ratios)
