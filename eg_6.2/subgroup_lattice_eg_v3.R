library(rdist)
library(MASS)
library(np)

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
  ## first 9 components relate to \langle R_{2 pi / 11}^{u_i} \rangle, then SO(3), then SL(3) 
  if (sum( rejection_vec[1:9] ) == 0) {
    return(0)
  }
  if (prod( rejection_vec[1:9] ) * rejection_vec[10]  == 0) {
    accepted <- which( rejection_vec[1:6] == 1, arr.ind = TRUE)
    ico_ind <- accepted[sample( length(accepted), 1)]
    return(ico_ind)
  }
  if ( rejection_vec[11] == 0 ) {
    return(10)
  } 
  return(11)
}

sample_G <- function( index ) {
  if (index < 16) {
    # Finite order 11 rotation around icosahedral axis
    angle <- runif(1, min = 0, max = 2 * pi)
    axis <- icosahedron[index, ]
    return( rotation_matrix( angle, axis ) )
  }
  if (index == 16) {
    # Samples from SO(3)
    angle <- runif( 1, min = 0, max = pi)
    axis <- rnorm(3)
    axis <- axis / sqrt( sum( axis^2 ) )
    return( rotation_matrix( angle , axis ) )
  }
  if (index == 17) {
    # Samples from SL(3, R)
    A <- matrix( 0, nrow = 3, ncol = 3)
    while( det(A) == 0 ) {
      A <- matrix( rnorm(9), nrow = 3, ncol = 3)
    }
    a <- sign( det(A) ) * ( sign(det(A)) * det(A) )^(-1/3)
    return( a * A )
  }
}

random_transforms_X <- function ( m_X, G_ind, G_action ) {
  n <- dim(m_X)[1]
  X_out <- matrix( 0, nrow = n, ncol = dim(m_X)[2] ) 
  for (i in 1:n) {
    g_i <- sample_G( G_ind )
    X_out[i,] <- array( G_action(g_i, m_X[i,]) )
  }
  return( X_out )
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

find_MSE_S1u <- function (X, Y, X.test, Y.test, ico_ind, actual = FALSE) {
  proj_vec <- icosahedron[ico_ind, ]
  X.prime.1 <- X %*% proj_vec
  X.prime.2 <- sqrt( rowSums(X^2) - X.prime.1^2 )
  
  m_df <- as.data.frame( cbind( X.prime.1, X.prime.2 , Y ) )
  colnames(m_df) <- c( "X1", "X2", "Y" )
  
  X.prime.test.1 <- X.test %*% proj_vec
  X.prime.test.2 <- sqrt( rowSums(X.test^2) - X.prime.test.1^2 )
  
  m_df_T <- as.data.frame( cbind( X.prime.test.1, X.prime.test.2, Y.test ) )
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
  
  # print("-------------")
  # print( sum(I( m_V(dists) < difs - m_thresh) ) )
  # print( sum(I(difs>0)) )
  # print( m_p_t(m_thresh) )
  
  p_val <- pbinom(sum(I( m_V(dists) < difs - m_thresh) ) - 1, 
                  sum(I(difs>0)), m_p_t(m_thresh), lower.tail = FALSE )
   
  return(p_val)
}

PermVarTest_2 <- function ( m_X, m_Y, m_q, m_G_index, m_B, m_g_dot, plot_flag = FALSE) {
  n <- dim(m_X)[1]
  
  ## Baseline Stat From Original Data
  dist_mat_orig <- pdist( m_X )
  dif_mat_orig <- pdist( m_Y )
  orig_ratios <- dif_mat_orig[which(dist_mat_orig != 0)] / dist_mat_orig[which(dist_mat_orig != 0)]
  A_0 <- quantile(orig_ratios, probs = m_q, type = 1, names = FALSE)
  if (plot_flag) {
    plot( dist_mat_orig, dif_mat_orig, type = "p", pch = 4)
    lines( c(0, 100), c(0, 100 * A_0) )
  }
  ## Permutations, Interchangable under H_0
  As <- rep(0, m_B)
  
  for (i in 1:m_B) {
    new_X <- random_transforms_X( m_X, m_G_index, m_g_dot )
    dist_mat <- pdist( new_X )
    ratios <- dif_mat_orig[which(dist_mat != 0)] / dist_mat[which(dist_mat != 0)]
    As[i] <- quantile(ratios, probs = m_q, type = 1, names = FALSE)
    if (plot_flag) { lines( c(0, 100), c(0, 100 * As[i]), col = "red", lty = 2) }
  }
  if (plot_flag) { lines( c(0, 100), c(0, 100 * A_0) ) }
  
  ## Under the alternative we would expect larger ratios for the permutations
  p_val <- mean( I(As <= A_0) )
  return( p_val )
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



#### Sims

run_sims <- function( f, sigma,  Sigma_X, Sigma_X_test, ns, num_sims ) {
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
      # print(rejections[i,])
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
      
      if ( G_hats[i, k] == 0 ) {
        MSEs_Full[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats[i,k] < 10 ) {
        MSEs_Full[i,k] <- find_MSE_S1u(X, Y, X.test, Y.test, G_hats[i,k] )
      } else if ( G_hats[i, k] == 10 ) {
        MSEs_Full[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else {
        MSEs_Full[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test)
      }
      
      if ( G_hats_split[i, k] == 0 ) {
        MSEs_Split[i, k] <- find_MSE_np(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else if ( G_hats_split[i, k] < 10 ) {
        MSEs_Split[i, k] <- find_MSE_S1u(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test, G_hats_split[i, k] )
      } else if ( G_hats_split[i, k] == 10 ) {
        MSEs_Split[i, k] <- find_MSE_SO3(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else {
        MSEs_Split[i, k] <- find_MSE_SL3(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
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
      
      p_vals[i, 1] <- PermVarTest_2(X, Y, 1, 1, num_perms, g_dot)
      p_vals[i, 2] <- PermVarTest_2(X, Y, 1, 2, num_perms, g_dot)
      p_vals[i, 3] <- PermVarTest_2(X, Y, 1, 3, num_perms, g_dot)
      p_vals[i, 4] <- PermVarTest_2(X, Y, 1, 4, num_perms, g_dot)
      p_vals[i, 5] <- PermVarTest_2(X, Y, 1, 5, num_perms, g_dot)
      p_vals[i, 6] <- PermVarTest_2(X, Y, 1, 6, num_perms, g_dot)
      p_vals[i, 7] <- PermVarTest_2(X, Y, 1, 7, num_perms, g_dot)
      p_vals[i, 8] <- PermVarTest_2(X, Y, 1, 8, num_perms, g_dot)
      p_vals[i, 9] <- PermVarTest_2(X, Y, 1, 9, num_perms, g_dot)
      p_vals[i, 10] <- PermVarTest_2(X, Y, 1, 10, num_perms, g_dot)
      # p_vals[i, 11] <- PermVarTest_2(X, Y, 1, 8, num_perms, g_dot)
      
      
      rejections[i, ] <- I(p_vals[i, ] > alpha)
      G_hats[i, k] <- G_hat(rejections[i, ])
      
      p_vals_split[i, 1] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 1, num_perms, g_dot)
      p_vals_split[i, 2] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 2, num_perms, g_dot)
      p_vals_split[i, 3] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 3, num_perms, g_dot)
      p_vals_split[i, 4] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 4, num_perms, g_dot)
      p_vals_split[i, 5] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 5, num_perms, g_dot)
      p_vals_split[i, 6] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 6, num_perms, g_dot)
      p_vals_split[i, 7] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 7, num_perms, g_dot)
      p_vals_split[i, 8] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 8, num_perms, g_dot)
      p_vals_split[i, 9] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 9, num_perms, g_dot)
      p_vals_split[i, 10] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 10, num_perms, g_dot)
      # p_vals_split[i, 11] <- PermVarTest_2(X[1:n/2, ], t(t(Y[1:n/2, ])), 1, 8, num_perms, g_dot)
      
      rejections_split[i, ] <- I(p_vals_split[i, ] > alpha)
      G_hats_split[i, k] <- G_hat(rejections_split[i, ])    
      
      MSEs_Baseline[i, k] <- find_MSE_np(X, Y, X.test, Y.test)
      
      if ( G_hats[i, k] == 0 ) {
        MSEs_Full[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats[i,k] < 10 ) {
        MSEs_Full[i,k] <- find_MSE_S1u(X, Y, X.test, Y.test, G_hats[i,k] )
      } else if ( G_hats[i, k] == 10 ) {
        MSEs_Full[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else {
        MSEs_Full[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test) 
      }
      
      if ( G_hats_split[i, k] == 0 ) {
        MSEs_Split[i, k] <- find_MSE_np(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else if ( G_hats_split[i, k] < 10 ) {
        MSEs_Split[i, k] <- find_MSE_S1u(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test, G_hats_split[i, k] )
      } else if ( G_hats_split[i, k] == 10 ) {
        MSEs_Split[i, k] <- find_MSE_SO3(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
      } else {
        MSEs_Split[i, k] <- find_MSE_SL3(X[(n/2 + 1):n, ], Y[(n/2 + 1):n, ], X.test, Y.test)
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
write.csv(scenario_1_perm, file = "scenario_1_LCE_perms.csv")
set.seed(1010)
scenario_2_perm <- run_sims_perm_var(f2, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_2_perm, file = "scenario_2_LCE_perms.csv")
set.seed(1010)
scenario_3_perm <- run_sims_perm_var(f3, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_3_perm, file = "scenario_3_LCE_perms.csv")
set.seed(1010)
scenario_4_perm <- run_sims_perm_var(f4, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_4_perm, file = "scenario_4_LCE_perms.csv")

scenario_1_perm <- read.csv( "scenario_1_LCE_perms.csv", header = TRUE )[,-1]
scenario_2_perm <- read.csv( "scenario_2_LCE_perms.csv", header = TRUE )[,-1]
scenario_3_perm <- read.csv( "scenario_3_LCE_perms.csv", header = TRUE )[,-1]
scenario_4_perm <- read.csv( "scenario_4_LCE_perms.csv", header = TRUE )[,-1]

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
colMeans( scenario_1_perm[, inds] == 10)
colMeans( scenario_2_perm[, inds] == 7)
colMeans( scenario_3_perm[, inds] == 0)
colMeans( scenario_4_perm[, inds] == 0)




#### Effect of Changing K

icosahedron <- matrix( 0 , nrow = 15, ncol = 3)
icosahedron[1, ] <- c(0, 1, gold_ratio)
icosahedron[2, ] <- c(0, 1, -gold_ratio)
icosahedron[3, ] <- c(1, gold_ratio, 0)
icosahedron[4, ] <- c(1, -gold_ratio, 0)
icosahedron[5, ] <- c(gold_ratio, 0, 1)
icosahedron[6, ] <- c(gold_ratio, 0, -1)
icosahedron[7, ] <- c(1, 0, gold_ratio)
icosahedron[8, ] <- c(1, 0, -gold_ratio)
icosahedron[9, ] <- c(gold_ratio, 1, 0)
icosahedron[10, ] <- c( -gold_ratio, 1, 0)
icosahedron[11, ] <- c( 0, gold_ratio, 1)
icosahedron[12, ] <- c( 0, gold_ratio, -1)
icosahedron <- icosahedron / sqrt( 1 + gold_ratio^2 )
icosahedron[13,] <- c(1,0,0)
icosahedron[14,] <- c(0,1,0)
icosahedron[15,] <- c(0,0,1)


G_hat_full <- function ( rejection_vec ) {
  if (sum( rejection_vec[1:15] ) == 0) {
    return(0)
  }
  if (prod( rejection_vec[1:15] ) * rejection_vec[16]  == 0) {
    accepted <- which( rejection_vec[1:15] == 1, arr.ind = TRUE)
    ico_ind <- accepted[sample( length(accepted), 1)]
    return(ico_ind)
  }
  return( 16 )
}

G_hat_smaller_with_x <- function( rejection_vec ) {
  ## first 9 components relate to \langle R_{2 pi / 11}^{u_i} \rangle, then SO(3), then SL(3) 
  if (sum( rejection_vec[1:13] ) == 0) {
    return(0)
  }
  if (prod( rejection_vec[1:13] ) * rejection_vec[16]  == 0) {
    accepted <- which( rejection_vec[1:13] == 1, arr.ind = TRUE)
    ico_ind <- accepted[sample( length(accepted), 1)]
    return(ico_ind)
  }
  return( 16 )
}

G_hat_smaller <- function( rejection_vec ) {
  ## first 9 components relate to \langle R_{2 pi / 11}^{u_i} \rangle, then SO(3), then SL(3) 
  if (sum( rejection_vec[1:12] ) == 0) {
    return(0)
  }
  if (prod( rejection_vec[1:12] ) * rejection_vec[16]  == 0) {
    accepted <- which( rejection_vec[1:12] == 1, arr.ind = TRUE)
    ico_ind <- accepted[sample( length(accepted), 1)]
    return(ico_ind)
  }
  return( 16 )
}

G_hat_smallest <- function( rejection_vec ) {
  ## first 9 components relate to \langle R_{2 pi / 11}^{u_i} \rangle, then SO(3), then SL(3) 
  if (sum( rejection_vec[1:6] ) == 0) {
    return(0)
  }
  if (prod( rejection_vec[1:6] ) * rejection_vec[16]  == 0) {
    accepted <- which( rejection_vec[1:6] == 1, arr.ind = TRUE)
    ico_ind <- accepted[sample( length(accepted), 1)]
    return(ico_ind)
  }
  return ( 16 )
}


run_sims_perm_var_varying_K <- function( f, sigma,  Sigma_X, Sigma_X_test, ns, num_sims, num_perms ) {
  G_hats_full <- matrix(0, ncol = length(ns), nrow = num_sims)
  G_hats_smaller_with_x <- matrix(0, ncol = length(ns), nrow = num_sims)
  G_hats_smaller <- matrix(0, ncol = length(ns), nrow = num_sims)
  G_hats_smallest <- matrix(0, ncol = length(ns), nrow = num_sims)
  
  MSEs_Baseline <-  matrix(0, ncol = length(ns), nrow = num_sims)
  MSEs_Full <-  matrix(0, ncol = length(ns), nrow = num_sims)
  MSEs_Smaller_with_x <-  matrix(0, ncol = length(ns), nrow = num_sims)
  MSEs_Smaller <-  matrix(0, ncol = length(ns), nrow = num_sims)
  MSEs_Smallest <-  matrix(0, ncol = length(ns), nrow = num_sims)
  
  pb <- txtProgressBar(min = 0, max = sum(ns) * num_sims, style = 3, width = 50, char = "=") 
  counter <- 0
  
  for (k in 1:length(ns)) {
    n <- ns[k]
    m <- ns[k]
    p_vals <- matrix(0, nrow = num_sims, ncol = 16)
    rejections <- matrix(0, nrow = num_sims, ncol = 16)

    for (i in 1:num_sims) {
      X <- mvrnorm(n, mu_X, Sigma_X)
      Y <- as.matrix( f(X) + rnorm(n, 0, sigma) )
      X.test <- mvrnorm(n, mu_X, Sigma_X_test) 
      Y.test <- as.matrix( f(X.test) + rnorm(n, 0, sigma) )
      
      p_vals[i, 1] <- PermVarTest_2(X, Y, 1, 1, num_perms, g_dot)
      p_vals[i, 2] <- PermVarTest_2(X, Y, 1, 2, num_perms, g_dot)
      p_vals[i, 3] <- PermVarTest_2(X, Y, 1, 3, num_perms, g_dot)
      p_vals[i, 4] <- PermVarTest_2(X, Y, 1, 4, num_perms, g_dot)
      p_vals[i, 5] <- PermVarTest_2(X, Y, 1, 5, num_perms, g_dot)
      p_vals[i, 6] <- PermVarTest_2(X, Y, 1, 6, num_perms, g_dot)
      p_vals[i, 7] <- PermVarTest_2(X, Y, 1, 7, num_perms, g_dot)
      p_vals[i, 8] <- PermVarTest_2(X, Y, 1, 8, num_perms, g_dot)
      p_vals[i, 9] <- PermVarTest_2(X, Y, 1, 9, num_perms, g_dot)
      p_vals[i, 10] <- PermVarTest_2(X, Y, 1, 10, num_perms, g_dot)
      p_vals[i, 11] <- PermVarTest_2(X, Y, 1, 11, num_perms, g_dot)
      p_vals[i, 12] <- PermVarTest_2(X, Y, 1, 12, num_perms, g_dot)
      p_vals[i, 13] <- PermVarTest_2(X, Y, 1, 13, num_perms, g_dot)
      p_vals[i, 14] <- PermVarTest_2(X, Y, 1, 14, num_perms, g_dot)
      p_vals[i, 15] <- PermVarTest_2(X, Y, 1, 15, num_perms, g_dot)
      p_vals[i, 16] <- PermVarTest_2(X, Y, 1, 16, num_perms, g_dot)
      
      rejections[i, ] <- I(p_vals[i, ] > alpha)
      # print(  rejections[i, ]  )
      
      G_hats_full[i, k] <- G_hat_full(rejections[i, ])
      G_hats_smaller_with_x[i, k] <- G_hat_smaller_with_x(rejections[i, ])
      G_hats_smaller[i, k] <- G_hat_smaller(rejections[i, ])
      G_hats_smallest[i, k] <- G_hat_smallest(rejections[i, ])
      
      MSEs_Baseline[i, k] <- find_MSE_np(X, Y, X.test, Y.test)
      
      if ( G_hats_full[i, k] == 0 ) {
        MSEs_Full[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats_full[i,k] < 16 ) {
        MSEs_Full[i,k] <- find_MSE_S1u(X, Y, X.test, Y.test, G_hats_full[i,k] )
      } else if ( G_hats_full[i, k] == 16 ) {
        MSEs_Full[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else {
        MSEs_Full[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test) 
      }
      
      if ( G_hats_smaller_with_x[i, k] == 0 ) {
        MSEs_Smaller_with_x[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats_smaller_with_x[i,k] < 16 ) {
        MSEs_Smaller_with_x[i,k] <- find_MSE_S1u(X, Y, X.test, Y.test, G_hats_smaller_with_x[i,k] )
      } else if ( G_hats_smaller_with_x[i, k] == 16 ) {
        MSEs_Smaller_with_x[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else {
        MSEs_Smaller_with_x[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test) 
      }
      
      if ( G_hats_smaller[i, k] == 0 ) {
        MSEs_Smaller[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats_smaller[i,k] < 16 ) {
        MSEs_Smaller[i,k] <- find_MSE_S1u(X, Y, X.test, Y.test, G_hats_smaller[i,k] )
      } else if ( G_hats_smaller[i, k] == 16 ) {
        MSEs_Smaller[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else {
        MSEs_Smaller[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test) 
      }
      
      if ( G_hats_smallest[i, k] == 0 ) {
        MSEs_Smallest[i, k] <- MSEs_Baseline[i, k]
      } else if ( G_hats_smallest[i,k] < 16 ) {
        MSEs_Smallest[i,k] <- find_MSE_S1u(X, Y, X.test, Y.test, G_hats_smallest[i,k] )
      } else if ( G_hats_smallest[i, k] == 16 ) {
        MSEs_Smallest[i, k] <- find_MSE_SO3(X, Y, X.test, Y.test)
      } else {
        MSEs_Smallest[i, k] <- find_MSE_SL3(X, Y, X.test, Y.test) 
      }

      counter <- counter + n
      setTxtProgressBar(pb, counter)
    }
    
    print("----------")
    print( n )
    print( mean( MSEs_Baseline[, k] ) )
    print( mean( MSEs_Full[, k] ) )
    print( mean( MSEs_Smaller_with_x[, k] ) )
    print( mean( MSEs_Smaller[, k] ) )
    print( mean( MSEs_Smallest[, k] ) )
    print( table( G_hats_full[ ,k ]))
    print( table( G_hats_smaller_with_x[ ,k ]))
    print( table( G_hats_smaller[ ,k ]))
    print( table( G_hats_smallest[ ,k ]))
  }
  
  return( cbind( G_hats_full, G_hats_smaller_with_x, G_hats_smaller, G_hats_smallest, MSEs_Baseline, MSEs_Full, MSEs_Smaller_with_x, MSEs_Smaller, MSEs_Smallest ) )
}

set.seed(1010)
scenario_2_perm_varying_K <- run_sims_perm_var_varying_K(f2, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_2_perm_varying_K, file = "scenario_2_LCE_perms_varying_K_2.csv")
set.seed(1010)
scenario_4_perm_varying_K <- run_sims_perm_var_varying_K(f4, sigma, diag(c(2,2,2)), diag(c(2,2,2)), ns, num_sims, num_perms)
write.csv(scenario_4_perm_varying_K, file = "scenario_4_LCE_perms_varying_K_2.csv")

scenario_2_perm_varying_K <- read.csv("scenario_2_LCE_perms_varying_K_2.csv")[,-1]
scenario_4_perm_varying_K <- read.csv("scenario_4_LCE_perms_varying_K_2.csv")[,-1]


scenario <- scenario_2_perm_varying_K

m <- length(ns)
MSPE_baseline <- scenario[,(4*m+1):(5*m)]
MSPE_full <- scenario[,(5*m+1):(6*m)]
MSPE_smaller_with_x <- scenario[,(6*m+1):(7*m)]
MSPE_smaller <- scenario[,(7*m+1):(8*m)]
MSPE_smallest <- scenario[,(8*m+1):(9*m)]

plot( ns, colMeans(MSPE_baseline), type = "l", ylim = c(-sigma, 0.7), 
      ylab = "MSPE", xlab = "Sample Size",  lty = 1)
lines( ns - 1, colMeans(MSPE_full), col = "blue", lty = 1)
lines( ns - 2, colMeans(MSPE_smaller_with_x), col = "red", lty = 1)
lines( ns - 2, colMeans(MSPE_smaller), col = "green", lty = 1)
lines( ns - 2, colMeans(MSPE_smallest), col = "purple", lty = 1)

MSPE_baseline.sd <- 1.96 * apply(MSPE_baseline, 2, sd) / sqrt( num_sims )
MSPE_full.sd <-  1.96 * apply(MSPE_full, 2, sd) / sqrt( num_sims )
MSPE_smaller_with_x.sd <-  1.96 * apply(MSPE_smaller_with_x, 2, sd) / sqrt( num_sims )
MSPE_smaller.sd <-  1.96 * apply(MSPE_smaller, 2, sd) / sqrt( num_sims )
MSPE_smallest.sd <-  1.96 * apply(MSPE_smallest, 2, sd) / sqrt( num_sims )

arrows(x0=ns, y0= colMeans(MSPE_baseline) - MSPE_baseline.sd - 0.01, x1=ns, y1=colMeans(MSPE_baseline) + MSPE_baseline.sd + 0.01, code=3, angle=90, length=0.05)
arrows(x0=ns - 1, y0= colMeans(MSPE_full) - MSPE_full.sd - 0.01, x1=ns - 1, y1=colMeans(MSPE_full) + MSPE_full.sd + 0.01, code=3, angle=90, length=0.05, col = "blue")
arrows(x0=ns - 2, y0= colMeans(MSPE_smaller_with_x) - MSPE_smaller_with_x.sd - 0.01, x1=ns - 2, y1=colMeans(MSPE_smaller_with_x) + MSPE_smaller_with_x.sd + 0.01, code=3, angle=90, length=0.05, col = "red")
arrows(x0=ns - 2, y0= colMeans(MSPE_smaller) - MSPE_smaller.sd - 0.01, x1=ns - 2, y1=colMeans(MSPE_smaller) + MSPE_smaller.sd + 0.01, code=3, angle=90, length=0.05, col = "green")
arrows(x0=ns - 3, y0= colMeans(MSPE_smallest) - MSPE_smallest.sd - 0.01, x1=ns - 3, y1=colMeans(MSPE_smallest) + MSPE_smallest.sd + 0.01, code=3, angle=90, length=0.05, col = "purple")


