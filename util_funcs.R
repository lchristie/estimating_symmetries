library(rdist)
library(MASS)
library(np)

gold_ratio <- ( 1 + sqrt(5) ) / 2
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

#### Helper Functions Used in Symmetry Testing

true_sample <- function (a_vec) {
  return( a_vec[ sample(length(a_vec), 1)] )
}

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

random_transforms_X <- function ( m_X, m_Sample_Function, G_action ) {
  n <- dim(m_X)[1]
  X_out <- matrix( 0, nrow = n, ncol = dim(m_X)[2] ) 
  for (i in 1:n) {
    X_out[i,] <- array( G_action(m_Sample_Function(), m_X[i,]) )
  }
  return( X_out )
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


#### For Assessing Performance in Regression

find_MSE_np <- function (X, Y, X.test, Y.test, actual = FALSE) {
  m_df <- as.data.frame( cbind(X,Y) )
  colnames(m_df) <- c( "X1", "X2", "X3", "Y" )
  m_df_T <- as.data.frame( cbind(X.test,Y.test) )
  colnames(m_df_T) <- c( "X1", "X2", "X3", "Y" )
  
  bw <- npregbw(Y ~ X1 + X2 + X3, data = m_df)
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
  
  bw <- npregbw(Y ~ X1, data = m_df)
  mod <- npreg(bw, data = m_df)
  fits <- predict( mod, data = m_df, newdata = m_df_T )
  
  if (actual) { return( mean( (fits - f(X.test))^2 ) ) }
  else { return( mean( (fits - Y.test)^2 ) ) }
}


#### For Smoothing Data with High Noise

m_dist <- function ( a_vec, another_vec ) {
  return( sqrt( sum( (a_vec - another_vec)^2 ) ) )
}

local_average <- function( X, all_data, m_bandwidth ) {
  ## This function averages over the outputs in all_data that have inputs within $h$ of the points in $X$
  
  n <- dim(X)[1]
  m <- dim( all_data )[1]
  output <- rep(0, n)
  
  for (i in 1:n) {
    val <- 0
    weights <- 0
    for (j in 1:m) {
      this_dist <- m_dist( X[i,], all_data[j,-4] )
      if ( this_dist < m_bandwidth ) {
        val <- val + all_data[j, dim( all_data )[2]] 
        weights <- weights + 1
      }
    }
    output[i] <- val / weights
    # print( paste( output[i], " from ", count, " data points." ) )
  }
  return( output )
}


g_dot <- function(g, X) {
    return( g %*% X )
}
