

#### Testing Functions

source( "util_funcs.R" )

AsymVarTest <- function( m_X, m_Y, d_X, norm_Y, m_V, m_Sample_G, m_p_t, m_thresh, m_m, m_g_dot) {
  n <- dim(m_X)[1]
  dist_mat <- matrix(0, nrow = m_m , ncol = n)
  dif_mat <- matrix(0, nrow = m_m, ncol = n)
  
  for (i in 1:m_m) {
    g_i <- m_Sample_G()
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

PermVarTest <- function ( m_X, m_Y, m_Sample_Function, m_B, m_g_dot, plot_flag = FALSE, m_lambda = 0, d_X = "euclidean") {

  
  n <- dim(m_X)[1]
  output <- 0
  
  if (m_lambda < 0) { 
    m_lambda <- Select_Lambda(m_X, m_Y, d_X)
  }

  for (i in 1:m_B) {
    ## Splitting Data Indices
    all_inds <- as.vector( c( rep( TRUE, floor(n/2) ), rep( FALSE, n - floor(n/2) ) ) )
    shuffled <- sample(all_inds)
    
    ## Baseline Stat From Original Data
    dist_mat_orig <- pdist( m_X[shuffled,] )
    if (d_X == "Great Circle") {
      dist_mat_orig <- acos( 1 - 0.5 * ( dist_mat_orig^2 ) )
    }
    dif_mat_orig <- pdist( m_Y[shuffled] )
    orig_ratios <- dif_mat_orig[which(dist_mat_orig != 0)] / ( m_lambda + dist_mat_orig[which(dist_mat_orig != 0)] )
    A_0 <- quantile(orig_ratios, probs = 1, type = 1, names = FALSE)

    
    new_X <- random_transforms_X( m_X[!shuffled,], m_Sample_Function, m_g_dot )
    dist_mat <- pdist( new_X )
    if (d_X == "Great Circle") {
      dist_mat <- acos( 1 - 0.5 * ( dist_mat^2 ) )
    }
    dif_mat <- pdist( m_Y[!shuffled] )
    ratios <- dif_mat[which(dist_mat != 0)] / (m_lambda + dist_mat[which(dist_mat != 0)] )
    A_i <- quantile(ratios, probs = 1, type = 1, names = FALSE)
    # if (plot_flag) { lines( c(0, 100), c(0, 100 * A_i), col = "red", lty = 2) }
    if (A_i <= A_0) { output <- output + 1 }
  }
  # if (plot_flag) { lines( c(0, 100), c(0, 100 * A_0) ) }
  
  if (plot_flag) {
    plot( dist_mat_orig, dif_mat_orig, type = "p", pch = 4)
    lines( c(0, 100), c(0, 100 * A_0) )
    lines( c(0, 100), c(0, 100 * A_i), col = "red", lty = 2) 
  }
  
  ## Under the alternative we would expect larger ratios for the permutations
  p_val <- output / m_B 
  return( p_val )
}


Select_Lambda <- function (m_X, m_Y, d_X = "euclidean", plot_flag = FALSE) {
  ## This function automatically selects Lambda from the data
  ## This is called from PermVarTest by using the argument m_lambda = -1
  
  dist_mat_orig <- pdist( m_X )
  if (d_X == "Great Circle") {
    dist_mat_orig <- acos( 1 - 0.5 * ( dist_mat_orig^2 ) )
  }
  dif_mat_orig <- pdist( m_Y )
  
  if (plot_flag) {
    plot( dist_mat_orig, dif_mat_orig, pch = 4, type = "p" )
  }
  
  max_ind <- which.max( dif_mat_orig[which(dist_mat_orig != 0)] / dist_mat_orig[which(dist_mat_orig != 0)] )
  X_dist_at_max <- dist_mat_orig[which(dist_mat_orig != 0)][max_ind]
  Y_dist_at_max <-  dif_mat_orig[which(dist_mat_orig != 0)][max_ind]
  
  new_ind <- which.max( 
    ( dif_mat_orig[which(dist_mat_orig > X_dist_at_max)] - Y_dist_at_max ) / 
      ( dist_mat_orig[which(dist_mat_orig > X_dist_at_max)] - X_dist_at_max ) 
  )
  
  Y_dist_at_2nd_max <- dif_mat_orig[which(dist_mat_orig > X_dist_at_max)][new_ind]
  X_dist_at_2nd_max <- dist_mat_orig[which(dist_mat_orig > X_dist_at_max)][new_ind]
  
  ideal_slope <- (Y_dist_at_2nd_max - Y_dist_at_max) / (X_dist_at_2nd_max - X_dist_at_max)
  lambda = (Y_dist_at_max / ideal_slope) - X_dist_at_max
  
  if (plot_flag) {
    print( lambda )
    points( c(X_dist_at_max, X_dist_at_2nd_max), c(Y_dist_at_max, Y_dist_at_2nd_max), pch = 4, col = "red")
    lines( c(-lambda, max(dist_mat_orig)), c(0, ideal_slope * (max(dist_mat_orig) + lambda) ), col = "red", lty = 2) 
  }

  if (lambda < 0) { return(0) }
  return( lambda )
}




