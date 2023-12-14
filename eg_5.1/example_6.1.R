library(rdist)
library(MASS)

#### Parameters

dimension <- 10
alpha <- 0.05
ns <- c(20,30,40,50,60,70,80,90,100,125,150,200,300,400,500,750)

mu_X <- rep(0, dimension)
Sigma_X = diag(rep(2, dimension))

sigma <- 0.05

a_B <- 50

num_sims <- 200

avt.file.name <- paste0( "sim_results_avt_f", dimension, ".csv" )
pv.file.name <- paste0( "sim_results_pv_f", dimension, ".csv" )


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

sample_G <- function( G_ind ) {
  if (G_ind == 1) {
    power = sample(2,1)
    if (power == 2) {
      return( identity_mat )
    } else {
      return( R_h )
    }
  } 
  if (G_ind == 2) {
    power = sample(2,1)
    if (power == 2) {
      return( identity_mat )
    } else {
      return( R_v )
    }
  }
  if (G_ind == 3) {
    power = sample(2,1)
    if (power == 2) {
      return( identity_mat )
    } else {
      return( R_180 )
    }
  }
  if (G_ind == 4) {
    power = sample(2,1)
    if (power == 2) {
      return( identity_mat )
    } else {
      return( R_diag_forw )
    }
  }
  if (G_ind == 5) {
    power = sample(2,1)
    if (power == 2) {
      return( identity_mat )
    } else {
      return( R_diag_back )
    }
  }
  if (G_ind == 6) {
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
}

random_transforms_X <- function ( m_X, G_ind ) {
  n <- dim(m_X)[1]
  X_out <- matrix( 0, nrow = n, ncol = dim(m_X)[2] ) 
  for (i in 1:n) {
    g_i <- sample_G( G_ind )
    X_out[i,] <- array( g_i %*% m_X[i,] )
  }
  return( X_out )
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


#### Testing Functions

AsymVarTest <- function( m_X, m_Y, d_X, norm_Y, m_V, g_ind, m_p_t, m_thresh, m_m ) {
  n <- dim(m_X)[1]
  dist_mat <- matrix(0, nrow = m_m , ncol = n)
  dif_mat <- matrix(0, nrow = m_m, ncol = n)
  
  for (i in 1:m_m) {
    g_i <- sample_G( g_ind )
    ind <- sample(1:n, 1)
    X_i <- m_X[ind,]
    Y_i <- m_Y[ind,]
    g_cdot_X_i <- g_i %*% X_i
    dist_mat[i,] <- d_X( g_cdot_X_i, m_X )
    dif_mat[i,] <- norm_Y( Y_i, m_Y)
  }
  
  inds <- row_min_inds(dist_mat)
  dists <- select_from_ind(dist_mat, inds)
  difs <- select_from_ind(dif_mat, inds)
  
  p_val <- pbinom(sum(I( m_V(dists) < difs - m_thresh) ) - 1, 
                  sum(I(difs>0)), m_p_t(m_thresh), lower.tail = FALSE )
  
  return(p_val)
}

PermVarTest <- function( m_X, m_Y, m_G_index, m_B, plot_flag = FALSE ) {
  n <- dim(m_X)[1]
  m_q <- 1
  
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
    new_X <- random_transforms_X( m_X, m_G_index )
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


#### f_d Sims

run_sims <- function( num_sims ) {

  rejections_fd_avt <- matrix(0, ncol = num_sims * length(ns), nrow = 6)
  rejections_fd_pv <- matrix(0, ncol = num_sims * length(ns), nrow = 6)
  
  pb <- txtProgressBar(min = 0, max = sum(ns) * num_sims, style = 3, width = 50, char = "=") 
  counter <- 0
  
  for (k in 1:length(ns)) {
    n <- ns[k]
    m <- ns[k]
    
    for (i in 1:num_sims) {
      X <- mvrnorm(n, mu_X, Sigma_X)
      Y <- as.matrix( f_d(X) + rnorm(n, 0, sigma) )
      
      rejections_fd_avt[1, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, 1, a_p_t, 2*sigma, m)
      rejections_fd_avt[2, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, 2, a_p_t, 2*sigma, m)
      rejections_fd_avt[3, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, 3, a_p_t, 2*sigma, m)
      rejections_fd_avt[4, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, 4, a_p_t, 2*sigma, m)
      rejections_fd_avt[5, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, 5, a_p_t, 2*sigma, m)
      rejections_fd_avt[6, (k-1)*num_sims + i] <- AsymVarTest(X, Y, d_X, norm_Y, a_V_e, 6, a_p_t, 2*sigma, m)
      
      rejections_fd_pv[1, (k-1)*num_sims + i] <- PermVarTest(X, Y, 1, a_B)
      rejections_fd_pv[2, (k-1)*num_sims + i] <- PermVarTest(X, Y, 2, a_B)
      rejections_fd_pv[3, (k-1)*num_sims + i] <- PermVarTest(X, Y, 3, a_B)
      rejections_fd_pv[4, (k-1)*num_sims + i] <- PermVarTest(X, Y, 4, a_B)
      rejections_fd_pv[5, (k-1)*num_sims + i] <- PermVarTest(X, Y, 5, a_B)
      rejections_fd_pv[6, (k-1)*num_sims + i] <- PermVarTest(X, Y, 6, a_B)

      counter <- counter + n
      setTxtProgressBar(pb, counter)
    }
    
    print(paste("--Rejections for n=", n, "--"))
    
  }
  write.csv( rejections_fd_avt, file = avt.file.name)
  write.csv( rejections_fd_pv, file = pv.file.name )
}

set.seed(1010)
# run_sims(num_sims)

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

## Confidence Region Plots

k_hat_errs_avt <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_avt <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_errs_pv <- matrix( 0, ncol = length(ns), nrow = num_sims )
k_hat_props_contain_g_max_pv <- matrix( 0, ncol = length(ns), nrow = num_sims )

for (k in 1:length(ns)) {
  for (i in 1:num_sims) {
    ind <- (k-1)*num_sims + i
    k_hat_errs_avt[i, k] <- K_hat_sym_diff( rejections_fd_avt[,ind] > alpha )
    k_hat_errs_pv[i, k] <- K_hat_sym_diff( rejections_fd_pv[,ind] > alpha )
    k_hat_props_contain_g_max_avt[i, k] <- ( sum( ( rejections_fd_avt[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 )
    k_hat_props_contain_g_max_pv[i, k] <- ( sum( ( rejections_fd_pv[,ind] > alpha ) * c(1,1,1,0,0,0) ) == 3 ) 
  }
}

plot( ns, colMeans(k_hat_errs_avt), pch = 15, type = "p", ylim = c(0,6), ylab = "Confidence Region Excess", xlab = "Sample Size (n)", log = "x")
points( ns-1, colMeans(k_hat_errs_pv), pch = 17, col = "red")
lines( ns, colMeans(k_hat_errs_avt), lty = 2)
lines( ns-1, colMeans(k_hat_errs_pv), lty = 4, col = "red")

avt.sd <- 1.96 * apply(k_hat_errs_avt, 2, sd) / sqrt( num_sims )
pv.sd <- 1.96 * apply(k_hat_errs_pv, 2, sd) / sqrt( num_sims )

arrows(x0=ns, y0= colMeans(k_hat_errs_avt) - avt.sd - 0.01, x1=ns, y1=colMeans(k_hat_errs_avt) + avt.sd + 0.01, code=3, angle=90, length=0.1)
arrows(x0=ns-1, y0= colMeans(k_hat_errs_pv) - pv.sd - 0.01, x1=ns-1, y1=colMeans(k_hat_errs_pv) + pv.sd + 0.01, code=3, angle=90, length=0.1, col = "red")


## Containment Proportion Plots

plot( ns, colMeans(k_hat_props_contain_g_max_avt), pch = 15, type = "p", 
      ylim = c(0,1), ylab = "Proportion of Correct Containment", xlab = "Sample Size (n)", log = "x")
points( ns-1, colMeans(k_hat_props_contain_g_max_pv), pch = 17, col = "red")
lines( ns, colMeans(k_hat_props_contain_g_max_avt), lty = 2)
lines( ns-1, colMeans(k_hat_props_contain_g_max_pv), lty = 4, col = "red")

avt.sd <- 1.96 * apply(k_hat_props_contain_g_max_avt, 2, sd) / sqrt( num_sims )
pv.sd <- 1.96 * apply(k_hat_props_contain_g_max_pv, 2, sd) / sqrt( num_sims )

lines( ns, rep( 0.7, length(ns )), type = "l", lty = 3, col = "blue")
lines(ns, rep( 0.95, length(ns )), type = "l", lty = 5, col = "blue")


## Estimate Accuracy Plots

G_hats_avt <- matrix(0, nrow = num_sims * length(ns), ncol = 6)
G_hats_pv <- matrix(0, nrow = num_sims * length(ns), ncol = 6)

for (k in 1:(num_sims * length(ns))) {
  # test_vec <- rejections_fd_avt[,k]
  G_hats_avt[k,] <- G_hat_from_TR( rejections_fd_avt[,k] > alpha )
  G_hats_pv[k,] <- G_hat_from_TR( rejections_fd_pv[,k] > alpha )
}

G_hat_correct_prop_avt <- rep( 0, length(ns) )
G_hat_correct_prop_pv <- rep(0, length(ns) )

for (i in 1:length(ns)) {
  m_prop_avt <- 0
  m_prop_pv <- 0
  for (k in 1:num_sims) {
    ind <- (i - 1)*num_sims + k
    if ( sum( G_hats_avt[ind,] == c(1,1,1,0,0,0) ) == 6 ) {
      m_prop_avt <- m_prop_avt + 1
    }
    if ( sum( G_hats_pv[ind,] == c(1,1,1,0,0,0) ) == 6 ) {
      m_prop_pv <- m_prop_pv + 1
    }
  }
  G_hat_correct_prop_avt[i] <- m_prop_avt / num_sims
  G_hat_correct_prop_pv[i] <- m_prop_pv / num_sims
}

print( rbind( ns, G_hat_correct_prop_avt, G_hat_correct_prop_pv)[ , c(1,4,9,12,15,16)])





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






## Rejection Proportion Plots

rej_props_fd_avt <- matrix(0, nrow = 6, ncol = length(ns) )
rej_props_fd_pv <- matrix(0, nrow = 6, ncol = length(ns) )
# 
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
}

## Used in Table 3 and 4
rej_props_fd_avt[,c(1,4,9,12,15)]
rej_props_fd_pv[,c(1,4,9,12,15)]

