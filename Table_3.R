library(rdist)
library(MASS)
library(np)

source("util_funcs.R")
source("test_function.R")


#### Data

mnist <- read.csv("mnist_data.csv", header = TRUE)
mnist[1:10, 1:10]



#### Setup Functions

# We are enumerating the elements of D_4 as follows:
# 1: e
# 2: R_90
# 3: R_180
# 4: R_270
# 5: R_h
# 6: R_v
# 7: R_/
# 8: R_\

mnist_g_dot <- function( D_4_index, X ) {
  
  X_mat <- matrix( unlist(X), nrow = 28, ncol = 28, byrow = TRUE )
  if (D_4_index == 1) {
    output <- X_mat
  }
  if (D_4_index == 2) {
    output <- t( apply( X_mat, 2, rev ) )
  }
  if (D_4_index == 3) {
    output <- t( apply( apply( X_mat, 2, rev ), 1, rev ) )
  }
  if (D_4_index == 4) {
    output <- apply( t(X_mat), 2, rev )
  }
  if (D_4_index == 5) {
    output <- apply( X_mat, 2, rev )
  }
  if (D_4_index == 6) {
    output <- t(apply( X_mat, 1, rev ))
  }
  if (D_4_index == 7) {
    output <- t( X_mat )
  }
  if (D_4_index == 8) {
    output <-  apply( apply( X_mat, 2, rev ), 1, rev )
  }

  return( matrix( output, nrow = 1 ) )
}

# visualise = function(vec, ...) {
#   image(matrix(as.numeric(vec),nrow = 28), col=gray((255:0)/255), ...)
# }
# old_par <- par(mfrow=c(2,4))
# for (i in 1:8) visualise( g_dot( i, mnist[1, 1:784] ) )
# par(old_par)

mu_g_1 <-           function() { return( 1 ) }
mu_g_R_h <-         function() { return( sample( c(1,5), 1 ) ) }
mu_g_R_v <-         function() { return( sample( c(1,6), 1 ) ) }
mu_g_R_slash <-     function() { return( sample( c(1,7), 1 ) ) }
mu_g_R_backslash <- function() { return( sample( c(1,8), 1 ) ) }
mu_g_R_180 <-       function() { return( sample( c(1,3), 1 ) ) }
mu_g_R_90 <-        function() { return( sample( 1:4, 1 ) ) }
mu_g_R_90_non_U <-  function() { return( sample( c(1,2), 1 ) ) }
mu_g_D_4 <-         function() { return( sample( 1:8, 1 ) ) }
mu_g_R_90_non_U_2 <-  function() { return( sample( c(1,2,4), 1 ) ) }
mu_g_R_90_poisson_power <- function() { return( 1 + (2^rpois(1,1) %% 4) ) }



mnist_0 <- subset( mnist, is0 > 0, pixel.1:pixel.784 )
mnist_6 <- subset( mnist, is6 > 0, pixel.1:pixel.784 )
mnist_9 <- subset( mnist, is9 > 0, pixel.1:pixel.784 )
merged <- rbind(mnist_0, mnist_6, mnist_9)


is_0_vec <- c( rep(1, dim(mnist_0)[1]), 
               rep(0, dim(mnist_6)[1]), 
               rep(0, dim(mnist_9)[1]) )

is_6_vec <- c( rep(0, dim(mnist_0)[1]), 
               rep(1, dim(mnist_6)[1]), 
               rep(0, dim(mnist_9)[1]) )

is_9_vec <- c( rep(0, dim(mnist_0)[1]), 
               rep(0, dim(mnist_6)[1]), 
               rep(1, dim(mnist_9)[1]) )


#### Testing

set.seed(2025)

p_vals <- matrix(0, nrow = 3, ncol = 6)

p_vals[1,1] <- PermVarTest( merged, is_0_vec, mu_g_1, 100, mnist_g_dot ) 
p_vals[1,2] <- PermVarTest( merged, is_0_vec, mu_g_R_90_non_U, 100, mnist_g_dot ) 
p_vals[1,3] <- PermVarTest( merged, is_0_vec, mu_g_R_90, 100, mnist_g_dot ) 
p_vals[1,4] <- PermVarTest( merged, is_0_vec, mu_g_R_90_non_U_2, 100, mnist_g_dot ) 
p_vals[1,5] <- PermVarTest( merged, is_0_vec, mu_g_D_4, 100, mnist_g_dot ) 
p_vals[1,6] <- PermVarTest( merged, is_0_vec, mu_g_R_90_poisson_power, 100, mnist_g_dot ) 

p_vals[2,1] <- PermVarTest( merged, is_6_vec, mu_g_1, 100, mnist_g_dot ) 
p_vals[2,2] <- PermVarTest( merged, is_6_vec, mu_g_R_90_non_U, 100, mnist_g_dot ) 
p_vals[2,3] <- PermVarTest( merged, is_6_vec, mu_g_R_90, 100, mnist_g_dot ) 
p_vals[2,4] <- PermVarTest( merged, is_6_vec, mu_g_R_90_non_U_2, 100, mnist_g_dot ) 
p_vals[2,5] <- PermVarTest( merged, is_6_vec, mu_g_D_4, 100, mnist_g_dot ) 
p_vals[2,6] <- PermVarTest( merged, is_6_vec, mu_g_R_90_poisson_power, 100, mnist_g_dot ) 

p_vals[3,1] <- PermVarTest( merged, is_9_vec, mu_g_1, 100, mnist_g_dot ) 
p_vals[3,2] <- PermVarTest( merged, is_9_vec, mu_g_R_90_non_U, 100, mnist_g_dot ) 
p_vals[3,3] <- PermVarTest( merged, is_9_vec, mu_g_R_90, 100, mnist_g_dot ) 
p_vals[3,4] <- PermVarTest( merged, is_9_vec, mu_g_R_90_non_U_2, 100, mnist_g_dot ) 
p_vals[3,5] <- PermVarTest( merged, is_9_vec, mu_g_D_4, 100, mnist_g_dot ) 
p_vals[3,6] <- PermVarTest( merged, is_9_vec, mu_g_R_90_poisson_power, 100, mnist_g_dot ) 

p_vals

