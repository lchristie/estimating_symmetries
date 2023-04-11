library(rdist)
library(MASS)
library(np)
library(ds4psy)
library(ggplot2)
world <- map_data("world")

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

sample_G <- function( index ) {
  if (index < 18) {
    # Uniform angle rotations around icosahedral axis
    angle <- rnorm(1, 0, 0.2)
    axis <- icosahedron[index, ]
    return( rotation_matrix( 2 * pi * angle , axis ) )
  }
  if (index == 18) {
    # Samples from SO(3)
    angle <- runif( 1, min = 0, max = pi)
    axis <- rnorm(3)
    axis <- axis / sqrt( sum( axis^2 ) )
    return( rotation_matrix( angle , axis ) )
  }
}

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

find_MSE_SL3 <- function (X, Y, X.test, Y.test, actual = FALSE) {
  fits <- mean( Y )
  if (actual) { return( mean( (fits - f(X.test))^2 ) ) }
  else { return( mean( (fits - Y.test)^2 ) ) }
}


#### Testing Functions

PermVarTest <- function( m_X, m_Y, d_X, norm_Y, m_V_Cal, m_q, m_G_index, m_m, m_B, m_g_dot) {
  n <- dim(m_X)[1] 
  ratios <- rep(0, m_B)
  
  pb <- txtProgressBar(min = 0, max = m_B * m_m, style = 3, width = 50, char = "=") 
  counter <- 0
  
  for (j in 1:m_B) {
    dist_mat <- matrix(0, nrow = m_m , ncol = n)
    dif_mat <- matrix(0, nrow = m_m, ncol = n)
    for (k in 1:m_m) {
      g_i <- sample_G( m_G_index )
      ind <- sample(1:n, 1)
      X_i <- m_X[ind,]
      Y_i <- m_Y[ind,]
      g_dot_X_i <- m_g_dot( g_i, X_i )
      dist_mat[k,] <- d_X( g_dot_X_i, m_X )
      dif_mat[k,] <- norm_Y(Y_i, m_Y)
      
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
    
    rats <- dif_mat[which(dist_mat != 0)] / dist_mat[which(dist_mat != 0)]
    ratios[j] <- quantile(rats, probs = m_q, type = 1, names = FALSE)  
  }
  
  dist_mat_orig <- matrix(0, nrow = m_m , ncol = n - 1)
  dif_mat_orig <- matrix(0, nrow = m_m, ncol = n - 1)
  
  for (k in 1:m_m) {
    ind <- sample(1:n, 1)
    X_i <- m_X[ind,]
    Y_i <- m_Y[ind,]
    dist_mat_orig[k,] <- d_X( X_i, m_X[-ind, ] )   # Important to exclude the original data point so we don't just see 0/0
    dif_mat_orig[k,] <- norm_Y( Y_i, m_Y[-ind,] )
  }
  
  orig_rats <- dif_mat_orig[which(dist_mat_orig != 0)] / dist_mat_orig[which(dist_mat_orig != 0)]
  orig_ratio <- quantile(orig_rats, probs = m_q, type = 1, names = FALSE)
  
  p_val <- mean( I(ratios <= orig_ratio) )
  return(p_val)
}


#### Set-up Pars

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

a_p_t <- function ( m_thres ) {
  return( 2 * exp( - (m_thres / sigma) ^2 / 4) / ( (m_thres / sigma) * sqrt(2 * pi )) )
}

## Source British Geological Survey
geo_dip_pole <- c(86.49, 162.76)
geo_mag_pole <- c(80.65, -72.68)


gold_ratio <- ( 1 + sqrt(5) ) / 2
icosahedron <- matrix( 0 , nrow = 17, ncol = 3)
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
icosahedron[13,] <- c(0,0,1)
icosahedron[14,] <- c(0,1,0)
icosahedron[15,] <- c(1,0,0)

#### Testing For Symmetries of Magnetosphere Data

## Load data
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

m_data_dense <- m_data[which(m_data$LAT > 0 & m_data$LONG > 0), ] ## Data over Eurasia
m_data_sparse <- m_data[which(m_data$LAT > 0 & m_data$LONG < 0), ] ## Data over America

set.seed(2021)
num_samples <- 200
inds <- sample( 1:(dim(m_data_dense)[1] * 1), num_samples, replace = FALSE)
inds_sparse <- sample( 1:(dim(m_data_sparse)[1] * 1), num_samples / 10, replace = FALSE)

m_data_used <- rbind( m_data_dense[inds,], m_data_sparse[inds_sparse, ])

m <- ggplot(m_data_used, aes(LONG, LAT, colour = F_NT)  ) +
  labs( title = "Locations of Geomagnetic Intensity Measurements" ) +
  scale_x_continuous( "Longitude", breaks = c(-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180 ), limits = c(-180,180) ) +
  scale_y_continuous( "Latitude", breaks = c(-90, -60, -30, 0, 30, 60, 90), limits = c(-75,90) ) +
  theme_bw()
m + geom_map( data = world, map = world,
              aes(x = long, y = lat, map_id = region),
              color = "black", fill = "lightgray", size = 0.1) + 
  geom_point()


#### Symmetry Tests

a_q <- 0.95
m <- 3000
a_B <- 200
alpha <- 0.05
set.seed(2020)

X <- as.matrix( cbind( m_data_used["x_coord"], m_data_used["y_coord"], m_data_used["z_coord"] ) )
F_NT_data <- t(t(as.matrix( m_data_used["F_NT"])))

p_vals_F <- matrix(0, nrow = 1, ncol = 18)
for (i in 1:18) {
  p_vals_F[i] <- PermVarTest(X, F_NT_data, d_X, norm_Y, a_V, a_q, i, m, a_B, g_dot)
  print( p_vals_F[i] )
}

I(p_vals_F > alpha)


#### Modelling

ref_grid <- read.csv("reference_EMF.csv", header = TRUE)
ref_grid <- cbind( ref_grid, cos( ref_grid[,"LAT"] * pi / 180 ) * cos (ref_grid[,"LONG"] * pi / 180 ) )
ref_grid <- cbind( ref_grid, cos( ref_grid[,"LAT"] * pi / 180 ) * sin (ref_grid[,"LONG"] * pi / 180 ) )
ref_grid <- cbind( ref_grid, sin( ref_grid[,"LAT"] * pi / 180 )  )
colnames(ref_grid) <- c(  colnames(ref_grid)[1:6],  "x_coord", "y_coord", "z_coord")
for (i in 1:15) {
  ref_grid <- cbind( ref_grid, acos( icosahedron[i,1] * ref_grid["x_coord"] + 
                                       icosahedron[i,2] * ref_grid["y_coord"] + 
                                       icosahedron[i,3] * ref_grid["z_coord"] ) ) 
  colnames(ref_grid)[9+i] <- paste0( "proj_coord", i) 
}

bw_F <- npregbw( F_NT ~ LAT + LONG, data = m_data_used, regtype = "ll", bwtype = "adaptive_nn")
f_F_hat <- npreg( bw_F, data = m_data_used )
fits_F <- predict( f_F_hat, data = m_data_used, newdata = ref_grid )

bw_F_S13 <- npregbw( F_NT ~ proj_coord13, data = m_data_used, regtype = "ll", bwtype = "adaptive_nn" )
f_F_S_hat13 <- npreg( bw_F_S13, data = m_data_used )
fits_F_S13 <- predict( f_F_S_hat13, data = m_data_used, newdata = ref_grid )

ref_grid <- cbind( ref_grid, fits_F, fits_F_S13 )

errors_LLE_eurasia <- mean( (ref_grid[which( (ref_grid$LAT > 0 ) & (ref_grid$LONG > 0 ) ), "F_NT"] - 
                       ref_grid[which( (ref_grid$LAT > 0 ) &  (ref_grid$LONG > 0 )), "fits_F"] )^2 ) 
errors_LLE_eurasia
errors_SLLE_eurasia <- mean( (ref_grid[which( (ref_grid$LAT > 0 ) & (ref_grid$LONG > 0 ) ), "F_NT"] - 
                       ref_grid[which( (ref_grid$LAT > 0 ) &  (ref_grid$LONG > 0 )), "fits_F_S13"] )^2 ) 
errors_SLLE_eurasia
errors_SLLE_eurasia / errors_LLE_eurasia

errors_LLE_north <- mean( (ref_grid[which( (ref_grid$LAT > 0 ) ), "F_NT"] - 
                               ref_grid[which( (ref_grid$LAT > 0 ) ), "fits_F"] )^2 ) 
errors_LLE_north
errors_SLLE_north <- mean( (ref_grid[which( (ref_grid$LAT > 0 )  ), "F_NT"] - 
                                ref_grid[which( (ref_grid$LAT > 0 )), "fits_F_S13"] )^2 ) 
errors_SLLE_north
errors_SLLE_north / errors_LLE_north

errors_LLE <- mean( (ref_grid[, "F_NT"] - ref_grid[, "fits_F"] )^2 ) 
errors_LLE

errors_SLLE13 <- mean( (ref_grid[, "F_NT"] - ref_grid[, "fits_F_S13"] )^2 ) 
errors_SLLE13


m <- ggplot(ref_grid[which(ref_grid$LAT > 0 ), ] , aes(LONG, LAT, z = fits_F) ) +
  # labs( title = "Estimated Intensity Contours" ) +
  scale_x_continuous( "Longitude", breaks = c(-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180 ), limits = c(-180,180) ) +
  scale_y_continuous( "Latitude", breaks = c(-90, -60, -30, 0, 30, 60, 90), limits = c(0,90) ) +
  theme_bw()
m + geom_map( data = world, map = world,
              aes(long, lat, z = 1, map_id = region),
              color = "black", fill = "lightgray", size = 0.1) +
  # geom_point() 
  # scale_colour_gradient2()
  geom_contour()
  
