library(rdist)
library(MASS)
library(np)
library(ds4psy)
library(ggplot2)
library(sphereplot)
world <- map_data("world")
library(mapproj)
library(rotasym)
library( dgof )
library(fields)


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
  if (index < 16) {
    # Uniform angle rotations around icosahedral axis
    # angle <- rnorm( 1, 0, 0.05 )
    angle <- runif( 1, 0, 2 * pi )
    axis <- icosahedron[index, ]
    return( rotation_matrix( angle , axis ) )
  }
  if (index == 16) {
    # Samples from SO(3)
    angle <- runif( 1, min = 0, max = pi)
    axis <- rnorm(3)
    axis <- axis / sqrt( sum( axis^2 ) )
    return( rotation_matrix( angle , axis ) )
  }
  if (index > 16) {
    axis_ind <- (index - 1) %% 15
    if (axis_ind <- 0) { axis_ind <- 15 }
    angle_order <- 1 + floor( (index - 1) / 15)
    axis <- icosahedron[axis_ind, ]
    angle <- sample( angle_order, 1 )
    return( rotation_matrix( 2 * pi * angle / angle_order , axis ) )
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
    count <- 0
    for (j in 1:m) {
      if ( m_dist( X[i,], all_data[j,-4] ) < m_bandwidth ) {
        val <- val + all_data[j, dim( all_data )[2]]
        count <- count + 1
      }
    }
    output[i] <- val / count
    # print( paste( output[i], " from ", count, " data points." ) )
  }
  return( output )
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

a_q <- 1
a_B <- 200
alpha <- 0.05
set.seed(2021)

X <- as.matrix( cbind( m_data_used["x_coord"], m_data_used["y_coord"], m_data_used["z_coord"] ) )
F_NT_data <- t(t(as.matrix( m_data_used["F_NT"])))

p_vals_F <- matrix(0, nrow = 1, ncol = 15)
for (i in 1:15) {
  p_vals_F[i] <- PermVarTest_2(X, F_NT_data, a_q, i, a_B, g_dot, plot_flag = FALSE)
  print( paste( "Group index:", ind, "p value:", p_vals_F[i] ) )
}

I(p_vals_F > alpha)


## QQ Plots

orig_dist_mat <- pdist( X )

set.seed(1010)
par( mfrow = c(3,5))

for (ind in 1:15) {
  new_X <- random_transforms_X( X, ind, g_dot )
  new_dist_mat <- pdist( new_X )
  
  title.as <- paste( "Group Number:", ind) 
  
  plot( sort(orig_dist_mat), sort(new_dist_mat), type = "l", xlab = "", ylab = "", main = title.as )
  lines( c(0,4), c(0,4) , col = "red")
  lines( sort(orig_dist_mat), sort(new_dist_mat) )
  
  a <- ks.test( sort(orig_dist_mat), sort(new_dist_mat))
  
  mtext(paste( "p value = ", round( a$p.value, digits = 4 ) ), side=3)  
}

## Plots of the Magnetic Field

data_sphere <- car2sph( X[,1], X[,2], X[,3], deg = FALSE)
nColors <- 64
colindex <- as.integer(cut(m_data_used[,7],breaks=nColors))
# col_index_sample <- as.integer(cut(all_data[m_inds,4],breaks=nColors))
# cols <- two.colors(n=nColors, start="green", end="midnightblue", middle="lightseagreen", alpha=1.0)
rgl.sphgrid(radaxis = FALSE)
rgl.sphpoints(data_sphere[,"long"], data_sphere[,"lat"], radius = 1, deg = FALSE, col = tim.colors(nColors)[colindex] )


m <- ggplot(m_data_used, aes(LONG, LAT, colour = F_NT)  ) +
  labs( title = "Locations of Geomagnetic Intensity Measurements" ) +
  scale_x_continuous( "Longitude", breaks = c(-180, -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180 ), limits = c(-180,180) ) +
  scale_y_continuous( "Latitude", breaks = c(-90, -60, -30, 0, 30, 60, 90), limits = c(-75,90) ) +
  theme_bw()
m + geom_map( data = world, map = world,
              aes(x = long, y = lat, map_id = region),
              color = "black", fill = "lightgray", size = 0.1) + 
  geom_point()






#### Symmetries of Sunspots

set.seed(1010)
data("sunspots_births")
sunspots_births$X <-
  cbind(cos(sunspots_births$phi) * cos(sunspots_births$theta),
        cos(sunspots_births$phi) * sin(sunspots_births$theta),
        sin(sunspots_births$phi))

# head(sunspots_births$X)
m_cycle = 23
num_samps <- 250
inds <- which( sunspots_births$cycle == m_cycle, arr.ind = TRUE )
all_data <- cbind( sunspots_births$X[inds, ], sunspots_births$date[inds] )
m_inds <- sample( dim(all_data)[1], num_samps, replace = FALSE )
m_data <- all_data[m_inds, ]

m_data <- cbind( m_data, local_average( m_data[, -4], all_data, 0.2 ) )
plot( data_sphere[m_inds,"lat"], m_data[,5],  #col = tim.colors(nColors)[col_index_sample],
      pch = 4, xlab = "Latitude", ylab = "Time of Sunspot Occurance", yaxt="n")


a_q <- 1
a_B <- 100
alpha <- 0.05
num_syms <- 15

set.seed(2020)

p_vals_F_2 <- matrix(0, nrow = 1, ncol = num_syms)

for (ind in 1:num_syms) {
  p_vals_F_2[ind] <- PermVarTest_2( m_data[,1:3] , m_data[,5], a_q, ind, a_B, g_dot, plot_flag = FALSE)
}

p_vals_F_2
I(p_vals_F_2 >= alpha)

## Plots

data_sphere <- car2sph( m_data[,1], m_data[,2], m_data[,3], deg = FALSE)
nColors <- 64
colindex <- as.integer(cut(m_data[,5],breaks=nColors))
col_index_sample <- as.integer(cut(all_data[m_inds,4],breaks=nColors))
# cols <- two.colors(n=nColors, start="green", end="midnightblue", middle="lightseagreen", alpha=1.0)
rgl.sphgrid(radaxis = FALSE)
rgl.sphpoints(data_sphere[,"long"], data_sphere[,"lat"], radius = 1, deg = FALSE, col = tim.colors(nColors)[colindex] )

plot( data_sphere[m_inds,"lat"], all_data[m_inds,4],  col = tim.colors(nColors)[col_index_sample], pch = 4, xlab = "Latitude", ylab = "Time of Sunspot Occurance", yaxt="n")



#### QQ Plots

orig_X <- m_data[,1:3]
orig_dist_mat <- pdist( orig_X )

set.seed(1010)
par( mfrow = c(3,5))

for (ind in 1:15) {
  new_X <- random_transforms_X( orig_X, ind, g_dot )
  new_dist_mat <- pdist( new_X )
  
  title.as <- paste( "Group Number:", ind) 
  
  plot( sort(orig_dist_mat), sort(new_dist_mat), type = "l", xlab = "", ylab = "", main = title.as )
  lines( c(0,4), c(0,4) , col = "red")
  lines( sort(orig_dist_mat), sort(new_dist_mat) )

  a <- ks.test( sort(orig_dist_mat), sort(new_dist_mat))

  mtext(paste( "p value = ", round( a$p.value, digits = 4 ) ), side=3)  
}
