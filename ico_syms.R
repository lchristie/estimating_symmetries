

#### Icosahedral Symmetries

source("util_funcs.R")

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


sample_ico_unif_axis_1 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[1, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_2 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[2, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_3 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[3, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_4 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[4, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_5 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[5, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_6 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[6, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_7 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[7, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_8 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[8, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_9 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[9, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_10 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[10, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_11 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[11, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_12 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[12, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_13 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[13, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_14 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[14, ]
  return( rotation_matrix( angle , axis ) )
}

sample_ico_unif_axis_15 <- function() {
  angle <- runif( 1, 0, 2 * pi )
  axis <- icosahedron[15, ]
  return( rotation_matrix( angle , axis ) )
}

sample_so3 <- function() {
  angle <- runif( 1, min = 0, max = pi)
  axis <- rnorm(3)
  axis <- axis / sqrt( sum( axis^2 ) )
  return( rotation_matrix( angle , axis ) )
}


ico_syms_by_ind <- c( sample_ico_unif_axis_1, 
                      sample_ico_unif_axis_2, 
                      sample_ico_unif_axis_3, 
                      sample_ico_unif_axis_4, 
                      sample_ico_unif_axis_5, 
                      sample_ico_unif_axis_6, 
                      sample_ico_unif_axis_7, 
                      sample_ico_unif_axis_8, 
                      sample_ico_unif_axis_9, 
                      sample_ico_unif_axis_10, 
                      sample_ico_unif_axis_11, 
                      sample_ico_unif_axis_12, 
                      sample_ico_unif_axis_13, 
                      sample_ico_unif_axis_14, 
                      sample_ico_unif_axis_15, 
                      sample_so3
                      )


