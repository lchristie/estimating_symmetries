library(rdist)
library(MASS)
library(np)
library(ds4psy)
library(ggplot2)
# library(sphereplot)
world <- map_data("world")
library(mapproj)
library(rotasym)
library( dgof )
library(fields)

source("util_funcs.R")
source("test_function.R")
source("ico_syms.R")

#### Setup Pars

d_X <- function ( x_i, X ) {
  return( sqrt( colSums( ( t(X) - as.numeric(x_i) )^2 ) ) )
}

norm_Y <- function ( y_i, Y ) {
  return( sqrt( colSums( (y_i - t(Y) )^2 ) ) )
}

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

m_data_used <- as.matrix( m_data[inds,8:10] )

## QQ Plots

orig_dist_mat <- pdist( m_data_used )

set.seed(1010)
par( mfrow = c(3,5))


for (ind in 1:15) {
  new_X <- random_transforms_X( m_data_used, ico_syms_by_ind[[ind]], g_dot )
  new_dist_mat <- pdist( new_X )
  
  title.as <- paste( "Group Number:", ind) 
  
  plot( sort(orig_dist_mat), sort(new_dist_mat), type = "l", xlab = "", ylab = "", main = title.as )
  lines( c(0,4), c(0,4) , col = "red")
  lines( sort(orig_dist_mat), sort(new_dist_mat) )
}




