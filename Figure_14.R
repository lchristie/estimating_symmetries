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

source("util_funcs.R")
source("test_function.R")
source("ico_syms.R")


set.seed(1010)
data("sunspots_births")
sunspots_births$X <-
  cbind(cos(sunspots_births$phi) * cos(sunspots_births$theta),
        cos(sunspots_births$phi) * sin(sunspots_births$theta),
        sin(sunspots_births$phi))

set.seed(1011)

num_samps <- 2000
a_B <- 400
a_lambda <- -0.1

m_cycle <- 23
inds <- which( sunspots_births$cycle == m_cycle, arr.ind = TRUE )
all_data <- cbind( sunspots_births$X[inds, ], sunspots_births$date[inds] )
m_inds <- sample( dim(all_data)[1], num_samps, replace = FALSE )
m_data <- all_data[m_inds, ]

m_data <- cbind( m_data, local_average( m_data[, -4], all_data, 0.1 ) )


orig_X <- m_data[,1:3]
orig_dist_mat <- acos(1 - 0.5 * pdist( orig_X )^2 )

set.seed(1010)
par( mfrow = c(3,5))


for (ind in 1:15) {
  new_X <- random_transforms_X( orig_X, ico_syms_by_ind[[ind]] , g_dot )
  new_dist_mat <- acos( 1 - 0.5 * pdist( new_X ) )
  
  title.as <- paste( "Group Number:", ind) 
  
  plot( sort(orig_dist_mat), sort(new_dist_mat), type = "l", xlab = "", ylab = "", main = title.as )
  lines( c(0,4), c(0,4) , col = "red")
  lines( sort(orig_dist_mat), sort(new_dist_mat) )
}

