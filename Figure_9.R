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

data_sphere <- car2sph( m_data[,1], m_data[,2], m_data[,3], deg = FALSE)
nColors <- 64
colindex <- as.integer(cut(m_data[,5],breaks=nColors))
col_index_sample <- as.integer(cut(all_data[m_inds,4],breaks=nColors))
rgl.sphgrid(radaxis = FALSE)
rgl.sphpoints(data_sphere[,"long"], data_sphere[,"lat"], radius = 1, deg = FALSE, col = tim.colors(nColors)[colindex] )

plot( data_sphere[,"lat"], all_data[m_inds,4],  col = tim.colors(nColors)[col_index_sample], pch = 4, xlab = "Latitude", ylab = "Time of Sunspot Occurance", yaxt="n")




