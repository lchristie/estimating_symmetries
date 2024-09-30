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

X <- as.matrix( cbind( m_data_used["x_coord"], m_data_used["y_coord"], m_data_used["z_coord"] ) )



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



