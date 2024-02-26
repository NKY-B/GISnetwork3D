##GISnetwork3D=group
##line=vector line
##ped=boolean FALSE
##edges=output vector
##vertices=output vector

library(sf)
library(igraph)
library(dplyr)
library(purrr)
library(furrr)
library(qgisprocess)
library(GISnetwork3D)

result = build3DNET(line, ped = ped)

edges = result$edges
vertices = result$vertices


