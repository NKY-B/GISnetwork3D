##GISnetwork3D=group
##pt=vector point
##vertices=vector point
##ptWithV=output vector

library(sf)
library(igraph)
library(dplyr)
library(purrr)
library(furrr)
library(qgisprocess)
library(GISnetwork3D)

result = findVertices3D(pt = pt, vertices = vertices)

ptWithV = result


