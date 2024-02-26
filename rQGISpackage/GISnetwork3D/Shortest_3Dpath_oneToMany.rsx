##GISnetwork3D=group

##NETedges=vector line
##NETvertices=vector point
##directed=boolean FALSE

##origin=vector point
##oV=field origin

##des=vector point
##dV=field origin

##mode=string 'all'
##weight=field NETedges

##paths=output vector

library(sf)
library(igraph)
library(dplyr)
library(purrr)
library(furrr)
library(qgisprocess)
library(GISnetwork3D)

oID = origin[[oID]]
oV = origin[[oV]]

dID = des[[dID]]
dV = des[[dV]]

weight = NETedges[[weight]]

if( ("fid" %in% names(NETedges)) == T ){
NETedges = NETedges %>% dplyr::select(-fid)
}

if( ("fid" %in% names(NETvertices)) == T ){
NETvertices = NETvertices %>% dplyr::select(-fid)
}

graph = GISigraph3D(edges = NETedges, vertices = NETvertices, directed = directed)

paths = Shortest_3Dpath_oneToMany(graph = graph, edges = NETedges, oV = oV, dV = dV, 
weight = weight, mode = mode)












