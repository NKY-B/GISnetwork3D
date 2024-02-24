##GISnetwork3D=group

##NETedges=vector line
##NETvertices=vector point
##directed=boolean FALSE
##ODMsummary=vector point

##weight=field NETedges
##mode= string 'all'

##paths=output vector

library(sf)
library(igraph)
library(dplyr)
library(purrr)
library(furrr)
library(qgisprocess)
library(GISnetwork3D)

oID = ODMsummary$oID
oV = ODMsummary$oV
dID = ODMsummary$dID
dV = ODMsummary$dV
weight = NETedges[[weight]]

if( ("fid" %in% names(NETedges)) == T ){
NETedges = NETedges %>% dplyr::select(-fid)
}

if( ("fid" %in% names(NETvertices)) == T ){
NETvertices = NETvertices %>% dplyr::select(-fid)
}

graph = GISigraph3D(edges = NETedges, vertices = NETvertices, directed = directed)

paths = Shortest_3Dpaths_for_pairs(graph = graph, edges = NETedges, oID = oID, oV = oV, dID = dID, dV = dV, weight = weight, mode = mode)


