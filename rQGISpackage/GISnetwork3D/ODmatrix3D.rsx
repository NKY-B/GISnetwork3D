##GISnetwork3D=group

##NETedges=vector line
##NETvertices=vector point
##directed=boolean FALSE

##origin=vector point
##des=vector point
##oID=string ID
##oV=string vID
##dID=string ID
##dV=string vID
##weight=field NETedges
##mode= string 'all'

##threshold=number 500
##avaiW=field des

##ODM=output table
##ODMsummary=output vector

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
avaiW = des[[avaiW]]


if( ("fid" %in% names(NETedges)) == T ){
NETedges = NETedges %>% dplyr::select(-fid)
}

if( ("fid" %in% names(NETvertices)) == T ){
NETvertices = NETvertices %>% dplyr::select(-fid)
}

graph = GISigraph3D(edges = NETedges, vertices = NETvertices, directed = directed)
ODM = ODmatrix3D(graph = graph, oID = oID, oV = oV, dID = dID, dV = dV, weight = weight, mode = mode)

ODMsummary = ODmatrix3D_summary(frame = ODM, oID = oID, oV = oV, dID = dID, dV = dV, threshold = threshold, avaiW = avaiW)
ODMsummary = origin %>% cbind(ODMsummary)




