##GISnetwork3D=group

##NETedges=vector line
##NETvertices=vector point
##directed=boolean FALSE

##origin=vector point
##oID=field origin
##oV=field origin

##mode= string 'all'
##weight=field NETedges

##search=number 500
##costThreshold=number 500

##edgesReached=output vector
##ListWithinV=output table

library(sf)
library(igraph)
library(dplyr)
library(purrr)
library(furrr)
library(qgisprocess)
library(GISnetwork3D)

oID = origin[[oID]]
oV = origin[[oV]]

weight = NETedges[[weight]]

if( ("fid" %in% names(NETedges)) == T ){
NETedges = NETedges %>% dplyr::select(-fid)
}

if( ("fid" %in% names(NETvertices)) == T ){
NETvertices = NETvertices %>% dplyr::select(-fid)
}

graph = GISigraph3D(edges = NETedges, vertices = NETvertices, directed = directed)

ServiceArea = ServiceArea3D(graph = graph, edges = NETedges, vertices = NETvertices, oID = oID, oV = oV, mode = mode,
 weight = weight, search = search, costThreshold = costThreshold)

edgesReached = ServiceArea$edgesReached
ListWithinV = ServiceArea$ListWithinV

