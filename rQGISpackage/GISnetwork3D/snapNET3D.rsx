##GISnetwork3D=group
##line=vector line
##threshold=number 0.5
##snapped=output vector

library(sf)
library(igraph)
library(dplyr)
library(purrr)
library(furrr)
library(qgisprocess)
library(GISnetwork3D)

result = line %>% snapNET3D(threshold = threshold)

snapped = result


