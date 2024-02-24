##GISnetwork3D=group
##line=vector line
##filter=number -999
##lineComp=output vector

library(sf)
library(igraph)
library(dplyr)
library(purrr)
library(furrr)
library(qgisprocess)
library(GISnetwork3D)

# Identify component for each edges
result = line %>% Identify_components()

if(filter != -999){
result = result %>% subset(component <= filter)
}

lineComp = result


