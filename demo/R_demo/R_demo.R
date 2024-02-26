
### >>> This R script demonstrate the usage of "GISnetwork3D" package for 3D spatial network analysis
### >>> This demonstration selected a small community in the Eastern District on Hong Kong Island (in Hong Kong, China) as a study area and used the "GISnetwork3D" package to model the spatial accessibility to primary school from residential buildings.
### >>> In this demonstration, we will use 3D building data, a Digital Terrain Model (DTM), the location of the primary school, 3D pedestrian network, and they are preprocessed in advance. Please note that the data is for demonstration only and is not intended for formal investigation! These demonstration data are in the data folder >> the relative path to the directory is "./data/".

# This demonstration requires several R packages
library(sf)
library(raster)
library(GISnetwork3D)
library(furrr)
library(purrr)
library(qgisprocess)
library(igraph)
library(tmap)
library(scales)
###### Work flow
### ## 1 ~ Data import
### ## 2 ~ Data preprocessing
### ## 3 ~ Analysis: Origin-Destination matrix and summary
### ## 4 ~ Analysis: Deriving shortest path for each O-D pair
### ## 5 ~ Analysis: Deriving service area

#_______________________________________________________________________________________________________________________
### ## 1 ~ Data import

# Objective: Import all necessary data for the analysis
## Processes
## 1 ~ Import boundary of the study area
## 2 ~ Import building data
## 3 ~ Import 3D pedestrian network data
## 4 ~ Import the location of primary schools
## 5 ~ Import digital terrain model (DTM)
## 6 ~ Plot the geography of the study area

#_________________________________________
## 1 ~ Import boundary of the study area

# We have two boundaries. First, the studyArea is the boundary of the study area of which the buildings (or the demands) are located.
# Second, the studyAreaBuf is the boundary of the 2000m buffer of the study area. This buffer area accounts for the cross-boundary movement of the demand near the fringe of the study area. All data, except for the buildings, were delineated by this boundary. The raw data were retrieved from the District Council Constituency Areas (DCCA) census map, downloaded from an open data portal by the Hong Kong government (https://www.csdi.gov.hk/).

studyArea = st_read("data/studyArea.gpkg")
studyAreaBuf = st_read("data/studyAreaBuf.gpkg")

#_________________________________________
## 2 ~ Import building data
# We have two building data. First is the resBuild. This includes residential buildings. Second, it is otherBuild. This includes buildings other than residential buildings, such as podium structures.
# The raw data were retrieved from the iB1000 dataset, which was downloaded from an open data portal managed by the Hong Kong government (https://www.csdi.gov.hk/). The landuse of buildings is defined by the Land Utilization map 2022, downloaded from the CSDI portal (https://portal.csdi.gov.hk/geoportal/#metadataInfoPanel).

# The building data are in polygon format. Though this data contains the third dimension (z value), the z value is 0 for all polygons. The true z value is recorded in the attribute table with the field name called "BASELEVEL".

resBuild = st_read("data/resBuild.gpkg")
otherBuild = st_read("data/otherBuild.gpkg")

#_________________________________________
## 3 ~ Import 3D pedestrian network data
# This 3D pedestrian data recorded pedestrian paths in the study area in 3D. This dataset does not simply drape the network on the earth's terrain and record its surface height. This dataset records the actual location of each path. Therefore, this dataset also includes underground, footbridge, and paths inside some major shopping malls. The raw data were downloaded from an open data portal managed by the Hong Kong government (https://www.csdi.gov.hk/).
ped = st_read("./data/ped.gpkg")

#_________________________________________
## 4 ~ Import the location of primary schools
# The location of the primary school within the buffer of the study area was retrieved from the iGeoCom dataset from the open data portal managed by the Hong Kong government (https://www.csdi.gov.hk/). This data records the 2D coordinates of the points of interest. Therefore, we will need to estimate the z value from the digital terrain model (DTM).
PRiScH = st_read("./data/PRiScH.gpkg")

#_________________________________________
## 5 ~ Import digital terrain model (DTM)
# This data recorded the terrain of Hong Kong with 5m resolution in raster format. This raw dataset was obtained from the Hong Kong government (https://www.csdi.gov.hk/).
DTM = raster("./data/DTM.tif")

#_________________________________________
## 6 ~ Plot the geography of the study area
# Now, we can create a map to better understand the geography of the study area. The created map will be saved in the "Export > map".
png("./Export/map/Map01_Geography_of_study_area.png", width = 180, height = 100, units = "mm", res = 600)
tm_shape(DTM) + tm_raster(title = "DTM(m)") +
  tm_shape(studyAreaBuf %>% st_cast("MULTILINESTRING")) + tm_lines(col = "studyAreaBuf", lwd = 2, labels = "2,000m buffer of study area", title.col = "") +
  tm_shape(ped) + tm_lines(col = "ped", title.col = "", labels = "Pedestrian path", palette = "grey", lwd = 0.6 ) +
  tm_shape(studyArea %>% st_cast("MULTILINESTRING")) + tm_lines(col = "studyArea", lwd = 1.5, labels = "Study area",
                                                                title.col = "", palette = "grey40", alpha = 0.5) +
  tm_shape(otherBuild) + tm_fill(col = "otherBuild", labels = "Other buildings", title = "", palette = "grey70") +
  tm_shape(resBuild) + tm_fill(col = "resBuild", labels = "Residential buildings", title = "", palette = "blue") +
  tm_shape(PRiScH) + tm_dots(col = "PRiScH", size = 0.1, title = "", labels = "Primary school", palette = "purple", shape = 17) +
  tm_compass(position = c("RIGHT", "TOP")) + tm_scale_bar() +
  tm_layout(legend.outside = T)
dev.off()

#_______________________________________________________________________________________________________________________
### ## 2 ~ Data preprocessing

# Objective: Preprocess data for latter spatial accessibility analysis

## Processes
## 1 ~ Convert resBuild to centroid points and define the Z value
## 2 ~ Define the z value of primary school
## 3 ~ Topology correction of the pedestrian network data
## 4 ~ Create the pedestrian network by extracting the edges and vertices
## 5 ~ Identify the nearest network vertices for residential buildings and primary schools
## 6 ~ Create network graph

#_________________________________________
## 1 ~ Convert resBuild to centroid points and define the Z value
# In this study, the spatial accessibility measures the access to primary school from resBuild. We will convert the resBuild to centroid point for each polygon and then use this centroid as the origins for the spatial accessibility analysis. We will create a 3D centroid point layer.

# Convert resBuild to point and drop the z value
resBuildPT = resBuild %>% st_zm() %>% st_centroid()

# Get the coordinate of the resBuildPT and convert it to data frame
resBuildPT = resBuildPT %>% cbind( st_coordinates(resBuildPT) ) %>% st_drop_geometry()

# Create the 3D centroid point. The BASELEVEL defines the z value.
resBuildPT$Z = resBuildPT$BASELEVEL
resBuildPT = resBuildPT %>% st_as_sf(coords = c("X", "Y", "Z"), crs = 2326)

#_________________________________________
## 2 ~ Define the z value of primary school
# As the primary school data lack the Z value. Therefore, we will extract the Z value from the DTM and convert the PRiScH to 3D point features.

# Extract Z value
PRiScH$Z = raster::extract(DTM, PRiScH)

# Get the X Y coordinate of PRiScH and drop geometry
PRiScH = PRiScH %>% cbind( st_coordinates(PRiScH) ) %>% st_drop_geometry()

# Create 3D point layer
PRiScH = PRiScH %>% st_as_sf(coords = c("X", "Y", "Z"), crs = 2326)

#_________________________________________
## 3 ~ Topology correction of the pedestrian network data
# This section will show how to use the "GISnetwork3D" to conduct some geometry errors of the input 3D line layer.
# First, we will demonstrate the component identification. Second, we will demonstrate how to connect the disconnected network nodes.


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 1 ~ Component identification
# In a network, each component represents a group of connected network nodes. Each component is isolated, meaning components A and B are not connected.
#We will use the "Identify_components" function from the "GISnetwork3D" package to identify components for the ped layer. The input is the 3D line feature and the output is the input plus an additional field of "component". The component is a numeric vector. Each unique value represents the identifier of each isolated group of line segments.
ped = ped %>% Identify_components()

# Here, we use table function to count the frequency of each component and the total number of component
# frequency of each component
table(ped$component)

# total number of component. We can see that there are 47 components in the pedestrian path.
table(ped$component) %>% length()


# We can map the component to inspect what happened. We can define a group as the main component (mainComp). mainComp = 1 if the line is within the main component and 0 otherwise.
ped$mainComp = 1
ped$mainComp[ped$component > 1] = 0
tm_shape(ped) + tm_lines(col = "mainComp", style = "cat", palette = c("red", "black"))

# From the map, we can see that some components are isolated, which can be related to the clipping of the network by the buffer of the study area and also because of digitization errors, such as overshoot and undershoot.

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2 ~ Snapping disconnected segments

# Inaccurate routing results may emerge if the lines supposed to connect are disconnected. We can define a distance threshold. Then a 3D line is created to connect two disconnected vertices within this threshold.
# We can use the "snapNET3D" function from the "GISnetwork3D" package to do this task.
# We can use an example to better demonstrate the functionality

# We created two lines (LineA and LineB). Both lines have the same X and Y coordinates but the Z coordinates. The difference in Z value for the first points is 1 meter while the second points is 0.1 m. We can create a small line to connect the second points of the two lines.

LineA = data.frame( X = c(838150, 838140), Y = c(815228.3, 815227.4), Z = c(0, 0), L1 = 1) # Vertices of Line A
LineB = data.frame( X = c(838150, 838140), Y = c(815228.3, 815227.4), Z = c(1, 0.1), L1 = 2) # Vertices of Line B
lines = LineA %>% rbind(LineB) %>% st_as_sf(coords = c("X", "Y", "Z"), crs = 2326) # Combine vertices of Lines A and B
lines = qgis_run_algorithm("native:pointstopath", INPUT = lines, GROUP_EXPRESSION = "L1", OUTPUT_TEXT_DIR =  tempdir())$OUTPUT %>% st_read() # Create lines from vertices

# We can use the "snapNET3D" function to conenct vertices within a threshold 3D distance. The example belows connect vertices within 0.1m
snapLines = lines %>% snapNET3D(0.1)
print(snapLines) # We can inspect the snapLines and see two additional lines with the field of "snap" as 1. This indicates that the line is a snapping line.

snapLines %>% st_coordinates() # We can then inspect the coordinates of the lines. We can see that the two snapping lines are identical but in opposite directions. This is to avoid routing errors when we choose directional routing.


# Now we try to snap vertices within 0.5m for the pedestrian network
pedSnapped = ped %>% dplyr::select(-component, -mainComp) %>% snapNET3D(threshold = 0.5) # snap vertices within 0.5m
pedSnapped$snap %>% sum # We created 2190 snapping lines, equivalent to connecting 1095 pairs of vertices
pedSnapped = pedSnapped %>% Identify_components() # Identify components
length(table(ped$component)) - length(table(pedSnapped$component)) # Now we can compare the snapped and unsnapped ped and the snapping reduced 7 components
pedSnapped$mainComp = 1
pedSnapped$mainComp[pedSnapped$component > 1] = 0

# We can create a map to compare the distribution of disconnected components between snapped and unsnapped ped
png("./Export/map/Map02_disconnected_snapped_vs_unsnapped_comparison.png", width = 180, height = 100, units = "mm", res = 600)
tm_shape(ped %>% subset(mainComp == 0)) + tm_lines(col = "mainComp", style = "cat", palette = "red", lwd = 1.5, title.col = "", labels = "Disconnected(unsnapped)") +
  tm_shape(pedSnapped %>% subset(mainComp == 0)) + tm_lines(col = "mainComp", style = "cat", palette = "green", title.col = "", labels = "Disconnected(snapped)") +
  tm_compass(position = c("RIGHT", "TOP")) + tm_scale_bar() +
  tm_layout(legend.outside = T)
dev.off()

# >>>>> For smooth routing, we will only retain the main component for later routing
pedSnapped = pedSnapped %>% subset(component == 1)

#_________________________________________
## 4 ~ Create the pedestrian network by extracting the edges and vertices
# In network analysis, we need to have the information of the edges and nodes/vertices. Each edge connects two nodes. Edges record the weight/cost. The edges and vertices layers are the prerequisites for constructing the graph for routing and network analysis.
# Here, we create two sets of networks. One is a usual network, and the other is an anisotropic network.
# We can use the "build3DNET" function from the "GISnetwork3D" package for this task.
# The "build3DNET" will produce a list with edges and vertice layers. By default, the argument "ped" is FALSE. This means the network built has the direction of each edge based on the digitization direction. In later routing, you can still use directional or bi-directional routing. If "ped" is TRUE, a pedestrian network is created. This acknowledges the usual bidirectional nature of the pedestrian route and the direction-dependent weight, known as anisotropic movement. By enabling this argument, the output network will calculate the walking duration for each line segment, which depends on direction. Each line segment has two travel times depending on the direction. The calculation is based on Tobler's hiking function.

# Create the usual network
pedNET = pedSnapped %>% build3DNET()

# Create an anisotropic network
pedNETa = pedSnapped %>% build3DNET(ped = TRUE)


# We can print the edges layer of the output of pedNET and pedNETa. We can see that the first two columns are "from" and "to". They correspond to the direction of each line and the connected vertices. The number corresponds to the vID of the vertices layers of the output of pedNET and pedNETa. We can see that the pedNETa has two additional fields: the walking pace and the "second", which refers to the duration of walking through the edge. Other additional information produced by this function include the unique IF of each edge (eID), the slope in degree, and the length of each edge in 3D and 2D.
pedNET$edges
pedNETa$edges
pedNET$vertices
pedNETa$vertices


#_________________________________________
## 5 ~ Identify the nearest network vertices for residential buildings and primary schools

# The "findVertices3D" function from the "GISnetwork3D" package helps find the nearest network vertices for each input 3D point layer. The input is two 3D point layers; one is the point of interest, and the other is the vertice layer from the "build3DNET" function.

# Find nearest network vertices for resBuildPT
resBuildPT = resBuildPT %>% findVertices3D(pedNET$vertices)
# Find nearest network vertices for PRiScH
PRiScH = PRiScH %>% findVertices3D(pedNET$vertices)

# We can print the output, and we can see that there is an additional field "vID". This is the unique ID of the network vertex nearest to the input POI.
print(resBuildPT)
print(PRiScH)

#_________________________________________
## 6 ~ Create network graph
# The final step is to create an igraph for later routing using the "GISigraph3D" function from the "GISnetwork3D" package.
# The input is the network edges and vertices created previously.

# Create an undirected graph
pedGraph_iso = GISigraph3D(pedNET$edges, pedNET$vertices, directed = FALSE)

# Create an anisotropic graph
pedGraph_aniso = GISigraph3D(pedNETa$edges, pedNETa$vertices, directed = TRUE)

#_______________________________________________________________________________________________________________________
### ## 3 ~ Analysis: Origin-Destination matrix and summary

# This section shows how to use the "GISnetwork3D" to calculate the least-cost travel between each origin and each destination. Origin-Destination Matrix (ODM) is a data matrix recording the cost of travel between each origin and each destination. Each row represents each origin, and each column represents each destination.

# Before starting, we must give each POI a unique ID in ascending chronological order.
resBuild$ID = 1:nrow(resBuild)
resBuildPT$ID = 1:nrow(resBuild)
PRiScH$ID = 1:nrow(PRiScH)

## Processes
## 1 ~ Create OD matrix
## 2 ~ Summarize OD matrix

#_________________________________________
## 1 ~ Create OD matrix

# We will calculate the OD matrix for each residential building and each primary school for both isotropic and anisotropic routing using the "ODmatrix3D" function from the "GISnetwork3D" package.
# For better comparison, we need to calculate the travel time for the bi-directional pedestrian network. As it does not consider the direction, the slope will become useless. Therefore, we assume 0 degrees of slope and use the 3D length to calculate the travel time.

# According to Tobler's hiking function, the walking pace at 0-degree slope is 0.71m/s
pedNET$edges$second = pedNET$edges$eLength3D * 0.71


# Here, we will create three scenarios. First, the travel time between each O-D based on an undirected movement and assumes a 0-degree slope. Second is the travel time from each residential building to each school. Third is the travel time from each school to each residential building. This directional movement is defined by the argument "mode" as "in" or "out" - by default, it is "all", assuming bi-directional.

# Create ODM for undirected isotropic movement
ODMiso = ODmatrix3D(pedGraph_iso, oID = resBuildPT$ID, oV = resBuildPT$vID, dID = PRiScH$ID, dV = PRiScH$vID, weight = pedNET$edges$second)

# Create ODM for anisotropic movement
ODManiso_from =  ODmatrix3D(pedGraph_aniso, oID = resBuildPT$ID, oV = resBuildPT$vID, dID = PRiScH$ID, dV = PRiScH$vID, weight = pedNETa$edges$second, mode = "out") # ODM according to leaving from the building tower to the school

ODManiso_to =  ODmatrix3D(pedGraph_aniso, oID = resBuildPT$ID, oV = resBuildPT$vID, dID = PRiScH$ID, dV = PRiScH$vID, weight = pedNETa$edges$second, mode = "in") # ODM according to moving to the building tower from the school


#_________________________________________
## 2 ~ Summarize OD matrix

# Now, we need to summarize the OD matrix for each origin (here, it is the building tower). We can use the "ODmatrix3D_summary" function from the "GISnetwork3D" package to do this task.

# This function can calculate the mean, median travel cost, and availability within a cost threshold. In the example below, we count the number of schools within 600s walks with a weight of availability for each destination as 1.
ODMiso_summary = ODmatrix3D_summary(ODMiso, oID = resBuildPT$ID, oV = resBuildPT$vID, dID = PRiScH$ID, dV = PRiScH$vID, threshold = 600, avaiW = rep(1, nrow(PRiScH)))

ODManiso_from_summary = ODmatrix3D_summary(ODManiso_from, oID = resBuildPT$ID, oV = resBuildPT$vID, dID = PRiScH$ID, dV = PRiScH$vID, threshold = 600, avaiW = rep(1, nrow(PRiScH)))

ODManiso_to_summary = ODmatrix3D_summary(ODManiso_to, oID = resBuildPT$ID, oV = resBuildPT$vID, dID = PRiScH$ID, dV = PRiScH$vID, threshold = 600, avaiW = rep(1, nrow(PRiScH)))

# Here, we plot three scatter plots to compare travel time for the three scenarios mentioned above.
# We can see the variation in travel time between scenarios. Besides, we can also see that isotropic routing tends to have lower values than anisotropic routing. We can also see the variability in travel time between travel from and to residential buildings.
png("./Export/plot/Plot01_ODM_time_comparison.png", width = 180, height = 60, units = "mm", res = 600)

par(mfrow = c(1, 3), mar = c(4.1, 4.1, 0.5, 0.5))

plot(ODManiso_from_summary$minCost, ODManiso_to_summary$minCost, pch = 21, col = alpha("grey", 0.5), bg = alpha("red", 0.5), xlab = "Travel time (from)", ylab = "Travel time (to)", xlim = c(120, 800), ylim = c(0, 730))

plot(ODMiso_summary$minCost, ODManiso_to_summary$minCost, pch = 21, col = alpha("grey", 0.5), bg = alpha("red", 0.5), xlab = "Travel time (isotropic)", ylab = "Travel time (to)", xlim = c(120, 800), ylim = c(0, 730))

plot(ODMiso_summary$minCost, ODManiso_from_summary$minCost, pch = 21, col = alpha("grey", 0.5), bg = alpha("red", 0.5), xlab = "Travel time (isotropic)", ylab = "Travel time (from)", xlim = c(120, 800), ylim = c(0, 730))
dev.off()

#_______________________________________________________________________________________________________________________
### ## 4 ~ Analysis: Deriving shortest path for each O-D pair
# In the previous section, we summarized the ODM. The ODM summary recorded each origin's destination with the least-cost travel time. We can use this O-D pair and derive the path in the 3D line feature using the "Shortest_3Dpaths_for_pairs" function from the "GISnetwork3D" package.

## Processes
## 1 ~ Create paths
## 2 ~ Summary of paths

#_________________________________________
## 1 ~ Create paths

# Path of isotropic walking
paths_iso = Shortest_3Dpaths_for_pairs(pedGraph_iso, edges = pedNET$edges, oID = ODMiso_summary$oID, oV = ODMiso_summary$oV, dID = ODMiso_summary$dID, dV = ODMiso_summary$dV, weight = pedNET$edges$second)

# Path of walking from the origin
paths_from = Shortest_3Dpaths_for_pairs(pedGraph_aniso, edges = pedNETa$edges, oID = ODManiso_from_summary$oID, oV = ODManiso_from_summary$oV, dID = ODManiso_from_summary$dID, dV = ODManiso_from_summary$dV, weight = pedNETa$edges$second, mode = "out")

# Path of walking toward the origin
paths_to = Shortest_3Dpaths_for_pairs(pedGraph_aniso, edges = pedNETa$edges, oID = ODManiso_to_summary$oID, oV = ODManiso_to_summary$oV, dID = ODManiso_to_summary$dID, dV = ODManiso_to_summary$dV, weight = pedNETa$edges$second, mode = "in")

# We can print one of the examples. We can see that the exported data are the edges of the line layers with an additional field oID. The oID is the ID of the origin.
print(paths_to)


#_________________________________________
## 2 ~ Summary of paths

# We can summarize the exported path. For instance, we can calculate the mean walking pace or summarize environmental variables along the way.
# Here, we calculate the average walk pace.

# We use paths_to as example.
paths_toSummary = paths_to
paths_toSummary = aggregate(walkpace ~ oID, data = paths_toSummary %>% st_drop_geometry(), FUN = mean ) # average the walk pace

# We can see that some routes are having faster pace while some are slower.
summary(paths_toSummary$walkpace)

# The "Shortest_3Dpath_oneToMany" function from the "GISnetowrkGIS" is similar to the "Shortest_3Dpaths_for_pairs" while the former support one-to-one and one-to-many routing.

#_______________________________________________________________________________________________________________________
### ## 5 ~ Analysis: Deriving service area

# This section demonstrates how to use the "ServiceArea3D" function from the "GISnetwork3D" package to explore the reachable areas by an origin point.
# We will create a service area for a sampled building using both isotropic and anisotropic methods.


# 1 ~ Create sample building
sampleBuild = resBuild[1, ]
sampleBuildPT = resBuildPT[1, ]

# 2 ~ Create service area

# Path of isotropic walking
catchment_iso = ServiceArea3D(pedGraph_iso, edges = pedNET$edges, vertices = pedNET$vertices, oID = sampleBuildPT$ID, oV = sampleBuildPT$vID, weight = pedNET$edges$second, search = 2000, costThreshold = 600)

# Catchment of walking from the origin
catchment_from = ServiceArea3D(pedGraph_aniso, edges = pedNETa$edges, vertices = pedNETa$vertices, oID = sampleBuildPT$ID, oV = sampleBuildPT$vID, weight = pedNETa$edges$second, search = 2000, costThreshold = 600, mode = "out")

# Catchment of walking toward the origin
catchment_to = ServiceArea3D(pedGraph_aniso, edges = pedNETa$edges, vertices = pedNETa$vertices, oID = sampleBuildPT$ID, oV = sampleBuildPT$vID, weight = pedNETa$edges$second, search = 2000, costThreshold = 600, mode = "in")

class(catchment_to) # The output is a list.

# The output list consists of a line layer showing the reachable edges and a data frame listing all reachable vertices. This function first calculates the cost between the origin and all network vertices, and then subsets all vertices within reachable cost, and then extracts the edges containing these vertices.

# 3 ~ Plot the service area

png("./Export/map/Map03_serviceArea_sample.png", width = 180, height = 100, units = "mm", res = 600)
tm_shape(catchment_iso$edgesReached) + tm_lines() +
  tm_shape(DTM) + tm_raster(title = "DTM(m)", breaks = c(0, 20, 40, 80, 100, 150, 200, 1000) ) +
  tm_shape(ped) + tm_lines(col = "ped", title.col = "", labels = "Pedestrian path", palette = "grey", lwd = 0.6 ) +
  tm_shape(catchment_iso$edgesReached) + tm_lines(col = "oID", title.col = "", labels = "Bi-directional", palette = "grey20") +
  tm_shape(catchment_from$edgesReached) + tm_lines(col = "oID", title.col = "", labels = "Outward", palette = "red", lwd = 2) +
  tm_shape(catchment_to$edgesReached) + tm_lines(col = "oID", title.col = "", labels = "Toward", palette = "green", lwd = 0.85) +
  tm_shape(sampleBuild) + tm_fill(col = "res", title = "", labels = "Sample building", palette = "purple") +
  tm_compass(position = c("LEFT", "TOP")) + tm_scale_bar() +
  tm_layout(legend.outside = T)
dev.off()

# From the result, we can see that the elevation increases southward. The walking pace is faster for going downhill than going uphill. Therefore, the outward catchment extends northward while the inward catchment retreats southward. The isotropic travel (assuming 0 slopes) gets the greatest areal footprint.

