#### 1 ~ dist3D - calculate the straight-line 3D distance between two point
#### 2 ~ build3DNET - Build network for 3D line vector
#### 3 ~ findVertices3D - Find the nearest network vertices for each point
#### 4 ~ GISigraph3D - Create igraph from edges and vertices
#### 5 ~ Shortest_3Dpath_oneToMany - Create paths from one origin to many destinations
#### 6 ~ ODmatrix3D - Create OD matrix
#### 7 ~ ODmatrix3D_summary - Create a summary of OD matrix
#### 8 ~ Shortest_3Dpaths_for_pairs - Create one or more paths for each pair of origin and destination
#### 9 ~ Identify_components - Identify the components of the network
#### 10 ~ snapNET3D - Connect disconnected line
#### 11 ~ ServiceArea3D - Create 3D catchment
#### 12 ~ calculate_slope - Calculate the degree slope
#### 13 ~ tobr_hf_pace - Calculate the walking speed according to the Tobler's hiking function
#### 14 ~ reverseLine_by_condition - Reverse line by condition

#___________________________________________________________________________________________________________________________
#### 1 ~ dist3D - Calculate the straight-line 3D distance between two point

dist3D = function(x1, y1, z1, x2, y2, z2){
  distance = sqrt( (x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2 )
  return(distance)
}


#___________________________________________________________________________________________________________________________
#### 2 ~ build3DNET - Build network for 3D line vector.

build3DNET = function(line, ped = F){
  if( "fid" %in% names(line) == T ){
    line = line %>% dplyr::select(-fid) # remove the fid to avoid the saving error owing to the duplicated fid after exploding the line
  }
  edges = qgis_run_algorithm("native:explodelines", INPUT = line)$OUTPUT %>% st_read() # Explode line that each line segment only has two vertices
  vertices = edges %>% st_cast("POINT") %>% dplyr::select() # Extract the vertices of each edge
  vertices = vertices %>% st_coordinates() %>% as.data.frame() %>% unique() # Get the coordinate of vertices and remove the duplicated vertices
  vertices$vID = 1:nrow(vertices) # Assign ID to vertices

  edges$eID = 1:nrow(edges) # Assign ID to each segment
  edgesCoords = edges %>% st_coordinates() # Find the coordinates of each segment
  fromX = edgesCoords[, "X"][seq(1, nrow(edgesCoords), 2)] # X of the from node
  fromY = edgesCoords[, "Y"][seq(1, nrow(edgesCoords), 2)] # Y of the from node
  fromZ = edgesCoords[, "Z"][seq(1, nrow(edgesCoords), 2)] # Z of the from node
  toX = edgesCoords[, "X"][seq(2, nrow(edgesCoords), 2)] # X of the to node
  toY = edgesCoords[, "Y"][seq(2, nrow(edgesCoords), 2)] # Y of the to node
  toZ = edgesCoords[, "Z"][seq(2, nrow(edgesCoords), 2)] # Z of the to node

  # Find the from and to vertices for each edge
  # Identify the vertices (vID) that their coordinates match the coordinates of the from and to vertices of the edges.
  edges$from = ( data.frame( coords = paste(fromX, fromY, fromZ) ) %>% left_join( data.frame(coords = paste(vertices$X, vertices$Y, vertices$Z),
                                                                                             vID = vertices$vID) ))$vID
  edges$to = ( data.frame( coords = paste(toX, toY, toZ) ) %>% left_join( data.frame(coords = paste(vertices$X, vertices$Y, vertices$Z),
                                                                                     vID = vertices$vID) ))$vID

  # Calculate the 3D and 2D length of each edge
  edges$eLength3D = dist3D(fromX, fromY, fromZ, toX, toY, toZ) # Calculate the 3D length
  edges$eLength2D = st_length(edges) %>% as.numeric() # Calculate the 2D length
  edges$slopeDeg = calculate_slope(fromX, fromY, fromZ, toX, toY, toZ)

  # Calculate the slope in degree
  edges$slopeDeg[edges$slopeDeg > 89] = 89
  edges$slopeDeg[edges$slopeDeg < -89] = -89

  if(ped == T){

    if("snap" %in% names(edges) == T){
      addEdges = edges %>% subset(snap == 0)
    } else {
      addEdges = edges}
    if( ("directCond" %in% names(addEdges)) == TRUE ){
      addEdges = addEdges %>% subset(directCond == 0)
    }
    add = list(from = addEdges$to, to = addEdges$from)
    addEdges$from = add$from
    addEdges$to = add$to
    addEdges$slopeDeg = addEdges$slopeDeg * (-1)
    addEdges = addEdges %>% st_reverse()
    edges = edges %>% rbind(addEdges)

    # Calculate the seconds
    edges$walkpace = tobr_hf_pace(edges$slopeDeg)
    edges$second = edges$walkpace * edges$eLength3D

  }

  edges = edges$from %>% cbind(edges$to) %>% cbind(edges %>% dplyr::select(-from, -to)) # Relocate the from and to to the first two columes
  names(edges)[1:2] = c("from", "to")
  vertices = st_as_sf(vertices, coords = c("X", "Y", "Z"), crs = st_crs(line) ) # Convert the vertices to be spatial point layer
  edges = edges %>% subset(from != to)
  edges$eID = 1:nrow(edges)
  return(list(edges = edges, vertices = vertices))
}

#___________________________________________________________________________________________________________________________
#### 3 ~ findVertices3D - Find the nearest network vertex for each input point

findVertices3D = function(pt, vertices){
    from = st_coordinates(pt) %>% as.data.frame() # Find the coordinates of each point of the input "pt"
    to = st_coordinates(vertices) %>% as.data.frame() # Find the coordinates of each vertex of the input "vertices"
    future::plan(multisession) # Activate the multi-threading
    vID = future_pmap( list(from$X, from$Y, from$Z),
                       function(fromX, fromY, fromZ, toX, toY, toZ, vID){
                         dist = dist3D(fromX, fromY, fromZ, toX, toY, toZ) # Calculate the one-to-many 3D distance from each point to all vertices
                         vID = vID[ dist == min(dist) ][1] # index the vertex that the one-to-many 3D distance is the minimal
                         return(vID)},
                       toX = to$X,
                       toY = to$Y,
                       toZ = to$Z,
                       vID = vertices$vID) %>% unlist()
    result = pt
    result$vID = vID
return(result)
    future:::ClusterRegistry("stop")
}

#___________________________________________________________________________________________________________________________
#### 4 ~ GISigraph3D - Create igraph from edges and vertices

GISigraph3D = function(edges, vertices, directed = FALSE){
  if(directed == FALSE){
    graph = graph_from_data_frame(edges %>% st_drop_geometry(),
                                  directed = F,
                                  vertices = vertices %>% st_drop_geometry()) # Obtain the un-directed igraph network object
  } else if(directed == TRUE){
    graph = graph_from_data_frame(edges %>% st_drop_geometry(),
                                  directed = T,
                                  vertices = vertices %>% st_drop_geometry()) # Obtain the directed igraph network object
  }
  return(graph)
}

#___________________________________________________________________________________________________________________________
#### 5 ~ Shortest_3Dpath_oneToMany - Create paths from one origin to many destinations

Shortest_3Dpath_oneToMany = function(graph, edges, weight, oV, dV, mode = "all"){
  path = shortest_paths(graph, from = oV, to = dV, weights = weight, mode = mode, output = "epath")$epath # Get the OD paths list
  dV = rep(dV, map(path, ~.x %>% length) %>% unlist) # Get the dV for each paths
  path = edges[path %>% unlist(), ] # Subset edges by paths
  path = oV %>% cbind(dV) %>% cbind(path)
  names(path)[1:2] = c("oV", "dV")
  return(path)
}

#___________________________________________________________________________________________________________________________
#### 6 ~ ODmatrix3D - Create OD matrix

# This function produce an OD matrix that calculate the cost of the least-cost paths between each origin and destination.

ODmatrix3D = function(graph, oID, oV, dID, dV, weight, mode = "all"){
  # Some desination could share the same vertices and only unique destination vertices are taken in the routing. Therefore, dID_index is to find the position of each dID in the unique list of dID.
  indexFrame = data.frame(dID = dID, dV = dV) %>% left_join( data.frame( dV = unique(dV), order = 1:length(unique(dV)) ) )
  dID_index = indexFrame$order

  future::plan(multisession) # Activate multi-threading
  # The column of ODmatrix is the destination ID and the row is the origin ID. Because the output OD matrix will remove the duplication of destination vertices. Therefore the operation below will duplicate those removed column to include the duplicated destination vertices but with unique destination ID.

  ODmatrix = future_map(oV, distances, graph = graph, to = unique(dV), mode = mode, weights = weight) # Calculate the OD matrix
  ODmatrix = ODmatrix %>%
    future_map(.,
               function(x, order){ x = as.numeric(x)[order]
               return(x)
               }, order = dID_index) %>% unlist() %>% matrix( ncol = length(dID),
                                                              byrow = T,
                                                              dimnames = list(NULL, dID) ) %>% as.data.frame()
  future:::ClusterRegistry("stop")
  return(ODmatrix)
}

#___________________________________________________________________________________________________________________________
#### 7 ~ ODmatrix3D_summary - Create a summary of OD matrix

ODmatrix3D_summary = function(frame, oID, oV, dID, dV, threshold = 500, avaiW){
  input = frame
  input[sapply(input, is.infinite)] = NA # Replace those value with "Inf" with NA
  input = input %>% split(1:nrow(input)) # Split the ODm by row for row-based operation
  result = future_map( input,
                       function(cost, dID, dV, threshold, avaiW){
                         cost = cost %>% as.numeric()
                         minCost = min(cost, na.rm = T) # Get the minimum cost
                         minCost[is.infinite(minCost) == T] = NA
                         meanCost = mean(cost, na.rm = T) # Get the mean cost
                         medianCost = median(cost, na.rm = T) # Get the median cost
                         if(is.na(minCost) == F){
                           dV = (dV[cost == minCost] %>% na.omit())[1] # Get the destination vertex of the least-cost destination, only the first matching result is subseted
                           dID = (dID[cost == minCost] %>% na.omit())[1] # Get the destination ID of the least-cost destination, only the first matching result
                         } else {
                           dV = NA
                           dID = NA
                         }
                         avai = avaiW[cost <= threshold] %>% na.omit() %>% sum() # Calculate the availability
                         result = c(dID, dV, minCost, avai, meanCost, medianCost)
                         return(result)
                       }, dID = dID, dV = dV, threshold = threshold, avaiW = avaiW )
  result = result %>% unlist() %>% matrix( ncol = 6,
                                           byrow = T,
                                           dimnames = list(NULL, c("dID", "dV", "minCost", "avai", "meanCost", "medianCost")) ) %>% as.data.frame()
  result = data.frame( oID = oID, oV = oV ) %>% cbind(result)
  return(result)
}

#___________________________________________________________________________________________________________________________
#### 8 ~ Shortest_3Dpaths_for_pairs - Create one or more paths for each pair of origin and destination

Shortest_3Dpaths_for_pairs = function(graph, edges, oID, oV, dID, dV, weight, mode = "all"){
  future::plan(multisession) # Activate muti-threading
  path = future_pmap( list(oID, oV, dV),
                      function(fromID, fromV, toV, weight, mode, graph){
                        path = shortest_paths(graph, from = fromV, to = toV, weights = weight,
                                              mode = mode, output = "epath")$epath %>% unlist() # Get the ID of the graph's edges that the shortest path pass through
                        oID = rep(fromID, length(path)) # Get a vector of oID of which the length is the same as the number of edges of the path
                        return( list( oID = oID, path = path ) )
                      }, weight = weight, mode = mode, graph = graph)
  oID = map(path, ~.x$oID) %>% unlist() # Get all of the oID index
  path = map(path, ~.x$path) %>% unlist() # Get all of the path's edges inedx
  path = edges[path, ] # Extract edges in GIS format
  path$oID = oID # Index the origin ID for each paths' edges
  future:::ClusterRegistry("stop")
  return(path)
}

#___________________________________________________________________________________________________________________________
#### 9 ~ Identify_components - Identify the components of the network

Identify_components = function(line){
  net = line %>% build3DNET # Create network
  graph = GISigraph3D(net$edges, net$vertices) # Convert to igraph object
  compV = components(graph) # Find the component of each vertices in the graph
  edges = net$edges %>% left_join(data.frame(from = net$vertices$vID, component = compV$membership)) # Identify the component of each edge by matching the component of vertices to the from vertex of each edge
  edges = edges %>% subset(from != to)
  edges = edges %>% dplyr::select(-from, -to, -eID, -eLength3D, -eLength2D, -slopeDeg)
  return(edges)
}

#___________________________________________________________________________________________________________________________
#### 10 ~ snapNET3D - Connect disconnected line

snapNET3D = function(line, threshold = 0.5){
  net = build3DNET(line) # Convert line to edges and vertices
  within = net$vertices %>% st_intersection(net$vertices %>% st_buffer(threshold) ) %>% st_drop_geometry() # Identify vertices within threshold distance in 2D for each input vertex
  names(within)[1:2] = c("from", "to")
  within = within %>% subset(from != to) #remove those having the same from and to vertices ID

  # Remove those OD vertices pairs that are already connected in the edges
  within = within[(paste(within$from, within$to) %in% paste(net$edges$from, net$edges$to) | paste(within$from, within$to) %in% paste(net$edges$to, net$edges$from)) == F, ]

  # Get the coordinates of the origin and destination vertices for each pair
  from = data.frame(vID = within$from) %>% cbind(net$vertices[within$from,] %>% st_coordinates())
  to = data.frame(vID = within$to) %>% cbind(net$vertices[within$to,] %>% st_coordinates())

  # Calculate the 3D distance between from and to for each OD pair
  dist = dist3D(from$X, from$Y, from$Z, to$X, to$Y, to$Z)

  # Subset from and to pair of which their seperation is within the threshold distance
  from = from[dist <= threshold, ]
  from$segment = 1:nrow(from)
  from$seq = 1
  to = to[dist <= threshold, ]
  to$segment = 1:nrow(to)
  to$seq = 2

  within = from %>% rbind(to) %>% st_as_sf(., coords = c("X", "Y", "Z"), crs = st_crs(net$vertices))

  # convert the identified pair into line
  within = qgis_run_algorithm("native:pointstopath",
                              INPUT = within,
                              GROUP_EXPRESSION = "segment",
                              ORDER_EXPRESSION = "seq",
                              OUTPUT_TEXT_DIR =  tempdir())$OUTPUT %>% st_read()

  # Join the snapped line into the edges
  within = within %>% dplyr::select() %>% cbind(matrix(NA, nrow = nrow(within), ncol = ncol(net$edges)-1))
  names(within)[1:(ncol(within)-1)] = names(net$edges)[1:(ncol(net$edges)-1)]
  net$edges$snap = 0
  within$snap = 1 # Create a field of "snap" of which 1 is a snapped network
  net$edges = net$edges %>% rbind(within)
  net$edges = net$edges %>% dplyr::select(-from, -to, -eID, -eLength3D, -eLength2D, -slopeDeg)
  return(net$edges)
}


#___________________________________________________________________________________________________________________________
#### 11 ~ ServiceArea3D - Create 3D catchment

# This function identifies vertices reachable within a cost threshold from an origin vertex and also all edges linking to these vertices. Therefore, the output is a list of (1) sf point showing the reachable vertices, and (2) the edges within threshold cost from each origin location.

ServiceArea3D = function(graph, edges, vertices, oID, oV, mode = "all", weight, search = 500, costThreshold = 500){

  radiusV = (vertices[oV, ] %>% st_buffer(search) %>% st_union() %>% st_as_sf() %>% st_intersection(vertices))$vID %>% unique()  # Find all vertices within search radius
  index = list( oID = rep( oID, rep( length(radiusV),  length(oV)) ),
                oV = rep( oV, rep( length(radiusV),  length(oV)) ),
                dV = rep( radiusV, length(oV) )) # dataframe of origin ID , origin vertex and all potential reachable vertices

  dist = distances(graph = graph, v = oV, to = radiusV, mode = mode, weights = weight) %>% t() %>% as.numeric() # Calculate the 3D distance

  # Subset the ID and vertex ID of the origin and the destination vertex ID reachable within a network cost
  fromID = index$oID[dist <= costThreshold]
  fromV = index$oV[dist <= costThreshold]
  toV = index$dV[dist <= costThreshold]

  ListWithinV = data.frame(oID = fromID, oV = fromV, dV = toV)

  # Subset edges that their from and to vertices contain the reachable vertices
  edgesReached = edges[edges$from %in% ListWithinV$dV & edges$to %in% ListWithinV$dV, ] %>% st_drop_geometry() %>% dplyr::select(from, to, eID)

  # Subset edges that their from and to vertices contain the reachable vertices for each input origin
  edgesReached = map( ListWithinV %>% split(ListWithinV$oID),
                      function(input, df){
                        eID = df[ df$from %in% input$dV & df$to %in% input$dV, ]$eID
                        oID = rep( unique(input$oID), length(eID) )
                        return( list(oID = oID, eID = eID) )
                      }, df =  edgesReached)
  edgesReached = edges[(map(edgesReached, ~.x$eID) %>% unlist), ] %>% cbind( data.frame(oID = map(edgesReached, ~.x$oID) %>% unlist) )
  return(list(edgesReached = edgesReached, ListWithinV = ListWithinV)) # Export edges
}

#___________________________________________________________________________________________________________________________
#### 12 ~ calculate_slope - Calculate the degree slope

calculate_slope = function(x1, y1, z1, x2, y2, z2) {

  run = sqrt((x2 - x1)^2 + (y2 - y1)^2)
  rise = z2 - z1

  degree_slope = atan(rise/run) * ( 180.0 / pi )

  return(degree_slope)
}

#___________________________________________________________________________________________________________________________
#### 13 ~ tobr_hf_pace - Calculate the walking speed according to the Tobler's hiking function

tobr_hf_pace = function(Deg){
  S = tan(Deg * pi / 180)
  result = 0.6 * exp(3.5*abs(S+0.05))
  return(result)
}

#___________________________________________________________________________________________________________________________
#### 14 ~ reverseLine_by_condition - Reverse line by condition

reverseLine_by_condition = function(line, condition = -999){
  if(((-999) %in% condition) == TRUE){
    # If no condition is specified, all lines are reversed
    result = line %>% st_reverse()
    result$reversed = 1
  } else {
    # Reverse line if condition = 1
    revLine = line[condition == 1, ] %>% st_reverse() # Subset lines for reversing
    revLine$reversed = 1
    result = line[condition == 0, ] # Subset non-reversed lines
    result$reversed = 0
    result = result %>% rbind(revLine) # combine reversed and non-reversed lines
    }
  return(result)
}














