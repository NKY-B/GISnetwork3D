
##_____________________________________________________________________________
#' Identify and Filter Clusters in a 3D Linestring Network
#'
#' Identifies clusters (connected components) in a network of 3D linestrings and optionally filters to retain only specified clusters.
#'
#' @param l An \code{sf} object with 3D \code{LINESTRING} geometry. The input network of linestrings.
#' @param filter A positive \code{integer} specifying the indices of clusters to retain. If \code{NULL} (default), all clusters are retained.
#'
#' @details
#' This function identifies connected components (clusters) in a network of 3D linestrings, where linestrings sharing common nodes are grouped into the same cluster. It adds a \code{cluster} field to the output, indicating the \code{numeric} cluster membership for each linestring. If \code{filter} is provided, only linestrings from the specified clusters are returned.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, including a \code{cluster} field with \code{numeric} values indicating cluster membership. If \code{filter} is specified, only linestrings from the selected clusters are included.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Create two disconnected groups of linestrings (two clusters)
#' l1 = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 2,0,0), ncol=3, byrow=TRUE))
#' ))
#' l2 = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(5,5,5, 6,6,6), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(6,6,6, 7,5,5), ncol=3, byrow=TRUE))
#' ))
#' l = rbind(l1, l2)
#'
#' # Identify clusters
#' lc = net.cluster3d(l)
#'
#' # Plot clusters in different colors
#' open3d()
#' cols = rainbow(length(unique(lc$cluster)))
#' for(i in seq_len(nrow(lc))) {
#'   coords = st_coordinates(lc[i,])
#'   lines3d(coords[,1], coords[,2], coords[,3],
#'           color = cols[lc$cluster[i]], lwd=5)
#' }
#' legend3d("topright", legend=paste("Cluster", unique(lc$cluster)), col=cols, lwd=5)
#'
#' # Filter to only cluster 2 and plot
#' lc2 = net.cluster3d(l, filter=2)
#' open3d()
#' for(i in seq_len(nrow(lc2))) {
#'   coords = st_coordinates(lc2[i,])
#'   lines3d(coords[,1], coords[,2], coords[,3], color=cols[2], lwd=7)
#' }
#' legend3d("topright", legend="Cluster 2", col=cols[2], lwd=7)
#' }
#' @export

net.cluster3d = function(l, filter = NULL){
  # --- Input validation ---
  if (!inherits(l, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(l))) {
    stop("Invalid input for l: Only accept <sf LINESTRING feature>.") }

  if (!is.null(filter)) {
    if (!is.numeric(filter) || length(filter) == 0 || any(filter <= 0) || any(filter != round(filter))) {
      stop("Invalid input for filter: Must be a vector of positive integers with length > 0.")
    }
  }
  if( "cluster" %in% names(l)){
    l = l %>% dplyr::select(-cluster)
  }

  # --- Main function implementation ---
  # Create network from the input l
  net = sf_to_net3d(l) # Generate edges and nodes
  fromNode = net$edges$FROM # Extract the start node
  net = igraph::graph_from_data_frame(net$edges %>% st_drop_geometry(),
                                      directed = F, vertices = net$nodes %>% st_drop_geometry()) # convert sf to igraph
  compV = igraph::components(net) # Find the component of each vertices in the graph
  l = cbind( l, data.frame(cluster = compV$membership[fromNode]) ) # Join the cluster number to l

  if(is.null(filter)){
    message(paste("Number of clusters:", compV$no))
    return(l)  # <---- OUTPUT: sf LINESTRING feature
  }
  else {
    message(paste("Number of clusters:", compV$no))
    message(paste("Clusters retained:", paste(filter, collapse = " ")))
    return(l[l$cluster %in% filter,])  # <---- OUTPUT: sf LINESTRING feature
  }
}

##_____________________________________________________________________________
#' Mitigate Overshoot and Undershoot in a 3D Linestring Network
#'
#' Corrects overshoot and undershoot errors in an \code{sf} object with 3D \code{LINESTRING} geometry using two strategies: (1) snapping end nodes to nearby start nodes within a tolerance distance in 3D, or (2) connecting nodes to their k-nearest neighbours within a tolerance distance in 3D.
#'
#' @param x An \code{sf} object with 3D \code{LINESTRING} geometry. The input linestring network to be corrected.
#' @param tolerance A \code{numeric} value specifying the maximum 3D distance (in the same units as the CRS of \code{x}) within which snapping or connecting is performed.
#' @param behaviour An \code{integer} value, either \code{1} or \code{2}, specifying the correction strategy:
#'   \itemize{
#'     \item \code{1}: Snaps each end node to its nearest start node if within \code{tolerance} distance.
#'     \item \code{2}: Creates new linestrings connecting each node (start or end) to its \code{k}-nearest neighbour nodes within \code{tolerance}, unless already connected.
#'   }
#' @param k A \code{numeric} \code{integer} value (default 1). The number of nearest neighbours to connect to each node when \code{behaviour = 2}. Ignored when \code{behaviour = 1}.
#'
#' @details
#' Overshoot and undershoot errors occur in digitised 3D linestring networks when lines do not meet exactly due to spatial misalignments. This function provides two correction strategies:
#' \itemize{
#'   \item \strong{Behaviour 1}: Snaps end nodes to the nearest start node within the \code{tolerance} distance, merging dangling ends while preserving Z coordinates.
#'   \item \strong{Behaviour 2}: Creates new 3D \code{LINESTRING} features to connect each node to its \code{k}-nearest neighbours within \code{tolerance}, if not already connected, to bridge minor gaps.
#' }
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, with corrected geometry (in \code{behaviour = 1}) or additional linestrings (in \code{behaviour = 2}). For \code{behaviour = 2}, a \code{snapped} column (\code{numeric} values: 1 = added linestring, 0 = original linestring) indicates new connections.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Create example: two nearly-touching lines
#' l = st_sf(geometry = st_sfc(st_linestring(matrix(c(0,0,0, 3,3,0), ncol=3, byrow=TRUE)),
#'                             st_linestring(matrix(c(3,3,1, 6,6,1), ncol=3, byrow=TRUE)) ))
#'
#' # Snap ends within 1 unit (Behaviour 1)
#' fixed1 = fix.line3d(l, tolerance = 1, behaviour = 1)
#'
#' # Explicitly connect endpoints within 1 unit (Behaviour 2, k=1)
#' fixed2 = fix.line3d(l, tolerance = 1, behaviour = 2, k = 1)
#'
#' # 3D Plot
#' open3d()
#'
#' # Plot original lines (gray)
#' lines_coords = st_coordinates(l)
#' for(i in unique(lines_coords[,"L1"])){
#'   seg = lines_coords[lines_coords[,"L1"] == i, ]
#'   lines3d(seg[,1], seg[,2], seg[,3], col="gray", lwd=2)
#' }
#'
#' # Plot snapped lines (blue)
#' fixed1_coords = st_coordinates(fixed1)
#' for(i in unique(fixed1_coords[,"L1"])){
#'   seg = fixed1_coords[fixed1_coords[,"L1"] == i, ]
#'   lines3d(seg[,1], seg[,2], seg[,3], col="blue", lwd=4)
#' }
#'
#' # Plot connections from behaviour 2 (red)
#' fixed2_coords = st_coordinates(fixed2)
#' for(i in unique(fixed2_coords[,"L1"])){
#'   seg = fixed2_coords[fixed2_coords[,"L1"] == i, ]
#'   lines3d(seg[,1], seg[,2], seg[,3], col="red", lwd=6)
#' }
#'
#'
#' # Add legend
#' legend3d("topright",
#'   legend = c("Original lines", "Snapped lines", "Connected lines"),
#'   lwd = c(2, 2, 2),
#'   col = c("gray", "blue", "red"),
#'   pch = c(NA, NA, NA)
#' )
#' }
#'
#' @export
fix.line3d = function(x, tolerance, behaviour = 1, k = 1) {

  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(x))) {
    stop("Invalid input for x: Only accept <sf LINESTRING feature>.") }

  if (!is.numeric(tolerance) || length(tolerance) > 1 || tolerance < 0) {
    stop("Invalid inputs for tolerance: Only accept ONE non-negative numeric value")}

  if(length(behaviour) != 1){stop("Invalid inputs for behaviour: Only accept ONE <integer> from c(1, 2)")}
  if(!behaviour %in% c(1,2)){stop("Invalid inputs for behaviour: Only accept ONE <integer> from c(1, 2)")}

  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != round(k)) {
    stop("Invalid input for k: Must be a positive integer.")}

  # --- Main function implementation ---
  #->>> Behaviour 1
  if (behaviour == 1) {
    # Move end node to the nearest start node if they are within the tolerance distance
    Vlist = x %>% st_coordinates() %>% as.data.frame() # Extract the vertex list
    Vlist$NodeIndex = getNodeIndex_cpp(Vlist) # Get the node index: 1 = first node; 2 = last node
    result = snap3d_p2p_cpp(Vlist[Vlist$NodeIndex == 2,], Vlist[Vlist$NodeIndex == 1,], tolerance) # Snap end nodes to start nodes
    message(paste("Number of end nodes snapped: ",
                  as.integer(result$dist != 0 & result$dist <= tolerance) %>% sum))
    Vlist[Vlist$NodeIndex == 2,][c("X","Y","Z")] = result[c("X","Y","Z")] # Replace the coordinates of nodes by the snapped nodes
    result = sfheaders::sf_linestring(Vlist, x = "X", y = "Y", z = "Z", linestring_id = "L1") # Convert dataframe to sf linestring
    st_crs(result) = st_crs(x) # Define coordinate
    st_geometry(x) = st_geometry(result) # Define coordinate
    return(x)  # <---- OUTPUT: sf LINESTRING feature
  }
  #->>> Behaviour 2
  else if (behaviour == 2) {
    # Convert x to network
    net = sf_to_net3d(x)
    # Find the k neighbour
    snappedL = knn3d(net$nodes, k = k) %>% subset(dist <= tolerance)

    # Filter out nodes that are already connected
    if(nrow(snappedL) > 0){
      id1_Node = net$nodes$nID[snappedL$id1]
      id2_Node = net$nodes$nID[snappedL$id2]
      snappedL =  snappedL[!paste(id1_Node, id2_Node) %in% paste(net$edges$FROM, net$edges$TO) &
                             !paste(id1_Node, id2_Node) %in% paste(net$edges$TO, net$edges$FROM) &
                             !paste(id2_Node, id1_Node) %in% paste(net$edges$FROM, net$edges$TO) &
                             !paste(id2_Node, id1_Node) %in% paste(net$edges$TO, net$edges$FROM), ]}

    # Join snapped line to x
    if(nrow(snappedL) > 0){
      # Create linestring to connect the nodes
      snappedL = connect3d.pair(net$nodes[snappedL$id1,], net$nodes[snappedL$id2,])

      # Create index for the new rows
      newRow = nrow(x) + 1
      newRow = newRow:(newRow - 1 + nrow(snappedL))

      x$snapped = 0 # Add field to indicate whether the line is a snapped line
      x = x %>% rbind(x[newRow,]) # Add empty rows to store result
      st_geometry(x)[newRow] = st_geometry(snappedL) # Define coordinate
      x$snapped[is.na(x$snapped)] = 1 # Join the snapped result to the original line network
      message(paste("Number of connections created: ", length(newRow)))
    } else{
      message(paste("Number of connections created: ", 0))
    }

    return(x)  # <---- OUTPUT: sf LINESTRING feature
  }
}


##_____________________________________________________________________________
#' Split 3D Linestrings at Nearest Vertices to 3D Points
#'
#' Splits 3D linestrings at the vertices closest to a set of 3D points, with an optional restriction to prohibit splitting specific linestrings.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry. The points used to identify the nearest vertices for splitting.
#' @param y An \code{sf} object with 3D \code{LINESTRING} geometry. The linestrings to be split.
#' @param restriction A \code{character} value specifying a field in \code{y} that indicates whether each linestring can be split.
#' The field must contain only \code{integer} values 0 (not restricted) or 1 (restricted). If \code{NULL} (default), no restrictions apply.
#'
#' @details
#' For each 3D point in \code{x}, the nearest vertex on the 3D linestrings in \code{y} is identified, and the linestrings are split at these vertices. Splitting occurs only at unique nearest vertices; if multiple points identify the same vertex, it is split once. If \code{restriction} is specified, linestrings with a restriction value of 1 are not split.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, containing the linestrings from \code{y} with splits at the vertices closest to the points in \code{x}, respecting any \code{restriction} field.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#' set.seed(123)
#'
#' # Example 3D linestring (zig-zag)
#' y = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1, 2,0,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(2,0,2, 3,1,3, 4,2,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(4,2,3, 4,3,3, 4,4,3, 5,4,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(5,4,3, 6,5,3, 7,6,3, 8,7,4), ncol=3, byrow=TRUE)))
#' )
#' y$restriction = c(1,0,0,1) # Assign restriction
#'
#' # Example points near the line
#' x = st_sf(geometry = st_sfc(st_point(c(1,1,2)), st_point(c(6,6,3)), st_point(c(4,4.1,3.5)), st_point(c(6,6,3))))
#'
#' # Split linestring at nearest vertices
#' splitLine = splitNearestV3d(x, y) # Without restriction
#' splitLine.restriction = splitNearestV3d(x, y, "restriction") # With restriction
#'
#' # Plot the original linestring
#' open3d()
#' y_coords = st_coordinates(y)
#' seg_col = c("#00798c", "#d1495b", "#8d96a3", "blue")
#'
#' for(i in 1:4) {
#' seg = y_coords[y_coords[,"L1"] == i, ]
#' lines3d(seg[,1], seg[,2], seg[,3], color = seg_col[i], lwd=4)
#' }
#'
#' legend3d("topright",
#'           legend = c("Line 1", "Line 2", "Line 3", "Line 4"),
#'           lwd = c(4,4,4,4),
#'           col = seg_col)
#'
#'
#' # Plot split line without restriction
#' splitLine_coords = st_coordinates(splitLine)
#' Lid = unique(splitLine_coords[,"L1"])
#' seg_col = rainbow(length(Lid))
#'
#' open3d()
#'
#' # Plot split lines (distinct colors)
#' for(i in seq_along(Lid)) {
#'   seg = splitLine_coords[splitLine_coords[,"L1"] == Lid[i], ]
#'   lines3d(seg[,1], seg[,2], seg[,3], color = seg_col[i], lwd=4)
#' }
#'
#' # Plot points (red)
#' x_coords = st_coordinates(x)
#' points3d(x_coords[,1], x_coords[,2], x_coords[,3], color="red", size=12)
#'
#' # Add a legend (showing only a few segment colors for clarity)
#' legend3d("topright",
#'   legend = c(paste0("Line ", seq_along(Lid)), "Input points"),
#'   lwd = c(rep(4, length(Lid)), NA),
#'   col = c(seg_col, "red"),
#'   pch = c(rep(NA, length(Lid)), 16)
#' )
#'
#' # Plot split line with restriction
#' splitLine_coords = st_coordinates(splitLine.restriction)
#' Lid = unique(splitLine_coords[,"L1"])
#' seg_col = rainbow(length(Lid))
#'
#' open3d()
#'
#' # Plot split lines (distinct colors)
#' for(i in seq_along(Lid)) {
#'   seg = splitLine_coords[splitLine_coords[,"L1"] == Lid[i], ]
#'   lines3d(seg[,1], seg[,2], seg[,3], color = seg_col[i], lwd=4)
#' }
#'
#' # Plot points (red)
#' x_coords = st_coordinates(x)
#' points3d(x_coords[,1], x_coords[,2], x_coords[,3], color="red", size=12)
#'
#' # Add a legend (showing only a few segment colors for clarity)
#' legend3d("topright",
#'   legend = c(paste0("Line ", seq_along(Lid)), "Input points"),
#'   lwd = c(rep(4, length(Lid)), NA),
#'   col = c(seg_col, "red"),
#'   pch = c(rep(NA, length(Lid)), 16)
#' )
#'
#'
#' }
#' @export

splitNearestV3d = function(x, y, restriction = NULL){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("Invalid input for x: Only accept <sf POINT feature>.") }

  if (!inherits(y, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(y))) {
    stop("Invalid input for y: Only accept <sf LINESTRING feature>.") }

  if(!is.null(restriction)){
    if (!is.character(restriction) || length(restriction) != 1) {
      stop("Invalid input for restriction: Must be a positive character.")}

    if (!restriction %in% names(y)) {
      stop("Invalid input for restriction: Field name is not available in y.")}

    restriction = y[["restriction"]]

    if (!all(sort(unique(restriction)) %in% c(0, 1))) {
      stop("Invalid input for restriction: The values from y's restriction field must only be 0 and 1.")} }


  # --- Main function implementation ---
  # Find the nearest line for each point
  y$L1 = 1:nrow(y) # Assign index for each linestring

  # Remove linestrings with restriction
  if(!is.null(restriction)){
    subLine = y[restriction == 0,] # Remove restricted linestrings
    subLine = subLine[nearest3d(x,subLine) %>% unique(),] # Get the unique index of the nearest linestring from subLine for each point
    subLine = y[subLine$L1 %>% unique() , ] # Get the geometry from y
  }
  else if(is.null(restriction)){
    subLine = nearest3d(x,y) # Get the index of the nearest linestring for each point
    subLine = y[subLine %>% unique(), ] # Subset those linestrings that are indexed by points
  }

  # For the subLine, find the nearest vertex for each point
  subLine = subLine %>% st_cast("POINT") # Convert subLine from line to point
  subLine$split = 0 # Assign index of split point, 0 = not to split, 1 = point to split
  subLine$split[nearest3d(x, subLine) %>% unique()] = 1 # Identify split point by indexing the nearest vertex for each point x
  subLine = subLine %>% sfheaders::sf_to_df(fill = T) # Convert the subLine from sf POINT to dataframe

  # Split the subLine by vertices
  newLine = splitLineV3d_cpp(subLine$x, subLine$y, subLine$z, subLine$L1, subLine$split)

  # Convert data.frame to sf LINESTRING grouped by the NewL1
  newLine = sfheaders::sf_linestring(obj = newLine, x = "X", y = "Y", z = "Z",
                                     linestring_id = "NewL1", keep = T) %>% dplyr::select(-NewL1)
  st_crs(newLine) = st_crs(y) # Standardise the coordinate
  # Join table

  # Join the new linestring to the original linestring
  result = y[newLine$L1,]
  st_geometry(result) = st_geometry(newLine)
  result = y[y$L1 %in% newLine$L1 == F,] %>% rbind(result)
  result = result[order(result$L1),]
  return(result)  # <---- OUTPUT: sf LINESTRING feature
}

##_____________________________________________________________________________
#' Blend 3D Points into 3D Linestrings
#'
#' Projects each 3D point onto its nearest 3D linestring and incorporates the projected points as new vertices, with an optional restriction to prohibit blending into specific linestrings.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry. The points to be blended into the linestrings.
#' @param y An \code{sf} object with 3D \code{LINESTRING} geometry. The linestrings to receive the projected points as new vertices.
#' @param restriction A \code{character} value specifying a field in \code{y} that indicates whether each linestring can be blended. The field must contain only \code{integer} values 0 (not restricted) or 1 (restricted). If \code{NULL} (default), no restrictions apply.
#'
#' @details
#' For each 3D point in \code{x}, the function identifies its nearest location (via 3D orthogonal projection) on the linestrings in \code{y} and adds these projected points as new vertices. If \code{restriction} is specified, linestrings with a restriction value of 1 are not modified.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, where projected points from \code{x} are incorporated as additional vertices in \code{y}, respecting any \code{restriction} field.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Create a 3D linestrings
#' y = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,2,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,2,1, 3,3,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(3,3,2, 5,2,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(5,2,3, 7,0,4), ncol=3, byrow=TRUE)))
#' )
#'
#' y$restriction = c(1,0,0,1)
#'
#' # Create some points near the line
#' x = st_sf(geometry = st_sfc(
#'   st_point(c(1,1,1)),
#'   st_point(c(4,2.5,3.5)),
#'   st_point(c(6,0.5,3.5))
#' )) %>% st_as_sf()
#'
#' # Blend points into linestring
#' blended = blendPT3d(x, y) # Without restriction
#' blended.restriction = blendPT3d(x, y, "restriction") # With restriction
#'
#' # Print result
#' st_coordinates(blended)
#' st_coordinates(blended.restriction)
#'
#' }
#' @export

# This is an inner function in the blendPT3d to perform blending without restriction
blendPT3d = function(x,y, restriction = NULL){

  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("Invalid input for x: Only accept <sf POINT feature>.") }

  if (!inherits(y, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(y))) {
    stop("Invalid input for y: Only accept <sf LINESTRING feature>.") }

  if(!is.null(restriction)){
    if (!is.character(restriction) || length(restriction) != 1) {
      stop("Invalid input for restriction: Must be a positive character.")}

    if (!restriction %in% names(y)) {
      stop("Invalid input for restriction: Field name is not available in y.")}

    restriction = y[["restriction"]]

    if (!all(sort(unique(restriction)) %in% c(0, 1))) {
      stop("Invalid input for restriction: The values from y's restriction field must only be 0 and 1.")}
  }

  # --- Main function implementation ---
  if(!is.null(restriction)){
    y.index = 1:nrow(y) # Get the index of y
    y.index.NoRST = y.index[restriction == 0] # Get the index of y without restriction
    index = GISnetwork3D::nearest3d(x, y[restriction == 0, ]) %>% unique() # Get the indices of the unique linestring that are the nearest from x
    index = y.index.NoRST[index] # Get the original y index from the indices
  }
  else if(is.null(restriction)){
    index = nearest3d(x, y) %>% unique() # Index the nearest linestring for point
  }

  # Subtract indexed linestirng and segmentise these lines
  segLine = y[index,] %>% dplyr::select()
  segLine$ID = index # Assign ID <--- Index of linestring from y
  segLine = segLine %>% GISnetwork3D::segmentise3d(keep = T) # Segmentise line

  # Find the nearest segment for each point
  nearestSegIndex = GISnetwork3D::nearest3d(x, segLine) # Extract the segment index
  nearestPt = GISnetwork3D::closestPT3d(x, segLine[nearestSegIndex,]) # Project x to the nearest segment and find the coordinates
  nearestPt$segID = nearestSegIndex # Assign segment ID
  nearestPt = nearestPt[order(nearestPt$segID),] %>% unique() # Remove duplicated projected points
  firstN_seg = GISnetwork3D::getNode3d(segLine[unique(nearestPt$segID),], sf = F, position = "first")[,1:3] # Get the coordinate of the first node for each segment
  lastN_seg = GISnetwork3D::getNode3d(segLine[unique(nearestPt$segID),], sf = F, position = "last")[,1:3] # Get the coordinate of the last node for each segment

  # Blend point to segment
  blendedSeg = blendPT3d_cpp(nearestPt, firstN_seg, lastN_seg, unique(nearestPt$segID))
  blendedSeg = blendedSeg[order(blendedSeg$sid, blendedSeg$dist),] # Order segment by segment ID and and distance from the first node
  blendedSeg = sfheaders::sf_linestring(blendedSeg, x = "X", y = "Y", z = "Z", linestring_id = "sid") # Create linestring by segment ID
  st_crs(blendedSeg) = st_crs(x)

  # Replace the indexed segments by the blended segments
  st_geometry(segLine)[blendedSeg$sid] = st_geometry(blendedSeg)
  segLine = segLine %>% sfheaders::sf_to_df(fill = T) # Convert back to data.frame

  # Group segments by original linestring ID
  segLine = sfheaders::sf_linestring(segLine, x = "x", y = "y", z = "z", linestring_id = "ID")
  st_crs(segLine) = st_crs(x) # Set crs

  # Replace the selected original linestring by the new blended linestrings
  st_geometry(y)[segLine$ID] = st_geometry(segLine)
  return(y)  # <---- OUTPUT: sf LINESTRING object
}

##_____________________________________________________________________________
#' Calculate Travel Time and Energy for 3D Linestrings
#'
#' Calculates the total travel time, energy expenditure, 3D length, and slope-specific metrics for each 3D linestring.
#'
#' @param l An \code{sf} object with 3D \code{LINESTRING} geometry. The input linestrings for which travel costs are calculated.
#' @param CF A \code{character} value specifying the cost function for travel time. Must be one of \code{"Tobler"} (Tobler's Hiking Function) or \code{"Campbell"} (Campbell et al., 2022).
#' @param pSlope A \code{numeric} value specifying the percent slope threshold for calculating slope-specific metrics.
#' @param a,b,c,d,e \code{numeric} values specifying parameters for Campbell's hiking function. Defaults are the 50th percentile values from Campbell et al. (2022).
#'
#' @details
#' For each 3D linestring in \code{l}, the function segments it into constituent segments (each pair of consecutive vertices). It calculates travel time, energy expenditure, and 3D length for each segment using the specified cost function:
#' \itemize{
#'   \item \code{"Tobler"}: Uses Tobler's Hiking Function for travel time.
#'   \item \code{"Campbell"}: Uses the Campbell et al. (2022) model for travel time, parameterised by \code{a}, \code{b}, \code{c}, \code{d}, and \code{e}.
#' }
#' Energy expenditure is computed using the LCDA function from Looney et al. (2019).
#' Total metrics (\code{time}, \code{energy}, \code{len3d}) are summed across all segments for each linestring.
#' The \code{pSlope} argument is used to identify segments with an absolute percent slope at or above the specified threshold, and their travel time, energy expenditure, and 3D length are summed to produce \code{slope.time}, \code{slope.energy}, and \code{slope.len3d}, respectively.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, identical to \code{l}, with six additional \code{numeric} columns:
#' \describe{
#'   \item{\code{time}}{Total travel time (in seconds) for each linestring.}
#'   \item{\code{energy}}{Total energy expenditure (in kcal/kg) for each linestring.}
#'   \item{\code{len3d}}{Total 3D length (in CRS units) for each linestring.}
#'   \item{\code{slope.time}}{Total travel time (in seconds) for segments with percent slope at or above \code{pSlope}.}
#'   \item{\code{slope.energy}}{Total energy expenditure (in kcal/kg) for segments with percent slope at or above \code{pSlope}.}
#'   \item{\code{slope.len3d}}{Total 3D length (in CRS units) for segments with percent slope at or above \code{pSlope}.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' l <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,0,0, 1,1,0), ncol=3, byrow=TRUE))
#' ))
#' # Using Tobler's Hiking Function
#' calCost.line3d(l, CF = "Tobler")
#' # Using Campbell et al. (2022) function
#' calCost.line3d(l, CF = "Campbell")
#' }
#' @export

calCost.line3d = function(l, CF = "Campbell", pSlope = 8,
                          a =  -1.4579, b = 22.0787, c = 76.3271, d = 0.0525, e = -3.2002 * 10^(-4)){

  # --- Input validation ---
  if (!inherits(l, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(l))) {
    stop("Invalid input for l: Only accept <sf LINESTRING feature>.") }

  if(length(CF) != 1){stop("Invalid input for CF: Only accept a <character> from one of `Tobler`, `Campbell`.")}
  if(!CF %in% c("Tobler", "Campbell")){stop("Invalid input for CF: Only accept a <character> from one of `Tobler`, `Campbell.")}

  if(length(pSlope)!= 1){
    stop("Invalid input for pSlope: Only accept a numeric value")
  }

  if(!is.numeric(pSlope)){
    stop("Invalid input for pSlope: Only accept a numeric value")
  }

  if(length(a)!= 1 || length(b)!= 1 || length(c)!= 1 || length(d)!= 1 || length(e)!= 1){
    stop("Invalid input for a, b, c, d, or e: Only accept one <numeric> value for each.")}
  if(!is.numeric(a) || !is.numeric(b) || !is.numeric(c) || !is.numeric(d) || !is.numeric(e)){
    stop("Invalid input for a, b, c, d, or e: Only accept one <numeric> value for each.")}

  # --- Main function implementation ---
  ## Segmentise linestring
  segLine = l %>% st_coordinates() %>% segmentise_cpp() # Segmentise line
  fromIndex = seq(1, length(segLine$L1), 2) # Index of the startNode
  toIndex = seq(2, length(segLine$L1), 2) # Index of the endNode
  segLine = data.frame(from_x = segLine$X[fromIndex], from_y = segLine$Y[fromIndex], from_z = segLine$Z[fromIndex],
                       to_x = segLine$X[toIndex], to_y = segLine$Y[toIndex], to_z = segLine$Z[toIndex],
                       L1 = segLine$L1[fromIndex], L2 = segLine$L2[fromIndex])

  ## Calculate costs
  if(CF == "Tobler"){
    return(cbind(l, calCost_net3d_THF_cpp(segLine, 1:nrow(l), pSlope)[,-1]))  # <---- OUTPUT
  } else if(CF == "Campbell"){
    return(cbind( l, calCost_net3d_CHF_cpp(segLine, 1:nrow(l), pSlope, a, b, c, d, e)[,-1]))  # <---- OUTPUT
  } else {
    stop("Invalid inputs for CF: Only accept ONE character of either <Tobler> OR <Campbell>")
  }
}

##_____________________________________________________________________________
#' Define and Manipulate Directionality of 3D Linestrings
#'
#' Modifies the directionality of 3D linestrings according to a user-specified vector to support spatial network creation.
#'
#' @param l An \code{sf} object with 3D \code{LINESTRING} geometry. The input linestrings whose directionality is to be defined or modified.
#' @param direction An \code{integer} vector of length \code{nrow(l)} with values in \code{c(-1, 0, 1)}:
#'   \itemize{
#'     \item \code{-1}: Reverse the linestring (against the digitised direction).
#'     \item \code{0}: Bi-directional; duplicates and reverses the linestring (both directions).
#'     \item \code{1}: Keep the original direction.
#'   }
#'
#' @details
#' This function enables directional control for spatial network analysis by modifying each 3D linestring in \code{l} based on the \code{direction} vector:
#' \itemize{
#'   \item \code{-1}: Reverses the linestring geometry and adds a \code{reversed} column (\code{integer}: 1 = reversed, 0 = not reversed).
#'   \item \code{0}: Duplicates the linestring, reverses the duplicate, and adds a \code{pathseq} column (\code{integer}: 1 = original, 2 = duplicated and reversed).
#'   \item \code{1}: Leaves the linestring unchanged.
#' }
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, potentially with more rows (if duplication occurs), and additional \code{integer} columns:
#' \describe{
#'   \item{\code{reversed}}{1 = reversed, 0 = not reversed.}
#'   \item{\code{pathseq}}{1 = original, 2 = duplicated and reversed (present only if any duplication was performed).}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' l <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE))
#' ))
#' # Reverse first, duplicate/reverse second
#' l2 <- dir.line3d(l, direction = c(-1, 0))
#' }
#' @export
dir.line3d = function(l, direction){

  # --- Input validation ---
  if (!inherits(l, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(l))) {
    stop("Invalid input for l: Only accept <sf LINESTRING feature>.") }

  if (length(direction) != nrow(l) || (!is.numeric(direction) & !is.integer(direction)) ) {
    stop("Invalid input for direction: It must be a <integer vector> and have the same length as l") }

  # --- Main function implementation ---
  ## Reverse linestring if -1 is present in direction
  if(-1 %in% direction){
    # Reverse linestrings where direction is -1
    l$reversed = 0 # Assign a new filed to denote whether line is reversed (0= not, 1 = reversed)
    l[direction == -1 & is.na(direction) == F,] = st_reverse(l[direction == -1 & is.na(direction) == F,]) # Reverse selected linestrings
    l$reversed[direction == -1 & is.na(direction) == F] = 1 # Indicate reversed linestrings
    message(paste("Linestrings reversed:", sum(l$reversed))) }
  else {
    message(paste("Linestrings reversed:", 0))
  }

  ## Duplicate linestring if 0 is present in direction
  if(0 %in% direction){
    # Duplicates linsteings and reverse its direction for linestring with direction = 0
    id = seq_len(nrow(l)) # Assign id for linestring
    pathseq = rep(1, nrow(l)) # Indicate the path sequence. 1 = Original linestring. 2 = Duplicated linestrings with reversed direction
    l = l %>% rbind( l[direction == 0 & is.na(direction) == F,] %>% st_reverse() ) # Duplicates linestring and reverse its direction and adds to l
    nDup = id[direction == 0 & is.na(direction) == F] %>% length() # Find the numbers of linestring duplicated
    pathseq = c(pathseq, rep(2, nDup)) # Indicates the duplicated linestring
    id = c(id, id[direction == 0 & is.na(direction) == F]) # Assign id for the added linestring
    l = l %>% cbind(pathseq)
    l = l[order(id),] # Order the new linestring by the original sequence of linestrings
    message(paste("Linestrings duplicated and reversed for di-directional linestrings:", nDup))
  } else {
    message(paste("Linestrings duplicated and reversed for di-directional linestrings:", 0))
  }
  return(l)  # <---- OUTPUT: sf LINESTRING feature
}

##_____________________________________________________________________________
#' Convert 3D Linestrings to Network Representation
#'
#' Converts a 3D linestring object into a network structure with unique nodes (as 3D points) and edges (as 3D linestrings).
#'
#' @param l An \code{sf} object with 3D \code{LINESTRING} geometry. Each feature represents an edge in the network.
#'
#' @details
#' This function extracts the start and end nodes of each 3D linestring in \code{l}, assigns unique identifiers to nodes based on their 3D coordinates (X, Y, Z), and constructs a network representation:
#' \itemize{
#'   \item \code{edges}: The original linestrings, with added columns \code{FROM} and \code{TO} (node IDs for start and end) and column \code{eID} (edge ID).
#'   \item \code{nodes}: Unique 3D points, each with a \code{character} node ID (\code{nID}).
#' }
#' The function is useful for preparing spatial network data (e.g., street or river networks) for network analysis or graph algorithms.
#'
#' @return A \code{list} with two elements:
#' \describe{
#'   \item{\code{edges}}{An \code{sf} object with 3D \code{LINESTRING} geometry, including columns \code{FROM} and \code{TO} (node IDs), a column \code{eID} (edge ID), and original attributes.}
#'   \item{\code{nodes}}{An \code{sf} object with 3D \code{POINT} geometry, including a column \code{nID} (node ID).}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' l <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(0,0,0, 2,2,2), ncol=3, byrow=TRUE))
#' ))
#' net <- sf_to_net3d(l)
#' net$edges
#' net$nodes
#' }
#'
#' @export

sf_to_net3d = function(l){
  # --- Input validation ---
  if (!inherits(l, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(l))) {
    stop("Invalid input for l: Only accept <sf LINESTRING feature>.") }

  # --- Main function implementation ---
  ## Extract nodes from each linestring
  firstN = getNode3d(l, sf = F, position = "first") # Extract first node
  endN = getNode3d(l, sf = F, position = "last") # Extract last node

  ## Get the unique node list
  uniqueN = data.frame(X = c(rbind(firstN$X, endN$X)),
                       Y = c(rbind(firstN$Y, endN$Y)),
                       Z = c(rbind(firstN$Z, endN$Z))) %>% unique()
  uniqueN$nID = paste("n", 1:nrow(uniqueN), sep = "") # Assign node ID
  firstN = firstN %>% dplyr::left_join(uniqueN, by = c("X", "Y", "Z"))
  endN = endN %>% dplyr::left_join(uniqueN, by = c("X", "Y", "Z"))
  l = cbind(firstN$nID, endN$nID, paste("e",1:nrow(l), sep = ""), l) # Find the from and to nodes, and edge ID
  names(l)[1:3] = c("FROM", "TO", "eID")
  uniqueN = sfheaders::sf_point(uniqueN, x = "X", y = "Y", z = "Z", keep = T) # Convert from data.frame to sf POINT for the node list
  st_crs(uniqueN) = st_crs(l) # Define crs
  return(list(edges = l, nodes = uniqueN))  # <---- OUTPUT: list
}

##_____________________________________________________________________________
#' Convert 3D Linestrings to Directed Network with Travel Cost Attributes
#'
#' Converts a set of 3D linestrings into a directed network structure with unique nodes (3D points) and edges (3D linestrings), including travel time and energy cost attributes.
#'
#' @param l An \code{sf} object with 3D \code{LINESTRING} geometry. Each feature represents an edge in the network.
#' @param direction \code{NULL} or a \code{character} value specifying a field in \code{l} with \code{integer} values in \code{c(-1, 0, 1)}:
#'   \itemize{
#'     \item \code{NULL}: Edges follow the digitised direction.
#'     \item \code{1}: Edge follows the digitised direction.
#'     \item \code{-1}: Edge is reversed (against digitised direction).
#'     \item \code{0}: Edge is bi-directional (duplicated and reversed).
#'   }
#' @param CF A \code{character} value specifying the cost function for travel time. Must be one of \code{"Tobler"} (Tobler's Hiking Function) or \code{"Campbell"} (Campbell et al., 2022).
#' @param force.pace \code{NULL} or a \code{character} value specifying a field in \code{l} with \code{numeric} values for travel pace.
#' If not \code{NA}, travel time is calculated as 3D length multiplied by pace, and energy is updated using standing energy expenditure.
#' @param add.time \code{NULL} or a \code{character} value specifying a field in \code{l} with \code{numeric} values for additional time (e.g., waiting time).
#' Added to travel time, with energy adjusted for the extended period.
#' @param pSlope A \code{numeric} value specifying the percent slope threshold for calculating slope-specific metrics: Travel costs spending on path with degree slope at or grater than pSlope.
#' @param a,b,c,d,e \code{numeric} values specifying parameters for Campbell's hiking function. Defaults are the 50th percentile values from Campbell et al. (2022).
#'
#' @details
#' This function converts 3D linestrings into a directed network by:
#' \itemize{
#'   \item Assigning directionality to edges based on \code{direction}, using \code{dir.line3d()} if specified.
#'   \item Extracting unique 3D nodes (start and end points) and assigning \code{character} node IDs (\code{nID}).
#'   \item Adding edge attributes: \code{FROM} and \code{TO} (node IDs), \code{eID} (edge ID), and travel cost metrics.
#'   \item Calculating travel costs by segmenting each linestring into individual segments (pairs of consecutive vertices) and applying the specified cost function (\code{Tobler} or \code{Campbell}) for travel time and the LCDA function (Looney et al., 2019) for energy expenditure. The \code{pSlope} argument is used to identify segments with an absolute percent slope at or above the specified threshold, and their travel time, energy expenditure, and 3D length are summed to produce \code{slope.time}, \code{slope.energy}, and \code{slope.len3d}, respectively.
#'   \item Adjusting travel time and energy if \code{force.pace} or \code{add.time} are specified.
#' }
#'
#' @return A \code{list} with two elements:
#' \describe{
#'   \item{\code{edges}}{An \code{sf} object with 3D \code{LINESTRING} geometry, including columns \code{FROM} and \code{TO} (node IDs), \code{eID} (edge ID), \code{time} (seconds), \code{energy} (kcal/kg), \code{len3d} (CRS units), \code{slope.time} (seconds), \code{slope.energy} (kcal/kg), and \code{slope.len3d} (CRS units), and original attributes.}
#'   \item{\code{nodes}}{An \code{sf} object with 3D \code{POINT} geometry, including column \code{nID}.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Example: simple 3D lines
#' l <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE))
#' ))
#' net <- sf_to_net3d.full(l, direction = c(1, 0), CF = "Tobler")
#' net$edges
#' net$nodes
#' }
#'
#' @export

sf_to_net3d.full = function(l, direction = NULL, CF = "Campbell", force.pace = NULL, add.time = NULL, pSlope = 8,
                            a =  -1.4579, b = 22.0787, c = 76.3271, d = 0.0525, e = -3.2002 * 10^(-4)){

  # --- Input validation ---
  if (!inherits(l, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(l))) {
    stop("Invalid input for l: Only accept <sf LINESTRING feature>.") }

  if(!is.null(direction)){
    if(!length(direction) == 1 || !is.character(direction)) {
      stop("Invalid input for direction: It must be a <character> with a length of 1") }
    if(!direction %in% names(l)){
      stop("Invalid input for direction: Missing field name from l") }

    direction = l[[direction]] }

  if(!is.null(direction)){
    if (length(direction) != nrow(l) || (!is.numeric(direction) & !is.integer(direction)) ) {
      stop("Invalid input for direction: It extracted values must be a <integer vector> and have the same length as l") }}


  if(length(CF) != 1){stop("Invalid input for CF: Only accept a <character> from one of `Tobler`, `Campbell`.")}
  if(!CF %in% c("Tobler", "Campbell")){stop("Invalid input for CF: Only accept a <character> from one of `Tobler`, `Campbell.")}

  if(!is.null(force.pace)){
    if(!length(force.pace) == 1 || !is.character(force.pace)) {
      stop("Invalid input for force.pace: It must be a <character> with a length of 1") }
    if(!force.pace %in% names(l)){
      stop("Invalid input for force.pace: Missing field name from l") }}

  if(!is.null(add.time)){
    if(!length(add.time) == 1 || !is.character(add.time)) {
      stop("Invalid input for add.time: It must be a <character> with a length of 1") }
    if(!add.time %in% names(l)){
      stop("Invalid input for add.time: Missing field name from l") }}


  if(length(a)!= 1 || length(b)!= 1 || length(c)!= 1 || length(d)!= 1 || length(e)!= 1){
    stop("Invalid input for a, b, c, d, or e: Only accept one <numeric> value for each.")}
  if(!is.numeric(a) || !is.numeric(b) || !is.numeric(c) || !is.numeric(d) || !is.numeric(e)){
    stop("Invalid input for a, b, c, d, or e: Only accept one <numeric> value for each.")}

  # --- Main function implementation ---
  ## Define the direction if the direction is not null
  if(is.null(direction) == F){
    l = dir.line3d(l = l, direction = direction)}

  ## Extract nodes from each linestring
  firstN = getNode3d(l, sf = F, position = "first") # Extract first node
  endN = getNode3d(l, sf = F, position = "last") # Extract last node

  ## Get the unique node list
  uniqueN = data.frame(X = c(rbind(firstN$X, endN$X)),
                       Y = c(rbind(firstN$Y, endN$Y)),
                       Z = c(rbind(firstN$Z, endN$Z))) %>% unique()
  uniqueN$nID = paste("n", 1:nrow(uniqueN), sep = "") # Assign node ID
  firstN = firstN %>% dplyr::left_join(uniqueN, by = c("X", "Y", "Z"))
  endN = endN %>% dplyr::left_join(uniqueN, by = c("X", "Y", "Z"))
  l = cbind(firstN$nID, endN$nID, paste("e",1:nrow(l), sep = ""), l) # Find the from and to nodes, and edge ID
  names(l)[1:3] = c("FROM", "TO", "eID")
  uniqueN = sfheaders::sf_point(uniqueN, x = "X", y = "Y", z = "Z", keep = T)
  st_crs(uniqueN) = st_crs(l)

  ## Calculate travel costs
  l = calCost.line3d(l = l, CF = CF, a = a, b = b, c = c, d = d, e = e)

  # Re-calculate travel time if pace is specified
  if(!is.null(force.pace)){
    force.pace = l[[force.pace]]
    l$time[!is.na(force.pace)] = GISnetwork3D::len3d(l[!is.na(force.pace),]) # Calculate the 3D length
    l$time[!is.na(force.pace)] = l$time[!is.na(force.pace)] * force.pace[!is.na(force.pace)] # Multiply the 3D length by the force.pace
    l$energy[!is.na(force.pace)] = l$time[!is.na(force.pace)] * (1.44 * 2.39 * 10^-4) # Update Energy
  }

  # Re-calculate travel time if add.time is specified
  if(!is.null(add.time)){
    add.time = l[[add.time]] # Extract the add time
    l$time = l$time + add.time # Add time
    l$energy = l$energy + (add.time * (1.44 * 2.39 * 10^-4)) # Update Energy
  }

  return(list(edges = l, nodes = uniqueN))  # <---- OUTPUT: list
}


##_____________________________________________________________________________
#' Subset a 3D Network by Node Index
#'
#' Subsets a 3D network by node indices, returning selected nodes and edges in a specified format.
#'
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically output from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param v A \code{character} vector of node indices (\code{nID}) to subset by.
#' @param exact A \code{logical} value.
#' If \code{TRUE}, only edges where both \code{FROM} and \code{TO} nodes are in \code{v} are kept.
#' If \code{FALSE} (default), edges are kept if either \code{FROM} or \code{TO} is in \code{v}.
#' @param output A \code{character} value specifying the output format:
#' \code{"sf"} (default, \code{sf} objects), \code{"df"} (\code{data.frame}s), or \code{"graph"} (\code{igraph} object).
#'
#' @details
#' This function subsets a 3D network (nodes as 3D points, edges as 3D linestrings) by selecting nodes specified in \code{v}:
#' \itemize{
#'   \item If \code{exact = FALSE} (default), retains edges where either \code{FROM} or \code{TO} is in \code{v}, and includes all nodes referenced by those edges.
#'   \item If \code{exact = TRUE}, retains only edges where both \code{FROM} and \code{TO} are in \code{v}, and includes only those nodes.
#' }
#' The output format is determined by \code{output}. For \code{"df"} or \code{"graph"}, edge geometries are dropped.
#'
#' @return Depends on \code{output}:
#' \describe{
#'   \item{\code{"sf"}}{A \code{list} with \code{nodes} (3D \code{sf} \code{POINT} objects with \code{nID}) and \code{edges} (3D \code{sf} \code{LINESTRING} objects with \code{FROM}, \code{TO}, \code{eID}, and original attributes).}
#'   \item{\code{"df"}}{A \code{list} with \code{nodes} (\code{data.frame} with \code{nID}) and \code{edges} (\code{data.frame} with \code{FROM}, \code{TO}, \code{eID}, and original attributes, without geometry).}
#'   \item{\code{"graph"}}{An \code{igraph} object representing the directed subnetwork.}
#' }
#'
#' @examples
#' \dontrun{
#' # Create a line network
#' l <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(2,2,2, 3,3,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(3,3,3, 4,4,4), ncol=3, byrow=TRUE))
#' ))
#' # Create network
#' net = GISnetwork3D::sf_to_net3d(l)
#' net_sub = subset.net3d(net, v = c("n2", "n3", "n4"), exact = TRUE, output = "sf")
#' # Return as igraph object for network analysis:
#' g = subset.net3d(net, v = c("n2", "n3", "n4"), exact = FALSE, output = "graph")
#' }
#' @export
subset.net3d = function(net, v, exact = F, output = "sf"){

  # --- Input validation ---
  if (!inherits(net, "list") || !"nodes" %in% names(net) || !"edges" %in% names(net)) {
    stop("Invalid input for net: Only accepts <list> of network containing 'nodes' and 'edges'.")}

  if(!is.character(v)){stop("Invalid input for v: Only accept a <character vector>")}

  if(length(exact) != 1 & !is.logical(exact)){stop("Invalid input for exact: Only accept <TRUE> OR <FALSE>")}

  if(length(output) != 1 || class(output) != "character"){ stop("Invalid input for output: Only accepts one of 'sf', 'df', or 'graph'.") }
  if (!output %in% c("sf", "df", "graph")) {stop("Invalid input for output: Only accepts one of 'sf', 'df', or 'graph'.")}


  # --- Main function implementation ---
  if(output == "df" | output == "graph"){
    # Drop the geometry for edges if the output is df or graph
    net$edges = net$edges %>% st_drop_geometry()}

  ## Subset edges when the indices <v> are available in either one of the FROM and TO nodes.
  if(exact == F){
    net$edges = net$edges[net$edges$FROM %in% v | net$edges$TO %in% v,]
    net$nodes = net$nodes[net$nodes$nID %in% net$edges$FROM | net$nodes$nID %in% net$edges$TO,]}
  ## Subset edges only when the indices <v> are available in both FROM and TO nodes
  else if(exact == T){
    net$edges = net$edges[net$edges$FROM %in% v & net$edges$TO %in% v,]
    net$nodes = net$nodes[net$nodes$nID %in% v,]}
  ## Return empty output if no edges are subset.
  if(nrow(net$edges) == 0){
    warning("Empty edge is returned. Please consider to use <exact = F> OR check if both nodes of edge is covered by the node indexs.")
  }

  ## Define outputs
  if(output == "graph"){
    # Return igraph
    return(graph_from_data_frame(net$edges, directed = T, vertices = net$nodes))}  # <---- OUTPUT: igraph
  else if(output == "df"){
    # Return a list of edges and nodes in data.frame format
    net$nodes = data.frame(nID = net$nodes$nID)
    return(net)}  # <---- OUTPUT: list
  else{
    # Return a list of edges and nodes in sf format
    return(net)}  # <---- OUTPUT: list
}

##_____________________________________________________________________________
#' Find the Nearest Network Node for Each 3D Point
#'
#' Identifies the nearest 3D network node for each input 3D point and attaches its index (\code{nID}).
#'
#' @param node_sf An \code{sf} object with 3D \code{POINT} geometry, typically from \code{\link[GISnetwork3D]{sf_to_net3d}}, containing a column \code{nID}.
#' @param pt An \code{sf} object with 3D \code{POINT} geometry, representing the input points for which to find the nearest network node.
#'
#' @details
#' For each 3D point in \code{pt}, the function identifies the closest 3D node in \code{node_sf} based on 3D Euclidean distance and appends the corresponding \code{nID} as a new column.
#'
#' @return An \code{sf} object with 3D \code{POINT} geometry, identical to \code{pt}, with an additional column \code{nID} specifying the index of the nearest network node for each point.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Create example network nodes
#' nodes = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 0)),
#'   st_point(c(1, 1, 1)),
#'   st_point(c(2, 2, 2))
#' ))
#' nodes$nID = 1:nrow(nodes)  # Add node IDs
#'
#' # Create example points to find nearest nodes
#' pts = st_sf(geometry = st_sfc(
#'   st_point(c(0.5, 0.5, 0.5)),
#'   st_point(c(0, 0, 1.5)),
#'   st_point(c(3, 3, 3))
#' ))
#'
#' # Find nearest nodes
#' nearest_pts = nearestNode3d(nodes, pts)
#' }
#'
#' @export
nearestNode3d = function(node_sf, pt){

  # --- Input validation ---
  if (!inherits(node_sf, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(node_sf))) {
    stop("Invalid input for node_sf: Only accept <sf POINT feature>.") }

  if (!inherits(pt, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(pt))) {
    stop("Invalid input for pt: Only accept <sf POINT feature>.") }

  # --- Main function implementation ---
  pt = cbind(pt, data.frame(nID = node_sf$nID[nearest3d(pt, node_sf)])) # Find nearest nodes for the input point
  return(pt)  # <---- OUTPUT: sf POINT feature
}


##_____________________________________________________________________________
#' Find the Nearest Network Node for Each 3D Point with Restriction
#'
#' Identifies the nearest 3D network node for each input 3D point, excluding nodes connected to restricted edges, and attaches its index (\code{nID}).
#'
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects with \code{character} \code{nID} column) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param pt An \code{sf} object with 3D \code{POINT} geometry, representing the input points for which to find the nearest network node.
#' @param restriction A \code{character} value specifying a field in \code{net$edges} with \code{integer} values (0 = not restricted, 1 = restricted).
#'
#' @details
#' For each 3D point in \code{pt}, the function identifies the closest 3D node in \code{net$nodes} based on 3D Euclidean distance, excluding nodes connected to edges where the \code{restriction} field is 1. The \code{restriction} parameter ensures that only nodes from non-restricted edges (value 0) are considered. The corresponding \code{nID} is appended as a new column to the output.
#'
#' @return An \code{sf} object with 3D \code{POINT} geometry, identical to \code{pt}, with an additional column \code{nID} specifying the index of the nearest non-restricted network node for each point.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' pts = st_sf(geometry = st_sfc(
#'   st_point(c(1,1,2)),
#'   st_point(c(6,6,3)),
#'   st_point(c(11, 11, 7))
#' ))
#'
#' # Create a line network
#' l <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 2,2,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(2,2,2, 5,5,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(5,5,2, 9,9,5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(9,9,5, 11,11,6), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(11,11,6, 13,13,6), ncol=3, byrow=TRUE))
#' ))
#' l$restriction = c(0,1,1,0,0) # Assign restriction
#'
#' # Create network
#' net <- sf_to_net3d.full(l, direction = rep(0, nrow(l)), CF = "Tobler") # Assumes bi-directional
#'
#' # Find nearest nodes
#' pts.noRestriction = nearestNode3d(net$nodes, pts) # Without restriction
#' pts.restriction = nearestNode3d.rst(net, pts, "restriction") # With restriction
#' }
#'
#' @export
nearestNode3d.rst = function(net, pt, restriction){

  # --- Input validation ---
  if (!inherits(net, "list") || !"nodes" %in% names(net) || !"edges" %in% names(net)) {
    stop("Invalid input for net: Only accepts <list> of network containing 'nodes' and 'edges'.")}

  if (!inherits(pt, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(pt))) {
    stop("Invalid input for pt: Only accept <sf POINT feature>.") }

  if (!is.character(restriction) || length(restriction) != 1) {
    stop("Invalid input for restriction: Must be a positive character.")}

  if (!restriction %in% names(net$edges)) {
    stop("Invalid input for restriction: Field name is not available in the net's edges.")}

  restriction = net$edges[["restriction"]]

  if (!all(sort(unique(restriction)) %in% c(0, 1))) {
    stop("Invalid input for restriction: The values from y's restriction field must only be 0 and 1.")
  }

  # --- Main function implementation ---
  # subset nodes
  nodes = net$edges[restriction == 0,] # Retain edges without restriction
  nodes = c(nodes$FROM, nodes$TO) %>% unique() # Get the nID from edges without restriction
  nodes = net$nodes[net$nodes$nID %in% nodes, ] # Get nodes in sf format
  pt = cbind(pt, data.frame(nID = nodes$nID[nearest3d(pt, nodes)])) # Find nearest nodes for the input point
  return(pt)  # <---- OUTPUT: sf POINT feature
}


