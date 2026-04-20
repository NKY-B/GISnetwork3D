##_____________________________________________________________________________
#' Add Z (Elevation) Values to 2D POINT or LINESTRING Geometries
#'
#' Adds Z (elevation) values to each point or vertex of an \code{sf} object with 2D \code{POINT} or 2D \code{LINESTRING} geometry, using values extracted from a raster layer. If no raster is provided, Z values are set to 0.
#'
#' @param x An \code{sf} object with 2D \code{POINT} or 2D \code{LINESTRING} geometry. The input features to which Z values will be added.
#' @param y \code{NULL} or a single-layer \code{RasterLayer} object. If provided, raster values at each point or vertex are used as Z values. If \code{NULL}, Z values are set to 0.
#'
#' @details
#' If \code{y} is a \code{RasterLayer}, Z values are extracted using \code{raster::extract()}. Missing values (\code{NA}) in the raster result in \code{NA} Z values in the output, which may cause issues in downstream processing. Ensure the raster covers the extent of \code{x} to avoid \code{NA} values. The output is a 3D \code{sf} object with Z coordinates added to each point or vertex.
#'
#' @return A 3D \code{sf} object of the same geometry type as \code{x}, with Z values added as the third coordinate.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(raster)
#' # Example raster
#' r = raster::raster(matrix(1:100, 10, 10), xmn=0, xmx=10, ymn=0, ymx=10)
#' # Example sf point
#' pt = st_sf(id=1, geometry=st_sfc(st_point(c(5,5))))
#' # Add Z from raster
#' pt_z = addZ(pt, r)
#' # Add Z=0
#' pt_z0 = addZ(pt)
#' }
#' @export

addZ = function(x, y = NULL){
  # --- Input validation ---
  if (!inherits(x, "sf") || (!"sfc_LINESTRING" %in% class(sf::st_geometry(x)) & !"sfc_POINT" %in% class(sf::st_geometry(x)))) {
    stop("x must be an <sf> object with POINT or LINESTRING geometry.") }

  if(!is.null(y)){
    if(!inherits(y, "RasterLayer")){stop("y must be a <RasterLayer>.")}
  }

  # --- Main function implementation ---
  ## Assign Z = 0 if y is NULL
  if(is.null(y)){
    # POINT: Convert sf object to data frame, add Z = 0, convert back to sf POINT
    if( "sfc_POINT" %in% class(sf::st_geometry(x)) ){
      pt = x %>% sfheaders::sf_to_df(fill = T) # Convert sf to data frame
      pt$z = 0 # Assign 0 to Z coordinate
      pt = sfheaders::sf_point(pt, x = "x", y = "y", z = "z", keep = T) # Convert data frame back to sf POINT
      st_crs(pt) = st_crs(x) # Assign original CRS
      return(pt)}  # <----  Output: sf POINT feature

    # LINESTRING: Convert sf object to data frame, add Z = 0, convert back to sf LINESTRING
    if( "sfc_LINESTRING" %in% class(sf::st_geometry(x)) ){
      l = x %>% sfheaders::sf_to_df(fill = T) # Convert sf to data frame
      l$z = 0 # Assign 0 to Z coordinate
      l = sfheaders::sf_linestring(l, x = "x", y = "y", z = "z", linestring_id = "linestring_id", keep = T) # Convert data frame back to sf LINESTRING
      st_crs(l) = st_crs(x) # Assign original CRS
      return(l)}  # <---- Output: sf LINESTRING feature
  }

  # --- Assign Z value from raster 'y' to 'x' if 'y' is not NULL ---
  if(!is.null(y)){
    # POINT: Extract Z values from raster at point locations
    if( "sfc_POINT" %in% class(sf::st_geometry(x)) ){
      pt = x
      pt$z = raster::extract(y, pt) # Extract Z values from raster
      pt = pt %>% sfheaders::sf_to_df(fill = T) # Convert sf to data frame
      pt = sfheaders::sf_point(pt, x = "x", y = "y", z = "z", keep = T) # Convert back to sf POINT with Z
      pt = pt %>% dplyr::select(-sfg_id, -point_id)
      st_crs(pt) = st_crs(x) # Assign CRS
      return(pt)}  # <---- Output: sf POINT feature

    # LINESTRING: Extract Z values for each vertex of the linestring
    if( "sfc_LINESTRING" %in% class(sf::st_geometry(x)) ){
      OUTPUTcrs = st_crs(x) # Store original CRS
      l = x %>% sfheaders::sf_to_df(fill = T)
      # Convert sf linestring to points for extraction
      l = sfheaders::sf_point(l, x = "x", y = "y", keep = T)
      l = l %>% dplyr::select(-sfg_id) # Remove geometry group identifier
      st_crs(l) = st_crs(x) # Assign CRS
      l$z = raster::extract(y, l) # Extract Z values from raster
      # Convert points with Z back to linestring
      l = l %>% sfheaders::sf_to_df(fill = T) %>% dplyr::select(-sfg_id, -point_id)
      l = sfheaders::sf_linestring(l, x = "x", y = "y", z = "z", linestring_id = "linestring_id", keep = T)
      l = l %>% dplyr::select(-linestring_id)
      st_crs(l) = st_crs(x) # Assign CRS
      return(l)}  # <----  Output: sf LINESTRING feature
  }
}

##_____________________________________________________________________________
#' Segmentise 3D Linestrings into Individual Line Segments
#'
#' Converts an \code{sf} object with 3D \code{LINESTRING} geometry into individual segments, where each segment consists of exactly two consecutive points from the original linestring.
#'
#' @param x An \code{sf} object with 3D \code{LINESTRING} geometry. The input linestrings to be segmented.
#' @param keep A \code{logical} value indicating whether to retain attributes:
#' \itemize{
#'   \item \code{TRUE}: Retains and joins all attributes from the original \code{x} to each segment.
#'   \item \code{FALSE}: Removes all attributes from the original \code{x}, resulting in segments with only geometry.
#' }
#'
#' @details
#' This function splits each 3D linestring feature in the input into a series of two-point line segments, each represented as a 3D \code{LINESTRING} geometry with Z coordinates.
#' If \code{keep = TRUE}, all original attributes are retained; otherwise, only geometry is retained.
#'
#' @return A 3D \code{sf} object with \code{LINESTRING} geometry, where each feature represents a two-point segment from the original linestrings.
#' Includes columns \code{L1} (unique ID of the original linestring) and \code{L2} (unique ID for the segmented linestring).
#' If \code{keep = TRUE}, all original attributes are included;
#' if \code{keep = FALSE}, only geometry and \code{L1}, \code{L2} columns are returned.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Example 3D linestring
#' y = sf::st_sf(geometry = sf::st_sfc(sf::st_linestring(matrix(
#'   c(0, 0, 0, 1, 2, 1, 2, 0, 2, 3, 2, 3, 4, 0, 4), ncol = 3, byrow = TRUE))))
#'
#' # Segmentise the linestring
#' segs = segmentise3d(y, keep = FALSE)
#'
#' # Extract segment coordinates and plot
#' segs_coords = sf::st_coordinates(segs)
#' seg_ids = unique(segs_coords[, "L1"])
#' seg_col = rainbow(length(seg_ids))
#'
#' # 3D plot
#' rgl::open3d()
#' for (i in seq_along(seg_ids)) {
#'   seg = segs_coords[segs_coords[, "L1"] == seg_ids[i], ]
#'   rgl::lines3d(seg[, 1], seg[, 2], seg[, 3], color = seg_col[i], lwd = 2)
#' }
#'
#' # Plot vertices
#' y_coords = sf::st_coordinates(y)
#' rgl::points3d(y_coords[, 1], y_coords[, 2], y_coords[, 3], color = "black", size = 8)
#'
#' # Plot original linestring (grey)
#' rgl::lines3d(y_coords[, 1], y_coords[, 2], y_coords[, 3], color = "grey", lwd = 6)
#'
#' # Add legend
#' rgl::legend3d("topright",
#'   legend = c(paste("Segment", seq_along(seg_ids)), "Vertices", "Linestring"),
#'   col = c(seg_col, "black", "grey"),
#'   lwd = c(rep(2, length(seg_ids)), NA, 6),
#'   pch = c(rep(NA, length(seg_ids)), 16, NA)
#' )
#' }
#' @export

segmentise3d = function(x, keep = FALSE){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with LINESTRING geometry.") }

  if (!is.logical(keep) || length(keep) != 1){
    stop("keep must be a single logical value (TRUE/FALSE).")}

  # --- Main function implementation ---
  # Extract points from line and then convert to data.frame
  result = x %>% st_coordinates() %>% as.data.frame() # Convert sf to data.frame
  result = result %>% segmentise_cpp() # Segmentise each linestring
  result = data.frame(X = result$X, Y = result$Y, Z = result$Z, L1 = result$L1, L2 = result$L2) # Create data.frame

  # Prepare the output geometry
  if(keep == T){
    # Join attribute to the output
    x$L1 = 1:nrow(x)
    result = dplyr::left_join(result, st_drop_geometry(x), by = "L1")
    # Convert df to sf
    result = sfheaders::sf_linestring(obj = result, x = "X", y = "Y", z = "Z", linestring_id = "L2", keep = T)
    st_crs(result) = st_crs(x) # Define crs
  }
  else if(keep == F){
    # Convert df to sf
    result = sfheaders::sf_linestring(obj = result, x = "X", y = "Y", z = "Z", linestring_id = "L2", keep = T)
    st_crs(result) = st_crs(x) } # Define crs

  return(result) # <---- OUTPUT: sf LINESTRING feature
}



##_____________________________________________________________________________
#' Extract Start or End Nodes from 3D Linestrings
#'
#' Extracts the start (first) or end (last) node of each linestring in an \code{sf} object with 3D \code{LINESTRING} geometry. Nodes can be returned as an \code{sf} object with 3D \code{POINT} geometry or as a \code{data.frame} of coordinates.
#'
#' @param x An \code{sf} object with 3D \code{LINESTRING} geometry. The input linestrings.
#' @param sf A \code{logical} value indicating the output format. Default: \code{TRUE}.
#' \itemize{
#'   \item \code{TRUE}: Returns nodes as an \code{sf} object with 3D \code{POINT} geometry.
#'   \item \code{FALSE}: Returns nodes as a \code{data.frame}.
#' }
#' @param position A \code{character} string specifying which nodes to extract. Must be one of \code{"first"} or \code{"last"}.
#' \itemize{
#'   \item \code{"first"}: Returns the start node of each linestring.
#'   \item \code{"last"}: Returns the end node of each linestring.
#' }
#' @details
#' This function extracts either the first or last node of each 3D linestring in the input.
#'
#' @return
#' If \code{sf = TRUE}, returns an \code{sf} object with 3D \code{POINT} geometry containing the extracted nodes.
#' If \code{sf = FALSE}, returns a \code{data.frame} with columns \code{X}, \code{Y}, \code{Z}, and \code{ID}, where \code{ID} is a \code{numeric} value identifying the source linestring.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Create an sf linestring object
#' l = st_linestring(matrix(c(0,0,0, 1,1,1, 2,2,2), ncol=3, byrow=TRUE))
#' sfl = st_sf(id = 1, geometry = st_sfc(l))
#' # Extract start nodes as sf POINTS:
#' start_pts = getNode3d(sfl, sf = TRUE, position = "first")
#' # Extract end nodes as data.frame:
#' end_pts = getNode3d(sfl, sf = FALSE, position = "last")
#' }
#' @export

getNode3d = function(x, sf = TRUE, position = "first"){

  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with LINESTRING geometry..") }

  if (!is.logical(sf) || length(sf) != 1){
    stop("sf must be a single logical value (TRUE/FALSE).")}

  if(length(position) != 1){stop("position must be a character of either `first` OR `last`")}
  if(!position %in% c("first", "last")){stop("position must be a character of either `first` OR `last`")}

  # --- Main function implementation ---
  outputCRS = st_crs(x) # Get the CRS of the input
  x = x %>% st_coordinates() %>% as.data.frame() # Convert from sf to data.frame

  # Extract nodes
  # Extract start nodes when position = "first"
  if(position == "first"){
    x = getFirstNode_cpp(x) %>% as.data.frame() }

  # Extract end nodes when position = "last"
  if(position == "last"){
    x = getLastNode_cpp(x) %>% as.data.frame() }

  names(x) = c("X","Y","Z")
  x$ID = 1:nrow(x) # Assign ID

  if(sf == T){
    x = sfheaders::sf_point(obj = x, x = "X", y = "Y", z = "Z", keep = T)
    st_crs(x) = outputCRS
    return(x)}  # <---- OUTPUT: sf POINT feature (sf = TRUE)

  if(sf == F) {
    return(x)}  # <---- OUTPUT: data.frame (sf = FALSE)
  }

##_____________________________________________________________________________
#' Connect Pairs of 3D Points with 3D Linestring Segments
#'
#' Creates a 3D linestring for each pair of 3D points from two \code{sf} objects with 3D \code{POINT} geometry, connecting the i-th point in \code{x} to the i-th point in \code{y}.
#' Both input objects must have the same number of geometries.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry. The starting points of the segments.
#' @param y An \code{sf} object with 3D \code{POINT} geometry. The ending points of the segments.
#'
#' @details
#' This function generates a 3D linestring for each corresponding pair of 3D points in \code{x} and \code{y}, preserving the Z coordinates. The order of the linestring features in the output matches the order of the input pairs.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, where each feature connects a pair of 3D points from \code{x} and \code{y}, in the order of the input pairs.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Create two sets of 3D points
#' x = st_sf(geometry = st_sfc(st_point(c(0, 0, 1)), st_point(c(1, 1, 2))))
#' y = st_sf(geometry = st_sfc(st_point(c(1, 0, 3)), st_point(c(2, 2, 4))))
#'
#' # Generate connecting 3D line segments
#' segs = connect3d.pair(x, y)
#'
#' # Extract coordinates for plotting
#' x_coords = st_coordinates(x)
#' y_coords = st_coordinates(y)
#' segs_coords = st_coordinates(segs)
#'
#' # Open a 3D plot
#' open3d()
#' # Plot start points (red)
#' points3d(x_coords[, 1], x_coords[, 2], x_coords[, 3], color = "red", size = 8)
#' # Plot end points (blue)
#' points3d(y_coords[, 1], y_coords[, 2], y_coords[, 3], color = "blue", size = 8)
#' # Plot line segments (black)
#' for (i in unique(segs_coords[, "L1"])) {
#'   seg = segs_coords[segs_coords[, "L1"] == i, ]
#'   lines3d(seg[, 1], seg[, 2], seg[, 3], color = "black", lwd = 3)
#' }
#'
#' # Add legend
#' legend3d("topright", legend = c("x (start)", "y (end)", "segment"),
#'          pch = c(16, 16, NA), col = c("red", "blue", "black"), lwd = c(NA, NA, 3))
#' }
#' @export

connect3d.pair = function(x, y){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with POINT geometry.") }
  if (!inherits(y, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(y))) {
    stop("y must be an <sf> object with POINT geometry.") }

  if (nrow(x) != nrow(y)) {
    stop("x and y must have the same number of features.") }

  # --- Main function implementation ---
  OutputCRS = st_crs(x) # Define the output coordinate reference system

  # Connect each pair of points
  result = connect3d_pair_cpp(x %>% st_coordinates() %>% as.data.frame(),
                              y %>% st_coordinates() %>% as.data.frame())

  # Convert data frame to sf LINESTRING
  result = result %>% sfheaders::sf_linestring(x = "X", y = "Y", z = "Z", linestring_id = "id")
  st_crs(result) = OutputCRS # Assign the original CRS
  return(result)  # <---- OUTPUT: sf LINESTRING feature
}

##_____________________________________________________________________________
#' Connect All Pairs of 3D Points Between Two Sets (Many-to-Many)
#'
#' Creates 3D linestrings connecting every 3D point in \code{x} to every 3D point in \code{y}. The output includes attributes \code{id1} and \code{id2} indicating the \code{numeric} indices of the start and end points from \code{x} and \code{y}, respectively.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry. The start points.
#' @param y An \code{sf} object with 3D \code{POINT} geometry. The end points.
#'
#' @details
#' This function generates a 3D linstring for each possible pair of 3D points between \code{x} and \code{y}. Each linestring feature includes attributes \code{id1} and \code{id2}, which are \code{numeric} values identifying the indices of the start point (from \code{x}) and end point (from \code{y}).
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, where each feature represents a segment connecting one point in \code{x} to one point in \code{y}, with \code{numeric} attributes \code{id1} and \code{id2} indicating the source point indices.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#' # Create two sets of 3D points
#' x = st_sf(geometry = st_sfc(st_point(c(0, 0, 0)), st_point(c(1, 0, 1))))
#' y = st_sf(geometry = st_sfc(st_point(c(0, 1, 2)), st_point(c(1, 1, 3))))
#'
#' # Generate all connecting 3D line segments
#' segs = connect3d.ManyToMany(x, y)
#'
#' # Extract coordinates for plotting
#' x_coords = st_coordinates(x)
#' y_coords = st_coordinates(y)
#' segs_coords = st_coordinates(segs)
#'
#' open3d()
#' # Plot start points (red)
#' points3d(x_coords[, 1], x_coords[, 2], x_coords[, 3], color = "red", size = 10)
#' # Plot end points (blue)
#' points3d(y_coords[, 1], y_coords[, 2], y_coords[, 3], color = "blue", size = 10)
#' # Plot line segments (black)
#' for (i in unique(segs_coords[, "L1"])) {
#'   seg = segs_coords[segs_coords[, "L1"] == i, ]
#'   lines3d(seg[, 1], seg[, 2], seg[, 3], color = "black", lwd = 3)
#' }
#' legend3d("topright", legend = c("x (start)", "y (end)", "segment"),
#'          pch = c(16, 16, NA), col = c("red", "blue", "black"), lwd = c(NA, NA, 3))
#' }
#' @export

connect3d.ManyToMany = function(x, y){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with POINT geometry.") }
  if (!inherits(y, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(y))) {
    stop("y must be an <sf> object with POINT geometry.") }

  # --- Main function implementation ---
  OutputCRS = st_crs(x) # Define the output coordinate reference system
  # Generate all pairwise connections
  result = connect3d_ManyToMany_cpp(
    x %>% st_coordinates() %>% as.data.frame(),
    y %>% st_coordinates() %>% as.data.frame()
  )
  # Assign identifiers for each segment
  result$L1 = rep(1:(nrow(result)/2), each = 2)
  # Convert data frame to sf LINESTRING
  result = result %>% sfheaders::sf_linestring(
    x = "X", y = "Y", z = "Z",
    linestring_id = "L1",
    keep = TRUE
  )
  # Set the coordinate reference system
  st_crs(result) = OutputCRS
  return(result)  # <---- OUTPUT: sf LINESTRING feature
}


##_____________________________________________________________________________
#' Snap 3D Points to 3D Linestrings by Node, Vertex, or Projection
#'
#' Snaps each point in an \code{sf} object with 3D \code{POINT} geometry to its nearest location on an \code{sf} object with 3D \code{LINESTRING} geometry, according to the specified snapping rule.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry. The points to snap to \code{y}.
#' @param y An \code{sf} object with 3D \code{LINESTRING} geometry. The reference linestrings for snapping.
#' @param tolerance A \code{numeric} value specifying the maximum 3D distance within which snapping is performed.
#' @param rule A \code{character} string specifying the snapping rule. Must be one of \code{"node"}, \code{"vertex"}, or \code{"project"}.
#' \itemize{
#'   \item \code{"node"}: Snaps to the nearest node (start or end) of a linestring.
#'   \item \code{"vertex"}: Snaps to the nearest vertex along the linestring.
#'   \item \code{"project"}: Snaps to the 3D orthogonal projection onto the linestring segment.
#' }
#'
#' @details
#' This function snaps each 3D point in \code{x} to the nearest location on the 3D linestrings in \code{y} within the specified \code{tolerance}, based on the chosen \code{rule}. Points beyond the \code{tolerance} distance remain unchanged.
#'
#' @return An \code{sf} object with 3D \code{POINT} geometry, where points within \code{tolerance} are replaced by their snapped positions, and others remain unchanged.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # 3D linestring (wiggle)
#' y = st_sf(geometry = st_sfc(st_linestring(matrix(
#'   c(0,0,0, 2,1,1, 4,0,2, 6,1,3, 8,0,4), ncol=3, byrow=TRUE))))
#'
#' # Create sample points
#' x = st_sf(geometry = st_sfc(
#'   st_point(c(2,1.2,1)),
#'   st_point(c(4,0,2.5)),
#'   st_point(c(0,2,0)),
#'   st_point(c(8,0,4))
#' ))
#'
#'# Set tolerance to 1 unit
#' tolerance = 1.0
#'
#' # Snap points using each rule
#' snap_node   = snap3d.p2l(x, y, tolerance, rule="node")
#' snap_vertex = snap3d.p2l(x, y, tolerance, rule="vertex")
#' snap_proj   = snap3d.p2l(x, y, tolerance, rule="project")
#'
#' # Plot all results for comparison
#' open3d()
#' y_coords = st_coordinates(y)
#' x_coords = st_coordinates(x)
#' node_coords = st_coordinates(snap_node)
#' vertex_coords = st_coordinates(snap_vertex)
#' proj_coords = st_coordinates(snap_proj)
#'
#' # Reference linestring
#' lines3d(y_coords[,1], y_coords[,2], y_coords[,3], color="black", lwd=4)
#'
#' # Original points (red)
#' points3d(x_coords[,1], x_coords[,2], x_coords[,3], color="red", size=10)
#'
#' # Snapped points (different colors for each rule)
#' points3d(node_coords[,1], node_coords[,2], node_coords[,3], color="blue", size=12)
#' points3d(vertex_coords[,1], vertex_coords[,2], vertex_coords[,3], color="green", size=12)
#' points3d(proj_coords[,1], proj_coords[,2], proj_coords[,3], color="orange", size=12)
#'
#' # Add legend
#' legend3d("topright",
#'   legend = c("Reference line", "Original points", "Snapped (node)", "Snapped (vertex)", "Snapped (project)"),
#'   col = c("black", "red", "blue", "green", "orange"),
#'   pch = c(NA, 16, 16, 16, 16), lwd = c(4, NA, NA, NA, NA)
#' )
#' }
#' @export

snap3d.p2l = function(x, y, tolerance, rule = "node"){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with POINT geometry.") }

  if (!inherits(y, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(y))) {
    stop("y must be an <sf> object with LINESTRING geometry.") }

  if (!rule %in% c("node", "vertex", "project")) {
    stop("rule must be a character of either 'node', 'vertex', or 'project'.")}

  # --- Main function implementation ---
  if(rule == "node"){
    # Define the output coordinate
    OutputCRS = st_crs(x)

    # Snap the point to line if it is within tolerance distance
    result = snap3d_p2p_cpp( x %>% st_coordinates() %>% as.data.frame(),
                             rbind(getNode3d(y, sf = F, position = "first"),
                                   getNode3d(y, sf = F, position = "last")) %>% unique() %>% as.data.frame(),
                             tolerance)
    result = result %>% sfheaders::sf_point(x="X", y="Y", z="Z")
    st_crs(result) = OutputCRS
    st_geometry(x) = st_geometry(result)
    return(x)}  # <---- OUTPUT: sf POINT feature
  else if(rule == "vertex"){
    # Define the output coordinate
    OutputCRS = st_crs(x)

    # Snap the point to line if it is within tolerance distance
    result = snap3d_p2p_cpp( x %>% st_coordinates() %>% as.data.frame(),
                             y %>% st_coordinates() %>% unique() %>% as.data.frame(),
                             tolerance)
    result = result %>% sfheaders::sf_point(x="X", y="Y", z="Z")
    st_crs(result) = OutputCRS
    st_geometry(x) = st_geometry(result)
    return(x)}  # <---- OUTPUT: sf POINT feature
  else if(rule == "project"){
    # Define the output coordinate
    OutputCRS = st_crs(x)

    # Segmentise line
    y = y %>% st_coordinates() %>% segmentise_cpp()
    fromIndex = seq(1, length(y$L1), 2)
    toIndex = seq(2, length(y$L1), 2)
    y = data.frame(from_x = y$X[fromIndex], from_y = y$Y[fromIndex], from_z = y$Z[fromIndex],
                   to_x = y$X[toIndex], to_y = y$Y[toIndex], to_z = y$Z[toIndex],
                   L1 = y$L1[fromIndex], L2 = y$L2[fromIndex])

    # Snap the point to line if it is within tolerance distance
    result = snap3d_p2l_cpp(x %>% st_coordinates() %>% as.data.frame(), y, tolerance)
    result = result %>% sfheaders::sf_point(x="X", y="Y", z="Z")
    st_crs(result) = OutputCRS
    st_geometry(x) = st_geometry(result)
    return(x)  # <---- OUTPUT: sf POINT feature
  }
}


##_____________________________________________________________________________
#' Snap 3D Linestring Vertices to Reference 3D Linestring
#'
#' Snaps each vertex of an \code{sf} object with 3D \code{LINESTRING} geometry to the nearest location on a reference \code{sf} object with 3D \code{LINESTRING} geometry, according to the specified snapping rule.
#'
#' @param x An \code{sf} object with 3D \code{LINESTRING} geometry. The linestring whose vertices will be snapped.
#' @param y An \code{sf} object with 3D \code{LINESTRING} geometry. The reference linestring for snapping.
#' @param tolerance A \code{numeric} value specifying the maximum 3D distance within which snapping is performed.
#' @param rule A \code{character} string specifying the snapping rule. Must be one of \code{"node"}, \code{"vertex"}, or \code{"project"}.
#' \itemize{
#'   \item \code{"node"}: Snaps each vertex of \code{x} to the nearest node (start or end) of \code{y}.
#'   \item \code{"vertex"}: Snaps each vertex of \code{x} to the nearest vertex of \code{y}.
#'   \item \code{"project"}: Snaps each vertex of \code{x} to its 3D orthogonal projection onto \code{y}.
#' }
#'
#' @details
#' This function snaps each vertex of the 3D linestring in \code{x} to the nearest location on the 3D linestring in \code{y} within the specified \code{tolerance}, based on the chosen \code{rule}. Vertices beyond the \code{tolerance} distance remain unchanged.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry, where vertices within \code{tolerance} are snapped to the reference linestring according to the selected \code{rule}.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Reference 3D line (wiggle)
#' y = st_sf(geometry = st_sfc(st_linestring(matrix(
#'   c(0,0,0, 2,1,1, 4,0,2, 6,1,3, 8,0,4), ncol=3, byrow=TRUE))))
#'
#' # Line to be snapped (offset above)
#' x = st_sf(geometry = st_sfc(st_linestring(matrix(
#'   c(0,0,1, 2,1,2, 4,0,3, 6,1,4, 8,0,5), ncol=3, byrow=TRUE))))
#'
#' tolerance = 1.5
#'
#' # Snap with each rule
#' snap_node   = snap3d.l2l(x, y, tolerance, rule = "node")
#' snap_vertex = snap3d.l2l(x, y, tolerance, rule = "vertex")
#' snap_proj   = snap3d.l2l(x, y, tolerance, rule = "project")
#'
#' # Extract coordinates for plotting
#' y_coords    = st_coordinates(y)
#' x_coords    = st_coordinates(x)
#' node_coords = st_coordinates(snap_node)
#' vertex_coords = st_coordinates(snap_vertex)
#' proj_coords   = st_coordinates(snap_proj)
#'
#' # 3D plot
#' open3d()
#' # Reference line (black)
#' lines3d(y_coords[,1], y_coords[,2], y_coords[,3], color="black", lwd=4)
#' # Original line (red, dashed)
#' lines3d(x_coords[,1], x_coords[,2], x_coords[,3], color="red", lwd=3, lty=2)
#' # Snapped lines (distinct colors)
#' lines3d(node_coords[,1], node_coords[,2], node_coords[,3], color="blue", lwd=4)
#' lines3d(vertex_coords[,1], vertex_coords[,2], vertex_coords[,3], color="green", lwd=10)
#' lines3d(proj_coords[,1], proj_coords[,2], proj_coords[,3], color="orange", lwd=4)
#'
#' # Add legend
#' legend3d("topright",
#'   legend = c("Reference line", "Original line", "Snapped (node)", "Snapped (vertex)", "Snapped (project)"),
#'   col = c("black", "red", "blue", "green", "orange"),
#'   lwd = c(4, 3, 4, 4, 4),
#'   lty = c(1, 2, 1, 1, 1)
#' )
#' }
#' @export

snap3d.l2l = function(x, y, tolerance, rule = "node"){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with LINESTRING geometry.") }

  if (!inherits(y, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(y))) {
    stop("y must be an <sf> object with LINESTRING geometry.") }

  if (!is.numeric(tolerance) || length(tolerance) != 1 || tolerance < 0) {
    stop("tolerance must be a non-negative numeric value.")}

  if(is.character(rule) == F | length(rule) > 1){
    stop("rule must be a character of either node', 'vertex', or 'project'.")}

  if (!rule %in% c("node", "vertex", "project")) {
    stop("rule must be a character of either node', 'vertex', or 'project'.")}

  # --- Main function implementation ---
  if(rule == "node"){
    # Define the output coordinate
    OutputCRS = st_crs(x)

    # Convert x to data.frame
    result = x %>% st_coordinates() %>% as.data.frame()
    L1 = result$L1

    # Snap the vertices of x to y if it is within tolerance distance
    result = snap3d_p2p_cpp( result,
                             rbind(getNode3d(y, sf = F, position = "first"),
                                   getNode3d(y, sf = F, position = "last")) %>% unique() %>% as.data.frame(),
                             tolerance) %>% cbind(L1)
    result = result %>% sfheaders::sf_linestring(x="X", y="Y", z="Z", linestring_id = "L1")
    st_crs(result) = OutputCRS
    st_geometry(x) = st_geometry(result)
    return(x)}  # <---- OUTPUT: sf LINESTRING feature
  else if(rule == "vertex"){
    # Define the output coordinate
    OutputCRS = st_crs(x)

    # Convert x to data.frame
    result = x %>% st_coordinates() %>% as.data.frame()
    L1 = result$L1

    # Snap the point to line if it is within tolerance distance
    result = snap3d_p2p_cpp( result,
                             y %>% st_coordinates() %>% unique() %>% as.data.frame(),
                             tolerance) %>% cbind(L1)
    result = result %>% sfheaders::sf_linestring(x="X", y="Y", z="Z", linestring_id = "L1")
    st_crs(result) = OutputCRS
    st_geometry(x) = st_geometry(result)
    return(x)}  # <---- OUTPUT: sf LINESTRING feature
  else if(rule == "project"){
    # Define the output coordinate
    OutputCRS = st_crs(x)

    # Convert x to data.frame
    result = x %>% st_coordinates() %>% as.data.frame()
    L1 = result$L1

    # Segmentise line
    y = y %>% st_coordinates() %>% segmentise_cpp()
    fromIndex = seq(1, length(y$L1), 2)
    toIndex = seq(2, length(y$L1), 2)
    y = data.frame(from_x = y$X[fromIndex], from_y = y$Y[fromIndex], from_z = y$Z[fromIndex],
                   to_x = y$X[toIndex], to_y = y$Y[toIndex], to_z = y$Z[toIndex],
                   L1 = y$L1[fromIndex], L2 = y$L2[fromIndex])

    # Snap the point to line if it is within tolerance distance
    result = snap3d_p2l_cpp(result, y, tolerance) %>% dplyr::select(-L1, -L2) %>% cbind(L1)
    result = result %>% sfheaders::sf_linestring(x="X", y="Y", z="Z", linestring_id = "L1")
    st_crs(result) = OutputCRS
    st_geometry(x) = st_geometry(result)
    return(result)  # <---- OUTPUT: sf LINESTRING feature
  }
}
