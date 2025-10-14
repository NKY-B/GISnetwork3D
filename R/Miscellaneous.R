# -----> This function subdivides a bounding box of the input sf feature into smaller boxes
# --> x = <sf POINT> feature: A spatial feature object of class `sf`. It must be a POINT object.
# --> nrows = <Integer>: Number of rows to divide the bounding box into.
# --> ncols = <Integer>: Number of columns to divide the bounding box into.
# --> cover = <Logical>: If `TRUE`, ensures that the tiles "cover" points by filtering out tiles with no points.

tile.create = function(x, nrows = 4, ncols = 4, cover = T){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("Invalid input for x: Only accept ONE <sf POINT> object.") }

  if (length(nrows) != 1 || !is.numeric(nrows) || nrows <= 0 || nrows %% 1 != 0) {
    stop("Invalid input for nrows: Only accept ONE positive <Integer>.")}

  if (length(ncols) != 1 || !is.numeric(ncols) || ncols <= 0 || ncols %% 1 != 0) {
    stop("Invalid input for ncols: Only accept ONE positive <Integer>.")}

  if (!is.logical(cover) || length(cover) != 1) {
    stop("Invalid input for cover: Only accept ONE <Logical> (TRUE or FALSE).")}

  # --- Main function implementation ---
  # Get the bounding box of the input x
  bbox = sf::st_bbox(x)

  # Calculate the height and width for each subdivided box
  height = (bbox["ymax"] - bbox["ymin"]) / nrows
  width = (bbox["xmax"] - bbox["xmin"]) / ncols
  outputCRS = sf::st_crs(x) # Get the coordinate reference system of the input
  tiles = list()  # Initialize an empty list to store subdivided tiles
  ibox = 0 # Initialize the box counter

  # Loop through each column and row to calculate the xmin, xmax, ymin, and ymax values
  for (i in 0:(ncols - 1)) {
    for (j in 0:(nrows - 1)) {
      x_min = bbox["xmin"] + i * width
      x_max = bbox["xmin"] + (i + 1) * width
      y_min = bbox["ymin"] + j * height
      y_max = bbox["ymin"] + (j + 1) * height
      ibox = ibox + 1
      BBOX = sf::st_bbox(c(xmin = as.numeric(x_min), xmax = as.numeric(x_max), ymax = as.numeric(y_max), ymin = as.numeric(y_min)),
                         crs = outputCRS)
      # Convert the bounding box into an sf POLYGON object
      BBOX = BBOX %>% sf::st_as_sfc() %>% sf::st_as_sf()
      tiles[[ibox]] = BBOX # Store the tile in the list
    }
  }
  # Combine the list of tiles into a single sf POLYGON object
  tiles = tiles %>% purrr::reduce(rbind)
  tiles$gridID = 1:nrow(tiles) # Assign unique grid IDs to each tile

  # --- Handle additional attributes based on `cover` argument ---
  # If the input `x` has a "nID" column
  if("nID" %in% names(x)){
    if(isTRUE(cover)){
      pt = tileNode_cpp(tiles %>% sf::st_coordinates(), x %>% sf::st_coordinates(), x$nID) # Intersect point by tile
      pt = pt[pt$id != "-999",] # Remove rows without points identified
      pt = pt[duplicated(pt$id) == F,] # Remove duplicated points
      tiles = tiles[unique(pt$gridID),] # Subset tiles that intersect with points
      pt$gridID = match(pt$gridID, tiles$gridID) # Find the gridID for each point
      tiles$gridID = seq_len(nrow(tiles))
      pt = split(pt[c("id", "nID")], pt$gridID) # Convert data.frame to point
      result = list(tile = tiles, pt = pt) }
    if(isFALSE(cover)){
      pt = tileNode_cpp(tiles %>% sf::st_coordinates(), x %>% sf::st_coordinates(), x$nID) # Intersect point by tile
      pt = split(pt[c("id", "nID")], pt$gridID)
      result = list(tile = tiles, pt = pt)} }

  # If the input `x` does NOT have an "nID" column
  if("nID" %in% names(x) == F){
    if(isTRUE(cover)){
      pt = tilePoint_cpp(tiles %>% sf::st_coordinates(), x %>% sf::st_coordinates()) # Intersect point by tile
      pt = pt[pt$id != "-999",] # Remove rows without points identified
      pt = pt[duplicated(pt$id) == F,] # Remove duplicated points
      tiles = tiles[unique(pt$gridID),] # Subset tiles that intersect with points
      pt$gridID = match(pt$gridID, tiles$gridID) # Find the gridID for each point
      tiles$gridID = seq_len(nrow(tiles))
      pt = split(pt$id, pt$gridID) # Convert data.frame to point
      result = list(tile = tiles, pt = pt) }
    if(isFALSE(cover)){
      pt = tilePoint_cpp(tiles %>% sf::st_coordinates(), x %>% sf::st_coordinates()) # Intersect point by tile
      pt = split(pt$id, pt$gridID) # Convert data.frame to point
      result = list(tile = tiles, pt = pt)} }

  return(result)  # <---- OUTPUT
}





# -----> This function subdivides sf points by tiles
# --> tiles = <sf POLYGON> feature: A spatial feature object of class `sf`. It must be a POLYGON object.
# --> x = <sf POINT> feature: A spatial feature object of class `sf`. It must be a POINT object.
# --> buf = <Numeric> or <NULL>: A buffer distance for expanding tiles (default is NULL, meaning no buffer).


tile.point = function(tiles, x, buf = NULL){
  # --- Input validation ---
  if (!inherits(tiles, "sf") || !"sfc_POLYGON" %in% class(sf::st_geometry(tiles))) {
    stop("Invalid input for tiles: Only accept an <sf POLYGON> object.") }

  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("Invalid input for x: Only accept an <sf POINT> object.") }

  if (!is.null(buf) && (!is.numeric(buf) || length(buf) != 1 || buf < 0)) {
    stop("Invalid input for buf: Only accept a <Numeric> value or <NULL>.") }


  # --- Main function implementation ---
  if(!"nID" %in% names(x)){
    if(is.null(buf)){
      # If no buffer is specified, intersect points by tiles to identify the GridID for each point `x`
      result = tilePoint_cpp(sf::st_coordinates(tiles), sf::st_coordinates(x)) }
    else {
      # If buffer is specified, expand each tile using the specified buffer distance
      expandTile = sf::st_buffer(x = tiles, dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2)

      for(i in 1:nrow(expandTile)){
        st_geometry(expandTile)[i] = GISnetwork3D::tile.expand.untill.cover(expandTile[i,], pt = x, buf = buf, coef = 1.2) }

      result = tilePoint_cpp(sf::st_coordinates(expandTile), sf::st_coordinates(x)) } # Identify tile GridID for each point after expanding tiles
  }

  if("nID" %in% names(x)){
    if(is.null(buf)){
      # If no buffer is specified, intersect points by tiles to identify the GridID for each point `x`
      result = tileNode_cpp(sf::st_coordinates(tiles), sf::st_coordinates(x), x$nID) }
    else {
      # If buffer is specified, expand each tile using the specified buffer distance
      expandTile = sf::st_buffer(x = tiles, dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2)

      for(i in 1:nrow(expandTile)){
        st_geometry(expandTile)[i] = GISnetwork3D::tile.expand.untill.cover(expandTile[i,], pt = x, buf = buf, coef = 1.2) }

      result = tileNode_cpp(sf::st_coordinates(expandTile), sf::st_coordinates(x), x$nID) } # Identify tile GridID for each point after expanding tiles
  }

  result = split(result %>% dplyr::select(-gridID), result$gridID)

  return(result)
}


# -----> This function tiles the origin and destination points
# --> o = <sf POINT feature>: A spatial feature object of class `sf` to indicate the origin points.
# --> d = <sf POINT feature>: A spatial feature object of class `sf` to indicate the destination points.
# --> buf = <numeric>: A buffer distance to expand the tiles for subtracting destination points.
# --> nrows = <Integer>: Number of rows to divide the bounding box of the origin points.
# --> ncols = <Integer>: Number of columns to divide the bounding box of the origin points.
# --> This function first creates a bounding box for the origin points <o>. Then, it creates non-overlapping tiles inside this box based on the specified number of rows <nrows> and columns <ncols>. Next, each tile is expanded by the distance defined by <buf>. These expanded tiles are used to subtract the destination points <d>. The function returns a list of stratified origin and destination points based on the tiles.

tile.OD = function(o, d, buf, nrows = 4, ncols = 4){
  # --- Input validation ---
  if (!inherits(o, "sf") || !"sfc_POINT" %in% class(st_geometry(o))) {
    stop("Invalid input for o: Only accept <sf POINT feature>.") }

  if (!inherits(d, "sf") || !"sfc_POINT" %in% class(st_geometry(d))) {
    stop("Invalid input for d: Only accept <sf POINT feature>.") }

  if (!is.null(buf) && (!is.numeric(buf) || length(buf) != 1 || buf < 0)) {
    stop("Invalid input for buf: Only accept a positive <Numeric> value or <NULL>.") }

  if (length(nrows) != 1 || !is.numeric(nrows) || nrows <= 0 || nrows %% 1 != 0) {
    stop("Invalid input for nrows: Only accept ONE positive <Integer>.") }

  if (length(ncols) != 1 || !is.numeric(ncols) || ncols <= 0 || ncols %% 1 != 0) {
    stop("Invalid input for ncols: Only accept ONE positive <Integer>.") }

  # --- Main function implementation ---

  # Create tiles from the bounding box of the origin points `o`
  tile = GISnetwork3D::tile.create(o, nrows = nrows, ncols = ncols)

  # Subtract destination points for each tile by expanding the tiles using the buffer distance `buf`
  destination = GISnetwork3D::tile.point(tiles = tile$tile, x = d, buf = buf)

  # Return a list containing stratified origin and destination points
  return( list(o = tile$pt, d = destination) )  # <---- OUTPUT
}



# This is an inner function of tile.point
# This function checks if a tile contains at least one point. If not, it iteratively expands the tile's geometry by applying a buffer until at least one point is encompassed. The function returns the geometry of the original or expanded tile as an sf POLYGON feature.
# Arguments:
#   tile: An sf POLYGON object representing the initial tile geometry (single feature).
#   pt: An sf POINT object.
#   buf: Numeric; initial buffer distance to expand the tile (in units of the CRS).
#   coef: Numeric; coefficient to multiply the buffer distance in each iteration (default: 1.2).
# Details:
#   1. It checks if the tile contains at least one point.
#   2. If the tile already contains point(s), its original geometry is returned.
#   3. If not, the tile is iteratively expanded by multiplying the buffer distance by coef and applying a buffer with square end caps and mitered joins until it includes at least one point.
#
# Returns:
#   An sf POLYGON object representing the geometry of the original or expanded tile.

tile.expand.untill.cover = function(tile, pt, buf, coef){

  pt.subset = tilePoint_cpp(sf::st_coordinates(tile), sf::st_coordinates(pt)) # Subset the pt by the tile
  pt.subset = pt.subset[pt.subset$id != -999,]
  pt.subset = nrow(pt.subset) # Count the number of point from the pt.subset

  if(pt.subset >= 1){
    # If the tile already contains point(s), return its geometry
    return(sf::st_geometry(tile))}  # <---- OUTPUT: sf POLYGON feature

  # If the tile does not contain any point, expand it iteratively
  if(!pt.subset >= 1){
    # Loop until one point covered
    while (TRUE) {
      buf = buf * coef # Increase the buffer distance by the coefficient

      # Expand the tile by applying a buffer with specified parameters
      # endCapStyle = "SQUARE" ensures square ends for linear features
      # joinStyle = "MITRE" ensures sharp corners
      # nQuadSegs = 5 controls the number of segments for curves
      # mitreLimit = 2 limits the extension of mitre joins
      tile = sf::st_buffer(x = tile, dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2) # Expand the tile

      # Subset points within the expanded tile's geometry
      pt.subset = tilePoint_cpp(sf::st_coordinates(tile), sf::st_coordinates(pt))

      # Count the number of point in the new subset
      pt.subset = pt.subset[pt.subset$id != -999,]
      pt.subset = nrow(pt.subset) # Count the number of class from the pt.subset

      # If the tile contains point, return its geometry and exit the loop
      if (pt.subset >= 1) {
        break }
    }
    # Return the geometry of the expanded tile
    return(sf::st_geometry(tile))  # <---- OUTPUT: sf POLYGON feature
  }

}


# This is an inner function of access3d.multitasks
# tile.multitasks.OD divides the bounding box of origin points (o) into a grid of tiles based on specified rows and columns, identifies origin points within each tile, and finds destination points (d) within a buffer distance (buf) from each tile. If a tile does not contain at least one destination point from each unique class, it expands the tile's geometry by multiplying the buffer distance by coef until all classes are included. The function returns a list containing the origin and destination points associated with each tile.
#
# Arguments:
#   o: An sf POINT object representing origin locations.
#   d: An sf POINT object representing destination locations, with a required "class" field indicating point categories.
#   buf: Numeric; initial buffer distance to expand each tile (in units of the CRS).
#   nrows: Numeric; number of rows to divide the bounding box (default: 4).
#   ncols: Numeric; number of columns to divide the bounding box (default: 4).
#   coef: Numeric; coefficient to multiply the buffer distance in each iteration (default: 1.2).
#
# Details:
#   1. The function uses GISnetwork3D::tile.create to divide the bounding box of origin points (o) into a grid of tiles with specified rows (nrows) and columns (ncols).
#   2. For each tile, it extracts the origin points (o) that fall within the tile's geometry.
#   3. It calls tile.multitasks.point to identify destination points (d) within the buffer distance of each tile, expanding the tile if necessary to include at least one point from each unique class in d$class.
#   4. The output is a list with two components: 'o' (a list of sf POINT objects for origins in each tile) and 'd' (a list of sf POINT objects for destinations in each tile).
#
# Returns:
#   A list with two elements:
#     - o: A list of sf POINT objects, each containing the origin points within the corresponding tile.
#     - d: A list of sf POINT objects, each containing the destination points within the buffered (and possibly expanded) tile.

tile.multitasks.OD = function(o, d, buf, nrows = 4, ncols = 4, coef = 1.2){
  # Create a tile that divide the bounding box of o into predeficed nuimber of rows and columns
  tiles = GISnetwork3D::tile.create(o, nrows = nrows, ncols = ncols, cover = T) # A list with: 1. tile in sf POLYGON and 2. a list of the data.frame indicating the covered point for each tile

  o = purrr::map(tiles$pt, ~ o[.x$id,]) # Extract the points in sf format for each tile
  tiles = tiles$tile # Extract the tile in sf POLYGON geometry
  d = GISnetwork3D::tile.multitasks.point(tiles = tiles, x = d, buf = buf, coef = coef) # Find the destination for each tile in list.
  return(list(o = o, d = d$points, otiles = tiles, dtiles = d$tiles))  # <---- OUTPUT: List

}

# This is an inner function of access3d.multitasks
# This function identifies intercepting points within each tile, ensuring tiles encompass points from all specified classes by iteratively expanding their geometry.
# Arguments:
#   tiles: An sf POLYGON object representing the initial tile geometries.
#   x: An sf POINT object representing points, with a required "class" field indicating point categories.
#   buf: Numeric; initial buffer distance to expand each tile (in units of the CRS).
#   coef: Numeric; coefficient to multiply the buffer distance in each iteration (default: 1.2).
# Details:
#   1. The function buffers the input tiles using the specified buffer distance (buf) with square end caps and mitered joins for precise geometry expansion.
#   2. For each tile, it checks if the tile contains at least one point from each unique class in x. If not, it iteratively expands the tile's geometry by multiplying the buffer distance by coef until all classes are represented or expansion criteria are met.
#   3. The function returns a list of sf POINT objects, where each element corresponds to a tile and contains the points (with their attributes) that intersect with that tile.
#
# Returns:
#   A list of sf POINT objects, where each element contains the points (with attributes) covered by the corresponding tile.

tile.multitasks.point = function(tiles, x, buf = NULL, coef = 1.2){

  # Buffering the tiles
  tiles = sf::st_buffer(x = tiles, dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2)

  # For each tile, check if it contain at least one point from each class. If not, it expands the tile's geometry by applying a buffer (multiply the buf by coef) until it encompasses points from all unique classes.
  for(i in 1:nrow(tiles)){
    st_geometry(tiles)[i] = GISnetwork3D::tile.multitasks.expand.untill.cover(tiles[i,], pt = x, buf = buf, coef = coef) }

  # Return a list of sf POINT objects covered by each tile
  return(list(points = purrr::map(GISnetwork3D::tile.point(tiles, x), ~ sf::st_drop_geometry(x[.x$id,])),
              tiles = tiles))  # <---- OUTPUT: List

}


# This is an inner function of access3d.multitasks
# This function checks if a tile contains points from all unique classes in the provided point dataset. If not, it iteratively expands the tile's geometry by applying a buffer until points from all unique classes are encompassed. The function returns the geometry of the original or expanded tile as an sf POLYGON feature.
# Arguments:
#   tile: An sf POLYGON object representing the initial tile geometry (single feature).
#   pt: An sf POINT object with a required 'Class' column indicating class labels for each point.
#   buf: Numeric; initial buffer distance to expand the tile (in units of the CRS).
#   coef: Numeric; coefficient to multiply the buffer distance in each iteration (default: 1.2).
# Details:
#   1. The function counts the number of unique classes in the point dataset (pt$Class).
#   2. It checks if the tile contains points from all unique classes using tilePoint_cpp to subset points within the tile's geometry.
#   3. If the tile already contains points from all classes, its original geometry is returned.
#   4. If not, the tile is iteratively expanded by multiplying the buffer distance by coef and applying a buffer with square end caps and mitered joins until all classes are included.
#
# Returns:
#   An sf POLYGON object representing the geometry of the original or expanded tile.

tile.multitasks.expand.untill.cover = function(tile, pt, buf, coef){

  n.class = unique(pt$Class) %>% length() # Count the number of class
  pt.subset = tilePoint_cpp(sf::st_coordinates(tile), sf::st_coordinates(pt)) # Subset the pt by the tile
  pt.subset = unique(pt$Class[pt.subset$id]) %>% length() # Count the number of class from the pt.subset

  if(pt.subset == n.class){
    # If the tile already contains points from all unique classes, return its geometry
    return(sf::st_geometry(tile))}  # <---- OUTPUT: sf POLYGON feature

  # If the tile does not contain all unique classes, expand it iteratively
  if(pt.subset != n.class){
    # Loop until all classes are covered
    while (TRUE) {
      buf = buf * coef # Increase the buffer distance by the coefficient

      # Expand the tile by applying a buffer with specified parameters
      # endCapStyle = "SQUARE" ensures square ends for linear features
      # joinStyle = "MITRE" ensures sharp corners
      # nQuadSegs = 5 controls the number of segments for curves
      # mitreLimit = 2 limits the extension of mitre joins
      tile = sf::st_buffer(x = tile, dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2) # Expand the tile

      # Subset points within the expanded tile's geometry
      pt.subset = tilePoint_cpp(sf::st_coordinates(tile), sf::st_coordinates(pt))

      # Count the number of unique classes in the new subset
      pt.subset = unique(pt$Class[pt.subset$id]) %>% length()

      # If all unique classes are now included, exit the loop
      if (pt.subset == n.class) {
        break }
    }
    # Return the geometry of the expanded tile
    return(sf::st_geometry(tile))  # <---- OUTPUT: sf POLYGON feature
  }

}



##_____________________________________________________________________________
#' Split Vector, Dataframe or sf feature into n groups
#'
#' This function partitions \code{data.frame}, \code{vector}, or \code{sf} feature into n groups
#' @param x  \code{sf} object or \code{vector} or \code{data.frame}. The data to be split.
#' @param n  \code{integer}. Non-negative value with a length of 1 to indicate the number of groups the data is to be split.
#' @return Return a list with a length of n, where each element contains a partition of the input data.
#' @examples
#' # Example 1: Split a numeric vector into 3 groups
#' v = 1:10
#' split_v = split.data(v, 3)
#' print(split_v)
#'
#' # Example 2: Split a data.frame into 4 groups
#' df = data.frame(a = 1:8, b = letters[1:8])
#' split_df = split.data(df, 4)
#' print(lapply(split_df, nrow)) # Each group's size
#'
#' # Example 3: Split an sf object into 2 groups and plot
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   library(sf)
#'   nc = st_read(system.file("shape/nc.shp", package="sf"), quiet=TRUE)
#'   split_nc = split.data(nc, 2)
#'   plot(st_geometry(split_nc[[1]]), col='red')
#'   plot(st_geometry(split_nc[[2]]), col='blue', add=TRUE)
#' }
#' @export

split.data = function(x, n){
  # --- Input validation ---
  if (!inherits(x, "sf") && !inherits(x, "data.frame") && !is.vector(x)) {
    stop("Invalid input for x: Only accept <sf feature>, <vector>, or <data.frame>.") }

  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("Invalid input for n: n must be a single positive integer.")
  }

  # --- Main function implementation ---
  # Get the length of the x
  if(inherits(x, "sf") || inherits(x, "data.frame")){len = nrow(x)}
  if(is.vector(x)){len = length(x)}

  # Find the maximum number of each group
  n_per_split = ceiling(len / n)

  # Split the data
  return( split(x, (seq_len(len) - 1) %/% n_per_split) )  # <---- OUTPUT: list
}






















