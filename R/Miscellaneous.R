
# ___________________________________________________________________________________________________________________________________

tile.create = function(x, nrows = 4, ncols = 4, returnPT = FALSE){

  # Create tile
  tiles = x %>% sf::st_make_grid(n = c(ncols, nrows)) %>% sf::st_as_sf()
  tiles$grid = 1:nrow(tiles)

  # Only retain tiles that containing points
  intersectList = sf::st_intersects(tiles, x)
  intersectList = purrr::map2(intersectList, 1:length(intersectList),
                              function(x,y){
                                if(length(x)==0){ return(NULL)
                                } else { x = data.frame(id = x, grid = y) }
                                return(x) } ) %>% purrr::reduce(rbind)
  intersectList = intersectList[duplicated(intersectList$id) == F,]
  tiles = tiles[intersectList$grid %>% unique(),]

  # Assign Grid ID
  tiles$gridID = 1:nrow(tiles)
  intersectList = intersectList %>% dplyr::left_join(data.frame(grid = tiles$grid, gridID = tiles$gridID), by = "grid")

  # Export Result
  if(isTRUE(returnPT)){
    # Extract Geometry
    intersectList = split(x[intersectList$id,], intersectList$gridID)
  } else {
    if("nID" %in% names(x)){
      # Find nID
      intersectList$nID = x$nID[intersectList$id]
      intersectList = split(intersectList[c("id", "nID")], intersectList$gridID)
    } else {
      # Only retain id
      intersectList = split(intersectList$id, intersectList$gridID) }
  }

  return(list(tiles = tiles["gridID"], pt = intersectList)) # <---- Output
}


# ___________________________________________________________________________________________________________________________________

tile.point = function(tiles, pt, buf = NULL, coef = 1.2, returnPT = FALSE, returnTiles = FALSE){

  # Set buf = 0 if NULL
  if(is.null(buf)){buf = 0}

  # Expand tiles
  tiles = sf::st_buffer(x = tiles, dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2)

  # Find points within tiles
  intersectList = sf::st_intersects(tiles, pt)

  # Identify Tiles requiring expansion
  needExpand = which(unlist(purrr::map(intersectList, length))==0)

  # Find at least one points
  while (length(needExpand) > 0){
    buf = buf * coef # Increase the buffer size

    # Expand Tiles
    sf::st_geometry(tiles)[needExpand] = sf::st_geometry(tiles)[needExpand] %>% sf::st_buffer(dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2)
    intersectList.expanded = sf::st_intersects(sf::st_geometry(tiles)[needExpand], pt) # Find points within tiles
    intersectList.expanded = which(unlist(purrr::map(intersectList.expanded, length))>0) # Identify Tiles requiring expansion

    if(length(intersectList.expanded) != 0){
      needExpand = needExpand[-intersectList.expanded] } }

  # Final: Find points within tiles
  intersectList = sf::st_intersects(tiles, pt)

  # Add nID if it exist in the original pt
  if("nID" %in% names(pt) & isFALSE(returnPT)){ intersectList = purrr::map(intersectList, ~ data.frame(id = .x, nID = pt$nID[.x])) }

  # Add Geometry
  if(isTRUE(returnPT)){ intersectList = purrr::map(intersectList, ~ pt[.x,] ) }

  # Export data
  if(isFALSE(returnTiles)){return(intersectList)} # <--- Output
  if(isTRUE(returnTiles)){return(list(tiles = tiles, pt = intersectList))}  # <--- Output
}

# ___________________________________________________________________________________________________________________________________
tile.OD = function(o, d, buf, nrows = 4, ncols = 4, coef = 1.2, returnPT = FALSE, returnTiles = FALSE){

  # Create tiles from the bounding box of the origin points `o`
  tiles = GISnetwork3D::tile.create(o, nrows = nrows, ncols = ncols, returnPT = returnPT)

  # Subtract destination points for each tile by expanding the tiles using the buffer distance `buf`
  destination = GISnetwork3D::tile.point(tiles = tiles$tiles, pt = d, buf = buf, coef = coef, returnPT = returnPT, returnTiles = returnTiles)

  # Return a list containing stratified origin and destination points
  if(isTRUE(returnTiles)){ return(list(o = tiles$pt, d = destination$pt, o.tiles = tiles$tiles, d.tiles = destination$tiles)) } # <---- OUTPUT
  if(isFALSE(returnTiles)){ return(list(o = tiles$pt, d = destination)) } # <---- OUTPUT
}

# ___________________________________________________________________________________________________________________________________
tile.ODN = function(o, d, n, buf, nrows = 4, ncols = 4, coef = 1.2, returnPT = FALSE, returnTiles = FALSE){

  # Create tiles from the bounding box of the origin points `o`
  tiles = GISnetwork3D::tile.create(o, nrows = nrows, ncols = ncols, returnPT = returnPT)

  # Subtract destination points for each tile by expanding the tiles using the buffer distance `buf`
  destination = GISnetwork3D::tile.point(tiles = tiles$tiles, pt = d, buf = buf, coef = coef, returnPT = returnPT, returnTiles = TRUE)

  # Extract Network Within Tile
  Nodes = GISnetwork3D::tile.point(tiles = destination$tiles, pt = n$nodes, buf = buf, coef = coef)
  n = purrr::map(Nodes, ~ GISnetwork3D::subsetNet3d(net = n, v = .x$nID, exact = F, output = "sf"))

  # Return a list containing stratified origin and destination points
  if(isTRUE(returnTiles)){ return(list(o = tiles$pt, d = destination$pt, n = n, o.tiles = tiles$tiles, d.tiles = destination$tiles)) } # <---- OUTPUT
  if(isFALSE(returnTiles)){ return(list(o = tiles$pt, d = destination$pt, n = n)) } # <---- OUTPUT
}

# ___________________________________________________________________________________________________________________________________
tile.multitasks.point = function(tiles, pt, buf = NULL, coef = 1.2, returnPT = FALSE, returnTiles = FALSE){

  # Set buf = 0 if NULL
  if(is.null(buf)){buf = 0}

  # Expand tiles
  tiles = sf::st_buffer(x = tiles, dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2)

  # Count number of Class
  nClass = pt$Class %>% unique() %>% length()

  # Find points within tiles
  intersectList = sf::st_intersects(tiles, pt)
  intersectList = purrr::map(intersectList, ~ unique(pt$Class[.x]) %>% length)

  # Identify Tiles requiring expansion
  needExpand = which(unlist(intersectList)!=nClass)

  # Find at least one points
  while (length(needExpand) > 0){
    buf = buf * coef
    sf::st_geometry(tiles)[needExpand] = sf::st_geometry(tiles)[needExpand] %>% sf::st_buffer(dist = buf, endCapStyle = "SQUARE", joinStyle = "MITRE", nQuadSegs = 5, mitreLimit = 2)
    intersectList.expanded = sf::st_intersects(sf::st_geometry(tiles)[needExpand], pt)
    intersectList.expanded = purrr::map(intersectList.expanded, ~ unique(pt$Class[.x]) %>% length)

    intersectList.expanded = which(unlist(intersectList.expanded)==nClass)

    if(length(intersectList.expanded) != 0){
      needExpand = needExpand[-intersectList.expanded] }
  }

  # Final: Find points within tiles
  intersectList = sf::st_intersects(tiles, pt)

  # Add nID if it exist in the original pt
  if("nID" %in% names(pt) & isFALSE(returnPT)){ intersectList = purrr::map(intersectList, ~ data.frame(id = .x, nID = pt$nID[.x])) }

  # Add Geometry
  if(isTRUE(returnPT)){ intersectList = purrr::map(intersectList, ~ pt[.x,] ) }

  # Export data
  if(isFALSE(returnTiles)){return(intersectList)} # <--- Output
  if(isTRUE(returnTiles)){return(list(tiles = tiles, pt = intersectList))}  # <--- Output
}


# ___________________________________________________________________________________________________________________________________
tile.multitasks.OD = function(o, d, buf, nrows = 4, ncols = 4, coef = 1.2, returnTiles = FALSE, returnPT = FALSE){

  # Create tiles from the bounding box of the origin points `o`
  tiles = GISnetwork3D::tile.create(o, nrows = nrows, ncols = ncols, returnPT = returnPT)

  # Subtract destination points for each tile by expanding the tiles using the buffer distance `buf`
  destination = GISnetwork3D::tile.multitasks.point(tiles = tiles$tiles, pt = d, buf = buf, coef = coef, returnPT = returnPT, returnTiles = returnTiles)

  # Return a list containing stratified origin and destination points
  if(isTRUE(returnTiles)){ return(list(o = tiles$pt, d = destination$pt, o.tiles = tiles$tiles, d.tiles = destination$tiles)) } # <---- OUTPUT
  if(isFALSE(returnTiles)){ return(list(o = tiles$pt, d = destination)) }
}

# ___________________________________________________________________________________________________________________________________
tile.multitasks.ODN = function(o, d, n, buf, nrows = 4, ncols = 4, coef = 1.2, returnPT = FALSE, returnTiles = FALSE){

  # Create tiles from the bounding box of the origin points `o`
  tiles = GISnetwork3D::tile.create(o, nrows = nrows, ncols = ncols, returnPT = returnPT)

  # Subtract destination points for each tile by expanding the tiles using the buffer distance `buf`
  destination = GISnetwork3D::tile.multitasks.point(tiles = tiles$tiles, pt = d, buf = buf, coef = coef, returnPT = returnPT, returnTiles = TRUE)

  # Extract Network Within Tile
  Nodes = GISnetwork3D::tile.point(tiles = destination$tiles, pt = n$nodes, buf = buf, coef = coef)
  n = purrr::map(Nodes, ~ GISnetwork3D::subsetNet3d(net = n, v = .x$nID, exact = F, output = "sf"))

  # Return a list containing stratified origin and destination points
  if(isTRUE(returnTiles)){ return(list(o = tiles$pt, d = destination$pt, n = n, o.tiles = tiles$tiles, d.tiles = destination$tiles)) } # <---- OUTPUT
  if(isFALSE(returnTiles)){ return(list(o = tiles$pt, d = destination$pt, n = n)) } # <---- OUTPUT
}

# ___________________________________________________________________________________________________________________________________
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






















