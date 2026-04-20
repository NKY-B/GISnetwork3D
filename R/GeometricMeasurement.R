##_____________________________________________________________________________
#' Calculate 3D Distance Matrix Between 3D Points and 3D Points/Lines
#'
#' This function computes the 3D Euclidean distance from each point in \code{x} to each feature in \code{y}. It supports both point-to-point and point-to-linestring distances (many-to-many).
#'
#' @param x An \code{sf} object with \code{POINT} geometry. The origin points.
#' @param y An \code{sf} object with either \code{POINT} or \code{LINESTRING} geometry. The destination points or linestrings.
#'
#' @details
#' When \code{y} is a set of points, the function calculates the pairwise 3D Euclidean distance between every point in \code{x} and every point in \code{y}.
#' When \code{y} is a set of linestrings, the function calculates the shortest 3D distance from each point in \code{x} to each linestring in \code{y}.
#' The returned matrix has \code{nrow(x)} rows and \code{nrow(y)} columns.
#'
#' @return A numeric matrix of 3D distances, where:
#'   - Rows correspond to each geometry in \code{x}.
#'   - Columns correspond to each geometry in \code{y}.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Example 1: 3D point-to-point distances
#' x = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 0)),
#'   st_point(c(2, 0, 2)),
#'   st_point(c(1, 2, 1))
#' ))
#' y = st_sf(geometry = st_sfc(
#'   st_point(c(0, 1, 0)),
#'   st_point(c(2, 2, 2))
#' ))
#' dmat = dist3d(x, y)
#' print(dmat)
#'
#' # Example 2: 3D point-to-line distances
#' y_line = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 2,2,2, 4,0,2), ncol=3, byrow=TRUE))
#' ))
#' dmat2 = dist3d(x, y_line)
#' print(dmat2)
#'}
#' @export

dist3d = function(x,y){
  # --- Main function implementation ---
  if("sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_POINT" %in% class(sf::st_geometry(y))){
    return(dist3d_pt2pt_cpp(sf::st_coordinates(x),st_coordinates(y)))}  # <---- OUTPUT: matrix
  else if("sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_LINESTRING" %in% class(sf::st_geometry(y))){
    x = st_coordinates(x)
    y = y %>% st_as_sf() %>% st_coordinates() %>% as.data.frame()
    y = y %>% dplyr::select(X,Y,Z,L1)
    return(dist3d_pt2line_cpp(x,y))}  # <---- OUTPUT: matrix
  else{
    stop("Only accept an <sf> object with POINT geometry for x and an <sf> object with LINESTRING geometry for y.")
  }
}

##_____________________________________________________________________________
#' Calculate 3D Distance Between Corresponding Pairs of 3D Points and 3D Points or Linestrings
#'
#' Computes the 3D Euclidean distance between each corresponding pair of geometries in \code{x} and \code{y}. Supports both point-to-point and point-to-linestring distances.
#' Both \code{x} and \code{y} must have the same number of features; the first geometry in \code{x} is paired with the first geometry in \code{y}, the second with the second, and so on.
#'
#' @param x An \code{sf} object with \code{POINT} geometry. The origin points.
#' @param y An \code{sf} object with either \code{POINT} or \code{LINESTRING} geometry. The destination points or linestrings. Must have the same number of features as \code{x}.
#'
#' @details
#' This function calculates the 3D Euclidean distance between corresponding pairs of geometries in \code{x} and \code{y}.
#' For point-to-point distances, it computes the straight-line distance between each pair of points.
#' For point-to-linestring distances, it computes the shortest distance from each point in \code{x} to the corresponding linestring in \code{y}.
#' The inputs \code{x} and \code{y} must have the same number of rows.
#'
#' @return A numeric vector of 3D distances, where each element represents the distance between the corresponding pair of geometries in \code{x} and \code{y}.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Example 1: Point-to-point distances
#' x = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 0)),
#'   st_point(c(2, 0, 1)),
#'   st_point(c(1, 1, 3))
#' ))
#' y = st_sf(geometry = st_sfc(
#'   st_point(c(1, 1, 1)),
#'   st_point(c(2, 1, 2)),
#'   st_point(c(1, 2, 4))
#' ))
#' dists = dist3d.pair(x, y)
#' print(dists)
#'
#' # Example 2: Point-to-linestring distances
#' x = st_sf(geometry = st_sfc(
#'   st_point(c(1, 1, 2.5)),
#'   st_point(c(2, 0, 4)),
#'   st_point(c(4, 1, 3))
#' ))
#'
#' y_line = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 5,5,0), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(2,0,1, 2,2,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,3, 2,2,4), ncol=3, byrow=TRUE))
#' ))
#' dists2 <- dist3d.pair(x, y_line)
#' print(dists2)
#'
#' }
#' @export

dist3d.pair = function(x,y){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("l must be an <sf> object with POINT geometry>.") }

  if ( !inherits(y, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(y)) & !"sfc_LINESTRING" %in% class(sf::st_geometry(y)) ) {
    stop("y must be an <sf> object with POINT OR LINESTRING geometry.") }

  if(nrow(x) != nrow(y)){
    stop("x and y must have same length.")
  }

  # --- Main function implementation ---
  ## POINT to POINT operation
  if( "sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_POINT" %in% class(sf::st_geometry(y)) ){
    x = st_coordinates(x)
    y = st_coordinates(y)
    return(dist3d_pt2pt_pair_cpp(x,y))}  # <---- OUTPUT: numeric vector

  ## POINT to LINE operation
  if( "sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_LINESTRING" %in% class(sf::st_geometry(y)) ){
    x = st_coordinates(x)
    y = as.data.frame(sf::st_coordinates(y))
    return(dist3d_pt2line_pair_cpp(x,y))}  # <---- OUTPUT: numeric vector
}



##_____________________________________________________________________________
#' Calculate 3D Length of 3D Linestrings
#'
#' This function calculates the true 3D length for each linestring in an \code{sf} object with \code{LINESTRING} geometry, using Euclidean distance in 3D space.
#'
#' @param x An \code{sf} object with 3D \code{LINESTRING} geometry. The input linestrings.
#'
#' @details
#' The 3D length is computed as the sum of straight-line Euclidean distances between consecutive vertices of each 3D linestring, considering all three dimensions (X, Y, Z). The input \code{x} must contain 3D \code{LINESTRING} geometries with Z coordinates.
#'
#' @return A numeric vector containing the 3D length of each linestring in \code{x}.
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Create 3D linestrings
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 0,0,1, 0,0,5.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 3,1,2, 5,2,4), ncol=3, byrow=TRUE))
#' ))
#'
#' # Calculate their 3D lengths
#' len = len3d(l)
#' print(len)
#'
#' # Visualize linestrings in 3D, colored by length
#' open3d()
#' cols = c("blue", "red")
#' for(i in seq_len(nrow(l))) {
#'   coords = st_coordinates(l[i, ])
#'   lines3d(coords[,1], coords[,2], coords[,3], color=cols[i], lwd=5)
#'   text3d(mean(coords[,1]), mean(coords[,2]), mean(coords[,3]),
#'          paste0("Length = ", round(len[i],2)), color=cols[i], adj=1)
#' }
#' legend3d("topright", legend=paste("Line", 1:2), col=cols, lwd=5)
#' }
#' @export

len3d = function(x){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with LINESTRING geometry.") }

  # --- Main function implementation ---
  return(len3d_cpp(as.data.frame(sf::st_coordinates(x)))) # <---- OUTPUT: numeric vector
}




##_____________________________________________________________________________
#' Calculate the Slope of a 3D sf LINESTRING object
#'
#' This function calculates the slope between the first and last nodes of each 3D linestring in an \code{sf} object with 3D \code{LINESTRING} geometry. The slope can be expressed in different units: percentage, degrees, or gradient.
#'
#' @param l An \code{sf} object with 3D \code{LINESTRING} geometry. The input linestrings.
#' @param unit A \code{character} specifying the unit for the slope. Acceptable values are "percent", "degrees", or "gradient".
#'
#' @details
#' The slope is calculated based on the straight-line distance between the start and end nodes of each 3D linestring, considering all three dimensions (X, Y, Z).
#' The input \code{l} must contain 3D \code{LINESTRING} geometries with Z coordinates. This function works best for linestrings with only two points (start and end nodes) without intermediate vertices.
#' The function uses \code{GISnetwork3D::getNode3d} to extract the 3D coordinates of the first and last nodes.
#'
#' @return A numeric vector containing the slope between the first and last nodes of each linestring in \code{l}, expressed in the specified unit.
#'
#'
#' @examples
#' library(sf)
#' library(GISnetwork3D)
#'
#' # Create a simple sf LINESTRING object
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 0,0,1, 0,0,5.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 3,1,2, 5,2,4), ncol=3, byrow=TRUE))
#' ))
#' # Calculate slope in percent
#' slope_percent = calSlope(l, unit = "percent")
#' print(slope_percent)
#'
#' # Calculate slope in degrees
#' slope_degrees = calSlope(l, unit = "degrees")
#' print(slope_degrees)
#'
#' # Calculate slope in gradient
#' slope_gradient = calSlope(l, unit = "gradient")
#' print(slope_gradient)
#'
#' @export
calSlope = function(l, unit){

  # --- Input validation ---
  if (!inherits(l, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(l))) {
    stop("l must be an <sf> object with LINESTRING geometry>.") }

  if (length(unit) != 1) {
    stop("unit must be a character of either 'percent', 'degrees' or 'gradient'.")
  }

  if (!unit %in% c("percent", "degrees", "gradient")) {
    stop("unit must be a character of either 'percent', 'degrees' or 'gradient'.")
  }

  # --- Main function implementation ---
  return( calSlope_net3d_cpp(sf::st_coordinates( GISnetwork3D::getNode3d(l, position="first") ),
                             sf::st_coordinates( GISnetwork3D::getNode3d(l, position="last") ),
                             1:nrow(sf::st_coordinates( GISnetwork3D::getNode3d(l, position="last") )),
                             unit) ) # <---- OUTPUT: numeric vector
}






