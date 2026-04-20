##_____________________________________________________________________________
#' Find the Nearest Feature in 3D Space
#'
#' Identifies the index of the nearest 3D point or linestring for each input 3D point.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry, representing the query points.
#' @param y An \code{sf} object with 3D \code{POINT} or \code{LINESTRING} geometry, representing the target features.
#'
#' @details
#' The function finds the nearest feature in \code{y} for each 3D point in \code{x} using 3D Euclidean distance, supporting two scenarios:
#' \itemize{
#'   \item \code{POINT}-to-\code{POINT}: Identifies the index of the nearest point in \code{y} for each point in \code{x}.
#'   \item \code{POINT}-to-\code{LINESTRING}: Identifies the index of the nearest linestring in \code{y} for each point in \code{x}, based on the closest point on the linestring.
#' }
#'
#' @return A \code{numeric} vector of length \code{nrow(x)}, where each value is the index of the nearest feature in \code{y} for the corresponding point in \code{x}.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Example: point-to-point
#' pt1 = st_sf(geometry = st_sfc(
#'   st_point(c(0,0,0)),
#'   st_point(c(1,0,0)),
#'   st_point(c(0,1,1))
#' ))
#'
#' pt2 = st_sf(geometry = st_sfc(
#'   st_point(c(0,0,1)),
#'   st_point(c(1,1,2)),
#'   st_point(c(0,1,4.5))
#' ))
#'
#' print(nearest3d(pt1, pt2))
#'
#' # Example: point-to-line
#' pt = st_sf(geometry = st_sfc(
#'   st_point(c(0,0,1)),
#'   st_point(c(1,1,2)),
#'   st_point(c(0,1,4.5)) )
#'   )
#' lines = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(2,0,0, 2,2,2), ncol=3, byrow=TRUE))
#' ))
#'
#' print(nearest3d(pt, lines))
#' }
#'
#' @export

nearest3d = function(x, y){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with POINT geometry.") }

  if ( !inherits(y, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(y)) & !"sfc_LINESTRING" %in% class(sf::st_geometry(y)) ) {
    stop("y must be an <sf> object with POINT or LINESTRING geometry.") }

  # --- Main function implementation ---
  ## POINT to POINT operation
  if( "sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_POINT" %in% class(sf::st_geometry(y)) ){
    x = st_coordinates(x)
    y = st_coordinates(y)
    return(nearest3d_p2p_cpp(x,y))}  # <---- OUTPUT: numeric vector

  ## POINT to LINE operation
  if( "sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_LINESTRING" %in% class(sf::st_geometry(y)) ){
    y = y %>% st_coordinates() %>% segmentise_cpp()
    fromIndex = seq(1, length(y$L1), 2)
    toIndex = seq(2, length(y$L1), 2)
    y = data.frame(from_x = y$X[fromIndex], from_y = y$Y[fromIndex], from_z = y$Z[fromIndex],
                   to_x = y$X[toIndex], to_y = y$Y[toIndex], to_z = y$Z[toIndex],
                   L1 = y$L1[fromIndex], L2 = y$L2[fromIndex])
    x = x %>% st_coordinates()
    return(nearest3d_p2l_cpp(x,y))}  # <---- OUTPUT: numeric vector

}

##_____________________________________________________________________________
#' Find the Closest Projected Point on a 3D Linestring for Each 3D Point
#'
#' Computes the closest projected point on a 3D linestring for each corresponding 3D input point, returning results as coordinates or an \code{sf} object.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry, representing the query points.
#' @param y An \code{sf} object with 3D \code{LINESTRING} geometry, with the same number of rows as \code{x}, each corresponding to a point in \code{x}.
#' @param output A \code{character} value specifying the output format: \code{"df"} (default, coordinates as a \code{data.frame}) or \code{"sf"} (3D \code{sf} \code{POINT} object).
#'
#' @details
#' For each pair of 3D point in \code{x} and 3D linestring in \code{y} (matched by row), the function computes the point on the linestring closest to the input point using 3D Euclidean distance (X, Y, Z). The number of rows in \code{x} and \code{y} must be equal for one-to-one pairing. The output preserves the original coordinate reference system (CRS) and Z coordinates for \code{"sf"} output.
#'
#' @return Depends on \code{output}:
#' \describe{
#'   \item{\code{"df"}}{A \code{data.frame} with columns \code{X}, \code{Y}, \code{Z} containing the coordinates of the projected points.}
#'   \item{\code{"sf"}}{An \code{sf} object with 3D \code{POINT} geometry containing the projected points.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(rgl)
#'
#' # Create example points
#' pts <- st_sf(geometry = st_sfc(
#'   st_point(c(1, 1, 1)),
#'   st_point(c(7, 7, 2.2))
#' ))
#'
#' # Create example lines
#' lines <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 2,0,0, 2,2,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(5,5,0, 6,6,0, 8,8,0), ncol=3, byrow=TRUE))
#' ))
#'
#' # Output as data.frame
#' projected_df <- closestPT3d(pts, lines, output = "df")
#' print(projected_df)
#'
#' # Output as sf POINT
#' projected_sf <- closestPT3d(pts, lines, output = "sf")
#' print(projected_sf)
#'
#' # 3D plot
#' open3d()
#'
#' # Plot the original lines (dashed gray)
#' lines_coords = st_coordinates(lines)
#' for (i in seq(1, nrow(lines_coords), by=3)) {
#'   lines3d(lines_coords[i:(i+2), 1], lines_coords[i:(i+2), 2], lines_coords[i:(i+2), 3], col="gray", lwd=2, lty=2)
#' }
#'
#' # Plot the original points (red)
#' pts_coords = st_coordinates(pts)
#' points3d(pts_coords[,1], pts_coords[,2], pts_coords[,3], color="red", size=12)
#'
#' # Plot the projected points (blue)
#' proj_coords = st_coordinates(projected_sf)
#' points3d(proj_coords[,1], proj_coords[,2], proj_coords[,3], color="blue", size=12)
#'
#' # Add a legend
#' legend3d("topright",
#'   legend = c("Original lines", "Input points", "Projected points"),
#'   lwd = c(2, NA, NA),
#'   col = c("gray", "red", "blue"),
#'   pch = c(NA, 16, 16)
#' )
#' }
#'
#' @export

closestPT3d = function(x, y, output = "df"){
  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with POINT geometry.") }

  if (!inherits(y, "sf") || !"sfc_LINESTRING" %in% class(sf::st_geometry(y))) {
    stop("y must be an <sf> object with LINESTRING geometry.") }

  if(nrow(x) != nrow(y)){
    stop("x and y must have same length.")}

  if(is.character(output) == F | length(output) > 1){
    stop("output must be a character of either 'df' or 'sf'.")}

  if (!output %in% c("df", "sf")) {
    stop("output must be a character of either 'df' or 'sf'.")}

  # --- Main function implementation ---
  OutputCRS = st_crs(x) # Define crs
  y = y %>% st_coordinates() %>% as.data.frame() # Extract the coordinates from y
  x = x %>% st_coordinates() %>% as.data.frame() # Extract the coordinates from x
  if(output == "df"){
    return(closestPT3d_cpp(x,y))}  # <---- OUTPUT: data.frame
  if(output == "sf"){
    result = closestPT3d_cpp(x,y) %>% sfheaders::sf_point(x = "X", y = "Y", z = "Z")
    st_crs(result) = OutputCRS
    return(result)  # <---- OUTPUT: sf POINT feature
  }
}


##_____________________________________________________________________________
#' Detect Geometries Within a 3D Distance Threshold
#'
#' Identifies pairs of 3D geometries within a specified 3D Euclidean distance threshold.
#'
#' @param x An \code{sf} object with 3D \code{POINT} or \code{LINESTRING} geometry, representing the query geometries.
#' @param y \code{NULL} or an \code{sf} object with 3D \code{POINT} or \code{LINESTRING} geometry, representing the target geometries. If \code{NULL}, compares all pairs within \code{x}, excluding self-pairs.
#' @param dist A non-negative \code{numeric} value specifying the maximum 3D distance threshold.
#'
#' @details
#' The function identifies pairs of 3D geometries within \code{dist} using 3D Euclidean distance, supporting:
#' \itemize{
#'   \item \code{POINT}-to-\code{POINT} (\code{y = NULL}): Finds all pairs of points in \code{x} within \code{dist}, excluding self-matches.
#'   \item \code{POINT}-to-\code{POINT}: Finds all points in \code{y} within \code{dist} of each point in \code{x}.
#'   \item \code{LINESTRING}-to-\code{POINT}: Finds all points in \code{y} within \code{dist} of the closest point on each linestring in \code{x}.
#'   \item \code{POINT}-to-\code{LINESTRING}: Finds all linestrings in \code{y} within \code{dist} of the closest point to each point in \code{x}.
#' }
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{\code{id1}}{\code{integer}, index of the feature in \code{x}.}
#'   \item{\code{id2}}{\code{integer}, index of the feature in \code{y} (or \code{x} if \code{y = NULL}).}
#'   \item{\code{dist}}{\code{numeric}, 3D Euclidean distance between the pair.}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' # Create some 3D points
#' pts1 = st_sf(geometry = st_sfc(st_point(c(0,0,0)), st_point(c(1,0,0)), st_point(c(0,1,1))))
#' pts2 = st_sf(geometry = st_sfc(st_point(c(0,0,0)), st_point(c(2,0,0)), st_point(c(0,2,2))))
#' # Find all pairs within 1.5 units
#' withinDist3d(pts1, pts2, dist = 1.5)
#' # Find all pairs within 2 units in one set
#' withinDist3d(pts1, dist = 2)
#' }
#'
#' @export

withinDist3d = function(x, y = NULL, dist){
  # --- Input validation ---
  if (!inherits(x, "sf") ||
      (!"sfc_POINT" %in% class(sf::st_geometry(x)) & (!"sfc_LINESTRING" %in% class(sf::st_geometry(x)))) ) {
    stop("x must be an <sf> object with POINT or LINESTRING geometry.") }

  if(!is.null(y)){
    if (!inherits(y, "sf") ||
        (!"sfc_POINT" %in% class(sf::st_geometry(y)) & (!"sfc_LINESTRING" %in% class(sf::st_geometry(y)))) ) {
      stop("y must be NULL or an <sf> object with POINT or LINESTRING geometry.")}

    if( "sfc_LINESTRING" %in% class(sf::st_geometry(x)) & "sfc_LINESTRING" %in% class(sf::st_geometry(y)) ){
      stop("Currently does not support line-to-line operation.")}
  }

  if(is.null(y)){
    if(!"sfc_POINT" %in% class(sf::st_geometry(x))){
      stop("x must be an <sf> object with POINT geometry if y is null.")
    }
  }

  if (!is.numeric(dist) || length(dist) != 1 || dist < 0) {
    stop("dist must be a non-negative numeric value.")}

  # --- Main function implementation ---
  #->>> x = POINT & y = NULL
  if(is.null(y)){
    if( "sfc_POINT" %in% class(sf::st_geometry(x))){
      # Remove self and find all points within predefined distance from the input points
      return(withinDist3d_p_self_cpp(sf::st_coordinates(x), dist = dist))}   # <---- OUTPUT: data.frame
  } else {
    #->>> x = POINT & y = POINT
    if( "sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_POINT" %in% class(sf::st_geometry(y)) ){
      # Find all points within predefined distance from the input points
      return(withinDist3d_p2p_cpp(sf::st_coordinates(x), sf::st_coordinates(y), dist = dist))}   # <---- OUTPUT: data.frame

    #->>> x = LINESTRING & y = POINT
    else if( "sfc_LINESTRING" %in% class(sf::st_geometry(x)) & "sfc_POINT" %in% class(sf::st_geometry(y)) ){
      # Find all points within predefined distance from the input lines
      result = withinDist3d_l2p_cpp(sf::st_coordinates(x), sf::st_coordinates(y), dist = dist)
      names(result) = c("id1", "id2", "dist")
      return(result)}   # <---- OUTPUT: data.frame

    #->>> x = POINT & y = LINESTRING
    else if( "sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_LINESTRING" %in% class(sf::st_geometry(y)) ){
      result = withinDist3d_l2p_cpp(sf::st_coordinates(y), sf::st_coordinates(x), dist = dist)
      result = result[order(result$PTid),]
      return( data.frame(id1 = result$PTid, id2 = result$L1, dist = result$dist)) }  # <---- OUTPUT: data.frame
  }

}

##_____________________________________________________________________________
#' Find the K-Nearest Neighbours in 3D Space
#'
#' Identifies the \code{k} nearest 3D points or linestrings for each input 3D point.
#'
#' @param x An \code{sf} object with 3D \code{POINT} geometry, representing the query points.
#' @param y \code{NULL} or an \code{sf} object with 3D \code{POINT} or \code{LINESTRING} geometry, representing the target features. If \code{NULL}, finds neighbours within \code{x}.
#' @param k A positive \code{integer} specifying the number of nearest neighbours to retrieve for each point in \code{x}.
#'
#' @details
#' The function identifies the \code{k} nearest neighbours for each 3D point in \code{x} using 3D Euclidean distance, supporting:
#' \itemize{
#'   \item \code{POINT}-to-\code{POINT} (\code{y = NULL}): Finds the \code{k} nearest points in \code{x}, excluding self-matches.
#'   \item \code{POINT}-to-\code{POINT}: Finds the \code{k} nearest points in \code{y} for each point in \code{x}.
#'   \item \code{POINT}-to-\code{LINESTRING}: Finds the \code{k} nearest linestrings in \code{y} for each point in \code{x}, based on the closest point on each linestring.
#' }
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{\code{id1}}{\code{integer}, index of the query point in \code{x}.}
#'   \item{\code{id2}}{\code{integer}, index of the matched feature in \code{y} (or \code{x} if \code{y = NULL}).}
#'   \item{\code{dist}}{\code{numeric}, 3D Euclidean distance between the features.}
#' }
#' Rows are ordered by \code{id1} and \code{dist}.
#'
#' @examples
#' \dontrun{
#' library(sf)
#' pts <- st_sf(geometry = st_sfc(
#'   st_point(c(0,0,0)),
#'   st_point(c(1,0,0)),
#'   st_point(c(0,1,1)),
#'   st_point(c(0,2,0))
#' ))
#' lines <- st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(2,0,0, 2,2,2), ncol=3, byrow=TRUE))
#' ))
#' # 2 nearest neighbours within pts
#' knn3d(pts, k = 2)
#' # 2 nearest neighbours from pts to lines
#' knn3d(pts, lines, k = 2)
#' }
#'
#' @export

knn3d = function(x, y = NULL, k){

  # --- Input validation ---
  if (!inherits(x, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(x))) {
    stop("x must be an <sf> object with POINT geometry.") }

  if(!is.null(y)){
    if (!inherits(y, "sf") ||
        (!"sfc_POINT" %in% class(sf::st_geometry(y)) & (!"sfc_LINESTRING" %in% class(sf::st_geometry(y)))) ) {
      stop("y must be NULL or an <sf> object with POINT or LINESTRING geometry.")}}

  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != round(k)) {
    stop("k must be a positive integer.")}

  # --- Main function implementation ---
  #->>> x = POINT & y = NULL
  if(is.null(y)){
    if("sfc_POINT" %in% class(sf::st_geometry(x))){
      result = knn3d_p2p_cpp(sf::st_coordinates(x), sf::st_coordinates(x), k+1) # Find the k+1 nearest neighbours for each point
      result = result %>% subset(id1 != id2) # Remove the self from the nearest neighbours
      result = result[order(result$id1, result$dist), ] # Order by ID and distance
      return(result) }   # <---- OUTPUT: data.frame
  } else {
    #->>> x = POINT & y = POINT
    if( "sfc_POINT" %in% class(sf::st_geometry(x)) & "sfc_POINT" %in% class(sf::st_geometry(y)) ){
      result = knn3d_p2p_cpp(sf::st_coordinates(x), sf::st_coordinates(y), k)
      result = result[order(result$id1, result$dist), ] # Order by ID and distance
      return(result)}  # <---- OUTPUT: data.frame
    #->>> x = POINT & y = LINES
    else {
      return(knn3d.noIndex(x,y,k))   # <---- OUTPUT: data.frame
    }
  }
}

# -----> An inner function to find the indices with the k least values.
# --> value = A <numeric> vector
# --> index = A <numeric> vector to indicate the index for the value
# --> k = An <integer> to indicate the number of least values to be considered

vector.knn = function(value, index, k){
  index = index[order(value)][1:k]
  value = value[index]
  return(list(value = value, index = index))
}

# -----> An inner function to compute the knn for POINT-to-LINE
# --> x = <sf POINT feature>.
# --> y = <sf POINT feature> OR <sf LINESTRING feature> OR <NULL>. If null, find the k nearest points for each point in x.
# --> k = An <integer> to indicate the number of least values to be considered

knn3d.noIndex = function(x,y,k){

  Xindex = 1:nrow(x) # Get the index for x
  Yindex = 1:nrow(y) # Get the index for y
  distMatrix = dist3d(x,y) %>% as.data.frame() # Calculate the distance matrix

  # Create container to store date
  id1 = integer(0) # Index for the x
  id2 = integer(0) # Index for the y
  dist = numeric(0) # Distance

  for (i in Xindex) {
    result = vector.knn(distMatrix[i,] %>% as.numeric(), Yindex, k) # For each row of the OD matrix, find the indices of the k least values
    id1 = c(id1, rep(i, k)) # Index for point
    id2 = c(id2, result$index) # Index for line
    dist = c(dist, result$value) # Distance
  }
  result = data.frame(id1 = id1, id2 = id2, dist = dist) # Create data.frame
  result = result[order(result$id1, result$dist), ] # Order by ID and distance
  return(result)   # <---- OUTPUT: data.frame
}























