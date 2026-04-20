##_____________________________________________________________________________
#' Calculate Least-Cost Origin-Destination Matrix on a 3D Network
#'
#' Computes a least-cost origin-destination (OD) matrix between 3D points on a 3D network, using specified edge weights and directionality.
#'
#' @param net An \code{igraph} object or a \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param o Origins: 3D \code{sf POINT} object with \code{nID} column or a \code{character} vector of node IDs.
#' @param d Destinations: 3D \code{sf POINT} object with \code{nID} column or a \code{character} vector of node IDs.
#' @param mode A \code{character} specifying directionality: \code{"out"} (origin to destination), \code{"in"} (destination to origin), or \code{"both"} (sum of both directions).
#' @param weight A \code{character} specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#'
#' @details
#' The function calculates the least-cost paths between each pair of origins (\code{o}) and destinations (\code{d}) on a 3D network, using the specified \code{weight} attribute (e.g., travel time). Origins and destinations must have their nearest network node IDs in a \code{nID} column (if \code{sf}) or be provided as node IDs (if \code{character}). The network can be an \code{igraph} object or a \code{list} with 3D \code{sf} nodes and edges. The output matrix reflects the \code{mode} of directionality.
#'
#' @return A \code{matrix} of least-cost values, with rows corresponding to origins (\code{o}) and columns corresponding to destinations (\code{d}), containing \code{numeric} cost values based on \code{weight}.
#'
#' @examples
#' \dontrun{
#'
#' #' # Create origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' d = st_sf(geometry = st_sfc(
#'   st_point(c(10, 10, 6)),
#'   st_point(c(5, 5, 5)),
#'   st_point(c(2, 2, 0))
#' ))
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 1,1,1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(1,1,1, 2,2,2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(2,2,2, 3,3,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(3,3,3, 4,4,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(4,4,4, 8,8,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(8,8,4, 12,13,5), ncol=3, byrow=TRUE))
#' ))
#' l = blendPT3d(rbind(o,d), l) # Blend points to the line network
#' l = splitNearestV3d(rbind(o,d), l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create network
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for o and d
#' o = nearestNode3d(net$nodes, o)
#' d = nearestNode3d(net$nodes, d)
#'
#' # Compute the OD matrix
#' ODM = ODM3d(net, o, d, mode = "out", "time")
#' print(ODM)
#' }
#' @export

ODM3d = function(net, o, d, mode = "out", weight){
  # --- Input validation ---
  if (!inherits(net, c("igraph", "list"))) {
    stop("net must be an igraph object or a list of edges and nodes.") }

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a character indicating the column of net's edges for weights.") }

  if (!mode %in% c("out", "in", "both")) { stop("mode must be one of 'out', 'in' and 'both'") }

  # Validate `o` and `d` are either `sf` features or character vectors
  if (inherits(o, "sf") && inherits(d, "sf")) {
    # If both `o` and `d` are `sf` features
    if (!"nID" %in% colnames(o) || !"nID" %in% colnames(d)) {
      stop("o and d must be an <sf> with POINT geometry and contain a column named 'nID' or a character vector of node indices.")
    }
    oID = seq_len(nrow(o))  # Assign IDs for the origin points
    dID = seq_len(nrow(d))  # Assign IDs for the destination points
    o = o$nID  # Reduce `o` to a vector of node IDs
    d = d$nID  # Reduce `d` to a vector of node IDs
    unique.d = unique(d)  # Retain unique node IDs for destinations
  } else if (is.character(o) && is.character(d)) {
    # If both `o` and `d` are character vectors
    oID = seq_len(length(o))  # Assign IDs for the origin points
    dID = seq_len(length(d))  # Assign IDs for the destination points
    unique.d = unique(d)  # Retain unique node IDs for destinations
  } else {
    stop("o and d must be an <sf> with POINT geometry and contain a column named 'nID' or a character vector of node indices.")
  }

  # --- Network preparation ---

  if (inherits(net, "list")) {
    # If `net` is a list, convert it into an `igraph` object
    if (!("edges" %in% names(net)) || !("nodes" %in% names(net))) {
      stop("Invalid input for net: When using a <list>, it must contain 'edges' and 'nodes'.") }

    net$edges = net$edges[c("FROM", "TO", "eID", weight)] %>% st_drop_geometry() # Drop geometry from edges
    weight = net$edges[[weight]] # Extract weight column
    net = igraph::graph_from_data_frame(d = net$edges, directed = T, vertices = net$nodes %>% sf::st_drop_geometry()) # Create igraph
  } else {
    # If `net` is already an `igraph` object
    weight = igraph::as_data_frame(net)[[weight]] # Extract weight column
  }

  # --- OD matrix calculation ---
  if(mode == "out"){
    # Calculate distances from origin to destination (outward)
    result = igraph::distances(net, v = o, to = unique.d, mode = "out", weights = weight)
  } else if(mode == "in"){
    # Calculate distances from destination to origin (inward)
    result = igraph::distances(net, v = o, to = unique.d, mode = "in", weights = weight)
  } else if(mode == "both"){
    # Calculate distances in both directions and sum them
    result = igraph::distances(net, v = o, to = unique.d, mode = "out", weights = weight)
    result = result + igraph::distances(net, v = o, to = unique.d, mode = "in", weights = weight)
  } else {stop("mode must be one of 'out', 'in', or 'both'.")}

  # --- Post-processing ---
  # Reintroduce destinations with shared node IDs as columns in the result matrix
  result = result[,match(d, unique.d), drop = F]
  colnames(result) = dID # Rename matrix columns as IDs of destinations
  rownames(result) = 1:nrow(result)
  return(result)  # <---- OUTPUT: OD matrix
}


##_____________________________________________________________________________
#' Compute Spatial Accessibility on a 3D Network
#'
#' Computes spatial accessibility indicators between 3D origins and destinations on a 3D network, including least-cost paths and destination availability within a cost threshold.
#'
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param o An \code{sf} object with 3D \code{POINT} geometry representing origins, optionally including a column \code{nID} for the nearest network node.
#' @param d An \code{sf} object with 3D \code{POINT} geometry representing destinations, optionally including a column \code{nID} for the nearest network node.
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination), \code{"in"} (destination to origin), or \code{"both"} (sum of both directions, with nearest destination based on total cost).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#' @param threshold A non-negative \code{numeric} value specifying the cost threshold for destination availability (default: 600).
#' @param avaiW A \code{numeric} vector of weights for destinations, matching the number of rows in \code{d}, or \code{NULL} (default: each destination weighted as 1).
#' @param buf A non-negative \code{numeric} value specifying the buffer distance (in map units) for tiling (default: 2000).
#' @param nrows A non-negative \code{integer} value specifying the number of rows for spatial tiling (default: 1).
#' @param ncols A non-negative \code{integer} value specifying the number of columns for spatial tiling (default: 2).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 2).
#' @param splitNET A \code{logical} value specifying whether to split the network by tiles (default: FALSE).
#'
#' @details
#' The function computes accessibility for each 3D origin in \code{o}, including the nearest destination and the weighted count of destinations within \code{threshold}. If \code{nID} is missing in \code{o} or \code{d}, the nearest node is assigned using \code{\link[GISnetwork3D]{nearestNode3d}}. Origins are divided into \code{nrows} x \code{ncols} spatial tiles, and destinations within \code{buf} distance of each tile are processed in parallel using \code{ncores} workers. If \code{nrows} and \code{ncols} are 0, no tiling is applied and the results are computed with single core. User can decide whether to also tile the network with the splitNET argument. Tiling network may speed up the process if the spatial extent and the network data are VERY LARGE. However, it can also slow down the process if the network and spatial extent are moderate to small.
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{\code{oID}}{index of the origin in \code{o}.}
#'   \item{\code{oNode}}{node ID (\code{nID}) of the origin.}
#'   \item{\code{dID}}{index of the nearest destination in \code{d}.}
#'   \item{\code{dNode}}{node ID (\code{nID}) of the nearest destination.}
#'   \item{\code{costs}}{cost to the nearest destination.}
#'   \item{\code{avai}}{weighted count of destinations within \code{threshold}, using \code{avaiW}.}
#' }
#'
#' @examples
#' \dontrun{
#'
#' #' # Create origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' d = st_sf(geometry = st_sfc(
#'   st_point(c(410, 550, 5.1)),
#'   st_point(c(181, 281, 4.1)),
#'   st_point(c(91, 98, 02.9))
#' ))
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 30,500,1.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(30,500,1.5, 60,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(60,100,3, 90,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(90,100,3, 100,120,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,120,4, 180,280,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(180,280,4, 400,560,5), ncol=3, byrow=TRUE))
#' ))
#' l = blendPT3d(rbind(o,d), l) # Blend points to the line network
#' l = splitNearestV3d(rbind(o,d), l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create network
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for o and d
#' o = nearestNode3d(net$nodes, o)
#' d = nearestNode3d(net$nodes, d)
#'
#' # Compute accessibility
#' result = access3d.pairOD(net, o, d, mode = "out", "time", threshold = 900, avaiW = c(1,1,5))
#' print(result)
#' }
#'
#' @export

access3d.pairOD = function(net, o, d, mode = "out", weight, threshold = 600, avaiW = NULL,
                           buf = 2000, nrows = 1, ncols = 2, ncores = 2, splitNET = FALSE){

  # --- Input validation ---
  if (!inherits(net, "list") || !all(c("nodes", "edges") %in% names(net))) {
    stop("net must be a list with 'nodes' and 'edges'.")}

  if (!inherits(o, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(o))) {
    stop("o must be an <sf> object with POINT geometry.") }

  if (!inherits(d, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(d))) {
    stop("d must be an <sf> object with POINT geometry.") }

  if (!mode %in% c("out", "in", "both")) {
    stop("mode must be one of 'out', 'in', or 'both'.")}

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a character indicating the column from net's edges for weights.")}

  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a non-negative numeric value.")}

  if (!is.null(avaiW)) {
    if (!is.numeric(avaiW) || length(avaiW) != nrow(d)) {
      stop("avaiW must be a numeric vector and match the number of destinations.") }
    } else { avaiW = rep(1, nrow(d)) }

  if (!is.numeric(buf) || length(buf) != 1 || buf < 0) {
    stop("buf must be a non-negative numeric value.")}

  if (!is.numeric(nrows) || length(nrows) != 1 || nrows < 0 || nrows %% 1 != 0) {
    stop("nrows must be a non-negative integer value.")}

  if (!is.numeric(ncols) || length(ncols) != 1 || ncols < 0 || ncols %% 1 != 0) {
    stop("ncols must be a non-negative integer value.")}

  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores %% 1 != 0) {
    stop("ncores must be a non-negative integer value.")}

  if (!is.logical(splitNET)) {
    stop("splitNET must be a <logical> value (TRUE/FALSE).")}

  # --- Assign nearest nodes if missing ---
  # If "nID" is missing from `o`, assign nearest network nodes
  if (!"nID" %in% colnames(o)) {o = GISnetwork3D::nearestNode3d(net$nodes, o)}

  # If "nID" is missing from `d`, assign nearest network nodes
  if (!"nID" %in% colnames(d)) {d = GISnetwork3D::nearestNode3d(net$nodes, d)}

  # --- Main function implementation ---
  if(nrows == 0 & ncols == 0){
    # Compute accessibility without tiling
    net = igraph::graph_from_data_frame(d = sf::st_drop_geometry(net$edges), directed = T, vertices = sf::st_drop_geometry(net$nodes))
    result = GISnetwork3D::access3d.pairOD.sTask(net, o$nID, d$nID, mode = mode, weight = weight, threshold = threshold, avaiW = avaiW)
    return(result) }  # <---- OUTPUT: data.frame

  if(nrows > 0 & ncols > 0){

    # Split Origin and Destination by Tiles
    input = GISnetwork3D::tile.OD(o = o, d = d, buf = buf, nrows = nrows, ncols = ncols, returnTiles = TRUE)
    for(i in 1:length(input$d)){ input$d[[i]]$avaiW = avaiW[input$d[[i]]$id] } # Add Weight for each destination

    ### --- Compute Accessibility without spiting network ---
    if(isFALSE(splitNET)){

      net = igraph::graph_from_data_frame(d = sf::st_drop_geometry(net$edges), directed = T, vertices = sf::st_drop_geometry(net$nodes))

      future::plan("future::multisession", workers = ncores)
      result = furrr::future_map2(input$o, input$d,
                                  function(x, y, net, mode, weight, threshold){
                                    result = GISnetwork3D::access3d.pairOD.sTask(net = net, o = x$nID, d = y$nID, mode = mode, weight = weight, threshold = threshold, avaiW = x$avaiW)
                                    result$oID = x$id[result$oID]
                                    result$dID = y$id[result$dID]
                                    return(result) },
                                  net = net, mode = mode, weight = weight, threshold = threshold,
                                  .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf")))
      future::plan("future::sequential") }

    ### --- Compute Accessibility with spiting network ---
    if(isTRUE(splitNET)){

      net = purrr::map(GISnetwork3D::tile.point(tiles = input$d.tiles, pt = net$nodes, buf = buf),
                       ~ subsetNet3d(net = net, v = .x$nID, exact = F, output = "graph"))

      future::plan("future::multisession", workers = ncores)
      result = furrr::future_pmap(list(net, input$o, input$d),
                                  function(net, x, y, mode, weight, threshold){
                                    result = GISnetwork3D::access3d.pairOD.sTask(net = net, o = x$nID, d = y$nID, mode = mode, weight = weight, threshold = threshold, avaiW = x$avaiW)
                                    result$oID = x$id[result$oID]
                                    result$dID = y$id[result$dID]
                                    return(result)},
                                  mode = mode, weight = weight, threshold = threshold,
                                  .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf")))
      future::plan("future::sequential") }

    result = purrr::reduce(result, rbind)
    result = data.frame(oID = 1:nrow(o)) %>% dplyr::left_join(result, by = "oID")

    return(result) }
  }

##_____________________________________________________________________________
#' Compute Spatial Accessibility and Edge Attribute Summaries on a 3D Network
#'
#' Computes accessibility indicators and edge attribute summaries along shortest paths between 3D origins and destinations on a 3D network, supporting parallel computation and optional geometry output.
#'
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param o An \code{sf} object with 3D \code{POINT} geometry representing origins, optionally including a column \code{nID} for the nearest network node.
#' @param d An \code{sf} object with 3D \code{POINT} geometry representing destinations, optionally including a column \code{nID} for the nearest network node.
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination), \code{"in"} (destination to origin), or \code{"both"} (sum of both directions, with nearest destination based on total cost).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#' @param threshold A non-negative \code{numeric} value specifying the cost threshold for destination availability (default: 600).
#' @param avaiW A \code{numeric} vector of weights for destinations, matching the number of rows in \code{d}, or \code{NULL} (default: each destination weighted as 1).
#' @param buf A non-negative \code{numeric} value specifying the buffer distance (in map units) for tiling (default: 2000).
#' @param nrows A non-negative \code{integer} value specifying the number of rows for spatial tiling (default: 1).
#' @param ncols A non-negative \code{integer} value specifying the number of columns for spatial tiling (default: 2).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 2).
#' @param geom A \code{logical} value indicating whether to include route geometry in the output (default: \code{TRUE}).
#' @param sumV A \code{character} vector of edge column names for simple sum along each route, or \code{NULL} (default).
#' @param w.sumV A \code{character} vector of edge column names for weighted sum, or \code{NULL} (default).
#' @param w.sumW A \code{character} value specifying the edge column for weighting \code{w.sumV}, or \code{NULL} (default).
#' @param p.sumV A \code{character} vector of edge column names for normalised sum, or \code{NULL} (default).
#' @param pw.sumV A \code{character} vector of edge column names for weighted-normalised sum, or \code{NULL} (default).
#' @param p.norm A \code{character} value specifying the edge column for normalisation of \code{p.sumV} or \code{pw.sumV}, or \code{NULL} (default).
#' @param splitNET A \code{logical} value specifying whether to split the network by tiles (default: FALSE).
#' @param both.split A \code{logical} value specifying whether to split summary into outward and return trips when mode = 'both' (default: FALSE).
#'
#' @details
#' The function computes accessibility indicators for each 3D origin in \code{o}, including the nearest destination and the weighted count of destinations within \code{threshold}, using \code{avaiW}. It also summarises edge attributes along each shortest path, supporting:
#' \itemize{
#'   \item Simple sum (\code{sumV}): Sum of each column’s values across route edges.
#'   \item Weighted sum (\code{w.sumV}, \code{w.sumW}): Sum of (\code{w.sumV} * \code{w.sumW}) across route edges.
#'   \item Normalised sum (\code{p.sumV}, \code{p.norm}): Sum of \code{p.sumV} divided by sum of \code{p.norm} across route edges.
#'   \item Weighted-normalised sum (\code{pw.sumV}, \code{p.norm}): Sum of (\code{pw.sumV} * \code{p.norm}) divided by sum of \code{p.norm} across route edges.
#' }
#' If \code{nID} is missing in \code{o} or \code{d}, the nearest node is assigned using \code{\link[GISnetwork3D]{nearestNode3d}}. Origins are divided into \code{nrows} x \code{ncols} spatial tiles, processing destinations within \code{buf} distance in parallel using \code{ncores} workers. If \code{nrows = 0} and \code{ncols = 0}, tiling is disabled. User can decide whether to also tile the network with the \code{splitNET} argument. Tiling network may speed up the process if the spatial extent and the network data are VERY LARGE. However, it can also slow down the process if the network and spatial extent are moderate to small.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry (if \code{geom = TRUE}) or a \code{data.frame} (if \code{geom = FALSE}) with columns:
#' \describe{
#'   \item{\code{oID}}{index of the origin in \code{o}.}
#'   \item{\code{oNode}}{node ID of the origin.}
#'   \item{\code{dID}}{index of the nearest destination in \code{d}.}
#'   \item{\code{dNode}}{node ID of the nearest destination.}
#'   \item{\code{costs}}{cost to the nearest destination.}
#'   \item{\code{avai}}{weighted count of destinations within \code{threshold}, using \code{avaiW}.}
#'   \item{\code{...}}{columns for each requested edge attribute summary.}
#'   \item{\code{geometry}}{3D route geometry (if \code{geom = TRUE}).}
#' }
#'
#' @examples
#' \dontrun{
#'
#' # Create sample origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' d = st_sf(geometry = st_sfc(
#'   st_point(c(410, 550, 5.1)),
#'   st_point(c(181, 281, 4.1)),
#'   st_point(c(91, 98, 02.9))
#' ))
#'
#' # Create a 3D line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 30,500,1.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(30,500,1.5, 60,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(60,100,3, 90,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(90,100,3, 100,120,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,120,4, 180,280,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(180,280,4, 400,560,5), ncol=3, byrow=TRUE))
#' ))
#' # Assign random edge attributes for summarisation
#' l$v1 = sample(1:5, nrow(l), replace = T)
#' l$v2 = sample(1:5, nrow(l), replace = T)
#' l$v3 = sample(1:5, nrow(l), replace = T)
#' l$v4 = sample(1:5, nrow(l), replace = T)
#' l$v5 = sample(1:5, nrow(l), replace = T)
#'
#' # Blend points to the network and split at nearest vertices
#' l = blendPT3d(rbind(o,d), l) # Blend points to the line network
#' l = splitNearestV3d(rbind(o,d), l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create 3D network with Tobler's hiking function for cost
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify nearest network nodes for origins and destinations
#' o = nearestNode3d(net$nodes, o)
#' d = nearestNode3d(net$nodes, d)
#'
#' # Compute accessibility and summarize edge attributes
#' result = path3d.pairOD.summary(net, o, d, mode = "out", threshold = 900, avaiW = c(1,1,5), weight = "time",
#'                                sumV = c("time", "energy"),
#'                                w.sumV = c("v1", "v2"), w.sumW = "time",
#'                                p.sumV = c("v1", "v2", "v3"), pw.sumV = c("v4", "v5"), p.norm = "energy")
#' print(result)
#' }
#'
#' @export

access3d.pairOD.summary = function(net, o, d, mode = "out", weight, threshold, avaiW = NULL,
                                   buf = 2000, nrows = 2, ncols = 2, ncores = 4, geom = TRUE,
                                   sumV = NULL,
                                   w.sumV = NULL, w.sumW = NULL,
                                   p.sumV = NULL, pw.sumV = NULL, p.norm = NULL, splitNET = FALSE, both.split = FALSE){

  # --- Input validation ---
  # net
  if (!is.list(net) || !all(c("nodes", "edges") %in% names(net))) {
    stop("net must be a list containing 'nodes' and 'edges'.") }

  # Origin
  if (!inherits(o, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(o))) { stop("o must be an <sf> object with POINT geometry.") }

  # Destination
  if (!inherits(d, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(d))) { stop("d must be an <sf> object with POINT geometry.") }

  # Scalar parameters
  if (!mode %in% c("out", "in", "both")) {stop("mode must be one of 'out', 'in', or 'both'.")}

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a single column name present in net$edges.")}

  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a non-negative numeric value.")}

  if (!is.null(avaiW)) {
    if (!is.numeric(avaiW) || length(avaiW) != nrow(d)) {
      stop("avaiW must be a numeric vector and match the number of destinations.") }
  } else { avaiW = rep(1, nrow(d)) }

  if (!is.numeric(buf) || length(buf) != 1 || buf < 0) {stop("buf must be a non-negative numeric value.")}

  if (!is.numeric(nrows) || length(nrows) != 1 || nrows < 0 || nrows %% 1 != 0) {stop("nrows must be a non-negative integer value.")}
  if (!is.numeric(ncols) || length(ncols) != 1 || ncols < 0 || ncols %% 1 != 0) {stop("ncols must be a non-negative integer value.")}
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores %% 1 != 0) {stop("ncores must be a non-negative integer value.")}
  if (!is.logical(splitNET)) {stop("splitNET must be a <logical> value (TRUE/FALSE).")}

  # Optional named parameters (character checks)
  if (!is.null(sumV) && !is.character(sumV)) { stop("sumV must be NULL or character vector.") }
  if (!is.null(w.sumV) && !is.character(sumV)) { stop("w.sumV must be NULL or a character vector.") }
  if (!is.null(p.sumV) && !is.character(p.sumV)) { stop("p.sumV must be NULL or a character vector.") }
  if (!is.null(pw.sumV) && !is.character(pw.sumV)) { stop("pw.sumV must be NULL or a character vector.") }

  if (!is.null(w.sumW) && (!is.character(w.sumW) || length(w.sumW) > 1)) {stop("w.sumW must be NULL or a single character string.")}
  if (!is.null(p.norm) && (!is.character(p.norm) || length(p.norm) > 1)) {stop("p.norm must be NULL or a single character string.")}

  # Dependency rules
  if (!is.null(p.sumV)  && is.null(p.norm)) stop("p.norm is required when p.sumV is provided.")
  if (!is.null(pw.sumV) && is.null(p.norm)) stop("p.norm is required when pw.sumV is provided.")
  if (!is.null(p.norm)  && is.null(p.sumV) && is.null(pw.sumV)) {stop("At least one of p.sumV or pw.sumV is required when p.norm is provided.")}
  if (!is.null(w.sumV)  != !is.null(w.sumW)) {stop("w.sumV and w.sumW must both be provided or both be NULL.")}


  # --- Assign nearest nodes if missing ---
  # If "nID" is missing from `o`, assign nearest network nodes
  if (!"nID" %in% colnames(o)) {o = nearestNode3d(net$nodes, o)}

  # If "nID" is missing from `d`, assign nearest network nodes
  if (!"nID" %in% colnames(d)) {d = nearestNode3d(net$nodes, d)}

  # --- Main function implementation ---
  # Compute without tiling
  if(nrows == 0 & ncols == 0){
    result = GISnetwork3D::access3d.pairOD.summary.sTask(net = net, o = o, d = d, mode = mode, weight = weight, threshold = threshold, avaiW = avaiW,
                                                         geom = geom,
                                                         sumV = sumV,
                                                         w.sumV = w.sumV, w.sumW = w.sumW,
                                                         p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm, both.split = both.split) }

  # Compute with O and D being tiled
  if(nrows >= 1 & ncols >= 1 & isFALSE(splitNET)){
    # Split Origin and Destination by Tiles
    input = GISnetwork3D::tile.OD(o = o, d = d, buf = buf, nrows = nrows, ncols = ncols, returnTiles = TRUE)
    for(i in 1:length(input$d)){ input$d[[i]]$avaiW = avaiW[input$d[[i]]$id] } # Add Weight for each destination

    future::plan("future::multisession", workers = ncores)
    result = furrr::future_map2(input$o, input$d,
                                function(x, y, net, mode, weight, threshold, geom = geom, sumV, w.sumV, w.sumW, p.sumV, pw.sumV, p.norm, both.split){
                                  result = GISnetwork3D::access3d.pairOD.summary.sTask(net = net, o = x$nID, d = y$nID,
                                                                                       mode = mode, weight = weight, threshold = threshold, avaiW = y$avaiW,
                                                                                       geom = geom,
                                                                                       sumV = sumV,
                                                                                       w.sumV = w.sumV, w.sumW = w.sumW,
                                                                                       p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm, both.split = both.split)
                                  result$oID = x$id[result$oID]
                                  result$dID = y$id[result$dID]
                                  return(result) },
                                net = net, mode = mode, weight = weight, threshold = threshold, geom = geom, sumV = sumV,
                                w.sumV = w.sumV, w.sumW = w.sumW, p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm, both.split = both.split,
                                .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf"))) %>% purrr::reduce(rbind)
    future::plan("future::sequential")

    result = result[order(result$oID),] }


  # Compute with O, D and net being tiled
  if(nrows >= 1 & ncols >= 1 & isTRUE(splitNET)){

    # Split Origin, Destination and network by Tiles
    input = GISnetwork3D::tile.OD(o = o, d = d, buf = buf, nrows = nrows, ncols = ncols, returnTiles = TRUE)
    for(i in 1:length(input$d)){ input$d[[i]]$avaiW = avaiW[input$d[[i]]$id] } # Add Weight for each destination
    net = purrr::map(GISnetwork3D::tile.point(tiles = input$d.tiles, pt = net$nodes, buf = buf),
                     ~ subsetNet3d(net = net, v = .x$nID, exact = F, output = "sf"))


    future::plan("future::multisession", workers = ncores)
    result = furrr::future_pmap(list(net, input$o, input$d),
                                function(net, x, y, mode, weight, threshold, geom = geom, sumV, w.sumV, w.sumW, p.sumV, pw.sumV, p.norm, both.split){
                                  result = GISnetwork3D::access3d.pairOD.summary.sTask(net = net, o = x$nID, d = y$nID,
                                                                                       mode = mode, weight = weight, threshold = threshold, avaiW = y$avaiW,
                                                                                       geom = geom,
                                                                                       sumV = sumV,
                                                                                       w.sumV = w.sumV, w.sumW = w.sumW,
                                                                                       p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm, both.split = both.split)
                                  result$oID = x$id[result$oID]
                                  result$dID = y$id[result$dID]
                                  return(result) },
                                mode = mode, weight = weight, threshold = threshold, geom = geom, sumV = sumV,
                                w.sumV = w.sumV, w.sumW = w.sumW, p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm, both.split = both.split,
                                .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf"))) %>% purrr::reduce(rbind)
    future::plan("future::sequential")

    result = result[order(result$oID),] }

  rownames(result) = NULL
  return(result) # <--- Output: sf or data.frame
  }


##_____________________________________________________________________________
#' Extract Edge Indices for Shortest Paths in a 3D Network
#'
#' Computes the shortest path between paired origin and destination nodes on a 3D network, returning the edge indices (\code{eID}) for each path.
#'
#' @param net An \code{igraph} object or a \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param oNode A \code{character} vector of node IDs representing origins in the 3D network.
#' @param dNode A \code{character} vector of node IDs representing destinations in the 3D network, matching the length of \code{oNode}.
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination) or \code{"in"} (destination to origin).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net} for travel cost (e.g., \code{"time"}).
#' @param output A \code{character} value specifying the output format: \code{"list"} (list of edge indices) or \code{"df"} (data frame with \code{oID} and \code{eID}).
#' @param oID An \code{integer}, \code{numeric}, or \code{character} vector of unique origin identifiers, matching the length of \code{oNode}, or \code{NULL} (default: \code{1:length(oNode)}).
#'
#' @details
#' For each pair of origin (\code{oNode}) and destination (\code{dNode}) nodes in a 3D network, the function computes the shortest path using the \code{weight} attribute and returns the sequence of edge indices (\code{eID}). Node IDs must correspond to nodes in the 3D network. The output is either a \code{list} (each element named by \code{oID}, containing edge indices) or a \code{data.frame} (with \code{oID} and \code{eID} columns).
#'
#' @return Depends on \code{output}:
#' \describe{
#'   \item{\code{list}}{A \code{list} where each element, named by \code{oID}, contains a vector of edge indices (\code{eID}) for the shortest path.}
#'   \item{\code{df}}{A \code{data.frame} with columns: \code{oID} (\code{integer}, \code{numeric}, or \code{character}, origin identifier) and \code{eID} (\code{numeric}, edge index).}
#' }
#'
#' @examples
#' \dontrun{
#'
#' #' # Create origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' d = st_sf(geometry = st_sfc(
#'   st_point(c(410, 550, 5.1)),
#'   st_point(c(181, 281, 4.1)),
#'   st_point(c(91, 98, 02.9))
#' ))
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 30,500,1.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(30,500,1.5, 60,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(60,100,3, 90,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(90,100,3, 100,120,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,120,4, 180,280,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(180,280,4, 400,560,5), ncol=3, byrow=TRUE))
#' ))
#' l = blendPT3d(rbind(o,d), l) # Blend points to the line network
#' l = splitNearestV3d(rbind(o,d), l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create 3D network with Tobler's hiking function for cost
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for o and d
#' o = nearestNode3d(net$nodes, o)
#' d = nearestNode3d(net$nodes, d)
#'
#' # Compute accessibility
#' result = access3d.pairOD(net, o, d, mode = "out", "time", threshold = 900, avaiW = c(1,1,5))
#'
#' # Get the path between each pair of origin and destination
#' result = path3d.pairOD.get_eID(net, oNode = result$oNode, dNode = result$dNode, mode = "out", weight = "time", output = "list")
#' print(result)
#' }
#'
#' @export

path3d.pairOD.get_eID = function(net, oNode, dNode, mode = "out", weight, output = "list", oID = NULL){
  # --- Input validation ---
  if (!inherits(net, c("igraph", "list"))) {stop("net must be an igraph a list with 'nodes' and 'edges'.")}

  if (!mode %in% c("out", "in")) {stop("mode must be one of 'out' and 'in'.")}

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a character indicating the column from net's edges for weights.")}

  if (!is.character(oNode)) {stop("oNode must be a character vector specifying node indices.")}

  if (!is.character(dNode)) {stop("dNode must be a character vector specifying node indices.")}

  if (length(oNode) != length(dNode)) {stop("oNode and dNode must have the same length.")}

  if (!is.character(output) || length(output) != 1){
    stop("output must be a character of 'list' or 'df'.")}

  if (!output %in% c("list", "df")) {
    stop("output must be a character of 'list' or 'df'.")}

  if (!is.null(oID)) {
    if (length(oID) != length(oNode)) {
      stop("oID must have the same length as oNode.")}
    if (!is.numeric(oID) && !is.character(oID) && !is.integer(oID)) {
      stop("oID must be a vector of either integer, numeric, or character.")
    }
  }

  # --- Network preparation ---
  if (inherits(net, "list")) {
    # If `net` is a list, extract the weight column and convert to an igraph object
    weight = net$edges[[weight]] # Extract the weight
    # Convert to igraph
    net = igraph::graph_from_data_frame(d = net$edges, directed = T, vertices = net$nodes %>% st_drop_geometry())
  } else {
    # If `net` is already an igraph object, extract the weight column
    weight = igraph::as_data_frame(net)[[weight]] # Extract the weight
  }

  # Define `oID` if not provided
  if (is.null(oID)) { oID = 1:length(oNode)} # Generate default IDs for origins

  # --- Shortest path computation ---

  # Initialise a list to store the edges for each shortest path
  edges = list()

  # Loop through each pair of origin and destination nodes
  if(mode == "out"){
    for(i in 1:length(oNode)){
      edges[[i]] = igraph::shortest_paths(net, from = oNode[i], to = dNode[i], mode = "out", weights = weight,
                                          output = "epath")$epath[[1]]$eID } }

  if(mode == "in"){
    for(i in 1:length(oNode)){
      edges[[i]] = igraph::shortest_paths(net, from = dNode[i], to = oNode[i], mode = "out", weights = weight,
                                          output = "epath")$epath[[1]]$eID } }

  # Assign names to the edges list based on `oID`
  names(edges) = oID

  # --- Output processing ---

  if(output == "list"){
    # Return edge sequence objects
    return(edges) # <---- OUTPUT: list
    }  else if(output == "df"){
      # Extract edge indices (`eID`) from the edge sequence objects
      edges = data.frame(oID = rep(oID, base::unlist(purrr::map(edges, length), use.names = F)),
                         eID = base::unlist(edges, use.names = F))
      return(edges)}  # <---- OUTPUT: data.frame

}

##_____________________________________________________________________________
#' Extract Network Edges for Shortest Paths in a 3D Network
#'
#' Computes shortest paths between paired origin and destination nodes in a 3D network, returning the sequence of network edges as \code{sf} \code{LINESTRING} objects or \code{data.frame}s.
#'
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param oNode A \code{character} vector of node IDs representing origins in the 3D network.
#' @param dNode A \code{character} vector of node IDs representing destinations in the 3D network, matching the length of \code{oNode}.
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination), \code{"in"} (destination to origin), or \code{"both"} (origin to destination and back).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net} for travel cost (e.g., \code{"time"}).
#' @param oID An \code{integer}, \code{numeric}, or \code{character} vector of unique origin identifiers, matching the length of \code{oNode}, or \code{NULL} (default: \code{1:length(oNode)}).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 4).
#' @param tasks A positive \code{integer} value specifying the number of task chunks for parallel computation (default: 4).
#' @param output A \code{character} value specifying the output format: \code{"eID"} (edge IDs), \code{"df"} (edge attributes), \code{"raw"} (edge geometries), or \code{"dissolve"} (merged geometries by \code{oID}) (default: \code{"raw"}).
#'
#' @details
#' For each pair of origin (\code{oNode}) and destination (\code{dNode}) nodes in a 3D network, the function computes the shortest path using the \code{weight} attribute and extracts the traversed edges. Outputs include edge indices (\code{eID}), edge attributes, or 3D \code{LINESTRING} geometries, with all edge attributes preserved for \code{"df"} and \code{"raw"} outputs. The \code{"dissolve"} output merges edges by \code{oID} into single linestrings. Parallel computation splits the workload into \code{tasks} chunks across \code{ncores} workers. The function preserves the coordinate reference system (CRS) and Z coordinates of the input network and relies on \code{path3d.pairOD.get_eID}, \code{igraph}, \code{furrr}, \code{sf}, and \code{sfheaders} for path and geometry processing.
#'
#' @return Depends on \code{output}:
#' \describe{
#'   \item{\code{eID}}{A \code{data.frame} with columns: \code{oID} (\code{integer}, \code{numeric}, or \code{character}, origin identifier) and \code{eID} (\code{numeric}, edge index).}
#'   \item{\code{df}}{A \code{data.frame} with attributes of traversed edges and an \code{oID} column.}
#'   \item{\code{raw}}{An \code{sf} object with 3D \code{LINESTRING} geometries, edge attributes, and an \code{oID} column.}
#'   \item{\code{dissolve}}{An \code{sf} object with 3D \code{LINESTRING} geometries, merged by \code{oID}.}
#' }
#'
#' @examples
#' \dontrun{
#'
#' #' # Create origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' d = st_sf(geometry = st_sfc(
#'   st_point(c(410, 550, 5.1)),
#'   st_point(c(181, 281, 4.1)),
#'   st_point(c(91, 98, 02.9))
#' ))
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 30,500,1.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(30,500,1.5, 60,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(60,100,3, 90,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(90,100,3, 100,120,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,120,4, 180,280,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(180,280,4, 400,560,5), ncol=3, byrow=TRUE))
#' ))
#' # Prepare network
#' l = blendPT3d(rbind(o,d), l) # Blend points to the line network
#' l = splitNearestV3d(rbind(o,d), l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create 3D network with Tobler's hiking function for cost
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for o and d
#' o = nearestNode3d(net$nodes, o)
#' d = nearestNode3d(net$nodes, d)
#'
#' # Compute accessibility
#' result = access3d.pairOD(net, o, d, mode = "out", "time", threshold = 900, avaiW = c(1,1,5))
#'
#' # Get the path between each pair of origin and destination
#' result = path3d.pairOD.getEdges(net, oNode = result$oNode, dNode = result$dNode, mode = "out", weight = "time", output = "raw")
#' print(result)
#' }
#' @export

path3d.pairOD.getEdges = function(net, oNode, dNode, mode = "out", weight, oID = NULL, ncores = 4, tasks = 4, output = "raw"){

  # --- Input validation ---
  if (!inherits(net, c("igraph", "list"))) {stop("net must be a list with 'nodes' and 'edges'.")}

  if (!is.character(oNode)) {stop("oNode must be a character vector specifying node indices.")}

  if (!is.character(dNode)) {stop("dNode must be a character vector specifying node indices.")}

  if (length(oNode) != length(dNode)) {stop("oNode and dNode must have the same length.")}

  if (!mode %in% c("out", "in", "both")) {stop("mode must be one of 'out', 'in' and 'both'.")}

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a character indicating the column from net's edges for weights.")}

  if (!is.null(oID)) {
    if (length(oID) != length(oNode)) {
      stop("oID must have the same length as oNode.")}
    if (!is.numeric(oID) && !is.character(oID) && !is.integer(oID)) {
      stop("oID must be a vector of either integer, numeric, or character.")
    }
  }

  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores != round(ncores)) {
    stop("ncores must be a positive integer.")}

  if (!is.numeric(tasks) || length(tasks) != 1 || tasks <= 0 || tasks != round(tasks)) {
    stop("tasks must be a positive integer.")}

  if(is.null(oID)){oID = 1:length(oNode)}

  if(!is.character(output) || length(output) != 1){
    stop("output must be a character of either 'eID', 'df', 'raw' or 'dissolve'.")
  }
  if (!output %in% c("eID", "df", "raw", "dissolve")) {
    stop("output must be a character of either 'eID', 'df', 'raw' or 'dissolve'.")}

  # --- Prepare data for multithreading---
  # Split input vectors into chunks for parallel processing
  oNode = GISnetwork3D::split.data(oNode, tasks)
  dNode = GISnetwork3D::split.data(dNode, tasks)
  oID = GISnetwork3D::split.data(oID, tasks)

  # --- Extract network edges using multiple cores ---
  if(mode != "both"){

    future::plan("future::multisession", workers = ncores) # Set up parallel plan

    # Map over the split data, executing path computation for each chunk
    result = furrr::future_pmap( list(oNode, dNode, oID),
                                 function(o, d, id, n, m, w){
                                   result = GISnetwork3D::path3d.pairOD.get_eID(net = n, oNode = o, dNode = d, mode = m,
                                                                                weight = w, output = "df", oID = id)
                                   return(result)
                                 }, n = net, m = mode, w = weight,
                                 .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "sf", "GISnetwork3D")) )
    future::plan("future::sequential") # Stop parallel cluster
  }

  if(mode == "both"){

    future::plan("future::multisession", workers = ncores) # Set up parallel plan

    # Map over the split data, executing path computation for each chunk
    result = furrr::future_pmap( list(oNode, dNode, oID),
                                 function(o, d, id, n, m, w){
                                   result.out = GISnetwork3D::path3d.pairOD.get_eID(net = n, oNode = o, dNode = d, mode = "out",
                                                                                   weight = w, output = "list", oID = id)
                                   result.in = GISnetwork3D::path3d.pairOD.get_eID(net = n, oNode = o, dNode = d, mode = "in",
                                                                                   weight = w, output = "list", oID = id)
                                   result = purrr::map2(result.out, result.in, base::c)

                                   result = data.frame(oID = rep(names(result), base::unlist(purrr::map(result, length), use.names = F)),
                                                       eID = base::unlist(result, use.names = F))
                                   class(result$oID) = class(id)
                                   return(result)
                                 }, n = net, w = weight,
                                 .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "sf", "GISnetwork3D")) )
    future::plan("future::sequential") # Stop parallel cluster
  }

  # --- Prepare output based on the specified format ---
  if(output == "eID"){
    # Combine all results into one data.frame with oID and eID columns
    return(purrr::reduce(result, rbind)) }  # <---- OUTPUT: data.frame

  # If output = "df", export the traveled edges with all their original fields with an additional field called "oID" in data.frame format
  if(output == "df"){
    # Combine results and match with network edges to get full attribute data
    result = purrr::reduce(result, rbind)
    # Match eID with network edges to get full attributes
    result = net$edges[ base::match(result$eID, net$edges$eID), ] %>% cbind( data.frame(oID = result$oID) ) # Append oID to the data
    # Return as data.frame with attributes
    return(result %>% sf::st_drop_geometry()) }  # <---- OUTPUT: data.frame

  if(output == "raw"){
    result = purrr::reduce(result, rbind)
    # Match eID with network edges to get full attributes
    result = net$edges[ base::match(result$eID, net$edges$eID), ] %>% cbind( data.frame(oID = result$oID) ) # Append oID
    # Return as sf LINESTRING with attributes
    return(result) }  # <---- OUTPUT: sf LINESTRING object

  if(output == "dissolve"){
    # Combine results and match with network edges
    result = purrr::reduce(result, rbind)
    # Keep only 'eID' for dissolving
    net$edges = net$edges %>% dplyr::select("eID")
    # Append oID
    result = net$edges[ base::match(result$eID, net$edges$eID), ] %>% cbind( data.frame(oID = result$oID) ) %>% dplyr::select(oID)
    # Get the coordinates for each linestring
    result = sfheaders::sf_to_df(result, fill = T)
    # Group geometries by oID
    result = sfheaders::sf_linestring(result, x = "x", y = "y", z = "z", linestring_id = "oID", keep = T)
    sf::st_crs(result) = sf::st_crs(net$edges) # Define CRS
    # Return dissolved linestrings with oID
    return(result %>% dplyr::select(oID)) }  # <---- OUTPUT: sf LINESTRING object

}


##_____________________________________________________________________________
#' Summarise Edge Attributes Along Shortest Paths in a 3D Network
#'
#' Computes and summarises edge attribute values for shortest paths between paired origin and destination nodes in a 3D network, supporting simple, weighted, normalised, and weighted-normalised summaries.
#'
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param oNode A \code{character} vector of node IDs representing origins in the 3D network.
#' @param dNode A \code{character} vector of node IDs representing destinations in the 3D network, matching the length of \code{oNode}.
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination), \code{"in"} (destination to origin), or \code{"both"} (origin to destination and back).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#' @param oID An \code{integer}, \code{numeric}, or \code{character} vector of unique origin identifiers, matching the length of \code{oNode}, or \code{NULL} (default: \code{1:length(oNode)}).
#' @param geom A \code{logical} value indicating whether to include route geometry in the output (default: \code{TRUE}).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 4).
#' @param tasks A positive \code{integer} value specifying the number of task chunks for parallel computation (default: 4).
#' @param sumV A \code{character} vector of edge column names for simple sum along each route, or \code{NULL} (default).
#' @param w.sumV A \code{character} vector of edge column names for weighted sum, or \code{NULL} (default).
#' @param w.sumW A \code{character} value specifying the edge column for weighting \code{w.sumV}, or \code{NULL} (default).
#' @param p.sumV A \code{character} vector of edge column names for normalised sum, or \code{NULL} (default).
#' @param pw.sumV A \code{character} vector of edge column names for weighted-normalised sum, or \code{NULL} (default).
#' @param p.norm A \code{character} value specifying the edge column for normalisation of \code{p.sumV} or \code{pw.sumV}, or \code{NULL} (default).
#'
#' @details
#' The function computes shortest paths in a 3D network using the \code{weight} attribute and summarises edge attributes for each origin-destination pair, supporting:
#' \itemize{
#'   \item Simple sum (\code{sumV}): Sum of each column’s values across route edges.
#'   \item Weighted sum (\code{w.sumV}, \code{w.sumW}): Sum of (\code{w.sumV} * \code{w.sumW}) across route edges.
#'   \item Normalised sum (\code{p.sumV}, \code{p.norm}): Sum of \code{p.sumV} divided by sum of \code{p.norm} across route edges.
#'   \item Weighted-normalised sum (\code{pw.sumV}, \code{p.norm}): Sum of (\code{pw.sumV} * \code{p.norm}) divided by sum of \code{p.norm} across route edges.
#' }
#' Parallel computation splits the workload into \code{tasks} chunks across \code{ncores} workers. The output includes one row per route, with columns for \code{oID} and requested summaries. If \code{geom = TRUE}, 3D \code{LINESTRING} geometries are included.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry (if \code{geom = TRUE}) or a \code{data.frame} (if \code{geom = FALSE}) with columns:
#' \describe{
#'   \item{\code{oID}}{\code{integer}, \code{numeric}, or \code{character}, origin identifier.}
#'   \item{\code{...}}{\code{numeric}, columns for each requested edge attribute summary (e.g., \code{sumV}, \code{w.sumV}, \code{p.sumV}, \code{pw.sumV}).}
#'   \item{\code{geometry}}{3D \code{LINESTRING}, route geometry (if \code{geom = TRUE}).}
#' }
#'
#' @examples
#' \dontrun{
#'
#' #' # Create origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' d = st_sf(geometry = st_sfc(
#'   st_point(c(410, 550, 5.1)),
#'   st_point(c(181, 281, 4.1)),
#'   st_point(c(91, 98, 02.9))
#' ))
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 30,500,1.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(30,500,1.5, 60,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(60,100,3, 90,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(90,100,3, 100,120,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,120,4, 180,280,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(180,280,4, 400,560,5), ncol=3, byrow=TRUE))
#' ))
#' # Assign random variables for summary
#' l$v1 = sample(1:5, nrow(l), replace = T)
#' l$v2 = sample(1:5, nrow(l), replace = T)
#' l$v3 = sample(1:5, nrow(l), replace = T)
#' l$v4 = sample(1:5, nrow(l), replace = T)
#' l$v5 = sample(1:5, nrow(l), replace = T)
#'
#' l = blendPT3d(rbind(o,d), l) # Blend points to the line network
#' l = splitNearestV3d(rbind(o,d), l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create 3D network with Tobler's hiking function for cost
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for o and d
#' o = nearestNode3d(net$nodes, o)
#' d = nearestNode3d(net$nodes, d)
#'
#' # Compute accessibility
#' result = access3d.pairOD(net, o, d, mode = "out", "time", threshold = 900, avaiW = c(1,1,5))
#'
#' # Get the path between each pair of origin and destination
#' result = path3d.pairOD.summary(net, oNode = result$oNode, dNode = result$dNode, mode = "out", weight = "time",
#'                                sumV = c("time", "energy"),
#'                                w.sumV = c("v1", "v2"), w.sumW = "time",
#'                                p.sumV = c("v1", "v2", "v3"), pw.sumV = c("v4", "v5"), p.norm = "energy")
#' print(result)
#' }
#'
#' @export

path3d.pairOD.summary = function(net, oNode, dNode, mode = "out", weight, oID = NULL, geom = TRUE, ncores = 4, tasks = 4,
                                 sumV = NULL,
                                 w.sumV = NULL, w.sumW = NULL,
                                 p.sumV = NULL, pw.sumV = NULL, p.norm = NULL){

  # --- Input validation ---
  # net
  if (!is.list(net) || !all(c("nodes", "edges") %in% names(net))) { stop("net must be a list containing 'nodes' and 'edges'.") }

  # Origins and Destinations
  if (!is.character(oNode)) {stop("oNode must be a character vector specifying node indices.")}
  if (!is.character(dNode)) {stop("dNode must be a character vector specifying node indices.")}
  if (length(oNode) != length(dNode)) {stop("oNode and dNode must have the same length.")}

  # Scalar parameters
  if (!mode %in% c("out", "in", "both")) {stop("mode must be one of 'out', 'in', or 'both'.")}

  if (!is.character(weight) || length(weight) != 1) {stop("weight must be a character indicating the column name in net's edges for weight.")}

  if (!is.character(weight) || length(weight) != 1 || !(weight %in% names(net$edges))) {stop("weight must be a single column name present in net$edges.") }

  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores != round(ncores)) {stop("ncores must be a positive integer.")}

  # Validate tasks: must be a positive integer
  if (!is.numeric(tasks) || length(tasks) != 1 || tasks <= 0 || tasks != round(tasks)) {
    stop("tasks must be a positive integer.")}

  if (!is.null(oID)) {
    if (length(oID) != length(oNode)) {
      stop("oID must have the same length as `oNode`.")}
    if (!is.numeric(oID) && !is.character(oID) && !is.integer(oID)) {
      stop("oID must be either integer, numeric, or character.")}
  } else {
      oID = 1:length(oNode)
    }

  if(length(geom) != 1 & !is.logical(geom)){stop("Invalid input for geom: Only accept <TRUE> OR <FALSE>")}

  # Optional named parameters (character checks)
  if (!is.null(sumV) && !is.character(sumV)) { stop("sumV must be NULL or character vector.") }
  if (!is.null(w.sumV) && !is.character(sumV)) { stop("w.sumV must be NULL or a character vector.") }
  if (!is.null(p.sumV) && !is.character(p.sumV)) { stop("p.sumV must be NULL or a character vector.") }
  if (!is.null(pw.sumV) && !is.character(pw.sumV)) { stop("pw.sumV must be NULL or a character vector.") }

  if (!is.null(w.sumW) && (!is.character(w.sumW) || length(w.sumW) > 1)) {stop("w.sumW must be NULL or a single character string.")}
  if (!is.null(p.norm) && (!is.character(p.norm) || length(p.norm) > 1)) {stop("p.norm must be NULL or a single character string.")}

  # Dependency rules
  if (!is.null(p.sumV)  && is.null(p.norm)) stop("p.norm is required when p.sumV is provided.")
  if (!is.null(pw.sumV) && is.null(p.norm)) stop("p.norm is required when pw.sumV is provided.")
  if (!is.null(p.norm)  && is.null(p.sumV) && is.null(pw.sumV)) {stop("At least one of p.sumV or pw.sumV is required when p.norm is provided.")}
  if (!is.null(w.sumV)  != !is.null(w.sumW)) {stop("w.sumV and w.sumW must both be provided or both be NULL.")}


  # --- Compute shortest paths and extract edge attributes for each path ---
  if(isTRUE(geom)){
    vars = GISnetwork3D::path3d.pairOD.getEdges(net = net, oNode = oNode, dNode = dNode, mode = mode,
                                                weight = weight, oID = oID, ncores = ncores, tasks = tasks, output = "raw")
    output = data.frame(oID = oID)
    output.geom = vars %>% sfheaders::sf_to_df(fill = T) %>% dplyr::select(oID, x, y, z) # Convert sf to df
    vars = vars %>% sf::st_drop_geometry()
  }
  if(isFALSE(geom)){
    vars = GISnetwork3D::path3d.pairOD.getEdges(net = net, oNode = oNode, dNode = dNode, mode = mode,
                                                weight = weight, oID = oID, ncores = ncores, tasks = tasks, output = "df")
    output = data.frame(oID = oID)}

  # Only keep relevant columns for summarisation
  vars = vars[c("oID", sumV, w.sumV, w.sumW, p.sumV, pw.sumV, p.norm) %>% unique()]

  # --- Simple Summation ---
  if(!is.null(sumV)){
    sumV.df = vars[c("oID",sumV)] # Extract variables
    sumV.df = stats::aggregate(. ~ oID, FUN = sum, data = sumV.df) # Sum variables by oID
    output = dplyr::left_join(output, sumV.df, by = "oID") # Join to the output
  }

  # --- Weighted sums ---
  if(!is.null(w.sumV)){
    w.sumV.df = vars[c("oID",w.sumV)]  # Extract variables
    w.sumV.df[w.sumV] = w.sumV.df[w.sumV] * eval(parse(text = paste0("vars$", w.sumW)))  # Multiply variables by weight
    w.sumV.df = stats::aggregate(. ~ oID, FUN = sum, data = w.sumV.df) # Sum variables by oID
    names(w.sumV.df) = c("oID", paste("w.", w.sumV, sep = "")) # Add suffix for the variables
    output = output %>% dplyr::left_join(w.sumV.df, by = "oID") # Join to the output
  }

  # --- Normalised sums ---

  if(!is.null(p.sumV) & !is.null(p.norm)){
    p.sumV.df = vars[c("oID",p.sumV,p.norm)]  # Extract variables
    p.sumV.df = stats::aggregate(. ~ oID, FUN = sum, data = p.sumV.df) # Sum variables by oID
    p.sumV.df[p.sumV] = p.sumV.df[p.sumV] / eval(parse(text = paste0("p.sumV.df$", p.norm)))  # Normalise variables
    p.sumV.df = p.sumV.df[c("oID",p.sumV)]
    names(p.sumV.df) = c("oID", paste("p.", p.sumV, sep = ""))
    output = output %>% dplyr::left_join(p.sumV.df, by = "oID")
  }

  # --- Weighted normalised sums ---
  if(!is.null(pw.sumV) & !is.null(p.norm)){
    pw.sumV.df = vars[c("oID",pw.sumV,p.norm)]  # Extract variables
    pw.sumV.df[pw.sumV] = pw.sumV.df[pw.sumV] * eval(parse(text = paste0("pw.sumV.df$", p.norm)))  # Weight variables
    pw.sumV.df = stats::aggregate(. ~ oID, FUN = sum, data = pw.sumV.df) # Sum variables by oID
    pw.sumV.df[pw.sumV] = pw.sumV.df[pw.sumV] / eval(parse(text = paste0("pw.sumV.df$", p.norm)))  # Normalise variables
    pw.sumV.df = pw.sumV.df[c("oID",pw.sumV)]
    names(pw.sumV.df) = c("oID", paste("pW.", pw.sumV, sep = ""))
    output = output %>% dplyr::left_join(pw.sumV.df, by = "oID")
  }

  # Replace NA value by 0
  output[is.na(output)] = 0

  # Join geometry
  if(isTRUE(geom)){
    output = output %>% dplyr::left_join(output.geom, by = "oID")
    output = output %>% sfheaders::sf_linestring(x = "x", y = "y", z = "z", linestring_id = "oID", keep = T) # Group geometry by oID
    sf::st_crs(output) = sf::st_crs(net$edges)
  }

  return(output)  # <---- OUTPUT
}

##_____________________________________________________________________________
#' Optimal Location Selection for Multiple Tasks in a 3D Network
#'
#' Solves the Generalised Traveling Salesperson Problem (GTSP) in a 3D network using a Nearest Neighbour heuristic to select one location per task group, minimising total travel cost.
#'
#' @param startPT An \code{sf} object with 3D \code{POINT} geometry representing starting points, with a \code{character} or \code{numeric} column \code{ID} (unique identifier) and optional \code{character} column \code{nID} (nearest network node).
#' @param taskPT A \code{list} of \code{sf} objects with 3D \code{POINT} geometry, each representing a task group (e.g., grocery stores), with a \code{character} or \code{numeric} column \code{ID} and optional \code{character} column \code{nID}.
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param edges.output A \code{character} value specifying edge output format: \code{NULL} (no edges), \code{"dissolved"} (dissolved route geometry), \code{"raw"} (raw edge geometries), or \code{"raw.df"} (raw edge attributes without geometry) (default: \code{NULL}).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net} for travel cost (e.g., \code{"time"}).
#' @param buf A non-negative \code{numeric} value specifying the initial buffer distance (in map units) for spatial tiling (default: 2000).
#' @param nrows A non-negative \code{integer} value specifying the number of rows for spatial tiling (default: 4).
#' @param ncols A non-negative \code{integer} value specifying the number of columns for spatial tiling (default: 4).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 2).
#' @param coef A \code{numeric} value to multiply the buffer distance per iteration of spatial tiling (default: 1.2) untill the tile cover at least one point from each task group.
#' @param splitNET A \code{logical} value specifying whether to split the network by tiles (default: FALSE).
#'
#' @details
#' The function uses a Nearest Neighbour heuristic to solve the GTSP, selecting one location per task group in a 3D network to minimise total travel cost, starting and ending at \code{startPT}. If \code{nID} is missing in \code{startPT} or \code{taskPT}, the nearest node is assigned using \code{\link[GISnetwork3D]{nearestNode3d}}. Spatial tiling (\code{nrows}, \code{ncols}, \code{buf}) partitions data for parallel processing with \code{ncores} workers, with \code{coef} expanding the buffer iteratively to ensure all task groups are covered. If \code{nrows} and \code{ncols} are 0, no tiling is applied and the results are computed with single core. User can decide whether to also tile the network with the splitNET argument. Tiling network may speed up the process if the spatial extent and the network data are VERY LARGE. However, it can also slow down the process if the network and spatial extent are moderate to small.
#'
#' @return Depends on \code{edges.output}:
#' \describe{
#'   \item{\code{NULL}}{A \code{data.frame} with columns: \code{start.ID} (\code{character} or \code{numeric}, starting point ID), travel costs, selected location \code{ID} and \code{nID} per task.}
#'   \item{\code{dissolved}}{An \code{sf} object with 3D \code{LINESTRING} geometry (dissolved by \code{start.ID}), including travel costs, selected \code{ID} and \code{nID} per task.}
#'   \item{\code{raw}}{A \code{list} with: \code{access} (a \code{data.frame} as above) and \code{eID} (an \code{sf} object with 3D \code{LINESTRING} geometries of traversed edges, with \code{start.ID}).}
#'   \item{\code{raw.df}}{A \code{list} with: \code{access} (as above) and \code{eID} (a \code{data.frame} of traversed edge attributes, with \code{start.ID}).}
#' }
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(GISnetwork3D)
#'
#' # Create a network
#' coords = list(
#'   rbind(c(0, 0, 0), c(3, 3, 1)),      # Line 1
#'   rbind(c(0, 0, 0), c(3, -1, 1)),     # Line 2
#'   rbind(c(3, 3, 1), c(6, 3, 2)),      # Line 3
#'   rbind(c(3, -1, 1), c(6, -1, 2)),    # Line 4
#'   rbind(c(3, 3, 1), c(3, 6, 3)),      # Line 5
#'   rbind(c(6, 3, 2), c(6, 6, 3)),      # Line 6
#'   rbind(c(6, -1, 2), c(6, -4, 3)),    # Line 7
#'   rbind(c(3, -1, 1), c(0, -4, 2)),     # Line 8
#'   rbind(c(3, 6, 3), c(0, 6, 4)),      # Line 9
#'   rbind(c(0, 0, 0), c(0, 6, 4)),      # Line 10
#'   rbind(c(0, -4, 2), c(-3, -4, 1)),   # Line 11
#'   rbind(c(6, -4, 3), c(9, -4, 2))     # Line 12
#' )
#'
#' # Create sf LINESTRING
#' lines = lapply(coords, function(x) st_linestring(x))
#' road_network = st_sf(geometry = st_sfc(lines))
#' road_network$ID = 1:nrow(road_network)
#'
#' # Create points to visit
#' bbox = st_bbox(road_network) # Get the bounding box of the road network
#' n_points = 15 # Set the number of points
#'
#' # Generate random coordinates
#' x = runif(n_points, min = bbox["xmin"], max = bbox["xmax"])
#' y = runif(n_points, min = bbox["ymin"], max = bbox["ymax"])
#' z = runif(n_points, min = 0, max = 10)  # Z values between 0 and 10
#'
#' # Create sf POINT object
#' points_sf = data.frame(x, y, z)
#' points_sf = st_as_sf(points_sf, coords = c("x", "y", "z"))
#' points_sf$ID = 1:nrow(points_sf) # Assign ID
#' points_sf$tasks = sample(c("school", "transit", "grocery"), 15, replace = TRUE) # Stratify point by tasks
#'
#' # Create the origin points
#' n_points = 4 # Set the number of points
#' x = runif(n_points, min = bbox["xmin"], max = bbox["xmax"])
#' y = runif(n_points, min = bbox["ymin"], max = bbox["ymax"])
#' z = runif(n_points, min = 0, max = 10)  # Z values between 0 and 10
#' Origin_points_sf = data.frame(x, y, z)
#' Origin_points_sf = st_as_sf(Origin_points_sf, coords = c("x", "y", "z"))
#' Origin_points_sf$ID = 1:nrow(Origin_points_sf) # Assign ID
#'
#' # Prepare the network data
#' road_network = GISnetwork3D::blendPT3d(Origin_points_sf, road_network)
#' road_network = GISnetwork3D::blendPT3d(points_sf, road_network)
#' road_network = GISnetwork3D::splitNearestV3d(Origin_points_sf, road_network)
#' road_network = GISnetwork3D::splitNearestV3d(points_sf, road_network)
#' road_network$direction = 0
#' road_network = GISnetwork3D::sf_to_net3d.full(road_network, direction = "direction") # Activate network
#'
#' # Assign node IDs to points
#' Origin_points_sf = GISnetwork3D::nearestNode3d(road_network$nodes, Origin_points_sf)
#' points_sf = GISnetwork3D::nearestNode3d(road_network$nodes, points_sf)
#' points_sf = split(points_sf, points_sf$tasks) # Convert point to a list stratified by the tasks
#'
#' # Compute the access that travel through the optimised location for each tasks
#' access3d.multitasks(startPT = Origin_points_sf, taskPT = points_sf, net = road_network, edges.output = "dissolved", weight = "time", ncores = 1)
#' }
#' @export

access3d.multitasks = function(startPT, taskPT, net, weight, edges.output = NULL, buf = 2000, nrows = 4, ncols = 4, ncores = 2, coef = 1.2, splitNET = FALSE){
  # --- Input validation ---
  if (!inherits(startPT, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(startPT))) {
    stop("startPT must be an <sf> object with POINT geometry.") }

  if (!("ID" %in% names(startPT))) { stop("Invalid input for startPT: A field 'ID' is missing.") }

  if (!is.list(taskPT) || any(!sapply(taskPT, function(x) inherits(x, "sf") && "sfc_POINT" %in% class(sf::st_geometry(x))))) {
    stop("taskPT must be a list containing more than one <sf> objects with POINT geometry.")
  }

  if(any(!sapply(taskPT, function(x) "ID" %in% names(x)))){
    stop("Missing field 'ID' in the <sf> objects insite the taskPT")
  }

  if (!inherits(net, "list") || !"nodes" %in% names(net) || !"edges" %in% names(net)) {
    stop("net must be a list with 'nodes' and 'edges'.")}

  if (!(is.null(edges.output) || (is.character(edges.output) && length(edges.output) == 1 && edges.output %in% c("dissolved", "raw", "raw.df")))) {
    stop("edges.output must be NULL or either from c('dissolved', 'raw', 'raw.df').")
  }

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a character indicating the column name from net's edges for weight.")}

  if (!is.numeric(nrows) || length(nrows) != 1 || nrows < 0 || nrows %% 1 != 0) {
    stop("nrows must be a non-negative integer value.")}

  if (!is.numeric(ncols) || length(ncols) != 1 || ncols < 0 || ncols %% 1 != 0) {
    stop("ncols must be a non-negative integer value.")}

  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores != round(ncores)) {
    stop("ncores must be a positive integer.")}

  if (!is.numeric(coef) || length(coef) != 1) {
    stop("coef must be a numeric value.")}

  if (!is.logical(splitNET)) {
    stop("splitNET must be a <logical> value (TRUE/FALSE).")}

  # --- Prepare data ---
  # Search for nID for each StartPT
  if (!"nID" %in% colnames(startPT)) {startPT = GISnetwork3D::nearestNode3d(net$nodes, startPT)}
  sf::st_geometry(startPT) = "geom"

  # Search for nID for each taskPT if needed
  if(any(!sapply(taskPT, function(x) "nID" %in% names(x)))){
    for(i in 1:length(taskPT)){
      sf::st_geometry(taskPT[[i]]) = "geom"
      taskPT[[i]]$cluster = i
      taskPT[[i]] = GISnetwork3D::nearestNode3d(net$nodes, taskPT[[i]])
      taskPT[[i]] = taskPT[[i]][c("nID", "ID", "cluster")] }
  } else {
    for(i in 1:length(taskPT)){
      sf::st_geometry(taskPT[[i]]) = "geom"
      taskPT[[i]]$cluster = i
      taskPT[[i]] = taskPT[[i]][c("nID", "ID", "cluster")]}
  }

  # Combine taskPT into one sf POINT object
  taskPT = purrr::reduce(taskPT, rbind)

  # --- Main function implementation ---

  ### --- Compute accessibility without tiling ---
  if(nrows == 0 & ncols == 0){
    result = GISnetwork3D::access3d.multitasks.sTask(startPT = startPT, taskPT = taskPT, net = net, weight = weight,
                                                     export.IDs = T, export.nIDs = T, edges.output = edges.output)
    return(result) } # <---- OUTPUT: data.frame

  ### --- Compute Accessibility with tiling ---
  if(nrows >= 1 & ncols >= 1){

    # Compute with O and D being tiled
    if(isFALSE(splitNET)){
      # Split Origin and Destination by Tiles
      input = GISnetwork3D::tile.multitasks.OD(o = startPT, d = taskPT, buf = buf, nrows = nrows, ncols = ncols, returnTiles = TRUE, returnPT = T)

      future::plan("future::multisession", workers = ncores)
      result = furrr::future_map2(input$o, input$d,
                                  function(x, y, net, weight, export.IDs, export.nIDs, edges.output){
                                    if(!is.null(edges.output)){if(edges.output == "raw"){edges.output = "raw.df"} }
                                    result = GISnetwork3D::access3d.multitasks.sTask(startPT = x, taskPT = y, net = net, weight = weight,
                                                                                     export.IDs = T, export.nIDs = T, edges.output = edges.output)
                                    result$oID = x$id[result$oID]
                                    result$dID = y$id[result$dID]
                                    return(result) },
                                  net = net, weight = weight, edges.output = edges.output,
                                  .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf")))
      future::plan("future::sequential") }

    # Compute with O, D and net being tiled
    if(isTRUE(splitNET)){
      # Split Origin, Destination and net by Tiles
      input = GISnetwork3D::tile.multitasks.ODN(o = startPT, d = taskPT, n = net, buf = buf, nrows = nrows, ncols = ncols, returnTiles = TRUE, returnPT = T)

      future::plan("future::multisession", workers = ncores)
      result = furrr::future_pmap(list(input$o, input$d, input$n),
                                  function(x, y, net, weight, export.IDs, export.nIDs, edges.output){
                                    if(!is.null(edges.output)){if(edges.output == "raw"){edges.output = "raw.df"} }
                                    result = GISnetwork3D::access3d.multitasks.sTask(startPT = x, taskPT = y, net = net, weight = weight,
                                                                                     export.IDs = T, export.nIDs = T, edges.output = edges.output)
                                    result$oID = x$id[result$oID]
                                    result$dID = y$id[result$dID]
                                    return(result) },
                                  weight = weight, edges.output = edges.output,
                                  .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf")))
      future::plan("future::sequential") }

    if(is.null(edges.output)){
      # Export Result as a Data.frame
      result = result %>% purrr::reduce(rbind)
      result = result[match(startPT$ID, result$start.ID), ]
      rownames(result) = NULL
      return(result) # <---- OUTPUT: data.frame

      } else if(edges.output == "dissolved"){
        # Export Result as an SF
        result = result %>% purrr::reduce(rbind)
        result = result[match(startPT$ID, result$start.ID), ]
        rownames(result) = NULL
        return(result) # <---- OUTPUT: sf
      } else {
        # Export Result as a List
        result.access = purrr::map(result, ~.x$access) %>% purrr::reduce(rbind)
        result.access = result.access[match(startPT$ID, result.access$start.ID), ]
        result.eID = purrr::map(result, ~.x$eID) %>% purrr::reduce(rbind)
        result.eID = data.frame(start.ID = startPT$ID) %>% dplyr::left_join(result.eID, by = "start.ID")
        rownames(result.access) = NULL
        rownames(result.eID) = NULL

        if(edges.output == "raw"){
          result.eID = c(start.ID = result.eID$start.ID) %>% cbind(net$edges[match(result.eID$eID, net$edges$eID),])
          names(result.eID)[1] = "start.ID"
          rownames(result.eID) = NULL }

        return(list(access = result.access, eID = result.eID))} } # <---- OUTPUT: List
}


##_____________________________________________________________________________
#' Summarise Edge Attributes Along Shortest 3D Network Paths for Multiple Locations
#'
#' This function solves the Generalised Traveling Salesperson Problem (GTSP) in a 3D network using a Nearest Neighbour heuristic, selecting one location per task group to minimise travel cost, and summarises edge attributes along the resulting routes.
#'
#' @param startPT An \code{sf} object with 3D \code{POINT} geometry representing starting points, with a \code{character} or \code{numeric} column \code{ID} (unique identifier) and optional \code{character} column \code{nID} (nearest network node).
#' @param taskPT A \code{list} of \code{sf} objects with 3D \code{POINT} geometry, each representing a task group (e.g., grocery stores), with a \code{character} or \code{numeric} column \code{ID} and optional \code{character} column \code{nID}.
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param geom A \code{logical} value indicating whether to include route geometry in the output (default: \code{FALSE}).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#' @param buf A non-negative \code{numeric} value specifying the initial buffer distance (in map units) for spatial tiling (default: 2000).
#' @param nrows A non-negative \code{integer} value specifying the number of rows for spatial tiling (default: 1).
#' @param ncols A non-negative \code{integer} value specifying the number of columns for spatial tiling (default: 2).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 2).
#' @param coef A \code{numeric} value to multiply the buffer distance per iteration of spatial tiling (default: 1.2).
#' @param sumV A \code{character} vector of edge column names for simple sum along each route, or \code{NULL} (default).
#' @param w.sumV A \code{character} vector of edge column names for weighted sum, or \code{NULL} (default).
#' @param w.sumW A \code{character} value specifying the edge column for weighting \code{w.sumV}, or \code{NULL} (default).
#' @param p.sumV A \code{character} vector of edge column names for normalised sum, or \code{NULL} (default).
#' @param pw.sumV A \code{character} vector of edge column names for weighted-normalised sum, or \code{NULL} (default).
#' @param p.norm A \code{character} value specifying the edge column for normalisation of \code{p.sumV} or \code{pw.sumV}, or \code{NULL} (default).
#' @param splitNET A \code{logical} value specifying whether to split the network by tiles (default: FALSE).
#'
#' @details
#' The function solves the GTSP in a 3D network using a Nearest Neighbour heuristic, selecting one location per task group in \code{taskPT} to minimise travel cost from \code{startPT}, and computes edge attribute summaries along each route. Supported summaries include:
#' \itemize{
#'   \item Simple sum (\code{sumV}): Sum of each column’s values across route edges.
#'   \item Weighted sum (\code{w.sumV}, \code{w.sumW}): Sum of (\code{w.sumV} * \code{w.sumW}) across route edges.
#'   \item Normalised sum (\code{p.sumV}, \code{p.norm}): Sum of \code{p.sumV} divided by sum of \code{p.norm} across route edges.
#'   \item Weighted-normalised sum (\code{pw.sumV}, \code{p.norm}): Sum of (\code{pw.sumV} * \code{p.norm}) divided by sum of \code{p.norm} across route edges.
#' }
#' If \code{nID} is missing in \code{startPT} or \code{taskPT}, the nearest node is assigned using \code{\link[GISnetwork3D]{nearestNode3d}}. Spatial tiling (\code{nrows}, \code{ncols}, \code{buf}, \code{coef}) partitions data for parallel processing with \code{ncores} workers. If \code{nrows} and \code{ncols} are 0, no tiling is applied and the results are computed with single core. User can decide whether to also tile the network with the splitNET argument. Tiling network may speed up the process if the spatial extent and the network data are VERY LARGE. However, it can also slow down the process if the network and spatial extent are moderate to small. The output includes one row per route, with columns for \code{start.ID} and requested summaries. If \code{geom = TRUE}, 3D \code{LINESTRING} geometries are included.
#'
#' @return An \code{sf} object with 3D \code{LINESTRING} geometry (if \code{geom = TRUE}) or a \code{data.frame} (if \code{geom = FALSE}) with columns:
#' \describe{
#'   \item{\code{start.ID}}{\code{character} or \code{numeric}, starting point identifier.}
#'   \item{\code{...}}{\code{numeric}, columns for each requested edge attribute summary (e.g., \code{sumV}, \code{w.sumV}, \code{p.sumV}, \code{pw.sumV}).}
#'   \item{\code{geometry}}{3D \code{LINESTRING}, route geometry (if \code{geom = TRUE}).}
#' }
#'
#' @examples
#' \dontrun{
#'
#' # Create a network
#' coords = list(
#'   rbind(c(0, 0, 0), c(3, 3, 1)),      # Line 1
#'   rbind(c(0, 0, 0), c(3, -1, 1)),     # Line 2
#'   rbind(c(3, 3, 1), c(6, 3, 2)),      # Line 3
#'   rbind(c(3, -1, 1), c(6, -1, 2)),    # Line 4
#'   rbind(c(3, 3, 1), c(3, 6, 3)),      # Line 5
#'   rbind(c(6, 3, 2), c(6, 6, 3)),      # Line 6
#'   rbind(c(6, -1, 2), c(6, -4, 3)),    # Line 7
#'   rbind(c(3, -1, 1), c(0, -4, 2)),     # Line 8
#'   rbind(c(3, 6, 3), c(0, 6, 4)),      # Line 9
#'   rbind(c(0, 0, 0), c(0, 6, 4)),      # Line 10
#'   rbind(c(0, -4, 2), c(-3, -4, 1)),   # Line 11
#'   rbind(c(6, -4, 3), c(9, -4, 2))     # Line 12
#' )
#'
#' # Create sf LINESTRING
#' lines = lapply(coords, function(x) st_linestring(x))
#' road_network = st_sf(geometry = st_sfc(lines))
#' road_network$ID = 1:nrow(road_network)
#'
#' # Create points to visit
#' bbox = st_bbox(road_network) # Get the bounding box of the road network
#' n_points = 15 # Set the number of points
#'
#' # Generate random coordinates
#' x = runif(n_points, min = bbox["xmin"], max = bbox["xmax"])
#' y = runif(n_points, min = bbox["ymin"], max = bbox["ymax"])
#' z = runif(n_points, min = 0, max = 10)  # Z values between 0 and 10
#'
#' # Create sf POINT object
#' points_sf = data.frame(x, y, z)
#' points_sf = st_as_sf(points_sf, coords = c("x", "y", "z"))
#' points_sf$ID = 1:nrow(points_sf) # Assign ID
#' points_sf$tasks = sample(c("school", "transit", "grocery"), 15, replace = TRUE) # Stratify point by tasks
#'
#' # Create the origin points
#' n_points = 4 # Set the number of points
#' x = runif(n_points, min = bbox["xmin"], max = bbox["xmax"])
#' y = runif(n_points, min = bbox["ymin"], max = bbox["ymax"])
#' z = runif(n_points, min = 0, max = 10)  # Z values between 0 and 10
#' Origin_points_sf = data.frame(x, y, z)
#' Origin_points_sf = st_as_sf(Origin_points_sf, coords = c("x", "y", "z"))
#' Origin_points_sf$ID = 1:nrow(Origin_points_sf) # Assign ID
#'
#' # Prepare the network data
#' road_network = GISnetwork3D::blendPT3d(Origin_points_sf, road_network)
#' road_network = GISnetwork3D::blendPT3d(points_sf, road_network)
#' road_network = GISnetwork3D::splitNearestV3d(Origin_points_sf, road_network)
#' road_network = GISnetwork3D::splitNearestV3d(points_sf, road_network)
#' road_network$direction = 0
#' road_network = GISnetwork3D::sf_to_net3d.full(road_network, direction = "direction") # Activate network
#'
#' # Find nodes for each points
#' Origin_points_sf = GISnetwork3D::nearestNode3d(road_network$nodes, Origin_points_sf)
#' points_sf = GISnetwork3D::nearestNode3d(road_network$nodes, points_sf)
#' points_sf = split(points_sf, points_sf$tasks) # Convert point to a list stratified by the tasks
#'
#' # Compute the access that travel through the optimised location for each tasks
#' access3d.multitasks.summary(startPT = Origin_points_sf, taskPT = points_sf, net = road_network, geom = TRUE, weight = "time", nrows = 4, ncols = 4, ncores = 2, sumV = c("time", "energy"))
#' }
#'
#' @export
access3d.multitasks.summary = function(startPT, taskPT, net, geom = FALSE, weight, buf = 2000, nrows = 1, ncols = 2, ncores = 2, coef = 1.2,
                                       sumV = NULL,
                                       w.sumV = NULL, w.sumW = NULL,
                                       p.sumV = NULL, pw.sumV = NULL, p.norm = NULL, splitNET = FALSE){

  # --- Input validation ---
  # startPT
  if (!inherits(startPT, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(startPT))) {
    stop("startPT must be an <sf> object with POINT geometry.") }
  if (!"ID" %in% names(startPT)) {
    stop("startPT is missing required column 'ID'.") }

  # taskPT
  if (!is.list(taskPT) || length(taskPT) == 0 ||
      any(!sapply(taskPT, function(x) inherits(x, "sf") && "sfc_POINT" %in% class(sf::st_geometry(x))))) {
    stop("taskPT must be a non-empty list of <sf POINT> objects.") }
  if (any(!sapply(taskPT, function(x) "ID" %in% names(x)))) {
    stop("One or more objects in taskPT are missing required column 'ID'.") }

  # net
  if (!is.list(net) || !all(c("nodes", "edges") %in% names(net))) {
    stop("net must be a list containing 'nodes' and 'edges'.") }

  # Scalar parameters
  if (!is.logical(geom) || length(geom) != 1) { stop("geom must be a single logical value (TRUE/FALSE).") }

  if (!is.character(weight) || length(weight) != 1 || !(weight %in% names(net$edges))) {
    stop("weight must be a single column name present in net$edges.") }

  if (!is.numeric(nrows) || length(nrows) != 1 || nrows < 0 || nrows != round(nrows)) {stop("nrows must be a non-negative integer.")}
  if (!is.numeric(ncols) || length(ncols) != 1 || ncols < 0 || ncols != round(ncols)) {stop("ncols must be a non-negative integer.")}
  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores != round(ncores)) {stop("ncores must be a positive integer.")}
  if (!is.numeric(coef) || length(coef) != 1) {stop("coef must be a single numeric value.")}

  # Optional named parameters (character checks)
  if (!is.null(sumV) && !is.character(sumV)) { stop("sumV must be NULL or character vector.") }
  if (!is.null(w.sumV) && !is.character(sumV)) { stop("w.sumV must be NULL or a character vector.") }
  if (!is.null(p.sumV) && !is.character(p.sumV)) { stop("p.sumV must be NULL or a character vector.") }
  if (!is.null(pw.sumV) && !is.character(pw.sumV)) { stop("pw.sumV must be NULL or a character vector.") }

  if (!is.null(w.sumW) && (!is.character(w.sumW) || length(w.sumW) > 1)) {stop("w.sumW must be NULL or a single character string.")}
  if (!is.null(p.norm) && (!is.character(p.norm) || length(p.norm) > 1)) {stop("p.norm must be NULL or a single character string.")}

  if (!is.logical(splitNET) || length(splitNET) != 1) { stop("splitNET must be a single logical value (TRUE/FALSE).") }

  # Dependency rules
  if (!is.null(p.sumV)  && is.null(p.norm)) stop("p.norm is required when p.sumV is provided.")
  if (!is.null(pw.sumV) && is.null(p.norm)) stop("p.norm is required when pw.sumV is provided.")
  if (!is.null(p.norm)  && is.null(p.sumV) && is.null(pw.sumV)) {stop("At least one of p.sumV or pw.sumV is required when p.norm is provided.")}
  if (!is.null(w.sumV)  != !is.null(w.sumW)) {stop("w.sumV and w.sumW must both be provided or both be NULL.")}


  # --- Prepare data ---
  # Search for nID for each StartPT
  if (!"nID" %in% colnames(startPT)) {startPT = GISnetwork3D::nearestNode3d(net$nodes, startPT)}
  sf::st_geometry(startPT) = "geom"

  # Search for nID for each taskPT if needed
  if(any(!sapply(taskPT, function(x) "nID" %in% names(x)))){
    for(i in 1:length(taskPT)){
      sf::st_geometry(taskPT[[i]]) = "geom"
      taskPT[[i]]$cluster = i
      taskPT[[i]] = GISnetwork3D::nearestNode3d(net$nodes, taskPT[[i]])
      taskPT[[i]] = taskPT[[i]][c("nID", "ID", "cluster")] }
  } else {
    for(i in 1:length(taskPT)){
      sf::st_geometry(taskPT[[i]]) = "geom"
      taskPT[[i]]$cluster = i
      taskPT[[i]] = taskPT[[i]][c("nID", "ID", "cluster")]}
  }

  # Combine taskPT into one sf POINT object
  taskPT = purrr::reduce(taskPT, rbind)

  # --- Main function implementation ---
  # Compute without tiling
  if(nrows == 0 & ncols == 0){
    result = GISnetwork3D::access3d.multitasks.summary.sTask(startPT = startPT, taskPT = taskPT, net = net, geom = geom, weight = weight,
                                                             sumV = sumV, w.sumV = w.sumV, w.sumW = w.sumW, p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm)
    return(result)} # <--- Output: sf or df

  ### --- Compute with tiling ---
  if(nrows >= 1 & ncols >= 1){

    # Compute with O and D being tiled
    if(isFALSE(splitNET)){
      # Split Origin and Destination by Tiles
      input = GISnetwork3D::tile.multitasks.OD(o = startPT, d = taskPT, buf = buf, nrows = nrows, ncols = ncols, returnTiles = TRUE, returnPT = T)

      future::plan("future::multisession", workers = ncores)
      result = furrr::future_map2(input$o, input$d,
                                  function(x, y, net, geom, weight, sumV, w.sumV, w.sumW, p.sumV, pw.sumV, p.norm){

                                    result = GISnetwork3D::access3d.multitasks.summary.sTask(startPT = x, taskPT = y, net = net, geom = geom, weight = weight,
                                                                                             sumV = sumV, w.sumV = w.sumV, w.sumW = w.sumW, p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm)
                                    return(result) },
                                  net = net, geom = geom, weight = weight, sumV = sumV, w.sumV = w.sumV, w.sumW = w.sumW, p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm,
                                  .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf")))
      future::plan("future::sequential") }


    if(isFALSE(splitNET)){
      # Compute with O, D and net being tiled
      input = GISnetwork3D::tile.multitasks.ODN(o = startPT, d = taskPT, n = net, buf = buf, nrows = nrows, ncols = ncols, returnTiles = TRUE, returnPT = T)

      future::plan("future::multisession", workers = ncores)
      result = furrr::future_pmap(list(input$o, input$d, input$n),
                                  function(x, y, n, geom, weight, sumV, w.sumV, w.sumW, p.sumV, pw.sumV, p.norm){

                                    result = GISnetwork3D::access3d.multitasks.summary.sTask(startPT = x, taskPT = y, net = n, geom = geom, weight = weight,
                                                                                             sumV = sumV, w.sumV = w.sumV, w.sumW = w.sumW, p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm)
                                    return(result) },
                                  geom = geom, weight = weight, sumV = sumV, w.sumV = w.sumV, w.sumW = w.sumW, p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm,
                                  .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "GISnetwork3D", "sf")))
      future::plan("future::sequential") }

    # Export Results
    result = result %>% purrr::reduce(rbind)
    result = result[match(startPT$ID, result$start.ID),]
    rownames(result) = NULL
    return(result) }  # <--- Output: sf or df
}


##_____________________________________________________________________________
#' Extract Edge Indices for Shortest Multi-Stop Paths in a 3D Network
#'
#' Computes shortest paths for sequences of nodes (multi-stop routes) in a 3D network and returns the edge indices (\code{eID}) for each path.
#'
#' @param net An \code{igraph} object or a \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param nodes A \code{character} vector (node IDs), or a \code{data.frame} or \code{sf} object with 3D \code{POINT} geometry, containing columns \code{ID} and \code{nID} (\code{character}, node identifier).
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination) or \code{"in"} (destination to origin).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#' @param output A \code{character} value specifying the output format: \code{"list"} (list of edge indices) or \code{"df"} (data frame with edge indices) (default: \code{"list"}).
#'
#' @details
#' For each sequence of nodes (grouped by \code{ID}), the function computes the shortest path between consecutive node pairs in a 3D network using the \code{weight} attribute, extracting the indices (\code{eID}) of traversed edges. If \code{nodes} is a \code{character} vector, it is treated as a single sequence with \code{ID = 1}. The function supports both \code{igraph} and \code{sf_to_net3d} network inputs, converting the latter to \code{igraph} internally.
#'
#' @return Depends on \code{output}:
#' \describe{
#'   \item{\code{list}}{A named \code{list} of \code{integer} vectors containing edge indices (\code{eID}), with names corresponding to \code{ID} (group identifier).}
#'   \item{\code{df}}{A \code{data.frame} with columns: \code{path.ID} (\code{character} or \code{numeric}, group identifier) and \code{eID} (\code{numeric}, edge index).}
#' }
#'
#' @examples
#' \dontrun{
#'
#' # Create points to traverse
#' pt = st_sf(geometry = st_sfc(
#'        st_point(c(0, 0, 1)),
#'        st_point(c(3, 4, 3.5)),
#'        st_point(c(410, 550, 5.1)),
#'        st_point(c(181, 281, 4.1)),
#'        st_point(c(91, 98, 02.9)) ))
#'
#' pt$ID = c(1,1,2,2,2)
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 30,500,1.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(30,500,1.5, 60,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(60,100,3, 90,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(90,100,3, 100,120,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,120,4, 180,280,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(180,280,4, 400,560,5), ncol=3, byrow=TRUE))
#' ))
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create 3D network with Tobler's hiking function for cost
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for pt
#' l = GISnetwork3D::blendPT3d(pt, l) # Blend point to network
#' pt = nearestNode3d(net$nodes, pt)
#'
#' # Compute paths
#' result = path3d.multistops.get_eID(net, pt, mode = "out", "time")
#'
#' print(result)
#' }
#'
#' @export

path3d.multistops.get_eID = function(net, nodes, mode = "out", weight, output = "list"){
  # --- Input validation ---
  if (!inherits(net, c("igraph", "list"))) {
    stop("net must be an igraph or a list with 'nodes' and 'edges'.")}

  if (!is.character(nodes) && !is.data.frame(nodes)) {
    stop("nodes must be a character vector, <sf> object with POINT geometry or data.frame.")}

  if(is.data.frame(nodes)){
    if(!"ID" %in% names(nodes) || !"nID" %in% names(nodes)){
      stop("nodes must contain fields named 'ID' and 'nID'.") } }

  if (!is.character(mode) || length(mode) != 1 || !mode %in% c("out", "in")) {
    stop("mode must be one of 'out' or 'in'.")}

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a single character string indicating the cost column from net's edges.")}

  if (inherits(net, "list") && !(weight %in% names(net$edges))) {
    stop(paste0("Weight column '", weight, "' not found in net$edges."))
  }
  if (inherits(net, "igraph") && !(weight %in% igraph::edge_attr_names(net))) {
    stop(paste0("Weight attribute '", weight, "' not found in igraph edges."))
  }

  if (!is.character(output) || length(output) != 1 || !output %in% c("list", "df")) {
    stop("output must be either 'list' or 'df'.") }

  # --- Network preparation ---
  # Convert 'net' to igraph
  if (inherits(net, "list")) {
    # Convert to igraph
    net = igraph::graph_from_data_frame(d = net$edges, directed = T, vertices = net$nodes %>% st_drop_geometry()) }

  # Split the nodes by group ID
  if(is.vector(nodes)){
    nodes = data.frame(ID = 1, nID = nodes)
    ID.type = class(nodes$ID) # Get the data type of the input ID from the nodes
    nodes = split(nodes, nodes$ID)
  } else {
    ID.type = class(nodes$ID) # Get the data type of the input ID from the nodes
    nodes = nodes %>% dplyr::left_join( data.frame(ID = unique(nodes$ID), split.Order = 1:length(unique(nodes$ID))), by = "ID" )
    nodes = split(nodes, nodes$split.Order)
  }

  # --- Output processing ---
  # Initialise a list to store the edges for each shortest path
  edges = list()

  # Loop through each group of nodes
  for(i in 1:length(nodes)){
    nodes_len = length(nodes[[i]]$nID) # Get number of nodes in the group
    # Get edges ID for each pair of nodes (origin-destination) in the group
    edges[[i]] = GISnetwork3D::path3d.pairOD.get_eID(net = net,
                                                     oNode = nodes[[i]]$nID[seq(1, nodes_len - 1, 1)],
                                                     dNode = nodes[[i]]$nID[seq(2, nodes_len, 1)],
                                                     mode = mode, weight = weight, output = "list") %>% unlist(use.names = F) %>% unlist() }
  names(edges) = unlist(purrr::map(nodes, ~ .x$ID[1]))

  # --- Output processing ---
  if(output == "list"){
    # Return eID as a list; names represent unique IDs for each path
    return(edges) # <---- OUTPUT: list
  }  else if(output == "df"){
    # Return eID as a data frame; path.ID represents unique IDs for each path
    edges = data.frame(path.ID = rep(names(edges), base::unlist(purrr::map(edges, length), use.names = F)),
                       eID = base::unlist(edges, use.names = F))
    class(edges$path.ID) = ID.type
    return(edges)}  # <---- OUTPUT: data.frame
}



##_____________________________________________________________________________
#' Extract Edges for Shortest Multi-Stop Paths in a 3D Network
#'
#' Computes shortest paths for sequences of nodes (multi-stop routes) in a 3D network and extracts the network edges as \code{sf} \code{LINESTRING} objects or \code{data.frame}s, with support for parallel processing.
#'
#' @param net An \code{igraph} object or a \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param nodes A \code{character} vector (node IDs), or a \code{data.frame} or \code{sf} object with 3D \code{POINT} geometry, containing columns \code{ID} and \code{nID} (\code{character}, node identifier).
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination) or \code{"in"} (destination to origin).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 4).
#' @param tasks A positive \code{integer} value specifying the number of task chunks for parallel computation (default: 4).
#' @param output A \code{character} value specifying the output format: \code{"eID"} (edge IDs), \code{"df"} (edge attributes), \code{"raw"} (edge geometries), or \code{"dissolve"} (dissolved route geometries) (default: \code{"raw"}).
#'
#' @details
#' For each sequence of nodes (grouped by \code{ID}), the function computes shortest paths between consecutive node pairs in a 3D network using the \code{weight} attribute, extracting traversed edges. If \code{nodes} is a \code{character} vector, it is treated as a single sequence with \code{ID = 1}. The function supports \code{igraph} or \code{sf_to_net3d} network inputs, converting the latter to \code{igraph} internally. Parallel processing splits node sequences into \code{tasks} chunks across \code{ncores} workers. Outputs include edge indices, attributes, or 3D \code{LINESTRING} geometries.
#'
#' @return Depends on \code{output}:
#' \describe{
#'   \item{\code{eID}}{A \code{data.frame} with columns: \code{path.ID} (\code{character} or \code{numeric}, group identifier) and \code{eID} (\code{numeric}, edge index).}
#'   \item{\code{df}}{A \code{data.frame} with edge attributes and \code{path.ID} (\code{character} or \code{numeric}).}
#'   \item{\code{raw}}{An \code{sf} object with 3D \code{LINESTRING} geometries, edge attributes, and \code{path.ID} (\code{character} or \code{numeric}).}
#'   \item{\code{dissolve}}{An \code{sf} object with dissolved 3D \code{LINESTRING} geometries grouped by \code{path.ID} (\code{character} or \code{numeric}).}
#' }
#'
#' @examples
#' \dontrun{
#'
#' #' # Create origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' d = st_sf(geometry = st_sfc(
#'   st_point(c(410, 550, 5.1)),
#'   st_point(c(181, 281, 4.1)),
#'   st_point(c(91, 98, 02.9))
#' ))
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 30,500,1.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(30,500,1.5, 60,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(60,100,3, 90,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(90,100,3, 100,120,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,120,4, 180,280,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(180,280,4, 400,560,5), ncol=3, byrow=TRUE))
#' ))
#' l = blendPT3d(rbind(o,d), l) # Blend points to the line network
#' l = splitNearestV3d(rbind(o,d), l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create 3D network with Tobler's hiking function for cost
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for o and d
#' o = nearestNode3d(net$nodes, o)
#' d = nearestNode3d(net$nodes, d)
#'
#' # Compute accessibility
#' result = access3d.pairOD(net, o, d, mode = "out", "time", threshold = 900, avaiW = c(1,1,5))
#'
#' # Get the path between each pair of origin and destination
#' result = path3d.multistops.getEdges(net, pt, "out", "time", ncores = 2, tasks = 2)
#' print(result)
#' }
#'
#' @export
path3d.multistops.getEdges = function(net, nodes, mode = "out", weight, ncores = 4, tasks = 4, output = "raw"){
  # --- Input validation ---
  if (!inherits(net, "list") || !"nodes" %in% names(net) || !"edges" %in% names(net)) {
    stop("net must be a list with 'nodes' and 'edges'.")}

  if (!is.character(nodes) && !is.data.frame(nodes)) {
    stop("nodes must be a character vector, <sf> object with POINT geometry or data.frame.")}

  if(is.data.frame(nodes)){
    if(!"ID" %in% names(nodes) || !"nID" %in% names(nodes)){
      stop("nodes must contain fields named 'ID' and 'nID'.") } }

  if (!mode %in% c("out", "in")) {stop("mode must be one of 'out' or 'in'.")}

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a character indicating the column name in net's edges for weights.")}

  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores != round(ncores)) {
    stop("ncores must be a positive integer.")}

  if (!is.numeric(tasks) || length(tasks) != 1 || tasks <= 0 || tasks != round(tasks)) {
    stop("tasks must be a positive integer.")}

  if(!is.character(output) || length(output) != 1){
    stop("output must be one of 'eID', 'df', 'raw' or 'dissolve'.")
  }
  if (!output %in% c("eID", "df", "raw", "dissolve")) {
    stop("output must be one of 'eID', 'df', 'raw' or 'dissolve'.")}

  # Split input nodes to n number of tasks
  if(is.vector(nodes)){
    nodes = data.frame(ID = 1, nID = nodes) }

  unique.ID = unique(nodes$ID)
  unique.ID = GISnetwork3D::split.data(unique.ID, tasks)
  nodes = purrr::map(unique.ID,
                     function(x,y){
                       y = y[y$ID %in% x,]
                       }, y = nodes)

  # --- Network preparation ---
  # Convert 'net' to igraph
  net.edges = net$edges
  net = igraph::graph_from_data_frame(d = net$edges, directed = T, vertices = net$nodes %>% st_drop_geometry())

  # --- Output processing ---
  # --- Extract network edges ID using multiple cores for each set of nodes---
  future::plan("future::multisession", workers = ncores) # Set up parallel plan

  # Map over the split nodes, executing path computation for each chunk
  edges = furrr::future_map( nodes, ~ GISnetwork3D::path3d.multistops.get_eID(net = net, nodes = .x, mode = mode, weight = weight, output = "df"),
                             .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "sf", "GISnetwork3D")) ) %>% purrr::reduce(rbind)
  future::plan("future::sequential") # Stop parallel cluster

  if(output == "eID"){
    # If output = "eID", return data.frame of the eID
    return(edges) }

  if(output == "df"){
    return(cbind(net.edges[match(edges$eID, net.edges$eID),], data.frame(path.ID = edges$path.ID)) %>% st_drop_geometry()) }

  if(output == "raw"){
    return(cbind(net.edges[match(edges$eID, net.edges$eID),], data.frame(path.ID = edges$path.ID))) }

  if(output == "dissolve"){
    edges = cbind(net.edges[match(edges$eID, net.edges$eID),], data.frame(path.ID = edges$path.ID))
    edges = sfheaders::sf_to_df(edges, fill = T)
    # Group geometries by edges.ID
    edges = sfheaders::sf_linestring(edges, x = "x", y = "y", z = "z", linestring_id = "path.ID", keep = T)
    sf::st_crs(edges) = sf::st_crs(net.edges) # Define CRS
    return(edges)
  }
}

##_____________________________________________________________________________
#' Compute 3D Network Catchments for Multiple Origins
#'
#' Identifies the network catchment (reachable nodes or edges) within a specified cost threshold for each origin point in a 3D network, with support for parallel computation.
#'
#' @param o An \code{sf} object with 3D \code{POINT} geometry representing origin points, with an optional \code{character} column \code{nID} (nearest network node).
#' @param net A \code{list} containing \code{nodes} (3D \code{sf} \code{POINT} objects) and \code{edges} (3D \code{sf} \code{LINESTRING} objects), typically from \code{\link[GISnetwork3D]{sf_to_net3d}}.
#' @param mode A \code{character} value specifying directionality: \code{"out"} (origin to destination), \code{"in"} (destination to origin), or \code{"both"} (round trip, minimising sum of out and in costs).
#' @param weight A \code{character} value specifying the edge attribute column in \code{net$edges} for travel cost (e.g., \code{"time"}).
#' @param threshold A non-negative \code{numeric} value specifying the maximum cumulative cost (e.g., travel time) for the catchment.
#' @param searchDist A non-negative \code{numeric} value specifying the maximum Euclidean distance to search for reachable nodes (spatial cutoff).
#' @param output A \code{character} value specifying the output format: \code{"df"} (node IDs), \code{"point"} (reachable nodes as \code{sf} points), \code{"ConvexHull"} (convex hull polygons), or \code{"line"} (reachable edges as \code{sf} linestrings) (default: \code{"ConvexHull"}).
#' @param ncores A positive \code{integer} value specifying the number of parallel workers (default: 4).
#' @param nrows A non-negative \code{integer} value specifying the number of rows for spatial tiling (default: 2).
#' @param ncols A non-negative \code{integer} value specifying the number of columns for spatial tiling (default: 2).
#'
#' @details
#' For each origin point in \code{o}, the function computes the set of nodes or edges reachable within the \code{threshold} cost along a 3D network, constrained by a \code{searchDist} Euclidean cutoff. If \code{nID} is missing in \code{o}, the nearest node is assigned using \code{\link[GISnetwork3D]{nearestNode3d}}. The \code{mode} parameter determines directionality, with \code{"both"} minimising the sum of out-and-in costs. Outputs include node IDs, 3D point geometries, convex hull polygons, or 3D linestring geometries. Origins are divided into \code{nrows} x \code{ncols} spatial tiles, and reachable nodes within \code{buf} distance of each tile are identified in parallel using \code{ncores} workers.
#'
#' @return Depends on \code{output}:
#' \describe{
#'   \item{\code{df}}{A \code{data.frame} with columns: \code{oID} (\code{integer}, \code{numeric}, or \code{character}, origin identifier) and reachable node IDs.}
#'   \item{\code{point}}{An \code{sf} object with 3D \code{POINT} geometries of reachable nodes, including \code{oID}.}
#'   \item{\code{ConvexHull}}{An \code{sf} object with 3D \code{POLYGON} geometries (convex hulls of catchments), including \code{oID}.}
#'   \item{\code{line}}{An \code{sf} object with 3D \code{LINESTRING} or \code{MULTILINESTRING} geometries of reachable edges, including \code{oID}.}
#' }
#'
#' @examples
#' \dontrun{
#'
#' #' # Create origin and destination points
#' o = st_sf(geometry = st_sfc(
#'   st_point(c(0, 0, 1)),
#'   st_point(c(3, 4, 3.5))
#' ))
#'
#' # Create a line network
#' l = st_sf(geometry = st_sfc(
#'   st_linestring(matrix(c(0,0,0, 50, 50, 1), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(50, 50, 1, 100, 100, 2), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(0,0,0, 0,100,3), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(0,100,3, 0,150,3.5), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(0,150,3.5, 10,170,4), ncol=3, byrow=TRUE)),
#'   st_linestring(matrix(c(100,100,4, 130,180,5), ncol=3, byrow=TRUE))
#' ))
#' l = blendPT3d(o, l) # Blend points to the line network
#' l = splitNearestV3d(o, l) # Split network by the nearest vertices for the origins and destinations
#' l$direction = rep(0, nrow(l)) # Assumes bi-directional
#'
#' # Create 3D network with Tobler's hiking function for cost
#' net = sf_to_net3d.full(l, direction = "direction", CF = "Tobler")
#'
#' # Identify network nodes for o and d
#' o = nearestNode3d(net$nodes, o)
#'
#' # Compute accessibility
#' result = catchment3d(o, net, mode = "out", weight = "time", threshold = 200, searchDist = 300)
#' print(result)
#' }
#'
#' @export
catchment3d = function(o, net, mode, weight, threshold, searchDist, output = "ConvexHull", ncores = 4, nrows = 2, ncols = 2){

  # --- Input validation ---

  if (!inherits(o, "sf") || !"sfc_POINT" %in% class(sf::st_geometry(o))) {
    stop("o must be an <sf> object with POINT geometry.") }

  if (!inherits(net, "list") || !"nodes" %in% names(net) || !"edges" %in% names(net)) {
    stop("net must be a list with 'nodes' and 'edges'.")}

  if (!mode %in% c("out", "in", "both")) {stop("mode must be one of 'out', 'in', or 'both'.")}

  if (!is.character(weight) || length(weight) != 1) {
    stop("weight must be a character indicating the column name in net's edges for weights.")}

  if (!is.numeric(threshold) || length(threshold) != 1 || threshold < 0) {
    stop("threshold must be a non-negative <numeric> value.")}

  if (!is.numeric(searchDist) || length(searchDist) != 1 || searchDist < 0) {
    stop("searchDist must be a non-negative numeric value.")}

  if (!output %in% c("df", "point", "ConvexHull", "line")) {stop("output must be one of 'df', 'point', 'ConvexHull' or 'line'.")}

  if (!is.numeric(ncores) || length(ncores) != 1 || ncores <= 0 || ncores %% 1 != 0) {
    stop("ncores must be a positive integer value.")}

  if (!is.numeric(ncols) || length(ncols) != 1 || ncols <= 0 || ncols %% 1 != 0) {
    stop("ncols must be a positive integer value.")}

  if (!is.numeric(nrows) || length(nrows) != 1 || ncols <= 0 || ncols %% 1 != 0) {
    stop("ncols must be a positive integer value.")}


  # --- Prepare the origin point data ---

  # If "nID" is missing from `o`, assign nearest network nodes
  if (!"nID" %in% colnames(o)) {o = nearestNode3d(net$nodes, o)}

  # Generate default IDs for origins
  oID = 1:nrow(o)

  # Split Data
  inputData = tile.OD(o = o, d = net$nodes, buf = searchDist, nrows = nrows, ncols = ncols)
  SubsetNET = purrr::map(inputData$d, ~ subsetNet3d(net = net, v = .x$nID, exact = F, output = "sf"))

  # Compute Catchments
  future::plan("future::multisession", workers = ncores) # Activate cores
  result = furrr::future_pmap(list(inputData$o, inputData$d, SubsetNET), catchment3d.sTask, mode = mode, weight = weight, threshold = threshold, output = output,
                              .options = furrr::furrr_options(seed = TRUE, packages = c("igraph", "sf", "sfheaders")) ) %>% purrr::reduce(rbind)
  future::plan("future::sequential")

  # Export Result
  result = result[order(result$id),]
  names(result)[1] = "oID"
  return(result)   # <---- OUTPUT
  }


