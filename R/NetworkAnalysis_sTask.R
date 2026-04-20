
##_____________________________________________________________________________
# -----> This function computes spatial accessibility with a single core.

access3d.pairOD.sTask = function(net, o, d, mode = "out", weight, threshold = 600, avaiW = NULL){
  # o and d must be nID
  # mode: c("in", "out", "both")
  # weight: A character
  # avaiW: A numeric vector with the same length as the row number of d or NULL

  # --- Main function implementation ---
  # Calculate the ODM (Origin-Destination Matrix) using the ODM3d function
  if(mode == "out"){ODM = GISnetwork3D::ODM3d(net = net, o = o, d = d, mode = "out", weight = weight)}
  if(mode == "in"){ODM = t(GISnetwork3D::ODM3d(net = net, o = d, d = o, mode = "out", weight = weight))}
  if(mode == "both"){ODM =  GISnetwork3D::ODM3d(net = net, o = o, d = d, mode = "out", weight = weight) +
    t(GISnetwork3D::ODM3d(net = net, o = d, d = o, mode = "out", weight = weight))}

  # Assign weights for availability computation if `avaiW` is NULL
  if (is.null(avaiW)) { avaiW = rep(1, length(d)) }  # Default weight is 1 for each destination

  # --- Summarise the ODM ---
  result = data.frame(oID = 1:length(o), oNode = o, dID = NA, dNode = NA, costs = NA)
  result$costs = apply(ODM, 1, FUN = min, na.rm = TRUE)
  result$dID = apply(ODM, 1, which.min)
  result$dNode = d[result$dID]
  avai = data.frame(oID = rep(result$oID, each = length(d)), avai = avaiW[rep(1:length(d), length(o))],
                    costs = c(t(ODM))) %>% subset(costs <= threshold)
  avai = avai %>% dplyr::group_by(oID) %>% dplyr::summarise(avai = sum(avai))
  result = result %>% dplyr::left_join(avai, by = "oID")
  result$avai[is.na(result$avai)] = 0
  # Return a data frame summarising the results
  return(result)  # <---- OUTPUT: data.frame
}

##_____________________________________________________________________________
# Compute Spatial Accessibility and Edge Attribute Summaries on a 3D Network with Single core

access3d.pairOD.summary.sTask = function(net, o, d, mode = "out", weight, threshold, avaiW = NULL,
                                         geom = TRUE,
                                         sumV = NULL,
                                         w.sumV = NULL, w.sumW = NULL,
                                         p.sumV = NULL, pw.sumV = NULL, p.norm = NULL, both.split = FALSE){

  # Extract nID for o and d
  if(inherits(o, "sf")){o = o$nID}
  if(inherits(d, "sf")){d = d$nID}


  # --- Compute access ---
  ODpair = GISnetwork3D::access3d.pairOD.sTask(net = net, o = o, d = d, mode = mode, weight = weight, threshold = threshold, avaiW = avaiW)

  # --- Summarise Path Attributes ---
  path.summary = GISnetwork3D::path3d.pairOD.summary.sTask(net = net, oNode = ODpair$oNode, dNode = ODpair$dNode, mode = mode, weight = weight, oID = ODpair$oID, geom = geom,
                                                           sumV = sumV,
                                                           w.sumV = w.sumV, w.sumW = w.sumW,
                                                           p.sumV = p.sumV, pw.sumV = pw.sumV, p.norm = p.norm, both.split = both.split)

  return(dplyr::left_join(path.summary, ODpair, by = "oID")[unique(c(names(ODpair), names(path.summary)))]) # <--- Output: sf or data.frame
}


##_____________________________________________________________________________
# ----> This function finds the one-to-many least-costs routes

get_paths_o2m = function(graph, oNode, dNode, mode = "out", weight, E.seq = FALSE){

  # Find shortest paths
  paths = igraph::shortest_paths(graph, from = oNode, to = dNode, mode = mode, weights = weight, output = "epath")$epath
  path_lengths = lengths(paths) # Get path length
  dNode_rep = rep(dNode, path_lengths) # Get the dNode for each pair of OD
  oNode_rep = rep(oNode, length(dNode_rep)) # Get the oNode for each pair of OD
  path_ids = unlist(paths) # Get the edge sequence ID

  # Remove those cannot be reached
  oNode = oNode[oNode %in% oNode_rep]
  dNode = dNode[dNode %in% dNode_rep]
  path_lengths = path_lengths[path_lengths!=0]

  if(length(path_lengths)>0){
    # Reverse the edge sequence if inbound movement is specified
    if(mode == "in"){
      # Define sequence by original sequence of dNode and reversed order of the edge sequence
      seq_order = order(rep(seq_along(dNode), path_lengths), unlist(Map(seq, path_lengths, 1, -1)))
      path_ids = path_ids[seq_order] # Reorder the path_id
      dNode_rep = dNode_rep[seq_order] # Reorder the dNode
      oNode_rep = oNode_rep[seq_order] } # Reorder the oNode

    # Define if export eID as edge name or the sequence of edges
    eID = if (E.seq) path_ids else igraph::edge_attr(graph)$eID[path_ids]

    return(list(oNode = oNode_rep, dNode = dNode_rep, eID = eID)) # <---- Output: List
  } else {
    # Return empty if not path can be sorted
    return(list(oNode = NULL, dNode = NULL, eID = NULL)) # <---- Output: List
  }
}


##_____________________________________________________________________________
# ----> This function finds the shortest path for each pair of OD with optimised computation
# This is very efficient if there are many duplicates in dNode

path3d.pairOD.get_eID.opt = function(net, oNode, dNode, mode = "out", weight, output = "list", oID = NULL, E.seq = FALSE){

  # --- Network preparation ---
  if (inherits(net, "list")) {
    # If `net` is a list, extract the weight column and convert to an igraph object
    weight = net$edges[[weight]] # Extract the weight
    # Convert to igraph
    net = igraph::graph_from_data_frame(d = net$edges, directed = T, vertices = net$nodes %>% sf::st_drop_geometry())
  } else {
    # If `net` is already an igraph object, extract the weight column
    weight = igraph::as_data_frame(net)[[weight]] # Extract the weight
  }

  # Define `oID` if not provided
  if (is.null(oID)) { oID = 1:length(oNode)} # Generate default IDs for origins

  # --- Optimise OD pairs ---
  uniqueODpairs = data.frame(oNode = oNode, dNode = dNode) %>% dplyr::distinct()
  uniqueODpairs = split(uniqueODpairs, uniqueODpairs$dNode)

  # --- Find Shortest Paths ---
  if(mode == "out"){paths = purrr::map(uniqueODpairs, ~ GISnetwork3D::get_paths_o2m(graph = net,
                                                                                    oNode = unique(.x$dNode), dNode = .x$oNode,
                                                                                    mode = "in", weight = weight, E.seq = E.seq)) }

  if(mode == "in"){paths = purrr::map(uniqueODpairs, ~ GISnetwork3D::get_paths_o2m(graph = net,
                                                                                   oNode = unique(.x$dNode), dNode = .x$oNode,
                                                                                   mode = "out", weight = weight, E.seq = E.seq)) }

  paths = data.frame(oNode = unlist(purrr::map(paths, ~.x$dNode)), dNode = unlist(purrr::map(paths, ~.x$oNode)), eID = unlist(purrr::map(paths, ~.x$eID)))

  # Prepare Result
  if(nrow(paths)>0){
    result = data.frame(oID = oID, oNode = oNode, dNode = dNode) %>% dplyr::left_join(paths, by = c("oNode", "dNode"), relationship = "many-to-many")
  } else {
    result = data.frame(oID = oID, oNode = oNode, dNode = dNode, eID = NA)
  }

  if(output == "df"){return(result[c("oID", "eID")])} # <--- Output: data.frame
  if(output == "list"){return( split(result$eID, result$oID)[oID] )} # <--- Output: list
}


##_____________________________________________________________________________
# ----> This function summarises Edge Attributes Along Shortest Paths in a 3D Network using single core
path3d.pairOD.summary.sTask = function(net, oNode, dNode, mode = "out", weight, oID = NULL, geom = TRUE,
                                       sumV = NULL,
                                       w.sumV = NULL, w.sumW = NULL,
                                       p.sumV = NULL, pw.sumV = NULL, p.norm = NULL, both.split = FALSE){

  # Assign oID if not specified
  if(is.null(oID)) oID = 1:length(oNode)
  result = data.frame(oID = oID)

  # --- Get the edge sequences of the shortest paths ---
  if(mode == "both"){
    paths.out = GISnetwork3D::path3d.pairOD.get_eID.opt(net = net, oNode = oNode, dNode = dNode, mode = "out", weight = weight, output = "df", oID = oID, E.seq = TRUE)
    paths.in = GISnetwork3D::path3d.pairOD.get_eID.opt(net = net, oNode = oNode, dNode = dNode, mode = "in", weight = weight, output = "df", oID = oID, E.seq = TRUE)
    paths.out$Lseq = 1
    paths.in$Lseq = 2
    paths.out$Sseq = 1:nrow(paths.out)
    paths.in$Sseq = (1:nrow(paths.in)) + nrow(paths.out)
    paths = rbind(paths.out, paths.in)

    if(isFALSE(both.split)){ paths = paths[order(paths$oID, paths$Lseq, paths$Sseq),c("oID", "eID")] }
    if(isTRUE(both.split)){ paths = paths[order(paths$oID, paths$Lseq, paths$Sseq),c("oID", "eID", "Lseq")]} }
  if(mode != "both"){
    paths = GISnetwork3D::path3d.pairOD.get_eID.opt(net = net, oNode = oNode, dNode = dNode, mode = mode, weight = weight, output = "df", oID = oID, E.seq = TRUE)
  }

  # --->>> Path Summary <<<---
  # --- Simple Sum ---
  if(!is.null(sumV)){
    if("Lseq" %in% names(paths) == F){
      # Not Splitting
      tmp = cbind(oID = paths$oID, sf::st_drop_geometry(net$edges)[paths$eID, sumV, drop = F])
      tmp = aggregate(. ~ oID, data = tmp, FUN = sum)
      result = result %>% dplyr::left_join(tmp, by = "oID")}

    if("Lseq" %in% names(paths)){
      # Splitting
      tmp.out = cbind(paths, sf::st_drop_geometry(net$edges)[paths$eID, sumV, drop = F]) %>% base::subset(Lseq == 1) %>% dplyr::select(-eID, -Lseq)
      tmp.out = aggregate(. ~ oID, data = tmp.out, FUN = sum)
      names(tmp.out)[2:length(tmp.out)] = paste0(names(tmp.out)[2:length(tmp.out)], ".out")

      tmp.in = cbind(paths, sf::st_drop_geometry(net$edges)[paths$eID, sumV, drop = F]) %>% base::subset(Lseq == 2) %>% dplyr::select(-eID, -Lseq)
      tmp.in = aggregate(. ~ oID, data = tmp.in, FUN = sum)
      names(tmp.in)[2:length(tmp.in)] = paste0(names(tmp.in)[2:length(tmp.in)], ".return")

      result = result %>% dplyr::left_join(tmp.out, by = "oID") %>% dplyr::left_join(tmp.in, by = "oID")} }

  # --- Weighted Sum ---
  if(!is.null(w.sumV)){
    # Not Splitting
    if("Lseq" %in% names(paths) == F){
      tmp = sf::st_drop_geometry(net$edges)[paths$eID, w.sumV, drop = F]
      names(tmp) = paste("w.", names(tmp), sep = "")
      tmp = tmp * sf::st_drop_geometry(net$edges)[paths$eID, w.sumW]
      tmp = aggregate(. ~ oID, data = cbind(oID = paths$oID, tmp), FUN = sum)
      result = result %>% dplyr::left_join(tmp, by = "oID")}

    if("Lseq" %in% names(paths)){
      # Splitting
      tmp.out = sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 1], w.sumV, drop = F]
      names(tmp.out) = paste0("w.", names(tmp.out), ".out")
      tmp.out = tmp.out * sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 1], w.sumW]
      tmp.out = aggregate(. ~ oID, data = cbind(oID = paths$oID[paths$Lseq == 1], tmp.out), FUN = sum)

      tmp.in = sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 2], w.sumV, drop = F]
      names(tmp.in) = paste0("w.", names(tmp.in), ".return")
      tmp.in = tmp.in * sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 2], w.sumW]
      tmp.in = aggregate(. ~ oID, data = cbind(oID = paths$oID[paths$Lseq == 2], tmp.in), FUN = sum)

      result = result %>% dplyr::left_join(tmp.out, by = "oID") %>% dplyr::left_join(tmp.in, by = "oID")} }

  # --- Normalised sums ---
  if(!is.null(p.sumV) & !is.null(p.norm)){
    if("Lseq" %in% names(paths) == F){
      # Not Splitting
      tmp = sf::st_drop_geometry(net$edges)[paths$eID, c(p.sumV, p.norm), drop = F]
      tmp = aggregate(. ~ oID, data = cbind(oID = paths$oID, tmp), FUN = sum)
      tmp[p.sumV] = tmp[p.sumV] / unlist(tmp[p.norm])
      names(tmp)[(names(tmp) %in% p.sumV)] = paste("p.", names(tmp)[(names(tmp) %in% p.sumV)], sep = "")
      result = result %>% dplyr::left_join(tmp[names(tmp) != p.norm], by = "oID")}

    if("Lseq" %in% names(paths)){
      # Splitting
      tmp.out = sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 1], c(p.sumV, p.norm), drop = F]
      tmp.out = aggregate(. ~ oID, data = cbind(oID = paths$oID[paths$Lseq == 1], tmp.out), FUN = sum)
      tmp.out[p.sumV] = tmp.out[p.sumV] / unlist(tmp.out[p.norm])
      names(tmp.out)[(names(tmp.out) %in% p.sumV)] = paste0("p.", names(tmp.out)[(names(tmp.out) %in% p.sumV)], ".out")

      tmp.in = sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 2], c(p.sumV, p.norm), drop = F]
      tmp.in = aggregate(. ~ oID, data = cbind(oID = paths$oID[paths$Lseq == 2], tmp.in), FUN = sum)
      tmp.in[p.sumV] = tmp.in[p.sumV] / unlist(tmp.in[p.norm])
      names(tmp.in)[(names(tmp.in) %in% p.sumV)] = paste0("p.", names(tmp.in)[(names(tmp.in) %in% p.sumV)], ".return")

      result = result %>% dplyr::left_join(tmp.out[names(tmp.out) != p.norm], by = "oID") %>%
        dplyr::left_join(tmp.in[names(tmp.in) != p.norm], by = "oID")} }

  # --- Weighted Normalised sums ---
  if(!is.null(pw.sumV) & !is.null(p.norm)){
    if("Lseq" %in% names(paths) == F){
      # Not Splitting
      tmp = sf::st_drop_geometry(net$edges)[paths$eID, c(pw.sumV, p.norm), drop = F]
      tmp[pw.sumV] = tmp[pw.sumV] * unlist(tmp[p.norm])
      tmp = aggregate(. ~ oID, data = cbind(oID = paths$oID, tmp), FUN = sum)
      tmp[pw.sumV] = tmp[pw.sumV] / unlist(tmp[p.norm])
      names(tmp)[(names(tmp) %in% pw.sumV)] = paste("pW.", names(tmp)[(names(tmp) %in% pw.sumV)], sep = "")
      result = result %>% dplyr::left_join(tmp[names(tmp) != p.norm], by = "oID")}

    if("Lseq" %in% names(paths)){
      # Splitting
      tmp.out = sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 1], c(pw.sumV, p.norm), drop = F]
      tmp.out[pw.sumV] = tmp.out[pw.sumV] * unlist(tmp.out[p.norm])
      tmp.out = aggregate(. ~ oID, data = cbind(oID = paths$oID[paths$Lseq == 1], tmp.out), FUN = sum)
      tmp.out[pw.sumV] = tmp.out[pw.sumV] / unlist(tmp.out[p.norm])
      names(tmp.out)[(names(tmp.out) %in% pw.sumV)] = paste0("pW.", names(tmp.out)[(names(tmp.out) %in% pw.sumV)], ".out")

      tmp.in = sf::st_drop_geometry(net$edges)[paths$eID[paths$Lseq == 2], c(pw.sumV, p.norm), drop = F]
      tmp.in[pw.sumV] = tmp.in[pw.sumV] * unlist(tmp.in[p.norm])
      tmp.in = aggregate(. ~ oID, data = cbind(oID = paths$oID[paths$Lseq == 2], tmp.in), FUN = sum)
      tmp.in[pw.sumV] = tmp.in[pw.sumV] / unlist(tmp.in[p.norm])
      names(tmp.in)[(names(tmp.in) %in% pw.sumV)] = paste0("pW.", names(tmp.in)[(names(tmp.in) %in% pw.sumV)], ".return")
      result = result %>% dplyr::left_join(tmp.out[names(tmp.out) != p.norm], by = "oID")  %>%
        dplyr::left_join(tmp.in[names(tmp.in) != p.norm], by = "oID")} }

  # Replace NA value by 0
  result[is.na(result)] = 0

  # Find Geometry
  if(isTRUE(geom)){
    tmp = cbind(oID = paths$oID, net$edges[paths$eID,"eID"])
    tmp = sfheaders::sf_to_df(tmp, fill = T)
    tmp = tmp %>% sfheaders::sf_linestring(x = "x", y = "y", z = "z", linestring_id = "oID", keep = T) # Group geometry by oID
    sf::st_crs(tmp) = sf::st_crs(net$edges)
    result = cbind(dplyr::select(net$edges[(nrow(net$edges)+1):(nrow(net$edges)+nrow(result)),]), result)
    sf::st_geometry(result) = sf::st_geometry(tmp)[match(result$oID, tmp$oID)] }

  # Export Data
  return(result) # <--- Output: sf or data.frame
}



##_____________________________________________________________________________
# This is an inner function for the access3d.multitasks
access3d.multitasks.sTask = function(startPT, taskPT, net, weight, export.IDs = FALSE, export.nIDs = FALSE, edges.output = NULL){

  # Assign Unique Names for POIs
  startPT = startPT[c("ID", "nID")]
  taskPT = taskPT[c("ID", "nID", "cluster")]
  startPT$uniqueNameID = paste("S", 1:nrow(startPT), sep = "")
  taskPT$uniqueNameID = paste("T", 1:nrow(taskPT), sep = "")

  # Combine all POIs
  allPT = rbind(startPT, taskPT[c(c("ID", "nID", "uniqueNameID"))])

  # --- Compute the full origin-destination matrix (all points to all points) ---
  OD.matrix = GISnetwork3D::ODM3d(net, o = allPT$nID, d = allPT$nID, mode = "out", weight = weight)
  rownames(OD.matrix) = c(startPT$uniqueNameID, taskPT$uniqueNameID)
  colnames(OD.matrix) = c(startPT$uniqueNameID, taskPT$uniqueNameID)
  diag(OD.matrix) = NA

  # --- Pre-compute least-cost links between tasks, stratified by cluster ---
  # Task-to-task submatrix
  OD.matrix.tasks = OD.matrix[taskPT$uniqueNameID, taskPT$uniqueNameID, drop = FALSE] # Index ODM only containing taskPT
  OD.matrix.tasks = data.frame(costs = c(t(OD.matrix.tasks)),
                               ODO = rep(taskPT$uniqueNameID, each = nrow(taskPT)),
                               ODD = rep(taskPT$uniqueNameID, nrow(taskPT)),
                               Ocluster = rep(taskPT$cluster, each = nrow(taskPT)),
                               Dcluster = rep(taskPT$cluster, nrow(taskPT)) )

  # For each origin task (ODO) and destination cluster, keep only least-cost pair
  OD.matrix.tasks = OD.matrix.tasks %>%
    dplyr::group_by(ODO, Dcluster) %>%
    dplyr::slice_min(costs, with_ties = FALSE, na_rm = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(ODO, Dcluster) %>%
    dplyr::filter(Ocluster != Dcluster)

  # --- Containers for results ---
  output = data.frame(start.ID = startPT$ID, costs = NA)
  cluster = matrix(nrow = nrow(startPT), ncol = length(unique(taskPT$cluster)))
  NameID = matrix(nrow = nrow(startPT), ncol = length(unique(taskPT$cluster)))

  # --- First-mile costs (startPT -> first taskPT) ---
  OD.matrix.firstMile = OD.matrix[startPT$uniqueNameID, taskPT$uniqueNameID, drop = FALSE] # Extract ODM
  output$costs = apply(OD.matrix.firstMile, 1, FUN = min, na.rm = TRUE) # Compute the least-costs
  minIndex = apply(OD.matrix.firstMile, 1, which.min) # Find the index of taskPT with the least costs

  cluster[, 1] = taskPT$cluster[minIndex]
  NameID[, 1] = taskPT$uniqueNameID[minIndex]
  rm(OD.matrix.firstMile)

  # --- Intermediate legs between task clusters ---
  for(i in 2:ncol(cluster)){
    # Start from previously chosen stop (NameID[, i-1]) for each start
    masterTable = data.frame(startNameID = startPT$uniqueNameID) # Get the Original startPT NameID
    masterTable$ODO = NameID[, i-1] # Get the start location's NameID

    # Join with pre-computed task-to-task least-cost table
    masterTable = masterTable %>% dplyr::left_join(OD.matrix.tasks, relationship = "many-to-many", by = "ODO")

    # Remove clusters already visited for each start point
    masterTable = masterTable %>%
      dplyr::left_join( data.frame(startNameID = rep(startPT$uniqueNameID, length(1:(i-1))), Dcluster = c(cluster[,1:(i-1)]), rmv = 1),
                        by = c("startNameID", "Dcluster"), relationship = "many-to-many")
    masterTable = masterTable[is.na(masterTable$rmv),]

    # Identify the least-costs end stops
    masterTable = masterTable %>%
      dplyr::group_by(startNameID) %>%
      dplyr::slice_min(costs, with_ties = FALSE, na_rm = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(startNameID, Dcluster)

    # Reorder the result
    masterTable = data.frame(startNameID = startPT$uniqueNameID) %>% dplyr::left_join(masterTable, by = "startNameID")

    # Add results to the final product
    output$costs = output$costs + masterTable$costs
    cluster[,i] = masterTable$Dcluster
    NameID[,i] = masterTable$ODD

  }

  # --- Last-mile costs: from final task stop back to startPT ---
  output$costs = output$costs + diag(OD.matrix[NameID[, i], startPT$uniqueNameID, drop = F])

  # --- Optional: export stop IDs ---
  if(export.IDs == T){
    # Extract the unique ID of stops
    output.ID = NameID
    output.ID[] = taskPT$ID[match(NameID, taskPT$uniqueNameID)]
    colnames(output.ID) = paste("ID.stop_", 1:length(taskPT$cluster %>% unique), sep = "") # Define column names
    output = output %>% cbind(output.ID) # Join to the output
  }

  # --- Optional: export node IDs (nID) ---
  if(export.nIDs == T){
    # Extract the node ID (nID) of stops
    output.nID = NameID
    output.nID[] = taskPT$nID[match(NameID, taskPT$uniqueNameID)]
    output.nID = (nID.start = startPT$nID) %>% cbind(output.nID) %>% cbind(nID.end = startPT$nID)
    colnames(output.nID) = c("nID.start", paste("nID.stop_", 1:length(taskPT$cluster %>% unique), sep = ""), "nID.end") # Define column names
    output = output %>% cbind(output.nID) }

  # --- Optional: export edges ---
  if(is.null(edges.output) == F){

    # --- Get the edge index ---

    # Build unique (fromNode, toNode) pairs from nID matrix
    unique.paths = c(output.nID[,seq(1, ncol(output.nID)-1)])
    unique.paths = data.frame(fromNode = unique.paths, toNode = c(output.nID[,seq(2, ncol(output.nID))]))
    unique.paths$pathID = paste(unique.paths$fromNode, unique.paths$toNode, sep = "->")
    unique.paths = unique.paths %>% unique()

    # Get eID for each unique path segment
    paths = GISnetwork3D::path3d.pairOD.get_eID.opt(net = net, oNode = unique.paths$fromNode, dNode = unique.paths$toNode,
                                                    mode = "out", weight = weight, output = "df", oID = unique.paths$pathID)

    names(paths) = c("pathID", "eID")

    # Map back to each start.ID and segment
    edges = matrix(nrow = nrow(startPT), ncol = ncol(output.nID)-1)
    for(i in 1:ncol(edges)){
      edges[,i] = paste( output.nID[,i], output.nID[,i+1], sep = "->" )
    }
    edges = data.frame(start.ID = rep(startPT$ID, each = ncol(edges)), pathID = c(t(edges)))
    edges = edges %>% dplyr::left_join(paths, relationship = "many-to-many", by = "pathID")
    edges = edges[is.na(edges$eID) == F,]

    # --- Export raw edges as data.frame ---
    if(edges.output == "raw.df"){output = list(access = output, eID = edges)}

    # --- Export raw edges as sf ---
    if(edges.output == "raw"){
      edges = cbind(start.ID = edges$start.ID, net$edges[match(edges$eID, net$edges$eID),])
      output = list(access = output, eID = edges)}

    # --- Export Result with Dissolved Geometry ---
    if(edges.output == "dissolved"){
      edges = cbind(start.ID = edges$start.ID, net$edges[match(edges$eID, net$edges$eID),"eID"])
      edges = sfheaders::sf_to_df(edges, fill = T)
      edges = edges %>% sfheaders::sf_linestring(x = "x", y = "y", z = "z", linestring_id = "start.ID", keep = T) # Group geometry by oID
      sf::st_crs(edges) = sf::st_crs(net$edges)

      output = cbind(dplyr::select(net$edges[(nrow(net$edges)+1):(nrow(net$edges)+nrow(output)),]), output)
      sf::st_geometry(output) = sf::st_geometry(edges)[match(output$start.ID, edges$start.ID)] }
  }

  return(output)
}

##_____________________________________________________________________________
# ----> This function summarises Edge Attributes Along Multipurpose Shortest Paths in a 3D Network using single core

access3d.multitasks.summary.sTask = function(startPT, taskPT, net, geom = FALSE, weight,
                                             sumV = NULL,
                                             w.sumV = NULL, w.sumW = NULL,
                                             p.sumV = NULL, pw.sumV = NULL, p.norm = NULL){

  # --- Compute access and get the edges travelled ---
  result = access3d.multitasks.sTask(startPT = startPT, taskPT = taskPT, net = net, weight = weight, export.IDs = TRUE, export.nIDs = TRUE, edges.output = "raw.df")

  # --->>> Path Summary <<<---
  # --- Simple Sum ---
  if(!is.null(sumV)){
    tmp = cbind(start.ID = result$eID$start.ID, sf::st_drop_geometry(net$edges)[match(result$eID$eID, net$edges$eID), sumV, drop = F])
    tmp = aggregate(. ~ start.ID, data = tmp, FUN = sum)
    result$access = result$access %>% dplyr::left_join(tmp, by = "start.ID") }

  # --- Weighted Sum ---
  if(!is.null(w.sumV)){
    tmp = sf::st_drop_geometry(net$edges)[match(result$eID$eID, net$edges$eID), w.sumV, drop = F]
    names(tmp) = paste("w.", names(tmp), sep = "")
    tmp = tmp * sf::st_drop_geometry(net$edges)[match(result$eID$eID, net$edges$eID), w.sumW]
    tmp = aggregate(. ~ start.ID, data = cbind(start.ID = result$eID$start.ID, tmp), FUN = sum)
    result$access = result$access %>% dplyr::left_join(tmp, by = "start.ID") }

  # --- Normalised sums ---
  if(!is.null(p.sumV) & !is.null(p.norm)){
    tmp = sf::st_drop_geometry(net$edges)[match(result$eID$eID, net$edges$eID), c(p.sumV, p.norm), drop = F]
    tmp = aggregate(. ~ start.ID, data = cbind(start.ID = result$eID$start.ID, tmp), FUN = sum)
    tmp[p.sumV] = tmp[p.sumV] / unlist(tmp[p.norm])
    names(tmp)[(names(tmp) %in% p.sumV)] = paste("p.", names(tmp)[(names(tmp) %in% p.sumV)], sep = "")
    result$access = result$access %>% dplyr::left_join(tmp[names(tmp) != p.norm], by = "start.ID") }


  # --- Weighted Normalised sums ---
  if(!is.null(pw.sumV) & !is.null(p.norm)){
    tmp = sf::st_drop_geometry(net$edges)[match(result$eID$eID, net$edges$eID), c(pw.sumV, p.norm), drop = F]
    tmp[pw.sumV] = tmp[pw.sumV] * unlist(tmp[p.norm])
    tmp = aggregate(. ~ start.ID, data = cbind(start.ID = result$eID$start.ID, tmp), FUN = sum)
    tmp[pw.sumV] = tmp[pw.sumV] / unlist(tmp[p.norm])
    names(tmp)[(names(tmp) %in% pw.sumV)] = paste("pW.", names(tmp)[(names(tmp) %in% pw.sumV)], sep = "")
    result$access = result$access %>% dplyr::left_join(tmp[names(tmp) != p.norm], by = "start.ID") }

  # Replace NA value by 0
  result$access[is.na(result$access)] = 0

  # Find Geometry
  if(isTRUE(geom)){
    tmp = cbind(start.ID = result$eID$start.ID, net$edges[match(result$eID$eID, net$edges$eID),"eID"])
    tmp = sfheaders::sf_to_df(tmp, fill = T)
    tmp = tmp %>% sfheaders::sf_linestring(x = "x", y = "y", z = "z", linestring_id = "start.ID", keep = T) # Group geometry by start.ID
    sf::st_crs(tmp) = sf::st_crs(net$edges)
    result$access = cbind(dplyr::select(net$edges[(nrow(net$edges)+1):(nrow(net$edges)+nrow(result$access)),]), result$access)
    sf::st_geometry(result$access) = sf::st_geometry(tmp)[match(result$access$start.ID, tmp$start.ID)] }

  # Export Data
  return(result$access) # <--- Output: sf or data.frame
}


##_____________________________________________________________________________
# This is an inner function to compute the catchment with single core
catchment3d.sTask = function(oNode, dNode, net, mode, weight, threshold, output = "df"){

  # Get the origin and destinations
  oID = oNode$id
  oNode = oNode$nID
  dNode = dNode$nID %>% unique()

  # --- Compute the catchment analysis ---
  # Compute ODM of costs
  ODM = ODM3d(net = net, o = oNode, d = dNode, mode = mode, weight = weight)

  # Find the reachable nodes
  nodesReached = data.frame( id = rep(oID, each = length(dNode)), costs = c(t(ODM)), nID = rep(dNode, length(oNode)) ) %>% subset(costs <= threshold)
  nodesReached = net$nodes[match(nodesReached$nID, net$nodes$nID),] %>% cbind(nodesReached[c("id", "costs")])
  st_geometry(nodesReached) = "geom"

  # Export data
  if(output == "df"){
    # Export data.frame
    return(st_drop_geometry(nodesReached[c("id","nID")]))
  } else if(output == "point"){
    # Export sf POINT
    return(nodesReached[c("id","nID")])
  } else if(output == "ConvexHull"){
    # Export sf POLYGON
    catchPolygon = nodesReached %>% dplyr::group_by(id) %>% summarise(geometry=st_convex_hull(st_union(geom)))
  } else if(output == "line"){
    # Export LINESTRING
    catchLines = split(nodesReached, nodesReached$id)
    catchLines = purrr::map(catchLines,
                            function(x){
                              result = net$edges[net$edges$FROM %in% x$nID & net$edges$TO %in% x$nID,] %>% st_union() %>% st_as_sf()
                              result$id = unique(x$id)
                              return(result)
                            }) %>% purrr::reduce(rbind) }
  }

