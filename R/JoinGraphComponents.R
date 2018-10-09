#' In case the graph is not connected, join the components
#' This function updates the original adjacency matrix and returns a new object.
#' @param root_cell the root cell of the lineage tree
#' @param adjacency_matrix the graph embedding of the cells
#' @param flat_distances the flattened embedding of the cells
#' @param n_components the number of components
#' @param component_members a list of vectors containing the cells in each component
#' @return \code{adjacency_matrix} such that new edges between disconnected components have length 2
#' @examples my_matrix <- JoinGraphComponents(root_cell = RSoptSC::GuoPtime$Values$root_cell0,
#'     adjacency_matrix = RSoptSC::GuoPtime$Values$W_graph1,
#'     flat_distances = RSoptSC::GuoPtime$Values$low_dis,
#'     n_components = RSoptSC::GuoPtime$Values$nComponents,
#'     component_members = RSoptSC::GuoPtime$Values$members)
JoinGraphComponents <- function(root_cell, adjacency_matrix, flat_distances, n_components, component_members){
  visited_cells <- matrix(nrow = 0, ncol = 2)
  unvisited_cells <- matrix(nrow = 0, ncol = 2)
  #
  ## find the max-distance pairs within each graph component
  #
  for(i in 1:n_components){
    imembers <- as.vector(unlist(component_members[[i]]))
    component_flat_distances <- flat_distances[imembers, imembers]
    xy <- which(component_flat_distances == max(component_flat_distances), arr.ind = TRUE)
    unvisited_cells <- rbind(unvisited_cells,
                             imembers[which(component_flat_distances == max(component_flat_distances),
                                            arr.ind = TRUE)[1,]])
  }
  #
  ## find the start cell of the root component
  #
  for(i in 1:n_components){
    if(is.element(root_cell, component_members[[i]])){
      check_start <- unvisited_cells[i,]
      start_distance <- flat_distances[root_cell, check_start]
      start_cell <- unvisited_cells[i,which.max(start_distance)]
      non_start_cell <- setdiff(unvisited_cells[i,], start_cell) # need to exclude this from future MST search
      visited_cells <- rbind(visited_cells, unvisited_cells[i,])
      unvisited_cells <- unvisited_cells[-i,]
      break()
    }
  }

  #
  ## add an edge between the start cell and the closest component
  #
  distance_to_next <- matrix(flat_distances[start_cell, as.vector(unvisited_cells)], ncol = 2)
  next_index <- which(distance_to_next == min(distance_to_next), arr.ind = TRUE)
  next_cell <- unvisited_cells[next_index]
  adjacency_matrix[start_cell, next_cell] <- 2
  visited_cells <- rbind(visited_cells, unvisited_cells[next_index[1],])
  unvisited_cells <- unvisited_cells[-next_index[1],]
  vis <- sort(setdiff(as.vector(visited_cells), non_start_cell))
  unvis <- sort(as.vector(unvisited_cells))

  #
  ## get MST on remaining
  #
  while (length(unvis) >0){
    distance_to_next <- flat_distances[vis, unvis]
    new_edge <- which(distance_to_next == min(distance_to_next), arr.ind = TRUE)
    u <- vis[new_edge[1]]
    v<- unvis[new_edge[2]]
    adjacency_matrix[u, v] <- 2
    xy <- which(unvisited_cells == v, arr.ind = TRUE)[1]

    vv <- setdiff(unvisited_cells[xy,], v)
    vis <- sort(c(vis, v, vv))
    unvis <- sort(setdiff(unvis, c(v, vv)))

  }
  adjacency_matrix <- cbind(adjacency_matrix)
}
