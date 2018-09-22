#' Get the predecessor vector for a dominator tree encoded as an undirected tree and root
#'
#' @param minspantree the tree to be directionalized
#' @param root the id of the root
#'
#' @return a vector[] of node ids, where vector[i] is the predecessor of node i
#' @examples predecessors <- GetPredecessors(minspantree = RSoptSC::GuoPtimeFull$mst, root = 6)
#'
GetPredecessors <- function(minspantree, root){
  #browser()
  pathdata <- igraph::get.shortest.paths(graph = minspantree, from = root)
  pathlist <- prelist <- lapply(pathdata$vpath, function(x){
    y = as.vector(x)
    y[(length(y)-1)]
  })
  prelist[[root]] <- 0
  return(unlist(prelist))
}


#' Get a directed graph from a predecessor vector
#'
#' @param predecessors the tree to be directionalized
#'
#' @return a directed igraph object
#'
GetDominatorTree <- function(predecessors){
  edges <- c()
  for(i in 1:length(predecessors)){
    if(predecessors[i] != 0){
      edge <- c(predecessors[i], i)
      edges <- c(edges, edge)
    }
  }
  tree <- igraph::make_graph(edges = edges, directed = TRUE)
}

#' Get a weighted, directed graph from a table of edges and weights
#'
#' @param edge_table a numeric matrix whose rows are edges, col 1 is v1, col2 is v2, col3 is weight
#'
#' @return a directed igraph object
#'
GetDGFromTable <- function(directed_edge_table){
  n_vertices <- max(directed_edge_table[,1:2])
  n_edges <- nrow(directed_edge_table)
  adj_matrix <- matrix(0, n_vertices, n_vertices)
  for(i in 1:n_edges){
    adj_matrix[directed_edge_table[i,1], directed_edge_table[i,2]] <- directed_edge_table[i,3]
  }
  tree <- igraph::graph_from_adjacency_matrix(adjmatrix = adj_matrix,
                                              weighted = TRUE,
                                              mode = 'directed')
}


#' Convert a matlab edge table and predecessor list into a directed weighted edge table
#'
#' @param edge_table a numeric matrix whose rows are undirected edges of a tree,
#'     col 1 is v1, col2 is v2, col3 is weight
#' @param predecessors a vector of tree predecessors such that pred[i] = the predecessor of i
#'
#' @return a weighted, directed edge table
#'
ProcessMatlabDTree <- function(edge_table, predecessors){
  # make a directed edge list
  directed_edges <- cbind(c(predecessors), 1:length(predecessors), rep(0, length(predecessors)))
  # remove the null edge leading to root
  directed_edges <- directed_edges[-(which(directed_edges[,1] == 0, arr.ind = TRUE)),]

  n_edges <- nrow(edge_table)
  for(i in 1:n_edges){
    for(ii in 1:n_edges){
      if(sort(directed_edges[i,1:2])==sort(edge_table[ii, 1:2])){
        directed_edges[i,3] <- edge_table[ii,3]
      }
    }
  }
  return(directed_edges)
}

