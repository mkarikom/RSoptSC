#' Infer the root cell
#'
#' @param use_flat_dist use max path length heuristic on flat embedding
#' @param cluster_order_by cluster rank parameter for measuring cluster/cell rank correlation
#' @param cell_order_by cell-wise pseudotime parameter for measuring cluster/cell rank correlation
#' @param graph_cluster the weighted adjacency matrix for a graph of clusters
#' @param dist_graph a distance matrix of cells embedded in a graph
#' @param dist_flat the manifold embedding of the cells
#' @param cluster_labels the cluster label for each cell
#' @param root_cluster the id of the root cluster on the cluster-cluster graph
#'
#' @return integer index of the root cell
#' @examples root_cell <- FindRootCell(use_flat_dist = TRUE, <some_pca_embeddings>)
#'
#' Use either the primary manifold embedding of cell similarity to find the root cell by cluster/cell rank-correlation
#'  or use the flattened representation of this embedding to find the root cell by maximum separation heuristic
#'  possible values for \code{cluster_order_by} are "predecessor" and "distance"
#'  possible values for \code{cell_order_by} are "index" and "distance"
#'

FindRootCell <- function(use_flat_dist = TRUE,
                         cluster_order_by = "distance",
                         cell_order_by = "distance",
                         graph_cluster = NULL,
                         dist_graph  = NULL,
                         dist_flat,
                         cluster_labels = NULL,
                         root_cluster = NULL){
  if(use_flat_dist){
    browser()
    xy <- which(dist_flat == max(dist_flat), arr.ind = TRUE)
    max_cells <- sort(unique(as.vector(which(dist_flat == max(dist_flat), arr.ind = TRUE))))
    return(max_cells[length(max_cells)])
  } else {
    browser()
    n_clusters <- length(unique(cluster_labels))
    root_cluster_cells <- which(cluster_labels == root_cluster, arr.ind = TRUE)

    tau_score <- matrix(0, 1, length(root_cluster_cells))
    avg_ptime_to_clusters <- matrix(0, length(root_cluster_cells), n_clusters)
    for(i in 1:length(root_cluster_cells)){
      for(j in 1:n_clusters){
        avg_ptime_to_clusters[i,j] <- GetCellMetric(cell_order_by, i, j, cluster_labels, dist_graph)
      }
      avg_ptime_to_clusters[i,] <- avg_ptime_to_clusters[i,] / max(avg_ptime_to_clusters[i,])
      tau_score[i] <- GetCorrelation(cluster_order_by, avg_ptime_to_clusters[i,], graph_cluster, root_cluster)
    }
    winner <- which(tau_score == min(tau_score))
    if(length(winner > 1)){
      winner = winner[1]
    }
    root_cell <- root_cluster_cells[winner]
    return(root_cell)
  }
}

GetPredecessors <- function(minspantree, root){
  pathdata <- igraph::get.shortest.paths(graph = minspantree, from = root)
  pathlist <- prelist <- lapply(pathdata$vpath, function(x){
    y = as.vector(x)
    y[(length(y)-1)]
  })
  prelist[[root]] <- 0
  return(unlist(prelist))
}

GetCorrelation <- function(cluster_order_by, cell_metric, graph_cluster, root_cluster){
  if(cluster_order_by == "distance"){
    path_lengths <- igraph::shortest.paths(graph = minspantree, v = root_cluster) + 1
    tau <- cor(path_lengths/max(path_lengths), cell_metric, method = "kendall")
  } else if (cluster_order_by == "predecessor"){
    predecessors <- GetPredecessors(graph_cluster, root_cluster) + 1
    tau <- cor(predecessors/max(predecessors), cell_metric, method = "kendall")
  }
}

GetCellMetric <- function(cell_order_by, cell_id, cluster_id, cluster_labels, dist_graph){
  cells_in_cluster <- which(cluster_labels == cluster_id, arr.ind = TRUE)
  if(cell_order_by == "distance"){
    return(mean(dist_graph[cell_id, cells_in_cluster]))
  } else if (cell_order_by == "index"){
    return(mean(cells_in_cluster))
  }
}
