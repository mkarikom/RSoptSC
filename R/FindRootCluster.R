#' Find the root cluster, given a weighted adjacency matrix
#' Generate a cluser-to-cluster graph based on a similarity-matrix of connected cells
#' Generate the minimum spanning tree on the cluster-to-cluster graph
#' Use one of the max dist clusters from this graph, unless
#'     The cluster is given, then proceed to find the root cell
#'     The root cell is given, then use its cluster
#'
#' @param cluster_labels the cluster label for each cell
#' @param flat_embedding rows are cells and columns are coordinates in n-col space
#' @param dist_graph distances on a connected graph
#' @param dist_flat flat distance matrix for all cells
#' @param reverse a boolean variable whether to take the root cluster based on minimum dispersion on the flat embedding
#' @return a list containing:
#'     \item{cluster_adj_matrix}{a weighted upper-triangular adjacency matrix for the clusters based on avg pseudotime between their cells}
#'     \item{graph_cluster}{an igraph object on the clusters, with a minimum spanning tree}
#'     \item{root_cluster}{index of the root cluster}
#'     \item{cluster_mst}{an igraph mst object}
#'
#'
#' @examples FindRootCluster(cluster_labels = RSoptSC::GuoPtime$Params$cluster_label,
#'           flat_embedding = RSoptSC::GuoPtime$Params$latent,
#'           dist_graph = RSoptSC::GuoPtime$Values$Short_pathd,
#'           dist_flat = RSoptSC::GuoPtime$Values$low_dis,
#'           reverse = FALSE)
FindRootCluster <- function(cluster_labels, flat_embedding, dist_graph, dist_flat, reverse = FALSE){
  n_clusters <- length(unique(cluster_labels))
  cluster_adj_matrix <- matrix(0, nrow = n_clusters, ncol = n_clusters)
  for(i in 1:(n_clusters-1)){
    for(j in (i+1):n_clusters){
      i_j_distances <- dist_graph[which(cluster_labels == i), which(cluster_labels == j)]
      n_i_j <- prod(dim(i_j_distances))
      cluster_adj_matrix[i, j] <- sum(i_j_distances)/n_i_j
    }
  }
  graph_cluster <- igraph::graph_from_adjacency_matrix(adjmatrix = cluster_adj_matrix,
                                                       mode = "upper",
                                                        weighted = TRUE)

  cluster_mst <- igraph::mst(graph_cluster)
  root_cluster_candidates <- which(cluster_adj_matrix == max(cluster_adj_matrix), arr.ind = TRUE)
  cluster_variance <- FindVariance(n_clusters, cluster_labels, flat_embedding)
  root_cluster_variance <- cluster_variance[root_cluster_candidates]
  if(reverse){
    max <- which(root_cluster_variance == max(root_cluster_variance), arr.ind = TRUE)
    root_cluster <- root_cluster_candidates[max]
  }else{
    min <- which(root_cluster_variance == min(root_cluster_variance), arr.ind = TRUE)
    root_cluster <- root_cluster_candidates[min]
  }
  clusters <- list(cluster_adj_matrix = cluster_adj_matrix,
                   graph_cluster = graph_cluster,
                   root_cluster = root_cluster,
                   cluster_mst = cluster_mst)
}

FindVariance <- function(n_clusters, cluster_labels, flat_embedding){
  cluster_variance <- c()
  for(i in 1:n_clusters){
    c_embed <- flat_embedding[which(cluster_labels == i),]
    vec_norms <- apply(c_embed, 1, function(x){
      sqrt(sum(x^2))
    })
    cluster_variance <- c(cluster_variance, var(vec_norms))
  }
  return(cluster_variance)
}
