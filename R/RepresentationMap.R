#' Generate multiple representations
#' 
#' Generate representations of the data to be used in lineage analysis.
#'
#' @param similarity_matrix the graphical embedding of the cells
#' @param join_components boolean, whether or not to join disconnected components of the similarity matrix
#' @param flat_embedding optionally provided low dim embedding, if not then 2d tsne will be used
#' @param normalize_S whether or not to normalize the similarity matrix
#' @param flat_embedding_method what method to use for flat embedding.  Use umap by default.  Can be umap or tsne.
#' @param ... arguments to called functions
#' @return a list containing:
#'     \item{dist_flat}{distance matrix on flat embedding}
#'     \item{dist_graph}{distance matrix on graph}
#'     \item{adj_matrix}{unweighted adjacency matrix}
#'     \item{flat_embedding}{the low dimensional embedding of the cells}
#'     \item{similarity_graph}{igraph object on the unweighted adjacency matrix}
#'     \item{components}{an igraph components object based on \code{similarity_graph}}
#'     \item{n_components}{number of components in \code{components}}
#'     \item{sizes}{sizes of components in \code{components}}
#'     \item{members}{list of members of \code{components}, sorted by component size and member index}
#'
#' @importFrom igraph graph_from_adjacency_matrix components distances
#' @importFrom Rtsne Rtsne
#' @importFrom Matrix as.matrix which
#' @importFrom umap umap
#' @importFrom stats dist
#' 
#' @export
#'
RepresentationMap <- function(flat_embedding = NULL,
                               similarity_matrix, 
                               join_components = TRUE,
                               normalize_S = TRUE,
                               flat_embedding_method = 'umap',
                               ...){
  set.seed(1)
  if(normalize_S){
    similarity_matrix <- similarity_matrix / sum(similarity_matrix)
  }
  
  if(is.null(flat_embedding)){
    if(flat_embedding_method == 'tsne'){
      # set the seed for tsne
      # initialize with low dim embedding
      flat_embedding <- Rtsne(X = as.matrix(similarity_matrix), ...)$Y
      flat_embedding <- as.data.frame(flat_embedding)
      colnames(flat_embedding) <- c("C1", "C2")
    }else if(flat_embedding_method == 'umap'){
      # initialize with low dim embedding
      data.umap <- umap(d = as.matrix(similarity_matrix), ...)
      flat_embedding <- data.umap$layout
      flat_embedding <- as.data.frame(flat_embedding)
      colnames(flat_embedding) <- c("C1", "C2")
    }
  }
  
  dist_flat <- as.matrix(dist(flat_embedding, method="euclidean", diag = TRUE))
  
  numcells <- dim(similarity_matrix)[1]
  adj_matrix <- matrix(0, nrow = numcells, ncol = numcells)
  adj_matrix[which(similarity_matrix>0)] <- 1
  diag(adj_matrix) <- 0
  
  similarity_graph <- graph_from_adjacency_matrix(adj_matrix, mode = c("undirected"))
  components <- components(similarity_graph)
  n_components <- components$no
  
  ### Join the components if more than one
  # find an initial cell
  if(n_components > 1 & join_components){
    cell_init <- which(dist_flat == max(dist_flat), arr.ind = TRUE)[1]
    membership_list <- MembershipList(components$membership)
    adj_matrix <- JoinGraphComponents(root_cell = cell_init,
                                      adjacency_matrix = adj_matrix,
                                      flat_distances = dist_flat,
                                      n_components = n_components,
                                      component_members = membership_list)
    similarity_graph <- graph_from_adjacency_matrix(adj_matrix, mode = c("undirected"))
  }
  
  sizes <- sort(components$csize, decreasing = TRUE)
  members <- lapply(sizes, function(x){
    which(components$membership == unname(which(table(components$membership)==x)))
  })
  
  members <- lapply(members, function(x) sort(x))
  members <- members[order(sapply(members, length), decreasing=FALSE)]
  
  dist_graph <- distances(similarity_graph)
  
  mappings <- list(dist_flat = dist_flat,
                   dist_graph = dist_graph,
                   adj_matrix = adj_matrix,
                   flat_embedding = flat_embedding,
                   similarity_graph = similarity_graph,
                   components = components,
                   n_components = n_components,
                   sizes = sizes,
                   members = members)
}

#' Get a list of membership vectors
#' 
#' Get a list of membership vectors to facilitate component-joining
#'
#' @param membership_vector the component membership of cells
#' @return where list[[i]] = c(cells in component i)
#' 
MembershipList <- function(membership_vector){
  components <- unique(membership_vector)
  membership_list <- list()
  for(i in 1:length(components)){
    membership_list[[i]] <- which(membership_vector == i)
  }
  return(membership_list)
}




