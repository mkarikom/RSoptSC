#' generate convenient representations of the data
#' representations of the data are necessary for subsequent lineage analysis
#'
#' @param similarity_matrix the graphical embedding of the cells
#' @param join_components boolean, whether or not to join disconnected components of the similarity matrix
#' @param flat_embedding optionally provided low dim embedding, if not then 2d tsne will be used
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
#' @export
#'
RepresentationMap <- function(flat_embedding = NULL,
                              similarity_matrix, 
                              join_components = TRUE, 
                              ...){
  #browser()
  if(is.null(flat_embedding)){
    # set the seed for tsne
    set.seed(1)
    # initialize with low dim embedding
    print("computing flat embedding")
    lat_pc = prcomp(similarity_matrix, center = TRUE)
    # flat_embedding = Rtsne::Rtsne(X = lat_pc$x[,1:2],
    #                               perplexity = perplexity,
    #                               pca_center = TRUE,
    #                               pca_scale = TRUE,
    #                               dims = 2)$Y
    flat_embedding = Rtsne::Rtsne(X = similarity_matrix, ...)$Y
  }

  dist_flat <- as.matrix(dist(flat_embedding, method="euclidean", diag = TRUE))

  numcells <- dim(similarity_matrix)[1]
  adj_matrix <- matrix(0, nrow = numcells, ncol = numcells)
  adj_matrix[similarity_matrix>0] <- 1
  diag(adj_matrix) <- 0

  similarity_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = c("undirected"))
  components <- igraph::components(similarity_graph)
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
    similarity_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = c("undirected"))
  }
  
  sizes <- sort(components$csize, decreasing = TRUE)
  members <- lapply(sizes, function(x){
    which(components$membership == unname(which(table(components$membership)==x)))
  })

  members <- lapply(members, function(x) sort(x))
  members <- members[order(sapply(members, length), decreasing=FALSE)]

  dist_graph <- igraph::distances(similarity_graph)

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

# give a list of vectors where list[[i]] = c(cells in component i)
MembershipList <- function(membership_vector){
  components <- unique(membership_vector)
  membership_list <- list()
  for(i in 1:length(components)){
    membership_list[[i]] <- which(membership_vector == i)
  }
  return(membership_list)
}




