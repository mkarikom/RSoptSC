#' generate convenient representations of the data
#'
#' @param low_dim_embedding the flat embedding of the cells
#' @param similarity_matrix the graphical embedding of the cells
#' @param join_components boolean, whether or not to join disconnected components of the similarity matrix
#' @return a list containing:
#'     \item{dist_flat}{distance matrix on flat embedding}
#'     \item{dist_graph}{distance matrix on graph}
#'     \item{adj_matrix}{unweighted adjacency matrix}
#'     \item{similarity_graph}{igraph object on the unweighted adjacency matrix}
#'     \item{components}{an igraph components object based on \code{similarity_graph}}
#'     \item{n_components}{number of components in \code{components}}
#'     \item{sizes}{sizes of components in \code{components}}
#'     \item{members}{list of members of \code{components}, sorted by component size and member index}
#'
#'
#' @examples embeddings <- RepresentationMap(low_dim_embedding, adjacency_matrix)
#'
#' representations of the data are necessary for subsequent lineage analysis

RepresentationMap <- function(low_dim_embedding, similarity_matrix, join_components = TRUE){
  dist_flat <- as.matrix(dist(low_dim_embedding, method="euclidean", diag = TRUE))

  numcells <- dim(similarity_matrix)[1]
  adj_matrix <- matrix(0, nrow = numcells, ncol = numcells)
  adj_matrix[similarity_matrix>0] <- 1
  diag(adj_matrix) <- 0

  similarity_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = c("undirected"))
  components <- igraph::components(similarity_graph)
  n_components <- components$no
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
                   similarity_graph = similarity_graph,
                   components = components,
                   n_components = n_components,
                   sizes = sizes,
                   members = members)
}
