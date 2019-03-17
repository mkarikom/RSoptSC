#' Get the marker genes for each cluster
#' 
#' Get the marker genes for each cluster and report their 1-normalized expression value.
#'
#' @param counts_data a matrix of expression values for each cell (rows) and gene (columns)
#' @param cluster_labels a vector of cluster labels
#' @param H a nonnegative matrix such that W = H*t(H), H_{i,j} is the cells_weight
#'     by which cell i belongs to the jth cluster
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param n_features number of marker genes per cluster to retrieve
#'
#' @return a table of marker genes
#'
#' @export
#'
GetMarkerTable <- function(counts_data,
                           cluster_labels,
                           H,
                           gene_expression_threshold,
                           n_features){
  # filter the data according to desired samples (cells) and features (genes)

  M <- SelectData(counts_data,
                  gene_expression_threshold,
                  n_features)

  # normalize the expression cell by cell
  Mnorm <- apply(M$M_variable, 2, function(x){
    x / sum(x)
  })

  # for each gene (rows), for each cluster (columns)
  # get cell-weighted expression score
  gene_cluster_scores <- Mnorm %*% H

  # select max cluster score
  gene_assignments <- apply(gene_cluster_scores, 1, function(x){
    c(which(x == max(x)), max(x))
  })


  marker_table <- cbind(M$gene_use, t(gene_assignments))
  colnames(marker_table) <- c('geneID', 'clusterId', 'geneScore')
  return(marker_table)
}
