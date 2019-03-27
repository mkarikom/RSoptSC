#' Get the marker genes for each cluster
#' 
#' Get the marker genes for each cluster and report their 1-normalized expression value.
#'
#' @param counts_data a matrix of expression values for each cell (rows) and gene (columns)
#' @param cluster_labels a vector of cluster labels
#' @param H a nonnegative matrix such that W = H*t(H), H_{i,j} is the cells_weight
#'     by which cell i belongs to the jth cluster
#' @param n_sorted the number of markers to report in the n_sorted table
#' @param gene_names gene symbols corresponding to rows in counts_data
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param n_features number of marker genes per cluster to retrieve
#'
#' @return a list containing 
#'     \item{all}{all the markers listed by gene index, cluster, and score}
#'     \item{n_sorted}{a subset of the sorted marker table with gene symbols}
#'
#' @export
#'
GetMarkerTable <- function(counts_data,
                           cluster_labels,
                           H,
                           n_sorted = 50,
                           gene_names,
                           gene_expression_threshold,
                           n_features){
  # filter the data according to desired samples (cells) and features (genes)
  clusterId = geneScore = NULL # get r cmd check to pass
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
  
  sorted_table <- arrange(as.data.frame(marker_table), clusterId, desc(geneScore))
  sorted_table <- as_tibble(sorted_table) %>% group_by(clusterId) %>% top_n(n_sorted, geneScore)
  sorted_table$geneSymbol <- gene_names[sorted_table$geneID]
  
  output <- list()
  output$all <- marker_table
  output$n_sorted <- sorted_table
  return(output)
}
