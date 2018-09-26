#' Get the marker genes for each cluster
#'
#' @param M a matrix of expression values for each cell (rows) and gene (columns)
#' @param cluster_labels a vector of cluster labels
#' @param H a nonnegative matrix such that W = H*t(H), H_{i,j} is the cells_weight
#'     by which cell i belongs to the jth cluster
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param n_features number of marker genes per cluster to retrieve
#'
#' @return a table of marker genes
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
  M$M_variable <- apply(M$M_variable, 2, function(x){
    x / sum(x)
  })

  # for each gene (rows), for each cluster (columns)
  # get cell-weighted expression score
  gene_cluster_scores <- M$M_variable %*% H

  # select max cluster score
  gene_assignments <- apply(gene_cluster_scores, 1, function(x){
    c(which(x == max(x)), max(x))
  })


  marker_table <- cbind(M$gene_use, t(gene_assignments))
  return(marker_table)
}

#' Return a set of the most variable genes
#' First filter using the expression threshold
#' Then use the coefficient of the top variance PCA components
#' to determine the variability of the gene
#'
#' @param M a matrix of expression values for each cell (rows) and gene (columns)
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param n_features number of marker genes per cluster to retrieve
#'
#' @return a table of features (rows) and samples (columns)
#'
SelectData <- function(M, gene_expression_threshold, n_features){
  n_genes <- nrow(M)
  n_cells <- ncol(M)

  # remove genes expressed in min number of cells and not expressed in n - min cells
  if(gene_expression_threshold > 0){
    alpha_filter <- gene_expression_threshold / n_cells
    gene_nnz <- apply(M, 1, function(x){
      nnzero(x) / n_cells
    })
    gene_use <- intersect(which(gene_nnz > alpha_filter, arr.ind = TRUE),
                                which(gene_nnz < 1 - alpha_filter, arr.ind = TRUE))
  } else {
    gene_use = 1:n_genes
  }
  M_variable <- M[gene_use, ]

  # find top components based on max eigengap
  genes_pca <- prcomp(t(M_variable))
  eigengaps <-  abs(genes_pca$sdev[-c(1,length(genes_pca$sdev))] - genes_pca$sdev[-(1:2)])
  max_component <- which(eigengaps == max(eigengaps), arr.ind = TRUE) + 1

  # select most variable genes based on the max coefficient
  eigenvector_order <- order(genes_pca$sdev, decreasing = TRUE)
  eigenvectors <- genes_pca$rotation[,eigenvector_order]
  max_coeffs <- apply(eigenvectors[,1:max_component], 1, function(x){
    max(x)
  })
  ordered_max <- max_coeffs[order(max_coeffs, decreasing = TRUE)]

  # if the requested feature count is too high, take all available features
  adj_n <- min(nrow(M_variable), n_features)
  nth_highest <- ordered_max[adj_n]
  indices <- which(max_coeffs >= nth_highest, arr.ind = TRUE)
  return(list(gene_use = gene_use[indices], M_variable = M_variable[indices,]))
}
