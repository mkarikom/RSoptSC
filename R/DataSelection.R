#' Return a set of the most variable genes
#' First filter using the expression threshold
#' Then use the coefficient of the top variance PCA components
#' to determine the variability of the gene
#'
#' @param M a matrix of expression values for each gene (rows) and cell (columns)
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param n_features number of genes to retrieve
#'
#' @return a table of features (rows) and samples (columns)
#'
#' @export
#'
SelectData <- function(M, gene_expression_threshold, n_features){
  n_genes <- nrow(M)
  n_cells <- ncol(M)

  # remove genes expressed in min number of cells and not expressed in n - min cells
  if(gene_expression_threshold > 0){
    alpha_filter <- gene_expression_threshold / n_cells
    gene_nnz <- apply(M, 1, function(x){
      Matrix::nnzero(x) / n_cells
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

#' Scale and center a data matrix
#' We want a matrix whose values are standard deviations from the mean and which are centered at zero.
#'
#' @param M a matrix of expression values for each gene (rows) and cell (columns)
#'
#' @return a scaled and centered matrix of expression values
#'
ScaleCenterData <- function(M){
  # get a column vector of gene-wise means across cells
  means <- t(t(apply(M, 1, mean)))
  sds <- t(t(apply(M, 1, sd)))

  # make sure we dont divide by 0
  sds[which(sds == 0)] <- 1

  # scale and center the data
  adjM <- apply(M, 2, function(x){
      (x - means)/sds
  })
  adjM
}
