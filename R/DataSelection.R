#' Set of the most variable genes
#' 
#' First filter using the expression threshold.  Then use the coefficient of the top variance PCA components to determine the variability of the gene.
#'
#' @param M a matrix of expression values for each gene (rows) and cell (columns)
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param n_features number of genes to retrieve
#'
#' @return a table of features (rows) and samples (columns)
#'
#' @importFrom Matrix nnzero t
#' @importFrom stats prcomp
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
  browser()
  # select most variable genes based on the max coefficient
  eigenvector_order <- order(genes_pca$sdev, decreasing = TRUE)
  eigenvectors <- genes_pca$rotation[,eigenvector_order]
  max_coeffs <- apply(eigenvectors[,1:max_component], 1, function(x){
    max(abs(x))
  })
  ordered_max <- max_coeffs[order(max_coeffs, decreasing = TRUE)]

  # if the requested feature count is too high, take all available features
  adj_n <- min(nrow(M_variable), n_features)
  nth_highest <- ordered_max[adj_n]
  indices <- which(max_coeffs >= nth_highest, arr.ind = TRUE)
  return(list(gene_use = gene_use[indices], M_variable = M_variable[indices,]))
}

#' Scale and center a data matrix
#' 
#' We want a matrix whose values are standard deviations from the mean and which are centered at zero.
#'
#' @param M a matrix of expression values for each gene (rows) and cell (columns)
#' 
#' @importFrom stats sd
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

#' Remove rare and ubiquitous genes
#' 
#' This is based on Kiselev et al 2017 (SC3).  We remove genes where we observe at least countL copies in less than X% of cells or at least countU copies in more than 100 - X% of cells.  Adjust countL down and countU up in order to get more permissive filters.
#'
#' @param M a matrix of expression values for each gene (rows) and cell (columns)
#' @param countL the definition of "expressed" raw counts for lower threshold
#' @param countU the definition of "expressed" raw counts for upper threshold
#' @param X the threshold percentage: lower = X\%, upper = 100-X\%
#' 
#' @return indices of genes to keep
#' 
#' @export
#'
ApplyCountThreshold <- function(M, countL = 1, countU = 3, X = 3){
  n_cells <- ncol(M)
  
  upper <- apply(M, 1, function(x){
    !((length(which(x >= countU)) / n_cells) >= (100-X)*0.01)
  })
  
  lower <- apply(M, 1, function(x){
    (length(which(x >= countL)) / n_cells) >= X*0.01
  })
  
  keep <- which(as.logical(upper * lower))
}

#' Scale and center a data matrix
#' 
#' This function normalizes the data across cells, scales it and log transforms it
#'
#' @param M a matrix of expression values for each gene (rows) and cell (columns)
#' @param scale the scaling factor
#' @param normalize a boolean whether to normalize the data
#' 
#' @return a log normalized matrix of data
#' 
#' @export
#'
LogNormalize <- function(M, scale = 1e4, normalize = TRUE){
  if(normalize){
    MM <- apply(M, 2, function(x){
      if(max(x) > 0){
        x/max(x)
      }else{
        x
      }
    })
  } else {
    MM <- M
  }
  
  MM <- scale * MM
  
  MM <- log10(MM + 1)
}