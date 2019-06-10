#' Get the marker genes for each cluster
#' 
#' Get the marker genes for each cluster and report their 1-normalized expression value.
#'
#' @param counts_data a matrix of expression values for each cell (rows) and gene (columns) 
#' @param cluster_labels a vector of cluster labels
#' @param H a nonnegative matrix such that W = H*t(H), H_{i,j} is the cells_weight
#'     by which cell i belongs to the jth cluster
#' @param n_sorted the number of markers to report in the n_sorted table
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param use_H whether to use H as loading, default is false, instead using cumulative absolute difference to sort max-mean assigned labels
#' @param lognorm perform log transform, normalization, and scaling on the data before the computing the markers
#'
#' @return a list containing 
#'     \item{all}{all the markers listed by gene index, cluster, and score}
#'     \item{n_sorted}{a subset of the sorted marker table with gene symbols}
#'
#' @importFrom dplyr arrange top_n 
#' @importFrom magrittr %>%
#' @importFrom Matrix rowMeans
#'
#' @export
#'
GetMarkerTable <- function(counts_data,
                           cluster_labels,
                           H = NULL,
                           n_sorted = 50,
                           gene_expression_threshold = NULL,
                           use_H = FALSE,
                           lognorm = FALSE){
  # filter the data according to desired samples (cells) and features (genes)
  clusterId = geneScore = NULL # get r cmd check to pass
  genes <- matrix(rownames(counts_data), ncol = 1)
  cells <- colnames(counts_data)
  if(use_H){
    if(lognorm){
      Mnorm <- LogNormalize(counts_data)
    }else{
      Mnorm <- counts_data
    }
    
    # for each gene (rows), for each cluster (columns)
    # get cell-weighted expression score
    gene_cluster_scores <- Mnorm %*% H
    
    # select max cluster score
    gene_assignments <- matrix(apply(gene_cluster_scores, 2, max, 
                                     simplify = TRUE), ncol = 1)
    
    marker_table <- cbind(genes, gene_assignments)

    colnames(marker_table) <- c('geneID', 'clusterId', 'geneScore')

    sorted_table <- arrange(as.data.frame(marker_table), clusterId, desc(geneScore))
    sorted_table <- sorted_table %>% group_by(clusterId) %>% top_n(n_sorted, geneScore)
    sorted_table$geneSymbol <- genes[sorted_table$geneID]
    
    output <- list()
    output$all <- marker_table
    output$n_sorted <- sorted_table
  } else {
    # number of clusters
    # mean expression of each cluster
    n_clusters <- length(unique(cluster_labels))
    mean_exp <- matrix(0, nrow = nrow(counts_data), ncol = n_clusters)
    for (clust in 1:n_clusters){
      idx <- which(cluster_labels == clust)
      mean_exp[,clust] <- rowMeans(counts_data[,idx])
    }
    
    # find the cluster with the hightest mean expression of each gene
    cluster_assign <- apply(mean_exp, 1, function(x){
      sorted <- order(x, decreasing = TRUE)                  
      sorted[1]
    })
    cluster_assign <- unlist(cluster_assign)
    
    # for each gene x cluster
    # get difference betwen expression in the cluster and expression in all cells
    DE_exp <- matrix(0, nrow = nrow(counts_data), ncol = n_clusters)
    for (clust in 1:n_clusters){
      col_wise <- replicate(n_clusters, mean_exp[,clust])
      col_wise <- abs(col_wise - mean_exp)
      DE_exp[,clust] <- rowSums(col_wise)
    }
    
    # for each cluster, order the genes whose mean expression 
    # is max in that cluster by differential expression
    ccol <- c()
    scorecol <- c()
    genecol <- c()
    idxcol <- c()
    for (clust in 1:n_clusters){
      idx <- which(cluster_assign == clust) # find genes
      scores <- DE_exp[idx, clust] # get their scores
      ordering <- order(scores, decreasing = TRUE) # order the scores
      idx <- idx[ordering] # order the genes
      
      idxcol <- c(idxcol, idx)
      scorecol <- c(scorecol, DE_exp[idx, clust])
      ccol <- c(ccol, rep(clust, length(idx)))
      genecol <- c(genecol, genes[idx])
    }
    
    df <- data.frame(geneID = idxcol, clusterId = ccol, 
                            geneScore = scorecol, geneSymbol = genecol)
    top_n <- df %>% group_by(clusterId) %>% top_n(n_sorted, geneScore)
    allmarkers <- df %>% group_by(clusterId) %>% top_n(Inf, geneScore)
    output <- list(all = allmarkers, top_n =  top_n)
  }
  return(output)
}