#' Get the marker genes for each cluster
#' 
#' Get the marker genes for each cluster and report their 1-normalized expression value.
#'
#' @param counts_data a matrix of expression values for each cell (rows) and gene (columns) which will be normalized and 
#' @param cluster_labels a vector of cluster labels
#' @param H a nonnegative matrix such that W = H*t(H), H_{i,j} is the cells_weight
#'     by which cell i belongs to the jth cluster
#' @param n_sorted the number of markers to report in the n_sorted table
#' @param gene_names gene symbols corresponding to rows in counts_data (this should match rows of M, if provided)
#' @param gene_expression_threshold for n cells, for \code{gene_expression_threshold} = m, dont consider genes
#'     expressed in more than n-m cells or genes expressed in less than m cells
#' @param n_features number of marker genes per cluster to retrieve
#' @param use_H whether to use H as loading, default is false, instead using cumulative absolute difference to sort max-mean assigned labels
#'
#' @return a list containing 
#'     \item{all}{all the markers listed by gene index, cluster, and score}
#'     \item{n_sorted}{a subset of the sorted marker table with gene symbols}
#'
#' @importFrom dplyr arrange top_n 
#' @importFrom magrittr %>%
#'
#' @export
#'
GetMarkerTable <- function(counts_data,
                            cluster_labels,
                            H = NULL,
                            n_sorted = 50,
                            gene_names,
                            gene_expression_threshold = NULL,
                            n_features = NULL,
                            use_H = FALSE){
  # filter the data according to desired samples (cells) and features (genes)
  clusterId = geneScore = NULL # get r cmd check to pass
  
  if(use_H){
    M <- counts_data
    Mnorm <- LogNormalize(M$M_variable)
    
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
  } else {
    M <- counts_data
    gene_names <- gene_names[counts_data$gene_use]
    counts_data <- counts_data$M_variable
    n_clusters <- length(unique(cluster_labels))
    mean_exp <- matrix(0, nrow = nrow(counts_data), ncol = n_clusters)
    for (clust in 1:n_clusters){
      idx <- which(cluster_labels == clust)
      mean_exp[,clust] <- rowMeans(counts_data[,idx])
    }
    
    cluster_assign <- apply(mean_exp, 1, function(x){
      sorted <- order(x, decreasing = TRUE)                  
      sorted[1]
    })
    cluster_assign <- unlist(cluster_assign)
    
    DE_exp <- matrix(0, nrow = nrow(counts_data), ncol = n_clusters)
    for (clust in 1:n_clusters){
      col_wise <- replicate(n_clusters, mean_exp[,clust])
      col_wise <- abs(col_wise - mean_exp)
      DE_exp[,clust] <- rowSums(col_wise)
    }
    ccol <- c()
    scorecol <- c()
    genecol <- c()
    idxcol <- c()
    for (clust in 1:n_clusters){
      idx <- which(cluster_assign == clust)
      scores <- DE_exp[idx, clust]
      ordering <- order(scores, decreasing = TRUE)
      idx <- idx[ordering]
      
      idxcol <- c(idxcol, M$gene_use[idx])
      scorecol <- c(scorecol, DE_exp[idx, clust])
      ccol <- c(ccol, rep(clust, length(idx)))
      genecol <- c(genecol, gene_names[idx])
    }
    allmarkers <- cbind(geneID = idxcol, clusterId = ccol, geneScore = scorecol)
    rownames(allmarkers) <- genecol
    
    topsorted <- data.frame(geneID = idxcol, clusterId = ccol, geneScore = scorecol, geneSymbol = genecol)
    topsorted <- as_tibble(topsorted) %>% group_by(clusterId) %>% top_n(n_sorted, geneScore)
    output <- list()
    output$all <- allmarkers
    output$n_sorted <- topsorted
  }
  return(output)
}