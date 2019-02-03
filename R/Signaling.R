#' Compute Cell-cell interaction probability
#'
#' We can have a situation where ligand/receptor expression is low, but target genes (repressors and activators) are highly expressed.  In this case we will get a false positive for \eqn{P_{i,j}}.  In order to correct for this, we introduce normalizing coefficients for the relationship between \eqn{\alpha_{i,j}} and \eqn{\beta / \gamma}.  When \eqn{\alpha} is low, then K and D will decrease rapidly with increasing \eqn{\beta} and \eqn{\gamma}, penalizing the resulting increase in P (which is increasingly likely to be a false positive, given the disparity between \eqn{\alpha} and \eqn{\beta},\eqn{\gamma})
#' \eqn{K_{i,j} = \frac{\alpha_{i,j}}{\alpha_{i,j} + \beta_{i,j}}}
#' Note that activated genes are required under this model.
#'
#' @param M a matrix of expression values for each cell (rows) and gene (columns)
#' @param ids a vector of gene ids
#' @param ligand a named list:
#'     list(ligand1 = list(receptor1 = list(up = list(target1, target2),
#'                                          down = list(target3, target4))))
#' @export
#'
GetSignalingPartners <- function(M,
                                 ids,
                                 ligand,
                                 targets){
  n_cells <- ncol(M)
  n_genes <- nrow(M)

  # represent the pathway heirarchy as a data frame
  pathway <- (reshape2::melt(ligand)[,-2])[,c(4, 3, 2, 1)]
  colnames(pathway) <- c('ligand', 'receptor', 'logic', 'target')

  # ligand-receptor pairs
  LR_pairs <- unique(pathway[,1:2])

  P <- list()
  for(pair in 1:nrow(LR_pairs)){
    # find alpha
    lig_index <- which((LR_pairs[pair, 1]) == ids)
    lnsub <- NormalizeSubset(M, lig_index)

    rec_index <- which((LR_pairs[pair, 2]) == ids)
    rnsub <- NormalizeSubset(M, rec_index)
    alpha <- exp(-1/(matrix(lnsub,n_cells,1) %*%
                       matrix(rnsub,1,n_cells)))

    # find beta
    targ_up <- filter(pathway,
                      ligand == LR_pairs[pair,1] &
                        receptor == LR_pairs[pair,2] &
                        logic == 'up')
    up_index <- match(targ_up[,4],ids)
    if(length(up_index) > 0){
      avg_up <- TargetAvg(M, up_index)
      beta <- exp(-1/avg_up)
    }else{
      beta <- rep(1, n_cells)
    }

    # find gamma
    targ_down <- filter(pathway,
                        ligand == LR_pairs[pair,1] &
                          receptor == LR_pairs[pair,2] &
                          logic == 'down')
    down_index <-match(targ_down[,4],ids)
    if(length(down_index) > 0){
      avg_down <- TargetAvg(M, down_index)
      gamma <- exp(-avg_down)      
    }else{
      gamma <- rep(1, n_cells)
    }

    gammaM <- t(replicate(n_cells, gamma))
    betaM <- t(replicate(n_cells, beta))

    # find K
    if(length(up_index) > 0){
      K <- PenaltyCoeff(alpha = alpha, const = betaM, num = 'alpha', n_cells)
    }else{
      K <- matrix(rep(1, n_cells*n_cells), n_cells, n_cells)
    }
    # find D
    if(length(down_index) > 0){
      D <- PenaltyCoeff(alpha = alpha, const = gammaM, num = 'alpha', n_cells)
    }else{
      D <- matrix(rep(1, n_cells*n_cells), n_cells, n_cells)
    }

    P_num <- alpha*K*betaM*D*gammaM
    # compute the normalizing factor
    tempMat <- replicate(n_cells, apply(P_num, 1, sum))
    non_zero_entries <- which(tempMat > 0)
    denom_Mat <- matrix(0, n_cells, n_cells)
    updateNonZero <- 1/tempMat[non_zero_entries]
    denom_Mat[non_zero_entries] <- updateNonZero
    P[[pair]] <- P_num * denom_Mat
    P[[pair]][P[[pair]] <= 1e-6] <- 0
  }
  P_agg <- Reduce('+', P)/length(P)
  return(list(P = P, P_agg = P_agg))
}

#' Get average target expression for each cell
#'
#' @param M a matrix of expression values for each cell (rows) and gene (columns)
#' @param ind a vector of row indexes corresponding to genes
#'
TargetAvg <- function(M,
                      ids){
  # normalize expression values by max cell expression
  norm <- NormalizeSubset(M, ids)
  avg <- apply(norm, 2, function(x){
    mean(x)
  })
}

#' Normalize a Subset of the Data
#'
#' @param M a matrix of expression values for each cell (columns) and gene (rows)
#' @param ind a vector of row indexes corresponding to genes
#'
NormalizeSubset <- function(M,
                            ids){
  # normalize expression values by max cell expression
  M <- M[ids, ]
  if(is.null(nrow(M))){
    # convert to matrix so apply works
    M <- matrix(M, nrow = 1)
  }
  row_max <- apply(M, 1, function(x){
    if(max(x) > 0){
      max_adj <- max(x)
    } else {
      # if the gene is not expressed, we will divide by 1 so that the normalized expression is 0
      max_adj <- 1
    }
    max_adj
  })
  M_norm <- apply(M, 2, function(x){
    x/row_max
  })
  return(M_norm)
}

#' Generate the divergence penalty coefficient
#'
#' This will approach 1 as the difference between alpha and \code{const} approaches 0
#'
#' @param alpha a matrix of normalized LR expression values for each cell (columns) and gene (rows)
#' @param const either the normalized repression targets or the normalized expression targets
#' @param num either alpha or const
#' @param n_cells the number of cells (row-traversal of the alpha matrix)
#'
PenaltyCoeff <- function(alpha,
                         const,
                         num = 'alpha',
                         n_cells){
  proper <- matrix(0, n_cells, n_cells)
  ind <- which(alpha + const > 0)
  if(num == 'alpha'){
    numerator <- alpha
  } else {
    numerator <- const
  }
  nonzero <- numerator[ind]/(alpha[ind]+const[ind])
  proper[ind] <- nonzero
  proper
}

#' Average cluster to cluster gene expression
#'
#' @param P signaling probabilities cells x cells
#' @param cluster_labels labels of cells 1:n
#' @param remove_zeros boolean whether to remove zero rows from the P matrix
#'
#' @export
#'
ClusterSig <- function(P,
                       cluster_labels,
                       remove_zeros = TRUE){
  sorted_cell <- sort.int(cluster_labels, index.return = TRUE)
  cell_order <- sorted_cell$ix
  cell_labels <- sorted_cell$x

  n_clusters <- length(unique(cell_labels))
  n_cells <- length(cell_labels) 
  # for each cluster, get the location of the first and the last cell in the permutation of labels
  counts <- as.matrix(table(cell_labels))
  accumulator <- matrix(1, n_clusters, n_clusters)*lower.tri(matrix(1, n_clusters, n_clusters), diag = TRUE)
  last_cells <- accumulator %*% counts
  first_cells <- last_cells - counts + 1
  
  # matrix to sum all the cells into clusters by column
  summing_matrix <- NULL
  for(i in 1:n_clusters){
    if(is.null(summing_matrix)){
      summing_matrix <- c(rep(0, n_cells))
      summing_matrix[first_cells[i]:last_cells[i]] <- 1
    } else {
      new_col <- c(rep(0, n_cells))
      new_col[first_cells[i]:last_cells[i]] <- 1
      summing_matrix <- cbind(summing_matrix, new_col)
    }
  }
  
  # collapse the colums by summing, note that because each row is normalized,
  # the following procedure will give normalized rows
  sums <- P[cell_order, cell_order] %*% summing_matrix
  
  if(remove_zeros){
    # get the nonzero rows of P
    nzrow <- which(rowSums(sums) > 0)
    # get the nonzero labels
    nzlabel <- cell_labels[nzrow]
    # adjust the summing matrix to handle non-zero rows only
    nzsumming_matrix <- summing_matrix[nzrow,]
    # get the counts of nonzero cells from each cluster
    nzcounts <- as.matrix(table(nzlabel))
    
    nzsums <- sums[nzrow,]
    nzsums <- t(nzsums) %*% nzsumming_matrix
    avg_matrix <- apply(nzsums, 1, function(x){
      x/nzcounts
    })
  } else {
    sums <- t(sums) %*% summing_matrix
    avg_matrix <- apply(sums, 1, function(x){
      x/counts
    })
  }
  
  arclabs <- paste0("C", c(1:n_clusters))
  rownames(avg_matrix) <- arclabs
  colnames(avg_matrix) <- arclabs
  return(avg_matrix)
}



#' Create a heatmap with signaling markers over clusters
#'
#' @param M a matrix of genes x cells where each entry is normalized expression
#' @param gene_names the names of the genes in the rows of M
#' @param cell_order permutation of the labels of P, cells ordered by cluster
#' @param cell_labels cluster labels of the permuted cells
#' @param gene_list filter by a list of gene names
#'
#' @export
#'
ClusterAvg <- function(M,
                    gene_names,
                    cell_order,
                    cell_labels,
                    gene_list){
  n_clusters <- length(unique(cell_labels))
  n_cells <- length(cell_labels)
  
  # for each cluster, get the location of the first and the last cell in the permutation of labels
  counts <- as.matrix(table(cell_labels))
  accumulator <- matrix(1, n_clusters, n_clusters)*lower.tri(matrix(1, n_clusters, n_clusters), diag = TRUE)
  last_cells <- accumulator %*% counts
  first_cells <- last_cells - counts + 1
  
  summing_matrix <- NULL
  for(i in 1:n_clusters){
    if(is.null(summing_matrix)){
      # we must initialize
      summing_matrix <- c(rep(0, n_cells))
      summing_matrix[first_cells[i]:last_cells[i]] <- 1
    } else {
      new_col <- c(rep(0, n_cells))
      new_col[first_cells[i]:last_cells[i]] <- 1
      summing_matrix <- cbind(summing_matrix, new_col)
    }
  }
  
  sums <- M[, cell_order] %*% summing_matrix
  avg_matrix <- apply(sums, 1, function(x){
    x/counts
  })

  return(avg_matrix[,match(gene_list, gene_names)])
  ## insert code for plot
}
