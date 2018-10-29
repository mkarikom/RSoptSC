#' Compute Cell-cell interaction probability
#'
#' We can have a situation where ligand/receptor expression is low, but target genes (repressors and activators) are highly expressed.  In this case we will get a false positive for P_{i,j}.  In order to correct for this, we introduce normalizing coefficients for the relationship between \alpha_{i,j} and \beta / \gamma.  When \alpha is low, then K and D will decrease rapidly with increasing \beta and \gamma, penalizing the resulting increase in P (which is increasingly likely to be a false positive, given the disparity between \alpha and \beta/\gamma)
#' \eqn{K_{i,j} = \frac{\alpha_{i,j}}{\alpha_{i,j} + \beta_{i,j}}}
#'
#' @param M a matrix of expression values for each cell (rows) and gene (columns)
#' @param ids a vector of gene ids
#' @param ligand a named list:
#'     list(ligand1 = list(receptor1 = list(up = list(target1, target2),
#'                                          down = list(target3, target4))))
#'
GetSignalingPartners <- function(M = M,
                                 ids = ids,
                                 ligand = ligand,
                                 targets){
  n_cells <- nrow(M)

  # represent the pathway heirarchy as a data frame
  pathway <- (reshape2::melt(ligand)[,-2])[,c(4, 3, 2, 1)]
  colnames(pathway) <- c('ligand', 'receptor', 'logic', 'target')

  # ligand-receptor pairs
  LR_pairs <- unique(pathway[,1:2])

  P <- list()
  for(pair in 1:nrow(LR_pairs)){
    # find alpha
    lig_index <- which((LR_pairs[pair, 1]) == ids)
    lig_cells <- t(replicate(n_cells, NormalizeSubset(M, lig_index)))
    rec_index <- which((LR_pairs[pair, 1]) == ids)
    rec_cells <- t(replicate(n_cells, NormalizeSubset(M, lig_index)))
    alpha <- exp(-1/(lig_cells * rec_cells))

    # find beta
    targ_up <- filter(pathway,
                      ligand == LR_pairs[pair,1] &
                        receptor == LR_pairs[pair,2] &
                        logic == 'up')
    up_index <- match(targ_up[,4],ids)
    avg_up <- TargetAvg(M, up_index)
    beta <- exp(-1/avg_up)

    # find gamma
    targ_down <- filter(pathway,
                        ligand == LR_pairs[pair,1] &
                          receptor == LR_pairs[pair,2] &
                          logic == 'down')
    down_index <-match(targ_down[,4],ids)
    avg_down <- TargetAvg(M, down_index)
    gamma <- exp(-avg_down)
    print(paste0("a,b,g ", pair))


    # find K
    K <- PenaltyCoeff(alpha, beta, n_cells)
    print(paste0("K ", pair))

    # find D
    D <- PenaltyCoeff(alpha, gamma, n_cells)
    print(paste0("D ", pair))

    # compute P
    P_num <- alpha*K*beta*D*gamma
    P_denom <- apply(P, 2, function(x){
      sum(x)})
    P_ind <- apply(P_num, 1, function(x){
        x/P_denom
      })
    P[i] <- P_ind
    print(paste0("pair ", pair))
  }

  P_agg <- Reduce('+', P)/length(P)

  return(list(P, P_agg))
}


#' Get average target expression for each cell
#'
#' @param M a matrix of expression values for each cell (rows) and gene (columns)
#' @param ind a vector of row indexes corresponding to genes
#'
TargetAvg <- function(M = M,
                      ids = ids){
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
NormalizeSubset <- function(M = M,
                            ids = ids){
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
}

#' Generate the divergence penalty coefficient
#'
#' This will approach 1 as the difference between alpha and \code{const} approaches 0
#'
#' @param alpha a matrix of normalized LR expression values for each cell (columns) and gene (rows)
#' @param const either the normalized repression targets or the normalized expression targets
#'
PenaltyCoeff <- function(alpha = alpha,
                         const = const,
                         n_cells = n_cells){
  coeff <- apply(alpha, 1, function(x){
    ind <- which(x+const > 0, arr.ind = TRUE)
    nonzero <- as.vector(replicate(n_cells, 0))
    nonzero[ind] <- x[ind]/(x+const)[ind]
  })
  coeff
}
