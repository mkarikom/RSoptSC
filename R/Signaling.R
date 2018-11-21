#' Compute Cell-cell interaction probability
#'
#' We can have a situation where ligand/receptor expression is low, but target genes (repressors and activators) are highly expressed.  In this case we will get a false positive for \eqn{P_{i,j}}.  In order to correct for this, we introduce normalizing coefficients for the relationship between \eqn{\alpha_{i,j}} and \eqn{\beta / \gamma}.  When \eqn{\alpha} is low, then K and D will decrease rapidly with increasing \eqn{\beta} and \eqn{\gamma}, penalizing the resulting increase in P (which is increasingly likely to be a false positive, given the disparity between \eqn{\alpha} and \eqn{\beta},\eqn{\gamma})
#' \eqn{K_{i,j} = \frac{\alpha_{i,j}}{\alpha_{i,j} + \beta_{i,j}}}
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
    print(LR_pairs[pair,])
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

    gammaM <- t(replicate(n_cells, gamma))
    betaM <- t(replicate(n_cells, beta))

    # find K
    K <- PenaltyCoeff(alpha = alpha, const = betaM, num = 'alpha', n_cells)
    print(paste0("K ", pair))

    # find D
    D <- PenaltyCoeff(alpha = alpha, const = gammaM, num = 'alpha', n_cells)
    print(paste0("D ", pair))

    P_num <- alpha*K*betaM*D*gammaM
    #browser()
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
  #browser()
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
