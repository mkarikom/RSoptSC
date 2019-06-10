#' Compute Cell-cell interaction probability
#'
#' We can have a situation where ligand/receptor expression is low, but target genes (repressors and activators) are highly expressed.  In this case we will get a false positive for \eqn{P_{i,j}}.  In order to correct for this, we introduce normalizing coefficients for the relationship between \eqn{\alpha_{i,j}} and \eqn{\beta / \gamma}.  When \eqn{\alpha} is low, then K and D will decrease rapidly with increasing \eqn{\beta} and \eqn{\gamma}, penalizing the resulting increase in P (which is increasingly likely to be a false positive, given the disparity between \eqn{\alpha} and \eqn{\beta},\eqn{\gamma})
#' \eqn{K_{i,j} = \frac{\alpha_{i,j}}{\alpha_{i,j} + \beta_{i,j}}}
#' Note that activated genes are required under this model.
#'
#' @param M a matrix of expression values for each cell (rows) and gene (columns)
#' @param ids a vector of gene ids
#' @param pathway a data frame with "ligands", "receptors", "direction", and "targets"
#' @param normalize_aggregate whether or not to normalize the P_agg matrix by dividing each row by its sum default is true
#'                                          
#' @return a list containing:
#'     \item{P}{a list of the cell-cell signaling probabilities for all ligand/receptor pairs}
#'     \item{P_agg}{the aggregate matrix of all ligand/receptor pairs}
#'                                          
#' @importFrom reshape2 melt
#' @importFrom dplyr filter mutate_all
#' @importFrom Matrix as.matrix rowSums colMeans
#' 
#' @export
#'
GetSignalingPartners <- function(M,
                                 ids,
                                 pathway,
                                 normalize_aggregate = TRUE){
  receptor = ligand = direction = NULL   # r cmd check pass
  n_cells <- ncol(as.matrix(M))
  n_genes <- nrow(as.matrix(M))
  cellnames <- colnames(M)
  
  # fix the names
  ids <- tolower(ids)
  pathway <- mutate_all(pathway, .funs = tolower)
  
  # ligand-receptor pairs
  LR_pairs <- unique(pathway[c("receptor", "ligand")])
  
  # version of model
  target <- FALSE
  if('target' %in% colnames(pathway)){
    target <- TRUE
  }
  
  P <- list()
  for(pair in 1:nrow(LR_pairs)){
    P[[pair]] <- matrix(0, nrow = n_cells, ncol = n_cells)
    rownames(P[[pair]]) <- cellnames
    colnames(P[[pair]]) <- cellnames
    
    rec_index <- which(ids == LR_pairs[pair, 1])
    lig_index <- which(ids == LR_pairs[pair, 2])
    ldata <- NormalizeGene(M[lig_index,])
    rdata <- NormalizeGene(M[rec_index,])
    
    alpha <- matrix(0, nrow = n_cells, ncol = n_cells)
    for(i in 1:n_cells){
      for(j in 1:n_cells){
        alpha[i,j] <- exp(-1/(ldata[i] * rdata[j]))
      }
    }
    if(target){
      # note the following may require at least 2 targets of up/down, up, or down
      pathsubset <- pathway %>% filter(receptor == LR_pairs[pair, 1] &
                                         ligand == LR_pairs[pair, 2])
      targupdown <- targdown <- targup <- FALSE
      if('up' %in% pathsubset$direction & 'down' %in% pathsubset$direction){
        targupdown <- TRUE
      } else if('up' %in% pathsubset$direction){
        targup <- TRUE
      } else {
        targdown <- TRUE
      }
      targets <- pathsubset$target 
      if(targupdown){
        upsubset <- pathsubset %>% filter(direction == 'up')
        updata <- M[drop=FALSE,match(upsubset$target,ids),]
        for(targind in 1:nrow(as.matrix(updata))){
          updata[targind,] <- NormalizeGene(updata[targind,])
        }
        updata <- colMeans(updata)          
        beta <- exp(-1/updata)
        
        downsubset <- pathsubset %>% filter(direction == 'down')
        downdata <- M[drop=FALSE,match(downsubset$target,ids),]
        for(targind in 1:nrow(as.matrix(downdata))){
          downdata[targind,] <- NormalizeGene(downdata[targind,])
        }
        downdata <- colMeans(downdata)          
        gamma <- exp(-downdata)
        K <- D <- alpha
        for(i in 1:n_cells){
          for(j in 1:n_cells){
            if(beta[j] <= 0){
              K[i,j] <- 0
            } else {
              K[i,j] <- alpha[i,j] / (alpha[i,j] + beta[j])
            }
            if(gamma[j] <= 0){
              D[i,j] <- 0
            } else {
              D[i,j] <- alpha[i,j] / (alpha[i,j] + gamma[j])
            }
          }
        }
        for(ii in 1:n_cells){
          b <- sum(alpha[ii,] * beta * K[ii,] * gamma * D[ii,])
          for(kk in 1:n_cells){
            a <- alpha[ii,kk] * beta[kk] * K[ii,kk] * gamma[kk] * D[ii,kk]
            if(b == 0){
              P[[pair]][ii,kk] <- 0
            }else{
              P[[pair]][ii,kk] <- a/b
            }
          }
        }
      } else if(targup){
        upsubset <- pathsubset %>% filter(direction == 'up')
        updata <- M[drop=FALSE,match(upsubset$target,ids),]
        for(targind in 1:nrow(as.matrix(updata))){
          updata[targind,] <- NormalizeGene(updata[targind,])
        }
        updata <- base::colMeans(updata)          
        beta <- exp(-1/updata)
        K <- alpha
        for(i in 1:n_cells){
          for(j in 1:n_cells){
            if(beta[j] <= 0){
              K[i,j] <- 0
            } else {
              K[i,j] <- alpha[i,j] / (alpha[i,j] + beta[j])
            }
          }
        }
        for(ii in 1:n_cells){
          b <- sum(alpha[ii,] * beta * K[ii,])
          for(kk in 1:n_cells){
            a <- alpha[ii,kk] * beta[kk] * K[ii,kk]
            if(b == 0){
              P[[pair]][ii,kk] <- 0
            }else{
              P[[pair]][ii,kk] <- a/b
            }
          }
        }
      } else {
        downsubset <- pathsubset %>% filter(direction == 'down')
        downdata <- M[drop=FALSE,match(downsubset$target,ids),]
        for(targind in 1:nrow(as.matrix(downdata))){
          downdata[targind,] <- NormalizeGene(downdata[targind,])
        }
        downdata <- colMeans(downdata)          
        gamma <- exp(-downdata)
        D <- alpha
        for(i in 1:n_cells){
          for(j in 1:n_cells){
            if(gamma[j] <= 0){
              D[i,j] <- 0
            } else {
              D[i,j] <- alpha[i,j] / (alpha[i,j] + gamma[j])
            }
          }
        }
        for(ii in 1:n_cells){
          b <- sum(alpha[ii,] * gamma * K[ii,])
          for(kk in 1:n_cells){
            a <- alpha[ii,kk] * gamma[kk] * K[ii,kk]
            if(b == 0){
              P[[pair]][ii,kk] <- 0
            }else{
              P[[pair]][ii,kk] <- a/b
            }
          }
        }      
      }
    } else {
      # no targets were provided
      for(ii in 1:n_cells){
        b <- sum(alpha[ii,])
        for(kk in 1:n_cells){
          a <- alpha[ii,kk]
          if(b == 0){
            P[[pair]][ii,kk] <- 0
          }else{
            P[[pair]][ii,kk] <- a/b
          }
        }
      }      
    }
    P[[pair]][P[[pair]] <= 1e-6] <- 0
  }
  
  if(normalize_aggregate){
    P_agg <- Reduce('+', P)
    nnorm <- apply(P_agg, 1, function(x){
      if(max(x) > 0){
        x/sum(x)
      }else{
        x
      }
    })
    P_agg <- t(nnorm)
  }else{
    P_agg <- Reduce('+', P)/length(P)
  }
  rownames(P_agg) <- cellnames
  colnames(P_agg) <- cellnames
  return(list(P = P, P_agg = P_agg))
}

#' Normalize signal across cells
#'
#' Normalize signal across cells
#'
#' @param M a matrix of expression values for each cell (columns) and gene (rows)
#'
#' @return a normalized vector of expression values
#'
NormalizeGene <- function(M){
  # normalize expression values by max cell expression
  M_norm <- M
  if(max(M_norm) > 0){
    M_norm <- M_norm/max(M_norm)
  }
  return(M_norm)
}

#' Average cluster to cluster signaling
#'
#' Average cluster to cluster signaling.
#'
#' @param P signaling probabilities cells x cells
#' @param cluster_labels labels of cells 1:n
#' @param normalize_rows normalize the rows of the output matrix so that the rows sum to 1, default true
#'
#' @return a matrix of cluster to cluster signaling
#'
#' @export
#'
ClusterSig <- function(P,
                       cluster_labels,
                       normalize_rows = TRUE){
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
  
  nzrow <- which(rowSums(P[cell_order, cell_order]) > 0)
  nzlabel_row <- cell_labels[nzrow]
  nzcounts_row <- matrix(0, nrow = n_clusters, ncol = 1)
  rownames(nzcounts_row) <- rownames(counts)  
  cnums_row <- sort(as.integer(rownames(counts)))
  for(c in 1:n_clusters){
    clust <- cnums_row[c]
    nzcounts_row[c,1] <- length(which(nzlabel_row == clust))
  }
  
  nzcol <- which(colSums(P[cell_order,cell_order]) > 0)
  nzlabel_col <- cell_labels[nzcol]
  nzcounts_col <- matrix(0, nrow = n_clusters, ncol = 1)
  rownames(nzcounts_col) <- rownames(counts)  
  cnums_col <- sort(as.integer(rownames(counts)))
  for(c in 1:n_clusters){
    clust <- cnums_col[c]
    nzcounts_col[c,1] <- length(which(nzlabel_col == clust))
  }
  
  sums <- P[cell_order, cell_order] %*% summing_matrix
  sums <- t(sums) %*% summing_matrix
  sums <- t(sums)
  
  for(col in 1:ncol(sums)){
    for(row in 1:nrow(sums)){
      if(nzcounts_row[row,1] > 0 & nzcounts_col[col,1] > 0){
        sums[row,col] <- sums[row,col]/nzcounts_row[row,1]/nzcounts_col[col,1]
      }else{
        sums[row,col] <- 0
      }
    }
  }
  
  if(normalize_rows){
    normsums <- apply(sums, 1, function(x){
      if(max(x) > 0){
        x / sum(x)
      }else{
        x
      }
    })
    sums <- t(normsums)
  }
  
  arclabs <- paste0("C", c(1:n_clusters))
  rownames(sums) <- rownames(counts)
  colnames(sums) <- rownames(counts)
  return(sums)
}

#' Create a heatmap with signaling markers over clusters
#' 
#' Create a heatmap with signaling markers over clusters.
#'
#' @param M a matrix of genes x cells where each entry is normalized expression
#' @param gene_names the names of the genes in the rows of M
#' @param cell_order permutation of the labels of P, cells ordered by cluster
#' @param cell_labels cluster labels of the permuted cells
#' @param gene_list filter by a list of gene names
#' 
#' @return a matrix of avg expression values
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

#' Process siganling pathway data
#'
#' Given two tsv tables, one including ligand-receptor pairs, and the other including receptor-target pairs, and a data set, inner-join the two tables and make sure no ligands, receptors, or targets are missing from the data.  This function returns a summary table of all requested markers.  If the data is a subset of some larger study, genes may be reported (a column exists in the data matrix) for which zero copies were observed.  In this case, "dropout in subset" applies.  There are also known gene symbols which were not measured in the study (no corresponding row exists in the data matrix).  In this case, "missing from all subsets" applies.
#'
#' @param lig_table a data frame with "ligand" in column 1 and "receptor" in column 2
#' @param rec_table a data frame with "receptor" in column 1, "target" in column 2, and up/down "direction"
#' @param lig_table_path a tsv/csv file with "ligand" in column 1 and "receptor" in column 2
#' @param rec_table_path a tsv/csv file with "receptor" in column 1, "target" in column 2, and up/down "direction" in column 3
#' @param data a matrix of the data, with gene names as row names (required if no gene list provided)
#' @param gene_names an optional list of gene names (required if rownames of data are null)
#'                                           
#' @return a list containing:
#'     \item{lignames}{vector of all ligands}
#'     \item{recnames}{vector of all receptors}
#'     \item{targets}{vector of all targets}
#'     \item{pathway}{a table with combined pathway information}
#'     \item{pathway_removed}{a table with markers excluded}
#'     \item{removed}{a table with the excluded markers}
#'           
#' @importFrom Matrix rowSums
#' @importFrom utils read.delim
#' 
#' @export
#'
ImportPathway <- function(lig_table = NULL,
                          rec_table = NULL,
                          lig_table_path,
                          rec_table_path = NULL,
                          data,
                          gene_names = NULL){
  if(is.null(rownames(data))){
    genes <- gene_names
  } else {
    genes <- rownames(data)
  }
  
  if(is.null(lig_table)){
    tableLig <- read.delim(lig_table_path, stringsAsFactors=FALSE)
    colnames(tableLig) <- tolower(sub('\\.', tolower(colnames(tableLig)), replacement = "_"))
  } else {
    tableLig <- lig_table
    colnames(tableLig) <- tolower(sub('\\.', tolower(colnames(tableLig)), replacement = "_"))
  }
  if(is.null(rec_table)){
    if(is.null(rec_table_path)){
      targets = FALSE
    } else {
      targets = TRUE
      tableRec <- read.delim(rec_table_path, stringsAsFactors=FALSE)
      colnames(tableRec) <- tolower(sub('\\.', tolower(colnames(tableRec)), replacement = "_"))
    }
  } else {
    targets = TRUE
    tableRec <- rec_table
    colnames(tableRec) <- tolower(sub('\\.', tolower(colnames(tableRec)), replacement = "_"))
  }
  if(targets){
    # remove references and ambiguous targets
    tableRec <- tableRec[which(tableRec$direction != "?"),
                         which(names(tableRec) %in% c("receptor", "target", "direction"))]
    
    # make gene names systematic
    tableRec$receptor <- unname(sapply(tableRec$receptor, SimpleCap))
    tableRec$target <- unname(sapply(tableRec$target, SimpleCap))
    tableLig$ligand <- unname(sapply(tableLig$ligand, SimpleCap))
    tableLig$receptor <- unname(sapply(tableLig$receptor, SimpleCap))
    genes <- unname(sapply(genes, SimpleCap))
    
    # find all the symbols in the pathway
    lignames <- unique(tableLig$ligand)
    recnames <- unique(tableLig$receptor)
    targets <- unique(tableRec$target)
    
    # inner join the tables
    pathway <- merge(tableLig, tableRec)
    
    # pathway elements for which we have no data
    ms_lig_ind <- which(is.na(match(lignames, genes)))
    ms_lig <- lignames[ms_lig_ind]
    
    ms_rec_ind <- which(is.na(match(recnames, genes)))
    ms_rec <- recnames[ms_rec_ind]
    
    ms_targ_ind <- which(is.na(match(targets, genes)))
    ms_targ <- targets[ms_targ_ind]
    
    # zeros across all cells in the time point
    dropout_ind <- which(!(rowSums(data) > 0))
    dropout_genes <- genes[dropout_ind]
    
    pathway_elements <- unique(c(lignames, recnames, targets))
    
    type_vec <- c()
    dropvec <- c()
    missingvec <- c()
    for (i in 1:length(pathway_elements)){
      elem <- pathway_elements[i]
      if(elem %in% lignames & elem %in% targets){
        type_vec[i] <- "lig, targ"
      } else if(elem %in% recnames & elem %in% targets){
        type_vec[i] <- "rec, targ"
      } else if(elem %in% targets){
        type_vec[i] <- "targ"
      } else if(elem %in% recnames){
        type_vec[i] <- "rec"
      } else {
        type_vec[i] <- "lig"
      }
      
      if(elem %in% dropout_genes){
        dropvec[i] <- "dropout"
      } else {
        dropvec[i] <- "observed"
      }
      
      if(elem %in% genes){
        missingvec[i] <- "observed"
      }else {
        missingvec[i] <- "missing"
      }
    }
    
    
    # see function description for these table elements
    removed_table <- cbind(pathway_elements, type_vec, dropvec, missingvec)
    colnames(removed_table) <- c("gene name", "gene type", "dropout in subset", "missing from all subsets")
    
    # remove the excluded markers from the pathway table
    msl <- which(pathway$ligand %in% ms_lig)
    msr <- which(pathway$receptor %in% ms_rec)
    mst <- which(pathway$target %in% ms_targ)
    
    if(length(msl) > 0 |  length(msr) > 0 | length(mst) > 0){
      pathway_removed <- pathway[-c(msl, msr, mst),]
    } else {
      pathway_removed <- pathway
    }
  } else {
    # make gene names systematic
    tableLig$ligand <- unname(sapply(tableLig$ligand, SimpleCap))
    tableLig$receptor <- unname(sapply(tableLig$receptor, SimpleCap))
    genes <- unname(sapply(genes, SimpleCap))
    
    # find all the symbols in the pathway
    lignames <- unique(tableLig$ligand)
    recnames <- unique(tableLig$receptor)
    
    # inner join the tables
    pathway <- tableLig
    
    # pathway elements for which we have no data
    ms_lig_ind <- which(is.na(match(lignames, genes)))
    ms_lig <- lignames[ms_lig_ind]
    
    ms_rec_ind <- which(is.na(match(recnames, genes)))
    ms_rec <- recnames[ms_rec_ind]
    
    # zeros across all cells in the time point
    dropout_ind <- which(!(rowSums(data) > 0))
    dropout_genes <- genes[dropout_ind]
    
    pathway_elements <- unique(c(lignames, recnames))
    
    type_vec <- c()
    dropvec <- c()
    missingvec <- c()
    for (i in 1:length(pathway_elements)){
      elem <- pathway_elements[i]
      if(elem %in% recnames){
        type_vec[i] <- "rec"
      } else {
        type_vec[i] <- "lig"
      }
      if(elem %in% dropout_genes){
        dropvec[i] <- "dropout"
      } else {
        dropvec[i] <- "observed"
      }
      if(elem %in% genes){
        missingvec[i] <- "observed"
      }else {
        missingvec[i] <- "missing"
      }
    }
    # see function description for these table elements
    removed_table <- cbind(pathway_elements, type_vec, dropvec, missingvec)
    colnames(removed_table) <- c("gene name", "gene type", "dropout in subset", "missing from all subsets")
    
    # remove the excluded markers from the pathway table
    msl <- which(pathway$ligand %in% ms_lig)
    msr <- which(pathway$receptor %in% ms_rec)
    if(length(msl) > 0 |  length(msr) > 0){
      pathway_removed <- pathway[-c(msl, msr),]
    } else {
      pathway_removed <- pathway
    }
  }
  
  return(list(lignames = lignames, 
              recnames = recnames, 
              targets = NULL,
              pathway = pathway,
              pathway_removed = pathway_removed, 
              removed = removed_table))
}



#' Capitalize the first letter of a string
#' 
#' Capitalize the first letter of a string.
#'
#' @param x the string to run
#'
#' @return a string
#' 
SimpleCap <- function(x) {
  s <- strsplit(as.character(x), " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}