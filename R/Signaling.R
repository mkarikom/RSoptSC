getK <- function(a,b){
  if(a == 0 && b == 0){
    return(0)
  }else{
    return(a/(a+b))
  }
}

getD <- function(a,g){
  if(a == 0 && g == 0){
    return(0)
  }else{
    return(a/(a+g))
  }
}

#' compute lig and rec
#'
#' @param i the ligand
#' @param j the receptor
#'                                          
#' @return a vector of probs
#' 
getAlpha <- function(i,j){
  if(i == 0 || j == 0){
    return(0)
  }else{
    return(exp(-1/(i*j)))
  }
}

#' compute uptargets
#'
#' @param Yv the target expression
#'                                          
#' @return a vector of probs
#' 
getBetaj <- function(Yv){
  vsum = sum(Yv)
  if(vsum == 0){
    return(0)
  }else{
    return(exp(-length(Yv)/vsum))
  }
}

#' compute downtargets
#'
#' @param Ystarv the target expression
#'                                          
#' @return a vector of probs
#' 
getGammaj <- function(Ystarv){
  vsum = sum(Ystarv)
  if(vsum == 0){
    return(0)
  }else{
    return(exp(-vsum/length(Ystarv)))
  }
}

#' compute the unnormalized probability
#'
#' @param ligname the ligand
#' @param recname the receptor
#' @param beta the up targets
#' @param gamma the down targets
#' @param flat the flattened p matrix
#' @param data the expression
#' 
#'                                          
#' @return a vector of probs
#' 
#' @import doParallel foreach                                           
updateFlatBoth <- function(ligname,recname,beta,gamma,flat,data){
  cat("\n using parallel update")
  flatval = foreach(v=1:dim(flat)[1], .export=c("getAlpha","getK","getD"),  .combine = 'c') %dopar% {
    ligcelli = data[tolower(ligname),flat[v,1]]
    reccellj = data[tolower(recname),flat[v,2]]
    a = getAlpha(ligcelli,reccellj)
    b = beta[flat[v,2]]
    g = gamma[flat[v,2]]
    kappa = getK(a,b)
    delta = getD(a,g)
    a * kappa * delta * b * g
  }
  return(flatval)
}

#' compute the unnormalized probability
#'
#' @param ligname the ligand
#' @param recname the receptor
#' @param beta the up targets
#' @param flat the flattened p matrix
#' @param data the expression
#'                                          
#' @return a vector of probs
#' 
#' @import doParallel foreach                                           
updateFlatUp <- function(ligname,recname,beta,flat,data){
  cat("\n using parallel update")
  flatval = foreach(v=1:dim(flat)[1], .export=c("getAlpha","getK"),  .combine = 'c') %dopar% {
    
    ligcelli = data[tolower(ligname),flat[v,1]]
    reccellj = data[tolower(recname),flat[v,2]]
    
    a = getAlpha(ligcelli,reccellj)
    b = beta[flat[v,2]]
    kappa = getK(a,b)
    a * kappa * b
  }
  return(flatval)
}

#' compute the unnormalized probability
#'
#' @param ligname the ligand
#' @param recname the receptor
#' @param gamma the down targets
#' @param flat the flattened p matrix
#' @param data the expression
#'                                          
#' @return a vector of probs
#' 
#' @import doParallel foreach                                           
updateFlatDown <- function(ligname,recname,gamma,flat,data){
  cat("\n using parallel update")
  flatval = foreach(v=1:dim(flat)[1], .export=c("getAlpha","getD"),  .combine = 'c') %dopar% {
    
    ligcelli = data[tolower(ligname),flat[v,1]]
    reccellj = data[tolower(recname),flat[v,2]]
    
    a = getAlpha(ligcelli,reccellj)
    g = gamma[flat[v,2]]
    delta = getD(a,g)
    a * delta * g
  }
  return(flatval)
}

#' compute the unnormalized probability
#'
#' @param ligname the ligand
#' @param recname the receptor
#' @param flat the flattened p matrix
#' @param data the expression
#'                                          
#' @return a vector of probs
#' 
#' @import doParallel foreach                                   
updateFlatNone <- function(ligname,recname,flat,data){
  cat("\n using parallel update")
  flatval = foreach(v=1:dim(flat)[1], .export=c("getAlpha"),  .combine = 'c') %dopar% {
    
    ligcelli = data[tolower(ligname),flat[v,1]]
    reccellj = data[tolower(recname),flat[v,2]]
    
    a = getAlpha(ligcelli,reccellj)
    a
  }
  return(flatval)
}

#' Compute Cell-cell interaction probability
#'
#' We can have a situation where ligand/receptor expression is low, but target genes (repressors and activators) are highly expressed.  In this case we will get a false positive for \eqn{P_{i,j}}.  In order to correct for this, we introduce normalizing coefficients for the relationship between \eqn{\alpha_{i,j}} and \eqn{\beta / \gamma}.  When \eqn{\alpha} is low, then K and D will decrease rapidly with increasing \eqn{\beta} and \eqn{\gamma}, penalizing the resulting increase in P (which is increasingly likely to be a false positive, given the disparity between \eqn{\alpha} and \eqn{\beta},\eqn{\gamma})
#' \eqn{K_{i,j} = \frac{\alpha_{i,j}}{\alpha_{i,j} + \beta_{i,j}}}
#' Note that activated genes are required under this model.
#'
#' @param data a matrix of expression values for each cell (rows) and gene (columns)
#' @param ids a vector of gene ids
#' @param pathway a data frame with "ligands", "receptors", "direction", and "targets"
#' @param normalize_aggregate whether or not to normalize the P_agg matrix by dividing each row by its sum default is true
#'                                          
#' @return a list containing:
#'     \item{P}{a list of the cell-cell signaling probabilities for all ligand/receptor pairs}
#'     \item{P_agg}{the aggregate matrix of all ligand/receptor pairs}
#'                          
#' @import doParallel                                              
#' @importFrom reshape2 melt
#' @importFrom dplyr filter mutate_all
#' @importFrom Matrix as.matrix rowSums colMeans
#' 
#' @export
#'
GetSignalingPartners <- function(data,
                                 ids,
                                 pathway,
                                 normalize_aggregate = TRUE){
  colnames(data) = tolower(colnames(data))
  rownames(data) = tolower(rownames(data))
  RCpairs = dplyr::distinct(pathway$pathway_removed, receptor, ligand)
  P = list()
  LR = list()
  
  chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  
  if (nzchar(chk) && chk == "TRUE") {
    # use 2 cores in CRAN/Travis/AppVeyor
    ncores <- 1
  } else {
    # use all cores in devtools::test()
    ncores <- detectCores()
  }
  doParallel::registerDoParallel(ncores)
  
  
  
  # loop over the RCpairs list
  for(u in 1:dim(RCpairs)[1]){
    # add if/else to switch between these
    ligname = RCpairs$ligand[u]
    recname = RCpairs$receptor[u]
    dirs = dplyr::select(pathway$pathway_removed,receptor,ligand,target,direction) %>% 
      dplyr::filter(receptor==recname,ligand==ligname) %>% 
      .$direction %>% 
      unique() 
    ################################################
    # there are no targets
    ################################################
    if(length(dirs)==0){
      cat("\n no targets detected, directions = ", dirs, ", all(dirs=='up')=", all(dirs == "up"))
      flat = reshape2::melt(matrix(0,ncol(data),ncol(data)))
      flat = flat[,-3]
      flat$value = updateFlatNone(ligname,recname,flat,data)
    }else if(all(dirs=="up")){
      cat("\n only up-targets detected, directions = ", dirs, ", all(dirs=='up')=", all(dirs == "up"))
      ################################################
      # there are only up targets
      ################################################
      flat = reshape2::melt(matrix(0,ncol(data),ncol(data)))
      flat = flat[,-3]
      
      uptargs = dplyr::select(pathway$pathway_removed,receptor,target,direction) %>% 
        dplyr::filter(receptor == recname,direction=="up") %>% 
        dplyr::select(target) %>% 
        .$target %>% unique()
      
      bvec = rep(0,ncol(data))
      for(b in 1:length(bvec)){
        bvec[b] = getBetaj(data[tolower(uptargs),b]) 
      }
      flat$value = updateFlatUp(ligname,recname,bvec,flat,data)
      
    }else if(all(dirs == "down")){
      cat("\n only down-targets detected, directions = ", dirs, ", all(dirs=='up')=", all(dirs == "up"))
      ################################################
      # there are only down targets
      ################################################
      flat = reshape2::melt(matrix(0,ncol(data),ncol(data)))
      flat = flat[,-3]
      
      downtargs = dplyr::select(pathway$pathway_removed,receptor,target,direction) %>% 
        dplyr::filter(receptor == recname,direction=="down") %>% 
        dplyr::select(target) %>% 
        .$target %>% unique()
      
      gvec = rep(0,ncol(data))
      
      for(g in 1:length(gvec)){
        gvec[g] = getGammaj(data[tolower(downtargs),g]) 
      }
      flat$value = updateFlatDown(ligname,recname,gvec,flat,data)
    }else{
      cat("\n up and down targets detected, directions = ", dirs, ", all(dirs=='up')=", all(dirs == "up"))
      ################################################
      # there are uptargets and downtargets
      ################################################
      flat = reshape2::melt(matrix(0,ncol(data),ncol(data)))
      flat = flat[,-3]
      
      uptargs = dplyr::select(pathway$pathway_removed,receptor,target,direction) %>% 
        dplyr::filter(receptor == recname,direction=="up") %>% 
        dplyr::select(target) %>% 
        .$target %>% unique()
      
      bvec = rep(0,ncol(data))
      for(b in 1:length(bvec)){
        bvec[b] = getBetaj(data[tolower(uptargs),b]) 
      }
      
      downtargs = dplyr::select(pathway$pathway_removed,receptor,target,direction) %>% 
        dplyr::filter(receptor == recname,direction=="down") %>% 
        dplyr::select(target) %>% 
        .$target %>% unique()
      
      gvec = rep(0,ncol(data))
      for(g in 1:length(gvec)){
        gvec[g] = getGammaj(data[tolower(downtargs),g]) 
      }
      flat$value = updateFlatBoth(ligname,recname,bvec,gvec,flat,data)
    }
    
    # normalize and store the values
    LR[[u]] = list("lig"=ligname,"rec"=recname)
    P[[u]] = flat %>% dplyr::group_by(Var1) %>% dplyr::mutate(value = value/sum(value))
    P[[u]]$value[is.nan(P[[u]]$value)] = 0
    P[[u]]$value = as.numeric(P[[u]]$value) # convert the normalized values to lower precision
  }
  P_tot = P[[1]]
  P_tot$value = rep(0,length(P_tot$value))
  for(i in 1:length(P)){
    P_tot$value = P_tot$value + P[[i]]$value
  }
  P_tot = P_tot %>% dplyr::group_by(Var1) %>% dplyr::mutate(value = value/sum(value))
  P_tot$value[is.nan(P_tot$value)] = 0
  return(list("P"=P,"P_tot"=P_tot,"LR"=LR))
}


#' Average cluster to cluster signaling.
#'
#' @param Pcell signaling probabilities cells x cells
#' @param cluster_labels labels of cells 1:n
#' @param normalize_rows normalize the rows of the output matrix so that the rows sum to 1, default true
#'
#' @return a matrix of cluster to cluster signaling
#'
#' @export
#'
ClusterSig <- function(Pcell,
                       cluster_labels,
                       normalize_rows = FALSE){
  PClust = list()
  P = Pcell$P
  for(i in 1:length(P)){

    P_i = dplyr::inner_join(P[[i]],
                            data.frame(id=1:length(cluster_labels),
                                       cluster=cluster_labels), 
                            by = c("Var1" = "id")) %>% dplyr::rename(cluster.Var1 = cluster)
    P_i = dplyr::inner_join(P_i,
                            data.frame(id=1:length(cluster_labels),
                                       cluster=cluster_labels), 
                            by = c("Var2" = "id")) %>% dplyr::rename(cluster.Var2 = cluster)
    P_i = P_i %>% filter(value > 0)
    if(!normalize_rows){
      PClust[[i]] = dplyr::group_by(P_i,cluster.Var1, cluster.Var2) %>% dplyr::summarize(value = sum(value))
    }else{
      PClust[[i]] = dplyr::group_by(P_i,cluster.Var1, cluster.Var2) %>% dplyr::summarize(value = sum(value)/dplyr::n_distinct(Var1))
    }
    PClust[[i]]$cluster.Var1 = as.factor(PClust[[i]]$cluster.Var1)
    PClust[[i]]$cluster.Var2 = as.factor(PClust[[i]]$cluster.Var2)
  }
  P_tot = PClust[[1]]
  colnames(P_tot) = c("cluster.Var1","cluster.Var2","RC1")
  for(i in 2:length(PClust)){
    colnames(PClust[[i]]) = c("cluster.Var1","cluster.Var2",paste0("RC",i))
    P_tot = full_join(P_tot, PClust[[i]])
  }
  tmp = as.matrix(P_tot[3:dim(P_tot)[2]])
  ind = which(is.na(tmp))
  tmp[ind] = 0
  P_tot$value = rowSums(tmp)
  P_tot = P_tot[,c("cluster.Var1","cluster.Var2","value")]
  
  P_tot = P_tot %>% dplyr::group_by(cluster.Var1) %>% dplyr::mutate(value = value/sum(value))
  return(list("P"=PClust,"P_tot"=P_tot))
}

#' Average cluster to cluster signaling.
#'
#' @param Pcell signaling probabilities cells x cells
#' @param cluster_labels labels of cells 1:n
#' @param normalize_rows normalize the rows of the output matrix so that the rows sum to 1, default true
#'
#' @return a matrix of cluster to cluster signaling
#'
#' @export
#'
ClusterSig_single <- function(Pcell,
                       cluster_labels,
                       normalize_rows = FALSE){

    P_i = dplyr::inner_join(Pcell,
                            data.frame(id=1:length(cluster_labels),
                                       cluster=cluster_labels), 
                            by = c("Var1" = "id")) %>% dplyr::rename(cluster.Var1 = cluster)
    P_i = dplyr::inner_join(P_i,
                            data.frame(id=1:length(cluster_labels),
                                       cluster=cluster_labels), 
                            by = c("Var2" = "id")) %>% dplyr::rename(cluster.Var2 = cluster)
    P_i = P_i %>% filter(value > 0)
    if(!normalize_rows){
      P_i = dplyr::group_by(P_i,cluster.Var1, cluster.Var2) %>% dplyr::summarize(value = sum(value))
    }else{
      P_i = dplyr::group_by(P_i,cluster.Var1, cluster.Var2) %>% dplyr::summarize(value = sum(value)/dplyr::n_distinct(Var1))
    }
    P_i$cluster.Var1 = as.factor(P_i$cluster.Var1)
    P_i$cluster.Var2 = as.factor(P_i$cluster.Var2)
  return(P_i)
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