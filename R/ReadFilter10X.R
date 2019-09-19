#' Read 10X Data
#' 
#' Read data from 10X counts or per-molecule info.
#'
#' @param filepath full path of the data directory
#' @param inf_exp_rate the minimum mean expression rate of retained genes
#' @param inf_gene_count retain a gene if more than this number of cells express it
#' @param inf_cell_count retain a cell if more than this number of genes are expressed in it
#' @param annotation_db annotation package for the data, eg "org.Mm.eg.db"
#' @param feature_annotation_key the key for annotation.db lookups, eg "SYMBOL"
#' @param remove_spike_ERCC remove ERCC spike ins after they are used for normalization
#' @param pre_clust whether or not to precluster the cells before getting size factors
#' @param rescale whether to compute size factors and rescale counts
#' @param subsample randomly pick n cells
#' @param filtercells whether to do QC in scater
#' @param var_gene_method how to select variable genes, eg 'mean_var', 'pca'.  mean_var selects the genes whose variance is highest after removing the fitted spike-in-based technical variance. pca finds the genes with the max rotation coefficients 
#' @param ... arguments to called functions
#' @return a list containing:
#'     \item{c(seur, sce)}{a seurat or sce object}
#'
#' @importFrom SingleCellExperiment isSpike isSpike<- colData rowData SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment rowData<- colData<- assay assay<-
#' @importFrom AnnotationDbi mapIds
#' @importFrom grDevices cairo_ps
#' @importFrom scater calculateQCMetrics isOutlier calcAverage nexprs normalize
#' @importFrom scran computeSumFactors computeSpikeFactors trendVar decomposeVar quickCluster
#' @importFrom DropletUtils read10xCounts 
#'  
#' @export
#' 
ReadFilter10X <- function(filepath = NULL, 
                    inf_exp_rate = 0,
                    inf_gene_count = 0,
                    inf_cell_count = 0,
                    annotation_db,
                    feature_annotation_key,
                    remove_spike_ERCC = TRUE,
                    pre_clust = FALSE,
                    rescale = TRUE,
                    subsample = NULL,
                    filtercells = TRUE,
                    var_gene_method = 'mean_var',
                    ...) {
  x = NULL
  # get metadata
  sce <- read10xCounts(filepath)
  sce.orig <- sce
  
  if(!is.null(subsample)){
    sce <- sce[,sample.int(n = ncol(sce), size = subsample, replace=FALSE)]
  }

  if(feature_annotation_key == "SYMBOL"){
    rowData(sce)$SYMBOL <- rownames(sce)
    ENSEMBL <- mapIds(get(annotation_db), keys=rownames(sce),
                      column="ENSEMBL", keytype=feature_annotation_key)
    rowData(sce)$ENSEMBL <- ENSEMBL
    
    # find and add chromosomal location
    location <- mapIds(get(annotation_db), keys=rownames(sce),
                       column="CHR", keytype=feature_annotation_key)
    rowData(sce)$CHR <- location
    
    # add spike in metadata
    isSpike(sce, "ERCC") <- grepl("^ERCC", rowData(sce)$ENSEMBL)
  }else if(feature_annotation_key == "ENSEMBL"){
    rowData(sce)$ENSEMBL <- rownames(sce)
    
    # find and add symbols
    SYMBOL <- mapIds(get(annotation_db), keys=rownames(sce),
                      column="SYMBOL", keytype=feature_annotation_key)
    rowData(sce)$ENSEMBL <- rownames(sce)
    
    # for missing symbols, use the Ensembl id
    SYMBOL[which(is.na(SYMBOL))] <- names(which(is.na(SYMBOL)))
    
    # find and add chromosomal location
    location <- mapIds(get(annotation_db), keys=rownames(sce),
                       column="CHR", keytype=feature_annotation_key)
    rowData(sce)$CHR <- location
    
    rownames(sce) <- SYMBOL
    isSpike(sce, "ERCC") <- grepl("^ERCC", rowData(sce)$ENSEMBL)
  }


  if(filtercells){
    # find mito genes
    mito <- which(rowData(sce)$CHR=="MT")
    sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
    head(colnames(colData(sce)), 10)
    
    libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
                              log=TRUE)
    feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", 
                              log=TRUE)
    mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher")
    
    #keep <- !(libsize.drop | feature.drop | mito.drop)
    keep <- !(mito.drop)
    
    # finally remove all the cells that did not pass QC
    sce$PassQC <- keep
    sce <- sce[,keep]
  }
  
  # filter out lowly expressed genes
  # this removes unexpressed genes from sce, and creates a second filtered where avg expressio for all genes is > 1
  
  num.cells <- nexprs(sce, byrow=TRUE)
  num.genes <- nexprs(sce, byrow=FALSE)
  to.keep.gene <- num.cells > inf_gene_count
  to.keep.cell <- num.genes > inf_cell_count
  sce <- sce[to.keep.gene,to.keep.cell]
  
  ave.counts <- calcAverage(sce, use_size_factors=FALSE)
  to.keep <- ave.counts >= inf_exp_rate
  sce <- sce[to.keep,]


  
  # normalization
  if(pre_clust){
    clusters <- quickCluster(sce, use.ranks=FALSE)
  }else{
    clusters <- NULL
  }

  if(rescale){
    sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
    sce <- normalize(sce)
  }else{
    assay(sce, 'logcounts') <- log2(counts(sce) + 1)
  }

  var.fit <- trendVar(sce, parametric=TRUE,loess.args=list(span=0.3),use.spikes=FALSE)
  var.out <- decomposeVar(sce, var.fit)
  dec.bio <- order(var.out$bio, decreasing=TRUE)
  dec.bio.genes <- rownames(sce)[dec.bio]

  # remove spike ins prior to downstream analysis
  if(remove_spike_ERCC){
    ind <- which(isSpike(sce))
    if(length(ind) > 0){
      sce <- sce[-ind,]
    }
  }
  
  # variable gene selection
  if(var_gene_method == 'mean_var'){
    chosen.genes.names <- dec.bio.genes
    chosen.genes.ind <- match(dec.bio.genes, rownames(sce))
    chosen.genes.ind <- chosen.genes.ind[!is.na(chosen.genes.ind)]
    chosen.genes <- list(names = chosen.genes.names, ind = chosen.genes.ind)
    plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
         ylab="Variance of log-expression")
    curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
  }else if(var_gene_method == 'varimax'){
    pr <- prcomp(t(logcounts(sce)))
    coef <- pr$rotation[,1:5]
    max_coef <- apply(coef, 1, max)
    ord <- order(max_coef, decreasing = TRUE)
    chosen.genes.ind <- ord
    chosen.genes.names <- rownames(sce)[ord]
    chosen.genes <- list(names = chosen.genes.names, ind = chosen.genes.ind)
  }
  
  return(list(orig = sce.orig, qc = sce, 
              vargenes = chosen.genes))
}