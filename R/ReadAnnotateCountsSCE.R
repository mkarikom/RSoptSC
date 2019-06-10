#' Read Counts Data
#' 
#' Read counts data from file and annotate the generated Single Cell Experiment Object.  Right now this function assumes:
#' - mouse cells
#' - only spike ins are ERCC
#' - row and column names are provided in the data
#' - the row names are gene symbols
#' In the case of GEO data, it takes a GSE accession code and assumes that the supp files contain all the expression data
#'
#' @param filepath the full or relative (to the workspace) path of the data file
#' @param inf_exp_rate the minimum mean expression rate of retained genes
#' @param inf_gene_count retain a gene if more than this number of cells express it
#' @param inf_cell_count retain a cell if more than this number of genes are expressed in it
#' @param series_accession the series accession code, eg "GSE67602"
#' @param metad_title the series accession code, eg "GSE67602"
#' @param cell_id the series accession code, eg "GSE67602"
#' @param GEOfilecache location for caching series supplemental files, default is current directory
#' @param annotation_db annotation package for the data, eg "org.Mm.eg.db"
#' @param feature_annotation_key the key for annotation.db lookups, eg "SYMBOL"
#' @param pheno_filter a vector of strings to partially match cell metadata 
#' @param plot_title path of the saved mean/variance plot, default NULL, if NULL the plot will be printed instead
#' @param remove_spike_ERCC remove ERCC spike ins after they are used for normalization
#' @param pre_clust whether or not to precluster the cells before getting size factors
#' @param rescale whether to compute size factors and rescale counts
#' @param var_gene_method how to select variable genes, eg 'mean_var', 'pca'.  mean_var selects the genes whose variance is highest after removing the fitted spike-in-based technical variance. pca finds the genes with the max rotation coefficients 
#' @param ... arguments to called functions
#' @return a list containing:
#'     \item{c(seur, sce)}{a seurat or sce object}
#'
#' @importFrom SingleCellExperiment isSpike isSpike<- colData rowData SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment rowData<- colData<- assay assay<-
#' @importFrom AnnotationDbi mapIds
#' @importFrom GEOquery getGEO getGEOSuppFiles
#' @importFrom grDevices cairo_ps
#' @importFrom scater calculateQCMetrics isOutlier calcAverage nexprs normalize
#' @importFrom scran computeSumFactors computeSpikeFactors trendVar decomposeVar quickCluster
#'  
#' @export
#' 
ReadAnnotateCountsSCE <- function(filepath = NULL, 
                        inf_exp_rate = 0,
                        inf_gene_count = 0,
                        inf_cell_count = 0,
                        series_accession,
                        metad_title = 'characteristics_ch1',
                        cell_id = 'title',
                        pheno_filter = NULL,
                        plot_title = NULL,
                        GEOfilecache = getwd(),
                        annotation_db,
                        feature_annotation_key,
                        remove_spike_ERCC = FALSE,
                        pre_clust = FALSE,
                        rescale = TRUE,
                        var_gene_method = 'mean_var',
                        ...) {

  x = NULL
  # get metadata
  series <- getGEO(series_accession, GSEMatrix = TRUE) # 'GSE67602'
  pheno <- cbind(pData(series[[1]])[eval(cell_id)],
                 pData(series[[1]])[eval(metad_title)])
    
  # get the expression data
  filepaths = rownames(getGEOSuppFiles(series_accession, baseDir = GEOfilecache))
  counts <- read.csv(filepaths[1], sep="", row.names = 1)
  #strip leading X's from the cell names
  names <- unlist(lapply(strsplit(colnames(counts), split = 'X'), function(x){
    x[2]
  }))
  
  sce <- SingleCellExperiment(list(counts=as.matrix(counts)))
  sce$Phenotype <- pheno[,eval(metad_title)]
  if(!is.null(pheno_filter)){
    matches <- unique(grep(paste(pheno_filter,collapse="|"),
                           sce$Phenotype,value=FALSE))
    sce <- sce[,matches]
  }
  
  if(feature_annotation_key == "SYMBOL"){
    rowData(sce)$SYMBOL <- rownames(sce)
    ENSEMBL <- mapIds(get(annotation_db), keys=rownames(sce),
                      column="ENSEMBL", keytype=feature_annotation_key)
    rowData(sce)$ENSEMBL <- ENSEMBL
    
    # add spike in metadata
    isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
  }else if(feature_annotation_key == "ENSEMBL"){
    rowData(sce)$ENSEMBL <- rownames(sce)
    
    SYMBOL <- mapIds(get(annotation_db), keys=rownames(sce),
                      column="SYMBOL", keytype=feature_annotation_key)
    rowData(sce)$SYMBOL <- SYMBOL
    rownames(sce) <- SYMBOL
    isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
  }

  
  # find location of genes and anotate mito genes
  # location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rowData(sce)$SYMBOL,
  location <- mapIds(get(annotation_db), keys=rownames(sce),
                     column="CHR", keytype=feature_annotation_key)
  rowData(sce)$CHR <- location
  
  
  mito <- which(rowData(sce)$CHR=="MT")
  sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
  head(colnames(colData(sce)), 10)
  
  libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
                            log=TRUE, batch=sce$PlateOnco)
  feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", 
                            log=TRUE, batch=sce$PlateOnco)
  
  spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher",
                          batch=sce$PlateOnco)
  
  keep <- !(libsize.drop | feature.drop | spike.drop)
  
  # finally remove all the cells that did not pass QC
  sce$PassQC <- keep
  sce.orig <- sce
  sce <- sce[,keep]
  
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


  
  # spike in based normalization
  if(pre_clust){
    clusters <- quickCluster(sce, use.ranks=FALSE)
  }else{
    clusters <- NULL
  }

  if(rescale){
    sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
    sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
    sce <- normalize(sce)
  }else{
    assay(sce, 'logcounts') <- log2(counts(sce) + 1)
  }

  var.fit <- trendVar(sce, parametric=TRUE,loess.args=list(span=0.3))
  var.out <- decomposeVar(sce, var.fit)
  dec.bio <- order(var.out$bio, decreasing=TRUE)
  dec.bio.genes <- rownames(sce)[dec.bio]

  # remove spike ins prior to downstream analysis
  if(remove_spike_ERCC){
    ind <- which(rowData(sce)$is_feature_control_ERCC)
    sce <- sce[-ind,]
  }
  
  # variable gene selection
  if(var_gene_method == 'mean_var'){
    chosen.genes.names <- dec.bio.genes
    chosen.genes.ind <- match(dec.bio.genes, rownames(sce))
    chosen.genes.ind <- chosen.genes.ind[!is.na(chosen.genes.ind)]
    chosen.genes <- list(names = chosen.genes.names, ind = chosen.genes.ind)
    
    if(!is.null(plot_title)){
      cairo_ps(plot_title, width = 5, height = 4)
      plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
           ylab="Variance of log-expression")
      curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
      cur.spike <- isSpike(sce)
      points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
      dev.off()
    } else {
      plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
           ylab="Variance of log-expression")
      curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
      cur.spike <- isSpike(sce)
      points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
    }
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