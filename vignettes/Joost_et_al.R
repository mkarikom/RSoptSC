## ----setup, include = FALSE----------------------------------------------
rm(list = ls())
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(13)
par(mar=c(1,1,1,1))

## ---- results='asis'-----------------------------------------------------
library(RSoptSC)
data("GSE67602_Joost")
logdata <- log10(GSE67602_Joost$data + 1)

## ------------------------------------------------------------------------
# clean the gene names
gene_names <- GSE67602_Joost$gene_names
spikein <- grep('ERCC', gene_names)
gene_names <- gene_names[-spikein]

# clean the data
logdata <- logdata[-spikein,]

## ---- results='asis'-----------------------------------------------------
gene_expression_threshold <- 0.03
n_features <- 3000
filtered_data <- SelectData(logdata, gene_expression_threshold, n_features)

## ---- results='asis'-----------------------------------------------------
S <- SimilarityM(lambda = 0.05, data = filtered_data$M_variable)

## ---- results='asis'-----------------------------------------------------
low_dim_mapping <- RepresentationMap(similarity_matrix = S$W,
                            join_components = TRUE)

## ---- results='asis'-----------------------------------------------------
clusters_evalues <- CountClusters(data = S$W)
n_clusters <- clusters_evalues$n

## ---- results='asis'-----------------------------------------------------
par(mar=c(1,1,1,1))
plot(c(1:20), 
     clusters_evalues$eigs$val[1:20],
     xlab = NA,
     ylab = 'eigenvalues',
     main = 'Eigenvalues of the Graph Laplacian')

## ---- results='asis', warning=FALSE--------------------------------------
output_NMF <- NMF::nmf(x = S$W,
              rank = n_clusters,
              method = 'lee',
              seed = 'nndsvd',
              .options = 'nP');
H <- NMF::basis(output_NMF)

## ---- results='asis'-----------------------------------------------------
labels <- apply(H, 1, function(x){
  which(x == max(x))})

## ---- results='asis', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
# plot clusters
par(mar=c(1,1,1,1))
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                         feature = as.factor(labels),
                         title = "NMF Cluster Labeling",
                         subtitle = "t-SNE Embedding",
                         featurename = "Cluster ID")

## ---- results='asis', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
true_labels <- RSoptSC::GSE67602_Joost$annotation
# plot clusters
par(mar=c(1,1,1,1))
FeatureScatterPlot(flat_embedding = low_dim_mapping$flat_embedding,
                         feature = as.factor(true_labels),
                         title = "True Labeling",
                         subtitle = "t-SNE Embedding",
                         featurename = "Annotated Cell Types")

