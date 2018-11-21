## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- results='asis'-----------------------------------------------------
library(RSoptSC)
data("GSE67602_Joost")
logdata <- log10(GSE67602_Joost$data + 1)

## ---- results='asis'-----------------------------------------------------
gene_expression_threshold <- 0.03
n_features <- 3000
filtered_data <- SelectData(logdata, gene_expression_threshold, n_features)

## ---- results='asis'-----------------------------------------------------
sim <- SimilarityM(lambda = 0.05, data = filtered_data$M_variable)

## ---- results='asis'-----------------------------------------------------
n_clusters <- CountClusters(data = sim$W)

## ---- results='asis'-----------------------------------------------------
H <- NNLM::nnmf(sim$W, n_clusters, rel.tol = 1e-6);

## ---- results='asis'-----------------------------------------------------
labels <- apply(H$W, 1, function(x){
  which(x == max(x))})

## ---- results='asis'-----------------------------------------------------
markers <- GetMarkerTable(counts_data = logdata,
                          cluster_labels = labels,
                          H = H$W,
                          gene_expression_threshold = gene_expression_threshold,
                          n_features = n_features)
heatplot <- MarkerHeatmap(logdata,labels,markers,n_markers = 5)


## ---- results='asis'-----------------------------------------------------
mapped <- RepresentationMap(similarity_matrix = sim$W,
                            join_components = TRUE)

## ---- results='asis'-----------------------------------------------------
cluster_ptime <- FindRootCluster(cluster_labels = labels,
                                 flat_embedding = mapped$flat_embedding,
                                 dist_graph = mapped$dist_graph,
                                 dist_flat = mapped$dist_flat,
                                 reverse = FALSE)

## ---- results='asis'-----------------------------------------------------
root_cell <- FindRootCell(cluster_order_by = "distance",
                           cell_order_by = "distance",
                           graph_cluster_mst = cluster_ptime$cluster_mst,
                           dist_graph  = mapped$dist_graph,
                           dist_flat = mapped$dist_flat,
                           cluster_labels = labels,
                           root_cluster = cluster_ptime$root_cluster)

## ---- results='asis', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
cluster_predecessors <- GetPredecessors(cluster_ptime$cluster_mst, cluster_ptime$root_cluster)
cluster_dtree <- GetDominatorTree(cluster_predecessors, cluster_ptime$graph_cluster)
PlotLineage(cluster_dtree)

## ---- results='asis', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
ptime_distance <- cluster_ptime$cluster_adj_matrix[cluster_ptime$root_cluster,]
ptime_distance <- 1 + ptime_distance/max(ptime_distance)
PlotLineage(cluster_dtree, node_color = ptime_distance)

## ---- results='asis', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
pseudotime <- mapped$dist_graph[root_cell,]

# plot pseudotime
FeatureScatterPlot(flat_embedding = mapped$flat_embedding,
                         feature = pseudotime,
                         title = "Pseudotime Labeling",
                         subtitle = "t-SNE Embedding")

## ---- results='asis', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
gene_index <- which(GSE67602_Joost$gene_names == 'Krt10')
labeling <- logdata[gene_index,]

# plot features
FeatureScatterPlot(flat_embedding = mapped$flat_embedding,
                    feature = labeling,
                    title = "Krt10 Labeling",
                    subtitle = "t-SNE Embedding")

## ---- results='asis', out.width = '100%', out.height = '100%', fig.width = 7, fig.heigt = 7----
feature_avg <- AvgFeatureExpression(cluster_labels = labels,
                                    M = logdata,
                                    feature = "Krt10",
                                    feature_list = GSE67602_Joost$gene_names)
PlotLineage(cluster_dtree, node_color = feature_avg)

