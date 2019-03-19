context("pseudotime inference and lineage trajectories")
library(RSoptSC)
library(igraph)

test_that("lineage graph is correct", {
  low_dim_mapping <- RepresentationMap(similarity_matrix = joostTest$S,
                                       join_components = TRUE,
                                       input = 'data',
                                       random_state = 0,
                                       knn.repeat = 3,
                                       min_dist = 0.6,
                                       n_neighbors = 30)
  
  cluster_ptime <- FindRootCluster(cluster_labels = joostTest$labels,
                                   flat_embedding = low_dim_mapping$flat_embedding,
                                   dist_graph = low_dim_mapping$dist_graph,
                                   dist_flat = low_dim_mapping$dist_flat,
                                   reverse = TRUE)
  
  
  root_cell <- FindRootCell(use_flat_dist = FALSE,
                            cluster_order_by = "distance",
                            cell_order_by = "distance",
                            graph_cluster_mst = cluster_ptime$cluster_mst,
                            dist_graph  = low_dim_mapping$dist_graph,
                            dist_flat = low_dim_mapping$dist_flat,
                            cluster_labels = joostTest$labels,
                            root_cluster = cluster_ptime$root_cluster)
  
  
  cluster_predecessors <- GetPredecessors(cluster_ptime$cluster_mst, cluster_ptime$root_cluster)
  
  
  cluster_dtree <- GetDominatorTree(cluster_predecessors, cluster_ptime$graph_cluster)
  
  expect_true(identical_graphs(cluster_dtree, joostTest$lineage_graph))
})
