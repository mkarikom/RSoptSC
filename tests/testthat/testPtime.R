context("decomposition values")
library(RSoptSC)


test_that("root cell is detected when flat dist used", {
  root_cell <- FindRootCell(dist_flat = RSoptSC::GuoPtime$Values$low_dis)
  expect_equal(root_cell, RSoptSC::GuoPtime$Values$root_cell0)
})

test_that("root cell is detected when graph dist used", {
  check_cluster <- igraph::graph_from_adjacency_matrix(adjmatrix = RSoptSC::GuoPtime$Values$CC_adjacent,
                                                       mode = "upper",
                                                       weighted = TRUE)

  root_cell <- FindRootCell(use_flat_dist = FALSE,
                            cluster_order_by = "predecessor",
                            cell_order_by = "index",
                            graph_cluster = check_cluster,
                            dist_graph = RSoptSC::GuoPtime$Values$Short_pathd,
                            dist_flat = RSoptSC::GuoPtime$Values$low_dis,
                            cluster_labels = RSoptSC::GuoPtime$Params$cluster_label,
                            root_cluster = RSoptSC::GuoPtime$Values$root_cluster)
  browser()
  expect_equal(root_cell, RSoptSC::GuoPtime$Values$root_cell)
})

test_that("root cluster is detected", {
  root_cluster <- FindRootCluster(cluster_labels = RSoptSC::GuoPtime$Params$cluster_label,
                                  flat_embedding = RSoptSC::GuoPtime$Params$latent,
                                  dist_graph = RSoptSC::GuoPtime$Values$Short_pathd,
                                  dist_flat = RSoptSC::GuoPtime$Values$low_dis,
                                  reverse = FALSE)
  expect_equal(root_cluster$cluster_adj_matrix, RSoptSC::GuoPtime$Values$CC_adjacent)
  expect_equal(root_cluster$root_cluster, RSoptSC::GuoPtime$Values$root_cluster)
})

test_that("added inter-component edges correctly", {
  connected_graph <- JoinGraphComponents(root_cell = RSoptSC::GuoPtime$Values$root_cell0,
                                         adjacency_matrix = RSoptSC::GuoPtime$Values$W_graph1,
                                         flat_distances = RSoptSC::GuoPtime$Values$low_dis,
                                         n_components = RSoptSC::GuoPtime$Values$nComponents,
                                         component_members = RSoptSC::GuoPtime$Values$members)
  expect_equal(connected_graph, RSoptSC::GuoPtime$Values$W_graph)
})

test_that("mapping is correct, when components have been previously joined", {
  mapping <- RepresentationMap(low_dim_embedding = RSoptSC::GuoPtime$Params$latent,
                               similarity_matrix = RSoptSC::GuoPtime$Params$W,
                               join_components = FALSE)
  mem <- RSoptSC::GuoPtime$Values$members
  mem <- lapply(mem, function(x) sort(x))
  mem <- mem[order(sapply(mem, length), decreasing=FALSE)]

  expect_equal(mapping$dist_flat, RSoptSC::GuoPtime$Values$low_dis, check.attributes = FALSE)
  expect_equal(mapping$dist_graph, RSoptSC::GuoPtime$Values$Short_pathd1)
  expect_equal(mapping$adj_matrix, RSoptSC::GuoPtime$Values$W_graph1)
  expect_equal(mapping$n_components, RSoptSC::GuoPtime$Values$nComponents)
  expect_equal(mapping$sizes, RSoptSC::GuoPtime$Values$sizes)
  expect_equal(mapping$members, mem)

})
