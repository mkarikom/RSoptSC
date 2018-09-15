Lineage_Ptime <- function(W,No_cluster,cluster_label,root_cluster = 0,root_cell = 0,latent,reverse){

  # create all necessary representations
  # if the similarity matrix based graph embedding of the cells is not connected,
  #   then we make it connected by calling JoinGraphComponents() from inside
  #   RepresentationMap()
  mapping <- RepresentationMap(low_dim_embedding = latent, similarity_matrix = W)

  root_cell <- FindRootCell(flat_distances = mapping$dist_flat)


  ## connect the components of mapping$adj_matrix
  ## need to redo mapping components that are generated from adj_matrix
    ## run mapping again
  ##

  ## generate adjacency matrix for clusters
  ## generate graph based on this adjacency matrix

  ## generate MST on cluster graph
    ## use the provided root cluster
      ## infer the root cell
        ## get average similarity matrix graph path length to the cells in each cluster,
        ##    from each cell in the root cluster
        ## normalize each row in the cell/cluster table by dividing by the max value for that row
        ## compute the kendall tau for each set of pairs
          ## measure the concordance between pairs
        ## max tau gives the root cell
    ## infer the root cluster
      ## if the root cell is given, use the cluster of that cell
        ## generate MST on cluster graph
      ## if the root cell is not given
        ## first, infer the root cluster
          ## take the two clusters farthest apart on the similarity matrix graph
          ## if \VAR{Reverse} == 0, then take the cluster with the lowest dispersion as the root cluster
          ## PROCEED AS ON LINE 20

}
