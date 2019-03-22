#' Plot the pseudotime ordering of clusters
#' 
#' Plot the pseudotime ordering of clusters using igraph tools.
#'
#' @param network an igraph network with weighted edges corresponding to the pseudotime distance
#' @param node_color an optional vector of colors for node labeling.  If null then the nodes are all colored black
#' @param alpha_color the gradient of the nodes
#' @param exaggerate whether to exponentiate then linearly scale the edge weights, thus emphasizing their differences
#' @param weight_scaler the weight scaler for edges
#' @param arrow_scaler the arrow width scaler for edges
#' 
#' @return nothing
#'
#' @importFrom igraph edges E plot.igraph E<-
#'
#' @export
#'
PlotLineage <- function(network, node_color = NULL, alpha_color = 1, 
                        exaggerate = TRUE,
                        weight_scaler = 0.2,
                        arrow_scaler = 0.8){
  set.seed(1)
  if(is.null(node_color)){
    node_color <- ColorHue(n = length(edges(network)[[1]][1]))
    node_color <- node_color$hex.1.n.
  }
  if(exaggerate){
    E(network)$width <-
      (exp(E(network)$weight)) * weight_scaler
    E(network)$arrow.width <- arrow_scaler * E(network)$weight
  }else{
    E(network)$width <-
      E(network)$weight
  }

  vertex_order <- as.numeric(names(edges(network)[[1]][1]))
  node_color <- node_color[vertex_order]
  plot.igraph(network, vertex.color = node_color)
}

#' Produce a plot of matlab DTree data and return the object
#' 
#' If an output dir and filename are provided, a plot will be saved, otherwise the function will just return the graph.  This is used for internal analysis to compare with matlab implementation.
#'
#' @param edge_table a numeric matrix whose rows are directed edges of a tree:
#'     col 1 is v1, col2 is v2, col3 is weight
#' @param predecessors a vector of tree predecessors such that pred[i] = the predecessor of i
#' @param outputdir the output directory, relative to getwd()
#' @param outputfile the output file
#'
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @importFrom igraph plot.igraph layout_with_kk
#' 
#' @return an igraph representation of the tree
#'
PlotMatlabDtree <- function(edge_table, predecessors, outputdir = NULL, outputfile = NULL){
  directed_edge_table <- ProcessMatlabDTree(edge_table, predecessors)
  directed_graph <- GetDGFromTable(directed_edge_table)
  if(!is.null(outputdir) && !is.null(outputfile)){
    file_path <- paste0(getwd(),
                        .Platform$file.sep,
                        outputdir,
                        .Platform$file.sep,
                        outputfile)
    pdf(file_path)
    plot(directed_graph, layout=layout_with_kk)
    dev.off()
  }
  return(directed_graph)
}

#' Produce a scatter plot of the cells
#'
#' Create a scatter plot of the data on a 2-dim ebedding colored by feature.  If discrete (factor) data is passed, the colorscale parameter must be provided, otherwise the plot will default to gradient shading. 
#'
#' @param flat_embedding a low dim embedding of cells
#' @param feature a scalar representation of feature
#' @param title the title of the plot
#' @param subtitle the subtitle of the plot
#' @param featurename the name of the feature to title the legend
#' @param colorscale color scale for discrete data
#'
#' @return a ggplot2 object
#' 
#' @import ggplot2
#'
#' @export
#'
FeatureScatterPlot <- function(flat_embedding,
                               feature,
                               title,
                               subtitle,
                               featurename,
                               colorscale = NULL){
  n_features <- length(unique(feature))

  p <- ggplot(flat_embedding,
         aes_(x = as.name(colnames(flat_embedding)[1]),
              y = as.name(colnames(flat_embedding)[2]),
              color = ~`feature`))
    
  
  if(is.null(colorscale)){
    p +
      geom_point() +
      labs(title = title, subtitle = subtitle, color = featurename) +
      theme_minimal()
  }else{
    p +
    geom_point() +
    labs(title = title, subtitle = subtitle, color = featurename) +
    theme_minimal() +
    scale_color_manual(values = colorscale)
  }
}


#' Heatmap of the top n markers
#' 
#' Plot the heatmap of the top n inferred markers for each cluster.
#'
#' @param data expression data with genes x cells
#' @param gene_names a vector of symbolic gene names corresponding to the rows in the data matrix
#' @param cluster_labels a vector of cluster labels corresponding to the columns in the data matrix
#' @param markers a table of markers
#' @param n_features the top n features to plot
#' 
#' @return nothing
#'
#' @import dplyr
#' @import RColorBrewer
#' @import gplots
#' @importFrom tibble as_tibble
#' @importFrom Matrix as.matrix
#'
#' @export
#'
PlotTopN <- function(data,
                 gene_names,
                 cluster_labels, 
                 markers,
                 n_features){
  clusterId = geneScore = NULL # r cmd check pass
  n_cells <- ncol(data)
  n_clusters <- length(unique(cluster_labels))
  sorted_cell <- sort.int(cluster_labels, index.return = TRUE)

  markers_table <- as.data.frame(markers)
  sortedmarkers <- arrange(markers_table, clusterId, desc(geneScore))
  sorted_gene_table <- as_tibble(sortedmarkers) %>% group_by(clusterId) %>% top_n(n_features, geneScore)
  sorted_gene <- sorted_gene_table$geneID

  plot_data <- data[sorted_gene, sorted_cell$ix]


  rlabs <- gene_names[sorted_gene]

  clabs <- c(rep(NA, n_cells))
  counts <- as.matrix(table(cluster_labels))
  accumulator <- matrix(1, n_clusters, n_clusters)*lower.tri(matrix(1, n_clusters, n_clusters), diag = TRUE)
  last_cells <- accumulator %*% counts
  offset <- floor(counts/2)
  labeled_columns <- last_cells - offset
  clabs[labeled_columns] <- c(1:n_clusters)

  p <- heatmap.2(as.matrix(plot_data),
                 col = rev(brewer.pal(11,"RdBu")),
                 trace = 'none',
                 dendrogram='none',
                 Rowv=FALSE, Colv=FALSE,
                 labCol = clabs,
                 labRow = rlabs,
                 srtCol = 0,
                 #cexCol =1,
                 #cexRow = .7,
                 #lhei = c(1,3.5),
                 key.title = NA,
                 scale = "row")
}

#' Heatmap of specific genes
#' 
#' Given a set of markers, plot their expression across clusters on a heatmap.
#'
#' @param data expression data with genes x cells
#' @param gene_names a vector of symbolic gene names corresponding to the rows in the data matrix
#' @param cluster_labels a vector of cluster labels corresponding to the columns in the data matrix
#' @param markers a vector of markers
#' 
#' @return nothing
#'
#' @import RColorBrewer
#' @import gplots
#'
#' @export
#'
PlotClusterExpression <- function(data,
                                  gene_names,
                                  cluster_labels, 
                                  markers){

  n_clusters <- length(unique(cluster_labels))
  sorted_cell <- sort.int(cluster_labels, index.return = TRUE)
  
  filtered <- ClusterAvg(M = data,
                         gene_names = gene_names,
                         cell_order = sorted_cell$ix,
                         cell_labels = cluster_labels,
                         gene_list = markers)
  
  p <- heatmap.2(log10(t(filtered) + 1),
                 col = rev(brewer.pal(11,"RdBu")),
                 trace = 'none',
                 dendrogram='none',
                 Rowv=FALSE, Colv=FALSE,
                 labCol = c(1:n_clusters),
                 labRow = markers,
                 srtCol = 0,
                 cexCol =1,
                 cexRow = .7,
                 lhei = c(1,2),
                 key.title = NA,
                 scale = "row")
}


#' Violin plot of gene expression
#' 
#' Produce a violin plot of gene expression.
#'
#' @param data a matrix with genes in rows and cells in columns
#' @param gene_names a vector of names corresponding to the rows in data
#' @param labels a vector of cluster assignments for the cells
#' @param gene_name the name of the gene to plot
#' @param colorscale optionally provided color scale
#'
#' @return nothing
#'
#' @import ggplot2
#'
#' @export
#'
ViolinPlotExpression <- function(data,
                                 gene_names,
                                 labels,
                                 gene_name,
                                 colorscale = NULL){
  n_features <- length(unique(labels))
  if(is.null(colorscale)){
    colorscale <- ColorHue(n = n_features)
    colorscale <- colorscale$hex.1.n.
  }
  gene_index <- which(gene_names == gene_name)
  expression <- data[gene_index,]
  plot_data <- as.data.frame(cbind(exp = expression, cluster = labels))
  plot_data$cluster = with(plot_data, reorder(cluster, exp, mean))

  #reorder the color scale so that ascending-mean factor order is retained
  colorscale <- colorscale[as.numeric(levels(plot_data$cluster))]
  print(ggplot(plot_data, aes_string(x = 'cluster', y = 'exp', fill = 'cluster')) +
          geom_violin(alpha = 0.6) +
          theme(legend.position = 'none') +
          labs(title = paste0("Expression of ", gene_name), x = "Cluster", y = "Relative Target Expression") +
          geom_jitter(shape=16, position=position_jitter(0.2)) +
          scale_fill_manual(values = colorscale) + 
          theme_minimal())
}

#' Circos plot of signaling score
#' 
#' Produce a plot of a cell signaling network according to the cell-normalized expression score.
#'
#' @param P signaling probabilities cells x cells
#' @param rgb_gap real percent of the inter-cluster spectrum gap to encompass cluster highlighting
#' @param cluster_labels labels of cells ordered from 1 to n
#' @param lig_cells_per_cluster int 0 to keep all
#' @param lig_cells_total int 0 to keep remaining cells per cluster
#' @param lig_edge_per_cell int 0 to keep all
#' @param lig_edge_per_cluster int 0 to keep remaining edge per cell
#' @param lig_edge_total int 0 to keep keep remaining edge per cluster
#' @param rec_cells_per_cluster int 0 to keep all
#' @param rec_cells_total int 0 to keep remaining cells per cluster
#' @param zero_threshold real value for zero cutoff (0 to keep all edges)
#' @param cD_reduce ratio of the circlos plot grid element (cell or cluster) to the whole circle is less than this value, the item is omitted from the plot.  Default is zero, don't change unless you are inspecting an unhighlighted chord diagram, since it breaks highlighting if one or more pairs above the zero_threshold happen to have a low grid/circle ratio.
#' @param highlight_clusters whether to label the plot with clusters.  Default is true.
#' 
#' 
#' @import dplyr
#' @import circlize
#' @importFrom reshape2 melt
#'
#' @export
#'
SigPlot <- function(P,
                     rgb_gap = 0.2,
                     cluster_labels,
                     lig_cells_per_cluster = 10,
                     lig_cells_total = 0,
                     lig_edge_per_cell = 0,
                     lig_edge_per_cluster = 0,
                     lig_edge_total = 0,
                     rec_cells_per_cluster = 10,
                     rec_cells_total = 0,
                     zero_threshold = 0,
                     cD_reduce = 0,
                     highlight_clusters = TRUE){
  label = lig_cluster_number = lig_cell = NULL # r cmd check pass
  circos.clear()
  # compute
  # ordering: int vector of indices corresponding to the labels
  # labels: int vector cluster labels of the cells ordered from cluster 1 to n
  sorted_cell <- sort.int(cluster_labels, index.return = TRUE)
  ordering <- sorted_cell$ix
  labels <- sorted_cell$x
  
  n_clusters <- length(unique(labels))
  n_cells <- length(labels)
  
  # find nonzero rows and columns of ordered P
  P <- P[ordering, ordering]
  nzrow <- which(rowSums(P) > zero_threshold)
  nzcol <- which(colSums(P) > zero_threshold)
  
  if(length(nzrow) < n_clusters | length(nzcol) < n_clusters){
    print("number of nonzero cells is less than the number of clusters, reduce zero threshold")
    return()
  }
  
  # prune the ordering and labels for rows and cols
  nzordering_lig <- ordering[nzrow]
  nzordering_rec <- ordering[nzcol]
  nzlabel_lig <- labels[nzrow]
  nzlabel_rec <- labels[nzcol]
  
  # prune P
  P_ordered <- P[nzrow, nzcol]
  rownames(P_ordered) <- nzordering_lig
  colnames(P_ordered) <- nzordering_rec
  
  
  # for each cluster, get the location of the first and the last cell in the permutation of labels
  lig_counts <- as.matrix(table(nzlabel_lig))
  accumulator <- matrix(1, n_clusters, n_clusters)*lower.tri(matrix(1, n_clusters, n_clusters), diag = TRUE)
  lig_last_cells <- accumulator %*% lig_counts
  lig_first_cells <- lig_last_cells - lig_counts + 1
  
  rec_counts <- as.matrix(table(nzlabel_rec))
  accumulator <- matrix(1, n_clusters, n_clusters)*lower.tri(matrix(1, n_clusters, n_clusters), diag = TRUE)
  rec_last_cells <- accumulator %*% rec_counts
  rec_first_cells <- rec_last_cells - rec_counts + 1
  
  ### apply cell and edge thresholds ###
  
  # apply top cell per cluster lig
  if(lig_cells_per_cluster > 0){
    sums <- rowSums(P_ordered)
    tab <- data.frame(cell = nzordering_lig,
                      label = nzlabel_lig,
                      sum = sums)
    filtered <- tab %>%
      group_by(label) %>%
      top_n(n = lig_cells_per_cluster,
            wt = sum)
    P_ordered <- P_ordered[as.character(filtered$cell),]  
  }
  
  # apply top cell per cluster rec
  if(rec_cells_per_cluster > 0){
    sums <- colSums(P_ordered)
    tab <- data.frame(cell = nzordering_rec,
                      label = nzlabel_rec,
                      sum = sums)
    filtered <- tab %>%
      group_by(label) %>%
      top_n(n = rec_cells_per_cluster,
            wt = sum)
    P_ordered <- P_ordered[,as.character(filtered$cell)]  
  }
  
  
  
  # flatten the adjacency matrix:
  # col1 is the ligand cell
  # col2 is the receptor cell
  # col3 is the weight of the link
  P_table <- melt(P_ordered)
  colnames(P_table) <- c("lig_cell", "rec_cell", "link_weight")
  ind_lig <- match(P_table$lig_cell, ordering)
  lig_clust_num <- labels[ind_lig]
  P_table$lig_cluster_number <- lig_clust_num
  
  ind_rec <- match(P_table$rec_cell, ordering)
  rec_clust_num <- labels[ind_rec]
  P_table$rec_cluster_number <- rec_clust_num
  
  P_table <- arrange(P_table, lig_cluster_number, lig_cell)
  
  # get the ordering of the receptors
  lig_cell_order <- unique(P_table$lig_cell[which(P_table$link_weight > zero_threshold)])
  P_tab <- arrange(P_table[which(P_table$link_weight > zero_threshold),], rec_cluster_number, rec_cell)
  rec_cell_unique <- unique(P_tab$rec_cell)
  rec_cell_unique <- paste(rec_cell_unique, "R", sep = "_")
  chord_plot_sector_order <- c(lig_cell_order, rec_cell_unique)
  
  # add a prefix to make the rec cells unique for sector highlighting
  P_table$rec_cell <- paste(P_table$rec_cell, "R", sep = "_")
  
  # get rgb codes and 0-360 hue map for the clusters
  cluster_colors <- ColorHue(n_clusters)
  alt_cluster_colors <- ColorHue(n_clusters,
                                 luminance = 100,
                                 chroma = 90)
  
  # get the rgb codes for the sectors (cells), based on 20% of the spectrum starting from the cluster hue
  gap <- rgb_gap*(cluster_colors[2,1] - cluster_colors[1,1])
  
  # get the chordcolors for the ligands
  cols <- ChordColors(P_table, cluster_colors, gap)
  
  # plot the chords
  circos.clear()
  chordDiagram(P_table[which(P_table$link_weight > zero_threshold),1:3],
               order = chord_plot_sector_order,
               directional = TRUE, 
               direction.type = c("diffHeight", "arrows"), 
               link.arr.type = "big.arrow", 
               annotationTrack = "grid", 
               grid.col = cols, 
               preAllocateTracks = list(list(track.height = 0.05), list(track.height = 0.05)),
               reduce = cD_reduce)
  # apply highlighting to the ligand signaling cells
  
  # Circlize only plots the P_table connections that are non-zero
  # In case zero_threshold is <= 0, find which clusters are being plotted
  
  # find highlightable pairs
  if(length(which(P_table$link_weight <= zero_threshold)) > 0){
    nz_lig_clust <- unique(P_table$lig_cluster_number[-(which(P_table$link_weight <= zero_threshold))])
    nz_rec_clust <- unique(P_table$rec_cluster_number[-(which(P_table$link_weight <= zero_threshold))])
  } else {
    nz_lig_clust <- unique(P_table$lig_cluster_number)
    nz_rec_clust <- unique(P_table$rec_cluster_number)
  }
  
  if(highlight_clusters){
    for(i in nz_lig_clust){
      lig_cells <- unique(P_table$lig_cell[which(P_table$lig_cluster_number == i)])
      highlight_col <- cluster_colors$hex.1.n.[i]
      cluster_name <- paste0("C", i)
      highlight.sector(sector.index = lig_cells,
                       col = highlight_col, 
                       text = cluster_name, 
                       text.vjust = -1, 
                       niceFacing = TRUE, 
                       track.index = 2)
    }
    for(i in nz_rec_clust){
      rec_cells <- unique(P_table$rec_cell[which(P_table$rec_cluster_number == i)])
      highlight_col <- cluster_colors$hex.1.n.[i]
      cluster_name <- paste0("C", i)
      highlight.sector(sector.index = rec_cells,
                       col = highlight_col, 
                       text = cluster_name, 
                       text.vjust = -1, 
                       niceFacing = TRUE, 
                       track.index = 2)
    }  
  }
}

#' Get a vector of n equally spaced rgb colors
#' 
#' Get a vector of n equally spaced rgb colors.
#'
#' @param n integer number of hex codes to return
#' @param starthue real hue argument for the grDevices::hcl() generates value 1
#' @param endhue real hue argument
#' @param luminance the luminance of the hcl value
#' @param chroma the chroma of the hcl value
#' 
#' @importFrom grDevices hcl
#' 
#' @export
#'
ColorHue <- function(n, starthue = 15, endhue = 360,
                     luminance = 65, chroma = 100) {
  hues <- seq(starthue, endhue, length = n + 1)
  hex <- hcl(h = hues, l = luminance, c = chroma)[1:n]
  hue_color_map <- data.frame(hues[1:n], hex[1:n])
  hue_color_map[,2] <- as.character(hue_color_map[,2])
  return(hue_color_map)
}

#' Get chord colors for cells
#' 
#' Get chord colors for cells.
#'
#' @param edge_table a table with lig, rec, score, lig_cluster, rec_cluster
#' @param cluster_cols colors for the clusters
#' @param gap the hue gap
#'
ChordColors <- function(edge_table, cluster_cols, gap) {
  chords <- c()
  for(i in 1:length(unique(edge_table$lig_cluster_number))){
    cell_labels <- edge_table[which(edge_table$lig_cluster_number == i),]$lig_cell
    starthue <- cluster_cols$hues.1.n.[i]
    cols <- ColorHue(n = length(cell_labels), 
                     starthue = starthue,
                     endhue = starthue + gap)
    cols <- cols$hex.1.n.
    names(cols) <- cell_labels
    chords <- c(chords, cols)
  }
  cols_rec <- rep("grey", length(edge_table$rec_cell))
  names(cols_rec) <- edge_table$rec_cell
  chords <- c(chords, cols_rec)
}

