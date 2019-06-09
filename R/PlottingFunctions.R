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
#' @param pointsize the size of geom point
#' @param hide_legend hide the legend and return a list, containing the legendless plot and the legend
#' @param legsize the size of the legend
#' @param legt size of legend title
#' @param legtxt size of legend text
#' @param grad_lim the limits of the color gradient for continuous features, this should usually be in log-expression space, overrides autoscale_dropout
#' @param autoscale_dropout whether to set plots where the feature space is 0 to #102033, default min for scale_color_gradient auto scaling, overrides dropout_NA
#' @param dropout_NA whether to set plots where the feature space is 0 to NA, note: the legend element in a hide_legend output will be NULL
#'
#' @return a ggplot2 object
#' 
#' @import ggplot2
#' @importFrom cowplot get_legend
#'
#' @export
#'
FeatureScatterPlot <- function(flat_embedding,
                               feature,
                               title,
                               subtitle,
                               featurename,
                               colorscale = NULL,
                               pointsize = 1,
                               hide_legend = FALSE,
                               legsize = 5,
                               legt = 5,
                               legtxt = 5,
                               grad_lim = NULL,
                               autoscale_dropout = FALSE,
                               dropout_NA = TRUE){
  n_features <- length(unique(feature))
  
  p <- ggplot(flat_embedding,
              aes_(x = as.name(colnames(flat_embedding)[1]),
                   y = as.name(colnames(flat_embedding)[2]),
                   color = ~`feature`))
  
  
  if(is.null(colorscale)){
    p <- p +
      geom_point(size = pointsize) +
      labs(title = title, subtitle = subtitle, color = featurename) +
      theme_minimal() + 
      theme(legend.title = element_text(size=legt), legend.text = element_text(size = legtxt)) 
    if(!is.null(grad_lim)){
      p <- p + scale_colour_gradient(limits = grad_lim)
    } else if (autoscale_dropout) {
      if(length(unique(feature)) == 1){
        p <- p + scale_colour_gradient(low = "#102033", high = "#102033")
      }
    } else if (dropout_NA){
      if(length(unique(feature)) == 1){
        feature[which(feature == 0)] <- NA
        p <- p + scale_color_gradient(na.value = "grey50")
      }
    }
  }else{
    p <- p +
      geom_point(size = pointsize) +
      labs(title = title, subtitle = subtitle, color = featurename) +
      theme_minimal() +
      scale_color_manual(values = colorscale) +
      guides(colour = guide_legend(override.aes = list(size=legsize))) +
      theme(legend.title = element_text(size=legt), legend.text = element_text(size = legtxt))
  }
  if(hide_legend){
    if(length(unique(feature)) == 1){
      # NAs have replaced the dropout 0's
      p <- p + theme_minimal()
      p <- p + theme(legend.position = "none")
      list(legend = NULL, plot = p)
    } else {
      leg <- get_legend(p)
      p <- p + theme_minimal()
      p <- p + theme(legend.position = "none")
      list(legend = leg, plot = p)
    }
  }else{
    p
  }
}

#' Factor Grid heatmap of the top n markers
#' 
#' Plot the heatmap of the top n inferred markers for each cluster.
#'
#' @param data expression data with genes x cells
#' @param cluster_labels a vector of cluster labels corresponding to the columns in the data matrix
#' @param markers a table of markers
#' @param n_features the top n features to plot
#' @param y_lsize y axis label size
#' @param y_tsize y axis title size
#' @param x_lsize x axis label size
#' @param x_tsize x axis title size
#' @param use_z use z-score per gene instead of raw expression
#' @param spacing the number of lines between grids
#' 
#' @return a ggplot2 plot object
#'
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom Matrix as.matrix
#'
#' @export
#'

PlotTopN_Grid <- function(data,
                          cluster_labels, 
                          markers,
                          n_features,
                          y_lsize = ((length(unique(cluster_labels))*n_features)/120)*5,
                          y_tsize = ((length(unique(cluster_labels))*n_features)/120)*15,
                          x_lsize = ((length(unique(cluster_labels))*n_features)/120)*14,
                          x_tsize = ((length(unique(cluster_labels))*n_features)/120)*15,
                          use_z = TRUE,
                          spacing = 0.1){
  clusterId = geneScore = fac_barcode = fac_symbol = NULL
  names(cluster_labels) <- colnames(data)
  ind <- order(cluster_labels)
  barcode_order <- colnames(data)[ind]
  
  filtered_markers <- markers %>% group_by(clusterId) %>% top_n(n_features, geneScore) 
  if(use_z){
    data <- apply(data, 1, function(x){
      m <- mean(x)
      sd <- sd(x)
      z <- (x-m)/sd
    })
    data <- t(data)
    expr_label <- "Z-score"
  }else{
    expr_label <- "Log Expression"
  }
  
  melted <- melt(as.matrix(data))
  colnames(melted) <- c("geneSymbol", "barcode", "expression")
  
  melted$geneSymbol <- factor(melted$geneSymbol, levels=markers$geneSymbol)
  full_melted <- inner_join(melted, filtered_markers)
  full_melted <- arrange(full_melted, clusterId, desc(geneScore))
  full_melted$fac_barcode <- factor(full_melted$barcode, levels=barcode_order)
  full_melted$fac_symbol <- factor(full_melted$geneSymbol, levels=markers$geneSymbol)
  
  # get x axis cluster ticks
  counts <- unname(table(cluster_labels))
  names <- c()
  pos <- c()
  current <- 0
  for(i in 1:length(counts)){
    names[i] <- paste0(i)
    current <- current + counts[i]
    pos[i] <- floor(current - counts[i]/2)
  }
  p <- ggplot(full_melted, aes(x = fac_barcode, y = fac_symbol, fill = expression, color=expression)) +
    geom_tile() +
    scale_fill_gradient2(low = 'purple',
                         mid = 'black',
                         high = 'yellow',
                         midpoint = 0,
                         name = eval(expr_label)) +
    scale_color_gradient2(low = 'purple',
                          mid = 'black',
                          high = 'yellow',
                          midpoint = 0,
                          name = eval(expr_label)) +
    facet_grid(clusterId~., scales = "free", space = "free") + #facet by group
    ylab("Gene Symbol") +
    theme(strip.text.y = element_text(angle = 0),
          panel.border = element_rect(colour = "black", fill = NA), #add black border
          panel.spacing = unit(spacing, "lines"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size = y_tsize),
          axis.text.y = element_text(size = y_lsize),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  plot(p)
  return(p)
}

#' Heatmap of the top n markers
#' 
#' Plot the heatmap of the top n inferred markers for each cluster.
#'
#' @param data expression data with genes x cells
#' @param cluster_labels a vector of cluster labels corresponding to the columns in the data matrix
#' @param markers a table of markers
#' @param n_features the top n features to plot
#' @param y_lsize y axis label size
#' @param y_tsize y axis title size
#' @param x_lsize x axis label size
#' @param x_tsize x axis title size
#' @param use_z use z-score per gene instead of raw expression
#' 
#' @return a ggplot2 plot object
#'
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom Matrix as.matrix
#'
#' @export
#'
PlotTopN <- function(data,
                     cluster_labels, 
                     markers,
                     n_features,
                     y_lsize = ((length(unique(cluster_labels))*n_features)/120)*5,
                     y_tsize = ((length(unique(cluster_labels))*n_features)/120)*15,
                     x_lsize = ((length(unique(cluster_labels))*n_features)/120)*14,
                     x_tsize = ((length(unique(cluster_labels))*n_features)/120)*15,
                     use_z = TRUE){
  clusterId = geneScore = fac_barcode = fac_symbol = NULL
  names(cluster_labels) <- colnames(data)
  ind <- order(cluster_labels)
  barcode_order <- colnames(data)[ind]
  # # markers <- data.frame(symbol = gene_names[markers[,'geneID']], 
  # #                       cluster = markers[,'clusterId'], 
  # #                       score = markers[,'geneScore'])
  # markers <- markers %>% group_by(clusterId) 
  # markers$symbol <- factor(markers$geneSymbol, levels=markers$symbol)
  filtered_markers <- markers %>% group_by(clusterId) %>% top_n(n_features, geneScore) 
  if(use_z){
    data <- apply(data, 1, function(x){
      m <- mean(x)
      sd <- sd(x)
      z <- (x-m)/sd
    })
    data <- t(data)
    expr_label <- "Z-score"
  }else{
    expr_label <- "Log Expression"
  }
  
  melted <- melt(as.matrix(data))
  colnames(melted) <- c("geneSymbol", "barcode", "expression")
  
  melted$geneSymbol <- factor(melted$geneSymbol, levels=markers$geneSymbol)
  full_melted <- inner_join(melted, filtered_markers)
  full_melted <- arrange(full_melted, clusterId, desc(geneScore))
  full_melted$fac_barcode <- factor(full_melted$barcode, levels=barcode_order)
  full_melted$fac_symbol <- factor(full_melted$geneSymbol, levels=markers$geneSymbol)
  
  # get x axis cluster ticks
  counts <- unname(table(cluster_labels))
  names <- c()
  pos <- c()
  current <- 0
  for(i in 1:length(counts)){
    names[i] <- paste0(i)
    current <- current + counts[i]
    pos[i] <- floor(current - counts[i]/2)
  }
  p <- ggplot(full_melted, aes(x = fac_barcode, y = fac_symbol, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(low = 'purple',
                         mid = 'black',
                         high = 'yellow',
                         midpoint = 0,
                         name = eval(expr_label)) +
    xlab("Cluster") +
    ylab("Gene Symbol") +
    scale_x_discrete(breaks = barcode_order[pos], labels = names) +
    theme_bw() +
    theme(axis.title.x = element_text(size = x_tsize, angle = 0, vjust = 0.5),
          axis.title.y = element_text(size = y_tsize),
          axis.text.x = element_text(size = x_lsize),
          axis.text.y = element_text(size = y_lsize))
  plot(p)
  return(p)
  
}

#' Heatmap of specific genes
#' 
#' Given a set of markers, plot their expression across clusters on a heatmap.
#'
#' @param data expression data with genes x cells
#' @param cluster_labels a vector of cluster labels corresponding to the columns in the data matrix
#' @param y_lsize y axis label size
#' @param x_lsize x axis label size
#' @param y_tsize y axis title size
#' @param x_tsize x axis title size
#' @param lti_size size of legend title
#' @param lt_size size of legend text
#' @param markers a vector of markers
#' @param use_z use z-score per gene instead of raw expression
#' 
#' @return nothing
#'
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @export
#'
PlotClusterExpression <- function(data,
                                   cluster_labels,
                                   y_lsize = (length(unique(cluster_labels))/8)*7,
                                   y_tsize = (length(unique(cluster_labels))/8)*7,
                                   x_lsize = (length(unique(cluster_labels))/8)*7,
                                   x_tsize = (length(unique(cluster_labels))/8)*7,
                                   lti_size = (length(unique(cluster_labels))/8)*7,
                                   lt_size = (length(unique(cluster_labels))/8)*7,
                                   markers,
                                   use_z = TRUE){
  cluster = symbol = NULL
  n_clusters <- length(unique(cluster_labels))
  sorted_cell <- sort.int(cluster_labels, index.return = TRUE)
  gene_names <- rownames(data)
  ind <- match(markers, gene_names)
  naind <- which(is.na(ind))
  if(length(naind) > 0){
    markers <- markers[-naind]
    print("Omit missing genes:")
    print(markers[naind])
  }
  filtered <- ClusterAvg(M = data,
                         gene_names = gene_names,
                         cell_order = sorted_cell$ix,
                         cell_labels = cluster_labels,
                         gene_list = markers)
  
  if(use_z){
    filtered <- apply(filtered, 1, function(x){
      m <- mean(x)
      sd <- sd(x)
      z <- (x-m)/sd
    })
    filtered <- t(filtered)
    expr_label <- "Z-score"
  }else{
    expr_label <- "Log Expression"
  }
  melted <- melt(filtered)
  colnames(melted) <- c("cluster", "symbol", "expression")
  p <- ggplot(melted, aes(x = as.factor(cluster), y = symbol, fill = expression)) +
    geom_tile() +
    scale_fill_gradient2(low = 'purple',
                         mid = 'black',
                         high = 'yellow',
                         midpoint = 0,
                         name = eval(expr_label)) +
    xlab("Cluster") +
    ylab("Gene Symbol") +
    theme_bw() + 
    theme(axis.title.x = element_text(size = x_tsize, angle = 0, vjust = 0.5),
          axis.title.y = element_text(size = y_tsize),
          axis.text.x = element_text(size = x_lsize),
          axis.text.y = element_text(size = y_lsize),
          legend.title = element_text(size = lti_size),
          legend.text = element_text(size = lt_size),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  plot(p)
  return(p)
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
#' @param jitsize the jitter size parameter
#'
#' @return a violin plot
#'
#' @import ggplot2
#'
#' @export
#'
ViolinPlotExpression <- function(data,
                                 gene_names,
                                 labels,
                                 gene_name,
                                 colorscale = NULL,
                                 jitsize = 0.2){
  n_features <- length(unique(labels))
  if(is.null(colorscale)){
    colorscale <- ColorHue(n = n_features)
    colorscale <- colorscale$hex.1.n.
  }
  gene_index <- which(tolower(gene_names) == tolower(gene_name))
  expression <- data[gene_index,]
  plot_data <- as.data.frame(cbind(exp = expression, cluster = labels))
  plot_data$cluster <- as.factor(plot_data$cluster)
  colors <- colorscale[sort(unique(labels))];
  p <- ggplot(plot_data, aes_string(x = 'cluster', y = 'exp', fill = 'cluster')) +
    geom_violin(alpha = 0.6) +
    theme(legend.position = 'none') +
    labs(title = paste0("Expression of ", gene_name), x = "Cluster", y = "Relative Target Expression") +
    geom_jitter(shape=16, position=position_jitter(jitsize)) +
    scale_fill_manual(values = colors) + 
    theme_minimal()
}

#' Circos plot of signaling score
#' 
#' Produce a plot of a cell signaling network according to the cell-normalized expression score.
#'
#' @param P signaling probabilities cells x cells
#' @param rgb_gap real percent of the inter-cluster spectrum gap to encompass cluster highlighting
#' @param cluster_labels labels of cells ordered from 1 to n
#' @param lig_cells_per_cluster int 0 to keep all
#' @param rec_cells_per_cluster int 0 to keep all
#' @param zero_threshold real value for zero cutoff (0 to keep all edges)
#' @param cD_reduce ratio of the circlos plot grid element (cell or cluster) to the whole circle is less than this value, the item is omitted from the plot.  Default is zero, don't change unless you are inspecting an unhighlighted chord diagram, since it breaks highlighting if one or more pairs above the zero_threshold happen to have a low grid/circle ratio.
#' @param highlight_clusters whether to label the plot with clusters.  Default is true.
#' @param title_text the title of the plot
#' @param n_clusters if n is higher than the number of unique labels (colorscale can be used for meta plotting)
#' @param n_cluster_cols for plots with different subpopulations, use this to reserve dropout clusters in individual plots.  WARNING: this is strongly recommended if there are dropouts, ie unique(labels) != 1:length(labels)
#' @param receptor_chord_color color the chords according to the receptor-bearing cell they point to.  default is false
#' 
#' 
#' @import dplyr
#' @import circlize
#' @importFrom reshape2 melt
#' @importFrom graphics legend title
#'
#' @export
#'
SigPlot <- function(P,
                    rgb_gap = 0.2,
                    cluster_labels,
                    lig_cells_per_cluster = 10,
                    rec_cells_per_cluster = 10,
                    zero_threshold = 0,
                    cD_reduce = 0,
                    highlight_clusters = TRUE,
                    title_text = NULL,
                    n_clusters = NULL,
                    n_cluster_cols = NULL,
                    receptor_chord_color = FALSE){
  if(max(P) <= zero_threshold){
    print('no signaling in P')
    return()
  }
  label = lig_cluster_number = lig_cell = rec_cell = 
    rec_cluster_number = legend = title = link_weight = arbitrary_order = NULL # r cmd check pass
  circos.clear()
  # compute
  # ordering: int vector of indices corresponding to the labels
  # labels: int vector cluster labels of the cells ordered from cluster 1 to n
  sorted_cell <- sort.int(cluster_labels, index.return = TRUE)
  ordering <- sorted_cell$ix
  labels <- sorted_cell$x
  
  if(is.null(n_clusters)){
    n_clusters <- length(unique(labels))
  }
  if(is.null(n_cluster_cols)){
    n_cluster_cols <- max(labels)
  }
  
  n_cells <- length(labels)
  
  # find nonzero rows and columns of ordered P
  P <- P[drop=FALSE,ordering, ordering]
  nzrow <- which(rowSums(P) > zero_threshold)
  nzcol <- which(colSums(P) > zero_threshold)
  
  
  # prune the ordering and labels for rows and cols
  nzordering_lig <- ordering[nzrow]
  nzordering_rec <- ordering[nzcol]
  nzlabel_lig <- labels[nzrow]
  nzlabel_rec <- labels[nzcol]
  
  # prune P
  P_ordered <- P[drop=FALSE,nzrow, nzcol]
  rownames(P_ordered) <- nzordering_lig
  colnames(P_ordered) <- nzordering_rec
  
  # for each cluster, get the location of the first and the last cell in the permutation of labels
  lig_counts <- matrix(0, nrow = n_clusters, ncol = 1)
  rec_counts <- matrix(0, nrow = n_clusters, ncol = 1)
  rownames(lig_counts) <- rownames(rec_counts) <- as.character(c(1:n_clusters))
  for(c in 1:n_clusters){
    cc <- unique(labels)[c]
    lig_counts[c,1] <- length(which(nzlabel_lig == cc))
    rec_counts[c,1] <- length(which(nzlabel_rec == cc))
  }
  
  # for each cluster, get the location of the first and the last cell in the permutation of labels
  accumulator <- matrix(1, n_clusters, n_clusters)*lower.tri(matrix(1, n_clusters, n_clusters), diag = TRUE)
  lig_last_cells <- accumulator %*% lig_counts
  lig_first_cells <- lig_last_cells - lig_counts + 1
  
  accumulator <- matrix(1, n_clusters, n_clusters)*lower.tri(matrix(1, n_clusters, n_clusters), diag = TRUE)
  rec_last_cells <- accumulator %*% rec_counts
  rec_first_cells <- rec_last_cells - rec_counts + 1
  
  ### apply cell and edge thresholds ###
  # apply top cell per cluster lig
  if(lig_cells_per_cluster > 0){
    sums <- rowSums(P_ordered)
    tab <- data.frame(cell = nzordering_lig,
                      label = nzlabel_lig,
                      sum = sums,
                      arbitrary_order = 1:length(sums))
    grouped <- tab %>% group_by(label)
    filtered <- top_n(grouped, n = lig_cells_per_cluster,wt = arbitrary_order)
    P_ordered <- P_ordered[drop=FALSE,as.character(filtered$cell),]  
  }
  
  # apply top cell per cluster rec
  if(rec_cells_per_cluster > 0){
    sums <- colSums(P_ordered)
    tab <- data.frame(cell = nzordering_rec,
                      label = nzlabel_rec,
                      sum = sums,
                      arbitrary_order = 1:length(sums))
    grouped <- tab %>%
      group_by(label)
    filtered <- top_n(grouped, n = rec_cells_per_cluster, wt = arbitrary_order)
    P_ordered <- P_ordered[drop=FALSE,,as.character(filtered$cell)]  
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
  cluster_colors <- ColorHue(n_cluster_cols)
  alt_cluster_colors <- ColorHue(n_cluster_cols,
                                 luminance = 100,
                                 chroma = 90)
  # rownames(cluster_colors) <- sort(unique(labels))
  # rownames(alt_cluster_colors) <- sort(unique(labels))
  
  # cluster_colors <- cluster_colors[sort(unique(cluster_labels)),]
  # alt_cluster_colors <- alt_cluster_colors[sort(unique(cluster_labels)),]
  
  # get the rgb codes for the sectors (cells), based on 20% of the spectrum starting from the cluster hue
  gap <- rgb_gap*(cluster_colors[2,1] - cluster_colors[1,1])
  
  if(!receptor_chord_color){
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
        ii <- which(rownames(cluster_colors) == i)
        lig_cells <- unique(P_table$lig_cell[which(P_table$lig_cluster_number == i)])
        highlight_col <- cluster_colors$hex.1.n.[ii]
        cluster_name <- paste0("C", i)
        highlight.sector(sector.index = lig_cells,
                         col = highlight_col, 
                         #text = cluster_name, # these may be cramped
                         text.vjust = -1, 
                         niceFacing = TRUE, 
                         track.index = 2)
      }
      for(i in nz_rec_clust){
        ii <- which(rownames(cluster_colors) == i)
        rec_cells <- unique(P_table$rec_cell[which(P_table$rec_cluster_number == i)])
        highlight_col <- cluster_colors$hex.1.n.[ii]
        cluster_name <- paste0("C", i)
        highlight.sector(sector.index = rec_cells[which(rec_cells %in% get.all.sector.index())],
                         col = highlight_col, 
                         #text = cluster_name, # these may be cramped
                         text.vjust = -1, 
                         niceFacing = TRUE, 
                         track.index = 2)
      }
      subset_cols <- cluster_colors$hex.1.n.[sort(unique(labels))]
      legend("topleft", legend = paste0("C", sort(unique(cluster_labels))), pch=16, pt.cex=1.5, cex=1, bty='n',
             col = subset_cols)
      title(title_text, cex.main = 1.5, line = -0.5)
    }
  } else {
    P_table <- arrange(P_table, lig_cluster_number, rec_cluster_number, desc(link_weight), rec_cell)
    P_table <- P_table[which(P_table$link_weight > zero_threshold),]
    #colnames(P_table) <- c("lig_link", "rec_cell", "link_weight", "lig_cluster_number", "rec_cluster_number")
    P_table$lig_cell <- 1:nrow(P_table) # this is not actually cells anymore, more like cell-extracted links
    chord_plot_sector_order <- c(P_table$lig_cell, rec_cell_unique)
    # get the chordcolors for the ligands
    cols <- ChordColors(P_table, cluster_colors, gap, receptor_chord_color = TRUE)
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
        ii <- which(rownames(cluster_colors) == i)
        lig_cells <- unique(P_table$lig_cell[which(P_table$lig_cluster_number == i)])
        highlight_col <- cluster_colors$hex.1.n.[ii]
        cluster_name <- paste0("C", i)
        highlight.sector(sector.index = lig_cells,
                         col = highlight_col, 
                         #text = cluster_name, # these may be cramped
                         text.vjust = -1, 
                         niceFacing = TRUE, 
                         track.index = 2)
      }
      for(i in nz_rec_clust){
        ii <- which(rownames(cluster_colors) == i)
        rec_cells <- unique(P_table$rec_cell[which(P_table$rec_cluster_number == i)])
        highlight_col <- cluster_colors$hex.1.n.[ii]
        cluster_name <- paste0("C", i)
        highlight.sector(sector.index = rec_cells[which(rec_cells %in% get.all.sector.index())],
                         col = highlight_col, 
                         #text = cluster_name, # these may be cramped
                         text.vjust = -1, 
                         niceFacing = TRUE, 
                         track.index = 2)
      }
      subset_cols <- cluster_colors$hex.1.n.[sort(unique(labels))]
      legend("topleft", legend = paste0("C", sort(unique(cluster_labels))), pch=16, pt.cex=1.5, cex=1, bty='n',
             col = subset_cols)
      title(title_text, cex.main = 1.5, line = -0.5)
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
#' @param receptor_chord_color color the chords according to the receptor-bearing cell they point to.  default is false
#'
ChordColors <- function(edge_table, cluster_cols, gap,
                        receptor_chord_color = FALSE) {
  chords <- c()
  if(!receptor_chord_color){
    for(i in unique(edge_table$lig_cluster_number)){
      cell_labels <- edge_table[which(edge_table$lig_cluster_number == i),]$lig_cell
      starthue <- cluster_cols[as.character(i),1]
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
    
  } else {
    for(i in unique(edge_table$lig_cluster_number)){
      sub_table <- edge_table[which(edge_table$lig_cluster_number == i),]
      for(ii in unique(sub_table$rec_cluster_number)){
        sub_sub_table <- sub_table[which(sub_table$rec_cluster_number == ii),]
        starthue <- cluster_cols[as.character(ii),1]
        cols <- ColorHue(n = nrow(sub_sub_table), 
                         starthue = starthue,
                         endhue = starthue + gap)
        cols <- cols$hex.1.n.
        names(cols) <- sub_sub_table$lig_cell
        chords <- c(chords, cols)
      }
    }
    cols_rec <- rep("grey", length(edge_table$rec_cell))
    names(cols_rec) <- edge_table$rec_cell
    chords <- c(chords, cols_rec)
  }
}
