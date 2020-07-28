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
#' @return as heatmap
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
  clusterId = geneScore = fac_barcode = fac_symbol = cluster_label = NULL
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
  
  labelstab <- data.frame(fac_barcode = names((cluster_labels)), cluster_label = cluster_labels)
  full_melted <- inner_join(labelstab, full_melted) # get the cell cluster labels
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
    facet_grid(rows = vars(clusterId), cols = vars(cluster_label), scales = "free", space = "free") + #facet by group
    xlab("Cell Cluster") +
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
#' @return a heatmap
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
#' @return a grid plot
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
  return(p)
}

#' Heatmap plot of cell signaling
#' 
#' Produce a plot of a cell signaling network according to the cell-normalized expression score.
#'
#' @param p signaling probabilities cells x cells
#' @param labels labels of cells ordered from 1 to n
#' @param textsize the ggplot text size
#' @param plottitle the plot title
#' @param plotsubtitle the plot sub title
#' 
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom graphics legend title
#'
#' @export
#'
SigPlot <- function(p,
                    labels,
                    textsize = 10,
                    plottitle = NULL,
                    plotsubtitle = NULL){
  withlab = inner_join(p,data.frame(labelsVar2 = labels,Var2=1:length(labels)))
  withlab = inner_join(withlab,data.frame(labelsVar1 = labels,Var1=1:length(labels)))
  withlab$Var2 = factor(withlab$Var2, levels = sort(labels,index.return=1)$ix)
  withlab$Var1 = factor(withlab$Var1, levels = sort(labels,index.return=1)$ix)
  var1sum = withlab %>% group_by(Var1) %>% summarize(Var1.sum = sum(value))
  var2sum = withlab %>% group_by(Var2) %>% summarize(Var2.sum = sum(value))
  withlab_nonzero = inner_join(withlab,var1sum[which(var1sum$Var1.sum>0),])
  withlab_nonzero = inner_join(withlab_nonzero,var2sum[which(var2sum$Var2.sum>0),])
  
  
  pp = ggplot(withlab_nonzero,mapping=aes(x = Var2,y=Var1,fill=value)) + geom_tile()
  pp + facet_grid(labelsVar1 ~ labelsVar2,scales = "free", space='free',drop=T,switch="y") + xlab("Receptor Cell") + ylab("Ligand Cell") + ggtitle(plottitle)
}

#' Heatmap plot of cluster signaling
#' 
#' Produce a plot of a cluster signaling network according to the cell-normalized expression score.
#'
#' @param p signaling probabilities clusters x clusters
#' @param labels labels of cells ordered from 1 to n
#' @param textsize the ggplot text size
#' @param plottitle the plot title
#' @param plotsubtitle the plot sub title
#' 
#' @import dplyr
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom graphics legend title
#'
#' @export
#'
SigPlot_Cluster <- function(p,
                    labels,textsize=30,
                    plottitle = NULL,
                    plotsubtitle = NULL){
  pp = ggplot(p,mapping=aes(x = cluster.Var2,y=cluster.Var1,fill=value),color="gray") + geom_tile()+ guides(color = FALSE)
  pp  + xlab("Receptor Cluster") + ylab("Ligand Cluster") + 
    labs(title = plottitle,
         subtitle = plotsubtitle) + 
    theme(legend.title = element_blank()) + 
    theme(text = element_text(size=20)) 
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