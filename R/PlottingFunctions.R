#' Plot the pseudotime ordering of clusters
#'
#' @param network an igraph network with weighted edges corresponding to the pseudotime distance
#' @param d_thickness controls edge width.  c(1, m) => m/distance from L2R2, c(0,m) => m.
#' @param node_color an optional vector of colors for node labeling.  If null then the nodes are all colored black
#' @param alpha_color the gradient of the nodes
#'
#' @return nothing
#'
#' @import intergraph
#' @import network
#' @import sna
#' @importFrom igraph E
#' @import ggnetwork
#'
#' @export
#'
PlotLineage <- function(network, d_thickness = c(1,2), node_color = NULL, alpha_color = 1){
  if(is.null(node_color)){
    node_color = replicate(igraph::vcount(network), "black")
  }
  set.seed(1)
  if(d_thickness[1]){
    print(
      ggplot(network,
             arrow.gap = 0.05,
             aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(arrow = arrow(length = unit(1, "lines"),
                                 type = 'open'),
                   color = "black",
                   aes(size=d_thickness[2]/weight),
                   show.legend = FALSE) +
        geom_nodes(color = "gold", size = 8) +
        geom_nodetext(aes(color = "black",
                          label = vertex.names,
                          size = 2),
                      fontface = "bold",
                      show.legend = FALSE) +
        theme_blank())
  }
}

#' Produce a plot of matlab DTree data and return the object
#' If an output dir and filename are provided, a plot will
#' be saved, otherwise the function will just return the graph
#'
#' @param edge_table a numeric matrix whose rows are directed edges of a tree:
#'     col 1 is v1, col2 is v2, col3 is weight
#' @param predecessors a vector of tree predecessors such that pred[i] = the predecessor of i
#' @param outputdir the output directory, relative to getwd()
#' @param outputfile the output file
#'
#' @return an igraph representation of the tree
#'
#' @export
#'
PlotMatlabDtree <- function(edge_table, predecessors, outputdir = NULL, outputfile = NULL){
  directed_edge_table <- RSoptSC::ProcessMatlabDTree(edge_table, predecessors)
  directed_graph <- RSoptSC::GetDGFromTable(directed_edge_table)
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

#' Produce a scatter plot of the cells on selected 2-dim ebedding colored by feature
#'
#' If an output dir and filename are provided, a plot will
#' be saved
#'
#' @param flat_embedding a low dim embedding of cells
#' @param feature a scalar representation of feature
#' @param outputdir the output directory, relative to getwd()
#' @param outputfile the output file
#' @param title the title of the plot
#' @param subtitle the subtitle of the plot
#' @param featurename the name of the feature to title the legend
#'
#' @return a ggplot2 object
#'
#' @export
#'
FeatureScatterPlot <- function(flat_embedding,
                               feature,
                               outputdir = NULL,
                               outputfile = NULL,
                               title,
                               subtitle,
                               featurename){
  p <- ggplot2::ggplot(as.data.frame(flat_embedding), ggplot2::aes(x=V1, y=V2, color=feature))

  if(!is.null(outputdir) && !is.null(outputfile)){
    # check if the dir exists and if not then create it
    if (!file.exists(outputdir)){
      dir.create(file.path(getwd(), outputdir))
    }
    file_path <- file.path(getwd(), outputdir, outputfile)
    pdf(file_path)
    p +
    ggplot2::geom_point() +
    labs(title = title, subtitle = subtitle, color = featurename)
    dev.off()
  }
  p +
  ggplot2::geom_point() +
  labs(title = title, subtitle = subtitle, color = featurename)
}

#' Get the marker genes for each cluster
#'
#' @param counts_data a matrix of expression values for each cell (rows) and gene (columns)
#' @param cluster_labels a vector of cluster labels
#' @param markerTabele a matrix of expression values for each cell (rows) and gene (columns)
#' @param range the interval over which centered, normalized expression is displayed. expression values below range[1] are increased to range[1].  values above range[2] are decreased to range[2].
#' @param n_markers number of marker genes per cluster to retrieve (0 = all)
#'
#' @return a table of marker genes
#'
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#'
#' @export
#'
MarkerHeatmap <- function(counts_data,
                          cluster_labels,
                          markerTable,
                          range = c(-3,3),
                          n_markers = 0){
  # work with the data in tibble form
  tibbleData <- tibble::as.tibble(markerTable)
  # sort the genes by cluster
  if(n_markers > 0){
    byCluster <- tibbleData[order(tibbleData$clusterId),]
    clusterList <- split(byCluster, byCluster$clusterId)
    topN <- lapply(clusterList, function(x){
                  head(x, n_markers)})
    gene_order <- dplyr::pull(do.call(rbind, topN), geneID)
  } else {
    byCluster <- tibbleData[order(tibbleData$clusterId),]
    clusterList <- split(byCluster, byCluster$clusterId)
    topN <- lapply(clusterList, function(x){
      head(x, nrow(markerTable))})
    gene_order <- dplyr::pull(do.call(rbind, topN), geneID)
  }


  # flatten the data (matrix form -> table form) so that ggplot can work
  cell_order <- order(cluster_labels)
  plot_data <- counts_data[gene_order, cell_order]
  plot_data <- ScaleCenterData(plot_data)
  plot_data_table <- melt(plot_data)
  colnames(plot_data_table) <- c("geneID", "cellId", "expression")

  # reduce the dynamic range of the plot, focusing on variation close to the mean of the centered data
  plot_data_table[which(plot_data_table[,3] < range[1]),3] <- range[1]
  plot_data_table[which(plot_data_table[,3] > range[2]),3] <- range[2]

  print(ggplot(plot_data_table, aes(cellId, geneID )) +
    geom_tile(aes(fill = expression), color = "white") +
    scale_fill_gradient(low = "red", high = "green") +
    ylab("List of Genes ") +
    xlab("List of Cells") +
    theme(legend.title = element_text(size = 10),
    legend.text = element_text(size = 12),
    plot.title = element_text(size=16),
    axis.title=element_text(size=14,face="bold"),
    axis.text.x = element_text(angle = 90, hjust = 1)) +labs(fill = "expression"))
}

#' Produce a violin plot of gene expression
#'
#' @param data a matrix with genes in rows and cells in columns
#' @param gene_names a vector of names corresponding to the rows in data
#' @param labels a vector of cluster assignments for the cells
#' @param gene_name the name of the gene to plot
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
                                 gene_name){
  gene_index <- which(gene_names == gene_name)
  expression <- data[gene_index,]
  plot_data <- as.data.frame(cbind(exp = expression, cluster = labels))
  plot_data$cluster = with(plot_data, reorder(cluster, exp, mean))
  print(ggplot(plot_data, aes(x = cluster, y = exp, fill = cluster)) +
    geom_violin(alpha = 0.6) +
    scale_fill_discrete() +
    theme(legend.position = 'none') +
    labs(title = paste0("Expression of ", gene_name), x = "Cluster", y = "Relative Target Expression") +
    geom_jitter(shape=16, position=position_jitter(0.2)))
}
