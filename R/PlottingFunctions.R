#' Plot the pseudotime ordering of clusters
#'
#' @param network an igraph network with weighted edges corresponding to the pseudotime distance
#' @param d_thickness controls edge width.  c(1, m) => m/distance from L2R2, c(0,m) => m.
#'
#' @return nothing
#'
#' @import ggplot2
#' @import intergraph
#' @import ggnetwork
#' @import sna
#' @import network
#' @importFrom igraph E
#'
#' @export
#'
PlotLineage <- function(network, d_thickness = c(1,2)){
  browser()
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
        geom_nodes(color = "black",
                   size =10) +
        geom_nodetext(aes(color = 'red',
                          label = vertex.names,
                          size = 2),
                      fontface = "bold",
                      show.legend = FALSE) +
        theme_blank())
  }else{
    print(
      ggplot(network,
             arrow.gap = 0.05,
             aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(arrow = arrow(length = unit(2, "lines"),
                                 type = 'open'),
                   color = "black",
                   aes(size=d_thickness[2]),
                   show.legend = FALSE) +
        geom_nodes(color = "black",
                   size =15) +
        geom_nodetext(aes(color = 'red',
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

#' Produce a scatter plot of the cells on selected 2-dim ebedding colored by pseudotime
#' Here pseudotime is defined as the distance from the root cell according to
#' the pseudotime metric recorded in \code{pseudotime}
#'
#' If an output dir and filename are provided, a plot will
#' be saved, otherwise just return the plot
#'
#' @param flat_embedding a low dim embedding of cells
#' @param pseudotime a scalar representation of pseudotime
#' @param outputdir the output directory, relative to getwd()
#' @param outputfile the output file
#'
#' @return a ggplot2 object
#'
#' @export
#'
PseudotimeScatterPlot <- function(flat_embedding, pseudotime, outputdir = NULL, outputfile = NULL){
  p <- ggplot2::ggplot(as.data.frame(flat_embedding), ggplot2::aes(x=V1, y=V2, color=pseudotime))

  if(!is.null(outputdir) && !is.null(outputfile)){
    # check if the dir exists and if not then create it
    if (!file.exists(outputdir)){
      dir.create(file.path(getwd(), outputdir))
    }
    file_path <- file.path(getwd(), outputdir, outputfile)
    pdf(file_path)
    p + ggplot2::geom_point()
    dev.off()
  }
  return(p)
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
#' @import tibble
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
  tibbleData <- as.tibble(markerTable)
  # sort the genes by cluster
  if(n_markers > 0){
    byCluster <- tibbleData[order(tibbleData$clusterId),]
    clusterList <- split(byCluster, byCluster$clusterId)
    topN <- lapply(clusterList, function(x){
                  head(x, n_markers)})
    gene_order <- pull(do.call(rbind, topN), geneID)
  } else {
    byCluster <- tibbleData[order(tibbleData$clusterId),]
    clusterList <- split(byCluster, byCluster$clusterId)
    topN <- lapply(clusterList, function(x){
      head(x, nrow(markerTable))})
    gene_order <- pull(do.call(rbind, topN), geneID)
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
