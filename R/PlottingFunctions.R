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
