#' Load sequencing data and annotations
#' 
#' A matrix of counts with rows as genes and columns as cells must be provided.  Gene and cell ids are provided in separate text files.  For labeled data, the vector must be provided as a separate file.  The counts are loaded as a sparse matrix, the other values as character vectors.
#'
#' @param d_f path to the file containing counts
#' @param gene_f path the file containing unique gene ids 
#' @param cell_f path the file containing unique cell ids 
#' @param annotation_f path to the file containing labels, e.g. tissue type of the cells
#'
#' @return a list containing:
#'     \item{cell_names}{cell ids}
#'     \item{gene_names}{eg gene symbols}
#'     \item{data}{counts matrix}
#'     \item{annotation}{cell labels}
#'
#' @importFrom Matrix Matrix
#' @importFrom utils read.csv
#' 
#' @export
#'
LoadData <- function(d_f, 
                     gene_f, 
                     cell_f,
                     annotation_f){
  datas <- read.csv(d_f, row.names = 1)
  genes <- read.csv(gene_f,row.names = 1)
  cells <- read.csv(cell_f,row.names = 1)
  annotations <- read.csv(annotation_f,row.names = 1)
  rownames(datas) <- genes[,1]
  colnames(datas) <- cells[,1]
  list(cell_names = as.vector(cells[,1]),
       gene_names = as.vector(genes[,1]),
       annotation = as.vector(annotations[,1]),
       data = Matrix(data = as.matrix(datas), sparse = TRUE))
}
