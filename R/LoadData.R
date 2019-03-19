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
  datas <- read.csv(d_f)
  genes <- read.csv(gene_f)
  cells <- read.csv(cell_f)
  annotations <- read.csv(annotation_f)
    
  list(cell_names = as.vector(cells[,2]),
       gene_names = as.vector(genes[,2]),
       annotation = as.vector(annotations[,2]),
       data = Matrix(data = as.matrix(datas)[,-c(1)], sparse = TRUE))
}
