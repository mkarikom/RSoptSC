csv2matrix <- function(filename,
                       foldername = "testdata", thispackage = "RSoptSC" ){
  pathstring <- system.file(foldername,filename,thispackage)
  dframe <- read.csv(pathstring)
  dmatrix <- as.matrix(dframe[-1,])
  return(dmatrix)
}
