context("signaling partners")
library(RSoptSC)
library(Matrix)
library(dplyr)
library(reshape2)

expressionfn <- system.file("extdata", "test_gene_expression.txt", package = "RSoptSC")
ligrecfn <- system.file("extdata", "test_ligrec.txt", package = "RSoptSC")
rectargfn <- system.file("extdata", "test_rectarg_updown.txt", package = "RSoptSC")
test_that("files are found", {
  expect_true(nchar(expressionfn) > 0)
  expect_true(nchar(ligrecfn) > 0)
  expect_true(nchar(rectargfn) > 0)
})

test_that("up down both work", {
  expression = read.table(expressionfn, header = T, sep = ",")
  ligrec = read.table(ligrecfn, header = T, sep = ",", dec = ",")
  rectarg = read.table(rectargfn, header = T, sep = ",")
  labels = c(1,2,1)
  
  data.frame(id=c(1,2,3),cluster=c(1, 2, 1));
  
  pathway = RSoptSC::ImportPathway(lig_table = ligrec, 
                                   rec_table = rectarg, 
                                   data = expression)
  P = GetSignalingPartners(expression,rownames(expression),pathway,normalize_aggregate=T)

  PClust = ClusterSig(P,labels,normalize_rows = T)
  pLig1Rec1 = matrix(0,3,3)
  #lig1,rec1, up
  Lcell1_Rcell1_up = exp(-1/(1*2))/(exp(-1/(1*2))+exp(-2/3))*exp(-2/3) ## unnormalized
  Lcell1_Rcell2_up = exp(-1/(1*1))/(exp(-1/(1*1))+exp(-2/4))*exp(-2/4) ## unnormalized
  Lcell1_Rcell3_up = exp(-1/(1*2))/(exp(-1/(1*2))+exp(-2/3))*exp(-2/3) ## unnormalized
  
  Lcell2_Rcell1_up = exp(-1/(2*2))/(exp(-1/(2*2))+exp(-2/3))*exp(-2/3) ## unnormalized
  Lcell2_Rcell2_up = exp(-1/(2*1))/(exp(-1/(2*1))+exp(-2/4))*exp(-2/4) ## unnormalized
  Lcell2_Rcell3_up = exp(-1/(2*2))/(exp(-1/(2*2))+exp(-2/3))*exp(-2/3) ## unnormalized
  
  Lcell3_Rcell1_up = exp(-1/(3*2))/(exp(-1/(3*2))+exp(-2/3))*exp(-2/3) ## unnormalized
  Lcell3_Rcell2_up = exp(-1/(3*1))/(exp(-1/(3*1))+exp(-2/4))*exp(-2/4) ## unnormalized
  Lcell3_Rcell3_up = exp(-1/(3*2))/(exp(-1/(3*2))+exp(-2/3))*exp(-2/3) ## unnormalized
  
  # lig1, rec1, down
  Lcell1_Rcell1_down = exp(-1/(1*2))/(exp(-1/(1*2))+exp(-8/1))*exp(-8/1) ## unnormalized
  Lcell1_Rcell2_down = exp(-1/(1*1))/(exp(-1/(1*1))+exp(-9/1))*exp(-9/1) ## unnormalized
  Lcell1_Rcell3_down = exp(-1/(1*2))/(exp(-1/(1*2))+exp(-8/1))*exp(-8/1) ## unnormalized
  
  Lcell2_Rcell1_down = exp(-1/(2*2))/(exp(-1/(2*2))+exp(-8/1))*exp(-8/1) ## unnormalized
  Lcell2_Rcell2_down = exp(-1/(2*1))/(exp(-1/(2*1))+exp(-9/1))*exp(-9/1) ## unnormalized
  Lcell2_Rcell3_down = exp(-1/(2*2))/(exp(-1/(2*2))+exp(-8/1))*exp(-8/1) ## unnormalized
  
  Lcell3_Rcell1_down = exp(-1/(3*2))/(exp(-1/(3*2))+exp(-8/1))*exp(-8/1) ## unnormalized
  Lcell3_Rcell2_down = exp(-1/(3*1))/(exp(-1/(3*1))+exp(-9/1))*exp(-9/1) ## unnormalized
  Lcell3_Rcell3_down = exp(-1/(3*2))/(exp(-1/(3*2))+exp(-8/1))*exp(-8/1) ## unnormalized
  
  
  pLig1Rec1[1,1] = exp(-1/(1*2))*Lcell1_Rcell1_up*Lcell1_Rcell1_down
  pLig1Rec1[1,2] = exp(-1/(1*1))*Lcell1_Rcell2_up*Lcell1_Rcell2_down
  pLig1Rec1[1,3] = exp(-1/(1*2))*Lcell1_Rcell3_up*Lcell1_Rcell3_down
  
  pLig1Rec1[2,1] = exp(-1/(2*2))*Lcell2_Rcell1_up*Lcell2_Rcell1_down
  pLig1Rec1[2,2] = exp(-1/(2*1))*Lcell2_Rcell2_up*Lcell2_Rcell2_down
  pLig1Rec1[2,3] = exp(-1/(2*2))*Lcell2_Rcell3_up*Lcell2_Rcell3_down
  
  pLig1Rec1[3,1] = exp(-1/(3*2))*Lcell3_Rcell1_up*Lcell3_Rcell1_down
  pLig1Rec1[3,2] = exp(-1/(3*1))*Lcell3_Rcell2_up*Lcell3_Rcell2_down
  pLig1Rec1[3,3] = exp(-1/(3*2))*Lcell3_Rcell3_up*Lcell3_Rcell3_down
  
  nrml = rowSums(pLig1Rec1)
  pLig1Rec1[,1] = pLig1Rec1[,1] / nrml
  pLig1Rec1[,2] = pLig1Rec1[,2] / nrml
  pLig1Rec1[,3] = pLig1Rec1[,3] / nrml
  
  checkcell = all(signif(P$P[[1]]$value,2) == signif(reshape2::melt(pLig1Rec1)$value,digits=2))
  expect_true(checkcell)
  
  clusttotal = PClust$P[[1]] %>% dplyr::group_by(cluster.Var1) %>% dplyr::summarize(check = sum(value))
  checkclust = all(clusttotal$check == 1)
  expect_true(checkclust)
})
