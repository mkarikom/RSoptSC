#' Plot the signaling between individual cells
#'
#' @param P probability of signal between cell i and cell j
#' @param cluster_label the ordered vector of cluster labels for each cell
#' @param threshold the threshold for p_ij whether to report signal
#' @param sigonly default is false
#' @param plotclusternames boolean whether to plot the names
#' @param subsample boolean whether to subsample
#' @param nsample number of samples
#' @param plottitle the main title
#'
#' @return nothing
#'
#' @import circlize
#' 
#' @export
#'
plotSig <- function(P, 
                    cluster_label, 
                    threshold, 
                    sigonly=FALSE, 
                    plotclusternames=TRUE, 
                    subsample, 
                    nsample, 
                    outfile="Psig.pdf", 
                    plottitle){ 
  
  if(!sigonly){
    Pdata = thresholdP(P, threshold)
  }
  
  cell_id = seq(1,dim(Pdata)[1])
  cid_out = rep(0, length(cell_id))
  for(i in 1:length(cell_id)){
    x = which(Pdata[i,]!=0)
    if(length(x) > 0){
      cid_out[i] = 1
    }
  }
  cid_in = rep(0, length(cell_id))
  for(i in 1:length(cell_id)){
    x = which(Pdata[,i]!=0)
    if(length(x) > 0){
      cid_in[i] = 1
    }
  }
  cid_out_ids = which(cid_out==1)
  cid_in_ids = which(cid_in==1)
  
  edge_pairs = matrix(nrow=0, ncol=2) 
  for(i in 1:length(cid_out_ids)){
    ee = which( Pdata[cid_out_ids[i],] > 0)
    pairs = cbind(rep(cid_out_ids[i], length(ee)), ee)  
    edge_pairs = rbind(edge_pairs, pairs)
  }
  
  ## Create dataframe from signaling probabilities
  edges = data.frame(cluster_from = cluster_label[edge_pairs[,1]], cell_from = edge_pairs[,1], 
                     cluster_to = cluster_label[edge_pairs[,2]], cell_to = edge_pairs[,2]
  )
  weights = rep(0, dim(edges)[1])
  for(i in 1:dim(edges)[1]){
    weights[i] = Pdata[edges[i,"cell_from"], edges[i,"cell_to"]]
  }
  edges$weight = weights
  
  
  if(subsample){
    edgesample = sample(1:dim(edges)[1], nsample, replace=FALSE)
    df = data.frame(edges[edgesample,], stringsAsFactors=FALSE)
  }
  else{        
    df = data.frame(edges, stringsAsFactors=FALSE)
  }
  
  if(dim(df)[1] > 200){    
    warning('\n Warning: The signaling network has > 200 edges; this is too big to visualize. 
            Either increase the edge threshold or subsample fewer edges. \n\n')
    return(df)
  }
  
  
  ## Start plot
  circos.clear()  
  
  cluster = c(structure(df$cluster_from, names=df$cell_from), structure(df$cluster_to,names= df$cell_to))
  cluster = cluster[!duplicated(names(cluster))]
  cluster = cluster[order(cluster, names(cluster))]
  
  gap.after = do.call("c", lapply(table(cluster), function(i) c(rep(2, i-1), 5)))
  circos.par(start.degree = 230, gap.after = gap.after, cell.padding = c(0, 0, 0, 0))
  
  col_fun = colorRamp2(range(Pdata), c("#CCCCFF", "#4D4DFF"))
  
  chordDiagram(df[, c(2,4,5)], order = names(cluster), directional = 1,
               col = col_fun,
               direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", diffHeight = uh(5,"mm"),
               annotationTrack = "grid",
               preAllocateTracks = list(list(track.height = 0.05))
  )
  
  circos.track(track.index = 2, panel.fun = function(x, y){
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), mean(ylim), labels='', col = "white", cex = 0.2, facing = "inside", niceFacing = TRUE)
    highlight.sector(sector.index, col = "#E6E6E6", border="#FFFFFF") # cell track color (grey) 
  }, bg.border = NA
  )
  
  cluster_order = c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12')
  ctext = setnames(plotclusternames, cluster_order)
  
  for(b in unique(cluster)){
    cell = names(cluster[cluster == b])
    highlight.sector(sector.index = cell, track.index = 1,
                     col = clswitch(cluster_order[b]),
                     text = ctext[b],
                     text.vjust = -1, niceFacing = TRUE
    )
  }
  
  title(plottitle)
}

clswitch <- function(label) {
    transparency = ""
    switch(label,
           "C1" = paste0("#0072BD", transparency),
           "C2" = paste0("#D95319", transparency),
           "C3" = paste0("#EDB120", transparency),
           "C4" = paste0("#7E2F8E", transparency),
           "C5" = paste0("#77AC30", transparency),
           "C6" = paste0("#4DBEEE", transparency),
           "C7" = paste0("#006400", transparency),
           "C8" = paste0("#A2142F", transparency),
           "C9" = paste0("#282828", transparency),
           "C10" = paste0("#cc00cc", transparency),
           "C11" = paste0("#006500", transparency),
           "C12" = paste0("#000080", transparency)
          )
}

thresholdP <- function(P, threshold){
    data = P
    for(i in 1:dim(data)[1]){
        for(j in 1:dim(data)[2]){
        
            if(data[i,j] < threshold){
                data[i,j] = 0.0
            }
        }
    }
    return(data)
}

setnames <- function(plot_clusternames, cluster_order){
    if(plot_clusternames){
        return(cluster_order)
    }
    else{ 
        return('')
    }
}