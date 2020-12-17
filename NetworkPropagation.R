

LoadNumPyS_matrix <- function (matrixPath, nodesTablePath){
  library (RcppCNPy)
  S_matrix <- RcppCNPy::npyLoad(matrixPath)
  nodeNames <- fread (nodesTablePath)
  stopifnot (nrow(S_matrix) == nrow(nodeNames))
  dimnames(S_matrix) <- list (nodeNames$gene_name, nodeNames$gene_name)
  return(S_matrix)
}



# geneHeats : data.table with gene and heat columns
NetworkPropagateS_matrix <- function(S_matrix, geneHeats,numPermutations = 20000, networkHeatOnly = FALSE, permuteOnlyInObserved=FALSE){
  message (now(), " Checking S matrix")
  stopifnot(check_s_matrix(S_matrix))
  
  if(not(all(geneHeats$gene %in% rownames(S_matrix)))){
    numHeats <- nrow(geneHeats)
    geneHeats <- geneHeats[gene %in% rownames(S_matrix)]
    message (sprintf("Not all genes with heat are found in network; heats reduced from %d to %d", numHeats, nrow(geneHeats)))
  }
  
  message (now(), " Propagating initial heat")
  
  heat.0 <- makeHeat.0(geneHeats, fullGenes = rownames(S_matrix))
  propagatedHeat <- S_matrix %*% heat.0
  # confirm conservation of heat: make sure total heat before and after propagation matches
  stopifnot ( abs(sum(propagatedHeat) - sum(heat.0))  < 0.001)  # allow rounding errors
  
  heatToBeat <- propagatedHeat
  if (networkHeatOnly){
    #remove self-heat
    heatToBeat <- propagatedHeat - heat.0 * diag(S_matrix)
  }
  
  #permute analysis
  message (now(), " Building permutation matrix")
  permutedColumns <- buildPermutedColumns(geneHeats, rownames(S_matrix), numPermutations, limitToObserved = permuteOnlyInObserved)
  message (now(), " Propagating permutation matrix")
  permutedHeats <- S_matrix %*% permutedColumns
  
  # remove self-heat if requested to allow comparison of just network-propagated heats
  if (networkHeatOnly) permutedHeats <- permutedHeats - permutedColumns * diag(S_matrix)#sweep (permutedColumns, 1, diag(S_matrix), FUN = "*")
  
  message (now(), " Summarizing to p-values")
  # in next two lines, we count (permute-propagate) > 0 and divide by total number of permutes
  deltaPermutedHeats <- sweep (permutedHeats, 1, heatToBeat) # permute - propagate
  pValue <- rowSums(deltaPermutedHeats>0)/ncol(deltaPermutedHeats) # count(()>0) /total
  
  results <- data.table (gene = rownames(S_matrix),
                         heat.0 = heat.0,
                         prop.heat = propagatedHeat[,1],
                         pvalue = pValue,
                         adj.pvalue = p.adjust(pValue, method="BH"),
                         self.heat = heat.0 * diag(S_matrix))
  
  return(results)
}



now <- function(){
  format(Sys.time(), "%H:%M:%S")
} 

check_s_matrix <- function(S_matrix){
  # cols should sum to 1.  allow rounding error
  colSums <- colSums(S_matrix)
  error <- max(colSums-1)
  return (error < 0.001)
}

makeHeat.0 <- function (geneHeats, fullGenes){
  #expand heats to full set in the network, as indicated by names on inv_denom
  fullHeats <- geneHeats[fullGenes, .(gene, heat), on = "gene"]
  fullHeats[is.na(heat), heat := 0]
  heat.0 <- fullHeats$heat
  names(heat.0) <- fullHeats$gene
  return(heat.0)
}

buildPermutedColumns <- function(geneHeats, fullGenes, numPermutations, limitToObserved = FALSE){
  if(limitToObserved  & (length(unique(geneHeats$heat)) == 1)){
    message ("limitToObserved = TRUE requires a full set of observed heats, significant and not", 
             "\n\tgeneHeats currently only has a single heat value",
             "\n\trunning with limitToObserved is set to FALSE")
    limitToObserved <- FALSE
  }
  if (limitToObserved){
    heat0 <- makeHeat.0(geneHeats, fullGenes)
    observedIndex <- match(geneHeats$gene, fullGenes)
    x0_perm_mat <- do.call(cbind,lapply(1:numPermutations, FUN = function(i){heat0[observedIndex] <- sample(heat0[observedIndex]); return(heat0)}))
    rownames(x0_perm_mat) <- names(heat0)
  }else {
    heat0 <- makeHeat.0(geneHeats, fullGenes)
    x0_perm_mat <- do.call(cbind,lapply(1:numPermutations, FUN = function(i)sample(heat0)))
    rownames(x0_perm_mat) <- names(heat0)
  }
  return (x0_perm_mat)
}

