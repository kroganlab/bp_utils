library(data.table)
# also requires RcppCNpy for LoadNumPyS_matrix



LoadNumPyS_matrix <- function (matrixPath, nodesTablePath){
  library (RcppCNPy)
  S_matrix <- RcppCNPy::npyLoad(matrixPath)
  nodeNames <- fread (nodesTablePath)
  stopifnot (nrow(S_matrix) == nrow(nodeNames))
  dimnames(S_matrix) <- list (nodeNames$gene_name, nodeNames$gene_name)
  return(S_matrix)
}



# idea here is to propagate two or more datasets at a time.  
# Heat is the product of each heat propagated individually.
#
# Theory is that this will best identify the "linker" genes between two datasets
#' @param geneHeatsList list (ab = ab.geneHeats, ppi = ppi.geneHeats), for example
#' each ab.geneHeats etc is a data.table with columns `gene` and `heat`
NetworkPropagate.multiHeats.S_matrix <- function(S_matrix, geneHeatsList,numPermutations = 20000, networkHeatOnly = TRUE, 
                                                 permuteOnlyInObserved=TRUE, calculateContributions = FALSE, genesInContributionsTable = NULL){
  message (now(), " Checking S matrix")
  stopifnot(check_s_matrix(S_matrix))
  
  message (now(), " Rescaling input heat to size of S_matrix (expected heat per gene = 1.0)")
  # lapply for its data.table side effects by :=
  lapply(geneHeatsList, function(geneHeats)
    geneHeats[, heat := nrow(S_matrix) * heat/sum(heat, na.rm=TRUE)])
  
  # expected heat outs is row sums of S_matrix. Assumes all input genes are equally likely to have heat.
  # expectedHeatOut <- rowSums(S_matrix)
  # if (networkHeatOnly){
  #   expectedHeatOut <- expectedHeatOut - diag(S_matrix)
  # }
  # 
  
  # expected Heat will vary based on coverage of dataset, calculate per dataset:
  expectedHeatOut <- sapply(geneHeatsList, calculateExpectedHeat, S_matrix, permuteOnlyInObserved, networkHeatOnly)
  


  # check that all heats are in network
  for (heatName in names(geneHeatsList)){
    # rescale input heat to size of S_matrix, thus expected heat per gene is 1.0
    geneHeats <- geneHeatsList[[heatName]]
    if(!all(geneHeats$gene %in% rownames(S_matrix))){
      numHeats <- nrow(geneHeats)
      geneHeats <- geneHeats[gene %in% rownames(S_matrix)]
      message (sprintf("%s: Not all genes with heat are found in network; heats reduced from %d to %d", heatName, numHeats, nrow(geneHeats)))
      geneHeatsList[[heatName]] <- geneHeats
    }
  }
  
  # check for duplicate heats
  for (heatName in names(geneHeatsList)){
    geneHeats <- geneHeatsList[[heatName]]
    duplicates <- geneHeats[, .N, by = gene][N > 1, ]$gene
    if(length(duplicates) > 0){
      message (sprintf("%s: Multiple heats found for %d genes; will take the max; be sure this is expected:", heatName, length(duplicates)),
               paste0(head(duplicates, 6), collapse = ","), 
               ifelse(length(duplicates) > 6, ",...", ""))
      geneHeats <- geneHeats[, .(heat = max(heat, na.rm = TRUE)), by = gene]
      geneHeatsList[[heatName]] <- geneHeats
    }
  }
  
  
  # make initial heat matrix (each column is one source of heat)
  heat.0.mat <- sapply(geneHeatsList, makeHeat.0, fullGenes = rownames(S_matrix))
  
  message (now(), " Propagating initial heat")
  propagatedHeat.mat <- S_matrix %*% heat.0.mat
  

  if (networkHeatOnly){
    #remove self-heat before taking product across columns (by row)
    propagatedHeat.mat <- propagatedHeat.mat - sweep(heat.0.mat,1,diag(S_matrix), '*')
  }
  
  # transform heats by log2FC (cut off less than zero)
  stopifnot (all(dim(propagatedHeat.mat) == dim(expectedHeatOut)))
  log2PropHeat <- log2(propagatedHeat.mat) - log2(expectedHeatOut)
  log2PropHeat[log2PropHeat < 0.0] <- 0.0
  
  heatToBeat <- apply (log2PropHeat, 1, prod)
  
  
  # permute
  # define a list of vectors. Each vector is the items in each data source (columns) to permute
  heatIndcesToPermute <- lapply (geneHeatsList,
                                 function(gh){
                                   if(permuteOnlyInObserved == FALSE  | (length(unique(gh$heat)) == 1)){
                                     return (1:nrow(S_matrix))
                                   }else{
                                     return (which(rownames(S_matrix) %in% gh$gene))
                                   }
                                 })

  # all permutations at once, create a 3d Array, gene * permutation * dataType
  
  #permute analysis
  message (now(), " Building permutation matrix")
  permutation.mat.list <- pbapply::pblapply (geneHeatsList, buildPermutedColumns, fullGenes= rownames(S_matrix), numPermutations = numPermutations, limitToObserved = permuteOnlyInObserved)

  #permutedColumns <- buildPermutedColumns(geneHeats, rownames(S_matrix), numPermutations, limitToObserved = permuteOnlyInObserved)
  message (now(), " Propagating permutation matrix")
  permutedHeats.mat3d <- pbapply::pbsapply (permutation.mat.list, function(permutedColumns)S_matrix %*% permutedColumns, simplify = "array")
  #permutedHeats <- S_matrix %*% permutedColumns
  
  # remove self-heat if requested to allow comparison of just network-propagated heats
  if (networkHeatOnly) {
    permutedHeats.mat3d <- permutedHeats.mat3d - sapply(permutation.mat.list, function(permutedColumns)permutedColumns * diag(S_matrix), simplify = "array") #sweep (permutedColumns, 1, diag(S_matrix), FUN = "*")
  }
  
  # done with permutation.mat.list, free up its memory
  rm (permutation.mat.list)
  gc()

  message (now(), " Transforming permutations to log2FC vs expected")
  permutedHeats.mat3d <- sweep (log2(permutedHeats.mat3d), c(1,3), log2(expectedHeatOut))
  permutedHeats.mat3d[permutedHeats.mat3d < 0] <- 0.0
  
  # 3rd dimension of permutedHeats.mat3d is data type.  collapse/integrate across data types by taking the product
  message (now(), " Integrating permutations by product")
  
  #the obvious way of doing it is very slow.
  #all.permuted <- pbapply::pbapply(permutedHeats.mat3d, c(1,2), prod)
  
  #faster to do use rowSums on a log transformed array:
  log.h <- log2(permutedHeats.mat3d)
  log.ap <- rowSums(log.h, dims = 2)
  all.permuted <- exp(log.ap)
  
  # done with permutedHeats.mat3d, free up the memory
  rm(permutedHeats.mat3d, log.ap, log.h)
  gc()
  
  # count number of times we beat or exceed the heatToBeat (by rows)
  message (now(), " Converting to p values")
  countsAbove <- rowSums (sweep(all.permuted, 1, heatToBeat) >= 0)
  pValue <- countsAbove/ncol(all.permuted)
  
  results <- data.table (gene = rownames(S_matrix),
                         heat0 = heat.0.mat,
                         expectH = expectedHeatOut,  # a matrix, will get colnames
                         propH = propagatedHeat.mat,
                         avgLog2FC.PvsE = heatToBeat^(1/ncol(heat.0.mat)),  # undo the product, effectively the geometric mean now
                         pvalue = pValue,
                         adj.pvalue = p.adjust(pValue, method="BH") )

  #cbind (results, heat.0.mat, propagatedHeat.mat)
  return (results)
}



# geneHeats : data.table with gene and heat columns
NetworkPropagateS_matrix <- function(S_matrix, geneHeats,numPermutations = 20000, networkHeatOnly = FALSE, 
                                     permuteOnlyInObserved=FALSE, calculateContributions = FALSE, genesInContributionsTable = NULL){
  message (now(), " Checking S matrix")
  stopifnot(check_s_matrix(S_matrix))
  
  if(!all(geneHeats$gene %in% rownames(S_matrix))){
    numHeats <- nrow(geneHeats)
    geneHeats <- geneHeats[gene %in% rownames(S_matrix)]
    message (sprintf("Not all genes with heat are found in network; heats reduced from %d to %d", numHeats, nrow(geneHeats)))
  }
  
  duplicates <- geneHeats[, .N, by = gene][N > 1, ]$gene
  if(length(duplicates) > 0){
    message (sprintf("Multiple heats found for %d genes; will take the max; be sure this is expected:", length(duplicates)),
             paste0(head(duplicates, 6), collapse = ","), 
             ifelse(length(duplicates) > 6, ",...", ""))
    geneHeats <- geneHeats[, .(heat = max(heat, na.rm = TRUE)), by = gene]
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
  # deltaPermutedHeats <- sweep (permutedHeats, 1, heatToBeat) # permute - propagate
  # pValue <- rowSums(deltaPermutedHeats>0)/ncol(deltaPermutedHeats) # count(()>0) /total
  

  # maybe more memory friendly
  countsAbove <- rep(0, nrow(permutedHeats))
  for (i in 1:ncol(permutedHeats)){
    countsAbove <- countsAbove + as.integer(permutedHeats[,i] > heatToBeat)
  }
  pValue <- countsAbove/ncol(permutedHeats)
  
  
  results <- data.table (gene = rownames(S_matrix),
                         heat.0 = heat.0,
                         prop.heat = propagatedHeat[,1],
                         pvalue = pValue,
                         adj.pvalue = p.adjust(pValue, method="BH"),
                         self.heat = heat.0 * diag(S_matrix))
  
  if (calculateContributions == TRUE){
    if (!is.null(genesInContributionsTable)){
      sigGenes <-genesInContributionsTable 
      message (now(), sprintf(" Calculating per-gene contributions to final heat for %d requested genes", length(genesInContributionsTable)))
    }else{
      sigGenes <- results[pvalue < 0.05,]$gene
      message (now(), " Calculating per-gene contributions to final heat for genes with p.value < 0.05")
    }
    toFromTable <- CalculateContributions (S_matrix, heat.0, sigGenes, networkHeatOnly)
    return(list(results = results, contributions = toFromTable))
  }else{
    return(results)
  }
}


CalculateContributions <- function(S_matrix, startHeats, sigGenes, networkHeatOnly){
  if (length(startHeats) != nrow(S_matrix)){
    duplicates <- startHeats[, .N, by = gene][N > 1, ]$gene
    if(length(duplicates) > 0){
      message (sprintf("Multiple heats found for %d genes; will take the max; be sure this is expected:", length(duplicates)),
               paste0(head(duplicates, 6), collapse = ","), 
               ifelse(length(duplicates) > 6, ",...", ""))
      startHeats <- startHeats[, .(heat = max(heat, na.rm = TRUE)), by = gene]
    }
    
    
  }
    startHeats <- makeHeat.0(startHeats, fullGenes = rownames(S_matrix))
  
  
  stopifnot (all (names(startHeats) == rownames(S_matrix)))
  
  fromToMat <- t(S_matrix[sigGenes, ]) * startHeats  #this works based on recycling the startHeats vector over t(S), column-wise.  We can safely subset rows, but not columns.
  fromToTable <- melt(as.data.table(fromToMat, keep.rownames= TRUE), id.vars = "rn", variable.name = "to", value.name = "contributed.heat", )
  setnames(fromToTable, old = "rn", new = "from")
  if (networkHeatOnly == TRUE)
    fromToTable[to == from, contributed.heat := 0.0]
  fromToTable[, percent.contributed.heat := 100 * contributed.heat/sum(contributed.heat, na.rm = TRUE), by = to]
  setorder (fromToTable, -percent.contributed.heat, na.last = TRUE)
  fromToTable[, cumulative.percent.contributed.heat := cumsum(percent.contributed.heat), by = to]
  return (fromToTable[contributed.heat > 0.0]) # trim table to positive contributions only (no 0.0)
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


calculateExpectedHeat <- function (geneHeats, S_matrix, limitToObserved, networkHeatOnly){
  if(limitToObserved  & (length(unique(geneHeats$heat)) == 1)){
    # should be warned, but don't warn here, they'll get another warning below when buildPermutedColumns
    limitToObserved <- FALSE
  }
  
  rcIndex <- which (rownames(S_matrix) %in% geneHeats$gene)
  
  # all calculations below are based on an input heat of 1.0 per gene, on average.
  # if actual average (across whole network) input heat is higher or lower, need to adjust the final answers
  scaleFactor <- sum(geneHeats$heat)/nrow(S_matrix)

  if (limitToObserved){
    expectedHeatOut <- rowSums(S_matrix[, rcIndex])
    if (networkHeatOnly){
      expectedHeatOut[rcIndex] <- expectedHeatOut[rcIndex] - diag(S_matrix)[rcIndex]
    }
  }else{
    # expected heat outs is row sums of S_matrix. Assumes all input genes are equally likely to have heat.
    expectedHeatOut <- rowSums(S_matrix)
    if (networkHeatOnly){
      expectedHeatOut <- expectedHeatOut - diag(S_matrix)
    }
  }
  return(expectedHeatOut * scaleFactor)
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

