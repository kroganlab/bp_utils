

# Modularization of Gene Sets based on network connections


require(data.table)





# functions to choose best representative term based on enricher results
termScore <- function (GeneRatio, BgRatio){
  gr <- lapply(tstrsplit(GeneRatio, "/"), as.integer)
  br <- lapply(tstrsplit(BgRatio, "/"), as.integer)
  
  return (gr[[1]]**2/(gr[[2]] * br[[1]]))
  
}

enricherJaccard <- function (GeneRatio, BgRatio){
  gr <- lapply(tstrsplit(GeneRatio, "/"), as.integer)
  br <- lapply(tstrsplit(BgRatio, "/"), as.integer)
  
  # jaccard is intersection/union
  
  return (gr[[1]]/ # intersection is genes in group AND term
            (br[[1]]  # union is all genes that match the term...
             + gr[[2]] # plus the genes in group...
             - gr[[1]]) # minus those already counted that match the term 
  )
}






