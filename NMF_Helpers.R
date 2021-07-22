
#' Do a full series of NMFS with multiple iterations at each of several ranks. 
#' The best nmf result per rank is kept, where best is chosen based on minimum residuals
doParallelNMF <- function (data.mat, ranks = 1:20, numIterations=24, numProc=12){
  best.nmfs <- list()
  cl <- parallel::makeCluster(numProc)
  for (rank in ranks){
    print (sprintf("rank.%02d", rank))
    
    rep.nmf <- parallel::parLapply (cl, 1:numIterations, function(i, data, rank, .options)NMF::nmf(data, rank, .options = .options),
                                    data = data.mat, rank = rank, .options = list(parallel = FALSE))
    
    best.nmf <- rep.nmf[[which.min(sapply (rep.nmf, NMF::residuals))]]
    best.nmfs[[rank]] <- best.nmf
  }
  parallel::stopCluster(cl)
  return (best.nmfs)
}


#' Some times is useful to run doParallelNMF in the background so the R console remains interactive
#' This function helps with that.
#' Methods of the return object are useful for accessing results: res2$get.result()
#' 

doParallelNMF_BG <- function (data.mat, ranks = 1:20, numIterations=24, numProc=12){
  res2 <- callr::r_bg( doParallelNMF, args = list (data.mat, ranks = ranks, numIterations = numIterations, numProc = numProc))
  message ("NMF in parallel running in the background. Use return.value$get.result() to access the result when finsihed. return.value$read.output() to consume status lines")  
  return (res2)
}


#' Simple function to put both the basis and coefficients tables into an excel sheet
writeNMFMatrices2Excel <- function (nmf, fileName){
  coefficients <- NMF::coef(nmf)
  basis <- as.data.table (NMF::basis (nmf), keep.rownames = TRUE) 
  openxlsx::write.xlsx(list(coefficients = coefficients, basis = basis), fileName )
}



#' depends on my  ScriptAndDatedFileName function (currently in the boiler-plate of my notebooks)
writeNMFListToScriptAndDatedFileName <- function (best.nmfs, fileNameFormat = "nmf.rank%02d.xlsx"){
  lapply(best.nmfs,
         function (n){
           writeNMFMatrices2Excel(n, ScriptAndDatedFileName(sprintf (fileNameFormat, ncol(NMF::basis(n)))))
         }) %>% invisible()
}


loadNMFBasisVectors <- function (fileName) {
  table <- openxlsx::read.xlsx(fileName, sheet = "basis")
  setDT(table)
  mat <- as.matrix(table, rownames = "rn")
  colnames(mat) <- sprintf("basis.%02d", seq_len(ncol(mat)))
  return (mat)
}


# loadNMFCoefficients <- function (fileName) {
#   table <- openxlsx::read.xlsx(fileName, sheet = "basis")
#   setDT(table)
#   mat <- as.matrix(table, rownames = "rn")
#   colnames(mat) <- sprintf("basis.%02d", seq_len(ncol(mat)))
#   return (mat)
# }