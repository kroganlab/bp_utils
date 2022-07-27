



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


#' Sometimes is useful to run doParallelNMF in the background so the R console remains interactive
#' This function helps with that.
#' Methods of the return object are useful for accessing results: res2$get.result()
#' 


doParallelNMF_BG <- function (data.mat, ranks = 1:20, numIterations=24, numProc=12){
  res2 <- callr::r_bg( doParallelNMF, args = list (data.mat, ranks = ranks, numIterations = numIterations, numProc = numProc))
  message ("NMF in parallel running in the background. Use proc$get_result() to access the result when finsihed. proc$read_output() to consume status lines")  
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
           if(!is.null(n)){
             writeNMFMatrices2Excel(n, ScriptAndDatedFileName(sprintf (fileNameFormat, ncol(NMF::basis(n)))))
           }
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


#' Useful for taking an MSstats RunLevelData table and converting to a matrix suitable for NMF
#' Current strategy for missing values is simply to set to zero (after converting to 0-1 linear scale)
#' This is likely not always the best case.
#' @param protQuant data.table that holds the MSstats dataProcess()$RunLevelData table
#' @param requireComplete a number, if < 1 it is taken as the portion of total columns in matrix required to be not-missing.
#' if >= 1 it is taken as the number of columns to require.
#' @param dcast.formula a two sided formula describing the rows~columns of the matrix
#' @param logIntColumn a column name to pass to dcast for value.var. Expected to be log transformed intensity
PrepareProteinLinearMat <- 
  function (protQuant, requireComplete = 0.5, dcast.formula = Protein~SUBJECT_ORIGINAL, logIntColumn  = "LogIntensities"){
  prot.mat <- as.matrix(dcast(protQuant, dcast.formula, value.var = logIntColumn), rownames = as.character(dcast.formula[[2]]))
  numNotMissing <- apply (prot.mat, 1, function(x)sum(!is.na(x)))
  #hist(numNotMissing)
  if (requireComplete < 1.0)
    requireComplete <- ncol(prot.mat) * requireComplete
  halfComplete <- numNotMissing >= requireComplete 
  halfComplete.mat <- prot.mat[halfComplete,]
  linear.mat <- 2^halfComplete.mat
  linear.mat <- sweep (linear.mat, 
                       1, 
                       apply (linear.mat, 1, max, na.rm = TRUE),
                       "/")
  # missing values...
  linear.mat[is.na(linear.mat)] <- 0.0
  return (linear.mat)
}




#' Makes use of NNLM::nnmf to estimate background contributions to individual proteins in APEX datasets
#' It does this by first defining a H0 row matrix which matches the expected contribution of a background basis
#' This H0 row matrix is calculated by assuming that the endogenous biotin carboxylases are 100% the result of background
#' and that all background will correlate with these.
#' Normalization offsets are calculated by median polish on the background-subtracted intensity matrix. Thus the result
#' will be a dataset with the backgrounds not normalized (beware!) but with the non-backgrounds normalized.
#' https://github.com/linxihui/NNLM for NNLM info/code
#' @param protQuant data.table version of MSstats::dataProcess()$RunlevelData
#' @param k integer specifying number of ranks in addition to background. Most often should be set to approximately the number of conditions (probably)
#' @param biotin.carboxylases gene names or uniprot IDS that match Protein in protQuant. 
NormalizeToNonEndogenousBiotin <- function (protQuant, k = 5, 
                                            biotin.carboxylases = c("ACACA","PC","ACACB","PCCA","MCCC1"),
                                            biotin.carboxylases.up = c("O00763","P05165","P11498","Q13085","Q96RQ3")){
  linear.mat <- PrepareProteinLinearMat(protQuant, dcast.formula = Protein~GROUP_ORIGINAL+SUBJECT_ORIGINAL)
  bcRows <- intersect (rownames(linear.mat), c(biotin.carboxylases, biotin.carboxylases.up))
  if (length(bcRows) > 0){
    message (sprintf("Using %d known background proteins as pure background indicators: ", length(bcRows)), paste0(bcRows, collapse = ", "))
  }else{
    stop("No known background proteins in protQuant")
  }
  #row matrix to match background:
  bg <- matrix(apply (linear.mat[bcRows,],2, median),
               nrow =1 )
  colnames(bg) <- colnames(linear.mat)
  print (round(100*bg))
  
  print ("Calculating NMF with fixed background coefficients")
  nmf.res <- NNLM::nnmf(linear.mat, k, n.threads = 0,
                  init = list(H0 = bg) )
  
  tmp <- nmf.res$H
  tmp <- sweep (tmp, 1, apply(tmp,1, max), "/")
  draw(Heatmap (tmp, cluster_columns = FALSE, cluster_rows = FALSE, col = c("white", "firebrick")))
  draw(Heatmap (tmp, cluster_columns = TRUE, cluster_rows = FALSE, col = c("white", "firebrick")))
  
  # here we subtract the background from the NMF fitted values.  Alternatively we could subtract from the original linear.mat
  bg.subtracted <- (nmf.res$W %*% nmf.res$H)  - (nmf.res$W[,k+1, drop = FALSE] %*% nmf.res$H[k+1,, drop = FALSE])
  # close to zero is basically zero and not trustworthy once we convert to log space (2^-10 and 2^-15 are not significantly different)
  # here we say we don't trust any log2FC > 10, or that 0.001 is same as zero.
  bg.subtracted[bg.subtracted < 2^-10] <- NA
  
  offsets <-  medpolish(log2(bg.subtracted), na.rm = TRUE)$col
  print (round(offsets, 2))
  
  message ("Modifying protQuant.  Adding new columns normalizeOffset, preNormalize and over-writing values in LogIntensities")
  protQuant[, normalizeOffset := offsets[paste0(GROUP_ORIGINAL, "_", SUBJECT_ORIGINAL)]]
  protQuant[, preNormalize := LogIntensities]
  protQuant[, LogIntensities := preNormalize - normalizeOffset]
  
  invisible(protQuant[])
}


#' Given the output of NMF::NMF, this plots the H, W matrices as annotations on originalMatrix or fitted matrix
#' relies on ComplexHeatmap 
#' @param nmf.out an NMF object
#' @param originalMatrix a matrix that was used as input for NMF.  If NULL, it will plot NMF::fitted(nmf.out)
#' @param mainColors colors (vector) or color mapping function (like from colorRamp2) to draw originalMatrix
#' @param basisColors vector of colors to be used for H and W matrices
#' @param topAnno string that selects type of annotation for the coefficient matrix (H)
#' @param ... other arguments passed to Heatmap (split_column, show_row_names, etc...)
PlotNMFHeatmap <- function (nmf.out, 
                            originalMatrix = NULL,
                            name = "linearScale",
                            mainColors = c("white", "firebrick"),
                            basisColors = RColorBrewer::brewer.pal(8, "Dark2"),
                            topAnno = c("lines", "bars", "heatmap"),
                            ...){
  coef.mat <- t(NMF::coef(nmf.out))
  rank <- ncol(coef.mat)
  colnames(coef.mat) <-   sprintf("basis.%02d", 1:rank)
  
  basis.dt <- setnames(as.data.table(NMF::basis(nmf.out)), new = sprintf("basis.%02d", 1:rank))
  
  basisColors <- rep(basisColors, ceiling(rank/length(basisColors)))[1:rank] # repeat the requested colors as necessary
  names(basisColors) <- sprintf("basis.%02d", 1:rank)
  
  mainMatrix <- originalMatrix
  if(is.null(mainMatrix)){
    mainMatrix <- NMF::fitted(nmf.out)
  }
  
  # top anno choices
  if (topAnno[1] == "lines"){
    # line annotation
    topAnnotation <- HeatmapAnnotation( coef = anno_lines(coef.mat, gp = gpar (col = basisColors),
                                                          add_points = TRUE, pt_gp = gpar(col = basisColors)),  
                                        which = "column",
                                        height = unit(2, "cm"))
  } else if (topAnno[1] == "bars"){
    # or barplot annotation
    topAnnotation <- HeatmapAnnotation( coef = anno_barplot(coef.mat, gp = gpar (fill = basisColors, col = basisColors)),  
                                        which = "column",
                                        height = unit(2, "cm"))
  } else if (topAnno[1] =="heatmap"){
    # or heatmap annotation
    upperThreshold <- quantile(coef.mat, .95)
    basisColorRamps <- lapply(basisColors, function(c)circlize::colorRamp2(breaks = c(0, upperThreshold), colors = c("white", c)) )
    coef.dt <- as.data.table(coef.mat)
    topAnnotation <-  HeatmapAnnotation( df = coef.dt, which = "column", col = basisColorRamps)
  } else {
    stop("unrecognized top annotation type")
  }
  
  # for left anno
  
  # in a group or not.
  geneGroups <- geneGroupsFromBasis(nmf.out)
  basisMembership <- geneGroups[rownames(mainMatrix), basis, on = "gene"]
  
  upperThreshold <- quantile(as.matrix(basis.dt), .99)
  basisColorRamps <- lapply(basisColors, function(c)circlize::colorRamp2(breaks = c(0, upperThreshold), colors = c("white", c)) )
  fullColorList <- c(basisColorRamps, list(basisM = basisColors))
  leftAnnotation <- HeatmapAnnotation(basisM =  basisMembership, df = basis.dt, which = "row", col = fullColorList, show_legend = FALSE, na_col = gray (0.95))
  
  
  Heatmap(mainMatrix, name = name,
          cluster_columns = FALSE, col = mainColors,
          top_annotation = topAnnotation,
          left_annotation = leftAnnotation,
          ...) 
}



geneGroupsFromBasis <- function (nmf.out, rowSumThreshold = 1.0){
  basis.mat <- NMF::basis(nmf.out)
  rowsums <- apply (basis.mat, 1, sum)
  rowsums[rowsums < rowSumThreshold] <- rowSumThreshold
  basis.mat.rowPortions <- sweep (basis.mat, 1, rowsums, FUN = "/")
  
  acceptThreshold <- 1/ncol(basis.mat.rowPortions) * 1.5
  geneGroups <- data.table (basis = sprintf("basis.%02d", apply(basis.mat.rowPortions, 1, which.max)),
                            gene = rownames(basis.mat.rowPortions),
                            portion = apply(basis.mat.rowPortions, 1, max))
  geneGroups <- geneGroups[portion > acceptThreshold]
  return (geneGroups)
}


#' relies heavily on enrichmentTestFunctions.R


characterizeNMFBasesByEnrichment <- function (nmf.out, gmt, 
                                              rowSumThreshold = 1.0,
                                              universe = NULL,
                                              basisColors = RColorBrewer::brewer.pal(8, "Dark2"),
                                              ...){
  if (!exists("enrichmentOnGroupsPL")){
    stop("This depends on functions in bp_utils/enrichmentTestFunctions.R. Be sure to source that file first")
  }
  basis.mat <- NMF::basis(nmf.out)
  rank <- ncol(basis.mat)
  colnames(basis.mat) <- sprintf("basis.%02d", 1:rank)
  
  # prepare color bars to match NMF plots
  basisColors <- rep(basisColors, ceiling(rank/length(basisColors)))[1:rank] # repeat the requested colors as necessary
  names(basisColors) <- sprintf("basis.%02d", 1:rank)
  top_annotation <- HeatmapAnnotation(basis = colnames(basis.mat), col = list(basis = basisColors), show_legend = FALSE)
  
  # rowsums <- apply (basis.mat, 1, sum)
  # rowsums[rowsums < rowSumThreshold] <- rowSumThreshold
  # basis.mat.rowPortions <- sweep (basis.mat, 1, rowsums, FUN = "/")
  # 
  # acceptThreshold <- 1/ncol(basis.mat.rowPortions) * 1.5
  # geneGroups <- data.table (basis = sprintf("basis.%02d", apply(basis.mat.rowPortions, 1, which.max)),
  #                           gene = rownames(basis.mat.rowPortions),
  #                           portion = apply(basis.mat.rowPortions, 1, max))
  # geneGroups <- geneGroups[portion > acceptThreshold]
  geneGroups <- geneGroupsFromBasis(nmf.out, rowSumThreshold)
  if (is.null(universe)){
    message("Using only rows of nmf.out as universe. To set a larger universe pass it as universe = myUniverse,")
    universe <- rownames(basis.mat)
  }
  
  
  
  
  res <- enrichmentOnGroupsPL(geneGroups, "gene", "basis", gmt, universe = universe,
                              top_annotation = top_annotation,
                              cols = colnames(basis.mat),
                              ...)
  # numProcessors = 2, max_pAdjust = 0.1, topN = 10,
  #                                  reduceRedundantsAcrossGroups = FALSE
  res$geneGroups <- geneGroups
  return (res)
}



#' Given an nmf object and the original matrix, this performs linear models on the residuals (fit-original)
#' @param nmf.out the output of NMF::NMF(originalMatrix, ...)
#' @param originalMatrix matrix should match the input that was used to generate nmf.out
#' @param parseColumnNamesFunction a function that takes a data.table with column run (from column names of originalMatirx)
#'                                 and adds columns.  See default for idea
#' @param formulas list of formulas that will be passed to my linearModelsAllFunctions. See example which uses `-1` to remove the
#'                 the intercept. This lets all coefficients be tested for difference from zero in the coef table, which is a usual
#'                 test for residuals.

testNMFResiduals <- function (nmf.out, originalMatrix,
                              parseColumnNamesFunction = function (dt)dt[, c("exp", "cell", "virus", "time", "rep") :=
                                                                           tstrsplit(run, "-")],
                              formulas = list(allVs0 = residual~virus-1),
                              ...){
  fittedMat <- NMF::fitted(nmf.out)
  stopifnot (rownames(originalMatrix) == rownames(fittedMat) )
  stopifnot (colnames(originalMatrix) == colnames(fittedMat) )
  
  residuals.long <- (originalMatrix - fittedMat ) |>
    as.data.table(keep.rownames = TRUE) |>
    melt(id.vars = "rn", variable.name = "run", value.name = "residual")
  
  parseColumnNamesFunction(residuals.long)
  
  setnames(residuals.long, old = "rn", new = "gene")
  
  lm.out <- linearModelsAllProteins(residuals.long, formulaList = formulas,splitColumn = "gene", ...)
  
  return (lm.out)
}

