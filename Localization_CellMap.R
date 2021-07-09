library (data.table)




#' @param listOfSets boolean indicating to return the full data.table (FALSE, default) or a list of sets
#' @param identifiers columns from cellMap table to use as identifiers, only used if listOfSets == TRUE
loadCellMap <- function (listOfSets = FALSE, identifiers = c("symbol", "UniProt"),
                         fileOrUrl = "https://humancellmap.org/resources/downloads/preys-latest.txt"){
  message ("Reading cellMap info from ", fileOrUrl )
  cellMap <- fread (fileOrUrl)
  if(listOfSets ){
    return (cellMapTable2ListOfSets(cellMap, identifiers))
  }
  return (cellMap)
}

cellMapTable2ListOfSets <- function (cellMap, identifiers = c("symbol", "UniProt")){
  location2ID <- cellMap[, .(ids = unlist(.SD)), by = `MMF localization`, .SDcols = identifiers][(!is.na(ids)) & ids != "",] 
  sets <- split(location2ID$ids, location2ID$`MMF localization`)
  return (sets)
}


#' Given a matrix and a column name perform fgsea on a single column
#' sets a named list of genes/proteins vectors. Each vector is a single set. vector items must match rownames of mat
#' scoreType, usually "std" for positive and negative or "pos" for positive only
oneColumnFGSEA <- function (columnName, mat, sets, scoreType = c("std", "pos", "neg"), ...){
  print (columnName)
  log2FC <- mat[,columnName]
  # reorder randomly so ties are less likely to influence results
  log2FC <- sample(log2FC, length(log2FC))
  
  #nperm=1000, gseaWeightParam = 1, nproc=1
  seaRes <- fgsea::fgsea(pathways = sets, stats = log2FC, gseaParam=1, ...)
  seaRes[, sigScore := -log10(pval) * ifelse(ES < 0, -1, 1) ]
  return (seaRes)
}

#' Given a matrix  fgsea on all columns, return table of results with column "group" identifying the column of the matrix
#' @param sets a named list of genes/proteins vectors. Each vector is a single set. vector items must match rownames of mat
#' @param scoreType is chosen based on presence of negative or positive values in the matrix
#' @param ... arguments to pass to oneColumnFGSEA and ultimately fgsea::fgsea

matrixFGSEA <- function (mat, sets, ...){
  anyPositive <- any(mat > 0)
  anyNegative <- any(mat < 0)
  if (anyPositive & anyNegative)
    scoreType = "std"
  if (anyPositive & !anyNegative)
    scoreType = "pos"
  if (!anyPositive & anyNegative)
    scoreType= "neg"
  if(!anyPositive & !anyNegative)
    stop("No variance or no values detected in matrix")
  
  all.sea <- lapply (colnames(mat), oneColumnFGSEA, mat, sets, scoreType, ...)
  names(all.sea) <- colnames(mat)
  sea.dt <- rbindlist(all.sea, idcol = "group")
  
  return (sea.dt)
}


#' adds a column to scores.dt that contains cellMap localization
setLocateScoresTable <- function(scores.dt, cellMap.dt = NULL, idCol = "Protein", cellMap.idCol = "UniProt"){
  if (is.null(cellMap.dt))
    cellMap.dt <- loadCellMap()
  
  names(cellMap.idCol) <- idCol
  scores.dt[cellMap.dt, location := `MMF localization`, on = cellMap.idCol] # c(Protein = "UniProt")
  
  numLocated <- scores.dt[is.na(location), length(unique(.SD[[idCol]])) ]
  numNotLocated <- scores.dt[!is.na(location), length(unique(.SD[[idCol]])) ]
  sprintf("%d proteins located and %d not located in cellMap. The not located include mismatched identifiers (if too low, check that identifiers are correct)",
          numLocated, numNotLocated)
  
  invisible(scores.dt)
}



cellMapLocalizationScores <- function(scores.dt, groupCol = "Label", scoreCol = "log2FC", idCol = "Protein", type = c("UNIPROT", "SYMBOL")){
  cellMap.dt <- loadCellMap()
  cellMap.idCol <- c(UNIPROT = "UniProt",
                     SYMBOL = "symbol")[
                       toupper(type[1])]
  cellMap.sets <- cellMapTable2ListOfSets(cellMap.dt,
                                          identifiers = cellMap.idCol)
  
  # convert scores to matrix for fgsea analysis
  formula <- as.formula (sprintf ("%s~%s", idCol, groupCol))
  wide.dt <- dcast (scores.dt, formula, value.var = scoreCol)
  scores.mat <- as.matrix (wide.dt, rownames = idCol )
  
  #FGSEA
  print (sprintf ("Starting FGSEA on scores in %d groups, (%d values)", ncol(scores.mat), nrow(scores.mat)))
  sea.dt <- matrixFGSEA (scores.mat, cellMap.sets)
  
  #combine localization scores into scores.dt
  # per-protein localizaitons
  setLocateScoresTable (scores.dt, cellMap.dt, idCol,cellMap.idCol)
  # per group FGSEA scores
  scores.dt[sea.dt, sea.sigScore := sigScore , on = c(setNames("group", groupCol), location = "pathway")]

  invisible(list(scores.dt = scores.dt, sea.dt = sea.dt))
}

#' @param gridFormula optional, something like receptor~time to see receptors in rows and time in columns. if null then .~groupCol

violinsAndScatterLocations <- function (scores.dt, scoreCol = "log2FC", groupCol = "Label", xlimits = NULL, gridFormula = NULL, reorder = TRUE){
  if (reorder){
    sea.mat <- summarizesSEAMatrixFromScoresDT(scores.dt)
    orderIdx <- dist(sea.mat) %>% hclust %>% as.dendrogram %>% order.dendrogram
    scores.dt[, location := factor(location, levels = rownames(sea.mat)[orderIdx]) ]
  }
  
  
  
  if (is.null(gridFormula))
    gridFormula <- as.formula (sprintf (".~%s", groupCol))
  p <- ggplot (scores.dt[!is.na(location)],
               aes(y = location, x = .data[[scoreCol]], fill = sea.sigScore )) +
    geom_vline(xintercept = 0) +
    geom_jitter(alpha = 0.1, width = 0.0, height = 0.1) +
    geom_violin(lwd = 0.2) + 
    facet_grid(gridFormula ) +
    #scale_x_continuous(limits = c(-2,2)) +
    coord_cartesian(xlim = xlimits) + 
    scale_fill_gradient2(low = "blue", high = "red") +
    theme_bw()
  return (p)
}



#' First 
summarizesSEAMatrixFromScoresDT <- function (scores.dt, groupCol = "Label"){
  formula <- as.formula (sprintf ("location~%s", groupCol))
  as.matrix (dcast (scores.dt[!is.na(sea.sigScore)], formula, value.var = "sea.sigScore", fun.aggregate = max),
             rownames = "location")
}


usage <- function(results.txt){
  results <- fread (result.txt)
  source ("../../bp_utils/Localization_CellMap.R")
  result.list <- cellMapLocalizationScores(results)
  p <- violinsAndScatterLocations(results)
  
}

