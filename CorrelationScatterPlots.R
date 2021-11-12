



#' Some functions to help make pairwise correlation plots for features
#' in evidence and spectronaut output files.  Not recommended for a large number
#' of runs.  

quickUsage <- function(){
  pairsTable <- BuildPairsTable.spectronautFile("/Users/ben/Box/5HT2A data analysis/Scatter plot/HP_80uL_1mg_MNL/20210924_100025_SR01_HP_80uL_1mg_MNL_Report.xls")
  p <- PairsPlot(pairsTable)
  png(filename = ScriptAndDatedFileName("DIA_MNL_Pairwise_ScatterPlots.png"), width = 8, height = 8, units = "in", res = 200)
  p
  dev.off()
}



BuildPairsTable <- function (dt){
  # get all pairs of groups
  allByAll <- merge(dt[, .(xGroup = unique(runID), dummy = 1)],
                    dt[, .(yGroup = unique(runID), dummy = 1)],
                    allow.cartesian = TRUE)[
                      xGroup < yGroup
                    ][]
  #expand for all observed features:
  allFeatures <- unique(dt$featureID)
  allByAllFeatures <- merge (allByAll,
                             data.table (dummy = 1, featureID = allFeatures),
                             by = c("dummy"),
                             allow.cartesian = TRUE)
  
  
  # fill out the exes column
  exes <- merge (allByAllFeatures,
                 dt[, .(xGroup = runID, featureID, log2Int.x = log2Intensity)],
                 by = c("xGroup", "featureID"),
                 all.x = TRUE)
  
  #fill out they whys column
  exesAndWhys <- merge (exes,
                        dt[, .(yGroup = runID, featureID, log2Int.y = log2Intensity)],
                        allow.cartesian = TRUE,
                        by = c("yGroup", "featureID"),
                        all.x = TRUE)
  
  return (exesAndWhys[, dummy := NULL][])
}



PairsPlot <- function(pairsTable, missingXColor = "red", missingYColor = "orange"){
  minLog2Int <- min(c(pairsTable$log2Int.x, pairsTable$log2Int.y), na.rm = TRUE)
  maxLog2Int <- max(c(pairsTable$log2Int.x, pairsTable$log2Int.y), na.rm = TRUE)
  missingPlotLocation <- minLog2Int - 0.1 * (maxLog2Int - minLog2Int)
  
  p <- ggplot (pairsTable, aes(x = log2Int.x, y= log2Int.y)) +
    geom_point(alpha = 0.02, size = 1) +
    geom_abline(intercept = 0,  slope = 1, lty = "dotted", color = "gray") +
    geom_density2d() +
    facet_grid(yGroup~xGroup) +
    theme_bw() +
    ggpubr::stat_cor(r.digits = 3)
  
  if (any (is.na(pairsTable$log2Int.x))){
    #missing in one or the other
    p <- p + geom_jitter(data = pairsTable[is.na(log2Int.x)], 
                         mapping = aes(x = missingPlotLocation),
                         color = missingXColor, alpha = 0.01)
    
  }
  if (any (is.na(pairsTable$log2Int.y))){
    p <- p + geom_jitter(data = pairsTable[is.na(log2Int.y)], 
                         mapping = aes(y = missingPlotLocation),
                         color = missingYColor, alpha = 0.01)
  }
  
  return (p)
}


simplifyEvidence <- function(ev.dt){
  a <- ev.dt[, .(featureID = paste(`Leading razor protein`, `Modified sequence`, Charge, sep = "_"),
                 runID = `Raw file`,
                 log2Intensity = log2(Intensity))]
  a <- a[!grepl("(^CON)|(^REV)", featureID)]
  a <- a[, .(log2Intensity = max(log2Intensity)), by = .(runID, featureID)]
  return (a)
  
}

simplifySpectronaut <- function (spec.dt){
  a <- spec.dt[,.(featureID = paste(ProteinName, PeptideSequence, PrecursorCharge, sep = "_"),
                  runID = paste(Condition, BioReplicate, sep = "_"),
                  log2Intensity = log2(Intensity))]
  
  a <- a[!grepl("(^CON)|(^REV)", featureID)]
  a <- a[, .(log2Intensity = max(log2Intensity)), by = .(runID, featureID)]
  
}

BuildPairsTable.evidenceFile <- function(evidencePath){
  a <- simplifyEvidence(fread(evidencePath))
  return (BuildPairsTable(a))
}

BuildPairsTable.spectronautFile <- function(path){
  a <- simplifySpectronaut(fread(path))
  return (BuildPairsTable(a))
}


