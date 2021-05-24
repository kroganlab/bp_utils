

library (data.table)
library (ComplexHeatmap)
library (ggplot2)


# A generic function that should work on any table provided it has columns:
#   runID, proteinID, featureID, logIntensity, 
# standards :  a vector of proteinIDs to normalize with
# doPlots : make plots that show progress.  Not recommended with standards = NULL
 
normalizeByMedianPolish <- function(featureIntensitiesDT, standards = NULL, doPlots = TRUE ){
  if (!is.null(standards)){
    normSubset <- featureIntensitiesDT[proteinID %in% standards]
    if (nrow(normSubset) == 0){
      stop("No features match requested standard")
    }else{
      cat ( sprintf("Normalizing using the %d features in protein(s): %s\n", normSubset[,length(unique(featureID))], paste0(standards, collapse = ",")))
    }
  }else
    normSubset <- featureIntensitiesDT
  
  if (any (normSubset[, .N, by = .(featureID, runID)]$N > 1)){
    message ("Some features to normalize on appear multiple times per run. Will use the maximum intensity peak only.")
    normSubset <- normSubset[, .SD[which.max(logIntensity)], by = .(featureID, runID)]
  }
  
  
  featureInt.mat <- as.matrix (dcast (normSubset, featureID~runID, value.var = "logIntensity", na.rm = TRUE),
                               rownames = "featureID")
  
  # plotFeatures <- rownames(featureInt.mat)
  # if (length(plotFeatures) > plotMaxFeatures)
  #   plotFeatures <- sample(plotFeatures, plotMaxFeatures)
  
  if(doPlots){
    p <- ggplot(normSubset, aes(x = runID, y = logIntensity, color = featureID)) +
      geom_point(show.legend = FALSE) + 
      geom_line(aes(group = featureID), show.legend = FALSE) + 
      theme(axis.text.x = element_text(angle = 90)) +
      ggtitle(sprintf("Features before normalization, %s", paste0(standards, collapse = ",")))
    print(p)
  }
  
  cat ("median polish iteration: sum absolute residuals\n")
  mp <- medpolish(featureInt.mat,  na.rm = TRUE)
  
  if (doPlots){
    rowAnno <- rowAnnotation(pep = anno_barplot(mp$row))
    colAnno <- HeatmapAnnotation(run = anno_barplot(mp$col))
    
    mp.r <- mp$residuals
    mp.r.noNA <- mp.r
    mp.r.noNA[is.na(mp.r.noNA)] <- 0.0
    ddr <- as.dendrogram(hclust(dist(mp.r.noNA)))
    hm <- Heatmap (mp.r, cluster_columns = FALSE,
                   name = "residuals",
                   cluster_rows = ddr,
                   row_names_gp = gpar(fontsize = 7),
                   col = circlize::colorRamp2(breaks = c(-2,0,2), colors = c("blue", "white", "red")),
                   top_annotation = colAnno, right_annotation = rowAnno,
                   column_title = "MP Residuals (center) and pep/run effects (barplots)")
    draw (hm)
  }
  

  offsets.table <- data.table(runID = names(mp$col), offset = mp$col)
  
  # do a merge-modify into a new variable 
  featureIntensitiesDT[offsets.table, normLogIntensity := logIntensity - offset,on = "runID"]
  
  if (doPlots){
    if (!is.null(standards)) {
      normSubset <- featureIntensitiesDT[proteinID %in% standards]
    }else{
      normSubset <- featureIntensitiesDT
    }
    normSubset <- normSubset[, .SD[which.max(logIntensity)], by = .(featureID, runID)]
    
    p <- ggplot(normSubset, aes(x = runID, y = normLogIntensity, color = featureID)) +
      geom_point(show.legend = FALSE) + 
      geom_line(aes(group = featureID), show.legend = FALSE) + 
      theme(axis.text.x = element_text(angle = 90))+
      ggtitle(sprintf("Features after normalization, %s", paste0(standards, collapse = ",")))
    print(p)
    
    p <- ggplot(featureIntensitiesDT, aes(x = runID)) + 
      geom_boxplot (aes(y = logIntensity), color = "black") + 
      geom_boxplot(aes(y = normLogIntensity), color = "red", alpha = 0.5) + 
      theme(axis.text.x = element_text(angle = 90)) +
      ggtitle ("Before (black) and after(red) normalization, all features")
    print (p)
  }
  invisible (featureIntensitiesDT)
}



# Normalize a maxquant evidence file. 
# inFilePath : the path to the evidence file
# standards  :  a vector of proteinIDs to normalize with. Must match IDs in "Leading razor protein" column
# doPlots    : makes plots that show progress, to the default graphics device.  Not recommended with standards = NULL
# outFilePath : where to write the output.  If not included it will be a file with similar name as inFilePath
normalizeByMedianPolish.evidenceFile <- function(inFilePath, standards = NULL, doPlots = FALSE, outFilePath = NULL){
  if (is.null(outFilePath))
    outFilePath <- gsub ("(\\.txt)?$", ".MP_normalized.txt", inFilePath)

  ev <- fread (inFilePath, integer64= "double")
  colNames <- c("runID", "proteinID", "featureID", "logIntensity")
  # avoid over-writing pre-existing columns
  stopifnot (!any(colNames %in% colnames(ev)))
  # define columns that normalizeByMedianPolish depends on
  ev[, c("runID", "proteinID", "featureID", "logIntensity") :=
       .(`Raw file`,
         `Leading razor protein`,
         paste(`Leading razor protein`, `Modified sequence`, Charge, sep = "_"),
         log2(Intensity))
  ]
  
  ev.norm <- normalizeByMedianPolish (ev, standards = standards, doPlots = doPlots)
  
  # clean up
  ev.norm[, Intensity.prenormalization := Intensity]
  ev.norm[, Intensity := 2^normLogIntensity]
  ev.norm[, c(colNames, "normLogIntensity") := rep(NULL, length(colNames))]
  
  cat (sprintf("Writing normalized evidence file to %s\n", outFilePath))
  fwrite (ev.norm, outFilePath, sep = "\t")
  invisible (ev.norm)
}







