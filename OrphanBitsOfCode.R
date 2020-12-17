
# Random bits of code

# reusable bits of code will be dumped here before finding a better home



# take a data.table of msstats results and define significance
# useful for tabulating types of significance and coloring volcano plots

defineSig <- function(resultsDT, p.threshold = 0.05, log2FC.threshold = 1, pvalueCol = c("pvalue", "adj.pvalue")[1]){
  resultsDT[,sig := "not"]
  resultsDT[is.finite (log2FC) &
              abs(log2FC) > log2FC.threshold &
              resultsDT[[pvalueCol]] < p.threshold,
            sig := ifelse(log2FC > 0, "up", "down")]
  resultsDT[is.infinite(log2FC),
            sig := ifelse(log2FC > 0, "infiniteUp", "infiniteDown")]
  return (resultsDT)
}


rankVolcanoHits <- function(resultsDT, log2FC.t = 1, pvalue.t = 0.05){
  resultsDT[is.finite(log2FC),volcanoDistance := sqrt(log2FC^2 + log10(pvalue)^2)]
  #resultsDT[,rankVolcano := NA]
  resultsDT[,rankVolcano := NULL]
  resultsDT[pvalue < pvalue.t & abs(log2FC) > log2FC.t, rankVolcano := rank(-volcanoDistance), by = .(Label, sig)]
  invisible(resultsDT)
}


# not a working function at the moment, mostly a place to store the code...
ggplotVolcano <- function (resultsDT, labelPattern = "CD74|UBE2M|APOBEC3G|STAU1|MACROH2A1|UBE2N|RNF213"){
  
  
  p<-ggplot(defineSig(abResults, pvalueCol = "pvalue"), aes(x = log2FC, y = -log10(pvalue), col = sig, fill=sig, label= geneLabel)) + 
    geom_hline(yintercept = -log10(0.05), lty="dashed", lwd=0.5, col = gray(0.2)) + geom_vline(xintercept = c(-1,1), lty="dashed", lwd=0.5, col = gray(0.2)) +
    geom_point(show.legend = FALSE, alpha = 0.5, size = 1, shape = 21, stroke=0) +
    scale_color_manual(values = c(up = "red", down="blue", not = gray(0.2), infiniteUp = "firebrick", infiniteDown = "navy")) +
    scale_fill_manual(values = c(up = "red", down="blue", not = "gray", infiniteUp = "firebrick", infiniteDown = "navy")) +
    ggforce::facet_wrap_paginate(~Label, ncol = 1,nrow=1, page = i) +
    #facet_wrap(~Label, ncol=2) +
    ggrepel::geom_text_repel(data = rankVolcanoHits(abResults)[(grepl ("CD74|UBE2M|APOBEC3G|STAU1|MACROH2A1|UBE2N|RNF213", gene) & sig %in% c("up", "down")) |
                                                                 (rankVolcano <= numToLabel & sig != "not" & is.finite(log2FC)) ],
                             show.legend = FALSE,min.segment.length=0, force=4, size=2.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print (p)
  
}

# Most useful for understanding whats happening with the linear models in MSstats
#
# using the results of dataProcess (protInt) and a selection of genes, make plots
# of the actual data fitted value overall and fitted per SUBJECT_ORIGINAL
# will probably only work with bioreplicate-matched analysis
# see file 2020_09_22_Paige_Donor_Matched_Analylsis.Rmd for an example.
#
# 
plotLMERfit <- function(protInt, geneIDs = c("RAB11FIP5_K583", "CDK6_K257","MYO1G_K87","ENO1_K193", "ATP5F1B_K201", "ENO1_K60", "ENO1_K120", "WDR1_K81","PKM_K66","PKM_K270")){
  for (geneID in geneIDs){
    
    l <- lmer (LogIntensities~GROUP_ORIGINAL+(1|SUBJECT_ORIGINAL)-1, data = protInt[gene %in% c(geneID)])
    
    
    fitted <- data.table (LogIntensities = summary(l)$coef[,1],  GROUP_ORIGINAL = gsub("^GROUP_ORIGINAL", "", names(summary(l)$coef[,1])), SUBJECT_ORIGINAL = "average")
    donorEffects <- data.table (coef = ranef(l)$SUBJECT_ORIGINAL[,1], donor = rownames (ranef(l)$SUBJECT_ORIGINAL))
    donorEffects <- fitted[, .(LogIntensities = LogIntensities + donorEffects[,coef], SUBJECT_ORIGINAL = donorEffects$donor ), by = GROUP_ORIGINAL]
    
    print(ggplot(protInt[gene %in% c(geneID)], aes(x = GROUP_ORIGINAL, y = LogIntensities, color = SUBJECT_ORIGINAL,group = SUBJECT_ORIGINAL)) + 
            geom_point() + 
            #geom_line(lty="dotted") + 
            geom_line(data=fitted) +
            geom_line(data = donorEffects, lty = "dotted")+
            ggtitle(label = geneID))
  }
  
}


# t tests and some cool visualizations for pairwise comparisons
# see https://cran.r-project.org/web/packages/emmeans/vignettes/comparisons.html
contrastsViaEmmeans <- function( protInt, geneID){
  l <- lmer (LogIntensities~GROUP_ORIGINAL+(1|SUBJECT_ORIGINAL)-1, data = protInt[gene %in% c(geneID)])
  emmeans::contrast (emmeans::emmeans(l, "GROUP_ORIGINAL"), method="pairwise", adjust = NULL)
  
  p <- plot(emmeans::emmeans(l, "GROUP_ORIGINAL") ) + theme_bw() + coord_flip()
  print (p)
  
  p <- emmeans::pwpp(emmeans::emmeans(l, "GROUP_ORIGINAL")) + theme_bw()
  print (p)
}



# PCA####



completeMat.pepIntensity <- function(pepInt, normalize=TRUE){
  expectedCols <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "Condition", "Run", "Intensity")
  if (!all( expectedCols %in% colnames(pepInt))){
    stop("Expected a table of peptide Intensity in MSstats input format with columns: ", paste(expectedCols, collapse=", "))
  }
  # convert to Matrix
  pep.wide <- dcast (pepInt, paste(ProteinName, PeptideSequence, PrecursorCharge, sep= ".")~Condition+Run, sep="-", value.var = "Intensity", fun.aggregate = max)
  pep.mat <- as.matrix(pep.wide, rownames = "ProteinName")
  pep.mat[pep.mat == 0] <- NA
  completeMat <- pep.mat[complete.cases(pep.mat),]
  completeMat <- log2(completeMat)
  
  message (sprintf("Data set reduced to %d peptide ions of %d in original by requiring complete cases", nrow(completeMat), nrow(pep.wide)))
  
  if (normalize){
    #normalize by centering medians
    colMedians <- apply (completeMat, 2, median, na.rm=TRUE)
    globalMedian <- median(colMedians)
    message ("Normalizing with adjustments in range ", paste(round(range (colMedians-globalMedian), 2), collapse = " to "))
    completeMat <- scale(completeMat, center = colMedians-globalMedian, scale=FALSE)
  }
  return (completeMat)
}

PCA.pepIntensity <- function(pepInt, doPrint  = TRUE, normalize = TRUE, 
                             runLabels = NULL, conditionLabels = NULL, conditionShapes = NULL,
                             colorSet = "Set2", title = "",
                             brewerPalette = "Spectral"){
  
  completeMat <- completeMat.pepIntensity (pepInt, normalize)
  
  #PCA
  pcaOut <- prcomp(t(completeMat))
  pcaDT <- as.data.table(pcaOut$x, keep.rownames=TRUE)
  
  # handle annotations, using passed in labels if present
  rowInfo <- unique (pepInt[, .(Condition, Run, BioReplicate)])
  rowInfo[,rn := paste(Condition, Run, sep="-")]
  
  if (is.null(runLabels)){
    rowInfo[,runLabel := BioReplicate]
  } else {
    stopifnot (all (rowInfo$Run %in% runLabels$Run) )
    rowInfo[runLabels, runLabel := label,on= "Run"]
  }
  
  if (is.null(conditionLabels)){
    rowInfo[,conditionLabel := Condition]
  } else {
    stopifnot (all (rowInfo$Condition %in% conditionLabels$Condition) )
    rowInfo[conditionLabels, conditionLabel := label,on= "Condition"]
  }
  
  if (is.null(conditionShapes)){
    rowInfo[,shape := 1] 
  }else{
    stopifnot (all (rowInfo$Condition %in% conditionShapes$Condition) )
    rowInfo[conditionShapes, shape := i.shape,on= "Condition"]
  }
  
  pcaDT <- merge (pcaDT, rowInfo, by = c("rn"))
  
  pcaPercentVar <- round(100 * (pcaOut$sdev^2)/sum(pcaOut$sdev^2), 1)
  
  #define colors
  conditions <- sort(unique(rowInfo$conditionLabel))
  conditionColors <- RColorBrewer::brewer.pal(length(conditions), brewerPalette)
  names(conditionColors) <- conditions
  
  
  #plot first two components
  p <- ggplot (pcaDT, aes(x=PC1, y=PC2, col=conditionLabel, shape = shape)) + 
    geom_point(alpha=1.0, size=4) + 
    ggrepel::geom_text_repel(aes(label=runLabel), show.legend = FALSE) +
    theme_bw() + 
    xlab (sprintf ("PC1, %.1f%%", pcaPercentVar[1])) + 
    ylab (sprintf ("PC2, %.1f%%", pcaPercentVar[2])) + 
    ggtitle (sprintf ("PCA %s using %d peptides (log intensity)", title, nrow(completeMat))) +
    scale_color_manual(values = conditionColors, name = "Condition")
  
  if (doPrint) print(p)
  invisible(p)
}



#  ARGUMENTS
# evidenceDT: a data.table holding the contents of an evidence file
# globalStandards: a character vector of one or more protein IDs
# subsetRuns: a character vector of runs to include.  Useful if evidence file includes extra runs, in which case the result will exclude the extra runs. Leave NULL or empty to ignore. 

#
# USAGE :
# library (data.table)
# evidenceDT <- fread ("evidence.txt")
# globalStandards <- "P38398"
# cleanEv <- removeIncompleteGlobalStandardFeatures(evidenceDT, globalStandards)
# fwrite (cleanEv, "clean.evidence.txt")


removeIncompleteGlobalStandardFeatures <- function(evidenceDT, globalStandards, subsetRuns=NULL){
  # make a copy or subset as needed
  if (!is.null(subsetRuns) & length(subsetRuns) > 0){
    message ("Will only return rows in subsetRuns")
    evDT <- evidenceDT[`Raw file` %in% subsetRuns]
  }else{
    evDT <- copy (evidenceDT)
  }
  
  #the global standard features
  gsFeatures <- evDT[Proteins %in% globalStandards]
  
  # the set of full runs
  fullRuns <- unique( evDT$`Raw file` )
  
  #presence table
  presenceTable <- dcast (gsFeatures,`Modified sequence`+Charge+Proteins~`Raw file`,
                          value.var = "Intensity",
                          # insert NAs here for easy use of complete.casese below
                          fun.aggregate = function(x)(if(length(x) > 0) length(x) else NA) )
  
  # did we get all runs:
  if (length(setdiff(fullRuns, colnames(presenceTable))) > 0){
    stop("global standards completely missing from ", paste0(setdiff(fullRuns, colnames(presenceTable)), collapse=", "))
  }
  
  complete <- presenceTable[complete.cases(presenceTable)]
  if (nrow(complete) == 0){
    stop ("no global standard features are complete in all runs")
  }
  
  message (sprintf("%d features are complete in all %d runs", nrow(complete), length(fullRuns)))
  completeFeatureDT <- gsFeatures[complete[,.(`Modified sequence`, Charge, Proteins)],, on = c("Modified sequence", "Charge", "Proteins")]
  message (sprintf ("Of %d global standard rows in evidence file, keeping %d global standard rows", nrow(gsFeatures), nrow(completeFeatureDT)))
  
  return (rbind (completeFeatureDT, evDT[!Proteins %in% globalStandards]))
}







