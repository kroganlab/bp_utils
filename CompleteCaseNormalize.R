

# peptideDF is expected to be a table ready to pass to MSstats::dataProcess
# it will be modified in place using data.table functions

CompleteCaseNormalize.inPlace <- function(peptideDF){
  if (!"data.table" %in% class (peptideDF)){
    setDT(peptideDF)
    data.frame.only <- TRUE
  } else{
    data.frame.only <- FALSE
  }
  
  # make a dt of complete cases, by passing through a wide format
  pepDF.wide <- dcast(peptideDF, ProteinName+PeptideSequence+PrecursorCharge~BioReplicate+Condition+Run, value.var="Intensity")
  completeFeatures <- pepDF.wide[complete.cases(pepDF.wide), .(ProteinName, PeptideSequence, PrecursorCharge)]
  message ("Normalizing based on intensities in ", length(completeFeatures), " complete features")
  setkey (peptideDF, ProteinName, PeptideSequence, PrecursorCharge)
  completeLong <- peptideDF[completeFeatures]
  
  overallMedian <- median (log2(completeLong$Intensity))
  adjustments <- completeLong[,.(adjustment = overallMedian - median(log2(Intensity))),by = .(Condition,  BioReplicate)]
  print(adjustments)
  
  setkey(peptideDF, Condition, BioReplicate)
  
  #normalization:
  
  peptideDF[adjustments, Intensity := Intensity * 2^adjustment]
  
  if (data.frame.only)
    setDF(peptideDF)

  invisible(peptideDF)
}



NearlyCompleteCaseNormalize.inPlace <- function(peptideDF, require = 0.50){
  if (!"data.table" %in% class (peptideDF)){
    setDT(peptideDF)
    data.frame.only <- TRUE
  } else{
    data.frame.only <- FALSE
  }
  
  # make a dt of complete cases, by passing through a wide format
  totalRuns <- nrow(unique(peptideDF[,.(BioReplicate,Condition,Run)]))
  minRunsReq <- ceiling(totalRuns * require)
  passingFeatures <- peptideDF[!is.na(Intensity),.N, by = .(ProteinName,PeptideSequence,PrecursorCharge)][N >= minRunsReq,.(ProteinName,PeptideSequence,PrecursorCharge) ]

  message (nrow(passingFeatures), " features present in at least ", minRunsReq, " runs for normalization.")
  completeLong <- peptideDF[passingFeatures, , on = .(ProteinName, PeptideSequence, PrecursorCharge)]
  
  overallMedian <- median (log2(completeLong$Intensity), na.rm=TRUE)
  adjustments <- completeLong[,.(adjustment = overallMedian - median(log2(Intensity), na.rm=TRUE)),by = .(Condition,  BioReplicate)]
  print(adjustments)
  
  setkey(peptideDF, Condition, BioReplicate)
  
  #normalization:
  
  peptideDF[adjustments, Intensity := Intensity * 2^adjustment]
  
  if (data.frame.only)
    setDF(peptideDF)
  
  invisible(peptideDF)
  
  
}


