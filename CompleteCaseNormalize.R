

# peptideDF is expected to be a table ready to pass to MSstats::dataProcess
# it will be modified in place using data.table functions

CompleteCaseNormalize.inPlace <- function(peptideDF, minIntensityQuantile = 0){
  if (!"data.table" %in% class (peptideDF)){
    setDT(peptideDF)
    data.frame.only <- TRUE
  } else{
    data.frame.only <- FALSE
  }
  
  minIntensity <- quantile(peptideDF$Intensity, minIntensityQuantile)
  
  #make a dt of complete cases, by passing through a wide format
  pepDF.wide <- dcast(peptideDF[Intensity > minIntensity & !is.na(Intensity)],
                      ProteinName+PeptideSequence+PrecursorCharge~BioReplicate+Condition+Run, value.var="Intensity")
  completeFeatures <- pepDF.wide[complete.cases(pepDF.wide), .(ProteinName, PeptideSequence, PrecursorCharge)]
  message ("Normalizing based on intensities in ", nrow(completeFeatures), " complete features (having intensity > ", minIntensity, " ")
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


# as above, but for evidence files
CompleteCaseNormalize.MaxQuant.inPlace <- function(peptideDF, minIntensityQuantile = 0){
  if (!"data.table" %in% class (peptideDF)){
    setDT(peptideDF)
    data.frame.only <- TRUE
  } else{
    data.frame.only <- FALSE
  }
  
  minIntensity <- quantile(peptideDF$Intensity, minIntensityQuantile, na.rm=TRUE)
  
  #make a dt of complete cases, by passing through a wide format
  pepDF.wide <- dcast(peptideDF[Proteins != "" & !grepl("CON|REV", Proteins) & Intensity > minIntensity & !is.na(Intensity)],
                      Proteins+`Modified sequence`+Charge~`Raw file`, value.var="Intensity", fun.aggregate = sum)


  completeFeatures <- pepDF.wide[,.(anyZero = any (unlist(.SD) == 0)), by = .(Proteins, `Modified sequence`, Charge)
                                 ][anyZero == FALSE, .(Proteins, `Modified sequence`, Charge)]
  message ("Normalizing based on intensities in ", nrow(completeFeatures), " complete features (having intensity > ", minIntensity, " ")
  setkey (peptideDF, Proteins, `Modified sequence`, Charge)
  completeLong <- peptideDF[!is.na(Intensity)][completeFeatures]
  
  overallMedian <- median (log2(completeLong$Intensity))
  adjustments <- completeLong[,.(adjustment = overallMedian - median(log2(Intensity))),by = .(`Raw file`)]
  print(adjustments)
  
  setkey(peptideDF, `Raw file`)
  
  #normalization:
  
  peptideDF[adjustments, Intensity := Intensity * 2^adjustment]
  
  if (data.frame.only)
    setDF(peptideDF)
  
  invisible(peptideDF)
}




NearlyCompleteCaseNormalize.inPlace <- function(peptideDF, require = 0.50, minIntensityQuantile = 0){
  if (!"data.table" %in% class (peptideDF)){
    setDT(peptideDF)
    data.frame.only <- TRUE
  } else{
    data.frame.only <- FALSE
  }
  
  minIntensity <- quantile(peptideDF$Intensity, minIntensityQuantile)
  
  
  # make a dt of complete cases, by passing through a wide format
  totalRuns <- nrow(unique(peptideDF[,.(BioReplicate,Condition,Run)]))
  minRunsReq <- ceiling(totalRuns * require)
  passingFeatures <- peptideDF[!is.na(Intensity) & Intensity > minIntensity,
                               .N,
                               by = .(ProteinName,PeptideSequence,PrecursorCharge)
                               ][N >= minRunsReq,.(ProteinName,PeptideSequence,PrecursorCharge) ]

  message (nrow(passingFeatures), " features present in at least ", minRunsReq, " runs for normalization (having intensity > ", minIntensity, " ")
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


