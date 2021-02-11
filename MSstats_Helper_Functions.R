


subsetDataProcessOutput <- function(mssquant.full, groups, newRunLevelData = NULL){
  subset <- list()
  subset$ProcessedData <- data.table(mssquant.full$ProcessedData)[GROUP_ORIGINAL %in% groups,]
  if (is.null(newRunLevelData)){
    subset$RunlevelData <- data.table(mssquant.full$RunlevelData)[GROUP_ORIGINAL %in% groups,]
  } else{
    subset$RunlevelData <- data.table(newRunLevelData)[GROUP_ORIGINAL %in% groups,]
  }
  
  # # RunlevelData might need columns in a specific order:
  # subset$RunlevelData <- subset$RunlevelData[,c("RUN", "Protein", "LogIntensities", "NumMeasuredFeature", "MissingPercentage", 
  #                                               "more50missing", "NumImputedFeature", "originalRUN", "GROUP", "GROUP_ORIGINAL",
  #                                               "SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT"),
  #                                            with = FALSE
  #                                            ]
  
  notEntirelyMissing <- subset$RunlevelData[,.(numNotMissing = sum(!is.na(LogIntensities))), by = Protein][numNotMissing > 0,]$Protein
  subset$RunlevelData <- subset$RunlevelData[Protein %in% notEntirelyMissing,]
  subset$RunlevelData[,Protein := factor(Protein)]
  subset$ProcessedData <- subset$ProcessedData[PROTEIN %in% unique(subset$RunlevelData$Protein)]
  subset$ProcessedData[,PROTEIN := factor(PROTEIN)]
  
  subset$SummaryMethod <- mssquant.full$SummaryMethod
  subset$ModelQC <- mssquant.full$ModelQC
  subset$PredictBySurvival <- mssquant.full$PredictBySurvival
  
  .fixLevels <- function(dt){
    # fix levels in GROUP
    dt[,GROUP_ORIGINAL := factor(GROUP_ORIGINAL)]
    dt[,GROUP := factor(as.character(as.integer(GROUP_ORIGINAL)), 
                        levels = as.character(sort(unique(as.integer(GROUP_ORIGINAL)))))]
    
    # SUBJECT
    dt[,SUBJECT_ORIGINAL := factor(SUBJECT_ORIGINAL)]
    dt[,SUBJECT := factor(as.character(as.integer(SUBJECT_ORIGINAL)), 
                          levels = as.character(sort(unique(as.integer(SUBJECT_ORIGINAL)))))]
    #SUBJECT_NESTED
    dt[,SUBJECT_NESTED := factor(paste(GROUP, SUBJECT, sep="."))]
    
    # RUN
    dt[,originalRUN := factor(originalRUN)]
    dt[,RUN := factor(as.character(as.integer(originalRUN)), 
                      levels = as.character(sort(unique(as.integer(originalRUN)))))]
    
  }
  
  .fixLevels(subset$ProcessedData)
  .fixLevels(subset$RunlevelData)
  stopifnot (all(levels (subset$ProcessedData$GROUP_ORIGINAL) == levels (subset$RunlevelData$GROUP_ORIGINAL)))
  
  setDF(subset$ProcessedData)
  setDF(subset$RunlevelData)
  return (subset)
}




makeContrast.regEx <- function(mssQ, regEx){
  columnNames <- as.character(levels(mssQ$RunlevelData$GROUP_ORIGINAL))
  
  positives = c()
  negatives = c()
  for (i in columnNames){
    for (j in columnNames){
      if (i != j){
        positives <- c(positives,i)
        negatives <- c(negatives,j)
      }
    }
  }
  contrasts <- data.table (positives, negatives)
  contrasts[,name := paste(positives,negatives, sep="-")]
  
  contrastNames <- character(0)
  for (re in regEx){
    contrastNames <- union(contrastNames, grep(re, contrasts$name, value = TRUE))
  }
  print (contrastNames)
  contrasts <- contrasts[name %in% contrastNames]
  print (contrasts)
  
  contrastMat <- matrix(rep(0,length(columnNames) * nrow(contrasts)), nrow=nrow(contrasts), ncol=length(columnNames),
                        dimnames = list(contrasts$name, columnNames))
  
  for (i in seq_len(nrow(contrasts))){
    #rows are named according to the positive condition:
    contrastMat[contrasts$name[i],contrasts$negatives[i]] = -1
    contrastMat[contrasts$name[i],contrasts$positives[i]] =  1
  }
  return (contrastMat)  
}
