


subsetDataProcessOutput <- function(mssquant.full, groups = NULL, proteins = NULL, subjects = NULL, newRunLevelData = NULL){
  subset <- list()
  subset$ProcessedData <- data.table(mssquant.full$ProcessedData)[ (is.null(groups) | GROUP_ORIGINAL %in% groups) &
                                                                     (is.null(proteins) | PROTEIN %in% proteins) &
                                                                     (is.null(subjects) | SUBJECT_ORIGINAL %in% subjects),]
  if (is.null(newRunLevelData)){
    subset$RunlevelData <- data.table(mssquant.full$RunlevelData)[(is.null(groups) | GROUP_ORIGINAL %in% groups) &
                                                                    (is.null(proteins) | Protein %in% proteins)&
                                                                    (is.null(subjects) | SUBJECT_ORIGINAL %in% subjects),]
  } else{
    subset$RunlevelData <- data.table(newRunLevelData)[(is.null(groups) | GROUP_ORIGINAL %in% groups) &
                                                         (is.null(proteins) | Protein %in% proteins)&
                                                         (is.null(subjects) | SUBJECT_ORIGINAL %in% subjects),]
    
    #make the subjects match based on RUN to allow for renaming of subjects for different nested design
    rld <- unique(subset$RunlevelData[,.(RUN, originalRUN, GROUP, GROUP_ORIGINAL, SUBJECT_ORIGINAL, SUBJECT_NESTED, SUBJECT)])
    subset$ProcessedData[rld, c("SUBJECT_ORIGINAL", "SUBJECT_NESTED", "SUBJECT") := .(i.SUBJECT_ORIGINAL, i.SUBJECT_NESTED, i.SUBJECT), on = "RUN"]
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


# artMS writes the output of dataProcess in two files,  results-mss-normalized.txt and
# results_RunLevelData.txt. This function loads these files and formats/packages them in a list
# so that MSstats can continue with groupComparison
loadProcessedDataFromArtMS <- function (results.mss.normalized.txt, results_RunLevelData.txt, groupSubset = NULL){
  normalizedPep <- fread(results.mss.normalized.txt)
  proteinQuant <- fread (results_RunLevelData.txt)
  
  #"ProcessedData"     "RunlevelData"      "SummaryMethod"     "ModelQC"           "PredictBySurvival"
  mssquant <- list (ProcessedData  = normalizedPep,
                    RunlevelData = proteinQuant,
                    SummaryMethod = "TMP",  # the most likely guess...I don't think it matters for downstream analysis
                    ModelQC = NULL,
                    PredictBySurvival = NULL)
  
  if (is.null(groupSubset))
    groupSubset <- unique(mssquant$RunlevelData$GROUP_ORIGINAL) # all groups

  # regarless of subsetting, we use this function to make everything a factor that needs to be
  mssquant <- subsetDataProcessOutput (mssquant, groupSubset)
  return (mssquant)
}


#makes use of artMS internal function to load an artMS formatted contrasts file
makeContrast.artMSFile <- function (contrasts.txt){
  contrast.mat <- artMS:::.artms_writeContrast(contrasts.txt)
  return (contrast.mat)
}


# 
checkVersionDataProcessOutput <- function (mssQ){
  if (!"RunlevelData" %in% names(mssQ)){
    stop("Wrong version of MSstats. Source MSstats_V4_Functions.R and try again")
  }
}


makeContrast.regEx <- function(mssQ, regEx){
  checkVersionDataProcessOutput(mssQ)
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

makeContrast.AllByAll <- function(mssQ){
  checkVersionDataProcessOutput(mssQ)
  columnNames <- as.character(levels(mssQ$RunlevelData$GROUP_ORIGINAL))
  
  positives = c()
  negatives = c()
  for (i in columnNames){
    for (j in columnNames){
      if(i < j){
        positives <- c(positives,i)
        negatives <- c(negatives,j)
      }
    }
  }
  contrasts <- data.table (positives, negatives)
  contrasts[,name := paste(positives,negatives, sep="-")]
  
  contrastMat <- matrix(rep(0,length(columnNames) * nrow(contrasts)), nrow=nrow(contrasts), ncol=length(columnNames),
                        dimnames = list(contrasts$name, columnNames))
  
  for (i in seq_len(nrow(contrasts))){
    #rows are named according to the positive condition:
    contrastMat[contrasts$name[i],contrasts$negatives[i]] = -1
    contrastMat[contrasts$name[i],contrasts$positives[i]] =  1
  }
  return (contrastMat)  
}




#' This function takes the output of MSstats::dataProcess with SUBJECTS like GroupA.1 and GroupB.1 and converts
#' both to match batch.1.  This will tell MSstats to do a "paired" analysis, aka include a SUBJECT term in the model.
#' This will only work if the second field when splitting by "." denotes a meaningful batch identifier. 
#' @param dp.out The output list of MSstats::dataProcess
#' @param data.frame.only logical ,if TRUE will make return values data.frames (instead of data.tables)

batchifyMSQuant <- function(dp.out, data.frame.only = NULL){
  if(is.null(data.frame.only))
    data.frame.only <- !"data.table" %in% class(dp.out$RunlevelData)
  
  setDT(dp.out$ProcessedData)
  dp.out$ProcessedData[, SUBJECT_ORIGINAL := paste("batch", tstrsplit(SUBJECT_ORIGINAL, split = "\\.")[[2]], sep = ".")]
  dp.out$ProcessedData[, SUBJECT := as.integer(as.factor(SUBJECT_ORIGINAL))]
  dp.out$ProcessedData[, SUBJECT_NESTED := sprintf("%d.%d", GROUP,SUBJECT)]
  if (data.frame.only)
    setDF(dp.out$ProcessedData)
  
  setDT(dp.out$RunlevelData)
  dp.out$RunlevelData[, SUBJECT_ORIGINAL := paste("batch", tstrsplit(SUBJECT_ORIGINAL, split = "\\.")[[2]], sep = ".")]
  dp.out$RunlevelData[, SUBJECT := as.integer(as.factor(SUBJECT_ORIGINAL))]
  dp.out$RunlevelData[, SUBJECT_NESTED := sprintf("%d.%d", GROUP,SUBJECT)]
  if (data.frame.only)
    setDF(dp.out$RunlevelData)
  
  return (dp.out)
}





specFileToCompleteMSstats <- function(specFile){
  #can't have duploicate entries per peptide ion/run
  stopifnot (all(specFile[,.N, by  = .(ProteinName, PeptideSequence, PrecursorCharge, Run)]$N == 1))
  # define the two things to do a full cross join on, add dummy column to allow the full cross join
  full.peptide.ions <- unique (specFile[, .(ProteinName, 
                                            PeptideSequence,
                                            PrecursorCharge,
                                            dummyIndex = 1)])
  full.runs <- unique(specFile[,.(Condition, Run, BioReplicate,dummyIndex = 1)])
  
  # do the full cross join
  full.peptide.ions.runs <- merge (full.peptide.ions, full.runs, by = "dummyIndex",allow.cartesian=TRUE)[,dummyIndex := NULL][]
  message (nrow(full.peptide.ions), " peptide ions in at least one of ", nrow(full.runs), " runs")
  # merge back to the actual intensity data
  specFile.complete <- merge (full.peptide.ions.runs, specFile, all.x=TRUE, 
                              by = intersect(colnames(full.peptide.ions.runs),
                                             colnames(specFile)))

  print ("Not adding IsotopeLabelType to table. You likely need to add it like this: mss[, IsotopeLabelType := 'L']")
  
  return (specFile.complete)
}




# A few experimental functions for different volcano like statistics ---- 
#  

lowerConfidenceInterval <- function (log2FC, SE, DF, p = 0.025){
  # work with positive numbers as much as possible, so we want right side of distribution
  # p should be > 0.5
  if (p < 0.5) p <- 1-p
  
  direction <- sign(log2FC)
  tError <- qt ( p, DF)
  
  
  
  result <- direction *(abs(log2FC) - (tError * SE))
  result <- ifelse(sign(result) != direction, 0, result)
  result
}
# Usage: 
# results[, log2FC.95ci := lowerConfidenceInterval(log2FC, SE, DF)]


pFromTNonZero <- function (log2FC, SE, DF, threshold = 0.0){
  # work in the positive direction for simplicity
  log2FC <- abs(log2FC)
  t <- (log2FC - threshold)/SE
  pt(t, DF, lower.tail = FALSE)
}

# Usage
# results[, pGt1 := pFromTNonZero(log2FC, SE, DF, 1.0)]


pFromT.S0 <- function (log2FC, SE, DF, S0=0){
  t <- log2FC/(SE + S0)
  pt(abs(t), DF, lower.tail = FALSE) * 2
}

# Usage
# results[, testP := pFromT(log2FC, SE, DF, S0 = 0.05)]

# ---- 
