library (data.table)
library (nnls)
#also requires artms (for now)...


# dataPrep ####

#a function to take the various names for raw file columns and standardize them
.fixRawFileColumnName <- function(datTab, newName = "RawFile"){
  rawFilePattern = "Raw[. ]?(F|f)ile"
  oldName = grep(rawFilePattern, names(datTab), value=TRUE)
  stopifnot(length(oldName) == 1)
  if (oldName!=newName){
    message("Converting column name ", oldName, " to ", newName)
    setnames(datTab, oldName, newName)
  }
}
#prepares for processData 
prepareDataForMSStats = function(evidenceFile, keysFile, outfile=NULL){
  
  if (is.null(outfile)){
    outfile <- "mss.input.txt"
  }
  
  #read evidence and fix a column name
  ev = fread (evidenceFile)
  .fixRawFileColumnName(ev)
  #read keys and fix a column name
  keys = fread (keysFile)
  .fixRawFileColumnName(keys)
  
  #merge
  mergeColumns = c("RawFile")
  if ('IsotopeLabelType' %in% colnames(ev)){
    mergeColumns = c(mergeColumns, "IsotopeLabelType")
  }
  # see artMS::artMSMergeEvidenceAndKeys for more checks, but this suffices for most
  evAndKeys = merge (ev, keys, by = mergeColumns, suffixes = c("Ev", "Keys"))  #
  
  #filter
  before = nrow(evAndKeys)
  #filter evidence (empty protein names, groups, contaminants, reversed decoys)
  evAndKeys <- evAndKeys[!grepl("^$|;|CON__|REV__", Proteins),                       ]
  after = nrow(evAndKeys)
  message ("Filtered ", before-after, " rows from keyed evidence , new size = ", after, "  (empty protein names, groups, contaminants, reversed decoys)")
  evAndKeys <- evAndKeys[!is.na(Intensity)]
  evAndKeys <- evAndKeys[Intensity > 0]
  after2 = nrow(evAndKeys)
  if(after2 < after)
    message ("Filtered ", after-after2, " rows with NA intensity or intensity <= 0, new size = ", after2)
  
  # in cases of  concatenated evidence files, we might still have non-unique-to-protein peptides
  nonUniquePeptides <- evAndKeys[,.(numProteins = length(unique(Proteins))), by = `Modified sequence`][numProteins >1]$`Modified sequence`
  if (length(nonUniquePeptides) > 0){
    evAndKeys <- evAndKeys[!`Modified sequence` %in% nonUniquePeptides]
    after3 <- nrow(evAndKeys)
    message ("Filtered ", after2-after3, " rows with non-specific-to-protein peptides (", length(nonUniquePeptides),"), new size ", after3)
  }

  
  #convert to msstatsFormat
  #fileOut <- 'fullAnalysis/results_noGroups/evAndKeys.mss.txt'
  #for now I use the artMS function for this, but first have to fix a column name it expects
  # setnames(evAndKeys, "Modified sequence", "Modified.sequence")
  # evAndKeys.mss <- artMS:::.artms_getMSstatsFormat(evAndKeys, fraction=FALSE,
  #                                                  output_name=outfile,
  #                                                  data_object = TRUE)
  
  # main steps above
  # 1 sum up intensity for duplicate features
  # 2 make table complete: add rows with NA intensity for missing values.
  # 3 get columns and names suitable for MSstats

  #1
  evAndKeys.noDups <- evAndKeys[,.(Intensity = sum(Intensity), count.features = .N),
                                by = .(Proteins, 
                                        `Modified sequence`,
                                        Charge,
                                        IsotopeLabelType,
                                        Condition,
                                        BioReplicate,
                                        Run)]
  countSummed <- sum(evAndKeys.noDups$count.features > 1)
  if (countSummed > 0)
    message (countSummed, " rows (peptide_ion x run) had >1 feature; using summed intensity for these")
  
  #2
  # define the two things to do a full cross join on, add dummy column to allow the full cross join
  full.peptide.ions <- unique (evAndKeys.noDups[, .(Proteins, 
                                                   `Modified sequence`,
                                                   Charge,
                                                   dummyIndex = 1)])
  full.runs <- unique(evAndKeys.noDups[,.(IsotopeLabelType,
                                          Condition,
                                          BioReplicate,
                                          Run,
                                          dummyIndex  = 1)])
  # do the full cross join
  full.peptide.ions.runs <- merge (full.peptide.ions, full.runs, by = "dummyIndex",allow.cartesian=TRUE)[,dummyIndex := NULL][]
  message (nrow(full.peptide.ions), " peptide ions in at least one of ", nrow(full.runs), " runs")
  # merge back to the actual intensity data
  evAndKeys.noDups.complete <- merge (full.peptide.ions.runs, evAndKeys.noDups, all.x=TRUE, 
                                      by = intersect(colnames(full.peptide.ions.runs),
                                                     colnames(evAndKeys.noDups)))
  
  #3 in msstats format
  evAndKeys.mss <- evAndKeys.noDups.complete[,.(ProteinName = Proteins,
                                                PeptideSequence = `Modified sequence`,
                                                PrecursorCharge = Charge,
                                                FragmentIon = NA,
                                                ProductCharge  = NA,
                                                IsotopeLabelType,
                                                Condition,
                                                BioReplicate,
                                                Run,
                                                Intensity)]
  
  
  evAndKeys.mss  
}

# Fit BGs ####
# 
# spatialReferences <- c("Cyto","Endo_C19", 
#                        "Endo_C20",  
#                        "GalT",     "LAMT",     "PM" )

# simLocations.list <- list(Cyto = c(),
#                           Endo_C19 = c("LAMT", "Endo_C20"),
#                           Endo_C20 = c("Endo_C19", "LAMT"),
#                           GalT = c(),
#                           LAMT = c("Endo_C19", "Endo_C20"),
#                           PM = c())


locationSpecificProteins <- function (gc.results, n.equal = 0,  threshold.log2FC = 1, 
                                      threshold.pvalue = 0.005, similarLocations.list = NULL,
                                      spatialReferences = NULL){
  
  
  if (is.null(spatialReferences))
    spatialReferences <- unique(unlist(strsplit(gc.results$Label, "-")))
  
  spatialReferencePattern <- paste0(gsub ("^([-^^])", "^\\1", spatialReferences), # if first character is not ^, make it a ^
                                    collapse = "|")
  
  prots.list <- lapply (spatialReferences, locationSpecificProteins.singleLocation,
                        gc.results = gc.results, 
                        n.equal = n.equal, threshold.log2FC = threshold.log2FC, threshold.pvalue = 0.005,
                        similarLocations.list = similarLocations.list)
  names(prots.list) <- spatialReferences
  prots.list
}

locationSpecificProteins.singleLocation <- function(location, gc.results, n.equal = 0, threshold.log2FC = 1, threshold.pvalue = 0.005, similarLocations = c(), similarLocations.list = NULL){
  
  posPattern = sprintf("^%s-", location)
  negPattern = sprintf("-%s$", location)
  
  posSubset <- gc.results[grep (posPattern, Label)]
  negSubset <- gc.results[grep (negPattern, Label)]
  negSubset[, Label := gsub ("^(.*)-(.*)$", "\\2-\\1", Label)]
  negSubset[, log2FC := -log2FC]
  
  locationSubset <- rbind (posSubset, negSubset)
  locationSubset[, otherGroup  := tstrsplit(Label, "-", keep=2)]
  
  if (!is.null(similarLocations.list))
    similarLocations <- c(similarLocations, similarLocations.list[[location]]) |> unique()
  
  contrastingGroups <- setdiff(unique(locationSubset$otherGroup), similarLocations)
  
  print( sprintf("For location %s, allowing matches to %s and requiring contrast to %s",
                 location, 
                 paste0(similarLocations, collapse = ","),
                 paste0(contrastingGroups, collapse = ",")
  ))
  
  locationSpecific <- locationSubset[otherGroup %in% contrastingGroups & 
                                       log2FC > threshold.log2FC & 
                                       pvalue < threshold.pvalue, 
                                     .N, 
                                     by = Protein][ N >= length(contrastingGroups)-n.equal]$Protein
  
  return (locationSpecific)
}









#@ sampleQuant.long requires columns "Protein", "GROUP_ORIGINAL", "LogIntensities" , 
#                   names match the MSstats::dataProcess()$RunlevelData
selectProteinsToFit <- function(sampleQuant.long, groupsToCountIn = c("Endo","Lyso", "PM"), 
                                minLogIntensity = 21, minPassingReplicates=3){
  
  setDT(sampleQuant.long)
  #general idea is to find proteins that have at least minLogIntensity in at least minPassingReplicates per a single background
  
  #these proteins will be duplicated if they pass thresholds in more than 1 group:
  dupProteins <- sampleQuant.long[GROUP_ORIGINAL %in% groupsToCountIn & LogIntensities > minLogIntensity ,
                                  .N,
                                  by = .(Protein, GROUP_ORIGINAL)
                                  ][N>=minPassingReplicates,]$Protein 
  unique(dupProteins)
}




# processedData is mssstats::dataProcess()$ProcessedData

selectProteinsToFit.featureCount <- function (processedData, groupsToCountIn = c("Endo", "Lyso", "PM"), 
                                              minPassingReplicates = 3, minFeatureCountQuantile = 0.75){
  # summarize featureCount per protein/run

  featureCounts <- processedData[feature_quality == "Informative" & is_outlier == FALSE & !is.na(INTENSITY),
                                 .(numFeatures = length(unique(FEATURE))),
                                 by = .(PROTEIN, GROUP_ORIGINAL, SUBJECT_ORIGINAL)]
  
  minFeatureCount <- quantile (featureCounts[GROUP_ORIGINAL %in% groupsToCountIn, numFeatures], minFeatureCountQuantile)
  message ("Requiring ", minFeatureCount, " features")
  
  dupProteins <- featureCounts[GROUP_ORIGINAL %in% groupsToCountIn & numFeatures >= minFeatureCount,
                               .N,
                               by = .(PROTEIN, GROUP_ORIGINAL)
                               ][N >= minPassingReplicates,]$PROTEIN
    
  unique(dupProteins)
}





nnls.ftest <- function(simpleModel, fullModel){
  k <- length(fullModel$x)
  q <-  k - length(simpleModel$x) # difference in terms
  N <- nrow(fullModel$residuals)
  f <- ((simpleModel$deviance - fullModel$deviance)/q)/
    (fullModel$deviance/(N-k))
  
  pf (f, q, N-k)
}


nnls.testRestrictedModels <- function (inputMatrix, outputVector, fullModel){
  pValues <- lapply (seq_along(colnames(inputMatrix)),
          FUN=function(i){
            simpleModel <- nnls(inputMatrix[,-i], outputVector)
            return (nnls.ftest(simpleModel, fullModel))
          })
  names(pValues) <- colnames(inputMatrix)
  return (pValues)
}


#@ dataSubset :  a table of protein intensities, a subtable of MSstats::dataProcess()$RunlevelData
#@ naIntensityValue : what to replace the NA value
#
# if dataSubset has column scaledIntensity it uses that, otherwise it is scaled here
# this is important if you're passing in a subset, the scaling here won't know about extreme values in the superset
calculateFitsLong = function(dataSubset, naIntensityValue = NA, 
                             backgroundRuns = c("Endo", "Lyso", "PM"),
                             experimentRuns = NULL,
                             forceNonNegative = TRUE,
                             randomizedColumns = 0,
                             removeOutlierProteins = FALSE,
                             withIntercept = FALSE){
  
  ..handleNA <- function(x){
    ifelse(is.na(x), naIntensityValue, x)
  }
  
  setDT(dataSubset)
  if (!"GROUP_ORIGINAL" %in% colnames(dataSubset)){
    dataSubset[, GROUP_ORIGINAL := GROUP]
  }
  
  if ("scaledIntensity" %in% colnames(dataSubset)){
    summaryData = dataSubset[,.(normalizedMean = mean(scaledIntensity),
                                normalizedMedian = median(scaledIntensity)),
                             by = .(Protein, GROUP_ORIGINAL)]
    
  }else{
    message ("Doing protein intensity scaling inside of calculateFitsLong...You may want to do this ahead of time by providing a column 'scaledIntensity', especially if subsetting data")
    summaryData = dataSubset[,.(meanInt = mean(2^..handleNA(LogIntensities)), 
                                medianInt = median(2^..handleNA(LogIntensities))),
                             by = .(Protein, GROUP_ORIGINAL)]
    # TODO: doing this after taking the means seems slightly problematic... 
    #normalize (aka scale) across proteins
    summaryData[,normalizedMean := meanInt/max(meanInt, na.rm=TRUE), by = .(Protein)]
    summaryData[,normalizedMedian := medianInt/max(medianInt, na.rm=TRUE), by = .(Protein)]
  }

  
  #long to wide ( asrequired by nnls):
  normalizedMeansWide = dcast (summaryData, Protein~GROUP_ORIGINAL, value.var = "normalizedMean")
  setnafill(normalizedMeansWide, fill=naIntensityValue, cols = setdiff(colnames(normalizedMeansWide),  "Protein"))
  
  
  if (randomizedColumns >0){
    randCols <- lapply (1:randomizedColumns, function(i)sample (unlist (normalizedMeansWide[,backgroundRuns, with=FALSE]), nrow(normalizedMeansWide)))
    names(randCols) <- paste ("zzRand", 1:randomizedColumns, sep="_")
    normalizedMeansWide[,names(randCols) := randCols]
    backgroundRuns <- c(backgroundRuns, names(randCols))
  }
  
  out = list()
  pValues <- list()
  nnls.outs  = list()
  allResiduals <- data.table (Protein = unique(normalizedMeansWide$Protein))    
  
  
  nmMat <- as.matrix(normalizedMeansWide, rownames="Protein")
  nmMat[is.na(nmMat)] <- naIntensityValue
  
  if(removeOutlierProteins){
    threshold <- ifelse(is.numeric(removeOutlierProteins), removeOutlierProteins,  quantile(nmMat, 0.999, na.rm=TRUE))
    outliers <- which(apply (nmMat, 1, function(x)any(x >threshold, na.rm=TRUE)))
    nmMat <- nmMat[!rownames(nmMat) %in% names(outliers),]
    
    message ("Removed ", length(outliers), "proteins that had intensity in top 0.1%; >", threshold)
    message (paste0(names(outliers), collapse="; "))
  }
  
  
  inputMatrix <- nmMat[,backgroundRuns]
  if(is.null(experimentRuns)){
    experimentRuns <- colnames(nmMat)[!colnames(nmMat) %in% backgroundRuns]
  }
  

  if (forceNonNegative){
    if (withIntercept){
      inputMatrix <- cbind (inputMatrix, Intercept = rep(1.0, nrow(inputMatrix)))
    }
    
    
    for (group in experimentRuns){
      outputVector <-  nmMat[,group]
      missing <- !complete.cases(inputMatrix) | is.na(outputVector)
      nnls.out <- nnls(inputMatrix[!missing,], outputVector[!missing])
      nnls.outs[[group]] <- nnls.out
      coef <- nnls.out$x
      names(coef) <- colnames(inputMatrix)
      out[[group]] <- split (coef, names(coef))
      # pValues are currently garbage here...
      pValues[[group]] <-nnls.testRestrictedModels(inputMatrix[!missing,], outputVector[!missing], nnls.out) 
      #nnls.outs[[group]] = nnls.out

      residuals <- as.data.table(residuals(nnls.out), keep.rownames = TRUE)
      
      setnames(residuals, old = c("rn", "V1"), new = c("Protein", paste0("residual", ".", group)))
      allResiduals <- merge (allResiduals, residuals, all.x=TRUE, by = "Protein")
    }
    
    
  }else { #don't force negative, use basic lm
    for (group in experimentRuns){
      termString <- paste(backgroundRuns, collapse="+")
      formulaString <- paste0 (group, "~", termString)
      if(!withIntercept){
        formulaString <- paste0(formulaString, "-1")
      }
      l <- lm (as.formula(formulaString), data = normalizedMeansWide)
      nnls.outs[[group]] <- l
      lcs <- coef(summary(l))
      coef <- lcs[,"Estimate"]
      #names(coef) <- colnames(inputMatrix)
      out[[group]] <- split (coef, names(coef))
      p <- lcs[,"Pr(>|t|)"]
      pValues[[group]] <- split(p, names(p)) 
      
      
      residuals <- data.table (Protein = normalizedMeansWide$Protein, residuals = normalizedMeansWide[[group]] - predict (l, normalizedMeansWide))
      setnames(residuals, old = c("residuals"), new = paste0("residual", ".", group))
      allResiduals <- merge (allResiduals, residuals, all.x=TRUE, by = "Protein")

      #pValues[[group]] <-nnls.testRestrictedModels(inputMatrix, outputVector, nnls.out) 
      #nnls.outs[[group]] = nnls.out
    }

  }

  
  return (list(coefficients = rbindlist(out, idcol="group"), pValues = rbindlist(pValues, idcol="group"), models = nnls.outs, residuals = allResiduals))
}




locationPredict.oneSample <- function(run, location.mat, 
                                      randomizedLocations = 3, 
                                      numIterations = ifelse(randomizedLocations > 0, 100, 1),
                                      na.value = 0.0){
  

  
  
  run[is.na(run)] <- na.value
  location.mat[is.na(location.mat)] <- na.value
  
  
  .oneIteration <- function(i, run, location.mat, randomizedLocations){
    if (randomizedLocations > 0){
      randCols <- sapply (1:randomizedLocations, function(i)sample (location.mat, nrow(location.mat)))
      colnames(randCols) <- paste ("zzRand", 1:randomizedLocations, sep="_")
      withRandom <- cbind(location.mat, randCols)
    } else withRandom <- location.mat
    nnls.out <- nnls::nnls(withRandom, run)
    coef <- nnls.out$x
    return(coef)
  }
  
  # all coefficients as list
  coef.list <- lapply(1:numIterations, .oneIteration, run, location.mat, randomizedLocations)
  # in  matrix form ...
  coef.mat <- do.call(cbind, coef.list)
  # ... so I can easily take a row median
  coef <- apply(coef.mat, 1, median)
  
  # the random columns will get an NA
  names(coef) <- colnames(location.mat)
  
  return (coef)
}



#' locationPredict.matrix
#' a simpler form of the above `calculateFitsLong`
#' Given a matrix of scaled protein intensities across runs calculate a location contribution to each run
#' @param run.mat scaled (usually between 0 and 1) protein intensity. Columns are runs, rows are proteins.
#' @param location.mat as `run.mat` but columns are expected protein intensity from spatial references (Cyto, Endo, etc)
#'        rows must match rows of run.mat
#' @param randomizedLocations integer number of randomized locations to add to locationMat. These are constructed by
#'                            random sampling of `location.mat`
#' @param numIterations integer, each iteration is a new sampling of randomizedLocations. This should be set to 1 when randomizedLocations == 0
#' @param na.value  value to substitute for missing values in either input matrix

locationPredict.matrix <- function(run.mat, location.mat,
                                   randomizedLocations = 3L, 
                                   numIterations = ifelse(randomizedLocations > 0, 100, 1),
                                   na.value = 0.0){
  
  stopifnot (`mismatched number of matrix rows` = nrow(run.mat) == nrow(location.mat),
             `mismatched row names` = all(rownames(run.mat) == rownames(location.mat)))

  if(max(run.mat, na.rm = TRUE) > 1 | max(location.mat, na.rm = TRUE) > 1){
    message ("WARNING: Maximum values detected above 1.0. This works best with linearly scaled intensities between 0 and 1.0")
  }  
  
    apply(run.mat, 2, locationPredict.oneSample,
        location.mat,
        randomizedLocations , 
        numIterations ,
        na.value )
}

#' Builds a set of locaiton references using Cyto and a set of vs-Cyto contrasts.
#' @param bgVsCyto expected columns Protein, bg, log2FC
#' @param protQuant expected columns Protein, LogIntensities, GROUP, SUBJECT. GROUP must contain "Cyto"

makeLocation.mat <- function(bgVsCyto, protQuant){
  stopifnot(`missing Cyto from protQuant$GROUP` = "Cyto" %in% protQuant$GROUP)
  
  referenceProteins <- unique(bgVsCyto$Protein)
  overlapProteins <- intersect (referenceProteins, protQuant[GROUP == "Cyto", unique(Protein)])
  message ( sprintf("Of %d proteins in bgVsCyto, %d found in protQuant", length(referenceProteins), length(overlapProteins)))
  
  predictedBackgrounds <- merge (protQuant[GROUP == "Cyto"], bgVsCyto, by = "Protein", allow.cartesian = TRUE)
  predictedBackgrounds[, LogIntensities := LogIntensities + log2FC ]
  predictedBackgrounds[, GROUP := bg]
  
  # rbind with cyto
  predictedBackgrounds <- rbind (predictedBackgrounds, 
                                 protQuant[GROUP == "Cyto" & Protein %in% predictedBackgrounds$Protein],
                                 use.names = TRUE, fill = TRUE)
  

  loc.mat <- as.matrix(dcast(predictedBackgrounds, Protein~GROUP, value.var = "LogIntensities", fun.aggregate = mean, na.rm = TRUE),
            rownames = "Protein")
  loc.mat <- 2^loc.mat # convert from log2 intensities
  loc.mat <- sweep (loc.mat, 1, apply (loc.mat, 1, max, na.rm = TRUE), FUN = "/")
  return (loc.mat)
}





#' @param dataSubset MSstats dataProcess()$runLevelData output with scaledIntensity (linear range 0,1) column
#'                   or other table with columns Protein, GROUP_ORIGINAL, scaledIntensity
#' @param proteins a character vector of protein IDs to include in matrix, NULL(default) for all
#' @param naIntensityValue

locationSignatureMatrix <- function(dataSubset,
                                    proteins = NULL,
                                    naIntensityValue = NA, 
                                    backgroundRuns = NULL){
  if(!is.null(proteins))
    dataSubset <- dataSubset[Protein %in% proteins]
  if(!is.null(backgroundRuns))
    dataSubset <- dataSubset[GROUP_ORIGINAL %in% backgroundRuns]

  # take the mean over all replicates
  # this works best if there are no missing values and/or no batch effects
  # TODO: noise filter when high numbers of missing values...if you have 2/3 missing, scaledIntensity = 1.0 is not trustworthy
  summaryData = dataSubset[,
                           .(normalizedMean = mean(scaledIntensity),
                              normalizedMedian = median(scaledIntensity)), # not used
                           by = .(Protein, GROUP_ORIGINAL)]
  
  #long to wide 
  normalizedMeansWide = dcast (summaryData, Protein~GROUP_ORIGINAL, value.var = "normalizedMean")
  data.table::setnafill(normalizedMeansWide, fill=naIntensityValue, cols = setdiff(colnames(normalizedMeansWide),  "Protein"))

  nmMat <- as.matrix(normalizedMeansWide, rownames="Protein")
  return (nmMat)
}




# runLevelData is from mssquant$runLevelData
# coefficiients is from calculateFitsLong()$coefficients
# bgNames should be the names of the columns in coefficients, also names of GROUP_ORIGINAL that match backgrounds
# WARNING: bgNames should not have anything that matches replicate names in them  TODO: improve this

makeBackgrounds <- function (runLevelData, coefficients, bgNames = NULL, naValue = 0.0, makeZeroesNA = TRUE){
  proteinIntensities = runLevelData
  setDT(proteinIntensities)
  
  if (is.null(bgNames)){
    bgNames = colnames(coefficients)
  }
  
  #protein intensities only for Endo|Lyso|PM, convert to wide format to create and fill in missing values with NA
  proteinIntensities[,runName := paste(GROUP_ORIGINAL, SUBJECT_ORIGINAL, sep=".")]
  groupMatchString <- paste(bgNames, collapse="|") # "Endo|Lyso|PM" for example  TODO avoid regular expression special chars
  piWide = dcast(proteinIntensities[grep(groupMatchString,GROUP_ORIGINAL),], Protein~runName, value.var = "LogIntensities")
  
  #replace the NA'S and linearize
  #naValue = quantile(piWide[,2:10], 0.001, na.rm=TRUE)

  for (col in grep("Protein",colnames(piWide), invert=TRUE, value=TRUE)){
    piWide[,(col) := 2^piWide[[col]]]
    piWide[is.na(piWide[[col]]), (col) := naValue]
  }
  
  bgRepMats = list()
  bgCompositeTables = list()
  
  #TODO generalize!!!
  #process by replicate
  replicates <- unique(proteinIntensities$SUBJECT_ORIGINAL)
  for (i in replicates){
    #matrix'ize for math purposes
    bgRepMats[[i]] = as.matrix(piWide[,grep(as.character(i), colnames(piWide)), with=FALSE], 
                               rownames.value = as.character(piWide$Protein))
    
    #matching column order matters, so use the order in the data matrix to get the right fits column order
    matchString <- sprintf ("(%s)\\..", groupMatchString)
    cols =gsub(matchString, "\\1", colnames(bgRepMats[[i]]))
    
    # matrix multiplication: fits X backgrounds; 
    # as.table %>% as.data.table conversion leads to long format
    # protein names come from rownames of bgRepMats[[i]]
    bgCompositeTables[[i]] = as.data.table(as.table((as.matrix(coefficients[, cols, with=FALSE], rownames.value = coefficients$group) 
                              %*% t(bgRepMats[[i]]) )))
    setnames(bgCompositeTables[[i]], c("composite", "Protein", "intensity"))
    #append the _compBG. label to end of group so we know its a control
    bgCompositeTables[[i]][,composite:=paste(composite, "_compBG.", i, sep="")]
  }
  #reassemble all replicates into one list
  bgComposites = do.call(rbind, bgCompositeTables)
  if (makeZeroesNA) bgComposites[intensity==0, intensity:=NA]
  
  return (bgComposites)
}

buildMSStatsDataWithCompositeBackgrounds <- function(bgComposites, mssquant){
  bgComposites <- copy(bgComposites)
  maxRun = max(as.integer(as.character(mssquant$RunlevelData$RUN))) # %>% as.character() %>% as.integer() %>% max()
  bgComposites[,RUN:= as.factor(as.integer(as.factor(composite))+maxRun) ]
  bgComposites[,LogIntensities := log2(intensity)]
  bgComposites[,NumMeasuredFeature := NA] # I hope this doesn't matter
  bgComposites[,MissingPercentage := NA] # I hope this doesn't matter
  bgComposites[,more50missing := NA] # I hope this doesn't matter
  bgComposites[,NumImputedFeature := NA] # I hope this doesn't matter
  bgComposites[,originalRUN := RUN]
  maxGroup = max(as.integer(as.character(mssquant$RunlevelData$GROUP)))# %>% as.character() %>% as.integer() %>% max()
  bgComposites[,GROUP:=NA] # a placeholder to get column ordering correct.
  bgComposites[,GROUP_ORIGINAL:=gsub("\\.[1-9]", "", composite)]
  bgComposites[,GROUP:= as.factor(as.integer(as.factor(GROUP_ORIGINAL)) + maxGroup)]# %>% as.factor()]
  
  #TODO  the bit below works only because SUBJECT_ORIGINAL is 1,2 or 3.  If they were different strings or nonsequential, I'll probably have to do a better mapping between SUBJECT_ORIGINAL and SUBJECT
  bgComposites[,SUBJECT_ORIGINAL:=as.factor(gsub(".*\\.([1-9])", "\\1", composite))]# %>% as.factor()]
  stopifnot(all(levels(unique(bgComposites)) == as.character(1:3)))  # fail if assumption in TODO above is wrong
  bgComposites[,SUBJECT_NESTED:= as.factor(paste (GROUP,SUBJECT_ORIGINAL, sep="."))]# %>% as.factor()]
  bgComposites[,SUBJECT := SUBJECT_ORIGINAL]
  
  bgComposites[,composite:=NULL]
  bgComposites[,intensity:=NULL]
  
  mssquantComp <- copy (mssquant)
  if (!setequal(colnames(mssquantComp$RunlevelData), colnames(bgComposites))){
    message ("unexpected columns in mssquant$RunlevelData: ", paste(setdiff (colnames(mssquantComp$RunlevelData), colnames(bgComposites)), collapse = ", "))
  }
  mssquantComp$RunlevelData <- rbind(mssquantComp$RunlevelData, bgComposites, fill=TRUE)
  
  
  dummyRows <- mssquantComp$ProcessedData[1:length(unique(bgComposites$GROUP_ORIGINAL)),]
  setDT(dummyRows)
  dummyRows[,(colnames(dummyRows)) := rep(NA, nrow(dummyRows))]
  dummyRows[,GROUP_ORIGINAL :=unique(bgComposites$GROUP_ORIGINAL)]
  
  if (!setequal(colnames(mssquantComp$ProcessedData), colnames(dummyRows))){
    message ("unexpected columns in mssquant$ProcessedData: ", paste(setdiff (colnames(mssquantComp$ProcessedData), colnames(dummyRows)), collapse = ", "))
    
  }
  
  mssquantComp$ProcessedData <- rbind (mssquantComp$ProcessedData, dummyRows, fill=TRUE)
  
  
  # fix GROUP, GROUP_ORIGINAL so order of levels is as expected by MSstats
  #
  groupNamesOrdered <- sort (levels(mssquantComp$RunlevelData$GROUP_ORIGINAL))
  #this fixes the order of the levels without changing which row gets which level assigned to it.
  mssquantComp$RunlevelData[,GROUP_ORIGINAL := factor(mssquantComp$RunlevelData$GROUP_ORIGINAL,
                                                      levels = groupNamesOrdered)]
  #now assign GROUP according to the integer value of GROUP_ORIGINAL
  mssquantComp$RunlevelData[,GROUP := as.factor(as.integer(GROUP_ORIGINAL))]
  #repeat for processedData
  mssquantComp$ProcessedData[,GROUP_ORIGINAL := factor(mssquantComp$ProcessedData$GROUP_ORIGINAL,
                                                       levels = groupNamesOrdered)]
  mssquantComp$ProcessedData[,GROUP := as.factor(as.integer(GROUP_ORIGINAL))]
  
  #inspect the mappings and stop if they don't meet expectations:
  checkGroupMappings <- function(data){
    groupMappings <- data[,.(groupOriginalOrder = unique(as.integer(GROUP_ORIGINAL)), group = unique(GROUP), groupOrder = unique(as.integer(GROUP))), by=GROUP_ORIGINAL]
    setorder(groupMappings, GROUP_ORIGINAL)
    stopifnot (
      # does order match GROUP and GROUP_ORIGINAL
      groupMappings$groupOriginalOrder == groupMappings$groupOrder,
      # are there any missing values
      groupMappings$groupOrder == 1:(max(groupMappings$groupOrder)), 
      # does character group match the order "1" and 1
      as.integer(as.character(groupMappings$group)) == groupMappings$groupOrder) 
  }
  checkGroupMappings(mssquantComp$ProcessedData)
  checkGroupMappings(mssquantComp$RunlevelData)
  
  
  # drop down to data.frames just in case it causes issues
  mssquantComp$ProcessedData <- as.data.frame(mssquantComp$ProcessedData)
  mssquantComp$RunlevelData  <- as.data.frame(mssquantComp$RunlevelData)
  return (mssquantComp)
}




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



groupQuantificationUsingEmmeans <- function (groupComparisonOutput){
  #this fails in some cases "singular" ...
  allMeans <- lapply(groupComparisonOutput$fittedmodel,
                     FUN = function(l)if(is.null(l)){return(NULL)} else {return (as.data.table(emmeans::emmeans(l, specs = "GROUP")))}
                     )

  rbindlist(allMeans)  
}




groupQuantificationUsingLMER <- function(groupComparisonOutput){
  predData <- pbapply::pblapply (1:length(groupComparisonOutput$fittedmodel), function(i)predictedData(gc.backgrounds, i))
  names(predData) <- levels(groupComparisonOutput$ModelQC$PROTEIN)
  predictedLong <-  rbindlist(predData, idcol = "PROTEIN")
  actualData <- groupComparisonOutput$ModelQC[,.(PROTEIN, GROUP, GROUP_ORIGINAL, SUBJECT, SUBJECT_ORIGINAL, ABUNDANCE )]
  
  datMerge <- merge (actualData, predictedLong, by = c("PROTEIN", "GROUP", "SUBJECT"), all=TRUE)
  
}

  
predictedData <- function(compResults, modelID){
    groups <- unique(compResults$ModelQC[PROTEIN == levels(PROTEIN)[modelID],.(GROUP)])
    subjects <- unique(compResults$ModelQC[PROTEIN == levels(PROTEIN)[modelID],.(SUBJECT)])
    newData <- merge (data.frame(GROUP = groups), data.frame(SUBJECT=subjects))
    setDT(newData)
    setorder(newData, GROUP, SUBJECT)
    modelClass <- class(compResults$fittedmodel[[modelID]]) 
    if ("lmerMod" %in% modelClass){
      model <- lmerTest::as_lmerModLmerTest(compResults$fittedmodel[[modelID]])
    } else if("lm" %in% modelClass){
      model <- compResults$fittedmodel[[modelID]]
    } else{
      if (!is.null(modelClass) & modelClass != "NULL"){
        message ("Unexpected model class, ", modelClass)
      }
      return (NULL)
    }
    return (newData[,prediction := predict (model, newdata = newData)] )
  }



detrendNMFFits <- function (nmf.res, actual, columnInfo, spatialRefPattern = spatialReferencePattern){
  if(is.null(nmf.res))
    return (NULL)
  nmf.predicted <- basis(nmf.res) %*% coef(nmf.res)
  nmf.predicted <- nmf.predicted[, !grepl(spatialReferencePattern, colnames(nmf.predicted))]
  actual <- actual[, !grepl(spatialReferencePattern, colnames(actual))]
  allFits.nmf.detrended <- detrendedFits(actual, nmf.predicted, columnInfo)
  return (allFits.nmf.detrended)
}




#' Works on matrices of actual and predicted. Detrends and fits using fitPoly.MultiplePowers, available in bp_utils/
#' @param actual a matrix of actual values
#' @param predicted a matrix of predicted values.  The fit will be done on log2(actual/predicted). Any infinite values will be replaced with NA, thus zeros will be treated as missing.
#' @param columnInfo a data.table of information to attach to the columns.  For now requires column, time, receptor, SUBJECT_ORIGINAL columns
  
detrendedFits <- function(actual, predicted, columnInfo){
    # check for duplicates
    duplicatedNames <- rownames(actual) [duplicated(rownames(actual))]
    if (length(duplicatedNames) > 0){
      message ("matrices contain duplicate entries. These will be removed: ", paste0(duplicatedNames, collapse = ","))
      actual <- actual[! (rownames(actual) %in% duplicatedNames),]
    }
    # allow predicted too be 1, then we do no detrending, otherwise make sure rows and columns match
    if (!(length(predicted) == 1 & predicted[1]==1)){
      stopifnot (ncol(actual) == ncol(predicted)) # column names allowed to differ.
      predicted <- predicted[rownames(actual),]
    }
    
    
    detrended <- log2(actual/predicted)
    if (any (is.infinite(detrended))){
      message ("Infinite values detected after detrendeing and log transform. Passed in 'actual' matrix likely contains values = 0; these are set to NA for lm fitting")
      detrended[is.infinite(detrended)] <- NA
    }
    
    long <- melt (as.data.table(detrended, keep.rownames = TRUE), id.vars = "rn")
    setnames(long, c("gene", "column", "detrendLog2Int"))
    long <- long[!is.na(detrendLog2Int)]
    long <- merge (long, columnInfo, by = "column")
    long[,orderedTime := as.integer(as.factor(time))]
    
    # do the fitting, per receptor
    groups <- unique(long$receptor)
    allFits <- list()
    
    numCores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(numCores)
    for(g in groups){
      print (g)
      subTables <- split (long[ receptor == g], by = "gene")
      #geneFits <- lapply (subTables, fitPoly.MultiplePowers, polyColumn = "orderedTime", otherTerms = c("SUBJECT_ORIGINAL"), yColumn = "detrendLog2Int")
      geneFits <- parLapply (cl, subTables, fitPoly.MultiplePowers, polyColumn = "orderedTime", otherTerms = c("SUBJECT_ORIGINAL"), yColumn = "detrendLog2Int")
      fitScores <- rbindlist(geneFits, idcol = "Protein", fill = TRUE)
      allFits[[g]] <- fitScores
    }
    parallel::stopCluster(cl)
    
    allFitsDT <- rbindlist(allFits, idcol = "receptor", use.names=TRUE, fill = TRUE)
    allFitsDT[,adj.f.pvalue := p.adjust(best.pF.polyColumn, method = "BH"), by = "receptor"]
    return (allFitsDT)
  }  
  


