require (ggplot2)
require (data.table)
require (pbapply)
rotate.x.axis.text <- theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))



#' @param secLong.dt a data.table with columns sample (char), fraction (integer), protein (char), intensity (numeric) 
#' 

# QC ----
medPolishFractionIntensity <- function(secLong.dt){
  int.mat <- dcast (secLong.dt, protein~sample+fraction, value.var = "intensity") |>
    as.matrix(rownames = "protein")
  int.mat[int.mat == 0] <- NA
  int.mat <- log2(int.mat)
  mp.out <- medpolish(int.mat, na.rm = TRUE)
  
  colInfo <- secLong.dt[, .(columnID = paste0(sample, "_", fraction)), by= .(sample, fraction)]
  
  data.table(columnID = names(mp.out$col),
             medPolishIntensity = mp.out$overall + mp.out$col)[
               colInfo,c("sample", "fraction") := .(i.sample, i.fraction),
               on =  "columnID", ]
  
}


qcSummaryTable <- function(secLong.dt){
  qcSummary <- secLong.dt[, 
                   .(  numProteins = length(unique(protein[!is.na(intensity) & intensity > 0])),
                       medIntensity = median(intensity[!is.na(intensity) & intensity > 0])),
                   by = .(sample, fraction)]
  
  qcSummary[medPolishFractionIntensity(secLong.dt),
            medPolishIntensity := i.medPolishIntensity ,
            on = c("sample", "fraction")]
  
  
  if (all(c("treatment", "replicate")%in% colnames(secLong.dt))){
    replicateInfo <- secLong.dt[, .(dummy = 1), by = .(sample, treatment, replicate)][, dummy := NULL]
    qcSummary <- merge(qcSummary, replicateInfo, by = "sample")
  }
  qcSummary
}


qcPlotProteinCount <- function(qcSummary){
  
  p <- ggplot (qcSummary, aes(x = fraction, y = numProteins))  + 
    geom_point() +
    geom_line(alpha = 0.5) + 
    rotate.x.axis.text
  
  if (all(c("treatment", "replicate")%in% colnames(qcSummary))){
    p <- p + facet_grid(treatment~replicate)                        
  }else{
    p <- p + facet_wrap(~sample)
  }
  p
}


qcPlotMedianInt <- function(qcSummary){
  p <- ggplot (qcSummary, aes(x = fraction, y = log2(medIntensity)))  + 
    geom_line(alpha = 0.5) + 
    geom_point(aes(color = "median")) +
    geom_point(aes(y = medPolishIntensity, color = "median polish"), alpha = 0.5) +
    scale_color_manual(values = c(median = "black", `median polish` = "red")) +
    rotate.x.axis.text
  
  if (all(c("treatment", "replicate")%in% colnames(qcSummary))){
    p <- p + facet_grid(treatment~replicate)                        
  }else{
    p <- p + facet_wrap(~sample)
  }
  p
}



qcPlotsSecMS <- function(secLong.dt){
  qcSummary <- qcSummaryTable(secLong.dt)
  
  list ( proteinCount = qcPlotProteinCount(qcSummary),
         medianIntensity = qcPlotMedianIntensity(qcSummary)
  )
}


# scaling ----


scaleByTotalIntensity <- function(secLong.dt){
  secLong.dt[, intensity_totalScaled := intensity/(sum(intensity, na.rm = TRUE)), by= .(sample, protein)]
}


scaleByMaxIntensity <- function(secLong.dt){
  secLong.dt[, intensity_maxScaled := intensity/(max(intensity, na.rm = TRUE)), by= .(sample, protein)]
}


# matrices ----

#' @param scaleDenom total/max. What to use for the denominator for scaled intensity, total=sum(intensity) or max(intensity)

scaledIntensityMatrices <- function(secLong.dt, scaleDenom = "total", reorder = TRUE){
  if(scaleDenom == "total" & !"intensity_totalScaled" %in% colnames(secLong.dt))
    scaleByTotalIntensity(secLong.dt)
  if(scaleDenom == "max" & !"intensity_maxScaled" %in% colnames(secLong.dt))
    scaleByMaxIntensity(secLong.dt)
  
  allProteins <- unique(secLong.dt$protein)
  
  v_var <- sprintf("intensity_%sScaled", scaleDenom)
  
  .oneMatrix <- function(sub.dt){
    mat <- dcast(sub.dt, protein~fraction, value.var = v_var)[allProteins,, on = "protein"] |> as.matrix(rownames = "protein")
    mat[is.na(mat)] <- 0.0
    mat[order(rownames(mat)),]
  }
  
  mats <- lapply(split(secLong.dt, secLong.dt$sample), .oneMatrix )
  
  if (reorder){
    proteinMaxFractions <- secLong.dt[, .(maxFraction = fraction[which.max(intensity)]), by= .(sample, protein)][, .(meanMaxFraction = mean(maxFraction, na.rm = TRUE)), by = protein]
    rowOrdering <- proteinMaxFractions[rownames(mats[[1]]), order(meanMaxFraction), on = "protein"]
    
    mats <- lapply(mats, function(x)x[rowOrdering,])
    
  }
  
  return (mats)
  
}


hclustRowsMultiMatrix <- function(matrixList, check.names = TRUE, ...){
  if (check.names){
    name.mat <- sapply(matrixList, rownames)
    allSame <- all(apply(name.mat, 1, function(x)1==length(unique(x))))
    stopifnot(allSame)
  }
  hclust(dist(do.call(cbind, matrixList)), ...)
    
}

#' @param mats a list of matrices or a single matrix, scaled or not doesn't matter.  column-normalized is probably best
#' @description
#' Reorders matrices rows based on two statistics: the fraction with the highest intensity, and secondarily the center of mass
#' When passed a list of matrices, these stats are averaged across all matrices before ordering 

reorderMatricesByMaxFraction <- function(mats){
  if(!"list" %in% class(mats)){
    if ("matrix" %in% class(mats)){
      mats <- list(mats)
    }else{
      stop("Expect a list of matrices or a single matrix")
    }
  }
  
  rowMaxes <- sapply(mats, function(x)apply(x,1, which.max), simplify = "array")
  avgRowMaxes <- apply(rowMaxes, 1, mean, na.rm = TRUE)
  
  centerOfMass <- function(x) sum(x * 1:length(x))/sum(x)
  rowCenters <- sapply(mats, function(x)apply(x,1, centerOfMass), simplify = "array")
  avgRowCenters <- apply(rowCenters, 1, mean, na.rm = TRUE)
  
  rowOrder <- order (avgRowCenters, avgRowMaxes)
  
  mats <- lapply(mats, function(x)x[rowOrder, ])
  return (mats)
  
}
# Heatmaps ----

intensityHeatmaps <- function(intMats, intensityName = "Scaled Intensity"){
  colorFun <- circlize::colorRamp2(breaks = (0:50)/200, colors = viridis::magma(51,direction = -1))
  samples <- names(intMats)
  
  sample <- samples[1]
  hml <- Heatmap (intMats[[sample]] ,
                  name = intensityName,
                  cluster_rows = FALSE,
                  row_dend_reorder = FALSE,
                  cluster_columns = FALSE,
                  col = colorFun ,
                  show_row_names = FALSE, 
                  column_names_gp = gpar(fontsize = 5),
                  column_labels = ifelse(as.integer(colnames(intMats[[sample]])) %% 5 == 0, colnames(intMats[[sample]]), ""),
                  column_title = sample,
                  
                  #first only:
                  show_heatmap_legend = TRUE, 
                  row_title = sprintf ("%d Proteins", nrow(intMats[[sample]]))
  )
  
  # first one gets row title, also gets a legend
  
  if (length(intMats) > 1){
    for (sample in samples[2:length(samples)]){
      hml <- hml + Heatmap (intMats[[sample]] ,
                            name = intensityName,
                            cluster_rows = FALSE,
                            row_dend_reorder = FALSE,
                            cluster_columns = FALSE,
                            col = colorFun ,
                            show_row_names = FALSE, 
                            column_names_gp = gpar(fontsize = 5),
                            column_labels = ifelse(as.integer(colnames(intMats[[sample]])) %% 5 == 0, colnames(intMats[[sample]]), ""),
                            column_title = sample,
                            # 
                            show_heatmap_legend = FALSE)      
    }
    
  }
  return (hml)
}


# Normalization and Outlier detection by fitting local cubics ----

#' @param sampleTerm additive/interaction whether to fit a different curve per sample or allow only
#'                   an offset per sample.  This is per window, so different offsets will be used, 
#'                   and some amount of curve-unique-to-sample persists weven with "additive"

fitLocalCubics<- function (qcSummary,window =15,extend  = 3, sampleTerm = "additive"){
  additiveModel <-  as.formula ("medPolishIntensity~poly(fraction,3) + sample")
  interactionModel <-     as.formula ("medPolishIntensity~poly(fraction,3) * sample")
  model <- c(interaction = interactionModel, additive =additiveModel)[[sampleTerm]]
  
  .doOneLocalCubic <- function( startRange, fullData, rangeWidth){
    
    subData <- fullData[inrange(fraction, startRange, startRange + rangeWidth)]
    l.out <- MASS::rlm( model , data = subData  )
    subData[!is.na(medPolishIntensity), fitted := predict(l.out)]
    subData[!is.na(medPolishIntensity), residual := residuals(l.out)]
    
    subData[, .(fit = sprintf("localFit.%02d.%02d", rangeWidth, startRange), sample,  fraction, fitted, residual)]
  }
  
  start <- min(qcSummary$fraction - extend)
  stop <- max(qcSummary$fraction - (window-1) + extend)
  
  allFits <-  lapply(start:stop, .doOneLocalCubic, fullData = qcSummary , rangeWidth = window) |> rbindlist()
  allFits
}


#' @param threshold absolute fold change (in linear space) to define an outlier. A value with abs(median(residuals)) > log2(threshold) will be labeled an outlier 

labelOutliers <- function(qcSummary, localCubicFits, threshold = 1.5){
  qcSummary[allFits[, .(medianResidual = median(residual)), by = .(sample , fraction)],
            medianResidual := i.medianResidual,
            on = c("sample",  "fraction")
  ]
  qcSummary[, isOutlier := abs(medianResidual) > log2(threshold)]
}

#' @param secLong.dt SEC-MS data in long format.
#' A data.table with columns sample (char), fraction (integer), protein (char), intensity (numeric) 
#' 
#' @description
#' Given residuals computed by fitLocalCubics and labelOutliers, this will adjust intensity by median(residual).
#' Original intensity will be stored in new column originalIntensity
#' If any intensity_*scaled columns exist, they will be recalculated
#' 
normalizeByResiduals <- function(secLong.dt, qcSummary){
  stopifnot ("medianResidual" %in% colnames(qcSummary))
  stopifnot(!"originalIntensity" %in% colnames(secLong.dt))
  
  secLong.dt[, originalIntensity := intensity]
  secLong.dt[, intensity := NA_real_]
  secLong.dt[qcSummary, intensity := originalIntensity / (2^i.medianResidual) , on = c("sample", "fraction")]
  
  if("intensity_totalScaled" %in% colnames(secLong.dt)){
    scaleByTotalIntensity(secLong.dt)
  }
  if("intensity_maxScaled" %in% colnames(secLong.dt)){
    scaleByMaxIntensity(secLong.dt)
  }
}

#' @description
#' Mostly for developing/debugging, espeically when normalizeByResiduals fails after updating intensity
UNnormalize <- function(secLong.dt){
  secLong.dt[, intensity := originalIntensity]
  secLong.dt[, originalIntensity := NULL]
  secLong.dt[, itensity_totalScaled := NULL]
  secLong.dt[, itensity_maxScaled := NULL]
}


#' @description 
#' A line and scatterplot showing per-fraciton values (by median polish), outliers labeled, and fitted curves
plotNormAndOutlierFits <- function(qcSummary, allFits){
  p <- ggplot (allFits, aes(x = fraction)) + 
    geom_line(aes(y = fitted, group = fit), show.legend = FALSE, alpha = 0.2) + 
    geom_point(data = qcSummary, aes(y = medPolishIntensity, color = isOutlier), show.legend = FALSE) +
    theme_bw()
  
  if (all(c("treatment", "replicate")  %in% colnames(allFits)) &
      all(c("treatment", "replicate") %in% colnames(qcSummary))){
    p <- p +  facet_grid(treatment~replicate)
  }else{
    p <- p + facet_wrap(~sample)
  }
  p
}


# Cosine Similarity ----

#' 
cosineMatrix <- function(mat){
  # see : https://stats.stackexchange.com/questions/31565/compute-a-cosine-dissimilarity-matrix-in-r
  # Thanks Brad!
  sim <- mat / sqrt(rowSums(mat * mat))
  sim <- sim %*% t(sim)  
  return (sim)
} 


# windowedCosineSimilarity <- function (intensity.mat, window = 15){
#   intensity.mat[is.na(intensity.mat)] <- 0.0
#   start <- 1
#   stop <- ncol(intensity.mat) - (window - 1)
#   
#   .doOneWindow <- function (start, intensity.mat, window = 15){
#     subMat <- intensity.mat[, start:(start + window-1)]
#     rs <- rowSums(subMat)
#     rmax <- apply(subMat, 1, which.max)
#     # limit to those proteins with decent coverage in this window
#     goodRows <- which( apply( subMat, 1, function(x)sum(x > 0)) > window/2)
#     goodRows <- intersect(goodRows, which( rs > 0.01)) # more than 1% of protein in window
#     subMat <- subMat[goodRows,]
#     
#     cos.dt <- setDT(reshape2::melt(cosineMatrix(subMat),
#                                    varnames = c("protein1", "protein2"),
#                                    value.name = "cosineSim"))[]
# 
#     cos.dt <- cos.dt[cosineSim %between% c(0.9, 1.0) & protein1 != protein2] # limit to just the most promising
#     
#     rs.dt <- data.table(protein = names(rs), portion = rs)
#     cos.dt[rs.dt, prot1Portion := i.portion, on = c(protein1 = "protein")]
#     cos.dt[rs.dt, prot2Portion := i.portion, on = c(protein2 = "protein")]
#     
#     rmax.dt <- data.table(protein = names(rmax), max = rmax)
#     cos.dt[rmax.dt, prot1Max := i.max, on = c(protein1 = "protein")]
#     
#     cos.dt[,start := start]
#     return (cos.dt)
#   }
#   
#   pbapply::pblapply(start:stop, .doOneWindow, intensity.mat, window) |> 
#     rbindlist()
# }

#' @param intensity.mat a normalized, but probably not smoothed, SEC matrix
#' @param goodPeaks.mat a boolean matrix with TRUE identifying a "good peak"
#' @description
#' This version filters protein profiles to compare based on presence of a good peak within center +/- peakRadius
#' It then calculates all-by-all cosine similarities based on profile in center +/- outerRadius 
#' Only values above 

windowedCosineSimilarity <- function (intensity.mat, goodPeaks.mat, outerRadius = 6, peakRadius = 2){
  intensity.mat[is.na(intensity.mat)] <- 0.0
  start <- 1 + outerRadius
  stop <- ncol(intensity.mat) - outerRadius
  
  
  .doOneWindow <- function (start, intensity.mat, goodPeaks.mat, outerRadius, peakRadius){
    rowsWithCentralPeaks <- apply(goodPeaks.mat[, (start - peakRadius):(start + peakRadius)],1,any)
    subMat <- intensity.mat[rowsWithCentralPeaks, (start-outerRadius):(start + outerRadius)]
    
    rs <- rowSums(subMat)
    
    cos.dt <- setDT(reshape2::melt(cosineMatrix(subMat),
                                   varnames = c("protein1", "protein2"),
                                   value.name = "cosineSim")
    )[]
    
    cos.dt[, obsP := rank(-cosineSim)/.N]
    cos.dt[, Z := (cosineSim -mean(cosineSim))/sd(cosineSim) ]
    
    cos.dt <- cos.dt[cosineSim > 0.9 & protein1 != protein2] # limit to just the most promising
    
    rs.dt <- data.table(protein = names(rs), portion = rs)
    cos.dt[rs.dt, prot1Portion := i.portion, on = c(protein1 = "protein")]
    cos.dt[rs.dt, prot2Portion := i.portion, on = c(protein2 = "protein")]
    
    cos.dt[,start := start]
    
    
    
    
    return (cos.dt)
  }
  
  pbapply::pblapply(start:stop, .doOneWindow, intensity.mat, goodPeaks.mat, outerRadius, peakRadius) |> 
    rbindlist()
}

# Correlation ----

#' @param intensity.mat a normalized, but probably not smoothed, SEC matrix
#' @param goodPeaks.mat a boolean matrix with TRUE identifying a "good peak"
#' @param goldStandardInteractome a data table with protein1, protein2 columns listing known interactions.
#'        If this is not NULL, it will be used to create a per-fraction confidence score.
#' @description
#' This version filters protein profiles to compare based on presence of a good peak within center +/- peakRadius
#' It then calculates all-by-all correlation based on profile in center +/- outerRadius 
#' Only values above 

windowedCorrelation <- function (intensity.mat, goodPeaks.mat, outerRadius = 6, peakRadius = 2, goldStandardInteractome = NULL, priorOdds = 0.01){
  intensity.mat[is.na(intensity.mat)] <- 0.0
  start <- 1 + outerRadius
  stop <- ncol(intensity.mat) - outerRadius
  
  
  .doOneWindow <- function (start, intensity.mat, goodPeaks.mat, outerRadius, peakRadius){
    rowsWithCentralPeaks <- apply(goodPeaks.mat[, (start - peakRadius):(start + peakRadius)],1,any)
    subMat <- intensity.mat[rowsWithCentralPeaks, (start-outerRadius):(start + outerRadius)]
    
    rs <- rowSums(subMat)
    
    cor.dt <- setDT(reshape2::melt(cor(t(subMat), use = "pairwise"),
                                   varnames = c("protein1", "protein2"),
                                   value.name = "pearsonR")
    )[]
    
    
    # translate cor to a log distance from 1.0
    cor.dt[, corScore := -log10(1 - pearsonR)]
    corScoreRange <- range (cor.dt[is.finite (corScore)]$corScore, na.rm = TRUE)
    cor.dt[corScore == Inf, corScore := corScoreRange[2] * 1.01]

    if (!is.null(goldStandardInteractome)){
      # label the known, for interactions in  both directions
      cor.dt[, knownInteractor := FALSE]
      cor.dt[goldStandardInteractome, knownInteractor := TRUE, on = c("protein1", "protein2")]
      cor.dt[goldStandardInteractome, knownInteractor := TRUE, on = c(protein2 = "protein1", protein1 = "protein2")]
      
      interactorDensityFun <- density(cor.dt[knownInteractor == TRUE , corScore],  bw = 0.2) |> stats::approxfun()
      bgDensityFun <-    density(cor.dt[knownInteractor == FALSE, corScore],  bw = 0.2) |> stats::approxfun()
      
      cor.dt <- cor.dt[pearsonR > 0.9 & protein1 != protein2] # limit to just the most promising
      cor.dt[, c("interactorLL", "bgLL") := .(interactorDensityFun(corScore), bgDensityFun(corScore))]
      cor.dt[, llRatio := interactorLL/bgLL ]
      cor.dt[, postOdds := llRatio * priorOdds ]
      cor.dt[, postProb := postOdds/(postOdds + 1) ]
      cor.dt[, start := start]

      xRange <- c(-1, 16.5)
      stats.dt <- data.table(start = start, x = seq (from = xRange[1], to = xRange[2], length.out = 512))
      stats.dt[, c("interactorLL", "bgLL") := .(interactorDensityFun(x), bgDensityFun(x))]
      stats.dt[, llRatio := interactorLL/bgLL ]
      stats.dt[, postOdds := llRatio * priorOdds ]
      stats.dt[, postProb := postOdds/(postOdds + 1) ]
      
    }else{
      cor.dt <- cor.dt[pearsonR > 0.9 & protein1 != protein2] # limit to just the most promising
      stats.dt <- data.table()
    }
    
    
    rs.dt <- data.table(protein = names(rs), portion = rs)
    cor.dt[rs.dt, prot1Portion := i.portion, on = c(protein1 = "protein")]
    cor.dt[rs.dt, prot2Portion := i.portion, on = c(protein2 = "protein")]
    

    return ( list(cor = cor.dt, stats = stats.dt))
  }
  
  out.ls <- pbapply::pblapply(start:stop, .doOneWindow, intensity.mat, goodPeaks.mat, outerRadius, peakRadius)
  
  cor.dt <- lapply(out.ls, function(x)x$cor) |> rbindlist()
  stats.dt <- lapply(out.ls, function(x)x$stats) |> rbindlist()
  
  return (list (cor = cor.dt, stats = stats.dt))
  
}



# Matrix Smoothing ----

#' @description
#' Smooths the rows of a matrix using gaussian kernel smoothing, row-wise
#' @param a matrix to smooth, like a SEC MS matrix. Rows are proteins, columns are equally spaced fractions
#' @param ... parameters to pass on, most relevantly `bandwidth`. Default bandwidth is 2.68, or SD = 1, same as PCprophet
#' @value a matrix of same dimensions and as input
smoothRowsOfMatrix <- function(mat, ...){
  result <- t(apply(mat, 1, smoothGuassianKernel, ...))
  colnames(result) <- colnames(mat)
  return (result)
}


#' @description
#' Smooth an equally spaced vector of numbers, such as a single SEC profile, a single row of a matrix.
#' 
#' Bandwidth effects the amount of smoothing. Wider bandwidths average over more points.
#' According to ksmooth help, quartiles of the gaussian weights should be +/- bandwidth/4
#' So 2.68 converts to 0.67, +/- which corresponds to 1st and 3rd quartiles when SD = 1
#' Thus 2.68 is equivalent to the smoothing in PCprophet.
#' To verify compare `dnorm(c(-4:4))` # SD = 1, mean = 0
#' With `ksmooth(1:9, c(0,0,0,0,1,0,0,0,0), n.points = 9, kernel = "normal", bandwidth = 2.68)$y`
#' 
smoothGuassianKernel <- function(yValues, bandwidth = 2.68){
  xValues <- 1:length(yValues)
  
  
  ksmooth(xValues, yValues, n.points = length(yValues),
          kernel = "normal",
          bandwidth = bandwidth)$y
}




# Peak detection ----

#' @description
#' Not updated. Don't use
#' 

findAllPeaksInLongDT <- function(secLong.dt, intensityColumn = "intensitySmoothed"){
  message ("This is an old/slow version, probably not updated")
  
  allFractions <- sort(unique(secLong.dt$fraction))
  stopifnot ("integer" %in% class (allFractions))
  
  .replaceNA <- function(x, replace= 0){
    x[is.na(x)] <- replace
    return (x)
  }
  
  .findPeaksConsistentReturn <- function(...){
    value = pracma::findpeaks(...)
    if(is.null(value)) value <- matrix(c(NA_real_, NA_integer_, NA_integer_, NA_integer_), nrow = 1)
    return (value)
  }
  
  # a rather complex data.table call here.
  # basically, we process by looping over all  sample/protein.
  # we take the subtable (.SD) and expand to allFractions as needed.
  # replace the NAs with zeros
  # use pracma::findpeaks, a simple peak detectiong algorithm
  # convert output to data.table
  # set meaningful column names
  
  # data.table then recombines all .SD together with sample/protein
  
  
  secLong.dt[,.SD[data.table(fraction = allFractions),,on = "fraction"][[intensityColumn]] |> # expand .SD by allFractions, then get just the intensityColumn as specified
               .replaceNA(0.0) |> 
               #pracma::findpeaks(nups =2, ndowns = 2) |>
               .findPeaksConsistentReturn(nups = 2, ndowns = 2) |> 
               as.data.table() |>
               setnames(new =  c("peakHeight", "peakLocation", "peakStart", "peakEnd")),
             by = .(sample, protein)] # loop over all sample/protein
  
  
}

#' faster than the LongDT form

findAllPeaksInSingleMatrix <- function(intMat){
  if (any(is.na(intMat)))
    intMat[is.na(intMat)] <- 0.0
  message ("peak locating...")
  allPeaks <- apply(intMat, 1, pracma::findpeaks, nups = 2, ndowns = 2 )
  allPeaks <- rbindlist(lapply(allPeaks, as.data.table), idcol = "protein")
  setnames(allPeaks, old = c("V1", "V2", "V3", "V4"), new =  c("peakHeight", "peakLocation", "peakStart", "peakEnd"))
  
  renumberFractions <- FALSE
  # CHECK FRACTION NAMES for unexpected things.
  # I favor warnings here, though these are probalby errors most of the time
  # use while as an if that I can break out of
  while (!is.null(colnames(intMat))){
    #check for integer colnames
    fractionNumbers <- as.integer(colnames(intMat))
    if (any(is.na(fractionNumbers))){
      warning("Column/fraction names are not integers; will use indexes as fraction numbers")
      break
    }
    # check for integers in order
    if (!(all(fractionNumbers == sort(fractionNumbers)))){
      warning ("Column/fraction names as integers are out of order. Ignoring the fraction names and using fraction index(order) as fraction numbers")
      break
    }
    # check for missing fractions by name
    if (length(fractionNumbers) != 1 + max(fractionNumbers) - min(fractionNumbers) ){
      stop("Missing fraction(s) ", paste0(setdiff(min(fractionNumbers):max(fractionNumbers), fractionNumbers), collapse = ","))
    }
    
    renumberFractions <- TRUE

    break # only one time through the while-as-if block
  }

  if(renumberFractions){
    message ("Using fraction numbers given in intensity matrix ", min(fractionNumbers), "-", max(fractionNumbers)," instead of fraction index 1-", length(fractionNumbers) )
    allPeaks[, peakLocation := fractionNumbers[peakLocation]]
    allPeaks[, peakStart := fractionNumbers[peakStart]]
    allPeaks[, peakEnd := fractionNumbers[peakEnd]]
  }  
  
  # center of mass for peaks
  .cofmCenterN <- function (protein, peakCenter, peakWindow =2 ){
    fractions <- peakCenter + (-peakWindow:peakWindow)
    fractionNames <- as.character(fractions) # as.character to support matrix columns that don't start at 1
    sum(intMat[protein, fractionNames ] * fractions/sum(intMat[protein, fractionNames ]))
  }
  .cofmFullPeak <- function (protein, peakCenter, start, stop){
    # in case there is a "lopsided" peak, keep it centered 
    radius <- min(abs(peakCenter - c(start, stop)))
    .cofmCenterN(protein, peakCenter, peakWindow = radius)
  }
  
  message ("center of mass calc...")
  allPeaks[, cofmN := .cofmCenterN(protein, peakLocation), by = .I]
  allPeaks[, cofmFull := .cofmFullPeak(protein, peakLocation, peakStart, peakEnd), by = .I]
  message ("...done")
  


  return(allPeaks)
}

#' @description
#' Takes a single sec matrix, smoothed ahead of time, and locates good peaks
#' Good peaks are defined as:
#' 1) peak shaped up-up-center-down-down
#' 2) center > 0.01 intensity (1% of total protein in one fraction, after smoothing)
#' 3) TBD... a width or sharp increase....

goodPeaksMatrixFromIntensityMatrix <- function (intMat, minPeakHeight = 0.01){
  stopifnot(!is.null(rownames(intMat)),
            !is.null(colnames(intMat)))
  
  allPeaks <- findAllPeaksInSingleMatrix(intMat)
  allPeaks[, goodPeak := (peakHeight > minPeakHeight)]
  
  #expand allPeaks table to include allProteins and allFractions
  allProteins <- rownames(intMat)
  allFractions <- as.integer(colnames(intMat))
  dummyRows <- data.table(protein = allProteins,
                          peakLocation = allFractions,
                          goodPeak = FALSE) |>
    suppressWarnings() # to ignore warning about recycling with remainder shorter vector
  
  #allPeaks <-  rbind(dummyRows, allPeaks, use.names = TRUE, fill = TRUE)
  
  # allPeaks <- allPeaks[
  #   data.table(protein = allProteins)[, .(fraction = allFractions), by = protein],
  #   , on = c("protein", peakLocation = "fraction")]

  peakMat <- dcast(rbind(dummyRows, allPeaks, use.names = TRUE, fill = TRUE), 
                   protein~peakLocation, value.var = "goodPeak",
                   fun.aggregate = any
  ) |> as.matrix(rownames = "protein")
  

  peakMat[is.na(peakMat)] <- FALSE
  
  return(list(matrix = peakMat[rownames(intMat), colnames(intMat)],
              table = allPeaks))
}

# Prob Density Estimation

# based on ggplot's density functions which deal with boundaries nicely.

.reflect_density <- 
  function (dens, bounds, from, to,...) 
  {
    if (all(is.infinite(bounds))) {
      return(dens)
    }
    f_dens <- stats::approxfun(x = dens$x, y = dens$y, method = "linear", 
                               yleft = 0, yright = 0)
    left <- max(from, bounds[1])
    right <- min(to, bounds[2])
    out_x <- seq(from = left, to = right, length.out = length(dens$x))
    left_reflection <- f_dens(bounds[1] + (bounds[1] - out_x))
    right_reflection <- f_dens(bounds[2] + (bounds[2] - out_x))
    out_y <- f_dens(out_x) + left_reflection + right_reflection
    list(x = out_x, y = out_y)
  }

boundedDensity <- function(x, bounds = c(-Inf, Inf), bw  = 0.01, cut = 3,...){
  
  from = min(x) - bw * cut
  to = max(x) + bw *cut
  
  
  
  dens <- stats::density(x, bw = bw, cut = cut, ...)
  dens <- .reflect_density(dens = dens, bounds = bounds, from = from, to = to, 
                           ...)
  return(dens)
  
}

# Peak alignment between samples ----

#' @description
#' A utility function to find closest peaks between two peak tables.
#' 
matchTwoPeakTables <- function (peaks.dt, i.peaks.dt, matchColumn  = "peakLocation"){
  # define a column to join on to preserve the original columns
  peaks.dt[, peakJoin := peaks.dt[[matchColumn]] ]
  i.peaks.dt[, peakJoin := i.peaks.dt[[matchColumn]] ]
  
  # also define a peakID column. this helps to make sure we have a one-to-one mapping in both directions
  peaks.dt[, peakID := paste0(protein, ".", 1:nrow(.SD)), by = protein]
  i.peaks.dt[, peakID := paste0(protein, ".",1:nrow(.SD)), by = protein]
  
  setkey(peaks.dt, protein, peakJoin)
  setkey(i.peaks.dt, protein, peakJoin)
  
  # do the joins in both directions
  combined <- peaks.dt[i.peaks.dt, roll = "nearest"][!is.na(peakLocation)]
  revCombined <- i.peaks.dt[peaks.dt, roll = "nearest"][!is.na(peakLocation)]
  
  # limit to the matching pairs closest A-B and B-A. Avoids problems of A-BC. (A in one profile, BC in the other) 
  # This will have B-A and C-A in reverse, but only B-A will be kept. C will be ignored. 
  peakPairs <- fintersect (combined[, .(peakID, i.peakID)],
                           revCombined[, .(peakID = i.peakID, i.peakID = peakID)])
  combined <- combined[peakPairs, on = c("peakID", "i.peakID")]
  #revCombined <- revCombined[peakPairs, on = c(peakID = "i.peakID", i.peakID = "peakID")]
  
  combined[, peakJoin := NULL] # don't need it any longer
  combined[, deltaPeak := i.peakLocation - peakLocation] 
  #combined[, i.deltaPeak := peakLocation - i.peakLocation]
  #revCombined[, peakJoin := NULL] # don't need it any longer
  #revCombined[, deltaPeak := i.peakLocation - peakLocation] 
  
  firstPeakFraction <- min(combined[, c(peakLocation, i.peakLocation)], na.rm = TRUE)
  lastPeakFraction <- max(combined[, c(peakLocation, i.peakLocation)], na.rm = TRUE)
  
  # provide balance at the ends by removing deltaPeaks that are not possible in reverse direction
  # min peak can only possibly shift positive so there's a bias to large positive
  combined[deltaPeak > peakLocation-firstPeakFraction, deltaPeak := NA]
  #revCombined[deltaPeak > peakLocation-firstPeakFraction, deltaPeak := NA]
  #combined[i.deltaPeak > i.peakLocation-firstPeakFraction, i.deltaPeak := NA]
  # and max peak can only possibly shift negative so there's a bias to large negative
  combined[deltaPeak < peakLocation-lastPeakFraction, deltaPeak := NA]
  #revCombined[deltaPeak < peakLocation-lastPeakFraction, deltaPeak := NA]
  #combined[i.deltaPeak < i.peakLocation-lastPeakFraction, i.deltaPeak := NA]
  
  return (combined)
  #return(list (forward = combined[], reverse = revCombined[]))
}



#' @param otherPeaks data.table of peaks (output of `findAllPeaks`) that will adjust to the standard
#' @param otherPeaks data.table of peaks, as above, that will serve as the standard. These will not be adjusted
#' @param minPeaksPerFraction required number of good peaks to detect in a fraction to include in fraction alignment
#'                            Fractions with low counts of peaks may skew the alignment
#' @param fitPortion Per fraction in otherPeaks, this is the portion of points (around the median) to include in loess fitting.
#'                   The default, 0.75, means bottom 12.5% and  top 12.5% of values are ignored as they are likely full of outliers. 
standardizeOnePeakTableToStandard <- function (otherPeaks, standardPeaks, sec.dt, sampleName, minPeaksPerFraction = 50, firstPeak = 1, lastPeak = 72, doPlots = TRUE, fitPortion = 0.75){
  message ("Matching peaks between samples...")
  matchedPeaks <- matchTwoPeakTables (otherPeaks, standardPeaks)
  
  goodFractions <- matchedPeaks[goodPeak == TRUE & i.goodPeak == TRUE, .N, by = peakLocation][N > minPeaksPerFraction, peakLocation]
  
  subData <- matchedPeaks[goodPeak == TRUE & i.goodPeak == TRUE][
    peakLocation %in% goodFractions][
      ,.SD[i.cofmN %between% quantile(i.cofmN, c( (1-fitPortion)/2, fitPortion + (1-fitPortion)/2))], by = peakLocation] # middle 75% only, or fitPortion
      
      #abs(i.cofmN-cofmN) < maxPeakShift]
  
  if(doPlots){
    p <- ggplot (subData, 
                 aes(x = cofmN, y = i.cofmN)) + 
      geom_point(alpha = 0.2, shape = ".") + 
      coord_fixed() + 
      geom_abline(slope = 1) +
      geom_density_2d() +
      geom_smooth(method = "loess", span = 0.25, se = FALSE, color = "red", lwd = 0.5) +
      ggtitle(label = sampleName) + 
      theme_bw()
    print (p)
    
  }

  # what is the best i.cofmN(standard) given the cofmN(other)
  loess.model <- loess(formula = i.cofmN~cofmN, 
                       data = subData,
                       span = 0.25)
  
  message ("Adding/updating column cofmN.standardized with standardized peak cofmN")
  otherPeaks[, cofmN.standardized := predict (loess.model, cofmN)]
  
  
  # interpolate the lower tail
  # a function to do simple two-point linear interpolation
  .linearInterpolateY <- function (newX, 
                                   xRange = c(1,5.660036),
                                   yRange = c(1, 5.228412)){
    (newX-xRange[1])/(xRange[2] - xRange[1]) * (yRange[2] - yRange[1]) + yRange[1]  
  }
  # define the tail boundaries in both spaces:
  standardRange <- otherPeaks[!is.na(cofmN.standardized), range (cofmN.standardized, na.rm = TRUE)]
  otherRange <- otherPeaks[!is.na(cofmN.standardized), range (cofmN, na.rm = TRUE)]
  
  otherPeaks[cofmN < min(otherRange), cofmN.standardized := .linearInterpolateY(cofmN, xRange = c(firstPeak, min(otherRange)), y = c(firstPeak, min(standardRange))) ]
  otherPeaks[cofmN > max(otherRange), cofmN.standardized := .linearInterpolateY(cofmN, xRange = c(max(otherRange), lastPeak), y = c(max(standardRange),  lastPeak)) ]

  message ("Updating sec.dt with standardFraction for sample ", sampleName)
  sec.dt[ sample == sampleName, standardFraction := predict (loess.model, fraction)]
  sec.dt[fraction < min(otherRange) & sample == sampleName, standardFraction := .linearInterpolateY(fraction, xRange = c(firstPeak, min(otherRange)), y = c(firstPeak, min(standardRange)))]
  sec.dt[fraction > max(otherRange) & sample == sampleName, standardFraction := .linearInterpolateY(fraction, xRange = c(max(otherRange), lastPeak), y = c(max(standardRange),  lastPeak)) ]
    
  invisible(otherPeaks)
}




#' @param allPeakTables a list of peak tables. One table per SEC sample
#' @param sec.dt a long data.table of sec data. Required columns are
#'  * `sample` must match values in `names(allPeakTables)`
#'  * `protein` must match values in `allPeakTables$protein`
#'  * `fraction` must match values in `allPeakTables$peakLocation`
#'  sec.dt will be updated with a `standardFraction` column
#' @param standardIdx Either an integer or character name identifying which table of allPeakTables to use
#'                    as the standard that all other tables will be adjusted to.
#' @param ... arguments passed to `standardizeOnePeakTableToStandard` that does the pairwise adjusting
#' 
standardizeAllPeakTablesToStandard <- function (allPeakTables, sec.dt,  standardIdx = 1, ...){
  if("character" %in% class(standardIdx)){
    standardIdx <- which(names(allPeakTables) == standardIdx)
  }
  
  standardPeakTable <- allPeakTables[[standardIdx]]
  
  firstPeak <- sec.dt[, min(fraction)]
  lastPeak  <- sec.dt[, max(fraction)]
  
  for (i in 1:length(allPeakTables)){
    sampleName <-   names(allPeakTables)[i]
    if (i == standardIdx){
      # set standardFraction to actual fraction
      sec.dt[sample == sampleName, standardFraction := fraction]
      next
    }
    message ("Calculating standard peak cofmN for ", names(allPeakTables)[standardIdx])
    standardizeOnePeakTableToStandard (allPeakTables[[i]], standardPeakTable, sec.dt, sampleName,firstPeak = firstPeak, lastPeak = lastPeak, ...)
  }
  
  standardPeakTable[, cofmN.standardized := cofmN] # standard doesn't change
  
  .cleanPeakJoiningColumns <- function(peakTable){
    peakTable[, peakJoin := NULL]
    peakTable[, peakID := NULL]
  }
  for (peak.dt in allPeakTables){
    cleanPeakJoiningColumns(peak.dt)
  }
  invisible(allPeakTables)
}

# peak differential statistics ----

clusterPeaks <- function(peakCenters, maxDistance = 2){
  if (length(peakCenters) <= 1){
    return (rep(1, length(peakCenters)))
  }
  dist(peakCenters, method = "manhattan") |>
    hclust(method = "complete") |>
    cutree(h = maxDistance)
}



diffAnovaOnePeak <- function (onePeakSec.dt, doPlotFunction = NULL){
  lmout <- lm( log2(intensity_totalScaled)~ poly(standardFraction, 4)*treatment, onePeakSec.dt)
  
  #onePeakSec.dt[, fitted := 2^predict(lmout)]
  
  if (!is.null(doPlotFunction)){
    fittedCurve <- data.table(treatment = unique(onePeakSec.dt$treatment))[, .(standardFraction = seq(from = min(onePeakSec.dt$standardFraction), to = max(onePeakSec.dt$standardFraction), length = 100)), by = treatment]
    fittedCurve[, intensity_totalScaled := 2^predict (lmout, newdata = fittedCurve)] 
    
    p <- ggplot(onePeakSec.dt, aes(x =standardFraction, y = intensity_totalScaled,color = treatment)) +
      geom_point(aes(shape = as.factor(replicate))) + 
      geom_line(aes(color = treatment, group = sample)) +
      geom_line(data= fittedCurve)
    
    if ("function" %in% class(doPlotFunction)){doPlotFunction(p)}else{print(p)}
    
  }
  return (anova(lmout))
}



# do the subset per protein ahead of time here

anovaPeaksInOneProtein <- function(sec.dt, peakClusters, radius = 5, ... ){

  anova.tables <- lapply(peakClusters$proteinPeakCluster,
                         function(clusterID){
                           peakClusters[clusterID, subsetAndDiffAnova.oneProtein(sec.dt, center, radius = radius, ...)] 
                           }) |> suppressWarnings()

  allAnova <- rbindlist(anova.tables, idcol = "proteinPeakCluster", fill = TRUE, use.names = TRUE)
  return (allAnova)
}

subsetAndDiffAnova.oneProtein <- function(sec.dt, peakCenter, radius = 5, ...){
  stopifnot (length(unique(sec.dt$protein)) == 1) # only intended for one protein
  onePeak <-  sec.dt[ standardFraction %between% (c(-radius, radius) + peakCenter)]
  result <- tryCatch( as.data.table(diffAnovaOnePeak(onePeak, ...), keep.rownames = TRUE),
                      error = function(cond){
                        return (data.table(error = conditionMessage(cond)))
                      }
  )
  return (result)
}



# this is slow on full table, because it has to subset 40K times, or however many peakClusters there are
# 
# subsetAndDiffAnova <- function(sec.dt, proteinOI, peakCenter, radius = 5, ...){
#   onePeak <-  sec.dt[protein == proteinOI  & standardFraction %between% (c(-radius, radius) + peakCenter)]
#   
#   
#   result <- tryCatch( as.data.table(diffAnovaOnePeak(onePeak, ...), keep.rownames = TRUE),
#                       error = function(cond){
#                         return (data.table(error = conditionMessage(cond)))
#                       }
#   )
#   return (result)
# }


