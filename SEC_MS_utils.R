require (ggplot2)
require (data.table)
require (pbapply)
rotate.x.axis.text <- theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))



#' @param secLong.dt a data.table with columns sample (char), fraction (integer), protein (char), intensity (numeric) 
#' 

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
  qcSummary <- secLong.dt[!is.na(intensity) & intensity > 0, 
                   .(  numProteins = length(unique(protein)),
                       medIntensity = median(intensity)),
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



scaleByTotalIntensity <- function(secLong.dt){
  secLong.dt[, intensity_totalScaled := intensity/(sum(intensity, na.rm = TRUE)), by= .(sample, protein)]
}


scaleByMaxIntensity <- function(secLong.dt){
  secLong.dt[, intensity_maxScaled := intensity/(max(intensity, na.rm = TRUE)), by= .(sample, protein)]
}

#' @param scaleDenom total/max. What to use for the denominator for scaled intensity, total=sum(intensity) or max(intensity)

scaledIntensityMatrices <- function(secLong.dt, scaleDenom = "total"){
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
  
  return (mats)
  
}


hclustRowsMultiMatrix <- function(matrixList, check.names = TRUE, ...){
  if (check.names){
    name.mat <- sapply(matrixList, rownames)
    allSame <- all(apply(name.mat, 1, function(x)1==length(unique(x))))
    stopifnot(allSame)
  }
  hclust(dist(do.call(cbind, intMats)), ...)
    
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
    subData[, fitted := predict(l.out)]
    subData[, residual := residuals(l.out)]
    
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


windowedCosineSimilarity <- function (intensity.mat, window = 15){
  intensity.mat[is.na(intensity.mat)] <- 0.0
  start <- 1
  stop <- ncol(intensity.mat) - (window - 1)
  
  .doOneWindow <- function (start, intensity.mat, window = 15){
    subMat <- intensity.mat[, start:(start + window-1)]
    rs <- rowSums(subMat)
    rmax <- apply(subMat, 1, which.max)
    # limit to those proteins with decent coverage in this window
    goodRows <- which( apply( subMat, 1, function(x)sum(x > 0)) > window/2)
    goodRows <- intersect(goodRows, which( rs > 0.01)) # more than 1% of protein in window
    subMat <- subMat[goodRows,]
    
    cos.dt <- setDT(reshape2::melt(cosineMatrix(subMat),
                                   varnames = c("protein1", "protein2"),
                                   value.name = "cosineSim"))[]

    cos.dt <- cos.dt[cosineSim %between% c(0.9, 1.0) & protein1 != protein2] # limit to just the most promising
    
    rs.dt <- data.table(protein = names(rs), portion = rs)
    cos.dt[rs.dt, prot1Portion := i.portion, on = c(protein1 = "protein")]
    cos.dt[rs.dt, prot2Portion := i.portion, on = c(protein2 = "protein")]
    
    rmax.dt <- data.table(protein = names(rmax), max = rmax)
    cos.dt[rmax.dt, prot1Max := i.max, on = c(protein1 = "protein")]
    
    cos.dt[,start := start]
    return (cos.dt)
  }
  
  pbapply::pblapply(start:stop, .doOneWindow, intensity.mat, window) |> 
    rbindlist()
}
