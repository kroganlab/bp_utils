require (ggplot2)
require (data.table)
require (pbapply)
rotate.x.axis.text <- theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))



#' @param secLong.dt a data.table with columns sample (char), fraction (integer), protein (char), intensity (numeric) 
#' 

summarizeSECTable <- function(secLong.dt){
  secLong.dt[, .(numFractions = length(unique(fraction)), minFraction = min(fraction), maxFraction = max(fraction), nrow = .N), by = .(treatment, replicate)]
}



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
    geom_line(aes(y = medPolishIntensity, color = "median polish"), alpha = 0.5) + 
    
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


qcFractionRunsPerProtein <- function(matrices= NULL, sec.dt = NULL){
  if (is.null(matrices))
    matrices <- scaledIntensityMatrices(sec.dt)

  .oneRow <- function(intensity){
    r <- rle(intensity>0)
    if (!any(r$values == TRUE))
      return(0)
    max(r$lengths[r$values == TRUE])
  }  
  .oneMatrix <- function(mat){
    apply(mat, 1, .oneRow)
  }
  
  runs <- lapply(matrices, .oneMatrix)
  # runs is a list of protein-named runs
  
  do.call(cbind, runs) |> 
    as.data.table(keep.rownames = TRUE) |>
    setnames(old = "rn", new = "protein") |>
    melt(id.vars = "protein", value.name = "longestRun", variable.name = "sample")
}

qcFractionsPerProtein <- function(matrices = NULL, sec.dt = NULL){
  if (is.null(matrices))
    matrices <- scaledIntensityMatrices(sec.dt)
  
  .oneRow <- function(x)sum(x>0, na.rm = TRUE)
  .oneMatrix <- function(mat){
    apply(mat, 1, .oneRow)
  }
  fractionCounts <- lapply(matrices, .oneMatrix)
  do.call(cbind, fractionCounts) |> 
    as.data.table(keep.rownames = TRUE) |>
    setnames(old = "rn", new = "protein") |>
    melt(id.vars = "protein", value.name = "countFractions", variable.name = "sample")
}



## correlations/overlaps ----


#' simplest correlation, sample by sample
#' concatenates all fraction vectors for a sample into a nFraction X nProtein length vector
#' correlates these full sample vectors

qcFullSampleCorrelation <- function(sec.dt){
  long.sample.mat <- dcast(sec.dt, sprintf ("%02d.%s", fraction, protein)~sample, value.var = "intensity_totalScaled") |> 
    as.matrix(rownames = 1)
  
  long.sample.mat[is.na(long.sample.mat)] <- 0.0
  sample.cor <- cor(long.sample.mat)
  return (sample.cor)
}


#' similar to `qcFullSampleCorrelation`, but builds offset matrices for each matrix
#' Each offset matrix is padded with zero-columns at left and right, up to radius, and leftpad + rightpad = 10
#' @param startFraction limits the fractions used to a middle-subset of all fractions
#' @param endFraction limits the fractions used to a middle-subset of all fractions
#' @returns A data.table with correlation values in long format
crossCorrelationFullSample <- function(sec.dt, startFraction = 15, endFraction = 60, radius = 5){
  # Helper functions to build zero-padded matrices
  ..padOneMat <- function(mat, pre, post, fill = 0.0){
    prePad = NULL
    postPad = NULL
    if(pre > 0)
      prePad <- matrix(rep(fill, pre * nrow(mat)), nrow = nrow(mat), dimnames = list(rownames(mat), -(pre-1):0))
    if (post > 0)
      postPad <- matrix(rep(fill, post * nrow(mat)), nrow = nrow(mat), dimnames = list(rownames(mat), (ncol(mat)+1):(ncol(mat)+post)))
    paddedMat <- cbind(prePad, mat, postPad)
  }
  .padListMatrices <- function (matList, offset, radius){
    stopifnot(abs(offset) <= radius)
    pre = radius - offset
    post = radius + offset
    paddedMats <- lapply (matList, ..padOneMat, pre = pre, post = post)
    names(paddedMats) <- sprintf("%s_offset%d", names(matList), offset)
    return (paddedMats)
  }
  # /Helper functions
  
  all.mats <- scaledIntensityMatrices(sec.dt[fraction %inrange% c(startFraction, endFraction)])

  paddedMats <- lapply (-radius:radius, function(x).padListMatrices(all.mats, offset = x, radius = radius))
  # a list of lists, but all should be named....
  paddedVecs <- lapply( unlist(paddedMats, recursive =FALSE), as.vector)
  matOfOffsetMats <- do.call(cbind, paddedVecs)
  
  
  matOfOffsetMats[is.na(matOfOffsetMats)] <- 0.0
  allCor <- cor(matOfOffsetMats)
  
  cor.long <- melt(as.data.table(allCor,keep.rownames = TRUE),id.vars = "rn", variable.name = "cn", value.name = "pearsonR")
  cor.long <- cor.long[grepl("offset0", rn)]
  cor.long[, refSample := gsub ("_offset0", "", rn)]
  cor.long[, c("otherSample", "offset") := tstrsplit(cn, "_offset")]
  cor.long[, offset := as.integer(offset)]
  
  offset.mat <- dcast(cor.long[, .SD[which.max(pearsonR)], by= .(refSample, otherSample)],
                      refSample~otherSample, value.var = "offset") |> as.matrix(rownames = 1) 
  
  offset.hc.order <- hclust(dist(offset.mat))$order
  cor.long[, refSample := factor(refSample, levels = rownames(offset.mat)[offset.hc.order])]
  cor.long[, otherSample := factor(otherSample, levels = rownames(offset.mat)[offset.hc.order])]
  
  
}


#' correlation of all fractions to all fractions
#' @param returnDataTable returns data.table if TRUE, or correlation matrix if FALSE
qcFractionByFractionCorrelation <- function(sec.dt = NULL, matrices = NULL, returnDataTable = TRUE, method = "pearson"){
  if (is.null(matrices))
    matrices <- scaledIntensityMatrices(sec.dt)

  # rename columns to include sample name
  res <- purrr::map2(names(matrices), matrices, function(name, mat){colnames(mat) <- sprintf("%s.fraction%02d", name, as.integer(colnames(mat))); mat})
  fullMat <- do.call(cbind, res)
  fullMat[is.na(fullMat)] <- 0.0
  suppressWarnings(fullCor <- cor(fullMat, use = "everything", method = method))
  
  
  if (returnDataTable){
    fullCor.long <- melt(as.data.table(fullCor, keep.rownames = TRUE), id.vars = "rn", value.name = method, variable.name = "cn")
    fullCor.long[, c("ref.sample", "ref.fraction") := tstrsplit(rn, "\\.")]
    fullCor.long[, c("other.sample", "other.fraction") := tstrsplit(cn, "\\.")]
    fullCor.long[, ref.fraction := as.integer(gsub("fraction", "", ref.fraction))]
    fullCor.long[, other.fraction := as.integer(gsub("fraction", "", other.fraction))]
    fullCor.long[, ref.treatment := tstrsplit(ref.sample, "_")[[1]]]
    fullCor.long[, other.treatment := tstrsplit(other.sample, "_")[[1]]]
    fullCor.long[, ref.replicate := tstrsplit(ref.sample, "_")[[2]]]
    fullCor.long[, other.replicate := tstrsplit(other.sample, "_")[[2]]]    
    
    return (fullCor.long)
  }
  return (fullCor) # matrix

}


#' @param fullCor.long output of qcFractionByFractionCorrelation(..., returnDataTable = TRUE)
#' @param allowFractionOffsets 
#' @param splitByReplicate
#' @param colorFacetStrips boolean indicating if an attempt should be made to color the facet headers to match
#'                         the treatment colors. This is hack-ish, may break, and requires that you use
#'                         grid.draw to draw the output.

qcFractionCorrelationLinePlot <- function(fullCor.long,
                                          colorFacetStrips = TRUE,
                                          allowFractionOffsets = FALSE,
                                          splitByReplicate = FALSE,
                                          facetTextColor = "white", # to look good on below colors
                                          treatmentColors = c("#E41A1C", "#377EB8", "#4DAF4A", 
                                                              "#984EA3", "#FF7F00", # "#FFFFFF",
                                                              "#A65628", "#F781BF", "#999999")){
  
  
  if (allowFractionOffsets){
    fullCor.long.fractionMatch <- fullCor.long[!is.na(pearson), .SD[which.max(pearson)], by = .(ref.sample, ref.fraction, other.sample)]
  }else{
    fullCor.long.fractionMatch <- fullCor.long[ref.fraction == other.fraction]
  }

  p <- ggplot(fullCor.long.fractionMatch[ref.sample!= other.sample], aes(x = ref.fraction, y = pearson, color = other.treatment)) + 
    geom_line( aes(lty = other.treatment != ref.treatment,  alpha = other.treatment != ref.treatment, group = interaction(ref.sample, other.sample)))  +
    theme_bw() +
    scale_linetype(name = 'contrast', label = c("within treatment", "between treatments")) +
    scale_color_manual(values = treatmentColors) +
    scale_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 1.0)) +
    guides(alpha = "none") +
    theme(strip.text = element_text(color = facetTextColor))
  
  if (splitByReplicate){
    p <- p + facet_grid(ref.treatment~ref.replicate)
  }  else{
    p <- p + facet_grid(ref.treatment~.)
  }
  
  if (colorFacetStrips){
    # see https://github.com/tidyverse/ggplot2/issues/2096#issuecomment-389825118
    g <- ggplot_gtable(ggplot_build(p))
    strip_r <- which(grepl('strip-r', g$layout$name))
    fills <- treatmentColors
    k <- 1
    for (i in strip_r) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    message ("Returning a grid graphical object, g. Print using grid.draw(g). Or BackupAsPDF(function()grid.draw(g))")
    return(g)
    #grid.draw(g)
  }
  return (p)
}


qcProteinOverlapByFraction <- function (sec.dt = NULL, matrices= NULL, returnTable = TRUE){
  if (is.null(matrices))
    matrices <- scaledIntensityMatrices(sec.dt)
  
  # rename columns to include sample name
  res <- purrr::map2(names(matrices), matrices, function(name, mat){colnames(mat) <- sprintf("%s.fraction%02d", name, as.integer(colnames(mat))); mat})
  fullMat <- do.call(cbind, res)
  fullMat[is.na(fullMat)] <- 0.0
  
  # we use dist, method = binary, which does 1-jaccard between rows of a matrix
  # ceiling converts numeric to smallest integer >= x, so 0 when 0, and 1 when >1
  allDist <- dist(t(ceiling(fullMat)), method = "binary")
  allJaccard <- 1- as.matrix(allDist)
  if (!returnTable) return (allJaccard)
  
  allJaccard <- melt(as.data.table(allJaccard, keep.rownames = TRUE), id.vars = "rn", value.name = 'jaccard', variable.name = "cn")
  allJaccard[, c("ref.sample", "ref.fraction") := tstrsplit(rn, "\\.")]
  allJaccard[, c("other.sample", "other.fraction") := tstrsplit(cn, "\\.")]
  allJaccard[, ref.fraction := as.integer(gsub("fraction", "", ref.fraction))]
  allJaccard[, other.fraction := as.integer(gsub("fraction", "", other.fraction))]
  allJaccard[, ref.treatment := tstrsplit(ref.sample, "_")[[1]]]
  allJaccard[, other.treatment := tstrsplit(other.sample, "_")[[1]]]
  allJaccard[, ref.replicate := tstrsplit(ref.sample, "_")[[2]]]
  allJaccard[, other.replicate := tstrsplit(other.sample, "_")[[2]]]    
  
  return (allJaccard)
}

qcProteinOverlapByFractionPlot <- function (overlap.dt,
                                            colorFacetStrips = TRUE,
                                            facetTextColor = "white", # to look good on below colors
                                            treatmentColors = c("#E41A1C", "#377EB8", "#4DAF4A", 
                                                                "#984EA3", "#FF7F00", # "#FFFFFF",
                                                                "#A65628", "#F781BF", "#999999")){
  p <- ggplot(proteinOverlap[ref.sample != other.sample & ref.fraction == other.fraction],
              aes(x = ref.fraction, y = jaccard, color = other.treatment)) +
    geom_line( aes(lty = other.treatment != ref.treatment,  alpha = other.treatment != ref.treatment, group = interaction(ref.sample, other.sample)))  +
    theme_bw() +
    scale_linetype(name = 'contrast', label = c("within treatment", "between treatments")) +
    scale_color_manual(values = treatmentColors) +
    scale_alpha_manual(values = c(`TRUE` = 0.5, `FALSE` = 1.0)) +
    guides(alpha = "none") +
    theme(strip.text = element_text(color = facetTextColor)) +
    facet_grid(ref.treatment~.)
  
  if (colorFacetStrips){
    # see https://github.com/tidyverse/ggplot2/issues/2096#issuecomment-389825118
    g <- ggplot_gtable(ggplot_build(p))
    strip_r <- which(grepl('strip-r', g$layout$name))
    fills <- treatmentColors
    k <- 1
    for (i in strip_r) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    message ("Returning a grid graphical object, g. Print using grid.draw(g). Or BackupAsPDF(function()grid.draw(g))")
    return(g)
    #grid.draw(g)
  }
  return (p)
}



qcProteinCorrelationBetweenSamples <- function(sec.dt = NULL, matrices = NULL){
  message ("not advised; little reason to do all-by-all protein correlation as QC.")
  if (is.null(matrices))
    matrices <- scaledIntensityMatrices(sec.dt)
  
  # rename rows to include sample name
  res <- purrr::map2(names(matrices), matrices, function(name, mat){rownames(mat) <- sprintf("%s.protein.%s", name, rownames(mat)); mat})
  fullMat <- do.call(rbind, res)
  fullMat[is.na(fullMat)] <- 0.0
  
  fullCor <- cor(t(fullMat), use = "pairwise")
  return (fullCor)
}


qcSameProteinCorrelationBetweenSamples <- function(sec.dt = NULL, matrices = NULL){
  if (is.null(matrices))
    matrices <- scaledIntensityMatrices(sec.dt)
  
  
  pairSamples <- data.table(ref.sample = names(matrices))[
    , .(other.sample = names(matrices)), by = ref.sample][
      other.sample!= ref.sample]
  
  pairSamples[, pair := paste0(ref.sample, "_", other.sample)]
  
  # take the transpose to correlate by protein (proteins in columns)
  tmatrices <- lapply(matrices, t)
  # convert to data.table or data.frame, so mapply can walk across either. 
  tmatrices <- lapply(tmatrices, as.data.table)
  
  .onePair <- function(ref, other)suppressWarnings( mapply(cor, tmatrices[[ref]], tmatrices[[other]]))
  
  
  allCors <- purrr::map2(pairSamples$ref.sample, pairSamples$other.sample, .onePair, .progress = "Protein auto-correlation between samples...")
  names(allCors) <- pairSamples$pair
  
  cor.dt <- do.call(cbind, allCors) |>
    as.data.table(keep.rownames = TRUE) |>
    melt(id.vars = "rn", value.name = "pearsonR", variable.name = "pair")
  
  cor.dt[pairSamples, c("ref.sample", "other.sample") := .(i.ref.sample, i.other.sample), on = "pair"]
  cor.dt[, other.treatment := tstrsplit(other.sample, "_")[[1]]]
  cor.dt[, ref.treatment := tstrsplit(ref.sample, "_")[[1]]]
  cor.dt[, other.replicate := tstrsplit(other.sample, "_")[[2]]]
  cor.dt[, ref.replicate := tstrsplit(ref.sample, "_")[[2]]]
  
  return (cor.dt[])
}



qcSameProteinCorrelationBetweenSamplesPlot <- function (cor.dt,
                                                        colorFacetStrips = TRUE,
                                                        facetTextColor = "white", # to look good on below colors
                                                        treatmentColors = c("#E41A1C", "#377EB8", "#4DAF4A", 
                                                                            "#984EA3", "#FF7F00", # "#FFFFFF",
                                                                            "#A65628", "#F781BF", "#999999")){
  
  
  p <- ggplot (cor.dt, aes(x = other.sample, y = pearsonR, fill = other.treatment)) +
    geom_hline(yintercept = 0.0, lty = 'dashed', color = "grey") + 
    geom_violin(aes(lty = other.treatment != ref.treatment) ) + 
    facet_grid(ref.treatment~ref.replicate) +
    theme_bw() +
    rotate.x.axis.text +
    scale_linetype(name = 'contrast', label = c("within treatment", "between treatments")) +
    scale_fill_manual(values = treatmentColors)
  
  
  if (colorFacetStrips){
    # see https://github.com/tidyverse/ggplot2/issues/2096#issuecomment-389825118
    g <- ggplot_gtable(ggplot_build(p))
    strip_r <- which(grepl('strip-r', g$layout$name))
    fills <- treatmentColors
    k <- 1
    for (i in strip_r) {
      j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
    message ("Returning a grid graphical object, g. Print using grid.draw(g). Or BackupAsPDF(function()grid.draw(g))")
    return(g)
    #grid.draw(g)
  }
  return (p)
}



crossCorrelateRatioOneVector <- function(x, maxShift = 5, numerator = 1, denominator = 3:5){
  if (all(x == 0)){
    return(NA)
  }
  offsetMatrix <- sapply(1:(maxShift + 1), function(start)x[start:(length(x) -1 - maxShift + start)])
  cor.mat <- cor(offsetMatrix[,1], offsetMatrix[, 2:maxShift])
  mean(cor.mat[numerator])/mean(cor.mat[denominator])
  
}


## cluster and GO enrich heatmap views ----

clusterOneMat <- function(mat, corThreshold = 0.9, minClusterSize = 5){
  mat[is.na(mat)] <- 0.0
  emptyProteins <- names(which(rowSums(mat) == 0))
  # before correlation, add consistent noise per fraction so completely-empty proteins have perfect correlation and cluster together.
  mat <- sweep(mat, 2, STATS = runif(ncol(mat))/1000, FUN = "+")
  
  clustering <- cutree(hclust(as.dist(1-cor(t(mat))), method = "complete"),  h = 1-corThreshold)
  clust.dt <- data.table(protein = names(clustering), cluster = sprintf("cluster.%04d", clustering))
  clust.dt[, size := .N,by = cluster]
  clust.dt[size < minClusterSize, cluster := 'cluster.0000']
  clust.dt[size < minClusterSize, size := .N,by = cluster]
  
  if (length(emptyProteins) > 0){
    # with the most missingProteins, proportionately
    missingCluster <- clust.dt[protein %in% emptyProteins, .N/size, by= cluster][order(-V1)]$cluster[1]
    clust.dt[cluster == missingCluster, cluster := "cluster.xxxx"]
    
  }
  
  computeClusterCenters(clust.dt, mat)
  computeClusterMaxes(clust.dt, mat)
  computeClusterSimilarityOrder (clust.dt, mat)
  
  return (clust.dt) 
}

#  source ("../../bp_utils/enrichmentTestFunctions.R")
enrichClusterTable <- function(cluster.dt, gmt, numProcessors = 8,...){
  enrich.dt <- enricherOnGroups(cluster.dt, geneColumn = "protein",
                                groupColumns = "cluster", term2gene.gmt = gmt,
                                universe = unique(cluster.dt$protein),
                                minGSSize = 4, maxGSSize = 2000,
                                numProcessors = numProcessors, ...)
  return (enrich.dt)
}

fastEnrichClusterTable <- function (cluster.dt, gmt, minMatchSize = 3){
  mutualProteins <- intersect(cluster.dt$protein, unique(gmt$gene))
  
  # size of gene sets
  bgCount <- gmt[gene %in% mutualProteins, .(setSize = .N), by = ont]
  
  # gene set overlap with group
  inGroupCount <- cluster.dt[gmt, , on = c(protein = "gene")][!is.na(cluster)][, .(groupANDsetSize = .N, proteins = paste0(protein, collapse = "/")), by = .(cluster, ont)]
  
  # group size
  groupDenom <- cluster.dt[protein %in% mutualProteins, .(groupSize = .N), by= cluster]
  
  contingencyStats <- inGroupCount[groupDenom, , on = "cluster"][bgCount, , on = "ont"]
  contingencyStats[, universeSize := length(mutualProteins)]
  
  # black and white balls in urn language. 
  contingencyStats[, whiteBallsDrawn := groupANDsetSize]
  contingencyStats[, totalWhiteInUrn := setSize]
  contingencyStats[, totalBlackInUrn := universeSize - setSize]
  contingencyStats[, totalBallsDrawn := groupSize]
  
  contingencyStats[ , log10P := phyper(whiteBallsDrawn-1, totalWhiteInUrn, totalBlackInUrn, totalBallsDrawn, lower.tail = FALSE, log.p = TRUE)/log(10)]
  
  contingencyStats[, TT := groupANDsetSize]
  contingencyStats[, TF := groupSize - groupANDsetSize]
  contingencyStats[, FT := setSize - groupANDsetSize]
  contingencyStats[, FF := universeSize - FT - TF - TT]
  contingencyStats[, oddsRatio := TT * FF/(TF * FT)]
  
  
  return(contingencyStats[groupANDsetSize>=minMatchSize, .(cluster, ont, groupANDsetSize, groupSize, setSize, universeSize, oddsRatio, log10P, proteins)])
}

fastEnrichTableBestPerProtein <- function(enrich.dt){
  perProtein <- enrich.dt[, .(protein = unlist(strsplit(geneID, "/"))), by= .(cluster, ID, GeneRatio, BgRatio, pvalue, p.adjust)]
  setorder(perProtein, pvalue, na.last = TRUE)
  perProtein[, symbol := multiUniprots2multiGenes(protein)]
  bestPerProtein <- perProtein[, .SD[1], by = .(protein)]
  return (bestPerProtein[])
  
}


enrichTableBestPerProtein <- function(enrich.dt){
  perProtein <- enrich.dt[, .(protein = unlist(strsplit(geneID, "/"))), by= .(cluster, ID, GeneRatio, BgRatio, pvalue, p.adjust)]
  setorder(perProtein, pvalue, na.last = TRUE)
  perProtein[, symbol := multiUniprots2multiGenes(protein)]
  bestPerProtein <- perProtein[, .SD[1], by = .(protein)]
  return (bestPerProtein[])
}


# for ordering of clusters:
computeClusterCenters <- function(cluster.dt, mat){
  .centerMassOne <- function(subMat){
    sum(colSums(subMat) * 1:ncol(subMat))/sum(subMat)
  }
  
  .centerMassOne <- function(subMat){
    rowCenters <- apply(subMat, 1, function(x){sum((1:ncol(subMat)) * x)/sum(x)} )
    median(rowCenters)
  }
  
  
  cluster.dt[, fractionCenter := .centerMassOne(mat[protein,]), by= cluster]
}

computeClusterMaxes <- function(cluster.dt, mat){
  .peakClusterOne <- function(subMat){
    as.numeric(median(peaksPerProtein <- apply(subMat,1, which.max)))
  }
  cluster.dt[, fractionMax := .peakClusterOne(mat[protein,]), by= cluster ]  
}


computeClusterSimilarityOrder <- function (cluster.dt, mat){
  specialClusters <-  c("cluster.0000" , "cluster.xxxx")
  splits <-split(cluster.dt[!cluster %in% specialClusters], by = "cluster")
  cluster.mat <- do.call (rbind, lapply(splits, function(dt)colSums(mat[dt$protein,])))
  
  cor <- cor(t(cluster.mat))
  
  cor[is.na(cor)] <- -1
  
  ddr <- as.dendrogram( hclust(as.dist(1-cor)))
  # desiredOrder based on fractionMax
  desiredOrder <- cluster.dt[cluster %in% labels(ddr)][order (fractionMax), unique(cluster)]
  ddr <- dendextend::rotate(ddr, desiredOrder)
  
  
  
  ordering <- order.dendrogram(ddr)
  namesInOrder <- rownames(cor)[ordering]
  ordering <- setNames(1:length(namesInOrder), namesInOrder)
  
  cluster.dt[, similarityOrder := ordering[cluster]]
  for(specialCluster in specialClusters){
    maxOrder <- max(cluster.dt$similarityOrder, na.rm = TRUE)
    cluster.dt[cluster == specialCluster, similarityOrder := maxOrder + 1]
  }
}



clusterOrderProteins <- function (cluster.dt, clustersOI = NULL){
  if (is.null(clustersOI))clustersOI <- unique(cluster.dt$cluster)
  specialClusterOffsets <- c("cluster.0000" = 1000, "cluster.xxxx" = 2000)

  #cluster.dt[cluster %in% clustersOI][order(fractionMax)]$protein
  # cluster.dt[cluster %in% clustersOI][order(ifelse(cluster %in% names(specialClusterOffsets),
  #                                                  fractionCenter + specialClusterOffsets[cluster],
  #                                                  fractionCenter))]$protein

  cluster.dt[cluster %in% clustersOI][order(similarityOrder)]$protein
  
  }


secHeatmap <- function (mat, cluster.dt, clustersOI = NULL, maxIntensity = 0.3){
  csDenom <- 50/maxIntensity
  proteinOrder <- clusterOrderProteins(cluster.dt, clustersOI)
  sec.mat <- t(mat[proteinOrder,])
  
  hm.sec <- Heatmap(sec.mat,
                    name = "Intensity (total scaled)",
                    column_split = factor(cluster.dt[proteinOrder,, on = "protein"]$cluster,
                                          levels = unique(cluster.dt[proteinOrder,, on = "protein"]$cluster)),
                    column_gap = unit(0.1, "mm"),
                    row_labels = ifelse(1:nrow(sec.mat) %% 5, "", 1:nrow(sec.mat)),
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    col <- circlize::colorRamp2(breaks = (0:50)/csDenom, colors = viridis::magma(51,direction = -1)),
                    show_column_names = FALSE, 
                    column_title_rot = -90,
                    column_title_gp = gpar(fontsize = 5),
                    row_names_gp = gpar(fontsize = 6),
                    
                    #bottom_annotation = columnAnnotation(GOenrich = sqrt(cluster.dt[order(fractionCenter)]$bestPP)),
                    height = unit(2,"inch"),
                    use_raster = FALSE,
                    heatmap_legend_param = list(direction = "horizontal"))
  
  return (hm.sec)
}

enrichHeatmap <- function(enrich.dt, cluster.dt, minPValue, clustersOI = NULL,
                          enrichColors = c(viridis::viridis(12, direction = -1)[1:12]) ){
  proteinOrder <- clusterOrderProteins(cluster.dt, clustersOI)
  # best enriched term per cluster
  terms <- enrich.dt[order (log10P), .SD[1], by = cluster][10^log10P< minPValue, unique(ont)]
  
  perProteinEn.dt <- enrich.dt[ont %in% terms, 
                               .(protein = unlist(strsplit(proteins, "/"))),
                               by = .(cluster, ont, log10P)
                               ] # any enrichment per protein
  
  enPG.mat <- dcast(perProteinEn.dt[ proteinOrder, , on = "protein"],
                    protein~ont, value.var = "log10P") |>
    as.matrix(rownames = "protein") |> t()
  
  enPG.mat <- enPG.mat[terms,] # removes the "NA"
  
  enPG.mat <- - enPG.mat
  #enPG.mat[is.na(enPG.mat)] <- 0
  
  # order columns
  enPG.mat <- enPG.mat[, proteinOrder]
  
  #order rows  
  enPG.mat <- enPG.mat[order(apply(enPG.mat, 1, which.max)),]
  
  
  hm.en <- Heatmap (sqrt(enPG.mat), name = "sqrt Enrichment",
                    
                    col = enrichColors,
                    na_col = "grey95",
                    
                    
                    show_column_names = FALSE, cluster_rows = FALSE,
                    row_names_gp = gpar(fontsize = 7),
                    row_split = factor(rownames(enPG.mat), levels = rownames(enPG.mat)),
                    gap = unit(0.5, "mm"), row_title_rot = 0, 
                    row_title_gp = gpar(fontsize = 7), cluster_columns = FALSE, use_raster = FALSE, 
                    heatmap_legend_param = list(direction = "horizontal"))
  return (hm.en)
}






# scaling ----




scaleByTotalIntensity <- function(secLong.dt){
  secLong.dt[, intensity_totalScaled := intensity/(sum(intensity, na.rm = TRUE)), by= .(sample, protein)]
}


scaleByMaxIntensity <- function(secLong.dt){
  secLong.dt[, intensity_maxScaled := intensity/(max(intensity, na.rm = TRUE)), by= .(sample, protein)]
}
scaleByMaxIntensity_global <- function(secLong.dt){
  secLong.dt[, intensity_maxScaled := intensity/(max(intensity, na.rm = TRUE)), by= .( protein)]
}



# matrices ----

#' @param scaleDenom total/max. What to use for the denominator for scaled intensity, total=sum(intensity) or max(intensity)

scaledIntensityMatrices <- function(secLong.dt, scaleDenom = "total", reorder = TRUE){
  if(scaleDenom == "total" & !"intensity_totalScaled" %in% colnames(secLong.dt))
    scaleByTotalIntensity(secLong.dt)
  if(scaleDenom == "max" & !"intensity_maxScaled" %in% colnames(secLong.dt))
    scaleByMaxIntensity(secLong.dt)

  
  allProteins <- unique(secLong.dt$protein)
  
  v_var <- ifelse(scaleDenom != "none", sprintf("intensity_%sScaled", scaleDenom), 'intensity')
  
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


interpolateMissingFractions <- function(matInput){
  mat <- matInput
  missingFractions <- which (apply(mat, 2, sum, na.rm = TRUE) == 0)
  print (missingFractions)
  #ignore the boundaries
  missingFractions <- setdiff(missingFractions, c(1,ncol(mat)))
  if (length(missingFractions) == 1){
    missingFractions <- sort(missingFractions)
    if (any(diff(missingFractions) > 1)){
      stop("only single missing fractions implemented")
    }
    for (mf in missingFractions){
      mat[, mf] <- (mat[, mf-1] + mat[, mf-2])/2
    }
    
  }
  return(mat)
}




# Heatmaps ----

intensityHeatmaps <- function(intMats, intensityName = "Scaled Intensity", topOfColorRange = 0.3, showRowNames = NULL,...){
  denom = 50 / topOfColorRange
  colorFun <- circlize::colorRamp2(breaks = (0:50)/denom, colors = viridis::magma(51,direction = -1))
  samples <- names(intMats)
  
  sample <- samples[1]
  
  if (is.null(showRowNames)){
    showRowNames <- nrow(intMats[[sample]]) < 100
  }
  
  hml <- Heatmap (intMats[[sample]] ,
                  name = intensityName,
                  cluster_rows = FALSE,
                  row_dend_reorder = FALSE,
                  cluster_columns = FALSE,
                  col = colorFun ,
                  show_row_names = (is.null(showRowNames) & nrow(intMats[[sample]]) > 100) | showRowNames,
                  row_names_side = "left",
                  column_names_gp = gpar(fontsize = 5),
                  column_labels = ifelse(as.integer(colnames(intMats[[sample]])) %% 5 == 0, colnames(intMats[[sample]]), ""),
                  column_title = sample,
                  
                  #first only:
                  show_heatmap_legend = TRUE, 
                  row_title = sprintf ("%d Proteins", nrow(intMats[[sample]])),
                  ...
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
                            show_heatmap_legend = FALSE,
                            ...)      
    }
    
  }
  return (hml)
}


# Normalization and Outlier detection by fitting local cubics ----

#' @param sampleTerm additive/interaction whether to fit a different curve per sample or allow only
#'                   an offset per sample.  This is per window, so different offsets will be used, 
#'                   and some amount of curve-unique-to-sample persists weven with "additive"

fitLocalCubics<- function (qcSummary,window =15,extend  = 3,
                           #choose one:
                           sampleTerm = c("additive",
                                          "interaction",
                                          "interactionByTreatment",
                                          "singleSample")[1]){
  additiveModel <-  as.formula ("medPolishIntensity~poly(fraction,3) + sample")
  interactionModel <-     as.formula ("medPolishIntensity~poly(fraction,3) * sample")
  interactionByTreatmentModel <- as.formula ("medPolishIntensity~poly(fraction,3) + poly(fraction,3):treatment + sample")
  singleSampleModel <- as.formula("medPolishIntensity~poly(fraction,3)")
  
  model <- c(interaction = interactionModel,
             additive =additiveModel,
             interactionByTreatment = interactionByTreatmentModel,
             singleSample = singleSampleModel)[[sampleTerm]]
  
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
  qcSummary[localCubicFits[, .(medianResidual = median(residual)), by = .(sample , fraction)],
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
plotNormAndOutlierFits <- function(qcSummary, localCubicFits){
  p <- ggplot (localCubicFits, aes(x = fraction)) + 
    geom_line(aes(y = fitted, group = fit), show.legend = FALSE, alpha = 0.2) + 
    geom_point(data = qcSummary, aes(y = medPolishIntensity, color = isOutlier), show.legend = FALSE) +
    theme_bw()
  
  if ("medianResidual" %in% colnames(qcSummary)){
    p <- p + geom_point(data = qcSummary, aes(y = medPolishIntensity - medianResidual), color= "grey20", size = 1,show.legend = FALSE)
  }
  
  
  if (all(c("treatment", "replicate")  %in% colnames(localCubicFits)) &
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
  diag(sim) <- 
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

windowedCorrelation <- function (intensity.mat, goodPeaks.mat,
                                 outerRadius = 6, peakRadius = 1,
                                 goldStandardInteractome = NULL,
                                 priorOdds = 0.01, thresholdR = 0.90){
  intensity.mat[is.na(intensity.mat)] <- 0.0
  start <- 1 + outerRadius
  stop <- ncol(intensity.mat) - outerRadius
  
  
  .doOneWindow <- function (start, intensity.mat, goodPeaks.mat, outerRadius, peakRadius){
    rowsWithCentralPeaks <- apply(goodPeaks.mat[, (start - peakRadius):(start + peakRadius)],1,any)
    if (length(rowsWithCentralPeaks) <= 1){
      warning ("For fraction ", start, " only 1 or 0 peaks found.")
      return ( list(cor = data.table(), stats = data.table()))
    }
    subMat <- intensity.mat[rowsWithCentralPeaks, (start-outerRadius):(start + outerRadius), drop = FALSE]
    
    peakIdx <- apply(goodPeaks.mat[rowsWithCentralPeaks, (start-peakRadius):(start + peakRadius), drop = FALSE],
                     1,
                     function(x)match(TRUE, x) )
    rs <- rowSums(subMat)
    
    if (nrow(subMat) <= 1){
      return ( list(cor = data.table(), stats = data.table()))
    }
    
    cor.dt <- setDT(reshape2::melt(cor(t(subMat), use = "pairwise"),
                                   varnames = c("protein1", "protein2"),
                                   value.name = "pearsonR")
    )[]
    
    cor.dt[, start := start]
    
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
      
      cor.dt <- cor.dt[pearsonR > thresholdR & protein1 != protein2] # limit to just the most promising
      cor.dt[, c("interactorLL", "bgLL") := .(interactorDensityFun(corScore), bgDensityFun(corScore))]
      cor.dt[, llRatio := interactorLL/bgLL ]
      cor.dt[, postOdds := llRatio * priorOdds ]
      cor.dt[, postProb := postOdds/(postOdds + 1) ]

      xRange <- c(-1, 16.5)
      stats.dt <- data.table(start = start, x = seq (from = xRange[1], to = xRange[2], length.out = 512))
      stats.dt[, c("interactorLL", "bgLL") := .(interactorDensityFun(x), bgDensityFun(x))]
      stats.dt[, llRatio := interactorLL/bgLL ]
      stats.dt[, postOdds := llRatio * priorOdds ]
      stats.dt[, postProb := postOdds/(postOdds + 1) ]
      
    }else{
      cor.dt <- cor.dt[pearsonR > thresholdR & protein1 != protein2] # limit to just the most promising
      stats.dt <- data.table()
    }
    
    
    #
    rs.dt <- data.table(protein = names(rs), portion = rs)
    cor.dt[rs.dt, prot1Portion := i.portion, on = c(protein1 = "protein")]
    cor.dt[rs.dt, prot2Portion := i.portion, on = c(protein2 = "protein")]
    cor.dt[, log2Ratio := log2(prot1Portion/prot2Portion)]
    
    #peakIdx
    pi.dt <- data.table(protein = names(peakIdx), subIndex = peakIdx)
    cor.dt[pi.dt, prot1Peak := start - peakRadius -1 + subIndex, on = c(protein1 = "protein")]
    cor.dt[pi.dt, prot2Peak := start - peakRadius -1 + subIndex, on = c(protein2 = "protein")]
    

    return ( list(cor = cor.dt, stats = stats.dt))
  }
  
  out.ls <- pbapply::pblapply(start:stop, .doOneWindow, intensity.mat, goodPeaks.mat, outerRadius, peakRadius)
  
  cor.dt <- lapply(out.ls, function(x)x$cor) |> rbindlist()
  stats.dt <- lapply(out.ls, function(x)x$stats) |> rbindlist()
  
  return (list (cor = cor.dt, stats = stats.dt))
  
}

# Gold Standard Interactors/Decoys ----

#' @param genes vector of gene names to limit string to
#' @param links.path path to string links file, 
#'                   example https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz
#' @param info.path path to string file with gene name etc info
#'                  example https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
#'                  or a data.table with columns `#string_protein_id` and `preferred_name`                  
#' @param combinedScoreThreshold only edges with combined score above this will be considered
#' @param stringDistThreshold protein pairs greater distance than this will be called a decoy
#' @param geneAliasFunction a function that will convert a list of genes to the canonical alias
#'                          `function(charGeneVector){... return(charGeneAliasVector)}`
#'                          The default, `identity`, is no conversion
#'                          See function `github/kroganlab/bp_utils/UniprotIDMapping.R :: geneAlias2officialGeneSymbol` as an example/possible
decoysFromString <- function (genes, 
                              links.path = "~/Downloads/9606.protein.physical.links.detailed.v12.0.txt.gz",
                              info.path = "~/Downloads/9606.protein.info.v12.0.txt.gz",
                              combinedScoreThreshold = 600,
                              stringDistThreshold = 5,
                              geneAliasFunction = identity){
  genes <- unique(genes)
  # remove KRT contaminants
  genes <- grep("^KRT", genes, invert = TRUE, value = TRUE)
  
  
  string <- fread (links.path)
  string <- string[combined_score > combinedScoreThreshold  ]
  
  if ("data.table" %in% class(info.path)){
    proteinNames <- info.path
  }else if (is.null(info.path)){
    proteinNames <- data.table (`#string_protein_id` = genes, preferred_name = genes)
  }else{
    proteinNames <- fread (info.path)
  }
  string[proteinNames, gene1 := i.preferred_name , on = c(protein1 = "#string_protein_id")]
  string[proteinNames, gene2 := i.preferred_name , on = c(protein2 = "#string_protein_id")]
  
  string[, alias1 := geneAliasFunction(gene1)]
  string[, alias2 := geneAliasFunction(gene2)]
  g <- igraph::graph_from_data_frame(string[, .(alias1, alias2)], directed = FALSE)
  # find distant genes in string
  rm.na <- function(x)x[!is.na(x)]
  dists <- igraph::distances(g, 
                             rm.na (match( genes, names(igraph::V(g)))),
                             rm.na (match( genes, names(igraph::V(g)))))
  distantGenes <- which(dists > stringDistThreshold, arr.ind = TRUE) |> as.data.table(keep.rownames = TRUE)
  # distantGenes is a data.table with columns rn, row, col. 
  # row and col are indeces to dimensions of dists matrix
  setnames(distantGenes, old = "rn", new = "gene1")
  distantGenes[, gene2 := colnames(dists)[col]]
  
  decoys <- unique(distantGenes[gene1 < gene2, .(gene1, gene2)])
  return (decoys)
}

scoreByGS <- function (sub.dt, denomDecoy, denomInteractor, column = "meanLL", groupByVariable = "treatment", pseudoCount = 1){
  # we need an encoding of gs levels so we can make decoy or interactor first:
  gsOtherLevels <- setdiff(unique(sub.dt$gs), c("decoy", "interactor"))
  sub.dt[, gs := factor(gs, levels = c("decoy", gsOtherLevels, "interactor"))]
  
  message ("Counting and scoring decoys...")
  setorderv(sub.dt, c(column, "gs"), order = c(-1, 1))
  
  sub.dt[, decoyCount := 0L]
  sub.dt[gs == "decoy", decoyCount := frankv(.SD,  column, order = -1, ties.method = "max"), by = eval(groupByVariable)]
  sub.dt[, decoyCount := cummax(decoyCount), by = eval(groupByVariable)] # requires decoy before others at same score, to copy count to other labeled data at same score
  sub.dt[, decoyCount := decoyCount + pseudoCount]
  sub.dt[, log10DecoyRate := log10(decoyCount) - log10(denomDecoy)]
  
  # reorder the gs labels, interactor before decoy
  message ("Counting and scoring interactors...")
  setorderv(sub.dt, c(column, "gs"), order = c(-1, -1))
  sub.dt[, interactorCount := 0L]
  sub.dt[gs == "interactor", interactorCount := frankv(.SD, cols = column, order = -1, ties.method = "max"), by = eval(groupByVariable)]
  sub.dt[, interactorCount := cummax(interactorCount), by = eval(groupByVariable)] # # requires interactor before others at same score, to copy count to other labeled data at same score
  sub.dt[, interactorCount := pseudoCount + interactorCount]
  sub.dt[, log10IntRate := log10(interactorCount) - log10(denomInteractor)]
  
  sub.dt[is.infinite(log10DecoyRate), log10DecoyRate := NA]
  sub.dt[is.infinite(log10IntRate), log10IntRate := NA]
  
  
  message ("Calculating ratios per score ", column, " ...")
  
  # ratios
  sub.dt[, log10RateRatioID := max(log10IntRate, na.rm = TRUE) -  max(log10DecoyRate, na.rm = TRUE), 
         by= c(eval(column), eval(groupByVariable))] #highest recovery at each value of score

  # take care of likelihood ratio apparently shrinking due to low numbers of itneractors at high scores
  # ... we expect ratio likelihood only increase as threshold increases
  #reverse order to increasing by score column
  setorderv(sub.dt, column, 1)
  sub.dt[, log10RateRatioID := cummax(log10RateRatioID), by = eval(groupByVariable)]
  
  #reverse back
  setorderv(sub.dt, column, -1)
}


#' @param genes vector of gene names to limit results to
#' @param corum.path path to corum file
#' 
#' @param string.links.path path to string links file, 
#'                   example https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz
#' @param string.info.path path to string file with gene name etc info
#'                  example https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
#'                  or a data.table with columns `#string_protein_id` and `preferred_name`
#' @param combinedScoreThreshold only edges with combined score above this will be considered
#' @param geneAliasFunction a function that will convert a list of genes to the canonical alias
#'                          `function(charGeneVector){... return(charGeneAliasVector)}`
#'                          The default, `identity`, is no conversion
#'                          See function `github/kroganlab/bp_utils/UniprotIDMapping.R :: geneAlias2officialGeneSymbol` as an example/possible

goldStandardPairs <- function (genes,
                               corum.path = "~/Downloads/corum_humanComplexes.txt",
                               string.links.path = "~/Downloads/9606.protein.physical.links.detailed.v12.0.txt.gz",
                               string.info.path = "~/Downloads/9606.protein.info.v12.0.txt.gz",
                               stringCombinedScoreThreshold = 600,
                               geneAliasFunction = identity){
  
  genes <- unique(genes)
  # remove KRT contaminants
  genes <- grep("^KRT", genes, invert = TRUE, value = TRUE)
  
  #corum
  if (!is.null(corum.path)){
    corumPairs <- unique(corumPairs(corum.path)[, .(gene1, gene2)])
  }else{
    corumPairs <- data.table (gene1 = c(), gene2 = c())
  }

  # string
  string <- fread (string.links.path)
  string <- string[combined_score > stringCombinedScoreThreshold  ]
  
  
  if ("data.table" %in% class(string.info.path)){
    proteinNames <- string.info.path
  }else if (is.null(string.info.path)){
    proteinNames <- data.table (`#string_protein_id` = genes, preferred_name = genes)
  }else{
    proteinNames <- fread (string.info.path)
  }
  string[proteinNames, gene1 := i.preferred_name , on = c(protein1 = "#string_protein_id")]
  string[proteinNames, gene2 := i.preferred_name , on = c(protein2 = "#string_protein_id")]
  
  string[, alias1 := geneAliasFunction(gene1)]
  string[, alias2 := geneAliasFunction(gene2)]
  
  stringPairs <- string[alias1 < alias2, .(gene1 = alias1, gene2 = alias2)]
  
  combinedPairs <- rbindlist( list(corum = corumPairs,
                                   string = stringPairs),
                              idcol = 'source'
  )[
    , .(source = paste0(source, collapse = ";")),
    by = .(gene1, gene2)]
  
  return (combinedPairs[gene1 %in% genes & gene2 %in% genes])
  
}

#' returns complex id and name together with gene1, gene2
corumPairs <- function(corum.path = "~/Downloads/corum_humanComplexes.txt",
                       geneAliasFunction = identity){
  # corum
  corum <- fread (corum.path)
  .allByAll <- function(genes){
    data.table(gene1 = genes)[, .(gene2 = genes), by = gene1][]
  }
  corumPairs <- corum[, .allByAll(unlist(strsplit(subunits_gene_name, ";"))), by = .(complex_id, complex_name)]
  corumPairs[, alias1 := geneAliasFunction(gene1)]
  corumPairs[, alias2 := geneAliasFunction(gene2)]
  
  corumPairs <- unique(corumPairs[alias1 < alias2, .(gene1 = alias1, gene2 = alias2, complex_id, complex_name)])
  return(corumPairs)
}


#' returns complex id and name together with gene1, gene2
corumPairs <- function(corum.path = "~/Downloads/corum_humanComplexes.txt"){
  # corum
  corum <- fread (corum.path)
  .allByAll <- function(genes){
    data.table(gene1 = genes)[, .(gene2 = genes), by = gene1][]
  }
  corumPairs <- corum[, .allByAll(unlist(strsplit(subunits_gene_name, ";"))), by = .(complex_id, complex_name)]
  
  corumPairs <- unique(corumPairs[gene1 < gene2, .(gene1, gene2, complex_id, complex_name)])
  return(corumPairs)
}



# Matrix Smoothing ----

#' @description
#' Smooths the rows of a matrix using gaussian kernel smoothing, row-wise
#' @param a matrix to smooth, like a SEC MS matrix. Rows are proteins, columns are equally spaced fractions
#' @param ... parameters to pass on, most relevantly `bandwidth`. Default bandwidth is 2.68, or SD = 1, same as PCprophet
#' @value a matrix of same dimensions and as input
smoothRowsOfMatrix <- function(mat, ...){
  message ("Smoothing matrix (don't accidentally smooth twice)...")
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

findAllPeaksInSingleMatrix <- function(originalMatrix, ...){
  intMat <- smoothRowsOfMatrix(originalMatrix, ...)
  if (any(is.na(intMat)))
    intMat[is.na(intMat)] <- 0.0
  message ("peak locating...")
  allPeaks <- apply(intMat, 1, pracma::findpeaks, nups = 2, ndowns = 2 )
  allPeaks <- rbindlist(lapply(allPeaks, as.data.table), idcol = "protein")
  setnames(allPeaks, old = c("V1", "V2", "V3", "V4"), new =  c("peakHeight", "peakLocation", "peakStart", "peakEnd"))
  
  renumberFractions <- FALSE
  # CHECK FRACTION NAMES for unexpected things.
  # I favor warnings here, though these are probably errors most of the time
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
    
    if (!(all(fractionNumbers == 1:length(fractionNumbers))))
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

  # kernel-based peak center
  
  .cofmKernel <- function (protein, peakCenter, start, stop){
    radius <- min(abs(peakCenter - c(start, stop)))
    fractions <- peakCenter + (-radius:radius)
    fractionNames <- as.character(fractions) # as.character to support matrix columns that don't start at 1
    
    # we use half the bandwidth of the original smoothing
    # this prevents large peak shifts when a peak is a shoulder of another peak
    kernel <- ksmooth(x = fractions, y = intMat[protein, fractionNames ],
                      kernel = "normal", bandwidth =  1.34,
                      x.points = seq(from = peakCenter-radius,
                                     to = peakCenter + radius,
                                     by = 0.01))
    result <- kernel$x[which.max(kernel$y)]
    # # debugging:
    # proteinOI <- protein
    # if (abs(allPeaks[protein == proteinOI & peakStart == start & peakEnd == stop]$peakLocation - result) > 1){
    #   message ("trouble")
    # }
    # # end debugging
    return (result)
  }
  
  .cofmFullPeak <- function (protein, peakCenter, start, stop){
    # in case there is a "lopsided" peak, keep it centered 
    radius <- min(abs(peakCenter - c(start, stop)))
    .cofmCenterN(protein, peakCenter, peakWindow = radius)
  }
  
  message ("center of mass calc...")
  allPeaks[, cofmN := .cofmCenterN(protein, peakLocation), by = .I]
  
  # inaccurate when peak is lopsided:
  #allPeaks[, cofmFull := .cofmFullPeak(protein, peakLocation, peakStart, peakEnd), by = .I]
  
  # accurate but slow:
  #allPeaks[, cofmKernel := .cofmKernel(protein, peakLocation, peakStart, peakEnd), by = .I]
  message ("...done")
  


  return(allPeaks)
}





#' @description
#' Takes a single sec matrix, NOT smoothed ahead of time, and locates good peaks
#' First step is to create a smoothed matrix, but it uses the not-smoothed for detection of 'good'
#' Good peaks are defined as:
#' 1) peak shaped up-up-center-down-down
#' 2) center > 0.01 intensity (1% of total protein in one fraction, after smoothing)
#' 3) center N peaks, given by centerRadius, are non-zero and non-missing
#' 4) a peak with "contrast" determined by coefficient of variation in peaks in varRadius
#'     this avoids flat peaks but also may exclude those peaks wiht high bases
#' @return A data.table of peaks with a goodPeak column labeling "goodPeaks"

goodPeaksTableFromIntensityMatrix <- function (intMat, minPeakHeight = 0.01, centerRadius = 2, varRadius = 4, minCV = 0.2){
  message("This now expects a not-smoothed matrix. Make sure you are not smoothing ahead of time")
  stopifnot(!is.null(rownames(intMat)),
            !is.null(colnames(intMat)))
  
  allPeaks <- findAllPeaksInSingleMatrix(intMat) # this function does the smoothing
  
  # detect if fractions around peak are complete (>0, not missing)
  allPeaks[, 
           centerComplete := apply(intMat[protein, peakLocation + (-centerRadius:centerRadius),
                                           drop = FALSE],
                                   1, function(x)all(x>0)),
           by = peakLocation]
  
  # variance
  allPeaks[, 
           var := apply(intMat[protein, 
                               intersect(1:ncol(intMat), peakLocation + (-varRadius:varRadius)),
                               drop= FALSE],
                        1, var, na.rm = TRUE),
           by = peakLocation]  
  # mean
  allPeaks[, 
           mean := apply(intMat[protein,
                                intersect(1:ncol(intMat), peakLocation + (-varRadius:varRadius)),
                                drop = FALSE],
                         1, mean, na.rm = TRUE),
           by = peakLocation]  
  
  #cv
  allPeaks[, cv := sqrt(var)/mean]
  
  allPeaks[, goodCV := cv > minCV ]
  allPeaks[, goodHeight := peakHeight > minPeakHeight]
  
  allPeaks[, goodPeak := goodHeight & centerComplete & goodCV]
  return (allPeaks[])
}



goodPeaksMatFromPeaksTable <- function(intMat, allPeaks){ 
  #create dummy rows that we use to expand allPeaks table to include allProteins and allFractions
  allProteins <- rownames(intMat)
  allFractions <- as.integer(colnames(intMat))
  dummyRows <- data.table(protein = allProteins,
                          peakLocation = allFractions,
                          goodPeak = FALSE) |>
    suppressWarnings() # to ignore warning about recycling with remainder shorter vector
  
  peakMat <- dcast(rbind(dummyRows, allPeaks, use.names = TRUE, fill = TRUE), 
                   protein~peakLocation, value.var = "goodPeak",
                   fun.aggregate = any
  ) |> as.matrix(rownames = "protein")
  
  peakMat[is.na(peakMat)] <- FALSE
  
  return( peakMat[rownames(intMat), colnames(intMat)] )
}





#' @description
#' Takes a single sec matrix, NOT smoothed ahead of time, and locates good peaks
#' First step is to create a smoothed matrix, but it uses the not-smoothed for detection of 'good'
#' Good peaks are defined as:
#' 1) peak shaped up-up-center-down-down
#' 2) center > 0.01 intensity (1% of total protein in one fraction, after smoothing)
#' 3) center N peaks, given by centerRadius, are non-zero and non-missing
#' 4) a peak with "contrast" determined by coefficient of variation in peaks in varRadius
#'     this avoids flat peaks but also may exclude those peaks wiht high bases

goodPeaksMatrixFromIntensityMatrix <- function (intMat, minPeakHeight = 0.01, centerRadius = 2, varRadius = 4, minCV = 0.2){
  warning ("Obsolete function. See goodPeaksTableFromIntensityMatrix")
  message("This now expects a not-smoothed matrix. Make sure you are not smoothing ahead of time")
  stopifnot(!is.null(rownames(intMat)),
            !is.null(colnames(intMat)))
  
  allPeaks <- findAllPeaksInSingleMatrix(intMat) # this function does the smoothing
  
  # detect if fractions around peak are complete (>0, not missing)
  allPeaks[, 
           centerComplete := apply(intMat[protein, peakLocation + (-centerRadius:centerRadius)],
                                    1, function(x)all(x>0)),
           by = peakLocation]
  
  # variance
  allPeaks[, 
           var := apply(intMat[protein, 
                               intersect(1:ncol(intMat), peakLocation + (-varRadius:varRadius))],
                                    1, var, na.rm = TRUE),
           by = peakLocation]  
  # mean
  allPeaks[, 
           mean := apply(intMat[protein,
                                intersect(1:ncol(intMat), peakLocation + (-varRadius:varRadius))],
                        1, mean, na.rm = TRUE),
           by = peakLocation]  
  
  #cv
  allPeaks[, cv := sqrt(var)/mean]
  
  allPeaks[, goodPeak := (peakHeight > minPeakHeight & centerComplete & cv > minCV)]
  
  
  # create goodPeak matrix... should really be its own function
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
  #/ end createGoodPeakMatrix, probably its own function

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
standardizeOnePeakTableToStandard <- function (otherPeaks, standardPeaks, sec.dt, sampleName, minPeaksPerFraction = 50,
                                               firstPeak = 1, lastPeak = 72,
                                               doPlots = TRUE, fitPortion = 0.75, startFitAtPeak = 15){
  message ("Matching peaks between samples...")
  matchedPeaks <- matchTwoPeakTables (otherPeaks, i.peaks.dt = standardPeaks)
  matchedPeaks <- matchedPeaks[cofmN >= startFitAtPeak & i.cofmN >= startFitAtPeak]
  
  # # barplots of matchedPeaks per fraction
  # if (doPlots){
  #   p <- ggplot(matchedPeaks, aes(x = peakLocation)) + geom_bar() + ggtitle(label = sprintf("%s Matched Peaks ", sampleName))
  #   print (p)
  # }
  
  goodFractions <- matchedPeaks[goodPeak == TRUE & i.goodPeak == TRUE, .N, by = peakLocation][N > minPeaksPerFraction, peakLocation]
  #goodFractions <- goodFractions[goodFractions > startFitAtPeak]
  
  subData <- matchedPeaks[goodPeak == TRUE & i.goodPeak == TRUE][
    peakLocation %in% goodFractions][
      ,.SD[i.cofmN %between% quantile(i.cofmN, c( (1-fitPortion)/2, fitPortion + (1-fitPortion)/2))], by = peakLocation] # middle 75% only, or fitPortion
      
      #abs(i.cofmN-cofmN) < maxPeakShift]
  
  if(doPlots){
    p <- ggplot (subData, 
                 aes(x = i.cofmN, y = cofmN)) +
      scale_x_continuous(name = "Peak location in standard (cofmN)") + 
      scale_y_continuous(name = sprintf("Peak location in %s (cofmN)", sampleName)) + 
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
#' @param ... arguments passed to `standardizeOnePeakTableToStandard` that does the pairwise adjusting,
#'            relevant arguments include
#'    * `startFitAtPeak` the integer value at which to start fits. Useful when initial peaks are difficult to align
#' 
#' 
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
      sec.dt[sample == sampleName, standardFraction := as.double(fraction)]
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
    .cleanPeakJoiningColumns(peak.dt)
  }
  invisible(allPeakTables)
}

# molecular weights by proteins/fractions ----

#' @description
#' Expected input is two columns that represent fraction number and molecular weight in kildaltons
#' Regardless of input names, output is titled fraction,mw and `mw` is 1000 * input molecular weight
#' 
loadStandardsFractionToMW <- function(path = "~/Downloads/cal_biosep.txt", ...){
  standards  <- fread (path, ...)
  setnames(standards, c("fraction", "mw"))
  standards[, mw := mw * 1000] # expect input in kD, convert to D
  return(standards[])
}

#' @description
#' Receives a two column table, First column should be uniprot or other designator of proteins used in the dataset
#' Second column is a mw in Daltons
#' 
loadUniprotToMW <- function(path = "~/Downloads/mw_uniprot_Accession.txt"){
  proteinMW <- fread (path)
  setnames(proteinMW, c("protein", "mw"))
  return (proteinMW[])
}

#calculateFractionMassConverters(loadStandardsFractionToMW())

calculateFractionMassConverters <- function(st){
  # conversion is based on a simple line log10(mw)~fraction
  
  ab <- coefficients(lm(log10(mw)~fraction, data = st))
  
  fraction2Mass <- function (fraction){
    y = ab["(Intercept)"] + ab["fraction"] * fraction
    mw = 10^y
    return (mw)
  }
  
  mass2Fraction <- function (mass){
    fraction = (log10(mass) - ab["(Intercept)"])/ab["fraction"]
    return (fraction)
  }
  
  return (list(fraction2Mass = fraction2Mass, mass2Fraction = mass2Fraction))
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
  # we work in log space, so zeros need to get replaced with NA

  onePeakSec.dt[intensity_totalScaled > 0, log2_intensity_totalScaled := log2(intensity_totalScaled)]
  
  lmout <- lm( log2_intensity_totalScaled~ poly(standardFraction, 4)*treatment, onePeakSec.dt)
  lmout.noInteraction <- lm( log2_intensity_totalScaled~ poly(standardFraction, 4)+treatment, onePeakSec.dt)
 
  
  #onePeakSec.dt[, fitted := 2^predict(lmout)]
  
  if (!is.null(doPlotFunction)){
    fittedCurve <- data.table(treatment = unique(onePeakSec.dt$treatment))[, .(standardFraction = seq(from = min(onePeakSec.dt$standardFraction), to = max(onePeakSec.dt$standardFraction), length = 100)), by = treatment]
    fittedCurve[, intensity_totalScaled := 2^predict (lmout, newdata = fittedCurve)] 
    fittedCurve[, intensity_totalScaled.noInteraction := 2^predict (lmout.noInteraction, newdata = fittedCurve)] 
    
    p <- ggplot(onePeakSec.dt, aes(x =standardFraction, y = intensity_totalScaled,color = treatment)) +
      geom_point(aes(shape = as.factor(replicate))) + 
      geom_line(aes(color = treatment, group = sample)) +
      geom_line(data= fittedCurve, lwd = 1.5, alpha = .5) +
      geom_line(data = fittedCurve, lwd = 1, lty = "dashed", mapping = aes(y =intensity_totalScaled.noInteraction )) + 
      theme_bw()
    
    if ("function" %in% class(doPlotFunction)){doPlotFunction(p)}else{print(p)}
    
  }
  return(data.frame(anova(lmout)))
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
  stopifnot (length(unique(sec.dt$protein)) <= 1) # only intended for one protein
  if (length(unique(sec.dt$protein)) == 0){
    return(data.table(error = "No data for protein"))
  }
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



## MG functions

#' return column with max consecutive detections per protein
MaxConsecutiveDetections <- function(secLong.dt, idcol='peptideSequence', intsCol='intensity', detectionCutoff=0, plot=F){
  
  .oneMatrix <- function(sub.dt){
    mat <- dcast(sub.dt, sub.dt[[idcol]]~fraction, value.var = intsCol) %>% 
      as.matrix(rownames = 1)
    mat[is.na(mat)] <- 0.0
    mat[order(rownames(mat)),]
  }
  
    #' Identify the max run of consecutive measurmeents
  .maxConsecFractionMeasurments <- function(row, cutoff=detectionCutoff){
    detectVec <- rle(row > cutoff)
    return(max(detectVec$lengths[detectVec$values == TRUE]))
  }
  
  # one matrix per sample 
  mats <- lapply(split(secLong.dt, list(secLong.dt$sample)), .oneMatrix)
  
  consec.dt <- pbapply::pblapply(mats, function(x) {
    apply(x, 1, function(feature){ .maxConsecFractionMeasurments(row = feature, cutoff = detectionCutoff) }) %>% 
    as.data.table(keep.rownames=T)
  }) %>% 
  rbindlist(idcol='sample')

  setnames(consec.dt, new=c('sample', 'feature', 'consecutiveDetections'))
 
  if (plot){
    
    g <- ggplot(consec.dt, aes(x=consecutiveDetections, col=sample)) +
      stat_ecdf() +
      labs(title='Distribution of maximum consecutive detections per protein') +
      theme_bw()
    
    p <- ggplot(consec.dt, aes(x=consecutiveDetections, col=sample)) +
      geom_density() +
      labs(title='Distribution of maximum consecutive detections per protein') +
      theme_bw()
 
    print(p); print(g)
  }
  
  sec.dt <- merge(x=secLong.dt, y=consec.dt, by.x=c('sample', idcol), by.y=c('sample', 'feature'), all.x=T)
  stopifnot(nrow(sec.dt) == nrow(secLong.dt))
  return(sec.dt)
}


#' function to interpolate missing values for missing/outlier fractions
#' requires a sec.dt (longformat) with sample, treatment, replicate and intensity columns
#' requires a qc.dt generated from the qcSummaryTable() function
#' datatable join to replace values in problematic fractions with NA
#' maxGap; treshold for N consec missing values for interpolation default is 1
interpolateMissingAndOutlierFractions <- function(secLong.dt, qc.dt, fractions, maxGap=1){
  
  # create copy to avoid modifying ori DT
  sec.dt <- copy(secLong.dt)
  # zero out the intensities of problematic/missing fractions
  sec.dt[, ori.intensity := intensity]
  sec.dt[qc.dt[isOutlier == TRUE], on=.(sample, fraction), intensity := NA]
  
  # add NA cols for missing fractions
  .addMissingFractions <- function(subDT){
    subMat <- dcast(subDT, protein~fraction, value.var='intensity') %>% 
      as.matrix(rownames='protein')
    
    if (!all(colnames(subMat) == fractions)){
      message('Some fractions are missing. Adding missing fractions populated with NA...')
      message('Missing fractions:\n', setdiff(fractions, colnames(subMat)))
      subMat <- subMat[, match(fractions, colnames(subMat)), drop=FALSE] #missing fractions assigned an NA col. drop=FALSE to avoid collapsing to vector
      colnames(subMat) <- fractions
    }
    return(subMat)
  }
  
  .interpolateOutlierFractions <- function(subMat, sampleOI){
    # interpolate intensity values in matrix rows
    interpMat <- apply(subMat, 1, function(x) zoo::na.approx(x, na.rm=F, maxgap=maxGap)) %>% 
      t()
    # fractions to update; restrict to outlier and/or missing
    fractionsOI <- c(qc.dt[sample == sampleOI & isOutlier == TRUE, unique(fraction)],
                     setdiff(fractions, sec.dt[sample == sampleOI, unique(fraction)])
                     )
    
    message('Handling problematic fractions for ', sampleOI, ': ', paste0(fractionsOI, collapse=','))
    colsToupdate <- colnames(subMat) %in% fractionsOI
    
    # now apply the values from the interpolated matrix to the original
    subMat[, colsToupdate] <- interpMat[, colsToupdate]
    subdt <- setDT(reshape2::melt(subMat))
    setnames(subdt, new=c('protein', 'fraction', 'intensity'))
    
    # add a flag to identify if the value is interpolated
    subdt[, interpolated := FALSE]
    subdt[fraction %in% fractionsOI, interpolated := TRUE]
    return(subdt)
  }

  # handling missing/outlir fractions
  sec.list <- pbapply::pblapply(split(sec.dt, sec.dt$sample), .addMissingFractions)
  
  # interpolate missing values
  message('interpolating missing values...')
  interp.list <- pbapply::pblapply(names(sec.list), function(n){.interpolateOutlierFractions(subMat = sec.list[[n]], sampleOI = n)}) 
  names(interp.list) <- names(sec.list)
  interp.dt <- rbindlist(interp.list, idcol='sample')
  
  interp.dt <- merge(x=interp.dt, y=sec.dt[, -c('intensity')], by=c('protein', 'sample', 'fraction'), all.x=T)
  
  # santy checks
  #  nrows for non-missing fractions and intensity vals for non-interpolated should match between input and output
  stopifnot( nrow(interp.dt[qc.dt, , on=.(sample, fraction)]) == nrow(sec.dt))
  stopifnot( nrow(interp.dt[interpolated == FALSE & (intensity != ori.intensity), ]) == 0 )
  return(interp.dt[, -c('ori.intensity')]) 
}