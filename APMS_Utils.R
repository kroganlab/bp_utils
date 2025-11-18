


#' requires MASS
#' @param subInt a data.table of apms intensity, subsetted to a single prey. 
#'               required columns are bait,log2Intensity
#'               columns added/overwritten are numMissing,predicted,SD,replaceNA
#'               missing values must be explicit, rows with `TRUE == is.na(log2Intensity)`
#'               
#' @param minOverall the minimum log2Intensity for replaceNA when no other infomation available
#' @param minOffset the minimum offset for log2Intensity to offset from minimum observed in complete.cases
#'                  Only used when the can not predict missing log2Intensity based on partly missing cases
#' @param minSD the minimum standard deviation used for imputing missing values (stored in replaceNA)
#'              SD will be estimated from model residuals. If not availalbe or too low, this SD will be used. 



replaceNAByModelingOnPartlyMissing <- function(subInt, minOverall = 20, minOffset = 2, minSD = 0.33, debugMessage = FALSE){
  suppressWarnings(  subInt[, numMissing := NULL])

  messageStr <- "no missing"
  if (any(is.na(subInt$log2Intensity))){
    baitsByMissing <-  subInt[is.na(log2Intensity) , .(numMissing = .N), by = bait]
    intensityByMissing <- merge(baitsByMissing, subInt[!is.na(log2Intensity) ,], by = "bait")
    subInt[baitsByMissing, numMissing := i.numMissing, on = "bait"]
    
    # initiliaze, cause we check its nonNullNess for if we have a good model or not
    rlmOut <- NULL
    
    # first case, we have baits with intensity and with different missing amounts: 1 of 3 and 2of 3, e.g.
    if (length(unique(intensityByMissing$numMissing)) > 1){
      rlmOut <- MASS::rlm(log2Intensity~numMissing, data = intensityByMissing)
      if (!coefficients(rlmOut)["numMissing"] < 0){ # doesn't behave as expected, set to NULL so we recompute below without dependence on nummissing
        if(debugMessage)message ( "reversed model, ", appendLF = FALSE)
        rlmOut <- NULL
      } else{messageStr <- "good model"}
    }
    
    # compute rlmOut just as median of all values
    if (is.null(rlmOut) & nrow(intensityByMissing) > 1){
      rlmOut <- MASS::rlm(log2Intensity~1, data = intensityByMissing)
      messageStr <- "intercept only model"
    }
    
    # predict values
    if (!is.null(rlmOut)){
      # we have good model, use that for predicted and SD of residuals for SD:
      subInt[is.na(log2Intensity), predicted := predict(rlmOut, newdata = .SD)]
      subInt[is.na(log2Intensity), SD := max(minSD,
                                             min(sd(residuals(rlmOut)), mad(residuals(rlmOut)))
                                             # median(abs(residuals(rlmOut))) * 
                                             #   1.4826 # the adustment for mad to sd, see ?mad
                                             # the actual SD:
                                             # sqrt(sum(residuals(rlmOut)^2)/
                                             #        (length(residuals(rlmOut)) - length(coefficients(rlmOut)))
                                             #)
                                             )]
    }else{
      if(nrow(intensityByMissing) > 0){ # should be exactly one row, but just in case use medians...
        messageStr <- "singleMissing"
        subInt[is.na(log2Intensity), predicted := median(intensityByMissing$log2Intensity)]
        subInt[is.na(log2Intensity), SD := minSD]
      } else{ # everyhting is completely obserfved or completely missing, so no info for the missing. 
        # fall back to the defaults
        messageStr <- "complete cases only"
        suppressWarnings( # when log2Intensity is totally empty, this will throw a warning
          subInt[is.na(log2Intensity), predicted := 
                   min(c(min(subInt$log2Intensity, na.rm = TRUE) - minOffset,
                         minOverall))]
        )
        subInt[is.na(log2Intensity), SD := minSD]
      }
    }
  }
  
  
  if(debugMessage)message (messageStr)
  subInt[, replaceNA :=  log2Intensity]
  subInt[is.na(log2Intensity), replaceNA :=  rnorm(.N, 
                                                   # ensure no predicted values drop below minOverall
                                                   ifelse(predicted < minOverall, minOverall, predicted),
                                                   SD)]
}


#' @param int.full A long data.table of apms data.  Data is normalized and must be "full". Explicitly missing.
#'                 Required columns are
#'                  `run` a unique identifer for a single APMS run. Each run has exactly 1 row per combination of bait/preyProtein
#'                  `bait` unique bait identifier
#'                  `preyProtein` unique prey identifier (usually Uniprot ID)
#'                  `log2Intensity` log2Intensity, normalized, and with explicit NA. Each run must have a value, NA or not.
#'                  `excludeGroup` a character vector that groups runs to exclude in Zmad calculation for a given row (but usually bait-wise)
#'                                 there must be a `control` group to serve as background for p value distribution
#'                  
computeZmadex <- function (int.full) {
  preyTables <- split (int.full , by = "preyProtein")
  
  print ("imputing...")
  lowImputeValue <- quantile(int.full$log2Intensity, 0.01, na.rm = TRUE) 
  suppressWarnings(
    purrr::walk(preyTables, replaceNAByModelingOnPartlyMissing, .progress = "Imputing NAs...", minOverall = lowImputeValue)
  )
  
  print ("mad median and Zmadex calculations...")
  purrr::walk(preyTables, function(dt)dt[,mad.exc := mad(   dt$replaceNA[dt$excludeGroup != excludeGroup]), by = .(bait, excludeGroup) ], .progress = "mad ...")
  purrr::walk(preyTables, function(dt)dt[,median.exc := median(dt$replaceNA[dt$excludeGroup != excludeGroup]), by = .(bait, excludeGroup)], .progress = "median ..." )
  purrr::walk(preyTables, function(dt)dt[, Zmadex  := (log2Intensity - median.exc)/mad.exc], .progress = "Zmadex ...")

  int.dt <- rbindlist(preyTables, use.names = TRUE, fill = TRUE)

  int.full[int.dt, c('replaceNA', 'mad.exc', 'median.exc', 'Zmadex') := .(replaceNA, mad.exc, median.exc, Zmadex), on = .(run, bait, preyProtein)]
  
  #convert to p values
  .z2p <- function(int.dt){
    converterFunction <- ecdf(int.dt[excludeGroup == "control", Zmadex])
    int.dt[, p := 1 - converterFunction(Zmadex)]
    minP <- int.dt[p > 0, min(p, na.rm = TRUE)]
    int.dt[p == 0, p := minP/2]
  }
  
  if ("control" %in% int.dt$excludeGroup){
    .z2p(int.full)
  }else{
    message ("'control' not found in excludeGroup. Will not compute empirical p-values.")
  }
  print ("p values...")

  return (int.full)
}



roc.dt <- function(dt, scoreColumn = "p", negativeLabels = c("decoy", "unknown"), positiveLabels = "interactor", order= 1L, labelColumn = "label"){
  totalPositive <- dt[get(labelColumn) %in% positiveLabels, .N]
  totalNegative <- dt[get(labelColumn) %in% negativeLabels, .N]
  #setorderv(dt, scoreColumn, na.last = TRUE)
  dt[, scoreRank := frankv(dt, scoreColumn, order = order)]
  rank2score <- unique(dt[, .SD , .SDcols = c(scoreColumn, "scoreRank") ])
  
  
  roc.dt <- dt[, .( numPositive = sum(get(labelColumn) %in% positiveLabels), numNegative = sum(get(labelColumn) %in% negativeLabels)), keyby = .(scoreRank)]
  roc.dt[, tpr := cumsum(numPositive)/totalPositive]
  roc.dt[, fpr := cumsum(numNegative)/totalNegative]
  
  roc.dt <- merge (roc.dt, rank2score, by  = "scoreRank")
  
  return (roc.dt)
}



roc.dt.groups <- function(dt, scoreColumn = "p", 
                          negativeLabels = c("decoy", "unknown"), 
                          positiveLabels = "interactor", 
                          order= 1L, groupByColumns = c(""),
                          maxPositives = NULL, maxNegatives = NULL){
  
  # create group-specific totalPositive and totalNegative
  if (is.null(maxPositives)){
    totalPositive <- dt[label %in% positiveLabels, .(totalPositive = .N), by =groupByColumns]
  }else{
    totalPositive <- dt[, .(totalPositive = maxPositives), by = groupByColumns]
  }
  if (is.null(maxNegatives)){
    totalNegative <- dt[label %in% negativeLabels, .(totalNegative = .N), by = groupByColumns]
  }else{
    totalNegative <- dt[, .(totalNegative = maxNegatives), by = groupByColumns]
  }
  
  
  #setorderv(dt, scoreColumn, na.last = TRUE)
  dt[, scoreRank := frankv(.SD, scoreColumn, order = order), by = groupByColumns]
  rank2score <- unique(dt[, .SD , .SDcols = c(scoreColumn, "scoreRank"), by = groupByColumns])
  
  
  roc.dt <- dt[, .( numPositive = sum(label %in% positiveLabels), numNegative = sum(label %in% negativeLabels)), keyby = c(groupByColumns, "scoreRank")]
  roc.dt[, cumulPositives := cumsum(numPositive), by = groupByColumns]
  roc.dt[, cumulNegatives := cumsum(numNegative), by = groupByColumns]
  roc.dt[totalPositive, tpr := cumulPositives/i.totalPositive, on = groupByColumns]
  roc.dt[totalNegative, fpr := cumulNegatives/i.totalNegative, on = groupByColumns]
  
  roc.dt <- merge (roc.dt, rank2score, by = c(groupByColumns, "scoreRank"))
  
  return (roc.dt)
}



# Gold Standard Interactors/Decoys ----

#' @param genes vector of gene names to limit string to
#' @param baits vector of gene names to take as baits. 
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
baitPreyDecoysFromString <- function (baits,
    preys,
                              links.path = "~/Downloads/9606.protein.physical.links.detailed.v12.0.txt.gz",
                              info.path = "~/Downloads/9606.protein.info.v12.0.txt.gz",
                              combinedScoreThreshold = 600,
                              stringDistThreshold = 5,
                              geneAliasFunction = identity){
  genes <- unique(c(preys, baits))
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
  
  decoys <- unique(distantGenes[gene1 %in% baits & gene2 %in% preys, .(gene1, gene2)])
  return (decoys)
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

baitPreyGoldStandardPairs <- function (baits,
                                       preys,
                               corum.path = "~/Downloads/corum_humanComplexes.txt", # corum 5.0 downloaded 10/31/2024
                               string.links.path = "~/Downloads/9606.protein.physical.links.detailed.v12.0.txt.gz",
                               string.info.path = "~/Downloads/9606.protein.info.v12.0.txt.gz",
                               stringCombinedScoreThreshold = 600,
                               geneAliasFunction = identity){
  
  genes <- unique(c(baits, preys))
  # remove KRT contaminants
  genes <- grep("^KRT", genes, invert = TRUE, value = TRUE)
  
  #corum
  if (!is.null(corum.path)){
    corumPairs <- unique(corumPairs(corum.path, geneAliasFunction)[, .(gene1, gene2)])
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
  
  stringPairs <- string[alias1 %in% baits & alias2 %in% preys, .(gene1 = alias1, gene2 = alias2)]
  
  combinedPairs <- rbindlist( list(corum = corumPairs,
                                   string = stringPairs),
                              idcol = 'source'
  )[
    , .(source = paste0(source, collapse = ";")),
    by = .(gene1, gene2)]
  
  return (combinedPairs[gene1 %in% baits & gene2 %in% preys])
  
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
  
  corumPairs <- unique(corumPairs[ , .(gene1 = alias1, gene2 = alias2, complex_id, complex_name)])
  return(corumPairs)
}



