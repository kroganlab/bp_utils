
# this is the original, but current thinking is that OmniPath is probably the best
# resource moving forward.  See function loadKinaseDataOmniPath May 2021
loadKinaseData <- function (kinaseDataFile = NULL){
  if (is.null(kinaseDataFile)){
    message ("path to kinaseDataFile is required")
    message ("consider downloading from https://www.biorxiv.org/content/biorxiv/early/2019/11/06/822668.1/DC1/embed/media-1.zip?download=true")
    message ("Should also be available in bp_utils/data/kinaseSiteList_BachmanGyoriSorger2019.csv.gz")
    stop()
  }
  
  #download file from https://www.biorxiv.org/content/10.1101/822668v3.supplementary-material
  # https://www.biorxiv.org/content/biorxiv/early/2019/11/06/822668.1/DC1/embed/media-1.zip?download=true
  #  QB3-901LVCF-LT:flu ben$ cp ~/Downloads/media-1/export.csv ./data/kinaseSiteList_BachmanGyoriSorger2019.csv
  
  kinaseSubstrates <- fread (kinaseDataFile)
  kinaseSubstrates <- kinaseSubstrates[CTRL_IS_KINASE == TRUE]
  
  #for CTRL_NS == FPLX the CTRL_ID looks like a kinase name but the CTRL_GENE_NAME field is blank
  # for now copy them whole-hog to CTRL_GENE_NAME
  # it will treat all (except for SRC, the only overlap) as new kinases for now.
  kinaseSubstrates[CTRL_NS == "FPLX" & CTRL_GENE_NAME == "", CTRL_GENE_NAME := CTRL_ID]
  
  #remove kinase self-phosphorylations
  #return (kinaseSubstrates[TARGET_GENE_NAME != CTRL_GENE_NAME])
  return (kinaseSubstrates)
}


loadKinaseDataFromKSEAFile <- function (ksDataFile = "./data/PSP&NetworKIN_Kinase_Substrate_Dataset_July2016.csv.gz", networKIN.minScore = 5){
  ksdata <- fread(ksDataFile)
  ksdata <- ksdata[networkin_score > networKIN.minScore]
  
  # for consistency with loadKinaseData I need columns: "CTRL_GENE_NAME", "TARGET_GENE_NAME", "TARGET_RES", "TARGET_POS"
  
  setnames(ksdata, old = c("GENE", "SUB_GENE"), new = c("CTRL_GENE_NAME", "TARGET_GENE_NAME") )
  ksdata[, TARGET_RES := substr(SUB_MOD_RSD, 1, 1)]
  ksdata[, TARGET_POS := as.integer (gsub("^.", "", SUB_MOD_RSD))]
  return (ksdata[])
}




# ... arguments are passed to omnipath function. 
#     #organism (9606 = human, default; 10116 rat; or 10090 mouse)
#     resources to select source of enz-sub information see get_ptms_resources for options


loadKinaseDataOmniPath <- function(species = "HUMAN", removeNonKinases = TRUE, fixNames = TRUE, ...){
  species <- toupper(species)
  organismID <- c(HUMAN = 9606, RAT = 10116, MOUSE = 10090)[species]
  if(is.na(organismID))
    stop(species, " not a recognized organism")
  enzsub <- OmnipathR::import_omnipath_enzsub(organism = organismID,...)
  setDT (enzsub)
  # for consistency with loadKinaseData I need columns: "CTRL_GENE_NAME", "TARGET_GENE_NAME", "TARGET_RES", "TARGET_POS", "TARGET_UP_ID"
  
  if (any(grepl("ProtMapper", enzsub$sources))){
    message ("OmniPath includes ProtMapper data which has many non-kinases.")
    
    if(removeNonKinases == TRUE) {
      message ("Removing all enzymes that are not kinases according to org.Hs.eg.db GO:0016301")
      # get GO kinase list
      
      if (species == "HUMAN"){
        library(org.Hs.eg.db)
        kinaseTable <- AnnotationDbi::select (org.Hs.eg.db, get(c("GO:0016301"), org.Hs.egGO2ALLEGS), c("ENTREZID", "GENENAME", "SYMBOL")) 
      } else if (species == "MOUSE"){
        library(org.Mm.eg.db)
        kinaseTable <- AnnotationDbi::select (org.Mm.eg.db, get(c("GO:0016301"), org.Mm.egGO2ALLEGS), c("ENTREZID", "GENENAME", "SYMBOL")) 
      } else{
        stop ("I don't yet know how to find kinases for ", species)
      }
      
      kinases <- unique(kinaseTable$SYMBOL)
      #label and filter into a new table
      enzsub[, is_kinase := enzyme_genesymbol %in% kinases]
      enzsub <- enzsub[is_kinase == TRUE]
    }
  }
  
  if (fixNames){
    setnames(enzsub, 
             old = c("enzyme_genesymbol", "substrate_genesymbol", "residue_type", "residue_offset", "substrate"), 
             new = c("CTRL_GENE_NAME", "TARGET_GENE_NAME", "TARGET_RES", "TARGET_POS", "TARGET_UP_ID"))
  }
  
  return (enzsub[modification == "phosphorylation"])
}

kinaseDataToGMT <- function(kd){
  kd[, .(gene = sprintf ("%s_%s%d", TARGET_GENE_NAME, TARGET_RES, TARGET_POS),
         kinase = CTRL_GENE_NAME,
         protein = sprintf ("%s_%s%d", TARGET_UP_ID, TARGET_RES, TARGET_POS))]
}



# Compute enrichment scores based on gene set enrichment analysis.
# uses data table structure that is built and returned by kinaseActivity
kinaseActivity.sea <- 
  function(kinaseMapped, nperm="not_used", gseaWeightParam = 1, nproc=1, fullData = NULL, center = FALSE, ...){
  stopifnot(length(unique(kinaseMapped$Label)) < 2) # we should only be computing this over a single label
  
  #collapse log2FC down to a single value per phSiteCombo
  if(is.null(fullData)){
  log2FC.table <-  kinaseMapped[is.finite(log2FC) & representative == TRUE, 
                                   .(log2FC = unique(log2FC)), by = .(phSiteCombo)]
  } else{
    stopifnot(length(unique(fullData$Label)) < 2) # we should only be computing this over a single label
    log2FC.table <- fullData[is.finite(log2FC) & representative == TRUE, 
                             .(log2FC = unique(log2FC)), by = .(phSiteCombo)]
  }
  
  #convert to a named vector
  log2FC<-log2FC.table$log2FC
  names(log2FC)<-log2FC.table$phSiteCombo
  
  if (center)
    log2FC <- log2FC - median(log2FC)
  
  setsTable <- kinaseMapped[representative == TRUE, .(phSiteCombo=unique(phSiteCombo)), by = .(CTRL_GENE_NAME)]
  sets <- split(setsTable$phSiteCombo,
                setsTable$CTRL_GENE_NAME)
  
  seaRes <- fgsea::fgsea(pathways = sets, stats = log2FC,  gseaParam=gseaWeightParam, nproc = nproc, ...)  #nperm = nperm
  seaRes[, sigScore := -log10(pval) * ifelse(ES < 0, -1, 1) ]
  return (seaRes)
}

# Demonstrating expected usage of function kinaseActivity
#' @param results a data.table
#' @param kinaseData a data.table, output of loadKinaseDataOmniPath etc
kinaseActivityOnResultsFile <- function(results, kinaseData){
  results[,gene := multiUniprotSites2multiGeneSites(Protein)]
  singleSiteResults <- prepare_AMSS_ResultsFile(results, column = "gene")
  
  labels <- unique(singleSiteResults$Label)
  
  kinActList <- lapply (labels, FUN=function(lab){kinaseActivity(singleSiteResults[Label == lab & representative==TRUE],
                                                                 plots = FALSE,
                                                                 kinaseData = kinaseData)})
  names(kinActList) <- labels
  
  kinActFull.scores <- rbindlist(lapply(kinActList, FUN = function(x)x$scores), idcol="Label")
  kinActFull.mapped <- rbindlist(lapply(kinActList, FUN = function(x)x$kinaseMapped)) # Label is already in these tables
  list(kinActFull.scores, kinActFull.mapped)
}




# Meant to be called once per different contrast or "Label" in results file.
# see function above kinaseActivityOnResultsFile for a code example

kinaseActivity <- function (log2FCData, kinaseData=NULL, kinaseDataFile = "./data/kinaseSiteList_BachmanGyoriSorger2019.csv.gz",
                            plots=TRUE, outlierLog2FC = 10, requireAAMatch = TRUE, uniprot=FALSE, do.sea = FALSE,
                            sea.center = FALSE, sea.weightParam = 1, sea.nproc = 1, subsetColumn = NULL){
  
  
  # deal with columns ... this is ugly ... can be improved
  requiredColumns <- c("log2FC", "pos", "phSiteCombo")
  if (requireAAMatch) requiredColumns <- c(requiredColumns, "aa")
  if (uniprot){
    requiredColumns <- c(requiredColumns, "uniprot")
  }else{
    requiredColumns <- c(requiredColumns, "Gene_Name")
  }

  missingColumns <- setdiff(requiredColumns, colnames(log2FCData))
  if (length(missingColumns) > 0) stop (length(missingColumns), " missing column(s) ", paste0(missingColumns, collapse=", "), " in log2FCData  table")
  
  if (uniprot){
    mergeXCols <- c("uniprot")
    mergeYCols <- c("TARGET_UP_ID")
  }else{
    mergeXCols <- c("Gene_Name")
    mergeYCols <- c("TARGET_GENE_NAME")
  }
  if(requireAAMatch){
    mergeXCols <- c(mergeXCols, "aa")
    mergeYCols <- c(mergeYCols, "TARGET_RES")
  }
  mergeXCols <- c(mergeXCols, "pos")
  mergeYCols <- c(mergeYCols, "TARGET_POS")
  
  # done with columns
  
  if ("representative" %in% colnames(log2FCData) & any(log2FCData$representative == FALSE)){
    warning ("data.table includes a 'representative' column, but some are set to FALSE. Did you forget to pre-filter?")
  }
    
  if (is.null(kinaseData)){
    kinaseData <- loadKinaseData (kinaseDataFile)
    #message ("Loading kinase data from OmniPath")
    #kinaseData <- loadKinaseDataOmniPath()
  }

  # redundancy check the kinaseData based on mergeYCols
  nrowFull <- nrow(kinaseData)
  nrowUnique <- nrow(unique(kinaseData[,c("CTRL_GENE_NAME", mergeYCols), with = FALSE]))
  if(nrowUnique < nrowFull){
    message (nrowFull - nrowUnique, " redundant rows will be dropped from kinase knowledge table; selecting first of each redundant set")
    kinaseData <- kinaseData[, .SD[1], by  = c("CTRL_GENE_NAME", mergeYCols)]
  }
  
    
  # clean up data
  ## remove NA and infinite
  cleanData <- log2FCData[!is.na(log2FC) & is.finite(log2FC)]
  if (plots) boxplot(cleanData$log2FC)
  ## outliers
  if (!is.na(outlierLog2FC)){
    cleanData <- cleanData[abs(log2FC) < abs(outlierLog2FC)]
    if (plots) boxplot(cleanData$log2FC)
  }
  

 
  kinaseMapped <- merge (cleanData, kinaseData, by.x=mergeXCols, by.y=mergeYCols)
  
  cat (sprintf("%d sites quantified mapped to %d rows in kinase-target data file results in %d kinase-site-log2FC relationships from %d phSiteCombo and %d kinases\n",
               nrow(cleanData), nrow(kinaseData), nrow(kinaseMapped), kinaseMapped[,length(unique(phSiteCombo))], kinaseMapped[, length(unique(CTRL_GENE_NAME))]))


    
  if (uniprot) kinaseMapped[, Gene_Name := uniprot]
  if (! "aa" %in% colnames(kinaseMapped)) kinaseMapped[,aa := ""]
  
  if(plots){
    plot (density(cleanData[,log2FC]))
    lines (density(kinaseMapped[,log2FC]), col="red")
    #legend ...
  }
  
  # the denominator is SD of the full dataset, whether it maps to a kinase or not:
  fullSD <- sd(cleanData[,log2FC])
  #each kinase mean will be compared to global mean
  fullMean <- mean(cleanData[,log2FC])
  
  subsetSelector <-TRUE  #select all
  if (!is.null(subsetColumn)){
    subsetSelector <- kinaseMapped[[subsetColumn]] == TRUE
  }else{
    subsetSelector <- rep(TRUE, nrow(kinaseMapped))
  }
  
  cat (sprintf("Computing using %d rows in site-kinase-log2FC table\n",sum(subsetSelector)))
  
  
  scores <- siteEnrichZScores(kinaseMapped[subsetSelector], fullMean, fullSD)
  
   
  if (do.sea == TRUE){
    if (!is.null(subsetColumn)) message ("subsetColumn does not work with SEA. SEA scores will not be subsetted")
    seaScores <- kinaseActivity.sea(kinaseMapped, fullData = cleanData, center = sea.center, gseaWeightParam = sea.weightParam,  nproc = sea.nproc) #nperm = 10000,
    setnames(seaScores, old = c("pval", "padj"), new = paste0(c("pval", "padj"), ".sea"))
    scores <- merge (scores, seaScores, by.x = c("CTRL_GENE_NAME"), by.y = c("pathway"),suffixes = c("", ".sea"), all = TRUE)
  }
  
  setorder(scores, pValue)
  
  if (plots){
    hist (scores$pValue, breaks=100)
    p <- ggplot(scores, aes(x=Z, y = -log10(pValue), col = sqrt(N))) + geom_point()
    print (p)
    p <- ggplot(scores, aes(x=Z, y = -log10(fdr.BH), col = sqrt(N))) + geom_point()
    print (p)
    p <- ggplot(scores, aes(x=meanLog2FC, y = -log10(pValue), col = sqrt(N))) + geom_point()
    print (p)
    p <- ggplot(scores, aes(x=meanLog2FC, y = -log10(fdr.BH), col = sqrt(N))) + geom_point() #+ scale_color_steps(low="white", high="navy", breaks=c(0,1,3,10,50,100,200)) + theme_classic()
    print (p)
    
    p <- BarplotKinaseActivities(scores, kinaseMapped, max_pValue = 0.05)
    print (p)
    p <- BarplotKinaseActivities(scores, kinaseMapped, max_pValue = 1.0)
    print (p)
    
  }
  return (list(scores = scores, kinaseMapped = kinaseMapped))
}

# The calculation of scores once all the cleaning of data has been done. Works on an already merged results/annotation file
#
#' @param mergedResults a data.table, with log2FC, phSiteCombo, and CTRL_GENE_NAME or with other column specified by termName
#' @param termName character that is name of column in `mergedResults` to enrich within. (kinase, ubiquitin ligase, etc..) 
#' @param siteColumn column name in mergedResults to collect in a sites column in the output
siteEnrichZScores <- function(mergedResults, fullMean = NULL, fullSD = NULL, termName  = "CTRL_GENE_NAME", siteColumn = "singleSite"){
  if(is.null(fullMean) | is.null(fullSD)) stop("Require calculation of background 'fullMean' and 'fullSD' ahead of time")
  if ("termName" %in% colnames(mergedResults))
    warning("termName used as column name. termName column will be used regardless of any value that the termName argument is set to")
  # summarize per kinase, Z score based on standard error of the mean log2FC per kinase
  scores <- mergedResults[,.(Z = (unique(.SD)[,mean(log2FC)] - fullMean)*sqrt(length(unique(phSiteCombo)))/fullSD,
                            N= length(unique(phSiteCombo)),#.N,  # don't double count when two sites are in the same peptide
                            meanLog2FC  = mean(log2FC),
                            sites = paste0(sort(unique(.SD[[siteColumn]])), collapse = ";") ),  # 'unique' because some sites will occur multiple times in different phospho-combos
                         by = termName,
                         .SDcols = c("phSiteCombo", "log2FC", siteColumn)]

  # p values.  2*pnorm... makes it two-tailed. This is different from KSEApp, but I think it is more appropriate here where we include both positive and negative effects
  scores[,pValue := 2*pnorm(-abs(Z), lower.tail= TRUE)]
  scores[, sigScore := -log10(pValue) * ifelse(Z < 0, -1, 1) ]
  scores[, fdr.BH := p.adjust(pValue, method = "fdr")]  #method=fdr is alias for "BH" Benjamini & Hochberg
  scores[,c("bgMean", "bgSD") := .(fullMean, fullSD)]
  return (scores[])
}



BarplotKinaseActivities <- function(scores, kinaseMapped, 
                                    max_pValue = 1.0, max_fdr = 1.0, min_N = 2,
                                    sigKinases = NULL, reverse = FALSE,
                                    useMonoFont = FALSE, useViolin = FALSE,
                                    useSEA = FALSE, reorder = TRUE,
                                    ncol = 3, labelPoints = FALSE){
  by.col <- c("CTRL_GENE_NAME")
  
  if ("Label" %in% colnames(kinaseMapped)){
    #expected in kinaseMapped, even with just 1 label
    if ("Label" %in% colnames(scores)){
      # if  we get here, make sure we are joining by Label
      by.col <- c(by.col, "Label")
    } else if(length(unique(kinaseMapped$Label)) != 1){
      # we're going to ignore label, so make sure data is only from one Label
      warning("Multiple Labels detected in kinaseMapped, but no Label information in scores.  You might be combining log2FC from multiple contrasts inappropriately: ", paste0(unique(kinaseMapped$Label), collpase = ", "))
    }
  }

  b <- merge (scores, kinaseMapped, by = by.col)
  
  if (is.null(sigKinases)){
    if (useSEA){
      sigKinases <-  unique(scores[pval.sea< max_pValue & padj.sea < max_fdr & size >= min_N]$CTRL_GENE_NAME)
    }else{
      sigKinases <-  unique(scores[pValue< max_pValue & fdr.BH < max_fdr & N >= min_N]$CTRL_GENE_NAME)
    }
  }
  
  if (length(sigKinases) == 0){
    return ("No significant kinases")
  }
  
  if (reorder == TRUE){
    scoreSummary <- scores[,.(meanZ = mean(Z, na.rm = TRUE)), by = CTRL_GENE_NAME]
    setorder(scoreSummary, meanZ)
    if (useSEA){
      scoreSummary <- scores[,.(meanNES = mean(NES, na.rm = TRUE)), by = CTRL_GENE_NAME]
      setorder(scoreSummary, meanNES)
    }
    kinasesSorted <- scoreSummary[,CTRL_GENE_NAME]
  }else kinasesSorted <- sigKinases
  
  b[,CTRL_GENE_NAME := factor(CTRL_GENE_NAME, levels = kinasesSorted)]
  
  if (reverse){
    p <- ggplot (b[CTRL_GENE_NAME %in% sigKinases,], aes(x=log2FC, y = Label, fill = sigScore, col = sigScore, label = phSiteCombo)) + 
      facet_wrap(facets = ~CTRL_GENE_NAME, ncol = ncol)
  }else{
    p <- ggplot (b[CTRL_GENE_NAME %in% sigKinases,], aes(x=log2FC, y = CTRL_GENE_NAME, fill = sigScore, col = sigScore, label = phSiteCombo)) + 
      facet_wrap(facets = ~Label, ncol = ncol)
  }
  if (useSEA){
    p <- p + aes(fill = sigScore.sea, col = sigScore.sea)
  }
  p <- p +
    geom_vline(xintercept=0.0, lty="dotted", col = "black") + 
    geom_jitter(width=0.0, height=0.1, col="black", alpha=0.5)
  if (useViolin){
    p <- p + geom_violin(scale = "area",  alpha=0.7) 
  } else{
    p <-  p + geom_boxplot( varwidth=FALSE, alpha=0.7,outlier.shape = NA)
  }
  p <- p + 
    scale_color_gradient2(low = "blue", mid= "gray", high="red", midpoint=0.0) + 
    scale_fill_gradient2(low = "blue", mid= "gray", high="red", midpoint=0.0) + 
    theme_classic() 
  
  if (useMonoFont) p <- p + theme(axis.text.y = element_text( size = 10, family = "mono"))
  

  if (labelPoints){
    p <- p + ggrepel::geom_text_repel(data = b[CTRL_GENE_NAME %in% sigKinases & abs(sigScore) > 1.5,
                                               .SD[which.max(abs(log2FC))], by = .(Label, CTRL_GENE_NAME)], size = 2, min.segment.length = 0 )
  }
  
  

  return (p)
}


# chooseRepsInSiteLog2FCData
#
# does the job of dealing with multiple site-combos per site.
# nothing is actually removed, but those rows that should be kept are tagged in the representative column
#  required columns in log2FCData:
#     aaPos   -  a single site (usually in string format like S35, T109 or Y99, but any consistent format is allowed)
#     Gene_Name - an identifier for the protein. The combination of Gene_Name and aaPos uniquely identifies the site.
#                 This defines where redundancy is removed.
#     contrast_label - the treatment groups inside which the rows for a single Gene_Name/aaPos are grouped
#     pvalue
#     log2FC
#
#  choices: per single site, takes the Protein (aka site combo) with the lowest pvalue as the representative, favoring non-infinites.
#           if there is a tie, then log2FC is the tie breaker
#           if there is a conflict in direction for the infinites (+Inf and -Inf), no representative is chosen
#    Note that for most kinase analysis to-date, infinites are ignored
#
#  WARNING: This modifies the data in place.  It will change the order and add four columns:
#       is_infinite
#       absLog2FC
#       infiniteConflict TRUE/FALSE does this site with the contrast_label have infinite fold chnages in opposite directions
#       representative   TRUE/FALSE should this row be selected to represent the site in the contrast
#
#    columnRename = c(Gene_Name = "singleProtein")
chooseRepsInSiteLog2FCData <- function(log2FCData, columnRename = NULL){
  reqColumns <- c("aaPos", "Gene_Name", "contrast_label")
  
  if  (!is.null(columnRename)){
    if (any(!names(columnRename) %in% reqColumns)){
      warning("Asking to rename a column to a non required column: ", paste0(columnRename, collapse = ", "), "\n",
              "Required columns: ", paste0(reqColumns, collapse = ", "))
      message (columnRename)
    }
    setnames(log2FCData, old = columnRename, new = names(columnRename), skip_absent = TRUE)
  }
  
  missingColumns <- setdiff(reqColumns, colnames(log2FCData))
  if(length(missingColumns) > 0){
    stop("chooseRepsInSiteLog2FCData: required columns are missing: ", paste0(missingColumns, collapse = ", "), "\n",
         "Make use of columnRename to rename an existing column")
  }
    
  #strategy is to sort by relevant data, then choose the top-ranked row after grouping by site
  # to facilitate this sorting, I have to produce a few columns
  message ("Adding column is_infinite")
  log2FCData[, is_infinite := is.infinite(log2FC)]
  message ("Adding column absLog2FC")
  log2FCData[,absLog2FC := abs(log2FC)]
  # mark those infinites that have conflicting directions.
  message ("Adding column infiniteConflict")
  log2FCData[is.infinite(log2FC),infiniteConflict := Inf %in% log2FC & -Inf %in% log2FC, by = .(contrast_label, Gene_Name, aaPos)]
  
  #order
  message ("Re-ordering")
  setorder(log2FCData, contrast_label, Gene_Name, aaPos,  # the grouping variables, sorted here for easy viewing
           is_infinite, pvalue, -absLog2FC)  # the ranking variables
  #initialize to be safe
  message ("Adding column representative")
  log2FCData[,representative := FALSE]
  # mark the top-ranked per site
  message ("choosing representative per contrast_label, Gene_Name, aaPos")
  log2FCData[, representative := .I == min(.I), by = .(contrast_label, Gene_Name, aaPos)]
  # unmark those infinite representatives that have a conflict in direction
  message ("un-choosing those infinites with conflict, ", log2FCData[representative == TRUE & infiniteConflict ==TRUE, .N], " cases")
  log2FCData[representative == TRUE & infiniteConflict ==TRUE, representative := FALSE]
  
  return(log2FCData)
}


test__chooseRepsInSiteLog2FCData <- function(){
  test <-
    structure(list(contrast_label = c("H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", 
                                      "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", 
                                      "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", 
                                      "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", 
                                      "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04", 
                                      "H5N1_D04-MOCK_D04", "H5N1_D04-MOCK_D04"), 
                   Gene_Name = c("SPTBN1", 
                                 "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", 
                                 "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", "SPTBN1", 
                                 "SPTBN1"), 
                   aaPos = c("S2165", "S2165", "S2165", "S2165", "S2165", 
                             "S2165", "S2165", "S2165", "S2165", "S2165", "S2165", "S2165", 
                             "S2165", "S2165", "S2165", "S2165"), 
                   phSiteCombo = c("Q62261 _ 2159_2160_2164_2168", 
                                   "Q62261 _ 2163_2164_2168", "Q62261 _ 2163_2164", "Q62261 _ 2159_2160_2164", 
                                   "Q62261 _ 2160_2164_2168", "Q62261 _ 2160_2163_2164_2168", "Q62261 _ 2164_2168", 
                                   "Q62261 _ 2160_2163_2164_2170", "Q62261 _ 2158_2163_2164_2168", 
                                   "Q62261 _ 2159_2164", "Q62261 _ 2159_2163_2164_2168", "Q62261 _ 2160_2164", 
                                   "Q62261 _ 2159_2164_2168", "Q62261 _ 2158_2159_2164_2168", "Q62261 _ 2158_2164", 
                                   "Q62261 _ 2158_2164_2168"), 
                   log2FC = c(2.05984381332201, 1.45147784388535, 
                              0.85490651326672, 1.08966755548881, 1.72535298829343, 1.55535703678434, 
                              0.969295916735991, 0.658954471460672, 0.567850158492438, 0.535303035883381, 
                              0.354740581106709, -0.0721370957516222, -Inf, Inf, -Inf, -Inf
                   ), 
                   pvalue = c(8.91229718349429e-05, 0.000701537373813954, 0.00289882470155267, 
                              0.00417113961607929, 0.0169145583382315, 0.0335536780777239, 
                              0.113939381684379, 0.14389315172982, 0.189011868032981, 0.201585918040218, 
                              0.620544071882386, 0.859791455469243, NA, NA, NA, NA)),
              row.names = c(NA, 
                            -16L), class = "data.frame")
  setDT(test)
  chooseRepsInSiteLog2FCData (test)
  View(test)
}




# Next few functions prepare an artMS/Msstats PH results file for kinase activity mapping
# deals with multiple proteins per row, multiple sites per protein

# this handles splitting the site off of a protein, sepaarated by a single underscore
# it should handle the possibility that there are other underscores in a protein name
# it will not handle the case where there are two sites  (other than simply splitting off the last one)
# example P789ASDF_S123_T125 will not split off both sites
splitProtein_Site <- function (protein_site){
  fullSplit <- strsplit (protein_site, "_")
  rejoinProteinSplit <- lapply (fullSplit, 
                                FUN = function(ss)
                                {
                                  lastI <- length(ss)
                                  return(c(paste(ss[1:(lastI-1)], collapse="_"), ss[lastI]))
                                })
  return (data.table::transpose(rejoinProteinSplit))
}

singleSiteMapping <- function(phData, proteinColumn = "Protein"){
  singleSiteMapping <- data.table (Protein = unique(phData[[proteinColumn]]))
  singleSiteMapping <- singleSiteMapping[,.(singleSite = unlist(strsplit(Protein, split=";"))), by = Protein]
  singleSiteMapping[,c("singleProtein", "aaAndPos"):=splitProtein_Site(singleSite)]
  singleSiteMapping[,c("aa", "pos") := .(substr(aaAndPos, 1,1), as.integer(substr(aaAndPos, 2, stringr::str_length(aaAndPos))))]
  return(singleSiteMapping[])
}


expandProteinToSingleSites <- function(data, proteinColumn = "Protein"){
  singleSiteMap <- singleSiteMapping (data, proteinColumn)
  return (merge (data, singleSiteMap, by.x=proteinColumn, by.y = "Protein", allow.cartesian=TRUE))
}


prepare_AMSS_ResultsFile <- function(resultsDT, column = "Protein"){
  
  expanded2SingleSites <- expandProteinToSingleSites(resultsDT, proteinColumn = column)
  

  expanded2SingleSites[, aaPos := aaAndPos]
  if (!"contrast_label" %in% colnames(expanded2SingleSites)){
    expanded2SingleSites[, contrast_label := Label]
  }
  expanded2SingleSites[,phSiteCombo := expanded2SingleSites[[column]] ]
  
  chooseRepsInSiteLog2FCData(expanded2SingleSites, columnRename = c(Gene_Name = "singleProtein"))
  
  return (expanded2SingleSites[])
}



# Other Kinase related info ----


#strong functional scores from Beltrao paper
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7100915/
#Supplemental table s3 has the scores per uniprot and site
LoadFunctionalScores_PHSite_Beltrao <- function(path = "../../bp_utils/data/EMS84831-supplement-Table_S3.xlsx"){
  sfs <- openxlsx::read.xlsx (path)
  setDT(sfs)
  return(sfs[])
}

#https://www.phosphosite.org/downloads/Regulatory_sites.gz
# Kinases only!
LoadRegulatoryKinasePhosphoSites <- function(path = "../../bp_utils/data/Regulatory_sites.gz",
                                             species = "human"){
  reg <- fread(path, skip = 2, fill = TRUE)
  reg <- reg[ORGANISM == tolower(species) & 
               grepl ("kinase", PROT_TYPE, ignore.case = TRUE) &
               grepl ("-p$", MOD_RSD)]
  
  
  reg[, activating := grepl("activity, induced", ON_FUNCTION)]
  reg[, inhibiting := grepl("activity, inhibited", ON_FUNCTION)]
  reg[, interaction.regulating := grepl("molecular association, regulation", ON_FUNCTION)]
  
  reg[, uniprot_site := paste0(ACC_ID, "_", tstrsplit(MOD_RSD, "-")[[1]])]
  reg[, gene_site := paste0(GENE, "_", tstrsplit(MOD_RSD, "-")[[1]])]
  
  
  return (reg[])
}

# Alll phospho-proteins
LoadPhosphoSiteRegulationDirection <- function(path = "../../bp_utils/data/Regulatory_sites.gz",
                                               species = "human"){
  reg <- fread(path, skip = 2, fill = TRUE)
  reg <- reg[ORGANISM == tolower(species) & 
               grepl ("-p$", MOD_RSD)]
  
  
  reg[, activating := grepl("activity, induced", ON_FUNCTION)]
  reg[, inhibiting := grepl("activity, inhibited", ON_FUNCTION)]
  reg[, interaction.regulating := grepl("molecular association, regulation", ON_FUNCTION)]
  
  reg[, uniprot_site := paste0(ACC_ID, "_", tstrsplit(MOD_RSD, "-")[[1]])]
  reg[, gene_site := paste0(GENE, "_", tstrsplit(MOD_RSD, "-")[[1]])]
  
  dir <- reg[activating == !inhibiting, .(uniprot_site, gene_site, direction = ifelse(activating, 1, -1))]
  
  
  return (dir[])
}




