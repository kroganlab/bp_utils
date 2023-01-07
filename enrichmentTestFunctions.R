library (data.table)


#' @param groupTable data.table with columns for gene symbols and group designtations
#' @param geneColumn character name of column which has gene symbols (or uniprot depending on gmt)
#' @param groupColumns character vector (>1 allowed) of factor-like columns that defines groups
#' @param term2gene.gmt the term2gene table.  First two columns (regardless of names) are used as term, gene
#' @param universe a character vector of gene IDs that define the background universe
#' @param numProcessors integer. Values > 1 enable multiprocessing

enricherOnGroups <- function(groupTable, geneColumn = "gene", groupColumns = c("group"),
                           term2gene.gmt = NULL, universe = NULL, numProcessors = 1,
                           ...){
  
  if (is.null(term2gene.gmt)){
    stop ("A term2gene.gmt is required. 1st column term, 2nd column gene IDs (or uniprots etc)")
  }
  if (is.null(universe)){
    warning("No universe chosen. Using the full set of genes in term2gene.gmt.  This is likely inappropriate")
    universe <- unique(term2gene.gmt[,2])
  }
  
  subTables <- split (groupTable, by = groupColumns, flatten = TRUE, drop = TRUE)
  
  message ("Computing enrichments on ", length(subTables), " groups defined by ", paste0(groupColumns, collapse = ","))
  enrichList <-  pbapply::pblapply(subTables, 
                          function(subDT){
                            setDT(as.data.frame(clusterProfiler::enricher(unique(subDT[[geneColumn]]),
                                                         pAdjustMethod="fdr",
                                                         universe=universe,
                                                         TERM2GENE =  term2gene.gmt,
                                                         qvalueCutoff = 1.1,
                                                         pvalueCutoff = 1.1,
                                                         ...)))
                          },
                          cl=numProcessors)
  
  enrichTable <-  rbindlist(enrichList, idcol= paste0(groupColumns, collapse = "."))
  return (enrichTable)
}

#' Compute an odds ratio using the text columns GeneRatio and BgRatio from the typical clusterProfiler output
#' The odds ratio is the ratio of odds between the selected and not selected groups for matching the term or not.
#' Thus we have to compute three numbers from GeneRatio and BgRatio:
#'     selectNoTerm:      those in our selection that don't match the term
#'     notSelectWiTerm:   those not selected that match the term
#'     notSelectNoTerm:   those not selected that don't match the term.
#' @param GeneRatio text column with entries in format 12/200 where 200 is the number of selected genes and 12 the number of those that match the term
#' @param BgRatio  text column with entries in format 120/2000 where 2000 is the number genes in the universe and 120 the number of those that match the term 
#' 
#' @return  a numeric vector of odds ratios, >1 implies more term matches observed in selected set than expected
#' 
OddsRatioFromEnricherRatios <- function (GeneRatio, BgRatio){
  selectParts <- tstrsplit(GeneRatio, "/")
  selectWiTerm <- as.integer(selectParts[[1]])
  selectTotal <- as.integer(selectParts[[2]])
  selectNoTerm <- selectTotal - selectWiTerm  #
  
  bgParts <- tstrsplit(BgRatio, "/")
  bgWiTerm <- as.integer(bgParts[[1]])
  bgTotal <- as.integer (bgParts[[2]])
  bgNoTerm <- bgTotal - bgWiTerm
  
  notSelectWiTerm <- bgWiTerm - selectWiTerm  #
  notSelectNoTerm <- bgNoTerm - selectNoTerm  #
  
  oddsRatio <- (selectWiTerm/selectNoTerm)/
               (notSelectWiTerm/notSelectNoTerm)
  
  return (oddsRatio)
}


simplifyEnrichBySimilarUniverseMembership.general <- function (enrichResultsTable, gmt, groupColumn=NULL, 
                                                       cutHeight = 0.99, broadest=TRUE, max_pAdjust = 0.01,
                                                       termColumn = "pathway", pvalueColumn = "padj"){
  if (length(unique(enrichResultsTable[[termColumn]])) < 2){
    message ("Nothing to simplify")
    return (list (enrichResultsTable, data.frame()))
  }
  setDT(enrichResultsTable); setDT(gmt)
  
  ##Select Significant Terms
  target_overrep_sig <- enrichResultsTable[enrichResultsTable[[pvalueColumn]] < max(max_pAdjust,0.05),]#qvalue < 0.01]

  ##Prepare Significant GO Term Jaccard Similarity Matrix
  sig_go_terms <- unique(target_overrep_sig[[termColumn]])
  
  if (!is.null(groupColumn)){
    message ("Computing universal gene overlap between ", length(sig_go_terms), " significant terms from ", length(unique(enrichResultsTable[[groupColumn]])), " ", groupColumn, "(s)")
  }else{
    message ("Computing universal gene overlap between ", length(sig_go_terms), " significant terms")
  }

  #cluster and divide the gmt based
  gmt.subset <- gmt[ont %in% sig_go_terms, .(ont=factor(ont), gene=factor(gene))]
  termByGeneMat <-  Matrix::sparseMatrix(as.integer(gmt.subset$ont), as.integer(gmt.subset$gene), 
                                         dimnames=list(levels(gmt.subset$ont), 
                                                       levels(gmt.subset$gene)))
  
  #go_dist_mat <- dist(termByGeneMat, method="binary")
  # maybe this is faster?  I haven't clocked it
  go_dist_mat <- parallelDist::parDist(as.matrix(termByGeneMat), method="binary")
  
  hc <- hclust(go_dist_mat)
  
  clusters <- cutree(hc, h=cutHeight)
  clusters <- data.table (cluster = as.numeric(clusters), ID = attributes(clusters)$names )
  
  message ("GO terms clustered into ", max(clusters$cluster), " clusters")
  
  ## go gene set lengths
  gmt.setLengths <- gmt[,.(setSize = length(unique(gene))), by = ont]
  clusters <- merge (clusters, gmt.setLengths, by.x="ID", by.y="ont")
  
  
  ## local gene set lengths  -- this will only count those genes that make it into an enriched list
  # specific to enricher output...

  # id.gene.long <- enrichResultsTable[,.(ID, gene = unlist(strsplit(geneID, "/"))), by = seq_len(nrow(enrichResultsTable))]
  # genesPerTerm <- id.gene.long[,.(localSetLength=length(unique(gene))), by = ID]
  # clusters <- merge (clusters, genesPerTerm, by.x="ID", by.y="ID")
  #setorder(clusters, -localSetLength)  # the tie breaker
  
  clusterInfo <- merge (enrichResultsTable, clusters, by.x = termColumn, by.y = "ID")
  
  clusterInfo[,maxSet := max (setSize), by = c("cluster", groupColumn)]

  if (broadest){
    winners <- clusterInfo[clusterInfo[[pvalueColumn]] < max_pAdjust,.SD[which.max(setSize),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners[[termColumn]])), " representative terms choosing the BROADEST significant term per term-cluster per ", groupColumn)
  }else{
    # setorder(clusterInfo, cols = pvalueColumn)
    # winners <- clusterInfo[clusterInfo[[pvalueColumn]] < max_pAdjust,.SD[1],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    winners <- clusterInfo[clusterInfo[[pvalueColumn]] < max_pAdjust,.SD[which.min(.SD[[pvalueColumn]]),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners[[termColumn]])), " representative  terms choosing the MOST significant term per term-cluster per ", groupColumn)
  }
  result <- enrichResultsTable[enrichResultsTable[[termColumn]] %in% winners[[termColumn]],]
  result[clusterInfo, cluster.id := cluster, on = termColumn]
  list(simplified = result, clusterInfo = clusterInfo)
}





simplifyEnrichBySimilarUniverseMembership <- function(enrichResultsTable, gmt, groupColumn=NULL, 
                                                      cutHeight = 0.99, broadest=FALSE, max_pAdjust = 0.01,
                                                      hclustMethod = "complete"){
  if (length(unique(enrichResultsTable$ID)) < 2){
    message ("Nothing to simplify")
    return (list (enrichResultsTable, data.frame()))
  }
  setDT(enrichResultsTable); setDT(gmt)
  
  ##Select Significant Terms
  #target_overrep_sig <- enrichResultsTable[p.adjust < max(max_pAdjust,0.05),]#qvalue < 0.01]
  target_overrep_sig <- enrichResultsTable[p.adjust < max_pAdjust,]#qvalue < 0.01]
  
  ##Prepare Significant GO Term Jaccard Similarity Matrix
  sig_go_terms <- unique(target_overrep_sig$ID)
  
  if (!is.null(groupColumn)){
    message ("Computing universal gene overlap between ", length(sig_go_terms), " significant GO terms from ", length(unique(enrichResultsTable[[groupColumn]])), " ", groupColumn, "(s)")
  }else{
    message ("Computing universal gene overlap between ", length(sig_go_terms), " significant GO terms")
  }
  
  if(groupColumn == "cluster"){
    message ("When group column = cluster, it is renamed to cluster.x to avoid clashing with column names used for GO clusters")
    groupColumn <- "cluster.x"
    setnames(enrichResultsTable, old= "cluster", new = "cluster.x")
  }
  
  # expect gmt in term2gene first 2 column format
  gmt.subset <- gmt[gmt[[1]] %in% sig_go_terms, 1:2]
  setnames(gmt.subset, new = c("ont", "gene"))
  gmt.subset <- gmt[ont %in% sig_go_terms, .(ont=factor(ont), gene=factor(gene))]
  termByGeneMat <-  Matrix::sparseMatrix(as.integer(gmt.subset$ont), as.integer(gmt.subset$gene), 
                                         dimnames=list(levels(gmt.subset$ont), 
                                                       levels(gmt.subset$gene)))
  
  #go_dist_mat <- dist(termByGeneMat, method="binary")
  # maybe this is faster?  I haven't clocked it
  go_dist_mat <- parallelDist::parDist(as.matrix(termByGeneMat), method="binary")
  
  if (hclustMethod != "complete" & cutHeight == 0.99){
      message ("You requested to cluster using method ", hclustMethod, " but didn't change from default cutHeight from 0.99. It is recommended you adjust cutHeight")
  }
  hc <- hclust(go_dist_mat, method = hclustMethod)
  
  clusters <- cutree(hc, h=cutHeight)
  clusters <- data.table (cluster = as.numeric(clusters), ID = attributes(clusters)$names )
  
  message ("GO terms clustered into ", max(clusters$cluster), " clusters")
  
  ## go gene set lengths
  gmt.setLengths <- gmt.subset[,.(setSize = length(unique(gene))), by = ont]
  clusters <- merge (clusters, gmt.setLengths, by.x="ID", by.y="ont")
  
  ## local gene set lengths  -- this will only count those genes that make it into an enriched list
  id.gene.long <- enrichResultsTable[,.(ID, gene = unlist(strsplit(geneID, "/"))), by = seq_len(nrow(enrichResultsTable))]
  genesPerTerm <- id.gene.long[,.(localSetLength=length(unique(gene))), by = ID]
  clusters <- merge (clusters, genesPerTerm, by.x="ID", by.y="ID")
  #setorder(clusters, -localSetLength)  # the tie breaker
  
  clusterInfo <- merge (enrichResultsTable, clusters, by = "ID")
  
  clusterInfo[,maxSet := max (setSize), by = c("cluster", groupColumn)]
  #winners <- clusterInfo[,.SD[which(count == maxSet)], by = .(cluster, Bait)]  #keeps all ties...needs updating
  
  if (broadest){
    winners <- clusterInfo[p.adjust < max_pAdjust,.SD[which.max(setSize),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners$ID)), " representative GO terms choosing the BROADEST (in GO) significant term per GO-cluster per ", groupColumn)
  }else{
    winners <- clusterInfo[p.adjust < max_pAdjust,.SD[which.min(p.adjust),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
    message (length(unique(winners$ID)), " representative GO terms choosing the MOST significant term per GO-cluster per ", groupColumn)
  }
  result <- enrichResultsTable[ID %in% winners$ID,]
  result[clusterInfo, cluster.id := cluster, on = "ID"]
  list(simplified = result, clusterInfo = clusterInfo, tree = hc)
}



fixMsigdbGONames <- function(names){
  names <- gsub("^GO_","",names)
  names <- gsub("_"," ",names)
  names <- tolower(names)
  patterns <- c("\\bdna|dna\\b", 
                "\\brna|rna\\b", 
                "\\bgtp|gtp\\b",
                "\\batp|atp\\b"
                )
  replacements <- c("DNA", 
                    "RNA", 
                    "GTP",
                    "ATP"
                    ) 
  for (i in seq_along(patterns)){
    names <- gsub (patterns[i], replacements[i], names)
  }
  #capitalize first letter
  substr(names,1,1) <- toupper(substr(names,1,1))
  return(names)
}
test <- function(){
  testStrings <- c("external", "rna processing", "gtpase", "atpase")
  correct <- c("External", "RNA processing", "GTPase", "ATPase")
  if  (!all(fixMsigdbGONames(testStrings) == correct)){
    message("fixMsigdbGONames FAIL ", paste(fixMsigdbGONames(testStrings), " "))
  }
}
test()


heatmapNumbered <- function (main.mat, counts.mat, negCols = NULL, title="",
                             borderMatrix = NULL, borderColFun = NULL,
                             borderMM = 2,
                             brewerPalette = "Blues",
                             show_column_dend = FALSE,
                             show_row_dend = FALSE,
                             border = TRUE,
                             max_pAdjust = 0.01,
                             row_names_gp = gpar(fontsize = 10),
                             column_names_gp = gpar(fontsize = 10),
                             upperThreshold = NULL,
                             colors = NULL,
                             ...){

  heatmap_legend_param = list(legend_direction="horizontal", title = "-log10(adj.p)") #global enrich\n
  
  if (is.null(colors)){  # by default we pass a long vector of colors and let ComplexHeatmap define the ranges
    Blues = colorRampPalette(RColorBrewer::brewer.pal(9, brewerPalette))
    colors <- Blues(100)
    
    # if we get an upperTreshold, we define the limits in the colors object
    if(!is.null(upperThreshold))
      colors <- circlize::colorRamp2(breaks =seq(0, upperThreshold, length.out = 100), colors = colors)
    
    
    if (!is.null(negCols)){
      limit <- ifelse(is.null (upperThreshold), 4, upperThreshold)
      
      #colors <- circlize::colorRamp2 (breaks=seq(from=-max(main.mat), to = max(main.mat), length.out=101), colors =colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))
      colors <- circlize::colorRamp2 (breaks=seq(from=-limit, to = limit, length.out=101), colors =colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101))
      main.mat[,negCols] = -main.mat[, negCols]
      heatmap_legend_param = c (heatmap_legend_param, list(at=c(-limit,-limit/2,0,limit/2,limit), labels = c(limit,limit/2,0,limit/2,limit)) )
    }
  }  
  
  ##Plot main figure heatmap
  hm <- ComplexHeatmap::Heatmap(main.mat, col = colors, border = border, rect_gp = gpar(col = "grey", lwd = 1),
                                #cluster_rows = ddr,
                                column_title = title,
                                column_names_rot = 90, row_names_gp = row_names_gp, column_names_gp = column_names_gp,
                                show_row_dend = show_row_dend, show_column_dend = show_column_dend, heatmap_legend_param = heatmap_legend_param,
                                row_names_max_width = 2*max_text_width(rownames(main.mat), gp = gpar(fontsize = 12)),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if (!is.null(borderMatrix) & !is.null(borderColFun)){
                                    lwd <- unit(borderMM,"mm")
                                    grid.rect(x, y, width, height, gp = gpar(fill = borderColFun(borderMatrix[i,j]), col = NA))
                                    grid.rect(x, y, width-lwd, height-lwd, gp = gpar(fill = fill, col = NA))
                                  }
                                  if (!is.na(counts.mat[i,j])){
                                    color <- ifelse (abs(main.mat[i,j]) < -log10(max_pAdjust), "grey", "white") # "white" #
                                    grid.text(sprintf("%.0f", counts.mat[i, j]), x, y, gp = gpar(fontsize=10, col=color))
                                  }
                                }, ...) #+1  # this makes it a list!
  
  
  if (!is.null(borderMatrix) & !is.null(borderColFun)){
    legendList <-  list (Legend(col_fun = colorFun, title= "viral enrich\n-log10(p)"))
    hm <- hm + 1 # this makes it a list so I can add annotations to it
  } else{
    legendList <- list()
  }
  
  hm <- draw(hm,heatmap_legend_side="top", annotation_legend_list = legendList,
             annotation_legend_side = "top")
  
  invisible (hm)
  
}


library (ComplexHeatmap)
enrichHeatmapBestPerGroup <- function(simplifiedEnrichTable, fullEnrichTable, groupColumn="bait", topN = 1, title="", cols = NULL, 
                                      negCols = NULL, reduceRedundantsAcrossGroups=TRUE, max_pAdjust = 0.01, minCount = 1,
                                      annotatePossibleMatches = TRUE,  row_names_gp = gpar(fontsize = 10),
                                      upperThreshold  = NULL,
                                      pvalColumn = "p.adjust", ...){
  setorderv(simplifiedEnrichTable, cols = pvalColumn)
  bestTermPerBait <- simplifiedEnrichTable[simplifiedEnrichTable[[pvalColumn]]<max_pAdjust & Count >= minCount,.(ID=ID[1:topN]),by=groupColumn]

    if (is.null(fullEnrichTable)){
    fullEnrichTable <- simplifiedEnrichTable
    reduceRedundantsAcrossGroups <- FALSE
  }

  if(reduceRedundantsAcrossGroups){  
    #reduce redundant based on clusters in fullEnrichTable
    countsByID <- fullEnrichTable[ID %in% bestTermPerBait$ID, .(geneCount  = length(unique(unlist(strsplit(geneID, "/"))))), by = .(ID, cluster)]
    # get the term with most genes across whole dataset per term-cluster
    setorder(countsByID, -geneCount)
    bestTerms <- countsByID[,.SD[1],by=cluster]$ID
  } else bestTerms <- unique(bestTermPerBait$ID)
  
  if (!is.null(negCols) & !is.null(cols)){
    negOnly <- setdiff(negCols, cols)
    if (length(negOnly) > 0)
      message ("Columns specified in negCols not found in cols, these will be removed ", paste (negOnly, collapse = ", "))
    negCols <- intersect (cols, negCols)
  }
  if (!is.null(negCols) & length(negCols) == 0){
    warning ("negCols is set to an empty vector. Did you get the right names? Set to NULL(default) for no negative columns")
  }
    
  main.wide <- dcast (fullEnrichTable[ID %in% bestTerms], as.formula(paste("Description", groupColumn, sep="~")), value.var="p.adjust")
  for(col in unique(c(cols,negCols))){
    if (is.null(main.wide[[col]])) main.wide[[col]] <- NA
  }
  
  main.mat <- -log10(as.matrix(main.wide, rownames = "Description"))
  main.mat[is.na(main.mat)] <- 0
  if (!is.null(upperThreshold)){
    main.mat[main.mat > upperThreshold] <- upperThreshold
  }
  if (all(grepl("^GO_", rownames(main.mat)))){
    rownames(main.mat) <- fixMsigdbGONames(rownames(main.mat))
  }
  
  counts.wide <- dcast (fullEnrichTable[ID %in% bestTerms], as.formula(paste("Description", groupColumn, sep="~")), value.var="Count")
  for(col in unique(c(cols,negCols))){
    if (is.null(counts.wide[[col]])) counts.wide[[col]] <- NA
  }
  counts.mat <- as.matrix(counts.wide, rownames="Description")
  
  
  geneTable <- fullEnrichTable[ID %in% bestTerms, .(gene = unlist(strsplit(geneID, split="/"))),by = ID]
  geneTable[,cleanName := fixMsigdbGONames(ID)]
  
  if (!is.null(cols)){
    if (! all(cols %in% colnames(counts.mat) & cols %in% colnames(main.mat))){
      message ("Not all requested columns for heatmap found in data")
      message ("main.mat ", paste(colnames(main.mat), collapse=" "))
      message ("counts.mat ", paste(colnames(counts.mat), collapse=" "))
    }else{
      counts.mat<- counts.mat[,cols]
      main.mat<- main.mat[,cols]
    }
  }

  if (annotatePossibleMatches==TRUE){
    genesInUniverseCounts <- unique(fullEnrichTable[, .( geneCount = as.integer(gsub("[0-9]+/", "", GeneRatio))), by = c(groupColumn)])
    if (nrow(genesInUniverseCounts) != length(unique(genesInUniverseCounts[[groupColumn]]))){
      stop("non-unique gene counts per group. If you didn't combine multiple differently grouped enrichments, this is unexpected. If it is, set annotatePossibleMatches = FALSE")
    }
    cols <- colnames(main.mat)
    setkeyv(genesInUniverseCounts, groupColumn)
    topBars <- HeatmapAnnotation(`Group Sizes` = anno_barplot ( genesInUniverseCounts[cols, geneCount] ))
  #   if (!is.null(top_annotation)){
  #     warning("over writing non-null top annotation with possible matches")
  #     #top_annotation <- top_annotation %v% topBars
  #   }#else{
  #     top_annotation <- topBars
  #   #}
  }else{
    topBars <- NULL
  }
  
  hm <- heatmapNumbered (main.mat, counts.mat, negCols, title, max_pAdjust = max_pAdjust, bottom_annotation = topBars, row_names_gp = row_names_gp,
                         upperThreshold = upperThreshold,...)
  
  invisible(list(geneTable = geneTable, main.mat = main.mat, counts.mat = counts.mat, hmList = hm))
}

#' Runs the three main GO (etc) enrichment functions together
#' 

enrichmentOnGroupsPL <- function (groupTable, geneColumn, groupColumns, gmt,universe = NULL, numProcessors = 1,
                                  max_pAdjust = 0.1, broadest = FALSE,
                                  topN = 4, reduceRedundantsAcrossGroups = TRUE,
                                  ...){
  
  
  enrich.dt  <- enricherOnGroups(groupTable, geneColumn, groupColumns, gmt, universe, numProcessors)
  groupColumn <- paste(groupColumns, collapse = ".")
  simp <- simplifyEnrichBySimilarUniverseMembership(enrich.dt, gmt, groupColumn,
                                                    cutHeight = 0.99, broadest = FALSE,
                                                    max_pAdjust = max_pAdjust)
  enrichHM <- enrichHeatmapBestPerGroup(simp[[1]], simp[[2]], 
                                        groupColumn, topN,
                                        reduceRedundantsAcrossGroups = reduceRedundantsAcrossGroups,
                                        max_pAdjust = max_pAdjust,
                                        ...)
  
  return (list(enrich.dt = enrich.dt, simp = simp, hm =  enrichHM))
  
}









loadGmtFromBioconductor <- function (dbName = "org.Hs.eg.db", ontology = "BP", keyType = "UNIPROT"){ #c("UNIPROT", "SYMBOL"
  message ("Using package ", dbName, " version ", packageVersion(dbName))
  GO <- clusterProfiler:::get_GO_data(dbName, ontology, keyType)
  gmt <- rbindlist(lapply (GO$EXTID2PATHID, function(x) data.table(ont.id = x)), idcol="gene")
  gmt$ont <- GO$PATHID2NAME[gmt$ont.id]
  setcolorder(gmt, c("ont", "gene", "ont.id"))
  return(gmt)
}



loadKegg <- function (organism=c("hsa", "mmu")[1], keyType = c("uniprot", "kegg", "ncbi-geneid", "ncbi-proteinid")[1]){
  message ("Current version of KEGG (clusterProfiler might use its cache, look for download messages below the KEGG Info)\n", format(Sys.time(), "%Y_%m_%d"))
  # download and display current kegg info:
  f <- tempfile()
  utils::download.file("https://rest.kegg.jp/info/kegg", f)
  message(paste0(readLines(f), collapse = "\n"))
  
  
  KEGG <- clusterProfiler:::prepare_KEGG(organism, "KEGG", keyType)
  gmt <- rbindlist(lapply (KEGG$EXTID2PATHID, function(x) data.table(ont.id = x)), idcol="gene")
  gmt$ont <- KEGG$PATHID2NAME[gmt$ont.id]
  setcolorder(gmt, c("ont", "gene", "ont.id"))
  return(gmt)
  
}


#  splitCircleHeatMap
# Instead of shaded rectangles, this plots split circles of varying size and color in each cell.
# For plotting up to four values per cell, for example up/down size and statistical strength.
# Size is denoted by area of the half circle and statistical strength by color. 
# This function takes up to four numeric matrices as arguments.
# If only the *Color matrices are passed in, they will be used for size as well.
# Size matrices must be positive.  Color matrices are more flexible and the mapping
# can be controlled through colFunLeft and colFunRight
# 
# It currently does not print a legend showing size of semicircles.
# 
# This relies on ComplexHeatmap and many parameters can be successfully passed through to Heatmap
# 
# For examples, see test_splitCircleHeatMap


splitCircleHeatMap <- function (matLeftColor, matRightColor,
                                matLeftSize=NULL, matRightSize=NULL,
                                colFunLeft = NULL, colFunRight = NULL,
                                leftLegendTitle = "left", rightLegendTitle = "right",
                                sizeLegendTitle = "size",
                                legends = NULL, ...){
  
  #sqrt <- identity  # a hack to test if ~area or ~radius looks better
  if (is.null(colFunLeft))colFunLeft <- circlize::colorRamp2(c(0, max(matLeftColor)), c("#EEEEEE", "blue"))
  if (is.null(colFunRight))colFunRight <- circlize::colorRamp2(c(0, max(matRightColor)), c("#EEEEEE", "red"))
  
  if (is.null(matLeftSize))matLeftSize <- matLeftColor
  if (is.null(matRightSize))matRightSize <- matRightColor
  
  stopifnot(all(matLeftSize > 0, na.rm=TRUE), all(matRightSize > 0, na.rm=TRUE))
  maxSize <- max(c(matLeftSize, matRightSize), na.rm=TRUE)
  minSize <- min(c(matLeftSize, matRightSize), na.rm=TRUE)

  # size legend breaks
  desiredBreaks <- 4
  breaks = labeling::extended(minSize, maxSize, m=desiredBreaks)
  # 
  if (any(breaks == 0)){
    if(length(breaks)<=desiredBreaks)
      breaks[breaks==0] <- min(breaks[breaks!=0])/2
    else breaks <- breaks[breaks!=0]
  }
  maxSize <- max(breaks)
  
  
  
  matLeftSize <- sqrt(matLeftSize/maxSize)
  matRightSize <- sqrt(matRightSize/maxSize)

  
  
  
  
  # define the function that will do the work of drawing inside a cell
  cell_fun <- function(j, i, x, y, width, height, fill){
    
    # get actual cell sizes in mm so we can draw accurate legends later
    if (is.null(cellWidthMM)){
      cellWidthMM <<- convertWidth(width, "mm")
      cellHeightMM <<- convertHeight(height, "mm")
    }
    
    # trick to draw half circles is to define a viewport that will clip, and then draw
    # the circle at the edge of the viewport.
    # https://stackoverflow.com/questions/31538534/plotting-half-circles-in-r
    vp <- viewport (x-0.5*width,y, width = width, height= height, clip="on")
    grid.circle(1.0,0.5, r = abs(matLeftSize[i, j])/2, 
                gp = gpar(fill = colFunLeft(matLeftColor[i, j]), col = NA), vp = vp)
    
    vp <- viewport (x+0.5*width,y, width = width, height= height, clip="on")
    grid.circle(0,0.5, r = abs(matRightSize[i, j])/2, 
                gp = gpar(fill = colFunRight(matRightColor[i, j]), col = NA), vp = vp)
    
  }
  
  if (is.null(legends)){
    legends <- ComplexHeatmap::packLegend(Legend(col_fun = colFunLeft, title= leftLegendTitle),
                          Legend(col_fun = colFunRight, title= rightLegendTitle))
  }  
  hm<-ComplexHeatmap::Heatmap (matLeftColor,  #this matrix is what it will cluster and label with, otherwise this is ignored
               rect_gp = gpar(type = "none"), 
               cell_fun = cell_fun,
               show_heatmap_legend=FALSE,
               row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
               ...)
  
  # set up variables to catch actual cell size in MM as drawing happens
  cellWidthMM <- NULL
  cellHeightMM <- NULL
  draw(hm, heatmap_legend_list = legends)
  
  # size legend
  desiredBreaks <- 4
  breaks = labeling::extended(minSize, maxSize, m=desiredBreaks)
  # 
  if (any(breaks == 0)){
    if(length(breaks)<=desiredBreaks)
      breaks[breaks==0] <- min(breaks[breaks!=0])/2
    else breaks <- breaks[breaks!=0]
  }
  sizeLegend <- splitCircleLegend(sizeLegendTitle, cellHeightMM, cellWidthMM,  maxVal = maxSize, breaks=breaks)
  
  return (list(legend = sizeLegend, cellWidth = cellWidthMM, cellHeight = cellHeightMM))
}


# Constructs a legend modeled on ComplexHeatmap legends. This can not be accurately constructed
# until the splitCircle heatmap actually draws to a device, so it can't be printed (afaik) using
# standard ComplexHetmap legend functions.
#
# See test_splitCircleHeatMap for one way to use.  Editing location of legend in a drawing program
# may be easiest, but you can certainly play with viewport coordinates
#
# this does not resize in interactive graphics devices that resize the heatmap, and so will not be
# accurate in RStudio's Plots window which seems to draw off-screen before resizing to current window size

splitCircleLegend <- function(title, cellHeight, cellWidth, breaks=c(5,10,20,40), maxVal=NULL, col = function(x){"gray"}){
  # create a legend body
  radii <- sqrt(breaks)
  n <- length(breaks)
  if (is.null(maxVal)){
    maxVal <- max(breaks)
  }
  cellHeight <- min(cellHeight, cellWidth)
  viewPorts <- lapply (0:(n-1), FUN=function(i) {viewport(y= i * cellHeight, height=cellHeight, width=cellWidth, clip=TRUE)})
  stackedCircles <- lapply (1:n, FUN=function(i) circleGrob( x = 0, r = 0.5/max(sqrt(maxVal)) * radii[i], gp = gpar(fill=col(breaks[i]), col=NA), vp=viewPorts[[i]]))
  
  labels <- textGrob( label = breaks, x=0.5, y =  (0:(n-1)) * cellHeight, just = "left")
  withLabels <- c (stackedCircles, list(labels=labels))
  
  class(withLabels) = "gList"
  gt = gTree(children = withLabels, cl = "legend_body", vp = viewport(width = cellWidth, 
                                                                      height = cellHeight * n))
  attr(gt, "height") = cellHeight * n
  attr(gt, "width") = cellWidth
  
  legendBody <- gt
  
  
  #, gp = title_gp
  title_grob <- textGrob(title)
  title_height <- convertHeight(grobHeight(title_grob), "mm")
  title_width <- convertWidth(grobWidth(title_grob), "mm")
  title_y = unit(1, "npc")
  title_just = c("left", "top")
  
  total_height <- title_height + convertHeight(grobHeight(legendBody), "mm")
  total_width <- max(title_width, cellWidth*2)
  legendVP <- viewport (width = total_width, height = total_height)
  
  gf = grobTree(textGrob(title, x = unit(0, "npc"), y = title_y, 
                         just = title_just), (legendBody), vp = viewport(width = total_width, height = total_height), 
                cl = "legend")
  attr(gf, "width") = total_width
  attr(gf, "height") = total_height
  
  object = new("Legends")
  object@grob = gf
  object@type = "single_legend"
  object@n = 1
  
  
  return (object)
}



test_splitCircleHeatMap <- function(){
  library (ComplexHeatmap)
  library (circlize)
  
  rn <- c("Rahman, Taylor", "Hageman, Allie", "el-Saeed, Misfar", "Brown, Rohan", 
          "el-Shah, Habeeba", "Dougal, William", "Schubert, William", "Laut, Jessica", 
          "Wimberly, Elias", "Garcia-Espinoza, Christian", "Buccieri, Thaylor", 
          "Magno, Paige", "Garcia, Jessica", "al-Ayoub, Rasheeq", "Woodall, Taylor", 
          "Thatcher, Alex", "Lin, Kyle", "Walck, William", "el-Qazi, Mutlaq", 
          "el-Azzi, Zain")
  
  cn <- c("Nguyen Do, Michael", "Begay, David", "Martinez, Marta", "Hotchkiss, Sheela", 
          "Beightel, Ruben", "Trujillo, Cameron", "King, Mary", "Aquiningoc, Yu Sung", 
          "Hall, Michaela", "Maestas, Thalia", "Ross, Lily", "Smith, Faith", 
          "al-Tabet, Muhyddeen", "Gurule, Alicia", "Tardif, Jordan", "Erdenebat, Taylor", 
          "el-Radwan, Saalima", "el-Aslam, Shaheed", "el-Omer, Azeema", 
          "Verde, Brandon")
  
  r <- splitCircleHeatMap(matrix (runif(400), nrow=20, ncol=20, dimnames = list (rn,cn)),
                          matrix (runif(400), nrow=20, ncol=20),
                          matrix (runif(400) * 50, nrow=20, ncol=20),
                          matrix (runif(400) * 50, nrow=20, ncol=20),
                          leftLegendTitle = "decreasers",
                          rightLegendTitle = "increasers",
                          #these next will get passed right to ComplexHeatmap
                          cluster_columns=FALSE, 
                          column_split = sample(c(1,2), 20, replace=TRUE, prob=c(0.2,0.8)),
                          row_km = 4,
                          column_gap = unit(1, "cm"),
                          row_gap = unit(0.5,"cm"))

  #draw legend in upper right hand corner  
  pushViewport(viewport(1,1,just = c("right", "top"), height = grobHeight(r$legend@grob), width = grobWidth(r$legend@grob)))
  grid.draw(r$legend)
  
}



loadCORUMasGMT <- function (path = NULL, species = c("HUMAN", "MOUSE"), idType = c("symbol", "uniprot")){
  if (is.null(path))
    stop ("Download and unzip file from ", "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip")
  corumDT <- fread (path)[toupper(Organism) == toupper (species)]
  if(nrow(corumDT) == 0)
    stop("no rows read from file that match the species requested", path, species)
  # lookup the column name
  idColumn <- c(uniprot = "subunits(UniProt IDs)",
                symbol =  "subunits(Gene name)")[idType]
  # just hte columns desired
  subTable <- corumDT[, .SD, by = .(ComplexID, ComplexName), .SDcols = idColumn]
  # expand the genes/uniprots column
  expanded <- subTable[, .(gene = unlist(strsplit(.SD[[idType[1]]], ";"))), by = . (ComplexID, ComplexName) ]
  # reorder to expected term2gene format
  return(expanded[, .(ont = ComplexName, gene, ComplexID)])
}


limitGMT2PantherGOslim <- function (gmt.dt, fileOrURL = "http://data.pantherdb.org/PANTHER16.0/ontology/PANTHERGOslim.obo"){
  slimGO <- ontologyIndex::get_ontology(fileOrURL)
  gmt.dt[ont.id %in% names(slimGO$name)]
}


# fgsea #################

#' Given a matrix and a column name perform fgsea on a single column
#' @param sets a named list of genes/proteins vectors. Each vector is a single set. vector items must match rownames of mat
#' @param scoreType usually "std" for positive and negative values together, or "pos" for positive only
oneColumnFGSEA <- function (columnName, mat, sets, scoreType = c("std", "pos", "neg"), fgseaFunction = fgsea::fgsea, ...){
  #print (columnName)
  log2FC <- mat[,columnName]
  
  # fgsea fails with infinite values
  if (any (is.infinite(log2FC))){
    message (sprintf ("Removing %d infinite values from scores for %s", sum(is.infinite(log2FC)), columnName))
    log2FC <- log2FC[is.finite(log2FC)]
  }
  
  if (var(log2FC, na.rm = TRUE) == 0){
    message ("No variance detected in ", columnName, " returning NULL to avoid a fgsea crash")
    return (NULL)
  }
  
  if (anyDuplicated(names(log2FC))){
    message ("Duplicate row names found, keeping only the greatest (this favors positive over negative and favors those more likely to be duplicated (maybe larger proteins))")
    log2FC <- sort (log2FC, decreasing = TRUE)
    log2FC <- log2FC[!duplicated(names(log2FC))]
    
  }
  
  # reorder randomly so ties are less likely to influence results
  log2FC <- sample(log2FC, length(log2FC))
  
  # fgsea fails with missing values
  log2FC <- log2FC[!is.na(log2FC)]
  

  
  #nperm=1000, gseaWeightParam = 1, nproc=1
  seaRes <- fgseaFunction(pathways = sets, stats = log2FC,  scoreType = scoreType, ...)
  seaRes[, sigScore := -log10(pval) * ifelse(ES < 0, -1, 1) ]
  return (seaRes)
}

#' Given a matrix run fgsea on all columns, return table of results with column "group" identifying the column of the matrix
#' @param mat
#' @param sets a named list of genes/proteins vectors. Each vector is a single set. vector items must match rownames of mat
#'             can optionally pass a term2gene data.frame with terms in first column and genes/uniprots in second
#' @param scoreType is chosen based on presence of negative or positive values in the matrix
#' @param ... arguments to pass to oneColumnFGSEA and ultimately fgsea::fgsea
matrixFGSEA <- function (mat, sets, ...){
  anyPositive <- any(mat > 0, na.rm = TRUE)
  anyNegative <- any(mat < 0, na.rm = TRUE)
  if (anyPositive & anyNegative)
    scoreType = "std"
  if (anyPositive & !anyNegative)
    scoreType = "pos"
  if (!anyPositive & anyNegative)
    scoreType= "neg"
  if(!anyPositive & !anyNegative)
    stop("No variance or no values detected in matrix")
  
  
  if ("data.frame" %in% class(sets)){
    # assume its in term2gene format and split the genes (2nd col) by term(1st)
    sets <- split(sets[[2]], sets[[1]])
  }
  
  all.sea <- pbapply::pblapply (colnames(mat), oneColumnFGSEA, mat, sets, scoreType, ...)
  names(all.sea) <- colnames(mat)
  sea.dt <- rbindlist(all.sea, idcol = "group")
  
  return (sea.dt)
}


# other stuff ###############


enrichmentAnnotationHeatmap <- function(enrich.dt = data.table (group = c(), pvalue = c(), term = c(), geneID = c()),
                                        terms,
                                        geneGroups = data.table(group = c(), gene = c()),
                                        matrixRowOrder,
                                        enrichPColors = c("white", "#807DBA"),      # c("white", RColorBrewer::brewer.pal(9, "Purples")[c(6)])
                                        goMatchColor = "#3F007D",                   # RColorBrewer::brewer.pal(9, "Purples")[c(9)]
                                        ...){
  # this only works if the genes are divided into mutually exclusive groups.  So make sure the number of groups per gene is exactly 1
  stopifnot (all(geneGroups[, .N, by = gene]$N==1))
  
  # define the go binary matrix
  miniGMT <- unique(enrich.dt[term %in% terms & group %in% unique(geneGroups$group), 
                              .(gene = unlist(strsplit(geneID, "/"))),
                              by = .(term) ])
  go.mat <- dcast (miniGMT[matrixRowOrder, , on = "gene"], # get all genes in table, this will generate an NA column which I clean up below
                   gene~term, 
                   fun.aggregate = length)[, `NA` := NULL] |>   
    as.matrix(rownames = "gene")
  
  # go pvalue matrix (by gene)
  p.mat <- dcast (enrich.dt[term %in% terms
  ][geneGroups, , on = "group", allow.cartesian = TRUE], #
  gene~term, value.var = "pvalue") |> 
    as.matrix(rownames = "gene")
  
  # missing values
  if (any(is.na(p.mat))){
    message ("not all terms have a p.value in all groups, setting those to p = 1.0")
    p.mat[is.na(p.mat)] <- 1.0
    
  }

  
  # enforce good ordering
  go.mat <- go.mat[matrixRowOrder, terms]
  p.mat <- p.mat[matrixRowOrder, terms]
  
    
  Heatmap(-log10(p.mat),
          name = "enrichment\n-log10(pvalue)",
          col = enrichPColors,
          layer_fun = function (j, i, x, y, width, height, fill){
            color <- ifelse(pindex (go.mat,i,j) > 0, goMatchColor, "transparent")
            grid.rect(x,y,width, height, gp = gpar (fill = color, col = NA))
          },
          show_row_dend = FALSE,
          ...)
}

