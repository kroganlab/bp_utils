library (data.table)


simplifyEnrichBySimilarUniverseMembership <- function(enrichResultsTable, gmt, groupColumn="bait"){
  if (length(unique(enrichResultsTable$ID)) < 2){
    message ("Nothing to simplify")
    return (list (enrichResultsTable, data.frame()))
  }
  setDT(enrichResultsTable); setDT(gmt)
  
  ##Select Significant Terms
  target_overrep_sig <- enrichResultsTable[p.adjust < 0.05,]#qvalue < 0.01]
  
  ##Prepare Significant GO Term Jaccard Similarity Matrix
  sig_go_terms <- unique(target_overrep_sig$ID)
  
  message ("Computing universal gene overlap between ", length(sig_go_terms), " significant GO terms from ", length(unique(enrichResultsTable[[groupColumn]])), " ", groupColumn, "(s)")
  
  gmt.subset <- gmt[ont %in% sig_go_terms, .(ont=factor(ont), gene=factor(gene))]
  termByGeneMat <-  Matrix::sparseMatrix(as.integer(gmt.subset$ont), as.integer(gmt.subset$gene), 
                                         dimnames=list(levels(gmt.subset$ont), 
                                                       levels(gmt.subset$gene)))
  
  go_dist_mat <- dist(termByGeneMat, method="binary")
  hc <- hclust(go_dist_mat)
  
  clusters <- cutree(hc, h=0.99)
  clusters <- data.table (cluster = as.numeric(clusters), ID = attributes(clusters)$names )
  
  message ("GO terms clustered into ", max(clusters$cluster), " clusters")
  
  ## go gene set lengths
  gmt.setLengths <- gmt[,.(setSize = length(unique(gene))), by = ont]
  clusters <- merge (clusters, gmt.setLengths, by.x="ID", by.y="ont")
  
  ## local gene set lengths  -- this will only count those genes that make it into an enriched list
  id.gene.long <- enrichResultsTable[,.(ID, gene = unlist(strsplit(geneID, "/"))), by = seq_len(nrow(enrichResultsTable))]
  genesPerTerm <- id.gene.long[,.(localSetLength=length(unique(gene))), by = ID]
  clusters <- merge (clusters, genesPerTerm, by.x="ID", by.y="ID")
  #setorder(clusters, -localSetLength)  # the tie breaker
  
  clusterInfo <- merge (target_overrep_sig, clusters, by = "ID")
  
  clusterInfo[,maxSet := max (setSize), by = c("cluster", groupColumn)]
  #winners <- clusterInfo[,.SD[which(count == maxSet)], by = .(cluster, Bait)]  #keeps all ties...needs updating
  winners <- clusterInfo[,.SD[which.max(setSize),],by=c("cluster", groupColumn)]  #chooses the first in case of tie breakers
  
  message (length(unique(winners$ID)), " representative GO terms choosing the broadest significant term per GO-cluster per ", groupColumn)
  
  result <- enrichResultsTable[ID %in% winners$ID,]
  list(simplified = result, clusterInfo = clusterInfo)
}




fixMsigdbGONames <- function(names){
  names <- gsub("GO_","",names)
  rnames <- gsub("_"," ",names)
  names <- tolower(names)
}

enrichHeatmapBestPerGroup <- function(simplifiedEnrichTable, fullEnrichTable, groupColumn="bait", topN = 1, title=""){
  setorder(simplifiedEnrichTable, p.adjust)
  bestTermPerBait <- simplifiedEnrichTable[p.adjust<0.05,.(ID=ID[1:topN]),by=groupColumn]
  
  main.wide <- dcast (fullEnrichTable[ID %in% bestTermPerBait$ID], as.formula(paste("Description", groupColumn, sep="~")), value.var="p.adjust")
  main.mat <- -log10(as.matrix(main.wide, rownames = "Description"))
  main.mat[is.na(main.mat)] <- 0
  main.mat[main.mat > 5] <- 5
  rownames(main.mat) <- fixMsigdbGONames(rownames(main.mat))
  
  counts.wide <- dcast (fullEnrichTable[ID %in% bestTermPerBait$ID], as.formula(paste("Description", groupColumn, sep="~")), value.var="Count")
  counts.mat <- as.matrix(counts.wide, rownames="Description")
  #print(str(counts.mat))
  
  
  Blues = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  
  
  ##Plot main figure heatmap
  hm <- ComplexHeatmap::Heatmap(main.mat, col = Blues(100), border = TRUE, rect_gp = gpar(col = "grey", lwd = 1),
                                column_title = title,
                                column_names_rot = 90, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                                show_row_dend = FALSE, show_column_dend = FALSE, heatmap_legend_param = list(legend_direction="horizontal", title = "-log10(q-value)"),
                                row_names_max_width = max_text_width(rownames(main.mat), gp = gpar(fontsize = 12)),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if (!is.na(counts.mat[i,j])){
                                    grid.text(sprintf("%.0f", counts.mat[i, j]), x, y, gp = gpar(fontsize=10, col="red"))
                                  }
                                })
  
  draw(hm,heatmap_legend_side="top")
}
