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
  names <- gsub("^GO_","",names)
  names <- gsub("_"," ",names)
  names <- tolower(names)
}

enrichHeatmapBestPerGroup <- function(simplifiedEnrichTable, fullEnrichTable, groupColumn="bait", topN = 1, title="", cols = NULL, negCols = NULL,...){
  setorder(simplifiedEnrichTable, p.adjust)
  bestTermPerBait <- simplifiedEnrichTable[p.adjust<0.05,.(ID=ID[1:topN]),by=groupColumn]
  
  main.wide <- dcast (fullEnrichTable[ID %in% bestTermPerBait$ID], as.formula(paste("Description", groupColumn, sep="~")), value.var="p.adjust")
  main.mat <- -log10(as.matrix(main.wide, rownames = "Description"))
  main.mat[is.na(main.mat)] <- 0
  main.mat[main.mat > 5] <- 5
  rownames(main.mat) <- fixMsigdbGONames(rownames(main.mat))
  
  counts.wide <- dcast (fullEnrichTable[ID %in% bestTermPerBait$ID], as.formula(paste("Description", groupColumn, sep="~")), value.var="Count")
  counts.mat <- as.matrix(counts.wide, rownames="Description")
  
  
  geneTable <- fullEnrichTable[ID %in% bestTermPerBait$ID, .(gene = unlist(strsplit(geneID, split="/"))),by = ID]
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
  #print(str(counts.mat))
  
  
  Blues = colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  #Reds = colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))
  colors <- Blues(100)
  
  if (!is.null(negCols)){
    colors <- circlize::colorRamp2 (breaks=seq(from=-max(main.mat), to = max(main.mat), length.out=101), colors =rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(101)))
    main.mat[,negCols] = -main.mat[, negCols]
  }
  
  
  ##Plot main figure heatmap
  hm <- ComplexHeatmap::Heatmap(main.mat, col = colors, border = TRUE, rect_gp = gpar(col = "grey", lwd = 1),
                                column_title = title,
                                column_names_rot = 90, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10),
                                show_row_dend = FALSE, show_column_dend = FALSE, heatmap_legend_param = list(legend_direction="horizontal", title = "-log10(q-value)"),
                                row_names_max_width = max_text_width(rownames(main.mat), gp = gpar(fontsize = 12)),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if (!is.na(counts.mat[i,j])){
                                    grid.text(sprintf("%.0f", counts.mat[i, j]), x, y, gp = gpar(fontsize=10, col="white"))
                                  }
                                }, ...)
  
  draw(hm,heatmap_legend_side="top")
  invisible(geneTable)
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
  if (is.null(colFunLeft))colFunLeft <- colorRamp2(c(0, max(matLeftColor)), c("#EEEEEE", "blue"))
  if (is.null(colFunRight))colFunRight <- colorRamp2(c(0, max(matRightColor)), c("#EEEEEE", "red"))
  
  if (is.null(matLeftSize))matLeftSize <- matLeftColor
  if (is.null(matRightSize))matRightSize <- matRightColor
  
  stopifnot(all(matLeftSize > 0), all(matRightSize > 0))
  maxSize <- max(c(matLeftSize, matRightSize), na.rm=TRUE)
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
    legends <- packLegend(Legend(col_fun = colFunLeft, title= leftLegendTitle),
                          Legend(col_fun = colFunRight, title= rightLegendTitle))
  }  
  hm<-Heatmap (matLeftColor,  #this matrix is what it will cluster and label with, otherwise this is ignored
               rect_gp = gpar(type = "none"), 
               cell_fun = cell_fun,
               show_heatmap_legend=FALSE,
               ...)
  
  # set up variables to catch actual cell size in MM as drawing happens
  cellWidthMM <- NULL
  cellHeightMM <- NULL
  draw(hm, heatmap_legend_list = legends)
  
  # size legend
  desiredBreaks <- 4
  breaks = labeling::extended(min(c(matLeftSize, matRightSize)), maxSize, m=desiredBreaks)
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

