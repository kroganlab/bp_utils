require(data.table)
# Modularization of Gene Sets based on network connections


# ######################



# functions to choose best representative term based on enricher results
termScore <- function (GeneRatio, BgRatio){
  gr <- lapply(tstrsplit(GeneRatio, "/"), as.integer)
  br <- lapply(tstrsplit(BgRatio, "/"), as.integer)
  
  return (gr[[1]]**2/(gr[[2]] * br[[1]]))
  
}

enricherJaccard <- function (GeneRatio, BgRatio){
  gr <- lapply(tstrsplit(GeneRatio, "/"), as.integer)
  br <- lapply(tstrsplit(BgRatio, "/"), as.integer)
  
  # jaccard is intersection/union
  
  return (gr[[1]]/ # intersection is genes in group AND term
            (br[[1]]  # union is all genes that match the term...
             + gr[[2]] # plus the genes in group...
             - gr[[1]]) # minus those already counted that match the term 
  )
}






# #######################3

baitPrey2baitModule <- function (tsne.dt, baitPrey.dt){
  moduleMedoids.dt <- tsne.dt[nameMatch == TRUE & type == "prey",
                              .(x = median(x), y = median(y), .N),
                              by = .(clusterID, cluster, name, go, nodeColor)]
  
  # radius based on N
  moduleMedoids.dt[, radius := sqrt(N)]
  
  # count edges bait to module
  baitModule.dt <-baitPrey.dt[tsne.dt[nameMatch == TRUE & type == "prey",], on  = c(Prey = "gene")
                              ][, .N, by =.(Bait, Prey = cluster, cluster, name, go)]
  
  return(baitModule.dt)
}



shiftNodes2AvoidOverlaps <- function (tsne.dt, sumFunction = sum, overlapParam = 0.6,
                                      baitRadius = 2.5,
                                      preyRadius = function(dt)dt[, sqrt(baitPerPrey/4)],
                                      maxIterations = 100){
  
  .oneIteration <- function (subTable, sumFunction = sum, overlapParam = 1.0 ){
    subTable[, dummy := TRUE]
    # collision table
    ct <- subTable[subTable,, on = "dummy", allow.cartesian  = TRUE][gene != i.gene & sqrt((x-i.x)^2 + (y - i.y)^2) < overlapParam * (r + i.r) ]
    cat (nrow(ct), " rows in collision table\n")
    
    # when x,y match exactly, adjust them  based on 1% radius:
    ct[ x == i.x & y == i.y, c("x", "y") := .(ifelse(gene < i.gene, x + r/100, x - r/100),
                                              ifelse(gene < i.gene, y + r/100, y - r/100))]
    
    ct[, c("dx", "dy") := .(x-i.x, y-i.y)]
    ct[, distance := sqrt(dx^2 + dy^2)]
    ct[, overlap  := (i.r + r) - distance ]
    # repel amount = overlap*i.r/ (i.r + r) , most repulsion happens on the smaller circle
    ct[, repelLength := overlap*i.r/(i.r + r)]
    ct[, c("repel.x", "repel.y") := .( dx * repelLength/distance,
                                       dy * repelLength/distance) ]
    
    # means or sums???
    repelTotal <- ct[,  .(repel.x = sumFunction(repel.x), repel.y = sumFunction(repel.y)), by = gene]
    
    # update original table
    subTable[repelTotal, x := x + repel.x ,on = "gene"]
    subTable[repelTotal, y := y + repel.y ,on = "gene"]
    return (nrow(ct))  
  }
  
  # set radii in tsne.dt  
  if ( is.null(tsne.dt$radius) ){
    # bait radii, is only a fixed number
    tsne.dt[type == "bait", radius := baitRadius]
    
    # prey radii might be based on other columns in table
    if ("function" %in% class(preyRadius))
      tsne.dt[type == "prey", radius := preyRadius(.SD)]
    else
      tsne.dt[type == "prey", radius := preyRadius]
  }
  
  # build a table to work on
  subTable <- tsne.dt[, .(x = V1, y = V2, r = radius, gene)]
  
  for (i in 1:maxIterations){
    cat (i, " ") # function below will close the line with \n
    if (.oneIteration(subTable, sumFunction = sumFunction, overlapParam = overlapParam) == 0)
      break
  }
  
  # now update tsne.dt
  tsne.dt[subTable, c("x", "y") := .(i.x, i.y), on = "gene"]
}


moveLabelsToPerimeter <- function (moduleMedoids.dt, center = c(0,0), xlim = c(-60, 55), ylim = c(-75, 60)){
  
  # get current angles using some trigonetry
  angles <- moduleMedoids.dt[, atan2(y,x)] 
  spreadAngles <- seq(from = -pi, to = pi, length.out = length(angles) + 1)
  spreadAngles <- spreadAngles[rank(angles)]
  moduleMedoids.dt[, spreadAngles := spreadAngles]
  moduleMedoids.dt[, coss := cos(spreadAngles)]
  moduleMedoids.dt[, sins := sin(spreadAngles)]
  
  # y value intersect with x boundaries
  # for positive x
  moduleMedoids.dt[coss > 0, YatX.Boundary := sins * max(xlim)/coss ]
  # for negative x
  moduleMedoids.dt[coss < 0, YatX.Boundary := sins * min(xlim)/coss ]
  
  # x value intersect with y boundaries
  # for positive y
  moduleMedoids.dt[sins > 0, XatY.Boundary := coss * max(ylim)/sins ]
  # for negative y
  moduleMedoids.dt[sins < 0, XatY.Boundary := coss * min(ylim)/sins ]
  
  # in order so we can alternate the labels along the y boundaries
  setorder(moduleMedoids.dt, spreadAngles)
  yStepSize <- max(ylim) * 0.05
  
  # actual location depends on staying in bounds
  moduleMedoids.dt[XatY.Boundary %between% xlim,
                   c("perimeterX", "perimeterY") := .(XatY.Boundary, yStepSize * c(-1, 0,1)[(.I%%3)+1] + ifelse(sins < 0, ylim[1], ylim[2]))]
  
  moduleMedoids.dt[YatX.Boundary %between% ylim,
                   c("perimeterX", "perimeterY") := .( ifelse(coss < 0, xlim[1], xlim[2]), YatX.Boundary)]
  
  # cleanup
  moduleMedoids.dt[, c("coss", "sins", "YatX.Boundary", "XatY.Boundary") := NULL]
  
  invisible(moduleMedoids.dt)
  
}


cloudAsPolygon <- function (subData,  lims, color = NULL, gridLines = 200, level = 0.002){
  #subData = tsne.dt[go == "protein localization to nucleus" & nameMatch == TRUE, .(x,y)]
  k2d <- MASS::kde2d(subData$x, subData$y, lims = lims, n = gridLines)
  l <- contourLines(k2d)
  contour.dt <- lapply(l, function(ll) data.table(level = ll$level, x = ll$x, y= ll$y)) |> rbindlist(idcol = "contour")
  #contour <- contour.dt[level == 0.002]
  contour <- contour.dt[level == sort(unique(level))[2]]
  
  if (nrow(contour) == 0){
    contour <- contour.dt[level == sort(unique(level))[2]]
  }
  if (is.null(color) & !is.null(subData$nodeColor)){
    color <-  subData$nodeColor[1]
  }
  if (is.null(color))
    color <- "black"
  
  geom_polygon(data = contour, aes(x =x, y = y, group = contour), fill = NA, color = color) 
}

wrapStrings <- function (strings, width){
  lapply (stringi::stri_wrap (strings, width,
                              simplify = FALSE, whitespace_only = TRUE, cost_exponent = 2),
          paste0, collapse = "\n")
}


moduleNetworkView <- function(tsne.dt, edgeView.dt, clusterCloudExclude, moduleMedoids.dt, maxEdgeLength = 25, doClouds = TRUE, polygonCloud = FALSE,
                              emphasizeBaitBaitEdges = FALSE, baitLabels = TRUE){
  # MAKE the view
  p <- ggplot(tsne.dt, aes (x = x, y = y))
  
  
  # cluster clouds
  if (doClouds == TRUE){
    for (i in setdiff(unique(tsne.dt$clusterID), clusterCloudExclude)){
      subData <- tsne.dt[clusterID ==i & type == "prey" & nameMatch == TRUE]
      if(nrow(subData) >= 5){
        print (i)
        if (polygonCloud){
          p <- p + cloudAsPolygon(subData, lims = c(range(tsne.dt$x), range (tsne.dt$y)))
        }else {
          p <- p + stat_density_2d(geom = "polygon", data = subData, aes(alpha = ..level.., fill = nodeColor), show.legend = FALSE)
        }
      }
    }
  }else{
    p <- p + ggforce::geom_circle(aes(x0 = x, y0 = y, r = radius, color = nodeColor, fill = nodeColor), data = moduleMedoids.dt)
  }

  # setorder(edgeView.dt, -numPPI, na.last = TRUE)
  # edgeView.dt[, baitRank := 1:nrow(.SD), by  = Bait]
  # edgeView.dt[, moduleRank := 1:nrow(.SD), by  = Prey]
  # 
    
  # bait-prey edges
  p <- p +  
    geom_segment(aes(x,y, xend = xend,      yend = yend,      color = preyColor), alpha = 0.2, linewidth = 0.25, size= NULL,  data =  edgeView.dt[ ]) +
    
    #geom_segment(aes(x,y, xend = xend.stub, yend = yend.stub, color = preyColor), alpha = 0.5, linewidth = 0.3, size= NULL,  data =  edgeView.dt[ ]) +
    
    # bait module edges
    #geom_segment(aes(x,y, xend = xend,      yend = yend,      color = preyColor, linewidth = numPPI), alpha = 0.3,    data =  edgeView.dt[type == "module" & numPPI > 1 & (baitRank < 3 | moduleRank < 3)]) +
    scale_linewidth_binned(breaks = c(0,1,2,4,8,16,32,64), range = c(0.1, 4)) + 
    
    
    # bait - bait edges
    # ifelse(emphasizeBaitBaitEdges == TRUE,
            #geom_segment(aes(x,y, xend = xend, yend = yend, color = preyColor), alpha = 0.5, linewidth = 0.3, size= NULL,  data =  edgeView.dt[Prey %in% Bait]) +
    #        c()) +
    # 
    # baits
    geom_point(data = tsne.dt[type == "bait"], alpha = 0.9, size = 3, shape = 23, show.legend = FALSE,  mapping = aes(fill = nodeColor)) +
    scale_fill_identity( na.value= "grey30") + 
    
    # preys
    geom_point(data = tsne.dt[type == "prey"], alpha = 0.5, shape = 21, show.legend = FALSE, mapping = aes(fill = nodeColor, color = nodeColor, size = baitPerPrey)) +
    scale_color_identity(na.value = "grey70") +
    scale_size_area(max_size = 2) + 
    
    
    # module labels
    #ggrepel::geom_text_repel(aes(label  = wrapStrings(go, width = 22) ), data = moduleMedoids.dt[N >=5 & !clusterID %in% clusterCloudExclude], size = 2, max.overlaps = 100, show.legend = FALSE, fill = NA, bg.color = "white", bg.r = 0.15) +
    #, color = nodeColor 
    geom_text(aes(label  = wrapStrings(go, width = 18)), data = moduleMedoids.dt[N >=5 & !clusterID %in% clusterCloudExclude], size = 3,  show.legend = FALSE) +
    #geom_text(aes(label  = go, x = perimeterX, y = perimeterY, color = nodeColor), data = moduleMedoids.dt[N >=5 & !clusterID %in% clusterCloudExclude], size = 2, max.overlaps = 100, show.legend = FALSE, fill = NA, bg.color = "white", bg.r = 0.15) +
    #geom_segment(aes( x = perimeterX, y = perimeterY, xend = x, yend = y, color = nodeColor), data = moduleMedoids.dt[N >=5 & !clusterID %in% clusterCloudExclude], show.legend = FALSE, linetype = "dotted") +
    
    # bait labels
    #ifelse(baitLabels == TRUE,
           ggrepel::geom_text_repel(data = tsne.dt[type == "bait"], aes (label = gene ), size = 2, show.legend = FALSE, bg.color = "white", bg.r = 0.15, seed = 1)+

    coord_fixed(clip = "off", xlim = range (tsne.dt$V1) * 1.2, ylim = range(tsne.dt$V2) * 1.2) + # only works if range is negative/positive
    theme_void()
  
  return (p)
}


defineModulesFromDistance <- function (dsd.path, geneNamesOI = NULL,
                                       deepSplit = 0, minClusterSize = 5){
  dsd.dt <- fread (dsd.path)
  dsd.mat <- as.matrix(dsd.dt, rownames = "V1")
  
  if (!is.null(geneNamesOI)){
    # remove the 1 hop connectors, etc
    stringsInMat <- intersect(rownames(dsd.mat), geneNamesOI)
    dsd.mat <- dsd.mat[geneNamesOI, geneNamesOI]
  }
  
  dsd.dist <- as.dist(dsd.mat)
  hc.average <- hclust (dsd.dist, method= "average")
  dtc <- dynamicTreeCut::cutreeHybrid(hc.average, distM = dsd.mat, deepSplit = deepSplit, minClusterSize = 5)
  
  cluster.dt <- data.table (gene = hc.average$labels, clusterID = dtc$labels, core = dtc$cores)
  cluster.dt[clusterID == 0, clusterID := NA]
  cluster.dt[!is.na(clusterID), cluster := sprintf("clust.%03d", clusterID)]
  return (cluster.dt[])
}



labelModulesByEnrichment <- function (cluster.dt, gmt = NULL, numProcessors = 8){
  #source ("../../bp_utils/enrichmentTestFunctions.R")
  #gmt <- loadGmtFromBioconductor(ontology = "ALL", keyType = keyType)
  enrich.dt <- enricherOnGroups(cluster.dt, geneColumn = "gene", groupColumns = "cluster", numProcessors = numProcessors, term2gene.gmt = gmt, universe = unique(cluster.dt$gene))
  #fwrite (enrich.dt, ScriptAndDatedFileName("Enrichment_DTC_clusters.csv.gz"))

  enrich.dt[, jaccard := enricherJaccard(GeneRatio, BgRatio)]
  enrich.dt[, termScore := termScore(GeneRatio, BgRatio)]
  
  # 
  enrich.dt[, bgCount := as.integer (tstrsplit(BgRatio, "/")[[1]])]
  
  # top by p value
  setorder(enrich.dt, pvalue)
  # this takes the first wiht the greatest Count from among the top 8 by p value
  # clusterNames.dt <- enrich.dt[, .SD[1:8], by = cluster, .SDcols = c("ID", "pvalue", "GeneRatio", "Count", "bgCount")
  #                              ][, .SD[which(Count == max(Count))[1]], by = cluster
  #                                ][, .(cluster, name = paste0(ID, "; ", GeneRatio = GeneRatio, "(", bgCount, ")"), go = ID)]
  clusterNames.dt <- enrich.dt[, .SD[1], by = cluster
                            ][, .(cluster, name = paste0(ID, "; ", GeneRatio = GeneRatio, "(", bgCount, ")"), go = ID) ]  
  

  setorder (enrich.dt, -jaccard)
  clusterNamesJaccard <- enrich.dt[, .SD[1], by = cluster
                                   ][, .(cluster, name = paste0(ID, "; ", GeneRatio = GeneRatio, "(", bgCount, ")"), go = ID) ]  
  clusterNames.dt[clusterNamesJaccard, jaccard.name := i.name , on = "cluster"]
  clusterNames.dt[clusterNamesJaccard, jaccard.go := i.go , on = "cluster"]
  
  setorder (enrich.dt, -termScore)
  clusterNamesTermScore <- enrich.dt[, .SD[1], by = cluster
                                              ][, .(cluster, name = paste0(ID, "; ", GeneRatio = GeneRatio, "(", bgCount, ")"), go = ID) ]   
  clusterNames.dt[clusterNamesTermScore, termScore.name := i.name , on = "cluster"]
  clusterNames.dt[clusterNamesTermScore, termScore.go := i.go , on = "cluster"]

  
  return (list(clusterNames.dt = clusterNames.dt, enrich.dt = enrich.dt))  
}




buildModuleNetworkView <- function (data.list = NULL,
                                    baitPrey.dt = NULL,
                                    dsd.path = NULL,
                                    tsne.dt = NULL,
                                    clusterNames = NULL,
                                    enrich.dt = NULL,
                                    
                                    clusterColors = NULL,
                                    baitClusterColors = NULL,
                                    
                                    collapseModules = FALSE,
                                    
                                    maxEdgeLength = Inf, stubEdgeLength = 10,
                                    clusterCloudExclude = c( ),
                                    adjustModuleLabels = NULL,
                                    ...){
  
  # baitPrey.dt
  if(is.null(baitPrey.dt) && # short circuit is crucial
     is.null(baitPrey.dt <- data.list$baitPrey.dt)){
    stop("a baitPrey.dt is required")
  }
  
  #MAKE tsne.dt
  if(is.null(tsne.dt) && # short circuit is crucial
     is.null(tsne.dt <- data.list$tsne.dt)){
    stop("no code yet for generating a tsne")
  }
  ## baitsPerPrey
  if (is.null(tsne.dt$baitPerPrey)){
    bpp = baitPrey.dt[, .(baitPerPrey = .N), by = Prey]
    tsne.dt[bpp, baitPerPrey := baitPerPrey, on = c(gene = "Prey")]
  }
  
  # set x and y from V1/V2 if necessary
  if (!all(c("x", "y") %in% colnames(tsne.dt)) & collapseModules == FALSE){
    shiftNodes2AvoidOverlaps (tsne.dt)  # adds/changes the x,y column based on V1/V2
  } 
  
  
  
  # COLORS node colors
  if (is.null(tsne.dt$nodeColor)){
    if(is.null(clusterColors) && # short circuit is crucial
       is.null(clusterColors <- data.list$clusterColors)){
      
      clusterColors <- randomcoloR::distinctColorPalette(length(unique(tsne.dt$cluster)))
      names(clusterColors) <- sort (unique(tsne.dt$cluster))
    }
    tsne.dt[type == "prey", nodeColor  := clusterColors[cluster]]
    
    if (!is.null(tsne.dt$baitCluster)){
      if(is.null(baitClusterColors) && # short circuit is crucial
         is.null(baitClusterColors <- data.list$baitClusterColors)){
        baitClusterColors <- RColorBrewer::brewer.pal(length(unique(tsne.dt[!is.na(baitCluster)]$baitCluster)), "Dark2")
      }
      tsne.dt[type == "bait", nodeColor  := baitClusterColors[baitCluster]]
    }else tsne.dt[type == "bait", nodeColor = gray(0.2)]
    
  }

  # MODULES: center of modules for labels
  if (collapseModules == FALSE){
    moduleMedoids.dt <- tsne.dt[nameMatch == TRUE & type == "prey",
                                .(x = median(x), y = median(y), .N),
                                by = .(clusterID, cluster, name, go, nodeColor)]
    
  }else {
    moduleMedoids.dt <- tsne.dt[nameMatch == TRUE & type == "prey",
                                .(V1 = median(V1), V2 = median(V2), .N),
                                by = .(clusterID, cluster, name, go, nodeColor)]
    # radius based on N
    moduleMedoids.dt[, radius := sqrt(N)]
    
    # new nodes data.table with baits and modules
    tsneModules.dt <- rbind (tsne.dt[type == "bait"],
                             moduleMedoids.dt[, .(V1, V2, gene = cluster, nodeColor, radius, type = "module", go,clusterID, cluster)],
                             use.names = TRUE, fill = TRUE )
    original.tsne <- tsne.dt
    tsne.dt <- tsneModules.dt
    # now we can adjust to avoid overlaps
    shiftNodes2AvoidOverlaps (tsne.dt, overlapParam = 1.0)  # adds/changes the x,y column based on V1/V2
    # and reupdate moduleMedoids.dt with new x,y positions
    moduleMedoids.dt[tsne.dt[type == "module"], c("x", "y") := .(i.x, i.y), on = c(cluster = "gene")]
  }
  # compute a position on the perimeter
  moduleMedoids.dt <- moveLabelsToPerimeter(moduleMedoids.dt, xlim = range(tsne.dt$V1) * 1.1, ylim = range(tsne.dt$V2) * 1.2)
  
  if (!is.null(adjustModuleLabels)){
    adjust.dt <- data.table(go = names(adjustModuleLabels),
                            x.adjust = sapply(adjustModuleLabels, `[`, 1),
                            y.adjust = sapply(adjustModuleLabels, `[`, 2))
    moduleMedoids.dt[adjust.dt,
                     c("x", "y") := .(x + x.adjust, y + y.adjust)
                     , on = "go"]
    }

  
  # segment table (connect the dots in tsne.dt  )
  edgeView.dt <- tsne.dt[baitPrey.dt, .(Bait, Prey, x = x, y= y), on= c(gene = "Bait")
                         ][tsne.dt, .(Bait, Prey, x, y, xend = i.x, yend = i.y), on = c(Prey = "gene")
                           ][!is.na(Bait)]
  if (collapseModules){
    # edgeView.dt should only have bait-bait edges so far...
    # definte some bait to module edges
    baitModule.dt <- baitPrey2baitModule(original.tsne, baitPrey.dt)

    bm.edges <- tsne.dt[baitModule.dt, .(Bait, Prey, x = x, y= y, numPPI = i.N), on= c(gene = "Bait")
            ][tsne.dt, .(Bait, Prey, x, y, xend = i.x, yend = i.y, numPPI, type = "module"), on = c(Prey = "gene")
              ][!is.na(Bait)]
    edgeView.dt <- rbind (edgeView.dt, bm.edges, check.names = TRUE, fill= TRUE)
  }
  edgeView.dt[tsne.dt[type == "prey"], preyColor := nodeColor , on = c(Prey = "gene")]
  edgeView.dt[tsne.dt[type == "bait"], preyColor := "black",  on = c(Prey = "gene")]
  edgeView.dt[tsne.dt[type == "module"], preyColor := nodeColor,  on = c(Prey = "gene")]
  
  # define edge stubs so we can prune long edges   
  edgeView.dt[, length := sqrt( (xend-x)^2 + (yend-y)^2)]
  #, maxEdgeLength = 20, stubEdgeLength = 10
  #!Prey %in% Bait &
  edgeView.dt[length > maxEdgeLength, c("xend.stub", "yend.stub") := .(x + (xend-x)*stubEdgeLength/length,
                                                                       y + (yend-y)*stubEdgeLength/length)]
  edgeView.dt[length <= maxEdgeLength, c("xend.stub", "yend.stub") := .(xend, yend)]

  
  
  # turn off color for prey that don't match their cluster GO
  tsne.dt[type == "prey" & is.na(nameMatchCluster), nodeColor := "grey70"]
  
  
  
  # MAKE the view
  p <- moduleNetworkView(tsne.dt, edgeView.dt, clusterCloudExclude, moduleMedoids.dt, ...)
  print (p)
  
  # RETURN the passed-in and generated data
  if (collapseModules) tsne.dt <- original.tsne
  invisible(list (p = p,
                  tsne.dt = tsne.dt,
                  baitPrey.dt = baitPrey.dt,
                  
                  moduleMedoids.dt = moduleMedoids.dt,
                  
                  baitClusterColors = baitClusterColors,
                  clusterColors = clusterColors,
                  edgeView.dt = edgeView.dt))
}



