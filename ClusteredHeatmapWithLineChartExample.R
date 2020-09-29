library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)


# start with output from MSstats::DataProcess. In my case, mine is called 'dp.damgo.noBG' below

#####################
# prepare the data
#####################

# get a single quantity per protein/group
groupQuant <- MSstats::quantification(dp.morphine.noBG, type = "Group", format = "long")

setDT(groupQuant) # a memory efficient way to convert a data.frame to a data.table

# scaledIntensity has some advantages over LogIntensity:
#  1) Missing values can be assigned intensity = 0.0 with little impact on clustering
#  2) Linear scale is more intuitive for some people than a log scale
#  3) Low abundance proteins and high abundance proteins are treated equally
groupQuant[,scaledIntensity := 2^(LogIntensity - max(LogIntensity, na.rm=TRUE)), by = "Protein" ]

# covnert to wide format, then a matrix
quant.wide <- dcast (groupQuant, Protein~Group, value.var = "scaledIntensity")
quant.mat <- as.matrix(quant.wide, rownames = "Protein")

# this is important later for setting the missing values to gray in the heatmap
minNonMissing <- min (quant.mat, na.rm=TRUE)

# set missing values to 0.0
quant.mat[is.na(quant.mat)] <- 0.0



#################
# clustering
###################
# choose how many clusters you want, best is actually to try a few, here I just use 8

k <- 8

# kmeans depends on randomly choosing a starting point and iterating
# for repeatable results if/when I run this again, this initializes the random number generator 
set.seed(1)
km.clusters <- kmeans(quant.mat, centers=k, nstart=100, iter.max=20)


##################
#  Heatmap
#################

stopifnot(rownames(quant.mat) == names(km.clusters$cluster))
clusterColors <- brewer.pal(n = 12, name = "Paired")[km.clusters$cluster]
names(clusterColors) <- names(km.clusters$cluster)

cols <- brewer.pal(n = 12, name = "Paired")
names(cols) <- cols
clusterAnn <- rowAnnotation(cluster=clusterColors, col = list(cluster=cols))

topAnno <-HeatmapAnnotation(clusterCenters = anno_lines(t(km.clusters$centers), 
                                                        height=unit(1, "inches"), 
                                                        gp = gpar(col = cols)),
                            annotation_name_side = "left")

hms <- topAnno %v% Heatmap(quant.mat, name = "Scaled\nIntensity", 
                           cluster_columns = FALSE, show_row_dend = FALSE, show_row_names=FALSE,
                           row_split = km.clusters$cluster,left_annotation = clusterAnn,
                           col = circlize::colorRamp2(c(0,minNonMissing * 1e-3, 0.5, 1.0),
                                                      c("gray", "blue", "white", "red")))

draw(hms)



