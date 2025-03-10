# Functions that help measure significance of overlaps between the prey sets of baits within and between APMS experiments.
# Developed for the PCMI ASD project originally





#  File from Belinda, all HEK detectable proteins + PCMI prey
#scan("/Users/ben/Downloads/hek293tProteome_BekkerJenwen2017_ASDppiPrey.csv", what = character()) |> unique() |> length()
# 11169, use this number below as the background
# given two sets, how significant is their overlap
onePairwiseFisher <- function (bait1, bait2, matrix, bgSize = 11169){
  both <- sum(rowSums(matrix[, c(bait1, bait2)]) == 2)
  firstOnly  <- sum(matrix[, c(bait1)]) - both
  secondOnly <- sum(matrix[, c(bait2)]) - both
  neither <- bgSize - both - firstOnly  -secondOnly
  #neither <- nrow(matrix) - both- firstOnly- secondOnly
  
  fisher.test (matrix(c(both, firstOnly, secondOnly, neither), nrow = 2), alternative = "greater")$p.value
}


# within a single set of bait-preys (a single apms experiment)
# does all by all pairwise tests, calcluating a p value for significance overlap
doFisherTests <- function(baitSet, baitPreyTable){
  
  baitPreyMat <- dcast (baitPreyTable[bait %in% baitSet,], prey~bait, value.var = "bait", fun.aggregate = length) |> as.matrix(rownames= "prey")
  baitPreyMat[baitPreyMat > 1] <- 1
  allByAll <- data.table(bait1 = colnames(baitPreyMat))[, .(bait2 = colnames(baitPreyMat)), by = bait1][bait1 < bait2]
  allByAll[ , fisherP := onePairwiseFisher(bait1, bait2, baitPreyMat), by = .(bait1, bait2)]
  return (allByAll)
}

# like doFisherTests, but works only by comparing baits between two bait sets (i.e. two apms experiments)
doFisherTestsBetweenGroups <- function (baitSet1, baitSet2, baitPreyTable1, baitPreyTable2){
  if (length(intersect(baitSet1, baitSet2)) > 0){
    warning("Matching baits in baitSet1, baitSet2. Removing them from baitSet2")
    baitSet2 <- setdiff(baitSet2, baitSet1)
  }
  baitPreyTable2 <- baitPreyTable2[bait %in% baitSet2]
  baitPreyTable1 <- baitPreyTable1[bait %in% baitSet1]
  baitPreyTable <- rbind (baitPreyTable1, baitPreyTable2)
  
  baitPreyMat <- dcast (baitPreyTable[bait %in% c(baitSet1, baitSet2),],
                        prey~bait, value.var = "bait", fun.aggregate = length) |>
    as.matrix(rownames= "prey")
  baitPreyMat[baitPreyMat > 1] <- 1
  
  allByAll <- data.table(bait1 = intersect(baitSet1,colnames(baitPreyMat)))[, .(bait2 = intersect(baitSet2,colnames(baitPreyMat))), by = bait1][bait1 != bait2]
  allByAll[ , fisherP := onePairwiseFisher(bait1, bait2, baitPreyMat), by = .(bait1, bait2)]
  return (allByAll)
  
}


# on random background, do fisher tests for a random sampling of bioplex that matches the prey-size distribution of pcmiPreySizes
doAndCountFisherTests <- function (pcmiPreySizes, bpPreySizes, baitPreyTable){
  # todo: process the two intput bait-prey tables to make pcmiPreySizes and bpPreySizes
  doOneSample(pcmiPreySizes, bpPreySizes) # pcmiPreySizes is changed in place: column bpBaitSample has the sample
  aba <- doFisherTests(pcmiPreySizes$bpBaitSample, baitPreyTable )
  sum(aba$fisherP < 0.05)/length(aba$fisherP)
}

doAndCountFisherTestsTwoGroups <- function (set1PreySizes, set2PreySizes, bpPreySizes, baitPreyTable){
  doOneSample(set1PreySizes, bpPreySizes) # set1PreySizes is changed in place: column bpBaitSample has the sample
  doOneSample(set2PreySizes, bpPreySizes, excludeBaits = set1PreySizes$bpBaitSample) # set2PreySizes is changed in place: column bpBaitSample has the sample
  aba <- doFisherTestsBetweenGroups(set1PreySizes$bpBaitSample,set2PreySizes$bpBaitSample, baitPreyTable, baitPreyTable )
  sum(aba$fisherP < 0.05)/length(aba$fisherP)
}


# used by doAndCountFisherTests
# randomly samples bpPreySizes, one per row in pcmiPreySizes, based on the prey size and window radius columns in pcmiPreySizes
doOneSample <- function(preySizes, bpPreySizes, excludeBaits = c()){
  preySizes[, bpBaitSample := NA_character_]
  # setkey(preySizes, bait)
  # for (bb in preySizes$bait){
  #   r <- preySizes[bb, windowRadius]
  #   n <- preySizes[bb, countPrey]
  #   preySizes[bb, bpBaitSample := sample(setdiff(bpPreySizes[countPrey >= n - r & countPrey <= n+r, bait],
  #                                                   c(preySizes$bpBaitSample, excludeBaits)),1)]
  #   
  # }
  
  for (i in 1:nrow(preySizes)){
      r <- preySizes[i, windowRadius]
      n <- preySizes[i, countPrey]
      preySizes[i, bpBaitSample := sample(setdiff(bpPreySizes[countPrey >= n - r & countPrey <= n+r, bait],
                                                      c(preySizes$bpBaitSample, excludeBaits)),1)]
    
  }
}


