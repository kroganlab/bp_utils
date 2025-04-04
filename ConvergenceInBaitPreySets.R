
#' @returns a data.table of bait,prey,bait_prey

loadBioPlex <- function(path = "~/Downloads/BioPlex_3.0_293T_DirectedEdges.tsv"){
  bp <- fread (path)
  bp[, c("bait", "prey") := .(`Bait Symbol`, `Prey Symbol`)]
  bp[, bait_prey := paste0(bait, "_", prey)]
  return (bp)
}

#' @param desiredBaitPoolSize for selecting baits with same-sized prey sets, this 
#'        is the desired set of baits to sample among. For prey-set sizes that have few baits
#'        the windowSize determines the steps above and below to go to reach a number of baits
#'         at least this number
#' @param bp bioplex bait prey table
#' @returns a data.table of countPrey,numBaitsCP,windowSize

buildBioPlexPreySizeWindowTable <- function(bp, desiredBaitPoolSize = 100){
  bpPreySizes <- bp[, .(countPrey = length(unique(prey))), by = bait]
  bpPreySizeTable <- bpPreySizes[, .(numBaitsWithCP = .N),
                                 by = countPrey
                                 ][data.table(countPrey = 1:max(countPrey)),
                                   , on = "countPrey"][order(countPrey)]
  bpPreySizeTable[is.na(numBaitsWithCP), numBaitsWithCP := 0]
  
  # I want each size to have at least 100 baits to sample from, so window up and down until I get the 100
  .windowSize <- function(i, nn, threshold = 100){
    step = 0
    while(step <length(nn)){
      windowIndex <- (i-step):(i+step)
      if (sum(nn[windowIndex], na.rm = TRUE) > threshold){
        break
      }
      step <- step + 1
    }
    return (step)
  }
  
  bpPreySizeTable[, windowSize := sapply(1:.N, 
                                         .windowSize, 
                                         nn = bpPreySizeTable$numBaitsWithCP,
                                         threshold = desiredBaitPoolSize)]
  
  return(bpPreySizeTable[])
}



#' @param bait1 character name matching one of the columns of matrix
#' @param bait2 character name matching one of the columns of matrix
#' @param matrix numerical matrix, rows are preys and columns are baits
#'              expect a value of 0,1 or FALSE/TRUE indicating bait-prey relationships
#' @param proteomeSize integer giving the "background" size. What is the total set of 
#'                     proteins that we are selecting from? (This defines how large is the
#'                      "neither" region of hte contigency table.) The default number 11169
#'                      is from hek293tProteome_BekkerJenwen2017_ASDppiPrey.csv given to me
#'                      by Belinda Wang, as a count of HEK expressed proteins
#'                      `scan("/Users/ben/Downloads/hek293tProteome_BekkerJenwen2017_ASDppiPrey.csv", what = character()) |> unique() |> length()`
#'                      
#' @returns the p.value of a single fishers.test
#' 
onePairwiseFisher <- function (bait1, bait2, matrix, proteomeSize = 11169 ){
  both <- sum(rowSums(matrix[, c(bait1, bait2)]) == 2)
  firstOnly  <- sum(matrix[, c(bait1)]) - both
  secondOnly <- sum(matrix[, c(bait2)]) - both
  neither <- proteomeSize - both - firstOnly - secondOnly

  fisher.test (matrix(c(both, firstOnly, secondOnly, neither), nrow = 2), alternative = "greater")$p.value
}


#' @param returns a data.table bait1, bait2, fisherP of all nonredundant pairs of bait1,bait2
#'                the fisher test tests for non-random overlap between two prey sets
#'                background size matters...

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
    baitSet2 < setdiff(baitSet2, baitSet1)
  }
  baitPreyTable2 <- baitPreyTable2[bait %in% baitSet2]
  baitPreyTable1 <- baitPreyTable1[bait %in% baitSet1]
  baitPreyTable <- rbind (baitPreyTable1, baitPreyTable2)
  
  baitPreyMat <- dcast (baitPreyTable[bait %in% c(baitSet1, baitSet2),],
                        prey~bait, value.var = "bait", fun.aggregate = length) |>
    as.matrix(rownames= "prey")
  baitPreyMat[baitPreyMat > 1] <- 1
  
  allByAll <- data.table(bait1 = intersect(baitSet1,colnames(baitPreyMat)))[, .(bait2 = intersect(baitSet1,colnames(baitPreyMat))), by = bait1] 
  allByAll[ , fisherP := onePairwiseFisher(bait1, bait2, baitPreyMat), by = .(bait1, bait2)]
  return (allByAll)
  
}

#' @param baitPreySizes, a table of original baits, including columns bait,countPrey,windowRadius
#' @param bpPreySizeByBait, a table describing baits in bioplex, with columns bait,countPrey
#' @param excludeBaits, a list of baits in bioplex to exclude from sampling
#' @returns the input data.table baitPreySizes, changed in place, that includes a sampled value in (possibly new) column bpBaitSample
#'          sampled to match sizes given by countPrey and windowRadius
# used by doAndCountFisherTests
# randomly samples bpPreySizes, one per row in pcmiPreySizes, based on the prey size and window radius columns in pcmiPreySizes
doOneSample <- function(baitPreySizes, bpPreySizeByBait, excludeBaits = c()){
  baitPreySizes[, bpBaitSample := NA_character_]
  setkey(baitPreySizes, bait)
  for (b in baitPreySizes$bait){
    r <- baitPreySizes[b, windowRadius]
    n <- baitPreySizes[b, countPrey]
    baitPreySizes[b, bpBaitSample := sample(setdiff(bpPreySizeByBait[countPrey >= n - r & countPrey <= n+r, bait],
                                                    c(baitPreySizes$bpBaitSample, excludeBaits)),1)]
    
  }
}


#' on random background, do fisher tests for a random sampling of bioplex that matches the prey-size distribution of baitPreySizes
#' #' @returns A single number, range 0-1, the portion of fisher tests for all bait-bait pairwise comparisons that produce a p.value < 0.05
doAndCountFisherTests <- function (baitPreySizes, bpPreySizeByBait, baitPreyTable){
  doOneSample(baitPreySizes, bpPreySizeByBait) # baitPreySizes is changed in place: column bpBaitSample has the sample
  aba <- doFisherTests(baitPreySizes$bpBaitSample, baitPreyTable )
  sum(aba$fisherP < 0.05)/length(aba$fisherP)
}

# on random background, do fisher tests for a random sampling of bioplex that matches the prey-size distribution of pcmiPreySizes
#' #' @returns A single number, range 0-1, the portion of fisher tests for all pairwise comparisons that produce a p.value < 0.05
doAndCountFisherTestsTwoGroups <- function (set1PreySizes, set2PreySizes, bpPreySizes, baitPreyTable){
  doOneSample(set1PreySizes, bpPreySizes) # set1PreySizes is changed in place: column bpBaitSample has the sample
  doOneSample(set2PreySizes, bpPreySizes, excludeBaits = set1PreySizes$bpBaitSample) # set2PreySizes is changed in place: column bpBaitSample has the sample
  aba <- doFisherTestsBetweenGroups(set1PreySizes$bpBaitSample,set2PreySizes$bpBaitSample, baitPreyTable, baitPreyTable )
  sum(aba$fisherP < 0.05)/length(aba$fisherP)
}



#' @returns A vector of portion of p values < 0.05 on a permutation of bioplex to match baitPrey.dt
doPermutedFisherTests <- function(baitPrey.dt,  bioplex, desiredBaitPoolSize = 100, numProcessors = 1){
  kPreySizes <- baitPrey.dt[, .(countPrey = length(unique(prey))), by = bait]
  
  bpPreySizeByBait <- bioplex[, .(countPrey = length(unique(prey))), by = bait]
  bioplexPreySizeTable <- buildBioPlexPreySizeWindowTable(bioplex, desiredBaitPoolSize = desiredBaitPoolSize)

  kPreySizes[bioplexPreySizeTable, windowRadius := i.windowSize, on = "countPrey"]
  
  # doAndCountFisherTests <- function (baitPreySizes, bpPreySizeByBait, baitPreyTable){
  countPlt0.05.kaushikLike <- pbapply::pbsapply(1:1000,
                                                function(i, ...)doAndCountFisherTests(...),
                                                baitPreySizes = kPreySizes,
                                                bpPreySizeByBait = bpPreySizeByBait,
                                                baitPreyTable = bioplex,
                                                cl = numProcessors)
}



doConvergenceTestsOneDataset <- function(baitPrey.dt, bioplex, desiredBaitPoolSize = 100, numProcessors = 1){
  
  # actual
  baitPray.aba <- doFisherTests(unique(baitPrey.dt$bait), baitPrey.dt[, .(bait, prey)] )
  actualPortionPlt0.05 = sum(baitPray.aba$fisherP < 0.05)/length(baitPray.aba$fisherP)
  
  
  message ("Convergence rate in actual set among ", length(unique(baitPrey.dt$bait)), " baits is ", signif(actualPortionPlt0.05, 2))
  
  # random sets from bioplex
  
  randomPlt0.05 <- doPermutedFisherTests(baitPrey.dt, bioplex, desiredBaitPoolSize = desiredBaitPoolSize, numProcessors = numProcessors)
  
  return (list(actual = actualPortionPlt0.05, randoms = randomPlt0.05))
}





#' @param allRandom Output of ``, a data.table of CD from permutations of bioplex. Needs columns study, CD.05
#' @param actual Output of ``, a data.table of CD from actual prey-prey overlap fisher's tests.
#' 

doConvergencePlot <- function(allRandom, actual){
  allRandom[, source := "BioPlex sample"]
  actual[, source := "actual"]
  p <- ggplot(allRandom, aes(x = study, y= CD.05, color = source)) + 
    geom_violin(show.legend = FALSE) + 
    ggforce::geom_sina( alpha = 0.2, size  = 0.5, stroke = NA, shape = 16) +
    #geom_text(data = actual, size = 6, label = "\u2605", color = "blue",family = "Arial Unicode MS") +
    geom_point(data = actual, size = 3) +
    theme_classic()+
    theme(legend.position = "bottom") +
    theme(axis.title.x = element_blank()) + 
    scale_x_discrete(name = "") +
    #scale_y_continuous(limits = c(0, 0.15)) + 
    scale_color_manual(values = c(`BioPlex sample` = "grey", actual = "blue"), name = "") + 
    rotate.x.axis.text
  p
  
}