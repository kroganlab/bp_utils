



#' specificity.matrix 
#' Given a matrix of protein observation (intensity or spectral count), compute a MIST-style specificity
#' A MIST-style specificity is sum(prey-in-bait.i)/sum(prey_in_all_runs)
#' Currently, this assumes a balanced design and no correction is made for different number of replicates per bait
#' 
#' @param mat with row names (preys) and column names (runs). This is assumed normalized/cleaned etc. No additional processing is done
#' @param columnKeys a data.table with columns `bait` and `column` which matches column names of `mat`. Or with first two columns
#'                   being equivalent mapping (but with any name)
#' @param sep character. If no columnKeys is passed, attempt to extract bait information from column names by separating 
#'            column names with this character.  `data.table::tstrsplit` is used for the separation
#' 
#' @param nameparts character. This defines the  names to apply to the result of splitting columnn names by `sep`
#' 
#' @return a matrix of specificity avlues with preys in rows and baits in column.  
#' 

specificity.matrix <- function(mat, columnKeys= NULL, sep = "_", nameParts = c("batch", "rep", "bait")){
  columnKeys <- .getColumnKeys(mat, columnKeys, sep, nameParts)
  
  specificity.mat <- sapply( unique(columnKeys$bait),
                             function(b)rowSums(mat[, columnKeys[b, column]])
  )/rowSums(mat)
  
  return(specificity.mat)
}

#' a utilitity function for handling `columnKeys` as passed in table or from column names of `mat`

.getColumnKeys <- function(mat, columnKeys= NULL, sep = "_", nameParts = c("batch", "rep", "bait")){
  # no columnKeys, try to extract bait/batch/rep info from column names
  if (is.null(columnKeys)){
    columnKeys <- data.table(column = colnames(mat))
    columnKeys[, (nameParts) := tstrsplit(column, sep)]  # () on nameParts needed to evaluate it to its value
  }
  
  # check columns are named or ordered as expected:
  if (!all (c("bait", "column") %in% colnames(columnKeys))){
    message ("The columnKeys table did not have columns bait and column. Using first two columns, which are named: ", colnames(columnKeys)[1:2])
    message ("If this is not correct, reorder or rename columns in columnKeys")
    columnKeys <- columnKeys[, 1:2]
    setnames(columnKeys, new= c("bait", "column"))
  }
  
  setkey(columnKeys, bait)
  return (columnKeys)
}



# 

permute.specificity.matrix <- function(mat, numPermutations = 1000, columnKeys= NULL, sep = "_", nameParts = c("batch", "rep", "bait"), numProc = 1, replace = FALSE){
  columnKeys <- .getColumnKeys(mat, columnKeys, sep, nameParts)
  
  # first calculate the specificity
  s <- specificity.matrix(mat, columnKeys)
  
  # get all unique prey/specificity/numRep pairs.  Each of these will gets its own p.value
  # TODO: include the replicate counts.  With current study, everything had 3, perfectly balanced, so no need to adjust
  if (length(unique(columnKeys[, .N, by = bait]$N)) > 1)
    stop("Unbalanced design detected. I don't handle different number of reps per bait yet")
  preySpecValues <- unique(melt(as.data.table(s, keep.rownames = TRUE),
                                id.vars = "rn", value.name = "specificity",
                                variable.name = "bait"
                           )[, .(rn, specificity)])
  setnames(preySpecValues, old = "rn", new = "prey")
  setkey(preySpecValues, prey)
  
  listOfDT <- pbapply::pblapply(rownames(mat), .permute1Prey, mat, preySpecValues, columnKeys, numPermutations, seed = 1, cl = numProc)
  return (listOfDT)
  
}




.permute1Prey <- function (prey, mat, preySpecValues, columnKeys, numPermutations,  seed = 1, replace = FALSE){
  library (data.table)
  # a table of bait to scores. 
  # I'll permute this. This ensures that baits with different numbers of reps are handled separately (but that's not handled yet)
  x.dt <- cbind(columnKeys,  # bait (and other things)
                mat[prey,])  # scores for 1 prey, specral counts for example
  setnames(x.dt, old = "V2", new = "score")
  x.dt[, ogScore := score] # keep a copy of the scores that we sample with replacement (avoid random walk)
  
  actualSpecificities <- preySpecValues[prey]
  actualSpecificities[, countPermAsGood := 0]
  actualSpecificities[, totalPermCount := 0]
  
  
  set.seed(seed)
  for (i in 1:numPermutations){
    x.dt[, score := sample (ogScore, replace = replace)]
    rowTotal <- sum(x.dt$score)
    specificities <- x.dt[, .(specificity = sum(score)/rowTotal), by = bait]
    
    
    # per specificity actually observed, count how many of these random are as good as it.
    # TODO: make this num-rep specific  4 reps are more likely high scoring than 2, e.g....
    actualSpecificities[,countPermAsGood := countPermAsGood +  sum(specificities$specificity >= specificity), by =  specificity]
    actualSpecificities[,totalPermCount := totalPermCount + nrow(specificities)]
    
  }
  actualSpecificities[, pvalue := countPermAsGood/totalPermCount]
  return(actualSpecificities)
}





