






checkVersionDataProcessOutput <- function (mssQ){
  if (!"ProteinLevelData" %in% names(mssQ)  ){
    stop("Wrong version of MSstats. Source MSstats_Helper_Functions.R and try again")
  }
}


makeContrast.regEx <- function(mssQ, regEx){
  checkVersionDataProcessOutput(mssQ)
  columnNames <- as.character(levels(mssQ$ProteinLevelData$GROUP))
  
  positives = c()
  negatives = c()
  for (i in columnNames){
    for (j in columnNames){
      if (i != j){
        positives <- c(positives,i)
        negatives <- c(negatives,j)
      }
    }
  }
  contrasts <- data.table (positives, negatives)
  contrasts[,name := paste(positives,negatives, sep="-")]
  
  contrastNames <- character(0)
  for (re in regEx){
    contrastNames <- union(contrastNames, grep(re, contrasts$name, value = TRUE))
  }
  print (contrastNames)
  contrasts <- contrasts[name %in% contrastNames]
  print (contrasts)
  
  contrastMat <- matrix(rep(0,length(columnNames) * nrow(contrasts)), nrow=nrow(contrasts), ncol=length(columnNames),
                        dimnames = list(contrasts$name, columnNames))
  
  for (i in seq_len(nrow(contrasts))){
    #rows are named according to the positive condition:
    contrastMat[contrasts$name[i],contrasts$negatives[i]] = -1
    contrastMat[contrasts$name[i],contrasts$positives[i]] =  1
  }
  return (contrastMat)  
}




makeContrast.AllByAll <- function(mssQ){
  checkVersionDataProcessOutput(mssQ)
  columnNames <- as.character(levels(mssQ$ProteinLevelData$GROUP))
  
  positives = c()
  negatives = c()
  for (i in columnNames){
    for (j in columnNames){
      if(i < j){
        positives <- c(positives,i)
        negatives <- c(negatives,j)
      }
    }
  }
  contrasts <- data.table (positives, negatives)
  contrasts[,name := paste(positives,negatives, sep="-")]
  
  contrastMat <- matrix(rep(0,length(columnNames) * nrow(contrasts)), nrow=nrow(contrasts), ncol=length(columnNames),
                        dimnames = list(contrasts$name, columnNames))
  
  for (i in seq_len(nrow(contrasts))){
    #rows are named according to the positive condition:
    contrastMat[contrasts$name[i],contrasts$negatives[i]] = -1
    contrastMat[contrasts$name[i],contrasts$positives[i]] =  1
  }
  return (contrastMat)  
}


