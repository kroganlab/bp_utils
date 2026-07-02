

cleanUpIsoformLabels <- function(fullProteinNames, nameSep = ";"){
  mapper <- data.table (fullIsoProtein = unique(fullProteinNames)) [, .(singleIsoProtein  = unlist(strsplit(fullIsoProtein , nameSep))), by = "fullIsoProtein"]
  mapper[, singleCanonProtein := gsub ("-[0-9]+$", "", singleIsoProtein)]
  mapper[, isCanonical :=  singleCanonProtein == singleIsoProtein]
  mapper[, anyCanonicalPerSingle := any(isCanonical), by= .(fullIsoProtein, singleCanonProtein)]
  
  mapper[anyCanonicalPerSingle == TRUE, singleFixedProtein := singleCanonProtein[1], by= .(fullIsoProtein, singleCanonProtein)]
  mapper[anyCanonicalPerSingle == FALSE, singleFixedProtein := paste0(sort(unique(singleIsoProtein)), collapse = ";"), by= .(fullIsoProtein, singleCanonProtein)]
  
  mapper[, fixedProtein := paste0(sort(unique(singleFixedProtein)), collapse= ";"), by= .(fullIsoProtein)]

  unique(mapper[, .(fullIsoProtein, fixedProtein)])  
  
  
  
  # canonicalOnly <- mapper[isCanonical == TRUE, .(fullCanonProtein = paste0(sort(singleCanonProtein), collapse = ";")), by= .(fullIsoProtein)]
  # mapper[, anyCanonical := any(isCanonical), by = fullIsoProtein]
  # noCanonical <- mapper[anyCanonical == FALSE,  .(fullProtein = paste0(sort(unique(singleIsoProtein)), collapse= ";")), by = .(fullIsoProtein)]
  # 
  # cleanMapper <- rbindlist(list(canonical = canonicalOnly, nonCanonical = noCanonical), idcol = "type", use.names = FALSE)
  # setnames(cleanMapper, old = "fullCanonProtein", new = "isoCleanProtein")
  # return (cleanMapper)
}

