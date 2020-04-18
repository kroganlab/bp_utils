
# Not a finished script.  This was originally used for mapping monkey to human.
# with a little editing in phSiteToHumanSite_expanded is should be convertible easily to mouse or...

library (data.table)

# this handles splitting the site off of a protein, sepaarated by a single underscore
# it should handle the possibility that there are other underscores in a protein name
# it will not handle the case where there are two sites  (other than simply splitting off the last one)
# example P789ASDF_S123_T125 will not split off both sites
splitProtein_Site <- function (protein_site){
  fullSplit <- strsplit (protein_site, "_")
  rejoinProteinSplit <- lapply (fullSplit, 
                                FUN = function(ss)
                                {
                                  lastI <- length(ss)
                                  return(c(paste(ss[1:(lastI-1)], collapse="_"), ss[lastI]))
                                })
  return (transpose(rejoinProteinSplit))
}

#expects a data.table
#returns a new data.table, probably with less rows if there are any multi-site 
packMultipleSitesPerProtein <- function (siteTable, byColumn="Protein", protColumn = "human_uniprot", aaColumn = "human_aa", posColumn="human_pos"){
  setDT(siteTable)
  setnames(siteTable, old = c(byColumn, protColumn, aaColumn, posColumn), new = c("byColumn", "protColumn", "aaColumn", "posColumn"))
  tryCatch({
    res <-  siteTable[,fifelse(all(is.na(protColumn)),
                               NA_character_,
                               unique(.SD[,.(protColumn, aaColumn, posColumn)])[,paste(protColumn,
                                                                                       paste0(aaColumn, posColumn),
                                                                                       sep="_",
                                                                                       collapse=";")]),
                      by=.(byColumn, protColumn)]
    
  }, finally={
    setnames(siteTable, new = c(byColumn, protColumn, aaColumn, posColumn), old = c("byColumn", "protColumn", "aaColumn", "posColumn"))
  })
  setnames(res, new=c(byColumn, protColumn, "human"))
  return (res)
}

# removeNonHumanGaps  -- gaps get the number of  neighboring non-gap.  So a non-human site that is on the edge
# of a deletion will get mapped to two different human sites, but one will be labeled a gap.
# 
# 
phSiteToHumanSite_expanded <- function (Protein, mapTable, removeNonHumanGaps = TRUE, gaps=c("-")){
  #names in mapTable
  #TODO figure these out automatically or pass them in 
  otherCols <- list(pos = "vervet_pos", aa = "vervet_aa", uniprot="vervet_uniprot") #change these to mouse...
  humancols <- c(pos = "human_pos", aa = "human_aa", uniprot="human_uniprot")
  
  #prepare the Protein column by expanding  prot1_S12;prot1_T14;prot1_T16
  # results will be additional columns like prot1 S 12 
  protExp <- data.table(Protein = unique(Protein))
  protExp <- protExp[,.(singleSite = unlist(strsplit(Protein, split=";"))), by = Protein]
  protExp[,c("singleProtein", "aaAndPos"):=splitProtein_Site(singleSite)]
  protExp[,c("aa", "pos") := .(substr(aaAndPos, 1,1), as.integer(substr(aaAndPos, 2, str_length(aaAndPos))))]
  
  #now join to the map table on othercols
  mergedSingleSites <- merge (protExp, mapTable, by.x = c('singleProtein', 'pos'), by.y = c(otherCols$uniprot, otherCols$pos), all.x=TRUE)
  
  # we my have picked up some gaps, remove them if necessary, also note the non-symmetry here
  # we only have to deal with those that align after a site of interest  
  # S-V  align SSV  gets here, V-S align VSS does not because the gap gets the position number of the preceding AA not the following.
  if (removeNonHumanGaps){
    otherGaps <- mergedSingleSites[[otherCols$aa]] %in% gaps
    if (any(otherGaps)){
      message ("NOTE: There are ", sum(otherGaps), " human phospho sites inserted and neighboring after a site of interest. These will be ignored.")
      mergedSingleSites <- mergedSingleSites[!otherGaps]
    }
  }
  
  return (mergedSingleSites)
  
}

phSiteToHumanSite <- function (Protein, mapTable){
  mergedSingleSites <- phSiteToHumanSite_expanded(Protein, mapTable)
  packedHumanSites <- packMultipleSitesPerProtein (mergedSingleSites)
  if(nrow(packedHumanSites) > length(Protein)){
    message ("WARNING: phSiteToHumanSite will ignore all but one human protein mapped to your proteins of interest. ", nrow(packedHumanSites) - length(Protein), " cases ignored")
  }
  return (packedHumanSites[match(Protein, packedHumanSites$Protein), ]$human )
}


# results below is the output of a MSstats with Protein column like A134ASD_S123;A134ASD_S125 for example...
# requires data.table setDT(results)

sampleUsage <- function (results, siteMapperFile = "Phospho_Sites_Full_Mapping_GreenMonkey_Human.csv.gz", uniprotIDMapper = "HUMAN_9606_idmapping.dat.gz"){
  expandToHumanSites <- phSiteToHumanSite_expanded (results$Protein, mapTable = fread (siteMapperFile), )
  expandToHumanSites <- merge (results, expandToHumanSites, by = "Protein", all.x=TRUE, allow.cartesian=TRUE)
  
  packedHumanSites <- packMultipleSitesPerProtein (expandToHumanSites)
  
  expandToHumanSites <- merge (expandToHumanSites, packedHumanSites, by = c("Protein", "human_uniprot"))
  
  setnames(expandToHumanSites, old="Protein", new="msstats.Protein")

  expandToHumanSites <- merge (expandToHumanSites, humanIdMapper, by = "human_uniprot", all.x=TRUE)
  
  if (!is.null(uniprotIDMapper)){
    humanIdMapper <- fread (uniprotIDMapper, header=FALSE)[V2 == "Gene_Name", .(V3,V1)]
    names(humanIdMapper) <- c("Gene_Name", "human_uniprot")
    humanIdMapper <- humanIdMapper[,Gene_Name[1], by = human_uniprot] #get only 1 name per uniprot
    expandToHumanSites <- merge (expandToHumanSites, humanIdMapper, by = "human_uniprot", all.x=TRUE)
  }
  expandToHumanSites
}