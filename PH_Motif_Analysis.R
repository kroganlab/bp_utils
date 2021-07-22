




#' This will expand a results table to include multiple rows for proteins like A2A5R2_S218;A2A5R2_S227
#' @param results a data.table like fread (results.txt). Protein should be in format uniprot_S123. Only required column is Protein
#' @param fastaFilePath a character 
#' @param radius integer. subseq will be ph-radius:ph+radius, or length 2*radius +1
addSeqNeighborhoodToResults <- function (results, fastaFilePath, radius = 6){
  seqs <- Biostrings::readAAStringSet(fastaFilePath)
  #reduce names to just hte uniprots, second item split on "|"
  names(seqs) <- tstrsplit(names(seqs), split = "\\|")[[2]]
  
  # per unique protein expand to individual site
  expandProteins <- results[,.(singleSite = unlist(strsplit(Protein, split = ";"))), by = .(Protein)]
  expandProteins[, c("uniprot", "site") :=  tstrsplit(singleSite, split = "_")]
  expandProteins[, site := as.integer(gsub ("[STY]", "", site))]
  
  # lookup subseq per unique site
  sites <- unique(expandProteins[,.(uniprot,site)])
  sites[,seqLength := sapply (uniprot, function(up)length(seqs[[up]]))]
  sites[, siteCenteredSubSeq := as.character(Biostrings::subseq(seqs[uniprot],
                                                                start = min(max (1, site-radius), seqLength),
                                                                end = min (seqLength, site+radius))),
        by = .(uniprot, site)]
  
  #merge protein and site info
  expandProteins <- expandProteins[sites, , on= c("uniprot", "site")]
  
  # collapse common substrings within Protein (usually only in case of multi-protein peptides) and merge to results:
  results.withSites <- results[unique(expandProteins [, .(Protein, siteCenteredSubSeq)]),
                               , on = "Protein", allow.cartesian = TRUE]
  
  return (list(sites, results.withSites))
}






