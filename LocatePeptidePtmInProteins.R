library (data.table)
# requires stringr, RCurl, seqinr



convertMassModificationFormat <- function(specModSequence, mods=c("PH",  "CAM", "MOX", "NAC"), keepOnly = NULL, removeAll = FALSE){
  result <- specModSequence
  massFormats <- c(MOX = "(.)\\[15.9949\\]",
                   CAM = "(.)\\[57.0215\\]",
                   PH = "(.)\\[79.9663\\]",
                   NAC = "n\\[42.0106\\](.)",
                   ANY = "n?\\[[0-9.]+\\]")  # used for matching all or all but keepOnly
  
  
  artmsFormats <- list (PH='\\1\\(ph\\)',
                        UB='\\1\\(gl\\)',
                        CAM = '\\1\\(cam\\)',
                        MOX = '\\1\\(ox\\)',
                        NAC = '\\1\\(ac\\)',
                        ANY = '')  # for removal purposes
  
  
  if(!is.null(keepOnly) | removeAll == TRUE){
    mods <- c(keepOnly, "ANY")
  }
  
  # mass of UB is not yet known by me...
  for (mod in mods){
    if (mod %in% names(massFormats)){
      replace <- artmsFormats[[mod]]
      if (!is.null(keepOnly) && keepOnly != mod)
        replace <- ""
      result <- gsub(massFormats[[mod]], replace, result)
    }else (message("I don't yet know how to deal with mod: ", mod))
  }
  
  missedMods <- grep ("\\[.*\\]", result, value=TRUE)
  if (length(missedMods) > 1){
    #regex is: prefix([anything.not.`[`])suffix
    missedMasses <- unique(gsub (".*(\\[[^\\[]*\\]).*", "\\1", missedMods))
    message (length(missedMods), " peptides with unexpected modifications: ", paste0(missedMasses, collapse = ", "))
  }
  return (result)
}




mapSites.fragPipeFormat <- function(phosphoTableMSstats, site = "PH"){
  
  
}





#' this is unusually slow
#' it seems some connections to uniprot hang and fail
downloadUniprotSequences <- function(uniprots, numConnections = 1){
  fastaURL <- paste0('https://www.uniprot.org/uniprot/', uniprots, '.fasta')
  cat (sprintf("Downloading %d sequences directly from uniprot\n", length(uniprots)))
  cat ("This might hang  for an unusually long time, and maybe some uniprots will fail to download (with a message).")
  cat ("If so, and you need this functionality, bug or bribe me to fix it. -BP\n")
  
  
  #fastaTxt <- RCurl::getURLAsynchronous(url = fastaURL, perform = numConnections)
  fastaTxt <- RCurl::getURL(url = fastaURL, perform = numConnections)
  
  cat ("done\n")
  fasta.dt <- as.data.table(transpose(stringr::str_split(fastaTxt, n =2, pattern = "\n")))    
  setnames(fasta.dt, new = c("header", "sequence"))
  fasta.dt <- fasta.dt[grepl ("^>", header),]
  fasta.dt[, header := gsub("^>", "", header)]
  fasta.dt[, sequence := gsub ("\n", "", sequence)]
  fasta.dt[, c("db", "uniprot", "uniprotName") := tstrsplit (header, split = "[ |]", keep = 1:3)]
  numSuccess <- nrow(fasta.dt)
  if (numSuccess < length(uniprots) )
    message ("Failed to download ", length(uniprots) - numSuccess, " uniprots: ", setdiff(uniprots, fasta.dt$uniprot))
  return (fasta.dt[])
}



fastaFileToTable <- function(filePath){
  seqs <- seqinr::read.fasta(filePath, seqtype = "AA", 
                             as.string = TRUE, set.attributes = FALSE)
  dt <- data.table (header = names(seqs), sequence = unlist(seqs))
  dt[, c("db", "uniprot", "uniprotName") := tstrsplit(header, split = "\\|", keep = 1:3)]
  return (dt[])
}





#' @param proteinNames character vector, each in format like `Q9Y6Y0_SLS[79.9663]FEMQQDELIEKPMSPMQYAR`
#' @param ptmType  one of "PH" or in the future "UB", but for now only "PH"
#' @param fastaFile path to the fasta file, uniprot format, that should match the first field in protein names
#' @param downloadFromWeb logical TRUE = attempt to download uniprots not found in fasta file from web.

fragPipePTMFormat2SiteFormat <- function(proteinNames, ptmType = "PH", fastaFile = NULL, downloadFromWeb = TRUE){
  mapper <- data.table(fragPipeProteinName = unique(proteinNames))
  
  mapper[, c("uniprot", "modPeptide") := tstrsplit(fragPipeProteinName, "_")]
  
  mapper[, modPeptideLC.PTM := convertMassModificationFormat(modPeptide, keepOnly = ptmType)]
  
  formats <- c(PH = "(ph)") # depends on function convertMassModificationFormat
  ptmFormat <- formats[ptmType]
  
  # first locate ptm in peptide. Mulitple matches per peptide
  mapper[, ptmPositions := stringr::str_locate_all(modPeptideLC.PTM, stringr::fixed (ptmFormat))]
  
  # we have a list of matrices. one matrix per row. expand as many rows as there are sites in the row
  mapper.expanded <- mapper[, .(pepSitePTM = ptmPositions[[1]][,"start"] -
                                  # subtract 1 because the (ph) labels the preceding AA, and subtract an additional 4 per each preceding mod
                                  seq(from = 1, by =  stringr::str_length(ptmFormat), length.out =nrow(ptmPositions[[1]]))  ),
                            by = .(fragPipeProteinName, uniprot, modPeptide)]
  
  mapper.expanded[, cleanPeptide := convertMassModificationFormat (modPeptide, removeAll = TRUE)]
  mapper.expanded[, aaPTM := substr(cleanPeptide, pepSitePTM, pepSitePTM)]
  
  
  ## load uniprots from file and web
  fasta.dt <- fastaFileToTable (fastaFile)
  uniprots <- unique(mapper.expanded$uniprot)
  missingUniprots <- setdiff (uniprots, fasta.dt$uniprot)
  message (length(missingUniprots), " protein IDs not found in fasta file: ",
                  paste0(head(missingUniprots, 10), collapse = ","),
           ifelse(length(missingUniprots)>10, ",...", ""))
  
  if (downloadFromWeb){
    web.up.dt <- downloadUniprotSequences (missingUniprots)
    up.dt <- rbind (fasta.dt, web.up.dt)
  }
  
  
  # now combine protein sequences with peptide info
  # merge/modify with up.dt
  # first count matches in case the peptide is duplicated in the protein
  mapper.expanded[up.dt, peptideCopies := stringr::str_count(sequence, cleanPeptide),on = "uniprot"]
  if (any (mapper.expanded$peptideCopies == 0)){
    notMatched <- unique(mapper.expanded[peptideCopies == 0]$fragPipeProteinName)
    message (length(notMatched), " peptides could not be located within their protein. Is this the correct FASTA?")
    print(head(notMatched))
    message ("The unmatched peptides will be removed from futher analysis")
    mapper.expanded <- mapper.expanded[peptideCopies > 0]
  }
  if (any (mapper.expanded$peptideCopies > 1)){
    message ("Some peptides are duplicated within their protein.  In these cases only the first possible location will be labeled")
    print (mapper.expanded[peptideCopies > 1 , 
                           .(maxCopies = max(peptideCopies), countCopiedPeptides = length(unique(cleanPeptide))),
                           by = uniprot])
  }

  mapper.expanded[up.dt, peptideStart := stringr::str_locate(sequence, cleanPeptide)[,"start"],on = "uniprot"]
  mapper.expanded[, ptmPos := pepSitePTM + peptideStart -1 ]
  
  # now collapse back to single row per fragPipeProteinName
  mapper.collapsed <- mapper.expanded[,.(ptmSiteCombo = paste0(uniprot, "_", aaPTM, ptmPos, collapse = ";")), by = fragPipeProteinName]
  
  return (mapper.collapsed[proteinNames, , on =  "fragPipeProteinName"]$ptmSiteCombo)
}

