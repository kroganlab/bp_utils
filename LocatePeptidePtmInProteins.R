library (data.table)
# requires stringr, RCurl, seqinr

# next two files originally  in  FragPipe2MSstats.R see
# https://kroganucsf.slack.com/archives/CQ50RCEF7/p1629875481043800 
# (Slack accessible to Krogan Lab members only)

#' A wrapper written by Mehdi. Modifies ph (in place) and writes out contents to an output file
#' @param ph data.table that holds the processed contents of input_file, see FragPipe2MSstats.usage.R for processing
#' @input_file character string denoting the path of input--will be modified with _sitsmapped.csv for output file
#' @fasta_file character path to the fastafile holding sequences expected in ph
#' @protein_peptide_sep character string separator between proteinID and peptideID in ProteinKey column of ph
mehdi_map_sites = function(ph,input_file,fasta_file,protein_peptide_sep = "__"){
  
  ph[, ProteinSite := fragPipePTMFormat2SiteFormat(ProteinKey, fastaFile =  fasta_file, downloadFromWeb = FALSE, proteinPepSep = protein_peptide_sep)]
  #ph[,Protein := ProteinName]
  #ph[,ProteinName := newProteinName]
  #ph[,newProteinName := NULL]
  
  # remove those that didn't map...or that are don't have PH
  #noMissingPH <- ph[!is.na(ProteinName)]
  
  # Convert mass shifts to words
  ph$PeptideSequence = convertMassModificationFormat(ph$PeptideSequence)
  
  fwrite(ph, paste(gsub(".csv","",input_file),"_sitesmapped.csv",sep=""))
  
}

#' @param input_file character string
#' @param keys_file character string 
mehdi_fragpipe_addkeys = function(input_file, keys_file){

  # Load data
  keys = fread(keys_file)
  ph = fread (input_file)
  
  # Make Run column the Condition column
  ph$Run = ph$Condition
  
  # Replace the Condition column
  ph$Condition = keys$Condition[match(ph$Run,keys$Run)]
  
  # Replace BioReplicate column
  ph$BioReplicate = keys$BioReplicate[match(ph$Run,keys$Run)]
  
  # Remove any NA intensity values
  ph = ph[!is.na(ph$Intensity),]
  
  # Save out file
  fwrite(ph,paste(gsub('.csv','',input_file),"_wKeys.csv",sep=''),sep=',')
  
}
# /end new FragPipe2MSstats.R functions



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
  # expect uniprot format 'sp|Q13245|protein name'
  dt[grepl ("\\|", header), c("db", "uniprot", "uniprotName") := tstrsplit(header, split = "\\|", keep = 1:3)]
  # if no pipes then we use the whole header as the uniprot identifier
  dt[!grepl("\\|", header),  uniprot := header] 
  return (dt[])
}





#' @param proteinNames character vector, each in format like `Q9Y6Y0_SLS[79.9663]FEMQQDELIEKPMSPMQYAR`
#' @param ptmType  one of "PH" or in the future "UB", but for now only "PH"
#' @param fastaFile path to the fasta file, uniprot format, that should match the first field in protein names
#' @param downloadFromWeb logical TRUE = attempt to download uniprots not found in fasta file from web.

fragPipePTMFormat2SiteFormat <- function(proteinNames, ptmType = "PH", fastaFile = NULL, downloadFromWeb = TRUE, proteinPepSep = "_"){
  mapper <- data.table(fragPipeProteinName = unique(proteinNames))
  
  mapper[, c("uniprot", "modPeptide") := tstrsplit(fragPipeProteinName, proteinPepSep)]
  
  mapper[, modPeptideLC.PTM := convertMassModificationFormat(modPeptide, keepOnly = ptmType)]
  
  formats <- c(PH = "(ph)") # depends on function convertMassModificationFormat
  ptmFormat <- formats[ptmType]
  
  # first locate ptm in peptide. Mulitple matches per peptide
  mapper[, ptmPositions := stringr::str_locate_all(modPeptideLC.PTM, stringr::fixed (ptmFormat))]
  
  # we have a list of matrices. one matrix per row. expand as many rows as there are sites in the row
  mapper.expanded <- mapper[lapply(ptmPositions, nrow) > 0,
                            .(pepSitePTM = ptmPositions[[1]][,"start"] -
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
  }else{
    up.dt <- fasta.dt
  }
  
  
  # now combine protein sequences with peptide info
  # merge/modify with up.dt
  # first count matches in case the peptide is duplicated in the protein
  mapper.expanded[up.dt, peptideCopies := stringr::str_count(sequence, cleanPeptide),on = "uniprot"]
  if (any (mapper.expanded$peptideCopies == 0)){
    notMatched <- unique(mapper.expanded[peptideCopies == 0]$fragPipeProteinName)
    message (length(notMatched), " peptides could not be located within their protein. Is this the correct FASTA?","\n\t",
             "This message can be misleading if you have duplicate uniprot IDs in your FASTA, such as forward and reverse.")
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

