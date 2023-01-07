library (data.table)
# requires stringr, RCurl, seqinr

# next two functions originally  in  FragPipe2MSstats.R see
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


methylationPatterns <- c(single = "K[14.0156]",
                         double = "K[28.0313]",
                         triple = "K[42.0470]")


convertLysMethylationMassFormat <- function(specModSequence, mods=c("PH",  "CAM", "MOX", "NAC, KMET1, KMET2, KMET3"), keepOnly = NULL, removeAll = FALSE){
  newMassFormats = c(KMET1 = "(K)\\[14.0156\\]",
                     KMET2 = "(K)\\[28.0313\\]",
                     KMET3 = "(K)\\[42.0470\\]")
  newOutFormats = c(KMET1 = "\\1\\(kmeti\\)",
                    KMET2 = "\\1\\(kmetii\\)",
                    KMET3 = "\\1\\(kmetiii\\)")
  mods <- c(mods, names(newMassFormats))
  
  convertMassModificationFormat(specModSequence, mods, keepOnly, removeAll, newMassFormats, newOutFormats)
}
  

convertMassModificationFormat <- function(specModSequence, mods=c("PH",  "CAM", "MOX", "NAC"), keepOnly = NULL, removeAll = FALSE,
                                          newMassFormats = character(0), newOutFormats = character(0) ){
  result <- specModSequence
  massFormats <- c(MOX = "(.)\\[15.9949\\]",
                   CAM = "(.)\\[57.0215\\]",
                   PH = "(.)\\[79.9663\\]",
                   NAC = "n\\[42.0106\\](.)",
                   ANY = "n?\\[[0-9.]+\\]")  # used for matching all or all but keepOnly
  massFormats <- c(massFormats, newMassFormats)
  
  
  artmsFormats <- list (PH='\\1\\(ph\\)',
                        UB='\\1\\(gl\\)',
                        CAM = '\\1\\(cam\\)',
                        MOX = '\\1\\(ox\\)',
                        NAC = '\\1\\(ac\\)',
                        ANY = '')  # for removal purposes
  artmsFormats <- c(artmsFormats, newOutFormats)
  
  stopifnot (all(names(massFormats) %in% names(artmsFormats))) 
  
  
  if(!is.null(keepOnly) | removeAll == TRUE){
    mods <- c(keepOnly, "ANY")
  }
  
  # mass of UB is not yet known by me...
  # do the replacements one by one
  for (mod in mods){
    if (mod %in% names(massFormats)){
      replace <- artmsFormats[[mod]]
      if (!is.null(keepOnly) && !mod %in% keepOnly)
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

sitifyProteins_SpectronautFile <- function (specFile.dt, site = "PH", fastaFile = "~/UCSF/kroganlab/BenPolacco/data/human_all_proteins_canonical_uniprot-proteome_UP000005640.fasta.gz"){
  stopifnot (site == "PH")  # nothing else supported yet, I think. Def not tested at least
  shortPTM = tolower(site)
  
  
  mapper <- unique(specFile.dt[,.(ProteinName, specFormat = PeptideSequence)])
  mapper[, parenLowerCaseFormat := convertSpectronautModificationFormat (specFormat)]
  # spectronaut will put _ before and after seqeunce.  I don't want those
  mapper[, parenLowerCaseFormat := gsub("^_|_$", "", parenLowerCaseFormat)]
  mapper[, uniprotSite := parenLowerCaseToProteinSiteNames(ProteinName, 
                                                           parenLowerCaseFormat, 
                                                           shortPTM = shortPTM, 
                                                           fastaFile = fastaFile)]
  
  specFile.dt[mapper,
              c("ProteinName",   "parenLowerCaseFormat", "oldProteinName") :=
              .(i.uniprotSite,  i.parenLowerCaseFormat,     ProteinName),
              on = c("ProteinName", PeptideSequence = "specFormat")]
  
  invisible(specFile.dt[])

}







convertSpectronautModificationFormat <- function(specModSequence, mods=c("PH", "UB", "CAM", "MOX", "NAC"), keepOnly = NULL, removeAll = FALSE){
  result <- specModSequence
  # these are overly complicated because they match both S[Phospho (STY)] and S(Phospho (STY))
  specFormats <- list (PH='([STY])[[(]Phospho \\(STY\\)[])]',
                       UB='(K)[[(]GlyGly \\(K\\)[])]',
                       CAM = '([C])[[(]Carbamidomethyl \\(C\\)[])]',
                       MOX = '([M])[[(]Oxidation \\(M\\)[])]',
                       NAC =  '([A-Z_])[[(]Acetyl \\(Protein N-term\\)[])]',
                       ANY =          '[[(][^][)(]*\\([^][)(]*\\)[])]' )  # suggest regex101.com to parse this visually, paste and then change \\

  artmsFormats <- list (PH='\\1\\(ph\\)',
                        UB='\\1\\(gl\\)',
                        CAM = '\\1\\(cam\\)',
                        MOX = '\\1\\(ox\\)',
                        NAC = '\\1\\(ac\\)',
                        ANY = '')  # for removal purposes

  stopifnot(names(specFormats)==names(artmsFormats))
  
  if(!is.null(keepOnly) | removeAll == TRUE){
    mods <- c(keepOnly, "ANY")
  }
  
    
  for (mod in mods){
    if (mod %in% names(specFormats)){
      replace <- artmsFormats[[mod]]
      if (!is.null(keepOnly) && keepOnly != mod)
        replace <- ""
      result <- gsub(specFormats[[mod]], replace, result)
    }else (message("I don't yet know how to deal with mod: ", mod))
  }
  
  missedMods <- grep ("\\[.*\\]", result, value=TRUE)
  if (length(missedMods) > 1){
    message (length(missedMods), " peptides with unexpected modifications: ", paste0(head(missedMods), collapse = ", "))
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
  #check for empty strings:
  fastaTxt <- fastaTxt[fastaTxt != ""]
  if (length(fastaTxt) > 0){
    fasta.dt <- as.data.table(transpose(stringr::str_split(fastaTxt, n =2, pattern = "\n")))    
    setnames(fasta.dt, new = c("header", "sequence"))
    fasta.dt <- fasta.dt[grepl ("^>", header),]
    fasta.dt[, header := gsub("^>", "", header)]
    fasta.dt[, sequence := gsub ("\n", "", sequence)]
    fasta.dt[, c("db", "uniprot", "uniprotName") := tstrsplit (header, split = "[ |]", keep = 1:3)]
  }else{
    fasta.dt <- data.table()
  }
  
  
  numSuccess <- nrow(fasta.dt)
  if (numSuccess < length(uniprots) )
    message ("Failed to download ", length(uniprots) - numSuccess, " uniprots: ", paste0(setdiff(uniprots, fasta.dt$uniprot), collapse = "; ") )
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


loadUniprots <- function(ids, fastaFile = NULL, downloadFromWeb = TRUE){
  if (any(grepl(";", ids))){
    message ("Semicolons detected in uniprot IDs. These are likely multi-protein IDs and should be dealt with prior to looking them up in a file or online. These will be removed.")
    remove <- grep(";", ids, value = TRUE)
    ids <- setdiff(ids, remove)
    message ("\t", paste0(head(remove), collapse = "  "), ifelse(length(remove) > 6, "...", ""))
  }
  fasta.dt <- fastaFileToTable (fastaFile)
  missingUniprots <- setdiff (ids, fasta.dt$uniprot)
  if (length(missingUniprots) > 0)
    message (length(missingUniprots), " protein IDs not found in fasta file: ",
             paste0(head(missingUniprots, 10), collapse = ","),
             ifelse(length(missingUniprots)>10, ",...", ""))
  
  if (downloadFromWeb & length(missingUniprots) > 0){
    web.up.dt <- downloadUniprotSequences (missingUniprots)
    up.dt <- rbind (fasta.dt, web.up.dt)
  }else{
    up.dt <- fasta.dt
  }
  return (up.dt)
}

#' looks up S(ph) and similar (in the future) modifiaitons in a peptide and mathcing uniprots and returns  uniprot_AApos such as  A0AVK6_S417 and A0AVK6_S413;A0AVK6_S417 for doubly PH-ated
#' @param uniprots a vector of uniprot IDs
#' @param parenLowerCaseFormats a vector of peptides like "_INSAPSS[Phospho (STY)]PIK_"; should match `uniprots` (as two columns of a table)
#' @param shortPTM  "ph" or similar.  For now only "ph" is supported
#' @param fastaFile string path to the fastaFile from which to lookup the sequences that match `uniprots`
#' @param downloadFromWeb boolean
#' @param proteinSiteSep  The string to put between the uniprot and site ID
parenLowerCaseToProteinSiteNames <- function(uniprots, parenLowerCaseFormats, shortPTM = "ph", fastaFile = NULL, downloadFromWeb = TRUE, proteinSiteSep = "_"){
  stopifnot (shortPTM == "ph")  # Only ph is currently supported, I think.
  mapper <- data.table(uniprot = uniprots, parenLowerCaseFormat = parenLowerCaseFormats)
  
  # remove all but 1 type of PTM
  allButMyPTM.regex  <- sprintf("\\((?!%s)[a-z]*\\)", shortPTM) 
  mapper[, only1PTM := gsub ( allButMyPTM.regex, "",parenLowerCaseFormat, perl = TRUE)]

  
  # then locate ptm in peptide. Mulitple matches per peptide
  ptmString <- sprintf("(%s)", shortPTM)
  mapper[, ptmPositions := stringr::str_locate_all(only1PTM, stringr::fixed (ptmString))]
  
  # we have a list of matrices. one matrix per row. expand as many rows as there are sites in the row
  mapper.expanded <- mapper[lapply(ptmPositions, nrow) > 0,
                            .(pepSitePTM = ptmPositions[[1]][,"start"] -
                                # subtract 1 because the (ph) labels the preceding AA, and subtract an additional 4 per each preceding mod
                                seq(from = 1, by =  stringr::str_length(ptmString), length.out =nrow(ptmPositions[[1]]))  ),
                            by = .(uniprot, parenLowerCaseFormat, only1PTM)]
  
  mapper.expanded[, cleanPeptide := gsub(ptmString, "", only1PTM, fixed = TRUE)]
  mapper.expanded[, aaPTM := substr(cleanPeptide, pepSitePTM, pepSitePTM)]
  
  up.dt <- loadUniprots(unique(mapper.expanded$uniprot), fastaFile, downloadFromWeb)

  # now combine protein sequences with peptide info
  # merge/modify with up.dt
  # first count matches in case the peptide is duplicated in the protein
  mapper.expanded[up.dt, peptideCopies := stringr::str_count(sequence, cleanPeptide),on = "uniprot"]
  mapper.expanded[is.na(peptideCopies), peptideCopies := 0]
  
  if (any (mapper.expanded$peptideCopies == 0)){
    notMatched <- unique(mapper.expanded[peptideCopies == 0, .(up_pep = paste0(uniprot, "_", parenLowerCaseFormat))]$up_pep)
    message (length(notMatched), " peptides could not be located within their protein. Likely a few sequences not found in the FASTA. Is this the correct FASTA?","\n\t",
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
  
  mapper.collapsed <- mapper.expanded[,.(ptmSiteCombo = paste0(uniprot, proteinSiteSep, aaPTM, ptmPos, collapse = ";")), by = .(uniprot, parenLowerCaseFormat)]
  
  # order the table to match the inputs and return the ptmSiteCombo column
  return (mapper.collapsed[data.table(uniprot = uniprots, parenLowerCaseFormat = parenLowerCaseFormats),
                           ,
                           on = c("uniprot", "parenLowerCaseFormat")]$ptmSiteCombo)
  
}



#' @param proteinNames character vector, uniprots or fasta header like `sp|Q04726|TLE3_HUMAN`.
#'    Also possibly in format like `Q9Y6Y0_SLS[79.9663]FEMQQDELIEKPMSPMQYAR`
#'    this format may have been an artifact of a processed file I received
#' @param modPeptideSequence character vector, like SLS[79.9663]FEMQQDELIEKPMSPMQYAR.  
#'    If NULL, then will look in proteinNames for protein_peptide format
#' @param ptmType  one of "PH" or in the future "UB"...also KMET1, KMET2, KMET3
#' @param fastaFile path to the fasta file, uniprot format, that should match the first field in protein names
#' @param downloadFromWeb logical TRUE = attempt to download uniprots not found in fasta file from web.

fragPipePTMFormat2SiteFormat <- function(proteinNames, 
                                         modPeptideSequences = NULL, 
                                         ptmType = "PH", fastaFile = NULL, downloadFromWeb = TRUE, proteinPepSep = "_"){
  mapper <- unique(data.table (fragPipeProteinName = proteinNames, modPeptide = modPeptideSequences))
  
  # handle various proteinNames formats 
  if (is.null(modPeptideSequences)){
    mapper[, c("uniprot", "modPeptide") := tstrsplit(fragPipeProteinName, proteinPepSep)]
  } else{
    mapper[, uniprot := NA_character_]
    mapper[grep ("\\|", fragPipeProteinName), uniprot := tstrsplit(fragPipeProteinName, "\\|")[[2]]]
    mapper[is.na(uniprot), uniprot := fragPipeProteinName]
  }
  
  #  ?do anything about contam tags taht might be in `contam_sp|P00761|TRYP_PIG`  ?

  if (ptmType %in% c("KMET1", "KMET2", "KMET3")){
    mapper[, modPeptideLC.PTM := convertLysMethylationMassFormat(modPeptide, keepOnly = ptmType)]
  } else{
    mapper[, modPeptideLC.PTM := convertMassModificationFormat(modPeptide, keepOnly = ptmType)]
  }
  
  
  formats <- c(PH = "(ph)",
               UB = "(gl)",
               KMET1 = "(kmeti)",
               KMET2 = "(kmetii)",
               KMET3 = "(kmetiii)") # depends on output of function convertMassModificationFormat
  
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
  
  if (downloadFromWeb & length(missingUniprots) > 0){
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



#' @param msInput, a data.table in format suitable for MSstats input. 
#'                 expects columns ProteinName (e.g. a uniprot id); PeptideSequence (with modifications)
#' @param proteinSeqs, a data.table with full length seqeunce information for looking up the peptides in msInput
#'                expects columns `ProteinName` that matches format of `ProteinName` in `msInput`, and `sequence`
#'                The function `loadUniprots` is handy for this.
#'                Alternatively, this can be a path to a fasta file, in which case it will load and use uniprot ID as
#'                the `ProteinName``
#' @param ptms, a named character vector of exact PTM strings that are inserted between capital letter amino acids in teh modified peptide string
#'              Also, beware partial matches!  If one PTM format is a subset of another, this will fail.              
labelProteinsMultiPTMs <- function(msInput,
                                   proteinSeqs,
                                   ptms = c(ph = "[79.9663]",
                                            r = "[100.0]",
                                            rp = "[175.0]",
                                            rpp = "[250.0]")){
  
 
  ###################################################################
  # Warn and check about modification format: no capital letters allowed 
  ####################################################################
  for (ptm in ptms){
    if (sum(grepl(ptm, ptms, fixed = TRUE)) > 1){
      message (ptm, " is a substring of at least one other PTM")
      stop("PTMs can not be a substring of another PTM")
    }
  }
  
  peptides <- unique(msInput$PeptideSequence)
  detectedModifications <- strsplit(peptides, "[A-Z]")  |> # all things between capital letters
    unlist() |> unique() |>  # remove redundants
    setdiff ("")  # remove the empty string
  
  cat ("Notice: Any capital letter in a peptide sequence is interpreted as a peptide AA. Capital letters in the modification string, even within brackets, will cause this to fail. Detected modifications: \n")
  cat ("\t",  sprintf ('"%s" ', detectedModifications), "\n")
  cat ("If any of these seem wrong, it is likely because a capital letter was used. Fix and start over.", "\n")
  
  
  ####################################################
  # Locate each PTM of interest in its modified peptide
  # each PTM per peptide gets its own row.
  
  # helper function
  .locatePTMInPeptide <- function (format, peptides){
    dt <- data.table(modPeptide = peptides)
    dt[, ptmPositions := stringr::str_locate_all(modPeptide,  stringr::fixed  (format) )]
    expanded <- dt[lapply(ptmPositions, nrow) > 0,
                   .(modPepPos = ptmPositions[[1]][,"start"] - 1),  # subtract 1 because the (ph) labels the preceding AA,
                   by = .(modPeptide)]
    return (expanded)
  }
  
  modPepPos.ls <- lapply(ptms, .locatePTMInPeptide, peptides)
  modPepPos <- rbindlist(modPepPos.ls, idcol = "ptm")
  modPepPos[, modAA := substr(modPeptide, modPepPos, modPepPos)] # this will fail if somehow there are two modifications one after another.  I can't imagine that possibility for now...
  
  
  
  #######################################
  # locate peptides in proteins
  #########################################
  # do we have a data.table of sequences or a path to a fasta file:
  if (! "data.frame" %in% class(proteinSeqs)){ 
    cat  ("Loading seqeunces from fasta at ", proteinSeqs, "\n")
    proteinSeqs <- loadUniprots(unique(msInput$ProteinName), 
                                fastaFile = proteinSeqs,
                                downloadFromWeb = FALSE)
    proteinSeqs[, ProteinName := uniprot] # to match expectations below
    
  }

  # clean modified peptide for lookups
  modPepPos[, cleanPeptide := gsub ("[^A-Z]", "", modPeptide)]
  
  # map ProteinName back into table:
  modPepPos[unique(msInput[, .(modPeptide = PeptideSequence, ProteinName)]),
            ProteinName := i.ProteinName,
            on = "modPeptide"]
  
  # locate cleanSeq in parent sequence, but first a couple checks based on how many times each peptide occurs in the protein
  modPepPos[proteinSeqs, peptideCopies := stringr::str_count(sequence, cleanPeptide),on = "ProteinName"]
  if (any (modPepPos$peptideCopies == 0)){
    notMatched <- unique(modPepPos[peptideCopies == 0]$ProteinName)
    message (length(notMatched), " peptides could not be located within their protein. Is this the correct FASTA?","\n\t")
    print(head(notMatched))
    message ("The unmatched peptides will be removed from futher analysis")
    modPepPos <- modPepPos[peptideCopies > 0]
  }
  if (any (modPepPos$peptideCopies > 1)){
    message ("Some peptides are duplicated within their protein.  In these cases only the first possible location will be labeled")
    print (modPepPos[peptideCopies > 1 , 
                     .(maxCopyNumber = max(peptideCopies), countPeptidesWithDups = length(unique(cleanPeptide))),
                     by = ProteinName])
  }
  
  modPepPos[proteinSeqs, peptideStart := stringr::str_locate(sequence, cleanPeptide)[,"start"],on = "ProteinName" ]
  
  
  # ##############################################
  # translate modPepPos to cleanPepPos, by counting capital letters in modPeptide. 
  # the AA counts: 
  aaCounts <- lapply (strsplit(modPepPos$modPeptide, ""),  # split by each character
                      function(chars)cumsum(chars %in% LETTERS)) # use cumsum to count TRUEs
  
  
  # use counter variable to 'loop' over each row and map modPepPos to the aaCounts
  modPepPos[, cleanPepPos := sapply (1:length(aaCounts), function(i) aaCounts[[i]][modPepPos[i]] ) ]
  
  # sanity check that we get same AA in modPeptide as cleanPeptide
  modPepPos[, testModAA := substr(cleanPeptide, cleanPepPos, cleanPepPos)]
  stopifnot(modPepPos[, all(testModAA == modAA)])
  modPepPos[, testModAA := NULL]
  
  ################################
  # finally a protein-level position
  modPepPos[, pos := peptideStart + cleanPepPos - 1 ]
  
  
  ############################
  # package into names
  modPepPos[, siteProteinName := paste0(ProteinName, "_", modAA, pos, ptm)]
  # wrap up to per-peptide names
  setorder(modPepPos, ProteinName, pos)
  newProteinNames <- modPepPos[, .(ProteinName = paste0(siteProteinName, collapse = ";")),  by = .(modPeptide)]
  
  
  ###########################3
  # modify msInput 
  cat ("Modifying MSstats input file with new PTM-specific 'ProteinName's. Original ProteinName in column og.ProteinName\n")
  msInput[, og.ProteinName := ProteinName]
  msInput[newProteinNames, ProteinName := i.ProteinName, on = c(PeptideSequence = "modPeptide")]
  return (msInput)
}
