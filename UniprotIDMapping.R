

localDir <- getwd()


# uniprot data file ####
downloadHumanIDDatMap <- function(saveFile = file.path(localDir,"data/human.uniprot.idmap.dat.gz")){
  human.idmap.dat <- fread ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
                            header=FALSE)
  setnames(human.idmap.dat, c("uniprot", "idType", "id"))
  fwrite(human.idmap.dat, file=saveFile)
}

loadHumanIDDatMap <- function(reload=FALSE, path = NULL){
  if (is.null(path))
    path <- file.path(localDir,"data/human.uniprot.idmap.dat.gz")
  createTime <- file.info(path)$ctime
  if (is.na(createTime)){
    message ("Local ID mapping file doesn't exist at ", path, " reloading.")
    reload = TRUE
  }else {
    fileAge <- Sys.time() - createTime
    if (Sys.time() - file.info(path)$ctime > as.difftime(5, units="mins")){
      message ("Local file is ", as.integer(fileAge), attr(fileAge, "units"), " old. Consider forcing a reload")
    }
  }
  if (reload){
    downloadHumanIDDatMap (saveFile = path)
  }
  fread (path)
}


ensembl2gene <- function(){
  stop("not implemented yet, see comment below for a start")
}
# ids <- loadHumanIDDatMap()
# up2gene <- ids[idType == "Gene_Name",.(uniprot, gene = id)]
# up2ensembl <- ids[idType == "Ensembl", .(uniprot, ensembl = id)]
# 
# ensembl2gene <- merge (up2ensembl, up2gene, by = "uniprot")
# ensembl2gene <- ensembl2gene[,.(uniprots = paste0(uniprot, collapse = ";")), by = .(ensembl, gene)]


##########

# mouse uniprot data file ####
downloadMouseIDDatMap <- function(saveFile = file.path(localDir,"data/mouse.uniprot.idmap.dat.gz")){
  idmap.dat <- fread ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz",
                            header=FALSE)
  setnames(idmap.dat, c("uniprot", "idType", "id"))
  dir.create(dirname(saveFile), recursive=TRUE, showWarnings = FALSE)
  fwrite(idmap.dat, file=saveFile)
}

downloadRatIDDatMap<- function(saveFile = file.path(localDir,"data/mouse.uniprot.idmap.dat.gz")){
  
  idmap.dat <- fread ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/RAT_10116_idmapping.dat.gz",
                      header=FALSE)
  setnames(idmap.dat, c("uniprot", "idType", "id"))
  dir.create(dirname(saveFile), recursive=TRUE, showWarnings = FALSE)
  fwrite(idmap.dat, file=saveFile)
}

loadMouseIDDatMap <- function(reload=FALSE, path = NULL){
  if (is.null(path))
    path <- file.path(localDir,"data/mouse.uniprot.idmap.dat.gz")
  createTime <- file.info(path)$ctime
  if (is.na(createTime)){
    message ("Local ID mapping file doesn't exist at ", path, " reloading.")
    reload = TRUE
  }else {
    fileAge <- Sys.time() - createTime
    if (Sys.time() - file.info(path)$ctime > as.difftime(5, units="mins")){
      message ("Local file is ", as.integer(fileAge), attr(fileAge, "units"), " old. Consider forcing a reload")
    }
  }
  if (reload){
    downloadMouseIDDatMap (saveFile = path)
  }
  fread (path)
}
loadRatIDDatMap <- function(reload=FALSE, path = file.path(localDir,"data/rat.uniprot.idmap.dat.gz")){
  createTime <- file.info(path)$ctime
  if (is.na(createTime)){
    message ("Local ID mapping file doesn't exist at ", path, " reloading.")
    reload = TRUE
  }else {
    fileAge <- Sys.time() - createTime
    if (Sys.time() - file.info(path)$ctime > as.difftime(5, units="mins")){
      message ("Local file is ", as.integer(fileAge), attr(fileAge, "units"), " old. Consider forcing a reload")
    }
  }
  if (reload){
    downloadRatIDDatMap (saveFile = path)
  }
  fread (path)
}

translateUniprot2String <- function (uniprots, species = "MOUSE"){
  if (toupper(species) == "HUMAN"){
    idMapper <- loadHumanIDDatMap()[idType == "STRING",] #has columns uniprot,idType,id
  }else if (toupper(species) == "MOUSE"){
    idMapper <- loadMouseIDDatMap()[idType == "STRING",] #has columns uniprot,idType,id
  }
  setnames(idMapper, old=c("id"), new=c("geneName"))
  idMapper[match(uniprots, uniprot), geneName]
  
}


# string data from string ####

# ... arguments passed to fread (url)
readLocalVersionOfRemoteFile<- function(localPath, url, reload = FALSE, ...){
  createTime <- file.info(localPath)$ctime
  if (is.na(createTime)){
    message ("Local ID mapping file doesn't exist at ", localPath, " reloading.")
    reload = TRUE
  }else {
    fileAge <- Sys.time() - createTime
    if (Sys.time() - file.info(localPath)$ctime > as.difftime(5, units="mins")){
      message ("Local file is ", as.integer(fileAge), attr(fileAge, "units"), " old. Consider forcing a reload")
    }
  }
  if (reload){
    message ("Downloading ID mapping file to save locally")
    message ("remote: ", url)
    message ("local:  ", localPath)
    remoteDT <- fread (url, ...)
    dir.create(dirname(localPath), recursive=TRUE, showWarnings = FALSE)
    fwrite (remoteDT, localPath)
  }
  fread (localPath)
}


loadStringToUniprotFile <- function (species, dataDir = file.path(localDir, "data"), reload = FALSE){
  if(toupper(species) == "MOUSE"){
    url <- "https://string-db.org/mapping_files/uniprot/mouse.uniprot_2_string.2018.tsv.gz"
  } else if (toupper(species)  == "HUMAN"){
    url <- "https://string-db.org/mapping_files/uniprot/human.uniprot_2_string.2018.tsv.gz"
  }
  fileName <- basename(url)
  # try loading locally
  localFilePath <- file.path(dataDir, fileName)
  string2Uniprot.dt <- readLocalVersionOfRemoteFile (localFilePath, url, reload)
  
  setnames(string2Uniprot.dt, new = c("speciesID", "upAcc", "string", "pid", "length")) # those last 2 are guesses
  string2Uniprot.dt[, c("uniprot", "uniprotName") := tstrsplit(upAcc, split = "\\|")]
  return (string2Uniprot.dt)
    
}



# using string data
translateString2Uniprot <- function (stringIDs, species  = "MOUSE", fillMissing = TRUE, reload=FALSE){
  idMapper <- loadStringToUniprotFile (species, reload = reload)
  setkey(idMapper, string)
  matched <- idMapper[stringIDs,.(string, uniprot)]
  if(fillMissing == TRUE)
    matched[is.na(uniprot), uniprot := string]
  return (matched$uniprot)
}


loadString2EntrezFile <- function (species, dataDir = file.path(localDir, "data"), reload = FALSE){
  if(toupper(species) == "MOUSE"){
    url <- "https://string-db.org/mapping_files/entrez/mouse.entrez_2_string.2018.tsv.gz"
  } else if (toupper(species)  == "HUMAN"){
    url <- "https://string-db.org/mapping_files/entrez/human.entrez_2_string.2018.tsv.gz"
  }
  fileName <- basename(url)
  # try loading locally, and reload if necessary or requested
  localFilePath <- file.path(dataDir, fileName)
  string2Entrez.dt <- readLocalVersionOfRemoteFile (localFilePath, url, reload)
  
  setnames(string2Entrez.dt, new = c("speciesID",  "entrez", "string")) 
  return (string2Entrez.dt)
}


loadString2GeneFile <- function (species, dataDir = file.path(localDir, "data"), reload = FALSE){
  if(toupper(species) == "MOUSE"){
    url <- "https://string-db.org/mapping_files/STRING_display_names/mouse.name_2_string.tsv.gz"
  } else if (toupper(species)  == "HUMAN"){
    url <- "https://string-db.org/mapping_files/STRING_display_names/human.name_2_string.tsv.gz"
  }
  fileName <- basename(url)
  # try loading locally, and reload if necessary or requested
  localFilePath <- file.path(dataDir, fileName)
  string2Gene.dt <- readLocalVersionOfRemoteFile (localFilePath, url, reload)
  
  setnames(string2Gene.dt, new = c("speciesID", "gene", "string")) # those last 2 are guesses
  #string2Uniprot.dt[, c("uniprot", "uniprotName") := tstrsplit(upAcc, split = "\\|")]
  return (string2Gene.dt)
}



translateString2Gene <- function (stringIDs, species = "MOUSE", fillMissing = TRUE, reload = FALSE){
  idMapper <- loadString2GeneFile (species, reload = reload)
  setkey(idMapper, string)
  matched <- idMapper[stringIDs,.(string, gene)]
  if(fillMissing == TRUE)
    matched[is.na(gene), gene := string]
  return (matched$gene)
  
}

translateGene2String <- function (geneNames, species = "MOUSE", fillMissing = TRUE, reload = FALSE){
  idMapper <- loadString2GeneFile (species, reload = reload)
  setkey(idMapper, gene)
  matched <- idMapper[geneNames,.(string, gene)]
  if(fillMissing == TRUE)
    matched[is.na(string), string := gene]
  return (matched$string)
}


translateEntrez2String <- function (entrezNumbers, species = "MOUSE", fillMissing = TRUE, reload = FALSE){
  idMapper <- loadString2EntrezFile (species, reload = reload)
  setkey(idMapper, entrez)
  matched <- idMapper[entrezNumbers,.(string, entrez)]
  if(fillMissing == TRUE)
    matched[is.na(string), string := entrez]
  return (matched$string)
}


# . ####
translateString2Uniprot.usingUniprotData <- function (stringIDs, species = "MOUSE"){
  if (toupper(species) == "HUMAN"){
    idMapper <- loadHumanIDDatMap()[idType == "STRING",] #has columns uniprot,idType,id
  }else if (toupper(species) == "MOUSE"){
    idMapper <- loadMouseIDDatMap()[idType == "STRING",] #has columns uniprot,idType,id
  }
  setnames(idMapper, old=c("id"), new=c("string"))
  message ("Double check non-redundancy filter is required...")
  idMapper[,uniprot[1], by = string] # make 1 to 1
  idMapper[match(stringIDs, string), uniprot]
}




translateUniprot2GeneName.datFile <- function(uniprots, species="HUMAN", path = NULL){
  if (toupper(species) == "HUMAN"){
    idMapper <- loadHumanIDDatMap(path = path)[idType == "Gene_Name",] #has columns uniprot,idType,id
  }else if (toupper(species) == "MOUSE"){
    idMapper <- loadMouseIDDatMap(path = path)[idType == "Gene_Name",] #has columns uniprot,idType,id
  }else if (toupper(species) == "RAT"){
    idMapper <- loadRatIDDatMap(path = path)[idType == "Gene_Name",] #has columns uniprot,idType,id
  }
  setnames(idMapper, old=c("id"), new=c("geneName"))
  idMapper[match(uniprots, uniprot), geneName]
}

translateUniprot2EntrezGeneID <- function(uniprots, species="HUMAN"){
  if (toupper(species) == "HUMAN"){
    idMapper <- loadHumanIDDatMap()[idType == "GeneID",] #has columns uniprot,idType,id
  }else if (toupper(species) == "MOUSE"){
    idMapper <- loadMouseIDDatMap()[idType == "GeneID",] #has columns uniprot,idType,id
  }
  setnames(idMapper, old=c("id"), new=c("geneID"))
  idMapper[match(uniprots, uniprot), geneID]
}


translateEntrez2Uniprot <- function(entrez, species="HUMAN"){
  if (toupper(species) == "HUMAN"){
    idMapper <- loadHumanIDDatMap()[idType == "GeneID",] #has columns uniprot,idType,id
  }else if (toupper(species) == "MOUSE"){
    idMapper <- loadMouseIDDatMap()[idType == "GeneID",] #has columns uniprot,idType,id
  }
  setnames(idMapper, old=c("id"), new=c("geneID"))
  subset <- idMapper[geneID %in% entrez,]
  if (any(subset[,.N, by = geneID]$N > 1)){
    message ("Multiple uniprots per entrez, only returning a single one.")
  }
  subset[match(entrez, geneID), uniprot]
}






# using uniprot selected file resources ####

downloadHumanIDmap <- function(saveFile = file.path(localDir,"data/human.uniprot.idmap.gz")){
  human.idmap = fread ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping_selected.tab.gz")
  setnames(human.idmap, c("UniProtKB_AC", "UniProtKB_ID", "EntrezGene", "RefSeq", "GI", "PDB", "GO", "UniRef100", "UniRef90", "UniRef50", "UniParc", "PIR", "NCBI-taxon", "MIM", "UniGene", "PubMed", "EMBL", "EMBL-CDS", "Ensembl", "Ensembl_TRS", "Ensembl_PRO", "Additional_PubMed"))
  fwrite (human.idmap, file = saveFile)
}

  
loadHumanIDmap<- function(reload=FALSE, path = file.path(localDir,"data/human.uniprot.idmap.gz")){
  createTime <- file.info(path)$ctime
  if (is.na(createTime)){
    message ("Local ID mapping file doesn't exist at ", path, " reloading.")
    reload = TRUE
  }else {
    fileAge <- Sys.time() - createTime
    if (Sys.time() - file.info(path)$ctime > as.difftime(5, units="mins")){
      message ("Local file is ", as.integer(fileAge), attr(fileAge, "units"), " old. Consider forcing a reload")
    }
  }
  if (reload){
    downloadHumanIDmap (saveFile = path)
  }
  fread (path)
}

# using uniprot web services API ####

myQueryUniprot <- function (query = character(0L),
                            fields = c("accession", "id"), 
                            collapse = " OR ",
                            size= 25, ...) 
{
  #stopifnot(isCharacter(query), isCharacter(fields))
  if (!length(query)) 
    stop("<internal> 'qlist' must be populated with queries")
  resp <- GET(paste0(UNIPROT_REST_URL, "uniprotkb/search"), 
              query = list(query = paste(query, collapse = collapse), 
                           fields = paste(fields, collapse = ","),
                           format = "tsv",
                           size = size,
                           ...))
  fread(text = content(resp, encoding = "UTF-8"))
}



.query_uniprot_ws_for_sequence <- function(uniprots, fields = c("accession","sequence")){
  queries <- paste0("accession:", uniprots)
  res <- setDT(UniProt.ws::queryUniProt(queries, fields = fields))
  #print (res)
  res[uniprots,query := Entry , on = "Entry"]
  res[]
  
}

uniprotSequencesFromWeb <- function (uniprotIDs, chunkSize = 25, fields = c("accession","sequence")){
  uniqueUniprots <- unique(uniprotIDs)
  chunks <- split(uniqueUniprots, (1:length(uniqueUniprots))%/%chunkSize)
  
  infoMapList <- pbapply::pblapply(chunks, .query_uniprot_ws_for_sequence, fields = fields)

  infoMap <- rbindlist(infoMapList)
  return (infoMap)
  
}


.query_uniprot_ws_for_species <- function(uniprots){
  queries <- paste0("accession:", uniprots)
  res <- setDT(UniProt.ws::queryUniProt(queries, fields = c("accession", "id", "organism_name", "organism_id", "gene_primary", "gene_names", "protein_name")))
  #print (res)
  res[uniprots,query := Entry , on = "Entry"]
  res[]
}

# returns the table, not just species
uniprotInfoFromWeb <- function(uniprotIDs, chunkSize = 25){
  print (head(uniprotIDs))
  uniqueUniprots <- unique(uniprotIDs)
  
  chunks <- split(uniqueUniprots, (1:length(uniqueUniprots))%/%chunkSize)
  
  speciesMapList <- pbapply::pblapply(chunks, .query_uniprot_ws_for_species)
  # function(uniprotIDs)
  #   UniProt.ws::queryUniProt(uniprotIDs, fields = c("accession", "id", "organism_name", "organism_id")))
  
  speciesMap <- rbindlist(speciesMapList)
  return (speciesMap)
  
}

speciesFromUniprotID <- function (uniprotIDs, chunkSize = 25){
  #speciesMap <- UniProt.ws::queryUniProt(unique(uniprotIDs, fields = c("accession", "id", "organism_name", "organism_id")))
  #setDT[speciesMap]
  speciesMap <- uniprotInfoFromWS(uniprotIDs)
  return (speciesMap[uniprotIDs, Organism, on= "Entry"])
}



multiSpeciesFromMultiUniprots <- function(uniprots, sep = ";", simplify = TRUE){
  toSpecies <- data.table(uniprots = unique(uniprots))
  toSpecies <- toSpecies[,.(singleUniprot = unlist(strsplit(uniprots, sep))),by = uniprots]
  toSpecies[, singleSpecies :=  speciesFromUniprotID(singleUniprot)]
  if (simplify == TRUE){
    simplify = function(x)unique(sort(x))
  }else if (simplify == FALSE){
    simplify = identity # do nothing
  }else if (! "function" %in% class (simplify)){
    stop("unexpected simplify format")
  }
  
  toSpecies <- toSpecies[, .(species = paste(simplify(singleSpecies), collapse = sep)), by = uniprots]
  
  return (toSpecies[uniprots, species, by = "uniprots"])
}

#  using biomart ####
downloadHumanIDmap.biomart <- function(saveFile = file.path(localDir,"data/human.biomart.idmap.gz")){
  ensembl = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
  message ("Loading uniprot and hgnc identifiers from Ensembl Biomart.  This could take a few minutes.")
  idmap <- biomaRt::getBM(attributes = c("uniprotswissprot", "uniprot_gn_id",  "hgnc_symbol", "description"), mart=ensembl)
  fwrite(idmap, saveFile)

}

loadBioMartIDMap <- function(reload=FALSE, path = file.path(localDir,"data/human.biomart.idmap.gz")){
  createTime <- file.info(path)$ctime
  if (reload==TRUE | is.na(createTime)){
    message ("Local ID mapping file doesn't exist at ", path, " reloading.")
    reload = TRUE
  }else {
    fileAge <- Sys.time() - createTime
    if (Sys.time() - file.info(path)$ctime > as.difftime(10, units="days")){
      message ("Local file is ", as.integer(fileAge), units(fileAge), " old. Consider forcing a reload")
    }
  }
  if (reload){
    downloadHumanIDmap.biomart (saveFile = path)
  }
  fread (path)
}


translateUniprot2HGNC <- function(uniprots, dataPath = file.path(localDir,"data")){
  mapper <- loadBioMartIDMap(path=file.path(dataPath, "human.biomart.idmap.gz"))
  mapper[match(uniprots, uniprot_gn_id), hgnc_symbol]
}


# this one works a little different than the others
getDescription4Uniprots <- function (uniprots, dataPath = file.path(localDir,"data"), add_symbol = FALSE){
  mapper <- loadBioMartIDMap(path=file.path(dataPath, "human.biomart.idmap.gz"))
  if (add_symbol){
    mapper[,description := paste (hgnc_symbol, description, sep=" : ")]
  }
  merge (data.table(requested = uniprots),
         unique(mapper[,.(uniprot_gn_id, description)]),
         by.x="requested", by.y = "uniprot_gn_id", all.x=TRUE)
}





### Mouse gene to entrez from MGI ####

downloadRemoteTable<-function(url, saveFile){
  data <- fread (url)
  fwrite (data, saveFile)
  invisible(data)
}

loadFileWithDownloadOption <- function(url, path, reload=FALSE){
  createTime <- file.info(path)$ctime
  if (reload==TRUE | is.na(createTime)){
    message ("Local ID mapping file doesn't exist at ", path, " reloading.")
    reload = TRUE
  }else {
    fileAge <- Sys.time() - createTime
    if (Sys.time() - file.info(path)$ctime > as.difftime(10, units="days")){
      message ("Local file is ", as.integer(fileAge), units(fileAge), " old. Consider forcing a reload")
    }
  }
  if (reload){
    downloadRemoteTable(url, saveFile = path)
  }
  fread (path)
  
}


loadMouseGeneNameToEntrez.mgi <- function(path = file.path(localDir, "data", "mouse.mgi.gene2entrez.rpt.gz"),
                                          reload=FALSE){
  result <- loadFileWithDownloadOption ("http://www.informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt",
                                        path,
                                        reload)
  setnames(result, new = c("mgi.ac", "mgi.symbol", "status", "mgi.name", "cM_Position", "Chromosome",
                     "type", "secondary.ac", "entrez.geneID", "synonyms", "feature_types", "genome_start",
                     "genome_end", "strand", "biotypes"))
  result
}


translateGeneName2Entrez <- function (geneNames, species="MOUSE"){
  if (species=="MOUSE"){
    mapper <- loadMouseGeneNameToEntrez.mgi()
    return (mapper[match(geneNames, mgi.symbol), entrez.geneID])
  } else if (species == "HUMAN"){
    entrezs <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, geneNames, 'ENTREZID', 'SYMBOL', multiVals == "first") 
      stopifnot (all(names(entrezs) == geneNames, na.rm=TRUE))
    return (unlist(entrezs))
  }else{
    stop("unrecognized species", species)
  }
}


translateUniprot2GeneName <- function(uniprots, species = "HUMAN", useDatFile = FALSE, fillMissing = FALSE){
  translateUniprot2Something(uniprots, something = "SYMBOL", 
                             species = species, useDatFile = useDatFile, fillMissing = fillMissing)
}

translateUniprot2Something <- function (uniprots, something = "SYMBOL", species = "HUMAN", fillMissing = FALSE, useDatFile = FALSE) {
  
  if (useDatFile != FALSE){
    if (something != "SYMBOL"){
      stop("Currently only know how to use dat file for gene `SYMBOL`")
    }
    if ("character" %in% class(useDatFile))
      path = useDatFile
    else
      path = NULL
    if (fillMissing != TRUE)
      message("fillMissing is currently ignored when using datFile")
    return(translateUniprot2GeneName.datFile(uniprots, species, path = path))
  }
  
  uniprotsNoNA <- uniprots[!is.na(uniprots)]
  
  # do AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db) to look up allowed columns
  
  if (species == "HUMAN"){
    geneNames <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(uniprotsNoNA), something, 'UNIPROT', multiVals = "first")
  }else if (species == "MOUSE"){
    geneNames <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, unique(uniprotsNoNA), something, 'UNIPROT', multiVals = "first")
  }else if (species == "RAT"){
    geneNames <- AnnotationDbi::mapIds(org.Rn.eg.db::org.Rn.eg.db, unique(uniprotsNoNA), something, 'UNIPROT', multiVals = "first")
  } else {
    stop("unrecognized species", species)
  }
  
  # a change in AnnotationDbi::mapIds (updated Feb 2023) makes it return a list, even when multiVals = "first"
  # not due to a change in AnnodationDbi::mapIds.  
  # Some weird behavior with sapply (which I use here and also AnnotationDbi uses), not returning a vector when expected.
  # solved: it happens when there is a NULL in the list.  sapply(list("a", NULL), identity) does not return a vector
  # the NULL happens when any of uniprots above is NA.  Should be fixed, lets see...

  if ("list" %in% class(geneNames)){
    message ("debugMessage: list returned by AnnotationDbi::mapIds")
    geneNames <- unlist(sapply(geneNames, function(x)x[[1]]))
  }
  
  
  mapTable <- unique(data.table(uniprot = names(geneNames), gene = geneNames))
  setkey(mapTable, "uniprot")
  # this orders all and expands the missing cases
  mapTable <- mapTable[uniprots]
  if(fillMissing == TRUE){
    # where gene lookup failed, assign the original uniprot
    mapTable[is.na(gene), gene := uniprot]
  }
  return (mapTable$gene)
}

# do AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db) to look up allowed columns
translateSomething2SomethingElse <- function (somethings, originalType = "UNIPROT",  somethingElse = "SYMBOL", species = "HUMAN", fillMissing = FALSE ) {
  allowedColumns <- AnnotationDbi::columns(org.Hs.eg.db::org.Hs.eg.db)
  
  if (!all(c(originalType, somethingElse)%in% allowedColumns) ){
    message ("Requested unexpected types. Allowed types are: ")
    print (allowedColumns)
  }

  somethingsNoNA <- as.character(somethings[!is.na(somethings)])
  
  
  if (species == "HUMAN"){
    geneNames <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(somethingsNoNA), somethingElse, originalType, multiVals = "first")
  }else if (species == "MOUSE"){
    geneNames <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, unique(somethingsNoNA), somethingElse, originalType, multiVals = "first")
  }else if (species == "RAT"){
    geneNames <- AnnotationDbi::mapIds(org.Rn.eg.db::org.Rn.eg.db, unique(somethingsNoNA), somethingElse, originalType, multiVals = "first")
  } else {
    stop("unrecognized species", species)
  }
  
  # a change in AnnotationDbi::mapIds (updated Feb 2023) makes it return a list, even when multiVals = "first"
  # not due to a change in AnnodationDbi::mapIds.  
  # Some weird behavior with sapply (which I use here and also AnnotationDbi uses), not returning a vector when expected.
  # solved: it happens when there is a NULL in the list.  sapply(list("a", NULL), identity) does not return a vector
  # the NULL happens when any of uniprots above is NA.  Should be fixed, lets see...
  
  if ("list" %in% class(geneNames)){
    message ("debugMessage: list returned by AnnotationDbi::mapIds")
    geneNames <- unlist(sapply(geneNames, function(x)x[[1]]))
  }
  
  
  mapTable <- unique(data.table(something = names(geneNames), gene = geneNames))
  setkey(mapTable, something)
  # this orders all and expands the missing cases
  mapTable <- mapTable[as.character(somethings)]
  if(fillMissing == TRUE){
    # where gene lookup failed, assign the original uniprot
    mapTable[is.na(gene), gene := something]
  }
  return (mapTable$gene)
}



translateGeneName2Uniprot <- function(geneNames, species = "HUMAN", fillMissing = FALSE){

  if (species == "HUMAN"){
    dbName <- org.Hs.eg.db::org.Hs.eg.db
  }else if (species == "MOUSE"){
    dbName <- org.Mm.eg.db::org.Mm.eg.db
  }else if (species == "RAT"){
    dbName <- org.Rn.eg.db::org.Rn.eg.db
  } else {
    stop("unrecognized species", species)
  }
  uniprots <- AnnotationDbi::mapIds(dbName, unique(geneNames), 'UNIPROT', 'SYMBOL', multiVals == "first")
  mapTable <- unique(data.table(gene = names(uniprots), uniprot = uniprots))
  setkey(mapTable, "gene")
  # this orders all and expands the missing cases
  mapTable <- mapTable[geneNames]
  if (fillMissing == TRUE){
    # where gene lookup failed, assign the original uniprot
    mapTable[is.na(uniprot), uniprot := gene]
  }
  return (mapTable$uniprot)
}






multiUniprots2multiGenes <- function (uniprots, sep = ";", species = "HUMAN", simplify = FALSE, useDatFile = FALSE, allowDups = FALSE){
  toGenes <- data.table(uniprots = uniprots)
  toGenes <- toGenes[,.(singleUniprot = unlist(strsplit(uniprots, sep))),by = uniprots]
  toGenes[,singleGene := translateUniprot2GeneName(singleUniprot, species = species, useDatFile = useDatFile)]
  toGenes[is.na(singleGene), singleGene := singleUniprot]
  if (simplify == TRUE){
    simplify = function(x)unique(sort(x))
  }else if (simplify == FALSE){
    simplify = identity # do nothing
  }else if (! "function" %in% class (simplify)){
    stop("unexpected simplify format")
  }
  toGenes <- toGenes[, .(genes = paste(simplify (singleGene), collapse=sep)), by = uniprots]
  if(!allowDups){
    duplicatedGeneNames <- unique(toGenes[duplicated(genes)]$genes)
    toGenes[genes %in% duplicatedGeneNames, genes := paste0(genes, ".", uniprots)]
  }
  
  setkey(toGenes, uniprots)
  
  return(toGenes[uniprots, genes])
}


multiUniprotSites2multiGeneSites <- function (uniprotSites, sep = ";", siteSep = "_", species = "HUMAN", useDatFile = FALSE, uniqueOnly = FALSE){
  mapper <- data.table(uniprotSites = unique(as.character(uniprotSites)))
  #expand to singleSite, 1 per row
  mapper <- mapper[,.(singleSite = unlist(strsplit(uniprotSites, split = sep))), by = uniprotSites]
  
  #deal with protein IDs that contain siteSep
  idParts <-  strsplit(mapper$singleSite, split = siteSep)
  protID <- sapply (idParts, function(x)paste(x[1:(length(x)-1)], collapse = siteSep))  
  site <- sapply (idParts, function(x)x[length(x)])  
  
  mapper[,c("uniprot", "site") := .(protID, site)]
  mapper[,gene := translateUniprot2GeneName(uniprot, species, useDatFile = useDatFile)]
  if (any(is.na(mapper$gene))){
    message (length(mapper[is.na(gene), unique(uniprot)]), " uniprots could not be mapped to genes (using their uniprot ID in gene column)")
    mapper[is.na(gene), gene := uniprot]
  }
  
  mapper[, singleGeneSite := paste(gene,site, sep = siteSep)]
  
  if (uniqueOnly) f <- unique
  else f <- identity
  
  # collapse back to the combined uniprotSites
  mapper <- mapper[,.(geneSite  = paste0(f(singleGeneSite), collapse = sep)), by = uniprotSites]
  
  toGenes <- data.table (uniprotSites)
  
  return(mapper[toGenes,,on="uniprotSites"]$geneSite)
}



# gene aliases


# useful to make sure the same gene symbol is used when combining different datasets

geneAlias2officialGeneSymbol <- function(geneAliases, species = "HUMAN"){
# if (FALSE == require (limma)  )
#   return (geneAliases)
  if (species == "HUMAN")
    speciesCode <- "Hs"
  else if (species == "MOUSE")
    speciesCode <- "Mm"
  else
    speciesCode <- species
  
  unique.genes <- unique(geneAliases)
  aliasTable <- data.table (alias = unique.genes, symbol = limma::alias2SymbolTable(unique.genes, species = speciesCode) )
  message ("Of ", length(unique.genes), " unique genes, ",
           aliasTable[alias != symbol & !is.na(symbol), length(alias)],
           " will be replaced with official symbols, and ",
           aliasTable[is.na(symbol), length(alias)], " were not found in alias table") 
  
  print (aliasTable[is.na(symbol), unique(alias)])
  
  aliasTable[is.na(symbol), symbol := alias]
  setkey(aliasTable, alias)
  return (aliasTable[geneAliases,]$symbol)
}


library(data.table)
mouse2HumanUniprotsFullTable <- function(run = FALSE){
  if(run != TRUE)
    stop("This function downloads about 100MB of data to build the mouse2Human uniprots table.
         You've been warned--don't run this in a loop.  Try again and set run = TRUE as its only argument")

  # get the human to mouse mapping from ensembl bioMart, 
  # these are ensembl gene IDs like this ENSG00000000003 or ENSMUSG00000067377
  human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  human2mouse <- setDT(biomaRt::getBM(attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"), mart = human))

  
  #load uniprot mappings to map from ensembl genes to uniprot IDs
  human.idmap.dat <- fread ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
                            header=FALSE)
  setnames(human.idmap.dat, c("uniprot", "idType", "id"))
  mouse.idmap.dat <- fread ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping.dat.gz",
                      header=FALSE)
  setnames(mouse.idmap.dat, c("uniprot", "idType", "id"))
  # take relevant subset of idmap.dat
  ensemble2uniprot.mouse <- mouse.idmap.dat[idType == "Ensembl", .(uniprot.mouse = uniprot, ensembl = id)]
  ensemble2uniprot.human <- human.idmap.dat[idType == "Ensembl", .(uniprot.human = uniprot, ensembl = id)]

  # do the merges
  merge.mouse <- merge (human2mouse, ensemble2uniprot.mouse, by.x = "mmusculus_homolog_ensembl_gene", by.y = "ensembl" )
  merge.mouse.human <- merge (merge.mouse, ensemble2uniprot.human, by.x = "ensembl_gene_id", by.y = "ensembl", allow.cartesian=TRUE)
  
  return (merge.mouse.human)
}


