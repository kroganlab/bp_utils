

localDir <- getwd()


# uniprot data file ####
downloadHumanIDDatMap <- function(saveFile = file.path(localDir,"data/human.uniprot.idmap.dat.gz")){
  human.idmap.dat <- fread ("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz",
                            header=FALSE)
  setnames(human.idmap.dat, c("uniprot", "idType", "id"))
  fwrite(human.idmap.dat, file=saveFile)
}

loadHumanIDDatMap <- function(reload=FALSE, path = file.path(localDir,"data/human.uniprot.idmap.dat.gz")){
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

loadMouseIDDatMap <- function(reload=FALSE, path = file.path(localDir,"data/mouse.uniprot.idmap.dat.gz")){
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



translateString2Uniprot <- function (stringIDs, species = "MOUSE"){
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




translateUniprot2GeneName.datFile <- function(uniprots, species="HUMAN"){
  if (toupper(species) == "HUMAN"){
    idMapper <- loadHumanIDDatMap()[idType == "Gene_Name",] #has columns uniprot,idType,id
  }else if (toupper(species) == "MOUSE"){
    idMapper <- loadMouseIDDatMap()[idType == "Gene_Name",] #has columns uniprot,idType,id
  }else if (toupper(species) == "RAT"){
    idMapper <- loadRatIDDatMap()[idType == "Gene_Name",] #has columns uniprot,idType,id
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

translateUniprot2GeneName <- function(uniprots, species = "HUMAN", useDatFile = FALSE){
  if (useDatFile){
    return(translateUniprot2GeneName.datFile(uniprots, species))
  }
  if (species == "HUMAN"){
    geneNames <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, unique(uniprots), 'SYMBOL', 'UNIPROT', multiVals == "first")
  }else if (species == "MOUSE"){
    geneNames <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, unique(uniprots), 'SYMBOL', 'UNIPROT', multiVals == "first")
  }else if (species == "RAT"){
    geneNames <- AnnotationDbi::mapIds(org.Rn.eg.db::org.Rn.eg.db, unique(uniprots), 'SYMBOL', 'UNIPROT', multiVals == "first")
  } else {
    stop("unrecognized species", species)
  }
  mapTable <- unique(data.table(uniprot = names(geneNames), gene = geneNames))
  setkey(mapTable, "uniprot")
  # this orders all and expands the missing cases
  mapTable <- mapTable[uniprots]
  # where gene lookup failed, assign the original uniprot
  mapTable[is.na(gene), gene := uniprot]
  return (mapTable$gene)
}


multiUniprots2multiGenes <- function (uniprots, sep = ";", species = "HUMAN", simplify = FALSE){
  toGenes <- data.table(uniprots = uniprots)
  toGenes <- toGenes[,.(singleUniprot = unlist(strsplit(uniprots, sep))),by = uniprots]
  toGenes[,singleGene := translateUniprot2GeneName(singleUniprot, species = species)]
  toGenes[is.na(singleGene), singleGene := singleUniprot]
  if (simplify == TRUE){
    simplify = function(x)unique(sort(x))
  }else if (simplify == FALSE){
    simplify = identity # do nothing
  }else if (! "function" %in% class (simplify)){
    stop("unexpected simplify format")
  }
  toGenes <- toGenes[, .(genes = paste(simplify (singleGene), collapse=sep)), by = uniprots]
  setkey(toGenes, uniprots)
  
  return(toGenes[uniprots, genes])
}


multiUniprotSites2multiGeneSites <- function (uniprotSites, sep = ";", siteSep = "_", species = "HUMAN", useDatFile = FALSE, uniqueOnly = FALSE){
  mapper <- data.table(uniprotSites = unique(as.character(uniprotSites)))
  #expand to singleSite, 1 per row
  mapper <- mapper[,.(singleSite = unlist(strsplit(uniprotSites, split = sep))), by = uniprotSites]
  mapper[,c("uniprot", "site") := tstrsplit(singleSite, split = siteSep)]
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
