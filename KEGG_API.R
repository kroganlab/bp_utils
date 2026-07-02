
# if running as a script, modify this and then call it as last line in file.

MAIN <- function(){
  library (data.table)
  
  # define a vector of kegg compound ids, e.g. C00597
  inputKeggs <- fread ("unifiedMetaboliteToKeggToHmdbMap.tsv")$KEGG |> 
    trimws() |> 
    unique()
  
  c2r <- keggCpd2Reaction( unique(inputKeggs))
  r2ko <- keggReaction2KeggOrtholog( unique(c2r$reaction))
  ko2g <- keggOrtholog2genes(unique(r2ko$keggOrtholog))
  
  fwrite (c2r, "c2r.txt")
  fwrite (r2ko, "r2ko.txt")
  fwrite (ko2g, "ko2g.txt")
  
  
  merged <- merge (c2r[is.na(error)], r2ko[is.na(error)], by = "reaction", allow.cartesian = TRUE) |> 
    merge (ko2g[is.na(error)], by = "keggOrtholog", allow.cartesian = TRUE)
  
  fwrite (merged[, .N, by = .(compound, gene)], "c2g.merged.txt")
  
}




GLOBALSLEEP = 0.4 # "please limit to 3 times/second"

keggCpd2Reaction <- function (cpd){
  message (sprintf("Looking up %d queries at kegg will take at least %0.2f minutes", length(cpd), round(length(cpd) * GLOBALSLEEP / 60, 2)))
  
  .loadOne <- function (oneCpd){
    Sys.sleep(GLOBALSLEEP)    
    tryCatch({
      res <- fread (sprintf ("https://rest.kegg.jp/link/reaction/%s", oneCpd), header = FALSE)
      setnames(res, c("compound", "reaction"))
      res[, query := oneCpd]
      res
    },
    error = function(e) return(data.table(query = oneCpd, error = paste0(unlist(e), collapse = ";")))
    )
  }
  res <- lapply (unique(cpd[!is.na(cpd) & cpd != ""]), .loadOne) |> rbindlist(use.names = TRUE, fill = TRUE)
  if("reaction" %in% names(res)){
    
    res[, compound := gsub ("^cpd:", "", compound)]
    res[, reaction := gsub ("^rn:", "", reaction)]
  }
  return (res[])
}

keggReaction2KeggOrtholog <- function (reaction){
  message (sprintf("Looking up %d queries at kegg will take at least %0.2f minutes", length(reaction), round(length(reaction) * GLOBALSLEEP / 60, 2)))
  
  .loadOne <- function (one){
    Sys.sleep(GLOBALSLEEP)    
    tryCatch({
      res <- fread (sprintf ("https://rest.kegg.jp/link/ko/%s", one), header = FALSE)
      setnames(res, c("reaction", "keggOrtholog"))
      res[, query := one]
      res
    },
    error = function(e) return(data.table(query = one, error = paste0(unlist(e), collapse = ";")))
    )
  }
  res <-  lapply (unique(reaction[!is.na(reaction) & reaction != ""]), .loadOne) |>
    rbindlist(use.names = TRUE, fill = TRUE)
  if("keggOrtholog" %in% names(res)){
    res[, keggOrtholog := gsub ("^ko:", "", keggOrtholog)]
    res[, reaction := gsub ("^rn:", "", reaction)]
  }
  return (res[])
}

keggOrtholog2genes <- function (kegggOrtholog, species = "hsa"){
  message (sprintf("Looking up %d queries at kegg will take at least %0.2f minutes", length(kegggOrtholog), round(length(kegggOrtholog) * GLOBALSLEEP / 60, 2)))
  
  .loadOne <- function (one){
    Sys.sleep(GLOBALSLEEP)    
    tryCatch({
      res <- fread (sprintf ("https://rest.kegg.jp/link/%s/%s", species, one), header = FALSE)
      setnames(res, c("keggOrtholog", "gene"))
      res[, query := one]
      res
    },
    error = function(e) return(data.table(query = one, error = paste0(unlist(e), collapse = ";")))
    )
  }
  res <- lapply (unique(kegggOrtholog[!is.na(kegggOrtholog) & kegggOrtholog != ""]), .loadOne) |>
    rbindlist(use.names = TRUE, fill = TRUE)
  if("keggOrtholog" %in% names(res)){
    res[, keggOrtholog := gsub ("^ko:", "", keggOrtholog)]
    res[, gene := gsub (sprintf ("^%s:", species), "", gene)]
  }
  return (res[])
}


# MAIN()

