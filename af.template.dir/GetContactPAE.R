
pdbFile = character(0)
paeJsonFile = character(0)
contactResOutFile = character(0)

library (data.table)
library (bio3d)
library (jsonlite)



args = commandArgs(trailingOnly=TRUE)

pdbFile = grep("pdb$", args, value = TRUE)
stopifnot(`Pass only one pdb file` = length(pdbFile) <= 1)

paeJsonFile = grep("json$", args, value = TRUE)
stopifnot(`Pass only one PAE json file` = length(pdbFile) <= 1)

if (length(pdbFile) == 0 & length(paeJsonFile) == 1){
  # missing pdb
  # relaxed_model_1_multimer_v3_pred_0.pdb
  pdbFile = gsub("^(.*)result_(model_[0-9]+_multimer_v3_pred_[0-9]+).predicted_aligned_error.json$",
                 "\\1relaxed_\\2.pdb",
                 paeJsonFile)
}

if (length(pdbFile) == 1 & length(paeJsonFile) == 0){
  # missing json
  # result_model_1_multimer_v3_pred_0.predicted_aligned_error.json
  paeJsonFile = gsub("^(.*)relaxed_(model_[0-9]+_multimer_v3_pred_[0-9]+).pdb$",
                     "\\1result_\\2.predicted_aligned_error.json",
                     pdbFile)
}


contactResOutFile <- sprintf ("%s.contacts.csv",
                              gsub ("\\.pdb$", "",pdbFile)) # strip off the pdb


message ("Loading pdb from ", pdbFile)
message ("Loading pae from ", paeJsonFile)
message ("Writing contacts to ", contactResOutFile)



loadPDBAtoms <- function(path){
  pdb <- bio3d::read.pdb(path)
  atoms <- pdb$atom
  setDT(atoms)
  atoms[, idx := .I]
  return (atoms[])
}


interChainContacts <- function (pdbFile){
  atoms <- loadPDBAtoms(pdbFile)

  # all by all atom distance
  message (sprintf("All by all atom distance for %d atoms", nrow(atoms)))
  atomDistance <- as.matrix(dist(atoms[, .(x,y,z)]))
  
  # sweep out row and col radii
  message (sprintf("Done with distance, now calucating contacts"))
  vdw.radii <- c(H =  1.2, C =  1.7, N =  1.55, O =  1.52, S =  1.8)
  atomDistance <- sweep (atomDistance, 1, vdw.radii[atoms$elesy])
  atomDistance <- sweep (atomDistance, 2, vdw.radii[atoms$elesy])
  
  # if remaining distance is still less than 0.5, we declare it a contact 
  contactIndeces <- which (atomDistance < 0.5, arr.ind = TRUE) |> as.data.table()
  
  # label with chains from idx in atoms table
  contactIndeces[atoms, chainRow := i.chain , on =  c(row = "idx")]
  contactIndeces[atoms, chainCol := i.chain , on =  c(col = "idx")]
  
  # make crosschain only, and only in one direction:
  contactIndeces <- contactIndeces[chainRow < chainCol]
  
  # label with resno from atoms table
  contactIndeces[atoms, resnoRow := i.resno, on = c(row = "idx")]
  contactIndeces[atoms, resnoCol := i.resno, on = c(col = "idx")]
  
  # collapse from atoms to residues and sort
  contactRes <- setorder(unique(contactIndeces[, .(chainRow, resnoRow, chainCol, resnoCol)]))
  
  # translate per-chain resno to the multimer resno based on chain lengths (max(resno)) for all prior chains
  # assumptions!!!
  cl <- atoms[, .(l = max(resno)), by = chain]
  contactRes[, mmerResnoRow := resnoRow + sum(cl[chain < chainRow, sum(l)]), by = chainRow]
  contactRes[, mmerResnoCol := resnoCol + sum(cl[chain < chainCol, sum(l)]), by = chainCol]
  
  
  return (contactRes[])  
}



pae <- jsonlite::fromJSON( paeJsonFile)

contactRes <- interChainContacts(pdbFile)

contactRes[, pae := pae[mmerResnoRow, mmerResnoCol], by = .(mmerResnoRow, mmerResnoCol)]


fwrite (contactRes,file = contactResOutFile)

