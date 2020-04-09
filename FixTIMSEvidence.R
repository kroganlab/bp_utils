

library (data.table)



oxidizeMs <- function (Sequence, count){
  chars <- unlist(strsplit(Sequence, ""))
  ms <- which (chars == "M")
  chars[ms[1:count]] <- "M(ox)" 
  return (paste(chars, collapse=""))
}

modifySequence <- function (Sequence, Modifications){
  mods <- unlist(strsplit(Modifications, ","))
  for (mod in mods){
    words <- unlist(strsplit(mod, " "))
    if (!suppressWarnings(is.na(as.integer(words[1])))){
      count <- as.integer(words[1])
      mod <- paste(words[2:length(words)], collapse=' ')
    } else count <- 1
    if (mod == "Acetyl (Protein N-term)"){
      Sequence <- paste ("(ac)", Sequence, sep="", collapse="")
    } else if (mod == "Oxidation (M)"){
      Sequence <- oxidizeMs(Sequence, count)
    } else (message ("unrecognized mod type ", mod))
    
  }
  #finally put _ on ends
  paste (c("_", "_"), collapse = Sequence)
}


fixTIMSEvidenceFile <- function (fileIn, fileOut){
  
  tims <- fread(fileIn) # ("~/Box/Issues/Evidence_issue/evidence_TIMS.txt")
  
  message ("The input file contains ", sum(tims$`Modified sequence`==""), " empty Modified sequences")
  message ("This script will try to fill them in using information from Mod peptide ID column")
  message ("Where that fails, it will try to generate the Modified sequence based on the Sequence and Modifications column")
  message ("")
  
  # there are no empty Mod. 
  message ("Checking for empty Mod. peptide ID")
  stopifnot (nrow(tims[is.na(`Mod. peptide ID`),]) == 0)
  message ("...all rows have a Mod. peptide ID")
  message ("")

    
  
  # there are only 33 redundant mod peptide ids, ambiguous at Met-ox site
  redundantModPeptideIDs <- tims[`Modified sequence`!= "", .(length(unique(`Modified sequence`)), paste(unique(`Modified sequence`), collapse=" ")),by=`Mod. peptide ID`][V1>1,]
  message ("There are ", nrow(redundantModPeptideIDs), " redundant mod peptide IDs, examples below. The most commonly seen in this evidence file will be used")
  message (head (redundantModPeptideIDs$V2))
  message("")
  
  
  idMapping <- tims[`Modified sequence`!= "", .(seq = paste(unique(`Modified sequence`), collapse=" _or_ ")), by = `Mod. peptide ID`]
  
  ambiguousIDs <- idMapping[grep ("_or_", seq),]$`Mod. peptide ID`
  
  #choose most commonly seen oxidation position in our data for the ambiguous ones:
  mostCommon <- tims[`Modified sequence`!= "" & `Mod. peptide ID` %in% ambiguousIDs,
                     .( .N), 
                     by = .(`Modified sequence`, `Mod. peptide ID` )
                     ][,
                       .(seq = `Modified sequence`[which.max(N)]),
                       by=`Mod. peptide ID`]
  
  idMapping <- merge (idMapping, mostCommon, by = "Mod. peptide ID", all.x=TRUE, suffixes = c("", ".y"))
  idMapping[grep ("_or_", seq), seq := seq.y]
  idMapping$seq.y <- NULL
  
  
  
  #merge back to tims:
  
  tims2 <- merge (tims, idMapping, by = "Mod. peptide ID", all.x=TRUE)
  
  
  #some are still missing
  message ("After using Mod. peptide ID, ", nrow(tims2[is.na(seq),]), " rows still have missing Modified sequences, their modifications shown below. Now attempting to generate these Modified Sequences now")
  message(tims2[is.na(seq),]$Modifications %>% unique())
  message("")
  
  
  missingModPeptides <- tims2[is.na(seq), .(`Mod. peptide ID`, Sequence, Modifications)] %>% unique()
  seqs <-sapply (seq_len(nrow(missingModPeptides)), function(i){modifySequence (missingModPeptides[i,]$Sequence, missingModPeptides[i,]$Modifications)})
  missingModPeptides[,generatedSeq:=seqs]
  
  
  tims2 <- merge (tims2, missingModPeptides[,.(`Mod. peptide ID`, generatedSeq)], by = "Mod. peptide ID", all.x=TRUE)
  
  
  tims2[`Modified sequence`== "", `Modified sequence`:=seq]
  tims2[is.na(`Modified sequence`), `Modified sequence`:= generatedSeq]
  
  message ("There are now ", nrow(tims2[is.na('Modified sequence'),]), " rows with empty Modified Sequence columns")
  
  fwrite (tims2[,colnames(tims), with=FALSE], fileOut, sep="\t")
  message ("Complete!  File written to ", fileOut)
}




