#a function to take the various names for raw file columns and standardize them
.fixRawFileColumnName <- function(datTab, newName = "RawFile"){
  rawFilePattern = "Raw[. ]?(F|f)ile"
  oldName = grep(rawFilePattern, names(datTab), value=TRUE)
  stopifnot(length(oldName) == 1)
  if (oldName!=newName){
    message("Converting column name ", oldName, " to ", newName)
    setnames(datTab, oldName, newName)
  }
}
#prepares for processData 
prepareDataForMSStats = function(evidenceFile, keysFile, outfile=NULL){
  

  #read evidence and fix a column name
  if ("data.table" %in% class(evidenceFile)){
    ev <- evidenceFile
  }else{
    ev = fread (evidenceFile)
  }
  .fixRawFileColumnName(ev)
  
  #read keys and fix a column name
  if ("data.table" %in% class(keysFile)){
    keys <- keysFile
  } else{
    keys = fread (keysFile)
  }
  .fixRawFileColumnName(keys)
  
  #merge
  mergeColumns = c("RawFile")
  if ('IsotopeLabelType' %in% colnames(ev)){
    mergeColumns = c(mergeColumns, "IsotopeLabelType")
  }
  # see artMS::artMSMergeEvidenceAndKeys for more checks, but this suffices for most
  evAndKeys = merge (ev, keys, by = mergeColumns, suffixes = c("Ev", "Keys"))  #
  
  #filter
  before = nrow(evAndKeys)
  #filter evidence (empty protein names, groups, contaminants, reversed decoys)
  #evAndKeys <- evAndKeys[!grepl("^$|;|CON__|REV__", Proteins),                       ]
  evAndKeys <- evAndKeys[!grepl("^$|CON__|REV__", Proteins),]
  after = nrow(evAndKeys)
  message ("Filtered ", before-after, " rows from keyed evidence , new size = ", after, "  (empty protein names, contaminants, reversed decoys)")
  evAndKeys <- evAndKeys[!is.na(Intensity)]
  evAndKeys <- evAndKeys[Intensity > 0]
  after2 = nrow(evAndKeys)
  if(after2 < after)
    message ("Filtered ", after-after2, " rows with NA intensity or intensity <= 0, new size = ", after2)
  
  # in cases of  concatenated evidence files, we might still have non-unique-to-protein peptides
  nonUniquePeptides <- evAndKeys[,.(numProteins = length(unique(Proteins))), by = `Modified sequence`][numProteins >1]$`Modified sequence`
  if (length(nonUniquePeptides) > 0){
    evAndKeys <- evAndKeys[!`Modified sequence` %in% nonUniquePeptides]
    after3 <- nrow(evAndKeys)
    message ("Filtered ", after2-after3, " rows with non-specific-to-protein peptides (", length(nonUniquePeptides),"), new size ", after3)
  }

  
  #convert to msstatsFormat

  # main steps
  # 1 sum up intensity for duplicate features
  # 2 make table complete: add rows with NA intensity for missing values.
  # 3 get columns and names suitable for MSstats

  #1
  evAndKeys.noDups <- evAndKeys[,.(Intensity = sum(Intensity), count.features = .N),
                                by = .(Proteins, 
                                        `Modified sequence`,
                                        Charge,
                                        IsotopeLabelType,
                                        Condition,
                                        BioReplicate,
                                        Run)]
  countSummed <- sum(evAndKeys.noDups$count.features > 1)
  if (countSummed > 0)
    message (countSummed, " rows (peptide_ion x run) had >1 feature; using summed intensity for these")
  
  #2
  # define the two things to do a full cross join on, add dummy column to allow the full cross join
  full.peptide.ions <- unique (evAndKeys.noDups[, .(Proteins, 
                                                   `Modified sequence`,
                                                   Charge,
                                                   dummyIndex = 1)])
  full.runs <- unique(evAndKeys.noDups[,.(IsotopeLabelType,
                                          Condition,
                                          BioReplicate,
                                          Run,
                                          dummyIndex  = 1)])
  # do the full cross join
  full.peptide.ions.runs <- merge (full.peptide.ions, full.runs, by = "dummyIndex",allow.cartesian=TRUE)[,dummyIndex := NULL][]
  message (nrow(full.peptide.ions), " peptide ions in at least one of ", nrow(full.runs), " runs")
  # merge back to the actual intensity data
  evAndKeys.noDups.complete <- merge (full.peptide.ions.runs, evAndKeys.noDups, all.x=TRUE, 
                                      by = intersect(colnames(full.peptide.ions.runs),
                                                     colnames(evAndKeys.noDups)))
  
  #3 in msstats format
  evAndKeys.mss <- evAndKeys.noDups.complete[,.(ProteinName = Proteins,
                                                PeptideSequence = `Modified sequence`,
                                                PrecursorCharge = Charge,
                                                FragmentIon = NA,
                                                ProductCharge  = NA,
                                                IsotopeLabelType,
                                                Condition,
                                                BioReplicate,
                                                Run,
                                                Intensity)]
  
  if (!is.null(outfile)){
      fwrite (evAndKeys.mss, outfile, sep = "\t")
  }
  evAndKeys.mss  
}
