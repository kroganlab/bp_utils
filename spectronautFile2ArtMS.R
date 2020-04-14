
library (data.table)
library (artMS)
library (yaml)
library (stringr)

# spectroNautFile2ArtMS 
# inputs: 
# filepath: path to a file with columns: 
#          BioReplicate, Condition, Run, IsotopeLabelType, ProteinName, PeptideSequence, PrecursorCharge, Intensity
# outFilePrefix:  (optional) a partial path or directory where artMS input files will be written. Most
#                 useful for running artmsQuantification
# artmsConfig : (optional) a nested list object that can be used to set options in the output config.yaml file
#
#
# output: a list with two data.frames representing artMS expectation of evidence_file and keys_file
#
#
# expected usage:
# ######## QC usage ###############
#  
# artmsInput <- spectronautFile2ArtMS("spectronautOutputFile.tsv")
# artmsQualityControlEvidenceBasic(artmsInput$evidence_file, artmsInput$keys_file)
# 
# do.call ("artmsQualityControlEvidenceBasic", artmsInput)
#
# ##### artmsQuantification usage#######
# artmsInput <- spectronautFile2ArtMS("experiment1SpectronautOutputFile.tsv", outFilePrefix = "artms/exp1")
#
# # or to change some configurations from the default:
# presetConfig <- list()
# presetConfig$output_extras$annotate$species=FLY
# artmsInput <- spectronautFile2ArtMS("experiment1SpectronautOutputFile.tsv", outFilePrefix = "artms/exp1", artmsConfig = optionalConfig)
#
# # before actually running this step you'll probably want to edit the automatically-generated contrasts file
# artmsQuantification (artmsInput$config_file)
#
######### PH quantification usage ###############
# artmsInput <- spectronautFile2ArtMS("spectronautOutputFile.tsv", outFilePrefix = "artms.ph/exp1")
# artmsInput <- doSiteConversion (artmsInput, referenceProteome = "my_reference_proteome.fasta", site="PH")
# artmsQuantification (artmsInput$config_file)

# issues ----
# 1.) FIXED: dummy columns are inserted, some meaningless. 
#            --extended QC will fail, it relies on columns not available in spectronaut output.
#
# 2.) config$data$filters$protein_groups == "remove" will fail because artMS needs a Leading Razor Protein column
#    to replace a group with a single protein. No such column exists (in the spectronaut files I've seen)
#

spectronautFile2ArtMS <- function (filePath, outFilePrefix = NULL, artmsConfig = NULL){
  sn <- fread (filePath)
  
  #bioreplicate needs to be unique across all conditions
  # make sure this is the case by pasting them together:
  sn[,BioReplicate := paste (Condition, BioReplicate, sep=".")]
  
  # There's a bug in artMS where a numeric BioReplicate causes problems
  # because the bioreplicate is used as a column header and data[[br]] ends up selecting columns
  # by numeric index instead of by name.
  # I work around this  here by making it a character vector always
  if (is.numeric(sn$BioReplicate)){
    sn[,BioReplicate := paste ("BioRep", BioReplicate, sep=".")]
  }
  
  keys <- unique(sn[,.(Condition, 
                       RawFile = Run, 
                       BioReplicate, 
                       IsotopeLabelType)])
  keys[,Run:=.I]
  evidence <- sn[,.(RawFile=Run, 
                    Proteins=ProteinName, 
                    `Leading proteins` = ProteinName,# artmsProtein2SiteConversion requires with a space Leading.proteins = ProteinName, 
                    Modified.sequence = convertModificationFormat(PeptideSequence), 
                    Charge = PrecursorCharge, 
                    Intensity,
                    
                    sequence = cleanModifiedSequence(PeptideSequence),
                    oxidation..m. = str_count(PeptideSequence, coll("[Oxidation (M)]")),
                    type = "SPECTRONAUT", #should be MSMS, MULTI-MSMS, MULTI-SECPEP...
                    ms.ms.count = 99,
                    
                    retention.length = rnorm(.N, mean=1, sd=0.1),
                    uncalibrated.mass.error..ppm. = rnorm(.N, mean=5, sd=2),
                    m.z = length(PeptideSequence) * 100/PrecursorCharge
                    )]
  
  if (!is.null(outFilePrefix)){
    paths <- makeOutFilePaths(outFilePrefix)
    if (!dir.exists(dirname(paths["output/results.txt"]))){
      dir.create(dirname(paths["output/results.txt"]), recursive=TRUE)
    }
    fwrite (evidence, paths["evidence.txt"], sep="\t")
    fwrite (keys, paths["keys.txt"], sep="\t")
    configData <- writeConfigFile(paths, artmsConfig)
    writeContrastFile(paths, keys)
  } else paths <- c()
  
  invisible(list (evidence_file = evidence, 
                  keys_file = as.data.frame(keys), 
                  config_file = paths["config.yaml"],
                  config_data = configData))
}


doSiteConversion <- function(artmsInput, referenceProteome,site="PH"){
  # do artms2SiteConversion to generate a new evidence.txt file and then generate a new config.yaml file
  # evidence.txt becomes evidence.site.txt
  # config.yaml becomes config.site.yaml
  
  # set new paths and confirm that they are new (confirm that gsub found something to work on)
  newEvidence <- gsub (".txt$", paste0(".",site,".txt"), artmsInput$config_data$files$evidence)
  stopifnot (newEvidence != artmsInput$config_data$files$evidence)
  
  newConfigFile <- gsub (".yaml$", paste0(".", site, ".yaml"), artmsInput$config_file)
  stopifnot (newConfigFile != artmsInput$config_file)
  
  #for UB and PH, artMS will filter based on Modification column
  artmsInput$evidence_file <- setModificationsColumns (artmsInput$evidence_file,site)
  
  artMS::artmsProtein2SiteConversion(artmsInput$evidence_file,
                              ref_proteome_file = referenceProteome,
                              column_name = "Proteins",
                              output_file = newEvidence,
                              mod_type = site
                              )
  #this gets the default config
  newConfig <- artmsWriteConfigYamlFile (config_file_name = NULL, verbose=FALSE)
  #set all values passed in from artmsInput...
  newConfig <- modifyList(newConfig, artmsInput$config_data)
  #new values here:
  newConfig$files$evidence <- newEvidence
  newConfig$data$filters$modifications <- site
  
  write_yaml(x = newConfig, file = newConfigFile)
  message ("New config file for ", site, " analysis written to ", newConfigFile)
  
  invisible (modifyList (artmsInput, list(evidence_file = newEvidence,
                                          config_file  = newConfigFile,
                                          config_data = newConfig)))
  
}

setModificationsColumns <- function (dt, site){
  dt[,Modifications := ""]
  for (s in site){
    if (s == "PH"){
      dt[grep ("\\(ph\\)", Modified.sequence),Modifications := paste (Modifications, "x Phospho")]
    } else if (s=="UB"){
      stop ("setModificationsColumns, I still need to check how spectronaut identifies ubiquitination")
      dt[grep ("\\(ph\\)", Modified.sequence),Modifications := paste (Modifications, "x GlyGly")]
    } else stop ("setModificationsColumns, unknown site requested: ", site)
  }
  invisible(dt)
}

convertModificationFormat <- function(specModSequence, mods=c("PH", "CAM", "MOX", "NAC")){
  result <- specModSequence
  specFormats <- list (PH='([STY])\\[Phospho \\(STY\\)\\]',
                       CAM = '([C])\\[Carbamidomethyl \\(C\\)\\]',
                       MOX = '([M])\\[Oxidation \\(M\\)\\]',
                       NAC =  '([A-Z])\\[Acetyl \\(Protein N-term\\)\\]')
  artmsFormats <- list (PH='\\1\\(ph\\)',
                        CAM = '\\1\\(cam\\)',
                        MOX = '\\1\\(ox\\)',
                        NAC = '\\1\\(ac\\)')
  stopifnot(names(specFormats)==names(artmsFormats))
  for (mod in mods){
    if (mod %in% names(specFormats)){
      result <- gsub(specFormats[[mod]], artmsFormats[[mod]], result)
    }else (stop("I don't know how to deal with requested mod: ", mod))
  }
  return (result)
}


cleanModifiedSequence <- function (modSeq){
  modSeq <- gsub("_", "",modSeq)
  modSeq <- gsub("\\[[A-Za-z ()-]*\\]", "", modSeq)
  return(modSeq)
}

makeOutFilePaths <- function (outFilePrefix){
  suffixes = c("config.yaml", "evidence.txt", "keys.txt", "output/results.txt", "contrast.txt")
  
  if (dir.exists(outFilePrefix)){
    paths <- file.path(outFilePrefix, suffixes)    
  }else{
    sep = "_"
    paths <- paste (outFilePrefix, suffixes, sep = sep)
  }
  names(paths) <- suffixes
  return (paths)
}

# main purpose is to write paths into a new artms config file
# presetConfig lets a user pass in a (possibly partial) nested list of artms configurations
# that will override the default artMS settings.  Note that any of the paths this function writes
# will overwrite the corresponding path passed in.
writeConfigFile <- function (paths, presetConfig=NULL){
  #this gets the default config
  artConfig <- artmsWriteConfigYamlFile (config_file_name = NULL, verbose=FALSE)

  # we know some things won't work, so turn them off
  #artConfig$qc$extended <- as.integer(0)
  artConfig$data$filters$protein_groups <- "keep"
  
  #let the user define some presets that will over ride the defaults in artMS
  if (!is.null(presetConfig)){
    if ("character" %in% class(presetConfig)){
      #looks like a path, load settings as a yaml file
      presetConfig <- read_yaml(presetConfig)
    }
    artConfig <- modifyList (artConfig, presetConfig)
  }
  
  artConfig$files$keys <- paths["keys.txt"]
  artConfig$files$evidence <- paths["evidence.txt"]
  artConfig$files$contrasts <- paths["contrast.txt"]
  artConfig$files$output <- paths["output/results.txt"]
  
  write_yaml(x = artConfig, file = paths["config.yaml"])
  invisible(artConfig)
}

writeContrastFile <- function(paths, keys){
  conditions <- unique(keys$Condition)
  contrasts <- unlist(lapply (conditions[conditions != min(conditions)], FUN=function(cond)(paste(cond, conditions[conditions < cond], sep="-"))))
  write (contrasts, file=paths["contrast.txt"])
  message ("Wrote ", length(contrasts), " contrast to file ", paths["contrast.txt"], " for all by all contrasts:" )
  message (paste("\t", contrasts, collapse="\n"))
}

