



# A handful of functions that help organize output data and pdf images per notebook file in RStudio
# 
# Basic idea is each Rmd file (but should also work with .R etc) will have an associated filename_data/ directory
# relies on rstudioapi::getActiveDocumentContext() so will work everywhere that works


require(data.table)
require(devtools)
require(rstudioapi)



today <- function(){
  format(Sys.time(), "%Y_%m_%d")
}

DateFileName <- function(x){
  name <-   paste0(today(), "_", x)
  print (name)
  return (name)
}

#' Creates if needed and returns path to _data directory
#' 
ScriptNamedDir <- function(scriptName = NULL){
  if(is.null(scriptName))
    scriptName <- rstudioapi::getActiveDocumentContext()$path
  if (is.null (scriptName) || scriptName == "")
    stop("No script name found -- you may need to save this file first")
  outDir <- gsub(".R(md)?$", "_data", scriptName, ignore.case = TRUE)
  stopifnot( outDir != scriptName)
  if (!dir.exists(outDir)){
    message ("Creating directory associated with ", scriptName,", at ", outDir)
    dir.create(outDir)
  }
  return(outDir)
}

#' Given a simple file name, such as results.txt, this returns a path like
#' scriptName_data/todaysDate_results.txt
#' 
#' The date does not include time, so files will be over-written when generated on the same day,
#' but not between days
 
ScriptAndDatedFileName <- function(x, scriptName = NULL){
  dir <- ScriptNamedDir(scriptName)
  fileName <- DateFileName(x)
  path <- file.path(dir, fileName)
  cat (path, "\n")
  return (path)
}

#'  Given a simple file name, such as results.txt, search the data directory
#'  for the matching file name with the latest (most recent) matching name
#' @param x character file name
#' 
GetLatestScriptFile <- function(x, scriptName=NULL, include.dirs = FALSE){
  stopifnot (length(x) == 1)
  dir <- ScriptNamedDir(scriptName)
  filePattern <- paste0("^\\d{4}_\\d{2}_\\d{2}_", x, "$", collapse = "")
  filesFound <- list.files(dir, filePattern, include.dirs = include.dirs)
  stopifnot (length(filesFound) > 0)
  if (length(filesFound) > 1){
    message ("Multiple files  with matching names found.  Using the last one")
    print (filesFound)
  }
  return (file.path(dir, tail(filesFound, 1)))
} 

#' Quickly return a unique path to the scriptName_data/pdf/ directory
#' File names are guaranteed unique because they include date:time to the second
#' and when there is a clash a further counter suffix is added
#' @param prefix optional string to label pdf names
#' @subDir prefix an optional subdirectory name. Useful when a loop generates numerous related pdfs at a time
PDFBackupFileName <- function(prefix = "", subDir = ""){
  scriptDir <- ScriptNamedDir()
  imageDir <- file.path(scriptDir, "pdfs", subDir)
  if (!dir.exists(imageDir)) dir.create(imageDir, recursive = TRUE)
  now <- format(Sys.time(),  "%Y_%m_%d__%H_%M__%S")
  counter <- 0
  path <- file.path (imageDir, sprintf("%s%s.%02d.pdf", prefix, now, counter))
  while (file.exists(path)){
    counter <- counter + 1
    path <- file.path (imageDir, sprintf("%s%s.%02d.pdf", prefix, now, counter))
  }
  return (path)
}

BackupAsPDF <- function(graphics, prefix = "", subDir = "", dimensions = NULL){
  path <- PDFBackupFileName(prefix, subDir)
  if (is.null(dimensions))
    dimensions <- dev.size(units = "in")
  
  cat (sprintf("Writing image to:  \n%s\n", path))
  cairo_pdf(path, width = dimensions[1], height = dimensions[2])
  
  # handle functions, my enrichment heatmaps that are part of a list
  if ("function" %in% class(graphics)){
    graphics()
    g <- "finished" # something to print to console instead of graphics to device 
  }else if (! ("ggplot" %in% class(graphics) | "grob" %in% class(graphics) | "Heatmap" %in% class(graphics) | "HeatmapList" %in% class(graphics))){
    g <- graphics$hmList    
  }  else{
    g <- graphics
  }
  print (g)
  
  dev.off()
  return (graphics)
}


WriteSessionInfo <- function(path = NULL){
  if (is.null(path))
    path <- ScriptAndDatedFileName("SessionInfo.txt")
  si <- devtools::session_info()
  fileOut <- file(path, open = "wt")
  writeLines("─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────",
             con = fileOut)
  write.table(data.frame( value = unlist(si[[1]])), fileOut)
  #writeLines(capture.output(data.frame(value = unlist(si[[1]]))), con = fileOut) # for some mysterious reason, this always goes to the notebook output when run in a notebook. type = "message" is no help...[shrug]
  writeLines(capture.output(data.table(setting = names(si[[1]]), value = unlist(si[[1]]))), con = fileOut) # this too. Maybe an RStudio version issue
  writeLines("─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────",
             con  = fileOut)
  #write.table(si[[2]], fileOut)
  writeLines(capture.output(si[[2]]), con = fileOut)
  close(fileOut)
}

WriteInstalledPackages <- function (path = NULL){
  if (is.null(path))
    path <- ScriptAndDatedFileName("Installed.Packages.csv")
  package.mat <- installed.packages()
  fwrite (as.data.table(package.mat, keep.rownames = TRUE), path)
}
