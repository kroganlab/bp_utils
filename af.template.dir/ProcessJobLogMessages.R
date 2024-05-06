library (data.table)

# relies on R > 4.1 for native pipe operator: |>

main <- function(args = commandArgs(trailingOnly=TRUE)){
  if (length(args) == 0){
    f <- file("stdin")
    text <- paste0(readLines(f), collapse = "\n")
    close(f)
    
    dt <- parseJobMessagesToTable(text)  #reads from stdin
  }
  else dt <- parseJobMessagesToTable(args[1])
  outFile <- ifelse(is.na(args[2]), "./jobInfo.csv", args[2])
  cat (c("Writing output to ", outFile, "\n"))
  fwrite (dt, outFile) # writes to stdout, redirect as needed
}


parseJobMessagesToTable <- function (jobLogPathOrString = "2024_05_03_AF_Jyoti3_data/jobInfo.txt"){
  
  jobMessages <- fread (jobLogPathOrString,
                        sep = ']',
                        header = FALSE,
                        col.names = c("jobInfo", "message"))
  
  
  jobMessages[, logBaseName := tstrsplit(jobInfo, ".log")[[1]]]
  
  # relies on message being written in order
  jobMessages[, logMessageOrder := .SD[, .I], by = logBaseName]
  # relies on messages being written in order
  # it simply increments the ID every time it reads a new "monomer pipeline" message
  jobMessages[, msaID := cumsum(grepl("Running monomer pipeline on chain", message))]
  
  
  jobAndTaskFromLogFile <- function(logBaseName){
    parts = tstrsplit(logBaseName, "-")
    lapply(parts[(length(parts)-1):length(parts)],
           as.integer)
  }
  jobMessages[, c("jobID", "taskID"):= jobAndTaskFromLogFile(logBaseName)]
  
  # random seed
  seeds <-  jobMessages[grepl ("Using random seed", message), 
                        .(logBaseName, jobID, taskID, seedStr = tstrsplit(message, " " )[[4]])] # save seeds as strings, because they are too large for R ints, avoid accidents of range
  
  
  # chains and sequences
  chainsLong <- jobMessages[grepl ("Running monomer pipeline", message),
                            c(tstrsplit(message, ":? ", keep = 6:7, names = c("chainID", "sequenceName")), list(msaID = msaID)),
                            by = .(logBaseName, jobID, taskID)]
  
  chains <- chainsLong |> 
    dcast(logBaseName+jobID+taskID~chainID, value.var = "sequenceName")
  
  
  
  
  # MSA info
  
  msaSizeProcess <- function(message, msaID){
    msaType <- tstrsplit(message, " size: ")[[1]]
    msaType <- gsub ("[ )(]+","_",msaType)
    msaSize <-   as.integer( gsub ("^.* MSA size: ([0-9]+) sequences.*$", "\\1", message)) 
    list(msaID = msaID, msaType = msaType, msaSize = msaSize)
  }
  msaSizes <- jobMessages[grepl("MSA size", message), 
                          msaSizeProcess(message, msaID),
                          by = .(logBaseName, jobID, taskID)]
  
  msaSizes[chainsLong, c("chainID", "sequenceName") := .(i.chainID, i.sequenceName), on = "msaID" ]
  ## convert to wide
  msaSizes <- dcast (msaSizes, logBaseName+jobID+taskID~msaType+chainID, value.var = "msaSize")
  
  
  # template info
  
  templateCounts <- jobMessages[grepl("Total number of templates", message),
                                .(possibleTemplates = gsub ("^Total number of templates.* ([0-9]+).$", "\\1", message)),
                                by = .(logBaseName, jobID, taskID, msaID)]
  templateCounts[chainsLong, c("chainID", "sequenceName") := .(i.chainID, i.sequenceName), on = "msaID" ]
  templateCounts[, chainID := gsub ("^", "templates_", chainID)]
  templateCounts <- dcast (templateCounts, logBaseName+jobID+taskID~chainID, value.var = "possibleTemplates")
  
  
  # model start
  # Running model model_1_multimer_v3_pred_0 on Q9Y2L1__RaTG13_orf8
  
  modelStarts <- jobMessages[grepl("^Running model model_", message), 
                             tstrsplit( message, " ", keep = c(3,5), names = c("model", "pair")),
                             by = .(logBaseName, jobID, taskID)]
  
  
  
  # model complete
  # Total JAX model model_1_multimer_v3_pred_0 on Q9Y2L1__RaTG13_orf8 predict time (includes compilation time, see --benchmark): 440.9s
  modelEnds <- jobMessages[grepl("^Total JAX model model_", message), 
                           tstrsplit( message, " ", keep = c(4,6, 14), names = c("model", "pair", "predictTime")),
                           by = .(logBaseName, jobID, taskID)]
  
  modelEnds[, predictTime := as.numeric(gsub("s$", "", predictTime))]
  
  # merge all together
  # these first four tables (3 merges) keep only one row per job/task
  result <- merge (seeds, chains) |>
    merge(msaSizes, by = c("logBaseName", "jobID", "taskID")) |> 
    merge (templateCounts, by = c("logBaseName", "jobID", "taskID") )  |>
    
    # these next two expand to multiple models per task
    merge (modelStarts) |>
    merge ( modelEnds, by = c("logBaseName", "jobID", "taskID", "model", "pair"), all = TRUE)
  return (result)
  
}


main(args = commandArgs(trailingOnly=TRUE))

