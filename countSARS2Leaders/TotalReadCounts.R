library (data.table)

readCounts.map <- fread ("sgeOut/readCounts.map", header = FALSE)
fastq.map <- fread ("sgeOut/fastq.map", header = FALSE)

fastq.map[, c("taskFile", "fastqFile") := tstrsplit(V1, ":")]
readCounts.map[, taskFile := tstrsplit(V1, ":")[[1]]]
readCounts.map[, totalReads := V2]

totalReadCounts <- merge (readCounts.map[, .(taskFile, totalReads)],
                     fastq.map[, .(taskFile, fastqFile)],
                     all = TRUE)
totalReadCounts[, file := gsub ("\\.fastq(\\.gz)?", "", basename(fastqFile))]

fwrite (totalReadCounts, "totalReadCounts.csv")
