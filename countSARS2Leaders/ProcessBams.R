library (data.table)
library (Biostrings)
#library (pblapply)
library (Rsamtools)

#numProcessors = 2 # for genome alignments

# generate alignment table that maps positions between sequences
ref.genome.fa <- "./SARS-CoV-2.fa"
genomesTable <- fread ("fastqList.txt", header = FALSE)
setnames(genomesTable, new = c("fastq", "fasta"))
genomesTable[, fastqName := tstrsplit(basename(fastq), "\\.")[[1]]]

genomes <- unique(genomesTable$fasta)

mapPositions <- function(otherFasta, refFasta){
  refSeq <- Biostrings::readDNAStringSet(refFasta)[[1]]
  otherSeq <- Biostrings::readDNAStringSet(otherFasta)[[1]]
  aligned <- Biostrings::pairwiseAlignment(pattern = refSeq, subject = otherSeq)
  patternChars = strsplit(toString(Biostrings::alignedPattern(aligned)), split="")[[1]]
  subjectChars = strsplit(toString(Biostrings::alignedSubject(aligned)), split="")[[1]]
  
  alignedDT <- data.table (ref = patternChars, other = subjectChars)
  alignedDT[, ref.isAligned := ifelse(ref == "-", 0, 1)]
  alignedDT[, other.isAligned := ifelse(other == "-", 0, 1)]
  alignedDT[, ref.pos := cumsum(ref.isAligned)]
  alignedDT[, other.pos := cumsum(other.isAligned)]
  
  mapped <- alignedDT[, .(ref, other, ref.pos, other.pos)]  
  mapped[, aligned := ! (ref == "-" | other == "-")]
  return (mapped)
}



message ("Aligning ", length(genomes), " genomes to reference genome at ", ref.genome.fa)
#genomesMapped <- pbapply::pblapply(genomes, mapPositions, refFasta = ref.genome.fa, cl = numProcessors)
genomesMapped <- lapply(genomes, mapPositions, refFasta = ref.genome.fa )
names(genomesMapped) <- genomes


genomesMapped <- rbindlist (genomesMapped, idcol = "fasta")
#genomesMapped[, c("exp", "cell", "virus") := tstrsplit(fasta, "[-._]")[1:3]]

fwrite (genomesMapped, "genomesMapped.csv")
message ("genomesMapped table written to genomesMapped.csv")

# process bam files
dir <- "./withLeaders"
bamFiles <- list.files (dir, pattern = "\\.bam$")


bamToTable <- function(bamFile){
  bam <- Rsamtools::scanBam(bamFile)
  stopifnot (length(bam) == 1)
  bam <- bam[[1]]
  
  #store names of BAM fields
  bam_field <- names(bam)
  
  bam$seq <-  as.character(bam$seq)
  bam$qual <-  as.character(bam$qual)
  
  flags <- as.data.table (Rsamtools::bamFlagAsBitMatrix(bam$flag))
  
  
  #store as data table
  bam_df <- do.call("data.table", bam)
  setnames(bam_df, bam_field)
  bam_df <- cbind(bam_df, flags)
  
  return (bam_df)
}

startPosFromBam <- function(bamFile){
  bt <- bamToTable(bamFile)
  bt[strand == "+", .N, by = pos]  # POSITIVE STRANDS ONLY. Because TRIM LEFT by bbmap only works for positive strands. THIS IS LIKELY RNASEQ METHOD SPECIFIC...
}


names(bamFiles) <- tstrsplit(basename(bamFiles), "\\.")[[1]]

message ("Loading ", length(bamFiles), " bamFiles from ", dir)
startPosTables <- lapply (file.path(dir, bamFiles), startPosFromBam )


names(startPosTables) <- names(bamFiles)
startPos <- rbindlist(startPosTables, idcol = "file")
startPos[]

#startPos[, c("exp", "cell", "virus", "time", "rep", "r1or2") := tstrsplit(file, "[-_]")]


fwrite (startPos, "StartPositionCounts.csv")
message ("Counts of leader splice sites written to StartPositionCounts.csv")


# connect the startPos rows to the FASTA used...
startPos[genomesTable, genomeFile := i.fasta, on = c(file = "fastqName")]
# ...so we can get the right ref.pos from genomesMapped
startPos[genomesMapped[aligned == TRUE], ref.pos := i.ref.pos, on = c(genomeFile = "fasta", pos = "other.pos")]


# are any great number of htings not matched?
message ("There were ",          startPos[is.na(ref.pos), .N], 
         " positions found in ", startPos[is.na(ref.pos), sum(N)],
         " reads that could not map to reference genomes")

# label the sites according to reference fasta
spliceSites.ref <- c(genomic = 67,
                     #nonc_01 = 15777,
                     #nonc_02 = 21053,
                     S = 21553,
                     orf3 = 25382,
                     E  = 26237,
                     M = 26470,
                     orf6 = 27041,
                     orf7 = 27385,
                     #nonc_03_7b = 27674,
                     #nonc_04_7b = 27762,
                     orf8 = 27885,
                     N = 28257,
                     nonc_9b = 28280,
                     Nstar = 28878)

# above data in table form and create windows:
splices <- rbindlist(list(SARS2  = data.table (splice.name = names(spliceSites.ref),
                                               splice.center = spliceSites.ref)),idcol = "refDB")
splices[, windowRight := splice.center + 5]
splices[, windowLeft := splice.center -5 ]


# do the weird data.table joins using windows.  If windowLeft/windowRight straddle ref.pos then it is a match.
startPos.labeled <- splices[startPos,
                            ,
                            on = .( windowLeft < ref.pos, windowRight > ref.pos)]


# windowLeft and windowRight will take the value of what they are compared to.
# confirm behavior is as expected:
stopifnot (all(startPos.labeled$windowLeft == startPos.labeled$windoRight))

# fix some names, clean up unnecessary columns
startPos.labeled[, refPos := windowLeft]
startPos.labeled[, windowLeft := NULL]
startPos.labeled[, windowRight := NULL]
startPos.labeled[, refDB := NULL]

fwrite (startPos.labeled,  "Subgenomic_Starts_Labeled.csv")
message ("Counts of leader splice sites matched to known splice sites written to StartPositionCounts.csv")


