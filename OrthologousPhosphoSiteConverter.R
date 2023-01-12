# original file by Mehdi: do_convert_phosphorthologs_v2.R which wraps up some of my (Ben) code

# Load libraries  
library(Biostrings)
library(seqinr)
library(parallel)
library(tidyr)


do_convert_phosphorthologs = function(array_from, fasta_from, fasta_to, tab_orthologs, gaps=c("-"), siteAA = c("S", "T", "Y")){
  # array_from
  # Array of UNIPROT_SITE of sites you would like to convert  
  
  # fasta_from
  # Location of fasta file of species to convert From
  
  # fasta_to
  # Location of fasta file of species to convert To
  
  # table_orthologs, columns:
  #  From: uniprot ID in from species
  #  To: uniprot ID in to species
  # Can get from Biomart
  
  # siteAA
  # Amino acids you would like to systematically convert
  
  
  
  # Load fasta file
  message("LOADING FASTA FILES...")
  fas_from = Biostrings::readAAStringSet(fasta_from)  # this one gives more sequence coverage of the monkey sequences, with an additional about 28 seqs
  names(fas_from) = gsub("^(sp|tr)\\|(.*)\\|.*", "\\2", names(fas_from))
  
  fas_to = Biostrings::readAAStringSet(fasta_to)  # this one gives more sequence coverage of the monkey sequences, with an additional about 28 seqs
  names(fas_to) = gsub("^(sp|tr)\\|(.*)\\|.*", "\\2", names(fas_to))
  
  # Restrict table_ortholog to only consider proteins that you want converted now (to speed up algorithm)
  table_orthologs = table_orthologs[table_orthologs$From %in% gsub("_.*","",array_from),]
  
  # Create sequencePairList
  message("CREATING SEQUENCE PAIR LIST...")
  sequencePairList = lapply (1:nrow(table_orthologs), function(i)list(fas_from[[table_orthologs[i]$From]],
                                                                      fas_to[[table_orthologs[i]$To]]))
  
  # for (i in 1:length(sequencePairList)){
  #   pair = sequencePairList[[i]]
  #   Biostrings::pairwiseAlignment(pattern = pair[[1]], subject = pair[[2]])
  #   print(i)
  # }
  
  # Global alignment
  message("ALIGNING SEQUENCES PAIRWISE...")
  cl = makeCluster(detectCores() - 4)
  alignmentsGlobal = parallel::parLapply (cl,sequencePairList,
                                          fun=function(pair){
                                            if ( !is.null(pair[[1]]) & !is.null(pair[[2]]) ){
                                              a = Biostrings::pairwiseAlignment(pattern = pair[[1]], subject = pair[[2]])
                                              c(score=Biostrings::score(a), 
                                                patternChars = strsplit(toString(Biostrings::alignedPattern(a)), split=""),
                                                subjectChars = strsplit(toString(Biostrings::alignedSubject(a)), split="")
                                              )
                                            }
                                          })
  stopCluster(cl)
  names(alignmentsGlobal) <- paste(table_orthologs$From, table_orthologs$To, sep="_")
  
  # Make table that gets all sites of interest, annotates per position how it changes 
  message("CREATING FINAL OUTPUT TABLES...")
  allSites = lapply (alignmentsGlobal, getAllPhosphoSites, siteAA = siteAA) %>% rbindlist(idcol = "pair")
  allSites[, c("pattern", "subject") := tstrsplit(pair, split="_")]
  allSites[,pair:=NULL]
  colnames(allSites) = gsub("pattern","From",colnames(allSites))
  colnames(allSites) = gsub("subject","To",colnames(allSites))
  allSites$From_ProteinSite = paste0(allSites$From,"_",allSites$FromAA,allSites$FromPos)
  allSites$To_ProteinSite = paste0(allSites$To,"_",allSites$ToAA,allSites$ToPos)
  
  # Write out file
  fwrite(allSites, file="Table_Mapping_From-To.txt",sep="\t")
  
  # Extract sites that are specific to input
  table_out = allSites[allSites$From_ProteinSite %in% array_from, .(From_ProteinSite,To_ProteinSite)]
  fwrite(table_out, file="Table_MappedSites_From-To.txt",sep="\t")
  
  # Return
  message("ANALYSIS COMPLETE!")
  return(table_out)
  
  
}


getAllPhosphoSites <- function(alignment, gaps=c("-"), siteAA = c("S", "T", "Y")){
  #get a true position by counting cumulative gaps at each position
  patternPos <- ( seq_len(length(alignment$patternChars)) -
                    cumsum(alignment$patternChars %in% gaps) )
  subjectPos <- ( seq_len(length(alignment$subjectChars)) -
                    cumsum(alignment$subjectChars %in% gaps) )
  
  # this will make redundnat patternPos over the course of a gap
  # but you can check the chars to make sure we're not reading a gap
  
  data.table(patternChars = alignment$patternChars, 
             subjectChars=alignment$subjectChars, 
             patternPos, 
             subjectPos)[patternChars %in% siteAA | subjectChars %in% siteAA,
                         .(patternPos, patternAA = patternChars, subjectPos, subjectAA = subjectChars)]
}


# # # EXAMPLE USAGE
# # 
# # Define fasta files
# fasta_from = "~/github/mehdi/data/proteomes/uniprot-musmusculus_3AUP00000058-2023.01.10-20.03.40.85.fasta"
# fasta_to = "~/github/mehdi/data/proteomes/2023-01-11_HomoSapiens_uniprot-proteome_UP000005640.fasta"
# 
# # Load table_from
# library(data.table)
# D = fread("ph/dir/evidence.PH-mapping.txt")
# D.exp = D[!is.na(D$mod_sites), .(site = unlist(strsplit(mod_sites, ";"))), by = .(mod_seqs, Proteins, mod_sites)]
# array_from = unique(D.exp$site)
# 
# # Load tab_orthologs
# table_orthologs = fread("~/github/mehdi/data/mart/Table_Convert_Human-Mouse.txt")
# table_orthologs = table_orthologs[,.(uniprot_mouse,uniprot_human)]
# colnames(table_orthologs) = c("From","To")
# table_orthologs = unique(table_orthologs[table_orthologs$From!="" & table_orthologs$To!="", ])
# 
# table_out = do_convert_phosphorthologs(array_from, fasta_from, fasta_to, table_orthologs, gaps=c("-"), siteAA = c("S", "T", "Y"))



