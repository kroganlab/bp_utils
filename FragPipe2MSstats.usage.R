

## This script takes an MSstats.csv file from FragPipe and converts it to be compatible with MSstats and artMS workflows.
## 2021-08-14


# Load libraries and source files
source("~/github/mehdi/proteomics/FragPipe2MSstats.R")

# Define input files
input_file = "Ph_FragPipe_MSstats.csv"
fasta_file = "2021-08-12-decoys-contam-Human-canonical-20200228__RSV-AstrainA2__PIV-5strainW3__EGFP.fasta.fas"

# Load file and make sure protein names are correct
ph = fread(input_file)
ph$ProteinName = gsub("^[^|]*|\\s*||[^|]*$", "", ph$ProteinName)
ph$ProteinName = gsub("\\|", "", ph$ProteinName)
ph$ProteinKey = paste(ph$ProteinName,ph$PeptideSequence,sep="__") # sep here should match protein_peptide_sep below (and not used in any protein names)

# Map peptides to sites within proteins
map_sites(ph,input_file,fasta_file, protein_peptide_sep = "__")

# Add keys and conditions properly to file, remove any intensities with NA
input_file = "Ph_FragPipe_MSstats_sitesmapped.csv"
keys_file = "keys.txt"
fragpipe_addkeys(input_file,keys_file)

# Example Keys file
# SampleNumber	Run	Condition	BioReplicate
# 1	ROADS_P1_RB1_1	Mock_8h	1
# 2	ROADS_P2_RB2_1	Mock_8h	2
# 3	ROADS_P3_RB3_1	Mock_8h	3
# 4	ROADS_P4_RB4_1	RSV_8h	1
# 5	ROADS_P5_RB5_1	RSV_8h	2
# 6	ROADS_P6_RB6_1	RSV_8h	3


