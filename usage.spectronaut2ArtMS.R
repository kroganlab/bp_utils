

# define some settings that will over-write defaults in artMS config file:
cf<- list()
cf$output_extras$annotate$enabled <- as.integer(0) # turn off any attempt to annotate proteins (only works with HUMAN)
cf$qc$extended <- as.integer(1)  # extended QC plots

#cf$msstats$normalization_method <- FALSE # turn off normalization in MSstats...useful if I've already normalized the evidence file


#load the next couple of functions from my file on github...modify the path as appropriate for your working directory
source ("../../../bp_utils/spectronautFile2ArtMS.R")

globalInput <- spectronautFile2ArtMS(inputEvidenceFile,   # a path to a spectronaut "evidence" file
                                     outFilePrefix = "UB.full/2020_05_18.no24Hr",   # the prefix for file names of generated input (config,evidence,keys,contrast) files and output directory.
                                     artmsConfig = cf,  # the config object we created above with different-from-default settings
                                     # contrastPatterns (below) is a vector of contrasts and/or regular expressions that expand to contrasts
                                     # below I use slightly complex regular expressions that expand to the contrasts in the comment below
                                     contrastPatterns = c("^(..Hr)_(Inf|Mock)_MG132-\\1_\\2$",   # non-MG132 versus correspnoding MG132 case
                                                          "^(..Hr)_Inf((_MG132)?)-\\1_Mock\\2$")) # infected versus corresponding mock case

############################################
                # contrasts from above regular expressions
                # 06Hr_Mock_MG132-06Hr_Mock
                # 06Hr_Inf_MG132-06Hr_Inf
                # 12Hr_Mock_MG132-12Hr_Mock
                # 12Hr_Inf_MG132-12Hr_Inf
                # 
                # 06Hr_Inf-06Hr_Mock
                # 06Hr_Inf_MG132-06Hr_Mock_MG132
                # 12Hr_Inf-12Hr_Mock
                # 12Hr_Inf_MG132-12Hr_Mock_MG132
##########################################



# if this is an AB dataset, we're ready to run artMSQuantification (skip this for UB):
artmsQuantification(globalInput$config_file)

# its also possible to just run QC instead:
artmsQualityControlEvidenceBasic(globalInput$config_data$files$evidence, globalInput$config_data$files$keys)


# if this is a UB dataset, do the further work of grouping peptides by UB sites:
siteInput <- doSiteConversion(globalInput,   #the output of spectronautfile2ArtMS
                              referenceProteome = "data/monkeyAndVirus.prot.2.fasta", 
                              site="UB")  #alternatively "PH"

# now do the artmsQuantification for UB.  
# Note that we use siteInput here instead of globaInput
artmsQuantification(siteInput$config_file)

