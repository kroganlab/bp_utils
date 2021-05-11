





#' Compute linear models, one model per protein or gene, and for each formula.
#'
#' @param fullDataTable data.table.  Required columns are dictated by formulas and splitColumn
#' @param formulaList A list of formulas. They will be computed separately.
#' @param splitColumn  character indicating name of column in fullDataTable that defines the units to compute
#'                    the linear model on.  Usually "Protein" or  "gene", etc.

linearModelsAllProteins <- function (fullDataTable, formulaList, splitColumn = "Protein"){
  subsetTables <- split (fullDataTable, fullDataTable[[splitColumn]])
  
  # do al linear models, anova, etc. Result is a list of [[anova, coef]] lists
  coef.anova.list <- lapply (subsetTables, function(sDT){.linearModelsOneProtein(sDT, formulaList)})

  # get all anova in a  single table
  anovas.list <- lapply(coef.anova.list, `[[`, "anova")
  anovasTable <- rbindlist (anovas.list, idcol = splitColumn, fill = TRUE)
  setnames(anovasTable, old = c("rn", "Pr(>F)"), new = c("term", "p.value"))
  
  # get all coef in a single table
  coefs.list <- lapply(coef.anova.list, `[[`, "coef")
  coefTable <- rbindlist(coefs.list, idcol = splitColumn, fill = TRUE)
  
  return (list (anova = anovasTable, coef = coefTable))
} 





.linearModelsOneProtein <- function(subsetDT, formulas ){
  lms <-tryCatch({
    lapply (formulas, function(f)lm(f, data = subsetDT))
  }, error = function(err){
    # halt execution, but return the error message.
    print (err)
    return(err)
  })
  if("error" %in% class(lms)){
    return (data.table(error = as.character(lms)))
  }  
  
  coef.list <- lapply (lms, function(lm.out)as.data.table(coefficients(summary(lm.out)), keep.rownames = TRUE))
  coef.table <- rbindlist(coef.list, idcol = "model")
  setnames (coef.table, old = c("rn", "Pr(>|t|)"), new = c("term", "p.value"))
  
  #f.tests from anova function
  anova.list <- lapply (lms, function(l)as.data.table(anova(l), keep.rownames=TRUE))
  anovaTables <- rbindlist(anova.list, idcol = "model")
  
  anovaTables[, sigCode := dplyr::case_when(`Pr(>F)`< 0.001 ~ "***",
                                            `Pr(>F)`< 0.01  ~ "**",
                                            `Pr(>F)`< 0.05  ~ "*",
                                            `Pr(>F)`< 0.1  ~ ".",
                                            TRUE ~ "")]
  return (list (anova = anovaTables[], coef = coef.table))
}



#' this function holds a piece of sample code that demonstrate how to use LinearModelsAllProteins
#' not inteded to be run on its own, you should copy/paste and modify as needed in your own notebook/script
#' 
#' @param protQuant.txt here this is the output of MSstats::dataProcess()$RunLevelData
#' 
.linearModelUsageExample <- function (protQuant.txt = 
                                      "~/UCSF/kroganlab/BenPolacco/darpa/2021_01_25_Qiongyu_MOR01_PlasmaMembrane_data/2021_01_28_Full_DataProcess_RunlevelData.txt.gz"){
  
  protQuant <- fread (protQuant.txt)
  
  # define the columns you need based on the condition name stored in GROUP_ORIGINAL
  protQuant[, ptx := grepl("PTX", GROUP_ORIGINAL)]
  # this uses times as strings  -- no ordering information is used. Best when time points is < 5
  protQuant[, time := gsub ("D([0-9]+)(PTX)?", "\\1", GROUP_ORIGINAL)]
  protQuant[time %in% c("CTRL", "PTX"), time := "00"]
  
  formulas = list (full = LogIntensities~as.factor(time)*ptx + SUBJECT_ORIGINAL,
                   fixed.ptx = LogIntensities~as.factor(time)+ptx + SUBJECT_ORIGINAL,
                   no.subject = LogIntensities~as.factor(time)*ptx )
  
  
  
  res <- linearModelsAllProteins(protQuant, formulas, "Protein")
  
  fwrite (res$anova, "All_lm_ANOVA_results.csv")
  fwrite (res$coef,  "All_lm_coef_results.csv")
  
}
