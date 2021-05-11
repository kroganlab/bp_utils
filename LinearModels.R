





#' Compute linear models, one model per protein or gene, and for each formula.
#'
#' @param fullDataTable data.table.  Required columns are dictated by formulas and splitColumn
#' @param formulaList A list of formulas. They will be computed separately.
#' @param splitColumn  character indicating name of column in fullDataTable that defines the units to compute
#'                    the linear model on.  Usually "Protein" or  "gene", etc.

linearModelsAllProteins <- function (fullDataTable, formulaList, splitColumn = "Protein"){
  subsetTables <- split (fullDataTable, fullDataTable[[splitColumn]])
  
  # do al linear models, anova, etc. Result is a list of [[anova, coef]] lists
  coef.anova.list <- pbapply::pblapply (subsetTables, function(sDT){.linearModelsOneProtein(sDT, formulaList)})

  # get all anova in a  single table
  anovas.list <- lapply(coef.anova.list, `[[`, "anova")
  anovasTable <- rbindlist (anovas.list, idcol = splitColumn, fill = TRUE)
  setnames(anovasTable, old = c("rn", "Pr(>F)"), new = c("term", "p.value"), skip_absent = TRUE)
  
  # get all coef in a single table
  coefs.list <- lapply(coef.anova.list, `[[`, "coef")
  coefTable <- rbindlist(coefs.list, idcol = splitColumn, fill = TRUE)
  
  # get all errors and warnings in a single table
  errors.list <- lapply(coef.anova.list, `[[`, "errWarn")
  errWarnTable <- rbindlist(errors.list, idcol = splitColumn, fill = TRUE)
  errWarnTable[err == "NULL", err := NA]
  errWarnTable[warn == "NULL", warn := NA]
  
  return (list (anova = anovasTable, coef = coefTable, errWarn = errWarnTable))
} 



#' the one-protein linear model computation

.linearModelsOneProtein <- function(subsetDT, formulas ){
  lms.out <- lapply (formulas, .errorWarningCatcherFactory(lm), data = subsetDT )

  errorWarnTable <- data.table (model = names (lms.out),
                                err = as.character(lapply (lms.out, `[[`, "err")),
                                warn = as.character(lapply (lms.out, `[[`, "warn")))

  lms <- lapply (lms.out, `[[`, "value")
  lms <- lms[!sapply(lms, is.null)]
  
  if(length(lms) > 0){
    coef.list <- lapply (lms, function(lm.out)as.data.table(coefficients(summary(lm.out)), keep.rownames = TRUE))
    coef.table <- rbindlist(coef.list, idcol = "model")
    setnames (coef.table, old = c("rn", "Pr(>|t|)"), new = c("term", "p.value"), skip_absent = TRUE)
    
    #f.tests from anova function
    anova.list <- lapply (lms, function(l)as.data.table(anova(l), keep.rownames=TRUE))
    anovaTables <- rbindlist(anova.list, idcol = "model")
    
    anovaTables[, sigCode := dplyr::case_when(`Pr(>F)`< 0.001 ~ "***",
                                              `Pr(>F)`< 0.01  ~ "**",
                                              `Pr(>F)`< 0.05  ~ "*",
                                              `Pr(>F)`< 0.1  ~ ".",
                                              TRUE ~ "")]
  } else{
    anovaTables <- coef.table <-  NULL
  }
  return (list (anova = anovaTables, coef = coef.table, errWarn = errorWarnTable))
}


#' see https://stackoverflow.com/questions/4948361/how-do-i-save-warnings-and-errors-as-output-from-a-function/4952908#4952908
#' takes any function and runs it in an error/warning/message capturing environment. 
#' Returns a list of value, warnings (including messages) and errors
#' still can't get it to capture all warnings from anova ...
.errorWarningCatcherFactory <- function(fun)
  function(...) {
    warn <- err <- NULL
    res <- withCallingHandlers(
      tryCatch(fun(...), error=function(e) {
        #cat ("error handled...\n")
        err <<- conditionMessage(e)
        NULL
      }), warning=function(w) {
        warn <<- append(warn, conditionMessage(w))
        #cat ("warning handled...\n")
        invokeRestart("muffleWarning")
      },message=function(w) {
        warn <<- append(warn, conditionMessage(w))
        #cat ("message handled...\n")
        invokeRestart("muffleMessage")
      })
    list(value = res, warn=warn, err=err)
  }



#' this function holds a piece of sample code that demonstrate how to use LinearModelsAllProteins
#' not intended to be run on its own, you should copy/paste and modify as needed in your own notebook/script
#' 
#' @param protQuant.txt the input file here; in this case, the output of MSstats::dataProcess()$RunLevelData
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
