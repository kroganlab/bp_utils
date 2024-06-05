#

example.emmeans.contrastOfContrasts <- function (l, factorFormula = ~drug|tissue){  # how does the drug effect change per tissue
  emm <- emmeans(l, factorFormula)
  contrast1 <- pairs(emm)
  contrast2 <- pairs(contrast1, by = NULL, adjust = "none")
  return (as.data.table(contrast2))
}



#' Compute linear models, one model per protein or gene, and for each formula.
#'
#' @param fullDataTable data.table.  Required columns are dictated by formulas and splitColumn
#' @param formulaList A list of formulas. They will be computed separately.
#' @param splitColumn  character indicating name of column in fullDataTable that defines the units to compute
#'                    the linear model on.  Usually "Protein" or  "gene", etc.
#' @param emmeansFormula  a formula (or other object) passed as second argument to emmeans. It is expected to be
#'                        used int he context of contrasts so it must produce a $contrasts field in the output
#'                        examples: emmeansFormula = pairwise~GROUP_ORIGINAL or emmeansFormula = trt.vs.ctrl~GROUP_ORIGINAL or trt.vs.ctrl~timeStr|treatment
#' @param postProcessFunction  a function that receives an lm and returns a data.table. See example.emmeans.contrastOfContrasts for an example 

linearModelsAllProteins <- function (fullDataTable, formulaList, splitColumn = "Protein", cl = NULL,
                                     emmeansFormula = NULL, returnLMs = FALSE, postProcessFunction = NULL){
  subsetTables <- split (fullDataTable, fullDataTable[[splitColumn]])
  
  # do al linear models, anova, etc. Result is a list of [[anova, coef]] lists
  coef.anova.list <- pbapply::pblapply (subsetTables, function(sDT){.linearModelsOneProtein(sDT, formulaList, emmeansFormula, postProcessFunction)}, cl = cl)

  # get all anova in a  single table
  anovas.list <- lapply(coef.anova.list, `[[`, "anova")
  anovasTable <- rbindlist (anovas.list, idcol = splitColumn, fill = TRUE)
  setnames(anovasTable, old = c("rn", "Pr(>F)"), new = c("term", "p.value"), skip_absent = TRUE)
  
  # get all coef in a single table
  coefs.list <- lapply(coef.anova.list, `[[`, "coef")
  coefTable <- rbindlist(coefs.list, idcol = splitColumn, fill = TRUE)
  
  residuals.list <- lapply(coef.anova.list, `[[`, "residuals")
  residualsTable <- rbindlist(residuals.list, idcol = splitColumn,  fill = TRUE)
  
  # get all contrast in a single table, if any
  if (!is.null(emmeansFormula)){
    contrast.list <- lapply (coef.anova.list, `[[`, "contrast")
    contrastTable <- rbindlist(contrast.list, idcol = splitColumn, fill = TRUE)
  }else{
    contrastTable <- NULL
  }
    

  # get all postProcess in a single table, if any
  if (!is.null(postProcessFunction)){
    pp.list <- lapply (coef.anova.list, `[[`, "postProcess")
    postProcessTable <- rbindlist(pp.list, idcol = splitColumn, fill = TRUE)
  }else{
    postProcessTable <- NULL
  }
  
  
    
  if (returnLMs == TRUE){
    lms.list <- lapply(coef.anova.list, `[[`, "lms")
  } else{
    lms.list <- NULL
  }
  
  
  # get all errors and warnings in a single table
  errors.list <- lapply(coef.anova.list, `[[`, "errWarn")
  errWarnTable <- rbindlist(errors.list, idcol = splitColumn, fill = TRUE)
  errWarnTable[err == "NULL", err := NA]
  errWarnTable[warn == "NULL", warn := NA]
  
  return (list (anova = anovasTable, coef = coefTable, errWarn = errWarnTable, contrast = contrastTable, lms = lms.list, postProcess = postProcessTable, residuals = residualsTable))
} 



#' the one-protein linear model computation

.linearModelsOneProtein <- function(subsetDT, formulas, emmeansFormula, postProcessFunction){
  #print (subsetDT$Protein[1])
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
    
    # residuals
    residuals.list <- lapply(lms, function(l)cbind(as.data.table(l$model), data.table(residuals = residuals(l), fitted = fitted(l))))
    residuals.dt <- rbindlist(residuals.list, idcol = "model")
    
    # get all contrasts from emmeans
    contrastTable <- NULL
    if (!is.null(emmeansFormula)){
        contrasts.list <- lapply(lms, .errorWarningCatcherFactory(function(l)as.data.table(emmeans::emmeans(l, emmeansFormula)$contrasts))) |>
          lapply(`[[`, "value") # get values, discard errors/warnings
        contrasts.list <- contrasts.list[!sapply(contrasts.list, is.null)]
      
      if (length(contrasts.list) > 0){ # possible when function above fails
        contrastTable <- rbindlist(contrasts.list, idcol = "model")
        setnames(contrastTable, "p.value", "Tukey.p")
        if ("t.ratio" %in% colnames(contrastTable)){
          contrastTable[, p.t := pt(abs(t.ratio), df = df, lower.tail = FALSE) * 2]
        } else{
          contrastTable[, p.t := NA]
        }
      }
    }
    if(!is.null(postProcessFunction)){
      postProcess.out <- lapply (lms, .errorWarningCatcherFactory(postProcessFunction) )
      
      ..errorWarnTable.notUsed <- data.table (model = names (postProcess.out),
                                    err = as.character(lapply (postProcess.out, `[[`, "err")),
                                    warn = as.character(lapply (postProcess.out, `[[`, "warn")))
      
      postProcess.list <- lapply (postProcess.out, `[[`, "value")
      postProcess.list <- postProcess.list[!sapply(postProcess.list, is.null)]
      
      
      postProcess <- rbindlist(postProcess.list, idcol = "model")
    }else{
      postProcess <- NULL
    }
    
  } else{
    # 'initialize' the return values 
    anovaTables <- coef.table <- contrastTable <- postProcess <- residuals.dt <-  NULL
  }
  return (list (anova = anovaTables, coef = coef.table, errWarn = errorWarnTable, contrast = contrastTable, lms = lms, postProcess = postProcess, residuals = residuals.dt))
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


#' useful when there is a SUBJECT effect in MSstats protein quantities and you want to remove that from the scatter plot
#' @param protQuant a data.table that is the proteinLevelData from MSstats::dataProcess
#' @param formula leave it at default
#' @param splitColumn character, probably either "Protein" or "gene"
#' @param cl an integer (or other valid cl argument to pblapply) specifying number of processes to use, default NULL = no subprocesses
groupEmmeansAndResiduals <- function(protQuant, formula = LogIntensities ~ GROUP + SUBJECT, splitColumn = "gene", cl = NULL){
  
  .groupEMMeans <- function (l){
    as.data.table(emmeans::emmeans(l, c("GROUP" )))
  }
  
  # in case 
  protQuant[, SUBJECT := as.character(SUBJECT)]
  
  lm.out <- linearModelsAllProteins(protQuant, formulaList = list(base = formula), splitColumn = splitColumn,
                                    emmeansFormula = trt.vs.ctrl~GROUP, returnLMs = TRUE,
                                    postProcessFunction = .groupEMMeans, cl = cl)
  
  groupResiduals <- merge ( lm.out$postProcess, lm.out$residuals, by = c(splitColumn, "GROUP", "model"))
  groupResiduals[, groupResidual := emmean + residuals]
  setnames(groupResiduals, "emmean", "groupMean")
  
  return (groupResiduals)
}



