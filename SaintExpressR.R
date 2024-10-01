library (nloptr)
library (data.table)

# ### Main Function #########
#  Code to run SaintExpress using spectral counts
#  This has been coded to mimic the C code.
#  Main function is SaintExpressR.SPC(interactions, baits, preys)


#' SaintExpressR.SPC
#'  If run with only those three parameters, it should run (nearly) identically to SaintExpress C code fpr spectral counts
#'  interactors, baits, preys can be data.table or character paths to files that match the SaintExpress C code input
#'  optional parameters change from default behavior
#' @param interactors a character path to a file, or a data.table, as input to SaintExpress. Format requirements are
#'                    similar ot SaintExpress requires, though this may be more tolerant.  `baits` and `preys` (below)
#'                    are similarly the usual input to SaintExpress
#' @param baits a character path to a file, or a data.table, as input to SaintExpress
#' @param preys a character path to a file, or a data.table, as input to SaintExpress
#'
#' @param pseudoControls a data.table with columns "run", "prey", "mean", "omega" that define the background
#'                       distribution per run/prey.  Missing prey will use the experimental controls like typical saintExpress
#'                       omega is an over-dispersion parameter, and is 1 - sqrt(mean/variance)
#'                       Over-dispersion occurs in a Poisson-like distribution when variance is greater than mean
#' @param interactorFCthreshold SAINTExpress has this hard-coded as 5. A true interactor is assumed
#'                              to have a probability distribution with mean 5 x the control mean. If
#'                              observed spectral count for a bait are greater than 5x, the likelihood at 5x is
#'                              used for true-case so as to not penalize when spc are greater than 5x.
#' @param fixedBeta1 By default, SAINTExpress estimates a prior beta1 to optimize the overall likelihood. This prior
#'                   appears in some ways to be a log-odds of true:false in the whole dataset, though I haven't worked
#'                   out the math to confirm this. Instead of optimizing, you can set this to a fixed value.
#'                   log(1/3)  is approximately -1.10
#'                   log(1/9)  is approximately -2.20
#'                   log(1/19) is approximately -2.95
#'                   log(1/99) is approximately -4.60
SaintExpressR.SPC <- function(interactions, baits, preys, pseudoControls = NULL, interactorFCthreshold = 5, fixedBeta1 = NULL){
  # arguments can be data.frames or paths 
  .processInput <- function(input){
    if (!"data.frame" %in% class (input))
      input <- fread (input, header = FALSE)
    else setDT(input)
    return (input)
  }
  interactions <- .processInput(interactions)
  baits <- .processInput(baits)
  preys <- .processInput(preys)
  
  #if(!is.null(pseudoControls)) pseudoControls <- .processInput(pseudoControls)
  
  # set column names, only if all names are missing
  .colNames <- function (dt, cnames){
    if (all(!cnames %in% colnames(dt)))
      setnames(dt, cnames)
  }
  .colNames(baits, c("run", "bait", "treatment"))
  .colNames(interactions, c("run", "bait", "prey", "spc"))
  .colNames(preys, c("prey", "preyLength", "preyGene" ))
  
  
  # merge/expand to get all info in one table
  # expand interactions table so we can insert zeros
  # every run must have an entry for every prey
  # do a full cartesian join on a dummy column
  fullRows <- data.table(run = unique(interactions$run), dummy = 1)[
    data.table(prey = unique(interactions$prey), dummy = 1),
    ,on = "dummy", allow.cartesian = TRUE][
      , dummy := NULL][]
  
  interactions <- interactions[fullRows, , on = c("run", "prey")]
  interactions[is.na(spc), spc := 0]
  
  # merge the inputs into a single table, copying only the columns we need
  interactions[preys, preyGene := i.preyGene, on = "prey"]
  interactions[baits, c("bait", "treatment") := .(i.bait, i.treatment), on = "run"]
  
  # zero out the bait-bait interactions
  interactions[bait == preyGene, spc := 0]
  
  
  # ########## controls - basic stats ######
  preyStats.controls <- interactions[treatment == "C", .(mean= mean(spc), variance = var(spc), nreps= .N, spc = paste0(spc, collapse = "|")), by = prey]
  preyStats.controls[mean < 0.1, mean := 0.1] # when no or very low spectral counts detected in control, set a low non-zero value of 0.1
  preyStats.controls[, omega := ifelse(variance > mean, 1 - sqrt(mean/variance), 0.1)]
  # my addition:  otherwise omega can be less than 0.1 only if variance is slightly higher than mean.  (Meanwhile underdispersion gets 0.1...)
  preyStats.controls[omega < 0.1, omega := 0.1]
  
  if (is.null(pseudoControls)){
    message ("Using experimental controls as prey-specific background for all runs")
  } else{
    # pseudoControls holds the necessary info
    # requires columns "run", "prey", "mean", "omega"
    stopifnot(all(c("run", "prey", "mean", "omega") %in% colnames(pseudoControls)))
    message ("Using stats in pseudoControls as run-and-prey-specific background")
  }
  
  

  # ############# treatments - basic stats ###########
  if (is.null(pseudoControls)){
    # traditional way is to have a single mean/omega per prey
    preyStats.treatments <- interactions[treatment == "T"][
      preyStats.controls,
      c("controlMean", "controlOmega", "controlSPC") := .(i.mean, i.omega, i.spc),
      on = "prey"][]
  } else{
    # pseudoControl is specific to run and prey:
    preyStats.treatments <- interactions[treatment == "T"][
      pseudoControls,
      c("controlMean", "controlOmega", "controlSPC", "controlSource") := .(i.mean, i.omega, "viaPseudoControls", "pseudoControls"),
      on = c("run", "prey")][]
    
    
    if (any(is.na(preyStats.treatments$controlMean))){
      message (sprintf("%d of %d total preys in treatments are not found in pseudo-controls, using experimental controls for these",
                     length(preyStats.treatments[is.na(controlMean), unique(prey)]),
                     length(unique(preyStats.treatments$prey))))
      
      preyStats.treatments[preyStats.controls,
                           c("controlMean", "controlOmega", "controlSPC", "controlSource") :=
                             .(ifelse(is.na(controlMean),i.mean,controlMean),
                               ifelse(is.na(controlOmega), i.omega, controlOmega),
                               i.spc,
                               ifelse(is.na(controlSource), "experimentalControls", controlSource )
                               
                               #ifelse(is.na(controlSPC),    i.spc, controlSPC)
                               ),
                           on = "prey"]
      
    }
  }
  
  # when no or very low spectral counts detected in pseudo-control, set a low non-zero value of 0.1
  preyStats.treatments[controlMean < 0.1, controlMean := 0.1]
  
  preyStats.treatments[, trueOmega := 0.1] # this is hardcoded in SAINT
  
  
  preyStats.treatments[, trueMean := interactorFCthreshold * controlMean] # this factor of 5 is hardcoded in SAINT
  preyStats.treatments[, foldChange := spc/controlMean]
  preyStats.treatments[, Z := ifelse(foldChange > 2, 1,0) ] # a first guess at whether this is a true interactor or not; will update during beta1 optimization
  # Z is used by optimizeBeta1
  
  preyStats.treatments[, log.p.true := log_dgpois(ifelse(spc > trueMean, trueMean, spc), trueMean, trueOmega)] # the ~likelihood the spc came from a distirbution 5 * the control
  preyStats.treatments[, log.p.false := log_dgpois(ifelse(spc < controlMean, controlMean, spc), controlMean, controlOmega)] # the likelihood the spc came from the control distribution.
  
  
  
  # ######  optimizing beta1 and Z guesses   #############
  #  optimizing beta1 and Z guesses
  ################### #
  # Iterated Conditional Mode
  # iterative procedure here. Assign true/false guess based on likelihood at current beta1.
  # then optimize beta1 given current true/false guess
  # based on function icms in main.cpp
  if (is.null(fixedBeta1)){
    message ("Optimizing beta1 to maximize overall likelihood, this may take a minute...")
    beta1 <- 0.0
    newllik = dp.llikelihood(beta1, preyStats.treatments)
    for (i in 1:15){
      oldllik = newllik;
      #print (c(i, beta1, oldllik, newllik))
      updateZ(preyStats.treatments, beta1)
      beta1 <- optimizeBeta1(preyStats.treatments, beta1 = beta1)
      newllik = dp.llikelihood(beta1, preyStats.treatments)
      print (c(i = i, beta1 = beta1, oldllik = oldllik, newllik = newllik))
      if( (newllik >= oldllik) &
          (exp(newllik - oldllik) -1 < 1e-3) ) break
    }
  }else{
    message ("Fixing beta1 prior to ", fixedBeta1)
    beta1 <- fixedBeta1
    updateZ(preyStats.treatments, beta1)
  }
  
  message ("exp(beta1) prior = ", signif(exp(beta1), 2))

  
  # ###### saint score per prey in each run ##########
  # beta1 from above-optimized procedure:
  .unnorm_score_true <- mrf_true <- exp(beta1 + 0 * 0) # gamma * gsum, gamma = 0, gsum = 0
  preyStats.treatments[, unnorm_score_false := exp(log.p.false - log.p.true)] # this is likelihood ratio, false/true
  preyStats.treatments[, saintScore := .unnorm_score_true/(.unnorm_score_true + unnorm_score_false)]
  
  # ##### workup to final scores per prey/run ##########
  # when counts are lower than control or at 1 or zero, hard code the saintScore to 0.0
  preyStats.treatments[spc <=1 | spc < controlMean, saintScore := 0.0]
  preyStats.treatments[, avgScore := mean(saintScore), by = .(bait, prey)]
  
  # ##### summarize across reps/BFDR ###########
  bfdr.dt <- preyStats.treatments[, .(avgSpec = mean(spc),
                                      treatCounts = paste0(spc, collapse = "|"),
                                      avgScore = mean(saintScore),
                                      foldChange = mean(foldChange),
                                      controlMean = paste0(sprintf("%.1f", controlMean), collapse = "|"),
                                      controlOmega = unique(controlOmega),
                                      controlSource = unique(controlSource)
                                      ),
                                  by = .(bait, prey, preyGene)]
  
  # inefficient, but an attempt to match the c code which loops over whole table per each row...
  # now replaced with code below for a more efficient implementation.
  # bfdr.dt[,  c("numer", "denom") := .(sum((bfdr.dt$avgScore > avgScore) * bfdr.dt$avgScore),
  #                                     sum(bfdr.dt$avgScore > avgScore)),
  #         by = .I]
  
  setorder(bfdr.dt, -avgScore)
  # denom is basically a rank (how many score better than current score), starting at 0
  bfdr.dt[,denom := .N - frank(avgScore, ties.method = "max")]
  # numer is similar to a cumulative sum of avgScore, counting only those with lower rank.
  # to get only those with lower rank we jump through some hoops, and create a temporary table:
  numers <- bfdr.dt[, .(denom, avgScore, numer = cumsum(avgScore)) # create a table of cumsums
                    ][, .(numer = max(numer)), by = .(denom, avgScore) # summarize to max per rank aka denom. These should be applied to the next highest rank...
                      ][, numer := c(0.0,numer)[1:.N]] # ..so we offset the column by 1 position, inserting zero in first position
  bfdr.dt[numers, numer := i.numer , on = "denom"]

  bfdr.dt[, bfdr := ifelse(denom == 0, 0, 1-numer/denom)]
  
  
  
  return (list(controls = preyStats.controls, treatments = preyStats.treatments, bfdr = bfdr.dt, beta1 = beta1))
}

# Helper functions ####

# Generalized Poisson

# modeled loosely on the 
log_dgpois <- function(x, expect, omega){
  l1 <- expect * (1-omega)
  #if (x == 0) return(-l1)
  tmp <- l1 + x * omega
  return (log(l1) + (x-1) * log(tmp) - tmp - lfactorial(x))
}

# Model_data::llikelihood
# 
# why does htis have the exp(log.p.[trueOrfalse]) factor in the product while that is missing from the .llik_MRF_gamma_0 function?
# this is a copy of the SAINTexpress code.
dp.llikelihood <- function(beta1, treatments.dt){
  loglik = 0
  MRFtrue = exp(beta1 )
  MRFfalse = exp(0) #beta0
  sub.dt <- treatments.dt[, .(allZero = all(spc == 0),
                              m = .N,
                              paste0(Z, collapse = ","),
                              product = prod(ifelse(Z == 1, MRFtrue, MRFfalse)/(MRFtrue + MRFfalse) * exp(ifelse(Z == 1, log.p.true, log.p.false)))),
                          by = .(bait, prey)]
  loglik <- sub.dt[allZero != TRUE, sum(log(product)/m)]
  return(loglik)
}

optimizeBeta1 <- function (treatments.dt, beta1 = 0.0, lb = -15, ub = 15, print_level = 0){
  # the function to optimize, relies on local scoped treatments.dt
  .llik_MRF_gamma_0 <- function (x){
    beta1_ = x
    loglik = 0
    MRFtrue = exp(beta1_ )
    MRFfalse = exp(0) #beta0
    sub.dt <- treatments.dt[, .(allZero = all(spc == 0),
                                m = .N,
                                paste0(Z, collapse = ","),
                                product = prod(ifelse(Z == 1, MRFtrue, MRFfalse)/(MRFtrue + MRFfalse))),
                            by = .(bait, prey)]
    loglik <- sub.dt[allZero != TRUE, sum(log(product)/m)]
    return(loglik)
  }
  # nlopt R can only minimize but we want to maximize, so take negative
  .neg_llik_MRF_gamma_0 <- function(x){
    -.llik_MRF_gamma_0(x)
  }
  
  # options to match SaintExpress Stats.cpp::Model_data::wrt_MRF_gamma_0
  opts <- list(
    "algorithm"= "NLOPT_LN_COBYLA",
    #"ftol_abs"= 1.0e-4,  # takes way too long with very marginal gains
    "ftol_abs" = 1.0e-2,
    "maxeval"= 1e4,
    "print_level" = print_level
  )
  
  if (print_level > 0)
    cat ("Optimizing beta1 on treatment data")
  res <- nloptr (
    x0 = beta1,
    eval_f = .neg_llik_MRF_gamma_0,  # the R nloptr interface only minimizes, we want to maximize, so use the negative function
    lb = lb,
    ub = ub,
    #eval_g_eq = eval_g_eq,
    opts = opts
  )
  if (print_level > 0)
    cat (res$message)
  
  return(res$solution)
  
}


# icm_Z
updateZ <- function (treatments.dt, beta1){
  beta0 = 0
  treatments.dt[, loglikelihood_Ztrue := beta1 + log.p.true]
  treatments.dt[, loglikelihood_Zfalse := beta0 + log.p.false]
  treatments.dt[, Z := ifelse(loglikelihood_Ztrue > loglikelihood_Zfalse, 1, 0)]
  
}


# Main functions above
## #### test and setup functions below #####

.SaintExpressR.setupData <- function(){
  return(c("This function simply wraps some code I used to setup the data files. Not intended to be used as a function.",
           "This is here mainly as a record"))
  #psc <- fread ("../bp_utils/data/example/SAINTexpress/PC.Batch14.interactions.not.rounded.txt")
  psc <- fread ("/Users/ben/UCSF/kroganlab/SummerInterns2024/Miles/SAINTWithPseudoControls/PC.interactions.txt")
  # we need columns  "run", "prey", "mean", "omega"
  
  setnames(psc, new = c("run", "bait", "prey", "spc"))
  
  actualData <- psc[!grepl("-fake$", run)]
  pseudoControls   <- psc[grepl("-fake$", run)]
  
  
  # for now use actual controls as estimate of variance/omega
  controls <- actualData[grep ("EGFP|EmptyVector", bait)]

  # fill in zeros in spc:
  allControlRuns <- unique(controls[, .(run, bait)])
  allPreys <- unique(psc$prey)
  allByAll <- allControlRuns[, .(bait = bait, prey = allPreys), by = run]
  controls <- controls[allByAll, , on = c("run", "bait", "prey")]
  controls[is.na(spc), spc := 0]
  
  # mean, variance, omega
  control.stats <- controls[, .(m = mean(spc), v = var(spc), omega = 0.1), by = prey]
  control.stats[v > m, omega := 1 - sqrt(m/v)]
  control.stats[omega < 0.1, omega := 0.1]
  pseudoControls[control.stats, omega := i.omega, on = "prey"]
  setnames(pseudoControls, old= "spc", new = "mean")
  
  pseudoControls[, run := gsub ("-fake$", "", run)]
  pseudoControls[, bait := gsub ("-fake$", "", bait)]
  

  fwrite (pseudoControls, "../bp_utils/data/example/SAINTexpress/pseudoControls.txt", sep = "\t")
  

}


testSAINTR <- function(){
  interactions <- "../bp_utils/data/example/SAINTexpress/interactions.txt"
  baits <- "../bp_utils/data/example/SAINTexpress/baits.txt"
  preys <- "../bp_utils/data/example/SAINTexpress/preys.txt"
  pseudoControls <- fread ("../bp_utils/data/example/SAINTexpress/pseudoControls.txt")

  d <- SaintExpressR.SPC(interactions, baits, preys,
                    pseudoControls = pseudoControls,
                    # fixedBeta1 = -2.95,
                    # interactorFCthreshold = 3
                    )
    
}




