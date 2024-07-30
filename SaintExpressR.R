library (nloptr)
library (data.table)

# Model_data::llikelihood
# 
# why does htis have the exp(log.p*) factor in the product while that is missing from the .llik_MRF_gamma_0 function?
# this is a copy of the SAINTexpress code.
dp.llikelihood <- function(beta1, treatments.dt){
  loglik = 0
  MRFtrue = exp(beta1_ )
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
    "ftol_abs"= 1.0e-4,  
    "maxeval"= 1e4,
    "print_level" = print_level
  )
  
  if (print_level > 0)
    print ("Optimizing beta1 on treatment data")
  res <- nloptr (
    x0 = beta1,
    eval_f = .neg_llik_MRF_gamma_0,  # the R interface only minimizes, we want to maximize, so use the negative function
    lb = lb,
    ub = ub,
    #eval_g_eq = eval_g_eq,
    opts = opts
  )
  if (print_level > 0)
    print (res$message)
  
  return(res$solution)
  
}


# icm_Z
updateZ <- function (treatments.dt, beta1){
  beta0 = 0
  treatments.dt[, loglikelihood_Ztrue := beta1 + log.p.true]
  treatments.dt[, loglikelihood_Zfalse := beta0 + log.p.false]
  treatments.dt[, Z := ifelse(loglikelihood_Ztrue > loglikelihood_Zfalse, 1, 0)]
  
}


SaintExpressR <- function(interactions, baits, preys,   interactorFCthreshold = 5){
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
  
  
  ################
  # controls - get mean and variance and over-dispersion for the controls
  ################
  preyStats.controls <- interactions[treatment == "C", .(mean= mean(spc), variance = var(spc), nreps= .N, spc = paste0(spc, collapse = "|")), by = prey]
  preyStats.controls[mean < 0.1, mean := 0.1] # when no or very low spectral counts detected in control, set a low non-zero value of 0.1
  preyStats.controls[, omega := ifelse(variance > mean, 1 - sqrt(mean/variance), 0.1)]
  # my addition:  otherwise omega can be less than 0.1 only if variance is slightly higher than mean.  (Meanwhile underdispersion gets 0.1...)
  #preyStats.controls[omega < 0.1, omega := 0.1]

  #################
  # treatments - basic stats
  ################  
  preyStats.treatments <- interactions[treatment == "T"][preyStats.controls, c("controlMean", "controlOmega", "controlSPC") := .(i.mean, i.omega, i.spc), on = "prey"][]
  preyStats.treatments[, trueOmega := 0.1] # this is hardcoded in SAINT
  
  
  preyStats.treatments[, trueMean := interactorFCthreshold * controlMean] # this factor of 5 is hardcoded in SAINT
  preyStats.treatments[, foldChange := spc/controlMean]
  preyStats.treatments[, Z := ifelse(foldChange > 2, 1,0) ] # a first guess at whether this is a true interactor or not; will update during beta1 optimization
  # Z is used by optimizeBeta1
  
  preyStats.treatments[, log.p.true := log_dgpois(ifelse(spc > trueMean, trueMean, spc), trueMean, trueOmega)] # the ~likelihood the spc came from a distirbution 5 * the control
  preyStats.treatments[, log.p.false := log_dgpois(ifelse(spc < controlMean, controlMean, spc), controlMean, controlOmega)] # the likelihood the spc came from the control distribution.
  
  
  
  ####################
  #  optimizing beta1 and Z guesses
  ####################
  # iterative procedure here. Assign true/false guess based on likelihood at current beta1.
  # then optimize beta1 given current true/false guess
  # function icms in main.cpp
  beta1 <- 0.0
  newllik = dp.llikelihood(beta1, preyStats.treatments)
  for (i in 1:15){
    oldllik = newllik;
    print (c(i, beta1, oldllik, newllik))
    updateZ(preyStats.treatments, beta1)
    beta1 <- optimizeBeta1(preyStats.treatments, beta1 = beta1)
    newllik = dp.llikelihood(beta1, preyStats.treatments)
    print (c(i, beta1, oldllik, newllik))
    if( (newllik >= oldllik) &
        (exp(newllik - oldllik) -1 < 1e-3) ) break
  }

  
  #################
  # treatments - saint score per prey in each run
  ################
  # beta1 from above-optimized procedure:
  .unnorm_score_true <- mrf_true <- exp(beta1 + 0 * 0) # gamma * gsum, gamma = 0, gsum = 0
  preyStats.treatments[, unnorm_score_false := exp(log.p.false - log.p.true)]
  preyStats.treatments[, saintScore := .unnorm_score_true/(.unnorm_score_true + unnorm_score_false)]
  
  ################
  # treatments - workup to final scores
  ################  
  # when counts are lower than control or at 1 or zero, hard code the saintScore to 0.0
  preyStats.treatments[spc <=1 | spc < controlMean, saintScore := 0.0]
  preyStats.treatments[, avgScore := mean(saintScore), by = .(bait, prey)]
  
  #################
  # BFDR
  ################
  bfdr.dt <- unique(preyStats.treatments[, .(bait, prey, preyGene,avgScore)])
  setorder(bfdr.dt, -avgScore)
  
  # inefficient, but an attempt to match the c code which loops over whole table per each row...
  bfdr.dt[,  c("numer", "denom") := .(sum((bfdr.dt$avgScore > avgScore) * bfdr.dt$avgScore),
                                      sum(bfdr.dt$avgScore > avgScore)),
          by = .I]
  
  bfdr.dt[, bfdr := ifelse(denom == 0, 0, 1-numer/denom)]
  
  return (list(controls = preyStats.controls, treatments = preyStats.treatments, bfdr = bfdr.dt))
  
}