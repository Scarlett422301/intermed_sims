#! /usr/bin/env Rscript

# install.packages("~/scratch/intermed", type = "source", repos = NULL)

# get environment variables
MYSCRATCH <- Sys.getenv('MYSCRATCH')
RESULTDIR <- Sys.getenv('RESULTDIR')
STEPSIZE <- as.numeric(Sys.getenv('STEPSIZE'))
TASKID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# set defaults if nothing comes from environment variables
MYSCRATCH[is.na(MYSCRATCH)] <- '.'
RESULTDIR[is.na(RESULTDIR)] <- '.'
STEPSIZE[is.na(STEPSIZE)] <- 1
TASKID[is.na(TASKID)] <- 0

# get command lines arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Not enough arguments. Please use args 'listsize', 'prepare', 'run <itemsize>' or 'merge'")
}

# simulation parameters
ns <- c(2000, 4000)
# version <- 0:2
seed <- 1:1000
parm <- expand.grid(n = ns, seed = seed)

# save directory
path = "..."

save_dir <- paste0(path, "/output")
code_dir <- paste0(path, "/code")

# load packages
library(SuperLearner)
library(intermed)

make_data <- function(n = 1e2, 
                      A_success = "plogis(-1 + 0.125 * C$C1 + 0.25*C$C2)",
                      M1_success = "plogis(-0.25 + 0.45 * C$C1 + 0.125 * A + 0.5 * C$C1 * A)",
                      M2_success = "plogis(-0.25 + 0.2 * C$C2 - 0.2 * A + 0.25 * C$C2 * A)",
                      M1M2_threshold = 5,
                      Y_success = "-1 + C$C1 - C$C2 + 0.05 * M1 + 0.8 * (M1 - 3)^2 - 0.25 * M2 - 0.8*(M2 - 3)^2 + 0.125 * A + 0.15 * M1*A - 0.15*M2*A",
                      get_truth = FALSE){
  C <- data.frame(C1 = runif(n), C2 = runif(n), C5 = rbinom(n, 1, 0.5),
                  C4 = rbinom(n, 1, 0.25), C3 = runif(n))
  if(!get_truth){
    g0 <- eval(parse(text = A_success))
    A <- rbinom(n, 1, g0)
    success_prob_M1 <- eval(parse(text = M1_success))
    success_prob_M2 <- eval(parse(text = M2_success))
    M1 <- rbinom(n, 5, success_prob_M1)
    M2 <- rbinom(n, 5, success_prob_M2)
    M1[M1 > (M1M2_threshold - 1)] <- M1M2_threshold
    M2[M2 > (M1M2_threshold - 1)] <- M1M2_threshold
    Qbar0 <- eval(parse(text = Y_success))
    Y <- Qbar0 + rnorm(n, 1)
    return(list(C = C, A = A, M1 = M1, M2 = M2, Y = Y))
  }else{
    success_prob_M1_A1 <- eval(parse(text = gsub("A", "1", M1_success)))
    success_prob_M1_A0 <- eval(parse(text = gsub("A", "0", M1_success)))
    success_prob_M2_A1 <- eval(parse(text = gsub("A", "1", M2_success)))
    success_prob_M2_A0 <- eval(parse(text = gsub("A", "0", M2_success)))
    M1_A0 <- rbinom(n, 5, success_prob_M1_A0)
    M2_A0 <- rbinom(n, 5, success_prob_M2_A0)
    M1_A1 <- rbinom(n, 5, success_prob_M1_A1)
    M2_A1 <- rbinom(n, 5, success_prob_M2_A1)
    M1_A0[M1_A0 > (M1M2_threshold - 1)] <- M1M2_threshold
    M2_A0[M2_A0 > (M1M2_threshold - 1)] <- M1M2_threshold
    M1_A1[M1_A1 > (M1M2_threshold - 1)] <- M1M2_threshold
    M2_A1[M2_A1 > (M1M2_threshold - 1)] <- M1M2_threshold
    
    # total effect
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A1", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "0", Y_success)))))
    total_effect <- mean(Qbar0_A1 - Qbar0_A0)
    # direct effect
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "0", Y_success)))))
    direct_effect <- mean(Qbar0_A1 - Qbar0_A0)
    
    # indirect effect through M1
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A0", gsub("A", "1", Y_success)))))
    indirect_effect_M1 <- mean(Qbar0_A1 - Qbar0_A0)
    
    # indirect effect through M2
    Qbar0_A1 <- eval(parse(text = gsub("M2", "M2_A1", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    Qbar0_A0 <- eval(parse(text = gsub("M2", "M2_A0", gsub("M1", "M1_A1", gsub("A", "1", Y_success)))))
    indirect_effect_M2 <- mean(Qbar0_A1 - Qbar0_A0)
    
    return(list(total = total_effect,
                direct = direct_effect,
                indirect_M1 = indirect_effect_M1, 
                indirect_M2 = indirect_effect_M2,
                covar_M1M2 = total_effect - (direct_effect + indirect_effect_M1 + indirect_effect_M2
          )))
  }
}


# get the list size #########
if (args[1] == 'listsize') {
  cat(nrow(parm))
}

# execute prepare job ##################
if (args[1] == 'prepare') {
  print(paste0('initial datasets saved to: ~/drinf/scratch/dataList ... .RData'))
}

# execute parallel job #################################################
if (args[1] == 'run') {
  if (length(args) < 2) {
    stop("Not enough arguments. 'run' needs a second argument 'id'")
  }
  id <- as.numeric(args[2])
  print(paste(Sys.time(), "arrid:" , id, "TASKID:",
              TASKID, "STEPSIZE:", STEPSIZE))
  for (i in (id+TASKID):(id+TASKID+STEPSIZE-1)) {
    print(paste(Sys.time(), "i:" , i))
    print(parm[i,])
    
    # set seed
    set.seed(parm$seed[i])

    # make data
    data <- make_data(n = parm$n[i], M1M2_threshold = 5)
    table(data$M1, data$M2)

    rslt <- intermed(Y = data$Y, C = data$C, M1 = data$M1, M2 = data$M2, A = data$A, 
                     a = 1, 
                     a_star = 0,
                     n_SL = 1,
                     SL_Qbar = c("SL.glm", "SL.earth", "SL.ranger", "SL.gam"),
                     SL_g = c("SL.glm", "SL.earth", "SL.ranger"),
                     SL_Q_M = list(M1 = c("SL.glm", "SL.earth", "SL.ranger", "SL.gam"), 
                                   M2 = c("SL.glm", "SL.earth", "SL.ranger", "SL.gam")),
                     tolg = 1e-2, 
                     targeted_se = TRUE, 
                     return_models = FALSE,
                     verbose = FALSE,
                     stratify = FALSE,
                     max_iter = 2)

    rslt_parametric <- intermed(Y = data$Y, C = data$C, 
                                M1 = data$M1, M2 = data$M2, A = data$A, 
                               a = 1, 
                               a_star = 0,
                               n_SL = 1,
                               SL_Qbar = NULL,
                               glm_Qbar = "M1*A + I(M1^2) + M2*A + I(M2^2) + C1 + C2",
                               # SL_Qbar = "SL.glm",
                               # SL_g = "SL.glm",
                               glm_g = "C1 + C2",
                               SL_g = NULL,
                               glm_Q_M = list(M1 = "C1 + A", 
                                             M2 = "C2 + A"),
                               SL_Q_M = NULL,
                               tolg = 1e-2, 
                               targeted_se = TRUE, 
                               return_models = FALSE,
                               verbose = FALSE,
                               stratify = FALSE,
                               max_iter = 0)

    # parametric g-comp estimator
    get_gcomp <- function(input_data, version = 0){
      # fit linear outcome regression
      Cs = paste0('C', 1:5, collapse  = ' + ')
      if(version == 0){
        or_fit <- glm(paste0("Y ~ A*M1 + A*M2 +", Cs), 
                      family = gaussian(), data = input_data)
        # fit linear mediator regression
        M1_fit <- glm(paste0("M1 ~ A + ", Cs), family = gaussian(), data = input_data)
        M2_fit <- glm(paste0("M2 ~ A + ", Cs), family = gaussian(), data = input_data)
        
        # obtain coefficients
        thetas <- or_fit$coefficients
        beta1s <- M1_fit$coefficients
        beta2s <- M2_fit$coefficients
        # interventional direct effect
        direct_gcomp_est <- thetas["A"] + thetas["A:M1"]*(beta1s[1] + beta1s["A"]*0 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
          thetas["A:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean)))
        # interventional indirect effect via M1
        indirectM1_gcomp_est <- (thetas["M1"] +thetas["A:M1"])*beta1s["A"]
        # interventional indirect effect via M2
        indirectM2_gcomp_est <- (thetas["M2"] + thetas["A:M2"])*beta2s["A"]
        # covarM1M2 
        covarM1M2_gcomp_est <- 0
        # total effect
        total_gcomp_est <- direct_gcomp_est + indirectM1_gcomp_est + indirectM2_gcomp_est 
      }else if(version == 1){
        or_fit <- glm(paste0("Y ~ A*M1 + A*M2 + M1*M2 +", Cs), 
              family = gaussian(), data = input_data)
        # fit linear mediator regression
        M1_fit <- glm(paste0("M1 ~ A + ", Cs), family = gaussian(), data = input_data)
        M2_fit <- glm(paste0("M2 ~ A + ", Cs), family = gaussian(), data = input_data)
        
        # obtain coefficients
        thetas <- or_fit$coefficients
        beta1s <- M1_fit$coefficients
        beta2s <- M2_fit$coefficients
        # interventional direct effect
        direct_gcomp_est <- thetas["A"] + thetas["A:M1"]*(beta1s[1] + beta1s["A"]*0 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
          thetas["A:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean)))
        # interventional indirect effect via M1
        indirectM1_gcomp_est <- (thetas["M1"] + thetas["M1:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                                   thetas["A:M1"])*beta1s["A"]
        # interventional indirect effect via M2
        indirectM2_gcomp_est <- (thetas["M2"] + thetas["M1:M2"]*(beta1s[1] + beta1s["A"]*1 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                                   thetas["A:M2"])*beta2s["A"]
        # covarM1M2 
        covarM1M2_gcomp_est <- 0
        # total effect
        total_gcomp_est <- direct_gcomp_est + indirectM1_gcomp_est + indirectM2_gcomp_est 
      }else if(version == 2){
        or_fit <- lm(paste0("Y ~ A*M1 + A*M2 + M1*M2 +", Cs), data = input_data)
        # fit linear mediator regression
        M1_fit <- lm(paste0("M1 ~ A + ", Cs), data = input_data)
        M2_fit <- lm(paste0("M2 ~ A + ", Cs, "+ M1 + A*M1"), data = input_data)
        # obtain coefficients
        thetas <- or_fit$coefficients
        beta1s <- M1_fit$coefficients
        beta2s <- M2_fit$coefficients
        # interventional direct effect
        direct_gcomp_est <- thetas["A"] + thetas["A:M1"]*(beta1s[1] + beta1s["A"]*0 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
          thetas["A:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean)))
        # interventional indirect effect via M1
        indirectM1_gcomp_est <- (thetas["M1"] + thetas["M1:M2"]*(beta2s[1] + beta2s["A"]*0 + sum(beta2s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                                   thetas["A:M1"])*beta1s["A"]
        # interventional indirect effect via M2
        indirectM2_gcomp_est <- (thetas["M2"] + thetas["M1:M2"]*(beta1s[1] + beta1s["A"]*1 + sum(beta1s[3:7]*apply(input_data[,1:5], 2, mean))) + 
                                   thetas["A:M2"])*beta2s["A"]
        # covarM1M2 
        covarM1M2_gcomp_est <- summary(M1_fit)$sigma^2 * thetas["M1:M2"] * beta2s["A:M1"]
        # total effect
        total_gcomp_est <- direct_gcomp_est + indirectM1_gcomp_est + indirectM2_gcomp_est + covarM1M2_gcomp_est
        
      }
      return(list(total_gcomp_est = total_gcomp_est,
                  direct_gcomp_est = direct_gcomp_est,
                  indirectM1_gcomp_est = indirectM1_gcomp_est, 
                  indirectM2_gcomp_est = indirectM2_gcomp_est,
                  covarM1M2_gcomp_est = covarM1M2_gcomp_est))
    }
    
    # create a data frame
    full_data <- Reduce(cbind, data)
    colnames(full_data) <- c(colnames(data$C), names(data)[-1])
    full_data = full_data[,c(1:2, 5:3, 6:9)]
    # get g-comp on original data
    gcomp_est0 <- unlist(get_gcomp(full_data, version = 0))
    gcomp_est1 <- unlist(get_gcomp(full_data, version = 1))
    gcomp_est2 <- unlist(get_gcomp(full_data, version = 2))
    
    get_gcomp_boot <- function(full_data, version = 0){
      # g-comp confidence interval(bootstrap)
      M <- 1000
      # empty vectors to hold results
      boot_gcomp_direct <- boot_gcomp_indirectM1 <-boot_gcomp_indirectM2 <- boot_gcomp_total <- boot_gcomp_covarM1M2 <- rep(NA, M)
      # set a seed to ensure reproducibility
      set.seed(parm$seed[i])
      for(m in 1:M){
        # sample id's with replacement
        boot_ids <- sample(1:parm$n[i], replace = TRUE)
        boot_data <- full_data[boot_ids,]
        # get estimates based on bootstrap data
        boot_gcomp_est <- get_gcomp(boot_data, version = version)
        boot_gcomp_direct[m] <- boot_gcomp_est$direct_gcomp_est
        boot_gcomp_indirectM1[m] <- boot_gcomp_est$indirectM1_gcomp_est
        boot_gcomp_indirectM2[m] <- boot_gcomp_est$indirectM2_gcomp_est
        boot_gcomp_covarM1M2[m] <- boot_gcomp_est$covarM1M2_gcomp_est
        boot_gcomp_total[m] <- boot_gcomp_est$total_gcomp_est
      }
      
      gcomp_cil <- c(apply(cbind(boot_gcomp_total, boot_gcomp_direct, boot_gcomp_indirectM1, boot_gcomp_indirectM2, boot_gcomp_covarM1M2), 2, quantile, p = 0.025))
      gcomp_ciu <- c(apply(cbind(boot_gcomp_total, boot_gcomp_direct, boot_gcomp_indirectM1, boot_gcomp_indirectM2, boot_gcomp_covarM1M2), 2, quantile, p = 0.975))
      return(list(cil = gcomp_cil, ciu = gcomp_ciu))
    }

    gcomp_ci0 <- get_gcomp_boot(full_data, 0)
    gcomp_cil0 <- gcomp_ci0$cil
    gcomp_ciu0 <- gcomp_ci0$ciu

    gcomp_ci1 <- get_gcomp_boot(full_data, 1)
    gcomp_cil1 <- gcomp_ci1$cil
    gcomp_ciu1 <- gcomp_ci1$ciu

    gcomp_ci2 <- get_gcomp_boot(full_data, 2)
    gcomp_cil2 <- gcomp_ci2$cil
    gcomp_ciu2 <- gcomp_ci2$ciu

    # get truth
    set.seed(1234)
    truth <- c(unlist(make_data(n = 1e6, get_truth = TRUE, M1M2_threshold = 5)))

    # get confidence intervals
    all_ci <- ci(rslt, est = c("tmle", "aiptw"))
    all_ci_parametric <- ci(rslt_parametric, est = c("tmle", "aiptw"))

    # check for truth in CIs
    in_ci <- function(truth, ci){
      truth >= min(ci) & truth <= max(ci)
    }
    tmle_cover <- aiptw_cover <- rep(NA, 5)
    for(j in seq_len(5)){
      tmle_cover[j] <- in_ci(truth[j], all_ci$tmle[j, c(2,3)])
      aiptw_cover[j] <- in_ci(truth[j], all_ci$aiptw[j, c(2,3)])
    }
    tmle_cover_parametric <- aiptw_cover_parametric <- rep(NA, 5)
    for(j in seq_len(5)){
      tmle_cover_parametric[j] <- in_ci(truth[j], all_ci_parametric$tmle[j, c(2,3)])
      aiptw_cover_parametric[j] <- in_ci(truth[j], all_ci_parametric$aiptw[j, c(2,3)])
    }
    gcomp_cover0 <- rep(NA, 5)
    gcomp_cover1 <- rep(NA, 5)
    gcomp_cover2 <- rep(NA, 5)
    for(j in seq_len(5)){
      gcomp_cover0[j] <- in_ci(truth[j], c(gcomp_cil0[j], gcomp_ciu0[j]))
      gcomp_cover1[j] <- in_ci(truth[j], c(gcomp_cil1[j], gcomp_ciu1[j]))
      gcomp_cover2[j] <- in_ci(truth[j], c(gcomp_cil2[j], gcomp_ciu2[j]))
    }

    out <- c(parm$n[i], parm$seed[i], parm$version[i],
             rslt$plugin, 
             rslt_parametric$plugin,
             t(all_ci$tmle), 
             t(all_ci$aiptw), 
             tmle_cover, 
             aiptw_cover,
             t(all_ci_parametric$tmle), 
             t(all_ci_parametric$aiptw), 
             tmle_cover_parametric, 
             aiptw_cover_parametric,   
             gcomp_est0,
             gcomp_cil0, 
             gcomp_ciu0,
             gcomp_cover0,
             gcomp_est1,
             gcomp_cil1, 
             gcomp_ciu1,
             gcomp_cover1,
             gcomp_est2,
             gcomp_cil2, 
             gcomp_ciu2,
             gcomp_cover2,
             truth) 

    names(out) <- c("n", "seed",
                    "total_plugin_est", "direct_plugin_est", "indirectM1_plugin_est", "indirectM2_plugin_est", "covarM1M2_plugin_est",
                    "total_parametric_plugin_est", "direct_parametric_plugin_est", "indirectM1_parametric_plugin_est", "indirectM2_parametric_plugin_est", "covarM1M2_parametric_plugin_est",
                    "total_tmle_est", "total_tmle_cil", "total_tmle_ciu",
                    "direct_tmle_est","direct_tmle_cil","direct_tmle_ciu",
                    "indirectM1_tmle_est","indirectM1_tmle_cil","indirectM1_tmle_ciu",
                    "indirectM2_tmle_est","indirectM2_tmle_cil","indirectM2_tmle_ciu",
                    "covarM1M2_tmle_est","covarM1M2_tmle_cil","covarM1M2_tmle_ciu",
                    "total_aiptw_est", "total_aiptw_cil", "total_aiptw_ciu",
                    "direct_aiptw_est","direct_aiptw_cil","direct_aiptw_ciu",
                    "indirectM1_aiptw_est","indirectM1_aiptw_cil","indirectM1_aiptw_ciu",
                    "indirectM2_aiptw_est","indirectM2_aiptw_cil","indirectM2_aiptw_ciu",
                    "covarM1M2_aiptw_est","covarM1M2_aiptw_cil","covarM1M2_aiptw_ciu",
                    "total_tmle_cover", "direct_tmle_cover", "indirectM1_tmle_cover", "indirectM2_tmle_cover", "covarM1M2_tmle_cover",
                    "total_aiptw_cover", "direct_aiptw_cover", "indirectM1_aiptw_cover", "indirectM2_aiptw_cover", "covarM1M2_aiptw_cover",
                    "total_parametric_tmle_est", "total_parametric_tmle_cil", "total_parametric_tmle_ciu",
                    "direct_parametric_tmle_est","direct_parametric_tmle_cil","direct_parametric_tmle_ciu",
                    "indirectM1_parametric_tmle_est","indirectM1_parametric_tmle_cil","indirectM1_parametric_tmle_ciu",
                    "indirectM2_parametric_tmle_est","indirectM2_parametric_tmle_cil","indirectM2_parametric_tmle_ciu",
                    "covarM1M2_parametric_tmle_est","covarM1M2_parametric_tmle_cil","covarM1M2_parametric_tmle_ciu",
                    "total_parametric_aiptw_est", "total_parametric_aiptw_cil", "total_parametric_aiptw_ciu",
                    "direct_parametric_aiptw_est","direct_parametric_aiptw_cil","direct_parametric_aiptw_ciu",
                    "indirectM1_parametric_aiptw_est","indirectM1_parametric_aiptw_cil","indirectM1_parametric_aiptw_ciu",
                    "indirectM2_parametric_aiptw_est","indirectM2_parametric_aiptw_cil","indirectM2_parametric_aiptw_ciu",
                    "covarM1M2_parametric_aiptw_est","covarM1M2_parametric_aiptw_cil","covarM1M2_parametric_aiptw_ciu",
                    "total_parametric_tmle_cover", "direct_parametric_tmle_cover", "indirectM1_parametric_tmle_cover", "indirectM2_parametric_tmle_cover", "covarM1M2_parametric_tmle_cover",
                    "total_parametric_aiptw_cover", "direct_parametric_aiptw_cover", "indirectM1_parametric_aiptw_cover", "indirectM2_parametric_aiptw_cover", "covarM1M2_parametric_aiptw_cover",
                    "total_gcomp_est0", "direct_gcomp_est0", "indirectM1_gcomp_est0", "indirectM2_gcomp_est0", "covarM1M2_gcomp_est0",
                    "total_gcomp_cil0", "direct_gcomp_cil0", "indirectM1_gcomp_cil0", "indirectM2_gcomp_cil0", "covarM1M2_gcomp_cil0",
                    "total_gcomp_ciu0", "direct_gcomp_ciu0", "indirectM1_gcomp_ciu0", "indirectM2_gcomp_ciu0", "covarM1M2_gcomp_ciu0",
                    "total_gcomp_cover0", "direct_gcomp_cover0", "indirectM1_gcomp_cover0", "indirectM2_gcomp_cover0", "covarM1M2_gcomp_cover0",
                    "total_gcomp_est1", "direct_gcomp_est1", "indirectM1_gcomp_est1", "indirectM2_gcomp_est1", "covarM1M2_gcomp_est1",
                    "total_gcomp_cil1", "direct_gcomp_cil1", "indirectM1_gcomp_cil1", "indirectM2_gcomp_cil1", "covarM1M2_gcomp_cil1",
                    "total_gcomp_ciu1", "direct_gcomp_ciu1", "indirectM1_gcomp_ciu1", "indirectM2_gcomp_ciu1", "covarM1M2_gcomp_ciu1",
                    "total_gcomp_cover1", "direct_gcomp_cover1", "indirectM1_gcomp_cover1", "indirectM2_gcomp_cover1", "covarM1M2_gcomp_cover1",
                    "total_gcomp_est2", "direct_gcomp_est2", "indirectM1_gcomp_est2", "indirectM2_gcomp_est2", "covarM1M2_gcomp_est2",
                    "total_gcomp_cil2", "direct_gcomp_cil2", "indirectM1_gcomp_cil2", "indirectM2_gcomp_cil2", "covarM1M2_gcomp_cil2",
                    "total_gcomp_ciu2", "direct_gcomp_ciu2", "indirectM1_gcomp_ciu2", "indirectM2_gcomp_ciu2", "covarM1M2_gcomp_ciu2",
                    "total_gcomp_cover2", "direct_gcomp_cover2", "indirectM1_gcomp_cover2", "indirectM2_gcomp_cover2", "covarM1M2_gcomp_cover2",
                    "true_total", "true_direct", "true_indirectM1", "true_indirectM2", "true_covarM1M2"
    ) 
    save(out, file = paste0(save_dir, "fit_", i, ".RData"))
  }
}


