


# FUNCTIONS

# TRUE MODEL PARAMETERS ----------------------------------
drawParams <- function(pS1_S2 = 0.04, RRS1_D = 3, RRS2_D = 15) {
  pH_S1 <- 0.15
  pH_D <- 0.005
  pH_H <- 1 - (pH_S1 + pH_D)
  
  pS1_S2 <- pS1_S2 # runif(RUNS, 0.01, 0.48)    <- TRUE VALUE 1
  pS1_H <- 0.5
  RRS1_D <- RRS1_D # runif(RUNS, 1, 4.5)    <- TRUE VALUE 1
  pS1_D <- RRS1_D * pH_D
  pS1_S1 <- 1 - (pS1_H + pS1_D + pS1_S2)
  
  RRS2_D <- RRS2_D # runif(RUNS, 1, 20)    <- TRUE VALUE 1
  pS2_D <- RRS2_D * pH_D
  pS2_S2 <- 1 - pS2_D
  
  params <- cbind(pH_S1, pH_D, pH_H, pS1_S2, pS1_H, RRS1_D, pS1_D, pS1_S1, RRS2_D, pS2_D, pS2_S2)
  # params <- c(params[1,])
  
  return(params)
}

# helper function within run true markov return incremental NMB
computeCE <- function(markov_trace, DISC_RATE = 0.03, thresh = 100000){
  
  u_H <- 1.00
  u_S1 <- 0.75
  u_S1_tx <- 0.95
  u_S2 <- 0.5
  u_D <- 0
  
  c_H  <- 2000
  c_S1 <- 4000
  c_S2 <- 15000
  c_D  <- 0
  
  c_S1_tx <- 12000
  c_S2_tx <- 12000
  
  disc_rates <- 1 / (1 + DISC_RATE) ^ (0:(nrow(markov_trace)-1)) 
  
  costs_base <- markov_trace %*% c(c_H,c_S1, c_S2, c_D)
  costs_base <- sum(costs_base * disc_rates)
  
  costs_tx <- markov_trace %*% c(c_H,c_S1 + c_S1_tx, c_S2+c_S2_tx, c_D)
  costs_tx <- sum(costs_tx * disc_rates)
  
  qalys_base <- markov_trace %*% c(u_H,u_S1, u_S2, u_D)
  qalys_base <- sum(qalys_base * disc_rates)
  
  qalys_tx <- markov_trace %*% c(u_H,u_S1_tx, u_S2, u_D)
  qalys_tx <- sum(qalys_tx * disc_rates)
  
  i_Q <- mean(qalys_tx - qalys_base)
  i_C <- mean(costs_tx - costs_base)
  icer <- mean(costs_tx - costs_base)  / mean(qalys_tx - qalys_base)
  nmb <- mean( (qalys_tx - qalys_base)  * thresh  - (costs_tx - costs_base)  )
  
  return(c(i_Q,i_C,icer,nmb))
}

# returns incremental NMB
runTrueMarkov <- function(params, HORIZON = 30, DISC_RATE = 0.03, thresh = 100000, return_targets = F){
  
  # TRANSISTION MATRIX
  trans_mat <- matrix(
    nrow = 4, ncol = 4,
    # dimnames = list(c(STATE_NAMES),c(STATE_NAMES)),
    data = c(
      params["pH_H"],  params["pH_S1"], 0,  params["pH_D"],
      params["pS1_H"], params["pS1_S1"], params["pS1_S2"], params["pS1_D"],
      0,                              0, params["pS2_S2"], params["pS2_D"],
      0,                  0,          0,                1
    ),
    byrow = T
  )
  
  # INIT TRACE
  markov_trace <- matrix(nrow = HORIZON+1,ncol = 4)
  markov_trace[1,] <- c(1,0,0,0)
  
  # RUN MARKOV
  for (year in 2:(HORIZON+1)){
    markov_trace[year, ] <- markov_trace[year - 1, ] %*% trans_mat
  }
  
  # cost-effectiveness results
  ce_res <- computeCE(markov_trace = markov_trace,DISC_RATE = DISC_RATE,thresh = thresh)
  names(ce_res) <- c("incr_Q","incr_C","icer","nmb")
  
  
  # targets 
  if(return_targets){
    t1 <- 1 - markov_trace[6,4]  # survival at year 5
    t2 <- 1 - markov_trace[11,4] # survival at year 10
    t3 <- sum(markov_trace[6,2:3]) / sum(markov_trace[6,1:3]) # prevalence S1 and S2 in year 5
    # t4 <- log(markov_trace[6,2]/markov_trace[6,3])            # Ratio S1:S2 in year 5
    t4 <- markov_trace[6,2]/markov_trace[6,3]            # Ratio S1:S2 in year 5
    # additional targets
    t5 <- markov_trace[11,2] / sum(markov_trace[11,1:3])      # prevalence of S1 in year 10
    t6 <- markov_trace[11,3] / sum(markov_trace[11,1:3])      # prevalence of S2 in year 10
    res_t <- c(t1,t2,t3,t4,t5,t6)
    return(list(
      "ce_res" = ce_res,
      "targets" = res_t
    ))
    
  }
  
  # return res
  return(list(
    "ce_res" = ce_res
  ))
}

# runCalibration function: returns abslute different between estimated and observed targets 
runCalibration <- function(params, targets, RUNS, HORIZON = 30,  DISC_RATE = 0.035, thresh = 100000){
  
  #### CALIBRATION 
  res_mat <- matrix(
    nrow = RUNS, ncol = length(targets),
    dimnames = list(c(1:RUNS),paste0("t",1:length(targets)))
  )
  
  ce_mat <- matrix(
    nrow = RUNS, ncol = 4,
    dimnames = list(c(1:RUNS),c("incr_Q","incr_C","icer","nmb"))
  )
  
  for(r in 1:RUNS){
    res_r <- runTrueMarkov(params[r,],return_targets = T)
    res_mat[r,] <- abs(targets - res_r$targets)
    ce_mat[r,] <- res_r$ce_res
  }
  
  
  # # return res
  return(list(
    "res_mat" = res_mat,
    "ce_mat" = ce_mat
  ))
  
}  

# evaluates target differnces and returns sets of parameters that match criteria
evalTargets <- function(target_diff, indices = c("t1","t2","t3","t4")){
  
  strs <- paste(paste0("low(",indices[1],")"), paste0("* low(",indices[2:length(indices)],")", collapse = ""), collapse = "")
  
  
  # GOF 1: pareto
  # pareto_set <- psel.indices(as.data.frame(target_diff[,indices]), low(t1) * low(t2) * low(t3) * low(t4))
  pareto_set <- eval(parse(text = paste("psel.indices(as.data.frame(target_diff[,indices]),", strs, ")")))
  rmi <- !is.na(target_diff[pareto_set,1])
  pareto_set <- pareto_set[rmi]
  pareto_len <- length(pareto_set)
  
  # GOF 2: unweighted sum - pareto set size sized and 10%
  gof2 <- apply(target_diff[,indices],1, sum)
  sum_pareto <- which(rank(gof2) <= pareto_len)
  sum_percent <- which(rank(gof2) <= (length(gof2)/100))
  
  # GOF 3: weighted sum - pareto set size sized and 10%
  # target_diff[,"t4"] <- target_diff[,"t4"] / 100
  # gof3 <- apply(t(t(target_diff[,indicies]) * c(1,1,1,0.01)), 1, sum)
  # wsum_pareto <- which(rank(gof3) <= pareto_len)
  # wsum_percent <- which(rank(gof3) <= (length(gof3)/100))
  
  
  set_indices <- list(
    "pareto" = pareto_set,
    "sum_pareto" = sum_pareto,
    "sum_percent" = sum_percent# ,
    # "wsum_pareto" = wsum_pareto,
    # "wsum_percent" = wsum_percent
  )
  
  return(set_indices)
  
}

# simstudy helper function 1
iTransProbs <- function(states,i_trans_mat){
  
  individual_probs <- sapply(states,function(x) {
    switch(x,
           "H" = i_trans_mat[1,],
           "S1" = i_trans_mat[2,],
           "S2" = i_trans_mat[3,],
           "D" = i_trans_mat[4,])
  })
  return(individual_probs)
}

# simstudy helper function 2
nextCycle <- function(n, k, lev, i_trans_mat) {
  ran <- matrix(NA, ncol = 1, nrow = n)
  U <- i_trans_mat
  
  for (i in 2:k) { # start loop, from the 2nd health states
    U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual
  }
  if (any((U[k, ] - 1) > 1e-05)) { # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
    stop("error in multinom: probabilities do not sum to 1")
  }
  
  
  un <- rep(runif(n), rep(k, n)) # sample from a uniform distribution of length n*k
  ran <- lev[1 + colSums(un > U)] # store the health state at the jth column of the U matrix
  
  return(ran)
}

# simStudy function: simulate targets
simStudy <- function(n = 3, true_params, HORIZON = 10){
  
  
  trans_mat <- matrix(
    nrow = 4, ncol = 4,
    # dimnames = list(c(STATE_NAMES),c(STATE_NAMES)),
    data = c( 
      true_params["pH_H"],  true_params["pH_S1"],                  0,  true_params["pH_D"],
      true_params["pS1_H"], true_params["pS1_S1"], true_params["pS1_S2"], true_params["pS1_D"],
      0,                  0, true_params["pS2_S2"], true_params["pS2_D"],
      0,                  0,                  0,                1
    ),
    byrow = T
  )
  
  trans_mat_trimmed <- trans_mat[1:3,]
  k <- ncol(trans_mat)
  
  lev <- c("H","S1","S2","D")
  
  study_trace <- matrix(NA, ncol = HORIZON+1, nrow = n)
  study_trace[,1] <- "H"
  
  for (t in 2:(HORIZON + 1)) {
    individuals_probs <- iTransProbs(study_trace[, t - 1], trans_mat)
    next_states <- nextCycle(n = n, k = 4, lev = lev, i_trans_mat = individuals_probs)
    study_trace[, t] <- next_states
  }
  study_trace
  
  t1 <- 1 - (sum(study_trace[,6] == "D") / nrow(study_trace))
  t2 <- 1 - (sum(study_trace[,11] == "D") / nrow(study_trace))
  t3 <- sum(study_trace[,6] == "S1" | study_trace[,6] == "S2" ) / sum(study_trace[, 11] != "D")
  t4 <- log(sum(study_trace[, 6] == "S1") / sum(study_trace[, 6] == "S2"))
  t5 <- sum(study_trace[, 11] == "S1") / sum(study_trace[, 11] != "D")
  t6 <- sum(study_trace[, 11] == "S2") / sum(study_trace[, 11] != "D")
  
  if(t4 == "NaN"){
    t4 <- NA
  }
  
  res <- c(t1, t2, t3, t4, t5, t6)
  
  return(res)
}

# retrieve results over selected param sets
evalRes <- function(selected_sets, test_params, true_inmb){
  
  res_df <- melt(selected_sets)
  # res_df <- res_df[res_df$L1 %in% c("pareto","sum_pareto","sum_percent"),]
  res_df$incr_Q <- calib_res$ce_mat[res_df$value,"incr_Q"]
  res_df$incr_C <- calib_res$ce_mat[res_df$value,"incr_C"]
  res_df$inmb <- calib_res$ce_mat[res_df$value,"nmb"]
  res_df$inmb_error <- res_df$inmb - true_inmb
  res_df$inmb_error_abs = abs(res_df$inmb_error)
  
  res_df <- cbind(res_df, test_params[res_df$value,c("pS1_S2", "RRS1_D", "RRS2_D")])
  res_df$pS1_S2_diff <- res_df$pS1_S2 - true_params[1,"pS1_S2"] 
  res_df$RRS1_D_diff <- res_df$RRS1_D - true_params[1,"RRS1_D"] 
  res_df$RRS2_D_diff <- res_df$RRS2_D - true_params[1,"RRS2_D"] 
  
  inmb_df <- aggregate(cbind(inmb,inmb_error,inmb_error_abs) ~ L1, res_df, mean)
  inmb_df$pareto_diff <-  abs(inmb_df$inmb_error_abs) - abs(inmb_df$inmb_error_abs[inmb_df$L1 == "pareto"])
  
  # inmb_df
  
  params_df <- aggregate(cbind(pS1_S2_diff, RRS1_D_diff,RRS2_D_diff) ~ L1, res_df, mean)
  # params_df
  
  return(list(
    res_df = res_df,
    inmb_df = inmb_df,
    params_df = params_df
  ))
}



# aggregation function formatter
aggFormat <- function(x, digits = 0, show_median = F, show_range = F){
  y_raw <- c("mean" = mean(x), "sd" = sd(x), quantile(x, c(0.5, 0.25, 0.75, 0, 1)))
  y <- formatC(y_raw, digits = digits, format = "f", big.mark = ",")
  y_fmt <- c(
    "mean (SD)" = paste0(y[1], " (",y[2],")")
  )
  
  if(show_median){
    y_fmt <- c(y_fmt, "median (IQR)" = paste0(y[3], " (",y[4],"; ",y[5],")"))
  }
  
  if(show_range){
    y_fmt <- c(y_fmt, "range" = paste0("[",y[6], "; ",y[7],"]"))
  }
  y_fmt
}

# load object to file
load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}
