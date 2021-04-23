# Assess Calbration Performance I: the TRUE model

mySimulationStudyWrapper <- function(n_studies = 10, n_participants = 100, tresh = 100000){

  if(length(n_participants) < n_studies){
    n_participants <- rep(n_participants, length.out = n_studies)
  }
  
  # META PARAMS
  HORIZON <- 30
  DISC_RATE <- 0.03
  
  # TRUE MODEL FUNCTIONS ----------------------------------
  trueParams <- function(){
    
    pH_S1 <- 0.15
    pH_D <- 0.005
    pH_H <- 1-(pH_S1+pH_D)
    
    pS1_S2 <- 0.04 # runif(RUNS, 0.01, 0.48)    <- TRUE VALUE 1
    pS1_H <- 0.5
    RRS1_D <-  3 # runif(RUNS, 1, 4.5)    <- TRUE VALUE 1    
    pS1_D <- RRS1_D * pH_D
    pS1_S1 <- 1-(pS1_H+pS1_D+pS1_S2)
    
    RRS2_D <-15 #runif(RUNS, 1, 20)    <- TRUE VALUE 1   
    pS2_D <- RRS2_D * pH_D
    pS2_S2 <- 1-pS2_D 
    
    params <- cbind(pH_S1, pH_D, pH_H, pS1_S2, pS1_H, RRS1_D, pS1_D, pS1_S1, RRS2_D, pS2_D, pS2_S2 )
    params <- c(params[1,])
    
    return(params)
  }
  
  # normal computeCE
  computeCE <- function(markov_trace, DISC_RATE = 0.03, thresh = 20000){
    
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
  
  # run the TRUE Markov: what are the expecations?
  runTrueMarkov <- function(params, HORIZON, DISC_RATE = 0.03, thresh = 100000){
    
      # TRANSISTION MATRIX 
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
      
    # return res
    return(ce_res)
    
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
  nextCycle <- function(n, k, lev, i_trans_mat){
    
    ran <- matrix(NA, ncol = 1, nrow = n) 
    U <- i_trans_mat
    
    for(i in 2:k) {    # start loop, from the 2nd health states
      U[i, ] <- U[i, ] + U[i - 1, ] # start summing the probabilities of the different health states per individual
    }
    if (any((U[k, ] - 1) > 1e-05))  # sum of all probs per individual - 1 should be 0 (use 1e-05 for rounding issues), else print the error statement
      stop("error in multinom: probabilities do not sum to 1")
    
    
    un <- rep(runif(n), rep(k, n))       # sample from a uniform distribution of length n*k
    ran <- lev[1 + colSums(un > U)]      # store the health state at the jth column of the U matrix
    
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
    
    for(t in 2:(HORIZON+1)){
      individuals_probs <- iTransProbs(study_trace[,t-1],trans_mat)
      next_states <- nextCycle(n = n, k = 4, lev = lev, i_trans_mat = individuals_probs)
      study_trace[,t] <- next_states
    }
    
    t1 <- 1 - (sum(study_trace[,6] == "D") / nrow(study_trace))
    t2 <- 1 - (sum(study_trace[,11] == "D") / nrow(study_trace))
    t3 <- sum(study_trace[,6] == "S1" | study_trace[,6] == "S2" ) / nrow(study_trace)
    t4 <- sum(study_trace[,6] == "S1") / sum(study_trace[,6] == "S2")
    if(t4 == "NaN"){
      t4 <- 1
    }
    
    res <- c(
      "survival_at_5" = t1, 
      "survival_at_10" = t2, 
      "prev_at_5" = t3, 
      "ratio_s1_s2" = t4
      )
    return(res)
  }
  
  # -----------------------------
  
  # draw params
  true_params <- trueParams()
  
  # run true model
  true_res <- runTrueMarkov(true_params, HORIZON, DISC_RATE, tresh)
  
  target_sets <- matrix(ncol = n_studies, nrow = 4)
  i = 1
  while(i <= n_studies){
    res_i <- simStudy(n = n_participants[i], true_params)
    # avoid infinite values!
    if(!is.infinite(res_i[4])){
      target_sets[,i] <- res_i
      i <- i + 1
    }
  }
  rownames(target_sets) <- c("survival_at_5", "survival_at_10", "prev_at_5", "ratio_s1_s2")
  
  return(list("target_sets" = target_sets, "n_participants" = n_participants, "true_markov" = true_res))
}


