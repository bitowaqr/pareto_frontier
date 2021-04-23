# calibrate wrapper



calibWrapper <- function(
  targets = c("survival_at_5" = 0.98, 
              "survival_at_10" = 0.84, 
              "prev_at_5" = 0.15, 
              "ratio_s1_s2" = 8), 
  RUNS = 10000, 
  thresh = 100000){
  
  
  # META PARAMS ------
  HORIZON <- 30
  DISC_RATE <- 0.03
  weight_target_4 <- 0.001
  
  
  
  # targets ------
  
  
  
  
  # FUNCTIONS ----------------------------------------
  
  drawParams <- function(RUNS = 1){
    
    pH_S1 <- 0.15
    pH_D <- 0.005
    pH_H <- 1-(pH_S1+pH_D)
    
    pS1_S2 <- runif(RUNS, 0.01, 0.48)
    pS1_H <- 0.5
    RRS1_D <- runif(RUNS, 1, 4.5)    
    pS1_D <- RRS1_D * pH_D
    pS1_S1 <- 1-(pS1_H+pS1_D+pS1_S2)
    
    RRS2_D <- runif(RUNS, 1, 20)   
    pS2_D <- RRS2_D * pH_D
    pS2_S2 <- 1-pS2_D 
    
    params <- cbind(pH_S1, pH_D, pH_H, pS1_S2, pS1_H, RRS1_D, pS1_D, pS1_S1, RRS2_D, pS2_D, pS2_S2 )
    
    return(params)
  }
  
  
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
  
  
  runCalibration <- function(params, HORIZON, RUNS, targets, DISC_RATE){
    
    #### CALIBRATION 
    res_mat <- matrix(
      nrow = RUNS, ncol = 4,
      dimnames = list(c(1:RUNS),names(targets))
    )
    ce_mat <- matrix(
      nrow = RUNS, ncol = 4,
      dimnames = list(c(1:RUNS),c("incr_Q","incr_C","icer","nmb"))
    )
    
    for(r in 1:RUNS){
      
      # TRANSISTION MATRIX 
      trans_mat <- matrix(
        nrow = 4, ncol = 4,
        # dimnames = list(c(STATE_NAMES),c(STATE_NAMES)),
        data = c( 
          params[r,"pH_H"],  params[r,"pH_S1"],                  0,  params[r,"pH_D"],
          params[r,"pS1_H"], params[r,"pS1_S1"], params[r,"pS1_S2"], params[r,"pS1_D"],
          0,                  0, params[r,"pS2_S2"], params[r,"pS2_D"],
          0,                  0,                  0,                1
        ),
        byrow = T
      )
      
      # break if implausible transitions
      check1 <- sum(abs(rowSums(trans_mat)-1) < 0.0000001) < 4
      check2 <- min(trans_mat) < 0
      if(check2 | check1){
        # print(paste("ERROR in", r))
        res_mat[r,] <- NA
        next()
      }
      
      
      # INIT TRACE 
      markov_trace <- matrix(nrow = HORIZON+1,ncol = 4)
      markov_trace[1,] <- c(1,0,0,0)
      
      # RUN MARKOV
      for (year in 2:(HORIZON+1)){
        markov_trace[year, ] <- markov_trace[year - 1, ] %*% trans_mat
      }
      
      # cost-effectiveness results
      ce_res <- computeCE(markov_trace = markov_trace,DISC_RATE = DISC_RATE,thresh = thresh)
      
      # CHECK TARGETS 
      t1 <- abs(targets[1] - (1 - markov_trace[6,4]))
      t2 <- abs(targets[2] - (1 - markov_trace[11,4]))
      t3 <- abs(targets[3] - (markov_trace[6,2] + markov_trace[6,3]))
      t4 <- abs(targets[4] - (markov_trace[6,2] / markov_trace[6,3]))
      
      res_mat[r,] <- c(t1, t2, t3, t4)
      ce_mat[r,] <- ce_res
      
    }
    
    # return res
    return(list(
      "res_mat" = res_mat,
      "ce_mat" = ce_mat
    ))
    
  }
  
  
  # assess targets
  evalRes <- function(res_mat, target_weights = c(1,1,1,0.5)){
    
    require(rPref)
    # pareto frontier
    res_df <- data.frame(res_mat)
    pareto_set <- psel.indices(res_df, low(survival_at_5) * low(survival_at_10) * low(prev_at_5) * low(ratio_s1_s2))
    
    pareto_len <- length(pareto_set)
    
    # gof1
    gof1 <- apply(res_mat,1, sum)
    gof1_set <- which(rank(gof1) <= pareto_len)
    
    # gof2
    gof2 <- apply(res_mat,1, function(x){ sum(x * target_weights) })
    gof2_set <- which(rank(gof2) <= pareto_len)
    
    set_indices <- data.frame(
      "pareto" = pareto_set,
      "gof1_set" = gof1_set,
      "gof2_set" = gof2_set
    )
    
    rownames(set_indices) <- 1:nrow(set_indices)
    
    return(set_indices)
    
  }
  
  
  # get mean ICER and NMB
  showMeans <- function(res_mat = res_list$ce_mat, sets = accepted_sets, round.digits = NULL){
    
    res_out <- c()
    for(i in seq_along(sets)){
      c_i = mean(res_mat[sets[[i]],"incr_C"])
      q_i = mean(res_mat[sets[[i]],"incr_Q"])
      icer_i = mean(res_mat[sets[[i]],"icer"])
      nmb_i = mean(res_mat[sets[[i]],"nmb"])
      res_i = cbind(c_i, q_i, icer_i, nmb_i)
      res_out <- rbind(res_out, res_i)
    }
    res_out <- data.frame(res_out)
    
    if(!is.null(round.digits)){
      res_out <- round(res_out,round.digits)
    }
    res_out$methods = names(sets)
    
    return(res_out)
  }
  
  # RUN --------
  # draw params
  params <- drawParams(RUNS = RUNS)
  
  
  # run calibration
  res_list <- runCalibration(params, HORIZON, RUNS, targets, DISC_RATE)
  
  # remove NA's
  params <- params[!is.na(res_list$res_mat[,1]),]
  res_list$ce_mat <- res_list$ce_mat[!is.na(res_list$res_mat[,1]),]
  res_list$res_mat <- res_list$res_mat[!is.na(res_list$res_mat[,1]),]
  
  # eval target success
  accepted_sets <- evalRes(res_list$res_mat, target_weights = c(1,1,1,weight_target_4))
  
  # retrieve mean incremental costs, qalys, cer, and nmbs
  ce_res <- showMeans(res_mat = res_list$ce_mat, sets = accepted_sets, round.digits = 3)
  
  return(ce_res)
  
 
  
}


# calibWrapper(RUNS = 1000, thresh = 100000)