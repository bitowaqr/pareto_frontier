# ANALYIS 1

# rm(list=ls())
  
# LOAD PACKAGES AND FUNCTIONS
  library(rPref)
  library(reshape2)
  library(ggplot2)
  
  source("./R/functions.R")

# META PARAMS
  HORIZON <- 30
  DISC_RATE <- 0.03
  THRESH <- 1000
  
  set.seed(2021)

# Number of calibrations
  n_calib <- 10000

  
# TRUE MODEL  
  true_params <- drawParams(pS1_S2 = 0.02, RRS1_D = 3, RRS2_D = 15)
  true_model <- runTrueMarkov(params = true_params[1,], return_targets = T)
  true_inmb <- true_model$ce_res[4]

# SIMULATE TARGET SETS
  # target_set <- c(0.98, 0.84, 0.15, 8, 0, 0)   # from Enns et al
  # target_set <- true_model$targets             # True targets
  # target_set <- simStudy(n = 10000, true_params = true_params[1,])   # Estimated targets
  
  

# run calibration and estimate iNMB
  # draw random test params
  test_params <- drawParams(
    pS1_S2 = runif(n_calib, 0.01, 0.48),
    RRS1_D =  runif(n_calib, 1, 4.5),
    RRS2_D = runif(n_calib, 1, 20)
    )
  
  # run model for all random param sets and compre to targets
  calib_res <- runCalibration(params = test_params, targets = target_set, RUNS = n_calib)
  
  # evaluate target achievement and select winning param sets
  selected_sets <- evalTargets(target_diff = calib_res$res_mat)
    
  # retrieve results over selected param sets
  res_list <- evalRes(selected_sets, test_params, true_inmb)
  
  res_list$inmb_df
  
  ggplot() +
    geom_point(data = res_list$res_df, aes(x = incr_C, y = incr_Q, col = L1), alpha = 0.5, size = 1) +
    geom_point(aes(x = true_model$ce_res["incr_C"], y = true_model$ce_res["incr_Q"], col = "True"), alpha = 2) +
    xlim(c(45000,80000)) +
    ylim(c(0.35, 0.8)) +
    theme_minimal()
  
  
  
  table(res_list$res_df$L1)
  
  
  # 
  # ggplot() +
  #   geom_histogram(data = res_df, aes(inmb, fill = L1), alpha = 0.5) +
  #   geom_vline(xintercept = true_inmb) +
  #   facet_wrap(~L1) +
  #   theme_minimal()
  #   
  





  
  
  
  
  
  # dd <- replicate(4, runif(10000))
  # dd <- data.frame(dd)
  # names(dd) <- c("t1","t2","t3", "t4")
  # psel.indices(as.data.frame(dd[,1:4]), low(t1) * low(t2) * low(t3) * low(t4))
  # 