# RANDOMISTA 2

# EVALUATING PARETO FRONTIER PERFORMANCE ON UNCERTAIN TARGET SETS

# rm(list=ls())

# LOAD PACKAGES AND FUNCTIONS
library(rPref)
library(reshape2)
library(ggplot2)

# parallel
library(foreach)
library(parallel)

source("./R/functions.R")

# META PARAMS
HORIZON <- 30
DISC_RATE <- 0.03
THRESH <- 1000

set.seed(2021)

n_true_models <- 100000
n_calibration_runs <- c(10000)
set_targets_sets <- list(
  c("t1","t2","t3","t4"),
  c("t1","t5","t6"),
  c("t2","t5"),
  c("t1","t2","t3","t5","t6")
)


# SIMULATION
true_params <- drawParams(
  pS1_S2 = runif(n_true_models, 0.01, 0.48), 
  RRS1_D = runif(n_true_models, 1, 4.5), 
  RRS2_D = runif(n_true_models, 1, 20)
)

n_study_size = round(runif(n_true_models, min = 100, max = 2000))


# N_ROW <- 3 * n_true_models * length(n_calibration_runs) * length(set_targets_sets)
# res_mat <- matrix(data = NA, nrow = N_ROW, ncol = 14)

t1 <- Sys.time()

cl <- parallel::makeForkCluster(detectCores()-1)
doParallel::registerDoParallel(cl)

res_uncertain <- foreach(i = 1:n_true_models, .combine = 'rbind') %dopar% {
  
  res_mat_i <- matrix(
    data = NA, 
    nrow = 3 * length(n_calibration_runs) * length(set_targets_sets), 
    ncol = 13 # 14 - n_calib
  )
  
  true_model_i <- runTrueMarkov(params = true_params[i,], return_targets = F)
  true_inmb_i <- true_model_i$ce_res["nmb"]
  
  
  target_set_i <- simStudy(n = n_study_size[i], true_params[i,], HORIZON = 10)
  
  test_params <- drawParams(
    pS1_S2 = runif(max(n_calibration_runs), 0.01, 0.48),
    RRS1_D =  runif(max(n_calibration_runs), 1, 4.5),
    RRS2_D = runif(max(n_calibration_runs), 1, 20)
  )
  
  calib_res <- runCalibration(
    params = test_params, 
    targets = target_set_i, 
    RUNS = max(n_calibration_runs)
  )
  
    for(j in seq_along(set_targets_sets)){
      
      cat("\r   i =",i,";j=",j,"xxxxxx",sep = "")
      
      selected_sets_ij <- evalTargets(
        target_diff = calib_res$res_mat, 
        indices = set_targets_sets[[j]]
      )    
      
      res_list_ij <- evalRes(
        selected_sets_ij, 
        test_params, 
        true_inmb_i)
      
      res_ij <- suppressWarnings(
        as.matrix(cbind(
          "run" = i,
          res_list_ij$inmb_df,
          res_list_ij$params_df[,-1],
          target_set = j,
          # n_calib = n_calibration_runs,
          "pareto_n"  = length(selected_sets_ij$pareto),
          "n_study_size" = n_study_size[i],
          "true_inmb" = true_inmb_i
        ))
      )
      
      res_mat_i[which(is.na(res_mat_i[,1]))[1:3],] <- res_ij
      
    }
    
  
  
  return(res_mat_i)
}

Sys.time() - t1 
parallel::stopCluster(cl)  



# dim(res_uncertain)

saveRDS(res_uncertain, file = "res_uncertain.RDS")
