# RANDOMISTA 1

# EVALUATING PARETO FRONTIER PERFORMANCE ON TRUE TARGET SETS

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

n_true_models <- 10 # 1,0000
n_calibration_runs <- c(500, 1000, 5000, 10000, 20000, 50000)
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


# N_ROW <- 3 * n_true_models * length(n_calibration_runs) * length(set_targets_sets)
# res_mat <- matrix(data = NA, nrow = N_ROW, ncol = 14)

t1 <- Sys.time()

cl <- parallel::makeForkCluster(7)
doParallel::registerDoParallel(cl)

res_mat <- foreach(i = 1:n_true_models, .combine = 'rbind') %dopar% {
  
  res_mat_i <- matrix(
    data = NA, 
    nrow = 3 * length(n_calibration_runs) * length(set_targets_sets), 
    ncol = 14
    )
  
  true_model_i <- runTrueMarkov(params = true_params[i,], return_targets = T)
  true_inmb_i <- true_model_i$ce_res["nmb"]
  target_set_i <- true_model_i$targets

  
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
  
  for(k in seq_along(n_calibration_runs)){
    
    calib_subset_k <- 1:n_calibration_runs[k]
    
    for(j in seq_along(set_targets_sets)){
      
      cat("\r   i =",i,";j=",j,";k=",k,";xxxxxx",sep = "")
      
      selected_sets_ijk <- evalTargets(
        target_diff = calib_res$res_mat[calib_subset_k,], 
        indices = set_targets_sets[[j]]
        )    
      
      res_list_ijk <- evalRes(
        selected_sets_ijk, 
        test_params[calib_subset_k,], 
        true_inmb_i)
      
      res_ijk <- suppressWarnings(
        as.matrix(cbind(
        "run" = i,
        res_list_ijk$inmb_df,
        res_list_ijk$params_df[,-1],
        true_params_run = i,
        n_calib = n_calibration_runs[k],
        "pareto_n"  = length(selected_sets_ijk$pareto),
        "n_study_size" = NA,
        "true_inmb" = true_inmb_i
        ))
      )
      
      res_mat_i[which(is.na(res_mat_i[,1]))[1:3],] <- res_ijk
      # res_mat[which(is.na(res_mat[,1]))[1:3],] <- res_ijk
      
    }
    
  }
  
  return(res_mat_i)
}
  
t1 - Sys.time()  
parallel::stopCluster(cl)  



# dim(res_mat)

saveRDS(res_mat_i, file = "res_mat_i.RDS")