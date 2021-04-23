# simulation study 1

#  rm(list=ls())

  options(scipen = 99)

# load packages and functions ----------

  library(rPref)
  library(ggplot2)
  
  source("calib_true.R")
  source("calibWrapper.R")
  source("utils.R")


# SETUP ----------------------
  set.seed(1)
  n_target_sets <- 5000
  n_participants <- round(runif(n_target_sets, min = 50, max = 1000))

  
# Simulate true markov and targets -----------
  true_results <- mySimulationStudyWrapper(
    n_studies = n_target_sets, 
    n_participants = n_participants
  )

  # saveRDS(true_results, "true_results2.RDS")
  # true_results <- readRDS("true_results2.RDS")
  
  # true icer and nmb
  true_markov <- true_results$true_markov


  # sets of simulated targets
  target_sets <- true_results$target_sets


# CALIBRATION LOOP -------------------------------------------------

  # simulate calibration per target set
  sim_res <- data.frame()
  for(s in 1:n_target_sets){
    cat("\r  run:  ",s, "  of ", n_target_sets)
    runs_i <- runif(1, min = 1000, max = 50000)
    res_i <- calibWrapper(targets = target_sets[,s], RUNS = runs_i, thresh = 100000)
    res_i$tset_index <- s
    res_i$runs <- runs_i
    sim_res <- rbind(sim_res, res_i)
  }

# # SAVEGAME -----------------
  # saveRDS(sim_res,"save_sim_res2.RDS")
  # saveRDS(true_results,"save_true_results2.RDS")

  # sim_res <- readRDS("save_sim_res2.RDS")
  # true_results <- readRDS("save_true_results2.RDS")

# CHECK RESULTS ----------------------------------------------------
  res_table <- resTabler1(sim_res, true_markov, boot_iter =  1000)
  res_table
  
  diff_table <- resTabler2(sim_res, true_markov, boot_iter =  1000)
  diff_table
  
  by(sim_res$nmb_i,sim_res$methods, sd )
  by(sim_res$nmb_i,sim_res$methods, IQR )

  mean_df <- aggregate(cbind(c_i,q_i)~methods, sim_res, mean)
  
  ggplot() +
    geom_point(data = sim_res, aes(x = c_i, y= q_i, col = methods), size = 0.5, alpha = .5) +
    geom_point(data = mean_df, aes(x = c_i, y= q_i, col = methods)) +
    theme_minimal()+
    ylab("incremental QALYS") +
    xlab("incremental costs") +
    coord_cartesian(xlim = c(40000, 120000), ylim = c(0,1))
    
  
  
# # # ANY ASSOCIATIONS WITH TARGETS OR SETUP VARS? ------------
#   sim_res$nmb_diff <- sim_res$nmb_i - true_markov["nmb"]
#   
#   sim_res$t1 <- target_sets[1, sim_res$tset_index]
#   sim_res$t2 <- target_sets[2, sim_res$tset_index]
#   sim_res$t3 <- target_sets[3, sim_res$tset_index]
#   sim_res$t4 <- target_sets[4, sim_res$tset_index]
#   
#   sim_res$n_participants <- true_results$n_participants[sim_res$tset_index]
#   
#   ggplot(sim_res) +
#     geom_point(aes(x = n_participants, y = nmb_diff, col = methods), size = 0.2, alpha = 0.5) +
#     geom_smooth(aes(x = n_participants, y = nmb_diff, col = methods, fill = methods), se = T) +
#     theme_minimal()
#   
#   ggplot(sim_res) +
#     geom_point(aes(x = runs, y = nmb_diff, col = methods), size = 0.2, alpha = 0.5) +
#     geom_smooth(aes(x = runs, y = nmb_diff, col = methods), se = F) +
#     # geom_smooth(aes(x = runs, y = nmb_diff, col = methods, fill =methods), se = T) +
#     theme_minimal()
#   
#   ggplot(sim_res[sim_res$t1 > 0.94 & sim_res$t1 < 0.99,]) +
#     geom_point(aes(x = t1, y = nmb_diff, col = methods), size = 0.5, alpha = 0.5) +
#     geom_smooth(aes(x = t1, y = nmb_diff, col = methods), se = F) +
#     theme_minimal()
#     
#   ggplot(sim_res[sim_res$t2 > 0.87 & sim_res$t2 < 0.99,]) +
#     geom_point(aes(x = t2, y = nmb_diff, col = methods), size = 0.5, alpha = 0.5) +
#     geom_smooth(aes(x = t2, y = nmb_diff, col = methods), se = F) +
#     theme_minimal()
#   
#   ggplot(sim_res) +
#     geom_point(aes(x = t3, y = nmb_diff, col = methods), size = 0.5, alpha = 0.5) +
#     geom_smooth(aes(x = t3, y = nmb_diff, col = methods), se = F) +
#     theme_minimal()
#   
#   ggplot(sim_res[sim_res$t4 > 3 & sim_res$t4 < 15,]) +
#     geom_point(aes(x = t4, y = nmb_diff, col = methods), size = 0.5, alpha = 0.5) +
#     geom_smooth(aes(x = t4, y = nmb_diff, col = methods), se = F) +
#     theme_minimal()
# 
# 
