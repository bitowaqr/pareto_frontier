

# mean nmb and icer + 95% ci per calibration method

resTabler1 <- function(sim_res, true_markov, boot_iter = 1000){
  
  true_values = formatC(true_markov[c("nmb","icer")], format = "f", digits = 0, big.mark = ",")
  
  calib_res <- matrix(ncol = 4,nrow = 2)
  calib_res[,1] <- true_values 
  n_rows <- 1:(nrow(sim_res)/3)
  index <- 2
  name_m <- c()
  for(m in unique(sim_res$methods)){
    name_m <- c(name_m, m)
    df <- sim_res[sim_res$methods == m,]
    
    mean_nmb <- formatC(mean(df$nmb_i),format = "f",big.mark = ",", digits = 0)
    mean_icer <- formatC(mean(df$c_i) / mean(df$q_i),format = "f",big.mark = ",", digits = 0)
    
    nmb_ci95 <- icer_ci95 <- c()
    for(boot in 1:boot_iter){
      s_rows <- sample(x = n_rows,size = max(n_rows),replace = T)
      mean_nmb_i <- mean(df$nmb_i[s_rows])
      nmb_ci95 <- c(nmb_ci95, mean_nmb_i)
      mean_icer_i <- mean(df$c_i[s_rows]) / mean(df$q_i[s_rows])
      icer_ci95 <- c(icer_ci95, mean_icer_i)
    }
    nmb_ci95 <- formatC(round(quantile(nmb_ci95, c(0.025, 0.975))),format = "f",big.mark = ",", digits = 0)
    icer_ci95 <- formatC(round(quantile(icer_ci95, c(0.025, 0.975))),format = "f",big.mark = ",", digits = 0)
    
    res_df <- cbind(
      nmb = paste0(mean_nmb, " (",nmb_ci95[1],"; ", nmb_ci95[2],")"),
      icer = paste0(mean_icer, " (",icer_ci95[1],"; ", icer_ci95[2],")")
    )
    
    calib_res[,index] <- res_df
    index <- index +1
  }
  
  calib_res <- data.frame(calib_res)
  names(calib_res) <- c("true values", name_m)
  rownames(calib_res) <- c("NMB","ICER")
  return(calib_res)
}







resTabler2 <- function(sim_res, true_markov, boot_iter = 1000){
  
  true_values = true_markov[c("nmb","icer")]
  true_values_str = formatC(true_markov[c("nmb","icer")], format = "f", digits = 0, big.mark = ",")
  
  calib_res <- matrix(ncol = 4,nrow = 2)
  calib_res[,1] <- true_values_str 
  n_rows <- 1:(nrow(sim_res)/3)
  index <- 2
  name_m <- c()
  for(m in unique(sim_res$methods)){
    name_m <- c(name_m, m)
    df <- sim_res[sim_res$methods == m,]
    
    mean_nmb <- formatC(mean(df$nmb_i) - true_values[1],format = "f",flag = "+",big.mark = ",", digits = 0)
    mean_icer <- formatC(mean(df$c_i) / mean(df$q_i) - true_values[2],format = "f",flag = "+",big.mark = ",", digits = 0)
    
    nmb_ci95 <- icer_ci95 <- c()
    for(boot in 1:boot_iter){
      s_rows <- sample(x = n_rows,size = max(n_rows),replace = T)
      mean_nmb_i <- mean(df$nmb_i[s_rows]) - true_values[1]
      nmb_ci95 <- c(nmb_ci95, mean_nmb_i)
      mean_icer_i <- mean(df$c_i[s_rows]) / mean(df$q_i[s_rows]) - true_values[2]
      icer_ci95 <- c(icer_ci95, mean_icer_i)
    }
    nmb_ci95 <- formatC(round(quantile(nmb_ci95, c(0.025, 0.975))),format = "f",flag = "+",big.mark = ",", digits = 0)
    icer_ci95 <- formatC(round(quantile(icer_ci95, c(0.025, 0.975))),format = "f",flag = "+",big.mark = ",", digits = 0)
    
    res_df <- cbind(
      nmb = paste0(mean_nmb, " (",nmb_ci95[1],"; ", nmb_ci95[2],")"),
      icer = paste0(mean_icer, " (",icer_ci95[1],"; ", icer_ci95[2],")")
    )
    
    calib_res[,index] <- res_df
    index <- index +1
  }
  
  calib_res <- data.frame(calib_res)
  names(calib_res) <- c("true values", name_m)
  rownames(calib_res) <- c("NMB","ICER")
  return(calib_res)
}
