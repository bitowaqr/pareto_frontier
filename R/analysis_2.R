


res_mat <- data.frame(readRDS("res_mat.RDS"))
for(i in c(1,3:ncol(res_mat))){
  res_mat[,i] <- as.numeric(res_mat[,i])
}

names(res_mat) <- c(
  "run","method","inmb","inmb_error","inmb_error_abs","pareto_diff",
  "pS1_S2_diff", "RRS1_D_diff", "RRS2_D_diff","run2",
  "n_calib","pareto_n","n_study_size","true_inmb"
  )
res_mat$target_set <- rep(c(1,2,3,4), each = 3) # forgot to specify this previously
dim(res_mat)
head(res_mat)



# Pareto Frontier Calibration Method Evaluation

# 1. Reference case: n_calib = 10,000 
ref_df <- subset(res_mat, n_calib == 50000)



mean(ref_df$true_inmb); sd(ref_df$true_inmb)

# # Absolute error
# agg_inmb_error <- aggregate(inmb_error ~ method + target_set, ref_df,
#   function(x){
#     y <- c("mean" = mean(x), "sd" = sd(x), quantile(x, c(0.5, 0.25, 0.75, 0, 1)))
#     y <- formatC(y, digits = 0, format = "f", big.mark = ",")
#     y <- c(
#       "mean (SD)" = paste0(y[1], " (",y[2],")"),
#       "median (IQR)" = paste0(y[3], " (",y[4],"; ",y[5],")")#,
#       # "range" = paste0("[",y[6], "; ",y[7],"]")
#       )
#     y
#   })
# 
# agg_inmb_error <- reshape(agg_inmb_error, direction = "wide", timevar = "target_set", idvar = "method")
# colnames(agg_inmb_error) <- c("Method" , paste("Target",1:4))
# agg_inmb_error


# INMB ABS ERROR -----
agg_inmb_error_abs <- aggregate(inmb_error_abs ~ method + target_set, ref_df,
                            function(x){
                              y <- c("mean" = mean(x), "sd" = sd(x), quantile(x, c(0.5, 0.25, 0.75, 0, 1)))
                              y <- formatC(y, digits = 0, format = "f", big.mark = ",")
                              y <- c(
                                "mean (SD)" = paste0(y[1], " (",y[2],")"),
                                "median (IQR)" = paste0(y[3], " (",y[4],"; ",y[5],")")#,
                                # "range" = paste0("[",y[6], "; ",y[7],"]")
                              )
                              y
                            })

agg_inmb_error_abs <- reshape(agg_inmb_error_abs, direction = "wide", timevar = "target_set", idvar = "method")
colnames(agg_inmb_error_abs) <- c("Method" , paste("Target",1:4))
agg_inmb_error_abs

inmb_error_abs_means <- aggregate(inmb_error_abs ~method+target_set,ref_df,mean )
ggplot(ref_df) +
  geom_density(aes(x = inmb_error_abs,col = method,fill = method), alpha = 0.5, position = "identity") +
  geom_vline(data = inmb_error_abs_means, aes(xintercept = inmb_error_abs, col = method)) +
  facet_wrap(~target_set) +
  coord_cartesian(xlim = c(0, 20000)) +
  xlab("Absolute error in Incremental Net Monetary Benefit") +
  theme_minimal()


# Parameters -------
agg_params <- aggregate(cbind(pS1_S2_diff,RRS1_D_diff,RRS2_D_diff) ~ method + target_set, ref_df,
                                function(x){
                                  x <- abs(x)
                                  y <- c("mean" = mean(x), "sd" = sd(x), quantile(x, c(0.5, 0.25, 0.75, 0, 1)))
                                  y <- formatC(y, digits = 3, format = "f", big.mark = ",")
                                  y <- c(
                                    "mean (SD)" = paste0(y[1], " (",y[2],")")# ,
                                    # "median (IQR)" = paste0(y[3], " (",y[4],"; ",y[5],")")#,
                                    # "range" = paste0("[",y[6], "; ",y[7],"]")
                                  )
                                  y
                                })
agg_params


# calib n effect

# INMB ABS ERROR -----
agg_inmb_error_abs_all <- aggregate(inmb_error_abs ~ method  + n_calib +target_set, res_mat,
                                function(x){
                                  y <- c("mean" = mean(x), "sd" = sd(x), quantile(x, c(0.5, 0.25, 0.75, 0, 1)))
                                  y <- formatC(y, digits = 0, format = "f", big.mark = ",")
                                  y <- c(
                                    "mean (SD)" = paste0(y[1], " (",y[2],")")#,
                                    # "median (IQR)" = paste0(y[3], " (",y[4],"; ",y[5],")")#,
                                    # "range" = paste0("[",y[6], "; ",y[7],"]")
                                  )
                                  y
                                })

agg_inmb_error_abs_all

ggplot(res_mat[1:50000,]) +
  geom_boxplot(aes(y=inmb_error_abs, x=as.factor(n_calib), fill = method)) +
  facet_wrap(~target_set) +
  theme_minimal() 


# pareto n size -----

pareto_n_df <- subset(res_mat, method == "pareto")

pareto_n_df_agg <- aggregate(pareto_n ~ n_calib + target_set, pareto_n_df, 
          function(x){
            y <- c("mean" = mean(x), "sd" = sd(x), quantile(x, c(0.5, 0.25, 0.75, 0, 1)))
            y <- formatC(y, digits = 0, format = "f", big.mark = ",")
            y <- c(
              "mean (SD)" = paste0(y[1], " (",y[2],")"),
              # "median (IQR)" = paste0(y[3], " (",y[4],"; ",y[5],")")#,
              "range" = paste0("[",y[6], "; ",y[7],"]")
            )
            y})


pareto_n_df_agg <- reshape(pareto_n_df_agg, direction = "wide", timevar = "target_set", idvar = "n_calib")
colnames(pareto_n_df_agg) <- c("n_calib" , paste("Target",1:4))
pareto_n_df_agg