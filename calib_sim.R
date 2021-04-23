

# SETUP -------
library(rPref)
library(ggplot2)
library(plotly)


# META PARAMS ------
HORIZON <- 30
DISC_RATE <- 0.03
weight_target_4 <- 0.001


  
# targets ------
targets <- c(
  "survival_at_5" = 0.98, 
  "survival_at_10" = 0.84, 
  "prev_at_5" = 0.15, 
  "ratio_s1_s2" = 8
  )



# FUNCTIONS ----------------------------------------

source("calib_funs.R")

# RUN --------

# set iterations
RUNS = 10000

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


# 3D plot accepted params -----
repr_plot <- params[,c("pS1_S2","RRS1_D","RRS2_D")]

repr_plot_accepted <- rbind(
  repr_plot[accepted_sets$pareto,],
  repr_plot[accepted_sets$gof1_set,],
  repr_plot[accepted_sets$gof2_set,]
  )
repr_plot_accepted <- data.frame(repr_plot_accepted, method = rep(c("pareto","sum","wt.sum"), each = length(accepted_sets$pareto)))


fig1 <- plot_ly()
fig1 <- fig1 %>% add_markers(data = repr_plot_accepted, x = ~pS1_S2, y = ~RRS1_D, z = ~RRS2_D, size = 0.1)
fig1

fig_pareto <- fig1 %>% 
  add_mesh(
    x = repr_plot_accepted$pS1_S2[repr_plot_accepted$method == "pareto"],
    y = repr_plot_accepted$RRS1_D[repr_plot_accepted$method == "pareto"],
    z = repr_plot_accepted$RRS2_D[repr_plot_accepted$method == "pareto"],
    hoverinfo = 'text'
  ) 
fig_pareto

fig_gof1 <- fig1 %>% 
  add_markers(
    x = repr_plot_accepted$pS1_S2[repr_plot_accepted$method == "sum"],
    y = repr_plot_accepted$RRS1_D[repr_plot_accepted$method == "sum"],
    z = repr_plot_accepted$RRS2_D[repr_plot_accepted$method == "sum"],
    color = "red", size = 3,
    hoverinfo = 'text'
  ) 
fig_gof1

fig_gof2 <- fig1 %>% 
  add_markers(
    x = repr_plot_accepted$pS1_S2[repr_plot_accepted$method == "wt.sum"],
    y = repr_plot_accepted$RRS1_D[repr_plot_accepted$method == "wt.sum"],
    z = repr_plot_accepted$RRS2_D[repr_plot_accepted$method == "wt.sum"],
    color = "red", size = 3,
    hoverinfo = 'text'
  ) 
fig_gof2


# plot CE ----
df_plot <- data.frame(
  rbind(
    res_list$ce_mat[accepted_sets$pareto,1:3],
    res_list$ce_mat[accepted_sets$gof1_set,1:3],
    res_list$ce_mat[accepted_sets$gof2_set,1:3]
  ),
  method = rep(c("pareto","sum","wt.sum"), each = length(accepted_sets$pareto)),
  stringsAsFactors = F
)


fig2 <- plot_ly(type = 'scatter', 
                mode = 'markers',
                symbols = c("circle","triangle-up",3)
                ) %>%
    add_markers(data = df_plot, 
                x = ~incr_C, 
                y = ~incr_Q,
                symbol = ~method,
                mode = 'markers') %>% 
  layout(
    xaxis = list(title = "Incremental costs"), yaxis = list(title = "Incremental QALYs"))
fig2
    

# retrieve mean incremental costs, qalys, cer, and nmbs
showMeans(res_mat = res_list$ce_mat, sets = accepted_sets, round.digits = 3)


