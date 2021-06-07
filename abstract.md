
# Evaluation of the Pareto Frontier approach for model calibration

# Purpose
To compare the performance of the Pareto Frontier approach for model calibration, recently proposed by Enns et al. (2015), against a distance-based summary score (selecting the 1% of inputs sets that minimise the sum of absolute differences between model outputs and calibration targets).

# Methods
Pareto Frontier approach is a new method for identifying best-fitting inputs in health economic model calibration. Based on the concept of Pareto optimality, it does not consider trade-offs between calibration targets. Instead, inputs are selected if and only if their fit to one target can only be improved by worsening the fit to another. While this method may seem intuitively attractive, up until now, it has not been systematically evaluated.

To assess the performance of the Pareto Frontier approach, we conducted a simulation study, using the well-documented Sick-Sicker Model. As in the original study, three parameters were calibrated. We randomly generated 10,000 true parameter sets and computed corresponding incremental net monetary benefit values. Six different calibration targets were derived directly from the (true) state transition models. In addition, we also ran micro-simulations (with n = 100-2,000 agents), to generate stochastic targets. We then ran model calibration against those targets, under different scenarios, in which we varied the number of candidate parameter sets being evaluated (500 - 50,000), and the number and types of calibration targets (2-5).

Performance was assessed by comparing the estimated mean incremental net monetary benefit, and the estimated parameter input values, against their true counterparts.


# Results

The Pareto Frontier approach generally provided less accurate results compared to the distance-based method. The mean absolute error of the incremental net monetary benefit was higher in almost all scenarios, ranging from +1,809 (4,749 vs. 2,940) to +11,125 (13,316 vs. 2,191). We also observed a much greater variability of estimates within scenarios. However, for one specific scenario with stochastic targets, the Pareto Frontier approach performed better, i.e. had a lower absolute error of -7,091 (21,095 vs. 28,185). Results for the calibration parameters were more mixed, but by and large, the distance-based method tended to provide more accurate estimates for two of the three parameters.  
# Conclusion
The Pareto Frontier approach generally performed worse than a simple, distance-based objective function. However, further research is needed to validate our findings in other, more complex models.
