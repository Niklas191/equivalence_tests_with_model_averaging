# === This file is related to===================================================================================================================================
# Overcoming model uncertainty - how equivalence tests can benefit from model averaging 
# Niklas Hagemann and Kathrin Moellenhoff

# === File description ========================================================================================================================================
# This file contains the code for the case study. 

# === Setup (SET ALL THESE MANUALLY) ==========================================================================================================================

ms <- c("linear", "quadratic", "exponential", "emax", "sigEmax", "betaMod")
B <- 300

# === Packages, data and functions ============================================================================================================================

library(DoseFinding)
library(dplyr)

load("wdm_data_all.Rdata")

source("Testing_functions.R")

# === Analysis ================================================================================================================================================

set.seed(12345) 
seeds <- sample(1:1000000, 1000)

result <- list()

for(i in 1:1000){
  
  set.seed(seeds[i])
  
  y_sd <- sd_list[[i]]$value
  x_sd <- sd_list[[i]]$time
  y_wd <- wd_list[[i]]$value
  x_wd <- wd_list[[i]]$time

  ranges_i <- (max(c(y_sd, y_wd)) - min(c(y_sd, y_wd)))
  
  bt <- bootstrap_test_MA(x1 = x_sd, 
                          x2 = x_wd,
                          y1 = y_sd, 
                          y2 = y_wd, 
                          model_set1 = ms, 
                          model_set2 = ms, 
                          criterion = "AIC", 
                          alpha = 0.05, 
                          B = B, 
                          grid_width = 0.001, 
                          method = "both", 
                          relevance_theshold = NULL, 
                          bnds1 = NULL, 
                          bnds2 = NULL,
                          off = NULL,
                          scal = NULL)

  w1 <- as.data.frame(t(bt$weights_group_1))
  w2 <- as.data.frame(t(bt$weights_group_2))

  result <- list(result, list(bt$u_hybrid, w1, w2, ranges_i))
}
