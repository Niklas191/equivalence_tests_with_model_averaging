# === This file is related to===================================================================================================================================
# Overcoming model uncertainty - how equivalence tests can benefit from model averaging 
# Niklas Hagemann and Kathrin Moellenhoff

# === File description ========================================================================================================================================
# This file contains the code for the simulation of the first scenario. 

# === Setup (SET ALL THESE MANUALLY) ==========================================================================================================================

# set.seed()

# epsilon =
# beta_20 =
# sigma1 =
# sigma2 =
# n1 = n2 =


simu <- 1000
B <- 300

# === Packages and functions ==================================================================================================================================

library(DoseFinding)
library(dplyr)

source("Testing_functions.R")

# === Simulation ==============================================================================================================================================

  e0 <- beta_20

  dose1 <- rep(0:4, each = n1/5) 
  dose2 <- rep(0:4, each = n2/5) 
  
  ds_true <- c()
  ds_MA <- c()
  ds_miss1 <- c()
  ds_miss2 <- c()
  ds_miss3 <- c()
  
  u90percentil_true <- c()
  u95percentil_true <- c()
  u975percentil_true <- c()
  u99percentil_true <- c()
  
  u90hybrid_true <- c()
  u95hybrid_true <- c()
  u975hybrid_true <- c()
  u99hybrid_true <- c()
  
  se_hybrid_true <- c()
  
  u90percentil_MA <- c()
  u95percentil_MA <- c()
  u975percentil_MA <- c()
  u99percentil_MA <- c()
  
  u90hybrid_MA <- c()
  u95hybrid_MA <- c()
  u975hybrid_MA <- c()
  u99hybrid_MA <- c()
  
  se_hybrid_MA <- c()
  
  u90percentil_miss1 <- c()
  u95percentil_miss1 <- c()
  u975percentil_miss1 <- c()
  u99percentil_miss1 <- c()
  
  u90hybrid_miss1 <- c()
  u95hybrid_miss1 <- c()
  u975hybrid_miss1 <- c()
  u99hybrid_miss1 <- c()
  
  se_hybrid_miss1 <- c()
  
  u90percentil_miss2 <- c()
  u95percentil_miss2 <- c()
  u975percentil_miss2 <- c()
  u99percentil_miss2 <- c()
  
  u90hybrid_miss2 <- c()
  u95hybrid_miss2 <- c()
  u975hybrid_miss2 <- c()
  u99hybrid_miss2 <- c()
  
  se_hybrid_miss2 <- c()
  
  u90percentil_miss3 <- c()
  u95percentil_miss3 <- c()
  u975percentil_miss3 <- c()
  u99percentil_miss3 <- c()
  
  u90hybrid_miss3 <- c()
  u95hybrid_miss3 <- c()
  u975hybrid_miss3 <- c()
  u99hybrid_miss3 <- c()
  
  se_hybrid_miss3 <- c()
  
  for(i in 1:simu){
    resp1 <- DoseFinding::emax(dose = dose1, e0 = 1, eMax = 2, ed50 = 1) + rnorm(n = n1, mean = 0, sd = sigma1)
    resp2 <- DoseFinding::exponential(dose = dose2, e0 = e0, e1 = 2.2, delta = 8) + rnorm(n = n2, mean = 0, sd = sigma2)
    
    true_result <- CI_test(x1 = dose1,
                           x2 = dose2,
                           y1 = resp1,
                           y2 = resp2,
                           m1 = "emax",
                           m2 = "exponential",
                           B = B,
                           off = NULL,
                           scal = NULL,
                           bnds1 = defBnds(max(dose1, dose2))$emax,
                           bnds2 = defBnds(max(dose1, dose2))$exponential,
                           grid_width = 0.001)
    
    MA_result <- bootstrap_test_MA(x1 = dose1,
                                   x2 = dose2,
                                   y1 = resp1,
                                   y2 = resp2,
                                   model_set1 =  c("emax", "exponential"),
                                   model_set2 =  c("emax", "exponential"),
                                   B = B,
                                   criterion = "AIC", 
                                   alpha = c(0.1, 0.05, 0.025, 0.01), 
                                   grid_width = 0.001, 
                                   method = "both", 
                                   relevance_theshold = NULL, 
                                   bnds1 = NULL, 
                                   bnds2 = NULL,
                                   off = NULL,
                                   scal = NULL)
    
    miss_result1 <- CI_test(x1 = dose1,
                            x2 = dose2,
                            y1 = resp1,
                            y2 = resp2,
                            m1 = "emax",
                            m2 = "emax",
                            B = B,
                            off = NULL,
                            scal = NULL,
                            bnds1 = defBnds(max(dose1, dose2))$emax,
                            bnds2 = defBnds(max(dose1, dose2))$emax,
                            grid_width = 0.001)
    
    miss_result2 <- CI_test(x1 = dose1,
                            x2 = dose2,
                            y1 = resp1,
                            y2 = resp2,
                            m1 = "exponential",
                            m2 = "exponential",
                            B = B,
                            off = NULL,
                            scal = NULL,
                            bnds1 = defBnds(max(dose1, dose2))$exponential,
                            bnds2 = defBnds(max(dose1, dose2))$exponential,
                            grid_width = 0.001)
    
    miss_result3 <- CI_test(x1 = dose1,
                            x2 = dose2,
                            y1 = resp1,
                            y2 = resp2,
                            m1 = "exponential",
                            m2 = "emax",
                            B = B,
                            off = NULL,
                            scal = NULL,
                            bnds1 = defBnds(max(dose1, dose2))$exponential,
                            bnds2 = defBnds(max(dose1, dose2))$emax,
                            grid_width = 0.001)
    
    ds_true <- c(ds_true, true_result$d)
    ds_MA <- c(ds_MA, MA_result$d)
    ds_miss1 <- c(ds_miss1, miss_result1$d)
    ds_miss2 <- c(ds_miss2, miss_result2$d)
    ds_miss3 <- c(ds_miss3, miss_result3$d)
    
    u90percentil_true <- c(u90percentil_true, true_result$u_percentil[1])
    u95percentil_true <- c(u95percentil_true, true_result$u_percentil[2])
    u975percentil_true <- c(u975percentil_true, true_result$u_percentil[3])
    u99percentil_true <- c(u99percentil_true, true_result$u_percentil[4])
    
    u90hybrid_true <- c(u90hybrid_true, true_result$u_hybrid[1])
    u95hybrid_true <- c(u95hybrid_true, true_result$u_hybrid[2])
    u975hybrid_true <- c(u975hybrid_true, true_result$u_hybrid[3])
    u99hybrid_true <- c(u99hybrid_true, true_result$u_hybrid[4])
    
    se_hybrid_true <- c(se_hybrid_true, true_result$se_hybrid)

    u90percentil_MA <- c(u90percentil_MA, MA_result$u_percentil[1])
    u95percentil_MA <- c(u95percentil_MA, MA_result$u_percentil[2])
    u975percentil_MA <- c(u975percentil_MA, MA_result$u_percentil[3])
    u99percentil_MA <- c(u99percentil_MA, MA_result$u_percentil[4])
    
    u90hybrid_MA <- c(u90hybrid_MA, MA_result$u_hybrid[1])
    u95hybrid_MA <- c(u95hybrid_MA, MA_result$u_hybrid[2])
    u975hybrid_MA <- c(u975hybrid_MA, MA_result$u_hybrid[3])
    u99hybrid_MA <- c(u99hybrid_MA, MA_result$u_hybrid[4])
    
    se_hybrid_MA <- c(se_hybrid_MA, MA_result$se_hybrid)
    
    u90percentil_miss1 <- c(u90percentil_miss1, miss_result1$u_percentil[1])
    u95percentil_miss1 <- c(u95percentil_miss1, miss_result1$u_percentil[2])
    u975percentil_miss1 <- c(u975percentil_miss1, miss_result1$u_percentil[3])
    u99percentil_miss1 <- c(u99percentil_miss1, miss_result1$u_percentil[4])
    
    u90hybrid_miss1 <- c(u90hybrid_miss1, miss_result1$u_hybrid[1])
    u95hybrid_miss1 <- c(u95hybrid_miss1, miss_result1$u_hybrid[2])
    u975hybrid_miss1 <- c(u975hybrid_miss1, miss_result1$u_hybrid[3])
    u99hybrid_miss1 <- c(u99hybrid_miss1, miss_result1$u_hybrid[4])
    
    se_hybrid_miss1 <- c(se_hybrid_miss1, miss_result1$se_hybrid)

    u90percentil_miss2 <- c(u90percentil_miss2, miss_result2$u_percentil[1])
    u95percentil_miss2 <- c(u95percentil_miss2, miss_result2$u_percentil[2])
    u975percentil_miss2 <- c(u975percentil_miss2, miss_result2$u_percentil[3])
    u99percentil_miss2 <- c(u99percentil_miss2, miss_result2$u_percentil[4])
    
    u90hybrid_miss2 <- c(u90hybrid_miss2, miss_result2$u_hybrid[1])
    u95hybrid_miss2 <- c(u95hybrid_miss2, miss_result2$u_hybrid[2])
    u975hybrid_miss2 <- c(u975hybrid_miss2, miss_result2$u_hybrid[3])
    u99hybrid_miss2 <- c(u99hybrid_miss2, miss_result2$u_hybrid[4])
    
    se_hybrid_miss2 <- c(se_hybrid_miss2, miss_result2$se_hybrid)

    u90percentil_miss3 <- c(u90percentil_miss3, miss_result3$u_percentil[1])
    u95percentil_miss3 <- c(u95percentil_miss3, miss_result3$u_percentil[2])
    u975percentil_miss3 <- c(u975percentil_miss3, miss_result3$u_percentil[3])
    u99percentil_miss3 <- c(u99percentil_miss3, miss_result3$u_percentil[4])
    
    u90hybrid_miss3 <- c(u90hybrid_miss3, miss_result3$u_hybrid[1])
    u95hybrid_miss3 <- c(u95hybrid_miss3, miss_result3$u_hybrid[2])
    u975hybrid_miss3 <- c(u975hybrid_miss3, miss_result3$u_hybrid[3])
    u99hybrid_miss3 <- c(u99hybrid_miss3, miss_result3$u_hybrid[4])
    
    se_hybrid_miss3 <- c(se_hybrid_miss3, miss_result3$se_hybrid)
}