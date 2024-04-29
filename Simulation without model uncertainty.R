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

  u90percentil_true <- c()
  u95percentil_true <- c()
  u975percentil_true <- c()
  u99percentil_true <- c()
  
  u90hybrid_true <- c()
  u95hybrid_true <- c()
  u975hybrid_true <- c()
  u99hybrid_true <- c()
  
  se_hybrid_true <- c()

  ds_asymp <- as.numeric()
  
  u90asymp <- c()
  u95asymp <- c()
  u975asymp <- c()
  u99asymp <- c()
  se_asymp <- c()
  
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

    ds_true <- c(ds_true, true_result$d)
    
    u90percentil_true <- c(u90percentil_true, true_result$u_percentil[1])
    u95percentil_true <- c(u95percentil_true, true_result$u_percentil[2])
    u975percentil_true <- c(u975percentil_true, true_result$u_percentil[3])
    u99percentil_true <- c(u99percentil_true, true_result$u_percentil[4])
    
    u90hybrid_true <- c(u90hybrid_true, true_result$u_hybrid[1])
    u95hybrid_true <- c(u95hybrid_true, true_result$u_hybrid[2])
    u975hybrid_true <- c(u975hybrid_true, true_result$u_hybrid[3])
    u99hybrid_true <- c(u99hybrid_true, true_result$u_hybrid[4])
    
    se_hybrid_true <- c(se_hybrid_true, true_result$se_hybrid)
    
    test_asymp <- CI_test_aysmp_emax_exp(x1 = dose1,
                                          x2 = dose2,
                                          y1 = resp1,
                                          y2 = resp2,
                                          m1 = "emax",
                                          m2 = "exponential",
                                          bnds = defBnds(max(dose1, dose2)),
                                          grid_width = 0.001)
    
    ds_asymp <- c(ds_asymp, test_asymp$d)
    
    u90asymp <- c(u90asymp, test_asymp$u_asymptotic[1])
    u95asymp <- c(u95asymp, test_asymp$u_asymptotic[2])
    u975asymp <- c(u975asymp, test_asymp$u_asymptotic[3])
    u99asymp <- c(u99asymp, test_asymp$u_asymptotic[4])
    
    se_asymp <- c(se_asymp, test_asymp$se_asymptotic)
}