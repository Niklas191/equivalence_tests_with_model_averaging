# === This file is related to===================================================================================================================================
# Overcoming model uncertainty - how equivalence tests can benefit from model averaging 
# Niklas Hagemann and Kathrin Moellenhoff

# === File description =========================================================================================================================================
# This file contains the functions for the simulation study and the case study.

# === Functions ================================================================================================================================================

library(DoseFinding)

CI_test_aysmp_emax_exp <- function(x1, x2, y1, y2, m1, m2, bnds, grid_width = 0.001){
  
  mod1 <- fitMod(dose = x1, resp = y1, model = m1, bnds = bnds$emax)
  mod2 <- fitMod(dose = x2, resp = y2, model = m2, bnds = bnds$exponential)
  
  n1 <- length(y1)
  n2 <- length(y2)
  n <- length(y1) + length(y2)
  lambda <- (n1+n2)/n1
  
  coef1 <- mod1$coefs
  coef2 <- mod2$coefs
  
  sd1 <- sqrt(mod1$RSS/length(y1))
  sd2 <- sqrt(mod2$RSS/length(y2))
  
  dose_grid <- seq(min(x1,x2), max(x1,x2), grid_width)
  
  curve1 <- predict(mod1, predType = "ls-means", doseSeq = dose_grid)
  curve2 <- predict(mod2, predType = "ls-means", doseSeq = dose_grid)
  
  d <- max(abs(curve1 - curve2))
  x0 <- dose_grid[which.max(abs(curve1 - curve2))]
  
  G1=function(x){as.numeric(n1*(emaxGrad(dose = x0, eMax = coef1[2], ed50 = coef1[3]) %*% vcov(mod1) %*% t(emaxGrad(dose = x0, eMax = coef1[2], ed50 = coef1[3]))))}
  G2=function(x){as.numeric(n2*(exponentialGrad(dose = x0, e1 = coef2[2], delta = coef2[3]) %*% vcov(mod2) %*% t(exponentialGrad(dose = x0, e1 = coef2[2], delta = coef2[3]))))}
  ker=function(x){lambda*G1(x)+lambda/(lambda-1)*G2(x)}
  var_asymp=ker(x0)
  
  return_values_asymptotic <- d + sqrt(var_asymp/n) * qnorm(c(0.9, 0.95, 0.975, 0.99))
  
  names(return_values_asymptotic) <- c("alpha = 0.1", "alpha = 0.05", "alpha = 0.025", "alpha = 0.01")
  
  return(list(d = d, 
              u_asymptotic = return_values_asymptotic,
              se_asymptotic = sqrt(var_asymp/n)
  ))
}

CI_test <- function(x1, x2, y1, y2, m1, m2, alpha = c(0.1, 0.05, 0.025, 0.01), B = 300, bnds1, bnds2, off = NULL, scal = NULL, grid_width = 0.001, method = "both"){
  
  if(is.null(off)){
    off <- 0.01*max(c(x1, x2))
  }
  if(is.null(scal)){
    scal <- 1.2*max(c(x1, x2))
  }
  
  mod1 <- fitMod(dose = x1, resp = y1, model = m1, bnds = bnds1)
  mod2 <- fitMod(dose = x2, resp = y2, model = m2, bnds = bnds2)

  sd1 <- sqrt(mod1$RSS/length(y1))
  sd2 <- sqrt(mod2$RSS/length(y2))
  
  dose_grid <- seq(min(x1,x2), max(x1,x2), grid_width)
  
  curve1 <- predict(mod1, predType = "ls-means", doseSeq = dose_grid)
  curve2 <- predict(mod2, predType = "ls-means", doseSeq = dose_grid)
  
  d <- max(abs(curve1 - curve2))
  
  d_stars <- as.numeric()
  
  for (j in 1:B){
    y1_step <- rnorm(n = length(x1), mean = predict(mod1, predType = "ls-means"), sd = sd1)
    y2_step <- rnorm(n = length(x2), mean = predict(mod2, predType = "ls-means"), sd = sd2)
    
    mod1_step <- fitMod(dose = x1, resp = y1_step, model = m1, bnds = bnds1)
    mod2_step <- fitMod(dose = x2, resp = y2_step, model = m2, bnds = bnds2)
    
    curve1_step <- predict(mod1_step, predType = "ls-means", doseSeq = dose_grid)
    curve2_step <- predict(mod2_step, predType = "ls-means", doseSeq = dose_grid)
    
    d_star <- max(abs(curve1_step - curve2_step))
    
    d_stars <- c(d_stars, d_star)
  }
  
  sigma_hat <- sd(d_stars)
  
  return_values_percentil <- quantile(d_stars, probs = (1-alpha))
  return_values_hybrid <- d + sigma_hat * qnorm(p = (1-alpha))
  names(return_values_percentil) <- paste("alpha =", alpha)
  names(return_values_hybrid) <- names(return_values_percentil)
  
  if(method == "both"){
    return(list(d = d, 
                u_percentil = return_values_percentil,
                u_hybrid = return_values_hybrid,
                se_hybrid = sigma_hat))
  } else if(method == "percentil"){
    return(list(d = d, u_percentil = return_values_percentil))
  } else if(method == "hybrid"){
    return(list(d = d,
                u_hybrid = return_values_hybrid,
                se_hybrid = sigma_hat))
  } else {
    warning("Specification of argument method is not valid. method = \'both\' is used.")
    return(list(d = d, 
                u_percentil = return_values_percentil,
                u_hybrid = return_values_hybrid,
                se_hybrid = sigma_hat))
  }
}

set_up <- function(x1, x2, model_set1, model_set2, off , scal, bnds1, bnds2){
  builtIn <- c("linlog", "linear", "quadratic", "emax", "exponential", 
               "logistic", "betaMod", "sigEmax")

  mods <- c(linlog = function(d, e) {TestingSimilarity:::linlog(d, e, off = off)}, 
            TestingSimilarity:::linear, 
            TestingSimilarity:::quadratic, 
            TestingSimilarity:::emax, 
            TestingSimilarity:::exponential, 
            TestingSimilarity:::logistic, 
            betaMod = function(d,e) {TestingSimilarity:::betaMod(d, e, scal = scal)}, 
            TestingSimilarity:::sigEmax)
  
  models1 <- mods[builtIn %in% model_set1]
  models2 <- mods[builtIn %in% model_set2]
  
  models1char <- builtIn[builtIn %in% model_set1]
  models2char <- builtIn[builtIn %in% model_set2]

  if (is.null(bnds1)) {
    bnds1 = defBnds(max(x1))
  }
  if (is.null(bnds2)) {
    bnds2 = defBnds(max(x2))
  }
  
  return(list(models1 = models1, 
              models2 = models2, 
              models1char = models1char, 
              models2char = models2char, 
              bnds1 = bnds1, 
              bnds2 = bnds2))
}

use_fitMod <- function(x, y, model_set, bnds, off, scal){
  
  mods <- list()
  
  for(i in 1:length(model_set)){
    
    if(model_set[i] %in% names(bnds)){
      mod <- fitMod(x, y, model = model_set[i], bnds = bnds[[(model_set[i])]], addArgs = list(scal = scal, off = off))
    } else {
      mod <- fitMod(x, y, model = model_set[i], addArgs = list(scal = scal, off = off))
    }
    
    mods[[i]] <- mod
  }
  
  return(mods)  
}

fitted_values_MA <- function(x, y, mods, weights){
  n <- length(y)
  
  model_values <- matrix(0, nrow = n, ncol = length(mods))
  
  for(j in 1:length(mods)){
    model_values[,j] <- predict(mods[[j]], predType = "ls-means")
  }
  
  fitted_values <- drop(model_values %*% weights)
  sd <- sqrt((sum((y - fitted_values)^2))/n)
  
  return(list(fitted_values = fitted_values, sd = sd)) 
}

smooth_weights <- function(ICs, robust = TRUE){
  k <- length(ICs)
  
  if(robust == TRUE){ # The weights depend only on the difference of the AICs/BICs, not on their absolute value. For very large AIC/BIC values the denominator of the weighting function can become numerically zero. Therefore, we rescale the AIC/BIS. 
    ICs <- ICs - mean(ICs)
  } 
  
  denominator <- sum(exp(-1/2 * ICs))
  
  weights <- exp(-1/2* ICs) / denominator
  
  return(weights)
}

distance_and_weights <- function(mods1, mods2, dose_grid, criterion, n1, n2, round = FALSE, robust = TRUE) {
  
  ICs1 <- c()
  ICs2 <- c()
  
  if(criterion == "AIC"){
    for(j in 1:length(mods1)){
      ICs1[j] <- AIC(mods1[[j]])
    }
    
    for(j in 1:length(mods2)){
      ICs2[j] <- AIC(mods2[[j]])
    }
  } else if(criterion == "BIC") {
    for(j in 1:length(mods1)){
      ICs1[j] <- AIC(mods1[[j]], k = log(n1))
    }
    
    for(j in 1:length(mods2)){
      ICs2[j] <- AIC(mods2[[j]], k = log(n2))
    }
  } else {
    stop("Argument criterion must be \'AIC\' or \'BIC\'.")
  }
  
  w1 <- smooth_weights(ICs = ICs1, robust = robust)
  w2 <- smooth_weights(ICs = ICs2, robust = robust) 

  mod_fit1 <- matrix(nrow = length(dose_grid), ncol = length(mods1))
  mod_fit2 <- matrix(nrow = length(dose_grid), ncol = length(mods2))
  
  for(j in 1:length(mods1)){
    mod_fit1[,j] <- predict(mods1[[j]], predType = "ls-means", doseSeq = dose_grid)
  }
  
  for(j in 1:length(mods2)){
    mod_fit2[,j] <- predict(mods2[[j]], predType = "ls-means", doseSeq = dose_grid)
  }
  
  m1 <- drop(mod_fit1 %*% w1)
  m2 <- drop(mod_fit2 %*% w2)
  
  distance <- max(abs(m1 - m2))
  max_location <- dose_grid[which.max(abs(m1 - m2))]
  
  if(round == FALSE){
    return(list(distance = distance, max_location = max_location, weights_group_1 = w1, weights_group_2 = w2))
  } else {
    return(list(distance = round(distance, 3), 
                max_location = round(max_location, 3), 
                weights_group_1 = round(w1, 3), 
                weights_group_2 = round(w2, 3)))
  }
}

bootstrap_test_MA <- function(x1, x2, y1, y2, model_set1, model_set2, criterion, alpha, B, grid_width = 0.001, method = "both", relevance_theshold = NULL, bnds1 = NULL, bnds2 = NULL, scal = NULL, off = NULL) {
  
  # bnds1, bnds2 not being NULL is currently only supported if bnds1 an bnds2 are given in full and in correct format as it would be generated by defBnds().
  
  if(is.null(off)){
    off <- 0.01*max(c(x1, x2))
  }
  if(is.null(scal)){
    scal <- 1.2*max(c(x1, x2))
  }
  
  init <- set_up(x1 = x1, x2 = x2, model_set1 = model_set1, model_set2 = model_set2, off = off, scal = scal, bnds1 = bnds1, bnds2 = bnds2)

  models1char <- init$models1char
  models2char <- init$models2char
  bnds1 <- init$bnds1 
  bnds2 <- init$bnds2
  n1 <- length(y1)
  n2 <- length(y2)
  
  dose_grid <- seq(min(x1,x2), max(x1,x2), grid_width)
  
  fitMod_result1 <- use_fitMod(x = x1, y = y1, model_set = models1char, bnds = bnds1, scal = scal, off = off)
  fitMod_result2 <- use_fitMod(x = x2, y = y2, model_set = models2char, bnds = bnds2, scal = scal, off = off)

  if(! is.null(relevance_theshold)){
    if(relevance_theshold < 0) {stop("Model relevance theshold \'relevance_theshold\' cannot be negative!")}
    if(relevance_theshold > 1) {stop("Model relevance theshold \'relevance_theshold\' cannot be >1!")}
    if(relevance_theshold > 0.05) {warning("You have chosen a large model relevance theshold \'relevance_theshold\'. Usually this value is smaller than 0.05, e.g. 0.01 or 0.001.")}
    
    dw <- distance_and_weights(mods1 = fitMod_result1, mods2 = fitMod_result2, dose_grid = dose_grid, criterion = criterion, n1 = length(y1), n2 = length(y2))
    
    w1 <- dw$weights_group_1
    w2 <- dw$weights_group_2
    
    w1_indicator <- w1 > relevance_theshold
    w2_indicator <- w2 > relevance_theshold
    
    if(! all(w1_indicator, w2_indicator)){

      models1char <- models1char[w1_indicator]
      models2char <- models2char[w2_indicator]
      
      fitMod_result1 <- use_fitMod(x = x1, y = y1, model_set = models1char, bnds = bnds1, scal = scal, off = off)
      fitMod_result2 <- use_fitMod(x = x2, y = y2, model_set = models2char, bnds = bnds2, scal = scal, off = off)
    }
  }
  
  dw <- distance_and_weights(mods1 = fitMod_result1, mods2 = fitMod_result2, dose_grid = dose_grid, criterion = criterion, n1 = n1, n2 =n2)
  
  d <- dw$distance
  w1 <- dw$weights_group_1
  w2 <- dw$weights_group_2
  
  names(w1) <- models1char
  names(w2) <- models2char
  
  fitted1 <- fitted_values_MA(x = x1, y = y1, mods = fitMod_result1, weights = w1)
  fitted2 <- fitted_values_MA(x = x2, y = y2, mods = fitMod_result2, weights = w2)
  
  d_stars <- c()
  
  for (i in 1:B) {
    step_y1 <- rnorm(n = n1, mean = fitted1$fitted_values, sd = fitted1$sd)
    step_y2 <- rnorm(n = n2, mean = fitted2$fitted_values, sd = fitted2$sd)
    
    fitMod_result_step_1 <- use_fitMod(x = x1, y = step_y1, model_set = models1char, bnds = bnds1, scal = scal, off = off)
    fitMod_result_step_2 <- use_fitMod(x = x2, y = step_y2, model_set = models2char, bnds = bnds2, scal = scal, off = off)
    
    dw_step <- distance_and_weights(mods1 = fitMod_result_step_1, mods2 = fitMod_result_step_2, dose_grid = dose_grid, criterion = criterion, n1 = n1, n2 =n2)

    d_star <- dw_step$distance
    d_stars <- c(d_stars, d_star)
  }
  
  se_d <- sd(d_stars)
  
  return_values_percentil <- quantile(d_stars, probs = (1-alpha))
  return_values_hybrid <- d + se_d * qnorm(p = (1-alpha))
  names(return_values_percentil) <- paste("alpha =", alpha)
  names(return_values_hybrid) <- names(return_values_percentil)
  
  if(method == "both"){
    return(list(d = d, 
                u_percentil = return_values_percentil,
                u_hybrid = return_values_hybrid,
                se_hybrid = se_d, 
                weights_group_1 = w1, 
                weights_group_2 = w2))
  } else if(method == "percentil"){
    return(list(d = d, 
                u_percentil = return_values_percentil, 
                weights_group_1 = w1, 
                weights_group_2 = w2))
  } else if(method == "hybrid"){
    return(list(d = d,
                u_hybrid = return_values_hybrid,
                se_hybrid = se_d, 
                weights_group_1 = w1, 
                weights_group_2 = w2))
  } else {
    warning("Specification of argument method is not valid. method = \'both\' is used.")
    return(list(d = d, 
                u_percentil = return_values_percentil,
                u_hybrid = return_values_hybrid,
                se_hybrid = se_d, 
                weights_group_1 = w1, 
                weights_group_2 = w2))
  }
}



