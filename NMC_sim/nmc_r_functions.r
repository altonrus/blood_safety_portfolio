library(prevalence)
library(EnvStats)
library(data.table)
library(PearsonDS)
library(statmod)
##
#
# Net Monetary Cost
#
##

nmc_from_sim_output <- function(sim_output, transmissability_RBC, transmissability_PLT, transmissability_FFP,
                                RBC_per_don, PLT_per_don, FFP_per_don, WTP){
  sim_output <- data.table(sim_output)
  sim_output[ , cost_RBC_unit := cost_times_RBC/Tot_RBC]
  sim_output[ , cost_PLT_unit := cost_times_PLT/Tot_PLT]
  sim_output[ , cost_FFP_unit := cost_times_FFP/Tot_FFP]
  sim_output[ , QALYL_RBC_unit := QALY_times_RBC/Tot_RBC]
  sim_output[ , QALYL_PLT_unit := QALY_times_PLT/Tot_PLT]
  sim_output[ , QALYL_FFP_unit := QALY_times_FFP/Tot_FFP]
  sim_output[ , cost_RBC_donation := cost_RBC_unit*transmissability_RBC*RBC_per_don]
  sim_output[ , cost_PLT_donation := cost_PLT_unit*transmissability_PLT*PLT_per_don]
  sim_output[ , cost_FFP_donation := cost_FFP_unit*transmissability_FFP*FFP_per_don]
  sim_output[ , QALYL_RBC_donation := QALYL_RBC_unit*transmissability_RBC*RBC_per_don]
  sim_output[ , QALYL_PLT_donation := QALYL_PLT_unit*transmissability_PLT*PLT_per_don]
  sim_output[ , QALYL_FFP_donation := QALYL_FFP_unit*transmissability_FFP*FFP_per_don]
  sim_output[ , cost_per_donation := cost_RBC_donation+cost_PLT_donation+cost_FFP_donation]
  sim_output[ , QALYL_per_donation := QALYL_RBC_donation+QALYL_PLT_donation+QALYL_FFP_donation]
  sim_output[ , NMC_per_donation := cost_per_donation + WTP*QALYL_per_donation]
  return(sim_output)
}


##
#
# PROBABILISTIC SENSITIVITY ANALYSIS
#
##

gen_PSA_inputs <- function(params, n, two_index = FALSE){
  params <- as.data.table(params)
  inputs <- matrix(nrow = nrow(params), ncol = n)
  
  for (row in 1:nrow(params)){
    # Matrix with n cols, each row corresponds to
    
    
    
    if(is.na(params[row, Distribution])) {
      inputs[row, ] <- rep(params[row, Basecase], n)
    } else {
      inputs[row, ] <- samp_dist(
        dist = params[row, Distribution],
        param1 = params[row, Param1],
        param2 = params[row, Param2],
        n = n,
        mean = params[row, Basecase],
        low = params[row, Low],
        high = params[row, High])*ifelse(is.na(params[row, mult]), 1, params[row, mult])
    }
  }
  if (two_index == TRUE){
    return(cbind(params[, c("key", "row_idx", "col_idx")], inputs))
  } else {
    return(cbind(params[, c("key", "Index")], inputs))
  }
  
}


samp_dist <- function(dist, param1, param2, n, mean, low, high){
  #print(c(dist, param1, param2, n, mean, low, high))
  if(dist == "Beta"){
    if(is.na(param1) | is.na(param2)){
      dist <- betaExpert(best = mean, lower = low, upper = high, p = 0.95, method = "mean")
      return(rbeta(n, shape1 = dist$alpha, shape2 = dist$beta))
    } else {
      return (rbeta(n=n, shape1 = param1, shape2 = param2))
    }
  } else if (dist == "Gamma"){
    if(is.na(param1) | is.na(param2)){
      dist <- gammaExpert(best = mean, lower = low, upper = high, p = 0.95, method = "mean")
      return (rgamma(n=n, shape = dist$shape, scale = dist$scale))
    } else{
      return (rgamma(n=n, shape = param1, scale = param2))
    }
  } else if (dist == "Tri"){
    rtri(n = n, min = low, max = high, mode = mean)
  } else if (dist == "Expo"){
    return (rexp(n=n, rate = param1))
  } else if (dist == "Weibull"){
    return (rweibull(n=n, shape = param1, scale = param2))
  } else if (dist == "Pearson5"){
    return(rpearsonV(n=n, shape=param1, location=1, scale=param2))
  } else if (dist == "Invgauss"){
    return(rinvgauss(n=n, mean=param1, shape=param2))
  } else return(NA)
}


gammaExpert <- function(best, lower, upper, p = 0.95, method = "mode"){
    ## check presence
    if (missing(best))
      stop("'best' is missing")
    if (missing(lower) & missing(upper))
      stop("at least 'lower' or 'upper' must be specified")
    
    ## check input values: order
    if (!missing(lower))
      if (lower > best) stop("'lower' cannot be greater than 'best'")
    if (!missing(upper))
      if (upper < best) stop("'upper' cannot be smaller than 'best'")
    if (!missing(lower) & !missing(upper)) # useless??
      if (lower > upper) stop("'lower' cannot be greater than 'upper'")
    
    
    f_mean <-
      function(x, mean, p, target){
        return(
          sum(
            (qgamma(p = p,
                    shape = x,
                    scale = mean/x) -
               target) ^ 2
          ))
      }
    
    ## define 'target' and 'p'
    if (!missing(lower) & missing(upper)){
      target <- lower
      p <- 1 - p
    } else if (!missing(upper) & missing(lower)){
      target <- upper
    } else if (!missing(upper) & !missing(lower)){
      target <- c(lower, upper)
      p <- c(0, p) + (1 - p) / 2
    }
    
    ## derive a and b (=shape1 and shape2)
    
    if (method == "mean"){
      a <- optimize(f_mean, c(0, 1000),
                    mean = best, p = p, target = target)$minimum #a is shape parameter
      b <- best / a #b is scale parameter
    }
    
    ## create 'out' dataframe
    out <- list(shape = a, scale = b)
    #class(out) <- "betaExpert"
    
    ## return 'out'
    return(out)
  }
