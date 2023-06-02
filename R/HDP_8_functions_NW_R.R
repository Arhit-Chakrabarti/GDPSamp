# Define the target function for [beta| - ] in log scale with Metropolis ratio
# Include pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8 and alpha0 in pars
## NOTE: pars should be a named list with names pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8 and alpha0. 
# a correspons to gamma. dat is NULL here
Target_beta <- function(ycand, ycurrent, a, dat = NULL, pars) {
  L = length(dat)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha0 = as.numeric(pars$alpha0); # extract the alpha parameters
  pi1 = pars$pi1; pi2 = pars$pi2; pi3 = pars$pi3; pi4 = pars$pi4; pi5 = pars$pi5; pi6 = pars$pi6; pi7 = pars$pi7; pi8 = pars$pi8;  # extract the weight parameters. 
  # Numerator in the MH Ratio in the log scale
  proposed = sum((a/L - 1)*log(ycand) +  alpha0 * ycand * (log(pi1) + log(pi2) + log(pi3) + log(pi4) + log(pi5) + log(pi6) + log(pi7) + log(pi8)) - 8 * lgamma(alpha0 * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum((a/L - 1)*log(ycurrent) +  alpha0 * ycurrent * (log(pi1) + log(pi2) + log(pi3) + log(pi4) + log(pi5) + log(pi6) + log(pi7) + log(pi8)) - 8 * lgamma(alpha0 * ycurrent))
  
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 

# @pars argument should be a named list with a0, b0, beta, pi1, pi2, pi3, pi4, pi5, pi6, pi7, pi8
# @ycurrent corresponds to the previous value of alpha0
sample_alpha0 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta = as.numeric(pars$beta)
  a0 = pars$a0; b0 = pars$b0

  pi1 = pars$pi1; pi2 = pars$pi2; pi3 = pars$pi3; pi4 = pars$pi4; pi5 = pars$pi5; pi6 = pars$pi6; pi7 = pars$pi7; pi8 = pars$pi8;  # extract the weight parameters
  
  L = length(beta) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = L * lgamma(ycand) - 8 * sum(lgamma(ycand * beta)) + (a0 - 1)*log(ycand) - b0 * ycand + ycand *(sum(log(pi1)) + sum(log(pi2)) + sum(log(pi3)) + sum(log(pi4)) + sum(log(pi5)) + sum(log(pi6)) + sum(log(pi7)) + sum(log(pi8))) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = L * lgamma(ycurrent) - 8 * sum(lgamma(ycurrent * beta)) + (a0 - 1)*log(ycurrent) - b0 * ycurrent + ycurrent *(sum(log(pi1)) + sum(log(pi2)) + sum(log(pi3)) + sum(log(pi4)) + sum(log(pi5)) + sum(log(pi6)) + sum(log(pi7)) + sum(log(pi8))) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(list(out = out, accept = accept))
}
