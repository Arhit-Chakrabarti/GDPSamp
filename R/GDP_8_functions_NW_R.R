# @pars argument should be a named list with alpha0, alpha2, alpha3, alpha4 and beta1
# @ycurrent corresponds to the previous value of alpha1
sample_alpha1 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta1 = as.numeric(pars$beta1)
  alpha0 = pars$alpha0; alpha2 = pars$alpha2; alpha3 = pars$alpha3; alpha4 = pars$alpha4
  L = length(beta1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = alpha0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = -2*lgamma(ycand) - L*lgamma(ycand/L) + sum(ycand/L * log(beta1)) + (alpha0 - 1)*log(ycand) - (1 - log(alpha2*alpha3*alpha4))*ycand + dgamma(ycurrent, shape = alpha0, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = -2*lgamma(ycurrent) - L*lgamma(ycurrent/L) + sum(ycurrent/L * log(beta1)) + (alpha0 - 1)*log(ycurrent) - (1 - log(alpha2*alpha3*alpha4))*ycurrent + dgamma(ycand, shape = alpha0, log = TRUE)
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


# @pars argument should be a named list with alpha1, alpha3, alpha4, alpha5, alpha6,
# beta1, beta2, nu1, nu2, eta
# @ycurrent corresponds to the previous value of alpha2
sample_alpha2 <- function(ycurrent, pars){
  # Extract the parameters from the pars list
  beta1 = as.numeric(pars$beta1); beta2 = as.numeric(pars$beta2)
  nu1 = as.numeric(pars$nu1); nu2 = as.numeric(pars$nu2); eta = as.numeric(pars$eta)
  
  alpha1 = pars$alpha1; alpha3 = pars$alpha3; alpha4 = pars$alpha4; alpha5 = pars$alpha5; alpha6 = pars$alpha6   
  L = length(beta1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = alpha1) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) + lgamma(2*(ycand + alpha3 + alpha4)) - ((1 - log(alpha5*alpha6) - (sum(beta1*log(beta2) + beta1*log(nu1) + beta1*log(nu2) + 2*beta1*log(eta))))*ycand) + (alpha1-1)*log(ycand) - (sum(lgamma(ycand*beta1) + lgamma((ycand + alpha4)*beta1) + lgamma((ycand + alpha3)*beta1) + lgamma(2*(ycand + alpha3 + alpha4)*beta1))) + dgamma(ycurrent, shape = alpha1, log = TRUE) 
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) + lgamma(2*(ycurrent + alpha3 + alpha4)) - ((1 - log(alpha5*alpha6) - (sum(beta1*log(beta2) + beta1*log(nu1) + beta1*log(nu2) + 2*beta1*log(eta))))*ycurrent) + (alpha1-1)*log(ycurrent) - (sum(lgamma(ycurrent*beta1) + lgamma((ycurrent + alpha4)*beta1) + lgamma((ycurrent + alpha3)*beta1) + lgamma(2*(ycurrent + alpha3 + alpha4)*beta1))) + dgamma(ycand, shape = alpha1, log = TRUE) 
  
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

# @pars argument should be a named list with alpha1, alpha2, alpha4, alpha6, alpha7,
# beta1, beta3, nu2, nu3, eta
# @ycurrent corresponds to the previous value of alpha2
sample_alpha3 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta1 = as.numeric(pars$beta1); beta3 = as.numeric(pars$beta3)
  nu2 = as.numeric(pars$nu2); nu3 = as.numeric(pars$nu3); eta = as.numeric(pars$eta)
  
  alpha1 = pars$alpha1; alpha2 = pars$alpha2; alpha4 = pars$alpha4; alpha6 = pars$alpha6; alpha7 = pars$alpha7   
  L = length(beta1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = alpha1) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) + lgamma(2*(ycand + alpha2 + alpha4)) - ((1 - log(alpha6*alpha7) - (sum(beta1*log(beta3) + beta1*log(nu2) + beta1*log(nu3) + 2*beta1*log(eta))))*ycand) + (alpha1-1)*log(ycand) - (sum(lgamma(ycand*beta1) + lgamma((ycand + alpha2)*beta1) + lgamma((ycand + alpha4)*beta1) + lgamma(2*(ycand + alpha2 + alpha4)*beta1))) + dgamma(ycurrent, shape = alpha1, log = TRUE) 
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) + lgamma(2*(ycurrent + alpha2 + alpha4)) - ((1 - log(alpha6*alpha7) - (sum(beta1*log(beta3) + beta1*log(nu2) + beta1*log(nu3) + 2*beta1*log(eta))))*ycurrent) + (alpha1-1)*log(ycurrent) - (sum(lgamma(ycurrent*beta1) + lgamma((ycurrent + alpha2)*beta1) + lgamma((ycurrent + alpha4)*beta1) + lgamma(2*(ycurrent + alpha2 + alpha4)*beta1))) + dgamma(ycand, shape = alpha1, log = TRUE) 
  
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

# @pars argument should be a named list with alpha1, alpha2, alpha3, alpha5, alpha7,
# beta1, beta4, nu1, nu3, eta
# @ycurrent corresponds to the previous value of alpha2
sample_alpha4 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta1 = as.numeric(pars$beta1); beta4 = as.numeric(pars$beta4)
  nu1 = as.numeric(pars$nu1); nu3 = as.numeric(pars$nu3); eta = as.numeric(pars$eta)
  
  alpha1 = pars$alpha1; alpha2 = pars$alpha2; alpha3 = pars$alpha3; alpha5 = pars$alpha5; alpha7 = pars$alpha7   
  L = length(beta1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = alpha1) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) + lgamma(2*(ycand + alpha2 + alpha3)) - ((1 - log(alpha5*alpha7) - (sum(beta1*log(beta4) + beta1*log(nu1) + beta1*log(nu3) + 2*beta1*log(eta))))*ycand) + (alpha1-1)*log(ycand) - (sum(lgamma(ycand*beta1) + lgamma((ycand + alpha2)*beta1) + lgamma((ycand + alpha3)*beta1) + lgamma(2*(ycand + alpha2 + alpha3)*beta1))) + dgamma(ycurrent, shape = alpha1, log = TRUE) 
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) + lgamma(2*(ycurrent + alpha2 + alpha3)) - ((1 - log(alpha5*alpha7) - (sum(beta1*log(beta4) + beta1*log(nu1) + beta1*log(nu3) + 2*beta1*log(eta))))*ycurrent) + (alpha1-1)*log(ycurrent) - (sum(lgamma(ycurrent*beta1) + lgamma((ycurrent + alpha2)*beta1) + lgamma((ycurrent + alpha3)*beta1) + lgamma(2*(ycurrent + alpha2 + alpha3)*beta1))) + dgamma(ycand, shape = alpha1, log = TRUE) 
  
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

# @pars argument should be a named list with alpha2, alpha4, alpha6, alpha7, alpha8
# beta5, nu1, nu4, eta
# @ycurrent corresponds to the previous value of alpha2
sample_alpha5 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta5 = as.numeric(pars$beta5)
  nu1 = as.numeric(pars$nu1); nu4 = as.numeric(pars$nu4); eta = as.numeric(pars$eta)
  
  alpha2 = pars$alpha2; alpha4 = pars$alpha4; alpha6 = pars$alpha6; alpha7 = pars$alpha7; alpha8 = pars$alpha8;    
  L = length(beta5) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = (alpha2 + alpha4)) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) - ((1 - log(alpha8) - (sum(nu1*log(beta5) + eta*log(nu4))))*ycand) + ((alpha2 + alpha4 - 1)*log(ycand)) - (sum(lgamma(ycand*nu1) + lgamma((ycand + alpha6 + alpha7)*eta))) + dgamma(ycurrent, shape = (alpha2 + alpha4), log = TRUE) 
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) - ((1 - log(alpha8) - (sum(nu1*log(beta5) + eta*log(nu4))))*ycurrent) + ((alpha2 + alpha4 - 1)*log(ycurrent)) - (sum(lgamma(ycurrent*nu1) + lgamma((ycurrent + alpha6 + alpha7)*eta))) + dgamma(ycand, shape = (alpha2 + alpha4), log = TRUE) 
  
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

# @pars argument should be a named list with alpha2, alpha3, alpha5, alpha7, alpha8
# beta6, nu2, nu4, eta
# @ycurrent corresponds to the previous value of alpha2
sample_alpha6 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta6 = as.numeric(pars$beta6)
  nu2 = as.numeric(pars$nu2); nu4 = as.numeric(pars$nu4); eta = as.numeric(pars$eta)
  
  alpha2 = pars$alpha2; alpha3 = pars$alpha3; alpha5 = pars$alpha5; alpha7 = pars$alpha7; alpha8 = pars$alpha8;    
  L = length(beta6) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = (alpha2 + alpha3)) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) - ((1 - log(alpha8) - (sum(nu2*log(beta6) + eta*log(nu4))))*ycand) + ((alpha2 + alpha3 - 1)*log(ycand)) - (sum(lgamma(ycand*nu2) + lgamma((ycand + alpha5 + alpha7)*eta))) + dgamma(ycurrent, shape = (alpha2 + alpha3), log = TRUE) 
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) - ((1 - log(alpha8) - (sum(nu2*log(beta6) + eta*log(nu4))))*ycurrent) + ((alpha2 + alpha3 - 1)*log(ycurrent)) - (sum(lgamma(ycurrent*nu2) + lgamma((ycurrent + alpha5 + alpha7)*eta))) + dgamma(ycand, shape = (alpha2 + alpha3), log = TRUE) 
  
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

# @pars argument should be a named list with alpha3, alpha4, alpha5, alpha6, alpha8
# beta7, nu3, nu4, eta
# @ycurrent corresponds to the previous value of alpha2
sample_alpha7 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta7 = as.numeric(pars$beta7)
  nu3 = as.numeric(pars$nu3); nu4 = as.numeric(pars$nu4); eta = as.numeric(pars$eta)
  
  alpha3 = pars$alpha3; alpha4 = pars$alpha4; alpha5 = pars$alpha5; alpha6 = pars$alpha6; alpha8 = pars$alpha8;    
  L = length(beta7) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = (alpha3 + alpha4)) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) - ((1 - log(alpha8) - (sum(nu3*log(beta7) + eta*log(nu4))))*ycand) + ((alpha3 + alpha4 - 1)*log(ycand)) - (sum(lgamma(ycand*nu3) + lgamma((ycand + alpha5 + alpha6)*eta))) + dgamma(ycurrent, shape = (alpha3 + alpha4), log = TRUE) 
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) - ((1 - log(alpha8) - (sum(nu3*log(beta7) + eta*log(nu4))))*ycurrent) + ((alpha3 + alpha4 - 1)*log(ycurrent)) - (sum(lgamma(ycurrent*nu3) + lgamma((ycurrent + alpha5 + alpha6)*eta))) + dgamma(ycand, shape = (alpha3 + alpha4), log = TRUE) 
  
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

# @pars argument should be a named list with alpha5, alpha6, alpha7
# beta8, nu4
# @ycurrent corresponds to the previous value of alpha2
sample_alpha8 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta8 = as.numeric(pars$beta8)
  nu4 = as.numeric(pars$nu4)
  
  alpha5 = pars$alpha5; alpha6 = pars$alpha6; alpha7 = pars$alpha7;    
  L = length(beta8) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = (alpha5 + alpha6 + alpha7)) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) - ((1 - sum(nu4*log(beta8)))*ycand) + ((alpha5 + alpha6 + alpha7 - 1)*log(ycand)) - sum(lgamma(ycand*nu4)) + dgamma(ycurrent, shape = (alpha5 + alpha6 + alpha7), log = TRUE) 
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) - ((1 - sum(nu4*log(beta8)))*ycurrent) + ((alpha5 + alpha6 + alpha7 - 1)*log(ycurrent)) - sum(lgamma(ycurrent*nu4)) + dgamma(ycand, shape = (alpha5 + alpha6 + alpha7), log = TRUE) 
  
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

####################################################################################################
# SIMPLEX TARGETS
####################################################################################################
# Define the target function for [nu1| - ] in log scale with Metropolis ratio
# @pars includes beta1, beta5, alpha2, alpha4, alpha5
## NOTE: pars should be a named list with names beta1, beta5,  alpha2, alpha4 and alpha5. beta1, beta5 should be on a simplex 
# a corresponds to alpha2 + alpha4. dat is NULL for nu1 as there is no data corresponding to nu1
Target_nu1 <- function(ycand, ycurrent, a, dat = NULL, pars) {
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha2 = as.numeric(pars$alpha2); alpha4 = as.numeric(pars$alpha4); alpha5 = as.numeric(pars$alpha5) # extract the alpha parameters
  beta1 = pars$beta1; beta5 = pars$beta5 # extract the weight parameters. 
  L = length(beta1)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Numerator in the MH Ratio in the log scale
  proposed = sum(((a * beta1)  - 1)*log(ycand) + (alpha5 * ycand) * log(beta5) - lgamma(alpha5 * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum(((a * beta1)  - 1)*log(ycurrent) + (alpha5 * ycurrent) * log(beta5) - lgamma(alpha5 * ycurrent))
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 

# Define the target function for [nu2| - ] in log scale with Metropolis ratio
# @pars includes beta1, beta6, alpha2, alpha3, alpha6
## NOTE: pars should be a named list with names beta1, beta6,  alpha2, alpha3 and alpha6. beta1, beta6 should be on a simplex 
# a corresponds to alpha2 + alpha3. dat is NULL for nu2 as there is no data corresponding to nu2
Target_nu2 <- function(ycand, ycurrent, a, dat = NULL, pars) {
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha2 = as.numeric(pars$alpha2); alpha3 = as.numeric(pars$alpha3); alpha6 = as.numeric(pars$alpha6) # extract the alpha parameters
  beta1 = pars$beta1; beta6 = pars$beta6 # extract the weight parameters. 
  L = length(beta1)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Numerator in the MH Ratio in the log scale
  proposed = sum(((a * beta1)  - 1)*log(ycand) + (alpha6 * ycand) * log(beta6) - lgamma(alpha6 * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum(((a * beta1)  - 1)*log(ycurrent) + (alpha6 * ycurrent) * log(beta6) - lgamma(alpha6 * ycurrent))
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 

# Define the target function for [nu3| - ] in log scale with Metropolis ratio
# @pars includes beta1, beta7, alpha3, alpha4, alpha7
## NOTE: pars should be a named list with names beta1, beta7,  alpha3, alpha4 and alpha7. beta1, beta7 should be on a simplex 
# a corresponds to alpha3 + alpha4. dat is NULL for nu3 as there is no data corresponding to nu3
Target_nu3 <- function(ycand, ycurrent, a, dat = NULL, pars) {
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha3 = as.numeric(pars$alpha3); alpha4 = as.numeric(pars$alpha4); alpha7 = as.numeric(pars$alpha7) # extract the alpha parameters
  beta1 = pars$beta1; beta7 = pars$beta7 # extract the weight parameters. 
  L = length(beta1)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Numerator in the MH Ratio in the log scale
  proposed = sum(((a * beta1)  - 1)*log(ycand) + (alpha7 * ycand) * log(beta7) - lgamma(alpha7 * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum(((a * beta1)  - 1)*log(ycurrent) + (alpha7 * ycurrent) * log(beta7) - lgamma(alpha7 * ycurrent))
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 

# Define the target function for [nu4| - ] in log scale with Metropolis ratio
# @pars includes beta8, eta, alpha5, alpha6, alpha7, alpha8
## NOTE: pars should be a named list with names beta8, eta,  alpha5, alpha6, alpha7 and alpha8. beta8, eta should be on a simplex 
# a corresponds to alpha5 + alpha6 + alpha7. dat is NULL for nu4 as there is no data corresponding to nu4
Target_nu4 <- function(ycand, ycurrent, a, dat = NULL, pars) {
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha5 = as.numeric(pars$alpha5); alpha6 = as.numeric(pars$alpha6); alpha7 = as.numeric(pars$alpha7); alpha8 = as.numeric(pars$alpha8) # extract the alpha parameters
  beta8 = pars$beta8; eta = pars$eta # extract the weight parameters.
  L = length(beta8)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Numerator in the MH Ratio in the log scale
  proposed = sum(((a * eta)  - 1)*log(ycand) + (alpha8 * ycand) * log(beta8) - lgamma(alpha8 * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum(((a * eta)  - 1)*log(ycurrent) + (alpha8 * ycurrent) * log(beta8) - lgamma(alpha8 * ycurrent))
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 

# Define the target function for [eta| - ] in log scale with Metropolis ratio
# @pars includes beta1, nu4, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7
## NOTE: pars should be a named list with names beta1, nu4,  alpha2, alpha3, alpha4, alpha5, alpha6, alpha7. beta1, nu4 should be on a simplex 
# a corresponds to (alpha2 + alpha3 + alpha4). dat is NULL for eta as there is no data corresponding to eta
Target_eta <- function(ycand, ycurrent, a, dat = NULL, pars) {
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha2 = as.numeric(pars$alpha2); alpha3 = as.numeric(pars$alpha3); alpha4 = as.numeric(pars$alpha4)
  alpha5 = as.numeric(pars$alpha5); alpha6 = as.numeric(pars$alpha6); alpha7 = as.numeric(pars$alpha7) # extract the alpha parameters
  beta1 = pars$beta1; nu4 = pars$nu4 # extract the weight parameters. 
  L = length(beta1)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Numerator in the MH Ratio in the log scale
  proposed = sum(((2* a * beta1)  - 1)*log(ycand) + ((alpha5 + alpha6 + alpha7) * ycand) * log(nu4) - lgamma((alpha5 + alpha6 + alpha7)  * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum(((2* a * beta1)  - 1)*log(ycurrent) + ((alpha5 + alpha6 + alpha7) * ycurrent) * log(nu4) - lgamma((alpha5 + alpha6 + alpha7)  * ycurrent))
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 

# Define the target function for [beta1| - ] in log scale with Metropolis ratio
# Include beta2, beta3, beta4 nu1, nu2, nu3, eta, alpha2, alpha3, alpha4 in pars
## NOTE: pars should be a named list with names beta2, beta3, beta4 nu1, nu2, nu3, eta, alpha2, alpha3, alpha4. beta2, beta3, beta4 nu1, nu2, nu3, eta should be on a simplex 
# a correspons to alpha1. dat corresponds to m1
Target_beta1 <- function(ycand, ycurrent, a, dat, pars) {
  L = length(dat)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha2 = as.numeric(pars$alpha2); alpha3 = as.numeric(pars$alpha3); alpha4 = as.numeric(pars$alpha4) # extract the alpha parameters
  beta2 = pars$beta2; beta3 = pars$beta3; beta4 = pars$beta4; nu1 = pars$nu1; nu2 = pars$nu2; nu3 = pars$nu3; eta = pars$eta # extract the weight parameters. 
  # Numerator in the MH Ratio in the log scale
  proposed = sum((dat + a/L - 1)*log(ycand) + (alpha2*ycand)*log(beta2) + (alpha3*ycand)*log(beta3) + (alpha4*ycand)*log(beta4) + ((alpha2 + alpha4)*ycand)*log(nu1) + ((alpha2 + alpha3)*ycand)*log(nu2) + ((alpha3 + alpha4)*ycand)*log(nu3) + (2*(alpha2 + alpha3 + alpha4)*ycand)*log(eta) - lgamma(alpha2 * ycand) - lgamma(alpha3 * ycand) - lgamma(alpha4 * ycand) - lgamma((alpha2 + alpha4) * ycand) - lgamma((alpha2 + alpha3) * ycand) - lgamma((alpha3 + alpha4) * ycand) - lgamma(2*(alpha2 + alpha3 + alpha4) * ycand) )
  # Denominator in the MH ratio in the log scale
  previous = sum((dat + a/L - 1)*log(ycurrent) + (alpha2*ycurrent)*log(beta2) + (alpha3*ycurrent)*log(beta3) + (alpha4*ycurrent)*log(beta4) + ((alpha2 + alpha4)*ycurrent)*log(nu1) + ((alpha2 + alpha3)*ycurrent)*log(nu2) + ((alpha3 + alpha4)*ycurrent)*log(nu3) + (2*(alpha2 + alpha3 + alpha4)*ycurrent)*log(eta) - lgamma(alpha2 * ycurrent) - lgamma(alpha3 * ycurrent) - lgamma(alpha4 * ycurrent) - lgamma((alpha2 + alpha4) * ycurrent) - lgamma((alpha2 + alpha3) * ycurrent) - lgamma((alpha3 + alpha4) * ycurrent) - lgamma(2*(alpha2 + alpha3 + alpha4) * ycurrent) )
  
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 


