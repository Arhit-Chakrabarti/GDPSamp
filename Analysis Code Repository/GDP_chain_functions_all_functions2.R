if(!require("Rcpp")) install.packages("Rcpp")
if(!require("fossil")) install.packages("fossil")

sample_alpha <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta = as.numeric(pars$beta)
  a0 = pars$a0
  b0 = pars$b0
  L = length(beta) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  # num = lgamma(ycand) + sum(ycand/L * log(beta)) - L*lgamma(ycand/L) - ycand + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  num = lgamma(ycand) + sum(ycand/L * log(beta)) - L*lgamma(ycand/L) + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  # Denominator of MH ratio in log-scale
  # den = lgamma(ycurrent) + sum(ycurrent/L * log(beta)) - L*lgamma(ycurrent/L) - ycurrent + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  den = lgamma(ycurrent) + sum(ycurrent/L * log(beta)) - L*lgamma(ycurrent/L) + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(out)
}

library(Rcpp)
sourceCpp("GDP_chain_functions_c.cpp")

gibbs_MVN_dp_c <- function(init, data, num_iterations, burn, thin = 1, L.max, prior){
  # Define the lists to store the samples
  beta_samples <- list() # list to store samples of weights beta
  mu_samples <- list() # list to store Normal means samples 
  z_samples <- list() # list to store Normal latent indicator samples 
  n_samples <- list() # list to store latent cluster sizes
  probability_samples <- list() # list to store matrix of probabilities
  alpha_samples <- list() # list to store the concentration parameters alpha
  
  # Initialize the starting values of the different parameters
  beta_samples[[1]] <- init$beta.start
  mu_samples[[1]] <- init$mu.start 
  z_samples[[1]] <- init$z.start 
  alpha_samples[[1]] <- init$alpha.start
  
  prior_mean = prior$prior_mean
  prior_sigma = prior_sigma
  
  for(i in 2:num_iterations){
    # i = 2
    # Printing the iterations
    if(i == 2){
      cat(paste0("Iteration: ", (i-1), "\n"))
    }
    if(i %% floor((10/100)*(num_iterations + 1)) == 0) {
      cat(paste0("Iteration: ", i, "\n"))
    }
    # Sample the latent cluster indicators and matrix of probabilities
    z_samples[[i]] = sample_z_known_sigma_c(data = t(data), beta = beta_samples[[i-1]], mu = mu_samples[[i-1]], Sigma = Sigma.true, L = L.max) # Take out the latent cluster indicators
    
    n = ncol(data)
    # Sample the means of the normal dsitribution
    # sample_muC_tau <- function(data, z, L, nu, W, prior_mean, prior_prec)
    out = sample_mu_dp_c(data = data, z = z_samples[[i]], Sigma = Sigma.true, L = L.max, n = n, prior_mean = prior_mean,  prior_sigma = prior_sigma)
    
    mu_samples[[i]] =  t(out$mu_draw) # Means
    n_samples[[i]] =  out$n_all # Cluster sizes
    
    # Sample the weights
    beta_samples[[i]] = as.numeric(rdirichlet(n = 1, alpha = n_samples[[i]] + alpha_samples[[i-1]]/L.max))
    
    beta_samples[[i]] = beta_samples[[i]] + 1e-8 
    beta_samples[[i]] = beta_samples[[i]]/sum(beta_samples[[i]]) 
    
    # Sample the concentration parameter alpha's
    alpha_samples[[i]] = sample_alpha(ycurrent = alpha_samples[[i-1]], 
                                      pars = list(a0 = prior$a0,
                                                  b0 = prior$b0,
                                                  beta = beta_samples[[i]]) )
    
  }# End of Gibbs update
  
  thin_samples <- seq(from = (burn + 1), to = num_iterations, by = thin)
  
  return(list(n_samples = n_samples[thin_samples],
              z_samples = z_samples[thin_samples],
              alpha_samples = alpha_samples[thin_samples],
              beta_samp_post = beta_samples[thin_samples]))
}


# @ z.list is the list of indicator samples post burnin
getDahl <- function(z.list, z.true){
  membershipMatrices <- lapply(z.list, function(x){
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
  })
  
  membershipAverage <- Reduce("+", membershipMatrices)/length(z.list)
  
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  
  DahlIndex <- which.min(SqError) # Which index post burn-in which gives minimun Dahl index
  
  library(fossil)
  library(aricode)
  if(is.null(z.true)){
    return(list(DahlIndex = DahlIndex))
  }else{
    return(list(rand = aricode::ARI(z.list[[DahlIndex]], z.true), # Rand Index by Vanilla DP
                DahlIndex = DahlIndex))
  }
}

####################################################################################################
# SIMPLEX TARGETS
####################################################################################################
# Define the target function for [beta1| - ] in log scale with Metropolis ratio
# @pars includes beta2, alpha2
## NOTE: pars should be a named list with names beta2, alpha2. beta2, should be on a simplex 
#  dat is m_1h
# a is alpha1
Target_beta1 <- function(ycand, ycurrent, a, dat, pars) {
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha2 = as.numeric(pars$alpha2) # extract the alpha parameters
  beta2 = pars$beta2# extract the weight parameters. 
  L = length(beta2)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Numerator in the MH Ratio in the log scale
  proposed = sum((dat + a/L - 1)*log(ycand) + (alpha2 * ycand) * log(beta2) - lgamma(alpha2 * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum((dat + a/L - 1)*log(ycurrent) + (alpha2 * ycurrent) * log(beta2) - lgamma(alpha2 * ycurrent))
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 


####################################################################################################
# Define the target function for [beta1| - ] in log scale with Metropolis ratio
# @pars includes beta_j_minus_1, beta_j_plus_1, alpha_j_plus_1
## NOTE: pars should be a named list with names beta_j_minus_1, beta_j_plus_1, alpha_j_plus_1. beta_j_minus_1 and beta_j_plus_1 should be on a simplex 
#  dat is m_jh
# a is alpha_j
Target_betaj <- function(ycand, ycurrent, a, dat, pars) {
  # Extract all the constant parameters from pars. pars should be a names dataframe
  alpha_j_plus_1 = as.numeric(pars$alpha_j_plus_1) # extract the alpha parameters
  beta_j_minus_1 = pars$beta_j_minus_1; beta_j_plus_1 = pars$beta_j_plus_1; # extract the weight parameters. 
  L = length(beta_j_minus_1)
  ycand = rep(0, L); ycurrent = rep(0, L)
  ycand <- exp(LogPq(ycand)$logp); ycurrent <- exp(LogPq(ycurrent)$logp)
  # Numerator in the MH Ratio in the log scale
  proposed = sum((dat + (a * beta_j_minus_1) - 1)*log(ycand) + (alpha_j_plus_1 * ycand) * log(beta_j_plus_1) - lgamma(alpha_j_plus_1 * ycand))
  # Denominator in the MH ratio in the log scale
  previous = sum((dat + (a * beta_j_minus_1) - 1)*log(ycurrent) + (alpha_j_plus_1 * ycurrent) * log(beta_j_plus_1) - lgamma(alpha_j_plus_1 * ycurrent))
  out = proposed - previous # MH Ratio but due to log-scale becomes difference
  return(out)
} 



sample_alpha1 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta1 = as.numeric(pars$beta1)
  a0 = pars$a0
  b0 = pars$b0
  L = length(beta1) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  # num = lgamma(ycand) + sum(ycand/L * log(beta1)) - L*lgamma(ycand/L) - ycand + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  num = lgamma(ycand) + sum(ycand/L * log(beta1)) - L*lgamma(ycand/L) + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  # Denominator of MH ratio in log-scale
  # den = lgamma(ycurrent) + sum(ycurrent/L * log(beta1)) - L*lgamma(ycurrent/L) - ycurrent + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  den = lgamma(ycurrent) + sum(ycurrent/L * log(beta1)) - L*lgamma(ycurrent/L) + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(out)
}



sample_alphaj <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta_j = as.numeric(pars$beta_j)
  beta_j_minus_1 = as.numeric(pars$beta_j_minus_1)
  a0 = pars$a0
  b0 = pars$b0
  L = length(beta_j) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  # num = lgamma(ycand) + sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) - ycand + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  
  num = lgamma(ycand) + sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  # Denominator of MH ratio in log-scale
  # den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) - ycurrent + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(out)
}


# sample_alphaJ <- function(ycurrent, pars){
#   # Extract the parameters from the pars dataframe
#   beta_j = as.numeric(pars$beta_j)
#   beta_j_minus_1 = as.numeric(pars$beta_j_minus_1)
#   alpha_j_minus_1 = as.numeric(pars$alpha_j_minus_1)
#   a0 = pars$a0
#   b0 = pars$b0
#   L = length(beta_j) # Length of the simplex
#   
#   ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
#   # Numerator of MH ratio in log-scale
#   # num = lgamma(ycand) + sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) - ycand + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
#   
#   num = lgamma(ycand) + sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycand) - (ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
#   # Denominator of MH ratio in log-scale
#   # den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) - ycurrent + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
#   den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycurrent) - (ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
#   
#   MH = num - den # MH ratio which become difference in log-scale
#   if(log(runif(n = 1)) < MH){ 
#     out = ycand # Accept proposed value if log(U) < MH
#     accept = 1
#   }else{
#     out = ycurrent # Reject proposed value if log(U) < MH
#     accept = 0
#   }
#   return(out)
# }

sample_alphaJ <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta_j = as.numeric(pars$beta_j)
  beta_j_minus_1 = as.numeric(pars$beta_j_minus_1)
  alpha_j_minus_1 = as.numeric(pars$alpha_j_minus_1)
  a0 = pars$a0
  b0 = pars$b0
  L = length(beta_j) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) + sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) - ycand + (a0 - 1)*log(ycand) - (b0 * ycand) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  
  # num = lgamma(ycand) + sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycand) - (ycand) + dgamma(ycurrent, shape = alpha_j_minus_1, rate = 1, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) - ycurrent + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  
  # den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycurrent) - (ycurrent) + dgamma(ycand, shape = alpha_j_minus_1, rate = 1, log = TRUE)
  
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(out)
}

sample_alphaj_j_minus_1 <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta_j = as.numeric(pars$beta_j)
  beta_j_minus_1 = as.numeric(pars$beta_j_minus_1)
  alpha_j_minus_1 = as.numeric(pars$alpha_j_minus_1)
  alpha_j_plus_1 = as.numeric(pars$alpha_j_plus_1)

  a0 = pars$a0
  b0 = pars$b0
  
  L = length(beta_j) # Length of the simplex
  
  ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate

  num = sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycand) - (ycand) + (ycand * log(alpha_j_plus_1)) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
  
  # num = sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycand) - (ycand) + (ycand * log(alpha_j_plus_1)) + dgamma(ycurrent, shape = alpha_j_minus_1, rate = 1, log = TRUE)
  
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) - ycurrent + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
  
  # den = sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycurrent) - (ycurrent) + (ycurrent * log(alpha_j_plus_1))  + dgamma(ycand, shape = alpha_j_minus_1, rate = 1, log = TRUE)
  
  MH = num - den # MH ratio which become difference in log-scale
  if(log(runif(n = 1)) < MH){ 
    out = ycand # Accept proposed value if log(U) < MH
    accept = 1
  }else{
    out = ycurrent # Reject proposed value if log(U) < MH
    accept = 0
  }
  return(out)
}

# sample_alphaj_j_minus_1 <- function(ycurrent, pars){
#   # Extract the parameters from the pars dataframe
#   beta_j = as.numeric(pars$beta_j)
#   beta_j_minus_1 = as.numeric(pars$beta_j_minus_1)
#   alpha_j_minus_1 = as.numeric(pars$alpha_j_minus_1)
#   alpha_j_plus_1 = as.numeric(pars$alpha_j_plus_1)
#   
#   a0 = pars$a0
#   b0 = pars$b0
#   
#   L = length(beta_j) # Length of the simplex
#   
#   ycand = rgamma(n = 1, shape = a0, rate = b0) # Propose a candidate
#   
#   num = sum(ycand * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycand * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycand) - (ycand) + (ycand * log(alpha_j_plus_1)) + dgamma(ycurrent, shape = a0, rate = b0, log = TRUE)
#   # Denominator of MH ratio in log-scale
#   # den = lgamma(ycurrent) + sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) - ycurrent + (a0 - 1)*log(ycurrent) - (b0 * ycurrent) + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
#   den = sum(ycurrent * beta_j_minus_1 * log(beta_j)) - sum(lgamma(ycurrent * beta_j_minus_1)) + (alpha_j_minus_1 - 1)*log(ycurrent) - (ycurrent) + (ycurrent * log(alpha_j_plus_1))  + dgamma(ycand, shape = a0, rate = b0, log = TRUE)
#   
#   MH = num - den # MH ratio which become difference in log-scale
#   if(log(runif(n = 1)) < MH){ 
#     out = ycand # Accept proposed value if log(U) < MH
#     accept = 1
#   }else{
#     out = ycurrent # Reject proposed value if log(U) < MH
#     accept = 0
#   }
#   return(out)
# }
