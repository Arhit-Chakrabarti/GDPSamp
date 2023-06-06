sample_alpha <- function(ycurrent, pars){
  # Extract the parameters from the pars dataframe
  beta = as.numeric(pars$beta)
  alpha0 = pars$alpha0
  L = length(beta) # Length of the simplex
  
  ycand = stats::rgamma(n = 1, shape = alpha0) # Propose a candidate
  # Numerator of MH ratio in log-scale
  num = lgamma(ycand) + sum(ycand/L * log(beta)) - L*lgamma(ycand/L) - ycand + (alpha0 - 1)*log(ycand) + dgamma(ycurrent, shape = alpha0, log = TRUE)
  # Denominator of MH ratio in log-scale
  den = lgamma(ycurrent) + sum(ycurrent/L * log(beta)) - L*lgamma(ycurrent/L) - ycurrent + (alpha0 - 1)*log(ycurrent) + dgamma(ycand, shape = alpha0, log = TRUE)
  
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

# @prior is a named list with alpha0, a, b prior_mean, prior_prec
# @init is a named list with beta.start, mu.start, tau.start, z.start, alpha.start
#' @export
gibbs_MVNW_dp_c <- function(init, data, num_iterations, burn, thin = 1, L.max, prior){
  # Define the lists to store the samples
  beta_samples <- list() # list to store samples of weights beta
  mu_samples <- list() # list to store Normal means samples 
  tau_samples <- list() # list to store Normal precision samples 
  z_samples <- list() # list to store Normal latent indicator samples 
  n_samples <- list() # list to store latent cluster sizes
  alpha_samples <- list() # list to store the concentration parameters alpha
  
  nu = prior$nu; W = prior$W
  
  # Initialize the starting values of the different parameters
  beta_samples[[1]] <- init$beta.start
  mu_samples[[1]] <- init$mu.start 
  tau_samples[[1]] <- init$tau.start 
  z_samples[[1]] <- init$z.start 
  alpha_samples[[1]] <- init$alpha.start
  
  for(i in 2:num_iterations){
    # Printing the iterations
    if(i == 2){
      cat(paste0("Iteration: ", (i-1), "\n"))
    }
    if(i %% floor((10/100)*(num_iterations + 1)) == 0) {
      cat(paste0("Iteration: ", i, "\n"))
    }
    # Sample the latent cluster indicators and matrix of probabilities
    z_samples[[i]] = sample_z_c(data = t(data), beta = beta_samples[[i-1]], mu = mu_samples[[i-1]], tau = tau_samples[[i-1]], L = L.max) # Take out the latent cluster indicators
    
    # Sample the means of the normal dsitribution
    # sample_muC_tau <- function(data, z, L, nu, W, prior_mean, prior_prec)
    mu_tau_sampling = sample_mu_tau_c(data = t(data), z = z_samples[[i]], L = L.max, prior_mean = prior$prior_mean, prior_prec = prior$prior_prec, nu = nu, W = W)
    
    mu_samples[[i]] =  t(mu_tau_sampling$mu.draw) # Means
    tau_samples[[i]] =  mu_tau_sampling$tau.draw # Precisions
    n_samples[[i]] =  mu_tau_sampling$n_h # Cluster sizes
    
    # Sample the weights
    beta_samples[[i]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n_samples[[i]] + alpha_samples[[i-1]]/L.max))
    beta_samples[[i]] = beta_samples[[i]] + 1e-6
    beta_samples[[i]] = beta_samples[[i]]/sum(beta_samples[[i]])
    # Sample the concentration parameter alpha's
    alpha_samples[[i]] = sample_alpha(ycurrent = alpha_samples[[i-1]], 
                                      pars = list(alpha0 = prior$alpha0,
                                                  beta = beta_samples[[i]]) )
    
  }# End of Gibbs update
  
  thin_samples <- seq(from = (burn + 1), to = num_iterations, by = thin)
  
  return(list(n_samples = n_samples[thin_samples],
              z_samples = z_samples[thin_samples],
              alpha_samples = alpha_samples[thin_samples],
              beta_samp_post = beta_samples[thin_samples]))
  
}

# @ z.list is the list of indicator samples post burnin
#' @export
getDahl <- function(z.list, z.true){
  membershipMatrices <- lapply(z.list, function(x){
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
  })
  
  membershipAverage <- Reduce("+", membershipMatrices)/length(z.list)
  
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  
  DahlIndex <- which.min(SqError) # Which index post burn-in which gives minimun Dahl index
  
  
  if(is.null(z.true)){
    return(list(DahlIndex = DahlIndex))
  }else{
  return(list(rand = fossil::rand.index(z.list[[DahlIndex]], z.true), # Rand Index by Vanilla DP
              DahlIndex = DahlIndex))
  }
}
