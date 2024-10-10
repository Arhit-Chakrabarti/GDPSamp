# function to get posterior estimate of global cluster assignments 
# by pooling the data across populations using Dahl's method
getDahl = function(MCMCSample){
  
  niters = length(MCMCSample)  #number of MCMC Samples
  
  # extract posterior samples of cluster labels
  iters = lapply(seq_len(niters), function(i) unlist(MCMCSample[[i]]$Z))
  
  # calculate the membership matrices
  membershipMatrices = lapply(iters, function(x){
    clusterAssign = x
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  
  # average membership matrix
  membershipAverage = Reduce("+", membershipMatrices)/niters
  
  # evaluate the sum of squared errors
  SqError = sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                   av = membershipAverage)
  
  # extract the cluster assignment that minimizes the squared errors
  DahlIndex = which.min(SqError)
  DahlAns = MCMCSample[[DahlIndex]]$Z
  attr(DahlAns, "iterIndex") = DahlIndex
  DahlAns
  
  return(DahlAns)
}

# function to calculate the log likelihood of the normal-mixture model
# x: data
# Z: cluster labels
# phi: distinct cluster means
# tau: common precision of the normal distribution
log.lik.normal = function(x, Z, phi, tau = 1){
  M = length(phi)
  n = length(x)
  
  log.lik = rep(NA, M)
  for(m in 1:M){
    
    d = sapply(seq_len(n), function(i) 
      dnorm(x[i], mean = phi[[m]][(Z[[m]][i])], sd = sqrt(1/tau), log = TRUE))
    
    log.lik[m] = sum(d)
  }
  return(log.lik)
}


