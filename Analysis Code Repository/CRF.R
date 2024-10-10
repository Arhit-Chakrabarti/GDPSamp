# Code to implement Chinese restaurant franchise (CRF) based collapsed Gibbs sampler
# for a univariate conjugate Gaussian mixture model.
# Notations follow from Hierarchical Dirichlet Processes (Teh et al., 2006).

# required package
library(mvtnorm)

# Function to get marginal prior of each observation x_ji
# x_ji : observation corresponding to customer i in restaurant j.
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
f_k_new.xji = function(xji, xi, lambda, tau){
  
  res = dnorm(xji, mean = xi, sd = sqrt((1/tau) + (1/lambda)))
  return(res)
}

# Function to get conditional density of x_ji under mixture component k given all data items except x_ji
# x_ji : observation corresponding to customer i in restaurant j.
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
# b : #customers having dish k leaving out x_ji
# xbar_b : mean of x's having dish index k leaving out x_ji
f_k.xji = function(xji, xi, lambda, tau, b, xbar_b){
  r = lambda / ((b*tau)+lambda)
  mu = r * xi + (1-r) * xbar_b
  s = sqrt((1/tau) + (1/ ((b*tau)+lambda)))
  
  res = dnorm(xji, mean = mu, sd = s)
  return(res)
}

# Function to get marginal prior of each observation x_jt
# x_jt : observations corresponding to all customers in table t of restaurant j.
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
f_k_new.xjt = function(xjt, xi, lambda, tau){
  n = length(xjt)
  mu = rep(xi, n)
  S = ((1/tau) * diag(n)) +  matrix((1/lambda), nrow = n, ncol = n)
  
  res = mvtnorm::dmvnorm(x = xjt, mean = mu, sigma = S)
  return(res)
}

# Function to get conditional density of x_jt under mixture component k given all data items except x_jt
# x_jt : observations corresponding to all customers in table t of restaurant j.
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
# d : #customers having dish k leaving out x_jt
# xbar_d : mean of x's having dish index k leaving out x_jt
f_k.xjt = function(xjt, xi, lambda, tau, d, xbar_d){
  a = length(xjt)
  r = lambda / ((d*tau)+lambda)
  
  mu = r * xi + (1-r) * xbar_d
  mu = rep(mu, a)
  
  S = ((1/tau)*diag(a)) +  matrix((1/((d*tau)+lambda)), nrow = a, ncol = a)
  
  res = mvtnorm::dmvnorm(x = xjt, mean = mu, sigma = S)
  return(res)
}

# Function to get unique dish indices Z_ji = k_{jt_ji} for all i,j.
# k : list containing dish indices for each table of each restaurant
# t: list containing table indices for each customer in each restaurant
construct_Z = function(k, t){
  
  J = length(t); n = lengths(t)
  Z = vector(mode = "list", length = J)
  
  for(j in 1:J){
    Z[[j]] = sapply(1:n[j], function(i) k[[j]][t[[j]][i]])
  }
  return(Z)
}

# Function to get unique dish indices Z_ji = k_{jt_ji} for all i,j leaving out (j_,i_).
# k : list containing dish indices for each table of each restaurant
# t: list containing table indices for each customer in each restaurant
# j_, i_ : indices that should be dropped
construct_Z_ji = function(j_, i_, k, t){
  
  J = length(t); n = lengths(t)
  Z = vector(mode = "list", length = J)
  
  for(j in 1:J){
    Z[[j]] = sapply(1:n[j], function(i) k[[j]][t[[j]][i]])
  }
  
  Z[[j_]][i_] = NA
  return(Z)
}

# Function to remove all observations corresponding to table t of restaurant j.
# j,t : restaurant and table that should be removed
# t.list : list containing table indices for each customer in each restaurant
# y : list containing data for all populations
remove_jt = function(j, t, t.list, y){
  y_jt = y
  
  # get the indices corresponding to
  ind.yjt = which(t.list[[j]]== t)
  y_jt[[j]][ind.yjt] = NA
  
  return(y_jt)
}

# Function to sample the table indices
# x : list containing data for all populations
# k.list : list containing dish indices for each table in each restaurant
# t.list : list containing table indices for each customer in each restaurant
# alpha0 : shared concentration parameter of local DPs
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
# gam : concentration parameter of global DP
sample.t = function(x, k.list, t.list, alpha0, xi, lambda, tau, gam){
  
  J = length(x); n = lengths(x); K = max(unlist(k.list))
  
  t.out = t.list
  k.out = k.list
  
  # vector storing the number of tables in each restaurant
  mj0 = lengths(k.list)
  
  for(j in 1:J){
    
    Tj = mj0[j] # number of tables in restaurant j
    
    for(i in 1:n[j]){
      
      tji = t.out[[j]][i]
      
      Zji = k.out[[j]][tji]
      
      # number of customers in table t at restaurant j leaving out x_ji
      njt0_ji = sapply(seq_len(Tj), function(t) sum(t.out[[j]] == t))
      njt0_ji[tji] = njt0_ji[tji] - 1
      
      # if no other customer is present in table t_ji, then the table is dropped and 
      # number of tables in restaurant j is decreased by 1
      if(njt0_ji[tji] == 0){
        
        t.out[[j]][i] = NA
        t.ind = which(t.out[[j]] == Tj)
        t.out[[j]][t.ind] = tji
        
        k.out[[j]][tji] = k.out[[j]][Tj]
        k.out[[j]] = k.out[[j]][-Tj]
        
        njt0_ji[tji] = njt0_ji[Tj]
        njt0_ji = njt0_ji[-Tj]
        Tj = Tj - 1
      }
      
      # get the unique dish indices for every customer, Z_ji = k_{jt_ji}
      Z = construct_Z_ji(j, i, k.out, t.out)
      Z_ji.vec = unlist(Z)
      
      xji = x[[j]][i]
      x_ji = x; x_ji[[j]][i] = NA; x_ji = unlist(x_ji)
      
      # calculate the number of dishes and mean of observations leaving out x_ji
      n.dish = sapply(seq_len(K), function(k) sum(Z_ji.vec == k, na.rm = TRUE))
      xbar = sapply(seq_len(K), function(k) mean(x_ji[Z_ji.vec == k], na.rm = TRUE))
      
      # if no other table serves dish k_jt_ji, then the dish is dropped
      # and total number of dishes is decreased by 1
      if(n.dish[Zji] == 0){
        
        n.dish[Zji] = n.dish[K]
        xbar[Zji] = xbar[K]
        
        n.dish = n.dish[-K]
        xbar = xbar[-K]
        
        for(l in 1:J){
          k.ind = which(k.out[[l]] == K)
          k.out[[l]][k.ind] = Zji
          
          Z.ind = which(Z[[l]] == K)
          Z[[l]][Z.ind] = Zji
        }
        K = K-1
      }
      
      # update table indices according to the CRF scheme
      f_ji = sapply(seq_len(K), function(k) f_k.xji(xji = xji, xi = xi, lambda = lambda, tau = tau,
                                                    b = n.dish[k], xbar_b = xbar[k]))
      kj = k.out[[j]]
      
      f.kjt_ji = sapply(seq_len(Tj), function(t) f_ji[kj[t]])
      prob_t = log(njt0_ji) + log(f.kjt_ji)
      
      f_ji[K+1] = f_k_new.xji(xji = xji, xi = xi, lambda = lambda, tau = tau)
      
      m0k = sapply(seq_len(K), function(k) sum(unlist(k.out) == k, na.rm = TRUE))
      m_vec = c(m0k, gam)
      m_ratio = m_vec / sum(m_vec)
      
      lik_t_new = sum(m_ratio * f_ji)
      
      prob_t[Tj+1] = log(alpha0) + log(lik_t_new)
      
      prob_t  = prob_t - max(prob_t)
      prob_t = exp(prob_t)
      prob_t = prob_t/ sum(prob_t)
      
      t.out[[j]][i] = sample(seq_len(Tj+1), size = 1, prob = prob_t)
      
      if(t.out[[j]][i] == Tj+1){
        
        # if a new table is sampled, increase its count and sample a dish for that table
        Tj = Tj + 1
        prob_k = log(m_vec) + log(f_ji)
        
        prob_k  = prob_k - max(prob_k)
        prob_k = exp(prob_k)
        prob_k = prob_k/ sum(prob_k)
        
        samp.k_new = sample(seq_len(K+1), size = 1, prob = prob_k)
        
        # if a new dish is sampled, increase the total dish count
        if(samp.k_new == (K+1)){
          K = K+1
        }
        k.out[[j]][Tj] = samp.k_new
        
      }
    }
  }
  
  return(list("k.list" = k.out, "t.list" = t.out))
}

# Function to sample the dish indices
# x : list containing data for all populations
# k.list : list containing dish indices for each table in each restaurant
# t.list : list containing table indices for each customer in each restaurant
# alpha0 : shared concentration parameter of local DPs
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
# gam : concentration parameter of global DP
sample.k = function(x, k.list, t.list, alpha0, xi, lambda, tau, gam){
  
  J = length(x); K = max(unlist(k.list))
  
  k.out = k.list
  
  # vector storing the number of tables in each restaurant
  mj0 = lengths(k.list)
  n.tables = mj0
  
  for(j in 1:J){
    for(t in 1:n.tables[j]){
      
      # get the unique dish indices for every customer, Z_ji = k_{jt_ji}
      Z = construct_Z(k.out, t.list)
      
      # remove dish corresponding to table t in restaurant j
      kjt = k.out[[j]][t]
      k.out[[j]][t] = NA
      
      # extract observations corresponding to table t in restaurant j
      xjt = x[[j]][ t.list[[j]] == t]
      
      # remove x_jt from the list of data and cluster labels
      x_jt = unlist(remove_jt(j, t, t.list, x))
      Z_jt = unlist(remove_jt(j, t, t.list, Z))
      
      # calculate the number of dishes and mean of observations leaving out x_jt 
      n.dish = sapply(seq_len(K), function(k) sum(Z_jt == k, na.rm = TRUE))
      xbar = sapply(seq_len(K), function(k) mean(x_jt[Z_jt == k], na.rm = TRUE))
      
      # if no other table serves dish k_jt, 
      # then drop that dish and decrease total number of dishes K by 1
      if(n.dish[kjt] == 0){
        n.dish[kjt] = n.dish[K]
        xbar[kjt] = xbar[K]
        
        n.dish = n.dish[-K]
        xbar = xbar[-K]
        
        for(l in 1:J){
          k.ind = which(k.out[[l]]==K)
          k.out[[l]][k.ind] = kjt
        }
        K = K-1
      }
      
      # update dish indices according to the CRF scheme
      f.k_jt = sapply(seq_len(K), function(k) f_k.xjt(xjt = xjt, xi = xi, lambda = lambda, tau = tau,
                                                      d = n.dish[k], xbar_d = xbar[k]))
      f.k_jt[K+1] = f_k_new.xjt(xjt = xjt, xi = xi, lambda = lambda, tau = tau)
      
      m0k_jt = sapply(seq_len(K), function(k) sum(unlist(k.out) == k, na.rm = TRUE))
      m_vec = c(m0k_jt, gam)
      
      prob_k = log(m_vec) + log(f.k_jt) 
      prob_k  = prob_k - max(prob_k)
      prob_k = exp(prob_k)
      prob_k = prob_k/ sum(prob_k)
      
      samp.k = sample(seq_len(K+1), size = 1, prob = prob_k)
      k.out[[j]][t] = samp.k
      
      # if a new dish gets sampled, increase number of dishes by 1
      if(samp.k == (K+1)){
        K = K + 1
      }
    }
  }
  return(k.out)
  
}

# Function to sample alpha0 using the auxiliary variable scheme given in Appendix A of Teh et al.(2006).
# k.list : list containing dish indices for each table in each restaurant
# t.list : list containing table indices for each customer in each restaurant
# alpha0 : shared concentration parameter of local DPs
# a : shape parameter of the gamma prior on alpha0
# b : rate parameter of the gamma prior on alpha0
sample.alpha0 = function(k.list, t.list, alpha0, a, b){
  J = length(t.list); n = lengths(t.list)
  m = sum(lengths(k.list))
  
  # update auxiliary Beta random variables w 
  w = sapply(seq_len(J), function(j) rbeta(1, shape1 = (alpha0+1), shape2 = n[j]))
  
  # update auxiliary Bernoulli random variables s 
  p = n/(alpha0 + n)
  s = sapply(p, function(j) rbinom(1, size = 1, prob = j))
  
  # update alpha0 given w and s
  shape_alpha0 = a + m - sum(s)
  rate_alpha0 = b - sum(log(w))
  samp.alpha0 = rgamma(1, shape = shape_alpha0, rate = rate_alpha0)
  
  return(samp.alpha0)
  
}

# Function to sample the atoms corresponding to the occupied clustes
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
# x : list containing data for all populations
# Z :  list containing global cluster labels
sample.phi = function(xi, lambda, tau, x, Z){
  L = max(unlist(Z))
  
  # unlisting all observations and class labels
  x.all = unlist(x)
  Z.all = unlist(Z)
  
  # n = (n1,...,nL) where n_k = \sum_(j,i) I(Z_ji = k)
  n.group = sapply(seq_len(L), function(k) sum(Z.all == k))
  
  # xbar = (xbar_1,..., xbar_L), xbar_k =(1/n_k) * \sum_{(j,i): Z_ji = k} x_ji
  xbar = sapply(seq_len(L), function(k) mean(x.all[Z.all == k]))
  
  # posterior mean and standard deviation of phi
  mean_phi = (n.group*tau*xbar + xi*lambda)/((n.group*tau) + lambda)
  sd_phi = 1/ sqrt((n.group*tau) + lambda)
  
  phi = sapply(seq_len(L), function(k) rnorm(1, mean = mean_phi[k], sd = sd_phi[k]))
  
  return(phi)
  
}

# Function to get density estimates for each population using posterior samples from CRF scheme
# y.grid : grid points for estimating the density of the populations
# k.list : list containing dish indices for each table in each restaurant
# t.list : list containing table indices for each customer in each restaurant
# alpha0 : shared concentration parameter of local DPs
# phi : sampled atoms corresponding to each occupied cluster
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
get_density = function(y.grid, k.list, t.list, alpha0, phi, gam, xi, lambda, tau){
  
  J = length(t.list);  n = lengths(t.list); K = length(phi)
  
  # number of tables in restaurant j
  mj0 = lengths(k.list)
  # total number of tables occupied
  m00 = sum(mj0)
  # number of tables serving dish k
  m0k = sapply(seq_len(K), function(k) sum(unlist(k.list) == k))
  
  n.grid = length(y.grid)
  
  # J x n.grid matrix to store the densities along the rows.
  f = matrix(NA, nrow = J, ncol = n.grid)
  
  for(j in 1:J){
    
    # number of tables in restaurant j
    Tj = mj0[j]
    
    # number of customers in table t at restaurant j
    njt0 = sapply(seq_len(Tj), function(t) sum(t.list[[j]] == t))
    
    # evaluate the density using posterior samples from CRF scheme
    for(i in 1:n.grid){
      res1 = (alpha0/(alpha0 + n[j])) * (gam/(gam + m00)) * 
        dnorm(y.grid[i], mean = xi, sd = sqrt((1/tau) + (1/lambda)))
      
      res2 = (alpha0/(alpha0 + n[j])) * (1/(gam + m00)) *
        sum(m0k * dnorm(y.grid[i], mean = phi, sd = sqrt(1/tau)))
      
      phi_3 = sapply(seq_len(Tj), function(t) phi[k.list[[j]][t]])
      res3 = (1/(alpha0 + n[j])) * sum(njt0 * dnorm(y.grid[i], mean = phi_3, sd = sqrt(1/tau)))
      
      f[j, i] = res1 + res2 + res3
    }
  }
  return(f)
}

# CRF based collapsed Gibbs sampler
# x : list of length J, x[[j]] contains data in the jth population
# y.grid : grid points for estimating the density of the populations
# K.init : initial number of global clusters
# gam : concentration parameter of global DP
# a : shape parameter of the gamma prior on alpha0
# b : rate parameter of the gamma prior on alpha0
# xi : prior mean of phi
# lambda : prior precision of phi
# tau : common precision of the normal likelihood
# Burn.in : Burn in period of the MCMC chain
# M : number of MCMC iterations required
CRF_gibbs = function(x, y.grid, K.init, gam, a, b, xi, lambda, tau, Burn.in, M, thin){
  
  J = length(x)
  n = lengths(x)
  
  # set initial values for running gibbs sampler
  alpha0 = rgamma(1, shape = a, rate = b)
  
  k.list = vector(mode = "list", length = J)
  t.list = vector(mode = "list", length = J)
  for(j in 1:J){
    
    # randomly allocate table and dish indices with K.init as initial #dishes.
    if(n[j] >= K.init){
      k.list[[j]] = sample(1:K.init)
      t.list[[j]] = sample(rep(1:K.init, length.out = n[j]))
    }
    else {
      k.list[[j]] = sample(1:n[j])
      t.list[[j]] = sample(1:n[j])
    }
  }
  
  # list to store the posterior samples
  Iterates = vector(mode = "list", length = M)
  
  for(m in 1:M){
    
    # time at the beginning
    T1 = Sys.time()
    
    # update table indices
    res.t = sample.t(x = x, k.list = k.list, t.list = t.list, alpha0 = alpha0, 
                     xi = xi, lambda = lambda, tau = tau, gam = gam)
    
    k.list = res.t$k.list
    t.list = res.t$t.list
    
    # update dish indices
    k.list = sample.k(x = x, k.list = k.list, t.list = t.list, alpha0 = alpha0, 
                      xi = xi, lambda = lambda, tau = tau, gam = gam)
    
    # update alpha0
    alpha0 = sample.alpha0(k.list = k.list, t.list = t.list, alpha0 = alpha0, 
                           a = a, b = b)
    
    # get global cluster allocations
    Z = construct_Z(k = k.list, t = t.list)
    
    # update phi
    phi = sample.phi(xi = xi, lambda = lambda, tau = tau, x = x, Z = Z)
    
    # time at the end of all updates
    T2 = Sys.time()
    Tdiff =  difftime(T2, T1, units = "secs")
    
    # estimate the density for each population
    f = get_density(y.grid = y.grid, k.list = k.list, t.list = t.list, alpha0 = alpha0, 
                    phi = phi, gam = gam, xi = xi, lambda = lambda, tau = tau)
    
    # print every 200th iteration
    if(m %% 200 == 0){
      print(paste("iteration :", m))
    }
    # COMMENT THIS OUT TO STORE THINNED SAMPLES
    # # store samples after Burn in
    # if(m > Burn.in){
    #   Iterates[[m-Burn.in]] = list("k" = k.list, "t" = t.list, "Z" = Z, "phi" = phi, 
    #                                "alpha0" = alpha0, "density" = f, "time" = Tdiff)
    # }
    # store samples after Burn in
    
    Iterates[[m]] = list("k" = k.list, "t" = t.list, "Z" = Z, "phi" = phi,
                                   "alpha0" = alpha0, "density" = f, "time" = Tdiff)
    
  }
  # ADDED THIS PART TO STORE THINNED SAMPLES
  thin_samples <- seq(from = (Burn.in + 1), to = M, by = thin)
  Iterates_thinned = Iterates[thin_samples]
  return(Iterates_thinned)
  
}
