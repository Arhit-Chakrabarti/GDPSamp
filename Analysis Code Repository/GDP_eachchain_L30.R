rm(list = ls())
# Read the UMAP Data for different populations
X1.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP1_new4.rds")
X2.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP2_new4.rds")
X3.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP3_new4.rds")
X4.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP4_new4.rds")
X5.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP5_new4.rds")
X6.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP6_new4.rds")
X7.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP7_new4.rds")
X8.UMAP = readRDS(file = "~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/UMAP/Data/UMAP8_new4.rds")

# Change the observations according to the groups from the information provided by Ellen Ruth
X1 <- t(X3.UMAP)
X2 <- t(X1.UMAP)
X3 <- t(X4.UMAP)
X4 <- t(X7.UMAP)
X5 <- t(X5.UMAP)
X6 <- t(X2.UMAP)
X7 <- t(X8.UMAP)
X8 <- t(X6.UMAP)

################################################################################
# Run a Vanilla DP on each group
################################################################################
source("VDP_Gibbs.R")
if(!require("extraDistr")) install.packages("extraDistr"); library(extraDistr)
if(!require("MASS")) install.packages("MASS"); library(MASS)
if(!require("mvtnorm")) install.packages("mvtnorm"); library(mvtnorm)
if(!require("fossil")) install.packages("fossil"); library(fossil)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
################################################################################
# Prior specifications
################################################################################
alpha0 = 3
L.max = 30 # Truncation level of DP
alpha1.start = rgamma(n = 1, shape = alpha0) # Prior for Dirichlet distribution for group 1
alpha2.start = rgamma(n = 1, shape = alpha1.start) # Prior for Dirichlet distribution for group 2
alpha3.start = rgamma(n = 1, shape = alpha1.start) # Prior for Dirichlet distribution for group 3
alpha4.start = rgamma(n = 1, shape = alpha1.start) # Prior for Dirichlet distribution for group 4
alpha5.start = rgamma(n = 1, shape = (alpha2.start + alpha4.start)) # Prior for Dirichlet distribution for group 5
alpha6.start = rgamma(n = 1, shape = (alpha2.start + alpha3.start)) # Prior for Dirichlet distribution for group 6
alpha7.start = rgamma(n = 1, shape = (alpha3.start + alpha4.start)) # Prior for Dirichlet distribution for group 7
alpha8.start = rgamma(n = 1, shape = (alpha5.start + alpha5.start + alpha7.start)) # Prior for Dirichlet distribution for group 8
################################################################################
# SPECIFICATION
################################################################################
nu = nrow(X1)
W = diag(1, nrow = nu) # Rate of Gamma prior

mu0 = as.matrix(rep(0, nu), ncol = 1) # Prior mean of Normal
lambda = 0.01 # Prior precision of Normal
mu.start.list = list()
tau.start = rWishart(n = L.max, df = nu, Sigma = W)
for(i in 1:L.max){
  mu.start.list[[i]] = matrix(mvtnorm::rmvnorm(n = 1, mean = as.numeric(mu0), sigma = solve(lambda * tau.start[ , , i])), ncol = 1)
}

n1 = ncol(X1)
beta.start = rep(1/L.max, L.max)
z.start = rcat(n = n1, prob = beta.start)
alpha.start = rgamma(n = 1, shape = alpha0)
mu.start <- t(sapply(1:L.max, function(j){mu.start.list[[j]]}))

num_iterations <- 50000
burn = 35000

# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
# Initial values of different parameters
init = list(alpha.start = alpha.start, mu.start = mu.start, tau.start = tau.start, beta.start = beta.start, z.start = z.start)

gibbs = gibbs_MVNW_dp_c(init = init, data = X1, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)


# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 1
# cluster <- data.frame(x1 = X1[1, ],
#                       x2 = X1[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 1", x = "UMAP1", y = "UMAP2")

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta1 by running a Naive DP on group 1 data
beta1.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta1.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta1.start.small = beta1.start.small/sum(beta1.start.small)

# But for the next step again I would require beta1 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta1.start.long = c(beta1.start.small, rep(1e-8, abs(L.max - length(beta1.start.small))))
beta1.start = beta1.start.long/sum(beta1.start.long) # This becomes the starting value of beta1 in the Gibbs iteration
# Draw starting values of latent indicators for group 1
z1.start = rcat(n = n1, prob = beta1.start)

################################################################################
# GROUP 2
################################################################################
gibbs = gibbs_MVNW_dp_c(init = init, data = X2, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 2
# cluster <- data.frame(x1 = X2[1, ],
#                       x2 = X2[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 2", x = "UMAP1", y = "UMAP2")

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 2 data
beta2.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta2.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta2.start.small = beta2.start.small/sum(beta2.start.small)

# But for the next step again I would require beta2 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta2.start.long = c(beta2.start.small, rep(1e-8, abs(L.max - length(beta2.start.small))))
beta2.start = beta2.start.long/sum(beta2.start.long) # This becomes the starting value of beta2 in the Gibbs iteration
# Draw starting values of latent indicators for group 2
n2 = ncol(X2)
z2.start = rcat(n = n2, prob = beta2.start)
################################################################################
# GROUP 3
################################################################################
gibbs = gibbs_MVNW_dp_c(init = init, data = X3, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 3
# cluster <- data.frame(x1 = X3[1, ],
#                       x2 = X3[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 3", x = "UMAP1", y = "UMAP2")

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 3 data
beta3.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta2.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta3.start.small = beta3.start.small/sum(beta3.start.small)

# But for the next step again I would require beta3 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta3.start.long = c(beta3.start.small, rep(1e-8, abs(L.max - length(beta3.start.small))))
beta3.start = beta3.start.long/sum(beta3.start.long) # This becomes the starting value of beta3 in the Gibbs iteration
# Draw starting values of latent indicators for group 3
n3 = ncol(X3)
z3.start = rcat(n = n3, prob = beta3.start)

################################################################################
# GROUP 4
################################################################################
gibbs = gibbs_MVNW_dp_c(init = init, data = X4, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 4
# cluster <- data.frame(x1 = X4[1, ],
#                       x2 = X4[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 4", x = "UMAP1", y = "UMAP2")

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 4 data
beta4.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta2.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta4.start.small = beta4.start.small/sum(beta4.start.small)
# But for the next step again I would require beta4 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta4.start.long = c(beta4.start.small, rep(1e-8, abs(L.max - length(beta4.start.small))))
beta4.start = beta4.start.long/sum(beta4.start.long) # This becomes the starting value of beta4 in the Gibbs iteration
# Draw starting values of latent indicators for group 4
n4 = ncol(X4)
z4.start = rcat(n = n4, prob = beta4.start)
################################################################################
# GROUP 5
################################################################################
gibbs = gibbs_MVNW_dp_c(init = init, data = X5, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 5
# cluster <- data.frame(x1 = X5[1, ],
#                       x2 = X5[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 5", x = "UMAP1", y = "UMAP2")

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 5 data
beta5.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta2.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta5.start.small = beta5.start.small/sum(beta5.start.small)

# But for the next step again I would require beta5 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta5.start.long = c(beta5.start.small, rep(1e-8, abs(L.max - length(beta5.start.small))))
beta5.start = beta5.start.long/sum(beta5.start.long) # This becomes the starting value of beta5 in the Gibbs iteration
# Draw starting values of latent indicators for group 5
n5 = ncol(X5)
z5.start = rcat(n = n5, prob = beta5.start)

################################################################################
# GROUP 6
################################################################################
gibbs = gibbs_MVNW_dp_c(init = init, data = X6, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 6
# cluster <- data.frame(x1 = X6[1, ],
#                       x2 = X6[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 6", x = "UMAP1", y = "UMAP2")

DahlIndex = dp_nw.dahl$DahlIndex

z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 6 data
beta6.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta2.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta6.start.small = beta6.start.small/sum(beta6.start.small)

# But for the next step again I would require beta6 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta6.start.long = c(beta6.start.small, rep(1e-8, abs(L.max - length(beta6.start.small))))
beta6.start = beta6.start.long/sum(beta6.start.long) # This becomes the starting value of beta6 in the Gibbs iteration
# Draw starting values of latent indicators for group 6
n6 = ncol(X6)
z6.start = rcat(n = n6, prob = beta6.start)

################################################################################
# GROUP 7
################################################################################
gibbs = gibbs_MVNW_dp_c(init = init, data = X7, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 7
# cluster <- data.frame(x1 = X7[1, ],
#                       x2 = X7[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 7", x = "UMAP1", y = "UMAP2")


DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 7 data
beta7.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta2.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta7.start.small = beta7.start.small/sum(beta7.start.small)

# But for the next step again I would require beta7 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta7.start.long = c(beta7.start.small, rep(1e-8, abs(L.max - length(beta7.start.small))))
beta7.start = beta7.start.long/sum(beta7.start.long) # This becomes the starting value of beta7 in the Gibbs iteration
# Draw starting values of latent indicators for group 7
n7 = ncol(X7)
z7.start = rcat(n = n7, prob = beta7.start)

################################################################################
# GROUP 8
################################################################################
gibbs = gibbs_MVNW_dp_c(init = init, data = X8, num_iterations = num_iterations, burn = burn, thin = 15, L.max = L.max, prior = prior)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs$z_samples, z.true = NULL)

# # Population 8
# cluster <- data.frame(x1 = X8[1, ],
#                       x2 = X8[2, ],
#                       cluster = factor(gibbs$z_samples[[dp_nw.dahl$DahlIndex]]))
# 
# cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1) + labs(title = "Clustering for group 8", x = "UMAP1", y = "UMAP2")


DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]] # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 8 data
beta8.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
# beta2.start.small = sort(gibbs$beta_samp_post[DahlIndex, ], decreasing = TRUE)
beta8.start.small = beta8.start.small/sum(beta8.start.small)

# But for the next step again I would require beta8 starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta8.start.long = c(beta8.start.small, rep(1e-8, abs(L.max - length(beta8.start.small))))
beta8.start = beta8.start.long/sum(beta8.start.long) # This becomes the starting value of beta8 in the Gibbs iteration
# Draw starting values of latent indicators for group 8
n8 = ncol(X8)
z8.start = rcat(n = n8, prob = beta8.start)


################################################################################
# OTHER INITIALIZATIONS
################################################################################
## Starting values of nu1
nu1.start = beta1.start
nu1.start = nu1.start + 1e-8 # Add 10^-8
nu1.start = nu1.start/sum(nu1.start) # Re-normalize

## Starting values of nu2
nu2.start = beta1.start
nu2.start = nu2.start + 1e-8 # Add 10^-8
nu2.start = nu2.start/sum(nu2.start) # Re-normalize

## Starting values of nu3
nu3.start = beta1.start
nu3.start = nu3.start + 1e-8 # Add 10^-8
nu3.start = nu3.start/sum(nu3.start) # Re-normalize

## Starting values of eta and nu4
eta.start = beta1.start
eta.start = eta.start + 1e-8 # Add 10^-8
eta.start = eta.start/sum(eta.start) # Re-normalize

nu4.start = eta.start
nu4.start = nu4.start + 1e-8 # Add 10^-8
nu4.start = nu4.start/sum(nu4.start) # Re-normalize

#####################################################################################################
# RUNNING GDP
#####################################################################################################
#####################################################################################################
## Initialize storage
#####################################################################################################
# Latent indicators
z1_samples <- list(); z2_samples <- list(); z3_samples <- list(); z4_samples <- list()
z5_samples <- list(); z6_samples <- list(); z7_samples <- list(); z8_samples <- list()
# Weights
beta1_samples <- list(); beta2_samples <- list(); beta3_samples <- list();  beta4_samples <- list()
beta5_samples <- list(); beta6_samples <- list(); beta7_samples <- list();  beta8_samples <- list()
nu1_samples <- list(); nu2_samples <- list(); nu3_samples <- list(); nu4_samples <- list(); eta_samples <- list()
# To save cluster specific sample sizes
n1_samples <- list(); n2_samples <- list(); n3_samples <- list(); n4_samples <- list(); n5_samples <- list() 
n6_samples <- list(); n7_samples <- list(); n8_samples <- list() 
# Storing the means and precisions of the normal distribution
mu_samples <- list() ; tau_samples <- list()
# To store the concentration parameters of the GDP
alpha1_samples <- list(); alpha2_samples <- list(); alpha3_samples <- list(); alpha4_samples <- list()   
alpha5_samples <- list(); alpha6_samples <- list(); alpha7_samples <- list(); alpha8_samples <- list()  
# Initialize all the values at the starting values
beta1_samples[[1]] <- beta1.start; beta2_samples[[1]] <- beta2.start; beta3_samples[[1]] <- beta3.start; beta4_samples[[1]] <- beta4.start; beta5_samples[[1]] <- beta5.start; beta6_samples[[1]] <- beta6.start; beta7_samples[[1]] <- beta7.start; beta8_samples[[1]] <- beta8.start;

nu1_samples[[1]] <- nu1.start; nu2_samples[[1]] <- nu2.start; nu3_samples[[1]] <- nu3.start; nu4_samples[[1]] <- nu4.start
eta_samples[[1]] <- eta.start
mu_samples[[1]] <- mu.start; tau_samples[[1]] <- tau.start 
z1_samples[[1]] <- z1.start; z2_samples[[1]] <- z2.start; z3_samples[[1]] <- z3.start; z4_samples[[1]] <- z4.start 
z5_samples[[1]] <- z5.start; z6_samples[[1]] <- z6.start; z7_samples[[1]] <- z7.start; z8_samples[[1]] <- z8.start 
alpha1_samples[[1]] <- alpha1.start; alpha2_samples[[1]] <- alpha2.start; alpha3_samples[[1]] <- alpha3.start; alpha4_samples[[1]] <- alpha4.start; alpha5_samples[[1]] <- alpha5.start; alpha6_samples[[1]] <- alpha6.start; alpha7_samples[[1]] <- alpha7.start; alpha8_samples[[1]] <- alpha8.start 
# List to store the log likelihood values
log_likelihood1 <- list(); log_likelihood2 <- list(); log_likelihood3 <- list(); log_likelihood4 <- list()
log_likelihood5 <- list(); log_likelihood6 <- list(); log_likelihood7 <- list(); log_likelihood8 <- list()
################################################################################
# Running the Sampler
################################################################################
if(!require("SALTSampler")) install.packages("SALTSampler"); library(SALTSampler)
if(!require("Rcpp")) install.packages("Rcpp"); library(Rcpp)
source("RunMH_my.R")
source("GDP_8_functions_NW_R.R")
sourceCpp("GDP_8_functions_c.cpp")

L.obs = L.max
gibbs_iterations = 50000

nu = nrow(X1)
W = diag(1, nrow = nu) # Rate of Gamma prior

mu0 = as.matrix(rep(0, nu), ncol = 1) # Prior mean of Normal
lambda = 0.01 # Prior precision of Normal

prior_prec = lambda # Prior mean of Normal distribution of G0
prior_mean = mu0 # Prior precision of Normal distribution of G0

p = nrow(X1)

time.start = Sys.time()

for(b in 2:gibbs_iterations){
  # Printing the iterations
  if(b == 2){
    cat(paste0("Iteration: ", (b-1), "\n"))
  }
  if(b %% floor((10/100)*(gibbs_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", b, "\n"))
  }
  
  # Sample the latent indicators for population 1
  sample = sample_z_c(data = t(X1), beta = beta1_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z1_samples[[b]] = sample# Store the samples values of z1
  # Sample the latent indicators for population 2
  sample = sample_z_c(data = t(X2), beta = beta2_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z2_samples[[b]] = sample # Store the samples values of z2
  # Sample the latent indicators for population 2
  sample = sample_z_c(data = t(X3), beta = beta3_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z3_samples[[b]] = sample # Store the samples values of z3
  # Sample the latent indicators for population 4
  sample = sample_z_c(data = t(X4), beta = beta4_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z4_samples[[b]] = sample # Store the samples values of z4
  # Sample the latent indicators for population 5
  sample = sample_z_c(data = t(X5), beta = beta5_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z5_samples[[b]] = sample # Store the samples values of z5
  # Sample the latent indicators for population 6
  sample = sample_z_c(data = t(X6), beta = beta6_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z6_samples[[b]] = sample # Store the samples values of z6
  # Sample the latent indicators for population 7
  sample = sample_z_c(data = t(X7), beta = beta7_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z7_samples[[b]] = sample # Store the samples values of z7
  # Sample the latent indicators for population 8
  sample = sample_z_c(data = t(X8), beta = beta8_samples[[b-1]], mu = mu_samples[[b-1]], tau = tau_samples[[b-1]], L = L.obs)
  
  z8_samples[[b]] = sample # Store the samples values of z8
  
  
  mu_tau_sample = sample_mu_tau_gdp_c(t(X1), t(X2), t(X3), t(X4), t(X5), t(X6), t(X7), t(X8),
                                      z1 = z1_samples[[b]], z2 = z2_samples[[b]], z3 = z3_samples[[b]], 
                                      z4 = z4_samples[[b]], z5 = z5_samples[[b]], z6 = z6_samples[[b]], 
                                      z7 = z7_samples[[b]], z8 = z8_samples[[b]], L = L.obs,
                                      nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
  # Store the samples of mu, tau, n1, n2, n3, n4, n5, n6, n7 and n8
  mu_samples[[b]] =  t(mu_tau_sample$mu.draw);
  tau_samples[[b]] =  mu_tau_sample$tau.draw
  
  n1_samples[[b]] = mu_tau_sample$n1_h; n2_samples[[b]] = mu_tau_sample$n2_h; n3_samples[[b]] = mu_tau_sample$n3_h; n4_samples[[b]] = mu_tau_sample$n4_h
  n5_samples[[b]] = mu_tau_sample$n5_h; n6_samples[[b]] = mu_tau_sample$n6_h; n7_samples[[b]] = mu_tau_sample$n7_h; n8_samples[[b]] = mu_tau_sample$n8_h
  
  init = beta1.start
  
  # pars should include beta2, beta3, beta4, nu1, nu2, nu3, eta, alpha2, alpha3, alpha4
  # Run SALTSampler for Beta1
  beta1 <- RunMH_my(Target = Target_beta1, init = init, B = 1,
                    h = rep(2, L.obs), 
                    pars = list(alpha2 = alpha2_samples[[b-1]],
                                alpha3 = alpha3_samples[[b-1]],
                                alpha4 = alpha4_samples[[b-1]],
                                beta2 = beta2_samples[[b-1]],
                                beta3 = beta3_samples[[b-1]],
                                beta4 = beta4_samples[[b-1]],
                                nu1 = nu1_samples[[b-1]],
                                nu2 = nu2_samples[[b-1]],
                                nu3 = nu3_samples[[b-1]],
                                eta = eta_samples[[b-1]]),
                    dat = n1_samples[[b]], 
                    concentration = rep(alpha1_samples[[b-1]], length(init)))
  
  beta1_samples[[b]] = as.numeric(beta1$S)
  # Draw samples from the beta2
  new_concentration = alpha2_samples[[b-1]] * beta1_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta2_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n2_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta2_samples[[b]] = beta2_samples[[b]] + 1e-8 
  beta2_samples[[b]] = beta2_samples[[b]]/sum(beta2_samples[[b]]) # Re-normalize
  
  # Draw samples from the beta3
  new_concentration = alpha3_samples[[b-1]] * beta1_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta3_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n3_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta3_samples[[b]] = beta3_samples[[b]] + 1e-8 
  beta3_samples[[b]] = beta3_samples[[b]]/sum(beta3_samples[[b]]) # Re-normalize
  
  # Draw samples from the beta4
  new_concentration = alpha4_samples[[b-1]] * beta1_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta4_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n4_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta4_samples[[b]] = beta4_samples[[b]] + 1e-8 
  beta4_samples[[b]] = beta4_samples[[b]]/sum(beta4_samples[[b]]) # Re-normalize
  
  
  # Next sample nu1
  # pars should be a named list with names beta1, beta5,  alpha2, alpha4 and alpha5
  # a corresponds to alpha2 + alpha4
  # init = nu1_samples[[b-1]]
  init = nu1.start
  nu1 = RunMH_my(Target = Target_nu1, 
                 init = init, B = 1,
                 h = rep(2, L.obs), 
                 pars = list(alpha2 = alpha2_samples[[b-1]],
                             alpha4 = alpha4_samples[[b-1]],
                             alpha5 = alpha5_samples[[b-1]],
                             beta1 = beta1_samples[[b]],
                             beta5 = beta5_samples[[b-1]]), 
                 dat = NULL, 
                 concentration = rep((alpha2_samples[[b-1]] + alpha4_samples[[b-1]]), L.obs))
  
  nu1_samples[[b]] = as.numeric(nu1$S)
  
  # Draw samples from the beta5
  new_concentration = alpha5_samples[[b-1]] * nu1_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta5_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n5_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta5_samples[[b]] = beta5_samples[[b]] + 1e-8 
  beta5_samples[[b]] = beta5_samples[[b]]/sum(beta5_samples[[b]]) # Re-normalize
  
  
  # Next sample nu2
  # pars should be a named list with names beta1, beta6,  alpha2, alpha3 and alpha6.
  # a corresponds to alpha2 + alpha3
  # init = nu2_samples[[b-1]]
  init = nu2.start
  nu2 = RunMH_my(Target = Target_nu2, 
                 init = init, B = 1,
                 h = rep(2, L.obs), 
                 pars = list(alpha2 = alpha2_samples[[b-1]],
                             alpha3 = alpha3_samples[[b-1]],
                             alpha6 = alpha6_samples[[b-1]],
                             beta1 = beta1_samples[[b]],
                             beta6 = beta6_samples[[b-1]]), 
                 dat = NULL, 
                 concentration = rep((alpha2_samples[[b-1]] + alpha3_samples[[b-1]]), L.obs))
  
  nu2_samples[[b]] = as.numeric(nu2$S)
  
  # Draw samples from the beta6
  new_concentration = alpha6_samples[[b-1]] * nu2_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta6_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n6_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta6_samples[[b]] = beta6_samples[[b]] + 1e-8 
  beta6_samples[[b]] = beta6_samples[[b]]/sum(beta6_samples[[b]]) # Re-normalize
  
  # Next sample nu3
  # pars should be a named list with names beta1, beta7,  alpha3, alpha4 and alpha7
  # a corresponds to alpha3 + alpha4
  # init = nu3_samples[[b-1]]
  init = nu3.start
  nu3 = RunMH_my(Target = Target_nu3, 
                 init = nu3.start, B = 1,
                 h = rep(2, L.obs), 
                 pars = list(alpha3 = alpha3_samples[[b-1]],
                             alpha4 = alpha4_samples[[b-1]],
                             alpha7 = alpha7_samples[[b-1]],
                             beta1 = beta1_samples[[b]],
                             beta7 = beta7_samples[[b-1]]), 
                 dat = NULL, 
                 concentration = rep((alpha3_samples[[b-1]] + alpha4_samples[[b-1]]), L.obs))
  
  nu3_samples[[b]] = as.numeric(nu3$S)
  
  # Draw samples from the beta7
  new_concentration = alpha7_samples[[b-1]] * nu3_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta7_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n7_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta7_samples[[b]] = beta7_samples[[b]] + 1e-8 
  beta7_samples[[b]] = beta7_samples[[b]]/sum(beta7_samples[[b]]) # Re-normalize
  
  # Next sample eta
  # pars should be a named list with names beta1, nu4,  alpha2, alpha3, alpha4, alpha5, alpha6, alpha7
  # a corresponds to (alpha2 + alpha3 + alpha4)
  # init = eta_samples[[b-1]]
  init = eta.start
  eta = RunMH_my(Target = Target_eta, 
                 init = init, B = 1,
                 h = rep(2, L.obs), 
                 pars = list(alpha2 = alpha2_samples[[b-1]],
                             alpha3 = alpha3_samples[[b-1]],
                             alpha4 = alpha4_samples[[b-1]],
                             alpha5 = alpha5_samples[[b-1]],
                             alpha6 = alpha6_samples[[b-1]],
                             alpha7 = alpha7_samples[[b-1]],
                             beta1 = beta1_samples[[b]],
                             nu4 = nu4_samples[[b-1]]), 
                 dat = NULL, 
                 concentration = rep((alpha2_samples[[b-1]] + alpha3_samples[[b-1]] + alpha4_samples[[b-1]]), L.obs))
  
  eta_samples[[b]] = as.numeric(eta$S)
  
  # Next sample nu4
  # pars should be a named list with names beta8, eta,  alpha5, alpha6, alpha7 and alpha8.
  # a corresponds to alpha5 + alpha6 + alpha7
  # init = nu4_samples[[b-1]]
  init = nu4.start
  nu4 = RunMH_my(Target = Target_nu4, 
                 init = init, B = 1,
                 h = rep(2, L.obs), 
                 pars = list(alpha5 = alpha5_samples[[b-1]],
                             alpha6 = alpha6_samples[[b-1]],
                             alpha7 = alpha7_samples[[b-1]],
                             alpha8 = alpha8_samples[[b-1]],
                             beta8 = beta8_samples[[b-1]],
                             eta = eta_samples[[b]]), 
                 dat = NULL, 
                 concentration = rep((alpha5_samples[[b-1]] + alpha6_samples[[b-1]] + alpha7_samples[[b-1]]), L.obs))
  
  nu4_samples[[b]] = as.numeric(nu4$S)
  
  # Draw samples from the beta8
  new_concentration = alpha8_samples[[b-1]] * nu4_samples[[b]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta8_samples[[b]] = as.numeric(rdirichlet(n = 1, alpha = n8_samples[[b]] + new_concentration))
  # To avoid numerical issues I add a small quantity to each coordinate and re-normalize
  beta8_samples[[b]] = beta8_samples[[b]] + 1e-8 
  beta8_samples[[b]] = beta8_samples[[b]]/sum(beta8_samples[[b]]) # Re-normalize
  
  # Sample the alphas
  # pars argument should be a named list with alpha0, alpha2, alpha3, alpha4 and beta1
  alpha1_samples[[b]] = sample_alpha1(ycurrent = alpha1_samples[[b-1]], 
                                      pars = list(beta1 = beta1_samples[[b]], 
                                                  alpha0 = alpha0,
                                                  alpha2 = alpha2_samples[[b-1]],
                                                  alpha3 = alpha3_samples[[b-1]],
                                                  alpha4 = alpha4_samples[[b-1]]))$out
  # @pars argument should be a named list with alpha1, alpha3, alpha4, alpha5, alpha6,
  # beta1, beta2, nu1, nu2, eta
  alpha2_samples[[b]] = sample_alpha2(ycurrent = alpha2_samples[[b-1]], 
                                      pars = list(alpha1 = alpha1_samples[[b]], 
                                                  alpha3 = alpha3_samples[[b-1]], 
                                                  alpha4 = alpha4_samples[[b-1]],
                                                  alpha5 = alpha5_samples[[b-1]],
                                                  alpha6 = alpha6_samples[[b-1]],
                                                  beta1 = beta1_samples[[b]], 
                                                  beta2 = beta2_samples[[b]], 
                                                  nu1 = nu1_samples[[b]],
                                                  nu2 = nu2_samples[[b]],
                                                  eta = eta_samples[[b]]))$out
  
  # @pars argument should be a named list with alpha1, alpha2, alpha4, alpha6, alpha7,
  # beta1, beta3, nu2, nu3, eta
  alpha3_samples[[b]] = sample_alpha3(ycurrent = alpha3_samples[[b-1]], 
                                      pars = list(alpha1 = alpha1_samples[[b]], 
                                                  alpha2 = alpha2_samples[[b]], 
                                                  alpha4 = alpha4_samples[[b-1]], 
                                                  alpha6 = alpha6_samples[[b-1]],
                                                  alpha7 = alpha7_samples[[b-1]],
                                                  beta1 = beta1_samples[[b]], 
                                                  beta3 = beta3_samples[[b]], 
                                                  nu2 = nu2_samples[[b]],
                                                  nu3 = nu3_samples[[b]],
                                                  eta = eta_samples[[b]]))$out
  
  # @pars argument should be a named list with alpha1, alpha2, alpha3, alpha5, alpha7,
  # beta1, beta4, nu1, nu3, eta
  alpha4_samples[[b]] = sample_alpha4(ycurrent = alpha4_samples[[b-1]], 
                                      pars = list(alpha1 = alpha1_samples[[b]], 
                                                  alpha2 = alpha2_samples[[b]], 
                                                  alpha3 = alpha3_samples[[b]], 
                                                  alpha5 = alpha5_samples[[b-1]],
                                                  alpha7 = alpha7_samples[[b-1]],
                                                  beta1 = beta1_samples[[b]], 
                                                  beta4 = beta4_samples[[b]], 
                                                  nu1 = nu1_samples[[b]],
                                                  nu3 = nu3_samples[[b]],
                                                  eta = eta_samples[[b]]))$out
  
  # @pars argument should be a named list with alpha2, alpha4, alpha6, alpha7, alpha8
  # beta5, nu1, nu4, eta
  alpha5_samples[[b]] = sample_alpha5(ycurrent = alpha5_samples[[b-1]], 
                                      pars = list(alpha2 = alpha2_samples[[b]], 
                                                  alpha4 = alpha4_samples[[b]], 
                                                  alpha6 = alpha6_samples[[b-1]], 
                                                  alpha7 = alpha7_samples[[b-1]],
                                                  alpha8 = alpha8_samples[[b-1]],
                                                  beta5 = beta5_samples[[b]], 
                                                  nu1 = nu1_samples[[b]],
                                                  nu4 = nu4_samples[[b]],
                                                  eta = eta_samples[[b]]))$out
  
  # @pars argument should be a named list with alpha2, alpha3, alpha5, alpha7, alpha8
  # beta6, nu2, nu4, eta
  alpha6_samples[[b]] = sample_alpha6(ycurrent = alpha6_samples[[b-1]], 
                                      pars = list(alpha2 = alpha2_samples[[b]], 
                                                  alpha3 = alpha3_samples[[b]], 
                                                  alpha5 = alpha5_samples[[b]], 
                                                  alpha7 = alpha7_samples[[b-1]],
                                                  alpha8 = alpha8_samples[[b-1]],
                                                  beta6 = beta6_samples[[b]], 
                                                  nu2 = nu2_samples[[b]],
                                                  nu4 = nu4_samples[[b]],
                                                  eta = eta_samples[[b]]))$out
  
  # @pars argument should be a named list with alpha3, alpha4, alpha5, alpha6, alpha8
  # beta7, nu3, nu4, eta
  alpha7_samples[[b]] = sample_alpha7(ycurrent = alpha7_samples[[b-1]], 
                                      pars = list(alpha3 = alpha3_samples[[b]], 
                                                  alpha4 = alpha4_samples[[b]], 
                                                  alpha5 = alpha5_samples[[b]], 
                                                  alpha6 = alpha6_samples[[b]],
                                                  alpha8 = alpha8_samples[[b-1]],
                                                  beta7 = beta7_samples[[b]], 
                                                  nu3 = nu2_samples[[b]],
                                                  nu4 = nu4_samples[[b]],
                                                  eta = eta_samples[[b]]))$out
  # @pars argument should be a named list with alpha5, alpha6, alpha7
  # beta8, nu4
  alpha8_samples[[b]] = sample_alpha8(ycurrent = alpha8_samples[[b-1]], 
                                      pars = list(alpha5 = alpha5_samples[[b]], 
                                                  alpha6 = alpha6_samples[[b]],
                                                  alpha7 = alpha7_samples[[b]],
                                                  beta8 = beta8_samples[[b]], 
                                                  nu4 = nu4_samples[[b]]))$out
  
  ## Log likelihoods to monitor convergence of the chain
  # Calculate the log likelihood for group 1 for every iteration
  log_likelihood1[[b]] = log_l(data = t(X1), beta = beta1_samples[[b]], z = z1_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  # Calculate the log likelihood for group 2 for every iteration
  log_likelihood2[[b]] = log_l(data = t(X2), beta = beta2_samples[[b]], z = z2_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  # Calculate the log likelihood for group 3 for every iteration
  log_likelihood3[[b]] = log_l(data = t(X3), beta = beta3_samples[[b]], z = z3_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  # Calculate the log likelihood for group 4 for every iteration
  log_likelihood4[[b]] = log_l(data = t(X4), beta = beta4_samples[[b]], z = z4_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  # Calculate the log likelihood for group 5 for every iteration
  log_likelihood5[[b]] = log_l(data = t(X5), beta = beta5_samples[[b]], z = z5_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  # Calculate the log likelihood for group 6 for every iteration
  log_likelihood6[[b]] = log_l(data = t(X6), beta = beta6_samples[[b]], z = z6_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  # Calculate the log likelihood for group 7 for every iteration
  log_likelihood7[[b]] = log_l(data = t(X7), beta = beta7_samples[[b]], z = z7_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  # Calculate the log likelihood for group 8 for every iteration
  log_likelihood8[[b]] = log_l(data = t(X8), beta = beta8_samples[[b]], z = z8_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
  
}

time.end = Sys.time()
time.taken = time.end - time.start

gibbs_burn <- 35000 # Burn-in
thin = 15 # Thinning factor
thin_samples <- seq(from = (gibbs_burn + 1), to = gibbs_iterations, by = thin) # Thinned samples


#################################################################################################################
# PLOT OF TRACEPLOT OF LOG-LIKELIHOOD
#################################################################################################################
log_like1 = unlist(log_likelihood1[thin_samples])
log_like2 = unlist(log_likelihood2[thin_samples])
log_like3 = unlist(log_likelihood3[thin_samples])
log_like4 = unlist(log_likelihood4[thin_samples])
log_like5 = unlist(log_likelihood5[thin_samples])
log_like6 = unlist(log_likelihood6[thin_samples])
log_like7 = unlist(log_likelihood7[thin_samples])
log_like8 = unlist(log_likelihood8[thin_samples])

l1 = data.frame(x = 1:length(log_like1), y = log_like1) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 1", x = "Iteration", y = "LL") +
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l2 = data.frame(x = 1:length(log_like2), y = log_like2) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 2", x = "Iteration", y = "LL") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l3 = data.frame(x = 1:length(log_like3), y = log_like3) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 3", x = "Iteration", y = "LL") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l4 = data.frame(x = 1:length(log_like4), y = log_like4) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 4", x = "Iteration", y = "LL") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l5 = data.frame(x = 1:length(log_like5), y = log_like5) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 5", x = "Iteration", y = "LL") +
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l6 = data.frame(x = 1:length(log_like6), y = log_like6) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 6", x = "Iteration", y = "LL") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l7 = data.frame(x = 1:length(log_like7), y = log_like7) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 7", x = "Iteration", y = "LL") + theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l8 = data.frame(x = 1:length(log_like8), y = log_like8) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 8", x = "Iteration", y = "LL") + theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
gridExtra::grid.arrange(l1, l2, l3, l4, l5, l6, l7, l8, top=textGrob("Traceplot plot of log-likelihood for different groups", gp = gpar(fontface = "bold", fontsize = 16)))

#################################################################################################################
# EFFECTIVE SAMPLE SIZES
#################################################################################################################
if(!require("coda")) install.packages("coda"); library(coda)
effectiveSize(log_like1)
effectiveSize(log_like2)
effectiveSize(log_like3)
effectiveSize(log_like4)
effectiveSize(log_like5)
effectiveSize(log_like6)
effectiveSize(log_like7)
effectiveSize(log_like8)


#################################################################################################################
# IF THE MODEL SHOWS CONVERGENCE THEN SAVE THE INFO OF THE MCMC 
#################################################################################################################
# Extract the cluster indicators post burn-in for the thinned samples
z1_samples.post.burn <- z1_samples[thin_samples]
z2_samples.post.burn <- z2_samples[thin_samples]
z3_samples.post.burn <- z3_samples[thin_samples]
z4_samples.post.burn <- z4_samples[thin_samples]
z5_samples.post.burn <- z5_samples[thin_samples]
z6_samples.post.burn <- z6_samples[thin_samples]
z7_samples.post.burn <- z7_samples[thin_samples]
z8_samples.post.burn <- z8_samples[thin_samples]

z.chain <- list(z1_samples.post.burn,
                z2_samples.post.burn,
                z3_samples.post.burn,
                z4_samples.post.burn,
                z5_samples.post.burn, 
                z6_samples.post.burn,
                z7_samples.post.burn,
                z8_samples.post.burn)
#####################################################################################
# UNCOMMENT AND SAVE
#####################################################################################
# saveRDS(z.chain, "/Users/arhitchakrabarti/Desktop/New GDP UMAP/GDP_z_chain_l3_alpha_3_L_30.rds")

LL <- list(log_like1,
           log_like2,
           log_like3,
           log_like4,
           log_like5, 
           log_like6,
           log_like7,
           log_like8)

#####################################################################################
# UNCOMMENT AND SAVE
#####################################################################################
# saveRDS(LL, "/Users/arhitchakrabarti/Desktop/New GDP UMAP/GDP_LL_chain_l3_alpha_3_L_30.rds")
