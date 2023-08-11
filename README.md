Graphical Dirichlet process
================

## Installation

You can install GDPSamp R package from GitHub with:

``` r
# install.packages('devtools') install the package without building
# vignettes
devtools::install_github("Arhit-Chakrabarti/GDPSamp")
```

**Note**: The library *devtools* needs to be installed before installing
the R package from GitHub. Here is a detailed example of a simulation
study. The example may be used to understand the various functionalities
of the R package.

## Simulation Example

We first generate the data from the eight dependent groups. The
dependency between the groups can be described by the underlying DAG as
mentioned in the paper.

``` r
# Source all the required libraries
library(extraDistr)
library(MASS)
library(mvtnorm)
library(fossil)
library(tidyverse)
library(SALTSampler)
library(GDPSamp)

# Generate the data

L.true = 10  # true number of groups in population 1
alpha0 = 15

alpha1.true = 10  # True alpha1 for generating data
alpha2.true = 12  # True alpha2 for generating data
alpha3.true = 15  # True alpha3 for generating data
alpha4.true = 13  # True alpha4 for generating data
alpha5.true = 20  # True alpha5 for generating data
alpha6.true = 24  # True alpha6 for generating data
alpha7.true = 22  # True alpha7 for generating data
alpha8.true = 30  # True alpha8 for generating data


# True weights
beta1.true = as.numeric(rep(alpha1.true/L.true, L.true))  # True beta1
beta1.true[is.na(beta1.true)] <- 0
beta1.true = beta1.true + 1e-06
beta1.true = beta1.true/sum(beta1.true)
beta2.true = c(rep(0.15, 5), rep(0.05, 3), rep(0, 2))  # True beta2
beta2.true[is.na(beta2.true)] <- 0
beta2.true = beta2.true + 1e-06
beta2.true = beta2.true/sum(beta2.true)
beta3.true = c(rep(0.08, 3), rep(0, 2), rep(0.12, 5))  # True beta3
beta3.true[is.na(beta3.true)] <- 0
beta3.true = beta3.true + 1e-06
beta3.true = beta3.true/sum(beta3.true)
beta4.true = c(rep(0.03, 3), rep(0.18, 5), rep(0, 2))  # True beta4
beta4.true[is.na(beta4.true)] <- 0
beta4.true = beta4.true + 1e-06
beta4.true = beta4.true/sum(beta4.true)

nu1.true = as.numeric((alpha2.true + alpha4.true) * beta1.true)/L.true  #True nu1
nu1.true[is.na(nu1.true)] <- 0
nu1.true = nu1.true + 1e-06
nu1.true = nu1.true/sum(nu1.true)
beta5.true = beta2.true + beta4.true  # True beta5
beta5.true[is.na(beta5.true)] <- 0
beta5.true = beta5.true + 1e-06
beta5.true = beta5.true/sum(beta5.true)

nu2.true = as.numeric((alpha2.true + alpha3.true) * beta1.true)/L.true  #True nu2
nu2.true[is.na(nu2.true)] <- 0
nu2.true = nu2.true + 1e-06
nu2.true = nu2.true/sum(nu2.true)
beta6.true = beta2.true + beta3.true  # True beta6
beta6.true[is.na(beta6.true)] <- 0
beta6.true = beta6.true + 1e-06
beta6.true = beta6.true/sum(beta6.true)

nu3.true = as.numeric((alpha3.true + alpha4.true) * beta1.true)/L.true  #True nu3
nu3.true[is.na(nu3.true)] <- 0
nu3.true = nu3.true + 1e-06
nu3.true = nu3.true/sum(nu3.true)
beta7.true = beta3.true + beta4.true  # True beta7
beta7.true[is.na(beta7.true)] <- 0
beta7.true = beta7.true + 1e-06
beta7.true = beta7.true/sum(beta7.true)

eta.true = as.numeric((2 * (alpha2.true + alpha3.true + alpha4.true) *
    beta1.true)/L.true)  # True eta
eta.true[is.na(eta.true)] <- 0
eta.true = eta.true + 1e-06
eta.true = eta.true/sum(eta.true)

nu4.true = as.numeric((alpha5.true + alpha6.true + alpha7.true) * eta.true)/L.true  # True nu4
nu4.true[is.na(nu4.true)] <- 0
nu4.true = nu4.true + 1e-06
nu4.true = nu4.true/sum(nu4.true)
beta8.true = beta5.true + beta6.true + beta7.true  # True beta8
beta8.true[is.na(beta8.true)] <- 0
beta8.true = beta8.true + 1e-06
beta8.true = beta8.true/sum(beta8.true)

# Sample sizes
n1 = 80  # First population
n2 = 70  # Second population
n3 = 70  # Third population
n4 = 75  # Forth population
n5 = 83  # Fifth population
n6 = 88  # Sixth population
n7 = 92  # Seventh population
n8 = 88  # Eighth population

z1.true = sample(1:L.true, size = n1, prob = beta1.true, replace = TRUE)
z2.true = sample(1:L.true, size = n2, prob = beta2.true, replace = TRUE)
z3.true = sample(1:L.true, size = n3, prob = beta3.true, replace = TRUE)
z4.true = sample(1:L.true, size = n4, prob = beta4.true, replace = TRUE)
z5.true = sample(1:L.true, size = n5, prob = beta5.true, replace = TRUE)
z6.true = sample(1:L.true, size = n6, prob = beta6.true, replace = TRUE)
z7.true = sample(1:L.true, size = n7, prob = beta7.true, replace = TRUE)
z8.true = sample(1:L.true, size = n8, prob = beta8.true, replace = TRUE)
# Draw data from Normal populations True cluster specific means
mean.true = cbind(c(-2.5, 0), c(0, 0), c(2.5, 0), c(2.5, -2.5), c(-3, -3),
    c(2, 2), c(-2, 5), c(5, 8), c(-5, -8), c(8, -8))
# True covariance matrix for for different poulation
Sigma1.true = matrix(c(0.8, 0.3, 0.3, 0.8), nrow = 2)
Sigma2.true = matrix(c(0.85, 0.25, 0.25, 0.85), nrow = 2)
Sigma3.true = matrix(c(1, 0.1, 0.1, 1), nrow = 2)
Sigma4.true = matrix(c(0.8, -0.1, -0.1, 0.8), nrow = 2)
Sigma5.true = matrix(c(0.9, -0.2, -0.2, 0.9), nrow = 2)
Sigma6.true = matrix(c(0.8, 0, 0, 0.8), nrow = 2)
Sigma7.true = matrix(c(0.75, 0.2, 0.2, 0.75), nrow = 2)
Sigma8.true = matrix(c(1.1, 0.1, 0.1, 1.1), nrow = 2)

X1 = sapply(1:n1, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z1.true][, j], Sigma = Sigma1.true)
})  # Data from population 1
X2 = sapply(1:n2, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z2.true][, j], Sigma = Sigma2.true)
})  # Data from population 2 
X3 = sapply(1:n3, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z3.true][, j], Sigma = Sigma3.true)
})  # Data from population 3 
X4 = sapply(1:n4, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z4.true][, j], Sigma = Sigma4.true)
})  # Data from population 4 
X5 = sapply(1:n5, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z5.true][, j], Sigma = Sigma5.true)
})  # Data from population 5 
X6 = sapply(1:n6, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z6.true][, j], Sigma = Sigma6.true)
})  # Data from population 6 
X7 = sapply(1:n7, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z7.true][, j], Sigma = Sigma7.true)
})  # Data from population 7 
X8 = sapply(1:n8, function(j) {
    MASS::mvrnorm(n = 1, mu = mean.true[, z8.true][, j], Sigma = Sigma8.true)
})  # Data from population 8 
```

Once we have our data, we first run a vanilla DP for each group. We
estimate of the stick breaking weights for each group from this vanilla
DP. We will use this estimated stick-breaking weights as initial values
in our proposed GDP sampler. This ensures fast mixing of the proposed
MCMC sampler.

``` r
## Prior specifications

alpha0 = 2
L.max = 10  # Truncation level of DP
# Prior for Dirichlet distribution for group 1
alpha1.start = rgamma(n = 1, shape = alpha0)
# Prior for Dirichlet distribution for group 2
alpha2.start = rgamma(n = 1, shape = alpha1.start)
# Prior for Dirichlet distribution for group 3
alpha3.start = rgamma(n = 1, shape = alpha1.start)
# Prior for Dirichlet distribution for group 4
alpha4.start = rgamma(n = 1, shape = alpha1.start)
# Prior for Dirichlet distribution for group 5
alpha5.start = rgamma(n = 1, shape = (alpha2.start + alpha4.start))
# Prior for Dirichlet distribution for group 6
alpha6.start = rgamma(n = 1, shape = (alpha2.start + alpha3.start))
# Prior for Dirichlet distribution for group 7
alpha7.start = rgamma(n = 1, shape = (alpha3.start + alpha4.start))
# Prior for Dirichlet distribution for group 8
alpha8.start = rgamma(n = 1, shape = (alpha5.start + alpha5.start + alpha7.start))

## SPECIFICATION

nu <- base::nrow(X1)
W = diag(1, nrow = nu)  # Rate of Gamma prior
mu0 = as.matrix(rep(0, nu), ncol = 1)  # Prior mean of Normal
lambda = 0.01  # Prior precision of Normal
mu.start.list = list()
tau.start = rWishart(n = L.max, df = nu, Sigma = W)
for (i in 1:L.max) {
    mu.start.list[[i]] = matrix(mvtnorm::rmvnorm(n = 1, mean = as.numeric(mu0),
        sigma = solve(lambda * tau.start[, , i])), ncol = 1)
}


beta.start = rep(1/L.max, L.max)
z.start = extraDistr::rcat(n = n1, prob = beta.start)
alpha.start = rgamma(n = 1, shape = alpha0)
mu.start <- t(sapply(1:L.max, function(j) {
    mu.start.list[[j]]
}))

num_iterations <- 15000
burn = 5000

# RUN DP ON GROUP 1

alpha0 = 2
nu = nrow(X1)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
# Initial values of different parameters
init = list(alpha.start = alpha.start, mu.start = mu.start, tau.start = tau.start,
    beta.start = beta.start, z.start = z.start)

gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X1, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z1.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta1 by running a Naive DP on group 1 data
beta1.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta1.start.small = beta1.start.small/sum(beta1.start.small)

# But for the next step again I would require beta1 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta1.start.long = c(beta1.start.small, rep(1e-06, abs(L.max - length(beta1.start.small))))
# This becomes the starting value of beta1 in the Gibbs iteration
beta1.start = beta1.start.long/sum(beta1.start.long)
beta1.start[is.na(beta1.start)] <- 0
beta1.start = beta1.start/sum(beta1.start)
# Draw starting values of latent indicators for group 1
z1.start = extraDistr::rcat(n = n1, prob = beta1.start)


# RUN DP ON GROUP 2

alpha0 = 2
nu = nrow(X2)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X2, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z2.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 2 data
beta2.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta2.start.small = beta2.start.small/sum(beta2.start.small)

# But for the next step again I would require beta2 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta2.start.long = c(beta2.start.small, rep(1e-06, abs(L.max - length(beta2.start.small))))
# This becomes the starting value of beta2 in the Gibbs iteration
beta2.start = beta2.start.long/sum(beta2.start.long)
# Draw starting values of latent indicators for group 2
beta2.start[is.na(beta2.start)] <- 0
beta2.start = beta2.start/sum(beta2.start)
z2.start = extraDistr::rcat(n = n2, prob = beta2.start)

# RUN DP ON GROUP 3

alpha0 = 2
nu = nrow(X3)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X3, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z3.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 3 data
beta3.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta3.start.small = beta3.start.small/sum(beta3.start.small)

# But for the next step again I would require beta3 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta3.start.long = c(beta3.start.small, rep(1e-06, abs(L.max - length(beta3.start.small))))
# This becomes the starting value of beta3 in the Gibbs iteration
beta3.start = beta3.start.long/sum(beta3.start.long)
beta3.start[is.na(beta3.start)] <- 0
beta3.start = beta3.start/sum(beta3.start)
# Draw starting values of latent indicators for group 3
z3.start = extraDistr::rcat(n = n3, prob = beta3.start)

# RUN DP ON GROUP 4

alpha0 = 2
nu = nrow(X4)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X4, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z4.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 4 data
beta4.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta4.start.small = beta4.start.small/sum(beta4.start.small)
# But for the next step again I would require beta4 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta4.start.long = c(beta4.start.small, rep(1e-06, abs(L.max - length(beta4.start.small))))
# This becomes the starting value of beta4 in the Gibbs iteration
beta4.start = beta4.start.long/sum(beta4.start.long)
beta4.start[is.na(beta4.start)] <- 0
beta4.start = beta4.start/sum(beta4.start)
# Draw starting values of latent indicators for group 4
z4.start = extraDistr::rcat(n = n4, prob = beta4.start)

# RUN DP ON GROUP 5

alpha0 = 2
nu = nrow(X5)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X5, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z5.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 5 data
beta5.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta5.start.small = beta5.start.small/sum(beta5.start.small)

# But for the next step again I would require beta5 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta5.start.long = c(beta5.start.small, rep(1e-06, abs(L.max - length(beta5.start.small))))
# This becomes the starting value of beta5 in the Gibbs iteration
beta5.start = beta5.start.long/sum(beta5.start.long)
beta5.start[is.na(beta5.start)] <- 0
beta5.start = beta5.start/sum(beta5.start)
# Draw starting values of latent indicators for group 5
z5.start = extraDistr::rcat(n = n5, prob = beta5.start)

# RUN DP ON GROUP 6

alpha0 = 2
nu = nrow(X6)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X6, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z6.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 6 data
beta6.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta6.start.small = beta6.start.small/sum(beta6.start.small)

# But for the next step again I would require beta6 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta6.start.long = c(beta6.start.small, rep(1e-06, abs(L.max - length(beta6.start.small))))
# This becomes the starting value of beta6 in the Gibbs iteration
beta6.start = beta6.start.long/sum(beta6.start.long)
beta6.start[is.na(beta6.start)] <- 0
beta6.start = beta6.start/sum(beta6.start)
# Draw starting values of latent indicators for group 6
z6.start = extraDistr::rcat(n = n6, prob = beta6.start)

# RUN DP ON GROUP 7

alpha0 = 2
nu = nrow(X7)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X7, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z7.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 7 data
beta7.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta7.start.small = beta7.start.small/sum(beta7.start.small)

# But for the next step again I would require beta7 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta7.start.long = c(beta7.start.small, rep(1e-06, abs(L.max - length(beta7.start.small))))
# This becomes the starting value of beta7 in the Gibbs iteration
beta7.start = beta7.start.long/sum(beta7.start.long)
# Draw starting values of latent indicators for group 7
beta7.start[is.na(beta7.start)] <- 0
beta7.start = beta7.start/sum(beta7.start)
z7.start = extraDistr::rcat(n = n7, prob = beta7.start)

# RUN DP ON GROUP 8

alpha0 = 2
nu = nrow(X8)
W = diag(1, nrow = nu)
mu0 = as.matrix(rep(0, nu), ncol = 1)
lambda = 0.01
# Prior hyper-parameters
prior = list(alpha0 = alpha0, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)
gibbs <- GDPSamp::gibbs_MVNW_dp_c(init = init, data = X8, num_iterations = num_iterations,
    burn = burn, thin = 1, L.max = L.max, prior = prior)

# DAHL's METHOD
dp_nw.dahl = GDPSamp::getDahl(gibbs$z_samples, z8.true)

DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs$z_samples[[DahlIndex]]  # Indices with minimum Dahl Index

# Get the estimate of beta2 by running a Naive DP on group 8 data
beta8.start.small = gibbs$beta_samp_post[[DahlIndex]][unique(z.dahl)]
beta8.start.small = beta8.start.small/sum(beta8.start.small)

# But for the next step again I would require beta8 starting value to
# be of length L.max So I add a very small quantity to the remaining
# co-ordinates and re-normlaize
beta8.start.long = c(beta8.start.small, rep(1e-06, abs(L.max - length(beta8.start.small))))
# This becomes the starting value of beta8 in the Gibbs iteration
beta8.start = beta8.start.long/sum(beta8.start.long)
beta8.start[is.na(beta8.start)] <- 0
beta8.start = beta8.start/sum(beta8.start)
# Draw starting values of latent indicators for group 8
z8.start = extraDistr::rcat(n = n8, prob = beta8.start)
```

Once we have all the initilializations, we are ready to run our proposed
GDP Sampler. We consider our sampler here, wherein we consider 25,000
iterations of our sampler and 10,000 burn-in and thinning by a factor of
10.

``` r
# OTHER INITIALIZATIONS FOR RUNNING GDP SAMPLER

## Starting values of nu1
nu1.start = beta1.start
nu1.start = nu1.start + 1e-08  # Add 10^-8
nu1.start = nu1.start/sum(nu1.start)  # Re-normalize

## Starting values of nu2
nu2.start = beta1.start
nu2.start = nu2.start + 1e-08  # Add 10^-8
nu2.start = nu2.start/sum(nu2.start)  # Re-normalize

## Starting values of nu3
nu3.start = beta1.start
nu3.start = nu3.start + 1e-08  # Add 10^-8
nu3.start = nu3.start/sum(nu3.start)  # Re-normalize

## Starting values of eta and nu4
eta.start = beta1.start
eta.start = eta.start + 1e-08  # Add 10^-8
eta.start = eta.start/sum(eta.start)  # Re-normalize

nu4.start = eta.start
nu4.start = nu4.start + 1e-08  # Add 10^-8
nu4.start = nu4.start/sum(nu4.start)  # Re-normalize


### RUNNING GDP

## Initialize storage

# Latent indicators
z1_samples <- list()
z2_samples <- list()
z3_samples <- list()
z4_samples <- list()
z5_samples <- list()
z6_samples <- list()
z7_samples <- list()
z8_samples <- list()
# Weights
beta1_samples <- list()
beta2_samples <- list()
beta3_samples <- list()
beta4_samples <- list()
beta5_samples <- list()
beta6_samples <- list()
beta7_samples <- list()
beta8_samples <- list()
nu1_samples <- list()
nu2_samples <- list()
nu3_samples <- list()
nu4_samples <- list()
eta_samples <- list()
# To save cluster specific sample sizes
n1_samples <- list()
n2_samples <- list()
n3_samples <- list()
n4_samples <- list()
n5_samples <- list()
n6_samples <- list()
n7_samples <- list()
n8_samples <- list()
# Storing the means and precisions of the normal distribution
mu_samples <- list()
tau_samples <- list()
# To store the concentration parameters of the GDP
alpha1_samples <- list()
alpha2_samples <- list()
alpha3_samples <- list()
alpha4_samples <- list()
alpha5_samples <- list()
alpha6_samples <- list()
alpha7_samples <- list()
alpha8_samples <- list()
# Initialize all the values at the starting values
beta1_samples[[1]] <- rep(1/L.max, L.max)
beta2_samples[[1]] <- rep(1/L.max, L.max)
beta3_samples[[1]] <- rep(1/L.max, L.max)
beta4_samples[[1]] <- rep(1/L.max, L.max)
beta5_samples[[1]] <- rep(1/L.max, L.max)
beta6_samples[[1]] <- rep(1/L.max, L.max)
beta7_samples[[1]] <- rep(1/L.max, L.max)
beta8_samples[[1]] <- rep(1/L.max, L.max)
nu1_samples[[1]] <- rep(1/L.max, L.max)
nu2_samples[[1]] <- rep(1/L.max, L.max)
nu3_samples[[1]] <- rep(1/L.max, L.max)
nu4_samples[[1]] <- rep(1/L.max, L.max)
eta_samples[[1]] <- rep(1/L.max, L.max)


mu_samples[[1]] <- mu.start
tau_samples[[1]] <- tau.start
z1_samples[[1]] <- z1.start
z2_samples[[1]] <- z2.start
z3_samples[[1]] <- z3.start
z4_samples[[1]] <- z4.start
z5_samples[[1]] <- z5.start
z6_samples[[1]] <- z6.start
z7_samples[[1]] <- z7.start
z8_samples[[1]] <- z8.start
alpha1_samples[[1]] <- alpha1.start
alpha2_samples[[1]] <- alpha2.start
alpha3_samples[[1]] <- alpha3.start
alpha4_samples[[1]] <- alpha4.start
alpha5_samples[[1]] <- alpha5.start
alpha6_samples[[1]] <- alpha6.start
alpha7_samples[[1]] <- alpha7.start
alpha8_samples[[1]] <- alpha8.start
# List to store the log likelihood values
log_likelihood1 <- list()
log_likelihood2 <- list()
log_likelihood3 <- list()
log_likelihood4 <- list()
log_likelihood5 <- list()
log_likelihood6 <- list()
log_likelihood7 <- list()
log_likelihood8 <- list()
# Rand index
rand_gdp_1 <- list()
rand_gdp_2 <- list()
rand_gdp_3 <- list()
rand_gdp_4 <- list()
rand_gdp_5 <- list()
rand_gdp_6 <- list()
rand_gdp_7 <- list()
rand_gdp_8 <- list()

# Running the Sampler

L.obs = L.max
gibbs_iterations = 25000

nu = nrow(X1)

W = diag(1, nrow = nu)  # Rate of Gamma prior
mu0 = as.matrix(rep(0, nu), ncol = 1)  # Prior mean of Normal
lambda = 0.01  # Prior precision of Normal

prior_prec = lambda  # Prior mean of Normal distribution of G0
prior_mean = mu0  # Prior precision of Normal distribution of G0

p = nrow(X1)

for (b in 2:gibbs_iterations) {
    # Printing the iterations
    if (b == 2) {
        cat(paste0("Iteration: ", (b - 1), "\n"))
    }
    if (b%%floor((10/100) * (gibbs_iterations + 1)) == 0) {
        cat(paste0("Iteration: ", b, "\n"))
    }
    nu = base::nrow(X1)
    W = diag(1, nrow = nu)
    # Sample the latent indicators for population 1
    sample = sample_z_c(data = t(X1), beta = beta1_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z1_samples[[b]] = sample  # Store the samples values of z1
    # Sample the latent indicators for population 2
    sample = sample_z_c(data = t(X2), beta = beta2_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z2_samples[[b]] = sample  # Store the samples values of z2
    # Sample the latent indicators for population 2
    sample = sample_z_c(data = t(X3), beta = beta3_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z3_samples[[b]] = sample  # Store the samples values of z3
    # Sample the latent indicators for population 4
    sample = sample_z_c(data = t(X4), beta = beta4_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z4_samples[[b]] = sample  # Store the samples values of z4
    # Sample the latent indicators for population 5
    sample = sample_z_c(data = t(X5), beta = beta5_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z5_samples[[b]] = sample  # Store the samples values of z5
    # Sample the latent indicators for population 6
    sample = sample_z_c(data = t(X6), beta = beta6_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z6_samples[[b]] = sample  # Store the samples values of z6
    # Sample the latent indicators for population 7
    sample = sample_z_c(data = t(X7), beta = beta7_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z7_samples[[b]] = sample  # Store the samples values of z7
    # Sample the latent indicators for population 8
    sample = sample_z_c(data = t(X8), beta = beta8_samples[[b - 1]], mu = mu_samples[[b -
        1]], tau = tau_samples[[b - 1]], L = L.obs)

    z8_samples[[b]] = sample  # Store the samples values of z8

    mu_tau_sample = sample_mu_tau_gdp_c(t(X1), t(X2), t(X3), t(X4), t(X5),
        t(X6), t(X7), t(X8), z1 = z1_samples[[b]], z2 = z2_samples[[b]],
        z3 = z3_samples[[b]], z4 = z4_samples[[b]], z5 = z5_samples[[b]],
        z6 = z6_samples[[b]], z7 = z7_samples[[b]], z8 = z8_samples[[b]],
        L = L.obs, nu = nu, W = W, prior_mean = mu0, prior_prec = lambda)

    # Store the samples of mu, tau, n1, n2, n3, n4, n5, n6, n7 and n8
    mu_samples[[b]] = t(mu_tau_sample$mu.draw)
    tau_samples[[b]] = mu_tau_sample$tau.draw

    n1_samples[[b]] = mu_tau_sample$n1_h
    n2_samples[[b]] = mu_tau_sample$n2_h
    n3_samples[[b]] = mu_tau_sample$n3_h
    n4_samples[[b]] = mu_tau_sample$n4_h
    n5_samples[[b]] = mu_tau_sample$n5_h
    n6_samples[[b]] = mu_tau_sample$n6_h
    n7_samples[[b]] = mu_tau_sample$n7_h
    n8_samples[[b]] = mu_tau_sample$n8_h

    init = beta1.start
    # pars should include beta2, beta3, beta4, nu1, nu2, nu3, eta,
    # alpha2, alpha3, alpha4 Run SALTSampler for Beta1
    beta1 <- RunMH_my(Target = Target_beta1, init = init, B = 1, h = rep(2,
        L.obs), pars = list(alpha2 = alpha2_samples[[b - 1]], alpha3 = alpha3_samples[[b -
        1]], alpha4 = alpha4_samples[[b - 1]], beta2 = beta2_samples[[b -
        1]], beta3 = beta3_samples[[b - 1]], beta4 = beta4_samples[[b -
        1]], nu1 = nu1_samples[[b - 1]], nu2 = nu2_samples[[b - 1]], nu3 = nu3_samples[[b -
        1]], eta = eta_samples[[b - 1]]), dat = n1_samples[[b]], concentration = rep(alpha1_samples[[b -
        1]], length(init)))

    beta1_samples[[b]] = as.numeric(beta1$S)
    beta1_samples[[b]][is.na(beta1_samples[[b]])] <- 0

    beta1_samples[[b]] = beta1_samples[[b]] + 1e-06
    beta1_samples[[b]] = beta1_samples[[b]]/sum(beta1_samples[[b]])
    # Draw samples from the beta2
    new_concentration = alpha2_samples[[b - 1]] * beta1_samples[[b]]
    # Now we only need to consider components 1, 2,...,L.obs
    beta2_samples[[b]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n2_samples[[b]] +
        new_concentration))
    # To avoid numerical issues I add a small quantity to each
    # coordinate and re-normalize
    beta2_samples[[b]][is.na(beta2_samples[[b]])] <- 0

    beta2_samples[[b]] = beta2_samples[[b]] + 1e-06
    beta2_samples[[b]] = beta2_samples[[b]]/sum(beta2_samples[[b]])  # Re-normalize

    # Draw samples from the beta3
    new_concentration = alpha3_samples[[b - 1]] * beta1_samples[[b]]
    # Now we only need to consider components 1, 2,...,L.obs
    beta3_samples[[b]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n3_samples[[b]] +
        new_concentration))
    # To avoid numerical issues I add a small quantity to each
    # coordinate and re-normalize
    beta3_samples[[b]][is.na(beta3_samples[[b]])] <- 0

    beta3_samples[[b]] = beta3_samples[[b]] + 1e-06
    beta3_samples[[b]] = beta3_samples[[b]]/sum(beta3_samples[[b]])  # Re-normalize

    # Draw samples from the beta4
    new_concentration = alpha4_samples[[b - 1]] * beta1_samples[[b]]
    # Now we only need to consider components 1, 2,...,L.obs
    beta4_samples[[b]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n4_samples[[b]] +
        new_concentration))
    # To avoid numerical issues I add a small quantity to each
    # coordinate and re-normalize
    beta4_samples[[b]][is.na(beta4_samples[[b]])] <- 0

    beta4_samples[[b]] = beta4_samples[[b]] + 1e-06
    beta4_samples[[b]] = beta4_samples[[b]]/sum(beta4_samples[[b]])  # Re-normalize

    # Next sample nu1 pars should be a named list with names beta1,
    # beta5, alpha2, alpha4 and alpha5 a corresponds to alpha2 +
    # alpha4
    init = nu1.start
    nu1 = GDPSamp::RunMH_my(Target = GDPSamp::Target_nu1, init = init,
        B = 1, h = rep(2, L.obs), pars = list(alpha2 = alpha2_samples[[b -
            1]], alpha4 = alpha4_samples[[b - 1]], alpha5 = alpha5_samples[[b -
            1]], beta1 = beta1_samples[[b]], beta5 = beta5_samples[[b -
            1]]), dat = NULL, concentration = rep((alpha2_samples[[b -
            1]] + alpha4_samples[[b - 1]]), L.obs))

    nu1_samples[[b]] = as.numeric(nu1$S)
    nu1_samples[[b]][is.na(nu1_samples[[b]])] <- 0
    nu1_samples[[b]] = nu1_samples[[b]] + 1e-06
    nu1_samples[[b]] = nu1_samples[[b]]/sum(nu1_samples[[b]])
    # Draw samples from the beta5
    new_concentration = alpha5_samples[[b - 1]] * nu1_samples[[b]]
    # Now we only need to consider components 1, 2,...,L.obs
    beta5_samples[[b]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n5_samples[[b]] +
        new_concentration))
    # To avoid numerical issues I add a small quantity to each
    # coordinate and re-normalize
    beta5_samples[[b]][is.na(beta5_samples[[b]])] <- 0
    beta5_samples[[b]] = beta5_samples[[b]] + 1e-06
    beta5_samples[[b]] = beta5_samples[[b]]/sum(beta5_samples[[b]])  # Re-normalize

    # Next sample nu2 pars should be a named list with names beta1,
    # beta6, alpha2, alpha3 and alpha6.  a corresponds to alpha2 +
    # alpha3
    init = nu2.start
    nu2 = GDPSamp::RunMH_my(Target = GDPSamp::Target_nu2, init = init,
        B = 1, h = rep(2, L.obs), pars = list(alpha2 = alpha2_samples[[b -
            1]], alpha3 = alpha3_samples[[b - 1]], alpha6 = alpha6_samples[[b -
            1]], beta1 = beta1_samples[[b]], beta6 = beta6_samples[[b -
            1]]), dat = NULL, concentration = rep((alpha2_samples[[b -
            1]] + alpha3_samples[[b - 1]]), L.obs))

    nu2_samples[[b]] = as.numeric(nu2$S)
    nu2_samples[[b]][is.na(nu2_samples[[b]])] <- 0

    nu2_samples[[b]] = nu2_samples[[b]] + 1e-06
    nu2_samples[[b]] = nu2_samples[[b]]/sum(nu2_samples[[b]])
    # Draw samples from the beta6
    new_concentration = alpha6_samples[[b - 1]] * nu2_samples[[b]]
    # Now we only need to consider components 1, 2,...,L.obs
    beta6_samples[[b]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n6_samples[[b]] +
        new_concentration))
    # To avoid numerical issues I add a small quantity to each
    # coordinate and re-normalize
    beta6_samples[[b]][is.na(beta6_samples[[b]])] <- 0
    beta6_samples[[b]] = beta6_samples[[b]] + 1e-06
    beta6_samples[[b]] = beta6_samples[[b]]/sum(beta6_samples[[b]])  # Re-normalize

    # Next sample nu3 pars should be a named list with names beta1,
    # beta7, alpha3, alpha4 and alpha7 a corresponds to alpha3 +
    # alpha4
    init = nu3.start
    nu3 = GDPSamp::RunMH_my(Target = GDPSamp::Target_nu3, init = nu3.start,
        B = 1, h = rep(2, L.obs), pars = list(alpha3 = alpha3_samples[[b -
            1]], alpha4 = alpha4_samples[[b - 1]], alpha7 = alpha7_samples[[b -
            1]], beta1 = beta1_samples[[b]], beta7 = beta7_samples[[b -
            1]]), dat = NULL, concentration = rep((alpha3_samples[[b -
            1]] + alpha4_samples[[b - 1]]), L.obs))

    nu3_samples[[b]] = as.numeric(nu3$S)

    nu3_samples[[b]][is.na(nu3_samples[[b]])] <- 0
    nu3_samples[[b]] = nu3_samples[[b]] + 1e-06
    nu3_samples[[b]] = nu3_samples[[b]]/sum(nu3_samples[[b]])
    # Draw samples from the beta7
    new_concentration = alpha7_samples[[b - 1]] * nu3_samples[[b]]
    # Now we only need to consider components 1, 2,...,L.obs
    beta7_samples[[b]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n7_samples[[b]] +
        new_concentration))
    # To avoid numerical issues I add a small quantity to each
    # coordinate and re-normalize
    beta7_samples[[b]][is.na(beta7_samples[[b]])] <- 0
    beta7_samples[[b]] = beta7_samples[[b]] + 1e-06
    beta7_samples[[b]] = beta7_samples[[b]]/sum(beta7_samples[[b]])  # Re-normalize

    # Next sample eta pars should be a named list with names beta1,
    # nu4, alpha2, alpha3, alpha4, alpha5, alpha6, alpha7 a
    # corresponds to (alpha2 + alpha3 + alpha4)
    init = eta.start
    eta = GDPSamp::RunMH_my(Target = GDPSamp::Target_eta, init = init,
        B = 1, h = rep(2, L.obs), pars = list(alpha2 = alpha2_samples[[b -
            1]], alpha3 = alpha3_samples[[b - 1]], alpha4 = alpha4_samples[[b -
            1]], alpha5 = alpha5_samples[[b - 1]], alpha6 = alpha6_samples[[b -
            1]], alpha7 = alpha7_samples[[b - 1]], beta1 = beta1_samples[[b]],
            nu4 = nu4_samples[[b - 1]]), dat = NULL, concentration = rep((alpha2_samples[[b -
            1]] + alpha3_samples[[b - 1]] + alpha4_samples[[b - 1]]), L.obs))

    eta_samples[[b]] = as.numeric(eta$S)
    eta_samples[[b]][is.na(eta_samples[[b]])] <- 0

    eta_samples[[b]] = eta_samples[[b]] + 1e-06
    eta_samples[[b]] = eta_samples[[b]]/sum(eta_samples[[b]])
    # Next sample nu4 pars should be a named list with names beta8,
    # eta, alpha5, alpha6, alpha7 and alpha8.  a corresponds to
    # alpha5 + alpha6 + alpha7
    init = nu4.start
    nu4 = GDPSamp::RunMH_my(Target = GDPSamp::Target_nu4, init = init,
        B = 1, h = rep(2, L.obs), pars = list(alpha5 = alpha5_samples[[b -
            1]], alpha6 = alpha6_samples[[b - 1]], alpha7 = alpha7_samples[[b -
            1]], alpha8 = alpha8_samples[[b - 1]], beta8 = beta8_samples[[b -
            1]], eta = eta_samples[[b]]), dat = NULL, concentration = rep((alpha5_samples[[b -
            1]] + alpha6_samples[[b - 1]] + alpha7_samples[[b - 1]]), L.obs))

    nu4_samples[[b]] = as.numeric(nu4$S)
    nu4_samples[[b]][is.na(nu4_samples[[b]])] <- 0
    nu4_samples[[b]] = nu4_samples[[b]] + 1e-06
    nu4_samples[[b]] = nu4_samples[[b]]/sum(nu4_samples[[b]])
    # Draw samples from the beta8
    new_concentration = alpha8_samples[[b - 1]] * nu4_samples[[b]]
    # Now we only need to consider components 1, 2,...,L.obs
    beta8_samples[[b]] = as.numeric(extraDistr::rdirichlet(n = 1, alpha = n8_samples[[b]] +
        new_concentration))
    # To avoid numerical issues I add a small quantity to each
    # coordinate and re-normalize
    beta8_samples[[b]][is.na(beta8_samples[[b]])] <- 0
    beta8_samples[[b]] = beta8_samples[[b]] + 1e-06
    beta8_samples[[b]] = beta8_samples[[b]]/sum(beta8_samples[[b]])  # Re-normalize

    # Sample the alphas pars argument should be a named list with
    # alpha0, alpha2, alpha3, alpha4 and beta1
    alpha1_samples[[b]] = GDPSamp::sample_alpha1(ycurrent = alpha1_samples[[b -
        1]], pars = list(beta1 = beta1_samples[[b]], alpha0 = alpha0, alpha2 = alpha2_samples[[b -
        1]], alpha3 = alpha3_samples[[b - 1]], alpha4 = alpha4_samples[[b -
        1]]))$out
    # @pars argument should be a named list with alpha1, alpha3,
    # alpha4, alpha5, alpha6, beta1, beta2, nu1, nu2, eta
    alpha2_samples[[b]] = GDPSamp::sample_alpha2(ycurrent = alpha2_samples[[b -
        1]], pars = list(alpha1 = alpha1_samples[[b]], alpha3 = alpha3_samples[[b -
        1]], alpha4 = alpha4_samples[[b - 1]], alpha5 = alpha5_samples[[b -
        1]], alpha6 = alpha6_samples[[b - 1]], beta1 = beta1_samples[[b]],
        beta2 = beta2_samples[[b]], nu1 = nu1_samples[[b]], nu2 = nu2_samples[[b]],
        eta = eta_samples[[b]]))$out

    # @pars argument should be a named list with alpha1, alpha2,
    # alpha4, alpha6, alpha7, beta1, beta3, nu2, nu3, eta
    alpha3_samples[[b]] = GDPSamp::sample_alpha3(ycurrent = alpha3_samples[[b -
        1]], pars = list(alpha1 = alpha1_samples[[b]], alpha2 = alpha2_samples[[b]],
        alpha4 = alpha4_samples[[b - 1]], alpha6 = alpha6_samples[[b -
            1]], alpha7 = alpha7_samples[[b - 1]], beta1 = beta1_samples[[b]],
        beta3 = beta3_samples[[b]], nu2 = nu2_samples[[b]], nu3 = nu3_samples[[b]],
        eta = eta_samples[[b]]))$out

    # @pars argument should be a named list with alpha1, alpha2,
    # alpha3, alpha5, alpha7, beta1, beta4, nu1, nu3, eta
    alpha4_samples[[b]] = GDPSamp::sample_alpha4(ycurrent = alpha4_samples[[b -
        1]], pars = list(alpha1 = alpha1_samples[[b]], alpha2 = alpha2_samples[[b]],
        alpha3 = alpha3_samples[[b]], alpha5 = alpha5_samples[[b - 1]],
        alpha7 = alpha7_samples[[b - 1]], beta1 = beta1_samples[[b]], beta4 = beta4_samples[[b]],
        nu1 = nu1_samples[[b]], nu3 = nu3_samples[[b]], eta = eta_samples[[b]]))$out

    # @pars argument should be a named list with alpha2, alpha4,
    # alpha6, alpha7, alpha8 beta5, nu1, nu4, eta
    alpha5_samples[[b]] = GDPSamp::sample_alpha5(ycurrent = alpha5_samples[[b -
        1]], pars = list(alpha2 = alpha2_samples[[b]], alpha4 = alpha4_samples[[b]],
        alpha6 = alpha6_samples[[b - 1]], alpha7 = alpha7_samples[[b -
            1]], alpha8 = alpha8_samples[[b - 1]], beta5 = beta5_samples[[b]],
        nu1 = nu1_samples[[b]], nu4 = nu4_samples[[b]], eta = eta_samples[[b]]))$out

    # @pars argument should be a named list with alpha2, alpha3,
    # alpha5, alpha7, alpha8 beta6, nu2, nu4, eta
    alpha6_samples[[b]] = GDPSamp::sample_alpha6(ycurrent = alpha6_samples[[b -
        1]], pars = list(alpha2 = alpha2_samples[[b]], alpha3 = alpha3_samples[[b]],
        alpha5 = alpha5_samples[[b]], alpha7 = alpha7_samples[[b - 1]],
        alpha8 = alpha8_samples[[b - 1]], beta6 = beta6_samples[[b]], nu2 = nu2_samples[[b]],
        nu4 = nu4_samples[[b]], eta = eta_samples[[b]]))$out

    # @pars argument should be a named list with alpha3, alpha4,
    # alpha5, alpha6, alpha8 beta7, nu3, nu4, eta
    alpha7_samples[[b]] = GDPSamp::sample_alpha7(ycurrent = alpha7_samples[[b -
        1]], pars = list(alpha3 = alpha3_samples[[b]], alpha4 = alpha4_samples[[b]],
        alpha5 = alpha5_samples[[b]], alpha6 = alpha6_samples[[b]], alpha8 = alpha8_samples[[b -
            1]], beta7 = beta7_samples[[b]], nu3 = nu2_samples[[b]], nu4 = nu4_samples[[b]],
        eta = eta_samples[[b]]))$out
    # @pars argument should be a named list with alpha5, alpha6,
    # alpha7 beta8, nu4
    alpha8_samples[[b]] = GDPSamp::sample_alpha8(ycurrent = alpha8_samples[[b -
        1]], pars = list(alpha5 = alpha5_samples[[b]], alpha6 = alpha6_samples[[b]],
        alpha7 = alpha7_samples[[b]], beta8 = beta8_samples[[b]], nu4 = nu4_samples[[b]]))$out

    # Calculate Rand-index at every interation
    rand_gdp_1[[b]] = aricode::ARI(z1_samples[[b]], z1.true)
    rand_gdp_2[[b]] = aricode::ARI(z2_samples[[b]], z2.true)
    rand_gdp_3[[b]] = aricode::ARI(z3_samples[[b]], z3.true)
    rand_gdp_4[[b]] = aricode::ARI(z4_samples[[b]], z4.true)
    rand_gdp_5[[b]] = aricode::ARI(z5_samples[[b]], z5.true)
    rand_gdp_6[[b]] = aricode::ARI(z6_samples[[b]], z6.true)
    rand_gdp_7[[b]] = aricode::ARI(z7_samples[[b]], z7.true)
    rand_gdp_8[[b]] = aricode::ARI(z8_samples[[b]], z8.true)
    ## Log likelihoods to monitor convergence of the chain Calculate
    ## the log likelihood for group 1 for every iteration
    log_likelihood1[[b]] = log_l(data = t(X1), beta = beta1_samples[[b]],
        z = z1_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
    # Calculate the log likelihood for group 2 for every iteration
    log_likelihood2[[b]] = log_l(data = t(X2), beta = beta2_samples[[b]],
        z = z2_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
    # Calculate the log likelihood for group 3 for every iteration
    log_likelihood3[[b]] = log_l(data = t(X3), beta = beta3_samples[[b]],
        z = z3_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
    # Calculate the log likelihood for group 4 for every iteration
    log_likelihood4[[b]] = log_l(data = t(X4), beta = beta4_samples[[b]],
        z = z4_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
    # Calculate the log likelihood for group 5 for every iteration
    log_likelihood5[[b]] = log_l(data = t(X5), beta = beta5_samples[[b]],
        z = z5_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
    # Calculate the log likelihood for group 6 for every iteration
    log_likelihood6[[b]] = log_l(data = t(X6), beta = beta6_samples[[b]],
        z = z6_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
    # Calculate the log likelihood for group 7 for every iteration
    log_likelihood7[[b]] = log_l(data = t(X7), beta = beta7_samples[[b]],
        z = z7_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])
    # Calculate the log likelihood for group 8 for every iteration
    log_likelihood8[[b]] = log_l(data = t(X8), beta = beta8_samples[[b]],
        z = z8_samples[[b]], mu = mu_samples[[b]], tau = tau_samples[[b]])

}

gibbs_burn <- 10000
thin = 10
thin_samples <- seq(from = (gibbs_burn + 1), to = gibbs_iterations, by = thin)


## Clustering by DAHL's method

z1_samples.post.burn <- z1_samples[thin_samples]
z2_samples.post.burn <- z2_samples[thin_samples]
z3_samples.post.burn <- z3_samples[thin_samples]
z4_samples.post.burn <- z4_samples[thin_samples]
z5_samples.post.burn <- z5_samples[thin_samples]
z6_samples.post.burn <- z6_samples[thin_samples]
z7_samples.post.burn <- z7_samples[thin_samples]
z8_samples.post.burn <- z8_samples[thin_samples]

# Calculate the membership matrix for Group 1
M1 <- lapply(z1_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M1.mean <- Reduce("+", M1)/length(z1_samples.post.burn)
# Calculate the Frobenius norm of the differences
M1.Frobenius <- sapply(M1, function(x, av) sum((x - av)^2), av = M1.mean)
# Find out the minimums
k1.min <- which.min(M1.Frobenius)
# Remove the membership matrix to free up memory
rm(M1)
rm(M1.mean)

# Calculate the membership matrix for Group 2
M2 <- lapply(z2_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M2.mean <- Reduce("+", M2)/length(z2_samples.post.burn)
# Calculate the Frobenius norm of the differences
M2.Frobenius <- sapply(M2, function(x, av) sum((x - av)^2), av = M2.mean)
# Find out the minimums
k2.min <- which.min(M2.Frobenius)
# Remove the membership matrix to free up memory
rm(M2)
rm(M2.mean)

# Calculate the membership matrix for Group 3
M3 <- lapply(z3_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M3.mean <- Reduce("+", M3)/length(z3_samples.post.burn)
# Calculate the Frobenius norm of the differences
M3.Frobenius <- sapply(M3, function(x, av) sum((x - av)^2), av = M3.mean)
# Find out the minimums
k3.min <- which.min(M3.Frobenius)
# Remove the membership matrix to free up memory
rm(M3)
rm(M3.mean)

# Calculate the membership matrix for Group 4
M4 <- lapply(z4_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M4.mean <- Reduce("+", M4)/length(z4_samples.post.burn)
# Calculate the Frobenius norm of the differences
M4.Frobenius <- sapply(M4, function(x, av) sum((x - av)^2), av = M4.mean)
# Find out the minimums
k4.min <- which.min(M4.Frobenius)
# Remove the membership matrix to free up memory
rm(M4)
rm(M4.mean)

# Calculate the membership matrix for Group 5
M5 <- lapply(z5_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M5.mean <- Reduce("+", M5)/length(z5_samples.post.burn)
# Calculate the Frobenius norm of the differences
M5.Frobenius <- sapply(M5, function(x, av) sum((x - av)^2), av = M5.mean)
# Find out the minimums
k5.min <- which.min(M5.Frobenius)
# Remove the membership matrix to free up memory
rm(M5)
rm(M5.mean)


# Calculate the membership matrix for Group 6
M6 <- lapply(z6_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M6.mean <- Reduce("+", M6)/length(z6_samples.post.burn)
# Calculate the Frobenius norm of the differences
M6.Frobenius <- sapply(M6, function(x, av) sum((x - av)^2), av = M6.mean)
# Find out the minimums
k6.min <- which.min(M6.Frobenius)
# Remove the membership matrix to free up memory
rm(M6)
rm(M6.mean)

# Calculate the membership matrix for Group 7
M7 <- lapply(z7_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M7.mean <- Reduce("+", M7)/length(z7_samples.post.burn)
# Calculate the Frobenius norm of the differences
M7.Frobenius <- sapply(M7, function(x, av) sum((x - av)^2), av = M7.mean)
# Find out the minimums
k7.min <- which.min(M7.Frobenius)
# Remove the membership matrix to free up memory
rm(M7)
rm(M7.mean)

# Calculate the membership matrix for Group 8
M8 <- lapply(z8_samples.post.burn, function(x) {
    clusterAssign <- x
    1 * outer(clusterAssign, clusterAssign, FUN = "==")
})
# Mean membership matrix
M8.mean <- Reduce("+", M8)/length(z8_samples.post.burn)
# Calculate the Frobenius norm of the differences
M8.Frobenius <- sapply(M8, function(x, av) sum((x - av)^2), av = M8.mean)
# Find out the minimums
k8.min <- which.min(M8.Frobenius)
# Remove the membership matrix to free up memory
rm(M8)
rm(M8.mean)

# POST PROCESSING FOR CLUSTERING
singleton1 = as.numeric(names(which(table(z1_samples.post.burn[[k1.min]]) <=
    0.02 * n1)))
singleton.index1 = which(z1_samples.post.burn[[k1.min]] %in% singleton1)
if (length(singleton.index1) > 0) {
    z1.estimated = z1_samples.post.burn[[k1.min]][-singleton.index1]
    z1.true.subset = z1.true[-singleton.index1]
} else {
    z1.estimated = z1_samples.post.burn[[k1.min]]
    z1.true.subset = z1.true
}

rand.gdp1 = aricode::ARI(z1.estimated, z1.true.subset)

singleton2 = as.numeric(names(which(table(z2_samples.post.burn[[k2.min]]) <=
    0.02 * n2)))
singleton.index2 = which(z2_samples.post.burn[[k2.min]] %in% singleton2)
if (length(singleton.index2) > 0) {
    z2.estimated = z2_samples.post.burn[[k2.min]][-singleton.index2]
    z2.true.subset = z2.true[-singleton.index2]
} else {
    z2.estimated = z2_samples.post.burn[[k2.min]]
    z2.true.subset = z2.true
}
rand.gdp2 = aricode::ARI(z2.estimated, z2.true.subset)

singleton3 = as.numeric(names(which(table(z3_samples.post.burn[[k3.min]]) <=
    0.02 * n3)))
singleton.index3 = which(z3_samples.post.burn[[k3.min]] %in% singleton3)
if (length(singleton.index3) > 0) {
    z3.estimated = z3_samples.post.burn[[k3.min]][-singleton.index3]
    z3.true.subset = z3.true[-singleton.index3]
} else {
    z3.estimated = z3_samples.post.burn[[k3.min]]
    z3.true.subset = z3.true
}
rand.gdp3 = aricode::ARI(z3.estimated, z3.true.subset)

singleton4 = as.numeric(names(which(table(z4_samples.post.burn[[k4.min]]) <=
    0.02 * n4)))
singleton.index4 = which(z4_samples.post.burn[[k4.min]] %in% singleton4)
if (length(singleton.index4) > 0) {
    z4.estimated = z4_samples.post.burn[[k4.min]][-singleton.index4]
    z4.true.subset = z4.true[-singleton.index4]
} else {
    z4.estimated = z4_samples.post.burn[[k4.min]]
    z4.true.subset = z4.true
}
rand.gdp4 = aricode::ARI(z4.estimated, z4.true.subset)

singleton5 = as.numeric(names(which(table(z5_samples.post.burn[[k5.min]]) <=
    0.02 * n5)))
singleton.index5 = which(z5_samples.post.burn[[k5.min]] %in% singleton5)
if (length(singleton.index5) > 0) {
    z5.estimated = z5_samples.post.burn[[k5.min]][-singleton.index5]
    z5.true.subset = z5.true[-singleton.index5]
} else {
    z5.estimated = z5_samples.post.burn[[k5.min]]
    z5.true.subset = z5.true
}
rand.gdp5 = aricode::ARI(z5.estimated, z5.true.subset)

singleton6 = as.numeric(names(which(table(z6_samples.post.burn[[k6.min]]) <=
    0.02 * n6)))
singleton.index6 = which(z6_samples.post.burn[[k6.min]] %in% singleton6)
if (length(singleton.index6) > 0) {
    z6.estimated = z6_samples.post.burn[[k6.min]][-singleton.index6]
    z6.true.subset = z6.true[-singleton.index6]
} else {
    z6.estimated = z6_samples.post.burn[[k6.min]]
    z6.true.subset = z6.true
}
rand.gdp6 = aricode::ARI(z6.estimated, z6.true.subset)

singleton7 = as.numeric(names(which(table(z7_samples.post.burn[[k7.min]]) <=
    0.02 * n7)))
singleton.index7 = which(z7_samples.post.burn[[k7.min]] %in% singleton7)
if (length(singleton.index7) > 0) {
    z7.estimated = z7_samples.post.burn[[k7.min]][-singleton.index7]
    z7.true.subset = z7.true[-singleton.index7]
} else {
    z7.estimated = z7_samples.post.burn[[k7.min]]
    z7.true.subset = z7.true
}
rand.gdp7 = aricode::ARI(z7.estimated, z7.true.subset)

singleton8 = as.numeric(names(which(table(z8_samples.post.burn[[k8.min]]) <=
    0.02 * n8)))
singleton.index8 = which(z8_samples.post.burn[[k8.min]] %in% singleton8)
if (length(singleton.index8) > 0) {
    z8.estimated = z8_samples.post.burn[[k8.min]][-singleton.index8]
    z8.true.subset = z8.true[-singleton.index8]
} else {
    z8.estimated = z8_samples.post.burn[[k8.min]]
    z8.true.subset = z8.true
}
rand.gdp8 = aricode::ARI(z8.estimated, z8.true.subset)


log_like1 = unlist(log_likelihood1[thin_samples])
log_like2 = unlist(log_likelihood2[thin_samples])
log_like3 = unlist(log_likelihood3[thin_samples])
log_like4 = unlist(log_likelihood4[thin_samples])
log_like5 = unlist(log_likelihood5[thin_samples])
log_like6 = unlist(log_likelihood6[thin_samples])
log_like7 = unlist(log_likelihood7[thin_samples])
log_like8 = unlist(log_likelihood8[thin_samples])

ess.gdp1 = coda::effectiveSize(log_like1)
ess.gdp2 = coda::effectiveSize(log_like2)
ess.gdp3 = coda::effectiveSize(log_like3)
ess.gdp4 = coda::effectiveSize(log_like4)
ess.gdp5 = coda::effectiveSize(log_like5)
ess.gdp6 = coda::effectiveSize(log_like6)
ess.gdp7 = coda::effectiveSize(log_like7)
ess.gdp8 = coda::effectiveSize(log_like8)
```

Plot of the log-likelihood and clustering plots can be obtained as
below:

``` r
# PLOT OF CLUSTERING

if(!require("fossil")) install.packages("fossil"); library(fossil)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("latex2exp")) install.packages("latex2exp"); library(latex2exp)
library(aricode)

myvalues = c("1" = "#F8766D",
             "2" = "#00BA38",
             "3" = "#619CFF", 
             "4" = "blueviolet",
             "5" = "cyan4",
             "6" = "#E6AB02",
             "7" = "#E36EF6",
             "8" = "bisque4",
             "9" = "coral4",
             "10" = "darkslateblue")

factors = factor(z1_samples.post.burn[[k1.min]])

cluster <- data.frame(x1 = X1[1, ],
                      x2 = X1[2, ], 
                      cluster = factors)

p1.leg = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + 
  geom_point(size = 2) + 
  labs(title = "Clustering for group 1",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp1, 4)), 
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

p1 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + 
  geom_point() + 
  labs(title = "Clustering for group 1",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp1, 4)), 
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))


factors = factor(z2_samples.post.burn[[k2.min]])

cluster <- data.frame(x1 = X2[1, ],
                      x2 = X2[2, ], 
                      cluster = factors)

p2 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) +
  geom_point() + 
  labs(title = "Clustering for group 2",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp2, 4)), 
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) +
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

factors = factor(z3_samples.post.burn[[k3.min]])

cluster <- data.frame(x1 = X3[1, ],
                      x2 = X3[2, ], 
                      cluster = factors)


p3 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + 
  geom_point() + 
  labs(title = "Clustering for group 3",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp3, 4)),
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

factors = factor(z4_samples.post.burn[[k4.min]])

cluster <- data.frame(x1 = X4[1, ],
                      x2 = X4[2, ], 
                      cluster = factors)

p4 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + 
  geom_point() + 
  labs(title = "Clustering for group 4",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp4, 4)), 
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

factors = factor(z5_samples.post.burn[[k5.min]])

cluster <- data.frame(x1 = X5[1, ],
                      x2 = X5[2, ], 
                      cluster = factors)

p5 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + 
  geom_point() + 
  labs(title = "Clustering for group 5",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp5, 4)),
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

factors = factor(z6_samples.post.burn[[k6.min]])

cluster <- data.frame(x1 = X6[1, ],
                      x2 = X6[2, ], 
                      cluster = factors)

p6 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + 
  geom_point() +
  labs(title = "Clustering for group 6",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp6, 4)),
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))


factors = factor(z7_samples.post.burn[[k7.min]])

cluster <- data.frame(x1 = X7[1, ],
                      x2 = X7[2, ], 
                      cluster = factors)

p7 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) +
  geom_point() + 
  labs(title = "Clustering for group 7",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp7, 4)), 
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues)+ theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

factors = factor(z8_samples.post.burn[[k8.min]])

cluster <- data.frame(x1 = X8[1, ],
                      x2 = X8[2, ], 
                      cluster = factors)

p8 = cluster %>% ggplot(aes(x = x1, y = x2, col = cluster)) + 
  geom_point() + 
  labs(title = "Clustering for group 8",
       subtitle = paste0("Adjusted Rand index = ", round(rand.gdp8, 4)), 
       x = TeX(r'($X_1$)'), y = TeX(r'($X_2$)')) + 
  ylim(c(min(mean.true) - 4, max(mean.true) + 4)) +
  xlim(c(min(mean.true) - 4, max(mean.true) + 4)) + 
  scale_color_manual(values = myvalues) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
if(!require("cowplot")) install.packages("cowplot");library(cowplot)
p.len = get_legend(p1.leg)
gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p.len)


# PLOT OF TRACEPLOT OF LOG-LIKELIHOOD

l1 = data.frame(x = 1:length(log_like1), y = log_like1) %>% 
  ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 1", x = "Iteration", y = "") +
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l2 = data.frame(x = 1:length(log_like2), y = log_like2) %>%
  ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 2", x = "Iteration", y = "") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l3 = data.frame(x = 1:length(log_like3), y = log_like3) %>% 
  ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 3", x = "Iteration", y = "") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l4 = data.frame(x = 1:length(log_like4), y = log_like4) %>% 
  ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 4", x = "Iteration", y = "") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l5 = data.frame(x = 1:length(log_like5), y = log_like5) %>% 
  ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 5", x = "Iteration", y = "") +
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l6 = data.frame(x = 1:length(log_like6), y = log_like6) %>% 
  ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 6", x = "Iteration", y = "") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l7 = data.frame(x = 1:length(log_like7), y = log_like7) %>% 
  ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 7", x = "Iteration", y = "") + 
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))
l8 = data.frame(x = 1:length(log_like8), y = log_like8) 
%>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(title = "Log-likelihood for group 8", x = "Iteration", y = "") +
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
gridExtra::grid.arrange(l1, l2, l3, l4, l5, l6, l7, l8)
```
