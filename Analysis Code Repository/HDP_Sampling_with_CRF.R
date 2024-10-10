rm(list = ls())
# Generate data 
J = 3; L = 4; n = rep(300, J)
xi = 0; 
tau = 1

# # True cluster means
phi.true = c(-6, -2, 2, 6)
# True cluster means
# phi.true = c(-3, -1, 1, 3)

# True mixture weights
Pi1 = c(0.5, 0.5, 0, 0); Pi2 = rep(1/L, L); Pi3 = c(0, 0.1, 0.6, 0.3);
Pi.true = rbind(Pi1, Pi2, Pi3)

# True cluster labels
true.Z = vector(mode = "list", length = J)

# Observations in each group is stored as a list
# x[[j]] denotes observations in the jth group
x = vector(mode = "list", length = J)

set.seed(2024)
for(j in 1:J){
  true.Z[[j]] = sample(1:L, size = n[j], prob = Pi.true[j, ], replace = TRUE)
  x[[j]] = sapply(1:n[j], function(i) 
    rnorm(n = 1, mean = phi.true[true.Z[[j]][i]], sd = sqrt(1/tau)))
}

source("~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/CRF.R")
xmin = min(unlist(x)) - 1 ; xmax = max(unlist(x)) + 1
y.grid = seq(xmin, xmax, length = 100)
# Number of MCMC samples and Burn in period.
M = 10000; M.burn = 5000; thin = 5

# hyperparameter specifications
G = 1; a = 0.1; b = 0.1; lambda = 0.01

CRF_HDP = CRF_gibbs(x = x, y.grid = y.grid, K.init = 1, gam = G, a = a, b = b, xi = xi, tau = tau, lambda = lambda, Burn.in = M.burn, M = M, thin = thin)

source("/Users/arhitchakrabarti/Downloads/blockedHDP/postestimates.R")
# get BGS clusters using Dahl's method
Z.hat = getDahl(CRF_HDP)

# Adjusted Rand indices corresponding to group specific cluster labels
ARI = sapply(seq_len(J), function(j) mcclust::arandi(Z.hat[[j]], true.Z[[j]]))
ARI

# Adjusted Rand index corresponding to global cluster labels
ARI.global = mcclust::arandi(unlist(Z.hat), unlist(true.Z))
ARI.global



################################################################################
## Running a Vanilla DP to get an estimate of beta
################################################################################
library(Rcpp)
source("HDP_Functions_all.R")
# Prior specifications
p = nrow(matrix(x[[1]], nrow = 1))
prior_mean = as.matrix(rep(0, p), ncol = 1) # Prior mean of Normal
lambda.inv = 0.01
lambda = 1/lambda.inv
prior_sigma = diag(lambda, p)
mu.start.list = list()

L.max = 10

for(i in 1:L.max){
  mu.start.list[[i]] = matrix(mvtnorm::rmvnorm(n = 1, mean = as.numeric(prior_mean), sigma = prior_sigma), ncol = 1)
}

beta.start = rep(1/L.max, L.max)
library(extraDistr)
z.start = rcat(n = sum(n), prob = beta.start)
mu.start <- t(sapply(1:L.max, function(j){mu.start.list[[j]]}))
alpha.start <- 1
a0 = 0.1; b0 = 0.1
# Prior hyper-parameters
prior = list(prior_mean = prior_mean, prior_sigma = prior_sigma, a0 = a0, b0 = b0)
# Initial values of different parameters
init = list(mu.start = mu.start, beta.start = beta.start, z.start = z.start, alpha.start = alpha.start)
num_iterations <- 50000
burn = 25000
thin = 25


gibbs_out <- gibbs_N_dp_c(init = init, data = matrix(unlist(x)), num_iterations = num_iterations, burn = burn, thin = thin, L.max = L.max, prior = prior, Sigma.true = matrix(1/tau, 1, 1))
library(latex2exp)
library(tidyverse)
# DAHL's METHOD
dp_nw.dahl = getDahl(gibbs_out$z_samples, unlist(true.Z))
cluster <- data.frame(x = 1:sum(n),
                      y = unlist(x), 
                      z = factor(gibbs_out$z_samples[[dp_nw.dahl$DahlIndex]]))
plots <- cluster %>% ggplot(aes(x = x, y = y, col = z)) + geom_point(size = 2) + labs(subtitle = paste0("Adjusted Rand index = ", round(dp_nw.dahl$adj.rand, 4))) + labs(x = "Index", y = TeX("$x$")) + 
  theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  )

plots


DahlIndex = dp_nw.dahl$DahlIndex
z.dahl = gibbs_out$z_samples[[DahlIndex]]
# Get the estimate of beta1 by running a Naive DP on group 1 data
beta.start.small = gibbs_out$beta_samp_post[[DahlIndex]][unique(z.dahl)]

beta.start.small = beta.start.small/sum(beta.start.small)

# But for the next step again I would require beta starting value to be of length L.max
# So I add a very small quantity to the remaining co-ordinates and re-normlaize
beta.start.long = c(beta.start.small, rep(1e-8, abs(L.max - length(beta.start.small))))
beta.start_HDP = beta.start.long/sum(beta.start.long) # This becomes the starting value of beta1 in the Gibbs iteration


library(Rcpp)
################################################################################
source("GDP_chain_functions_all_functions.R")
################################################################################
num_iterations = 50000
#####################################################################################################
# RUNNING HDP
#####################################################################################################
#####################################################################################################
## Initialize storage
#####################################################################################################
# Latent indicators
z_samples <- replicate(num_iterations, list(replicate(J, list())))
# Weights
beta_samples <- list()
pi_samples <- replicate(num_iterations, list(replicate(J, list())))
# Concentration parameters
alpha_samples <- list()
gamma_samples <- list()
# To save cluster specific sample sizes
n_samples <- replicate(num_iterations, list(replicate(J, list())))
# Storing the means and precisions of the normal distribution
mu_samples <- replicate(num_iterations, list())
L.max = 10

prior_mean = matrix(xi, 1, 1)
prior_sigma = matrix(lambda, 1, 1)
mu.start_HDP.list = list()

for(i in 1:L.max){
  mu.start_HDP.list[[i]] = matrix(mvtnorm::rmvnorm(n = 1, mean = as.numeric(prior_mean), sigma = prior_sigma), ncol = 1)
}

mu.start_HDP <- t(sapply(1:L.max, function(j){mu.start_HDP.list[[j]]}))
beta_samples[[1]] <- beta.start_HDP

for(j in 1:J){
  pi_samples[[1]][[j]] <- rep(1/L.max, L.max)
  z_samples[[1]][[j]] <- sample(1:L.max, size = n[j], prob = pi_samples[[1]][[j]], replace = TRUE)
}

alpha_samples[[1]] <- 1
gamma_samples[[1]] <- 1

mu_samples[[1]] <- mu.start_HDP
sigma.true = matrix((1/tau), 1, 1)
################################################################################
# Running the Sampler
################################################################################
if(!require("SALTSampler")) install.packages("SALTSampler"); library(SALTSampler)
if(!require("Rcpp")) install.packages("Rcpp"); library(Rcpp)
if(!require("fossil")) install.packages("fossil"); library(fossil)
library(extraDistr)
source("HDP_Functions.R")
source("RunMH_my.R")
a0 = 0.1; b0 = 0.1; B_SALT = 1 # Number of MH Samples in the SALTSampler

for(b in 2:num_iterations){
  # b = 2
  # Printing the iterations
  if(b == 2){
    cat(paste0("Iteration: ", (b-1), "\n"))
  }
  if(b %% floor((10/100)*(num_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", b, "\n"))
  }
  
  for(j in 1:J){
    # Sample the latent indicators for population j
    sample = sample_z_known_sigma_c(data = matrix(x[[j]], ncol = 1), beta = pi_samples[[b-1]][[j]], mu = matrix(mu_samples[[b-1]]), Sigma = sigma.true, L = L.max)
    
    z_samples[[b]][[j]] = sample # Store the samples values of zj
  }
  
  
  out = sample_mu_gdp_c(data = lapply(x, matrix, nrow = 1), z = z_samples[[b]], Sigma = sigma.true, L = L.max, n = n, prior_mean = prior_mean,  prior_sigma = prior_sigma)
  # Store the samples of mu, and n_j
  mu_samples[[b]] =  out$mu_draw;
  
  for(j in 1:J){
    n_samples[[b]][[j]] = out$n_h[, 1][[j]]
  }
  
  
  init = beta.start_HDP
  # init = beta_samples[[b-1]]
  
  # Run SALTSampler for Beta
  ## NOTE: pars should be a named list with names pii which is a list and alpha. 
  # a correspons to gamma. dat is NULL here
  beta <- RunMH_my(Target = Target_beta, init = init, B = B_SALT,
                   h = rep(2, L.max), 
                   
                   pars = list(pii = pi_samples[[b-1]],
                               alpha = alpha_samples[[b-1]]),
                   dat = NULL, 
                   concentration = rep(gamma_samples[[b-1]], length(init)))
  
  beta_samples[[b]] = as.numeric(beta$S[B_SALT, ])
  beta_samples[[b]] =  beta_samples[[b]] + 1e-6
  beta_samples[[b]] =  beta_samples[[b]]/sum( beta_samples[[b]] )
  
  for(j in 1:J){
    new_concentration = alpha_samples[[b-1]] * beta_samples[[b]]
    pi_samples[[b]][[j]] = as.numeric(rdirichlet(n = 1, alpha = n_samples[[b]][[j]] + new_concentration))
    
    pi_samples[[b]][[j]] = pi_samples[[b]][[j]] + 1e-6
    pi_samples[[b]][[j]] = pi_samples[[b]][[j]]/sum(pi_samples[[b]][[j]])
  }
  
  alpha_samples[[b]] = sample_alpha(ycurrent = alpha_samples[[b-1]], 
                                    pars = list(a0 = a0,
                                                b0 = b0,
                                                pii = pi_samples[[b]],
                                                beta = beta_samples[[b]]) )$out
  
  gamma_samples[[b]] = sample_gamma(ycurrent = gamma_samples[[b-1]], 
                                    pars = list(a0 = a0,
                                                b0 = b0,
                                                beta = beta_samples[[b]]) )$out
  
  
}

burn = 25000
thin = 25
thin_samples <- seq(from = (burn + 1), to = num_iterations, by = thin)

plot(unlist(alpha_samples[thin_samples]), type = "l", ylab = "alpha", xlab = "Iteration post burn-in")
plot(unlist(gamma_samples[thin_samples]), type = "l", ylab = "gamma", xlab = "Iteration post burn-in")
#####################################################################################
## Clustering by VI
#####################################################################################
z_samples_VI <- list()

for(iter in 1:length(thin_samples)){
  z_samples_VI[[iter]] <- unlist(z_samples[[thin_samples[iter]]])
}

library(mcclust.ext)
posterior_samples <- matrix(0, nrow = length(z_samples_VI), ncol = length(z_samples_VI[[1]]))

for(i in 1:length(z_samples_VI)){
  posterior_samples[i, ] = z_samples_VI[[i]]
}

sim_mat = comp.psm(posterior_samples)
par(mfrow = c(1, 1))
plotpsm(sim_mat)


################################################################################
# GLOBAL LEVEL CLUSTERING
################################################################################
################################################################################
# CLUSTERING BY MINIMIZING VI
################################################################################
best_mcclust2 = minVI(sim_mat, posterior_samples, method = "all")

summary(best_mcclust2)

singleton_mcclust2 = as.numeric(names(which(table(best_mcclust2$cl[1, ]) <= 0.0 * length(best_mcclust2$cl[1, ]))))

singleton.index_mmclust2 = which(best_mcclust2$cl[1, ] %in% singleton_mcclust2)

z.estimated_mcclust2 = best_mcclust2$cl[1, ]
x1 = NULL
x2 = NULL
group = NULL
for(j in 1:J){
  x1 = c(x1, 1:length(x[[j]]))
  x2 = c(x2, x[[j]])
  group = c(group, rep(paste0("Population ", j), n[j]))
}

cluster.global_mcclust2 <- data.frame(x = x1,
                                      y = x2,
                                      cluster = factor(z.estimated_mcclust2),
                                      group = factor(group))


if(length(singleton.index_mmclust2) > 0){
  z.estimated_mcclust2 = z.estimated_mcclust2[-singleton.index_mmclust2]
  cluster.global_mcclust2 = cluster.global_mcclust2[-singleton.index_mmclust2, ]
}else{
  z.estimated_mcclust2 = z.estimated_mcclust2
  cluster.global_mcclust2 = cluster.global_mcclust2
}

library(tidyverse)

x.limit.lower <- min(x1)
x.limit.upper <- max(x1)

y.limit.lower <- min(x2)
y.limit.upper <- max(x2)

library(patchwork)
library(tidyverse)
library(latex2exp)
library(aricode)

myvalues_mcclust2 =  c("1" = "#F8766D",
                       "2" = "#00BA38",
                       "3" = "#619CFF", 
                       "4" = "blueviolet",
                       "5" = "cyan4",
                       "6" = "#E6AB02",
                       "7" = "#E36EF6",
                       "8" = "bisque4",
                       "9" = "coral4",
                       "10" = "darkslateblue",
                       
                       "11" = "lightseagreen",
                       "12" = "#E69F00", 
                       "13" = "#AA3377",
                       "14" = "sienna3",
                       "15" = "hotpink",
                       "16" = "sienna4",
                       "17" = "hotpink3",
                       "18" = "sienna1",
                       "19" = "dodgerblue4",
                       "20" = "bisque2",
                       
                       "21" = "darkgreen",
                       "22" = "orange", 
                       "23" = "maroon2",
                       "24" = "sienna2",
                       "25" = "hotpink4",
                       "26" = "sienna3",
                       "27" = "brown4",
                       "28" = "sienna1",
                       "29" = "dodgerblue3",
                       "30" = "bisque3")

ARI_population <- 0
for(j in 1:J){
  ARI_population[j] <- ARI(cluster.global_mcclust2 %>% filter(group == paste0("Population ", j)) %>% select(cluster) %>% pull(),
                           factor(true.Z[[j]]))
}

population = NULL
for(j in 1:J){
  population = c(population, rep(paste0("Population ", j, "\n ARI = ", round(ARI_population[j], 4)), n[j]))
}

cluster.global_mcclust2$population <- factor(population, 
                                             levels = factor(gtools::mixedsort(as.character(levels(factor(population))))))

cluster.global_mcclust2 %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 2) + labs(x = "Index", y = TeX("$x$")) + facet_wrap(~population, ncol = 5) + 
  scale_color_manual(values = myvalues_mcclust2) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() +  
  theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 


z1_samples <- list(); z2_samples <- list(); z3_samples <- list()
z_combined_samples <- list()

for(i in 1:length(thin_samples)){
  z1_samples[[i]] <- z_samples[[thin_samples[i]]][[1]]
  z2_samples[[i]] <- z_samples[[thin_samples[i]]][[2]]
  z3_samples[[i]] <- z_samples[[thin_samples[i]]][[3]]
  z_combined_samples[[i]] <- c(z1_samples[[i]],
                               z2_samples[[i]],
                               z3_samples[[i]])
}

getDahl(z1_samples, true.Z[[1]])
getDahl(z2_samples, true.Z[[2]])
getDahl(z3_samples, true.Z[[3]])
getDahl(z_combined_samples, unlist(true.Z))

n_samples[[num_iterations]]


ll <- 0
for(b in 2:num_iterations){
  # Printing the iterations
  if(b == 2){
    cat(paste0("Iteration: ", (b-1), "\n"))
  }
  if(b %% floor((10/100)*(num_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", b, "\n"))
  }
  ll[b] <- log_l2(data = matrix(x[[1]], ncol = 1), pi = pi_samples[[b]][[1]], z = z_samples[[b]][[1]], mu = matrix(mu_samples[[b]]), tau = array(tau, dim = c(1,1,L.max))) + log_l2(data = matrix(x[[2]], ncol = 1), pi = pi_samples[[b]][[2]], z = z_samples[[b]][[2]], mu = matrix(mu_samples[[b]]), tau = array(tau, dim = c(1,1,L.max))) + log_l2(data = matrix(x[[3]], ncol = 1), pi = pi_samples[[b]][[3]], z = z_samples[[b]][[3]], mu = matrix(mu_samples[[b]]), tau = array(tau, dim = c(1,1,L.max)))
}

ll_post_burn <- ll[thin_samples]
plot(ll_post_burn, type = "l", ylab = "", xlab = "Iteration post burn-in")
library(coda)
effectiveSize(ll_post_burn)
acf(ll_post_burn, lag = 50, main = "")

Pi1 = Pi1 + 1e-10; Pi2 = Pi2 + 1e-10; Pi3 = Pi3 + 1e-10; 
Pi1 = Pi1/sum(Pi1); Pi2 = Pi2/sum(Pi2); Pi3 = Pi3/sum(Pi3) 

ll.true <- log_l2(data = matrix(x[[1]], ncol = 1), pi = Pi1, z = true.Z[[1]], mu = matrix(phi.true), tau = array(tau, dim = c(1,1,L))) + log_l2(data = matrix(x[[2]], ncol = 1), pi = Pi2, z = true.Z[[2]], mu = matrix(phi.true), tau = array(tau, dim = c(1,1,L))) + log_l2(data = matrix(x[[3]], ncol = 1), pi = Pi3, z = true.Z[[3]], mu = matrix(phi.true), tau = array(tau, dim = c(1,1,L)))

plot(ll_post_burn, type = "l", ylab = "", xlab = "Iteration post burn-in")
abline(h = ll.true, col = "red")


# grid points for density estimation
xmin = min(unlist(x)) - 1 ; xmax = max(unlist(x)) + 1
y.grid = seq(xmin, xmax, length = 100)
# True density for each group
true.density = matrix(NA, nrow = J, ncol = n[1])
for(j in 1:J){
  
  # evaluate the density for each group
  true.density[j, ] = sapply(1:100, function(ii) 
    sum( Pi.true[j, ] * dnorm(y.grid[ii], mean = phi.true, sd = sqrt(1/tau)) ))
}


n.grid = length(y.grid)
density <- replicate(J, list())

for(iter in 1:length(thin_samples)){
  # J x n.grid matrix to store the densities along the rows.
  f = matrix(NA, nrow = J, ncol =  n.grid)
  phi = as.numeric(mu_samples[[thin_samples[iter]]])
  n.group = sapply(seq_len(L.max), function(j) sum(unlist(z_samples[[thin_samples[iter]]])==j))
  dat = data.frame(phi, n.group, ind = 1:L.max)
  dat = dat[order(-dat$n.group),]
  
  Pi <- matrix(0, J, L.max)
  for(j in 1:J) {
    Pi[j, ] <- pi_samples[[thin_samples[iter]]][[j]]
  }
  for(j in 1:J){
    
    # evaluate the density for each population
    f[j, ] = sapply(1:n.grid, function(ii) 
      sum( Pi[j, ] * dnorm(y.grid[ii], mean = phi, sd = sqrt(1/tau)) ))
    
    density[[j]][[iter]] <- f[j, ]
    
  }
  
}

# Estimated densities using posterior samples
est.density <- matrix(0, J, n[1])
for(j in 1:J){
  est.density[j, ] <- Reduce("+", density[[j]])/length(thin_samples)  
}


# data frames for plotting histograms and densities
dat_hist = data.frame(x = unlist(x), group = rep(c("Group 1", "Group 2", "Group 3"), each = n[1]))
dens = c(cbind(t(true.density), t(est.density)))
dat_dens = data.frame(grid = y.grid, density = dens, method = rep(c("True", "BGS with SALTSampler"), each = 3*n[1]), 
                      group = rep(c("Group 1", "Group 2", "Group 3"), each = n[1], times = 2))

# plot the histograms for each group
g.hist = ggplot(dat_hist, aes(x = x))  +
  geom_histogram(color="grey80", fill="grey90", aes(y=after_stat(density)), bins = 20)+
  facet_wrap(~as.factor(group), ncol = 1) + xlim(xmin, xmax) + theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

# overlay the histograms with density plots
g.hist + 
  geom_line(dat_dens, mapping = aes(x = grid, y = density, color = method, size = method)) +
  scale_size_manual(values = c("True" = 0.9, "BGS with SALTSampler" = 0.9)) +
  scale_color_manual(values=c("red2","blue4")) + labs(x = "", y = "") + 
  theme(legend.position = "top", legend.title=element_blank()) + theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 




density_CRF = replicate(J, list(replicate(length(CRF_HDP), list() )))
for(j in 1:J){
  for(i in 1:length(CRF_HDP)){
    density_CRF[[j]][[i]] = CRF_HDP[[i]]$density[j, ]
  }
}

# Estimated densities using posterior samples
est.density_CRF <- matrix(0, J, n[1])
for(j in 1:J){
  est.density_CRF[j, ] <- Reduce("+", density_CRF[[j]])/length(CRF_HDP)  
}


# data frames for plotting histograms and densities
dat_hist = data.frame(x = unlist(x), group = rep(c("Group 1", "Group 2", "Group 3"), each = n[1]))
dens = c(cbind(t(true.density), t(est.density), t(est.density_CRF)))
dat_dens = data.frame(grid = y.grid, density = dens, method = rep(c("True", "BGS with SALTSampler", "CRF"), each = 3*100), 
                      group = rep(c("Group 1", "Group 2", "Group 3"), each = 100, times = 3))

# plot the histograms for each group
g.hist = ggplot(dat_hist, aes(x = x))  +
  geom_histogram(color="grey80", fill="grey90", aes(y=after_stat(density)), bins = 20)+
  facet_wrap(~as.factor(group), ncol = 1) + xlim(xmin, xmax) + theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

# overlay the histograms with density plots
g.hist + 
  geom_line(dat_dens, mapping = aes(x = grid, y = density, color = method, size = method)) +
  scale_size_manual(values = c("True" = 0.9, "BGS with SALTSampler" = 0.9, "CRF" = 0.9)) +
  scale_color_manual(values=c("red2","blue4","darkgreen")) + labs(x = "", y = "") + 
  theme(legend.position = "top", legend.title=element_blank()) + theme(
    # LABLES APPEARANCE
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.title = element_text(hjust = 0.5, size=20, face= "bold", colour= "black" ),
    plot.subtitle = element_text(hjust = 0.5, size=16, face= "bold", colour= "black" ),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),    
    axis.title.y = element_text(size=16, face="bold", colour = "black"),    
    axis.text.x = element_text(size=16, face="bold", colour = "black"), 
    axis.text.y = element_text(size=16, face="bold", colour = "black"),
    strip.text.x = element_text(size = 14, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 14, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

# source("/Users/arhitchakrabarti/Downloads/blockedHDP/SS.R")
# 
# slice_HDP = hdp_slice_sampler(y = lapply(x, matrix), beta0=0.1, gam0=0.1, ITRmax=10000, Kcap=20, Tcap=20,
#                                           doubling_factor=1.5, categorical=F, W=10, cat_prior_alpha=NA,
#                                           randinit=F, prec2_y = tau, prec2_phi = lambda, y.grid = y.grid)
# 
# slice_HDP[[thin_samples[1]]]$Z
# slice_HDP[[10000]]
# source("/Users/arhitchakrabarti/Downloads/blockedHDP/postestimates.R")
# Z.hat.SS = getDahl(slice_HDP)
# 
# 
# # Adjusted Rand indices corresponding to group specific cluster labels
# ARI.SS = sapply(seq_len(J), function(j) mcclust::arandi(Z.hat.SS[[j]], true.Z[[j]]))
# ARI.SS
# 
# # Adjusted Rand index corresponding to global cluster labels
# ARI.global = mcclust::arandi(unlist(Z.hat), unlist(true.Z))
# ARI.global
