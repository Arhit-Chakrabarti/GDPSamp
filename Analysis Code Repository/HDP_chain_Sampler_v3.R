rm(list = ls())
################################################################################
# Generate the data
################################################################################
if(!require("extraDistr")) install.packages("extraDistr"); library(extraDistr)
if(!require("MASS")) install.packages("MASS"); library(MASS)
if(!require("mvtnorm")) install.packages("mvtnorm"); library(mvtnorm)
if(!require("fossil")) install.packages("fossil"); library(fossil)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

L.true = 5 # true number of groups in population 1
alpha0 = 15
J = 10
# True alpha for generating data
alpha.true <- list()

alpha.true[[1]] = rgamma(n = 1, shape = alpha0, rate = 1) 
for(j in 2:J){
  alpha.true[[j]] = rgamma(n = 1, shape = alpha.true[[j-1]], rate = 1)
}

# for(j in 2:J){
#   alpha.true[[j]] = rgamma(n = 1, shape = alpha0, rate = 1)
# }
# True weights
beta.true <- list()
beta.true[[1]] <- as.numeric(rdirichlet(n = 1, alpha = rep(alpha.true[[1]]/L.true, L.true))) 
beta.true[[1]] <- (beta.true[[1]]  + 1e-6)/sum((beta.true[[1]]  + 1e-6))

for(j in 2:J){
  beta.true[[j]] <- as.numeric(rdirichlet(n = 1, alpha = (1/L.true) + alpha.true[[j]] * beta.true[[j-1]])) 
  beta.true[[j]] <- (beta.true[[j]]  + 1e-6)/sum((beta.true[[j]]  + 1e-6))
}
# Sample sizes
n <- rep(100, J)

z.true <- list()
for(j in 1:J){
  z.true[[j]] <- sample(1:L.true, size = n[j], prob = beta.true[[j]], replace = TRUE)
}


# True cluster specific means
mean.true = cbind(c(-2, -5), c(0, 0), c(-3, 3), c(3,-3), c(8, 5)) # True means of the clusters

Sigma.true = matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2) # True covariance matrix for all populations

# Draw data from Normal populations
X <- list()
for(j in 1:J){
  X[[j]] = sapply(1:n[j], function(i){
    mvrnorm(n = 1, mu = mean.true[, z.true[[j]]][, i], Sigma = Sigma.true)
  })
}

x1 = NULL
x2 = NULL
group = NULL
for(j in 1:J){
  x1 = c(x1, X[[j]][1, ])
  x2 = c(x2, X[[j]][2, ])
  group = c(group, rep(paste0("Population ", j), n[j]))
}

x.limit.lower <- min(x1)
x.limit.upper <- max(x1)

y.limit.lower <- min(x2)
y.limit.upper <- max(x2)

library(tidyverse)
library(latex2exp)

data.frame(x = x1, y = x2, population = factor(group, levels = gtools::mixedsort(unique(group)))) %>%
  ggplot(aes(x = x, y = y)) + geom_point(size = 2) + labs(x = TeX("$X_1$"), y = TeX("$X_2$")) + facet_wrap(~population, ncol = 5)  +
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

library(Rcpp)
################################################################################
source("GDP_chain_functions_all_functions2.R")
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

# Prior specifications
p = nrow(X[[1]])
prior_mean = as.matrix(rep(0, p), ncol = 1) # Prior mean of Normal
lambda.inv = 0.01
lambda = 1/lambda.inv
prior_sigma = diag(lambda, p)

L.max = 10

mu.start_HDP.list = list()

for(i in 1:L.max){
  mu.start_HDP.list[[i]] = matrix(mvtnorm::rmvnorm(n = 1, mean = as.numeric(prior_mean), sigma = prior_sigma), ncol = 1)
}

mu.start_HDP <- t(sapply(1:L.max, function(j){mu.start_HDP.list[[j]]}))

beta.start_HDP = rep(1/L.max, L.max)

alpha.start <- 2
a0 = 0.1; b0 = 0.1


for(j in 1:J){
  pi_samples[[1]][[j]] <- rep(1/L.max, L.max)
  z_samples[[1]][[j]] <- sample(1:L.max, size = n[j], prob = pi_samples[[1]][[j]], replace = TRUE)
}

alpha_samples[[1]] <- 2
gamma_samples[[1]] <- 2

mu_samples[[1]] <- mu.start_HDP
sigma.true = Sigma.true


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
    sample = sample_z_known_sigma_c(data = t(X[[j]]), beta = pi_samples[[b-1]][[j]], mu = mu_samples[[b-1]], Sigma = Sigma.true, L = L.max)
    
    z_samples[[b]][[j]] = sample # Store the samples values of zj
  }
  
  
  out = sample_mu_gdp_c(data = X, z = z_samples[[b]], Sigma = Sigma.true, L = L.max, n = n, prior_mean = prior_mean,  prior_sigma = prior_sigma)
  # Store the samples of mu, and n_j
  mu_samples[[b]] =  t(out$mu_draw);
  
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

burn = 35000
thin = 15
thin_samples <- seq(from = (burn + 1), to = num_iterations, by = thin)

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
  x1 = c(x1, X[[j]][1, ])
  x2 = c(x2, X[[j]][2, ])
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

ARI_population <- 0
for(j in 1:J){
  ARI_population[j] <- aricode::ARI(cluster.global_mcclust2 %>% filter(group == paste0("Population ", j)) %>% dplyr::select(cluster) %>% pull(),
                                    factor(z.true[[j]]))
}

population = NULL
for(j in 1:J){
  population = c(population, rep(paste0("Population ", j, "\n ARI = ", round(ARI_population[j], 4)), n[j]))
}

cluster.global_mcclust2$population <- factor(population, 
                                             levels = factor(gtools::mixedsort(as.character(levels(factor(population))))))


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

cluster.global_mcclust2 %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 2) + labs(x = TeX("$X_1$"), y = TeX("$X_2$")) + facet_wrap(~population, ncol = 5) + 
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
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 12, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

################################################################################
# CLUSTERING BY LEAST SQUARES
################################################################################
# Calculate the membership matrix for Group 1
M <- lapply(z_samples_VI, function(x){
  clusterAssign <- x
  Matrix::Matrix(1 * outer(clusterAssign, clusterAssign, FUN = "=="), sparse = TRUE)
})
# Mean membership matrix
M.mean <- Reduce("+", M)/length(z_samples_VI)
# Calculate the Frobenius norm of the differences
M.Frobenius <- sapply(M, function(x, av) sum((x - av)^2),
                      av = M.mean)
# Find out the minimums
k.min <- which.min(M.Frobenius)
# Remove the membership matrix to free up memory
rm(M)
rm(M.mean)


z_samples_post_burn <- replicate(length(thin_samples), list(replicate(J, list())))

for(j in 1:J){
  for(iter in 1:length(thin_samples)){
  z_samples_post_burn[[iter]][[j]] <- z_samples[[thin_samples[iter]]][[j]]
  }
}

z.estimated_LS <- list()
for(j in 1:J){
  z.estimated_LS[[j]] <- z_samples_post_burn[[k.min]][[j]]
}


cluster.global_LS <- data.frame(x = x1,
                                        y = x2,
                                        cluster = factor(unlist(z.estimated_LS)),
                                        group = factor(group))
library(tidyverse)
library(latex2exp)
ARI_population_LS <- 0
for(j in 1:J){
  ARI_population_LS[j] <- aricode::ARI(cluster.global_LS %>% filter(group == paste0("Population ", j)) %>% dplyr::select(cluster) %>% pull(),
                                    factor(z.true[[j]]))
}

population = NULL
for(j in 1:J){
  population = c(population, rep(paste0("Population ", j, "\n ARI = ", round(ARI_population_LS[j], 4)), n[j]))
}

cluster.global_LS$population <- factor(population, 
                                             levels = factor(gtools::mixedsort(as.character(levels(factor(population))))))

cluster.global_LS %>%
  ggplot(aes(x = x, y = y, col = cluster)) + geom_point(size = 2) + labs(x = TeX("$X_1$"), y = TeX("$X_2$")) + facet_wrap(~population, ncol = 5) + 
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
    axis.text.x = element_text(size=12, face="bold", colour = "black"), 
    axis.text.y = element_text(size=12, face="bold", colour = "black"),
    strip.text.x = element_text(size = 12, face="bold", colour = "black" ),
    strip.text.y = element_text(size = 12, face="bold", colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill=NA, size=0.3),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14)
  ) 

# save.image("~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/SavedSimulationData/20NodesChainHDP_alpha0_15.RData")