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

################################################################################
# Run a Vanilla DP on each group
################################################################################
source("GDP_chain_functions_all_functions2.R")
################################################################################
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

# Prior specifications
p = nrow(X[[1]])
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
z.start = rcat(n = n[1], prob = beta.start)
mu.start <- t(sapply(1:L.max, function(j){mu.start.list[[j]]}))
alpha.start <- 2
a0 = 0.1; b0 = 0.1
# a0 = 1; b0 = 1
# Prior hyper-parameters
prior = list(prior_mean = prior_mean, prior_sigma = prior_sigma, a0 = a0, b0 = b0)
# Initial values of different parameters
init = list(mu.start = mu.start, beta.start = beta.start, z.start = z.start, alpha.start = alpha.start)
num_iterations <- 25000
burn = 15000
thin = 10

beta.start.small <- list()
beta.start.long <- list()
beta.start_GDP <- list()
z.start_GDP <- list()
plots <- list()

for(j in 1:J){
  cat(paste0("\nRunning a Vanilla DP for population ", j, "\n"))
  gibbs_out <- gibbs_MVN_dp_c(init = init, data = X[[j]], num_iterations = num_iterations, burn = burn, thin = thin, L.max = L.max, prior = prior)
  library(latex2exp)
  library(tidyverse)
  # DAHL's METHOD
  dp_nw.dahl = getDahl(gibbs_out$z_samples, z.true[[j]])
  cluster <- data.frame(x = X[[j]][1, ],
                        y = X[[j]][2, ], 
                        z = factor(gibbs_out$z_samples[[dp_nw.dahl$DahlIndex]]))
  plots[[j]] <- cluster %>% ggplot(aes(x = x, y = y, col = z)) + geom_point(size = 2) + labs(title = paste("Population ", j), subtitle = paste0("Adjusted Rand index = ", round(dp_nw.dahl$rand, 4))) + labs(x = TeX("$X_1$"), y = TeX("$X_2$")) + 
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
  
  
  DahlIndex = dp_nw.dahl$DahlIndex
  z.dahl = gibbs_out$z_samples[[DahlIndex]]
  # Get the estimate of beta1 by running a Naive DP on group 1 data
  beta.start.small[[j]] = gibbs_out$beta_samp_post[[DahlIndex]][unique(z.dahl)]
  
  beta.start.small[[j]] = beta.start.small[[j]]/sum(beta.start.small[[j]])
  
  # But for the next step again I would require beta starting value to be of length L.max
  # So I add a very small quantity to the remaining co-ordinates and re-normlaize
  beta.start.long[[j]] = c(beta.start.small[[j]], rep(1e-8, abs(L.max - length(beta.start.small[[j]]))))
  beta.start_GDP[[j]] = beta.start.long[[j]] /sum(beta.start.long[[j]] ) # This becomes the starting value of beta1 in the Gibbs iteration
  # Draw starting values of latent indicators for group 1
  z.start_GDP[[j]] = rcat(n = n[j], prob = beta.start_GDP[[j]])
}

library(gridExtra)
do.call(grid.arrange, plots)
#####################################################################################################
# RUNNING GDP
#####################################################################################################
#####################################################################################################
## Initialize storage
#####################################################################################################
num_iterations = 50000
# Latent indicators
z_samples <- replicate(num_iterations, list(replicate(J, list())))
# Weights
beta_samples <- replicate(num_iterations, list(replicate(J, list())))
# Concentration parameters
alpha_samples <- replicate(num_iterations, list(replicate(J, list())))
# To save cluster specific sample sizes
n_samples <- replicate(num_iterations, list(replicate(J, list())))
# Storing the means and precisions of the normal distribution
mu_samples <- replicate(num_iterations, list())
L.obs = L.max

mu.start_GDP.list = list()
for(i in 1:L.obs){
  mu.start_GDP.list[[i]] = matrix(mvtnorm::rmvnorm(n = 1, mean = as.numeric(prior_mean), sigma = prior_sigma), ncol = 1)
}

mu.start_GDP <- t(sapply(1:L.max, function(j){mu.start_GDP.list[[j]]}))

for(j in 1:J){
  beta_samples[[1]][[j]] <- beta.start_GDP[[j]]
  z_samples[[1]][[j]] <- z.start_GDP[[j]]
  alpha_samples[[1]][[j]] <- 2
}


mu_samples[[1]] <- mu.start_GDP
################################################################################
# Running the Sampler
################################################################################
if(!require("SALTSampler")) install.packages("SALTSampler"); library(SALTSampler)
if(!require("Rcpp")) install.packages("Rcpp"); library(Rcpp)
if(!require("fossil")) install.packages("fossil"); library(fossil)

source("RunMH_my.R")

time.start = Sys.time()
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
    sample = sample_z_known_sigma_c(data = t(X[[j]]), beta = beta_samples[[b-1]][[j]], mu = mu_samples[[b-1]], Sigma = Sigma.true, L = L.obs)
    
    z_samples[[b]][[j]] = sample # Store the samples values of zj
  }
  
  
  out = sample_mu_gdp_c(data = X, z = z_samples[[b]], Sigma = Sigma.true, L = L.obs, n = n, prior_mean = prior_mean,  prior_sigma = prior_sigma)
  # Store the samples of mu, and n_j
  mu_samples[[b]] =  t(out$mu_draw);
  
  for(j in 1:J){
    n_samples[[b]][[j]] = out$n_h[, 1][[j]]
  }
  
  
  init = beta.start_GDP[[1]]
  
  # pars should include beta2, beta3, beta4, nu1, nu2, nu3, eta, alpha2, alpha3, alpha4
  # Run SALTSampler for Beta1
  beta1 <- RunMH_my(Target = Target_beta1, init = init, B = 1,
                    h = rep(2, L.obs), 
                    pars = list(beta2 = beta_samples[[b-1]][[2]],
                                alpha2 = alpha_samples[[b-1]][[2]]),
                    dat = n_samples[[b]][[1]], 
                    concentration = rep(alpha_samples[[b-1]][[1]], length(init)))
  
  beta_samples[[b]][[1]] = as.numeric(beta1$S)
  
  for(j in 2:(J-1)){
    # j = 2
    init = beta.start_GDP[[j]]
    # Run SALTSampler for Beta1
    betaj <- RunMH_my(Target = Target_betaj, init = init, B = 1,
                      h = rep(2, L.obs), 
                      pars = list(beta_j_minus_1 = beta_samples[[b]][[(j - 1)]],
                                  beta_j_plus_1 = beta_samples[[b-1]][[(j + 1)]],
                                  alpha_j_plus_1 = alpha_samples[[b-1]][[(j + 1)]]),
                      dat = n_samples[[b]][[j]], 
                      concentration = rep(alpha_samples[[b-1]][[j]], length(init)))
    
    beta_samples[[b]][[j]] = as.numeric(beta1$S)
    beta_samples[[b]][[j]] = beta_samples[[b]][[j]] + 1e-8
    beta_samples[[b]][[j]] = beta_samples[[b]][[j]]/sum(beta_samples[[b]][[j]])
  }
  # Draw samples from the beta2
  new_concentration = alpha_samples[[b-1]][[J]] * beta_samples[[b]][[(J-1)]]
  # Now we only need to consider components 1, 2,...,L.obs 
  beta_samples[[b]][[J]] = as.numeric(rdirichlet(n = 1, alpha = n_samples[[b]][[J]] + new_concentration))
  
  beta_samples[[b]][[J]] = beta_samples[[b]][[J]] + 1e-8
  beta_samples[[b]][[J]] = beta_samples[[b]][[J]]/sum(beta_samples[[b]][[J]])
  
  alpha_samples[[b]][[1]] = sample_alpha1(ycurrent = alpha_samples[[b-1]][[1]], 
                                          pars = list(a0 = prior$a0,
                                                      b0 = prior$b0,
                                                      beta1 = beta_samples[[b]][[1]]) )
  
  
  for(j in 2:(J-1)){
    # Run MH Sampler for alpha_j
    alpha_samples[[b]][[j]] <- sample_alphaj_j_minus_1(ycurrent = alpha_samples[[b-1]][[j]], 
                                                       pars = list(a0 = prior$a0,
                                                                   b0 = prior$b0,
                                                                   beta_j = beta_samples[[b]][[j]],
                                                                   beta_j_minus_1 = beta_samples[[b]][[(j-1)]],
                                                                   alpha_j_minus_1 = alpha_samples[[b]][[j-1]],
                                                                   alpha_j_plus_1 = alpha_samples[[b-1]][[j+1]]) )
  }  
  
  alpha_samples[[b]][[J]] <- sample_alphaJ(ycurrent = alpha_samples[[b-1]][[J]], 
                                           pars = list(a0 = prior$a0,
                                                       b0 = prior$b0,
                                                       beta_j = beta_samples[[b]][[J]],
                                                       beta_j_minus_1 = beta_samples[[b]][[(J-1)]],
                                                       alpha_j_minus_1 = alpha_samples[[b]][[J-1]]) )
  
}
time.end = Sys.time()
time.end - time.start

log_likelihood = 0

for(b in 2:num_iterations){
  # Printing the iterations
  if(b == 2){
    cat(paste0("Iteration: ", (b-1), "\n"))
  }
  if(b %% floor((10/100)*(num_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", b, "\n"))
  }
  ll = 0
  for(j in 1:J){
    ll = ll + log_ll_chain(data = t(X[[j]]), pi = beta_samples[[b]][[j]], z = z_samples[[b]][[j]], mu = matrix(mu_samples[[b]], ncol = 2), tau = array(Sigma.true, dim = c(2, 2, L.max)))
  }
  log_likelihood[b] <- ll
}

burn = 35000
thin = 15
thin_samples <- seq(from = (burn + 1), to = num_iterations, by = thin)

log_like <- log_likelihood[thin_samples]

ll_plot <- data.frame(x = 1:length(log_like), y = log_like) %>% ggplot(aes(x = x, y = y)) + geom_line() +
  labs(
    # title = "Traceplot of log-likelihood",
    
    x = "Iteration post burn-in", y = "") +
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

library(forecast)
ACF_plot <-  ggAcf(x = log_like, lag.max = 40) +
  ggtitle("") +
  labs(y = "") +
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

ll_plot
ACF_plot
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


# save.image("~/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/SavedSimulationData/50NodesChainGDP_alpha0_15.RData")

alpha_samples_matrix <- matrix(unlist(alpha_samples),ncol=J,byrow=TRUE)[thin_samples, ]
library(latex2exp)
plots_alpha <- list()
for(j in 1:J){
  cluster <- data.frame(x = 1:length(alpha_samples_matrix[, j]),
                        y = alpha_samples_matrix[, j])
  new_text = paste0("Traceplot of $\\alpha_{",j,"}$")
  plots_alpha [[j]] <- cluster %>% ggplot(aes(x = x, y = y)) + geom_line() + theme_bw() + labs(x = "Iteration post burn-in", y = "", title = TeX(new_text)) + theme_classic() +  
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
}

library(gridExtra)
do.call(grid.arrange, plots_alpha)
