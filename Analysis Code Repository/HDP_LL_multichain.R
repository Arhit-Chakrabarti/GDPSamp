rm(list = ls())

LL1.mat <- matrix(NA, nrow = length(LL1[[1]]), ncol = length(LL1))
LL2.mat <- matrix(NA, nrow = length(LL2[[1]]), ncol = length(LL2))
LL3.mat <- matrix(NA, nrow = length(LL3[[1]]), ncol = length(LL3))
LL4.mat <- matrix(NA, nrow = length(LL4[[1]]), ncol = length(LL4))


for(i in 1:length(LL2)){
  LL1.mat[, i] <- LL1[[i]]
  LL2.mat[, i] <- LL2[[i]]
  LL3.mat[, i] <- LL3[[i]]
  LL4.mat[, i] <- LL4[[i]]
}


LL1.dataframe <- cbind(data.frame(LL1.mat), Chain = factor("Chain 1"), Iteration = rep(1:length(LL1[[1]]), 2))

LL2.dataframe <- cbind(data.frame(LL2.mat), Chain = factor("Chain 2"), Iteration = rep(1:length(LL2[[1]]), 2))

LL3.dataframe <- cbind(data.frame(LL3.mat), Chain = factor("Chain 3"), Iteration = rep(1:length(LL3[[1]]), 2))

LL4.dataframe <- cbind(data.frame(LL4.mat), Chain = factor("Chain 4"), Iteration = rep(1:length(LL4[[1]]), 2))

LL <- rbind(LL1.dataframe, LL2.dataframe, LL3.dataframe, LL4.dataframe)


# POPULATION 1
library(tidyverse)
ll1 = LL %>% ggplot(aes(x = Iteration, y = X1, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 1") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))


# POPULATION 2
library(tidyverse)
ll2 = LL %>% ggplot(aes(x = Iteration, y = X2, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 2") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# POPULATION 3
library(tidyverse)
ll3 = LL %>% ggplot(aes(x = Iteration, y = X3, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 3") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# POPULATION 4
library(tidyverse)
ll4 = LL %>% ggplot(aes(x = Iteration, y = X4, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 4") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# POPULATION 5
library(tidyverse)
ll5 = LL %>% ggplot(aes(x = Iteration, y = X5, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 5") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# POPULATION 6
library(tidyverse)
ll6 = LL %>% ggplot(aes(x = Iteration, y = X6, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 6") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# POPULATION 7
library(tidyverse)
ll7 = LL %>% ggplot(aes(x = Iteration, y = X7, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 7") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))


# POPULATION 8
library(tidyverse)
ll8 = LL %>% ggplot(aes(x = Iteration, y = X8, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 8") +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)
if(!require("grid")) install.packages("grid"); library(grid)
library(cowplot)
p.len = get_plot_component(LL %>% ggplot(aes(x = Iteration, y = X8, color = Chain)) + geom_line() + theme_bw() + labs(y = "", title = "Log-likelihood for group 8") +
                             theme(legend.title = element_text(face = "bold", size = 14),
                                   legend.text = element_text(size = 14),
                                   axis.text = element_text(face = "bold", size = 12),
                                   axis.title = element_text(face = "bold", size = 12),
                                   plot.title = element_text(face = "bold", size = 14)),
                           'guide-box-right', return_all = TRUE)
gridExtra::grid.arrange(ll1, ll2, ll3, ll4, ll5, ll6, ll7, ll8, p.len)
