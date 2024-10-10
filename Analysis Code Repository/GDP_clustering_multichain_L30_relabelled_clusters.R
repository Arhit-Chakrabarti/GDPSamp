if(!require("extraDistr")) install.packages("extraDistr"); library(extraDistr)
if(!require("MASS")) install.packages("MASS"); library(MASS)
if(!require("mvtnorm")) install.packages("mvtnorm"); library(mvtnorm)
if(!require("fossil")) install.packages("fossil"); library(fossil)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

# Change the observations according to the groups from the information provided by Ellen Ruth
X1 <- t(X3.UMAP)
X2 <- t(X1.UMAP)
X3 <- t(X4.UMAP)
X4 <- t(X7.UMAP)
X5 <- t(X5.UMAP)
X6 <- t(X2.UMAP)
X7 <- t(X8.UMAP)
X8 <- t(X6.UMAP)

n1 = ncol(X1)
n2 = ncol(X2)
n3 = ncol(X3)
n4 = ncol(X4)
n5 = ncol(X5)
n6 = ncol(X6)
n7 = ncol(X7)
n8 = ncol(X8)

z1.chain <- readRDS("/Users/arhitchakrabarti/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/RealDataMultiChain/GDP/GDP_z_chain_l1_alpha_3_L_30.rds")
z2.chain <- readRDS("/Users/arhitchakrabarti/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/RealDataMultiChain/GDP/GDP_z_chain_l2_alpha_3_L_30.rds")
z3.chain <- readRDS("/Users/arhitchakrabarti/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/RealDataMultiChain/GDP/GDP_z_chain_l3_alpha_3_L_30.rds")
z4.chain <- readRDS("/Users/arhitchakrabarti/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/RealDataMultiChain/GDP/GDP_z_chain_l4_alpha_3_L_30.rds")


z1_samples.post.burn <- c(z1.chain[[1]], z2.chain[[1]], z3.chain[[1]], z4.chain[[1]])
z2_samples.post.burn <- c(z1.chain[[2]], z2.chain[[2]], z3.chain[[2]], z4.chain[[2]])
z3_samples.post.burn <- c(z1.chain[[3]], z2.chain[[3]], z3.chain[[3]], z4.chain[[3]])
z4_samples.post.burn <- c(z1.chain[[4]], z2.chain[[4]], z3.chain[[4]], z4.chain[[4]])
z5_samples.post.burn <- c(z1.chain[[5]], z2.chain[[5]], z3.chain[[5]], z4.chain[[5]])
z6_samples.post.burn <- c(z1.chain[[6]], z2.chain[[6]], z3.chain[[6]], z4.chain[[6]])
z7_samples.post.burn <- c(z1.chain[[7]], z2.chain[[7]], z3.chain[[7]], z4.chain[[7]])
z8_samples.post.burn <- c(z1.chain[[8]], z2.chain[[8]], z3.chain[[8]], z4.chain[[8]])

z_samples <- list()

for(iter in 1:length(z1_samples.post.burn)){
  z_samples[[iter]] <- c(z1_samples.post.burn[[iter]],
                         z2_samples.post.burn[[iter]],
                         z3_samples.post.burn[[iter]],
                         z4_samples.post.burn[[iter]],
                         z5_samples.post.burn[[iter]],
                         z6_samples.post.burn[[iter]],
                         z7_samples.post.burn[[iter]],
                         z8_samples.post.burn[[iter]])
}


# Calculate the membership matrix for Group 1
M <- lapply(z_samples, function(x){
  clusterAssign <- x
  Matrix::Matrix(1 * outer(clusterAssign, clusterAssign, FUN = "=="), sparse = TRUE)
})
# Mean membership matrix
M.mean <- Reduce("+", M)/length(z_samples)
# Calculate the Frobenius norm of the differences
M.Frobenius <- sapply(M, function(x, av) sum((x - av)^2),
                      av = M.mean)
# Find out the minimums
k.min <- which.min(M.Frobenius)
# Remove the membership matrix to free up memory
rm(M)
rm(M.mean)

################################################################################
## POST PROCESSING
################################################################################
singleton1 = as.numeric(names(which(table(z1_samples.post.burn[[k.min]]) <= 0.05 * n1)))
singleton.index1 = which(z1_samples.post.burn[[k.min]] %in% singleton1)
if(length(singleton.index1) > 0){
  z1.estimated = z1_samples.post.burn[[k.min]][-singleton.index1]
}else{
  z1.estimated = z1_samples.post.burn[[k.min]]
}

singleton2 = as.numeric(names(which(table(z2_samples.post.burn[[k.min]]) <= 0.05 * n2)))
singleton.index2 = which(z2_samples.post.burn[[k.min]] %in% singleton2)
if(length(singleton.index2) > 0){
  z2.estimated = z2_samples.post.burn[[k.min]][-singleton.index2]
}else{
  z2.estimated = z2_samples.post.burn[[k.min]]
}


singleton3 = as.numeric(names(which(table(z3_samples.post.burn[[k.min]]) <= 0.05 * n3)))
singleton.index3 = which(z3_samples.post.burn[[k.min]] %in% singleton3)
if(length(singleton.index3) > 0){
  z3.estimated = z3_samples.post.burn[[k.min]][-singleton.index3]
}else{
  z3.estimated = z3_samples.post.burn[[k.min]]
}


singleton4 = as.numeric(names(which(table(z4_samples.post.burn[[k.min]]) <= 0.05 * n4)))
singleton.index4 = which(z4_samples.post.burn[[k.min]] %in% singleton4)
if(length(singleton.index4) > 0){
  z4.estimated = z4_samples.post.burn[[k.min]][-singleton.index4]
}else{
  z4.estimated = z4_samples.post.burn[[k.min]]
}

singleton5 = as.numeric(names(which(table(z5_samples.post.burn[[k.min]]) <= 0.05 * n5)))
singleton.index5 = which(z5_samples.post.burn[[k.min]] %in% singleton5)
if(length(singleton.index5) > 0){
  z5.estimated = z5_samples.post.burn[[k.min]][-singleton.index5]
}else{
  z5.estimated = z5_samples.post.burn[[k.min]]
}

singleton6 = as.numeric(names(which(table(z6_samples.post.burn[[k.min]]) <= 0.05 * n6)))
singleton.index6 = which(z6_samples.post.burn[[k.min]] %in% singleton6)
if(length(singleton.index6) > 0){
  z6.estimated = z6_samples.post.burn[[k.min]][-singleton.index6]
}else{
  z6.estimated = z6_samples.post.burn[[k.min]]
}

singleton7 = as.numeric(names(which(table(z7_samples.post.burn[[k.min]]) <= 0.05 * n7)))
singleton.index7 = which(z7_samples.post.burn[[k.min]] %in% singleton7)
if(length(singleton.index7) > 0){
  z7.estimated = z7_samples.post.burn[[k.min]][-singleton.index7]
}else{
  z7.estimated = z7_samples.post.burn[[k.min]]
}

singleton8 = as.numeric(names(which(table(z8_samples.post.burn[[k.min]]) <= 0.05 * n8)))
singleton.index8 = which(z8_samples.post.burn[[k.min]] %in% singleton8)
if(length(singleton.index8) > 0){
  z8.estimated = z8_samples.post.burn[[k.min]][-singleton.index8]
}else{
  z8.estimated = z8_samples.post.burn[[k.min]]
}

#################################################################################################################
# PLOT OF CLUSTERING
#################################################################################################################
library(tidyverse)

x.limit.lower <- min(X1[1, ], X2[1, ], X3[1, ], X4[1, ], X5[1, ], X6[1, ], X7[1, ], X8[1, ])
x.limit.upper <- max(X1[1, ], X2[1, ], X3[1, ], X4[1, ], X5[1, ], X6[1, ], X7[1, ], X8[1, ])

y.limit.lower <- min(X1[2, ], X2[2, ], X3[2, ], X4[2, ], X5[2, ], X6[2, ], X7[2, ], X8[2, ])
y.limit.upper <- max(X1[2, ], X2[2, ], X3[2, ], X4[2, ], X5[2, ], X6[2, ], X7[2, ], X8[2, ])

myvalues = c("1" = "#F8766D",
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

z1.estimated = factor(z1.estimated)
z1.estimated <- recode(z1.estimated, "1" = '1', 
                  "2" = '2',
                  "3" = '3', 
                  "4" = '4',
                  "13" = '5',
                  "16" = '6',
                  "20" = '7',
                  "26" = '8')
# Population 1
if(length(singleton.index1) > 0){
  cluster1 <- data.frame(x1 = X1[1, -singleton.index1],
                         x2 = X1[2, -singleton.index1], 
                         cluster = z1.estimated)
}else{
  cluster1 <- data.frame(x1 = X1[1, ],
                         x2 = X1[2, ], 
                         cluster = factor(z1.estimated))
}

library(latex2exp)

p1 = cluster1 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 1", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))   +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Population 2
z2.estimated = factor(z2.estimated)
z2.estimated <- recode(z2.estimated, "1" = '1', 
                       "2" = '2',
                       "3" = '3', 
                       "4" = '4',
                       "13" = '5',
                       "16" = '6',
                       "20" = '7',
                       "26" = '8')

if(length(singleton.index2) > 0){
  cluster2 <- data.frame(x1 = X2[1, -singleton.index2],
                         x2 = X2[2, -singleton.index2], 
                         cluster = factor(z2.estimated))
}else{
  cluster2 <- data.frame(x1 = X2[1, ],
                         x2 = X2[2, ], 
                         cluster = factor(z2.estimated))
}


p2 = cluster2 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 2", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))  +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper)  + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Population 3
z3.estimated = factor(z3.estimated)
z3.estimated <- recode(z3.estimated, "1" = '1', 
                       "2" = '2',
                       "3" = '3', 
                       "4" = '4',
                       "13" = '5',
                       "16" = '6',
                       "20" = '7',
                       "26" = '8')

if(length(singleton.index3) > 0){
  cluster3 <- data.frame(x1 = X3[1, -singleton.index3],
                         x2 = X3[2, -singleton.index3], 
                         cluster = factor(z3.estimated))
}else{
  cluster3 <- data.frame(x1 = X3[1, ],
                         x2 = X3[2, ], 
                         cluster = factor(z3.estimated))
}

p3 = cluster3 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 3", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))  +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Population 4
z4.estimated = factor(z4.estimated)
z4.estimated <- recode(z4.estimated, "1" = '1', 
                       "2" = '2',
                       "3" = '3', 
                       "4" = '4',
                       "13" = '5',
                       "16" = '6',
                       "20" = '7',
                       "26" = '8')

if(length(singleton.index4) > 0){
  cluster4 <- data.frame(x1 = X4[1, -singleton.index4],
                         x2 = X4[2, -singleton.index4], 
                         cluster = factor(z4.estimated))
}else{
  cluster4 <- data.frame(x1 = X4[1, ],
                         x2 = X4[2, ], 
                         cluster = factor(z4.estimated))
}

p4 = cluster4 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 4", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))   +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Population 5
z5.estimated = factor(z5.estimated)
z5.estimated <- recode(z5.estimated, "1" = '1', 
                       "2" = '2',
                       "3" = '3', 
                       "4" = '4',
                       "13" = '5',
                       "16" = '6',
                       "20" = '7',
                       "26" = '8')

if(length(singleton.index5) > 0){
  cluster5 <- data.frame(x1 = X5[1, -singleton.index5],
                         x2 = X5[2, -singleton.index5], 
                         cluster = factor(z5.estimated))
}else{
  cluster5 <- data.frame(x1 = X5[1, ],
                         x2 = X5[2, ], 
                         cluster = factor(z5.estimated))
}

p5 = cluster5 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 5", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))  +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Population 6
z6.estimated = factor(z6.estimated)
z6.estimated <- recode(z6.estimated, "1" = '1', 
                       "2" = '2',
                       "3" = '3', 
                       "4" = '4',
                       "13" = '5',
                       "16" = '6',
                       "20" = '7',
                       "26" = '8')

if(length(singleton.index6) > 0){
  cluster6 <- data.frame(x1 = X6[1, -singleton.index6],
                         x2 = X6[2, -singleton.index6], 
                         cluster = factor(z6.estimated))
}else{
  cluster6 <- data.frame(x1 = X6[1, ],
                         x2 = X6[2, ], 
                         cluster = factor(z6.estimated))
}

p6 = cluster6 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 6", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))   +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Population 7
z7.estimated = factor(z7.estimated)
z7.estimated <- recode(z7.estimated, "1" = '1', 
                       "2" = '2',
                       "3" = '3', 
                       "4" = '4',
                       "13" = '5',
                       "16" = '6',
                       "20" = '7',
                       "26" = '8')

if(length(singleton.index7) > 0){
  cluster7 <- data.frame(x1 = X7[1, -singleton.index7],
                         x2 = X7[2, -singleton.index7], 
                         cluster = factor(z7.estimated))
}else{
  cluster7 <- data.frame(x1 = X7[1, ],
                         x2 = X7[2, ], 
                         cluster = factor(z7.estimated))
}

p7 = cluster7 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 7", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))  +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

# Population 8
z8.estimated = factor(z8.estimated)
z8.estimated <- recode(z8.estimated, "1" = '1', 
                       "2" = '2',
                       "3" = '3', 
                       "4" = '4',
                       "13" = '5',
                       "16" = '6',
                       "20" = '7',
                       "26" = '8')

if(length(singleton.index8) > 0){
  cluster8 <- data.frame(x1 = X8[1, -singleton.index8],
                         x2 = X8[2, -singleton.index8], 
                         cluster = factor(z8.estimated))
}else{
  cluster8 <- data.frame(x1 = X8[1, ],
                         x2 = X8[2, ], 
                         cluster = factor(z8.estimated))
}

p8 = cluster8 %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(size = 2, alpha = 1) + labs(title = "Clustering for group 8", x = TeX(r'($UMAP_1$)'), y = TeX(r'($UMAP_2$)'))  +
  scale_color_manual(values = myvalues) +
  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) + theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))


cluster_merged = bind_rows(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7, cluster8)
cluster_merged$cluster <- factor(cluster_merged$cluster)
sorted_labels <- paste(sort(as.integer(levels(cluster_merged$cluster))))

plot_merged = cluster_merged %>% ggplot(aes(x = x1, y = x2, col = cluster)) + geom_point(alpha = 1, size = 3) + labs(title = "Clustering for group 8", x = "UMAP1", y = "UMAP2")  +  xlim(x.limit.lower, x.limit.upper) + ylim(y.limit.lower, y.limit.upper) +
  scale_color_manual(values = myvalues, breaks = sorted_labels)  +
  theme_classic() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

library(gridExtra)
library(grid)

library(cowplot)
p.len = get_plot_component(plot_merged, 'guide-box-bottom', return_all = TRUE)

gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p.len)
