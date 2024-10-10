out10 = readRDS("/Users/arhitchakrabarti/Documents/ServerGDPLS/GDP_L_10_final_LS.rds")
out20 = readRDS("/Users/arhitchakrabarti/Documents/ServerGDPLS/GDP_L_20_final_LS.rds")
out30 = readRDS("/Users/arhitchakrabarti/Documents/ServerGDPLS/GDP_L_30_final_LS.rds")
out50 = readRDS("/Users/arhitchakrabarti/Documents/ServerGDPLS/GDP_L_50_final_LS.rds")

cluster10 = 0
for(i in 1:length(out10)){
  cluster10[i] = out10[[i]][[1]]
}
sum(!is.na(cluster10))
# cluster10 = sort(cluster10, decreasing = T)[1:50]

cluster20 = 0
for(i in 1:length(out20)){
  cluster20[i] = out20[[i]][[1]]
}
sum(!is.na(cluster20))
# cluster20 = sort(cluster20, decreasing = T)[1:50]

cluster30 = 0
for(i in 1:length(out30)){
  cluster30[i] = out30[[i]][[1]]
}
sum(!is.na(cluster30))
# cluster30 = sort(cluster30, decreasing = T)[1:50]

cluster50 = 0
for(i in 1:length(out50)){
  cluster50[i] = out50[[i]][[1]]
}
sum(!is.na(cluster50))
# cluster50 = sort(cluster50, decreasing = T)[1:50]

library(tidyverse)
data.frame(N = c(cluster10, cluster20, cluster30, cluster50), 
           Truncation = c(rep("10", length(cluster10)),
                          rep("20", length(cluster20)),
                          rep("30", length(cluster30)),
                          rep("50", length(cluster50)))) %>% ggplot(aes(x = Truncation, y = N)) + geom_boxplot(width = 0.3, position = position_dodge(width = 0.9)) + labs(x = "Truncation level (L)", y = "Number of clusters", title = "Boxplot of number of clusters with truncation level of GDP", subtitle = "50 replicates") + theme_bw() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))

summary(cluster10)
summary(cluster20)
summary(cluster30)
summary(cluster50)