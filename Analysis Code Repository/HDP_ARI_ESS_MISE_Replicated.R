out = out_org[1:50]
length(out)

ARI_LS_Matrix = matrix(0, nrow = length(out), ncol = 3)
ARI_CRF_Matrix = matrix(0, nrow = length(out), ncol = 3)

for(b in 1:length(out)){
  ARI_LS_Matrix[b, ] <-  out[[b]]$ARI_LS
  ARI_CRF_Matrix[b, ] <-  out[[b]]$ARI_CRF
}

ARI_Data = NULL
for(j in 1:3){
  for(b in 1:length(out)){
    ARI_Data = rbind(ARI_Data, data.frame(BGS = ARI_LS_Matrix[b, j],
                                          CRF = ARI_CRF_Matrix[b, j],
                                          Group = paste0("Group", j)))
  }
}

ARI_Data_Long1 = unname(cbind(ARI_Data[, c(1,3)], "BGS with SALTSampler"))

names(ARI_Data_Long1) <- c("ARI", "Group", "Method")

ARI_Data_Long2 = unname(cbind(ARI_Data[, c(2,3)], "CRF"))
names(ARI_Data_Long2) <- c("ARI", "Group", "Method")

ARI_Data_Long = rbind(ARI_Data_Long1, ARI_Data_Long2)

library(tidyverse)
ARI_plot = ARI_Data_Long %>% group_by(Method) %>% ggplot(aes(x = Group, y = ARI, col = Method)) + geom_boxplot(width = 0.5, position = position_dodge(width = 0.8)) + labs(x = "Group", y = "ARI") + theme_bw() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))







ESS_LS_Matrix = matrix(0, nrow = length(out), ncol = 3)
ESS_CRF_Matrix = matrix(0, nrow = length(out), ncol = 3)

for(b in 1:length(out)){
  ESS_LS_Matrix[b, ] <-  out[[b]]$ESS_density
  ESS_CRF_Matrix[b, ] <-  out[[b]]$ESS_density_CRF
}

ESS_Data = NULL
for(j in 1:3){
  for(b in 1:length(out)){
    ESS_Data = rbind(ESS_Data, data.frame(BGS = ESS_LS_Matrix[b, j],
                                          CRF = ESS_CRF_Matrix[b, j],
                                          Group = paste0("Group", j)))
  }
}

ESS_Data_Long1 = unname(cbind(ESS_Data[, c(1,3)], "BGS with SALTSampler"))

names(ESS_Data_Long1) <- c("ESS", "Group", "Method")

ESS_Data_Long2 = unname(cbind(ESS_Data[, c(2,3)], "CRF"))
names(ESS_Data_Long2) <- c("ESS", "Group", "Method")

ESS_Data_Long = rbind(ESS_Data_Long1, ESS_Data_Long2)

library(tidyverse)
ESS_plot = ESS_Data_Long %>% group_by(Method) %>% ggplot(aes(x = Group, y = ESS, col = Method)) + geom_boxplot(width = 0.5, position = position_dodge(width = 0.8)) + labs(x = "Group", y = "ESS") + theme_bw() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))




MISE_LS_Matrix = matrix(0, nrow = length(out), ncol = 3)
MISE_CRF_Matrix = matrix(0, nrow = length(out), ncol = 3)

for(b in 1:length(out)){
  MISE_LS_Matrix[b, ] <-  out[[b]]$MISE
  MISE_CRF_Matrix[b, ] <-  out[[b]]$MISE_CRF
}

MISE_Data = NULL
for(j in 1:3){
  for(b in 1:length(out)){
    MISE_Data = rbind(MISE_Data, data.frame(BGS = MISE_LS_Matrix[b, j],
                                          CRF = MISE_CRF_Matrix[b, j],
                                          Group = paste0("Group", j)))
  }
}

MISE_Data_Long1 = unname(cbind(MISE_Data[, c(1,3)], "BGS with SALTSampler"))

names(MISE_Data_Long1) <- c("MISE", "Group", "Method")

MISE_Data_Long2 = unname(cbind(MISE_Data[, c(2,3)], "CRF"))
names(MISE_Data_Long2) <- c("MISE", "Group", "Method")

MISE_Data_Long = rbind(MISE_Data_Long1, MISE_Data_Long2)

library(tidyverse)
MISE_plot = MISE_Data_Long %>% group_by(Method) %>% ggplot(aes(x = Group, y = MISE, col = Method)) + geom_boxplot(width = 0.5, position = position_dodge(width = 0.8)) + labs(x = "Group", y = "MISE") + theme_bw() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14))


legend = cowplot::get_plot_component(MISE_Data_Long %>% group_by(Method) %>% ggplot(aes(x = Group, y = MISE, col = Method)) + geom_boxplot(width = 0.5, position = position_dodge(width = 0.8)) + labs(x = "Group", y = "MISE") + theme_bw() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.title.x=element_blank(),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14)),
  'guide-box-right', return_all = TRUE)

gridExtra::grid.arrange(ARI_plot, ESS_plot, MISE_plot, legend, nrow = 1)
