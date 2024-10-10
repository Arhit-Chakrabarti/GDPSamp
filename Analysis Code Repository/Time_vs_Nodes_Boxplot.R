out10 = readRDS("T_10_3.rds")
length(out10)
time10 = 0

for(i in 1:length(out10)){
  if(length(out10[[i]]) > 2){
    time10[i] = as.numeric(out10[[i]]$time.taken)
  }
}

time10 = sort(time10[!is.na(time10)], decreasing = T)[1:50]

out20 = readRDS("T_20_3.rds")
length(out20)
time20 = 0

for(i in 1:length(out20)){
  if(length(out20[[i]]) > 2){
    time20[i] = as.numeric(out20[[i]]$time.taken)
  }
}

time20 = sort(time20[!is.na(time20)], decreasing = T)[1:50]
time20 = time20*60

out50 = readRDS("T_50_3.rds")
length(out50)
time50 = 0

for(i in 1:length(out50)){
  if(length(out50[[i]]) > 2){
    time50[i] = as.numeric(out50[[i]]$time.taken)
  }
}

time50 = sort(time50[!is.na(time50)], decreasing = F)[1:50]
time50 = time50*60

output = data.frame(Time = c(time10,
                             time20,
                             time50),
                    Nodes = c(rep("10", length(time10)),
                              rep("20", length(time20)),
                               rep("50", length(time50))))

library(tidyverse)
time_taken <- output %>% ggplot(aes(x = Nodes, y = Time)) + geom_boxplot(width = 0.3, position = position_dodge(width = 0.9)) + labs(x = "Number of nodes (T)", y = "Run-time (minutes)") + theme_bw() + 
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 14),
        axis.text = element_text(face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14)) 

time_taken

mean(time10)
mean(time20)
mean(time50)
