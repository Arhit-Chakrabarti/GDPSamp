GDP_cluster <- readRDS("/Users/arhitchakrabarti/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/RealDataMultiChain/GDP/BestGDPCombinedCluster.rds")
HDP_cluster <- readRDS("/Users/arhitchakrabarti/Library/CloudStorage/OneDrive-TexasA&MUniversity/GDP/JMLR_Revision/RealDataMultiChain/HDP/BestHDPCombinedCluster.rds")
################################################################################################
# COMAPRE CLUSTERING PERFORMANCE
################################################################################################
# FOR EACH SAMPLE 1 to 4 CHANGE THE NUMBER IN THE list() accordingly
dataa.aug <- list()
# Comapre the Within SS
calc_SS <- function(df) sum(as.matrix(stats::dist(df)^2)) / (2 * nrow(df))

library(tidyverse)

WithinSS <- NULL

for(i in 1:length(GDP_cluster)){
  X <- cbind(GDP_cluster[[i]]$x1, GDP_cluster[[i]]$x2)
  dataa.aug[[i]] <- data.frame(X, 
                               GDP = as.character(GDP_cluster[[i]]$cluster), 
                               HDP = as.character(HDP_cluster[[i]]$cluster))
  
  
  # Function to compare clustering performance by looking at Total Within Sum of Squares
  WithinSS <- bind_rows(WithinSS, dataa.aug[[i]] %>%
                          gather(method, cluster, GDP, HDP) %>%
                          group_by(method, cluster) %>%
                          nest() %>%
                          transmute(
                            method,
                            cluster,
                            within_SS = map_dbl(data, ~calc_SS(.x))) %>%
                          group_by(method) %>%
                          summarise(total_within_SS = sum(within_SS)) %>%
                          spread(method, total_within_SS)
  )
}

colSums(WithinSS)

dataa.augC <- NULL
for(i in 1:length(GDP_cluster)){
  X <- cbind(GDP_cluster[[i]]$x1, GDP_cluster[[i]]$x2)
  dataa.augC <- rbind(dataa.augC, data.frame(X, 
                                             GDP = as.character(GDP_cluster[[i]]$cluster), 
                                             HDP = as.character(HDP_cluster[[i]]$cluster)))
}

dataa.augC %>%
  gather(method, cluster, GDP, HDP) %>%
  group_by(method, cluster) %>%
  nest() %>%
  transmute(
    method,
    cluster,
    within_SS = map_dbl(data, ~calc_SS(.x))) %>%
  group_by(method) %>%
  summarise(total_within_SS = sum(within_SS)) %>%
  spread(method, total_within_SS)

library(clusterSim)
GDP.cluster = factor(dataa.augC[, "GDP"])
levels(GDP.cluster)[match("7", levels(GDP.cluster))] <- "6"

HDP.cluster = factor(dataa.augC[, "HDP"])
# Calinski-Harabasz Index
CH.GDP = index.G1(dataa.augC[, c(1, 2)],cl = as.integer(GDP.cluster), centrotypes = "centroids")
CH.GDP
CH.HDP = index.G1(dataa.augC[, c(1, 2)],cl = as.integer(HDP.cluster), centrotypes = "centroids")
CH.HDP

# Davies-Bouldin Index
DB.GDP = index.DB(dataa.augC[, c(1, 2)],cl = as.integer(GDP.cluster), centrotypes = "centroids", p=2, q=2)

DB.HDP = index.DB(dataa.augC[, c(1, 2)],cl = as.integer(HDP.cluster), centrotypes = "centroids", p=2, q=2)

DB.GDP$DB
DB.HDP$DB

# Silhoutte Index
dataa <- dataa.augC[, c(1, 2)]
colnames(dataa) <- NULL
distance.data <- stats::dist(dataa)

sil.GDP <- index.S(distance.data, as.integer(GDP.cluster))
sil.HDP <- index.S(distance.data, as.integer(HDP.cluster))

sil.GDP
sil.HDP

