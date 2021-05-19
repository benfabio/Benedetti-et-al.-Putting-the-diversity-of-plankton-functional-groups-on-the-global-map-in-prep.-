
##### 27/08/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Assessing pairwise predictors correlations, from species-level to group level
#   - For total background data only, load the species data (1/0) and compute pairwise Spearman's rank corr coeff
 
### Last update: 27/08/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")
library("parallel")
library("corrgram")

world2 <- map_data("world2")
world <- map_data("world")
WD <- getwd() 

### Draw the correlation coeff heatmaps for each net
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
# Utiliser la corrélation entre les variables
  # comme mésure de distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <- cormat[hc$order, hc$order]
}


# --------------------------------------------------------------------------------------------------------------------------------

### A°) First, Zooplankton
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
files <- dir()[grep("data_",dir())]# ; files
f <- files[5]

res <- mclapply(files, function(f) {
    
            # Message
            data <- read.table(f, sep = ";", h = T)
            message(paste(unique(data[data$obs == 1,'species']), sep = ""))
            names <- colnames(data)[c(19,23:33,35,37:41)]
            mydata <- na.omit(data[,names]) #  head(mydata)
            cormat <- round(cor(mydata, method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = TRUE)
            melted_cormat <- melted_cormat[!(melted_cormat$Var1 == melted_cormat$Var2),]
            
            melted_cormat$species <- unique(data[data$obs == 1,'species'])
            melted_cormat$group <- unique(data[data$obs == 1,'group2'])
            return(melted_cormat)
                
    }, mc.cores = 30
    
) # eo mclapply
# Rbind
table <- bind_rows(res) ; rm(res) ; gc()
dim(table) ; head(table)
summary(table)

### Compute mean or median correlation
table$id <- paste(table$Var1, table$Var2, sep = "_")

ddf <- data.frame(table %>% group_by(id) %>% summarize(V1 = unique(Var1), V2 = unique(Var2), mean = round(mean(value),2), median = round(median(value),2)) )
head(ddf) ; summary(ddf)

# ggplot(ddf, aes(factor(V2), factor(V1), fill = mean)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Mean\nRho") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(V2), factor(V1), label = mean), color = "black", size = 3)
#
ggplot(ddf, aes(factor(V2), factor(V1), fill = median)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Median\nRho") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(V2), factor(V1), label = median), color = "black", size = 3)


### B°) Phytoplankton
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
files <- dir()[grep("data_",dir())] #; files
f <- files[100]

res <- mclapply(files, function(f) {
    
            # Message
            data <- read.table(f, sep = ";", h = T)
            message(paste(unique(data[data$obs == 1,'species']), sep = ""))
            names <- colnames(data)[c(18,22:32,34,36:40)]
            mydata <- na.omit(data[,names]) #  head(mydata)
            cormat <- round(cor(mydata, method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = TRUE)
            melted_cormat <- melted_cormat[!(melted_cormat$Var1 == melted_cormat$Var2),]
    
            melted_cormat$species <- unique(data[data$obs == 1,'species'])
            melted_cormat$group <- unique(data[data$obs == 1,'group2'])
            return(melted_cormat)
    
    }, mc.cores = 30
    
) # eo mclapply
table <- bind_rows(res) ; rm(res) ; gc()
dim(table) ; head(table)
table$id <- paste(table$Var1, table$Var2, sep = "")

### Compute mean or median correlation
table$id <- paste(table$Var1, table$Var2, sep = "_")

ddf <- data.frame(table %>% group_by(id) %>% summarize(V1 = unique(Var1), V2 = unique(Var2), mean = round(mean(value),2), median = round(median(value),2)) )
head(ddf) ; summary(ddf)

ggplot(ddf, aes(factor(V2), factor(V1), fill = mean)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Mean\nRho") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(V2), factor(V1), label = mean), color = "black", size = 3)
#
ggplot(ddf, aes(factor(V2), factor(V1), fill = median)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Median\nRho") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(V2), factor(V1), label = median), color = "black", size = 3)
    
### Use these to define pools of non collinear covariates

### Phytoplankton: discard O2, logPO4, SSS, Pstar
### Choose between NO3 and SiO2
### MLPAR vs MLD+PAR
# Set1: SST, logChla, pCO2, logNO3, Nstar, Sistar, MLPAR, logEKE, SLA, Wind
# Set2: SST, logChla, pCO2, logSiO2, Nstar, Sistar, MLD, PAR, logEKE, SLA, Wind

### Zooplankton: discard SSS, logPO4, Pstar
### Choose one nutrients field (SiO2 because less correlated to SST)
# Set1: SST, logChla, pCO2, logSiO2, Nstar, Sistar, MLPAR, logEKE, SLA, Wind
# Set2: SST, logChla, pCO2, logSiO2, Nstar, Sistar, MLD, PAR, logEKE, SLA, Wind



