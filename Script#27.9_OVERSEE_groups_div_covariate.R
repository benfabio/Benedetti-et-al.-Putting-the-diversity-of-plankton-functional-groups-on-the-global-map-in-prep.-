
### ================================================================================================================

library("tidyverse")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("parallel")

WD <- getwd()
world2 <- map_data("world2")
world <- map_data("world")

### ================================================================================================================

strategies <- c("group")
categories <- c("Phytoplankton","Zooplankton")
combin <- apply(expand.grid(categories, strategies), 1, paste, collapse = "_")
SDMs <- c("GLM","GAM","ANN")
ESMs <- c("CNRM-PISCES","IPSL-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")		

res <- mclapply(combin, function(c) {
            
            # Extract strat and cat from c
            cat <- do.call(cbind,strsplit(c, "_"))[1,1]
            strat <- do.call(cbind,strsplit(c, "_"))[2,1]
            
            # Go to the communities wd
            if(cat == "Phytoplankton" & strat == "total") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/Scores")
                data.wd.back <- getwd()
            } else if(cat == "Zooplankton" & strat == "total") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/Scores")
                data.wd.back <- getwd()
            } else if(cat == "Phytoplankton" & strat == "group") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/Scores")
                data.wd.back <- getwd()
            } else if(cat == "Zooplankton" & strat == "group") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/Scores")
                data.wd.back <- getwd()
            } # eo else if loop 
            
            files <- dir()[grep("for_projections",dir())]
            # Rbind all
            scores <- lapply(files, function(f) { d <- get((load(f))) ; return(d) } ) # eo lapply - f in files
            table.scores <- do.call(rbind, scores) ; rm(scores)
            table.scores$SDM <- factor(t(data.frame(str_split(as.character(rownames(table.scores)), pattern = "_", n = 2)))[,1])
            table.scores$cat <- cat
            table.scores$strat <- strat
            
            return(table.scores)
    
    }, mc.cores = 4 
    
) # eo mclapply 
# Rbind
table.classif <- bind_rows(res)
rm(res) ; gc()

### ================================================================================================================

### Load projections for both phyto and zoo
setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
base <- get(load("table_mon_rich_baseline_groups_22_10_20.Rdata"))
fut <- get(load("table_mon_rich_2100-2000_groups_22_10_20.Rdata"))
base.ens <- data.frame(base[base$rich > 0,] %>% group_by(cell_id, group, SDM) %>%
                summarize(x = unique(x), y = unique(y),
                rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T))
) # eo ddf

fut.ens <- data.frame(fut[fut$rich > 0,] %>% group_by(cell_id, group, SDM, ESM) %>%
                summarize(x = unique(x), y = unique(y),
                rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T))
) # eo ddf
groups <- unique(base.ens$group) ; groups
# g <- "Calanoida"
res <- mclapply(groups, function(g) {
            
            base.sub <- base.ens[base.ens$group == g,]
            nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
            message(paste("Compiling ensemble changes in diversity for ",g," (",nsp," species)", sep = ""))
            # Flip x coordinates for Pacific-centered map
            base.sub$x2 <- base.sub$x
            base.sub[base.sub$x < 0 ,"x2"] <- (base.sub[base.sub$x < 0 ,"x"]) + 360
            # New id and order
            base.sub$id <- factor(paste(base.sub$x2, base.sub$y, sep = "_"))
            base.sub <- base.sub[order(base.sub$id),]
            # Compute mean annual baseline richness
            ens <- data.frame(base.sub %>% group_by(id) %>% summarize(x = unique(x2), y = unique(y),
                    rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T)) )
            ens$group <- g 
            # head(ens) ; summary(ens)
            return(ens)
    
        }, mc.cores = length(groups)

) # eo lapply 
# Rbind
ddf.rich <- bind_rows(res)    
head(ddf.rich) ; dim(ddf.rich)
rm(res) ; gc()

### Dcast to put the groups' rich or per as columns
table.rich <- dcast(data = ddf.rich[,c(1:4,6)], formula = id + x + y ~ group, value.var = "rich")
head(table.rich) ; dim(table.rich)
summary(table.rich)
rm(ddf.rich) ; gc()

### Great, now combine with mean annual climatologies of environmental covariates
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
files <- dir()[grep("18_09_20",dir())] ; files
clims <- mclapply(files, function(f) {
            message(paste(f, sep = ""))
            c <- read.table(f, h = T, sep = ";")
            c <- c[-which(c$SSS < 20),] 
            c <- c[-which(c$Bathy > -175),] 
            c$x2 <- c$x
            c[c$x < 0 ,"x2"] <- (c[c$x < 0 ,"x"]) + 360
            c$id <- factor(paste(c$x2, c$y, sep = "_"))
            c <- c[order(c$id),] 
            return(c)
        }, mc.cores = 12
) # eo mclapply - f in files
# Rbind
clim.base <- bind_rows(clims)
rm(clims) ; gc()
head(clim.base) ; dim(clim.base) ; summary(clim.base$x2)

### Compute annual climatologies prior to combining with diversity estimates
ann.clim <- data.frame(clim.base %>% group_by(id) %>%
    summarize(x = unique(x2), y = unique(y), SST = mean(SST,na.rm=T), dSST = mean(dSST,na.rm=T),
    SSS = mean(SSS,na.rm=T), logEKE = mean(logEKE,na.rm=T), Wind = mean(Wind,na.rm=T), MLD = mean(MLD,na.rm=T), PAR = mean(PAR,na.rm=T), 
    pCO2 = mean(pCO2,na.rm=T), dO2 = mean(dO2,na.rm=T), logNO3 = mean(logNO3,na.rm=T), logSiO2 = mean(logSiO2,na.rm=T),
    Nstar = mean(Nstar,na.rm=T), Sistar = mean(Sistar,na.rm=T), logChl = mean(logChl,na.rm=T) )
) # eo ddf
summary(ann.clim)
rm(clim.base) ; gc()
# Quick map
# ggplot() + geom_raster(aes(x = x, y = y, fill = MLD), data = ann.clim) +
#     scale_fill_viridis(name = "MLD") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#        panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### Combine with table.rich using common cells
commons <- intersect(unique(table.rich$id), unique(ann.clim$id))
# length(commons)
ann.clim2 <- ann.clim[ann.clim$id %in% commons,]
table.rich2 <- table.rich[table.rich$id %in% commons,]
dim(ann.clim2) ; dim(table.rich2)
# cbind
table.rich3 <- cbind(table.rich2, ann.clim2[,c(4:length(ann.clim2))])
head(table.rich3)

### Save new data.frame because it take some time to compute
setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
save(table.rich3, file = "table_ann_rich_env_baseline_23_10_20.Rdata")
rm(ann.clim2,table.rich2,ddf.grp.div,fut,base) ; gc()

### ================================================================================================================

setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
table.rich <- get(load("table_ann_rich_env_baseline_23_10_20.Rdata"))

### Compute pairwise Spearman correlation coeff between groups' richness and env
library("corrplot")
library("corrgram")
library("ggcorrplot")
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

names <- colnames(table.rich)[c(4:length(table.rich))] ; names
mydata <- na.omit(table.rich[,names])
cormat <- round(cor(mydata, method = "spearman"),2)
p.mat <- cor_pmat(mydata, method = "spearman", conf.level = 0.95)
# Re-order?
cor_tri <- get_upper_tri(cormat)
cor_tri <- melt(cor_tri, na.rm = T)
pval_tri <- get_upper_tri(p.mat)
pval_tri <- melt(pval_tri, na.rm = T)
colnames(cor_tri)[3] <- "rho"
colnames(pval_tri)[3] <- "pval"
cor_tri$pval <- pval_tri$pval
rm(pval_tri,p.mat,cormat) ; gc()
# unique(cor_tri$Var1)
# unique(cor_tri$Var2)

### First heatmap: main zoo groups versus hydro
plankton <- names[c(1:14)]
env <- names[c(15:28)]
cormat <- cor_tri[which(cor_tri$Var1 %in% plankton & cor_tri$Var2 %in% env),]
cormat <- cormat[!(cormat$Var1 == cormat$Var2),]
head(cormat)

# Same for hydrobio data
levels(cormat$Var2)[levels(cormat$Var2) == "Sistar"] <- "Si*"
levels(cormat$Var2)[levels(cormat$Var2) == "Nstar"] <- "N*"
levels(cormat$Var2)[levels(cormat$Var2) == "dO2"] <- "O2 (200m)"
levels(cormat$Var2)[levels(cormat$Var2) == "logEKE"] <- "EKE (log)"
levels(cormat$Var2)[levels(cormat$Var2) == "dSST"] <- "SST range"
levels(cormat$Var2)[levels(cormat$Var2) == "logNO3"] <- "NO3 (log)"
levels(cormat$Var2)[levels(cormat$Var2) == "logSiO2"] <- "SiO2 (log)"
levels(cormat$Var2)[levels(cormat$Var2) == "logChl"] <- "Chl-a (log)"
 
# Add signif label
cormat$signif <- NA
cormat[cormat$pval <= 0.001,"signif"] <- "***"
cormat[cormat$pval > 0.001 & cormat$pval <= 0.01,"signif"] <- "**"
cormat[cormat$pval > 0.01 & cormat$pval < 0.05,"signif"] <- "*"
cormat[cormat$pval >= 0.05,"signif"] <- "-"

colnames(cormat)[c(1,2)] <- c("Group","Covariate")

plot <- ggplot(cormat, aes(factor(Covariate), factor(Group), fill = rho)) + geom_tile(color = "white") +
        scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), 
                name = paste("Spearman' rank\ncorrelation\ncoefficient", sep = ""), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), 
                guide = guide_colourbar(ticks = T, ticks.colour = "black", frame.colour = "black")) + 
        xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
                axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
        coord_fixed() + geom_text(aes(factor(Covariate), factor(Group), label = signif), color = "black", size = 4)
#
ggsave(plot = plot, filename = paste("heatmap_corr_ann_rich_env_hydro.jpg", sep = ""), dpi = 300, width = 8, height = 8)

### Issue: the global latitudinal SST-driven SR pattern hides some potentially interesting second-order patterns, derive 
### anomalies from linear models to analyze the part of the div gradients that do not mainly covary with SST

### OR: perform PCA on env and check correlation to PC1-3
require("FactoMineR")
data4pca <- na.omit(table.rich)
pca.env <- PCA(X = data4pca[,env], scale.unit = T, graph = F, ncp = 5)
summary(pca.env)
#plot(pca.env, choix = "var", axes = c(1,2))
#plot(pca.env, choix = "var", axes = c(3,4))
# Provide the PC coords to data4pca
data4pca[,paste("PC",c(1:4),sep="")] <- pca.env$ind$coord[,c(1:4)]
eig <- data.frame(perc = pca.env$eig[,"percentage of variance"], nb = c(1:nrow(pca.env$eig)) ) # eig
pca1 <- paste0("PC1 (",floor(eig$perc[1]*100)/100,"%)")
pca2 <- paste0("PC2 (",floor(eig$perc[2]*100)/100,"%)")
pca3 <- paste0("PC3 (",floor(eig$perc[3]*100)/100,"%)")
pca4 <- paste0("PC4 (",floor(eig$perc[4]*100)/100,"%)")
summary(data4pca)

# Map PCs quickly
map <- ggplot() + geom_raster(aes(x = x, y = y, fill = PC4), data = data4pca) +
    geom_contour(colour = "grey50", binwidth = 1, size = 0.25, aes(x = x, y = y, z = PC4), data = data4pca) +
    scale_fill_gradient2(name = pca4, mid = "white", low = "#2166ac", high = "#b2182b") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
	    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") + 
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

setwd(WD)
ggsave(plot = map, filename = "map_ann_env_PC4.jpg", dpi = 300, width = 7, height = 4)

### Interpretation 
# PC1: latitudinal SST and nutrients availability gradient (negative = tropics, positive = high latitudes)
# PC2 > 0 --> MLD/O2/Wind/Nstar/SSS ; < 0 --> S*, SiO2, Chl, pCO2
# PC3 > 0 --> PCO2, EKE ; < 0 --> SST range
# PC4 > 0 --> N*/Si*/SSS ; < 0 --> dSST/logChl/ EKE

### Very nice, and now re-compute spearman's correlations$
names <- colnames(data4pca)[c(4:17,32:35)] ; names
mydata <- data4pca[,names]
cormat <- round(cor(mydata, method = "spearman"),2)
p.mat <- cor_pmat(mydata, method = "spearman", conf.level = 0.95)
# Re-order?
cor_tri <- get_upper_tri(cormat)
cor_tri <- melt(cor_tri, na.rm = T)
pval_tri <- get_upper_tri(p.mat)
pval_tri <- melt(pval_tri, na.rm = T)
colnames(cor_tri)[3] <- "rho"
colnames(pval_tri)[3] <- "pval"
cor_tri$pval <- pval_tri$pval
rm(pval_tri,p.mat,cormat) ; gc()

### First heatmap: main zoo groups versus hydro
plankton <- names[c(1:14)]
env <- names[c(15:18)]
cormat <- cor_tri[which(cor_tri$Var1 %in% plankton & cor_tri$Var2 %in% env),]
cormat <- cormat[!(cormat$Var1 == cormat$Var2),]
head(cormat)

# Same for hydrobio data
levels(cormat$Var2)[levels(cormat$Var2) == "PC1"] <- pca1
levels(cormat$Var2)[levels(cormat$Var2) == "PC2"] <- pca2
levels(cormat$Var2)[levels(cormat$Var2) == "PC3"] <- pca3
levels(cormat$Var2)[levels(cormat$Var2) == "PC4"] <- pca4

# Add signif label
cormat$signif <- NA
cormat[cormat$pval <= 0.001,"signif"] <- "***"
cormat[cormat$pval > 0.001 & cormat$pval <= 0.01,"signif"] <- "**"
cormat[cormat$pval > 0.01 & cormat$pval < 0.05,"signif"] <- "*"
cormat[cormat$pval >= 0.05,"signif"] <- "-"

colnames(cormat)[c(1,2)] <- c("Group","Covariate")

plot <- ggplot(cormat, aes(factor(Covariate), factor(Group), fill = rho)) + geom_tile(color = "white") +
        scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), 
                name = paste("Spearman' rank\ncorrelation\ncoefficient", sep = ""), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), 
                guide = guide_colourbar(ticks = T, ticks.colour = "black", frame.colour = "black")) + 
        xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
                axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
        coord_fixed() + geom_text(aes(factor(Covariate), factor(Group), label = signif), color = "black", size = 4)
#
ggsave(plot = plot, filename = paste("heatmap_corr_ann_rich_env_PCs.jpg", sep = ""), dpi = 300, width = 8, height = 8)

### Interpretation: 
### PC1: The richness of all groups decerases with PC1 --> increase with SST/PAR/SSS but decreases with NO3/SiO2/Chl/Wind/MLD
### PC2: Amphipods/Append/Calanoida/Diatoms/Dinos/Forams/Pteropods/Salps increase with S*, SiO2, Chl, pCO2 but decrease with MLD/O2/Wind/Nstar
### PC2 Coccos/Euphausiids/Jellyfish/Oithona show the opposite pattern
### PC3: Poecilo/Oitho/Chaeto decrease with pCO2/EKE and increase with SST range
### PC3: Forams and Coccos decrease with SST range and increase with pCO2/EKE
### PC4: Calanoida/Chaeto/Cocco/Euphausiids/Jellyfish/Oitho/Poecilo/Ptero increase with N*/Si*/SSS but decerase with dSST/logChl/ EKE



