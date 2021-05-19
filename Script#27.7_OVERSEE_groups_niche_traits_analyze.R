
##### 08/10/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Loading the phytoplankton and zooplankton niche traits from Script#27.6
#   - Check their overall distribution and then per group
#   - Test variance in niche traits across groups, and phyto vs. zooplankton for common variables
#   - Ordinate phytoplankton and zooplankton specie sin niche traits space?
#   - Try to compute and plot a group-level mean response curves
 
### Last update: 21/10/2020 

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("biomod2")
library("viridis")
library("vegan")
library("FactoMineR")

WD <- getwd()
world2 <- map_data("world2")
world <- map_data("world")

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Phytoplankton 
# Go to working directories
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/resp_curves/niche_traits/")
files <- dir()[grep("niche_traits_GAM",dir())] # files
# Concatenate all niche traits
require("parallel")
res <- mclapply(files, function(f) {
            t <- get(load(f))
            return(t)
    }, mc.cores = 20
) # eo mclapply - f in files
# Rbind
traits_phy <- bind_rows(res)
rm(res,files) ; gc()
head(traits_phy) ; dim(traits_phy)

### Plot distrbution of niche traits (center and center2, width and width2)
ggplot(aes(x = factor(group), y = center, fill = factor(group)), data = traits_phy[traits_phy$var == "SST",]) +
     geom_boxplot(colour = "black") + scale_fill_brewer(name = "") +
     xlab("") + ylab("SST niche center (°C)") + theme_classic()

ggplot(aes(x = factor(group), y = width, fill = factor(group)), data = traits_phy[traits_phy$var == "SST",]) +
     geom_boxplot(colour = "black") + scale_fill_brewer(name = "") +
     xlab("") + ylab("SST niche width (°C)") + theme_classic()

ggplot(aes(x = center, y = width, fill = factor(group)), data = traits_phy[traits_phy$var == "SST",]) +
     geom_point(colour = "black", pch = 21) + scale_fill_brewer(name = "", palette = "Spectral") +
     ylab("SST niche width (°C)") + xlab("SST niche center (°C)") + theme_classic()

### Print those plots per variable
vars <- unique(traits_phy$var) ; vars
v <- "SST"
for(v in vars) {
    
    message(paste("Plotting niche traits ditsribution across groups for ",v, sep = ""))
    
    p1 <- ggplot(aes(x = factor(group), y = center2, fill = factor(group)), data = traits_phy[traits_phy$var == v,]) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(paste(v," niche center",sep="")) + scale_x_discrete(labels = NULL) +
        theme_classic()
    #
    p2 <- ggplot(aes(x = factor(group), y = width2, fill = factor(group)), data = traits_phy[traits_phy$var == v,]) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(paste(v," niche width",sep="")) + scale_x_discrete(labels = NULL) + 
        theme_classic()
    
    p3 <- ggplot(aes(x = center2, y = width2, fill = factor(group)), data = traits_phy[traits_phy$var == v,]) + 
        geom_point(colour = "black", pch = 21) + scale_fill_brewer(name = "", palette = "Paired") + 
        ylab(paste(v," niche width",sep="")) + xlab(paste(v," niche center",sep="")) + theme_classic()
    
    p4 <- ggplot(aes(x = factor(group), y = prob.range, fill = factor(group)), data = traits_phy[traits_phy$var == v,]) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(paste("HSI range (importance)",sep="")) + scale_x_discrete(labels = NULL) + 
        scale_y_continuous(limits = c(0,1)) + theme_classic()
    
    require("ggpubr") 
    panel <- ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, align = "hv")
    # panel
    setwd(WD)
    ggsave(plot = panel, filename = paste("panel_niche_traits_phyto_",v,".jpg", sep = ""), width = 8, height = 5, dpi = 300)
    
} # eo v in vars

### Summarize variations in those niche traits (let's keep center2 and width2) in a PCA
### But first...need to reformat the data.frame so columns look like: SST_center; SST_width; etc. etc.
### Need to use melt and dcast I guess
#head(traits_phy)
m_traits_phy <- melt(traits_phy, id.vars = c("group","species","var"))
unique(m_traits_phy$variable)
# Create every combination of var x variable in a factor column and dcast
m_traits_phy$new.col <- factor(paste(m_traits_phy$var ,m_traits_phy$variable, sep = "_"))
unique(m_traits_phy$new.col)
# Dcast to put those new.col as columns (return to wide format)
d_traits_phy <- dcast(m_traits_phy[,c(1,2,5,6)], group + species ~ new.col, value.var = "value") 
# Select the columns corresponding to: center2, width2 and weight.pca (weights to give to the columns, not the individuals)
cols2keep <- c(colnames(d_traits_phy)[grep("center2", colnames(d_traits_phy))],
                colnames(d_traits_phy)[grep("width2", colnames(d_traits_phy))],
                colnames(d_traits_phy)[grep("weight.pca", colnames(d_traits_phy))]
) # eo cols2keep
cols2keep
# Retai columns used for PCA
data4pca <- d_traits_phy[,c("group","species",cols2keep)]
dim(data4pca) ; head(data4pca)


### 21/10/2020: Use data4pca to test variations of niche traits across groups
colnames(data4pca) <- str_replace(as.character(colnames(data4pca)), "center2", "center")
colnames(data4pca) <- str_replace(as.character(colnames(data4pca)), "width2", "width")
# In a lapply, perform non parametric tests of variance (K-W) and return pvalue and Chi2 values
res.tests <- lapply(colnames(data4pca)[c(3:18)], function(v) {
            # v <- "SST_center"
            test <- kruskal.test(get(v) ~ factor(group), data = data4pca) # str(summary(test))
            # test$p.value
            return( data.frame(niche.trait = v, Chi2 = test$stat, pval = test$p.value) )
        } # eo fun
) # eo lapply
tests <- bind_rows(res.tests) ; rm(res.tests) ; gc()
rownames(tests) <- NULL
tests
### And select those with pval < 0.01
tests[tests$pval <= 0.01,]



### Check how to add weights properly in PCA
# ?PCA
# Determine the predictors mean weight.pca (from)
vars.weights <- data.frame(traits_phy %>% group_by(var) %>% summarize(w = mean(weight.pca), probs = mean(prob.range)) )
vars.weights
vars.weights$w2 <- vars.weights$w/max(vars.weights$w)
#vars.weights[order(vars.weights$w2, decreasing = T),] # OK 
### colnames in data4pca and variables in vars.weight$w2 follow the same order, so just double the vars.weight$w2 vector
col.weights <- c(vars.weights$w2,vars.weights$w2) ; col.weights

pca.phy <- PCA(data4pca[,c(3:18)], scale.unit = T, ncp = 10, graph = F, col.w = col.weights)
summary(pca.phy)
str(pca.phy)
# Provide the PC coords to data4pca
data4pca[,paste("PC",c(1:5),sep="")] <- pca.phy$ind$coord[,c(1:5)]
summary(data4pca)

eig <- data.frame(perc = pca.phy$eig[,"percentage of variance"], nb = c(1:nrow(pca.phy$eig)) ) # eig
pca1 <- paste0("PC1 (",floor(eig$perc[1]*100)/100,"%)")
pca2 <- paste0("PC2 (",floor(eig$perc[2]*100)/100,"%)")
pca3 <- paste0("PC3 (",floor(eig$perc[3]*100)/100,"%)")
pca4 <- paste0("PC4 (",floor(eig$perc[4]*100)/100,"%)")
pca5 <- paste0("PC5 (",floor(eig$perc[5]*100)/100,"%)")


ggplot(aes(x = factor(group), y = PC1, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(pca1) + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()
#
ggplot(aes(x = factor(group), y = PC2, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(pca2) + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()
#
ggplot(aes(x = factor(group), y = PC3, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(pca3) + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()
#
ggplot(aes(x = factor(group), y = PC4, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(pca4) + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()

### Compute the groups' centroid in PCA space
groups.coords <- data.frame(data4pca %>% group_by(group) %>% summarize(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), PC4 = mean(PC4)) )
groups.coords

### PC plot: plot the species coordinates in PC space and report the centroid of the 3 groups
ggplot() + geom_point(aes(x = PC1, y = PC2, fill = factor(group)), data = data4pca, pch = 21, colour = "black", alpha = .5) + 
    geom_point(aes(x = PC1, y = PC2, fill = factor(group)), data = groups.coords, pch = 21, colour = "black", size = 4) +
    scale_fill_brewer(name = "", palette = "Paired") + xlab(pca1) + ylab(pca2) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    theme_classic()
#
### Same with PC3 and PC4
ggplot() + geom_point(aes(x = PC3, y = PC4, fill = factor(group)), data = data4pca, pch = 21, colour = "black", alpha = .5) + 
    geom_point(aes(x = PC3, y = PC4, fill = factor(group)), data = groups.coords, pch = 21, colour = "black", size = 4) + 
    scale_fill_brewer(name = "", palette = "Paired") + xlab(pca3) + ylab(pca4) + 
    geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    theme_classic()

### 09/10/2020: Replot PCA space with variables
library("ggrepel")
library("broom")

### Make sure the dims arguments is as long as the nb of PCs retained

augment.PCA <- function(x, dims = c(1:6), which="col") {
  .get <- function(x, element, dims) {
    y <- as.data.frame(x[[element]]$coord[,dims])
    if (nrow(y) == 0) {
      y <- NULL
    } else {
      y$type <- element
    }
    return(y)
  }
  if (which == "col") {
    y <- rbind(.get(x, "var", dims), .get(x, "quanti.sup", dims))
  } else {
    y <- rbind(.get(x, "ind", dims), .get(x, "quali.sup", dims))
  }
  y$var <- row.names(y)
  row.names(y) <- NULL
  return(y)
}

pcad <- augment.PCA(pca.phy)
pcad$var # Re-name properly
pcad$var <- str_replace_all(pcad$var, "_", " ")
#pcad$var <- str_replace_all(pcad$var, "center2", "center")
#pcad$var <- str_replace_all(pcad$var, "width2", "width")

ggplot(pcad) +
    coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
    annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
    geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
    scale_colour_manual(name = "", values = c("#4d9221","#c51b7d"), guide = F) + 
    geom_text_repel(aes(x=Dim.1, y=Dim.2, label=var), 
            data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
    xlab(pca1) + ylab(pca2) + theme_bw()

ggplot(pcad) +
    coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
    annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
    geom_segment(aes(x=0, xend = Dim.2, y=0, yend = Dim.3), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
    scale_colour_manual(name = "", values = c("#4d9221","#c51b7d"), guide = F) + 
    geom_text_repel(aes(x=Dim.2, y=Dim.3, label = var), 
            data=filter(pcad, (Dim.2^2+Dim.3^2) > 0.2^2), segment.alpha=0.5) +
    xlab(pca3) + ylab(pca4) + theme_bw()


    
# ----------------------------------------------------------------

### 2°) Zooplankton 
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/resp_curves/niche_traits/")
files <- dir()[grep("niche_traits_GAM",dir())] ; files
require("parallel")
res <- mclapply(files, function(f) {
            t <- get(load(f))
            return(t)
    }, mc.cores = 20
) # eo mclapply - f in files
# Rbind
traits_zoo <- bind_rows(res)
rm(res,files) ; gc()
head(traits_zoo) ; dim(traits_zoo)

### Print those plots per variable
vars <- unique(traits_zoo$var) ; vars
v <- "SST"

for(v in vars) {
    
    message(paste("Plotting niche traits ditsribution for ",v, sep = ""))
    
    p1 <- ggplot(aes(x = factor(group), y = center2, fill = factor(group)), data = traits_zoo[traits_zoo$var == v,]) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(paste(v," niche center",sep="")) + scale_x_discrete(labels = NULL) +
        theme_classic()
    #
    p2 <- ggplot(aes(x = factor(group), y = width2, fill = factor(group)), data = traits_zoo[traits_zoo$var == v,]) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(paste(v," niche width",sep="")) + scale_x_discrete(labels = NULL) + 
        theme_classic()
    
    p3 <- ggplot(aes(x = center2, y = width2, fill = factor(group)), data = traits_zoo[traits_zoo$var == v,]) + 
        geom_point(colour = "black", pch = 21) + scale_fill_brewer(name = "", palette = "Paired") + 
        ylab(paste(v," niche width",sep="")) + xlab(paste(v," niche center",sep="")) + theme_classic()
    
    p4 <- ggplot(aes(x = factor(group), y = prob.range, fill = factor(group)), data = traits_zoo[traits_zoo$var == v,]) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab(paste("HSI range (importance)",sep="")) + scale_x_discrete(labels = NULL) + 
        scale_y_continuous(limits = c(0,1)) + theme_classic()
    
    require("ggpubr")
    panel <- ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, align = "hv")
    # panel
    setwd(WD)
    ggsave(plot = panel, filename = paste("panel_niche_traits_zoo_",v,".jpg", sep = ""), width = 10, height = 6, dpi = 300)
    
} # eo v in vars

### Plot range if HSI per var and group (facet per vars)
# head(traits_zoo)
# ggplot(aes(x = factor(group), y = weight.pca, fill = factor(group)), data = traits_zoo) +
#     geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Spectral") +
#      ylab("Importance in constraining HSI") + xlab("") + theme_classic() +
#      scale_x_discrete(labels = NULL) + scale_y_continuous(limits = c(0,1)) +
#      facet_wrap(~factor(var), nrow = 2, ncol = 4, scales = "fixed")

### Prepare data for PCA (analysis in multivariate env space)
m_traits_zoo <- melt(traits_zoo, id.vars = c("group","species","var"))
# Create every combination of var x variable in a factor column and dcast
m_traits_zoo$new.col <- factor(paste(m_traits_zoo$var ,m_traits_zoo$variable, sep = "_"))
unique(m_traits_zoo$new.col)
# Dcast to put those new.col as columns (return to wide format)
d_traits_zoo <- dcast(m_traits_zoo[,c(1,2,5,6)], group + species ~ new.col, value.var = "value")
 
# Select the columns corresponding to: center2, width2 and weight.pca (weights to give to the columns, not the individuals)
cols2keep <- c(colnames(d_traits_zoo)[grep("center2", colnames(d_traits_zoo))],
                colnames(d_traits_zoo)[grep("width2", colnames(d_traits_zoo))],
                colnames(d_traits_zoo)[grep("weight.pca", colnames(d_traits_zoo))]
) # eo cols2keep
cols2keep

data4pca <- d_traits_zoo[,c("group","species",cols2keep)]
dim(data4pca)
head(data4pca)

### 21/10/2020: Use 'data4pca' to test variantions of niche traits across groups of zooplankton
colnames(data4pca) <- str_replace(as.character(colnames(data4pca)), "center2", "center")
colnames(data4pca) <- str_replace(as.character(colnames(data4pca)), "width2", "width")
# In a lapply, perform non parametric tests of variance (K-W) and return pvalue and Chi2 values
res.tests <- lapply(colnames(data4pca)[c(3:18)], function(v) {
            # v <- "SST_center"
            test <- kruskal.test(get(v) ~ factor(group), data = data4pca) # str(summary(test))
            # test$p.value
            return( data.frame(niche.trait = v, Chi2 = test$stat, pval = test$p.value) )
        } # eo fun
) # eo lapply
tests <- bind_rows(res.tests) ; rm(res.tests) ; gc()
rownames(tests) <- NULL
tests
### And select those with pval < 0.01
tests[tests$pval <= 0.001,]


### Check how to add weights properly in PCA
vars.weights <- data.frame(traits_zoo %>% group_by(var) %>% summarize(w = mean(weight.pca)) )
vars.weights$w2 <- vars.weights$w/max(vars.weights$w)
#vars.weights[order(vars.weights$w2, decreasing = T),] # OK 
### colnames in data4pca and variables in vars.weight$w2 follow the same order, so just double the vars.weight$w2 vector
col.weights <- c(vars.weights$w2,vars.weights$w2) ; col.weights

pca.zoo <- PCA(data4pca[,c(3:18)], scale.unit = T, ncp = 10, graph = F, col.w = col.weights)
summary(pca.zoo)

eig <- data.frame(perc = pca.zoo$eig[,"percentage of variance"], nb = c(1:nrow(pca.zoo$eig)) ) # eig
pca1 <- paste0("PC1 (",floor(eig$perc[1]*100)/100,"%)")
pca2 <- paste0("PC2 (",floor(eig$perc[2]*100)/100,"%)")
pca3 <- paste0("PC3 (",floor(eig$perc[3]*100)/100,"%)")
pca4 <- paste0("PC4 (",floor(eig$perc[4]*100)/100,"%)")
pca5 <- paste0("PC5 (",floor(eig$perc[5]*100)/100,"%)")


# Provide the PC coords to data4pca
data4pca[,paste("PC",c(1:5),sep="")] <- pca.zoo$ind$coord[,c(1:5)]
summary(data4pca)

ggplot(aes(x = factor(group), y = PC1, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab("PC1") + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()
#
ggplot(aes(x = factor(group), y = PC2, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab("PC2") + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()
#
ggplot(aes(x = factor(group), y = PC3, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab("PC3") + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()
#
ggplot(aes(x = factor(group), y = PC4, fill = factor(group)), data = data4pca) + 
        geom_boxplot(colour = "black") + scale_fill_brewer(name = "", palette = "Paired") + 
        xlab("") + ylab("PC4") + scale_x_discrete(labels = NULL) +
        geom_hline(yintercept = 0, linetype = "dashed") + theme_classic()

### Compute the groups' centroid in PCA space
groups.coords <- data.frame(data4pca %>% group_by(group) %>% summarize(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), PC4 = mean(PC4)) )
groups.coords

### PC plot: plot the species coordinates in PC space and report the centroid of the 3 groups
ggplot() + geom_point(aes(x = PC1, y = PC2, fill = factor(group)), data = data4pca, pch = 21, colour = "black", alpha = .5) + 
    scale_fill_brewer(name = "", palette = "Paired") + 
    geom_point(aes(x = PC1, y = PC2, fill = factor(group)), data = groups.coords, pch = 21, colour = "black", size = 4) + 
    xlab("PC1") + ylab("PC2") + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    theme_classic()

### Same with PC3 and PC4
ggplot() + geom_point(aes(x = PC3, y = PC4, fill = factor(group)), data = data4pca, pch = 21, colour = "black", alpha = .5) + 
    scale_fill_brewer(name = "", palette = "Paired") + 
    geom_point(aes(x = PC3, y = PC4, fill = factor(group)), data = groups.coords, pch = 21, colour = "black", size = 4) + 
    xlab("PC3") + ylab("PC4") + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") + 
    theme_classic()
    
###
pcad <- augment.PCA(pca.zoo)
pcad$var # Re-name properly
pcad$var <- str_replace_all(pcad$var, "_", " ")
pcad$var <- str_replace_all(pcad$var, "center2", "center")
pcad$var <- str_replace_all(pcad$var, "width2", "width")

ggplot(pcad) +
    coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
    annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
    geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
    scale_colour_manual(name = "", values = c("#4d9221","#c51b7d"), guide = F) + 
    geom_text_repel(aes(x=Dim.1, y=Dim.2, label=var), 
            data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
    xlab(pca1) + ylab(pca2) + theme_bw()
#
ggplot(pcad) +
    coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
    annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
    geom_segment(aes(x=0, xend = Dim.2, y=0, yend = Dim.3), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
    scale_colour_manual(name = "", values = c("#4d9221","#c51b7d"), guide = F) + 
    geom_text_repel(aes(x=Dim.2, y=Dim.3, label = var), 
            data=filter(pcad, (Dim.2^2+Dim.3^2) > 0.2^2), segment.alpha=0.5) +
    xlab(pca3) + ylab(pca4) + theme_bw()


# --------------------------------------------------------------------------------------------------------------------------------

### 21/10/2020: 

setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/resp_curves/niche_traits/")
files <- dir()[grep("resp_curves_GAM",dir())] # files
# Concatenate all niche traits
require("parallel")
res <- mclapply(files, function(f) {
            t <- get(load(f))
            return(t)
    }, mc.cores = 20
) # eo mclapply - f in files
# Rbind
traits_phy <- bind_rows(res)
rm(res,files) ; gc()
head(traits_phy) ; dim(traits_phy)

### Make sur SST (or another variable) ranges are the same across species
unique(traits_phy$species)
summary(traits_phy[traits_phy$species == "Chaetoceros_laciniosus","SST"])
summary(traits_phy[traits_phy$species == "Tripos_pulchellus","SST"])


