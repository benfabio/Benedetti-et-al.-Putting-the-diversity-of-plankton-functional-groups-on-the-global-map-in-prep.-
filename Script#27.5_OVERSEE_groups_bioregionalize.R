
##### 05/10/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Loading monthly group-level (group-background) community composition estimates for all SDMs and derive mean annual comm
#   - Cbind all phyto and zoo estimates on annual scale
#   - Perform PCA to reduce dimensionality and summarize spatial variations in community composition
#   - Use clValid on PC scores to estimate best clustering method and N clusters (pam, hierarchical; internal, stability indices)
 
### Last update: 14/10/2020 

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("parallel")
library("clValid")
# ?clValid

WD <- getwd()
world2 <- map_data("world2")
world <- map_data("world")

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Load monthly SDM communities and derive mean annual probabilities for each species

### A) Phytoplankton 
# setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities")
# comm <- read.table(paste("table_mon_composition_","Phytoplankton","_","group",".txt", sep = ""), h = T, sep = "\t")
# head(comm) ; dim(comm) ; colnames(comm)
# # Melt priori to averaging probabilities
# m.comm <- melt(comm, id.vars = c("cell_id","x","y","SDM","month"))
# #head(m.comm)
# ann.comm <- data.frame(m.comm %>% group_by(cell_id, variable) %>% summarize(x = unique(x), y = unique(y), mean = mean(value, na.rm = T)) )
# head(ann.comm)
# # And dcast
# d.ann.comm <- dcast(ann.comm, cell_id + x + y ~ variable, fun.aggregate = mean, na.rm = T, value.var = "mean")
# dim(d.ann.comm)
# summary(d.ann.comm) ; head(d.ann.comm)
# # Save it somewhere
# write.table(d.ann.comm, file = "table_annual_composition_ensemble_Phytoplankton_group.txt", sep = "\t")


# ### B) Zooplankton
# setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities")
# comm <- read.table(paste("table_mon_composition_","Zooplankton","_","group",".txt", sep = ""), h = T, sep = "\t")
# head(comm) ; dim(comm)
# # Melt priori to averaging probabilities
# m.comm <- melt(comm, id.vars = c("cell_id","x","y","SDM","month"))
# # Derive mean annual probabilities
# ann.comm <- data.frame(m.comm %>% group_by(cell_id, variable) %>% summarize(x = unique(x), y = unique(y), mean = mean(value, na.rm = T)) )
# head(ann.comm) ; summary(ann.comm)
# # And dcast
# d.ann.comm <- dcast(ann.comm, cell_id + x + y ~ variable, fun.aggregate = mean, na.rm = T, value.var = "mean")
# dim(d.ann.comm)
# summary(d.ann.comm)
# write.table(d.ann.comm, file = "table_annual_composition_ensemble_Zooplankton_group.txt", sep = "\t")

# --------------------------------------------------------------------------------------------------------------------------------

### 2°) Load above-created tables of mean annual community composition and perform PCAs
# - total plankton
# - phytoplankton
# - zooplankton

setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities")
phyto <- read.table("table_annual_composition_ensemble_Phytoplankton_group.txt", h = T, sep = "\t")
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities")
zoo <- read.table("table_annual_composition_ensemble_Zooplankton_group.txt", h = T, sep = "\t")
#dim(phyto) ; dim(zoo)
phyto <- phyto[order(phyto$cell_id),]
zoo <- zoo[order(zoo$cell_id),]

library("FactoMineR")
data.pca.phyto <- na.omit(phyto)
pca.phyto <- PCA(X = data.pca.phyto[,c(4:length(data.pca.phyto))], scale.unit = T, ncp = 10, graph = F)
#summary(pca.phyto)
# eig <- data.frame(perc = pca.phyto$eig[,"eigenvalue"], nb = c(1:nrow(pca.phyto$eig)) )
# head(eig) ; dim(eig)
# avg <- mean(eig$perc) # avg
# ggplot(data = eig) + geom_bar(aes(x=nb, y=perc), stat="identity") + geom_line(aes(x=nb, y= mean(perc) )) +
#     xlab("PC") +  ylab("% of variance explained") + theme_classic()

data.pca.zoo <- na.omit(zoo)
pca.zoo <- PCA(X = na.omit(data.pca.zoo[,c(4:length(data.pca.zoo))]), scale.unit = T, ncp = 10, graph = F)
# summary(pca.zoo)
# head(eig)
# eig <- data.frame(perc = pca.zoo$eig[,"percentage of variance"], nb = c(1:nrow(pca.zoo$eig)) )
# ggplot(data = eig) + geom_bar(aes(x=nb, y=perc), stat="identity") + geom_line(aes(x=nb, y=mean(perc) )) +
#     xlab("PC") +  ylab("% of variance explained") + theme_classic()

### For both, keep the first 5 dimensions
coords.phyto <- data.frame(pca.phyto$ind$coord[,1:5])
coords.zoo <- data.frame(pca.zoo$ind$coord[,1:5])
# dim(coords.phyto) ; dim(coords.zoo)

### Perform clValid
# Phytoplankton first
methods <- c("average","ward")

# for(m in methods) {
#
#     valid.phyto <- clValid(obj = coords.phyto, nClust = c(4:13), clMethods = c("kmeans"),
#                     validation = "internal", metric = "euclidean",
#                     verbose = T, maxitems = nrow(coords.phyto))
#     #
#     setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/clustering/clValid")
#     save(valid.phyto, file = paste("clValid_kmeans_phyto_annual_ens_PC1-5_",m,".Rdata", sep = ""))
#
#     # Zooplankton second
#     valid.zoo <- clValid(obj = coords.zoo, nClust = c(4:13), clMethods = c("kmeans"),
#                     validation = "internal", metric = "euclidean",
#                     verbose = T, maxitems = nrow(coords.zoo))
#
#     ### Save
#     setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/clustering/clValid")
#     save(valid.zoo, file = paste("clValid_kmeans_zoo_annual_ens_PC1-5_",m,".Rdata", sep = ""))
#
# } # eo for loop



### 3°) Combine data.pca.phyto and data.pca.zoo by common cells and perform PCA on all plankton species 
# dim(data.pca.phyto) ; dim(data.pca.zoo)
commons <- intersect(unique(data.pca.phyto$cell_id),unique(data.pca.zoo$cell_id))
# length(commons)
data.pca.phyto2 <- data.pca.phyto[data.pca.phyto$cell_id %in% commons,]
data.pca.zoo2 <- data.pca.zoo[data.pca.zoo$cell_id %in% commons,]
data.pca.all <- cbind(data.pca.phyto2,data.pca.zoo2[,c(4:length(data.pca.zoo2))])
# dim(data.pca.all) ; head(data.pca.all)
rm(data.pca.phyto2,data.pca.zoo2,commons) ; gc()

# Perform PCA
pca.plankton <- PCA(X = data.pca.all[,c(4:length(data.pca.all))], scale.unit = T, ncp = 10, graph = F)
# summary(pca.plankton)
coords.all <- data.frame(pca.plankton$ind$coord[,1:5])

for(m in methods) {
    
    valid.all <- clValid(obj = coords.all, nClust = c(4:13), clMethods = c("hierarchical"),
                    method = m, validation = "internal", metric = "euclidean", 
                    verbose = T, maxitems = nrow(coords.all))

    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/clustering/clValid")
    save(valid.all, file = paste("clValid_plankton_annual_ens_PC1-5_",m,".Rdata", sep = ""))
    
} # eo for m in methods


# --------------------------------------------------------------------------------------------------------------------------------

### 06/10/2020: While clValid is still running, examine some bioregionalizations:
# - compute euclidean distance matrix from the PC1-5 coordinates
# - draw dendrogram (print on computer, no display) 
# - plot the outputs from various cutting levels

# # Distance matrices
# dist.phyto <- dist(coords.phyto, "euclidean")
# dist.zoo <- dist(coords.zoo, "euclidean")
# dist.all <- dist(coords.all, "euclidean")

# # Dendrogram (keep Ward's link as default)
# clust.phy <- hclust(dist.phyto, "ward.D2")
# clust.zoo <- hclust(dist.zoo, "ward.D2")
# clust.all <- hclust(dist.all, "ward.D2")
# # class(clust.phy)

# # Basic plotting (NOTE: ggdendro too slow because of the amount of tips)
# # library("ggplot2")
# # library("ggdendro")
# # dendroplot <- ggdendrogram(clust.zoo)
# # class(dendroplot)

# setwd(WD)

# pdf(file = "dendro_phyto_annual_ens_PC1-5.pdf", width = 10, height = 8)
# plot(clust.phy, hang = -1, labels = F)
# dev.off()

# pdf(file = "dendro_zoo_annual_ens_PC1-5.pdf", width = 10, height = 8)
# plot(clust.zoo, hang = -1, labels = F)
# dev.off()

# pdf(file = "dendro_all_annual_ens_PC1-5.pdf", width = 10, height = 8)
# plot(clust.all, hang = -1, labels = F)
# dev.off()

# ### Cut at k = 7 just to check
# # length(cutree(clust.all, k = 7)) ; colnames(data.pca.all)
# data.pca.phyto$k7 <- cutree(clust.phy, k = 10) # kind of pointless beyond 10
# data.pca.zoo$k7 <- cutree(clust.zoo, k = 13) # kind of pointless beyond 11
# data.pca.all$k7 <- cutree(clust.all, k = 10)
# # Check n per k
# # summary(factor(data.pca.phyto$k7)) ; summary(factor(data.pca.zoo$k7)) ; summary(factor(data.pca.all$k7))
# # Map
# ggplot() + geom_raster(aes(x = x, y = y, fill = factor(k7)), data = data.pca.zoo) +
#     scale_fill_brewer(name = "Bioregions", palette = "Paired") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "black", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#         panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) +
#     scale_x_continuous(name = "Longitude", limits = c(-180,180), expand = c(0,0), labels = NULL, breaks = NULL) +
#     scale_y_continuous(name = "Latitude", limits = c(-85,85), expand = c(0,0), labels = NULL, breaks = NULL)
#
# ### Check out potential indicator species values
# require("labdsv")
# ?indval
# indval.all <- indval(x = data.pca.all[,c(4:(length(data.pca.all)-1))], clustering = data.pca.all$k7)
# #class(indval.all) ; str(indval.all) ; summary(indval.all)
# #indval.all$indval
# inds.all <- melt(indval.all$maxcls)
# inds.all$species <- rownames(inds.all)
# colnames(inds.all)[1] <- "Bioregion" ; rownames(inds.all) <- NULL
# inds.all[inds.all$Bioregion == 5,]


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
