
##### 05/11/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Analyzing the drivers of changes in SR across groups by assessing their correlations with ens annual changes in ∆ predictors
#   - For SST only scenario
#   - For constant SST scenario
 
### Last update: 09/11/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("parallel")

WD <- getwd()
world2 <- map_data("world2")

SDMs <- c("GLM","GAM","ANN")
ESMs <- c("CNRM-PISCES","IPSL-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")		

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Load the table containing the 
setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
table.diff <- get(load("table_ann_diff_ens_constantSST.vs.SSTonly_groups.Rdata"))
head(table.diff)
d.table.diff <- dcast(data = table.diff[,c(1:4,6:7)], formula = cell_id + x + y ~ group + scenario, value.var = "perc")
groups <- unique(table.diff$group)


### 2°) Load the monthly climatologies of model deltas for all 5 ESMs and compute ensemble average deltas (annual)
setwd("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/mon_clims")
files <- dir()[grep("_diff_",dir())][c(6:60)]  # remove the first 5 because not used in the groups' projections (annual range of SST)
files

# Load one for test
#d <- get(load(files[43]))
#dim(d)
#head(d) # 1 column per month, easily derive mean annual diff with rowMeans

# Define vector or preds names as they are written in the file names
vars <- c("sst","o2","par","logchl","logno3","logsio2","nstar","sistar")
# v <- "o2"
ann.clims <- mclapply(vars, function(v) {
    
             message(paste("Computing mean annual delta for ", v, sep = ""))
             files2 <- files[grep(paste("diff_",v,"_rcp85",sep=""),files)]
             
             if(length(files2) == 5) {
                 # f <- files2[3]
                 res.esm <- lapply(files2, function(f) {
                             d <- get(load(f))
                             # Derive annual mean
                             d$Annual <- rowMeans(as.matrix(d[,c(4:length(d))]), na.rm = T)
                             d$variable <- v
                             return(d[,c("id","x","y","Annual","variable")])    
                         }
                 ) # eo lapply
                 # Rbind
                 ddf.esm <- bind_rows(res.esm)
                 rm(res.esm) ; gc()
                 return(ddf.esm)
                 
             } else {
                 
                 message(paste("!ERROR OwO", sep = ""))
                 
             }
    
    }, mc.cores = length(vars)   
    
) # eo mclapply - v in vars
# Rbind
ddf.clims <- bind_rows(ann.clims)
dim(ddf.clims)
head(ddf.clims)
summary(ddf.clims)
rm(ann.clims) ; gc()

### Derive mean annual ensemble across the 5 ESM with dplyr
ddf.ens <- data.frame(ddf.clims %>% group_by(id,variable) %>% summarize(x = unique(x)[1], y = unique(y)[1], mean.delta = mean(Annual, na.rm = T)) )
dim(ddf.ens)
head(ddf.ens)
summary(ddf.ens)
rm(ddf.clims); gc()

### And dcast
table.ens <- dcast(data = ddf.ens, formula = id + x + y ~ variable, value.var = "mean.delta")
dim(table.ens)
head(table.ens)
summary(table.ens)
# Re-order according to id
table.ens <- table.ens[order(table.ens$id),]
# And rename the cols so they match previous names
colnames(table.ens)[(4:11)] <- c("logChl","logNO3","logSiO2","N*","O2","PAR","Si*","SST")


### 3°) Great ! Now, for each group, extract the 'All scenario' projections of ∆perc in SR from 'd.table.diff', cbind with 'table.ens'
### and compute spearman's corr coeff. And geom_point 


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

# head(d.table.diff) ; head(table.ens)

g <- "Diatoms"

cormats <- lapply(groups, function(g) {
    
    if(g %in% c("Diatoms","Dinoflagellates","Coccolithophores")) {
        preds <- c("PAR","logSiO2","logChl","N*","Si*")
    } else {
        preds <- c("O2","logSiO2","logChl","PAR","Si*")
    }
    ### NOTE: removed "SST" for constant SST projections in order not to confuse readerand myself
    
    names <- colnames(d.table.diff)[grep(g,colnames(d.table.diff))]
    sub <- na.omit(d.table.diff[,c("cell_id","x","y",names[1])])
    
    # Find common cells
    commons <- intersect(unique(sub$cell_id), unique(table.ens$id)) # length(commons)
    sb2 <- sub[sub$cell_id %in% commons,]
    table.ens2 <- table.ens[table.ens$id %in% commons,]
    sb2 <- sb2[order(sb2$cell_id),]
    table.ens2 <- table.ens2[order(table.ens2$id),]
    # dim(sb2) ; dim(table.ens2)
    
    # Cbind
    sb3 <- cbind(sb2, table.ens2[,preds])
    # dim(sb3) ; head(sb3) ; summary(sb3)
    rm(sb2, table.ens2) ; gc()
    
    # Compute pairwise correlation heatmap
    names <- colnames(sb3)[c(4:length(sb3))] #; names
    mydata <- na.omit(sb3[,names])
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
    cormat <- cor_tri[which(cor_tri$Var1 %in% names[1] & cor_tri$Var2 %in% preds),]
    cormat <- cormat[!(cormat$Var1 == cormat$Var2),]
    # head(cormat)

    # Add signif label
    cormat$signif <- NA
    cormat[cormat$pval <= 0.001,"signif"] <- "***"
    cormat[cormat$pval > 0.001 & cormat$pval <= 0.01,"signif"] <- "**"
    cormat[cormat$pval > 0.01 & cormat$pval < 0.05,"signif"] <- "*"
    cormat[cormat$pval >= 0.05,"signif"] <- "-"
    colnames(cormat)[c(1,2)] <- c("Group","Covariate")
    cormat$Group <- str_replace_all(cormat$Group, "_Constant SST", "")
    
    return(cormat)
    
    }
    
) # eo lapply
# Rbind
cormat <- bind_rows(cormats)
rm(cormats) ; gc()

plot <- ggplot(cormat, aes(factor(Covariate), factor(Group), fill = rho)) + geom_tile(color = "white") +
        scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), 
                name = paste("Spearman' rank\ncorrelation\ncoefficient", sep = ""), breaks = c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1), 
                guide = guide_colourbar(ticks = T, ticks.colour = "black", frame.colour = "black")) + 
        xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1), 
                axis.text.y = element_text(vjust = 1, size = 10, hjust = 1))+
        coord_fixed() + geom_text(aes(factor(Covariate), factor(Group), label = rho), color = "black", size = 4)
#
setwd(WD)
ggsave(plot = plot, filename = paste("heatmap_corr_ann_ens_perc_deltas_groups_constantSST.jpg", sep = ""), dpi = 300, width = 7, height = 9)
### Not the greatest test indeed --> need to launch simualtions wiht constant SST and THEN look at correlations

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
