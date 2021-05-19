
##### 21/08/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	 - Examining diversity monthly and annual patterns (base and fut % changes) after aggregating per group/ order/ family
#    - Examine seasonality across the various groups?
#    - Examine (per lat band for instance) the median changes in HSI, at group/order/family levels
#    - Examine the groups that correlate best with NPP/ Cexp/ Plankton size (baseline conditions) 
 
### Last update: 24/08/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("sp")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")
library("viridis")
library("scales")
library("maps")
library("FactoMineR")
library("parallel")
library("ggpubr")

world2 <- map_data("world2")
WD <- getwd() 

# --------------------------------------------------------------------------------------------------------------------------------

### Load phyto and zoo classification
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology") ; dir()
classif <- read.csv("classif_all_modelled_taxa.csv", h = T, sep = ";")
head(classif) ; str(classif)
# 324 genera
# 147 families
# 71 orders
# 24 classes
unique(classif$class)
# One Larvacea should be an Appendicularia
levels(classif$class)[levels(classif$class)=="Larvacea"] <- "Appendicularia"


### Examine contribution of groups etc. with donut plots after tallying number of species
nsp.phy <- length(unique(classif[classif$Category == "Phytoplankton","species"]))
nsp.zoo <- length(unique(classif[classif$Category == "Zooplankton","species"]))
nsp.phy ; nsp.zoo

groups.phyto <- data.frame(classif[classif$Category == "Phytoplankton",] %>% group_by(group) %>% summarize(Nsp = n(), Psp = (n()/nsp.phy)*100))
groups.phyto[order(groups.phyto$Psp, decreasing = T),]
plot.phyto <- ggplot(groups.phyto[groups.phyto$Psp > 1,], aes(x = 2, y = Psp, fill = factor(group))) +
  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(name = "group", palette = "Paired") + theme_void() + xlim(0.5,2.5) 
#
groups.zoo <- data.frame(classif[classif$Category == "Zooplankton",] %>% group_by(group) %>% summarize(Nsp = n(), Psp = (n()/nsp.zoo)*100))
groups.zoo[order(groups.zoo$Psp, decreasing = T),]
plot.zoo <- ggplot(groups.zoo[groups.zoo$Psp > 1,], aes(x = 2, y = Psp, fill = factor(group))) +
  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(name = "group", palette = "Paired") + theme_void() + xlim(0.5,2.5) 

# Same at class level
class.phyto <- data.frame(classif[classif$Category == "Phytoplankton",] %>% group_by(class) %>% summarize(Nsp = n(), Psp = (n()/nsp.phy)*100))
class.zoo <- data.frame(classif[classif$Category == "Zooplankton",] %>% group_by(class) %>% summarize(Nsp = n(), Psp = (n()/nsp.zoo)*100))
class.phyto[order(class.phyto$Psp, decreasing = T),]
class.zoo[order(class.zoo$Psp, decreasing = T),]

# Same at order level
order.phyto <- data.frame(classif[classif$Category == "Phytoplankton",] %>% group_by(order) %>% summarize(Nsp = n(), Psp = (n()/nsp.phy)*100))
order.zoo <- data.frame(classif[classif$Category == "Zooplankton",] %>% group_by(order) %>% summarize(Nsp = n(), Psp = (n()/nsp.zoo)*100))
order.phyto[order(order.phyto$Psp, decreasing = T),]
order.zoo[order(order.zoo$Psp, decreasing = T),]

# Family level?
fam.phyto <- data.frame(classif[classif$Category == "Phytoplankton",] %>% group_by(family) %>% summarize(Nsp = n(), Psp = (n()/nsp.phy)*100))
fam.zoo <- data.frame(classif[classif$Category == "Zooplankton",] %>% group_by(family) %>% summarize(Nsp = n(), Psp = (n()/nsp.zoo)*100))
fam.phyto[order(fam.phyto$Psp, decreasing = T),]
fam.zoo[order(fam.zoo$Psp, decreasing = T),]

### Need to find an approprate way to break down those groups into meaningful functional units
### (groups of taxa that broadly perform the same function in the ecosystem)

### A) Phytoplankton: harder than for zooplankton, but essentially it seems like you have the following main groups
# Diatoms 
# Dinoflagellates
# Coccos
# Phaeocystis?

### And then break down per main orders/ main families?

### B) Zooplankton: easier
# Copepods --> need to break down (knowing that a dedicated study is also undergoing)
nsp.cops <- length(unique(classif[classif$group == "Copepoda","species"]))
data.frame(classif[classif$group == "Copepoda",] %>% group_by(order) %>% summarize(Nsp = n(), Psp = (n()/nsp.cops)*100))
# Calanoida vs. Cyclopoida (Poecilostomatoida vs. Oithonida)
### --> requires more thought, also depends on the nb of species modelled


# Malacostraca --> need to break down 
nsp.malaco <- length(unique(classif[classif$group == "Malacostraca","species"]))
data.frame(classif[classif$group == "Malacostraca",] %>% group_by(order) %>% summarize(Nsp = n(), Psp = (n()/nsp.malaco)*100))
# So basically Amphipoda vs. Euphausiids, good

# Jellyfish --> need to break down ? (gelatinous passive ambush feeders)
nsp.jel <- length(unique(classif[classif$group == "Jellyfish","species"]))
data.frame(classif[classif$group == "Jellyfish",] %>% group_by(order) %>% summarize(Nsp = n(), Psp = (n()/nsp.jel)*100))
# Hard to tell...keep as such

# Chaetognatha --> keep as such, coherent from a functional point of view (gelatinous active ambush predators)
# Foraminifera --> keep as such, coherent from a functional point of view (passive ambush feeders, photosymbionts and calcareous test)
# Pteropoda --> restrict to Thecosomata, coherent from a functional point of view (calcareous flux feeders)
# Chordata --> large gelatinous filter feeders that need to be broken down?
nsp.chord <- length(unique(classif[classif$group == "Chordata","species"]))
data.frame(classif[classif$group == "Chordata",] %>% group_by(order) %>% summarize(Nsp = n(), Psp = (n()/nsp.chord)*100))
# Doliod+Salp versus Appendicularia/Copelata

### So at least 10 to 11 groups
### Add a secondary group column in the classif data.frame
classif$group2 <- "Other"
classif[classif$group == "Bacillariophyceae","group2"] <- "Diatoms"
classif[classif$group == "Dinoflagellata","group2"] <- "Dinoflagellates"
classif[classif$class %in% c("Prymnesiophyceae","Coccolithophyceae"),"group2"] <- "Coccolithophores" 
# na.omit(classif[classif$group2 == "Coccolithophores","species"])
classif[classif$group == "Jellyfish","group2"] <- "Jellyfish"
classif[classif$group == "Chaetognatha","group2"] <- "Chaetognatha"
classif[classif$group == "Foraminifera","group2"] <- "Foraminifera"
classif[classif$order %in% c("Salpida","Doliolida"),"group2"] <- "Salps"
classif[classif$order %in% c("Copelata"),"group2"] <- "Appendicularians"
classif[classif$order == "Euphausiacea","group2"] <- "Euphausiids"
classif[classif$order == "Amphipoda","group2"] <- "Amphipods"
classif[classif$order == "Thecosomata","group2"] <- "Pteropods"
# And divide the copepods into 3 groups?
classif[classif$order == "Calanoida","group2"] <- "Calanoids"
classif[classif$family == "Oithonidae","group2"] <- "Oithonids"
classif[classif$family %in% c("Oncaeidae","Corycaeidae","Sapphirinidae"),"group2"] <- "Poecilostomatoids"
classif$group2 <- factor(classif$group2)
summary(classif$group2)
classif[is.na(classif$group2),] # can be drops 

classif2 <- classif[!is.na(classif$group2),]

nsp.phy <- length(unique(classif2[classif2$Category == "Phytoplankton","species"]))
nsp.zoo <- length(unique(classif2[classif2$Category == "Zooplankton","species"]))
nsp.phy ; nsp.zoo

groups.phyto <- data.frame(classif2[classif2$Category == "Phytoplankton",] %>% group_by(group2) %>% summarize(Nsp = n(), Psp = (n()/nsp.phy)*100))
groups.phyto[order(groups.phyto$Psp, decreasing = T),]
plot.phyto <- ggplot(groups.phyto[groups.phyto$Psp > 1,], aes(x = 2, y = Psp, fill = factor(group2))) +
  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(name = "Group", palette = "Paired") + theme_void() + xlim(0.5,2.5) 
#
groups.zoo <- data.frame(classif2[classif2$Category == "Zooplankton",] %>% group_by(group2) %>% summarize(Nsp = n(), Psp = (n()/nsp.zoo)*100))
groups.zoo[order(groups.zoo$Psp, decreasing = T),]
plot.zoo <- ggplot(groups.zoo[groups.zoo$Psp > 1,], aes(x = 2, y = Psp, fill = factor(group2))) +
  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(name = "Group", palette = "Paired") + theme_void() + xlim(0.5,2.5) 


# --------------------------------------------------------------------------------------------------------------------------------

### For each group in classif2, retreive the the corresponding monthly species composition and derive mean annual diversity
### for the baseline conditions first, then do the same for future and derive % change

groups <- unique(classif2$group2) ; groups # 14 groups, ok
g <- "Diatoms"

for(g in groups) {
    
        species <- classif2[classif2$group2 == g,"species"]
        category <- unique(classif2[classif2$group2 == g,"Category"])
        message(paste("Retrieving projections for ",g," (",length(species)," species) - ",category, sep = ""))
        setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections")
        # If else loop
        if(category == "Phytoplankton") {
            files <- dir()[grep("table_ann_compo_phyto_baseline_",dir())][c(1:12)]
        } else {
            files <- dir()[grep("table_ann_compo_zoo_baseline",dir())][c(1:12)]
        } # eo if else loop 
        
        require("parallel")
        res <- mclapply(files, function(f) {
                    # f <- files[1]
                    d <- get(load(f))
                    # Change some colnames
                    colnames(d) <- gsub(".", "", colnames(d), fixed = T)
                    # And remove bracket from 'species'
                    species <- gsub("\\(|\\)", "", species)
                    # restrict to species of interest # intersect(species,colnames(d))
                    d2 <- d[,c("cell_id","x","y",species)]
                    # rowSums
                    d2$SR <- rowSums(as.matrix(d2[,species])) # summary(d2$SR)
                    # divide by nb of modelled species
                    d2$percSR <- (d2$SR)/length(species) # summary(d2$percSR)
                    # retrieve predictor pool and SDM before returning data.frame
                    terms <- data.frame(t = unlist(str_split(as.character(f), "_"))) # terms
                    sdm <- terms[6,]
                    pool <- str_replace_all(as.character(terms[7,]),".Rdata","")
                    
                    ### !!! need to re-scame a bit for RF models: *3
                    if(sdm == "RF") {
                        d2$SR <- (d2$SR)*3.2
                        d2$percSR <- (d2$SR)/(length(species)*3.2) # summary(d2$percSR)
                    } # eo if loop
                    
                    # Return
                    ddf <- data.frame(group = g, pool = pool, SDM = sdm, d2[,c("cell_id","x","y","SR","percSR")]) # dim(ddf); head(ddf)
                    return(ddf)
            }, mc.cores = 20
        ) # eo lapply 
        # Rbind
        ddf <- bind_rows(res) ; rm(res) ; gc() # summary(ddf)
        
        setwd("/net/kryo/work/fabioben/OVERSEE/data/tables_composition_ensemble_rcp85/Individual_projections/maps")
        bin <- ceiling((ceiling(max(ddf$SR,na.rm=T)))/3)
        map <- ggplot(data=ddf) + geom_tile(aes(x = x, y = y, fill = SR) ) +
                    scale_fill_viridis(name = "Mean annual\nSR", limits = c(0,max(ddf$SR,na.rm=T))) +
                    geom_contour(colour = "grey75", binwidth = bin, size = 0.25, aes(x = x, y = y, z = SR) ) +
 	                geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
 	                coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
                        labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
 	                scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
 		                labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
   	                theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
 		                panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "none") + 
                    facet_grid(factor(pool) ~ factor(SDM)) + theme(legend.position = "bottom")
                    
        # Save maps
        ggsave(plot = map, filename = paste("map_facets_ann_SR_baseline_",g,"_noRF.jpg", sep = ""), dpi = 300, width = 14, height = 13)
        
    
} # eo 1st for loop


# --------------------------------------------------------------------------------------------------------------------------------

### Ok, plot distribution of the groups' TSS scores per SDM and predictor pool 

setwd( paste(WD,"/biology/species_data_v9v3.1/total_background/niche.modelling_future_rcp85_ensemble/", sep = "") )
zoo.wd <- getwd()
setwd( paste(WD,"/biology/phytoplankton_15_01_19/total_background/species_data/niche.modelling_future_rcp85_ensemble/", sep = "") )
phyto.wd <- getwd()
setwd(WD)

SDMs <- c('GLM','GAM','RF','ANN')
pools <- c("p1","p2","p3","p4")

### 2°) Plot distrbution of models'TSS values for zoo and phyto separately, using boxplots per SDM and facet per pool
setwd(zoo.wd)
res <- mclapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(zoo.wd,"/eval_score_",p,"/", sep = "") )
					zoo_spp <- str_replace_all(dir(), "eval_scores_", "")
					zoo_spp <- str_replace_all(zoo_spp, ".Rdata", "")
					scores <- lapply(zoo_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- do.call(rbind, scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}, mc.cores = 30
) # eo lapply
# Rbind
table_scores_zoo <- do.call(rbind,res)
head(table_scores_zoo) # No SDM indicator...
table_scores_zoo$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_zoo)), pattern = "_", n = 2)))[,1])
table_scores_zoo$kingdom <- "Zooplankton"
dim(table_scores_zoo); head(table_scores_zoo)
summary(table_scores_zoo)
rm(res)

### And for phyto
setwd(phyto.wd)
res <- mclapply(pools, function(p) {
  					message(paste("Loading scores for pool ",p, sep = ""))
  					setwd( paste(phyto.wd,"/eval_score_",p,"/", sep = "") )
					phyto_spp <- str_replace_all(dir(), "eval_scores_", "")
					phyto_spp <- str_replace_all(phyto_spp, ".Rdata", "")
					scores <- lapply(phyto_spp, function(sp) {
								message(paste("Loading scores for ",sp, sep = ""))
								s <- get(load(paste("eval_scores_",sp,".Rdata", sep = "")))
								s$species <- sp
								return(s)
							}
					) # eo lapply
					table <- do.call(rbind, scores)
					table$pool <- p
					rm(scores)
  					return(table)
 		   		}, mc.cores = 30
) # eo lapply
# Rbind
table_scores_phyto <- do.call(rbind,res)
table_scores_phyto$SDM <- factor(t(data.frame(str_split(as.character(rownames(table_scores_phyto)), pattern = "_", n = 2)))[,1])
table_scores_phyto$kingdom <- "Phytoplankton"
dim(table_scores_phyto); head(table_scores_phyto)
summary(table_scores_phyto)
rm(res)

### Join both tables and make plot
dim(table_scores_phyto) ; dim(table_scores_zoo)
colnames(table_scores_phyto); colnames(table_scores_zoo)
scores <- rbind(table_scores_phyto,table_scores_zoo)
summary(scores)
rm(table_scores_phyto,table_scores_zoo);gc()

### Provide classif to scores ddf and then plot per group2
unique(scores[scores$kingdom == "Zooplankton","species"]) # remove brackets
scores$species <- gsub("\\(|\\)", "", scores$species)
intersect(unique(scores$species),classif2$species)

# Provide group with for loop
scores$group <- NA
for(sp in unique(classif2$species)) {
        message(paste(sp, sep = ""))
        scores[scores$species == sp,"group"] <- unique(as.character(classif2[classif2$species == sp,"group2"]))
}
# Check
summary(factor(scores$group)) # Looks ok
head(scores)

plot <- ggplot(data = na.omit(scores), aes(x = factor(group), y = TSS)) + geom_violin(colour = "black", fill = "grey55") + 
        geom_boxplot(colour = "black", fill = "white", width = 0.2) + facet_grid(factor(pool) ~ factor(SDM)) + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Save plot
setwd(WD)
ggsave(plot = plot, filename = paste("plot_TSS_all_groups.jpg", sep = ""), dpi = 300, width = 14, height = 13)


# --------------------------------------------------------------------------------------------------------------------------------

### Go back to the raw occurrence data and plot sampling effort AND sampled species richness in space

setwd(paste(WD,"/biology/occurence_data_groups/v9/v9v8v5.1v3.1", sep = ""))
# Load all obs
require("parallel")
res <- mclapply(dir()[grep("30_04_19",dir())], function(f) {
			d <- get(load(f))
			return(d[,c("species","genus","family","order","class","phylum","month","year","xbin_1d","ybin_1d","source")])
	}, mc.cores = 30
) # eo lapply
ddf <- dplyr::bind_rows(res) 
dim(ddf);str(ddf)
rm(res);gc()
setwd(WD)
head(ddf)

### Compute and sampling effort 
# Add cell id
ddf$id <- factor(paste(ddf$xbin, ddf$ybin, sep = "_"))
ddf$x2 <- ddf$xbin_1d 
ddf[ddf$xbin_1d < 0 ,"x2"] <- (ddf[ddf$xbin_1d < 0 ,"xbin_1d"]) + 360

# Chck categories
unique(ddf$phylum)
unique(ddf$class)
unique(ddf$order)
unique(ddf$family)

### Provide groups2 like in classif2
ddf$group2 <- "Other_zoo"
# Jellyfish
ddf[ddf$phylum %in% c("Cnidaria","Ctenophora"),"group2"] <- "Jellyfish"
# Chaetognaths
ddf[ddf$phylum == "Chaetognatha","group2"] <- "Chaetognatha"
# Foraminifera
ddf[ddf$phylum == "Foraminifera","group2"] <- "Foraminifera"
# Salps & Doliolids
ddf[ddf$order %in% c("Salpida","Doliolida"),"group2"] <- "Salps"
# Appendicularia
ddf[ddf$class %in% c("Appendicularia","Larvacea"),"group2"] <- "Appendicularians"
# Pteropods
ddf[ddf$order %in% c("Thecosomata"),"group2"] <- "Pteropods"
# Amphipods
ddf[ddf$order %in% c("Amphipoda"),"group2"] <- "Amphipods"
# Euphausiids
ddf[ddf$order %in% c("Euphausiacea"),"group2"] <- "Euphausiids"
# And divide the copepods into 3 groups?
ddf[ddf$order %in% c("Calanoida"),"group2"] <- "Calanoids"
ddf[ddf$family %in% c("Oithonidae"),"group2"] <- "Oithonids"
ddf[ddf$family %in% c("Oncaeidae","Corycaeidae","Sapphirinidae"),"group2"] <- "Poecilostomatoids"
ddf$group2 <- factor(ddf$group2)
summary(ddf$group2)

### Tally observations and plot % occ per group2
Nocc <- nrow(ddf) # Nocc
groups.contrib <- data.frame(ddf %>% group_by(group2) %>% summarize(N = n(), Pocc = (n()/Nocc)*100))
groups.contrib[order(groups.contrib$Pocc, decreasing = T),]

ggplot(groups.contrib, aes(x = 2, y = Pocc, fill = factor(group2))) +
  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
  scale_fill_brewer(name = "Group", palette = "Paired") + theme_void() + xlim(0.5,2.5) 

### OK. Now plot sampling effort and richness sampled per 1°x1° cell bin
g <- "Oithonids"
for(g in unique(ddf$group2) ) {
    
    setwd(WD)
    subset <- ddf[ddf$group2 == g,]
    n <- nrow(subset)
    nsp <- length(unique(subset$species))
    message(paste("Mapping sampling efforts for ",g," (",nsp," species)", sep = ""))
    
    effort <- data.frame(subset %>% group_by(id) %>% summarize(x = unique(x2), y = unique(ybin_1d), n = n(), rich = length(unique(species))) ) # eo ddf
    # dim(effort) ; summary(effort)
    effort$logn <- log(effort$n)

    ### Make bins for nicer color palette
    effort$logn_bin <- factor(cut_interval(effort$logn,9))
    # levels(effort$logn_bin)
    levels <- str_replace_all(levels(effort$logn_bin), ",", "-")
    levels <- gsub("\\[|\\]", "", levels)
    levels <- gsub("\\(|\\)", "", levels)
    levels(effort$logn_bin) <- levels

    map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(logn_bin)), data = effort) + 
    			scale_fill_brewer(name = "Sampling effort\n(log)", palette = "YlOrRd") +
    			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "white", colour = "white", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
                    labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
                scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
	                labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "black"),legend.key = element_rect(fill = "black"),
	                panel.grid.major = element_line(colour = "black",linetype = "dashed"), legend.position = "right")

    ggsave(plot = map1, filename = paste("map_sampling_effort_",g,".jpg", sep = ""), dpi = 300, width = 7, height = 4)

    ### Plot sampled richness
    effort$rich_bin <- factor(cut_interval(effort$rich,7))
    # levels(effort$logn_bin)
    levels <- str_replace_all(levels(effort$rich_bin), ",", "-")
    levels <- gsub("\\[|\\]", "", levels)
    levels <- gsub("\\(|\\)", "", levels)
    levels(effort$rich_bin) <- levels
    
    map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(rich_bin)), data = effort) + 
    			scale_fill_brewer(name = "Richness\neffort", palette = "YlOrRd") +
    			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "white", colour = "white", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
                    labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
                scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
	                labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "black"),legend.key = element_rect(fill = "black"),
	                panel.grid.major = element_line(colour = "black",linetype = "dashed"), legend.position = "right")

    ggsave(plot = map2, filename = paste("map_richness_effort_",g,".jpg", sep = ""), dpi = 300, width = 7, height = 4)
        
}


### Same with phytoplankton
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/total_background/species_data")
files <- dir()[grep("data_total_",dir())] ; files
# f <- files[11]
res <- mclapply(files, function(f) {
			d <- read.table(f, h = T, sep = ";")
			return(d[,c("species","genus","family","order","class","phylum","month","av.year","cell_id","x","y")])
	}, mc.cores = 30
) # eo lapply
ddf <- dplyr::bind_rows(res) 
dim(ddf);str(ddf)
rm(res);gc()
head(ddf)

# Add cell id
ddf$id <- factor(paste(ddf$x, ddf$y, sep = "_"))
ddf$x2 <- ddf$x 
ddf[ddf$x < 0 ,"x2"] <- (ddf[ddf$x < 0 ,"x"]) + 360

unique(ddf$phylum)
unique(ddf$class)
unique(ddf$order)
unique(ddf$family)

ddf$group2 <- "Other_phyto"
ddf[ddf$class %in% c("Bacillariophyceae"),"group2"] <- "Diatoms"
ddf[ddf$class %in% c("Dinophyceae","Peridinea","Dinoflagellata incertae sedis"),"group2"] <- "Dinoflagellates"
ddf[ddf$class %in% c("Prymnesiophyceae","Coccolithophyceae","Haptophyta incertae sedis"),"group2"] <- "Coccolithophores"
ddf$group2 <- factor(ddf$group2)
summary(ddf$group2)

Nocc <- nrow(ddf) # Nocc
groups.contrib <- data.frame(ddf %>% group_by(group2) %>% summarize(N = n(), Pocc = (n()/Nocc)*100))
groups.contrib[order(groups.contrib$Pocc, decreasing = T),]

for(g in unique(ddf$group2) ) {
    
    setwd(WD)
    subset <- ddf[ddf$group2 == g,]
    n <- nrow(subset)
    nsp <- length(unique(subset$species))
    message(paste("Mapping sampling efforts for ",g," (",nsp," species)", sep = ""))
    
    effort <- data.frame(subset %>% group_by(id) %>% summarize(x = unique(x2), y = unique(y), n = n(), rich = length(unique(species))) ) # eo ddf
    # dim(effort) ; summary(effort)
    effort$logn <- log(effort$n)

    ### Make bins for nicer color palette
    effort$logn_bin <- factor(cut_interval(effort$logn,9))
    # levels(effort$logn_bin)
    levels <- str_replace_all(levels(effort$logn_bin), ",", "-")
    levels <- gsub("\\[|\\]", "", levels)
    levels <- gsub("\\(|\\)", "", levels)
    levels(effort$logn_bin) <- levels

    map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(logn_bin)), data = effort) + 
    			scale_fill_brewer(name = "Sampling effort\n(log)", palette = "YlOrRd") +
    			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "white", colour = "white", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
                    labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
                scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
	                labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "black"),legend.key = element_rect(fill = "black"),
	                panel.grid.major = element_line(colour = "black",linetype = "dashed"), legend.position = "right")

    ggsave(plot = map1, filename = paste("map_sampling_effort_",g,".jpg", sep = ""), dpi = 300, width = 7, height = 4)

    ### Plot sampled richness
    effort$rich_bin <- factor(cut_interval(effort$rich,7))
    # levels(effort$logn_bin)
    levels <- str_replace_all(levels(effort$rich_bin), ",", "-")
    levels <- gsub("\\[|\\]", "", levels)
    levels <- gsub("\\(|\\)", "", levels)
    levels(effort$rich_bin) <- levels
    
    map2 <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(rich_bin)), data = effort) + 
    			scale_fill_brewer(name = "Richness\neffort", palette = "YlOrRd") +
    			geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "white", colour = "white", size = 0.3) +
                coord_quickmap() + scale_x_continuous(name = "", breaks = c(0,60,120,180,240,300,360),
                    labels = c("GM","60°E","120°E","180°","120°W","60°W","GM"), expand = c(0,0)) +
                scale_y_continuous(name = "", breaks = c(-90,-60,-30,0,30,60,90), limits = c(-90,90), 
	                labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N"), expand = c(0,0)) +
                theme(panel.background = element_rect(fill = "black"),legend.key = element_rect(fill = "black"),
	                panel.grid.major = element_line(colour = "black",linetype = "dashed"), legend.position = "right")

    ggsave(plot = map2, filename = paste("map_richness_effort_",g,".jpg", sep = ""), dpi = 300, width = 7, height = 4)
        
}

### CONCLUSIONS:
# RF models are overfitting, sometimes ANN could be doing so too
# Re-run some tests with total background and group2 background when appropriate
# For groups with not enough data for target group-background, use sites (monthly 1x1 grid cells) that have at least 3 groups described in it?  
 

### Go to: /net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.1/total_background
### and print into:  /net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/no_background

setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.1/total_background")
files <- dir()[grep("data_total_",dir())] ; files
f <- files[5]
names <- colnames(read.table(f, h = T, sep = ";"))[c(1:15,17,18)] ; names

# f <- "data_total_Atrophia_glacialis.txt"

for(f in files) {
       
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/species_data_v9v3.1/total_background")
        message(paste(f, sep = ""))
        data <- read.table(f, h = T, sep = ";")
        subset <- data[data$obs == 1,names]
        
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/no_background")
        message(paste("Printing ",f, sep = ""))
        write.table(x = subset, file = str_replace_all(f,"_total","_occurrence"), sep = ";")
        rm(data, subset) ; gc()
               
} 


### Same with phytoplankton
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/total_background/species_data")
files <- dir()[grep("data_total_",dir())] ; files
f <- files[15]
names <- colnames(read.table(f, h = T, sep = ";"))[c(78,1:4,63,76,54:61)] ; names

for(f in files) {
       
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/total_background/species_data")
        message(paste(f, sep = ""))
        data <- read.table(f, h = T, sep = ";")
        subset <- data[data$obs == 1,names]
        
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/no_background")
        message(paste("Printing ",f, sep = ""))
        write.table(x = subset, file = str_replace_all(f,"_total","_occurrence"), sep = ";")
        rm(data, subset) ; gc()
               
} 


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


