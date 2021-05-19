
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

### 1°) Get a table associating the species to their groups (needed for computing group diversity)
### Simplest way is to get the TSS tables

# c <- "Phytoplankton_total"
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
table <- bind_rows(res)
rm(res) ; gc()
# dim(table) ; head(table)

# ### You may use this table to test score variations too
# summary(table[table$strat == "total","TSS"])
# summary(table[table$strat == "group","TSS"])
# summary(table[table$cat == "Phytoplankton","TSS"])
# summary(table[table$cat == "Zooplankton","TSS"])

# ggplot(data = table, aes(x = factor(cat), y = TSS, fill = factor(cat))) + geom_violin(colour = "black") +
#     geom_boxplot(colour = "black", fill = "white", width = 0.2) + xlab("") + ylab("TSS") + theme_bw()
# # Clear diff ^^

# ggplot(data = table, aes(x = factor(strat), y = TSS, fill = factor(strat))) + geom_violin(colour = "black") +
#     geom_boxplot(colour = "black", fill = "white", width = 0.2) + xlab("") + ylab("TSS") + theme_bw()
# # No diff

# ggplot(data = table, aes(x = factor(strat), y = TSS, fill = factor(strat))) + geom_violin(colour = "black") +
#     geom_boxplot(colour = "black", fill = "white", width = 0.2) + xlab("") + ylab("TSS") + theme_bw() +
#     facet_grid(.~ factor(cat))
# INTERESTING.
# kruskal.test(TSS ~ factor(strat), data = na.omit(table[table$cat == "Phytoplankton",]))
# chi-squared = 99.215, df = 1, p-value < 2.2e-16

### Anyways, you may now use this table to derive group-level diversity patterns
groups <- unique(table$group) # groups

### 2°) Extarct communities and map group-level diversity patterns for both background strategies
# c <- "Zooplankton_group"
groups.mapper <- function(c = combin) {
    
    message(paste("", sep = ""))
    message(paste("Mapping annual diversity for ",c, sep = ""))
    message(paste("", sep = ""))
    # Extract strat and cat from c
    cat <- do.call(cbind,strsplit(c, "_"))[1,1]
    strat <- do.call(cbind,strsplit(c, "_"))[2,1]
    
    # Go to the communities wd
    if(cat == "Phytoplankton" & strat == "total") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/communities")
        data.wd.back <- getwd()
    } else if(cat == "Zooplankton" & strat == "total") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/communities")
        data.wd.back <- getwd()
    } else if(cat == "Phytoplankton" & strat == "group") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities")
        data.wd.back <- getwd()
    } else if(cat == "Zooplankton" & strat == "group") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities")
        data.wd.back <- getwd()
    } # eo else if loop 
    
    comm <- read.table(paste("table_mon_composition_",cat,"_",strat,".txt", sep = ""), h = T, sep = "\t")
    
    # Define groups to map according to category
    if(cat == "Phytoplankton") {
        groups2 <- groups[c(1:3)]
    } else if(cat == "Zooplankton") {
        groups2 <- groups[c(4:14)]
    } # eo else if loop 
    
    for(g in groups2) {
        
        # g <- "Diatoms" ; g <- "Appendicularians"
    	message(paste(" ", sep = ""))	
    	message(paste("MAPPING ",g," DIVERSITY PATTERN BASED ON ",strat," BACKGROUND ", sep = ""))
    	message(paste(" ", sep = ""))
        
        species2keep <- unique(table[table$group == g,"species"])
        nsp <- length(species2keep)
        # Filter species n comm
        if( length(grep(pattern = '-', x = species2keep, fixed = T)) > 0 ) {
             species2keep <- str_replace_all(species2keep,'-','')   
        }
        cols2keep <- intersect(species2keep, colnames(comm))
        # Flip x coordinates for baseline projections
        comm$x2 <- comm$x
        comm[comm$x < 0 ,"x2"] <- (comm[comm$x < 0 ,"x"]) + 360
        comm2 <- comm[,c("cell_id","x2","y",cols2keep,"SDM","month")]
        # Compute rowSums 
        comm2$rich <- rowSums(as.matrix(comm2[,cols2keep]), na.rm = T)
        comm2$perc <- (comm2$rich)/nsp
        # summary(comm2)
        # sdm <- "GAM"
        for(sdm in unique(comm2$SDM)) {
            
            # Compute ensemble
            comm3 <- comm2[comm2$SDM == sdm,]
            ens2 <- data.frame(comm3[comm3$rich > 0,] %>% group_by(cell_id) %>% summarize(x = unique(x2), y = unique(y), perc = mean(perc, na.rm = T))) 
            # Derive zonal pattern and longitudinal variation (sd)   
            zonal <- data.frame(ens2 %>% group_by(y) %>% summarize(avg_perc = mean(perc, na.rm = T), sd = sd(perc, na.rm=T)) )
        
            map <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = ens2) +
                  scale_fill_viridis(name = paste("% spp modelled for\n",g," (",nsp,")", sep="") ) +
                  geom_contour(colour = "grey60", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = perc), data = ens2) +
                  geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                  coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    	            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
                  scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
                  scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
              
            z.plot <- ggplot() + geom_ribbon(aes(y = y, xmin = avg_perc-sd, xmax = avg_perc+sd), fill = "grey75", data = zonal) +
                geom_path(aes(y = y, x = avg_perc), data = zonal, colour = "black", linetype = "solid") +
                scale_x_continuous(name = paste(g," SR (%spp)", sep="")) +
                scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
                breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
                theme_classic()
      
            # Assemble in panel
            require("ggpubr"); require("grDevices")
            panel <- ggarrange(map,z.plot, ncol = 2, nrow = 1, widths = c(2.8,1), align = "v")
            setwd("/net/kryo/work/fabioben/OVERSEE/data/")
            ggsave(plot = panel, filename = paste("map_zonal_ann_baseline_SR_",g,"_",sdm,".pdf", sep = ""), dpi = 300, width = 8, height = 3)     
            
        } # eo for loop - sdm in SDMs
                
    } # eo for loop - g in groups2
    
} # eo FUN - groups.mapper

### Apply in parallel
require("parallel")
mclapply(X = combin, FUN = groups.mapper, mc.cores = 3)


# c <- combin[1]

all.plankton.mapper <- function(c = combin) {
    
    message(paste("Mapping annual diversity for ",c, sep = ""))
    message(paste("", sep = ""))
    # Extract strat and cat from c
    cat <- do.call(cbind,strsplit(c, "_"))[1,1]
    strat <- do.call(cbind,strsplit(c, "_"))[2,1]
    
    # Go to the communities wd
    if(cat == "Phytoplankton" & strat == "total") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/communities")
        data.wd.back <- getwd()
    } else if(cat == "Zooplankton" & strat == "total") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/communities")
        data.wd.back <- getwd()
    } else if(cat == "Phytoplankton" & strat == "group") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities")
        data.wd.back <- getwd()
    } else if(cat == "Zooplankton" & strat == "group") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities")
        data.wd.back <- getwd()
    } # eo else if loop 
    
    comm <- read.table(paste("table_mon_composition_",cat,"_",strat,".txt", sep = ""), h = T, sep = "\t")
    comm$x2 <- comm$x
    comm[comm$x < 0 ,"x2"] <- (comm[comm$x < 0 ,"x"]) + 360
   
    index <- (length(comm))-3
    nsp <- length(colnames(comm)[c(4:index)])

    # Compute rowSums and H' index
    comm$rich <- rowSums(as.matrix(comm[,c(4:index)]), na.rm = T)
    comm$perc <- (comm$rich)/nsp
        
    # # Compute ensemble
  #   ens <- data.frame(comm %>% group_by(cell_id,SDM) %>%
  #           summarize( x = unique(x), y = unique(y), H = mean(H, na.rm = T), perc = mean(perc, na.rm = T) )
  #   ) # eo ens
    # summary(ens)
    # map each month separately
    #months <- unique(ens$month)
    #SDMs <- unique(ens$SDM)
    # sdm <- "GAM"
        
    # for(sdm in SDMs) {
#
#             data2map <- ens[ens$SDM == sdm,] # dim(data2map) ; summary(data2map)
#
#             data2map[which(data2map$perc == 0 & is.na(data2map$H)),"perc"] <- NA
#
#             map1 <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = data2map) +
#                    scale_fill_viridis(name = paste("Richness\n(%spp modelled)\nn = ",nsp, sep = "") ) +
#                    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#                    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#                      panel.grid.major = element_line(colour = "grey70",linetype = "dashed") ) +
#                   scale_x_continuous(name = "Longitude", limits = c(-180,180), expand = c(0,0)) +
#                   scale_y_continuous(name = "Latitude", limits = c(-85,85), expand = c(0,0))
#
#             # Save
#             setwd("/net/kryo/work/fabioben/OVERSEE/data/")
#             ggsave(plot = map1, filename = paste("map_ann_SR_",cat,"_",strat,"_",sdm,".jpg", sep = ""), width = 7, height = 5, dpi = 300)
#             #ggsave(plot = map2, filename = paste("map_mon_H_",cat,"_",strat,"_",m,".jpg", sep = ""), width = 7, height = 5, dpi = 300)
#             setwd(data.wd.back)
#
#     } # eo for loop - sdm in SDMs
    
        # Compute ensemble
        ens2 <- data.frame(comm[comm$rich > 0,] %>% group_by(cell_id) %>% summarize(x = unique(x2), y = unique(y), perc = mean(perc, na.rm = T))) 
        # Derive zonal pattern and longitudinal variation (sd)   
        zonal <- data.frame(ens2 %>% group_by(y) %>% summarize(avg_perc = mean(perc, na.rm = T), sd = sd(perc, na.rm=T)) )
    
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = ens2) +
              scale_fill_viridis(name = paste("% spp modelled for\n",cat," (",nsp,")", sep="") ) +
              geom_contour(colour = "grey60", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = perc), data = ens2) +
              geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
              coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
	            panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
              scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
              scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
          
        z.plot <- ggplot() + geom_ribbon(aes(y = y, xmin = avg_perc-sd, xmax = avg_perc+sd), fill = "grey75", data = zonal) +
            geom_path(aes(y = y, x = avg_perc), data = zonal, colour = "black", linetype = "solid") +
            scale_x_continuous(name = paste(cat," SR (%spp)", sep="")) +
            scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
            breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
            theme_classic()
  
        # Assemble in panel
        require("ggpubr"); require("grDevices")
        panel <- ggarrange(map,z.plot, ncol = 2, nrow = 1, widths = c(2.8,1), align = "v")
        setwd("/net/kryo/work/fabioben/OVERSEE/data/")
        ggsave(plot = panel, filename = paste("map_zonal_ann_baseline_SR_",cat,"_ens.pdf", sep = ""), dpi = 300, width = 8, height = 3)     
    
} # eo FUN - groups.mapper

### Apply in parallel
mclapply(X = combin, FUN = all.plankton.mapper, mc.cores = 3)


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
