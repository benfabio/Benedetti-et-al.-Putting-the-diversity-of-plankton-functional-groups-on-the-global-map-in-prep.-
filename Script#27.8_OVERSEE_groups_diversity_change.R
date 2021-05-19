
##### 15/10/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Loading the monthly projections of group level diversity (SR and % species modelled)
#   - Make sure baseline projections match previous ones (those in RSCRIPTBATCH_group_diversity_mappe.R)
#   - Compute and map future changes in diversity for each group
#   - Analyze similarity between groups' ensemble projections (PCA, covariations in Script#27.9)
#   - Quantify similarity between baseline richness projections based on reduced predictors set and full predictors set
#   - Quantify similarity between scenario where evrything changes versus the one where only SST changes
#   - Plot distribution of changes in SR per groups and between the 2 scenarios
#   - Also for the 'constant SST' scenario
 
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
world <- map_data("world")

SDMs <- c("GLM","GAM","ANN")
ESMs <- c("CNRM-PISCES","IPSL-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")		

# --------------------------------------------------------------------------------------------------------------------------------

### Load classif table if needed further down
strategies <- c("group")
categories <- c("Phytoplankton","Zooplankton")
combin <- apply(expand.grid(categories, strategies), 1, paste, collapse = "_")

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


### Load projections for both phyto and zoo
setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
base <- get(load("table_mon_rich_baseline_groups_22_10_20.Rdata"))
fut <- get(load("table_mon_rich_2100-2000_groups_09_11_20_constant_SST.Rdata"))
# dim(base) ; dim(fut)
# head(base) ; head(fut)
# str(base) ; str(fut)
# Check some maps
# subset <- base[which(base$month == "sep" & base$group == "Appendicularians" & base$SDM == "ANN"),]
# ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = subset) + scale_fill_viridis(name = "Richness") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#     scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# subset <- fut[which(fut$month == "apr" & fut$group == "Calanoida" & fut$SDM == "GLM" & fut$ESM == "CNRM-PISCES"),]
# dim(subset)
# ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = subset) + scale_fill_viridis(name = "Richness") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#         panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

### For each group, compute mean annual rich and perc per SDM/ESM
base.ens <- data.frame(base[base$rich > 0,] %>% group_by(cell_id, group, SDM) %>%
                summarize(x = unique(x), y = unique(y),
                rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T))
) # eo ddf
# head(base.ens) ; dim(base.ens) ; summary(base.ens)

#subset <- base.ens[which(base.ens$group == "Appendicularians" & base.ens$SDM == "ANN"),]
# dim(subset)
# ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = subset) + scale_fill_viridis(name = "Richness") +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "grey70", colour = "black", size = 0.3) +
#      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#          panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#      scale_x_continuous(name = "", limits = c(-180,180), expand = c(0,0), labels = NULL) +
#      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

fut.ens <- data.frame(fut[fut$rich > 0,] %>% group_by(cell_id, group, SDM, ESM) %>%
                summarize(x = unique(x), y = unique(y),
                rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T))
) # eo ddf
# summary(fut.ens)

# subset2 <- fut.ens[which(fut.ens$group == "Appendicularians" & fut.ens$SDM == "ANN" & fut.ens$ESM == "CNRM-PISCES"),]
# # dim(subset2)
# ggplot() + geom_raster(aes(x = x, y = y, fill = rich), data = subset2) + scale_fill_viridis(name = "Richness") +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#      coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#          panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#      scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


### 16/10/2020: For each group, compute differences in mean annual diversity for each SDM*ESM cmbination (n=15)
### Map and plot zonal gradient like in your first paper
groups <- unique(base.ens$group)

# require("parallel")
g <- "Diatoms"
sdm <- "GLM"
esm <- "GFDL-TOPAZ"

for(g in groups) {
    
        nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
        message(paste("Plotting changes in diversity for ",g," (",nsp," species)", sep = ""))
        base.sub <- base.ens[base.ens$group == g,]
        fut.sub <- fut.ens[fut.ens$group == g,]
        
        for(sdm in SDMs) {
            
            message(paste("based on ",sdm, sep = ""))
            base.sub2 <- base.sub[base.sub$SDM == sdm,]
            fut.sub2 <- fut.sub[fut.sub$SDM == sdm,]
           
            for(esm in ESMs) {
                
                message(paste("based on ",esm, sep = ""))
                base.sub3 <- base.sub2
                fut.sub3 <- fut.sub2[fut.sub2$ESM == esm,]
                # head(base.sub) ; head(fut.sub)
                # summary(base.sub$x2) ; summary(fut.sub$x)
                # summary(base.sub$y) ; summary(fut.sub$y)
                # Flip x coordinates for baseline projections
                base.sub3$x2 <- base.sub3$x
                base.sub3[base.sub3$x < 0 ,"x2"] <- (base.sub3[base.sub3$x < 0 ,"x"]) + 360
                # New id and order
                base.sub3$id <- factor(paste(base.sub3$x2,base.sub3$y, sep = "_"))
                base.sub3 <- base.sub3[order(base.sub3$id),] ; fut.sub3 <- fut.sub3[order(fut.sub3$cell_id),]
                # Find common cells
                commons <- intersect(unique(base.sub3$id), unique(fut.sub3$cell_id)) 
                # length(commons)
                base.sub3 <- base.sub3[base.sub3$id %in% commons,]
                fut.sub3 <- fut.sub3[fut.sub3$cell_id %in% commons,]
                # head(base.sub3) ; head(fut.sub3) ; dim(base.sub3) ; dim(fut.sub3)
                # Compute difference
                fut.sub3$rich_base <- base.sub3$rich
                fut.sub3$diff <- (fut.sub3$rich)-(fut.sub3$rich_base) # summary(fut.sub2)
                # Compute corresponding diff in number of species
                fut.sub3$diff_perc <- (fut.sub3$diff/fut.sub3$rich_base)*100
                # Compute zonal plot
                zonal <- data.frame(fut.sub3 %>% group_by(y) %>% summarize(avg_perc = mean(diff_perc, na.rm = T), sd = sd(diff_perc, na.rm=T)) )
                # summary(zonal) ; summary(fut.sub3)
                
                # Map and Zonal plot
                map <- ggplot() + geom_raster(aes(x = x, y = y, fill = diff_perc), data = fut.sub3[fut.sub3$diff_perc < 75,]) +
                    geom_raster(aes(x = x, y = y), data = fut.sub3[fut.sub3$diff_perc > 50,], fill = "#9e0142") + 
                    geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = diff_perc), data = fut.sub3[fut.sub3$diff_perc < 75,]) +
                    scale_fill_gradient2(name = paste("%SR change for\n",g, sep=""), mid = "white", low = "#2166ac", high = "#b2182b") +
                    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
                    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
	                    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
                    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
                    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

                z.plot <- ggplot() + geom_ribbon(aes(y = y, xmin = avg_perc-sd, xmax = avg_perc+sd), fill = "grey75", data = zonal) +
                     geom_path(aes(y = y, x = avg_perc), data = zonal, colour = "black", linetype = "solid") +
                     geom_vline(xintercept = 0, linetype = "dashed") +
                     scale_x_continuous(name = paste("Difference in ",g," SR (%)", sep="")) +
                     scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
                         breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
                     theme_classic()
      
                # Assemble in panel
                require("ggpubr")
                require("grDevices")
                panel <- ggarrange(map,z.plot, ncol = 2, nrow = 1, widths = c(2.8,1), align = "v")   
                ggsave(plot = panel, filename = paste("map_zonal_ann_diff_",g,"_",sdm,"_",esm,"_SSTonly.pdf", sep = ""), dpi = 300, width = 8, height = 3)    
                                
            } # eo for loop - esm in ESMs
            
        } # eo for loop - sdm in SDMs
    
        rm(panel,map,z.plot,zonal,base.sub3,fut.sub3,base.sub2,fut.sub2,nsp)
        gc()
    
} # eo for loop - g in groups
         
         
         
### And now ensemble difference (compute the 15 (3 SDMs x 5 ESMs) diff possible and average)
cobins <- expand.grid(SDMs,ESMs)
cobins <- paste(cobins$Var1, cobins$Var2, sep = "_")
g <- "Diatoms"

for(g in groups) {
    
        nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
        message(paste("Plotting changes in diversity for ",g," (",nsp," species)", sep = ""))
        base.sub <- base.ens[base.ens$group == g,]
        fut.sub <- fut.ens[fut.ens$group == g,]
        
        # i <- 5
        require("parallel")
        res <- mclapply(c(1:15), function(i) {
                    
                    message(paste("based on ",cobins[i], sep = ""))
                    sdm <- unlist(strsplit(cobins[i], split = "_"))[1]
                    esm <- unlist(strsplit(cobins[i], split = "_"))[2]
                    
                    message(paste("based on ",sdm, sep = ""))
                    base.sub2 <- base.sub[base.sub$SDM == sdm,]
                    fut.sub2 <- fut.sub[fut.sub$SDM == sdm & fut.sub$ESM == esm,]
                    base.sub2$x2 <- base.sub2$x
                    base.sub2[base.sub2$x < 0 ,"x2"] <- (base.sub2[base.sub2$x < 0 ,"x"]) + 360
                    # New id and order
                    base.sub2$id <- factor(paste(base.sub2$x2,base.sub2$y, sep = "_"))
                    base.sub2 <- base.sub2[order(base.sub2$id),] ; fut.sub2 <- fut.sub2[order(fut.sub2$cell_id),]
                    # Find common cells
                    commons <- intersect(unique(base.sub2$id), unique(fut.sub2$cell_id)) 
                    base.sub3 <- base.sub2[base.sub2$id %in% commons,]
                    fut.sub3 <- fut.sub2[fut.sub2$cell_id %in% commons,]
                    # Compute difference
                    fut.sub3$rich_base <- base.sub3$rich
                    fut.sub3$diff <- (fut.sub3$rich)-(fut.sub3$rich_base) # summary(fut.sub2)
                    # Compute corresponding diff in number of species
                    fut.sub3$diff_perc <- (fut.sub3$diff/fut.sub3$rich_base)*100
                    fut.sub3$model <- cobins[i]
                                                         
                    return(fut.sub3[,c("cell_id","x","y","diff","diff_perc","model")])
                    
                }, mc.cores = 15
            
        ) # eo mclapply
        # Rbind and compute model ensemble average
        ddf <- bind_rows(res)
        rm(res) ; gc()
        ens <- data.frame(ddf %>% group_by(cell_id) %>% summarize(x = unique(x), y = unique(y), diff = mean(diff, na.rm = T), perc = mean(diff_perc, na.rm = T)))
        #summary(ens)
        
        zonal <- data.frame(ens %>% group_by(y) %>% summarize(avg_perc = mean(perc, na.rm = T), sd = sd(perc, na.rm=T)) )              
                        
        # Map and Zonal plot
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = ens[ens$perc < 75,]) +
            geom_raster(aes(x = x, y = y), data = ens[ens$perc > 50,], fill = "#9e0142") + 
            geom_contour(colour = "grey60", binwidth = 25, size = 0.25, aes(x = x, y = y, z = perc), data = ens[ens$perc < 75,]) +
            scale_fill_gradient2(name = paste("%SR change for\n",g, sep=""), mid = "white", low = "#2166ac", high = "#b2182b") +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
            scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
            scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
        
        z.plot <- ggplot() + geom_ribbon(aes(y = y, xmin = avg_perc-sd, xmax = avg_perc+sd), fill = "grey75", data = zonal) +
             geom_path(aes(y = y, x = avg_perc), data = zonal, colour = "black", linetype = "solid") +
             geom_vline(xintercept = 0, linetype = "dashed") +
             scale_x_continuous(name = paste("Difference in ",g," SR (%)", sep="")) +
             scale_y_continuous(position = "right", name = "Latitude", limits = c(-90,90),
                 breaks = c(-90,-60,-30,0,30,60,90), labels = c("90°S","60°S","30°S","Eq","30°N","60°N","90°N")) +
             theme_classic()

        # Assemble in panel
        require("ggpubr")
        require("grDevices")
        panel <- ggarrange(map,z.plot, ncol = 2, nrow = 1, widths = c(2.8,1), align = "v")   
        
        ggsave(plot = panel, filename = paste("map_zonal_ann_diff_",g,"_ens_constant_SST.pdf", sep = ""), dpi = 300, width = 8, height = 3)
                                
} # eo for loop - g in groups


# --------------------------------------------------------------------------------------------------------------------------------


### 20/10/2020: Use base.ens to re-plot the baseline mean annual richness patterns for each SDM
for(g in groups) {
    
        nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
        message(paste("Plotting baseline SR for ",g," (",nsp," species)", sep = ""))
        base.sub <- base.ens[base.ens$group == g,]
      
        for(sdm in SDMs) {
            
            message(paste("based on ",sdm, sep = ""))
            base.sub2 <- base.sub[base.sub$SDM == sdm,]
                           
            base.sub3 <- base.sub2
            # Flip x coordinates for baseline projections
            base.sub3$x2 <- base.sub3$x
            base.sub3[base.sub3$x < 0 ,"x2"] <- (base.sub3[base.sub3$x < 0 ,"x"]) + 360
            # New id and order
            base.sub3$id <- factor(paste(base.sub3$x2,base.sub3$y, sep = "_"))
            base.sub3 <- base.sub3[order(base.sub3$id),]
            # head(base.sub3)
            # Compute zonal plot
            zonal <- data.frame(base.sub3 %>% group_by(y) %>% summarize(avg_perc = mean(perc, na.rm = T), sd = sd(perc, na.rm=T)) )
            # summary(zonal) 
            nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
                
            # Map and Zonal plot
            map <- ggplot() + geom_raster(aes(x = x2, y = y, fill = perc), data = base.sub3) +
                    geom_contour(colour = "grey60", binwidth = 0.1, size = 0.25, aes(x = x2, y = y, z = perc), data = base.sub3) +
                    scale_fill_viridis(name = paste("% spp modelled for\n",g," (",nsp,")", sep="") ) +
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
            require("ggpubr")
            require("grDevices")
            panel <- ggarrange(map,z.plot, ncol = 2, nrow = 1, widths = c(2.8,1), align = "v")   
            ggsave(plot = panel, filename = paste("map_zonal_ann_baseline_SR_",g,"_",sdm,".pdf", sep = ""), dpi = 300, width = 8, height = 3)
                                     
        } # eo for loop - sdm in SDMs
    
        rm(panel,map,z.plot,zonal,base.sub3,fut.sub3,base.sub2,fut.sub2,nsp)
        gc()
    
} # eo for loop - g in groups
         
         
         
### And now ensemble perc in SR modelled
for(g in groups) {
    
        nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
        message(paste("Mapping ensemble SR for ",g," (",nsp," species)", sep = ""))
        base.sub <- base.ens[base.ens$group == g,]        
        
        require("parallel")
        res <- mclapply(SDMs, function(sdm) {
                message(paste("based on ",sdm, sep = ""))
                base.sub2 <- base.sub[base.sub$SDM == sdm,]
                base.sub2$x2 <- base.sub2$x
                base.sub2[base.sub2$x < 0 ,"x2"] <- (base.sub2[base.sub2$x < 0 ,"x"]) + 360
                base.sub2$id <- factor(paste(base.sub2$x2,base.sub2$y, sep = "_"))
                base.sub2 <- base.sub2[order(base.sub2$id),]                                               
                return(base.sub2[,c("cell_id","x2","y","perc","rich","SDM")])
            }, mc.cores = 5
        ) # eo mclapply
        # Rbind and compute model ensemble average
        ddf <- bind_rows(res)
        rm(res) ; gc()
        ens <- data.frame(ddf %>% group_by(cell_id) %>% summarize(x = unique(x2), y = unique(y), rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T)))
        # summary(ens)
        
        zonal <- data.frame(ens %>% group_by(y) %>% summarize(avg_perc = mean(perc, na.rm = T), sd = sd(perc, na.rm=T)) )              
                        
        # Map and Zonal plot
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = perc), data = ens) +
                geom_contour(colour = "grey60", binwidth = 0.1, size = 0.25, aes(x = x, y = y, z = perc), data = ens) +
                scale_fill_viridis(name = paste("% spp modelled for\n",g," (",nsp,")", sep="") ) +
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
        require("ggpubr")
        require("grDevices")
        panel <- ggarrange(map,z.plot, ncol = 2, nrow = 1, widths = c(2.8,1), align = "v")   
        ggsave(plot = panel, filename = paste("map_zonal_ann_baseline_SR_",g,"_ens.pdf", sep = ""), dpi = 300, width = 8, height = 3)
                                
} # eo for loop - g in groups


# --------------------------------------------------------------------------------------------------------------------------------

### 22/10/2020: Evaluate similarities in %SR changes across groups (ensembles)
### - PCA (classic) + clustering on the %SR changes maybe? To regionalize the global patterns in changes
### - Distribution of %losses and %gains across groups based on ensmebles
### - Quantify variability in projections of ∆SR across groups and relate to sampling effort? (x = Noccurrences, y = Var)

# First, to implement the ideas above, you need to concatenate all ens proj of diff in a data.frame with groups' projections as columns
res <- lapply(groups, function(g) {
            
            message(paste("Compiling ensemble changes in diversity for ",g," (",nsp," species)", sep = ""))
            base.sub <- base.ens[base.ens$group == g,]
            fut.sub <- fut.ens[fut.ens$group == g,]
            nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
        
            require("parallel")
            res <- mclapply(c(1:15), function(i) {
                    
                        message(paste("based on ",cobins[i], sep = ""))
                        sdm <- unlist(strsplit(cobins[i], split = "_"))[1]
                        esm <- unlist(strsplit(cobins[i], split = "_"))[2]
                    
                        message(paste("based on ",sdm, sep = ""))
                        base.sub2 <- base.sub[base.sub$SDM == sdm,]
                        fut.sub2 <- fut.sub[fut.sub$SDM == sdm & fut.sub$ESM == esm,]
                        base.sub2$x2 <- base.sub2$x
                        base.sub2[base.sub2$x < 0 ,"x2"] <- (base.sub2[base.sub2$x < 0 ,"x"]) + 360
                        # New id and order
                        base.sub2$id <- factor(paste(base.sub2$x2,base.sub2$y, sep = "_"))
                        base.sub2 <- base.sub2[order(base.sub2$id),] ; fut.sub2 <- fut.sub2[order(fut.sub2$cell_id),]
                        # Find common cells
                        commons <- intersect(unique(base.sub2$id), unique(fut.sub2$cell_id)) 
                        base.sub3 <- base.sub2[base.sub2$id %in% commons,]
                        fut.sub3 <- fut.sub2[fut.sub2$cell_id %in% commons,]
                        # Compute difference
                        fut.sub3$rich_base <- base.sub3$rich
                        fut.sub3$diff <- (fut.sub3$rich)-(fut.sub3$rich_base) # summary(fut.sub2)
                        # Compute corresponding diff in number of species
                        fut.sub3$diff_perc <- (fut.sub3$diff/fut.sub3$rich_base)*100
                        fut.sub3$model <- cobins[i]
                                                         
                        return(fut.sub3[,c("cell_id","x","y","diff","diff_perc","model")])
                    
                    }, mc.cores = 15
            
            ) # eo mclapply
            # Rbind and compute model ensemble average
            ddf <- bind_rows(res)
            rm(res) ; gc()
            ens <- data.frame(ddf %>% group_by(cell_id) %>%
                    summarize(x = unique(x), y = unique(y),
                    diff = mean(diff, na.rm = T), perc = mean(diff_perc, na.rm = T)) )
            #
            ens$group <- g
        
            return(ens)
    
        } # eo fun 

) # eo lapply 
# Rbind
table.diff <- bind_rows(res)    
rm(res) ; gc()
head(table.diff) ; dim(table.diff)
# Dcast to have groups' perc as columns
d.table.diff <- dcast(data = table.diff[,c(1:3,5,6)], formula = cell_id + x + y ~ group, value.var = "perc")
dim(d.table.diff)
summary(d.table.diff)
# Nice.


### A°) PCA to assess similarities in change sin diversity
require("FactoMineR") ; require("vegan")
data4PCA <- na.omit(d.table.diff)
res.pca <- PCA(data4PCA[,c(4:length(data4PCA))], scale.unit = T, graph = F)
summary(res.pca)

eig <- data.frame(perc = res.pca$eig[,"percentage of variance"], nb = c(1:nrow(res.pca$eig)) ) # mean(eig$perc)
pca1 <- paste0("PC1 (",floor(eig$perc[1]*100)/100,"%)")
pca2 <- paste0("PC2 (",floor(eig$perc[2]*100)/100,"%)")
pca3 <- paste0("PC3 (",floor(eig$perc[3]*100)/100,"%)")
pca4 <- paste0("PC4 (",floor(eig$perc[4]*100)/100,"%)")

# For better PCA plot
library("ggrepel")
library("broom")
augment.PCA <- function(x, dims = c(1:4), which="col") {
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
pcad <- augment.PCA(res.pca)
pcad$type # 5 - 7 --> phyto
pcad[c(1:4,8:14),"type"] <- "Zooplankton"
pcad[c(5:7),"type"] <- "Phytoplankton"

ggplot(pcad) +
    coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
    annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
    geom_segment(aes(x=0, xend = Dim.1, y=0, yend = Dim.2, colour = factor(type)), arrow=arrow(angle=20, length=unit(0.01, "npc"))) +
    scale_colour_manual(name = "", values = c("#4d9221","#c51b7d"), guide = F) + 
    geom_text_repel(aes(x=Dim.1, y=Dim.2, label=var), 
            data=filter(pcad, (Dim.1^2+Dim.2^2) > 0.2^2), segment.alpha=0.5) +
    xlab(pca1) + ylab(pca2) + theme_bw()

ggplot(pcad) +
    coord_fixed() + scale_x_continuous(breaks=0) + scale_y_continuous(breaks=0) +
    annotate(geom="path", colour="black", x=cos(seq(0, 2*pi, length.out=100)), y=sin(seq(0, 2*pi, length.out=100)), colour = "grey30") +
    geom_segment(aes(x=0, xend = Dim.2, y=0, yend = Dim.3, colour = factor(type)), arrow=arrow(angle=20, length=unit(0.01,"npc"))) +
    scale_colour_manual(name = "", values = c("#4d9221","#c51b7d"), guide = F) + 
    geom_text_repel(aes(x = Dim.2, y = Dim.3, label = var), 
            data=filter(pcad, (Dim.2^2+Dim.3^2) > 0.2^2), segment.alpha=0.5) +
    xlab(pca2) + ylab(pca3) + theme_bw()

### Retrieve coords along PC1-3 and map
data4PCA[,paste("PC",c(1:3),sep="")] <- res.pca$ind$coord[,c(1:3)]
summary(data4PCA)
### Map PC1
ggplot() + geom_raster(aes(x = x, y = y, fill = PC1), data = data4PCA) +
    geom_contour(colour = "grey60", binwidth = 3, size = 0.25, aes(x = x, y = y, z = PC1), data = data4PCA) +
    scale_fill_gradient2(name = "PC1", mid = "white", low = "#2166ac", high = "#b2182b") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
# Map PC2
ggplot() + geom_raster(aes(x = x, y = y, fill = PC2), data = data4PCA) +
    geom_contour(colour = "grey60", binwidth = 3, size = 0.25, aes(x = x, y = y, z = PC2), data = data4PCA) +
    scale_fill_gradient2(name = "PC2", mid = "white", low = "#2166ac", high = "#b2182b") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
# Map PC3
ggplot() + geom_raster(aes(x = x, y = y, fill = PC3), data = data4PCA) +
    geom_contour(colour = "grey60", binwidth = 3, size = 0.25, aes(x = x, y = y, z = PC3), data = data4PCA) +
    scale_fill_gradient2(name = "PC3", mid = "white", low = "#2166ac", high = "#b2182b") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
### OK. As before.


### Now, try to represent the differences in % changes across groups? From 'table.diff'
ggplot(data = table.diff, aes(x = factor(group), y = perc)) +
    geom_boxplot(aes(fill = factor(group)), colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    xlab("") + ylab("% change in diversity") + theme_classic()


### Same as above, but instead of showing mean projection, retrieve its variations across models
cobins <- expand.grid(SDMs,ESMs)
cobins <- paste(cobins$Var1, cobins$Var2, sep = "_")

res <- lapply(groups, function(g) {
            
            base.sub <- base.ens[base.ens$group == g,]
            fut.sub <- fut.ens[fut.ens$group == g,]
            nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
            message(paste("Compiling ensemble changes in diversity for ",g," (",nsp," species)", sep = ""))
        
            require("parallel")
            res <- mclapply(c(1:15), function(i) {
                    
                        message(paste("based on ",cobins[i], sep = ""))
                        sdm <- unlist(strsplit(cobins[i], split = "_"))[1]
                        esm <- unlist(strsplit(cobins[i], split = "_"))[2]
                    
                        message(paste("based on ",sdm, sep = ""))
                        base.sub2 <- base.sub[base.sub$SDM == sdm,]
                        fut.sub2 <- fut.sub[fut.sub$SDM == sdm & fut.sub$ESM == esm,]
                        base.sub2$x2 <- base.sub2$x
                        base.sub2[base.sub2$x < 0 ,"x2"] <- (base.sub2[base.sub2$x < 0 ,"x"]) + 360
                        # New id and order
                        base.sub2$id <- factor(paste(base.sub2$x2,base.sub2$y, sep = "_"))
                        base.sub2 <- base.sub2[order(base.sub2$id),] ; fut.sub2 <- fut.sub2[order(fut.sub2$cell_id),]
                        # Find common cells
                        commons <- intersect(unique(base.sub2$id), unique(fut.sub2$cell_id)) 
                        base.sub3 <- base.sub2[base.sub2$id %in% commons,]
                        fut.sub3 <- fut.sub2[fut.sub2$cell_id %in% commons,]
                        # Compute difference
                        fut.sub3$rich_base <- base.sub3$rich
                        fut.sub3$diff <- (fut.sub3$rich)-(fut.sub3$rich_base) # summary(fut.sub2)
                        # Compute corresponding diff in number of species
                        fut.sub3$diff_perc <- (fut.sub3$diff/fut.sub3$rich_base)*100
                        fut.sub3$model <- cobins[i]
                                                         
                        return(fut.sub3[,c("cell_id","x","y","diff","diff_perc","model")])
                    
                    }, mc.cores = 15
            
            ) # eo mclapply
            # Rbind and compute model ensemble average
            ddf <- bind_rows(res)
            rm(res) ; gc()
            ens <- data.frame(ddf %>% group_by(cell_id) %>%
                    summarize(x = unique(x), y = unique(y),
                        perc = mean(diff_perc, na.rm = T),
                        sd = sd(diff_perc, na.rm = T),
                        var = var(diff_perc, na.rm = T) ) 
            )
            #
            ens$group <- g
        
            return(ens)
    
        } # eo fun 

) # eo lapply 
# Rbind
table.var <- bind_rows(res)    
rm(res) ; gc()
head(table.var) ; dim(table.var)
summary(table.var)

### First, simply compute mean variance in projs per group
vars.group <- data.frame(table.var %>% group_by(group) %>% summarize(mean.sd = mean(sd), mean.var = mean(var)))
# vars.group[order(vars.group$mean.var, decreasing = T),]
vars.group[order(vars.group$mean.sd, decreasing = T),]
# Nice, really shows how the least best sampled groups are more variable in projetcions ! 

### And now, plot and map variances in projections per groups
for(g in unique(table.var$group)) {
    
        message(paste("Mapping variance in projections for ", g, sep = "")) 
        # g <- "Salps"
        sub <- table.var[table.var$group == g,]
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = log(sd)), data = sub) +
            scale_fill_viridis(name = paste("Projections stdev\n for ",g,sep=""), option = "B", limits = c(0,8)) +
            geom_contour(colour = "grey60", binwidth = 2, size = 0.25, aes(x = x, y = y, z = log(sd)), data = sub) +
            geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
            coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
            scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
            scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
                
        ggsave(plot = map, filename = paste("map_variance_ens_2100-2000_",g,".pdf", sep = ""), dpi = 300, width = 6, height = 3)
    
} # eo for loop - g in groups

### Very nce as well, now combine in a ddf:
# - 'vars.group'
# - 'nb of occurrences' (from table.classif)
# - 'nb species modelled'
vars.group$Nsp <- data.frame(table.classif %>% group_by(group) %>% summarize(nsp = length(unique(species))))$nsp
# And sampling effort?

### For zoo
setwd(paste(WD,"/biology/occurence_data_groups/v9/v9v8v5.1v3.1", sep = ""))
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

### Same with phytoplankton
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/phytoplankton_15_01_19/total_background/species_data")
files <- dir()[grep("data_total_",dir())] ; files
# f <- files[11]
res <- mclapply(files, function(f) {
			d <- read.table(f, h = T, sep = ";")
			return(d[,c("species","genus","family","order","class","phylum","month","av.year","cell_id","x","y")])
	}, mc.cores = 30
) # eo lapply
ddf2 <- dplyr::bind_rows(res) 
dim(ddf2);str(ddf2)
rm(res);gc()
head(ddf2)
ddf2$group2 <- "Other_phyto"
ddf2[ddf2$class %in% c("Bacillariophyceae"),"group2"] <- "Diatoms"
ddf2[ddf2$class %in% c("Dinophyceae","Peridinea","Dinoflagellata incertae sedis"),"group2"] <- "Dinoflagellates"
ddf2[ddf2$class %in% c("Prymnesiophyceae","Coccolithophyceae","Haptophyta incertae sedis"),"group2"] <- "Coccolithophores"
ddf2$group2 <- factor(ddf2$group2)
summary(ddf2$group2)

# Rbind 
names <- c("species","genus","family","order","group2")
ddf3 <- rbind(ddf[,names], ddf2[,names])
#dim(ddf3)
# And derive effort from 'ddf3'
effort <- data.frame(ddf3 %>% group_by(group2) %>% summarize(n = n()))
# Change some levels in groups name of 'effort'
levels(effort$group2)[levels(effort$group2) == "Calanoids"] <- "Calanoida"
levels(effort$group2)[levels(effort$group2) == "Oithonids"] <- "Oithonida"
levels(effort$group2)[levels(effort$group2) == "Poecilostomatoids"] <- "Poecilostomatoida"

vars.group$Nocc <- NA
for(g in unique(vars.group$group)) {
    vars.group[vars.group$group == g,"Nocc"] <- unique(effort[effort$group2 == g,"n"])
} # eo for loop
vars.group

### Finally, plot mean.sd or mean.var as a function of Nocc and change points size by Nsp
summary(lm(mean.sd ~ log(Nocc), data = vars.group)) # Adjusted R-squared: 0.4254; p-value: 0.0068
summary(lm(mean.var ~ log(Nocc), data = vars.group)) # Adjusted R-squared: 0.2406; p-value: 0.043
#summary(lm(log(mean.var) ~ log(Nocc), data = vars.group)) # Adjusted R-squared: 0.2406; p-value: 0.043

require("ggrepel")
qplot <- ggplot(data = vars.group) +
    geom_point(aes(x = log(Nocc), y = mean.sd, size = Nsp), pch = 21, colour = "black", fill = "grey55") + 
    geom_text_repel(aes(x = log(Nocc), y = mean.sd, label = group), segment.alpha = 0.5) +
    xlab("Sampling effort - log(#occurrences)") + ylab("Uncertainty in ensemble projections\n(average stdev)") + 
    scale_size_continuous(name = "Number of species\nmodelled") +
    theme_classic()
#
setwd(WD)
ggsave(plot = qplot, filename = "plot_uncertainty_vs_effort_groups.jpg", dpi = 300, width = 7, height = 5)


# --------------------------------------------------------------------------------------------------------------------------------

### 26/10/2020: Quantify similarity between baseline richness projections based on reduced predictors set and full predictors set
setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
rich.reduced <- get(load("table_ann_rich_env_baseline_23_10_20.Rdata"))
# dim(rich.reduced); summary(rich.reduced)

### And now get the same table but for the full set
setwd(paste(WD,"/biology/data_for_group_studies/zooplankton/group_background/communities", sep = ""))
zoo <- read.table("table_mon_composition_Zooplankton_group.txt", sep = "\t", h = T)
# head(zoo) ; dim(zoo)
setwd(paste(WD,"/biology/data_for_group_studies/phytoplankton/group_background/communities", sep = ""))
phy <- read.table("table_mon_composition_Phytoplankton_group.txt", sep = "\t", h = T)
# head(phy) ; dim(phy)

### From phy, zoo and table.classif, get ensemble baseline richness of eahc groups
groups <- unique(table.classif$group) ; groups
# g <- "Calanoida"
res <- mclapply(groups, function(g) {
            
            if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
                comm <- phy
            } else {
                comm <- zoo
            } # eo if else loop - comm
            
            nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
            species2keep <- unique(table.classif[table.classif$group == g,"species"])
            if( length(grep(pattern = '-', x = species2keep, fixed = T)) > 0 ) {
                 species2keep <- str_replace_all(species2keep,'-','')
            } # eo if loop
            message(paste("Compiling ensemble changes in diversity for ",g," (",nsp," species)", sep = ""))
            
            # Derive mean ensemble richness 
            cols2keep <- intersect(species2keep, colnames(comm))
            subset <- comm[,c("cell_id","x","y",cols2keep,"SDM","month")]
            subset$rich <- rowSums(subset[,cols2keep], na.rm = T)
            grp.div <- subset[,c("cell_id","x","y","month","SDM","rich")]
            grp.div$perc <- (grp.div$rich)/nsp
            grp.div$group <- g
            # Flip x coordinates for Pacific-centered map
            grp.div$x2 <- grp.div$x
            grp.div[grp.div$x < 0 ,"x2"] <- (grp.div[grp.div$x < 0 ,"x"]) + 360
            # New id and order
            grp.div$id <- factor(paste(grp.div$x2, grp.div$y, sep = "_"))
            grp.div <- grp.div[order(grp.div$id),]
            # Compute mean annual baseline richness
            ens <- data.frame(grp.div[grp.div$rich > 0,] %>% group_by(id) %>% summarize(x = unique(x2), y = unique(y),
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
# Nice
rich.full <- ddf.rich
rm(ddf.rich, zoo, phy, comm, ens) ; gc()

# Dcast rich.full
d.rich.full <- dcast(data = rich.full[,c(1:4,6)], formula = id + x + y ~ group, value.var = "rich")
head(rich.full) ; dim(rich.full)

# Quick map to check if 'rich.full' was computed properly
ggplot() + geom_raster(aes(x = x, y = y, fill = Pteropods), data = rich.full) +
    scale_fill_viridis(name = paste("Annual richness\n of ","Pteropods", sep = ""), option = "B") +
    geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
    coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
    panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") + 
    scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) + 
    scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


### Now, for each group, combine baselie richness pattern from reduced set to full set and quantify covariance with lm()
g <- "Diatoms"
for(g in groups) {
    
        message(paste("Quantifying similarity in baseline richness for ",g, sep = ""))
        message(paste("", sep = ""))
        sub.red <- rich.reduced[,c("id","x","y",g)]
        sub.full <- rich.full[,c("id","x","y",g)]
        # head(sub.red) ; head(sub.full)
        # Find common cells and cbind
        commons <- intersect(unique(sub.red$id), unique(sub.full$id)) # length(commons)
        sub.red2 <- sub.red[sub.red$id %in% commons,]
        sub.full2 <- sub.full[sub.full$id %in% commons,]               
        # head(sub.red2) ; head(sub.full2)
        # summary(sub.red2) ; summary(sub.full2)
        colnames(sub.red2) <- c("id","x","y","reduced")
        colnames(sub.full2) <- c("id","x","y","full")
        
        sub.red2 <- sub.red2[order(sub.red2$id),]
        sub.full2 <- sub.full2[order(sub.full2$id),]
        
        sub.full2$reduced <- sub.red2$reduced     
        
        lm <- lm(reduced ~ full, data = sub.full2) #; summary(lm) ; str(summary(lm))
        radj <- round(summary(lm)$adj.r.squared,5)
        plot <- ggplot(aes(x = full, y = reduced), data = sub.full2) + 
            geom_point(colour = "grey55", alpha = .5) + geom_smooth(method = "lm", colour = "#d53e4f") + 
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
            xlab("Species richness (all predictors)") + ylab("Species richness (reduced predictors)") + theme_classic() +
            ggtitle(paste(g," (r2=",radj,")", sep = ""))
            
        # Save
        setwd(WD)
        ggsave(plot = plot, filename = paste("plot_lm_ann_rich_ens_fullvsreduced_",g,".jpg", sep = ""), dpi = 300, width = 4, height = 3)
        
} # eo g in groups

### 26/10/2020: After discussing with Luke G., decide to share the ensemble annual SR patterns for the various groups (table prepared below)
### and ave them match with annual conditions

### First, dcast 'rich.full'
d.rich.full <- dcast(data = rich.full[,c(1:4,6)], formula = id + x + y ~ group, value.var = "rich")
head(d.rich.full) ; dim(d.rich.full)

d.perc.full <- dcast(data = rich.full[,c(1:3,5,6)], formula = id + x + y ~ group, value.var = "perc")
head(d.perc.full) ; dim(d.perc.full)

# ggplot() + geom_raster(aes(x = x, y = y, fill = Diatoms), data = d.perc.full) +
#     scale_fill_viridis(name = paste("Annual richness\n of ","Diatoms", sep = ""), option = "B") +
#     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
#     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

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

### Combine with table.perc using common cells
commons <- intersect(unique(d.perc.full$id), unique(ann.clim$id))
# length(commons)
ann.clim2 <- ann.clim[ann.clim$id %in% commons,]
d.perc.full2 <- d.perc.full[d.perc.full$id %in% commons,]
# Re-order
ann.clim2 <- ann.clim2[order(ann.clim2$id),]
d.perc.full2 <- d.perc.full2[order(d.perc.full2$id),]
dim(ann.clim2) ; dim(d.perc.full2)
# cbind
d.perc.full3 <- cbind(d.perc.full2, ann.clim2[,c(4:length(ann.clim2))])
head(d.perc.full3)

ggplot() + geom_raster(aes(x = x, y = y, fill = Coccolithophores), data = d.perc.full3) +
     scale_fill_viridis(name = 'Coccolithophores', option = "B") +
     geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
     coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "left") +
     scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
     scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)


### Save new data.frame because it take ssome time to compute
setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
write.table(d.perc.full3, file = "table_ann_perc_groups+env_baseline_27_10_20.txt", sep = "\t")
#rm(ann.clim2,table.perc2,ddf.grp.div,fut,base) ; gc()


# --------------------------------------------------------------------------------------------------------------------------------

### 03/11/2020: Compare (e.g. linear regressions) projections based on all predictors vs. the scenario where only SST is projected to change
setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
base <- get(load("table_mon_rich_baseline_groups_22_10_20.Rdata"))
fut1 <- get(load("table_mon_rich_2100-2000_groups_03_11_20_SSTonly.Rdata"))
fut2 <- get(load("table_mon_rich_2100-2000_groups_09_11_20_constant_SST.Rdata"))

### For each group, compute mean annual rich and perc per SDM/ESM
base.ens <- data.frame(base[base$rich > 0,] %>% group_by(cell_id, group, SDM) %>%
                summarize(x = unique(x), y = unique(y),
                rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T))
) # eo ddf

fut1.ens <- data.frame(fut1[fut1$rich > 0,] %>% group_by(cell_id, group, SDM, ESM) %>%
                summarize(x = unique(x), y = unique(y),
                rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T))
) # eo ddf

fut2.ens <- data.frame(fut2[fut2$rich > 0,] %>% group_by(cell_id, group, SDM, ESM) %>%
                summarize(x = unique(x), y = unique(y),
                rich = mean(rich, na.rm = T), perc = mean(perc, na.rm = T))
) # eo ddf
# dim(base.ens) ; dim(fut1.ens) ; dim(fut2.ens)


### For both scenarios (SST change sonly vs. all changes) and for every group, compute % change in SR and return ensemble perc/diff for both scanerios
cobins <- expand.grid(SDMs,ESMs)
cobins <- paste(cobins$Var1, cobins$Var2, sep = "_")
groups <- unique(base.ens$group)
# g <- "Coccolithophores"

require("parallel")
res <- mclapply(groups, function(g) {
            
            nsp <- length(unique(table.classif[table.classif$group == g,"species"]))
            message(paste("Compiling ensemble changes in diversity for ",g," (",nsp," species)", sep = ""))
            base.sub <- base.ens[base.ens$group == g,]
            fut.sub1 <- fut1.ens[fut1.ens$group == g,]
            fut.sub2 <- fut2.ens[fut2.ens$group == g,]
        
            res <- lapply(c(1:15), function(i) {
                    
                        message(paste("based on ",cobins[i], sep = ""))
                        sdm <- unlist(strsplit(cobins[i], split = "_"))[1]
                        esm <- unlist(strsplit(cobins[i], split = "_"))[2]
                        
                        ### Compute ensemble
                        base.sub2 <- base.sub[base.sub$SDM == sdm,]
                        base.sub2$x2 <- base.sub2$x
                        base.sub2[base.sub2$x < 0 ,"x2"] <- (base.sub2[base.sub2$x < 0 ,"x"]) + 360
                        base.sub2$id <- factor(paste(base.sub2$x2,base.sub2$y, sep = "_"))
                        base.sub2 <- base.sub2[order(base.sub2$id),]
                        
                        # Find common cells
                        fut.sub2.1 <- fut.sub1[fut.sub1$SDM == sdm & fut.sub1$ESM == esm,]
                        fut.sub2.1 <- fut.sub2.1[order(fut.sub2.1$cell_id),]
                        fut.sub2.2 <- fut.sub2[fut.sub2$SDM == sdm & fut.sub2$ESM == esm,]
                        fut.sub2.2 <- fut.sub2.2[order(fut.sub2.2$cell_id),]
                        
                        commons1 <- intersect(unique(base.sub2$id), unique(fut.sub2.1$cell_id)) 
                        commons2 <- intersect(unique(base.sub2$id), unique(fut.sub2.2$cell_id)) 
                        
                        base.sub3.1 <- base.sub2[base.sub2$id %in% commons1,]
                        fut.sub3.1 <- fut.sub2.1[fut.sub2.1$cell_id %in% commons1,]
                        base.sub3.2 <- base.sub2[base.sub2$id %in% commons2,]
                        fut.sub3.2 <- fut.sub2.2[fut.sub2.2$cell_id %in% commons2,]
                        
                        # Compute difference
                        fut.sub3.1$rich_base <- base.sub3.1$rich
                        fut.sub3.2$rich_base <- base.sub3.2$rich
                        fut.sub3.1$diff <- (fut.sub3.1$rich)-(fut.sub3.1$rich_base) 
                        fut.sub3.2$diff <- (fut.sub3.2$rich)-(fut.sub3.2$rich_base) 
                        
                        # Compute corresponding diff in number of species
                        fut.sub3.1$diff_perc <- (fut.sub3.1$diff/fut.sub3.1$rich_base)*100
                        fut.sub3.2$diff_perc <- (fut.sub3.2$diff/fut.sub3.2$rich_base)*100
                        
                        fut.sub3.1$model <- cobins[i]
                        fut.sub3.2$model <- cobins[i]
                        
                        fut.sub3.1$scenario <- "SST only"
                        fut.sub3.2$scenario <- "Constant SST"
                        
                        ### Rbind and return
                        ddf <- rbind(fut.sub3.1,fut.sub3.2)                                 
                        return(ddf[,c("cell_id","x","y","diff","diff_perc","model","scenario")])
                    
                    } #, mc.cores = 15
            
            ) # eo lapply
            
            # Rbind and compute model ensemble average
            ddf <- bind_rows(res)
            rm(res) ; gc()
            
            ens <- data.frame(ddf %>% group_by(cell_id,scenario) %>%
                    summarize(x = unique(x), y = unique(y),
                    diff = mean(diff, na.rm = T), perc = mean(diff_perc, na.rm = T))
            ) # eo data.frame

            ens$group <- g
        
            return(ens)
    
        }, mc.cores = length(groups)

) # eo lapply 
# Rbind
table.diff <- bind_rows(res)    
rm(res) ; gc()
head(table.diff) ; dim(table.diff)

# Save this table
save(table.diff, file = "table_ann_diff_ens_constantSST.vs.SSTonly_groups.Rdata")
table.diff <- get(load("table_ann_diff_ens_constantSST.vs.SSTonly_groups.Rdata"))

# Dcast to have groups' perc as columns
d.table.diff <- dcast(data = table.diff[,c(1:4,6:7)], formula = cell_id + x + y ~ group + scenario, value.var = "perc")
dim(d.table.diff)
head(d.table.diff)
summary(d.table.diff)
# Nice.

### For each group, select corresponding columns, remove NAs and quantify similarty with lm's r2 , cor and plot
g <- "Poecilostomatoida"
for(g in groups) {
    
        message(paste(g, sep = ""))
        names <- colnames(d.table.diff)[grep(g,colnames(d.table.diff))]
        sub <- na.omit(d.table.diff[,c("cell_id","x","y",names)])
        
        lm <- lm(sub[,names[1]] ~ sub[,names[2]], data = sub) 
        # summary(lm)$adj.r.squared
        radj <- round(summary(lm)$adj.r.squared, 3)
        cor.spear <- round(cor(sub[,names[1]], sub[,names[2]], method = "spearman"), 3)

        plot <- ggplot(data = sub, aes(x = get(names[2]), y = get(names[1]), colour = y)) +
            scale_colour_distiller(name = "Latitude", palette = "Spectral") + 
            geom_point(alpha = .5) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            geom_hline(yintercept = 0, linetype = "dotted") + geom_vline(xintercept = 0, linetype = "dotted") + 
            ylab(paste("% difference in SR based on 'constant SST' scenario ","(",g,")", sep = "")) + xlab("% difference in SR based on 'SST only' scenario") +
            theme_classic() + ggtitle(paste("r2 = ",radj," ; rho = ",cor.spear, sep = ""))
            
        # ggsave
        ggsave(plot = plot, filename = paste("plot_SSTonly_vs_All_ann_perc_ens_",g,".pdf", sep = ""), dpi = 300, width = 6, height = 5)
    
} # eo g in groups

### From 'table.diff', you should also compare the distribution of the 2 scenarios by faceting per groups ;-)
head(table.diff)

plot <- ggplot(data = table.diff[which(table.diff$perc < 200),], aes(x = factor(scenario), y = perc, fill = factor(scenario))) +
    geom_violin(colour = "black") + geom_boxplot(fill = "white", colour = "black", width = 0.1) + 
    scale_fill_manual(name = "Scenario", values = c("#3B9AB2","#F21A00")) + xlab("") + ylab("% difference in SR") + theme_classic() +
    facet_wrap(~ factor(group), scales = "free")
#
ggsave(plot = plot, filename = paste("plots_distrib_SSTonly_vs_All_ann_perc_ens_groups.pdf", sep = ""), dpi = 300, width = 9, height = 11)



### 04/11/2020: CONCEPT: separet cells according to their changes in SR (per groups) across the 2 scenarios and illustrate them 
### following the concept of "valleys" and "sides".
# - All > 0 & SST only > 0 and above 1:1 line --> valley of positive agreement, amplification side
# - All > 0 & SST only > 0 and below 1:1 line --> valley of positive agreement, attenuation side
# - All < 0 & SST only < 0 and above 1:1 line --> valley of negative agreement, attenuation side
# - All < 0 & SST only < 0 and below 1:1 line --> valley of negative agreement, amplification side
# - All > 0 & SST only < 0 --> valley of rescue !
# - All < 0 & SST only > 0 --> valley of peril !

setwd(paste(WD,"/biology/data_for_group_studies", sep = ""))
table.diff <- get(load("table_ann_diff_ens_constantSST.vs.SSTonly_groups.Rdata"))
d.table.diff <- dcast(data = table.diff[,c(1:4,6:7)], formula = cell_id + x + y ~ group + scenario, value.var = "perc")
groups <- unique(table.diff$group)

### Redraw plot biplot above and map regions based on their valleys
 g <- "Diatoms"

for(g in groups) {
    
        message(paste(g, sep = ""))
        names <- colnames(d.table.diff)[grep(g,colnames(d.table.diff))]
        sub <- na.omit(d.table.diff[,c("cell_id","x","y",names)])
        
        # NOTE:
        # - names#1 == constant SST
        # - names#2 == SST only
        
        # Define valleys and sides
        sub$valley <- NA
        sub[which(sub[,names[1]] > 0 & sub[,names[2]] > 0),"valley"] <- "Agreement"
        sub[which(sub[,names[1]] < 0 & sub[,names[2]] < 0),"valley"] <- "Agreement"
        sub[which(sub[,names[1]] < 0 & sub[,names[2]] > 0),"valley"] <- "Negative_overide"
        sub[which(sub[,names[1]] > 0 & sub[,names[2]] < 0),"valley"] <- "Positive_overide"
        # summary(factor(sub$valley))
        
        # Add whether the cell is above or below the 1:1 line
        sub$sign <- NA
        sub[which(sub$valley == "Negative_overide"),"sign"] <- " "
        sub[which(sub$valley == "Positive_overide"),"sign"] <- " "
        
        sub[which( sub[,names[1]] > 0 & sub[,names[2]] > 0 & sub[,names[1]] > sub[,names[2]] ),"sign"] <- "Positive_amplification"
        sub[which( sub[,names[1]] > 0 & sub[,names[2]] > 0 & sub[,names[1]] < sub[,names[2]] ),"sign"] <- "Positive_attenuation"
        
        sub[which( sub[,names[1]] < 0 & sub[,names[2]] < 0 & sub[,names[1]] > sub[,names[2]] ),"sign"] <- "Negative_attenuation"
        sub[which( sub[,names[1]] < 0 & sub[,names[2]] < 0 & sub[,names[1]] < sub[,names[2]] ),"sign"] <- "Negative_amplification"
        # summary(factor(sub$sign))
        
        # And derive combination of both valley and sign
        sub$category <- paste(sub$valley, sub$sign, sep = "_")
        # summary(factor(sub$category))
        # And correct slightly
        sub[which(sub$category == "Negative_overide_ "),"category"] <- "Rescuing losses not due to SST"
        sub[which(sub$category == "Positive_overide_ "),"category"] <- "Rescuing losses due to SST"
        
        sub[which(sub$category == "Agreement_Negative_amplification"),"category"] <- "Amplifying losses"
        sub[which(sub$category == "Agreement_Negative_attenuation"),"category"] <- "Attenuating losses"
        sub[which(sub$category == "Agreement_Positive_amplification"),"category"] <- "Amplifying gains"
        sub[which(sub$category == "Agreement_Positive_attenuation"),"category"] <- "Attenuating gains"
        
        cor.spear <- round(cor(sub[,names[1]], sub[,names[2]], method = "spearman"), 3)
        
        # Define color palette to choose for unique(sub$category)
        cols <- c("Rescuing losses not due to SST" = "#e31a1c", "Rescuing losses due to SST" = "#33a02c",
                "Amplifying losses" = "#1f78b4", "Attenuating losses" = "#a6cee3",
                "Amplifying gains" = "#ff7f00", "Attenuating gains" = "#fdbf6f")

        plot <- ggplot(data = sub[which(sub[,names[1]] < 250),], aes(x = get(names[2]), y = get(names[1]), colour = factor(category))) +
            scale_colour_manual(name = "", values = cols) + 
            geom_point() + geom_hline(yintercept = 0, linetype = "dotted") +
            geom_vline(xintercept = 0, linetype = "dotted") + 
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            ylab(paste("% difference in SR based on\n'constant SST' scenario", sep = "")) +
            xlab("% difference in SR based on\n'SST only' scenario") +
            theme_classic() + ggtitle(paste(g," (Spearman's rho = ",cor.spear,")", sep = "")) + guides(colour = F)
    
        # And map
        map <- ggplot() + geom_raster(aes(x = x, y = y, fill = factor(category)), data = sub) +
                 scale_fill_manual(name = "", values = cols) + 
                 geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "black", colour = "black", size = 0.3) +
                 coord_quickmap() + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
                 panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
                 scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
                 scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)
        #
        # ggsave
        #ggsave(plot = map, filename = paste("map_SSTonly_vs_All_ann_perc_ens_valleys_",g,".pdf", sep = ""), dpi = 300, width = 4, height = 7)
        require("ggpubr")
        panel <- ggarrange(plot, map, ncol = 2, nrow = 1, align = "hv", widths = c(1.5,2)) #; panel
        
        ggsave(plot = panel, filename = paste("panel_constantSST_vs_SSTonly_ann_perc_ens_valleys_",g,".pdf", sep = ""), dpi = 300, width = 10, height = 4)
    
} # eo g in groups


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
