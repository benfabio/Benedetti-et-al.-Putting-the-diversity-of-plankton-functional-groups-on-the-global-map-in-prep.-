
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

strat <- "group"
categories <- c("phytoplankton","zooplankton")

### 1°) Get the table associating the species to their groups (needed for computing group diversity)
res <- mclapply(categories, function(cat) {
            
            # Go to the communities wd
            if(cat == "phytoplankton" & strat == "total") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/Scores")
                data.wd.back <- getwd()
            } else if(cat == "zooplankton" & strat == "total") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/Scores")
                data.wd.back <- getwd()
            } else if(cat == "phytoplankton" & strat == "group") {
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/Scores")
                data.wd.back <- getwd()
            } else if(cat == "zooplankton" & strat == "group") {
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
    
    }, mc.cores = 2
    
) # eo mclapply 
table.classif <- bind_rows(res)
rm(res) ; gc()

table.classif$species <- str_replace_all(table.classif$species, "-", "")

# Use this table to derive group-level diversity patterns
groups <- unique(table.classif$group) # groups

### 2°) Load baseline monthly projections and derive mean annual SR for each group
# cat <- "phytoplankton"
projs <- mclapply(categories, function(cat) {

              message(paste(" ", sep = ""))
              message(paste("CONCATENATING ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))

              # Go to the communities wd
              if(cat == "phytoplankton" & strat == "total") {
                      setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/communities/climate_change_projections")
                      data.wd.back <- getwd()
              } else if(cat == "zooplankton" & strat == "total") {
                      setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/communities/climate_change_projections")
                      data.wd.back <- getwd()
              } else if(cat == "phytoplankton" & strat == "group") {
                      setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities/climate_change_projections")
                      data.wd.back <- getwd()
              } else if(cat == "zooplankton" & strat == "group") {
                      setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities/climate_change_projections")
                      data.wd.back <- getwd()
              } # eo else if loop

              comm <- read.table(paste("table_mon_composition_baseline_all_",cat,"_",strat,".txt", sep = ""), sep = "\t", h = T)

              # Define groups to map according to category
              if(cat == "phytoplankton") {
                  groups2 <- groups[c(1:3)]
              } else if(cat == "zooplankton") {
                  groups2 <- groups[c(4:14)]
              } # eo else if loop

              # In a lapply, subset species per group and compute diversity after melting and summarizing
              # g <- "Calanoida"
              res <- lapply(groups2, function(g) {

                           message(paste("Computing SR (sum of HSI) estimates for ",g, sep = ""))
                           species2keep <- unique(table.classif[table.classif$group == g,"species"])
                           nsp <- length(species2keep)
                           # Filter species n comm
                           if( length(grep(pattern = '-', x = species2keep, fixed = T)) > 0 ) {
                                species2keep <- str_replace_all(species2keep,'-','')
                           } # eo if loop
                           
                           cols2keep <- intersect(species2keep, colnames(comm))
                           
                           subset <- comm[,c("cell_id","x","y",cols2keep,"SDM","month")]
                           # Compute diversity
                           # dim(subset) ; head(subset) ; summary(subset[subset$SDM == "GLM","rich"])
                           subset$rich <- rowSums(subset[,cols2keep], na.rm = T)
                           # Drop col2keeps
                           grp.div <- subset[,c("cell_id","x","y","month","SDM","rich")]
                           # m.subset <- melt(subset, id.vars = c("cell_id","x","y","SDM","month"))
                           # colnames(m.subset)[c(6,7)] <- c("species","HSI")
                           # grp.div <- data.frame(m.subset %>% group_by(cell_id,month,SDM) %>% summarize(rich = sum(HSI,na.rm=T)))
                           # rm(subset,m.subset) ; gc()
                           # head(grp.div) ; summary(grp.div)
                           # Derive % species modelled
                           grp.div$perc <- (grp.div$rich)/nsp
                           grp.div$group <- g
                           grp.div$category <- cat
                           # summary(grp.div[grp.div$month == "apr",'rich']) ; summary(grp.div[grp.div$month == "dec",'rich'])
                           return(grp.div)

                       } # eo FUN

              ) # eo lapply
              # Rbind
              div <- dplyr::bind_rows(res)
              # head(div) ; summary(div)
              rm(res,comm) ; gc()

              # Return
              return(div)

      }, mc.cores = 2

) # eo mclapply
# Rbind
ddf.grp.div <- bind_rows(projs)
rm(projs) ; gc()

### Save it somewhere because the script above is long enough
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/")
save(ddf.grp.div, file = "table_mon_rich_baseline_groups_22_10_20.Rdata")
# logout and re-load it in a new session

rm(ddf.grp.div) ; gc()

### ================================================================================================================

### 14/10/2020: Same as above but for future projections
 cat <- "zooplankton"

projs <- mclapply(categories, function(cat) {
            
            message(paste(" ", sep = ""))	
            message(paste("CONCATENATING ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
  
            # Go to the communities wd
            if(cat == "phytoplankton" & strat == "total") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/communities/climate_change_projections")
                    data.wd.back <- getwd()
            } else if(cat == "zooplankton" & strat == "total") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/communities/climate_change_projections")
                    data.wd.back <- getwd()
            } else if(cat == "phytoplankton" & strat == "group") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities/climate_change_projections")
                    data.wd.back <- getwd()
            } else if(cat == "zooplankton" & strat == "group") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities/climate_change_projections")
                    data.wd.back <- getwd()
            } # eo else if loop 
            
            comm <- read.table(paste("table_mon_composition_2100-2000_all_",cat,"_",strat,"_constant_SST.txt", sep = ""), sep = "\t", h = T)
           
            # Define groups to map according to category
            if(cat == "phytoplankton") {
                groups2 <- groups[c(1:3)]
            } else if(cat == "zooplankton") {
                groups2 <- groups[c(4:14)]
            } # eo else if loop 
            
            # In a lapply, subset species per group and compute diversity after melting and summarizing
            # g <- "Poecilostomatoida"
            res <- lapply(groups2, function(g) {
                
                         message(paste("Computing SR (sum of HSI) estimates for ",g, sep = ""))
                         species2keep <- unique(table.classif[table.classif$group == g,"species"])
                         nsp <- length(species2keep)
                         if( length(grep(pattern = '-', x = species2keep, fixed = T)) > 0 ) {
                              species2keep <- str_replace_all(species2keep,'-','')   
                         } # eo if loop
                         
                         cols2keep <- intersect(species2keep, colnames(comm))
                         subset <- comm[,c("cell_id","x","y",cols2keep,"SDM","ESM","month")]
                         
                         # Compute diversity
                         # head(subset[subset$SDM == "GAM",]) ; summary(subset[subset$SDM == "GLM","rich"])
                         subset$rich <- rowSums(subset[,cols2keep], na.rm = T) 
                         # Drop col2keeps 
                         grp.div <- subset[,c("cell_id","x","y","month","SDM","ESM","rich")]
                         grp.div$perc <- (grp.div$rich)/nsp
                         grp.div$group <- g
                         grp.div$category <- cat
                         # dim(grp.div) ; head(grp.div) ; summary(grp.div)
                         # summary(grp.div[grp.div$SDM == "GAM",])
                         
                         return(grp.div)
                         
                     } # eo FUN
                     
            ) # eo lapply 
            # Rbind
            div <- dplyr::bind_rows(res)   
            
            rm(res,comm) ; gc()           
            # head(div) ; summary(div)
            
            # Return                     
            return(div)
    
    }, mc.cores = 2 
    
) # eo mclapply 
# Rbind
ddf.grp.div <- bind_rows(projs)
rm(projs) ; gc()

### Save it somewhere because the script above is long enough
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/")
save(ddf.grp.div, file = "table_mon_rich_2100-2000_groups_09_11_20_constant_SST.Rdata")


### ================================================================================================================
### ================================================================================================================
### ================================================================================================================
