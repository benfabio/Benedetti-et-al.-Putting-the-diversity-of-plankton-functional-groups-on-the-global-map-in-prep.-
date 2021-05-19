
### ================================================================================================================

library("raster")
library("sp")
library("dplyr")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")
library("viridis")
library("scales")
library("maps")

WD <- getwd()

### ================================================================================================================

firstup <- function(x) { substr(x,1,1) <- toupper(substr(x,1,1)); x } # eo firstup fun

# Vector of SDMs
SDMs <- c('GAM','GLM','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 
# Vector of months
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Set the working directories, vectors etc.

strategies <- c("group")
categories <- c("Phytoplankton","Zooplankton")

cat <- "Phytoplankton"
strat <- "group"

for(cat in categories) {
        
    	message(paste(" ", sep = ""))	
    	message(paste("EXTRACTING ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
    	message(paste(" ", sep = ""))
        
        # Set the working directories
        if(cat == "Phytoplankton" & strat == "total") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
            data.wd.back <- getwd()
        } else if(cat == "Zooplankton" & strat == "total") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
            data.wd.back <- getwd()
        } else if(cat == "Phytoplankton" & strat == "group") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background")
            data.wd.back <- getwd()
        } else if(cat == "Zooplankton" & strat == "group") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
            data.wd.back <- getwd()
        } 
        
        # Retrieve scores
        setwd(paste(data.wd.back,"/Scores/", sep = ""))
        files <- dir()[grep("scores_for_projections",dir())] #  files
		scores <- lapply(files, function(f) { s <- get(load(f)) ; return(s) } ) # eo lapply
        table <- do.call(rbind,scores)
        table$SDM <- factor(t(data.frame(str_split(as.character(rownames(table)), pattern = "_", n = 2)))[,1])
		rm(scores) ; gc()
        
        scores <- data.frame(table %>% group_by(species) %>% summarise(avg_TSS = mean(TSS, na.rm = T), sd_TSS = sd(TSS, na.rm = T)) ) # eo ddf
        # summary(scores)
        species2keep <- unique(scores[scores$avg_TSS > 0.3,"species"])
        
        # # While you're at it, plot distrbution per group and SDM
 #        plot <- ggplot(data = table, aes(x = factor(group), y = TSS)) + geom_boxplot(colour = "black", fill = "grey65") +
 #                facet_grid(. ~ factor(SDM)) + xlab("") + ylab("True Skill Statistics") +
 #                theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
     
        # setwd("/net/kryo/work/fabioben/OVERSEE/data/")
 #        if(cat == "Zooplankton") {
 #            ggsave(plot = plot, filename = paste("plot_scores_group_",strat,"_",cat,".jpg", sep = ""), dpi = 300, width = 8, height = 5)
 #        } else if(cat == "Phytoplankton") {
 #            ggsave(plot = plot, filename = paste("plot_scores_group_",strat,"_",cat,".jpg", sep = ""), dpi = 300, width = 6, height = 4)
 #        }
 
        setwd(data.wd.back)
        
        # m <- "mar"
        community.extracter <- function(m = months) {
    
                    message(paste(" ", sep = ""))
                    message(paste("Retrieving probabilities for ", m, sep = ""))
                    message(paste(" ", sep = ""))
                    
                    # Load env variables
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
                    # glob_stack_month_apr_18_09_20.txt
                    env <- read.table(paste("glob_stack_month_",m,"_18_09_20.txt", sep = ""), h = T, sep = ";")
                    env <- env[-which(env$SSS < 20),]
                    env <- env[-which(env$Bathy > -175),]
                    #env$x2 <- env$x
                    #env[env$x < 0 ,"x2"] <- (env[env$x < 0 ,"x"]) + 360
                    #env$id <- paste(env$x2, env$y, sep = "_")
            
                    setwd(data.wd.back)
            
                    # And for each SDM
                    # sdm <- "GLM"
                    for(sdm in SDMs) {
                    
                            message(paste("SDM == ",sdm, sep = ""))
                            # sp <- str_replace_all(species2keep, "_", ".")[401]
                            
                            probas <- lapply(X = str_replace_all(species2keep, "_", "."), FUN = function(sp) {

                                    # Got to species dir
                                    if( grepl(pattern = '-', x = sp, fixed = T) ) {
                                         sp <- str_replace_all(sp,'-','')   
                                    }
                                    
                                    message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
                                    setwd(paste(data.wd.back,"/niche.modelling/projections/",sp, sep = ""))
                                    
                                    # # UNSURE: Need to modify sp when there are 2 names and add brackets?
 #                                    if( length(strsplit(sp, ".", fixed = T)[[1]]) == 3 ) {
 #                                        sp <- paste(strsplit(sp, ".", fixed = T)[[1]][1],".(",
 #                                            strsplit(sp, ".", fixed = T)[[1]][2],").",
 #                                            strsplit(sp, ".", fixed = T)[[1]][3], sep = "" )
 #                                    } # eo if loop

                                    # If the 4 seasonal projections are done
                                    if( sum(grepl("proj_projection_", dir())) == 12 ) {
                                        
                                        # Load projections for each SDM
                                        setwd( paste(paste("proj_projection_",sp,"_",m, sep = ""),"/", sep = "") )
                                        d <- get(load( paste("proj_projection_",sp,"_",m,"_",sp, ".RData", sep = "") ))
                                        resModel <- d[,sdm,,]
                                        resModel <- apply(resModel, 1, mean, na.rm = F)
                                        resModel <- (resModel/1000)
                                        # Return
                                        return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y,
                                                species = gsub("\\.","_",sp), HSI = resModel) )
                                                
                                    } else {
                                        
                                        message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

                                    } # eo if else loop

                                } # eo FUN
                
                            ) # eo lapply
                            
                            tbl <- dplyr::bind_rows(probas)
                            # dim(tbl) ; colnames(tbl)
                            rm(probas); gc()
                            
                            # Dcast to put species as columns 
                            d <- dcast(tbl, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
                            
                            ### Add richness and shannon-weaver index
                            #require("vegan")
                            #d$rich <- rowSums(as.matrix( d[,c(4:length(d))]) )
                            #d$H <- diversity(d[,c(4:length(d))], "shannon")
                            
                            ### Save !
                            message(paste("Saving community table for ",sdm," | ",m, sep = ""))
                            message(paste(" ", sep = ""))
                            setwd(paste(data.wd.back,"/","communities","/", sep = ""))
                            write.table(d, file = paste("table_mon_composition+div_baseline_",cat,"_",strat,"_",sdm,"_",m,".txt", sep = ""), sep = "\t")
                    
                        } # eo for loop - sdm in SDMs
    
        } # eo FUN - community.extracter

        # Apply fun in mclapply
        require("parallel")
        mclapply(X = months, FUN = community.extracter, mc.cores = 15)

} # eo cat in categories


# --------------------------------------------------------------------------------------------------------------------------------

strategies <- c("group")
categories <- c("phytoplankton","zooplankton")
combin <- apply(expand.grid(categories, strategies), 1, paste, collapse = "_")

### Concatenate all community tables, and derive diversity patterns per group
# c <- "Zooplankton_group"
for(c in combin) {

        # Extract strat and cat from c
        cat <- do.call(cbind,strsplit(c, "_"))[1,1]
        strat <- do.call(cbind,strsplit(c, "_"))[2,1]
        
        message(paste(" ", sep = ""))	
	    message(paste("CONCATENATING ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
        
        # Go to the communities wd
        if(cat == "phytoplankton" & strat == "total") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/communities")
            data.wd.back <- getwd()
        } else if(cat == "zooplankton" & strat == "total") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/communities")
            data.wd.back <- getwd()
        } else if(cat == "phytoplankton" & strat == "group") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities")
            data.wd.back <- getwd()
        } else if(cat == "zooplankton" & strat == "group") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities")
            data.wd.back <- getwd()
        } # eo else if loop 
        
        files <- dir()[grep("div_baseline_",dir())]
        # f <- files[3]
        require("parallel")
        comms <- mclapply(files, function(f) {
                    d <- read.table(f, h = T, sep = "\t")
                    terms <- do.call(cbind,strsplit(f,"_"))
                    d$SDM <- terms[7,1]
                    d$month <- str_replace_all(terms[8,1],".txt","")
                    return(d)
            } , mc.cores = 30
        ) # eo lapply - f in files   
        # Rbind
        table <- bind_rows(comms)
        rm(comms) ; gc()
        # dim(table)
        
        # Save table
	    message(paste("SAVING ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
	    message(paste(" ", sep = ""))
        write.table(x = table, paste("table_mon_composition_",cat,"_",strat,".txt", sep = ""), sep = "\t")
    
} # eo for loop - c in combin


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------

