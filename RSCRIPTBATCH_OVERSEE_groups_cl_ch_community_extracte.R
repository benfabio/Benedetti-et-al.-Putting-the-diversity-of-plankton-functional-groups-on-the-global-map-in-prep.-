
# --------------------------------------------------------------------------------------------------------------------------------

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

# --------------------------------------------------------------------------------------------------------------------------------

firstup <- function(x) {
  substr(x,1,1) <- toupper(substr(x,1,1))
  x
} # eo firstup fun

# Vector of SDMs etc.
SDMs <- c('GAM','GLM','ANN')
# eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 
months <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
ESMs <- c("CESM-BEC","CNRM-PISCES","GFDL-TOPAZ","IPSL-PISCES","MRI-NEMURO")
rcp <- "rcp85"

# --------------------------------------------------------------------------------------------------------------------------------

# m <- "apr"
# sdm <- "GAM"
 
# baseline.community.extracter <- function(m = months) {
#
#              message(paste(" ", sep = ""))
#              message(paste("Retrieving baseline zoo probabilities for ", m, sep = ""))
#              message(paste(" ", sep = ""))
#              # Load env variables
#              setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
#              env <- read.table(paste("glob_stack_month_",m,"_18_09_20.txt", sep = ""), h = T, sep = ";")
#              env <- env[-which(env$SSS < 20),]
#              env <- env[-which(env$Bathy > -175),]
#
#              setwd(paste(data.wd.back,"/niche.modelling/future_projections/", sep = ""))
#              spp <- dir()
#              spp <- gsub("\\(|\\)", "", spp)
#
#              # And for each SDM
#              for(sdm in SDMs) {
#
#                      message(paste("SDM == ",sdm, sep = ""))
#                      require("parallel")
#                      # sp <- spp[5]
#
#                      probas <- lapply(X = spp, FUN = function(sp) {
#
#                                  # Got to species dir
#                                  setwd(paste(data.wd.back,"/niche.modelling/future_projections/",sp,"/", sep = ""))
#                                  message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))
#
#                                  # If the 4 seasonal projections are done
#                                  if( sum(grepl("proj_projection_", dir())) == 72 ) {
#
#                                          # Load projections for each SDM
#                                          if(sp %in% c("Coscinodiscus.oculusiridis","Pseudonitzschia.delicatissima","Pseudonitzschia.lineola",
#                                                      "Pseudonitzschia.pungens","Pseudonitzschia.seriata","Pseudosolenia.calcaravis") ) {
#
#                                              if(sp == "Coscinodiscus.oculusiridis") {
#                                                  sp2 <- "Coscinodiscus_oculus-iridis"
#                                              } else if(sp == "Pseudonitzschia.delicatissima") {
#                                                  sp2 <- "Pseudo-nitzschia.delicatissima"
#                                              } else if(sp == "Pseudonitzschia.lineola") {
#                                                  sp2 <- "Pseudo-nitzschia.lineola"
#                                              } else if(sp == "Pseudonitzschia.pungens") {
#                                                  sp2 <- "Pseudo-nitzschia.pungens"
#                                              } else if(sp == "Pseudonitzschia.seriata") {
#                                                  sp2 <- "Pseudo-nitzschia.seriata"
#                                              } else if(sp == "Pseudosolenia.calcaravis") {
#                                                  sp2 <- "Pseudosolenia.calcar-avis"
#                                              } # eo else if loop
#
#                                              setwd( paste(paste("proj_projection_",gsub("\\.","_",sp2),"_",m, sep = ""),"/", sep = "") )
#                                              d <- get(load( paste("proj_projection_", gsub("\\.","_",sp2),"_",m,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
#                                              resModel <- d[,sdm,,]
#                                              resModel <- apply(resModel, 1, mean, na.rm = F)
#                                              resModel <- (resModel/1000)
#
#                                          } else {
#
#                                              setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m, sep = ""),"/", sep = "") )
#                                              d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_",sp, ".RData", sep = "") ))
#                                              resModel <- d[,sdm,,]
#                                              resModel <- apply(resModel, 1, mean, na.rm = F)
#                                              resModel <- (resModel/1000)
#                                          }
#
#                                          # Return
#                                          return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y, species = gsub("\\.","_",sp), HSI = resModel) )
#
#                                      } else {
#                                          message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))
#                                      } # eo if else loop
#
#                              } # eo FUN
#
#                      ) # eo lapply
#                      tbl <- dplyr::bind_rows(probas)
#                      # dim(tbl) ; colnames(tbl)
#                      rm(probas); gc()
#                      # Dcast to put species as columns
#                      d_tbl <- dcast(tbl, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
#
#                      ### Save !
#                      message(paste("Saving community table for | ",sdm," | ",m, sep = ""))
#                      message(paste(" ", sep = ""))
#                      setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/",cat,"/group_background/communities/climate_change_projections", sep = ""))
#                      write.table(d_tbl, file = paste("table_",cat,"_mon_composition_baseline_",sdm,"_",m,".txt", sep = ""), sep = "\t")
#
#                      setwd(data.wd.back)
#                      rm(d_tbl,tbl) ; gc()
#
#              } # eo for loop - sdm in SDMs
#
# } # eo FUN - baseline.community.extracter
#
#
# #categories <- c("phytoplankton","zooplankton")

# cat <- "zooplankton"
# strat <- "group"

# #for(cat in categories) {
#
#      message(paste(" ", sep = ""))
#      message(paste("EXTRACTING ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
#      message(paste(" ", sep = ""))
#
#      # Set the working directories
#      if(cat == "phytoplankton" & strat == "total") {
#          setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
#          data.wd.back <- getwd()
#      } else if(cat == "zooplankton" & strat == "total") {
#          setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
#          data.wd.back <- getwd()
#      } else if(cat == "phytoplankton" & strat == "group") {
#          setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background")
#          data.wd.back <- getwd()
#      } else if(cat == "zooplankton" & strat == "group") {
#          setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
#          data.wd.back <- getwd()
#      }
#
#      # Apply fun in mclapply
#      require("parallel")
#      mclapply(X = months, FUN = baseline.community.extracter, mc.cores = 15)
#
# #} # eo

### --------------------------------------------------------------------------------------------------------------------------------

### Same but per each ESM now 
# esm <- "GFDL-TOPAZ"
# m <- "apr"
# sdm <- "GAM"

future.community.extracter <- function(m = months) {
    
            message(paste(" ", sep = ""))
            message(paste("Retrieving baseline zoo probabilities for ", m, sep = ""))
            message(paste(" ", sep = ""))
             
            # And for each SDM
            for(sdm in SDMs) {
                    
                    message(paste("SDM == ",sdm, sep = ""))
                    
                    for(esm in ESMs) {
                    
                            message(paste("ESM == ",esm, sep = ""))
                            setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims_constant_SST/", sep = ""))
                            mm <- firstup(m)
                            env <- read.table(paste("clims_mon_",mm,"_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")    
                            
                            setwd(paste(data.wd.back,"/niche.modelling/future_projections_constant_SST/", sep = ""))
                            spp <- dir()
                            spp <- gsub("\\(|\\)", "", spp)
                            
                            require("parallel")
                            # sp <- spp[1]
                            probas <- lapply(X = spp, FUN = function(sp) {

                                    # Got to species dir
                                    setwd(paste(data.wd.back,"/niche.modelling/future_projections_constant_SST/",sp,"/", sep = ""))
                                    message(paste("Loading projections for ", sp, "  ================================  ", sep = ""))

                                    # If the 12 future monthly projections are done
                                    if( sum(grepl("proj_projection_", dir())) == 72 ) {
                                            
                                        # Load projections
                                        if(sp %in% c("Coscinodiscus.oculusiridis","Pseudonitzschia.delicatissima","Pseudonitzschia.lineola",
                                                    "Pseudonitzschia.pungens","Pseudonitzschia.seriata","Pseudosolenia.calcaravis") ) {
                                                
                                            if(sp == "Coscinodiscus.oculusiridis") {
                                                sp2 <- "Coscinodiscus_oculus-iridis"
                                            } else if(sp == "Pseudonitzschia.delicatissima") {
                                                sp2 <- "Pseudo-nitzschia.delicatissima"
                                            } else if(sp == "Pseudonitzschia.lineola") {
                                                sp2 <- "Pseudo-nitzschia.lineola"
                                            } else if(sp == "Pseudonitzschia.pungens") {
                                                sp2 <- "Pseudo-nitzschia.pungens"
                                            } else if(sp == "Pseudonitzschia.seriata") {
                                                sp2 <- "Pseudo-nitzschia.seriata"
                                            } else if(sp == "Pseudosolenia.calcaravis") {
                                                sp2 <- "Pseudosolenia.calcar-avis"
                                            } # eo else if loop
                                         
                                            setwd( paste(paste("proj_projection_",gsub("\\.","_",sp2),"_",m,"_2100_",esm, sep = ""), sep = "") )
                                            d <- get(load( paste("proj_projection_", gsub("\\.","_",sp2),"_",m,"_2100_",esm,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
                                            resModel <- d[,sdm,,]
                                            resModel <- apply(resModel, 1, mean, na.rm = F)
                                            resModel <- (resModel/1000)
                                                
                                        } else {
                                                
                                            setwd( paste(paste("proj_projection_",gsub("\\.","_",sp),"_",m,"_2100_",esm, sep = ""), sep = "") )
                                            d <- get(load( paste("proj_projection_", gsub("\\.","_",sp),"_",m,"_2100_",esm,"_", gsub("\\(|)","",sp), ".RData", sep = "") ))
                                            resModel <- d[,sdm,,]
                                            resModel <- apply(resModel, 1, mean, na.rm = F)
                                            resModel <- (resModel/1000)
                                                
                                        }
                                             
                                            # Return
                                            return( data.frame(cell_id = paste(env$x, env$y, sep = "_"), x = env$x, y = env$y,species = gsub("\\.","_",sp), HSI = resModel) )
                                                    
                                    } else {

                                        message(paste("Skipping for not projection (yet)", sp, "  ================================", sep = ""))

                                    } # eo if else loop

                                } # eo FUN
                
                            ) # eo lapply
                            tbl <- dplyr::bind_rows(probas)
                            # dim(tbl) ; colnames(tbl)
                            rm(probas); gc()
                            # Dcast to put species as columns 
                            d_tbl <- dcast(tbl, cell_id + x + y ~ species, fun.aggregate = mean, na.rm = T, value.var = "HSI")
                            # dim(d_tbl) ; colnames(d_tbl)
                            
                            ### Save !
                            message(paste("Saving baseline community table for | ",sdm," | ",m," | ",esm, sep = ""))
                            message(paste(" ", sep = ""))
                            setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/",cat,"/group_background/communities/climate_change_projections", sep = ""))
                            write.table(d_tbl, file = paste("table_",cat,"_mon_composition_2100-2000_",esm,"_",sdm,"_",m,"_constant_SST.txt", sep = ""), sep = "\t")
                            
                            setwd(data.wd.back)
                            rm(d_tbl,tbl) ; gc()
                    
                   } # eo for loop - esm in ESMs
                
        } # eo for loop - sdm in SDMs
                        
} # eo FUN - future.community.extracter

categories <- c("phytoplankton","zooplankton")
#cat <- "zooplankton"
strat <- "group"

for(cat in categories) {
    
	message(paste(" ", sep = ""))	
	message(paste("EXTRACTING ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
	message(paste(" ", sep = ""))
    
    # Set the working directories
    if(cat == "phytoplankton" & strat == "total") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
        data.wd.back <- getwd()
    } else if(cat == "zooplankton" & strat == "total") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
        data.wd.back <- getwd()
    } else if(cat == "phytoplankton" & strat == "group") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background")
        data.wd.back <- getwd()
    } else if(cat == "zooplankton" & strat == "group") {
        setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
        data.wd.back <- getwd()
    } 
    
    # Apply fun in mclapply
    require("parallel")
    mclapply(X = months, FUN = future.community.extracter, mc.cores = 15)

} # eo - cat in categories

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


### A°) Concatenate all baseline mothly projections in one table (for diversity mapping after)
categories <- c("phytoplankton","zooplankton")
strat <- "group"

#cat <- "zooplankton"

# for(cat in categories) {
#
#         message(paste(" ", sep = ""))
#         message(paste("CONCATENATING ALL ",cat," BASELINE COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
#
#         # Go to the communities wd
#         if(cat == "phytoplankton" & strat == "total") {
#             setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/communities/climate_change_projections")
#             data.wd.back <- getwd()
#         } else if(cat == "zooplankton" & strat == "total") {
#             setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/communities/climate_change_projections")
#             data.wd.back <- getwd()
#         } else if(cat == "phytoplankton" & strat == "group") {
#             setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/communities/climate_change_projections")
#             data.wd.back <- getwd()
#         } else if(cat == "zooplankton" & strat == "group") {
#             setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/communities/climate_change_projections")
#             data.wd.back <- getwd()
#         } # eo else if loop
#
#         files <- dir()[grep("baseline",dir())]
#         # f <- files[3]
#         require("parallel")
#         comms <- mclapply(files, function(f) {
#                     d <- read.table(f, h = T, sep = "\t")
#                     terms <- do.call(cbind,strsplit(f,"_"))
#                     d$SDM <- terms[6,1]
#                     d$month <- str_replace_all(terms[7,1],".txt","")
#                     return(d)
#             } , mc.cores = 30
#         ) # eo lapply - f in files
#         table <- bind_rows(comms)
#         rm(comms) ; gc()
#
#         # Save table
#         message(paste("SAVING ALL ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
#         message(paste(" ", sep = ""))
#         write.table(x = table, paste("table_mon_composition_baseline_all_",cat,"_",strat,".txt", sep = ""), sep = "\t")
#
# } # eo for loop - c in combin


### B°) Do the same but for future monthly projections
#cat <- "phytoplankton"
strat <- "group"

for(cat in categories) {

        message(paste(" ", sep = ""))	
	    message(paste("CONCATENATING ALL ",cat," FUTURE COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
        
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
        
        files <- dir()[grep("_constant_SST",dir())] # files
        # f <- files[1]
        require("parallel")
        require("dplyr")
        comms <- mclapply(files, function(f) {
                    d <- read.table(f, h = T, sep = "\t")
                    terms <- do.call(cbind,strsplit(f,"_"))
                    d$SDM <- terms[7,1]
                    d$ESM <- terms[6,1]
                    d$month <- str_replace_all(terms[8,1],".txt","")
                    return(d)
            } , mc.cores = 30
        ) # eo lapply - f in files   
        table <- bind_rows(comms)
        rm(comms) ; gc()
        # dim(table) ; unique(table$month)
        # Save table
	    message(paste("SAVING ALL ",cat," COMMUNITIES BASED ON ",strat," BACKGROUND ", sep = ""))
	    message(paste(" ", sep = ""))
        write.table(x = table, paste("table_mon_composition_2100-2000_all_",cat,"_",strat,"_constant_SST.txt", sep = ""), sep = "\t")
    
} # eo for loop - c in combin


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


