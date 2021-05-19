
### ================================================================================================================

library("raster")
library("sp")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")
library("R.devices")

WD <- getwd()

### ================================================================================================================

# Vector of SDMs
SDMs <- c('GLM','GAM','ANN')
# Vector of eval_runs :
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 

### Load the env stacks (winter  and summer, like January & August) for projection
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")
apr <- read.table("glob_stack_month_apr_18_09_20.txt", h = T, sep = ";")
jul <- read.table("glob_stack_month_jul_18_09_20.txt", h = T, sep = ";")
oct <- read.table("glob_stack_month_oct_18_09_20.txt", h = T, sep = ";")
jan <- read.table("glob_stack_month_jan_18_09_20.txt", h = T, sep = ";")

apr <- apr[-which(apr$SSS < 20),] 
apr <- apr[-which(apr$Bathy > -175),] 
jul <- jul[-which(jul$SSS < 20),] 
jul <- jul[-which(jul$Bathy > -175),] 
oct <- oct[-which(oct$SSS < 20),] 
oct <- oct[-which(oct$Bathy > -175),] 
jan <- jan[-which(jan$SSS < 20),] 
jan <- jan[-which(jan$Bathy > -175),] 

### Load the remaining 8 climatologies for annual proj of diversity
annual <- TRUE
#annual <- FALSE

if(annual == TRUE) {
		
		feb <- read.table("glob_stack_month_feb_18_09_20.txt", h = T, sep = ";")
		mar <- read.table("glob_stack_month_mar_18_09_20.txt", h = T, sep = ";")
		may <- read.table("glob_stack_month_may_18_09_20.txt", h = T, sep = ";")
		jun <- read.table("glob_stack_month_jun_18_09_20.txt", h = T, sep = ";")
		aug <- read.table("glob_stack_month_aug_18_09_20.txt", h = T, sep = ";")
		sep <- read.table("glob_stack_month_sep_18_09_20.txt", h = T, sep = ";")
		nov <- read.table("glob_stack_month_nov_18_09_20.txt", h = T, sep = ";")
		dec <- read.table("glob_stack_month_dec_18_09_20.txt", h = T, sep = ";")
	
		feb <- feb[-which(feb$SSS < 20),] 
		feb <- feb[-which(feb$Bathy > -175),] 
		mar <- mar[-which(mar$SSS < 20),] 
		mar <- mar[-which(mar$Bathy > -175),] 
		may <- may[-which(may$SSS < 20),] 
		may <- may[-which(may$Bathy > -175),] 
		jun <- jun[-which(jun$SSS < 20),] 
		jun <- jun[-which(jun$Bathy > -175),] 
		aug <- aug[-which(aug$SSS < 20),] 
		aug <- aug[-which(aug$Bathy > -175),] 
		sep <- sep[-which(sep$SSS < 20),] 
		sep <- sep[-which(sep$Bathy > -175),] 
		nov <- nov[-which(nov$SSS < 20),] 
		nov <- nov[-which(nov$Bathy > -175),] 
		dec <- dec[-which(dec$SSS < 20),] 
		dec <- dec[-which(dec$Bathy > -175),]

} # eo if loop : annual switch


setwd(WD)


### Set modelling options
myBiomodOption <- BIOMOD_ModelingOptions(
	
						GLM = list( type = 'quadratic',
    								interaction.level = 0,
   			 						myFormula = NULL,
    								test = 'AIC',
    								family = binomial("logit"),
    								mustart = 0.5,
    								control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),

						GAM = list(algo = 'GAM_mgcv',
								type = 's_smoother',
								k = 5,
								interaction.level = 0,
								myFormula = NULL,
								family = binomial("logit"),
								method = 'GCV.Cp',
								optimizer = c('outer','newton'),
								select = FALSE,
								knots = NULL,
								paraPen = NULL,
								control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
								, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
								, rank.tol = 1.49011611938477e-08
								, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
								, optim = list(factr=1e+07)
								, newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE) ), 
                                
							 GBM = list( distribution = 'bernoulli',
							 	n.trees = 100,
							  	interaction.depth = 1,
							  	n.minobsinnode = 10,
							  	shrinkage = 0.001,
							  	bag.fraction = 0.5,
							  	train.fraction = 1,
							  	cv.folds = 2,
							  	keep.data = FALSE,
							  	verbose = FALSE,
							  	perf.method = "cv") 
								
) # eo modelling options

# f <- files[330]
# f <- files2[30]

SpeciesNicheModelling <- function(f = files) {
	
	                        # Get the data
	                        setwd(data.wd.back)
	                        message(paste("Modelling ",f, " =========================================================", sep = ""))
	                        data <- read.table(f, h = T, sep = ";")
                            g <- unique(data[data$obs == 1,"group2"])
                            
                            ### Change data colnames for: logChla --> logChl
                            colnames(data)[grep("logChla", colnames(data))] <- "logChl"
                            colnames(data)[grep("Chla", colnames(data))] <- "Chl"
                            
	                        data2 <- na.omit( data[,c("species","genus","x","y","obs")] )
	                        n <- nrow( data2[data2$obs == 1,] )
	                        rm(data2)
                            
                            ### Issue with some group names
                            if( g == "Coccolithophores" & unique(data[data$obs == 1,"group2"]) == "bacillariophyceae" ) {
                                g <- "Diatoms"
                                #levels(data[data$obs == 1,"group2"])[levels(data[data$obs == 1,"group2"]) == "Coccolithophores"] <- g
                            }
                            
                            ### Issue with some group names
                            if( g == "Coccolithophores" & unique(data[data$obs == 1,"group2"]) == "dinoflagellata" ) {
                                g <- "Dinoflagellates"
                                #levels(data[data$obs == 1,"group2"])[levels(data[data$obs == 1,"group2"]) == "Coccolithophores"] <- g
                            }
							
							# If n >= 85, continue, otherwise got to next species
							if( n >= 75 & g != "Other_zoo" & g != "Other_phyto" ) {
								
	  			  				### Initialisation: data formatting
								myRespName <- str_replace_all(unique(data[data$obs == 1,"species"]), "_", ".")
								myRespName <- gsub("\\(|\\)", "", myRespName)
                                
                                if( grepl(pattern = '-', x = myRespName, fixed = T) ) {
                                     myRespName <- str_replace_all(myRespName,'-','')   
                                }
                            
								# the presence/absences data for our species 
								myResp <- as.numeric(data$obs)
								# the XY coordinates of species data
								myRespXY <- data[,c("x","y")]
                                
                                # Environmental variables as a function of the group
                                if(g == "Diatoms" & strat == "total") {
                                    vars <- c("SST","Nstar","logChl","pCO2","PAR","logSiO2")
                                } else if(g == "Diatoms" & strat == "group") {
                                    vars <- c("SST","Nstar","logChl","Sistar","PAR","logSiO2")
                                } else if(g == "Dinoflagellates" & strat == "total") {
                                    vars <- c("SST","Nstar","Wind","PAR","MLD","logChl")
                                } else if(g == "Dinoflagellates" & strat == "group") {
                                    vars <- c("SST","Nstar","Wind","PAR","MLD","logChl")
                                } else if(g == "Coccolithophores" & strat == "total") {
                                    vars <- c("SST","Nstar","Wind","logChl","Sistar","pCO2")
                                } else if(g == "Coccolithophores" & strat == "group") {
                                    vars <- c("SST","Nstar","Wind","PAR","MLD","pCO2")
                                } else if(g == "Calanoida" & strat == "total") {
                                    vars <- c("SST","logSiO2","Sistar","dO2","MLD","PAR")
                                } else if(g == "Calanoida" & strat == "group") {
                                    vars <- c("SST","logSiO2","Sistar","dO2","MLD","PAR")
                                } else if(g == "Poecilostomatoida" & strat == "total") {
                                    vars <- c("SST","Sistar","logChl","PAR","MLD","logSiO2")
                                } else if(g == "Poecilostomatoida" & strat == "group") {
                                    vars <- c("SST","Sistar","Wind","PAR","MLD","logSiO2")
                                } else if(g == "Oithonida" & strat == "total") {
                                    vars <- c("SST","logSiO2","Sistar","dO2","MLD","logChl")
                                } else if(g == "Oithonida" & strat == "group") {
                                    vars <- c("SST","logSiO2","Sistar","dO2","MLD","PAR")
                                } else if(g == "Jellyfish" & strat == "total") {
                                    vars <- c("SST","logSiO2","Sistar","dO2","MLD","logChl")
                                } else if(g == "Jellyfish" & strat == "group") {
                                    vars <- c("SST","Nstar","MLD","PAR","dO2","logSiO2")
                                } else if(g == "Chaetognatha" & strat == "total") {
                                    vars <- c("SST","Sistar","logSiO2","logEKE","MLD","PAR")
                                } else if(g == "Chaetognatha" & strat == "group") {
                                    vars <- c("SST","Sistar","logSiO2","logEKE","MLD","PAR")
                                } else if(g == "Euphausiids" & strat == "total") {
                                    vars <- c("SST","Sistar","logSiO2","logEKE","MLD","PAR")
                                } else if(g == "Euphausiids" & strat == "group") {
                                    vars <- c("SST","Sistar","logSiO2","dO2","MLD","Nstar")
                                } else if(g == "Foraminifera" & strat == "total") {
                                    vars <- c("SST","dO2","PAR","logSiO2","Wind","Sistar")
                                } else if(g == "Foraminifera" & strat == "group") {
                                    vars <- c("SST","dO2","PAR","logSiO2","Wind","Sistar")
                                } else if(g == "Amphipods" & strat == "total") {
                                    vars <- c("SST","logSiO2","Sistar","PAR","dO2","MLD")
                                } else if(g == "Amphipods" & strat == "group") {
                                    vars <- c("SST","logSiO2","Sistar","PAR","dO2","logChl")
                                } else if(g == "Pteropods" & strat == "total") {
                                    vars <- c("SST","dO2","Sistar","logSiO2","PAR","logEKE")
                                } else if(g == "Pteropods" & strat == "group") {
                                    vars <- c("SST","dO2","Sistar","logSiO2","PAR","MLD")
                                } else if(g == "Salps" & strat == "total") {
                                    vars <- c("SST","dO2","Sistar","logSiO2","PAR","logEKE")
                                } else if(g == "Salps" & strat == "group") {
                                   vars <- c("SST","dO2","Sistar","logSiO2","PAR","MLD")
                                } else if(g == "Appendicularians" & strat == "total") {
                                    vars <- c("SST","dO2","Sistar","logSiO2","PAR","logEKE")
                                } else if(g == "Appendicularians" & strat == "group") {
                                    vars <- c("SST","Sistar","MLD","Wind","PAR","logSiO2")
                                } 
                                
								myExpl <- data[,vars]

								# weights vector
								weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
		
								# formatage des données
                                setwd(second.wd)
                                
	  			  				myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
							
	 			   				### Modelling
                                nulldev()
								myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = SDMs, models.options = myBiomodOption, NbRunEval = length(eval_runs), 
								   DataSplit = 80, Yweights = weights, VarImport = 0, models.eval.meth = c('ROC','TSS'), SaveObj = F, do.full.models = F
								) # eo modelling
                                dev.off()					
  
								### Get evaluation scores
								scores <- data.frame(get_evaluations(myBiomodModelOut))
								# Use lapply to summarize all scores (10 scores per 6 SDMs) AND PROBABILITY THRESH
								scores.ls <- lapply(SDMs, function(sdm) {
												# Retrieve all scores for the sdm 's' 
												tss <- scores["TSS",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
												tss_cutoffs <- scores["TSS",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
												colnames(tss) <- paste(sdm, eval_runs, sep = "_")
												colnames(tss) <- paste(sdm, eval_runs, sep = "_")
												tss <- t(tss)
												tss_cutoffs <- t(tss_cutoffs)
												colnames(tss_cutoffs) <- "Cutoff_TSS"
												return(cbind(tss, tss_cutoffs))
								} ) # eo lapply
								scores.tbl <- data.frame(do.call(rbind, scores.ls))
                                scores.tbl$species <- unique(data[data$obs == 1,"species"])
                                scores.tbl$group <- g
								rm(scores,scores.ls)
							
	  			  				### Make projections (they will be printed in the species folder)
								message(paste("Projecting ",myRespName, " =========================================================", sep = ""))
                                
                                # dim(apr[,vars]) ; summary(apr[,vars])                       
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = apr[,vars],
                                                     proj.name = paste("projection",myRespName,"apr", sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
    
                                 # Project niches in July conditions
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = jul[,vars],
                                                     proj.name = paste("projection",myRespName,"jul", sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
    
                                 # Project niches in October conditions
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = oct[,vars],
                                                     proj.name = paste("projection",myRespName,"oct", sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
    
                                 # Project niches in January conditions
                                 myBiomodProj <- BIOMOD_Projection(
                                                     modeling.output = myBiomodModelOut,
                                                     new.env = jan[,vars],
                                                     proj.name = paste("projection",myRespName,"jan", sep = "_"),
                                                     selected.models = 'all',
                                                     binary.meth = 'TSS',
                                                     compress = 'xz',
                                                     clamping.mask = FALSE
                                 ) # eo projection
							
								### Create a switch for annual projections (annual = T/F)
								if( annual == TRUE ) {
                                    
									# Feb
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = feb[,vars],
														proj.name = paste("projection",myRespName,"feb", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Mar
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = mar[,vars],
														proj.name = paste("projection",myRespName,"mar", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# May
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = may[,vars],
														proj.name = paste("projection",myRespName,"may", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Jun
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = jun[,vars],
														proj.name = paste("projection",myRespName,"jun", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Aug
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = aug[,vars],
														proj.name = paste("projection",myRespName,"aug", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection	
									# Sep
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = sep[,vars],
														proj.name = paste("projection",myRespName,"sep", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Nov
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = nov[,vars],
														proj.name = paste("projection",myRespName,"nov", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Dec
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = dec[,vars],
														proj.name = paste("projection",myRespName,"dec", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									
								} else { 
                                    
									message(paste("============== NOT PERFORMING ANNUAL PROJECTIONS ==============", sep = ""))
                                    
								} # eo if else loop

								### Save evaluation scores
								setwd(paste(data.wd.back,"/","Scores", sep = ""))
								save(scores.tbl, file = paste("eval_scores_for_projections_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
								setwd(WD)  			
								
							} else {
							
								message(paste("WRONG GROUP ('Other_') OR JUST NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
							}
  
} # eo NicheModelling


strategies <- c("group")
categories <- c("Phytoplankton")

cat <- "Phytoplankton"
strat <- "group"

for(cat in categories) {


    for(strat in strategies) {
	
    	message(paste(" ", sep = ""))	
    	message(paste("RUNNING ",cat," PROJECTIONS BASED ON ",strat," BACKGROUND ", sep = ""))
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
        
        setwd(paste(data.wd.back,"/","niche.modelling/projections", sep = ""))
        second.wd <- getwd()
        
        # List files of species that have been projected already
        modelled.spp <- dir()
        res <- lapply(modelled.spp, function(s) {
                    setwd(paste(data.wd.back,"/","niche.modelling/projections/",s,"/", sep = ""))
                    sum <- sum(grepl("proj_projection_", dir()))
                    if(annual == TRUE) {
                        nProjs <- 12
                    } else {
                        nProjs <- 4
                    } # eo if else loop
                    
                    if(sum == nProjs) {
                         message(paste(s," has ",sum," monthly projections", sep = ""))   
                         return(data.frame(s = "no"))
                    } else {
                        message(paste(s," NOT ENOUGH MONTHLY PROJECTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!", sep = ""))
                        return(data.frame(s = s))
                    } # eo if else loop
                    
                } # eo FUN
        ) # eo lapply
        table <- do.call(rbind, res); rm(res); gc()
        table <- table[table$s != "no",]
        #spp2model <- unique(table$s)
        # Replace "." by "_"
        spp2remodel <- gsub("\\.", "_", x = as.character(table))
        
        setwd(data.wd.back)
        
        #files <- dir()[grep("data_",dir())]
        #spp2model <- str_replace_all(files,paste("data_",strat,"_", sep = ""),"")
        #spp2model <- str_replace_all(spp2model,".txt","")
        
        ### Sometimes the server shuts off and interrupts the niche modelling...need to identify the species that are already modelled and projected and remove them from the 'species' vector
         #setwd(paste(data.wd.back,"/","Scores", sep = ""))
         #  dir() # just need to replace the dot by an underscore and there you go
         #donespp <- dir()[grep("for_projections", dir())] # But also need to add brackets for some annoying species names
         # remove "eval_scores_for_projections_" and ".Rdata"
         #donespp <- str_replace_all(donespp,"eval_scores_for_projections_","")
         #donespp <- str_replace_all(donespp,".Rdata","")

         # Identify species in spp2model that are not part of donespp and re-build their filenames
         if( length(spp2remodel) > 0 ) {
            #spp2remodel <- spp2model[!(spp2model %in% donespp)]
            files2 <- paste("data_",strat,"_",spp2remodel,".txt", sep = "")
         }

         if( exists("files2") & length(files2) > 0 ) {

             require("parallel")
             mclapply(X = files2, SpeciesNicheModelling, mc.cores = 30)

         } else if( length(donespp) == 0 ) {
            
            require("parallel")
    	    mclapply(X = files, SpeciesNicheModelling, mc.cores = 30)
            
        } # eo if else loop - files2 exists
        
    
    } # eo strat in strategies   
	
    
} # eo cat in categories


### ==============================================================================================================================
### ==============================================================================================================================
### ==============================================================================================================================
