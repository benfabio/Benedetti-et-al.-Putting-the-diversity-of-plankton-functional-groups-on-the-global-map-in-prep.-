
# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")
library("parallel")
library("corrgram")
library("biomod2")

world2 <- map_data("world2")
world <- map_data("world")
WD <- getwd() 

# --------------------------------------------------------------------------------------------------------------------------------

myBiomodOption <- BIOMOD_ModelingOptions(
	
						GLM = list(type = 'quadratic', interaction.level = 0, myFormula = NULL, test = 'AIC', family = binomial("logit"),
    								mustart = 0.5, control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),

						GAM = list(algo = 'GAM_mgcv', type = 's_smoother',
								k = 5,interaction.level = 0, myFormula = NULL,
								family = binomial("logit"),method = 'GCV.Cp',
								optimizer = c('outer','newton'), select = FALSE,
								knots = NULL, paraPen = NULL,
								control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
								, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
								, rank.tol = 1.49011611938477e-08
								, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
								, optim = list(factr=1e+07)
								, newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE) )
								
) # eo modelling options


VarsRanker <- function(f = files) {
	
							# Get the data
							setwd(data.wd.back)
							message(paste("Modelling ",f, " =========================================================", sep = ""))
							
                            data <- read.table(f, h = T, sep = ";")
							data2 <- na.omit( data[,c(vars,"obs")] )
							n <- nrow( data2[data2$obs == 1,] )
							rm(data2)
							
							# If n >= 85, continue, otherwise got to next species
							if(n >= 75) {
								
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

								# Environmental variables
								myExpl <- data[,vars]

								# weights vector
								weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
		
								# formatage des données
	  			  				myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
                                
                                setwd(paste(data.wd.back,"/","niche.modelling/",set, sep = ""))
								
                                myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = SDMs, models.options = myBiomodOption, 
								                          NbRunEval = length(eval_runs), DataSplit = 80, 
								                          Yweights = weights, VarImport = 10, models.eval.meth = c('TSS'),
								                          SaveObj = FALSE, do.full.models = FALSE
								) # eo modelling
                                
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
                                scores.tbl$group <- unique(data[data$obs == 1,"group2"])
								rm(scores,scores.ls)
							
	  			  				### ?variables_importance
                                ranks <- data.frame(get_variables_importance(myBiomodModelOut))
                                ranks$predictor <- rownames(ranks)
                                m.ranks <- melt(ranks, id.var = "predictor")
                                colnames(m.ranks) <- c("predictor","model","rank")
                                m.ranks$model <- gsub(".", "_", m.ranks$model, fixed = T)
                                m.ranks$SDM <- data.frame(do.call(rbind,strsplit(as.character(m.ranks$model), split = "_")))[,1]
                                m.ranks$RUN <- data.frame(do.call(rbind,strsplit(as.character(m.ranks$model), split = "_")))[,2]
                                m.ranks$species <- unique(data[data$obs == 1,"species"])
                                m.ranks$group <- unique(data[data$obs == 1,"group2"])
                                m.ranks$Set <- set
                                                                     
								### Save evaluation scores
								setwd(paste(data.wd.back,"/","Ranks", sep = ""))
                                save(m.ranks, file = paste("table_ranks_",set,"_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
                                
                                setwd(paste(data.wd.back,"/","Scores", sep = ""))
								save(scores.tbl, file = paste("eval_scores_",set,"_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
                                
								setwd(WD)  			
								
							} else {
							
								message(paste("NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
							}
  
} # eo NicheModelling

sets <- c("Set1","Set2")
categories <- c("Phytoplankton","Zooplankton")

strategies <- c("total","group")

SDMs <- c('GLM','GAM','ANN')
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5") 

# For testing
# cat <- "Phytoplankton"
# set <- "Set1"
# strategy <- "total"

for(strategy in strategies) {

    for(cat in categories) {
    
        message(paste(" ", sep = ""))
        message(paste(" ", sep = ""))
        message(paste("QUANTIFYING VARIABLES IMPORTANCE FOR ", cat, sep = ""))
        message(paste(" ", sep = ""))
        message(paste(" ", sep = ""))
    
        for(set in sets) {

            message(paste(" ", sep = ""))
            message(paste("BASED ON ", set," and ",strategy,"-background", sep = ""))
            message(paste(" ", sep = ""))

            if(cat == "Phytoplankton" & strategy == "total") {
    
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
                data.wd.back <- getwd()
                Set1 <- c("SST", "logChla", "pCO2", "logNO3", "Nstar", "Sistar", "MLPAR", "logEKE", "Wind")
                Set2 <- c("SST", "logChla", "pCO2", "logSiO2", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "Wind")
    
            } else if(cat == "Zooplankton" & strategy == "total") {
    
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
                data.wd.back <- getwd()
                Set1 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLPAR", "logEKE", "Wind")
                Set2 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "Wind")
    
            } else if(cat == "Phytoplankton" & strategy == "group") {
    
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background")
                data.wd.back <- getwd()
                Set1 <- c("SST", "logChla", "pCO2", "logNO3", "Nstar", "Sistar", "MLPAR", "logEKE", "Wind")
                Set2 <- c("SST", "logChla", "pCO2", "logSiO2", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "Wind")
    
            } else if(cat == "Zooplankton" & strategy == "group") {
    
                setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
                data.wd.back <- getwd()
                Set1 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLPAR", "logEKE", "Wind")
                Set2 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "Wind")
    
            } 
        
            ### 29/08/2020: Check which species/ files did not return any scores/ ranks table
            if(cat == "Phytoplankton" & strategy == "total") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/Ranks")
                    files <- dir()[grep(set,dir())]
            } else if(cat == "Zooplankton" & strategy == "total") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/Ranks")
                    files <- dir()[grep(set,dir())]
            } else if(cat == "Phytoplankton" & strategy == "group") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/Ranks")
                    files <- dir()[grep(set,dir())]
            } else if(cat == "Zooplankton" & strategy == "group") {
                    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/Ranks")
                    files <- dir()[grep(set,dir())]
            } # eo if else loop
        
            # Retrieve the names of the species modelled 
            res <- lapply(files, function(f) { d <- get(load(f)) ; return(data.frame(sp = unique(d$species)))} )
            spmodelled <- do.call(rbind,res) ; rm(res) ; gc() ; spmodelled <- unique(spmodelled$sp)
        
            if( length(spmodelled) < 3 ) {
            
                setwd(data.wd.back)
                files <- dir()[grep("data_",dir())]
            
            } else {
            
                # And retrieve all species that were supposed to be modelled 
                setwd(data.wd.back)
                files <- dir()[grep("data_",dir())]
                res <- lapply(files, function(f) { d <- read.table(f, h = T, sep = ";") ; return(data.frame(sp = unique(d[d$obs == 1,"species"]))) } )
                sp2model <- do.call(rbind,res) ; rm(res) ; gc() ; sp2model <- unique(sp2model$sp)
                # Find non overlapping species
                commons <- intersect(sp2model,spmodelled)
                missing <- sp2model[!(sp2model %in% commons)]
            
                files <- paste("data_",strategy,"_",missing,".txt", sep = "")   
                rm(missing, commons, sp2model) ; gc()
            
            } # eo if else loop 
    
	        # Set pool of variables
	        if(set == "Set1") {
		        vars <- Set1
	        } else if (set == "Set2") {
		        vars <- Set2
	        } # eo if else loop 
        
	        require("parallel")
            # f <- files[8]
	        mclapply(X = files, FUN = VarsRanker, mc.cores = 35)
	
        } # eo set in sets

    } # eo cat in categories

} # eo strat in strategies



# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
