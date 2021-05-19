
##### 06/10/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Loading the group-background species data for niche modelling
#   - Re-running GAMS to derive univarite response curves and env niche traits (center and width)
#   - Plot curves to make sure they make sense 
#   - Re-run to save the average resp curves as well
 
### Last update: 20/10/2020 

# --------------------------------------------------------------------------------------------------------------------------------

library("reshape2")
library("tidyverse")
library("biomod2")
library("R.devices")
library("parallel")

WD <- getwd()

# --------------------------------------------------------------------------------------------------------------------------------

# Vector of SDMs
#SDMs <- c('ANN','GLM')
SDMs <- "GAM"
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5","RUN6","RUN7","RUN8","RUN9","RUN10") 

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
								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE))
								
) # eo modelling options

#categories <- c("phytoplankton","zooplankton")

cat <- "zooplankton"
strat <- "group"

#for(cat in categories) {
    
	message(paste(" ", sep = ""))	
	message(paste("RUNNING ",cat," PROJECTIONS BASED ON ",strat," BACKGROUND ", sep = ""))
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
    
    files <- dir()[grep("data_group",dir())]
    
    setwd(paste(data.wd.back,"/","niche.modelling/response_curves", sep = ""))
    second.wd <- getwd()
    
    ### 08/10/2020: Re-start script by omitting Jaschnovia_tolli [242]
    # f <- files[5]
for(f in files[c(243:length(files))]) {
	
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
                }
                            
                ### Issue with some group names
                if( g == "Coccolithophores" & unique(data[data$obs == 1,"group2"]) == "dinoflagellata" ) {
                   g <- "Dinoflagellates"
                   #levels(data[data$obs == 1,"group2"])[levels(data[data$obs == 1,"group2"]) == "Coccolithophores"] <- g
                }
							
        		# If n >= 75, continue, otherwise got to next species
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
                                
                        if(g %in% c("Dinoflagellates","Diatoms","Coccolithophores")) {
                                    
                             vars <- c("SST","Nstar","logSiO2","logChl","Wind","PAR","Sistar","logEKE")
                                    
                        } else if(g %in% c("Calanoida","Poecilostomatoida","Oithonida","Jellyfish",
                                        "Chaetognatha","Euphausiids","Foraminifera","Amphipods",
                                        "Pteropods","Salps","Appendicularians")) {
                                                        
                             vars <- c("SST","logSiO2","dO2","Sistar","logChl","Wind","logEKE","PAR")
                                    
                        } # eo if else loop 
                                
        				myExpl <- data[,vars]

        				# weights vector
        				weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
		
        				# data formating
                        setwd(second.wd)
                                
        	  			myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
							
        	 			### Species Distribution Modelling
                        nulldev()
        				myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = SDMs, models.options = myBiomodOption, NbRunEval = length(eval_runs), 
        								   DataSplit = 80, Yweights = weights, VarImport = 0, models.eval.meth = c('ROC','TSS'), SaveObj = F, do.full.models = F
        				) # eo modelling
                        dev.off()
                                    
                        ### Get response curves
                        # sdm <- "GAM"
                        for(sdm in SDMs) {
                            
            				response.plot <- BIOMOD_LoadModels(myBiomodModelOut, models = sdm)
                                
        				
            				# Get response plot from biomod2's special function
                            nulldev()
            				myRespPlot2D <- response.plot2(models = response.plot,
            								Data = get_formal_data(myBiomodModelOut,'expl.var'), plot = F,
            								show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
            								do.bivariate = FALSE,fixed.var.metric= 'mean',legend = TRUE,
            								data_species = get_formal_data(myBiomodModelOut,'resp.var') 
                            ) # eo myRespPlot2D
                            dev.off()	

            				# Transform into a data.frame for ggploting							   
            				resp <- data.frame(myRespPlot2D)
                                
                            # Changes colnames (the order in 'vars' object dictates the sequenc eof response curves)
                            colnames(resp) <- paste(rep(vars, each = 11, times = 1), sdm, eval_runs, sep = "_")
                            colnames(resp)[c(1,12,23,34,45,56,67,78)] <- vars
                                
                            # Then for each predictor variable (n = 8), need to derive the mean response curve
                            avg.resp <- lapply(vars, function(v) {
                                    sub <- resp[,grep(v, colnames(resp))]
                                    sub[,paste(v,"_avg_resp",sep = "")] <- rowMeans(sub[,c(2:11)], na.rm = F) # average across the 10 CV runs
                                    return(sub[,c(1,12)])	    
                                }
                            ) # eo lapply - v in vars     
                            # cbind
                            table.resp <- do.call(cbind,avg.resp)     
                            rm(avg.resp,myRespPlot2D,resp) ; gc()
                            table.resp$species <- unique(data[data$obs == 1,"species"])
                            table.resp$group <- g
                            table.resp$category <- cat
           
                            ### Save mean response curves
                            setwd(paste(data.wd.back,"/resp_curves/niche_traits/", sep = ""))
            				save(table.resp, file = paste("table_mean_resp_curves_",sdm,"_group_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
 
                            rm(table.resp); gc()
                            setwd(second.wd)		
                            
                        } # eo for loop sdm in SDMs
                        
    				    setwd(second.wd)
								
        			} else {
							
        				    message(paste("WRONG GROUP ('Other_') OR JUST NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
        			}
        
} # eo for loop - f in files
    
    
    #} # eo for loop - cat in categories
    
    
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
