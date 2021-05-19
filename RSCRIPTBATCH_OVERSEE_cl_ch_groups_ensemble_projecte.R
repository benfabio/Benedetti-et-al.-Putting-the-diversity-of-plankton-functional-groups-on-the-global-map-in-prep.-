
### ================================================================================================================

library("raster")
library("sp")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")

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


setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims_constant_SST", sep = ""))

ESMs <- c("CNRM-PISCES","IPSL-PISCES","GFDL-TOPAZ","CESM-BEC","MRI-NEMURO")	
future <- TRUE

future_clims_ESMs <- lapply(ESMs, function(esm) {
    
            message(paste("Loading clims for ", esm, sep = ""))
            setwd(paste("/net/kryo/work/fabioben/OVERSEE/data/future/MAREMIP_data/future_mon_clims_constant_SST/", sep = ""))
            listMonths <- list()
    		future.jan <- read.table(paste("clims_mon_Jan_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.feb <- read.table(paste("clims_mon_Feb_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.mar <- read.table(paste("clims_mon_Mar_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.apr <- read.table(paste("clims_mon_Apr_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.may <- read.table(paste("clims_mon_May_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.jun <- read.table(paste("clims_mon_Jun_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.jul <- read.table(paste("clims_mon_Jul_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.aug <- read.table(paste("clims_mon_Aug_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.sep <- read.table(paste("clims_mon_Sep_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.oct <- read.table(paste("clims_mon_Oct_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.nov <- read.table(paste("clims_mon_Nov_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
    		future.dec <- read.table(paste("clims_mon_Dec_",esm,"_rcp85_base+2100-2031_constant_SST.txt", sep = ""), h = T, sep = "\t")
            # Add a month id
            future.jan$month <- "Jan"
            future.feb$month <- "Feb"
            future.mar$month <- "Mar"
            future.apr$month <- "Apr"
            future.may$month <- "May"
            future.jun$month <- "Jun"
            future.jul$month <- "Jul"
            future.aug$month <- "Aug"
            future.sep$month <- "Sep"
            future.oct$month <- "Oct"
            future.nov$month <- "Nov"
            future.dec$month <- "Dec"
            # And supply to listMonths
            listMonths[[1]] <- future.jan
            listMonths[[2]] <- future.feb
            listMonths[[3]] <- future.mar
            listMonths[[4]] <- future.apr
            listMonths[[5]] <- future.may
            listMonths[[6]] <- future.jun
            listMonths[[7]] <- future.jul
            listMonths[[8]] <- future.aug
            listMonths[[9]] <- future.sep
            listMonths[[10]] <- future.oct
            listMonths[[11]] <- future.nov
            listMonths[[12]] <- future.dec
            # dimnames(listMonths) <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
            # And return
            return(listMonths)
    
        } # eo FUN

) # eo lapply
names(future_clims_ESMs) <- ESMs

### 25/09/18: Sometimes the server shuts off and interrupts the niche modelling...need to identify the species that are already modelled and projected and remove them from the 'species' vector
#  setwd(second.wd)
#  dir() # just need to replace the dot by an underscore and there you go
#  donespp <- gsub("\\.", "_", dir()) # But also need to add brackets for some annoying species names
#  for(i in 1:length(donespp)) {
#     	sp <- donespp[i]
#     	if( length( strsplit(sp, "_", fixed = TRUE)[[1]] ) == 3) {
#     		  # Then add brackets around the second piece
#     		sp <- paste(strsplit(sp, "_", fixed = TRUE)[[1]][1],"_(",
#     				strsplit(sp, "_", fixed = TRUE)[[1]][2],")_",
#     				strsplit(sp, "_", fixed = TRUE)[[1]][3], sep = "")
#
#     		donespp[i] <- sp
#    		}
#  }  #eo for loop
#  # OK Now remove the 'donespp' from 'species'
#  species2 <- species[!(species %in% donespp)]
# #
# # # 25/09/18: And SOMETIMES, the models wouldn't get projected even though they're done !
#  setwd(second.wd)
#  modelled_spp <- dir()
#  list <- lapply(modelled_spp, function(s) {
#     				 # Go to species dir and list the number of elements
#     				 # First if loop to add brackets when necessary...
#     				if( length( strsplit(s, "_", fixed = TRUE)[[1]] ) == 3) {
#     					 # Then add brackets around the second piece
#     					s <- paste(strsplit(s, "_", fixed = TRUE)[[1]][1],"_(",
#     							strsplit(s, "_", fixed = TRUE)[[1]][2],")_",
#     							strsplit(s, "_", fixed = TRUE)[[1]][3], sep = "")
#     				}
#
#     				setwd( paste(second.wd, "/", str_replace_all(s, "_", "."), "/", sep = "") )
#     				n <- length(dir())
#     				 # If there are more than 5 elements (n >= 5), then it means the the projections have already been carried out sucessfully
#     				if(n >= 5) {
#     					rm(n)
#     				} else {
#     					return(s)
#     				}  # eo seond if else loop
# } ) # eo lapply
#
#  unproj_spp <- do.call(rbind, list)[,1]
#  #  So species2model is a combination of species2 & unproj_spp
#  species2model <- gsub("\\.", "_", c(species2, unproj_spp))
# # species2model <- gsub("\\.", "_", unproj_spp[13:24])
# # species2model <- gsub("\\.", "_", unproj_spp)
#  rm(modelled_spp, donespp)

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
								, outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE) )	
								
) # eo modelling options

# f <- files[5]

SpeciesNicheModelling <- function(f = files) {
	
                    	    # Get the data
                    	    setwd(data.wd.back)
                    	    message(paste("Modelling ",f, " =========================================================", sep = ""))
                    	    data <- read.table(f, h = T, sep = ";")
                            g <- unique(data[data$obs == 1,"group2"])
                            sp <- unique(data[data$obs == 1,"species"])
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
	
							# If n >= 85, continue, otherwise got to next species
							if( n >= 75 & g != "Other_zoo" & g != "Other_phyto" ) {
								
                                if(g == "Appendicularians") {
                                    # remove dO2 for appendicularia
                                    vars <- c("SST","logSiO2","logChl","PAR","Sistar")
                                } # eo if loop
                                
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
                              	# explanatory variables
                                myExpl <- data[,vars]

                				# weights vector
                				weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
		
                				# data formating
                                setwd(second.wd)
           
								# formatage des données
	  			  				myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
	                                      	 			 	expl.var = myExpl,
	                                       			  	 	resp.xy = myRespXY,
	                                       			  	 	resp.name = myRespName 
								)
							
	 			   				### Modelling
								setwd(second.wd)
								myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
								                                     models = SDMs, 
								                                     models.options = myBiomodOption, 
								                                     NbRunEval = length(eval_runs), 
								                                     DataSplit = 80, 
								                                     Yweights = weights,
								                                     VarImport = 0, 
								                                     models.eval.meth = c('ROC','TSS'),
								                                     SaveObj = FALSE,
								                                     do.full.models = FALSE
								) # eo modelling					
  		
	  			  				### Make projections (they will be printed in the species folder)
								message(paste("Projecting ",sp, " =========================================================", sep = ""))
								# Project niches in April conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = apr[,vars],
													proj.name = paste("projection",sp,"apr", sep = "_"),
													selected.models = 'all',
													binary.meth = 'TSS',
													compress = 'xz',
													clamping.mask = FALSE 
								) # eo projection
								# Project niches in July conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = jul[,vars],
													proj.name = paste("projection",sp,"jul", sep = "_"),
													selected.models = 'all',
													binary.meth = 'TSS',
													compress = 'xz',
													clamping.mask = FALSE 
								) # eo projection			
								# Project niches in Octpber conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = oct[,vars],
													proj.name = paste("projection",sp,"oct", sep = "_"),
													selected.models = 'all',
													binary.meth = 'TSS',
													compress = 'xz',
													clamping.mask = FALSE 
								) # eo projection			
								# Project niches in January conditions
								myBiomodProj <- BIOMOD_Projection(
													modeling.output = myBiomodModelOut,
													new.env = jan[,vars],
													proj.name = paste("projection",sp,"jan", sep = "_"),
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
														proj.name = paste("projection",sp,"feb", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Mar
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = mar[,vars],
														proj.name = paste("projection",sp,"mar", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# May
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = may[,vars],
														proj.name = paste("projection",sp,"may", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Jun
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = jun[,vars],
														proj.name = paste("projection",sp,"jun", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection
									# Aug
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = aug[,vars],
														proj.name = paste("projection",sp,"aug", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									) # eo projection	
									# Sep
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = sep[,vars],
														proj.name = paste("projection",sp,"sep", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Nov
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = nov[,vars],
														proj.name = paste("projection",sp,"nov", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									# Dec
									myBiomodProj <- BIOMOD_Projection(
														modeling.output = myBiomodModelOut,
														new.env = dec[,vars],
														proj.name = paste("projection",sp,"dec", sep = "_"),
														selected.models = 'all',
														binary.meth = 'TSS',
														compress = 'xz',
														clamping.mask = FALSE 
									)  # eo projection
									
								} else { 
									message(paste("============== NOT PERFORMING ANNUAL PROJECTIONS ==============", sep = ""))
								} # eo if else loop
								
								if(future == TRUE) {
                                    
                                    for(esm in ESMs) {
                                        
                                        # esm <- "IPSL-PISCES"
                                        message(paste("PROJECTING NICHES IN FUTURE CONDITIONS OF ", esm, sep = ""))
                                        
                                        future.apr <- future_clims_ESMs[esm][[1]][[4]]
                                        # dim(future.apr); summary(future.apr)
    									# Project niches in future April conditions
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.apr[,vars],
    														proj.name = paste("projection",sp,"apr_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection
                                        # Project niches in future July conditions
                                        future.jul <- future_clims_ESMs[esm][[1]][[7]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.jul[,vars],
    														proj.name = paste("projection",sp,"jul_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection			
                                        # Project niches in future October conditions
                                        future.oct <- future_clims_ESMs[esm][[1]][[10]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.oct[,vars],
    														proj.name = paste("projection",sp,"oct_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection		
    									# Project niches in future January conditions
                                        future.jan <- future_clims_ESMs[esm][[1]][[1]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.jan[,vars],
    														proj.name = paste("projection",sp,"jan_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection
    									# Feb
                                        future.feb <- future_clims_ESMs[esm][[1]][[2]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.feb[,vars],
    														proj.name = paste("projection",sp,"feb_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection
    									# Mar
                                        future.mar <- future_clims_ESMs[esm][[1]][[3]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.mar[,vars],
    														proj.name = paste("projection",sp,"mar_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection
    									# May
                                        future.may <- future_clims_ESMs[esm][[1]][[5]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.may[,vars],
    														proj.name = paste("projection",sp,"may_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection
    									# Jun
                                        future.jun <- future_clims_ESMs[esm][[1]][[6]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.jun[,vars],
    														proj.name = paste("projection",sp,"jun_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection
    									# Aug
                                        future.aug <- future_clims_ESMs[esm][[1]][[8]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.aug[,vars],
    														proj.name = paste("projection",sp,"aug_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									) # eo projection	
    									# Sep
                                        future.sep <- future_clims_ESMs[esm][[1]][[9]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.sep[,vars],
    														proj.name = paste("projection",sp,"sep_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									)  # eo projection
    									# Nov
                                        future.nov <- future_clims_ESMs[esm][[1]][[11]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.nov[,vars],
    														proj.name = paste("projection",sp,"nov_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									)  # eo projection
    									# Dec
                                        future.dec <- future_clims_ESMs[esm][[1]][[12]]
    									myBiomodProj <- BIOMOD_Projection(
    														modeling.output = myBiomodModelOut,
    														new.env = future.dec[,vars],
    														proj.name = paste("projection",sp,"dec_2100",esm, sep = "_"),
    														selected.models = 'all',
    														binary.meth = 'TSS',
    														compress = 'xz',
    														clamping.mask = FALSE 
    									)  # eo projection
                                        
                                    } # eo for loop - esm in ESMs
                                    
								} else {
                                    
									message(paste("============== NOT PERFORMING FUTURE PROJECTIONS ==============", sep = ""))
								
                                }

								### Save evaluation scores
								#setwd( paste(WD,"/","niche.modelling_future_",rcp,"_ensemble/eval_score_",p,"/", sep = "") )
								#save(scores.tbl, file = paste("eval_scores_",sp,".Rdata", sep = ""))
								setwd(WD)  			
								
							} else {
							
								message(paste("NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
							}
  
} # eo NicheModelling

### Set the dataset and the background selection strategy
categories <- c("Phytoplankton","Zooplankton")
strat <- "group"
rcp <- "rcp85"

# cat <- "Zooplankton"

for(cat in categories) {
	
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
    
    files <- dir()[grep("data_group",dir())]
    
    setwd(paste(data.wd.back,"/","niche.modelling/future_projections_constant_SST", sep = ""))
    second.wd <- getwd()
    setwd(data.wd.back)

    if(cat == "Phytoplankton") {
        vars <- c("SST","PAR","logSiO2","logChl","Nstar","Sistar")
    } else if(cat == "Zooplankton") {
        vars <- c("SST","dO2","logSiO2","logChl","PAR","Sistar")
    }

	require("parallel")
	message(paste(" ", sep = ""))	
	message(paste("RUNNING ",cat," PROJECTIONS ", sep = ""))
	message(paste(" ", sep = ""))
	
	mclapply(X = files, SpeciesNicheModelling, mc.cores = 30)
	
} # eo for cat in categories


### ==============================================================================================================================
### ==============================================================================================================================
### ==============================================================================================================================
