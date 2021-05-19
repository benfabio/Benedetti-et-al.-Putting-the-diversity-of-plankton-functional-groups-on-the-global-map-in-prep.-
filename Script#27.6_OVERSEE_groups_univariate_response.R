
##### 06/10/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Loading the group-background species data for niche modelling
#   - Re-running GAMS to derive univarite response curves and env niche traits (center and width)
#   - Plot curves to make sure they make sense 
 
### Last update: 13/10/2020 

# --------------------------------------------------------------------------------------------------------------------------------

library("raster")
library("sp")
library("reshape2")
library("tidyverse")
library("biomod2")
library("R.devices")

WD <- getwd()
world2 <- map_data("world2")
world <- map_data("world")

# --------------------------------------------------------------------------------------------------------------------------------

# Vector of SDMs
SDMs <- c('ANN','GLM')
#SDMs <- "ANN"
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


# f <- files[5] ; f
#
# GetResponseCurves <- function(f = files) {
#
#                             # Get the data
#                             setwd(data.wd.back)
#                             message(paste("Modelling ",f, " =========================================================", sep = ""))
#                             data <- read.table(f, h = T, sep = ";")
#                             g <- unique(data[data$obs == 1,"group2"])
#
#                             ### Change data colnames for: logChla --> logChl
#                             colnames(data)[grep("logChla", colnames(data))] <- "logChl"
#                             colnames(data)[grep("Chla", colnames(data))] <- "Chl"
#
#                             data2 <- na.omit( data[,c("species","genus","x","y","obs")] )
#                             n <- nrow( data2[data2$obs == 1,] )
#                             rm(data2)
#
#                             ### Issue with some group names
#                             if( g == "Coccolithophores" & unique(data[data$obs == 1,"group2"]) == "bacillariophyceae" ) {
#                                 g <- "Diatoms"
#                                 #levels(data[data$obs == 1,"group2"])[levels(data[data$obs == 1,"group2"]) == "Coccolithophores"] <- g
#                             }
#
#                             ### Issue with some group names
#                             if( g == "Coccolithophores" & unique(data[data$obs == 1,"group2"]) == "dinoflagellata" ) {
#                                 g <- "Dinoflagellates"
#                                 #levels(data[data$obs == 1,"group2"])[levels(data[data$obs == 1,"group2"]) == "Coccolithophores"] <- g
#                             }
#
#                             # If n >= 75, continue, otherwise got to next species
#                             if( n >= 75 & g != "Other_zoo" & g != "Other_phyto" ) {
#
#                                     ### Initialisation: data formatting
#                                 myRespName <- str_replace_all(unique(data[data$obs == 1,"species"]), "_", ".")
#                                 myRespName <- gsub("\\(|\\)", "", myRespName)
#
#                                 if( grepl(pattern = '-', x = myRespName, fixed = T) ) {
#                                     myRespName <- str_replace_all(myRespName,'-','')
#                                 }
#
#                                 # the presence/absences data for our species
#                                 myResp <- as.numeric(data$obs)
#                                 # the XY coordinates of species data
#                                 myRespXY <- data[,c("x","y")]
#
#                                 if(g %in% c("Dinoflagellates","Diatoms","Coccolithophores")) {
#
#                                     vars <- c("SST","Nstar","logSiO2","logChl","Wind","PAR","Sistar","logEKE")
#
#                                 } else if(g %in% c("Calanoida","Poecilostomatoida","Oithonida","Jellyfish",
#                                                     "Chaetognatha","Euphausiids","Foraminifera","Amphipods",
#                                                     "Pteropods","Salps","Appendicularians")) {
#
#                                     vars <- c("SST","logSiO2","dO2","Sistar","logChl","Wind","logEKE","PAR")
#
#                                 } # eo if else loop
#
#                                 myExpl <- data[,vars]
#
#                                 # weights vector
#                                 weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
#
#                                 # data formating
#                                 setwd(second.wd)
#
#                                     myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
#
#                                     ### Species Distribution Modelling
#                                 nulldev()
#                                 myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, models = SDMs, models.options = myBiomodOption, NbRunEval = length(eval_runs),
#                                    DataSplit = 80, Yweights = weights, VarImport = 0, models.eval.meth = c('ROC','TSS'), SaveObj = F, do.full.models = F
#                                 ) # eo modelling
#                                 dev.off()
#
#                                 ### Get evaluation scores
#                                 scores <- data.frame(get_evaluations(myBiomodModelOut))
#                                 ### Get response curves
#                                 sdm <- "GAM"
#                                 response.plot <- BIOMOD_LoadModels(myBiomodModelOut, models = sdm)
#
#                                 # Use lapply to summarize all scores (10 scores per 6 SDMs) AND PROBABILITY THRESH
#                                 scores.ls <- lapply(SDMs, function(sdm) {
#                                                 # Retrieve all scores for the sdm 's'
#                                                 tss <- scores["TSS",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
#                                                 tss_cutoffs <- scores["TSS",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
#                                                 colnames(tss) <- paste(sdm, eval_runs, sep = "_")
#                                                 colnames(tss) <- paste(sdm, eval_runs, sep = "_")
#                                                 tss <- t(tss)
#                                                 tss_cutoffs <- t(tss_cutoffs)
#                                                 colnames(tss_cutoffs) <- "Cutoff_TSS"
#                                                 return(cbind(tss, tss_cutoffs))
#                                 } ) # eo lapply
#                                 scores.tbl <- data.frame(do.call(rbind, scores.ls))
#                                 scores.tbl$species <- unique(data[data$obs == 1,"species"])
#                                 scores.tbl$group <- g
#                                 rm(scores,scores.ls)
#
#                                 # Get response plot from biomod2's special function
#                                 nulldev()
#                                 myRespPlot2D <- response.plot2(models = response.plot,
#                                             Data = get_formal_data(myBiomodModelOut,'expl.var'), plot = F,
#                                             show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
#                                             do.bivariate = FALSE,fixed.var.metric= 'mean',legend = TRUE,
#                                             data_species = get_formal_data(myBiomodModelOut,'resp.var')
#                                 ) # eo myRespPlot2D
#                                 dev.off()
#
#                                 # Transform into a data.frame for ggploting
#                                 resp <- data.frame(myRespPlot2D)
#
#                                 # Changes colnames (the order in 'vars' object dictates the sequenc eof response curves)
#                                 colnames(resp) <- paste(rep(vars, each = 11, times = 1), sdm, eval_runs, sep = "_")
#                                 colnames(resp)[c(1,12,23,34,45,56,67,78)] <- vars
#
#                                 # Then for each predictor variable (n = 8), need to derive the mean response curve
#                                 avg.resp <- lapply(vars, function(v) {
#                                         # v <- "Sistar"
#                                         sub <- resp[,grep(v, colnames(resp))]
#                                         sub[,paste(v,"_avg_resp",sep="")] <- rowMeans(sub[,c(2:11)], na.rm = F) # average across the 10 CV runs
#                                         return(sub[,c(1,12)])
#                                         }
#                                 ) # eo lapply - v in vars
#                                 # cbind
#                                 table.resp <- do.call(cbind,avg.resp)
#                                 rm(avg.resp,myRespPlot2D,resp) ; gc()
#
#                                 ### Derive univariate niche traits !
#                                 require("parallel")
#                                 niche.traits <- mclapply(vars, function(v) {
#
#                                         # v <- "Sistar" # for testing
#                                         sub <- table.resp[,grep(v, colnames(table.resp))]
#                                         colnames(sub) <- c('x','y')
#
#                                         ### Niche center = median
#                                         sub$area <- cumsum(sub$y)
#                                         medarea <- (sub[length(sub$area),"area"] / 2)
#                                         sub$dist1 <- abs(sub$area - medarea)
#                                         center <- sub[which(sub$dist1 == min(sub$dist1)),'x']
#                                         if( length(center) > 1) {
#                                                 center <- mean(center)
#                                         } # eo if loop
#
#                                         ### Find x value where y is maximal
#                                         center2 <- sub[which(sub$y == max(sub$y)),'x']
#                                         if( length(center2) > 1) {
#                                                 center2 <- mean(center2)
#                                         } # eo if loop
#
#                                         ### Niche breadth = interquantile range between the 10th and the 90th
#                                         q1 <- (sub[length(sub$area),"area"] / 10)
#                                         q3 <- (9*(sub[length(sub$area),"area"]) / 10)
#                                         sub$dist2 <- abs(sub$area - q1)
#                                         sub$dist3 <- abs(sub$area - q3)
#                                         lower <- sub[which(sub$dist2 == min(sub$dist2)),'x']
#                                         if( length(lower) > 1) {
#                                                 lower <- mean(lower)
#                                         } # eo if loop
#                                         upper <- sub[which(sub$dist3 == min(sub$dist3)),'x']
#                                         if( length(upper) > 1) {
#                                                 upper <- mean(upper)
#                                         } # eo if loop
#                                         breadth <- upper - lower
#
#                                         # And try also weighted mean and sd to estimate center and breadth
#                                         require("Hmisc")
#                                         wtd.quantiles <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)))
#                                         center3 <- wtd.quantiles[3]  # weighted median
#                                         lower2 <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)), probs = seq(from = 0, to = 1, by = 0.1) )[2]
#                                         upper2 <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)), probs = seq(from = 0, to = 1, by = 0.1) )[9]
#                                         breadth2 <- upper2 - lower2
#
#                                         # And compute the range of 'y' to idetify which covariates have flat responses or not
#                                         # (could be used to weight in the PCA, or whatever analysis of niche traits variations)
#                                         range <- max(sub$y) - min(sub$y)
#
#                                         ### Summarize in a data.frame that you will return
#                                         # center = niche center based on median derived from cumsums, like in Benedetti et al., 2018
#                                         # center2 = value where 'y' (mean HSI) is the highest (rename 'optim')
#                                         # center3 = niche center based on weighted median (weights = normalized HSI) (rename center2)
#                                         # lower = 10th quantile based on cumsum area
#                                         # upper = 90th quantile based on cumsum area
#                                         # breadth = upper - lower
#                                         # lower2 = weighted quantiles
#                                         # upper2 = weighted quantiles
#                                         # breadth2 = upper2 - lower2
#
#                                         traits <- data.frame(species = unique(data[data$obs == 1,"species"]), group = g,
#                                             var = v, center = center, width = breadth, lower = lower, upper = upper,
#                                             optim = center2, center2 = center3, width2 = breadth2, lower2 = lower2,
#                                             upper2 = upper2, prob.range = range)
#
#                                         # Return niche traits ; not plot yet
#                                         return(traits)
#
#                                     }, mc.cores = 10
#
#                                 ) # eo lapply
#
#                                 ### Bind species straits
#                                 table.niche.traits <- bind_rows(niche.traits)
#                                 # Idea: if you plan to choose the niche traits in a PCA-like analysis, you might want to weight in their contribution according to how much HSI varies (so prob.range)
#                                 table.niche.traits$weight.pca <- (table.niche.traits$prob.range)/ max(table.niche.traits$prob.range)
#
#                                 rm(niche.traits) ; gc()
#
#                                 ### Make a panel plot with all of species resp curves (n = 8 per species) and highlighting their niche traits
#                                 ### Make use of 'table.resp' and 'table.niche.traits' in a lapply and return a 'list' that will contain the 8 plots,
#                                 ### arrange them in a panel with ggarrange
#                                 niche.plots <- lapply(vars, function(v) {
#
#                                             sub <- table.resp[,grep(v, colnames(table.resp))]
#                                             colnames(sub) <- c('x','y')
#                                             sub.traits <- table.niche.traits[table.niche.traits$var == v,]
#
#                                             plot <- ggplot() + geom_path(aes(x = x, y = y), linetype = "solid", data = sub, group = 1) +
#                                                       geom_vline(xintercept = sub.traits$optim, colour = "#a50026", linetype = "solid") +
#                                                       geom_vline(xintercept = sub.traits$center, colour = "#d73027", linetype = "solid") +
#                                                       geom_vline(xintercept = sub.traits$center2, colour = "#313695", linetype = "solid") +
#                                                       geom_vline(xintercept = sub.traits$lower, colour = "#f46d43", linetype = "dashed") +
#                                                       geom_vline(xintercept = sub.traits$upper, colour = "#f46d43", linetype = "dashed") +
#                                                       geom_vline(xintercept = sub.traits$lower2, colour = "#74add1", linetype = "dashed") +
#                                                       geom_vline(xintercept = sub.traits$upper2, colour = "#74add1", linetype = "dashed") +
#                                                       xlab(v) + ylab("Mean HSI") + scale_y_continuous(limits = c(0,1)) +
#                                                       theme_classic()
#
#                                             # Return
#                                             return(plot)
#
#                                     } # eo FUN
#
#                                 ) # eo lapply - v in vars
#                                 # class(niche.plots) ; str(niche.plots)
#                                 #niche.plots[[1]]
#
#                                 require("ggpubr")
#                                 panel <- ggarrange(niche.plots[[1]],niche.plots[[2]],niche.plots[[3]],niche.plots[[4]],
#                                                   niche.plots[[5]],niche.plots[[6]],niche.plots[[7]],niche.plots[[8]],
#                                          ncol = 4, nrow = 2, align = "hv")
#
#                                 message(paste("", sep = ""))
#                                 message(paste("Saving plots and tables for ",unique(data[data$obs == 1,"species"]), sep = ""))
#                                 message(paste("", sep = ""))
#
#                                 ### Print panel in appropriate dir
#                                 setwd(paste(data.wd.back,"/resp_curves/plots/", sep = ""))
#                                 ggsave(plot = panel, filename = paste("panel_respcurv_",sdm,"_group_",unique(data[data$obs == 1,"species"]),".jpg", sep = ""), dpi = 300, height = 4, width = 10)
#
#                                 ### Save niche traits table
#                                 setwd(paste(data.wd.back,"/resp_curves/niche_traits/", sep = ""))
#                                 save(table.niche.traits, file = paste("table_niche_traits_",sdm,"_group_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
#
#                                 ### Save evaluation scores
#                                 setwd(paste(data.wd.back,"/","Scores", sep = ""))
#                                 save(scores.tbl, file = paste("eval_scores_for_respcurv_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
#
#                                 rm(niche.plots, table.niche.traits, scores.tbl, table.resp)
#                                 gc()
#                                 setwd(second.wd)
#
#                             } else {
#
#                                 message(paste("WRONG GROUP ('Other_') OR JUST NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
#
#                             }
#
# } # eo GetResponseCurves


categories <- c("phytoplankton","zooplankton")

cat <- "zooplankton"
strat <- "group"

#require("parallel")
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
    # f <- files[242]
    for(f in files[c(243:length(files))] ) {
	
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
  
        				### Get evaluation scores
        				scores <- data.frame(get_evaluations(myBiomodModelOut))
        				# Use lapply to summarize all scores 
        				scores.ls <- lapply(SDMs, function(sdm) {
        							# Retrieve all scores for the sdm 's' 
        							tss <- scores["TSS",grep(paste("Testing.data.",sdm, sep = ""), colnames(scores))]
        							tss_cutoffs <- scores["TSS",grep(paste("Cutoff.",sdm, sep = ""), colnames(scores))]
        							colnames(tss) <- paste(sdm, eval_runs, sep = "_")
        							tss <- t(tss)
        							tss_cutoffs <- t(tss_cutoffs)
        							colnames(tss_cutoffs) <- "Cutoff_TSS"
        							return(cbind(tss, tss_cutoffs))
        					} 
                        ) # eo lapply
        				scores.tbl <- data.frame(do.call(rbind, scores.ls))
                        scores.tbl$species <- unique(data[data$obs == 1,"species"])
                        scores.tbl$group <- g
        				rm(scores,scores.ls)
                        
                        ### Get response curves
                        # sdm <- "ANN"
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
                                        # v <- "Sistar"
                                        sub <- resp[,grep(v, colnames(resp))]
                                        sub[,paste(v,"_avg_resp",sep="")] <- rowMeans(sub[,c(2:11)], na.rm = F) # average across the 10 CV runs
                                        return(sub[,c(1,12)])	    
                                }
                            ) # eo lapply - v in vars     
                            # cbind
                            table.resp <- do.call(cbind,avg.resp)     
                            rm(avg.resp,myRespPlot2D,resp) ; gc()

                            ### Derive univariate niche traits !
                            require("parallel")
                            niche.traits <- mclapply(vars, function(v) {
                                        
                                                    # v <- "Sistar" # for testing
                                                    sub <- table.resp[,grep(v, colnames(table.resp))]
                                                    colnames(sub) <- c('x','y')
        								
                                                    ### Niche center = median
                    								sub$area <- cumsum(sub$y)
                    								medarea <- (sub[length(sub$area),"area"] / 2)
                    								sub$dist1 <- abs(sub$area - medarea)
                    								center <- sub[which(sub$dist1 == min(sub$dist1)),'x']
                    								if( length(center) > 1) {
                    										center <- mean(center)
                    								} # eo if loop
                                        
                                                    ### Find x value where y is maximal 
                                                    center2 <- sub[which(sub$y == max(sub$y)),'x']
                    								if( length(center2) > 1) {
                    										center2 <- mean(center2)
                    								} # eo if loop
                                        
                    								### Niche breadth = interquantile range between the 10th and the 90th
                    								q1 <- (sub[length(sub$area),"area"] / 10)
                    								q3 <- (9*(sub[length(sub$area),"area"]) / 10)
                    								sub$dist2 <- abs(sub$area - q1)
                    								sub$dist3 <- abs(sub$area - q3)
                    								lower <- sub[which(sub$dist2 == min(sub$dist2)),'x']
                    								if( length(lower) > 1) {
                    										lower <- mean(lower)
                    								} # eo if loop
                    								upper <- sub[which(sub$dist3 == min(sub$dist3)),'x']
                    								if( length(upper) > 1) {
                    										upper <- mean(upper)
                    								} # eo if loop
                    								breadth <- upper - lower
                                        
                                                    # And try also weighted mean and sd to estimate center and breadth
                                                    require("Hmisc")
                                                    wtd.quantiles <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)))
                                                    center3 <- wtd.quantiles[3]  # weighted median
                                                    lower2 <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)), probs = seq(from = 0, to = 1, by = 0.1) )[2]
                                                    upper2 <- wtd.quantile(x = sub$x, weights = ((sub$y)/max(sub$y)), probs = seq(from = 0, to = 1, by = 0.1) )[9]
                                                    breadth2 <- upper2 - lower2
                                        
                                                    # And compute the range of 'y' to idetify which covariates have flat responses or not
                                                    # (could be used to weight in the PCA, or whatever analysis of niche traits variations)
                                                    range <- max(sub$y) - min(sub$y)
                                        
                                                    ### Summarize in a data.frame that you will return
                                                    # center = niche center based on median derived from cumsums, like in Benedetti et al., 2018
                                                    # center2 = value where 'y' (mean HSI) is the highest (rename 'optim')
                                                    # center3 = niche center based on weighted median (weights = normalized HSI) (rename center2)
                                                    # lower = 10th quantile based on cumsum area
                                                    # upper = 90th quantile based on cumsum area
                                                    # breadth = upper - lower
                                                    # lower2 = weighted quantiles
                                                    # upper2 = weighted quantiles
                                                    # breadth2 = upper2 - lower2
                                                                                
                                                    traits <- data.frame(species = unique(data[data$obs == 1,"species"]), group = g,
                                                        var = v, center = center, width = breadth, lower = lower, upper = upper, 
                                                        optim = center2, center2 = center3, width2 = breadth2, lower2 = lower2, 
                                                        upper2 = upper2, prob.range = range)
                                            
                                                    # Return niche traits ; not plot yet 
                                                    return(traits)
        								
                                                }, mc.cores = 10
                                    
                                ) # eo lapply
                                
                                ### Bind species straits
                                table.niche.traits <- bind_rows(niche.traits)
                                # Idea: if you plan to choose the niche traits in a PCA-like analysis, you might want to weight in their contribution according to how much HSI varies (so prob.range)
                                table.niche.traits$weight.pca <- (table.niche.traits$prob.range)/ max(table.niche.traits$prob.range)    
                                         
                                rm(niche.traits) ; gc()	
                                
                                ### Make a panel plot with all of species resp curves (n = 8 per species) and highlighting their niche traits
                                ### Make use of 'table.resp' and 'table.niche.traits' in a lapply and return a 'list' that will contain the 8 plots, 
                                ### arrange them in a panel with ggarrange
                                niche.plots <- lapply(vars, function(v) { 
                                
                                                sub <- table.resp[,grep(v, colnames(table.resp))]
                                                colnames(sub) <- c('x','y')
                                                sub.traits <- table.niche.traits[table.niche.traits$var == v,]
                                            
                                                plot <- ggplot() + geom_path(aes(x = x, y = y), linetype = "solid", data = sub, group = 1) +
                                                            geom_vline(xintercept = sub.traits$optim, colour = "#a50026", linetype = "solid") +
                                                            geom_vline(xintercept = sub.traits$center, colour = "#d73027", linetype = "solid") +
                                                            geom_vline(xintercept = sub.traits$center2, colour = "#313695", linetype = "solid") +
                                                            geom_vline(xintercept = sub.traits$lower, colour = "#f46d43", linetype = "dashed") +
                                                            geom_vline(xintercept = sub.traits$upper, colour = "#f46d43", linetype = "dashed") +
                                                            geom_vline(xintercept = sub.traits$lower2, colour = "#74add1", linetype = "dashed") +
                                                            geom_vline(xintercept = sub.traits$upper2, colour = "#74add1", linetype = "dashed") +
                                                            xlab(v) + ylab("Mean HSI") + scale_y_continuous(limits = c(0,1)) +
                                                            theme_classic()
                                                      
                                                # Return 
                                                return(plot)
                                
                                        } # eo FUN 
                                    
                                ) # eo lapply - v in vars
                                
                                require("ggpubr")
                                panel <- ggarrange(niche.plots[[1]],niche.plots[[2]],niche.plots[[3]],niche.plots[[4]],
                                            niche.plots[[5]],niche.plots[[6]],niche.plots[[7]],niche.plots[[8]],
                                            ncol = 4, nrow = 2, align = "hv")
                                         
                                message(paste("", sep = ""))
                                message(paste("Saving plots and tables for ",unique(data[data$obs == 1,"species"]), sep = ""))
                                message(paste("", sep = ""))
                                         
                                ### Print panel in appropriate dir
                                setwd(paste(data.wd.back,"/resp_curves/plots/", sep = ""))
                                ggsave(plot = panel, filename = paste("panel_respcurv_",sdm,"_group_",unique(data[data$obs == 1,"species"]),".jpg", sep = ""), dpi = 300, height = 4, width = 10)
                                
                                ### Save niche traits table
                                setwd(paste(data.wd.back,"/resp_curves/niche_traits/", sep = ""))
            					save(table.niche.traits, file = paste("table_niche_traits_",sdm,"_group_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
 
                                rm(niche.plots, table.niche.traits, table.resp)
                                gc()
                                setwd(second.wd)		
                            
                        } # eo for loop sdm in SDMs
                        
    					### Save evaluation scores
    					setwd(paste(data.wd.back,"/","Scores", sep = ""))
    					save(scores.tbl, file = paste("eval_scores_for_respcurv_",sdm,"_",unique(data[data$obs == 1,"species"]),".Rdata", sep = ""))
                        rm(scores.tbl)
                        gc()
                        setwd(second.wd)
								
        			} else {
							
        				    message(paste("WRONG GROUP ('Other_') OR JUST NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
        			}
        
    } # eo for loop - f in files
    
    
# } # eo for loop - cat in categories
    
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
