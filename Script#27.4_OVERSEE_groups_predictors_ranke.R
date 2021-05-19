
##### 27/08/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Assessing pairwise predictors correlations, from species-level to group level
#   - For total background data only, load the species data (1/0) and compute pairwise Spearman's rank corr coeff
#     between env covariates
#   - Examine results from predictors rankings (total background and group-background)
 
### Last update: 10/11/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")
library("parallel")
library("corrgram")

world2 <- map_data("world2")
world <- map_data("world")
WD <- getwd() 

### Draw the correlation coeff heatmaps for each net
get_lower_tri <- function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
# Utiliser la corrélation entre les variables
  # comme mésure de distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
cormat <- cormat[hc$order, hc$order]
}

# --------------------------------------------------------------------------------------------------------------------------------

### A°) First, Zooplankton
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
files <- dir()[grep("data_",dir())]# ; files
f <- files[5]

res <- mclapply(files, function(f) {
    
            # Message
            data <- read.table(f, sep = ";", h = T)
            message(paste(unique(data[data$obs == 1,'species']), sep = ""))
            names <- colnames(data)[c(19,23:33,35,37:41)]
            mydata <- na.omit(data[,names]) #  head(mydata)
            cormat <- round(cor(mydata, method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = TRUE)
            melted_cormat <- melted_cormat[!(melted_cormat$Var1 == melted_cormat$Var2),]
            
            melted_cormat$species <- unique(data[data$obs == 1,'species'])
            melted_cormat$group <- unique(data[data$obs == 1,'group2'])
            return(melted_cormat)
                
    }, mc.cores = 30
    
) # eo mclapply
# Rbind
table <- bind_rows(res) ; rm(res) ; gc()
dim(table) ; head(table)
summary(table)

### Compute mean or median correlation
table$id <- paste(table$Var1, table$Var2, sep = "_")

ddf <- data.frame(table %>% group_by(id) %>% summarize(V1 = unique(Var1), V2 = unique(Var2), mean = round(mean(value),2), median = round(median(value),2)) )
head(ddf) ; summary(ddf)

# ggplot(ddf, aes(factor(V2), factor(V1), fill = mean)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Mean\nRho") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(V2), factor(V1), label = mean), color = "black", size = 3)
#
ggplot(ddf, aes(factor(V2), factor(V1), fill = mean)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Median\nRho") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(V2), factor(V1), label = mean), color = "black", size = 3)


### B°) Phytoplankton
setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
files <- dir()[grep("data_",dir())] #; files
f <- files[100]

res <- mclapply(files, function(f) {
    
            # Message
            data <- read.table(f, sep = ";", h = T)
            message(paste(unique(data[data$obs == 1,'species']), sep = ""))
            names <- colnames(data)[c(18,22:32,34,36:40)]
            mydata <- na.omit(data[,names]) #  head(mydata)
            cormat <- round(cor(mydata, method = "spearman"),2)
            upper_tri <- get_upper_tri(cormat)
            melted_cormat <- melt(upper_tri, na.rm = TRUE)
            melted_cormat <- melted_cormat[!(melted_cormat$Var1 == melted_cormat$Var2),]
    
            melted_cormat$species <- unique(data[data$obs == 1,'species'])
            melted_cormat$group <- unique(data[data$obs == 1,'group2'])
            return(melted_cormat)
    
    }, mc.cores = 30
    
) # eo mclapply
table <- bind_rows(res) ; rm(res) ; gc()
dim(table) ; head(table)
table$id <- paste(table$Var1, table$Var2, sep = "")

### Compute mean or median correlation
table$id <- paste(table$Var1, table$Var2, sep = "_")

ddf <- data.frame(table %>% group_by(id) %>% summarize(V1 = unique(Var1), V2 = unique(Var2), mean = round(mean(value),2), median = round(median(value),2)) )
head(ddf) ; summary(ddf)

# ggplot(ddf, aes(factor(V2), factor(V1), fill = mean)) + geom_tile(color = "white") +
#     scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Mean\nRho") +
#     xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
#         axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
#     geom_text(aes(factor(V2), factor(V1), label = mean), color = "black", size = 3)

ggplot(ddf, aes(factor(V2), factor(V1), fill = median)) + geom_tile(color = "white") +
    scale_fill_gradient2(low = "#3288bd", high = "#d53e4f", mid = "white", midpoint = 0, limits = c(-1,1), name = "Median\nRho") +
    xlab("") + ylab("") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 10, hjust = 1),
        axis.text.y = element_text(vjust = 1, size = 10, hjust = 1)) + coord_fixed() +
    geom_text(aes(factor(V2), factor(V1), label = median), color = "black", size = 3)
    
### Use these to define pools of non collinear covariates

### Phytoplankton: discard O2, logPO4, SSS, Pstar
### Choose between NO3 and SiO2
### MLPAR vs MLD+PAR
# Set1: SST, logChla, pCO2, logNO3, Nstar, Sistar, MLPAR, logEKE, SLA, Wind
# Set2: SST, logChla, pCO2, logSiO2, Nstar, Sistar, MLD, PAR, logEKE, SLA, Wind

### Zooplankton: discard SSS, logPO4, Pstar
### Choose one nutrients field (SiO2 because less correlated to SST)
# Set1: SST, logChla, dO2, pCO2, logSiO2, Nstar, Sistar, MLPAR, logEKE, SLA, Wind
# Set2: SST, logChla, dO2, pCO2, logSiO2, Nstar, Sistar, MLD, PAR, logEKE, SLA, Wind

# --------------------------------------------------------------------------------------------------------------------------------

### Use the variable_importance function of biomod2 for GLM, GAMs and ANN to rank variables of set 1 vs 2 and examine ranks per groups
library("biomod2")
# ?variables_importance

setwd(WD)

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
							
								# the presence/absences data for our species 
								myResp <- as.numeric(data$obs)

								# the XY coordinates of species data
								myRespXY <- data[,c("x","y")]

								# Environmental variables
								myExpl <- data[,vars]

								# weights vector
								weights <- na.omit(data[,c(vars,"weights")])[,"weights"]
		
								# formatage des données
	  			  				myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,expl.var = myExpl, resp.xy = myRespXY, resp.name = myRespName)
                                
                                setwd(paste(data.wd.back,"/","niche.modelling/",set, sep = ""))
								
                                myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
								                                     models = SDMs, 
								                                     models.options = myBiomodOption, 
								                                     NbRunEval = length(eval_runs), 
								                                     DataSplit = 80, 
								                                     Yweights = weights,
								                                     VarImport = 10, 
								                                     models.eval.meth = c('ROC','TSS','KAPPA'),
								                                     SaveObj = FALSE,
								                                     do.full.models = FALSE
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
                                save(m.ranks, file = paste("table_ranks_",set,"_",sp,".Rdata", sep = ""))
                                
                                setwd(paste(data.wd.back,"/","Scores", sep = ""))
								save(scores.tbl, file = paste("eval_scores_",set,"_",sp,".Rdata", sep = ""))
                                
								setwd(WD)  			
								
							} else {
							
								message(paste("NOT ENOUGH PRESENCES WITH ENV PREDICTORS VALUES || n = ", n, sep = ""))
								
							}
  
} # eo NicheModelling

sets <- c("Set1","Set2")
categories <- c("Phytoplankton","Zooplankton")

SDMs <- c('GLM','GAM','ANN')
eval_runs <- c("RUN1","RUN2","RUN3","RUN4","RUN5") 

# For testing
strategy <- "total"
cat <- "Zooplankton"
set <- "Set1"

for(cat in categories) {
    
    message(paste("QUANTIFYING VARIABLES IMPORTANCE FOR ", cat, sep = ""))
    message(paste(" ", sep = ""))
    message(paste(" ", sep = ""))
    
    for(set in sets) {

        message(paste("BASED ON ", set, sep = ""))
        message(paste(" ", sep = ""))

        if(cat == "Phytoplankton" & strategy == "total") {
    
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
            data.wd.back <- getwd()
            Set1 <- c("SST", "logChla", "pCO2", "logNO3", "Nstar", "Sistar", "MLPAR", "logEKE", "SLA", "Wind")
            Set2 <- c("SST", "logChla", "pCO2", "logSiO2", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "SLA", "Wind")
    
        } else if(cat == "Zooplankton" & strategy == "total") {
    
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
            data.wd.back <- getwd()
            Set1 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLPAR", "logEKE", "SLA", "Wind")
            Set2 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "SLA", "Wind")
    
        } else if(cat == "Phytoplankton" & strategy == "group") {
    
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background")
            data.wd.back <- getwd()
            Set1 <- c("SST", "logChla", "pCO2", "logNO3", "Nstar", "Sistar", "MLPAR", "logEKE", "SLA", "Wind")
            Set2 <- c("SST", "logChla", "pCO2", "logNO3", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "SLA", "Wind")
    
        } else if(cat == "Zooplankton" & strategy == "group") {
    
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
            data.wd.back <- getwd()
            Set1 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLPAR", "logEKE", "SLA", "Wind")
            Set2 <- c("SST", "logChla", "dO2", "pCO2", "logSiO2", "Nstar", "Sistar", "MLD", "PAR", "logEKE", "SLA", "Wind")
    
        } 

        files <- dir()[grep("data_",dir())]
	
	    # Set pool of variables
	    if(set == "Set1") {
		    vars <- Set1
	    } else if (set == "Set2") {
		    vars <- Set2
	    }
        
	    require("parallel")
	    mclapply(X = files, FUN = VarsRanker, mc.cores = 30)
	
    } # eo set in sets

}




### ------------------------------------------------------------------------------------------------------------------------

### 28/08/2020: Code above was ran with R CMD BATCH. Below develop script for examining results in predictors rankings (and TSS) per: 
# - SDMs
# - Trophic level
# - Groups
# (and Sets)

strategy <- "group"
cat <- "Phytoplankton"
set <- "Set2"

if(cat == "Phytoplankton" & strategy == "total") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/Ranks")
    data.wd.back <- getwd()
} else if(cat == "Zooplankton" & strategy == "total") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/Ranks")
    data.wd.back <- getwd()
} else if(cat == "Phytoplankton" & strategy == "group") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/Ranks")
    data.wd.back <- getwd()
} else if(cat == "Zooplankton" & strategy == "group") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/Ranks")
    data.wd.back <- getwd()
} 

files <- dir()[grep(set,dir())] ; files

res <- lapply(files, function(f) { data <- get(load(f)) ; return(data) } ) # eo lapply
table <- bind_rows(res) ; gc() ; rm(res)
head(table) 

# ggplot(data = table[table$group != "Other_phyto",], aes(x = factor(predictor), y = rank, fill = factor(predictor))) + geom_boxplot(colour = "black") +
#     scale_fill_brewer(palette = "Spectral", name = "") + xlab("") + ylab("Predictors ranking") +
#     theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# #
# ggplot(data = table[table$group != "Other_phyto",], aes(x = factor(predictor), y = rank, fill = factor(predictor))) + geom_boxplot(colour = "black") +
#     scale_fill_brewer(palette = "Spectral", name = "") + xlab("") + ylab("Predictors ranking") +
#     theme_bw() + facet_grid( ~ factor(SDM)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# #
# ggplot(data = table[table$group != "Other_phyto",], aes(x = factor(predictor), y = rank, fill = factor(predictor))) + geom_boxplot(colour = "black") +
#     scale_fill_brewer(palette = "Spectral", name = "") + xlab("") + ylab("Predictors ranking") +
#     theme_bw() + facet_grid( ~ factor(group)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#
# ggplot(data = table[table$group != "Other_phyto",], aes(x = factor(predictor), y = rank, fill = factor(predictor))) + geom_boxplot(colour = "black") +
#     scale_fill_brewer(palette = "Spectral", name = "") + xlab("") + ylab("Predictors ranking") +
#     theme_bw() + facet_grid(factor(SDM) ~ factor(group)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    


### 23/10/2020: From 'table', plot distrbution of ranks per species and then per group
ranks <- data.frame(table[table$group != "Other_phyto",] %>% group_by(group,species,predictor) %>% summarize(rank = mean(rank,na.rm=T)) )
dim(ranks)
head(ranks)

all.ranks <- data.frame(table[table$group != "Other_phyto",] %>% group_by(predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
all.ranks[order(all.ranks$avg.rank, decreasing = T),]
# Fix this order for the predictors in the y axis

plot <- ggplot(data = ranks, aes(x = rank, y = reorder(factor(predictor), rank) )) + geom_boxplot(aes(fill = factor(group)),colour = "black") + 
    scale_fill_brewer(palette = "Paired", name = "") + xlab("Mean covariate importance (% in permutations)") + ylab("") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + facet_wrap(factor(group)~., scales = "fixed")

setwd(WD)
ggsave(plot = plot, filename = "plot_vars_ranks_phyto_groups.jpg", dpi = 300, width = 12, height = 3)
    
    
### Tally numbers 
ranks <- data.frame(table[table$group != "Other_phyto",] %>% group_by(predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
ranks[order(ranks$avg.rank, decreasing = T),]

ranks <- data.frame(table[table$group != "Other_phyto",] %>% group_by(group,predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
r <- ranks[ranks$group == "Diatoms",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Dinoflagellates",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Coccolithophores",]; r[order(r$avg.rank, decreasing = T),]

ranks <- data.frame(table[table$group != "Other_phyto",] %>% group_by(SDM,predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
r <- ranks[ranks$SDM == "GLM",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "GAM",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "ANN",]; r[order(r$avg.rank, decreasing = T),]

ranks <- data.frame(table[table$group != "Other_phyto",] %>% group_by(group,SDM,predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
r <- ranks[ranks$SDM == "GLM" & ranks$group == "Coccolithophores",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "GAM" & ranks$group == "Coccolithophores",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "ANN" & ranks$group == "Coccolithophores",]; r[order(r$avg.rank, decreasing = T),]



### For zooplankton groups
ranks <- data.frame(table[table$group != "Other_zoo",] %>% group_by(predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
ranks[order(ranks$avg.rank, decreasing = T),]

ranks <- data.frame(table[table$group != "Other_zoo",] %>% group_by(group,predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
r <- ranks[ranks$group == "Calanoida",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Poecilostomatoida",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Oithonida",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Jellyfish",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Chaetognatha",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Euphausiids",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Foraminifera",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Amphipods",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Pteropods",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Salps",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$group == "Appendicularians",]; r[order(r$avg.rank, decreasing = T),]


ranks <- data.frame(table[table$group != "Other_zoo",] %>% group_by(SDM,predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
r <- ranks[ranks$SDM == "GLM",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "GAM",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "ANN",]; r[order(r$avg.rank, decreasing = T),]


unique(table$group)
v <- c("Calanoida","Poecilostomatoida","Oithonida","Jellyfish","Chaetognatha",
       "Euphausiids","Foraminifera","Amphipods","Pteropods","Salps","Appendicularians")
        
ranks <- data.frame(table[table$group != "Other_zoo",] %>% group_by(group,SDM,predictor) %>% summarize(avg.rank = mean(rank, na.rm = T)))
r <- ranks[ranks$SDM == "GLM" & ranks$group == "Appendicularians",] ; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "GAM" & ranks$group == "Appendicularians",]; r[order(r$avg.rank, decreasing = T),]
r <- ranks[ranks$SDM == "ANN" & ranks$group == "Appendicularians",]; r[order(r$avg.rank, decreasing = T),]


    
### Same with TSS scores: 
cat <- "Zooplankton"
set <- "Set2"
strategy <- "group"

if(cat == "Phytoplankton" & strategy == "total") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background/Scores")
    data.wd.back <- getwd()
} else if(cat == "Zooplankton" & strategy == "total") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background/Scores")
    data.wd.back <- getwd()
} else if(cat == "Phytoplankton" & strategy == "group") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background/Scores")
    data.wd.back <- getwd()
} else if(cat == "Zooplankton" & strategy == "group") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background/Scores")
    data.wd.back <- getwd()
} 
setwd(data.wd.back)
files <- dir()[grep(set,dir())] ; files

# f <- files[2]
res <- lapply(files, function(f) { data <- get(load(f)) ; return(data) } ) # eo lapply
table.scores <- do.call(rbind,res) ; gc() ; rm(res)
# Return SDM and species names 
terms <- data.frame(do.call(rbind,strsplit(as.character(rownames(table.scores)), split = "_")))
colnames(terms) <- c("SDM","RUN")
table.scores$SDM <- data.frame(do.call(rbind,strsplit(as.character(rownames(table.scores)), split = "_")))[,1]
table.scores$RUN <- data.frame(do.call(rbind,strsplit(as.character(rownames(table.scores)), split = "_")))[,2]
head(table.scores) ; dim(table.scores) ; summary(table.scores)

ranks <- data.frame(table.scores %>% group_by(group,SDM) %>% summarize(avg.rank = mean(TSS, na.rm = T)))
#ranks[order(ranks$avg.rank, decreasing = TRUE),]
ranks

ggplot(data = table.scores, aes(x = factor(SDM), y = TSS, fill = factor(SDM))) + geom_boxplot(colour = "black") + 
    scale_fill_brewer(palette = "Spectral", name = "") + xlab("") + ylab("TSS") + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data = table.scores, aes(x = factor(SDM), y = TSS, fill = factor(SDM))) + geom_boxplot(colour = "black") + 
    scale_fill_brewer(palette = "Spectral", name = "") + xlab("") + ylab("TSS") + 
    theme_bw() + facet_grid( ~ factor(group)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



### 17/09/2020: CONCLUSIONS FROM THE EXPERIMENTS ABOVE: 
# - no diff in models skills between Set1 and Set2 (keep one only)
# - no diff in models skills between total-background and group-background (test both for patterns)
# - clear variability in top ranking variables depending on PFT (for both phyto- and zoo-), so need to define group-level sets
#   (for total-background and group-background, Set2 only)

### Diatoms
# total-background = SST	Nstar	logChla	pCO2	PAR	logSiO2
# group-background = SST	Nstar	logSiO2	Sistar	logChla	PAR

### Dinoflagellates
# total-background = SST	Nstar	Wind	PAR	MLD	logChla
# group-background = SST	Nstar	Wind	PAR	MLD	logChla

### Coccolithophores
# total-background = Nstar	SST	Wind	logChla	Sistar	pCO2
# group-background = Nstar	SST	Wind	PAR     pCO2	MLD

### Calanoida
# total-background = SST	logSiO2	Sistar dO2	MLD	PAR
# group-background = SST	logSiO2	Sistar dO2	MLD	PAR

### Poecilostomatoida
# total-background = SST	Sistar	logChla	PAR	MLD	logSiO2
# group-background = SST	Sistar	MLD	PAR	logSiO2	Wind

### Oithonida
# total-background = SST	logSiO2	Sistar	dO2	MLD	logChla
# group-background = SST	logSiO2	Sistar	dO2	MLD	PAR

### Jellyfish
# total-background = SST	logSiO2	Sistar	dO2	MLD	logChla
# group-background = SST	Nstar	MLD	dO2	PAR	logSiO2

### Chaetognatha
# total-background = SST	Sistar	logSiO2	logEKE	MLD	PAR
# group-background = SST	Sistar	logEKE	logSiO2	MLD	PAR

### Euphausiids
# total-background = SST	Sistar	logSiO2	logEKE	MLD	PAR
# group-background = SST	Sistar	logSiO2 dO2	    MLD	Nstar

### Foraminifera
# total-background = SST	dO2	PAR	logSiO2	Wind	Sistar
# group-background = SST	dO2	PAR	logSiO2	Wind	Sistar

### Amphipods
# total-background = SST	logSiO2	Sistar	PAR	dO2	MLD
# group-background = SST	logSiO2	Sistar	dO2	PAR	logChla

### Pteropods
# total-background = SST	dO2	Sistar	logSiO2	PAR	logEKE
# group-background = dO2	SST	logSiO2	Sistar	PAR	MLD

### Salps
# total-background = SST	dO2	Sistar	logSiO2	PAR	logEKE
# group-background = dO2	SST	PAR	Wind	Sistar	MLD

### Appendicularians
# total-background = SST	dO2	Sistar	logSiO2	PAR	logEKE
# group-background = SST	Sistar	MLD	Wind	PAR	logSiO2


