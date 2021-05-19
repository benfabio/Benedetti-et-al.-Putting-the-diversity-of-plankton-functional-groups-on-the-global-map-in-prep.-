
##### 27/10/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Loading the annual estimates of the plankton groups' species richness and their annual covariates
#   - train GBM (boosted regression trees) and Random Forest models to model re-predict annual SR from the covariates
#   - evaluate sensitivity to various parameters (interaction depth, learning rate, min obs per tree, nb trees, nb of params to separet trees)
#   - randomly split the training data set 25/30 times into a 75/25 training/testing set
#   - evaluate RMSE/r2 for each model and keep the bets ones
#   - use iBreakDown's break_down() to evaluate the predictors contributions to the models (like in Busseni et al. 2020)
#   - summarize the information for eahc grid cell and evaluate regional/local covariates of annual SR for each group
 
### Last update: 01/11/2020 

# --------------------------------------------------------------------------------------------------------------------------------

# install.packages("iBreakDown")
library("tidyverse")
library("reshape2")
library("viridis")
library("scales")
library("maps")
library("parallel")
library("DALEX")
library("iBreakDown")
library("gbm")
library("ranger")
library("randomForest")

WD <- getwd()
world2 <- map_data("world2")

RMSE <- function(m,o) { sqrt(mean((m - o)^2)) } # m = predicted values by the model ; # o = observed values used to train the model

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Get the table with the annual baseline SR estimates based on the full predictors 
setwd(paste(WD,"/","biology/data_for_group_studies", sep = "")) ; dir()
table <- read.table("table_ann_rich_groups+env_baseline_27_10_20.txt", sep = "\t", h = T)
#dim(table) ; head(table)
summary(table)
colnames(table)
groups <- colnames(table)[c(4:17)] ; groups

### Reminder of the predictors used to model the groups' SR: 
# Environmental variables as a function of the group
# - Diatoms:            c("SST","Nstar","logChl","Sistar","PAR","logSiO2")
# - Dinos:              c("SST","Nstar","Wind","PAR","MLD","logChl")
# - Coccos:             c("SST","Nstar","Wind","PAR","MLD","pCO2")

### Vars to use for Phytoplankton: c("SST","Nstar","Wind","PAR","MLD","pCO2","logSiO2","logChl","Sistar")

# - Calanoida:          c("SST","logSiO2","Sistar","dO2","MLD","PAR")
# - Poecilostomatoida:  c("SST","Sistar","Wind","PAR","MLD","logSiO2")
# - Oithonida:          c("SST","logSiO2","Sistar","dO2","MLD","PAR")
# - Jellyfish:          c("SST","Nstar","MLD","PAR","dO2","logSiO2")
# - Chaetognatha:       c("SST","Sistar","logSiO2","logEKE","MLD","PAR")
# - Euphausiids:        c("SST","Sistar","logSiO2","dO2","MLD","Nstar")
# - Foraminifera:       c("SST","dO2","PAR","logSiO2","Wind","Sistar")
# - Amphipods:          c("SST","logSiO2","Sistar","PAR","dO2","logChl")
# - Pteropods:          c("SST","dO2","Sistar","logSiO2","PAR","MLD")
# - Salps:              c("SST","dO2","Sistar","logSiO2","PAR","MLD")
# - Appendicularians:   c("SST","Sistar","MLD","Wind","PAR","logSiO2")

### Vars to use for Zooplankton: c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")


### Implement some tests before defining the function to be ran in parallel
# g <- "Coccolithophores"
#
# if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
#     vars <- c("SST","pCO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
# } else {
#     vars <- c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
# } # eo if loop
#
#
# subset <- na.omit(table[,c("id","x","y",g,vars)])
# dim(subset)
# n <- nrow(subset)
#
# ### Split susbet into XX% training set and XX% testing set randomly
# ind <- sample(c(TRUE,FALSE), size = n, replace = T, prob = c(0.5,0.5))
# train <- subset[ind,] ; test <- subset[!ind,]
# # dim(train) ; dim(test) # Looks good
# # head(train) ; head(test)
#
# ### Derive formula
# form <- as.formula(paste(g, paste(vars, collapse = " + "), sep = " ~ "))
# # form
#
# ### Test one GBM
# # ?gbm
# # ?gbm.perf
# gbm.test <- gbm(formula = form, data = train, distribution = "gaussian", n.minobsinnode = 10, interaction.depth = 3,
#             shrinkage = 0.1, cv.folds = 0, verbose = T, n.cores = 5, n.trees = 100)
#
# summary.gbm(gbm.test)
# #gbm.perf(gbm.test, method = "cv")
#
# ### Predict SR fpr the testing set and derive RMSE
# test$fit <- predict(object = gbm.test, newdata = test[,vars], n.trees = 100)
# summary(test)
#
# ### Compute RMSE (and correlation between fitted and observed values)
# rmse <- round(RMSE(m = test$fit, o = test[,4]),2) ; rmse
# cor.spear <- round(cor(test$fit, test[,4], method = "spearman"),3) ; cor.spear
#
# # And a plot for looking at pred vs obs
# ggplot(data = test, aes(x = get(g), y = fit)) + geom_point(colour = "grey65", alpha = .5) +
#     geom_abline(slop = 1, intercept = 0, linetype = "dashed") +
#     #geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+1.5), label = paste("RMSE = ",rmse, sep = "")) +
#     #geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+0.5), label = paste("n = ", nrow(test), sep = "")) +
#     #geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+1), label = paste("rho = ", cor.spear, sep = "")) +
#     xlab(paste("Observed mean annual SR ","(",g,")", sep = "")) + ylab("Predicted mean annual SR from GBM") +
#     theme_classic() + ggtitle(paste("RMSE = ",rmse," ; rho = ",cor.spear, sep = ""))
#
# # --------------------------------------------------------------
#
# ### Test one RF with ranger
# # ?ranger
# rf.test <- ranger(formula = form, data = train, num.trees = 250, mtry = 3, importance = 'none',
#             write.forest = T, min.node.size = 50, verbose = T, save.memory = F, classification = F)
# #rf.test
# ### Predict SR fpr the testing set and derive RMSE
# p <- predict(object = rf.test, data = test[,vars], type = "response")
# test$fit <- p$predictions
# #summary(test)
# ### Compute RMSE (and correlation between fitted and observed values)
# rmse <- round(RMSE(m = test$fit, o = test[,4]),2) ; rmse
# cor.spear <- round(cor(test$fit, test[,4], method = "spearman"),3) ; cor.spear
#
# ### And a plot for looking at pred vs obs
# ggplot(data = test, aes(x = get(g), y = fit)) + geom_point(colour = "grey65", alpha = .5) +
#     geom_abline(slop = 1, intercept = 0, linetype = "dashed") +
#     #geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+1.5), label = paste("RMSE = ",rmse, sep = "")) +
#     #geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+0.5), label = paste("n = ", nrow(test), sep = "")) +
#     #geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+1), label = paste("rho = ", cor.spear, sep = "")) +
#     xlab(paste("Observed mean annual SR ","(",g,")", sep = "")) + ylab("Predicted mean annual SR from RF") +
#     theme_classic() + ggtitle(paste("RMSE = ",rmse," ; rho = ",cor.spear, sep = ""))


# --------------------------------------------------------------

### 28/10/2020: This seems to work nicely :) Now check the result of explain and break_down for some random cells
#  expl_gbm <- explain(gbm.test, data = train[,c(g,vars)], y = train[,g])
#  expl_rf <- explain(rf.test, data = train[,c(g,vars)], y = train[,g])
#bd_gbm <- break_down(x = expl_gbm, new_observation = test[1:100,c(g,vars)], keep_distributions = F, interactions = F)
#bd_rf <- break_down(x = expl_rf, new_observation = test[1:100,c(g,vars)], keep_distributions = F, interactions = F)
# Takes some time though
# plot(bd_gbm)
# plot(bd_rf)
#bd_gbm$variable ; bd_gbm$contribution
#bd_rf$variable ; bd_rf$contribution

# (abs(bd_gbm$contribution[2:10]) / sum(abs(bd_gbm$contribution[2:10])))*100
# sum(bd_gbm$contribution[1:10]) ; bd_gbm$contribution[12]

### To get basins codes and look at features importace regionally (https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18documentation.pdf)
# setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors")
# basins <- read.table("basinmask_01.msk", h = T, sep = ",")[,1:3]
# colnames(basins) <- c("y","x","basin")
# head(basins) ; dim(basins)
# summary(factor(basins$basin))
# basins$x2 <- basins$x
# basins[basins$x < 0 ,"x2"] <- (basins[basins$x < 0 ,"x"]) + 360
# basins$id <- factor(paste(basins$x2, basins$y, sep = "_"))
# summary(basins)
# # Translation fo the codes:
# # 1 = Atlantic
# # 2 = Pacific
# # 3 = Indian
# # 4 = Med Sea
# # 10 = Southern Ocean
# # 11 = Arctic Ocean
# ids2keep <- basins[basins$basin %in% c(1:4,10,11),"id"]
# commons <- intersect(unique(table$id), unique(ids2keep)) ; length(commons)
# basins2 <- basins[basins$id %in% commons,]
# table2 <- table[table$id %in% commons,]
# basins2 <- basins2[order(basins2$id),]
# table2 <- table2[order(table2$id),]
# dim(basins2) ; dim(table2)
# table2$basin <- basins2$basin
# head(table2) ; dim(table2)
#
# bd_gbm <- break_down(x = expl_gbm, new_observation = table2[which(table2$basin == 3),c(vars)], keep_distributions = F, interactions = F)
# bd_rf <- break_down(x = expl_rf, new_observation = table2[which(table2$basin == 3),c(vars)], keep_distributions = F, interactions = F)
#
# plot(bd_gbm)
# plot(bd_rf)


# --------------------------------------------------------------------------------------------------------------------------------

### 30/10/2020: Next, define a suite of range for the parameters of GBMs and RF to test (n trees, n.minobsinnode, shrinkage, mtry etc.)
### And train models for all combiantions posssible and for each group. Find the combiantion with the lowest RMSE. Then, train full model
### based on the "best" parameters suite for each group and examine features importance globally, or per main basin

### Range of parameters to define in GBMs:
# n.minobsinnode (Integer specifying the minimum number of observations in the terminal nodes of the trees) --> 5,10,25,50
# interaction.depth (Integer specifying the maximum depth of each tree (i.e. the highest level of variable interactions allowed. value == 1 implies an additive model) --> 1,3,5
# shrinkage (learning rate or step-size reduction): usually between 0.01 and 0.1 --> 0.01 & 0.1
# n.trees (This is equivalent to the number of iterations and the number of basis functions in the additive expansion) --> 50, 75, 100, 500, 750
### according to LG: interaction.depth & n.trees are the 2 main parameters controlling model complexity and overfitting, so 

min_obs_node <- c(5,10,25,50)
inter_depth <- c(1,3,5)
shrinkage <- c(0.01,0.1)
Ntree_gbm <- c(50,75,100,500,750)
GBM_params <- expand.grid(min_obs_node, inter_depth, shrinkage, Ntree_gbm)
colnames(GBM_params) <- c("min_obs_node","inter_depth","shrinkage","Ntree_gbm")
GBM_params$combination <- apply(expand.grid(min_obs_node, inter_depth, shrinkage, Ntree_gbm), 1, paste, collapse = "_")
GBM_params

### Range of parameters to define in Random Forests:
# num.trees --> 100, 500, 750, 1000, 2000
# mtry (# of variables to possibly split at in each node) --> 1,3,5
# min.node.size (default = 5 for regression) --> 5, 10, 30, 50 
# max.depth (Maximal tree depth. A value of NULL or 0 (the default) corresponds to unlimited depth --> 0 ; keep default
Ntree_rf <- c(100,500,750,1000,2000)
mtry <- c(1,3,5)
min.node <- c(5,10,30,50)
RF_params <- expand.grid(Ntree_rf, mtry, min.node)
colnames(RF_params) <- c("Ntree_rf","mtry","min.node")
RF_params$combination <- apply(expand.grid(Ntree_rf, mtry, min.node), 1, paste, collapse = "_")
RF_params

### For each of the 14 groups, split the global dataset into a 50/50 training/testing set 20 times and test all parameters combinations. Retrieve rmse and cor.spear. Check distributions of RMSE per parameters and groups and make final choice.
splits <- c(1:25)

### In a mclapply, create 25 random 50/50 splits and test each combination for each group and return a ddf with rmse, rho and parameters
# g <- "Diatoms"
# s <- 3
# i <- 15 # GBM_params[i,]

res.gbms <- mclapply(splits, function(s) {
            
            res.groups <- lapply(groups, function(g) {
                    
                    ### Testing GBM parameters for group g , split s
                    message(paste("Testing GBM parameters for ",g," || split #",s, sep = ""))
                    message(paste("", sep = ""))
                    
                    if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
                        vars <- c("SST","pCO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
                    } else {
                        vars <- c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
                    } # eo if loop 
                    
                    # Split the global susbet into XX% training set and XX% testing set randomly
                    subset <- na.omit(table[,c("id","x","y",g,vars)])
                    n <- nrow(subset)
                    ind <- sample(c(TRUE,FALSE), size = n, replace = T, prob = c(0.5,0.5))
                    train <- subset[ind,] ; test <- subset[!ind,]
                    # Derive formula
                    form <- as.formula(paste(g, paste(vars, collapse = " + "), sep = " ~ "))
                    
                    # And yet again, another lapply where you check every line of GBM_params
                    res.params <- lapply(c(1:nrow(GBM_params)), function(i) {
                        
                            # Perform GBM
                            message(paste("Running GBM based on parameters from row ", i, sep = ""))
                             
                            gbm.test <- gbm(formula = form, data = train, distribution = "gaussian",
                                        n.minobsinnode = GBM_params[i,"min_obs_node"], interaction.depth = GBM_params[i,"inter_depth"],
                                        n.trees = GBM_params[i,"Ntree_gbm"], shrinkage = GBM_params[i,"shrinkage"],
                                        cv.folds = 0, verbose = F, n.cores = 1)

                            # Predict SR fpr the testing set and derive RMSE
                            test$fit <- predict(object = gbm.test, newdata = test[,vars], n.trees = GBM_params[i,"Ntree_gbm"])

                            # Compute RMSE (and correlation between fitted and observed values)
                            rmse <- round(RMSE(m = test$fit, o = test[,4]),3)
                            cor.spear <- round(cor(test$fit, test[,4], method = "spearman"),3)
                            
                            # Return a table with all the info and results
                            table.res <- GBM_params[i,]
                            table.res$group <- g
                            table.res$RMSE <- rmse
                            table.res$rho <- cor.spear
                            table.res$slit <- s
                            
                            return(table.res)
                        
                        } # eo FUN - i in nrows(GBM_params)
                    
                    ) # eo lapply - res.params
                    # Rbind
                    table.params <- bind_rows(res.params)
                    rm(res.params) ; gc()
                    
                    return(table.params)  
                
                } # eo FUN g in groups 
                
            ) # eo lapply - res.groups
            # Rbind and return 
            table.group <- bind_rows(res.groups)
            rm(res.groups) ; gc()
            
            return(table.group)
    
    }, mc.cores = length(splits) 
    
    
) # eo mclapply
### Rbind
table_gbms <- bind_rows(res.gbms)
head(table_gbms) ; dim(table_gbms) ; summary(table_gbms)
# Make room
rm(res.gbms) ; gc()

save(table_gbms, file = "table_tests_parameters_GBMs_30_10_2020.Rdata")

### Same as above but for random forests
res.rfs <- mclapply(splits, function(s) {
            
            res.groups <- lapply(groups, function(g) {
                    
                    ### Testing GBM parameters for group g , split s
                    message(paste("", sep = ""))
                    message(paste("Testing RF parameters for ",g," || split #",s, sep = ""))
                    message(paste("", sep = ""))
                    
                    if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
                        vars <- c("SST","pCO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
                    } else {
                        vars <- c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
                    } # eo if loop 
                    
                    # Split the global susbet into XX% training set and XX% testing set randomly
                    subset <- na.omit(table[,c("id","x","y",g,vars)])
                    n <- nrow(subset)
                    ind <- sample(c(TRUE,FALSE), size = n, replace = T, prob = c(0.5,0.5))
                    train <- subset[ind,] ; test <- subset[!ind,]
                    # Derive formula
                    form <- as.formula(paste(g, paste(vars, collapse = " + "), sep = " ~ "))
                    
                    # And yet again, another lapply where you check every line of GBM_params
                    res.params <- lapply(c(1:nrow(RF_params)), function(i) {
                        
                            # Perform random forest
                            message(paste("Running RF based on parameters from row ", i, sep = ""))
                             
                            rf.test <- ranger(formula = form, data = train, num.trees = RF_params[i,"Ntree_rf"],
                                        mtry = RF_params[i,"mtry"], min.node.size = RF_params[i,"min.node"],
                                        write.forest = T, importance = 'none', verbose = F, save.memory = F, 
                                        classification = F)
                                        
                            # Predict SR fpr the testing set and derive RMSE
                            p <- predict(object = rf.test, data = test[,vars], type = "response")
                            test$fit <- p$predictions
                            rmse <- round(RMSE(m = test$fit, o = test[,4]),2) ; rmse
                            cor.spear <- round(cor(test$fit, test[,4], method = "spearman"),3) ; cor.spear

                            # Return a table with all the info and results
                            table.res <- RF_params[i,]
                            table.res$group <- g
                            table.res$RMSE <- rmse
                            table.res$rho <- cor.spear
                            table.res$slit <- s
                            
                            return(table.res)
                        
                        } # eo FUN - i in nrows(GBM_params)
                    
                    ) # eo lapply - res.params
                    # Rbind
                    table.params <- bind_rows(res.params)
                    rm(res.params) ; gc()
                    
                    return(table.params)  
                
                } # eo FUN g in groups 
                
            ) # eo lapply - res.groups
            # Rbind and return 
            table.group <- bind_rows(res.groups)
            rm(res.groups) ; gc()
            
            return(table.group)
    
    }, mc.cores = length(splits) 
    
) # eo mclapply
### Rbind 
table_rf <- bind_rows(res.rfs)
head(table_rf) ; dim(table_rf) ; summary(table_rf)
# Make room
rm(res.rfs) ; gc()

save(table_rf, file = "table_tests_parameters_RFs_30_10_2020.Rdata")

rm(table_rf, table_gbms) ; gc()

# --------------------------------------------------------------------------------------------------------------------------------

### 01/11/2020: Plot results from above:
# - plot RMSE per group and per model type (GBM vs. RF)
# - plot RMSE per parameters (n trees etc.)

### A°) GBM skills across parameters combinations    -----------------------------------------------------------------------------

table <- get(load("table_tests_parameters_GBMs_30_10_2020.Rdata"))
head(table) ; dim(table)

groups <- data.frame(table %>% group_by(group) %>% summarize(RMSE = mean(RMSE, na.rm = T)))
groups[order(groups$RMSE, decreasing = T),]

groups <- data.frame(table %>% group_by(group) %>% summarize(rho = mean(rho, na.rm = T)))
groups[order(groups$rho, decreasing = T),]

ggplot(aes(x = factor(group), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("Group") + ylab("RMSE (GBM models)") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#
ggplot(aes(x = factor(group), y = rho), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("Group") + ylab("rho (GBM models)") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Some plankton groups display high variance in RMSE/rho (e.g. Calanoida has higher variance in RMSE than Coccos or Amphipods)
### This means that not all parameters combinations are equally at predicting 50% of the global richness estimates. So it's worth
### examining the various parameters separately. It also means some groups hare harder to predict than other.

### A.1) RMSE across min_obs_node
ggplot(aes(x = factor(min_obs_node), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("min_obs_node") + ylab("RMSE (GBM models)") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# + Facet per group
ggplot(aes(x = factor(min_obs_node), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("min_obs_node") + ylab("RMSE (GBM models)") + theme_classic() +
    facet_wrap(~ factor(group))

### --> RMSE doesn't really change across min_obs_node...but across groups
groups <- data.frame(table %>% group_by(min_obs_node) %>% summarize(rmse = mean(RMSE, na.rm = T), sd = sd(RMSE, na.rm = T)))
groups[order(groups$rmse, decreasing = T),]
groups <- data.frame(table %>% group_by(min_obs_node) %>% summarize(r = mean(rho, na.rm = T), sd = sd(rho, na.rm = T)))
groups[order(groups$r, decreasing = T),]
### no impact

### A.2) RMSE across inter_depth
ggplot(aes(x = factor(inter_depth), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("inter_depth") + ylab("RMSE (GBM models)") + theme_classic() +
    scale_y_continuous(limits = c(0,5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(aes(x = factor(inter_depth), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("inter_depth") + ylab("RMSE (GBM models)") + theme_classic() +
    facet_wrap(~ factor(group), scale = "free")

groups <- data.frame(table %>% group_by(inter_depth) %>% summarize(rmse = mean(RMSE, na.rm = T), sd = sd(RMSE, na.rm = T)))
groups[order(groups$rmse, decreasing = T),]
groups <- data.frame(table %>% group_by(inter_depth) %>% summarize(r = mean(rho, na.rm = T), sd = sd(rho, na.rm = T)))
groups[order(groups$r, decreasing = T),]

### --> RMSE does decrease with increasing inter_depth, especially for some groups (e.g. chaetognaths, coccos)
### Choose inter_depth = 5 for all groups


### A.3) RMSE across shrinkage
ggplot(aes(x = factor(shrinkage), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("shrinkage") + ylab("RMSE (GBM models)") + theme_classic() +
    scale_y_continuous(limits = c(0,5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(aes(x = factor(shrinkage), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("shrinkage") + ylab("RMSE (GBM models)") + theme_classic() +
    facet_wrap(~ factor(group), scale = "free")

groups <- data.frame(table %>% group_by(shrinkage) %>% summarize(rmse = mean(RMSE, na.rm = T), sd = sd(RMSE, na.rm = T)))
groups[order(groups$rmse, decreasing = T),]
groups <- data.frame(table %>% group_by(shrinkage) %>% summarize(r = mean(rho, na.rm = T), sd = sd(rho, na.rm = T)))
groups[order(groups$r, decreasing = T),]

### --> RMSE decreases as skrinkage increases to 0.1. For all groups. Choose shrinkage == 0.1


### A.4) RMSE across Ntree
ggplot(aes(x = factor(Ntree_gbm), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("N trees") + ylab("RMSE (GBM models)") + theme_classic() +
    scale_y_continuous(limits = c(0,5)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(aes(x = factor(Ntree_gbm), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("N trees") + ylab("RMSE (GBM models)") + theme_classic() +
    facet_wrap(~ factor(group), scale = "free")

groups <- data.frame(table %>% group_by(Ntree_gbm) %>% summarize(rmse = mean(RMSE, na.rm = T), sd = sd(RMSE, na.rm = T)))
groups[order(groups$rmse, decreasing = T),]
groups <- data.frame(table %>% group_by(Ntree_gbm) %>% summarize(r = mean(rho, na.rm = T), sd = sd(rho, na.rm = T)))
groups[order(groups$r, decreasing = T),]

### --> RMSE decreases with increasing N trees, for all groups. No sure it really changes after N tree > 500

### CONCLUSIONS: 
# - min_obs_node = whatever (use 50)
# - inter_depth = 5
# - shrinkage = 0.1
# - Ntree_gbm = 750


### B°) RF skills across parameters combinations    -----------------------------------------------------------------------------

table <- get(load("table_tests_parameters_RFs_30_10_2020.Rdata"))
head(table) ; dim(table)

groups <- data.frame(table %>% group_by(group) %>% summarize(RMSE = mean(RMSE, na.rm = T)))
groups[order(groups$RMSE, decreasing = T),]

groups <- data.frame(table %>% group_by(group) %>% summarize(rho = mean(rho, na.rm = T)))
groups[order(groups$rho, decreasing = T),]

ggplot(aes(x = factor(group), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("Group") + ylab("RMSE (GBM models)") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Again: plankton groups display high variance in RMSE/rho
### This means that not all parameters combinations are equally at predicting 50% of the global richness estimates. So it's worth
### examining the various parameters separately. It also means some groups hare harder to predict than other.


### B.1) RF across Ntree
ggplot(aes(x = factor(Ntree_rf), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("N trees") + ylab("RMSE (GBM models)") + theme_classic() +
    scale_y_continuous(limits = c(0,2)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(aes(x = factor(Ntree_rf), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("N trees") + ylab("RMSE (GBM models)") + theme_classic() +
    facet_wrap(~ factor(group), scale = "free")

groups <- data.frame(table %>% group_by(Ntree_rf) %>% summarize(rmse = mean(RMSE, na.rm = T), sd = sd(RMSE, na.rm = T)))
groups[order(groups$rmse, decreasing = T),]
groups <- data.frame(table %>% group_by(Ntree_rf) %>% summarize(r = mean(rho, na.rm = T), sd = sd(rho, na.rm = T)))
groups[order(groups$r, decreasing = T),]

### --> RMSE does not vary with Ntree_rf

### B.2) RF across mtry
ggplot(aes(x = factor(mtry), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("mtry") + ylab("RMSE (GBM models)") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(aes(x = factor(mtry), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("mtry") + ylab("RMSE (GBM models)") + theme_classic() +
    facet_wrap(~ factor(group), scale = "free")

groups <- data.frame(table %>% group_by(mtry) %>% summarize(rmse = mean(RMSE, na.rm = T), sd = sd(RMSE, na.rm = T)))
groups[order(groups$rmse, decreasing = T),]
groups <- data.frame(table %>% group_by(mtry) %>% summarize(r = mean(rho, na.rm = T), sd = sd(rho, na.rm = T)))
groups[order(groups$r, decreasing = T),]

### --> RMSE decreases with mtry going for 1 to 3


### B.3) RF across min.node
ggplot(aes(x = factor(min.node), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("min.node") + ylab("RMSE (GBM models)") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(aes(x = factor(min.node), y = RMSE), data = table) +
    geom_boxplot(fill = "grey65", colour = "black") +
    xlab("min.node") + ylab("RMSE (GBM models)") + theme_classic() +
    facet_wrap(~ factor(group), scale = "free")

groups <- data.frame(table %>% group_by(min.node) %>% summarize(rmse = mean(RMSE, na.rm = T), sd = sd(RMSE, na.rm = T)))
groups[order(groups$rmse, decreasing = T),]
groups <- data.frame(table %>% group_by(min.node) %>% summarize(r = mean(rho, na.rm = T), sd = sd(rho, na.rm = T)))
groups[order(groups$r, decreasing = T),]

### --> RMSE increases with min.node, choose min.node = 5


### Make some plots to check your final GBM parameters choice
# - min_obs_node = whatever (use 50)
# - inter_depth = 5
# - shrinkage = 0.1
# - Ntree_gbm = 750

g <- "Oithonida"

if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
    vars <- c("SST","pCO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
} else {
    vars <- c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar")
} # eo if loop 

# Split the global susbet into XX% training set and XX% testing set randomly
subset <- na.omit(table[,c("id","x","y",g,vars)])
n <- nrow(subset)
ind <- sample(c(TRUE,FALSE), size = n, replace = T, prob = c(0.5,0.5))
train <- subset[ind,] ; test <- subset[!ind,]
form <- as.formula(paste(g, paste(vars, collapse = " + "), sep = " ~ "))

gbm.test <- gbm(formula = form, data = train, distribution = "gaussian",
            n.minobsinnode = 50, interaction.depth = 5, n.trees = 500, shrinkage = 0.1,
            cv.folds = 0, verbose = T, n.cores = 1)

# Predict SR fpr the testing set and derive RMSE
test$fit <- predict(object = gbm.test, newdata = test[,vars], n.trees = 500)

# Compute RMSE (and correlation between fitted and observed values)
rmse <- round(RMSE(m = test$fit, o = test[,4]),3)
cor.spear <- round(cor(test$fit, test[,4], method = "spearman"),3)

ggplot(data = test, aes(x = get(g), y = fit)) + geom_point(colour = "grey65", alpha = .5) +
     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
     # geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+1.5), label = paste("RMSE = ",rmse, sep = "")) +
     # geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+0.5), label = paste("n = ", nrow(test), sep = "")) +
     # geom_text(aes(x = max(test[,g])-1.5, y = min(test$fit)+1), label = paste("rho = ", cor.spear, sep = "")) +
     xlab(paste("Observed mean annual SR ","(",g,")", sep = "")) + ylab("Predicted mean annual SR from optimal GBM") +
     theme_classic() + ggtitle(paste("RMSE = ",rmse," ; rho = ",cor.spear, sep = ""))

### Check full model
subset <- na.omit(table[,c("id","x","y",g,vars)])
form <- as.formula(paste(g, paste(vars, collapse = " + "), sep = " ~ "))
gbm.test <- gbm(formula = form, data = subset, distribution = "gaussian",
            n.minobsinnode = 50, interaction.depth = 5, n.trees = 500, shrinkage = 0.1,
            cv.folds = 0, verbose = T, n.cores = 1)
#
summary(gbm.test)

# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
