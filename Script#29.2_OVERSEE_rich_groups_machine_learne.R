
##### 27/10/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Loading the annual estimates of the plankton groups' species richness and their annual covariates
#   - Train GBM based on the selected set of GBM parameters defined in the previous script
#   - Draw some plots of CV skill with RMSE for each group
#   - Use iBreakDown() to assess regional (main basins) drivers of groups' species richness
 
### Last update: 11/11/2020 

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

setwd(paste(WD,"/","biology/data_for_group_studies", sep = "")) ; dir()
table <- read.table("table_ann_rich_groups+env_baseline_27_10_20.txt", sep = "\t", h = T)
groups <- colnames(table)[c(4:17)] ; groups
# g <- "Chaetognatha"

for(g in groups) {
    
    message(paste("Plotting CV skill for ",g, sep = ""))
    
    if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
        vars <- c("SST","pCO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar","logEKE")
    } else {
        vars <- c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar","logEKE")
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

    plot <- ggplot(data = test, aes(x = get(g), y = fit)) + geom_point(colour = "grey65", alpha = .5) +
         geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
         xlab(paste("Observed mean annual SR ","(",g,")", sep = "")) + ylab("Predicted mean annual SR from optimal GBM") +
         theme_classic() + ggtitle(paste("RMSE = ",rmse," ; rho = ",cor.spear, sep = ""))
    #
    ggsave(plot = plot, filename = paste("plot_cv_opt_GBM_",g,".jpg", sep=""), dpi = 300, width = 4, height = 4)
    
} # eo for loop - g in groups


# ---------------------------------------------------------------

### Load basins code and make some modifications
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors")
basins <- read.table("basinmask_01.msk", h = T, sep = ",")[,1:3]
colnames(basins) <- c("y","x","basin")
basins$x2 <- basins$x
basins[basins$x < 0 ,"x2"] <- (basins[basins$x < 0 ,"x"]) + 360
basins$id <- factor(paste(basins$x2, basins$y, sep = "_"))
basins <- basins[order(basins$id),]
summary(basins)
str(basins)
# Translation fo the codes:
# 1 = Atlantic Ocean
# 2 = Pacific Ocean
# 3 = Indian Ocean
# 4 = Med Sea
# 10 = Southern Ocean
# 11 = Arctic Ocean
# Re-level
basins[which(basins$basin == 1 & basins$y > 0),"basin"] <- "North Atlantic"
basins[which(basins$basin == 1 & basins$y < 0),"basin"] <- "South Atlantic"
basins[which(basins$basin == 2 & basins$y > 0),"basin"] <- "North Pacific"
basins[which(basins$basin == 2 & basins$y < 0),"basin"] <- "South Pacific"
basins[which(basins$basin == 3),"basin"] <- "Indian Ocean"
basins[which(basins$basin == 4),"basin"] <- "Mediterranean Sea"
basins[which(basins$basin == 10),"basin"] <- "Southern Ocean"
basins[which(basins$basin == 11),"basin"] <- "Arctic Ocean"
unique(basins$basin)
summary(factor(basins$basin))

ids2keep <- basins[basins$basin %in% c("North Atlantic","South Atlantic","North Pacific",
            "South Pacific","Mediterranean Sea","Southern Ocean","Arctic Ocean","Indian Ocean"),"id"]
# length(ids2keep)
basins <- basins[basins$id %in% ids2keep,]
dim(basins) ; head(basins)
# Quickmap
# ggplot() + geom_raster(aes(x = x2, y = y, fill = factor(basin)), data = basins) +
#      scale_fill_brewer(name = 'Basins', palette = "RdYlBu") +
#      geom_polygon(aes(x = long, y = lat, group = group), data = world2, fill = "grey70", colour = "black", size = 0.3) +
#      coord_map(projection = "moll") + theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
#      panel.grid.major = element_line(colour = "grey70",linetype = "dashed"), legend.position = "bottom") +
#      scale_x_continuous(name = "", limits = c(0,360), expand = c(0,0), labels = NULL) +
#      scale_y_continuous(name = "", limits = c(-90,90), expand = c(0,0), labels = NULL)

# ---------------------------------------------------------------

### For all plankton groups, run GBM with optimal parameters and derive regional features importance
g <- "Diatoms"
require("parallel")

res.groups <- mclapply(groups, function(g) {
            
            message(paste("", sep = ""))
            message(paste("Deriving features importance for ", g, sep = ""))
            
            
            if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
                vars <- c("SST","pCO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar","logEKE")
            } else {
                vars <- c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar","logEKE")
            } # eo if loop 

            # Split the global susbet into XX% training set and XX% testing set randomly
            subset <- na.omit(table[,c("id","x","y",g,vars)])
            n <- nrow(subset)

            form <- as.formula(paste(g, paste(vars, collapse = " + "), sep = " ~ "))

            gbm.mod <- gbm(formula = form, data = subset, distribution = "gaussian",
                        n.minobsinnode = 50, interaction.depth = 5, n.trees = 500, shrinkage = 0.1,
                        cv.folds = 0, verbose = F, n.cores = 1)        
            
            # Convert to an explainer object
            expl_gbm <- explain(gbm.mod, data = subset[,c(g,vars)], y = subset[,g])
            
            # Now, for each basin, use break_down based on the "expl_gbm" to extract regional features importance
            # b <- "Southern Ocean"
            res.basins <- lapply(unique(basins$basin), function(b) {
                
                        # Select cells belonging to basin
                        basin2keep <- unique(basins[which(basins$basin == b),"id"])
                        commons <- intersect(unique(subset$id), unique(basin2keep))
                        basins2 <- basins[basins$id %in% commons,]
                        table2 <- table[table$id %in% commons,]
                        basins2 <- basins2[order(basins2$id),]
                        table2 <- table2[order(table2$id),]
                        table2$basin <- basins2$basin
                            
                        # perform break_down
                        message(paste(b, sep = ""))
                        bd_gbm <- break_down(x = expl_gbm, new_observation = table2[,vars], keep_distributions = F, interactions = F)
                        # Extract features relative importance (name, importance, sign)
                        # str(bd_gbm)
                        # bd_gbm$variable_name
                        # bd_gbm$contribution  # dont retrieve last one though
                        contribs <- data.frame(variable = bd_gbm$variable_name, sign = bd_gbm$sign, contrib = bd_gbm$contribution)
                        # Get rid of last row
                        contribs <- contribs[2:10,]
                        contribs$relative <- (abs(contribs$contrib) / sum(abs(contribs$contrib)))*100
                        # Add group and basin
                        contribs$basin <- b
                        contribs$group <- g
                        # Return
                        rm(table2,basins2,basin2keep) ; gc()
                        return(contribs)
                        
                } # eo b in basins   
            
            ) # eo lapply
            # Rbind and return
            ddf.basins <- bind_rows(res.basins)            
            rm(res.basins) ; gc()
            # Change some levels in vars names
            levels(ddf.basins$variable)[levels(ddf.basins$variable) == "Sistar"] <- "Si*"
            levels(ddf.basins$variable)[levels(ddf.basins$variable) == "Nstar"] <- "N*"
            levels(ddf.basins$variable)[levels(ddf.basins$variable) == "logSiO2"] <- "Silicates"
            levels(ddf.basins$variable)[levels(ddf.basins$variable) == "logChl"] <- "Chlorophyll"
            levels(ddf.basins$variable)[levels(ddf.basins$variable) == "dO2"] <- "Oxygen (200m)"
            levels(ddf.basins$variable)[levels(ddf.basins$variable) == "logEKE"] <- "EKE"
            
             # ggplot(ddf.basins, aes(x = 2, y = relative, fill = factor(variable))) +
#                  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
#                  scale_fill_brewer(name = "", palette = "Paired") + theme_void() +
#                  xlim(0.5,2.5) + facet_wrap(~ factor(basin)) + theme(legend.position = "bottom")
            
            return(ddf.basins)
    
    }, mc.cores = length(groups)
     
) # eo mclapply - g in groups 

# Rbind and plot results
ddf <- bind_rows(res.groups)
head(ddf) ; dim(ddf) ; summary(ddf)
rm(res.groups) ; gc()

# Plots (select sensible colours for variables)
cols <- c("SST" = "#e31a1c", "Oxygen (200m)" = "#1f78b4", "pCO2" = "#ffff99", "Wind" = "#a6cee3",
        "PAR" = "#ff7f00", "MLD" = "#d9d9d9", "N*" = "#6a3d9a", "Si*" = "#cab2d6",
        "Chlorophyll" = "#33a02c", "Silicates" = "#b2df8a", "EKE" = "#fb9a99")

ggplot(data = ddf[ddf$group == "Chaetognatha",], aes(x = 2, y = relative, fill = factor(variable))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_manual(name = "", values = cols) + theme_void() +
    xlim(0.5,2.5) + facet_wrap(~ factor(basin)) + theme(legend.position = "bottom")

### OR: compute average relative contribution across regions and plot per group as facet
ddf2 <- data.frame(ddf %>% group_by(group, variable) %>% summarize(mean = mean(relative), sd = sd(relative)))
# ddf2
plot <- ggplot(data = ddf2, aes(x = 2, y = mean, fill = factor(variable))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_manual(name = "", values = cols) + theme_void() +
    xlim(0.5,2.5) + facet_wrap(~ factor(group)) + theme(legend.position = "bottom")
#
setwd(WD)
ggsave(plot = plot, filename = "donuts_mean_contrib_vars_groups_opt_GBM.pdf", dpi = 300, width = 10, height = 15)

for(g in groups) {
    
    plot <- ggplot(data = ddf[ddf$group == g,], aes(x = 2, y = relative, fill = factor(variable))) +
        geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
        scale_fill_manual(name = "", values = cols) + theme_void() +
        xlim(0.5,2.5) + facet_wrap(~ factor(basin)) + theme(legend.position = "bottom")

    ggsave(plot = plot, filename = paste("donuts_mean_contrib_vars_opt_GBM_",g,".pdf", sep = ""), dpi = 300, width = 6, height = 13)
    
}


### 04/11/2020: As above but globally for each group
g <- "Coccolithophores"

res.groups <- mclapply(groups, function(g) {
            
            message(paste("", sep = ""))
            message(paste("Deriving global features importance for ", g, sep = ""))
            
            
            if(g %in% c("Coccolithophores","Diatoms","Dinoflagellates")) {
                vars <- c("SST","pCO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar","logEKE")
            } else {
                vars <- c("SST","dO2","Wind","PAR","MLD","Nstar","logSiO2","logChl","Sistar","logEKE")
            } # eo if loop 

            # Split the global susbet into XX% training set and XX% testing set randomly
            subset <- na.omit(table[,c("id","x","y",g,vars)])
            n <- nrow(subset)

            form <- as.formula(paste(g, paste(vars, collapse = " + "), sep = " ~ "))

            gbm.mod <- gbm(formula = form, data = subset, distribution = "gaussian",
                        n.minobsinnode = 50, interaction.depth = 5, n.trees = 500, shrinkage = 0.1,
                        cv.folds = 0, verbose = F, n.cores = 1)        
            
            # Convert to an explainer object
            expl_gbm <- explain(gbm.mod, data = subset[,c(g,vars)], y = subset[,g])
            
            # Perform break_down
            bd_gbm <- break_down(x = expl_gbm, new_observation = subset[,vars], keep_distributions = F, interactions = F)
            # Extract features relative importance (name, importance, sign)
            contribs <- data.frame(variable = bd_gbm$variable_name, sign = bd_gbm$sign, contrib = bd_gbm$contribution)
            # Get rid of last row
            contribs <- contribs[2:10,]
            contribs$relative <- (abs(contribs$contrib) / sum(abs(contribs$contrib)))*100
            # Add group and basin
            contribs$group <- g
            
            # Change some levels in vars names
            levels(contribs$variable)[levels(contribs$variable) == "Sistar"] <- "Si*"
            levels(contribs$variable)[levels(contribs$variable) == "Nstar"] <- "N*"
            levels(contribs$variable)[levels(contribs$variable) == "logSiO2"] <- "Silicates"
            levels(contribs$variable)[levels(contribs$variable) == "logChl"] <- "Chlorophyll"
            levels(contribs$variable)[levels(contribs$variable) == "dO2"] <- "Oxygen (200m)"
            levels(contribs$variable)[levels(contribs$variable) == "logEKE"] <- "EKE"
            
              #ggplot(contribs, aes(x = 2, y = relative, fill = factor(variable))) +
                 #  geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
                  # scale_fill_brewer(name = "", palette = "Paired") + theme_void() +
                  # xlim(0.5,2.5) + theme(legend.position = "bottom")
            
            return(contribs)
    
    }, mc.cores = length(groups)
     
) # eo mclapply - g in groups 

# Rbind and plot results
ddf <- bind_rows(res.groups)
head(ddf) ; dim(ddf) ; summary(ddf)
rm(res.groups) ; gc()

# Plots (select sensible colours for variables)
cols <- c("SST" = "#e31a1c", "Oxygen (200m)" = "#1f78b4", "pCO2" = "#ffff99", "Wind" = "#a6cee3",
        "PAR" = "#ff7f00", "MLD" = "#d9d9d9", "N*" = "#6a3d9a", "Si*" = "#cab2d6", 
        "Chlorophyll" = "#33a02c", "Silicates" = "#b2df8a", "EKE" = "#fb9a99")

plot <- ggplot(data = ddf, aes(x = 2, y = relative, fill = factor(variable))) +
    geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
    scale_fill_manual(name = "", values = cols) + theme_void() +
    xlim(0.5,2.5) + facet_wrap(~ factor(group)) + theme(legend.position = "bottom")

ggsave(plot = plot, filename = paste("donuts_global_contrib_vars_opt_GBM_groups.pdf", sep = ""), dpi = 300, width = 9, height = 13)


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------