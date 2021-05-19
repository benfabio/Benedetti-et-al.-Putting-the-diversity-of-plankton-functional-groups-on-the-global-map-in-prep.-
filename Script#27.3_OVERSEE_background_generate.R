
##### 25/08/2020 - ETHZ - Fabio Benedetti © UP Group, IBP, ETH Zürich
##### Script for : 
#	- Fitting total-target and group-target background for the species-level occurrences
#   - First, inform the group from classif
#   - define which groups present enough data for a group-target background
#   - define the total background
 
### Last update: 27/08/2020

# --------------------------------------------------------------------------------------------------------------------------------

library("rgeos")
library("raster")
library("maptools")
library("rgdal")
library("tidyverse")
library("stringr")
library("reshape2")
library("geosphere")
library("ncdf4")
library("classInt")
library("parallel")

world2 <- map_data("world2")
world <- map_data("world")
WD <- getwd() 

# --------------------------------------------------------------------------------------------------------------------------------

### Load phyto and zoo classification
# setwd("/net/kryo/work/fabioben/OVERSEE/data/biology") ; dir()
# classif <- read.csv("classif_all_modelled_taxa.csv", h = T, sep = ";")
# levels(classif$class)[levels(classif$class) == "Larvacea"] <- "Appendicularia"
# classif$group2 <- "Other"
# classif[classif$Category == "Phytoplankton","group2"] <- "Other_phyto"
# classif[classif$Category == "Zooplankton","group2"] <- "Other_zoo"
#
# classif[classif$group == "Bacillariophyceae","group2"] <- "Diatoms"
# classif[classif$group == "Dinoflagellata","group2"] <- "Dinoflagellates"
# classif[classif$class %in% c("Prymnesiophyceae","Coccolithophyceae"),"group2"] <- "Coccolithophores"
# classif[classif$group == "Jellyfish","group2"] <- "Jellyfish"
# classif[classif$group == "Chaetognatha","group2"] <- "Chaetognatha"
# classif[classif$group == "Foraminifera","group2"] <- "Foraminifera"
# classif[classif$order %in% c("Salpida","Doliolida"),"group2"] <- "Salps"
# classif[classif$order %in% c("Copelata"),"group2"] <- "Appendicularians"
# classif[classif$order == "Euphausiacea","group2"] <- "Euphausiids"
# classif[classif$order == "Amphipoda","group2"] <- "Amphipods"
# classif[classif$order == "Thecosomata","group2"] <- "Pteropods"
# classif[classif$order == "Calanoida","group2"] <- "Calanoids"
# classif[classif$family == "Oithonidae","group2"] <- "Oithonids"
# classif[classif$family %in% c("Oncaeidae","Corycaeidae","Sapphirinidae"),"group2"] <- "Poecilostomatoids"
# classif$group2 <- factor(classif$group2)
# summary(classif$group2)
# classif[is.na(classif$group2),] # should be zero
# setwd(WD)

# --------------------------------------------------------------------------------------------------------------------------------

### 1°) Inform groups of interest
#cat <- "Zooplankton"
cat <- "Phytoplankton"

# Define the background sampling strategy (either "total" or "target_group") and draw psAbs 
strategy <- "group" 
#strategy <- "total" 

# Choose diretcory as a FUN of cat
if(cat == "Phytoplankton") {
    
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/no_background")
    wd2 <- getwd()
    files <- dir()[grep("data_fitted",dir())]
    
} else if(cat == "Zooplankton") {
   
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/no_background")
    wd2 <- getwd()
    files <- dir()[grep("data_fitted",dir())]
    
} # eo if else loop

### Rbind all and provide 'group2'
#f <- files[5]
res <- mclapply(files, function(f) {
            message(paste("Retrieving ",f, sep = ""))
            message(paste("", sep = ""))
            data <- read.table(f, h = T, sep = ";")
            return(data)
        }, mc.cores = 30
    
) # eo mclapply
# Rbind
table <- bind_rows(res); rm(res); gc()
# head(table); dim(table)
# summary(table)

# Choose diretcory as a FUN of cat
if(cat == "Phytoplankton") {
    
    table$group2 <- "Other_phyto"
    table[table$group == "bacillariophyceae","group2"] <- "Diatoms"
    table[table$group == "dinoflagellata","group2"] <- "Dinoflagellates"
    table[table$class %in% c("Prymnesiophyceae","Coccolithophyceae","Haptophyta incertae sedis",NA),"group2"] <- "Coccolithophores" 
    # summary(factor(table$group2))
    
    # Drop rotated x coord
    table <- subset(table, select = -c(obs))
      
} else if(cat == "Zooplankton") {
   
    table$group2 <- "Other_zoo"
    
    table[table$phylum %in% c("Cnidaria","Ctenophora"),"group2"] <- "Jellyfish"
    
    table[table$phylum %in% c("Chaetognatha"),"group2"] <- "Chaetognatha"
    
    table[table$phylum %in% c("Foraminifera"),"group2"] <- "Foraminifera"
    
    table[table$class %in% c("Thaliacea"),"group2"] <- "Salps"
    
    table[table$class %in% c("Appendicularia","Larvacea"),"group2"] <- "Appendicularians"
    
    table[table$order %in% c("Euphausiacea"),"group2"] <- "Euphausiids"
    
    table[table$order %in% c("Amphipoda"),"group2"] <- "Amphipods"
    
    table[table$order %in% c("Thecosomata"),"group2"] <- "Pteropods"
    
    table[table$order %in% c("Calanoida"),"group2"] <- "Calanoida"
    
    table[table$family %in% c("Oithonidae"),"group2"] <- "Oithonida"
    
    table[table$family %in% c("Oncaeidae","Corycaeidae","Sapphirinidae"),"group2"] <- "Poecilostomatoida"
    
    # summary(factor(table$group2)) 
    #unique(table[table$group2 == "Other_zoo","species"])
    
} # eo if else loop

# nsp <- length(unique(table[,"species"]))
# nsp
#
# groups <- data.frame(table %>% group_by(group2) %>% summarize(Nocc = n(), Pocc = (n()/nrow(table))*100))
# groups[order(groups$Pocc, decreasing = T),]
#
# ggplot(groups, aes(x = 2, y = Pocc, fill = factor(group2))) +
#   geom_bar(stat = "identity", color = "black") + coord_polar(theta = "y", start = 0) +
#   scale_fill_brewer(name = "Group", palette = "Paired") + theme_void() + xlim(0.5,2.5)

### DECISION:
# Salps --> not enough, total only
# Appendicularia --> not enough, total only
# Forams --> not enough, total only, though likely not reliable
# For group-target, put Oithonida and Cyclopoida together
### Threshold for groups who cannot be tested with target-group background: 20'000 occurrences (Other_zoo)


if(cat == "Phytoplankton" & strategy == "total") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
    wd.back <- getwd()
} else if(cat == "Zooplankton" & strategy == "total") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
    wd.back <- getwd()
} else if(cat == "Phytoplankton" & strategy == "group") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background")
    wd.back <- getwd()
} else if(cat == "Zooplankton" & strategy == "group") {
    setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
    wd.back <- getwd()
} 

# Specify parameters to stratify the sampled environment:
vec.strat <- c("SST")

# Remove brackets of species names
if(cat == "Zooplankton") {
    table$species <- gsub("\\(|\\)", "", table$species)
    species <- unique(table$species)# ; species   
} else if(cat == "Phytoplankton") {
    table$species <- str_replace_all(table$species," ","_")
    species <- unique(table$species)# ; species   
}

#sp <- species[420]; sp  # for testing

### Sites id
table$id2 <- paste(table$xbin_1d, table$ybin_1d, table$month, sep = "_")
 
for(sp in species) {
	
		 message(paste("Drawing psAbs for ",sp, " ============================================ ", sep = ""))
		 grp <- unique(table[table$species == sp,"group2"])
								
		 # Get species data
		 subset <- table[table$species == sp,] ; n <- nrow(subset)
         # Identify sites where the species if found as present
         sites.with.pres <- unique(subset$id2)
			
		 # Choice of background data depending on the "strategy"
		 if( strategy == "total" ) {	
				
			    # Compute n obs per species
			    riches <- data.frame(table[!(table$id2 %in% sites.with.pres),] %>% group_by(id2) %>% summarise(n = n(), group.div = length(unique(group2))))
                # summary(riches)
            
                if(cat == "Zooplankton") {
                
                        cells2keep <- riches[riches$group.div >= 3,"id2"] # length(cells2keep)
    		            bckgrnd <- table[table$id2 %in% cells2keep,]
                
                } else if(cat == "Phytoplankton") {
                
                        cells2keep <- riches[riches$group.div >= 3,"id2"] # length(cells2keep)
    			        bckgrnd <- table[table$id2 %in% cells2keep,] # dim(bckgrnd)
                
                }
										
	      } else if( strategy == "group" ) {
					
			    # Then use the groups' data as bckgrnd
			    if(grp %in% c("Calanoida","Chaetognatha","Pteropods","Euphausiids","Jellyfish","Diatoms","Dinoflagellates","Coccolithophores")) {
						
                         bckgrnd <- table[table$group2 == grp,]	
                         bckgrnd <- bckgrnd[!(bckgrnd$id2 %in% sites.with.pres),] 
                    
			     } else if(grp %in% c("Oithonida","Poecilostomatoida")) {
						
                         bckgrnd <- table[table$group2 %in% c(grp,"Oithonida"),]	
                         bckgrnd <- bckgrnd[!(bckgrnd$id2 %in% sites.with.pres),] 
                        
			     } else {
					    
    			         riches <- data.frame(table[!(table$id2 %in% sites.with.pres),] %>% group_by(id2) %>% summarise(n = n(), group.div = length(unique(group2))))
    			         # Keep only those occ of table that belong to these monthly resolved 1° cells
                         if(cat == "Zooplankton") {
                             cells2keep <- riches[riches$group.div >= 3,"id2"] # length(cells2keep)
        			         bckgrnd <- table[table$id2 %in% cells2keep,]
                         } else if(cat == "Phytoplankton") {
                             cells2keep <- riches[riches$group.div >= 3,"id2"] # length(cells2keep)
        			         bckgrnd <- table[table$id2 %in% cells2keep,]
                         }
                       
			      }
							
		  } # eo if else loop
          
          # intersect(unique(bckgrnd$id2), sites.with.pres) # should be 0
								
		  ### Specify range of SST into which the values fall to drive sampling of absences proportionally to the overall presences points
		  x_envir <- bckgrnd[,c(vec.strat[1])] # vector of SST in the bckgrnd data
		  #y_envir <- bckgrnd[,c(vec.strat[2])]
		  # Split ranges into 9  strata
		  breaks <- 9
		  ### Create a matrix that divides range into 9 equal parts; with two variables we get a maximum of 81 strata
		  x_breaks <- classIntervals(na.omit(x_envir), breaks, style = "equal")
		  x_matrix <- cbind(x_breaks$brks[1:breaks], x_breaks$brks[2:(breaks+1)], ID = 1:breaks)
		  colnames(x_matrix) <- c("low","up","ID")	   
          #y_breaks <- classIntervals(na.omit(y_envir), breaks, style = "equal")
		  #y_matrix <- cbind(y_breaks$brks[1:breaks], y_breaks$brks[2:(breaks + 1)], ID = 1:breaks )
		  #colnames(y_matrix) <- c("low","up","ID")  
          # Define vector of length of total points of environmental variable
		  x_reclass <- c(1:length(x_envir))
		  #y_reclass <- c(1:length(y_envir))					
		  # Allocate points from full data to one of the nine environmental strata per variable
		  for(i in 1:breaks) {	
		      x_reclass[which(x_envir >= x_matrix[i,"low"] & x_envir <= x_matrix[i,"up"] )] <- x_matrix[i,"ID"]
			  #y_reclass[which(y_envir >= y_matrix[i,"low"] & y_envir <= y_matrix[i,"up"] )] <- y_matrix[i,"ID"] 		
		  } # eo for loop

		  ### Create an ID indicating the stratum (unique combination of variables) into which each point falls in full data-frame
		  bckgrnd$x_rcls <- x_reclass
		  #bckgrnd$y_rcls <- y_reclass
		  #bckgrnd$xy_rcls <- x_reclass+10*y_reclass
								
		  print( paste0(sp,", ",grp," | n = ",n, " | drawing psAbs")) # eo print
								
		  ### Extract frequencies by which points/sites of the target group fall into environmental strata. 
		  # Then, derive the number of desired absences for the focal model species per stratum. 
		  x_rcls_freq <- data.frame( table(bckgrnd$x_rcls) / length(bckgrnd$x_rcls) ) 		
		  # Give name to column
		  colnames(x_rcls_freq)[1] <- "x_rcls" 
		  # Convert to numeric
		  x_rcls_freq$x_rcls <- as.numeric(as.character(x_rcls_freq$x_rcls)) 
		  ### Add desired background points to be produced per stratum: generally 10 x more absences than presences
		  x_rcls_freq$prop_abs <- (n*10)*x_rcls_freq$Freq
		  # To round desired absences to integer: adds column difference between smaller closest integer and desired number
		  x_rcls_freq$prop_abs_0 <- ( x_rcls_freq$prop_abs - floor(x_rcls_freq$prop_abs) )
		  # To add column with random number between 0 and 1 (with steps of 0.01)
		  x_rcls_freq$prob <- sample(seq(0, 1, 0.01), nrow(x_rcls_freq), replace = T)
		  # To add column with "1"
		  x_rcls_freq$absences <- 1
		  # Round up absences for random subset
		  x_rcls_freq$absences[which(x_rcls_freq$prop_abs_0 > x_rcls_freq$prob)] <- ceiling(x_rcls_freq$prop_abs[which(x_rcls_freq$prop_abs_0 > x_rcls_freq$prob)])
		  # Round absences down for random subset
		  x_rcls_freq$absences[which(x_rcls_freq$prop_abs_0 < x_rcls_freq$prob)] <- floor(x_rcls_freq$prop_abs[which(x_rcls_freq$prop_abs_0 < x_rcls_freq$prob)]) 
		  
          # Skip strata without presences
		  absence_groups <- x_rcls_freq[x_rcls_freq$absences > 0,] 
		  # Select backround data, here including the points/sites of the focal species ('overlapping background')
		  absence_table <- bckgrnd ; gc()
			
		  # Randomly select background pts
		 # nnn <- nrow(absence_groups) 
          
          nnn <- 9
          
		  #require("parallel")
		  psAbs <- lapply(X = c(1:nnn), FUN = function(i) {
	
						# Select available absences within stratum in question
						message(paste(i, sep = ""))
						grp_abs_table <- absence_table[absence_table$x_rcls == absence_groups[i,"x_rcls"],]

						# Define the max nb of absences that can be drawn 
						absence_num <- ifelse(
								# Test if the number of desired background pts is bigger than the available background points
								absence_groups[i,"absences"] > nrow( absence_table[absence_table$x_rcls == absence_groups[i,"x_rcls"],]),		
								# if TRUE the potential points are insufficient - however, save the number of available points as absence_num
								nrow(absence_table[absence_table$x_rcls == absence_groups[i,"x_rcls"],]),		
								# ELSE: save the number of desired background points as absence_num   		
								absence_groups[i,"absences"]
						) # eo if else loop
						
						# Randomly sample the background points from the table containing all possible absences for the stratum in question
						sampled_grp_abs_table <- grp_abs_table[sample(1:nrow(grp_abs_table), size = absence_num),]
						rm(absence_num, grp_abs_table)
				
						return(sampled_grp_abs_table)
				
					} # eo fun
			
		   ) # eo lapply
			
		   ### Merge presences (obs = 1) with absences (obs = 0)
		   pseudoabs <- data.frame(do.call("rbind", psAbs), obs = 0)
		   occ_table <- rbind(data.frame(subset, obs = 1),  pseudoabs[,c(colnames(subset),"obs")])
				
           rm(psAbs, nnn) ; gc()
           
		   ### Create column with weights = 1; weights are associated with presences and absences for modelling
		   occ_table$weights <- 1
		   # Compute ratio of presences to absences
		   abs_ratio <- nrow(occ_table[occ_table$obs == 1,]) / nrow(occ_table[occ_table$obs == 0,])
			
		   # Add the ratio as weight for the psAbs
		   occ_table$weights[occ_table$obs == 0] <- abs_ratio # For observation that are absences we add the ratio
		   row.names(occ_table) <- c(1:nrow(occ_table)) # 
                
           # Final checks
           # dim(occ_table) ; colnames(occ_table) ; str(occ_table)
           # summary(occ_table)
           if( length(unique(occ_table$obs)) == 2 ) {
               
    		   ### Save the data to train some ENMs later
    		   setwd(wd.back)
    		   message(paste("Saving species dataset for ",sp, " ============================================ ", sep = ""))
    		   write.table(occ_table, paste("data_",strategy,"_",sp,".txt", sep = ""), sep = ";")
			
    		   ### Clean some stuff 
    		   rm(occ_table, abs_ratio, pseudoabs, absence_table, x_rcls_freq, x_reclass, x_matrix, x_envir)
    		   gc()
				
               setwd(wd2)
               
           } else {
               
               message(paste("  ", sep = ""))
               message(paste(" |||  ISSUE, NO PSEUDO-ABSENCES WERE GENERATED  ||| ", sep = ""))
               message(paste("  ", sep = ""))
               
           }
            	
} # eo for loop

#require("ggthemes")
#ggplot() + geom_point(aes(x = x, y = y, colour = factor(obs)), data = occ_table, alpha = .5) +
 #    geom_polygon(aes(x = long, y = lat, group = group), data = world,
  #       fill = "grey75", colour = "grey25", size = 0.2) +
   #  coord_quickmap() + theme_map()


### 27/08/2020: Making sanity checks for the various datasets created: look at proprotions of 1 and 0 for each file
categories <- c("Phytoplankton","Zooplankton")
strategies <- c("total","group")
c <- "Zooplankton"
s <- "total"

for(c in categories) {
    
    message(paste("Checking background data for ", c, "  ", sep = ""))
    message(paste("", sep = ""))
    message(paste("", sep = ""))
    
    for(s in strategies) {
        
        message(paste("Checking background data for ", s, sep = ""))
        message(paste("", sep = ""))
   
        if(c == "Phytoplankton" & s == "total") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/total_background")
            wd.back <- getwd()
        } else if(c == "Zooplankton" & s == "total") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/total_background")
            wd.back <- getwd()
        } else if(c == "Phytoplankton" & s == "group") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/phytoplankton/group_background")
            wd.back <- getwd()
        } else if(c == "Zooplankton" & s == "group") {
            setwd("/net/kryo/work/fabioben/OVERSEE/data/biology/data_for_group_studies/zooplankton/group_background")
            wd.back <- getwd()
        } 
        
        files <- dir()[grep("data",dir())]
        
        for(f in files) {
            
            # f <- files[5]
            d <- read.table(f, sep = ";", h = T)
            avg <- round(mean(d$obs, na.rm = T),4)
            
            if(avg < 0.1) {
                message(paste("", sep = ""))
            } else {
                message(paste("Mean proportion of 1:0 for ",unique(d[d$obs == 1,"species"])," ||  ",avg, sep = ""))
                require("ggthemes")
                map <- ggplot() + geom_point(aes(x = x, y = y, colour = factor(obs)), data = d, alpha = .5) +
                     geom_polygon(aes(x = long, y = lat, group = group), data = world,
                         fill = "grey75", colour = "black", size = 0.2) +
                     coord_quickmap() + theme_map()
               # save
               setwd(WD)
               ggsave(plot = map, filename = paste("map_bckgrnd_",s,"_",unique(d[d$obs == 1,"species"]),".jpg", sep = ""), dpi = 300, width = 4, height = 7)
               setwd(wd.back)
                
            }
            
        } # eo for loop 
        
    } # eo for loop 
    
} # eo for loop 


# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------
