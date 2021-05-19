
### ================================================================================================================

library("raster")
library("sp")
library("stringr")
library("reshape2")
library("tidyverse")
library("biomod2")

WD <- getwd()

### ================================================================================================================


### Load the env stacks (winter  and summer, like January & August) for projection
setwd("/net/kryo/work/fabioben/OVERSEE/data/env_predictors/global_monthly_clims_1d/")

apr <- read.table("glob_stack_month_apr_21_02_19.txt", h = T, sep = ";")
jul <- read.table("glob_stack_month_jul_21_02_19.txt", h = T, sep = ";")
oct <- read.table("glob_stack_month_oct_21_02_19.txt", h = T, sep = ";")
jan <- read.table("glob_stack_month_jan_21_02_19.txt", h = T, sep = ";")
feb <- read.table("glob_stack_month_feb_21_02_19.txt", h = T, sep = ";")
mar <- read.table("glob_stack_month_mar_21_02_19.txt", h = T, sep = ";")
may <- read.table("glob_stack_month_may_21_02_19.txt", h = T, sep = ";")
jun <- read.table("glob_stack_month_jun_21_02_19.txt", h = T, sep = ";")
aug <- read.table("glob_stack_month_aug_21_02_19.txt", h = T, sep = ";")
sep <- read.table("glob_stack_month_sep_21_02_19.txt", h = T, sep = ";")
nov <- read.table("glob_stack_month_nov_21_02_19.txt", h = T, sep = ";")
dec <- read.table("glob_stack_month_dec_21_02_19.txt", h = T, sep = ";")

### First change some colnames
colnames(apr)[7] <- "dSST"
colnames(jul)[7] <- "dSST"
colnames(oct)[7] <- "dSST"
colnames(jan)[7] <- "dSST"
colnames(feb)[7] <- "dSST"
colnames(mar)[7] <- "dSST"
colnames(may)[7] <- "dSST"
colnames(jun)[7] <- "dSST"
colnames(aug)[7] <- "dSST"
colnames(sep)[7] <- "dSST"
colnames(nov)[7] <- "dSST"
colnames(dec)[7] <- "dSST"

colnames(apr)[13] <- "MLD"
colnames(jul)[13] <- "MLD"
colnames(oct)[13] <- "MLD"
colnames(jan)[13] <- "MLD"
colnames(feb)[13] <- "MLD"
colnames(mar)[13] <- "MLD"
colnames(may)[13] <- "MLD"
colnames(jun)[13] <- "MLD"
colnames(aug)[13] <- "MLD"
colnames(sep)[13] <- "MLD"
colnames(nov)[13] <- "MLD"
colnames(dec)[13] <- "MLD"

### Reorder both to make sure they follow the same ID order
jan$id <- factor(paste(jan$x, jan$y, sep = "_")) ; jan <- jan[order(jan$id),]
feb$id <- factor(paste(feb$x, feb$y, sep = "_")) ; feb <- feb[order(feb$id),]
mar$id <- factor(paste(mar$x, mar$y, sep = "_")) ; mar <- mar[order(mar$id),]
apr$id <- factor(paste(apr$x, apr$y, sep = "_")) ; apr <- apr[order(apr$id),]
may$id <- factor(paste(may$x, may$y, sep = "_")) ; may <- may[order(may$id),]
jun$id <- factor(paste(jun$x, jun$y, sep = "_")) ; jun <- jun[order(jun$id),]
jul$id <- factor(paste(jul$x, jul$y, sep = "_")) ; jul <- jul[order(jul$id),]
aug$id <- factor(paste(aug$x, aug$y, sep = "_")) ; aug <- aug[order(aug$id),]
sep$id <- factor(paste(sep$x, sep$y, sep = "_")) ; sep <- sep[order(sep$id),]
oct$id <- factor(paste(oct$x, oct$y, sep = "_")) ; oct <- oct[order(oct$id),]
nov$id <- factor(paste(nov$x, nov$y, sep = "_")) ; nov <- nov[order(nov$id),]
dec$id <- factor(paste(dec$x, dec$y, sep = "_")) ; dec <- dec[order(dec$id),]

### Then load the new satellite products to add ("27_11_19" in filename)
# dir()[grep("27_11_19", dir())]
new.chla <- read.table("clim_month_Chla_27_11_19.txt", h = T, sep = "\t")
new.mld <- read.table("clim_month_MLD_27_11_19.txt", h = T, sep = "\t")
new.par <- read.table("clim_month_PAR_27_11_19.txt", h = T, sep = "\t")
new.logchla <- read.table("clim_month_logChla_27_11_19.txt", h = T, sep = "\t")
new.mlpar <- read.table("clim_month_MLPAR_27_11_19.txt", h = T, sep = "\t")

# ### Check dimensions
# dim(new.chla) ; dim(mar)
# colnames(mar) ; colnames(new.chla)
# head(mar)
# head(new.chla)
# Check if longitudes match 
#summary(new.chla$x)
#summary(oct$x)
# Rotate new clims

new.mld$x2 <- new.mld$x
new.mld[new.mld$x > 180 ,"x2"] <- (new.mld[new.mld$x > 180 ,"x"]) - 360
# summary(new.mld$x2)
# ggplot() + geom_raster(aes(x = x2, y = y, fill = Dec), data = new.mld) + coord_quickmap() + scale_fill_viridis() + theme_bw()
new.chla$x2 <- new.chla$x
new.chla[new.chla$x > 180 ,"x2"] <- (new.chla[new.chla$x > 180 ,"x"]) - 360

new.par$x2 <- new.par$x
new.par[new.par$x > 180 ,"x2"] <- (new.par[new.par$x > 180 ,"x"]) - 360

new.logchla$x2 <- new.logchla$x
new.logchla[new.logchla$x > 180 ,"x2"] <- (new.logchla[new.logchla$x > 180 ,"x"]) - 360

new.mlpar$x2 <- new.mlpar$x
new.mlpar[new.mlpar$x > 180 ,"x2"] <- (new.mlpar[new.mlpar$x > 180 ,"x"]) - 360
# ggplot() + geom_raster(aes(x = x2, y = y, fill = Jul), data = new.logchla) + coord_quickmap() + scale_fill_viridis() + theme_bw()

# Re-order based on new id
new.mld$id <- factor(paste(new.mld$x2, new.mld$y, sep = "_")) ; new.mld <- new.mld[order(new.mld$id),]
new.par$id <- factor(paste(new.par$x2, new.par$y, sep = "_")) ; new.par <- new.par[order(new.par$id),]
new.mlpar$id <- factor(paste(new.mlpar$x2, new.mlpar$y, sep = "_")) ; new.mlpar <- new.mlpar[order(new.mlpar$id),]
new.chla$id <- factor(paste(new.chla$x2, new.chla$y, sep = "_")) ; new.chla <- new.chla[order(new.chla$id),]
new.logchla$id <- factor(paste(new.logchla$x2, new.logchla$y, sep = "_")) ; new.logchla <- new.logchla[order(new.logchla$id),]

# Check intersect length between ids
length(intersect(unique(new.chla$id), unique(jan$id))) # should be 64800
unique(new.chla$id)[1:250] ; unique(jan$id)[1:250]
# Checks out
 
# Replace the variables 
jan$Chl <- new.chla$Jan
jan$logChl <- new.logchla$Jan
jan$MLD <- new.mld$Jan
jan$PAR <- new.par$Jan
jan$MLPAR <- new.mlpar$Jan

 
feb$Chl <- new.chla$Feb
feb$logChl <- new.logchla$Feb
feb$MLD <- new.mld$Feb
feb$PAR <- new.par$Feb
feb$MLPAR <- new.mlpar$Feb

 
mar$Chl <- new.chla$Mar
mar$logChl <- new.logchla$Mar
mar$MLD <- new.mld$Mar
mar$PAR <- new.par$Mar
mar$MLPAR <- new.mlpar$Mar

 
apr$Chl <- new.chla$Apr
apr$logChl <- new.logchla$Apr
apr$MLD <- new.mld$Apr
apr$PAR <- new.par$Apr
apr$MLPAR <- new.mlpar$Apr

 
may$Chl <- new.chla$May
may$logChl <- new.logchla$May
may$MLD <- new.mld$May
may$PAR <- new.par$May
may$MLPAR <- new.mlpar$May

 
jun$Chl <- new.chla$Jun
jun$logChl <- new.logchla$Jun
jun$MLD <- new.mld$Jun
jun$PAR <- new.par$Jun
jun$MLPAR <- new.mlpar$Jun

 
jul$Chl <- new.chla$Jul
jul$logChl <- new.logchla$Jul
jul$MLD <- new.mld$Jul
jul$PAR <- new.par$Jul
jul$MLPAR <- new.mlpar$Jul

 
aug$Chl <- new.chla$Aug
aug$logChl <- new.logchla$Aug
aug$MLD <- new.mld$Aug
aug$PAR <- new.par$Aug
aug$MLPAR <- new.mlpar$Aug

 
sep$Chl <- new.chla$Sep
sep$logChl <- new.logchla$Sep
sep$MLD <- new.mld$Sep
sep$PAR <- new.par$Sep
sep$MLPAR <- new.mlpar$Sep

 
oct$Chl <- new.chla$Oct
oct$logChl <- new.logchla$Oct
oct$MLD <- new.mld$Oct
oct$PAR <- new.par$Oct
oct$MLPAR <- new.mlpar$Oct

 
nov$Chl <- new.chla$Nov
nov$logChl <- new.logchla$Nov
nov$MLD <- new.mld$Nov
nov$PAR <- new.par$Nov
nov$MLPAR <- new.mlpar$Nov

 
dec$Chl <- new.chla$Dec
dec$logChl <- new.logchla$Dec
dec$MLD <- new.mld$Dec
dec$PAR <- new.par$Dec
dec$MLPAR <- new.mlpar$Dec


summary(dec)

### Quck map
ggplot() + geom_raster(aes(x = x, y = y, fill = PAR), data = may) + coord_quickmap() + scale_fill_viridis() + theme_bw()


### Ok, now also update pCO2 - fit new pCO2 product by P. Landschüzter (ESSD)
dir()
pCO2 <- raster::stack("MPI-ULB-SOM_FFN_clim.nc")
pCO2[pCO2 < 0] <- 50 # replacing wrong values
# Disaggregate resolution from 1/4 to 1°
pCO2v2 <- aggregate(pCO2, fact = 4) # pCO2v2
rm(pCO2)
class(pCO2v2)

pco2 <- as.data.frame(pCO2v2, xy = T)
class(pco2)
dim(pco2)
colnames(pco2)[c(3:14)] <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
# Make sure they match coordinates
pco2$id <- factor(paste(pco2$x, pco2$y, sep = "_"))
pco2 <- pco2[order(pco2$id),]
# head(pco2) ; summary(pco2)
ggplot() + geom_raster(aes(x = x, y = y, fill = Jun), data = pco2) + coord_quickmap() + scale_fill_viridis() + theme_bw() 

jan$pCO2 <- pco2$Jan
feb$pCO2 <- pco2$Feb
mar$pCO2 <- pco2$Mar
apr$pCO2 <- pco2$Apr
may$pCO2 <- pco2$May
jun$pCO2 <- pco2$Jun
jul$pCO2 <- pco2$Jul
aug$pCO2 <- pco2$Aug
sep$pCO2 <- pco2$Sep
oct$pCO2 <- pco2$Oct
nov$pCO2 <- pco2$Nov
dec$pCO2 <- pco2$Dec

world <- map_data("world")
ggplot() + geom_raster(aes(x = x, y = y, fill = logNO3), data = feb) + scale_fill_viridis() + 
    geom_polygon(aes(x = long, y = lat, group = group), data = world, fill = "white", colour = "black", size = 0.3) +
    coord_quickmap() + theme_bw() + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
       	theme(panel.background = element_rect(fill = "white"),legend.key = element_rect(fill = "grey50"),
     		panel.grid.major = element_line(colour = "grey70",linetype = "dashed") )
# Looks ok, save

write.table(x = jan, "glob_stack_month_jan_18_09_20.txt", sep = ";")
write.table(x = feb, "glob_stack_month_feb_18_09_20.txt", sep = ";")
write.table(x = mar, "glob_stack_month_mar_18_09_20.txt", sep = ";")
write.table(x = apr, "glob_stack_month_apr_18_09_20.txt", sep = ";")
write.table(x = may, "glob_stack_month_may_18_09_20.txt", sep = ";")
write.table(x = jun, "glob_stack_month_jun_18_09_20.txt", sep = ";")
write.table(x = jul, "glob_stack_month_jul_18_09_20.txt", sep = ";")
write.table(x = aug, "glob_stack_month_aug_18_09_20.txt", sep = ";")
write.table(x = sep, "glob_stack_month_sep_18_09_20.txt", sep = ";")
write.table(x = oct, "glob_stack_month_oct_18_09_20.txt", sep = ";")
write.table(x = nov, "glob_stack_month_nov_18_09_20.txt", sep = ";")
write.table(x = dec, "glob_stack_month_dec_18_09_20.txt", sep = ";")





