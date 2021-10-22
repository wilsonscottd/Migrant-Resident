### R code for analysis of alignment between migratory species and resident vertebrates in the Neotropics
##-------------------------------------------------------------------------------------------------------------------------------------------------------
## The code below is for the initial step to download and stack eBird nonbreeding distribution maps for the focal species
## The supplemental file with extension "ebird_list.xlsx" has the full list of species available in eBird, 
#  the column "focal" has the listed species used in this analysis

rmlist=((ls))
getwd()
library(tidyverse)
library(raster)
library(readxl)
library(ebirdst)
library(RColorBrewer)
library(viridis)
library(here)
library(GADMTools)
library(maptools)
library(rgdal)

setwd("C:/Users/scott/Documents/R/Spatial_Analyses/Data")
## Process eBird nonbreeding season data
ebird_df<-read_xlsx("ebird_list.xlsx") # after moving to working directory and renaming
ebird_list<-list()

for (ii in 1:nrow(ebird_df)){
  sp_path<-ebirdst_download(species=ebird_df$species_code[ii],path=here("E:/data/Neotrop"),force=TRUE)
  nb_layer<-load_raster("abundance_seasonal", path=sp_path)
  if(any(names(nb_layer) == "nonbreeding")&(ebird_df$focal==1)) { # select the abundance_seasonal rasters with "nonbreeding" and the focal species 
    tmp_r<-raster::subset(nb_layer, grep('nonbreeding', names(nb_layer), value = T)) 
    tmp_r[tmp_r==0]<-NA 
    
    # select areas with 95% abundance to avoid including peripheral areas with very low abundance
    tmp_r <- tmp_r >= quantile(tmp_r, probs = c(0.05))
    tmp_r[tmp_r == 0] <- NA
    ebird_list[[ii]] <- tmp_r
  }else {
    ebird_list[[ii]]<-nb_layer[[1]]%>%raster()%>%setValues(NA)
  } # for all others that don't meet criteria set raster to 0
  
  names(ebird_list[[ii]]) <- ebird_df$species_code[ii]
}

nb_stack<-stack(ebird_list)
richness<-sum(nb_stack, na.rm=TRUE) # create the focal richness raster, next line sets NA for all those outside it
richness[richness == 0] <- NA
writeRaster(richness, file="C:/Users/scott/Documents/R/Spatial_Analyses/Results/eBird_Nonbreed_richness_focal.tif", overwrite = TRUE) # save the nonbreeding stacked raster, this was used for Figure S1

## Create the migrant focal area raster based on minimum species cutoff
pal<-colorRampPalette(plasma(10,direction=-1))
plot(richness,col=pal(256),maxpixels=ncell(richness)) # quick plot to visualize

## Raster math on richness raster
cellStats(richness,stat='mean',na.rm=TRUE) # average richness of cells with values is 1.90 (won't work with na.rm=FALSE)
cellStats(richness,stat='min',na.rm=TRUE) # min richness is 1
cellStats(richness,stat='max',na.rm=TRUE) # max richness is 10
quantile(richness,probs=c(seq(0.05,1.0,by=0.05)),na.rm=TRUE) # find number of species corresponding to upper 10% of pixels with highest richness (= 4species)

richness.10 <- (richness>=4) # at least 4 species corresponding to 10% of pixels
Richness.10 <- richness.10 %>% na_if(0) # this raster will be used with resident data after resampling to same resolution (see next step)
##-----------------------------------------------------------------------------------------------------------------------------------------------------------