### R code for analysis of alignment between migratory species and resident vertebrates in the Neotropics

## STEP2 - PROCESSING THE TERRESTRIAL AREA OF HABITAT DATA

rmlist=((ls))
getwd()
library(tidyverse)
library(raster)
library(readxl)
library(ebirdst)
library(RColorBrewer)
library(viridis)
library(sf)
library(here)
library(GADMTools)
library(maptools)
library(rgdal)

memory.limit()
memory.limit(size=50000)
raster::rasterOptions(maxmemory=1e+11,chunksize=1e+10) 

substrRight<-function(x,n){
  substr(x,nchar(x)-n+1,nchar(x))
}

## The Terrestrial Area of Habitat maps need to be in the working directory.
## These maps are not yet on a public repository but are available from the authors, contact scottd.wilson@ec.gc.ca 
## Mammals are shown as an example for the process of species selection and masking to the focal area but same process exactly for other taxa

## The Terrestrial Area of Habitat maps need to be in the working directory. 
## These maps are not yet on a public repository but are available from the authors, contact scottd.wilson@ec.gc.ca 
## Mammals are shown as an example for the process of species selection and masking to the focal area but same process exactly for other taxa

setwd("C:/Users/scott/Documents/R/Spatial_Analyses/Data")

## Move species distribution raster maps for threatened species to a new folder

IUCN_mamm_list<-read_xlsx("Final Mammal List.xlsx","mammals")%>%
  filter(concern1==1) # concern 1 refers to species listed as Near Threatened, Vulnerable, Endangered or Critically Endangered as of March, 2020
IUCN_mamm_list # 543 species NT or higher, 527 species after processing below (14 missing, 2 subspecies excluded)

## Pull files from the list of all AOH files 
IUCN_mamm_files<-list.files(path=here("E:/data/AOH Data/Mammals_1km_AOH/Mammals"),pattern=".tif$") 
IUCN_mamm_files_df<-tibble(string=gsub(".tif","",IUCN_mamm_files),
                           filename=list.files(path=here("E:/data/AOH Data/Mammals_1km_AOH/Mammals"),pattern=".tif$",full.names=TRUE)
)%>%
  mutate(name=substr(string,1,nchar(string)-4)) 

IUCN_mamm_files_df_red<-IUCN_mamm_files_df[IUCN_mamm_files_df$name %in% IUCN_mamm_list$Name,]
IUCN_mamm_files_df_red
table(IUCN_mamm_files_df_red$name)
write.csv(IUCN_mamm_files_df_red,file="mammals_output.csv")

## Place the new set of reduced files with threatened species only into a separate folder
file.copy(IUCN_mamm_files_df_red$filename,here("E:/data/AOH Data/Mammals_1km_AOH/Mammals_IUCN")) 

## The next step is to mask the resident maps to include only those species that fall within the migrant focal area
## Before doing this we need to resample the focal area richness raster so that it has the same resolution as the AOH maps
## An example AOH map is used to do this

mammal_ex<-raster("E:/data/AOH Data/Mammals_1km_AOH/Mammals_IUCN/Alouatta_pigra_aoh.tif")
focalArea_1k<-resample(Richness.10,mammal_ex,method="ngb")
writeRaster(focalArea_1k, file="E:/data/intermediate/eBird_migrant_focal_rp.tif", overwrite = TRUE)

## Next step is to run the mask with gdalwarp, the code below will retain those resident AOH rasters that overlap the migrant focal area (i.e. "focalArea_1k")
## All threatened species AOH files were moved into a new folder, listed in next line, for this step since they are removed if they do not overlap the focal area

mammals<-list.files(here("E:/data/intermediate/AOH/Mammals"),full.name=TRUE)
bound_rast<-raster("E:/data/intermediate/eBird_migrant_focal_rp.tif") # this is the migrant focal area 1km raster
mammal_out <- list.files(here("E:/data/intermediate/AOH/Mammals"), full.name = TRUE)

# need to install gdal for this next step, see: https://trac.osgeo.org/osgeo4w/
# then set the path variable to know where to look.

for(ii in 1:length(mammals)) {
  system(paste("gdalwarp -r near",
               mammals[ii],
               mammal_out[ii],
               sep=" "))           
}  

for(ii in 1:length(mammal_out)){
  print(ii)
  flush.console()
  
  tmp <- mammal_out[ii] %>%
    raster() %>%
    mask(bound_rast)
  
  if(cellStats(tmp, sum) > 0){
    tmp %>%
      writeRaster(mammal_out[ii], overwrite = TRUE)
  } else {
    file.remove(mammal_out[ii])
  }
  
  rm(tmp)
}

## The remaining rasters in the folder (i.e. "E:/data/intermediate/AOH/Mammals"in this example) will be those that overlap the migrant focal area

mammal_out <- list.files(here("E:/data/intermediate/AOH/Mammals"), full.name = TRUE)
mammal_stack<-stack(mammal_out)
mammal_richness<-sum(mammal_stack, na.rm=TRUE)
writeRaster(mammal_richness, file="C:/Users/scott/Documents/R/Spatial_Analyses/Results/resident_mammal_richness.tif", overwrite = TRUE) # resident mammal richness raster, used in Figure 2
mammal_sp_count<-cellStats(mammal_stack,stat='sum',na.rm=TRUE) # calculates the number of pixels for each species within the focal area
write.csv(mammal_sp_count,file="E:/data/Result/mamm_sp_count.csv")

mamm_outF <- list.files(here("E:/data/Mammals_IUCN"), full.name = TRUE)
mamm_stackF<-stack(mamm_outF)
xF<-cellStats(mamm_stackF,stat='sum',na.rm=TRUE)
write.csv(xF,file="E:/data/Result/mamm_outputF.csv")
