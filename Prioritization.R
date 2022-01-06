### Code below is to run the prioritization analyses using prioritizr package
rmlist=((ls))
library(gurobi)
library(slam)
library(prioritizr)
library(prioritizrdata)
library(raster)
library(tidyverse)
library(rgdal)
library(here)

memory.limit()
memory.limit(size=50000)
raster::rasterOptions(maxmemory=1e+12,chunksize=1e+10) #100GB max memory, 10GB chunksize
raster::rasterOptions(tmpdir="E:/data/temp/storage/")

# load and view planning unit data
setwd("C:/Users/scott/Documents/R/Spatial_Analyses/Data/Mig_Res_Prior")
units<-raster("eBird_migrant_focal_rp.tif") # the migrant focal area raster is used for the planning units with each unit having equal cost
print(units)
PA_ext<-raster("PA_1km.tif") # protected areas raster
PA_loss<-raster("PA_ForestLoss.tif") # raster containing protected areas and areas of projected forest loss by 2050 (joint layer needed for prioritization1), see text for how forest loss was defined

# lock in the existing protected areas into the prioritization, these account for approximately 22% of the migrant focal area, therefore for the
# goal of 30% protected area, the prioritization will focus on identifying 8% that is unprotected, projected to lose forest and maximizes the number
# of resident species included
PA_locked <- ifelse(!is.na(PA_ext[][!is.na(units[])]), TRUE, FALSE)

# create the resident raster stack,this will include all terrestrial area of habitat (AOH) maps for resident species that overlapped the focal area 
# (note the AOH maps are not yet on a public repository, contact the authors for access)

# load and stack features
all.taxa.M<-list.files(here("C:/Users/scott/Documents/R/Spatial_Analyses/Data/Mammals_IUCN/"),full.name=TRUE)
all.taxa.B<-list.files(here("C:/Users/scott/Documents/R/Spatial_Analyses/Data/Birds_IUCN/"),full.name=TRUE)
all.taxa.R<-list.files(here("C:/Users/scott/Documents/R/Spatial_Analyses/Data/Reptiles_IUCN"),full.name=TRUE)
all.taxa.A<-list.files(here("C:/Users/scott/Documents/R/Spatial_Analyses/Data/Amphibians_IUCN"),full.name=TRUE)

all_stack.M<-stack(all.taxa.M)
all_stack.B<-stack(all.taxa.B)
all_stack.R<-stack(all.taxa.R)
all_stack.A<-stack(all.taxa.A)

all_stack<-stack(all_stack.M,all_stack.B,all_stack.R,all_stack.A)

compareRaster(units,all_stack)
print(all_stack)

### First prioritization focused on maximizing resident species richness within projected forest loss areas in the migrant focal layer. 

# Creation of the rij matrix below helps speed up processing during the prioritization
rij <- rij_matrix(PA_loss, all_stack)

feat <- data.frame(id = 1:nlayers(all_stack),
                   name = names(all_stack),
                   stringsAsFactors = FALSE) # dataframe for biodiversity features to be prioritized (i.e. Neotropical residents)

mig_df <- data.frame(id = 1:ncell(PA_loss),
                     cost = PA_loss[]) %>% 
                      drop_na() # dataframe for the units layer representing areas of projected forest loss and protected areas within the migrant focal area

# Prioritization based on a maximum utility objective, for more detail see prioritizr package and https://prioritizr.net/

PA_locked <- ifelse(!is.na(PA_ext[][!is.na(PA_loss[])]), TRUE, FALSE)
area_amt<-round((freq(units,digits=0,value=1))*0.30)-(length(PA_locked[PA_locked==TRUE])) # how many new pixels needed to reach 30% including already protected areas
budget<-area_amt+(length(PA_locked[PA_locked==TRUE])) # total budget is the 30% total area and equals protected areas plus the area we want to prioritize to reach 30%

p1 <- problem(mig_df$cost, features=feat, rij_matrix=rij, run_checks=FALSE) %>%
  add_max_utility_objective(budget) %>% # set the objective to minimize the cost of the solution while making sure targets are met
  #add_feature_weights(weights) %>% # no weights assigned to species in this analysis but a possibility (e.g. ranking by IUCN status)
  add_locked_in_constraints(PA_locked) %>% # protected areas locked in
  add_binary_decisions() %>% # set a binary decision to the problem, default is binary = include the planning unit when it meets the criteria
  add_default_solver() # identify the best solver on the system (gurobi)
print(p1)
s1<-solve(p1)

# reconnect spatial information
r1 <- PA_loss
r1[][!is.na(r1[])] <- s1
plot(r1, main="", breaks=c(0,0.5,1),col=c("grey70", "darkgreen"))

# check feature representation - may be too large for memory, use function below instead
# f1<-eval_feature_representation_summary(p1,s1)
# summary(f1$relative_held)

res_processing <- function(rij_in = NULL, solution = NULL){
  
  tmp <- rij_in %*% solution
  nms <- tmp@Dimnames[[1]]
  tmp <- as.vector(tmp)
  names(tmp) <- nms
  
  return(tmp)
}

feat_rep <- res_processing (rij, s1)
rij_stat <- rowSums(rij, na.rm =T)
output_p1<-data.frame(feat_rep,rij_stat) # data frame where for each species shows the number of pixels in the focal area and the number selected by the prioritization
colnames(output_p1)<-c("selected","total")
setwd("C:/Users/scott/Documents/R/Spatial_Analyses/Results")
write.csv(output_p1, file="maxutility_lossonly_Aug17.csv")
writeRaster(r1, filename="max_utility_lossonly_Aug17.tif", overwrite=TRUE) # output raster showing the prioritized area, used to create Figure 3 in the manuscript

### Second prioritization focused on maximizing resident species richness within the migrant focal layer without considering forest loss


rij <- rij_matrix(units, all_stack) # rij matrix

feat <- data.frame(id = 1:nlayers(all_stack),
                   name = names(all_stack),
                   stringsAsFactors = FALSE) # dataframe for biodiversity features to be prioritized (i.e. Neotropical residents)

mig_df <- data.frame(id = 1:ncell(units),
                     cost = units[]) %>% 
                      drop_na() # dataframe for the units layer representing the migrant focal area

###----- maximum utility objective - all species -------

area_amt<-round((freq(units,digits=0,value=1))*0.30)-(length(PA_locked[PA_locked==TRUE])) # how many new pixels needed to reach 30% including already protected areas
budget<-area_amt+(length(PA_locked[PA_locked==TRUE])) # total budget is the 30% total area and equals protected areas plus the area we want to prioritize to reach 30%

p2 <- problem(mig_df$cost, features=feat, rij_matrix=rij, run_checks=FALSE) %>%
  add_max_utility_objective(budget) %>% # set the objective to minimize the cost of the solution while making sure targets are met
  #add_feature_weights(weights) %>% # no weights assigned to species in this analysis but a possibility (e.g. ranking by IUCN status)
  add_locked_in_constraints(PA_locked) %>% # protected areas locked in
  add_binary_decisions() %>% # set a binary decision to the problem, default is binary = include the planning unit when it meets the criteria
  add_default_solver() # identify the best solver on the system (gurobi)
print(p2)
s2<-solve(p2)

# reconnect spatial information 
r2 <- units
r2[][!is.na(r1[])] <- s2
plot(r2, main="", breaks=c(0,0.5,1),col=c("grey70", "darkgreen"))

# check feature representation but may be too large for memory and so can use "res_processing" function below instead
# f2<-eval_feature_representation_summary(p2,s2)
# summary(f2$relative_held)

res_processing <- function(rij_in = NULL, solution = NULL){
  
  tmp <- rij_in %*% solution
  nms <- tmp@Dimnames[[1]]
  tmp <- as.vector(tmp)
  names(tmp) <- nms
  
  return(tmp)
}

feat_rep <- res_processing (rij, s2)
rij_stat <- rowSums(rij, na.rm =T)
output_p2<-data.frame(feat_rep,rij_stat) # data frame where for each species shows the number of pixels in the focal area and the number selected by the prioritization
colnames(output_p2)<-c("selected","total")
setwd("C:/Users/scott/Documents/R/Spatial_Analyses/Results")
write.csv(output_p2, file="maxutility_all_Aug17.csv")
writeRaster(r2, filename="max_utility_all_Aug17.tif", overwrite=TRUE) # output raster showing the prioritized area, used to create Figure S2 in the manuscript