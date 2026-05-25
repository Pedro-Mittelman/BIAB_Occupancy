#Single Season Single species ODer Delta old grid 

#Capreolus capreolus (roe deer) as an example

#packages-----
library(tidyverse)
library(lubridate)
library(dplyr)
library(tidyr)
library(unmarked)
library(broom)

#loading, checcking and preparing data----
deployments <- read_csv("Data/2023/deployments.csv")
observations <- read_csv("Data/2023/observations.csv")

#adding more info to the observations csv
deploymentlocsinfo <- deployments %>% dplyr::select(deploymentID,locationID,locationName)
observations <- observations %>% left_join(deploymentlocsinfo)

# Ensure Date column is in Date format
deployments$deploymentStart <- as.Date(deployments$deploymentStart)
class(deployments$deploymentStart)
deployments$deploymentEnd <- as.Date(deployments$deploymentEnd)
class(deployments$deploymentEnd)

# season selection-----
#filtering for one season 
#here I am defining a two month period to maintain the assumption of closure
observationsf <- observations %>% filter(eventStart>as.Date("2023-06-14")) %>% 
  filter(eventStart<as.Date("2023-08-15")) %>% 
  filter(observationType=="animal")

# Defining the start and end date for the survey period
#here I am defining a two month period to maintain the assumption of closure
start_date <- "2023-06-15" %>% as.Date("%Y-%m-%d")
end_date <- "2023-08-14" %>% as.Date( "%Y-%m-%d")

#exploring data with all species-----
observationsf %>% group_by(scientificName) %>% summarise(n())

# species selection----
#only Capreolus
data.Capreolus.occu <-  filter(observationsf, scientificName =="Capreolus capreolus") %>% 
  left_join(dplyr::select(deployments,deploymentID,locationID,locationName))

#creating the detection matrix----
# Creating occasions
# using one week as standard
survey_occasions <- seq.Date(start_date, end_date, by = "7 days")

# Add a column for survey occasion
data.Capreolus.occu <- data.Capreolus.occu %>%
  mutate(Occasion = findInterval(as.Date(eventStart), survey_occasions))

# Get unique sites 
all.sites <- as.factor(sort(unique(observationsf$deploymentID)))

# Create an empty matrix that will be populated with detection/non-detection data
detection_data.Capreolus.occu <- matrix(0, nrow = length(all.sites), ncol = length(survey_occasions))

rownames(detection_data.Capreolus.occu) <- all.sites

#records per site per occasion
# here we get how many individuals were seen (for the target species) in each occasion
data.Capreolus.occu1 <- data.Capreolus.occu %>% mutate(day=substr(eventStart,1,7))%>%
  group_by(deploymentID,scientificName,Occasion,day) %>% summarize(max_animalsday=max(count))

#max counts per occasion only ( the minimum count of individuals seen per occasion)
data.Capreolus.occu2 <- data.Capreolus.occu1 %>% group_by(deploymentID, scientificName, Occasion)%>% 
  summarise(min_n_idiv=max(max_animalsday))

#we need to write an if statement for when there is no count data and just use 1 and 0s for detection
#already having count data in the data matrix does not affect occupancy and will help the integration
#with future models (e.g. abundance estimation)

# Fill the matrix with detections
for (i in 1:nrow(data.Capreolus.occu2)) {
  site <- data.Capreolus.occu2$deploymentID[i]
  occasion <- data.Capreolus.occu2$Occasion[i]
  min_n_idiv <- filter(data.Capreolus.occu2, deploymentID==site & Occasion==occasion)$min_n_idiv
  detection_data.Capreolus.occu[as.character(site), occasion] <- min_n_idiv
  print(i)
}

# Convert matrix to data frame
detection_data.Capreolus.occu <- as.data.frame(detection_data.Capreolus.occu)

#make sure the order of deploymentID is the same as the covriable matrix inserted later (alphabetical)
# important for when joining with site covariable
detection_data.Capreolus.occu <-detection_data.Capreolus.occu[order(rownames(detection_data.Capreolus.occu)), ]


## Until now we cannot different 0s (no detection) from NAs (camera not active)
#let's insert the NAs
## adding start and end occasions for each tree/camera so the rest is filled with NAs
deployments <- deployments %>% mutate(Occasionstart = findInterval(deploymentStart, survey_occasions))
deployments <- deployments %>% mutate(Occasionend = findInterval(deploymentEnd, survey_occasions))

#max and minimal occasion per camera location/deployment ID
location_timerang <- deployments %>% group_by(deploymentID) %>% 
  summarise(minlocOccasion=min(Occasionstart),maxlocOccasion=max(Occasionend)) %>% 
  filter(minlocOccasion!=0&maxlocOccasion!=0)

location_timerang2 <- deployments %>% group_by(deploymentID) %>% 
  summarise(minlocOccasion=min(Occasionstart),maxlocOccasion=max(Occasionend)) %>% 
  filter(minlocOccasion!=0&maxlocOccasion!=0)

detection_data.Capreolus.occu2 <- detection_data.Capreolus.occu

#Start
for (i in 1:length(location_timerang$deploymentID)) {
  ifelse(location_timerang$minlocOccasion[i]==min(location_timerang$minlocOccasion),"no NAs - camera set on first occasion",
         detection_data.Capreolus.occu2[location_timerang$deploymentID[i],1:(location_timerang$minlocOccasion[i])]  <- NA)
}

#End
for (i in 1:length(location_timerang$deploymentID)) {
  ifelse(location_timerang$maxlocOccasion[i]==max(location_timerang$maxlocOccasion),"no NAs - camera set until last occasion",
         detection_data.Capreolus.occu2[location_timerang$deploymentID[i],(location_timerang$maxlocOccasion[i]):max(location_timerang$maxlocOccasion)]  <- NA)
}
detection_data.Capreolus.occu2
#NAs inserted 
###### add problems within the sample period if this information is available-----
#if we have information about the dates that the camera were not working WITHIN
#each deployment the user need to provide a cameraOpertaiontable as set per
#camtrapR
#most people do not use this and fix the deployment end or start date of each
#camera deployment
# I encourage this approach as then the cameratrap_dp is the only file needed

#Adding camera operation table----
# 
# camop_problem_lubridate2 <- camtrapR::cameraOperation(CTtable      = deploymentst,
#                                                       stationCol   = "locationName",
#                                                       sessionCol   = "sessionperloc",
#                                                       setupCol     = "deploymentStart",
#                                                       retrievalCol = "deploymentEnd",
#                                                       writecsv     = FALSE,
#                                                       byCamera     = FALSE,
#                                                       hasProblems  = FALSE,
#                                                       dateFormat   = "%Y-%m-%d"
# )
# 
# cam_op_table <- as_tibble((camop_problem_lubridate2),rownames=NA)
# cam_op_table <- cam_op_table %>%  
#   mutate(locationName=sub("_.*", "", row.names(cam_op_table)))
# 
# cam_op_table[is.na(cam_op_table)] <- 0
# 
# cam_op_table <- cam_op_table %>%
#   group_by(locationName) %>%
#   summarise(across(where(is.numeric), sum, na.rm = TRUE))
# 
# cam_op_table2 <- as_tibble(t(cam_op_table %>% dplyr::select(-locationName)),rownames=NA)
# 
# #adding occasion column
# cam_op_table2 <- cam_op_table2 %>% mutate(date=row.names(cam_op_table2))
# #adding one more survey occasion after the period so everythin
# cam_op_table2 <- cam_op_table2 %>% 
#   mutate(Occasion = findInterval(as.Date(date), survey_occasions))
# 
# #mutate(Occasion = findInterval(as.Date(date), c(survey_occasions, as.Date("2030-01-04"))))
# #was each camera active in each occasion?
# cam_op_tablecam <- cam_op_table2 %>%
#   dplyr::select(-date) %>% 
#   group_by(Occasion) %>% 
#   summarise(across(everything(), sum, na.rm = TRUE))
# colnames(cam_op_tablecam) <- c("Occasion",cam_op_table$locationName)
# 
# 
# cam_op_tablecamfinal <- as_tibble(ifelse(cam_op_tablecam %>% dplyr::select(-Occasion)>3,1,NA)) %>% 
#   mutate(Occasion=cam_op_tablecam$Occasion) %>% 
#   dplyr::select(Occasion,everything())
# # 
# # #adding NAs back to the matrix
# # 
# # #first filtering only the stations that are in the Capreolus data.frame
# # cam_op_tablecamfinal2 <- cam_op_tablecamfinal %>% 
# #   dplyr::select(Occasion, row.names(detection_data.Capreolus.occu2)) %>% 
# #   filter(Occasion!=0,Occasion!=158) %>% #removing dates that fall out of the survey period
# #   dplyr::select(-Occasion)
# # 
# # #reansforming 0s to NAs, dayswere cameras were not active
# # cam_op_tablecamfinal2[cam_op_tablecamfinal2==0] <-NA
# 
# #checking sizes for matrix multiplication
# #should be the same
# cam_op_tablecamfinal2 <- cam_op_tablecamfinal %>%
#   dplyr::select(-Occasion)
# length(as.matrix(cam_op_tablecamfinal2))==length(as.matrix(detection_data.Capreolus.occu2))
# if_else(length(as.matrix(cam_op_tablecamfinal2))==length(as.matrix(detection_data.Capreolus.occu2)),
#         "WOW, the matrixes have exact the same size, you can proceed!!", "ohoh, probably something wrong")
# #check in case itehy do not match 
# # rownames(t(cam_op_tablecamfinal))  
# # rownames(detection_data.Capreolus.occu2)  
# # 
# # setdiff(rownames(t(cam_op_tablecamfinal)) ,rownames(detection_data.Capreolus.occu2))
# # setdiff(rownames(detection_data.Capreolus.occu2), rownames(t(cam_op_tablecamfinal)) )
# 
# activeoccasions<- t(as.matrix(cam_op_tablecamfinal2))
# observationsocassion_Capreolus <- as.matrix(detection_data.Capreolus.occu2)
# 
# #multiplyng the matrices 
# detection_data.Capreolus.occu3 <- as_tibble(activeoccasions * observationsocassion_Capreolus) 
# row.names(detection_data.Capreolus.occu3) <- row.names(detection_data.Capreolus.occu2)
# row.names(detection_data.Capreolus.occu3)



#removing any ocassion that contain NAs only (optional)-------
detection_data.Capreolus.occu3  <- detection_data.Capreolus.occu2[, colSums(is.na(detection_data.Capreolus.occu2)) < nrow(detection_data.Capreolus.occu2)]
print("detection matrix ready")

#Null model----
#creating unmarked data type
umf.Capreolus.occu <- unmarked::unmarkedFrameOccu(
  y = detection_data.Capreolus.occu3,
  siteCovs = tibble(deploymentID=rownames(detection_data.Capreolus.occu3), var= 1)
)
summary(umf.Capreolus.occu)
plot(umf.Capreolus.occu)


Capreolus_single_null<- unmarked::occu(~1 ~ 1,
                                  ,data=umf.Capreolus.occu) 
summary(Capreolus_single_null)

#obtaining posterior estimation per site-------
ranef_Capreolus<- ranef(Capreolus_single_null)
#probability of occupancy per site
ranef_Capreolus@post
bup(ranef_Capreolus)

#script to estimate posterior distribution errors (per site/ camera location) with null model
# is being devloped-----


#Inserting environmental covariables-----
# if they are in a table already
# this tablemust have latitude and longitude or any other info that match the camera locations
site_covs <-  read_csv("Data/all_years_together/OD_oldgrid_vars2.bilinearex.csv")
str(site_covs)
print ( "environmental covariables table loaded succesfully")

#add extraction from a geotiff or any other raster data

#scaling covariables (I am a fan of scaling, for model fitting and comparing
#covariables but it might make interpretaion harder )

## we can give this as an option for the user
site_covs_scaled <- site_covs
cols_to_scale <- !names(site_covs) %in% c("longitude", "latitude", "deploymentID", "locationName", "locationID")
site_covs_scaled[ , cols_to_scale] <- scale(site_covs[ , cols_to_scale])
#putting in the same order as the detection data (alphabetical by depolyment ID)
site_covs_scaled <-site_covs_scaled[order(rownames(site_covs_scaled)), ]
print ("covariables sucessfully scaled")

# adding information about deployments
site_covs_scaled <- deployments %>% dplyr::select(deploymentID, locationName) %>% left_join(site_covs_scaled)
#putting in the same order as the detection data (alphabetical by depolyment ID)
site_covs_scaled <-site_covs_scaled[order(rownames(site_covs_scaled)), ]
#filtering for only the sites present in the detection matrix

site_covs_scaled <- site_covs_scaled %>%
  filter(deploymentID %in% rownames( detection_data.Capreolus.occu3))

#Creating unmarked frame occu-------
umf.Capreolus.occu_cov<- unmarked::unmarkedFrameOccu(
  y = detection_data.Capreolus.occu3,
  siteCovs = site_covs_scaled
)
print("Unmarked frame successfully created")
summary(umf.Capreolus.occu_cov)
plot(umf.Capreolus.occu_cov)

#creating model with predictors------
# some variables selected as example
Capreolus_single_season_cov<- unmarked::occu(~trees  ~
                                           trees + human_footprint,
                                   data=umf.Capreolus.occu_cov)
summary(Capreolus_single_season_cov)

#later we can add a model selection here

#graphs------

covs <- Capreolus_single_season_cov@data@siteCovs
#getting mean of all variable to plot
typical_value <- function(x) {
  if (is.numeric(x)) mean(x, na.rm=TRUE) else levels(x)[1]
}

#creating data frame to predict
trees_seq <- seq(min(site_covs_scaled$trees[!is.na(site_covs_scaled$trees)]), max(site_covs_scaled$trees[!is.na(site_covs_scaled$trees)]), length.out=101)
nd <- as.data.frame(lapply(covs, typical_value))
nd<- nd[rep(1, length(trees_seq)), , drop = FALSE]
nd$trees <- trees_seq
occ_tree_Capreolus <- unmarked::predict(Capreolus_single_season_cov, type="state", newdata=nd, re.form=NA, level=0.95)
occ_tree_Capreolus$Species <- "Capreolus capreolus"
occ_tree_Capreolus$trees <- trees_seq 

plotCapreolustrees <- ggplot(data= occ_tree_Capreolus, aes(x=trees, y=Predicted))+
  geom_line(linewidth=2)+
  ylim(0,1)+
  geom_ribbon( aes(ymin = lower, ymax = upper, fill=Species), fill="darkgreen", alpha = 0.1)+
  labs(x="Forested area in one kilometer buffer (scaled)", y="Occupancy")+
  theme_minimal();plotCapreolustrees

