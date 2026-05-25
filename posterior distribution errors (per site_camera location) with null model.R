#script to estiamte posterior distribution errors (per site/ camera location) with null model-----

#obtainig errors by bootstrapping------
#only cervus
#put this into a four 
raneflist <- list()
for (j in 1:1000) {
  
  
  data.cervus.occub <-  sample_n(data.cervus.occu, size=nrow(data.cervus.occu), replace = T)
  
  # Add a column for survey occasion
  data.cervus.occub <- data.cervus.occub %>%
    mutate(Occasion = findInterval(as.Date(eventStart), survey_occasions))
  
  # Create an empty matrix that will be populated with detection/non-detection data
  detection_data.cervus.occub <- matrix(0, nrow = length(all.sites), ncol = length(survey_occasions))
  rownames(detection_data.cervus.occub) <- all.sites
  
  #records per site per occasion
  # here we get how many individuals were seen (for the target species) in each occasion
  data.cervus.occub1 <- data.cervus.occub %>% mutate(day=substr(eventStart,1,7))%>%
    group_by(deploymentID,scientificName,Occasion,day) %>% summarize(max_animalsday=max(count))
  
  #max counts per occasion only ( the minimum count of individuals seen per occasion)
  data.cervus.occub2 <- data.cervus.occub1 %>% group_by(deploymentID, scientificName, Occasion)%>% 
    summarise(min_n_idiv=max(max_animalsday))
  
  # Fill the matrix with detections
  for (i in 1:nrow(data.cervus.occub2)) {
    site <- data.cervus.occub2$deploymentID[i]
    occasion <- data.cervus.occub2$Occasion[i]
    min_n_idiv <- filter(data.cervus.occub2, deploymentID==site & Occasion==occasion)$min_n_idiv
    detection_data.cervus.occub[as.character(site), occasion] <- min_n_idiv
    #print(i)
  }
  
  # Convert matrix to data frame
  detection_data.cervus.occu <- as.data.frame(detection_data.cervus.occu)
  
  #make sure the order of deploymentID is the same as the covriable matrix inserted later (alphabetical)
  # important for when joining with site covariable
  detection_data.cervus.occub <-detection_data.cervus.occub[order(rownames(detection_data.cervus.occub)), ]
  
  
  #Start
  for (i in 1:length(location_timerang$deploymentID)) {
    ifelse(location_timerang$minlocOccasion[i]==min(location_timerang$minlocOccasion),"no NAs - camera set on first occasion",
           detection_data.cervus.occub[location_timerang$deploymentID[i],1:(location_timerang$minlocOccasion[i])]  <- NA)
  }
  
  #End
  for (i in 1:length(location_timerang$deploymentID)) {
    ifelse(location_timerang$maxlocOccasion[i]==max(location_timerang$maxlocOccasion),"no NAs - camera set until last occasion",
           detection_data.cervus.occub[location_timerang$deploymentID[i],(location_timerang$maxlocOccasion[i]):max(location_timerang$maxlocOccasion)]  <- NA)
  }
  detection_data.cervus.occub
  
  
  #removing any ocassion that contain NAs only (optional)-------
  detection_data.cervus.occub2  <- detection_data.cervus.occub[, colSums(is.na(detection_data.cervus.occub)) < nrow(detection_data.cervus.occub)]
  print("detection matrix ready")
  
  #Null model
  #creating unmarked data type
  umf.cervus.occub <- unmarked::unmarkedFrameOccu(
    y = detection_data.cervus.occub2,
    siteCovs = tibble(deploymentID=rownames(detection_data.cervus.occu3), var= 1)
  )
  
  cervus_single_nullb<- unmarked::occu(~1 ~ 1,
                                       ,data=umf.cervus.occub) 
  raneflist[[j]] <- bup(ranef(cervus_single_nullb))
  print(j)
  
}


ranef_estimatesdf <- as.data.frame(raneflist)
colnames(ranef_estimatesdf) <- paste0("boot",1:1000)

ranef_estimatesdf2 <- as_tibble(t(ranef_estimatesdf))
ci95up <- vector()
ci95low <- vector()
for (i in 1:length(ranef_estimatesdf2)) {
  vec <- sort(unlist(as.vector(ranef_estimatesdf2[,i])))
  ci95up[i] <- vec[975]
  ci95low[i] <- vec[24]
}

boot_cols <- paste0("boot",1:1000)

ranef_estimatesdf_errors <- ranef_estimatesdf %>%
  rowwise() %>%
  mutate(
    #mean_boot_psi_hat = mean(c_across(all_of(boot_cols)), na.rm = TRUE),
    #estimate_median = median(c_across(all_of(boot_cols)), na.rm = TRUE),
    SE = sd(c_across(all_of(boot_cols)), na.rm = TRUE),
    #ci_lower_95 = quantile(c_across(all_of(boot_cols)), 0.025, na.rm = TRUE),
    #ci_upper_95 = quantile(c_across(all_of(boot_cols)), 0.975, na.rm = TRUE)
  ) %>% ungroup() %>%
  mutate(psi_hat=bup(ranef_cervus), 
         ci95up=ci95up,
         ci95low=ci95low,
         ci95low=if_else(psi_hat==1,1,ci95low),
         SE=if_else(psi_hat==1,NA,SE),
         ci95low=if_else(psi_hat<ci95low,psi_hat,ci95low))%>% 
  dplyr::select(-all_of(boot_cols)) %>% 
  dplyr::select(psi_hat, dplyr::everything())