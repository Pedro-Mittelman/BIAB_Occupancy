library(camtrapR)
library(unmarked)
library(terra)
library(elevatr)
library(sf)

# ── User Defined ───────────────────────────────────────────────────────────────

# number of occasions per sample occasion
n_occ = 5

# set date window for inference
custom_start = "2023-06-15"
custom_end = "2023-08-14"

# manual set a single species... for now
species = "Sus scrofa"

#variable of interest to generate a marginal effects plot
name_var <- "TRI"

##possibility for the user to upload their own covariates per site
#site.covs_user <- read.csv("site.covs_user.csv")


# ── Data ───────────────────────────────────────────────────────────────────────

# camtrapR function that reads camtrapDP to generate camOp and deteciton histories
# User input? Function needs the file paths

# camtrapdp = readcamtrapDP(
#   deployments_file  = "C:\\Users\\wong\\Desktop\\Bon-in-a-box\\deployments.csv",
#   media_file        = "C:\\Users\\wong\\Desktop\\Bon-in-a-box\\media.csv",
#   observations_file = "C:\\Users\\wong\\Desktop\\Bon-in-a-box\\observations.csv",
#   datapackage_file  = "C:\\Users\\wong\\Desktop\\Bon-in-a-box\\datapackage.json"
# )
# I think it is better if this is a relative path 
# need to check how it works in bon in a box

camtrapdp = readcamtrapDP(
  deployments_file  = "Data/camtrapDP_Donana/camtrapDP/deployments.csv",
  media_file        = "Data/camtrapDP_Donana/camtrapDP/media.csv",
  observations_file = "Data/camtrapDP_Donana/camtrapDP/observations.csv",
  datapackage_file  = "Data/camtrapDP_Donana/camtrapDP/datapackage.json"
)
print("camtrapdp read successfully")

# ── Camera operation and detection history ─────────────────────────────────────

camOp = cameraOperation(
  CTtable           = camtrapdp$CTtable,
  stationCol        = "Station",
  setupCol          = "Setup_date",
  retrievalCol      = "Retrieval_date",
  dateFormat        = "%Y-%m-%d",
  hasProblems       = TRUE,
  occasionStartTime = 0,
  writecsv          = FALSE
)

dh = as.data.frame(detectionHistory(
  recordTable              = camtrapdp$recordTable,
  species                  = species,
  camOp                    = camOp,
  output                   = "count",
  stationCol               = "locationID.x",
  speciesCol               = "scientificName",
  recordDateTimeCol        = "DateTimeOriginal",
  recordDateTimeFormat     = "%Y-%m-%d %H:%M:%S",
  occasionLength           = n_occ,
  minActiveDaysPerOccasion = 1,
  day1                     = "survey",
  includeEffort            = FALSE,
  scaleEffort              = FALSE,
  datesAsOccasionNames     = TRUE,
  writecsv                 = FALSE
))

# set the column names to dates so it can be restricted based on custom survey start/end dates
colnames(dh) = colnames(camOp)[seq(1, ncol(camOp), n_occ)]

# format column names to dates to be able to draw out the specified survey window
col_dates   = as.Date(colnames(dh))
in_window   = col_dates >= as.Date(custom_start) & col_dates <= as.Date(custom_end)
dh_filtered = dh[, in_window, drop = FALSE]

print("detection history created")


### Point Estimates

# ── Null model ─────────────────────────────────────────────────────────────────
umf_null <- unmarkedFrameOccu(y = dh_filtered)
fit_null  <- occu(~1 ~1, umf_null)
summary(fit_null)

print("Null model fitted")


  # ── Null model posterior estimates ────────────────────────────────────────────
psi_hat_null <- bup(ranef(fit_null))

results_null <- data.frame(
  station = rownames(dh_filtered),
  psi_hat = psi_hat_null,
  row.names = NULL
)
print("Posterior estimate of occupancy per station")
print(results_null)

### Spatial Predictions/Rasters 

# How much of the buffers for covariate extraction and area boundary should/can be user defined?

# ── 1. Station coordinates ─────────────────────────────────────────────────────

coords = unique(camtrapdp$CTtable[, c("Station", "latitude", "longitude")])

# SpatVector for terra package; for creating polygon around outer points in next step
stations_wgs84 = vect(
  coords,
  geom = c("longitude", "latitude"),
  crs  = "EPSG:4326"
)

stations_utm = project(stations_wgs84, "EPSG:32630")

# ── 2. Concave hull + 1km buffer = survey boundary ─────────────────────────────
survey_hull     = hull(stations_utm, type = "concave_ratio", param = 0.3)
survey_boundary = buffer(survey_hull, width = 1000)

# ── 3. Download DEM for survey boundary ───────────────────────────────────────
boundary_sf = st_as_sf(project(survey_boundary, "EPSG:4326"))

# z is the zoom level (between 1-14), z = 11 is 30m resolution
dem_raw = get_elev_raster(boundary_sf, z = 11, clip = "locations")
dem     = project(rast(dem_raw), "EPSG:32630")

# ── 4. Derive terrain indices at native resolution ────────────────────────────
TRI_native = terrain(dem, v = "TRI")

# Clip native resolution rasters to boundary for extraction
# Maybe redundant here since we already clip the raw DEM by the boundary...
elevation_native = crop(mask(dem, survey_boundary), survey_boundary)
TRI_native       = crop(mask(TRI_native, survey_boundary), survey_boundary)

names(elevation_native) = "elevation"
names(TRI_native)       = "TRI"

# Keep the finer resolution for extraction of covariate values
terrain_stack_native = c(elevation_native, TRI_native)

# ── 5. Extract mean within 500m buffer per station (from native resolution) ───
station_buffers   = buffer(stations_utm, width = 500)
terrain_extracted = terra::extract(terrain_stack_native, station_buffers, fun = mean, na.rm = TRUE)

terrain_df           = terrain_extracted[, -1, drop = FALSE]
rownames(terrain_df) = coords$Station
terrain_df           = terrain_df[rownames(dh_filtered), , drop = FALSE]

# ── 6. Scale covariates — store mean and SD for later use in prediction ────────
elev_mean = mean(terrain_df$elevation, na.rm = TRUE)
elev_sd   = sd(terrain_df$elevation,   na.rm = TRUE)
tri_mean  = mean(terrain_df$TRI,       na.rm = TRUE)
tri_sd    = sd(terrain_df$TRI,         na.rm = TRUE)

terrain_df_scaled = terrain_df
terrain_df_scaled$elevation = (terrain_df$elevation - elev_mean) / elev_sd
terrain_df_scaled$TRI       = (terrain_df$TRI - tri_mean)  / tri_sd

# ── 4b. Resample to target resolution for prediction raster only ──────────────

# resolution here set to 500m to help predictions go a little bit faster
target_res = 500

template     = rast(ext(dem), resolution = target_res, crs = crs(dem))
elevation_rs = resample(dem,        template, method = "bilinear")
TRI_rs       = resample(TRI_native, template, method = "bilinear")

elevation_rs = crop(mask(elevation_rs, survey_boundary), survey_boundary)
TRI_rs       = crop(mask(TRI_rs,       survey_boundary), survey_boundary)

names(elevation_rs) = "elevation"
names(TRI_rs)       = "TRI"

terrain_stack = c(elevation_rs, TRI_rs)  # this is used only for prediction

# ── 6. Fit occupancy model ─────────────────────────────────────────────────────
umf = unmarkedFrameOccu(y = dh_filtered, siteCovs = terrain_df)

fit = occu(~1 ~TRI + elevation, umf)
summary(fit)

print("Model with covariates sucessfully fitted")

# ── 7. Predict occupancy across survey boundary raster ────────────────────────
pred_df           = as.data.frame(terrain_stack, xy = FALSE, na.rm = FALSE)
colnames(pred_df) = c("elevation", "TRI")

psi_pred          = predict(fit, type = "state", newdata = pred_df)

occ_rast          = terrain_stack[[1]]  # use elevation layer as template
values(occ_rast)  = psi_pred$Predicted
names(occ_rast)   = "occupancy"

plot(occ_rast)

print("Prediction raster created")


# ── 8. Create rasters for SE and confidence intervals ─────────────────────────

# SE raster
se_rast         = terrain_stack[[1]]
values(se_rast) = psi_pred$SE
names(se_rast)  = "occupancy_SE"

# Lower 95% CI raster
lower_rast          = terrain_stack[[1]]
values(lower_rast)  = psi_pred$lower
names(lower_rast)   = "occupancy_lower_CI"

# Upper 95% CI raster
upper_rast          = terrain_stack[[1]]
values(upper_rast)  = psi_pred$upper
names(upper_rast)   = "occupancy_upper_CI"

print("Prediction error rasters created")


# ── 9. Stack all four rasters together ────────────────────────────────────────
occ_stack <- c(occ_rast, se_rast, lower_rast, upper_rast)
names(occ_stack) <- c("Predicted", "SE", "Lower_CI", "Upper_CI")

plot(occ_stack)


# ── 10. Save raster ───────────────────────────────────────────────────────────
# writeRaster(
#   occ_rast,
#   filename  = paste("C:\\Users\\wong\\Desktop\\Bon-in-a-box\\", 
#                      species, 
#                      "_occupancy_predicted.tif", sep = ""),
#   overwrite = TRUE
# )

writeRaster(occ_stack, paste("Outputs\\", 
                             species, 
                             "_occupancy_predictions.tif", sep = ""), overwrite = TRUE)

#11. marginal effects graph 


covs <- fit@data@siteCovs
#getting mean of all variable to plot
typical_value <- function(x) {
  if (is.numeric(x)) mean(x, na.rm=TRUE) else levels(x)[1]
}

#select variable to plot the effect on occupancy
#example
name_var <- "TRI"
var_plot <- covs[,name_var]

#creating data frame to predict
seq_var <- seq(min(var_plot), max(var_plot), length.out=101)
nd <- as.data.frame(lapply(covs, typical_value))
nd<- nd[rep(1, length(seq_var)), , drop = FALSE]
nd[,name_var] <- seq_var
occ_pred <- unmarked::predict(fit, type="state", newdata=nd, re.form=NA, level=0.95)
occ_pred$Species <- species
occ_pred$seq_var <- seq_var
occ_pred$seq_var <- seq_var
names(occ_pred)[names(occ_pred) == "seq_var"] <- name_var

library(ggplot2)
plotocc <- ggplot2::ggplot(data= occ_pred, aes(x=.data[[name_var]], y=Predicted))+
  geom_line(linewidth=2)+
  ylim(0,1)+
  geom_ribbon( aes(ymin = lower, ymax = upper, fill=Species), fill="darkgreen", alpha = 0.1)+
  labs(x=name_var, y="Occupancy")+
  theme_minimal();plotocc

print("Prediction raster created")
