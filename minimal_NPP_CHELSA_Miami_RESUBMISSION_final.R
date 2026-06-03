library(pastclim)
library(tidyr)
library(ggplot2)
library(dplyr)
library(binford)
library(readr)
library(readxl)
library(rstudioapi)
library(mgcv)

script_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(script_path))

####EXTERNAL DATA####
download_dataset(dataset = "CHELSA_trace21k_1.0_0.5m_vsi", bio_variables = c("bio01", "bio12"))

t <- get_time_bp_steps("CHELSA_trace21k_1.0_0.5m_vsi")

t4 <- read_csv("fishing_combo.csv") #"https://zenodo.org/records/1167852/files/Tallavaara_Dataset_4.xls?download=1"; updated with % fishing presented by Binford 2001 and Kelly 2013, average where present in both
t4_n <- subset(t4, subpop != "x") #exclude groups with potential farming contact for subsistence

sea_level_Lambeck_PNAS <- readxl::read_excel("Lambeck et al. 2014 PNAS Table S3.xlsx") #the SI
bath <- rast("C:/Users/jblin/OneDrive - The University of Liverpool/Malta Population/malta_bath_proj.tif")#this is GEBCO2019 that has been cropped in ESRI ArcMAP 10.5 to E12-16 N35-39 and projected to UTM33N

####GENERATE MIAMI ####

t4_n$time_bp <- 0
t4_n <- location_slice(x=t4_n, bio_variables = c("bio01", "bio12"), dataset = "CHELSA_trace21k_1.0_0.5m_vsi", nn_interpol = F)
t4_n$CHELSA_miami <- NA

for(i in 1:nrow(t4_n)){ # of site locations
  t4_n$NPPt[i] <- 3000/(1+ exp(1.315-0.119*t4_n$bio01[i]))
  t4_n$NPPp[i] <- 3000*(1- exp(-0.000664*t4_n$bio12[i]))
  t4_n$temp_control[i]  <- t4_n$NPPp[i]>t4_n$NPPt[i]
  t4_n$CHELSA_miami[i] <- min(t4_n$NPPt[i], t4_n$NPPp[i])
}

chelsa <- region_slice(time_bp=0, bio_variables = "bio01", dataset = "CHELSA_trace21k_1.0_0.5m_vsi")

download_etopo()
relief_rast <- load_etopo()
#relief_rast <- terra::resample(relief_rast, chelsa)
relief_proj <- project(relief_rast, "EPSG:3857")

sea_model_df <- data.frame(time_bp = sea_level_Lambeck_PNAS$ka*-1000, esl = sea_level_Lambeck_PNAS$esl)#model_df, converts from ka and has sea level
interpolated_esl <- approx(sea_model_df$time_bp, sea_model_df$esl, xout = t)$y #interpolates to regular timing of CHELSA
interpolated_esl[1] <- 0 #forces esl to be 0
interpolated_esl[2] <- 0
df_interpolated_average <- data.frame(time_bp = t, esl = interpolated_esl) 

contour_esl <- as.contour(relief_proj, levels = df_interpolated_average$esl[1])

sites <- vect(t4_n, geom = c("longitude", "latitude"), crs = "EPSG:4326")
sites_proj <- project(sites, "EPSG:3857")

dist_to_coast <- distance(sites_proj, contour_esl)

t4_n$dist_to_coast_m <- as.numeric(dist_to_coast)
t4_n$dist_to_coast_km <- t4_n$dist_to_coast_m / 1000

#subset based on fishing %

filtered50 <- subset(t4_n, fishingC <= 50)
filtered25 <- subset(t4_n, fishingC <= 25)
filtered10 <- subset(t4_n, fishingC <= 10)

precip_control <- subset(t4_n, temp_control ==F)
precip_control50 <- subset(precip_control, fishingC <= 50)
precip_control25 <- subset(precip_control, fishingC <= 25)
precip_control10 <- subset(precip_control, fishingC <= 10)

#models

CHELSA_models <- list()

CHELSA_models[[1]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = t4_n) # /100 as density data is at #/100km2 whereas Miami estimated for 1km2
CHELSA_models[[2]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = filtered50)
CHELSA_models[[3]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = filtered25)
CHELSA_models[[4]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = filtered10)

CHELSA_models[[5]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = t4_n) # /100 as density data is at #/100km2 whereas Miami estimated for 1km2
CHELSA_models[[6]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = filtered50)
CHELSA_models[[7]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = filtered25)
CHELSA_models[[8]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = filtered10)

CHELSA_models[[9]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = precip_control) # /100 as density data is at #/100km2 whereas Miami estimated for 1km2
CHELSA_models[[10]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = precip_control50)
CHELSA_models[[11]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = precip_control25)
CHELSA_models[[12]] <- gam(log10(densityC/100) ~ s(log10(CHELSA_miami)) + s(dist_to_coast_km), data = precip_control10)

CHELSA_models[[13]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = precip_control) # /100 as density data is at #/100km2 whereas Miami estimated for 1km2
CHELSA_models[[14]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = precip_control50)
CHELSA_models[[15]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = precip_control25)
CHELSA_models[[16]] <- gam(log10(densityC/100) ~ ts(log10(CHELSA_miami)) + ts(dist_to_coast_km), data = precip_control10)


library(tidyverse)

model_eval <- tibble( #set up tibble
  n = numeric(),
  dev_exp = numeric(),
  r2 = numeric(),
  edf_chelsa = numeric(),
  F_chelsa = numeric(),
  edf_dist = numeric(),
  F_dist = numeric()
)

for (i in seq_along(CHELSA_models)) { #extract model fit stats
  s <- summary(CHELSA_models[[i]])
  
  model_eval <- add_row(
    model_eval,
    n = s$n,
    dev_exp = s$dev.expl,
    r2 = if (!is.null(s$r.sq)) s$r.sq else NA,
    edf_chelsa = s$s.table[1, "edf"],
    F_chelsa = s$s.table[1, "F"],
    edf_dist = s$s.table[2, "edf"],
    F_dist = s$s.table[2, "F"]
  )
}

model_eval #model fit statistics

CHELSA_models <- CHELSA_models[c(1,2,7,8, 13:16)] #final model selection 

####CREATE THE MIAMI NPP RASTERS - SKIP TO LOAD FILE####
interpolated_list <- list()

sea_model_df <- data.frame(time_bp = sea_level_Lambeck_PNAS$ka*-1000, esl = sea_level_Lambeck_PNAS$esl)#model_df, converts from ka and has sea level
interpolated_esl <- approx(sea_model_df$time_bp, sea_model_df$esl, xout = t)$y #interpolates to regular timing of CHELSA
interpolated_esl[1] <- 0 #forces esl to be 0
interpolated_esl[2] <- 0
df_interpolated_average <- data.frame(time_bp = t, esl = interpolated_esl) #produces the time and esl interpolated to match CHELSA

sea_model_df <- data.frame(time_bp = sea_level_Lambeck_PNAS$ka*-1000, esl = sea_level_Lambeck_PNAS$`esl minus sigma`)#model_df, converts from ka and has sea level
interpolated_esl <- approx(sea_model_df$time_bp, sea_model_df$esl, xout = t)$y #interpolates to regular timing of CHELSA
interpolated_esl[1] <- 0 #forces esl to be 0
interpolated_esl[2] <- 0
df_interpolated_low <- data.frame(time_bp = t, esl = interpolated_esl) #produces the time and esl interpolated to match CHELSA

sea_model_df <- data.frame(time_bp = sea_level_Lambeck_PNAS$ka*-1000, esl = sea_level_Lambeck_PNAS$`esl plus sigma`)#model_df, converts from ka and has sea level
interpolated_esl <- approx(sea_model_df$time_bp, sea_model_df$esl, xout = t)$y #interpolates to regular timing of CHELSA
interpolated_esl[1] <- 0 #forces esl to be 0
interpolated_esl[2] <- 0
df_interpolated_high <- data.frame(time_bp = t, esl = interpolated_esl) #produces the time and esl interpolated to match CHELSA

interpolated_list <- list(df_interpolated_average, df_interpolated_low, df_interpolated_high)
names <- c("average", "low", "high")
dataset <- c("CHELSA_trace21k_1.0_0.5m_vsi")

for(j in 1:3){
NPP_list <- list()
cropped_list1 <- list()
cropped_list2 <- list()
dist_list <- list()
for(i in 1:length(t)){#for every timeslice
    cropped <- region_series(#crop
    time_bp = t[[i]], ##this timeslice
    bio_variables = c("bio01", "bio12"),#these vars
    dataset = dataset, #this data
    ext = c(12, 16, 35, 39))#this extent
  
  binary_raster <- project(bath, cropped[[1]]) #project to match crs
  binary_raster <- (binary_raster > interpolated_list[[j]][i, 2])#a raster that matches the correct sea level by timeslice
  
  #calculate NPP using MIAMI
  NPPt <- 3000/(1+ exp(1.315-0.119*(mask(cropped[[1]], binary_raster, maskvalue = FALSE)))) 
  NPPp <- 3000*(1- exp(-0.000664*(mask(cropped[[2]], binary_raster, maskvalue = FALSE)))) 
  
  NPP_list[[i]] <- min(NPPt, NPPp) #https://rdrr.io/github/Mavbegg/MIAMI/src/R/NPP.R
  cropped_list1[[i]] <- mask(cropped[[1]], binary_raster, maskvalue = FALSE)
  cropped_list2[[i]] <- mask(cropped[[2]], binary_raster, maskvalue = FALSE)
  
  
  relief <- project(relief_proj, cropped[[1]])
  relief <- mask(relief, NPP_list[[i]])
  flipped <- ifel(is.na(relief), 1, NA)
  dist_to_coast <- distance(flipped)
  dist_list[[i]] <- mask(dist_to_coast/1000, NPP_list[[i]])
  print(i)
  }

names(NPP_list) <- t #set time slices as names
names(cropped_list1) <- t
names(cropped_list2) <- t
names(dist_list) <- t
NPP <- rast(Filter(Negate(is.null), NPP_list)) #stack as rasters of NPP cropped to correct bathymetry
cropped1 <- rast(Filter(Negate(is.null), cropped_list1)) #stack as rasters of NPP cropped to correct bathymetry
cropped2 <- rast(Filter(Negate(is.null), cropped_list2)) #stack as rasters of NPP cropped to correct bathymetry
Dist <- rast(Filter(Negate(is.null), dist_list))
writeCDF(NPP, paste0("Malta_Miami_NPP_CHELSA_", names[[j]], ".nc"), overwrite = TRUE) #save it
writeCDF(cropped1, paste0("Malta_Miami_bio01_CHELSA_", names[[j]], ".nc"), overwrite = TRUE) #save it
writeCDF(cropped2, paste0("Malta_Miami_bio12_CHELSA_", names[[j]], ".nc"), overwrite = TRUE) #save it
writeCDF(Dist, paste0("Malta_Miami_Dist_CHELSA_", names[[j]], ".nc"), overwrite = TRUE)
}

####identify which parameter limits NPP via Miami####

bio01 <- rast("Malta_CHELSA_Miami_bio01_average.nc")
bio12 <- rast("Malta_CHELSA_Miami_bio12_average.nc")

NPPt <- 3000/(1+ exp(1.315-0.119*(bio01))) 
NPPp <- 3000*(1- exp(-0.000664*(bio12))) 

test <- NPPp>NPPt
test2 <- app(test, mean, na.rm=T)
plot(test2) #this shows that our Miami estimates are precipitation limited

####LOAD PROCESSED NPP RASTERS####

NPP_rasts <- list(
NPP <- rast("Malta_Miami_NPP_CHELSA_average.nc"),
bio01 <- rast("Malta_Miami_bio01_CHELSA_average.nc"),
bio12 <- rast("Malta_Miami_bio12_CHELSA_average.nc"),
Dist <- rast("Malta_Miami_Dist_CHELSA_average.nc"),
NPP_hi <- rast("Malta_Miami_NPP_CHELSA_high.nc"),
bio01_hi <- rast("Malta_Miami_bio01_CHELSA_high.nc"),
bio12_hi <- rast("Malta_Miami_bio12_CHELSA_high.nc"),
Dist_hi <- rast("Malta_Miami_Dist_CHELSA_high.nc"),
NPP_lo <- rast("Malta_Miami_NPP_CHELSA_low.nc"),
bio01_lo <- rast("Malta_Miami_bio01_CHELSA_low.nc"),
bio12_lo <- rast("Malta_Miami_bio12_CHELSA_low.nc"),
Dist_lo <- rast("Malta_Miami_Dist_CHELSA_low.nc")
)

names(NPP_rasts) <- c("NPP", "bio01", "bio12", "Dist", "NPP_hi", "bio01_hi", "bio12_hi","Dist_hi", "NPP_lo", "bio01_lo", "bio12_lo", "Dist_lo")

plot(NPP_rasts[[4]][[181]])

####GENERATE POPULATION DENSITY ESTIMATES####

cell_size_list <- list((111.32 * res(NPP_rasts[[1]])[2]) * (111.32 * cos(36 * pi/180)*res(NPP_rasts[[1]])[1]))# first is N/S, second is E/W for 36 north (converting from deg to rad), where 111.32 is km/1deg at equator

CHELSA_density_rasters <- lapply(seq_along(CHELSA_models), function(k) {
  mod <- CHELSA_models[[k]]

  # predict each layer separately
  r_list <- lapply(seq_len(nlyr(NPP_rasts[[1]])), function(j) {

    predictors <- c(NPP_rasts[[1]][[j]], NPP_rasts[[4]][[j]])
    names(predictors) <- c("CHELSA_miami", "dist_to_coast_km")
    
    10^terra::predict(predictors, mod, type="response") * #back from log10 popn density
      cell_size_list[[1]] #density to count; transforming back from log10
 
  })
    # stack back into a multilayer SpatRaster
  rast(r_list)
})

CHELSA_density_rasters_hi <- lapply(seq_along(CHELSA_models), function(k) {
  
  mod <- CHELSA_models[[k]]
  
  # predict each layer separately
  r_list <- lapply(seq_len(nlyr(NPP_rasts[[1]])), function(j) {
    
    predictors <- c(NPP_rasts[[5]][[j]], NPP_rasts[[8]][[j]])
    names(predictors) <- c("CHELSA_miami", "dist_to_coast_km")
    
    10^terra::predict(predictors, mod, type="response") * #back from log10 popn density
      cell_size_list[[1]]
  })
  # stack back into a multilayer SpatRaster
  rast(r_list)
  
})

CHELSA_density_rasters_lo <- lapply(seq_along(CHELSA_models), function(k) {
  
  mod <- CHELSA_models[[k]]
  
  # predict each layer separately
  r_list <- lapply(seq_len(nlyr(NPP_rasts[[1]])), function(j) {
    
    predictors <- c(NPP_rasts[[9]][[j]], NPP_rasts[[12]][[j]])
    names(predictors) <- c("CHELSA_miami", "dist_to_coast_km")
    
    10^terra::predict(predictors, mod, type="response") * #back from log10 popn density
      cell_size_list[[1]]
  })
  # stack back into a multilayer SpatRaster
  rast(r_list)
  
})

# Compute mean and sd for cropped1

r_mean1 <- app(NPP_rasts$bio01, mean, na.rm=TRUE)
r_sd1   <- app(NPP_rasts$bio01, sd,   na.rm=TRUE)
r_mean2 <- app(NPP_rasts$bio12, mean, na.rm=TRUE)
r_sd2   <- app(NPP_rasts$bio12, sd,   na.rm=TRUE)
NPP_mean1 <- app(NPP_rasts$NPP, mean, na.rm=TRUE)
NPP_sd1   <- app(NPP_rasts$NPP, sd,   na.rm=TRUE)

####Figure 2a####
library(terra)
library(viridis)

par(mfrow = c(2,3), mar = c(1,1,1,1), oma = c(3,3,8,1), ps = 14, cex.main = 0.8)

plot(r_mean1, main="Bio01 (°C)",
     col = colorRampPalette(c("white","yellow","orange","red","darkred","black"))(50))
plot(r_mean2, main="Bio12 (mm)",
     col = colorRampPalette(c("white","skyblue","darkblue"))(50))
plot(NPP_mean1, main="NPP (g DM m² yr)",
     col = colorRampPalette(c("yellow","green","forestgreen","darkgreen"))(50))

plot(r_sd1, col = colorRampPalette(c("white","yellow","orange","red","darkred","black"))(50))
plot(r_sd2, col = colorRampPalette(c("white","skyblue","darkblue"))(50))
plot(NPP_sd1,  col = colorRampPalette(c("yellow","green","forestgreen","darkgreen"))(50))

mtext("Mean", side = 2, line = 1, outer = TRUE, at = 0.75, cex = 1.2)
mtext("Standard Deviation", side = 2, line = 1, outer = TRUE, at = 0.25, cex = 1.2)

####Figure 2b####

density_mean_list <- list()
density_sd_list <- list()
for(i in 1:length(CHELSA_density_rasters)){
density_mean_list[[i]] <- app(CHELSA_density_rasters[[i]], mean, na.rm=T)
density_sd_list[[i]] <- app(CHELSA_density_rasters[[i]], sd, na.rm=T)
}

library(ggplot2)
library(scico)
library(dplyr)
library(tidyr)

mean_r <- rast(density_mean_list)
names(mean_r) <- paste0("mean_", 1:8)

mean_mm <- minmax(mean_r)
mean_global_min <- min(mean_mm[1, ])
mean_global_max <- max(mean_mm[2, ])

sd_r <- rast(density_sd_list)
names(sd_r) <- paste0("sd_", 1:8)

sd_mm <- minmax(sd_r)
sd_global_min <- min(sd_mm[1, ])
sd_global_max <- max(sd_mm[2, ])

rasters <- c(
  mean_r[[1]], sd_r[[1]], mean_r[[5]], sd_r[[5]],   
  mean_r[[2]], sd_r[[2]], mean_r[[6]], sd_r[[6]],   
  mean_r[[3]], sd_r[[3]], mean_r[[7]], sd_r[[7]],  
  mean_r[[4]], sd_r[[4]], mean_r[[8]], sd_r[[8]]  
)

titles <- c(
  "Mean 1", "SD 1", "Mean 5", "SD 5",
  "Mean 2", "SD 2", "Mean 6", "SD 6",
  "Mean 3", "SD 3", "Mean 7", "SD 7",
  "Mean 4", "SD 4", "Mean 8", "SD 8"
)

par(mfrow = c(4,4),  mar = c(0.5,0.5,0.5,0.5),  oma = c(3, 4, 7, 2))

for (i in 1:16) {
  ttl <- titles[i]
  r <- rasters[[i]]
  
  if (startsWith(ttl, "Mean")) {
    plot(r, col = scico(50, palette = "vik"), range = c(mean_global_min, mean_global_max))
  } else {
    plot(r, col = scico(50, palette = "cork"), range = c(sd_global_min, sd_global_max))}
  
  
  if (i %in% c(1, 3)) {mtext("Mean", side = 3, line = 0.5, cex = 1.1, font = 2)}
  if (i %in% c(2, 4)) {mtext("SD",   side = 3, line = 0.5, cex = 1.1, font = 2)}
  if (i %in% c(1, 5, 9, 13)) {row_lab <- c("Full", "<50%", "<25%", "<10%")[(i - 1) / 4 + 1]
    mtext(row_lab, side = 2, line = 2, cex = 1.1, font = 2)}
}

par(xpd = NA)
mtext("Precipitation & Temperature",side  = 3,outer = TRUE, line  = 2.5,  adj   = 0.08,  # left adjust
      cex   = 1.3,       
      font  = 2)
mtext("Precipitation Only", side  = 3,outer = TRUE, line  = 2.5,adj   = 0.82, cex   = 1.3,font  = 2)
par(xpd = FALSE)
#exp as Fig 2b 8.27 x 8.27

####some useful stats####
max(density_mean_list[[2]])

max(CHELSA_density_rasters[[2]]@pntr$range_max)
min(CHELSA_density_rasters[[2]]@pntr$range_min)

max(CHELSA_density_rasters[[5]]@pntr$range_max)
min(CHELSA_density_rasters[[5]]@pntr$range_min)

#####SUM DENSITIES BY LANDMASS TO GENERATE POPULATION SIZE####

for(i in 1:8){
patchy_rasters <- lapply(CHELSA_density_rasters[[i]], patches)
patchy_rasters_lo <- lapply(CHELSA_density_rasters_lo[[i]], patches)
patchy_rasters_hi <- lapply(CHELSA_density_rasters_hi[[i]], patches)

writeCDF(rast(patchy_rasters), paste0("Malta_CHELSA_patchy_rasters_model",i, ".nc"), overwrite = TRUE) #s
writeCDF(rast(patchy_rasters_lo), paste0("Malta_CHELSA_patchy_rasters_model",i, "_lo.nc"), overwrite = TRUE) #s
writeCDF(rast(patchy_rasters_hi), paste0("Malta_CHELSA_patchy_rasters_model",i, "_hi.nc"), overwrite = TRUE) #s
}

CHELSA_patchy_rasters <- list(rast("Malta_CHELSA_patchy_rasters_model1.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model2.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model3.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model4.nc"),
                              rast("Malta_CHELSA_patchy_rasters_model5.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model6.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model7.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model8.nc"))

CHELSA_patchy_rasters_lo <- list(rast("Malta_CHELSA_patchy_rasters_model1_lo.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model2_lo.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model3_lo.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model4_lo.nc"),
                              rast("Malta_CHELSA_patchy_rasters_model5_lo.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model6_lo.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model7_lo.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model8_lo.nc"))

CHELSA_patchy_rasters_hi <- list(rast("Malta_CHELSA_patchy_rasters_model1_hi.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model2_hi.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model3_hi.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model4_hi.nc"),
                              rast("Malta_CHELSA_patchy_rasters_model5_hi.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model6_hi.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model7_hi.nc"), 
                              rast("Malta_CHELSA_patchy_rasters_model8_hi.nc"))

density_sums <- list()
for(k in 1:length(CHELSA_density_rasters)){
  
  # Initialize a list to store results
  density_sums_list1 <- list() 
  
  # Loop through each raster layer
  for (i in 1:nlyr(CHELSA_patchy_rasters[[k]])) {
    density_sums_list1[[i]] <- zonal(CHELSA_density_rasters[[k]][[i]], CHELSA_patchy_rasters[[k]][[i]], "sum") #calculates sum of density by patch
    
    # Ensure consistent column names
    colnames(density_sums_list1[[i]]) <- c("Patch_ID", paste0("Sum_Layer_", i))
  }
  
density_sums[[k]] <- density_sums_list1}

density_sums_hi <- list()
for(k in 1:length(CHELSA_density_rasters_hi)){
  
  # Initialize a list to store results
  density_sums_list1 <- list() 
  
  # Loop through each raster layer
  for (i in 1:nlyr(CHELSA_patchy_rasters_hi[[k]])) {
    density_sums_list1[[i]] <- zonal(CHELSA_density_rasters_hi[[k]][[i]], CHELSA_patchy_rasters_hi[[k]][[i]], "sum") #calculates sum of density by patch
    
    # Ensure consistent column names
    colnames(density_sums_list1[[i]]) <- c("Patch_ID", paste0("Sum_Layer_", i))
  }
  
  density_sums_hi[[k]] <- density_sums_list1}

density_sums_lo <- list()
for(k in 1:length(CHELSA_density_rasters_lo)){
  
  # Initialize a list to store results
  density_sums_list1 <- list() 
  
  # Loop through each raster layer
  for (i in 1:nlyr(CHELSA_patchy_rasters_lo[[k]])) {
    density_sums_list1[[i]] <- zonal(CHELSA_density_rasters_lo[[k]][[i]], CHELSA_patchy_rasters_lo[[k]][[i]], "sum") #calculates sum of density by patch
    
    # Ensure consistent column names
    colnames(density_sums_list1[[i]]) <- c("Patch_ID", paste0("Sum_Layer_", i))
  }
  
  density_sums_lo[[k]] <- density_sums_list1}

density_rasters_x <- list()
density_rasters_x_hi <- list()
density_rasters_x_lo <- list()

for(k in 1:length(CHELSA_patchy_rasters)){

  density_rasters1x <- CHELSA_patchy_rasters[[k]]  # Copy structure
  
  # Loop through layers
  for (i in 1:nlyr(CHELSA_patchy_rasters[[k]])) {
    # Extract current layer
    layer <- CHELSA_patchy_rasters[[k]][[i]]
    density_sums_df <- density_sums[[k]][[i]]
    
    # Create a lookup table from density_sums_df
    lookup_table <- data.frame(
      patch_value = density_sums_df[[1]],    # Patch identifiers
      density_value = density_sums_df[[2]]# Corresponding density sums
    )
    
    # Apply classification to swap patch values with density values
    density_rasters1x[[i]] <- classify(layer, lookup_table, right = FALSE) #this lets you plot maps showing density by patch, rather than by cell values
  }
  writeCDF(rast(density_rasters1x), paste0("Malta_CHELSA_sum_density_rasters_model",k, ".nc"), overwrite = TRUE) #s
  density_rasters_x[[k]] <- density_rasters1x
}

for(k in 1:length(CHELSA_patchy_rasters_lo)){
  # Create an updated raster stack
  density_rasters1x <- CHELSA_patchy_rasters_lo[[k]]  # Copy structure
  
  # Loop through layers
  for (i in 1:nlyr(CHELSA_patchy_rasters_lo[[k]])) {
    # Extract current layer
    layer <- CHELSA_patchy_rasters_lo[[k]][[i]]
    density_sums_df <- density_sums_lo[[k]][[i]]
    
    # Create a lookup table from density_sums_df
    lookup_table <- data.frame(
      patch_value = density_sums_df[[1]],    # Patch identifiers
      density_value = density_sums_df[[2]]# Corresponding density sums
    )
    
    # Apply classification to swap patch values with density values
    density_rasters1x[[i]] <- classify(layer, lookup_table, right = FALSE) #this lets you plot maps showing density by patch, rather than by cell values
  }
  writeCDF(rast(density_rasters1x), paste0("Malta_CHELSA_sum_density_rasters_model",k, "_lo.nc"), overwrite = TRUE) #s
  density_rasters_x_lo[[k]] <- density_rasters1x
}

for(k in 1:length(CHELSA_patchy_rasters_hi)){
  # Create an updated raster stack
  density_rasters1x <- CHELSA_patchy_rasters_hi[[k]]  # Copy structure
  
  # Loop through layers
  for (i in 1:nlyr(CHELSA_patchy_rasters_hi[[k]])) {
    # Extract current layer
    layer <- CHELSA_patchy_rasters_hi[[k]][[i]]
    density_sums_df <- density_sums_hi[[k]][[i]]
    
    # Create a lookup table from density_sums_df
    lookup_table <- data.frame(
      patch_value = density_sums_df[[1]],    # Patch identifiers
      density_value = density_sums_df[[2]]# Corresponding density sums
    )
    
    # Apply classification to swap patch values with density values
    density_rasters1x[[i]] <- classify(layer, lookup_table, right = FALSE) #this lets you plot maps showing density by patch, rather than by cell values
  }
  writeCDF(rast(density_rasters1x), paste0("Malta_CHELSA_sum_density_rasters_model",k, "_hi.nc"), overwrite = TRUE) #s
  density_rasters_x_hi[[k]] <- density_rasters1x
}

####extract population size for landmass based on patch matching specific locations####

library(terra)
library(ggplot2)

density_df_list <- list(data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                         hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                         lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA),
                        data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                                   hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                                   lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA))
                        

patch_values_list <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

target_coord <- matrix(c(14.439837, 35.892727), ncol=2)  # Longitude first #target location is the Bizzle

for(k in 1:length(CHELSA_patchy_rasters)){
#k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
  # i <- 2
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[6]] <- patch_values1
  density_df_list[[k]][,6] <- unlist(density_values1)
  }

target_coord <- matrix(c(14.336429, 36.011046), ncol=2)  # middle of Comino

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
    # i <- 2
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[7]] <- patch_values1
  density_df_list[[k]][,7] <- unlist(density_values1)
}

target_coord <- matrix(c(14.2394401, 36.0463615), ncol=2)  # gozo

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
    # i <- 2
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[8]] <- patch_values1
  density_df_list[[k]][,8] <- unlist(density_values1)
}

target_coord <- matrix(c(14.996843, 37.754805), ncol=2) #the crater of Etna

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
    # i <- 2
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[9]] <- patch_values1
  density_df_list[[k]][,9] <- unlist(density_values1)
}

target_coord <- matrix(c(14.439837, 35.892727), ncol=2)  # Longitude first #target location is the Bizzle

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    # i <- 2
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[10]] <- patch_values1
  density_df_list[[k]][,10] <- unlist(density_values1)
}

target_coord <- matrix(c(14.336429, 36.011046), ncol=2)  # middle of Comino

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    # i <- 2
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[11]] <- patch_values1
  density_df_list[[k]][,11] <- unlist(density_values1)
}

target_coord <- matrix(c(14.2394401, 36.0463615), ncol=2)  # gozo

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    # i <- 2
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[12]] <- patch_values1
  density_df_list[[k]][,12] <- unlist(density_values1)
}

target_coord <- matrix(c(14.996843, 37.754805), ncol=2) #the crater of Etna

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    # i <- 2
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[13]] <- patch_values1
  density_df_list[[k]][,13] <- unlist(density_values1)
}

target_coord <- matrix(c(14.439837, 35.892727), ncol=2)  # Longitude first #target location is the Bizzle

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    # i <- 2
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[2]] <- patch_values1
  density_df_list[[k]][,2] <- unlist(density_values1)
}

target_coord <- matrix(c(14.336429, 36.011046), ncol=2)  # middle of Comino

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    # i <- 2
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[3]] <- patch_values1
  density_df_list[[k]][,3] <- unlist(density_values1)
}

target_coord <- matrix(c(14.2394401, 36.0463615), ncol=2)  # gozo

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    # i <- 2
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[4]] <- patch_values1
  density_df_list[[k]][,4] <- unlist(density_values1)
}

target_coord <- matrix(c(14.996843, 37.754805), ncol=2) #the crater of Etna

for(k in 1:length(CHELSA_patchy_rasters)){
  #k <- 1  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  density_values1 <- numeric(nlyr(CHELSA_patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    # i <- 2
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- CHELSA_patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[k]][[5]] <- patch_values1
  density_df_list[[k]][,5] <- unlist(density_values1)
}

for(k in 1:8){write.csv(density_df_list[[k]], paste0("density_df_model",k,".csv"))
patch_values_df <- do.call(cbind, patch_values_list[[k]][-1])
write.csv(patch_values_df, paste0("patch_values_list_model", k,".csv"))}

view(density_df_list[[1]])

####AREA CALCULATIONS####

#patchy_rasters2 <- rast("Malta_CHELSA_patchy_rasters_model1.nc")
patchy_rasters2 <- CHELSA_patchy_rasters[[1]]
patch_values_list <- read.csv("patch_values_list_model1.csv")

min_dist <- numeric()
total_true_area_malta <- numeric()
total_true_area_comino <- numeric()
total_true_area_gozo <- numeric()
total_true_area_sicily <- numeric()

for(i in 1:221){
id_malta <- patch_values_list[[2]]
id_comino <- patch_values_list[[3]]
id_gozo <-patch_values_list[[4]]
id_sicily <- patch_values_list[[5]]
  
p1 <- patchy_rasters2[[i]] == id_malta[[i]]
p1a <- patchy_rasters2[[i]] == id_comino[[i]]
p1b <- patchy_rasters2[[i]] == id_gozo[[i]]
p2 <- patchy_rasters2[[i]] == id_sicily[[i]]

# Replace FALSE with NA
p1_clean <- ifel(p1, 1, NA)
poly1 <- as.polygons(p1_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)

p1a_clean <- ifel(p1a, 1, NA)
poly1a <- as.polygons(p1a_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)

p1b_clean <- ifel(p1b, 1, NA)
poly1b <- as.polygons(p1b_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)

p2_clean <- ifel(p2, 1, NA)
poly2 <- as.polygons(p2_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)

plot(poly1)

d <- distance(poly1, poly2)
min_dist[[i]] <- min(d, na.rm = TRUE)

patches_1 <- project(p1, "EPSG:32633")
patches_1a <- project(p1a, "EPSG:32633")
patches_1b <- project(p1b, "EPSG:32633")
patches_2 <- project(p2, "EPSG:32633")

p1_bin <- classify(patches_1, rcl = matrix(c(-Inf, 0.5, 0,
                                             0.5, Inf, 1), ncol=3, byrow=TRUE))
area_1 <- cellSize(p1_bin, unit="m")
true_area_raster_1 <- area_1 * p1_bin
total_true_area_malta[i] <- global(true_area_raster_1, "sum", na.rm=TRUE)

p1a_bin <- classify(patches_1a, rcl = matrix(c(-Inf, 0.5, 0,
                                             0.5, Inf, 1), ncol=3, byrow=TRUE))
area_1a <- cellSize(p1a_bin, unit="m")
true_area_raster_1a <- area_1a * p1a_bin
total_true_area_comino[i] <- global(true_area_raster_1a, "sum", na.rm=TRUE)

p1b_bin <- classify(patches_1b, rcl = matrix(c(-Inf, 0.5, 0,
                                             0.5, Inf, 1), ncol=3, byrow=TRUE))
area_1b <- cellSize(p1b_bin, unit="m")
true_area_raster_1b <- area_1b * p1b_bin
total_true_area_gozo[i] <- global(true_area_raster_1b, "sum", na.rm=TRUE)

p2_bin <- classify(patches_2, rcl = matrix(c(-Inf, 0.5, 0,
                                             0.5, Inf, 1), ncol=3, byrow=TRUE))
area_2 <- cellSize(p2_bin, unit="m")
true_area_raster_2 <- area_2 * p2_bin
total_true_area_sicily[i] <- global(true_area_raster_2, "sum", na.rm=TRUE)
}

write.csv(unlist(total_true_area_malta), "malta_area.csv")
write.csv(unlist(total_true_area_comino), "comino_area.csv")
write.csv(unlist(total_true_area_gozo), "gozo_area.csv")
write.csv(unlist(total_true_area_sicily), "sicily_area.csv")

patchy_rasters2 <- CHELSA_patchy_rasters_hi[[1]]
patch_values_list <- read.csv("patch_values_list_model1.csv")

min_dist <- numeric()
total_true_area_malta_hi <- numeric()
total_true_area_comino_hi <- numeric()
total_true_area_gozo_hi <- numeric()
total_true_area_sicily_hi <- numeric()

for(i in 1:221){
  id_malta <- patch_values_list[[6]]
  id_comino <- patch_values_list[[7]]
  id_gozo <-patch_values_list[[8]]
  id_sicily <- patch_values_list[[9]]
  
  p1 <- patchy_rasters2[[i]] == id_malta[[i]]
  p1a <- patchy_rasters2[[i]] == id_comino[[i]]
  p1b <- patchy_rasters2[[i]] == id_gozo[[i]]
  p2 <- patchy_rasters2[[i]] == id_sicily[[i]]
  
  # Replace FALSE with NA
  p1_clean <- ifel(p1, 1, NA)
  poly1 <- as.polygons(p1_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
  p1a_clean <- ifel(p1a, 1, NA)
  poly1a <- as.polygons(p1a_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
  p1b_clean <- ifel(p1b, 1, NA)
  poly1b <- as.polygons(p1b_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
  p2_clean <- ifel(p2, 1, NA)
  poly2 <- as.polygons(p2_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
polys <- rbind(poly1, poly1a, poly1b)
  
crs(p1)

  d <- distance(aggregate(polys), poly2)
  min_dist[[i]] <- min(d, na.rm = TRUE)
  
  patches_1 <- project(p1, "EPSG:32633")
  patches_1a <- project(p1a, "EPSG:32633")
  patches_1b <- project(p1b, "EPSG:32633")
  patches_2 <- project(p2, "EPSG:32633")
  
  p1_bin <- classify(patches_1, rcl = matrix(c(-Inf, 0.5, 0,
                                               0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_1 <- cellSize(p1_bin, unit="m")
  true_area_raster_1 <- area_1 * p1_bin
  total_true_area_malta_hi[i] <- global(true_area_raster_1, "sum", na.rm=TRUE)
  
  p1a_bin <- classify(patches_1a, rcl = matrix(c(-Inf, 0.5, 0,
                                                 0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_1a <- cellSize(p1a_bin, unit="m")
  true_area_raster_1a <- area_1a * p1a_bin
  total_true_area_comino_hi[i] <- global(true_area_raster_1a, "sum", na.rm=TRUE)
  
  p1b_bin <- classify(patches_1b, rcl = matrix(c(-Inf, 0.5, 0,
                                                 0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_1b <- cellSize(p1b_bin, unit="m")
  true_area_raster_1b <- area_1b * p1b_bin
  total_true_area_gozo_hi[i] <- global(true_area_raster_1b, "sum", na.rm=TRUE)
  
  p2_bin <- classify(patches_2, rcl = matrix(c(-Inf, 0.5, 0,
                                               0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_2 <- cellSize(p2_bin, unit="m")
  true_area_raster_2 <- area_2 * p2_bin
  total_true_area_sicily_hi[i] <- global(true_area_raster_2, "sum", na.rm=TRUE)
}

write.csv(unlist(total_true_area_malta_hi), "malta_area_hi.csv")
write.csv(unlist(total_true_area_comino_hi), "comino_area_hi.csv")
write.csv(unlist(total_true_area_gozo_hi), "gozo_area_hi.csv")
write.csv(unlist(total_true_area_sicily_hi), "sicily_area_hi.csv")

patchy_rasters2 <- CHELSA_patchy_rasters_lo[[1]]
patch_values_list <- read.csv("patch_values_list_model1.csv")

min_dist <- numeric()
total_true_area_malta_lo <- numeric()
total_true_area_comino_lo <- numeric()
total_true_area_gozo_lo <- numeric()
total_true_area_sicily_lo <- numeric()

for(i in 1:221){

  id_malta <- patch_values_list[[10]]
  id_comino <- patch_values_list[[11]]
  id_gozo <-patch_values_list[[12]]
  id_sicily <- patch_values_list[[13]]
  
  p1 <- patchy_rasters2[[i]] == id_malta[[i]]
  p1a <- patchy_rasters2[[i]] == id_comino[[i]]
  p1b <- patchy_rasters2[[i]] == id_gozo[[i]]
  p2 <- patchy_rasters2[[i]] == id_sicily[[i]]
  
  # Replace FALSE with NA
  p1_clean <- ifel(p1, 1, NA)
  poly1 <- as.polygons(p1_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
  p1a_clean <- ifel(p1a, 1, NA)
  poly1a <- as.polygons(p1a_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
  p1b_clean <- ifel(p1b, 1, NA)
  poly1b <- as.polygons(p1b_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
  p2_clean <- ifel(p2, 1, NA)
  poly2 <- as.polygons(p2_clean, dissolve = TRUE, values = FALSE, na.rm = TRUE)
  
  polys <- rbind(poly1, poly1a, poly1b)
  
  
  crs(p1)
  
  d <- distance(aggregate(polys), poly2)
  min_dist[[i]] <- min(d, na.rm = TRUE)
  
  patches_1 <- project(p1, "EPSG:32633")
  patches_1a <- project(p1a, "EPSG:32633")
  patches_1b <- project(p1b, "EPSG:32633")
  patches_2 <- project(p2, "EPSG:32633")
  
  p1_bin <- classify(patches_1, rcl = matrix(c(-Inf, 0.5, 0,
                                               0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_1 <- cellSize(p1_bin, unit="m")
  true_area_raster_1 <- area_1 * p1_bin
  total_true_area_malta_lo[i] <- global(true_area_raster_1, "sum", na.rm=TRUE)
  
  p1a_bin <- classify(patches_1a, rcl = matrix(c(-Inf, 0.5, 0,
                                                 0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_1a <- cellSize(p1a_bin, unit="m")
  true_area_raster_1a <- area_1a * p1a_bin
  total_true_area_comino_lo[i] <- global(true_area_raster_1a, "sum", na.rm=TRUE)
  
  p1b_bin <- classify(patches_1b, rcl = matrix(c(-Inf, 0.5, 0,
                                                 0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_1b <- cellSize(p1b_bin, unit="m")
  true_area_raster_1b <- area_1b * p1b_bin
  total_true_area_gozo_lo[i] <- global(true_area_raster_1b, "sum", na.rm=TRUE)
  
  p2_bin <- classify(patches_2, rcl = matrix(c(-Inf, 0.5, 0,
                                               0.5, Inf, 1), ncol=3, byrow=TRUE))
  area_2 <- cellSize(p2_bin, unit="m")
  true_area_raster_2 <- area_2 * p2_bin
  total_true_area_sicily_lo[i] <- global(true_area_raster_2, "sum", na.rm=TRUE)
}

write.csv(unlist(total_true_area_malta_lo), "malta_area_lo.csv")
write.csv(unlist(total_true_area_comino_lo), "comino_area_lo.csv")
write.csv(unlist(total_true_area_gozo_lo), "gozo_area_lo.csv")
write.csv(unlist(total_true_area_sicily_lo), "sicily_area_lo.csv")


area_df <- data.frame(t=t,
                      malta_average = read.csv("malta_area.csv")[[2]]/1000000,
                      malta_lo = read.csv("malta_area_lo.csv")[[2]]/1000000,
                      malta_hi = read.csv("malta_area_hi.csv")[[2]]/1000000,
                      sicily_average = read.csv("sicily_area.csv")[[2]]/1000000,
                      sicily_lo = read.csv("sicily_area_lo.csv")[[2]]/1000000,
                      sicily_hi = read.csv("sicily_area_hi.csv")[[2]]/1000000,
                      gozo_average = read.csv("gozo_area.csv")[[2]]/1000000,
                      gozo_lo = read.csv("gozo_area_lo.csv")[[2]]/1000000,
                      gozo_hi = read.csv("gozo_area_hi.csv")[[2]]/1000000,
                      comino_average = read.csv("comino_area.csv")[[2]]/1000000,
                      comino_lo = read.csv("comino_area_lo.csv")[[2]]/1000000,
                      comino_hi = read.csv("comino_area_hi.csv")[[2]]/1000000)

combine_values <- function(malta, gozo, comino) {
  gc <- ifelse(gozo == comino, gozo, gozo + comino)
  total <- ifelse(malta == gc, gc, malta + gc)
  return(total)
}

area_df <- area_df %>%
  mutate(
    average_maltese = mapply(combine_values, malta_average, gozo_average, comino_average),
    hi_maltese      = mapply(combine_values, malta_hi, gozo_hi, comino_hi),
    lo_maltese      = mapply(combine_values, malta_lo, gozo_lo, comino_lo)
  )

# Find all rows where Maltese != Sicily
distinct_rows <- area_df %>%
  filter(average_maltese != sicily_average |
           hi_maltese != sicily_hi |
           lo_maltese != sicily_lo)

# Get the last distinct t based on your data order
last_distinct_t <- tail(distinct_rows$t, 1)

area_df_full <- area_df %>%
  mutate(
    average_maltese_plot = ifelse(t > last_distinct_t, average_maltese, NA),
    hi_maltese_plot      = ifelse(t > last_distinct_t, hi_maltese, NA),
    lo_maltese_plot      = ifelse(t > last_distinct_t, lo_maltese, NA)
  )

library(ggplot2)

ggplot(area_df_full, aes(x = t)) +
  # Maltese Islands
  geom_ribbon(aes(ymin = lo_maltese_plot, ymax = hi_maltese_plot, fill = "Maltese Islands"), alpha = 0.2) +
  geom_line(aes(y = average_maltese_plot, color = "Maltese Islands"), size = 1) +
  
  # Sicily
  geom_ribbon(aes(ymin = sicily_lo, ymax = sicily_hi, fill = "Sicily"), alpha = 0.2) +
  geom_line(aes(y = sicily_average, color = "Sicily"), size = 1) +
  
  scale_color_manual(values = c("Maltese Islands" = "red", "Sicily" = "blue")) +
  scale_fill_manual(values = c("Maltese Islands" = "red", "Sicily" = "blue")) +
  labs(x = "Time", y = "Population size",
       #title = "Population size: Maltese Islands vs Sicily",
       color = "Region", fill = "Region") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(),   # keep major gridlines at 10, 100, 1000, 10000
    panel.grid.minor = element_blank(),  # remove minor gridlines
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line()
  )

####PLOT for FIG 1A AND 1B####

library(patchwork)

# Save each plot to an object
p1 <- ggplot(area_df_full, aes(x = t/1000)) +
  geom_ribbon(aes(ymin = lo_maltese_plot, ymax = hi_maltese_plot, fill = "Maltese Islands"), alpha = 0.5) +
  geom_line(aes(y = average_maltese_plot, color = "Maltese Islands"), linewidth = 0.8) +
  geom_ribbon(aes(ymin = sicily_lo, ymax = sicily_hi, fill = "Sicily"), alpha = 0.5) +
  geom_line(aes(y = sicily_average, color = "Sicily"), linewidth = 0.8) +
  scale_y_log10(
    limits = c(100, 100000),
    breaks = c(100, 1000, 10000, 100000),
    labels = scales::comma
  ) +
  scale_color_manual(values = c("Maltese Islands" = "red", "Sicily" = "blue"),
                     breaks = c("Maltese Islands", "Sicily")) +
  scale_fill_manual(values = c("Maltese Islands" = "red", "Sicily" = "blue"),
                    breaks = c("Maltese Islands", "Sicily")) +
  labs(x = "Time (ka BP)", y = "Area (km2)", color = "Region", fill = "Region",
       title = "Exposed land surface area") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),
    legend.position = c(0.05, 0.15),         # relative coordinates inside plot (x,y)
    legend.justification = c("left", "bottom"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA) # semi-transparent box
  )

plot_df <- sea_level_Lambeck_PNAS[1:240,]

p2 <- ggplot(plot_df) +
  geom_ribbon(
    aes(x = -ka, ymin = `esl-2sigma`, ymax = `esl+2sigma`),
    fill = "steelblue", alpha = 0.2
  ) +
  geom_line(
    aes(x = -ka, y = esl),
    color = "steelblue", linewidth = 0.8
  ) +
  labs(
    x = "Time (ka BP)",
    y = "Eustatic Sea Level (m)",
    title = "Sea Level Reconstruction"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(),   # keep major gridlines
    panel.grid.minor = element_blank(),  # remove minor gridlines
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line()
  )

# Arrange side by side
p1 + p2

####COSTPATH DISTANCE ANALYSIS####

library(terra)
library(gdistance)

line_list <- list()
cost_dist <- numeric()
max_cross <- numeric()
sea_list <- list()

patchy_rasters2 <- CHELSA_patchy_rasters[[1]]
# store largest sea crossing per raster
library(raster)
library(gdistance)


for (i in 1:nlyr(patchy_rasters2)) {
  
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)

  etna_ll   <- vect(matrix(c(14.996843, 37.754805), ncol=2), crs="EPSG:4326")
  bizzle_ll <- vect(matrix(c(14.439837, 35.892727), ncol=2), crs="EPSG:4326")
  etna   <- project(etna_ll, patches_r)
  bizzle <- project(bizzle_ll, patches_r)

  etna_sp   <- as(etna, "Spatial")
  bizzle_sp <- as(bizzle, "Spatial")
  
  tr <- transition(raster(cost_r), function(x) 1 / mean(x), directions = 8)
  tr_corr <- geoCorrection(tr, type = "c")
  
  path_line <- shortestPath(tr_corr, etna_sp, bizzle_sp, output = "SpatialLines")
  line_list[[i]] <- path_line

  cost_dist[[i]] <- costDistance(tr_corr, etna_sp, bizzle_sp)
}

#Max Crossing Distance#

max_dist_list <- list()
max_dist <- numeric()
total_list <- list()
total_dist <- numeric()

for (i in 1:nlyr(patchy_rasters2)) {

patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")

cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
path_vect <- vect(line_list[[i]])
crs(path_vect) <- crs(cost_r)

patch_ids <- patches(patches_r, directions = 8)
patch_polys <- as.polygons(patch_ids, dissolve = TRUE, values = TRUE)
patch_polys$cost_val <- terra::extract(cost_r, patch_polys, fun = mean, na.rm = TRUE)[,2]
path_segments <- terra::intersect(patch_polys, path_vect)
x <- path_vect-path_segments
max_dist_list[[i]] <- perim(disagg(x))
max_dist[i] <- max(perim(disagg(x)))
}

line_list_hi <- list()
cost_dist_hi <- numeric()
max_cross_hi <- numeric()
sea_list_hi <- list()


patchy_rasters2 <- CHELSA_patchy_rasters_hi[[1]]

# store largest sea crossing per raster
library(raster)
library(gdistance)

for (i in 1:nlyr(patchy_rasters2)) {

  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  etna_ll   <- vect(matrix(c(14.996843, 37.754805), ncol=2), crs="EPSG:4326")
  bizzle_ll <- vect(matrix(c(14.439837, 35.892727), ncol=2), crs="EPSG:4326")
  etna   <- project(etna_ll, patches_r)
  bizzle <- project(bizzle_ll, patches_r)
  etna_sp   <- as(etna, "Spatial")
  bizzle_sp <- as(bizzle, "Spatial")
  tr <- transition(raster(cost_r), function(x) 1 / mean(x), directions = 8)
  tr_corr <- geoCorrection(tr, type = "c")
  path_line <- shortestPath(tr_corr, etna_sp, bizzle_sp, output = "SpatialLines")
  line_list_hi[[i]] <- path_line
  cost_dist_hi[[i]] <- costDistance(tr_corr, etna_sp, bizzle_sp)
  }

#Max Crossing Distance#

max_dist_list_hi <- list()
max_dist_hi <- numeric()
total_list_hi <- list()
total_dist_hi <- numeric()

for (i in 1:nlyr(patchy_rasters2)) {
  
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  path_vect <- vect(line_list[[i]])
  crs(path_vect) <- crs(cost_r)
  patch_ids <- patches(patches_r, directions = 8)
  patch_polys <- as.polygons(patch_ids, dissolve = TRUE, values = TRUE)
  patch_polys$cost_val <- terra::extract(cost_r, patch_polys, fun = mean, na.rm = TRUE)[,2]
  path_segments <- terra::intersect(patch_polys, path_vect)
  x <- path_vect-path_segments
  max_dist_list_hi[[i]] <- perim(disagg(x))
  max_dist_hi[i] <- max(perim(disagg(x)))
}

line_list_lo <- list()
cost_dist_lo <- numeric()
max_cross_lo <- numeric()
sea_list_lo <- list()

patchy_rasters2 <- CHELSA_patchy_rasters_lo[[1]]

# store largest sea crossing per raster
library(raster)
library(gdistance)


for (i in 1:nlyr(patchy_rasters2)) {
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  etna_ll   <- vect(matrix(c(14.996843, 37.754805), ncol=2), crs="EPSG:4326")
  bizzle_ll <- vect(matrix(c(14.439837, 35.892727), ncol=2), crs="EPSG:4326")
  etna   <- project(etna_ll, patches_r)
  bizzle <- project(bizzle_ll, patches_r)
  etna_sp   <- as(etna, "Spatial")
  bizzle_sp <- as(bizzle, "Spatial")
  tr <- transition(raster(cost_r), function(x) 1 / mean(x), directions = 8)
  tr_corr <- geoCorrection(tr, type = "c")
  path_line <- shortestPath(tr_corr, etna_sp, bizzle_sp, output = "SpatialLines")
  line_list_lo[[i]] <- path_line
  cost_dist_lo[[i]] <- costDistance(tr_corr, etna_sp, bizzle_sp)
  }


#Max Crossing Distance#

max_dist_list_lo <- list()
max_dist_lo <- numeric()
total_list_lo <- list()
total_dist_lo <- numeric()

for (i in 1:nlyr(patchy_rasters2)) {
  
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  path_vect <- vect(line_list[[i]])
  crs(path_vect) <- crs(cost_r)
  patch_ids <- patches(patches_r, directions = 8)
  patch_polys <- as.polygons(patch_ids, dissolve = TRUE, values = TRUE)
  patch_polys$cost_val <- terra::extract(cost_r, patch_polys, fun = mean, na.rm = TRUE)[,2]
  path_segments <- terra::intersect(patch_polys, path_vect)
  x <- path_vect-path_segments
  max_dist_list_lo[[i]] <- perim(disagg(x))
  max_dist_lo[i] <- max(perim(disagg(x)))
}


####PLOT FOR FIGURE 3####

combine_values <- function(malta, gozo, comino) {
  gc <- ifelse(gozo == comino, gozo, gozo + comino)
  total <- ifelse(malta == gc, gc, malta + gc)
  return(total)
}

max_dist_lo_test <- numeric()
sum_dist_lo_test <- numeric()
for(i in 1:length(max_dist_list_lo)){
max_dist_lo_test[i] <- max(max_dist_list_lo[[i]][max_dist_list_lo[[i]]>1200])
sum_dist_lo_test[i] <- sum(max_dist_list_lo[[i]][max_dist_list_lo[[i]]>1200])
}

max_dist_hi_test <- numeric()
sum_dist_hi_test <- numeric()
for(i in 1:length(max_dist_list_hi)){
  max_dist_hi_test[i] <- max(max_dist_list_hi[[i]][max_dist_list_hi[[i]]>1200])
  sum_dist_hi_test[i] <- sum(max_dist_list_hi[[i]][max_dist_list_hi[[i]]>1200])
}

max_dist <- numeric()
sum_dist <- numeric()
for(i in 1:length(max_dist_list)){
  max_dist[i] <- max(max_dist_list[[i]][max_dist_list[[i]]>1200])
  sum_dist[i] <- sum(max_dist_list[[i]][max_dist_list[[i]]>1200])
}


plot_df2 <- data.frame(t=t,
  max_dist = max_dist,
  max_dist_hi_test = max_dist_hi_test,
  max_dist_lo_test = max_dist_lo_test,
  sum_dist = sum_dist,
  sum_dist_hi_test = sum_dist_hi_test,
  sum_dist_lo_test = sum_dist_lo_test
  )

plot_df2[plot_df2 == -Inf] <- 0

library(ggplot2)



color_values <- c(
  "P&T Maltese Islands" = "darkorange2",
  "P&T Sicily"          = "red3",
  "P Maltese Islands"   = "purple3",
  "P Sicily"            = "blue3"
)

linetypes <- c(
  "Full" = "solid",
  "50%"  = "longdash",
  "25%"  = "twodash",
  "10%"  = "dotted"
)

density_df_all <- density_df_all %>%
  mutate(
    model_family = ifelse(grepl("P&T", scenario), "P&T", "P"),
    fishing_level = case_when(
      grepl("Full", scenario) ~ "Full",
      grepl("50%", scenario)  ~ "50%",
      grepl("25%", scenario)  ~ "25%",
      grepl("10%", scenario)  ~ "10%"
    ),
    region_maltese = "Maltese Islands",
    region_sicily  = "Sicily"
  )


l1 <- ggplot(density_df_all, aes(x = t)) +
  
  # Maltese ribbon (masked before divergence)
  geom_ribbon(
    aes(
      ymin = lo_maltese_plot,
      ymax = hi_maltese_plot,
      fill = paste(model_family, region_maltese),
      linetype = fishing_level
    ),
    alpha = 0.12
  ) +
  
  # Sicily ribbon (always present)
  geom_ribbon(
    aes(
      ymin = lo_sicily,
      ymax = hi_sicily,
      fill = paste(model_family, region_sicily),
      linetype = fishing_level
    ),
    alpha = 0.12
  ) +
  
  # Maltese line
  geom_line(
    aes(
      y = average_maltese_plot,
      color = paste(model_family, region_maltese),
      linetype = fishing_level
    ),
    linewidth = 1
  ) +
  
  # Sicily line
  geom_line(
    aes(
      y = average_sicily,
      color = paste(model_family, region_sicily),
      linetype = fishing_level
    ),
    linewidth = 1
  ) +
  
  # Manual scales 
  scale_color_manual(values = color_values) +
  scale_fill_manual(values = color_values) +
  scale_linetype_manual(values = linetypes) +
  
  # Log scale 
  scale_y_log10(
    limits = c(10, 20000),
    breaks = c(10, 100, 1000, 10000),
    labels = scales::comma
  ) +
  
  labs(
    x = "Time",
    y = "Population size",
    color = "Model × Region",
    fill  = "Model × Region",
    linetype = "Fishing level",
    title = "Population size"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "left",
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(),
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.6, "cm")
  )


####l2 plot####
l2 <- ggplot(plot_df2, aes(x = t)) +
  # Line for max_dist
  geom_line(aes(y = max_dist/1000, color = "Largest single distance"), size = 1) +
  
  # Line for sum_dist
  geom_line(aes(y = sum_dist/1000, color = "Total distance"), size = 1) +
  
  # Styling
  labs(x = "t", y = "Distance (km)",
       title = "Sea crossing distance",
       color = "Metric") +   # legend title
  theme_minimal(base_size = 14) +
  theme(legend.position = c(0.95, 0.05),   # bottom right inside plot
        legend.justification = c("right", "bottom"))

l2 <- ggplot(plot_df2, aes(x = t)) +
  # Largest single distance
  geom_line(aes(y = max_dist/1000, color = "Largest single distance"), size = 1) +
  
  # Total distance
  geom_line(aes(y = sum_dist/1000, color = "Total distance"), size = 1) +
  
  # Manual colors (steelblue / darkorange)
  scale_color_manual(values = c("Largest single distance" = "steelblue",
                                "Total distance" = "darkorange")) +
  
  labs(x = "Time", y = "Distance (km)",
       title = "Sea crossing distance",
       color = "Metric") +
  
  theme_minimal(base_size = 11) +   # match l1’s default font size
  theme(
    panel.grid.major = element_line(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),
    legend.position = c(0.95, 0.05),
    legend.justification = c("right", "bottom")
  )

library(patchwork)
l1/l2

