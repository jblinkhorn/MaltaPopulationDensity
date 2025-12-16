library(pastclim)
library(tidyr)
library(ggplot2)
library(dplyr)
library(binford)
library(readr)
library(readxl)
library(rstudioapi)

script_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(script_path))

####EXTERNAL DATA####
download_dataset(dataset = "CHELSA_trace21k_1.0_0.5m_vsi", bio_variables = c("bio01", "bio12"))
t <- get_time_bp_steps("CHELSA_trace21k_1.0_0.5m_vsi")
t4 <- read_excel("Tallavaara_Dataset_4.xls") #"https://zenodo.org/records/1167852/files/Tallavaara_Dataset_4.xls?download=1"
t4_n <- subset(t4, subpop != "x") #exclude groups with potential farming contact for subsistence
sea_level_Lambeck_PNAS <- readxl::read_excel("Lambeck et al. 2014 PNAS Table S3.xlsx") #the SI
bath <- rast("malta_bath_proj.tif")#this is GEBCO2019 thats been cropped in ESRI ArcMAP 10.5 to E12-16 N35-39 and projected to UTM33N

####GENERATE MIAMI FROM CHELSA DATA####

t4_n$time_bp <- 0
t4_n <- location_slice(x=t4_n, bio_variables = c("bio01", "bio12"), dataset = "CHELSA_trace21k_1.0_0.5m_vsi", nn_interpol = F)
t4_n$CHELSA_miami <- NA

for(i in 1:nrow(t4_n)){ # of site locations
  NPPt <- 3000/(1+ exp(1.315-0.119*t4_n$bio01[i]))
  NPPp <- 3000*(1- exp(-0.000664*t4_n$bio12[i]))
  t4_n$CHELSA_miami[i] <- min(NPPt, NPPp)
}

CHELSA_miami_model <- lm(log10(t4_n$densityC/100) ~ t4_n$CHELSA_miami) # /100 as density data is at #/100km2 whereas Miami estimated for 1km2
CHELSA_miami_model
summary(CHELSA_miami_model)
model_df <- data.frame(model=NA, intercept=NA, slope=NA)

#Regression#

model_df[1,1] <- "CHELSA_miami_model"
model_df[1,2] <- coef(CHELSA_miami_model)[1] # 0.3050509
model_df[1,3] <- coef(CHELSA_miami_model)[2]

ggplot(t4_n, aes(x = CHELSA_miami, y = log10(densityC/100))) + #densityC is in individuals per 100km2
  geom_point(alpha = 0.6) +  # scatter plot of the data
  geom_smooth(method = "lm", se = TRUE, color = "blue") +  # regression line with confidence interval
  labs(
    x = "CHELSA_miami",
    y = "Log10 Density (per km)",
    title = "Linear Regression: log10(Density) ~ CHELSA_miami"
  ) +
  theme_minimal()

####CREATE THE MIAMI NPP RASTERS - SKIP TO LOAD FILE####

sea_model_df <- data.frame(time_bp = sea_level_Lambeck_PNAS$ka*-1000, esl = sea_level_Lambeck_PNAS$esl)#model_df, converts from ka and has sea level
interpolated_esl <- approx(sea_model_df$time_bp, sea_model_df$esl, xout = t)$y #interpolates to regular timing of CHELSA
interpolated_esl[1] <- 0 #forces esl to be 0
interpolated_esl[2] <- 0
df_interpolated_average <- data.frame(time_bp = t, esl = interpolated_esl) #produces the time and esl interpolated to match CHELSA

sea_model_df <- data.frame(time_bp = sea_level_Lambeck_PNAS$ka*-1000, esl = sea_level_Lambeck_PNAS$`esl-2sigma`)#model_df, converts from ka and has sea level
interpolated_esl <- approx(sea_model_df$time_bp, sea_model_df$esl, xout = t)$y #interpolates to regular timing of CHELSA
interpolated_esl[1] <- 0 #forces esl to be 0
interpolated_esl[2] <- 0
df_interpolated_low <- data.frame(time_bp = t, esl = interpolated_esl) #produces the time and esl interpolated to match CHELSA

sea_model_df <- data.frame(time_bp = sea_level_Lambeck_PNAS$ka*-1000, esl = sea_level_Lambeck_PNAS$`esl+2sigma`)#model_df, converts from ka and has sea level
interpolated_esl <- approx(sea_model_df$time_bp, sea_model_df$esl, xout = t)$y #interpolates to regular timing of CHELSA
interpolated_esl[1] <- 0 #forces esl to be 0
interpolated_esl[2] <- 0
df_interpolated_high <- data.frame(time_bp = t, esl = interpolated_esl) #produces the time and esl interpolated to match CHELSA

interpolated_list <- list(df_interpolated_average, df_interpolated_low, df_interpolated_high)
names <- c("average", "low", "high")

for(j in 1:3){
NPP_list <- list()
cropped_list1 <- list()
cropped_list2 <- list()

for(i in 1:221){#for every timeslice
  cropped <- region_series(#crop
    time_bp = t[[i]], ##this timeslice
    bio_variables = c("bio01", "bio12"),#these vars
    dataset = "CHELSA_trace21k_1.0_0.5m_vsi", #this data
    ext = c(12, 16, 35, 39))#this extent
  
  binary_raster <- project(bath, cropped[[1]]) #project to match crs
  binary_raster <- (binary_raster > interpolated_list[[j]][i, 2])#a raster that matches the correct sea level by timeslice
  
  #calculate NPP using MIAMI
  NPPt <- 3000/(1+ exp(1.315-0.119*(mask(cropped[[1]], binary_raster, maskvalue = FALSE)))) #https://rdrr.io/github/Mavbegg/MIAMI/src/R/NPP.R
  NPPp <- 3000*(1- exp(-0.000664*(mask(cropped[[2]], binary_raster, maskvalue = FALSE)))) #https://rdrr.io/github/Mavbegg/MIAMI/src/R/NPP.R
  
  NPP_list[[i]] <- min(NPPt, NPPp) #https://rdrr.io/github/Mavbegg/MIAMI/src/R/NPP.R
  cropped_list1[[i]] <- mask(cropped[[1]], binary_raster, maskvalue = FALSE)
  cropped_list2[[i]] <- mask(cropped[[2]], binary_raster, maskvalue = FALSE)
  }

names(NPP_list) <- t #set time slices as names
names(cropped_list1) <- t
names(cropped_list2) <- t
NPP <- rast(Filter(Negate(is.null), NPP_list)) #stack as rasters of NPP cropped to correct bathymetry
cropped1 <- rast(Filter(Negate(is.null), cropped_list1)) #stack as rasters of NPP cropped to correct bathymetry
cropped2 <- rast(Filter(Negate(is.null), cropped_list2)) #stack as rasters of NPP cropped to correct bathymetry
writeCDF(NPP, paste0("Malta_CHELSA_Miami_NPP_", names[[j]], ".nc"), overwrite = TRUE) #save it
writeCDF(cropped1, paste0("Malta_CHELSA_Miami_bio01_", names[[j]], ".nc"), overwrite = TRUE) #save it
writeCDF(cropped2, paste0("Malta_CHELSA_Miami_bio12_", names[[j]], ".nc"), overwrite = TRUE) #save it

}

####identify which parameter limits NPP via Miami####

bio01 <- rast("Malta_CHELSA_Miami_bio01_average.nc")
bio12 <- rast("Malta_CHELSA_Miami_bio12_average.nc")

NPPt <- 3000/(1+ exp(1.315-0.119*(bio01))) #https://rdrr.io/github/Mavbegg/MIAMI/src/R/NPP.R
NPPp <- 3000*(1- exp(-0.000664*(bio12))) #https://rdrr.io/github/Mavbegg/MIAMI/src/R/NPP.R

test <- NPPp>NPPt
test2 <- app(test, mean, na.rm=T)
plot(test2) #this shows that our Miami estimates are precipitation limited

####LOAD PROCESSED NPP RASTERS####

NPP <- rast("Malta_CHELSA_Miami_NPP_average.nc")
bio01 <- rast("Malta_CHELSA_Miami_bio01_average.nc")
bio12 <- rast("Malta_CHELSA_Miami_bio12_average.nc")

NPP_hi <- rast("Malta_CHELSA_Miami_NPP_high.nc")
bio01_hi <- rast("Malta_CHELSA_Miami_bio01_high.nc")
bio12_hi <- rast("Malta_CHELSA_cMiami_bio12_high.nc")

NPP_lo <- rast("Malta_CHELSA_Miami_NPP_low.nc")
bio01_lo <- rast("Malta_CHELSA_Miami_bio01_low.nc")
bio12_lo <- rast("Malta_CHELSA_cMiami_bio12_low.nc")

####GENERATE POPULATION DENSITY ESTIMATES####

cell_size_km2 <- (111.32 * res(NPP)[1]) * (111.32 * res(NPP)[2]) 

density_rasters <- lapply(seq_len(nrow(model_df)), function(k) {
  intercept <- model_df$intercept[k]
  slope     <- model_df$slope[k]
  
  app(NPP, function(x) {
    log10_density <- intercept + slope * x
    pmax(10^log10_density, 0) * cell_size_km2
  })
})

density_rasters_hi <- lapply(seq_len(nrow(model_df)), function(k) {
  intercept <- model_df$intercept[k]
  slope     <- model_df$slope[k]
  
  app(NPP_hi, function(x) {
    log10_density <- intercept + slope * x
    pmax(10^log10_density, 0) * cell_size_km2
  })
})

density_rasters_lo <- lapply(seq_len(nrow(model_df)), function(k) {
  intercept <- model_df$intercept[k]
  slope     <- model_df$slope[k]
  
  app(NPP_lo, function(x) {
    log10_density <- intercept + slope * x
    pmax(10^log10_density, 0) * cell_size_km2
  })
})

# Compute mean and sd for cropped1
r_mean1 <- app(bio01, mean, na.rm=TRUE)
r_sd1   <- app(bio01, sd,   na.rm=TRUE)
r_mean2 <- app(bio12, mean, na.rm=TRUE)
r_sd2   <- app(bio12, sd,   na.rm=TRUE)
NPP_mean1 <- app(NPP, mean, na.rm=TRUE)
NPP_sd1   <- app(NPP, sd,   na.rm=TRUE)
density_mean <- app(rast(density_rasters), mean, na.rm=T)
density_sd <- app(rast(density_rasters), sd, na.rm=T)

####CREATE PLOT FOR FIGURE 2####
# Arrange plots in a 2x2 grid
library(terra)
library(viridis)

par(mfrow = c(2,4),
    mar = c(1,1,1,1),   # shrink inner margins
    oma = c(3,3,8,1),   # shrink outer margins
    #mgp = c(1.5,0.5,0), # tighter axis spacing
    ps = 14,            # global font size
    cex.main = 0.8)

plot(r_mean1, main="Bio01 (°C)",
     col = colorRampPalette(c("white","yellow","orange","red","darkred","black"))(50))
plot(r_mean2, main="Bio12 (mm)",
     col = colorRampPalette(c("white","skyblue","darkblue"))(50))
plot(NPP_mean1, main="NPP (g DM m² yr)",
     col = colorRampPalette(c("yellow","green","forestgreen","darkgreen"))(50))
plot(density_mean, main="Predicted Forager Density (#/km²)",
     col = viridis(50, option="plasma"))

plot(r_sd1, col = colorRampPalette(c("white","yellow","orange","red","darkred","black"))(50))
plot(r_sd2, col = colorRampPalette(c("white","skyblue","darkblue"))(50))
plot(NPP_sd1,  col = colorRampPalette(c("yellow","green","forestgreen","darkgreen"))(50))
plot(density_sd, col = viridis(50, option="plasma"))

mtext("Mean", side = 2, line = 1, outer = TRUE, at = 0.75, cex = 1.2)
mtext("Standard Deviation", side = 2, line = 1, outer = TRUE, at = 0.25, cex = 1.2)


#some useful stats
max(NPP_mean1)
max2 <- rast(density_rasters)

max_vals <- global(max2, "range", na.rm = TRUE)
min(max_vals$min)
max(max_vals$max)

max_vals <- global(NPP, "range", na.rm = TRUE)
min(max_vals$min)
max(max_vals$max)

#####SUM DENSITIES BY LANDMASS TO GENERATE POPULATION SIZE####
patchy_rasters <- lapply(density_rasters, patches)
patchy_rasters_lo <- lapply(density_rasters_lo, patches)
patchy_rasters_hi <- lapply(density_rasters_hi, patches)

writeCDF(rast(patchy_rasters), paste0("Malta_CHELSA_patchy_rasters", ".nc"), overwrite = TRUE) #s
writeCDF(rast(patchy_rasters_lo), paste0("Malta_CHELSA_patchy_rasters", "_lo.nc"), overwrite = TRUE) #s
writeCDF(rast(patchy_rasters_hi), paste0("Malta_CHELSA_patchy_rasters", "_hi.nc"), overwrite = TRUE) #s

density_sums <- list()
for(k in 1:nrow(model_df)){
  
  # Initialize a list to store results
  density_sums_list1 <- list() 
  
  # Loop through each raster layer
  for (i in 1:nlyr(patchy_rasters[[k]])) {
    density_sums_list1[[i]] <- zonal(density_rasters[[k]][[i]], patchy_rasters[[k]][[i]], "sum") #calculates sum of density by patch
    
    # Ensure consistent column names
    colnames(density_sums_list1[[i]]) <- c("Patch_ID", paste0("Sum_Layer_", i))
  }
  
  density_sums[[k]] <- density_sums_list1}

density_sums_lo <- list()
for(k in 1:nrow(model_df)){
  
  # Initialize a list to store results
  density_sums_list1 <- list() 
  
  # Loop through each raster layer
  for (i in 1:nlyr(patchy_rasters_lo[[k]])) {
    density_sums_list1[[i]] <- zonal(density_rasters_lo[[k]][[i]], patchy_rasters_lo[[k]][[i]], "sum") #calculates sum of density by patch
    
    # Ensure consistent column names
    colnames(density_sums_list1[[i]]) <- c("Patch_ID", paste0("Sum_Layer_", i))
  }
  
  density_sums_lo[[k]] <- density_sums_list1}

density_sums_hi <- list()
for(k in 1:nrow(model_df)){
  
  # Initialize a list to store results
  density_sums_list1 <- list() 
  
  # Loop through each raster layer
  for (i in 1:nlyr(patchy_rasters_hi[[k]])) {
    density_sums_list1[[i]] <- zonal(density_rasters_hi[[k]][[i]], patchy_rasters_hi[[k]][[i]], "sum") #calculates sum of density by patch
    
    # Ensure consistent column names
    colnames(density_sums_list1[[i]]) <- c("Patch_ID", paste0("Sum_Layer_", i))
  }
  
  density_sums_hi[[k]] <- density_sums_list1}

density_rasters_x <- list()
density_rasters_x_hi <- list()
density_rasters_x_lo <- list()

for(k in 1:nrow(model_df)){
  # Create an updated raster stack
  density_rasters1x <- patchy_rasters[[k]]  # Copy structure
  
  # Loop through layers
  for (i in 1:nlyr(patchy_rasters[[k]])) {
    # Extract current layer
    layer <- patchy_rasters[[k]][[i]]
    density_sums_df <- density_sums[[k]][[i]]
    
    # Create a lookup table from density_sums_df
    lookup_table <- data.frame(
      patch_value = density_sums_df[[1]],    # Patch identifiers
      density_value = density_sums_df[[2]]# Corresponding density sums
    )
    
    # Apply classification to swap patch values with density values
    density_rasters1x[[i]] <- classify(layer, lookup_table, right = FALSE) #this lets you plot maps showing density by patch, rather than by cell values
  }
  
  density_rasters_x[[k]] <- density_rasters1x
}

for(k in 1:nrow(model_df)){
  # Create an updated raster stack
  density_rasters1x <- patchy_rasters_lo[[k]]  # Copy structure
  
  # Loop through layers
  for (i in 1:nlyr(patchy_rasters_lo[[k]])) {
    # Extract current layer
    layer <- patchy_rasters_lo[[k]][[i]]
    density_sums_df <- density_sums_lo[[k]][[i]]
    
    # Create a lookup table from density_sums_df
    lookup_table <- data.frame(
      patch_value = density_sums_df[[1]],    # Patch identifiers
      density_value = density_sums_df[[2]]# Corresponding density sums
    )
    
    # Apply classification to swap patch values with density values
    density_rasters1x[[i]] <- classify(layer, lookup_table, right = FALSE) #this lets you plot maps showing density by patch, rather than by cell values
  }
  
  density_rasters_x_lo[[k]] <- density_rasters1x
}

for(k in 1:nrow(model_df)){
  # Create an updated raster stack
  density_rasters1x <- patchy_rasters_hi[[k]]  # Copy structure
  
  # Loop through layers
  for (i in 1:nlyr(patchy_rasters_hi[[k]])) {
    # Extract current layer
    layer <- patchy_rasters_hi[[k]][[i]]
    density_sums_df <- density_sums_hi[[k]][[i]]
    
    # Create a lookup table from density_sums_df
    lookup_table <- data.frame(
      patch_value = density_sums_df[[1]],    # Patch identifiers
      density_value = density_sums_df[[2]]# Corresponding density sums
    )
    
    # Apply classification to swap patch values with density values
    density_rasters1x[[i]] <- classify(layer, lookup_table, right = FALSE) #this lets you plot maps showing density by patch, rather than by cell values
  }
  
  density_rasters_x_hi[[k]] <- density_rasters1x
}

library(terra)
library(ggplot2)

density_df <- data.frame(t = t, average_malta = NA, average_comino = NA, average_gozo = NA, average_sicily = NA, 
                         hi_malta = NA, hi_comino = NA, hi_gozo = NA, hi_sicily = NA,
                         lo_malta = NA, lo_comino = NA, lo_gozo = NA, lo_sicily = NA)

patch_values_list <- list()

target_coord <- matrix(c(14.439837, 35.892727), ncol=2)  # Longitude first #target location is the Bizzle

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(patchy_rasters_hi[[k]]))
  density_values1 <- numeric(nlyr(patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  
  patch_values_list[[6]] <- patch_values1
  density_df[,6] <- unlist(density_values1)}

target_coord <- matrix(c(14.336429, 36.011046), ncol=2)  # middle of Comino

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1a <- numeric(nlyr(patchy_rasters_hi[[k]]))
  density_values1a <- numeric(nlyr(patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1a[i] <- patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1a[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[7]] <- patch_values1a
  density_df[,7] <- unlist(density_values1a)}

target_coord <- matrix(c(14.2394401, 36.0463615), ncol=2)  # gozo

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1b <- numeric(nlyr(patchy_rasters_hi[[k]]))
  density_values1b <- numeric(nlyr(patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1b[i] <- patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1b[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[8]] <- patch_values1b
  density_df[,8] <- unlist(density_values1b)}

target_coord <- matrix(c(14.996843, 37.754805), ncol=2) #the crater of Etna

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values2 <- numeric(nlyr(patchy_rasters_hi[[k]]))
  density_values2 <- numeric(nlyr(patchy_rasters_hi[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_hi[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_hi[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values2[i] <- patchy_rasters_hi[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values2[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[9]] <- patch_values2
  density_df[,9] <- unlist(density_values2)}

target_coord <- matrix(c(14.439837, 35.892727), ncol=2)  # Longitude first #target location is the Bizzle

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(patchy_rasters_lo[[k]]))
  density_values1 <- numeric(nlyr(patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[10]] <- patch_values1
  density_df[,10] <- unlist(density_values1)}

target_coord <- matrix(c(14.336429, 36.011046), ncol=2)  # middle of Comino

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1a <- numeric(nlyr(patchy_rasters_lo[[k]]))
  density_values1a <- numeric(nlyr(patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1a[i] <- patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1a[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[11]] <- patch_values1a
  density_df[,11] <- unlist(density_values1a)}

target_coord <- matrix(c(14.2394401, 36.0463615), ncol=2)  # gozo

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1b <- numeric(nlyr(patchy_rasters_lo[[k]]))
  density_values1b <- numeric(nlyr(patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1b[i] <- patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1b[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[12]] <- patch_values1b
  density_df[,12] <- unlist(density_values1b)}

target_coord <- matrix(c(14.996843, 37.754805), ncol=2) #the crater of Etna

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values2 <- numeric(nlyr(patchy_rasters_lo[[k]]))
  density_values2 <- numeric(nlyr(patchy_rasters_lo[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x_lo[[k]])){
    #i <- 1
    density_layer <- density_rasters_x_lo[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values2[i] <- patchy_rasters_lo[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values2[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[13]] <- patch_values2
  density_df[,13] <- unlist(density_values2)}

target_coord <- matrix(c(14.439837, 35.892727), ncol=2)  # Longitude first #target location is the Bizzle

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1 <- numeric(nlyr(patchy_rasters[[k]]))
  density_values1 <- numeric(nlyr(patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    #i <- 1
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1[i] <- patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[2]] <- patch_values1
  density_df[,2] <- unlist(density_values1)}

target_coord <- matrix(c(14.336429, 36.011046), ncol=2)  # middle of Comino

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1a <- numeric(nlyr(patchy_rasters[[k]]))
  density_values1a <- numeric(nlyr(patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    #i <- 1
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1a[i] <- patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1a[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[3]] <- patch_values1a
  density_df[,3] <- unlist(density_values1a)}

target_coord <- matrix(c(14.2394401, 36.0463615), ncol=2)  # gozo

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values1b <- numeric(nlyr(patchy_rasters[[k]]))
  density_values1b <- numeric(nlyr(patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    #i <- 1
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values1b[i] <- patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values1b[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[4]] <- patch_values1b
  density_df[,4] <- unlist(density_values1b)}

target_coord <- matrix(c(14.996843, 37.754805), ncol=2) #the crater of Etna

for(k in 1:nrow(model_df)){
  
  # Initialize vectors
  patch_values2 <- numeric(nlyr(patchy_rasters[[k]]))
  density_values2 <- numeric(nlyr(patchy_rasters[[k]]))
  
  # Extract patch and corresponding density value for each layer
  for (i in 1:nlyr(density_rasters_x[[k]])){
    #i <- 1
    density_layer <- density_rasters_x[[k]][[i]]
    
    # Find nearest cell and patch value
    nearest_cell <- cellFromXY(layer, target_coord)
    patch_values2[i] <- patchy_rasters[[k]][[i]][nearest_cell]
    
    # Lookup corresponding density value in density_sums_df
    density_values2[i] <- density_layer[nearest_cell]
    
  }
  patch_values_list[[5]] <- patch_values2
  density_df[,5] <- unlist(density_values2)}

write.csv(density_df, "density_df_full.csv")

patch_values_df <- do.call(cbind, patch_values_list[-1])
write.csv(patch_values_df, "patch_values_list.csv")

density_df_full <- read.csv("density_df_full.csv")
density_df_full <- density_df_full[,-1]

####AREA CALCULATIONS####

patchy_rasters2 <- rast("Malta_CHELSA_patchy_rasters.nc")

patch_values_list <- read.csv("patch_values_list.csv") 

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
  
p1 <- patchy_rasters2[[i]] == id_malta
p1a <- patchy_rasters2[[i]] == id_comino
p1b <- patchy_rasters2[[i]] == id_gozo
p2 <- patchy_rasters2[[i]] == id_sicily

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

patches_1 <- project(rast(p1), "EPSG:32633")
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

patchy_rasters2 <- rast("Malta_CHELSA_patchy_rasters_hi.nc")
plot(patchy_rasters2[[1]])
crs(patchy_rasters2)

min_dist <- numeric()
total_true_area_malta_hi <- numeric()
total_true_area_comino_hi <- numeric()
total_true_area_gozo_hi <- numeric()
total_true_area_sicily_hi <- numeric()

for(i in 1:221){
  id_malta <- patch_values_list[[6]][[i]]
  id_comino <- patch_values_list[[7]][[i]]
  id_gozo <-patch_values_list[[8]][[i]]
  id_sicily <- patch_values_list[[9]][[i]]
  
  p1 <- patchy_rasters2[[i]] == id_malta
  p1a <- patchy_rasters2[[i]] == id_comino
  p1b <- patchy_rasters2[[i]] == id_gozo
  p2 <- patchy_rasters2[[i]] == id_sicily
  
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

  # terra version
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

patchy_rasters2 <- rast("Malta_CHELSA_patchy_rasters_lo.nc")
plot(patchy_rasters2[[1]])
crs(patchy_rasters2)

min_dist <- numeric()
total_true_area_malta_lo <- numeric()
total_true_area_comino_lo <- numeric()
total_true_area_gozo_lo <- numeric()
total_true_area_sicily_lo <- numeric()

for(i in 1:221){

  id_malta <- patch_values_list[[10]][[i]]
  id_comino <- patch_values_list[[11]][[i]]
  id_gozo <-patch_values_list[[12]][[i]]
  id_sicily <- patch_values_list[[13]][[i]]
  
  p1 <- patchy_rasters2[[i]] == id_malta
  p1a <- patchy_rasters2[[i]] == id_comino
  p1b <- patchy_rasters2[[i]] == id_gozo
  p2 <- patchy_rasters2[[i]] == id_sicily
  
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
  
  # terra version
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

patchy_rasters2 <- rast("Malta_CHELSA_patchy_rasters.nc")

# store largest sea crossing per raster
library(raster)
library(gdistance)


for (i in 1:nlyr(patchy_rasters2)) {
  
  # 1. Project raster
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  
  # 2. Build cost raster: 0 = patch, 1 = sea
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  
  # 3. Define points
  etna_ll   <- vect(matrix(c(14.996843, 37.754805), ncol=2), crs="EPSG:4326")
  bizzle_ll <- vect(matrix(c(14.439837, 35.892727), ncol=2), crs="EPSG:4326")
  etna   <- project(etna_ll, patches_r)
  bizzle <- project(bizzle_ll, patches_r)
  
  # Convert to SpatialPoints for gdistance
  etna_sp   <- as(etna, "Spatial")
  bizzle_sp <- as(bizzle, "Spatial")
  
  # 4. Build transition object with queen adjacency
  tr <- transition(raster(cost_r), function(x) 1 / mean(x), directions = 8)
  tr_corr <- geoCorrection(tr, type = "c")
  
  # 5. Compute shortest path
  path_line <- shortestPath(tr_corr, etna_sp, bizzle_sp, output = "SpatialLines")
  line_list[[i]] <- path_line
  
  # 6. Compute cost distance
  cost_dist[[i]] <- costDistance(tr_corr, etna_sp, bizzle_sp)

}

#Max Crossing Distance#

max_dist_list <- list()
max_dist <- numeric()
total_list <- list()
total_dist <- numeric()

for (i in 1:nlyr(patchy_rasters2)) {

patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")

# 2. Build cost raster: 0 = patch, 1 = sea
cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
path_vect <- vect(line_list[[i]])
crs(path_vect) <- crs(cost_r)

# Identify contiguous patches (queen adjacency by default)
patch_ids <- patches(patches_r, directions = 8)

# Convert to polygons, dissolving by patch ID
patch_polys <- as.polygons(patch_ids, dissolve = TRUE, values = TRUE)

patch_polys$cost_val <- terra::extract(cost_r, patch_polys, fun = mean, na.rm = TRUE)[,2]

path_segments <- terra::intersect(patch_polys, path_vect)

x <- path_vect-path_segments

#lines(x, col="blue")

max_dist_list[[i]] <- perim(disagg(x))
max_dist[i] <- max(perim(disagg(x)))
}

line_list_hi <- list()
cost_dist_hi <- numeric()
max_cross_hi <- numeric()
sea_list_hi <- list()


patchy_rasters2 <- rast("Malta_CHELSA_patchy_rasters_hi.nc")

# store largest sea crossing per raster
library(raster)
library(gdistance)

for (i in 1:nlyr(patchy_rasters2)) {
  
  # 1. Project raster
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  
  # 2. Build cost raster: 0 = patch, 1 = sea
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  
  # 3. Define points
  etna_ll   <- vect(matrix(c(14.996843, 37.754805), ncol=2), crs="EPSG:4326")
  bizzle_ll <- vect(matrix(c(14.439837, 35.892727), ncol=2), crs="EPSG:4326")
  etna   <- project(etna_ll, patches_r)
  bizzle <- project(bizzle_ll, patches_r)
  
  # Convert to SpatialPoints for gdistance
  etna_sp   <- as(etna, "Spatial")
  bizzle_sp <- as(bizzle, "Spatial")
  
  # 4. Build transition object with queen adjacency
  tr <- transition(raster(cost_r), function(x) 1 / mean(x), directions = 8)
  tr_corr <- geoCorrection(tr, type = "c")
  
  # 5. Compute shortest path
  path_line <- shortestPath(tr_corr, etna_sp, bizzle_sp, output = "SpatialLines")
  line_list_hi[[i]] <- path_line
  
  # 6. Compute cost distance
  cost_dist_hi[[i]] <- costDistance(tr_corr, etna_sp, bizzle_sp)
  
}

#Max Crossing Distance#

max_dist_list_hi <- list()
max_dist_hi <- numeric()
total_list_hi <- list()
total_dist_hi <- numeric()

for (i in 1:nlyr(patchy_rasters2)) {
  
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  
  # 2. Build cost raster: 0 = patch, 1 = sea
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  path_vect <- vect(line_list[[i]])
  crs(path_vect) <- crs(cost_r)
  
  # Identify contiguous patches (queen adjacency by default)
  patch_ids <- patches(patches_r, directions = 8)
  
  # Convert to polygons, dissolving by patch ID
  patch_polys <- as.polygons(patch_ids, dissolve = TRUE, values = TRUE)
  
  patch_polys$cost_val <- terra::extract(cost_r, patch_polys, fun = mean, na.rm = TRUE)[,2]
  
  path_segments <- terra::intersect(patch_polys, path_vect)
  
  x <- path_vect-path_segments
  
  #lines(x, col="blue")
  
  max_dist_list_hi[[i]] <- perim(disagg(x))
  max_dist_hi[i] <- max(perim(disagg(x)))
}

line_list_lo <- list()
cost_dist_lo <- numeric()
max_cross_lo <- numeric()
sea_list_lo <- list()

patchy_rasters2 <- rast("Malta_CHELSA_patchy_rasters_lo.nc")

# store largest sea crossing per raster
library(raster)
library(gdistance)


for (i in 1:nlyr(patchy_rasters2)) {
  
  # 1. Project raster
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  
  # 2. Build cost raster: 0 = patch, 1 = sea
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  
  # 3. Define points
  etna_ll   <- vect(matrix(c(14.996843, 37.754805), ncol=2), crs="EPSG:4326")
  bizzle_ll <- vect(matrix(c(14.439837, 35.892727), ncol=2), crs="EPSG:4326")
  etna   <- project(etna_ll, patches_r)
  bizzle <- project(bizzle_ll, patches_r)
  
  # Convert to SpatialPoints for gdistance
  etna_sp   <- as(etna, "Spatial")
  bizzle_sp <- as(bizzle, "Spatial")
  
  # 4. Build transition object with queen adjacency
  tr <- transition(raster(cost_r), function(x) 1 / mean(x), directions = 8)
  tr_corr <- geoCorrection(tr, type = "c")
  
  # 5. Compute shortest path
  path_line <- shortestPath(tr_corr, etna_sp, bizzle_sp, output = "SpatialLines")
  line_list_lo[[i]] <- path_line
  
  # 6. Compute cost distance
  cost_dist_lo[[i]] <- costDistance(tr_corr, etna_sp, bizzle_sp)
  
}


#Max Crossing Distance#

max_dist_list_lo <- list()
max_dist_lo <- numeric()
total_list_lo <- list()
total_dist_lo <- numeric()

for (i in 1:nlyr(patchy_rasters2)) {
  
  patches_r <- project(patchy_rasters2[[i]], "EPSG:32633")
  
  # 2. Build cost raster: 0 = patch, 1 = sea
  cost_r <- ifel(is.na(patches_r), 1, 0.00000001)
  path_vect <- vect(line_list[[i]])
  crs(path_vect) <- crs(cost_r)
  
  # Identify contiguous patches (queen adjacency by default)
  patch_ids <- patches(patches_r, directions = 8)
  
  # Convert to polygons, dissolving by patch ID
  patch_polys <- as.polygons(patch_ids, dissolve = TRUE, values = TRUE)
  
  patch_polys$cost_val <- terra::extract(cost_r, patch_polys, fun = mean, na.rm = TRUE)[,2]
  
  path_segments <- terra::intersect(patch_polys, path_vect)
  
  x <- path_vect-path_segments
  
  #lines(x, col="blue")
  
  max_dist_list_lo[[i]] <- perim(disagg(x))
  max_dist_lo[i] <- max(perim(disagg(x)))
}


####PLOT FOR FIGURE 3####

combine_values <- function(malta, gozo, comino) {
  gc <- ifelse(gozo == comino, gozo, gozo + comino)
  total <- ifelse(malta == gc, gc, malta + gc)
  return(total)
}

density_df_full <- density_df_full %>%
  mutate(
    average_maltese = mapply(combine_values, average_malta, average_gozo, average_comino),
    hi_maltese      = mapply(combine_values, hi_malta, hi_gozo, hi_comino),
    lo_maltese      = mapply(combine_values, lo_malta, lo_gozo, lo_comino)
  )

# Find all rows where Maltese != Sicily
distinct_rows <- density_df_full %>%
  filter(average_maltese != average_sicily |
           hi_maltese != hi_sicily |
           lo_maltese != lo_sicily)

# Get the last distinct t based on your data order
last_distinct_t <- tail(distinct_rows$t, 1)

density_df_full <- density_df_full %>%
  mutate(
    average_maltese_plot = ifelse(t > last_distinct_t, average_maltese, NA),
    hi_maltese_plot      = ifelse(t > last_distinct_t, hi_maltese, NA),
    lo_maltese_plot      = ifelse(t > last_distinct_t, lo_maltese, NA)
  )


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

l1 <- ggplot(density_df_full, aes(x = t)) +
  # Maltese Islands
  geom_ribbon(aes(ymin = lo_maltese_plot, ymax = hi_maltese_plot, fill = "Maltese Islands"), alpha = 0.2) +
  geom_line(aes(y = average_maltese_plot, color = "Maltese Islands"), size = 1) +
  
  # Sicily
  geom_ribbon(aes(ymin = lo_sicily, ymax = hi_sicily, fill = "Sicily"), alpha = 0.2) +
  geom_line(aes(y = average_sicily, color = "Sicily"), size = 1) +
  
  scale_y_log10(
    limits = c(10, 10000),                     # force axis range
    breaks = c(10, 100, 1000, 10000),          # show only these tick marks
    labels = scales::comma                     # format nicely (10, 100, 1,000, 10,000)
  ) +
  scale_color_manual(values = c("Maltese Islands" = "red", "Sicily" = "blue")) +
  scale_fill_manual(values = c("Maltese Islands" = "red", "Sicily" = "blue")) +
  labs(x = "Time", y = "Population size",
       title = "Population size",
       color = "Region", fill = "Region") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(),   # keep major gridlines at 10, 100, 1000, 10000
    panel.grid.minor = element_blank(),  # remove minor gridlines
    axis.ticks.y = element_line(),
    axis.ticks.x = element_line(),
    legend.position = c(0.95, 0.05),   # bottom right inside plot
    legend.justification = c("right", "bottom")
  )

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


l1+l2
