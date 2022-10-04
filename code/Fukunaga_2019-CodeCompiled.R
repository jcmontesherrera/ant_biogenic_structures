#################### R code for processing DEMs ####################
#
# Integrating Three-Dimensional Benthic Habitat Characterization Techniques
# into Ecological Monitoring of coral reefs
# A. Fukunaga, J.H.R. Burns, B.K. Craig and R.K. Kosaki 
# Journal of Marine Science and Engineering (2019)

####################################################################

# JCMH version (2022): Towards reproducible analysis of benthos structural complexity: 
# A case study on Antarctic polychaete reefs using action cameras and remotely operated vehicles 

# Aggregation values changed to 16 and 32 cm

library(raster)  # need v. 2.5-16 or ealier or uncomment the two lines to convert NaN to NA
library(rgeos)
library(ggplot2)

terrain_fun <- function(data, cell_size) {
  
  d0 <- as.matrix(data)
  if (ncol(d0) == 1 | nrow(d0) == 1) {
    terrain_list <- list(NA, NA, NA)
    names(terrain_list) <- c("mean_slope", "mean_profile_curvature", "mean_plan_curvature")
  } else {
    da <-  cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                 rbind(matrix(NA, nrow = 1, ncol = ncol(d0) - 1), 
                       matrix(d0[1 : (nrow(d0) - 1), 1 : (ncol(d0) - 1)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1)))
    db <- rbind(matrix(NA, nrow = 1, ncol = ncol(d0)), 
                matrix(d0[1 : (nrow(d0) - 1), ], nrow = nrow(d0) - 1, ncol = ncol(d0)))
    dc <- cbind(rbind(matrix(NA, nrow = 1, ncol = ncol(d0) - 1), 
                      matrix(d0[1 : (nrow(d0) - 1), 2 : ncol(d0)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1)),
                matrix(NA, nrow = nrow(d0), ncol = 1))
    dd <- cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                matrix(d0[, 1 : (ncol(d0) - 1)], nrow = nrow(d0), ncol = ncol(d0) -1))
    df <- cbind(matrix(d0[, 2 : ncol(d0)], nrow = nrow(d0), ncol = ncol(d0) - 1), 
                matrix(NA, nrow = nrow(d0), ncol = 1))
    dg <- cbind(matrix(NA, nrow = nrow(d0), ncol = 1), 
                rbind(matrix(d0[2 : nrow(d0), 1 : (ncol(d0) - 1)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1), 
                      matrix(NA, nrow = 1, ncol = ncol(d0) - 1)))
    dh <- rbind(matrix(d0[2 : nrow(d0), ], nrow = nrow(d0) - 1, ncol = ncol(d0)), 
                matrix(NA, nrow = 1, ncol = ncol(d0)))
    di <- cbind(rbind(matrix(d0[2 : nrow(d0), 2 : ncol(d0)], nrow = nrow(d0) - 1, ncol = ncol(d0) - 1), 
                      matrix(NA, nrow = 1, ncol = (ncol(d0) - 1))), 
                matrix(NA, nrow = nrow(d0), ncol = 1))
    
    ## slope
    
    x_rate <- ((dc + (2 * df) + di) - (da + (2 * dd) + dg)) / (8 * cell_size)
    y_rate <- ((dg + (2 * dh) + di) - (da + (2 * db) + dc)) / (8 * cell_size)
    
    slope_degrees <- atan(sqrt(x_rate ^ 2 + y_rate ^ 2)) * 57.29578
    mean_slope <- mean(slope_degrees, na.rm = TRUE)
    
    ## curvature
    
    coefD <- (((dd + df) / 2) - d0) / (cell_size ^ 2)
    coefE <- (((db + dh) / 2) - d0) / (cell_size ^ 2)
    coefF <- (dc + dg -da - di) / (4 * (cell_size ^ 2))
    coefG <- (df - dd) / 2 * cell_size
    coefH <- (db - dh) / 2 * cell_size
    
    prof_curv <- -2 * ((coefD * coefG ^ 2 + coefE * coefH ^ 2 + coefF * coefG * coefH) / 
                         (coefG ^ 2 + coefH ^ 2))
    mean_prof_curv <- mean(prof_curv, na.rm = TRUE)
    
    plan_curv <- 2 * ((coefD * coefH ^ 2 + coefE * coefG ^ 2 - coefF * coefG * coefH) / 
                        (coefG ^ 2 + coefH ^ 2))
    mean_plan_curv <- mean(plan_curv, na.rm = TRUE)
    
    ## output list
    
    terrain_list <- list(mean_slope, mean_prof_curv, mean_plan_curv)
    names(terrain_list) <- c("mean_slope", "mean_profile_curvature", "mean_plan_curvature")
  }
  
  return(terrain_list)
  
}

#source("Text_S2.R")  # or run scripts in Text_S2.R
# install.pacakges("rgdal")   # need to have the rgdal pacakge installed

########## analysis prep #########
resolution <- 0.01  # enter DEM resolution in meter

### initialize output data frame
hab_data <- data.frame(file_name = character(0), 
                       orig_area = numeric(0), 
                       surface16 = numeric(0), planerS32 = numeric(0), 
                       max_h = numeric(0), min_h = numeric(0), 
                       fd16 = numeric(0), fd32 = numeric(0), 
                       surface_complexity = numeric(0), 
                       mean_slope = numeric(0), 
                       mean_profile_curvature = numeric(0), 
                       mean_plan_curvature = numeric(0)
)

########## process files #########

dir = "ant_biogenic_structures\\data\\quadrats\\"
setwd(dir)

file_names <- list.files(path = ".", pattern = "DEM")

## enter file names
files <- c(file_names, sep = "")

for (k in 1:length(files)) {
  
  ras <- raster(files[k])
  ras_g <- as(ras, "SpatialGridDataFrame")
  orig_area <- sum(!is.na(ras_g@data)) * resolution ^ 2
  
  print(file_names[k])
  
  # Aggregate function - "Splits the data into subsets, computes summary statistics for each, and returns the result in a convenient form."
  
  #modified to 16 due to total size of raster
  extent16 <- aggregate(ras, fac = 16, fun = mean, expand = FALSE, na.rm = FALSE)
  extent16_uniform <- extent16 > -Inf
  extent16_polygon <- rasterToPolygons(extent16_uniform, dissolve = TRUE)
  
  #modified to 32
  extent32 <- aggregate(ras, fac = 32, fun = mean, expand = FALSE, na.rm = FALSE)
  extent32_uniform <- extent32 > -Inf
  extent32_polygon <- rasterToPolygons(extent32_uniform, dissolve = TRUE)
  
  ras_clip_16 <- mask(crop(ras, extent(extent16_polygon)), extent16_polygon)
  ras_clip_32 <- mask(crop(ras, extent(extent32_polygon)), extent32_polygon)
  
  # # uncomment to plot raster
  #plot(ras)
  #plot(ras_clip_16)
  #plot(ras_clip_32)
  
  dat16 <- data.frame(fac = c(1, 2, 4, 8, 16), 
                     s_area = NA, area = NA, 
                     max_height = NA, min_height = NA, 
                     mean_slope = NA,  
                     mean_profile_curvature = NA, mean_plan_curvature = NA)
  dat32 <- data.frame(fac = c(1, 2, 4, 8, 16, 32), 
                      s_area = NA, area = NA, 
                      max_height = NA, min_height = NA, 
                      mean_slope = NA, 
                      mean_profile_curvature = NA, mean_plan_curvature = NA)
  
  dat <- list(dat16, dat32)
  
  for (j in 1:length(dat)) {
    
    if (j == 1) {
      for (i in 1:nrow(dat[[j]])) {
        
        print(dat[[j]]$fac[i])
        
        if (dat[[j]]$fac[i] == 1) {
          reef <- ras_clip_16
        } else {
          reef <- aggregate(ras_clip_16, fac = dat[[j]]$fac[i], fun = mean, 
                            expand = FALSE, na.rm = FALSE)
        }
        reef_g <- as(reef, "SpatialGridDataFrame")
        reef_g@data[, 1][is.nan(reef_g@data[, 1])] <- NA   # uncomment if using raster 2.6-7 or later
        terrain_res <- terrain_fun(reef, resolution)
        dat[[j]]$s_area[i] <- surfaceArea(reef_g)
        dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (dat[[j]]$fac[i] * resolution) ^ 2
        dat[[j]]$max_height[i] <- max(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$min_height[i] <- min(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$mean_slope[i] <- terrain_res$mean_slope
        dat[[j]]$mean_profile_curvature[i] <- terrain_res$mean_profile_curvature
        dat[[j]]$mean_plan_curvature[i] <- terrain_res$mean_plan_curvature
      }
    } else if (j == 2) {
      for (i in 1:nrow(dat[[j]])) {
        
        print(dat[[j]]$fac[i])
        
        if (dat[[j]]$fac[i] == 1) {
          reef <- ras_clip_32
        } else {
          reef <- aggregate(ras_clip_32, fac = dat[[j]]$fac[i], fun = mean, 
                            expand = FALSE, na.rm = FALSE)
        }
        reef_g <- as(reef, "SpatialGridDataFrame")
        reef_g@data[, 1][is.nan(reef_g@data[, 1])] <- NA   # uncomment if using raster 2.6-7 or later
        terrain_res <- terrain_fun(reef, resolution)
        dat[[j]]$s_area[i] <- surfaceArea(reef_g)
        dat[[j]]$area[i] <- sum(!is.na(reef_g@data)) * (dat[[j]]$fac[i] * resolution) ^ 2
        dat[[j]]$max_height[i] <- max(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$min_height[i] <- min(reef_g@data[[1]], na.rm = TRUE)
        dat[[j]]$mean_slope[i] <- terrain_res$mean_slope
        dat[[j]]$mean_profile_curvature[i] <- terrain_res$mean_profile_curvature
        dat[[j]]$mean_plan_curvature[i] <- terrain_res$mean_plan_curvature  }
    } else {
      break
    }
  }
  
  p16 <- ggplot(dat[[1]], aes(x = log(fac * 0.01), y = log(s_area))) +
    geom_point() +
    geom_smooth(method = "lm", se = F) +
    xlim(c(-5, 1)) + ylim(c(4, 6)) + 
    labs(x = expression(paste("log(", delta, ")")),
         y = expression(paste("logS(", delta, ")"))) +
    ggtitle(file_names[k])
  theme_bw()
  
  p32 <- ggplot(dat[[2]], aes(x = log(fac * 0.01), y = log(s_area))) +
    geom_point() +
    geom_smooth(method = "lm", se = F) + 
    xlim(c(-5, 1)) + ylim(c(4, 6)) + 
    labs(x = expression(paste("log(", delta, ")")),
         y = expression(paste("logS(", delta, ")"))) +
    ggtitle(file_names[k])
  theme_bw()
  
  #plot(p16)
  #plot(p32)
  
  #d16 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[1]])
  #slope16 <- coef(d16)[[2]]
  #fd16 <- 2 - slope16
  
  d16 <- lm(log(s_area/area) ~ log(fac * resolution), data = dat[[1]])  # using log(s_area/area) to account for changes in are due to raster aggregate function
  s_coef <- coef(d16)[[2]]
  fd16 <- 2 - s_coef  # fractal dimension at 1-16 cm resolution range
  
  d32 <- lm(log(s_area/area) ~ log(fac * 0.01), data = dat[[2]])
  slope32 <- coef(d32)[[2]]
  fd32 <- 2 - slope32
  
  surface16 <- dat[[1]]$s_area[1]
  surface32 <- dat[[2]]$s_area[1]
  
  planerS16 <- dat[[1]]$area[1]
  planerS32 <- dat[[2]]$area[1]
  
  surface_complexity <- dat[[1]]$s_area[1]/dat[[1]]$area[1]
  surface_complexity2 <- dat[[1]]$s_area[2]/dat[[1]]$area[2]
  
  max_h <- dat[[1]]$max_height[1]
  min_h <- dat[[1]]$min_height[1]
  
  mean_slope <- dat[[1]]$mean_slope[1]
  
  mean_profile_curvature <- dat[[1]]$mean_profile_curvature[1]
  mean_plan_curvature <- dat[[1]]$mean_plan_curvature[1]
  
  temp <- data.frame(file_name = file_names[k], 
                     orig_area = orig_area, 
                     surface16 = surface16, planerS16 = planerS16,
                     surface32 = surface32, planerS32 = planerS32, 
                     max_h = max_h, min_h = min_h, 
                     fd16 = fd16, fd32 = fd32, 
                     surface_complexity = surface_complexity,
                     mean_slope = mean_slope, 
                     mean_profile_curvature = mean_profile_curvature, 
                     mean_plan_curvature = mean_plan_curvature)
  
  hab_data <- rbind(hab_data, temp)
  
}

hab_data

########## export output data frame #########
write.csv(hab_data, "habitat_complexity_16-32.csv", row.names = FALSE)

rm(list = ls())  # cleanup
