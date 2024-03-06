################################################################################
#######         Calculate regional rates of glacier lake growth           ######
#######                                                                   ######
#######                            by Georg Veh                           ######
#######                checked and comments added March 05, 2024          ######
################################################################################

# DISCLAIMER: This script is meant to make our findings reproducible. However,
# we do not redistribute the original lake inventories, which might be subject
# to differing licenses. Therefore, please contact the authors of the underlying
# studies to run this script.

# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

# Important: we run brms models using the cmdstan backend. This needs to be
# installed separately using the instructions given here: 
# https://mc-stan.org/cmdstanr/ 

require(tidyverse)
require(tidybayes)
require(readODS)
require(sf)
require(gfcanalysis)
require(pbapply)
require(parallel)
require(brms)

# Set YOUR working directory folder where to find all lake inventories, necessary 
# to run this script. Change the location appropriately.

setwd("D:/data/BoxUP/Published_lake_databases")

################################################################################
### Split lake inventories to the extent of the study regions ##################

# We split each lake inventory according the outlines of the study regions. 
# The lakes need to be within a 5-km buffer around glaciers. To this end,
# we collect the glacier buffer we generated previously.

buffer.dissolve <- readRDS("D:/data/BoxUP/Work 2022/GLOFsize/dissolved_buffer.RDS")

# We then create inventories where the individual inventories should be deposited.

sapply(buffer.dissolve$FULL_NAME, function (x) {
  
  dir.create(path = paste0("D:/data/BoxUP/Published_lake_databases/Lakes_split/", 
                           x)) 
})

# We list the name of the folders, as we will drop the inventories there.

out.folder.list <- paste0("D:/data/BoxUP/Published_lake_databases/Lakes_split/", 
                          buffer.dissolve$FULL_NAME)

# We read the table including all lake outlines, either as point or polygon
# features. The data come as ESRI shapefiles (*.shp) or geopackages (*.gpkg).

tab <- read_ods("Glacier_lakes_global.ods", na = "NA") %>% 
  as_tibble() %>%
  filter(Feature_type == "Polygon",
         Data_format == "shp" | Data_format == "gpkg",
         !is.na(Year_of_inventory),
         !is.na(File),
         !is.na(Underlying_data)) %>%
  group_by(Authors) %>%
  mutate(nobs = n()) %>%
  mutate(Year_of_inventory = as.numeric(Year_of_inventory))

# We generate a look up table of all sensorsand their spatial resolution 
# of a pixel used to map glacier lakes. This resolution will be used as the
# uncertainty in lake area.

sensor.lookup <- tribble(~sensor, ~resolution,
                         "Aerial images", 1,
                         "ASTER", 15,
                         "ALOS", 30,
                         "Corona KH-4", 2.75,
                         "Gaofen-1", 2,
                         "GLAD", 30,
                         "GoogleEarth", 1,
                         "Hexagon KH-9", 9,
                         "JRC", 30,           
                         "Landsat", 30, 
                         "Landsat 4", 30,  
                         "Landsat 5", 30,
                         "Landsat 7", 30,     
                         "Landsat 8", 30,
                         "Landsat MSS", 60,
                         "PlanetScope", 3.125,
                         "Pleiades 1", 2.8,
                         "Sentinel-1", 10,
                         "Sentinel-2", 10,
                         "LIS-III", 23.5,
                         "LISS-IV", 5,
                         "Cartosat-1", 2.5,
                         "Cartosat-2A", 0.8,
                         "SPOT 4", 20,
                         "SPOT 5", 10,
                         "ESRI World Imagery", 0.6,
                         "Rapid Eye", 6.5)

# We now iterate over each entry in the glacier lake database.

for (i in 1:nrow(tab)) {
  
  # We read the vector file, repair its geometry, and project the data to
  # WGS84 (EPSG code: 4326).
  
  fname <- paste0("./", tab$Region[i], "/", tab$Authors[i], "/", tab$File[i])
  
  Author <- tab$Authors[i]
  
  sf_use_s2(TRUE)
  
  shp <- st_read(fname) %>% 
    st_zm(drop = TRUE, what = "ZM") %>%
    st_make_valid() %>%
    st_transform(4326) 
  
  # We ensure that each file is a MULTIPOLYGON.
  
  if (st_geometry_type(shp, by_geometry = F) != "MULTIPOLYGON") {
    
    shp <- st_cast(x = shp , 
                   to = "MULTIPOLYGON")  
    
  }
  
  # We delete, if available, any column that include the year, area,
  # latitude and longitude of the lakes. Instead, we add new columns, and 
  # also add the minimum mapping unit.
  
  shp[ ,grepl(pattern = "year", colnames(shp), ignore.case = T)] <- NULL
  shp[ ,grepl(pattern = "area", colnames(shp), ignore.case = T)] <- NULL
  shp[ ,grepl(pattern = "lat", colnames(shp), ignore.case = T)] <- NULL
  shp[ ,grepl(pattern = "lon", colnames(shp), ignore.case = T)] <- NULL
  
  shp <- shp %>% 
    mutate(Year = tab$Year_of_inventory[i],
           MMU = tab$`Min_lake_area [m²]`[i])
  
  # We delete columns that have non-standard characters, as those are difficult
  # to be machine readable.
  
  shp <- shp[, grep("[^\x20-\x7F]", colnames(shp), invert = T)]
  
  # We then try to find out, which sensors and resolution were used to map the  
  # lakes in the original publication.
  
  sensors.used <- str_split(tab$Underlying_data[i], pattern = ", ")[[1]]
  
  res.used <- sapply(sensors.used, function (x) {
    
    res <- sensor.lookup$resolution[grep(paste0("^", x, "$"), sensor.lookup$sensor)]
    
  })
  
  # We use the minimum resolution to draw buffers of uncertainty around the lake.
  
  min.res <- min(unlist(res.used), na.rm = T)
  
  sf_use_s2(FALSE)
  
  # The next step is important: we find out, which lakes intersect with the RGI 
  # buffer. We only keep lakes that are within a distance of 5 km around the
  # glaciers in the RGI V6.0 to warrant comparability across inventories.
  
  int <- st_intersects(buffer.dissolve, shp) 
  
  # We loop over RGI buffers that intersect with the lake inventories. If
  # no lake intersects with the buffer, we stop the loop.
  
  rgi.int <- which(sapply(int, function (x) !identical(x, integer(0))))
  
  if(length(rgi.int) == 0) next
  
  for (j in rgi.int) {
    
    regional.int <- st_intersects(buffer.dissolve[j, ], shp, sparse = F)[1, ]
    
    # We obtain the coordinate from point inside the polygon. This is
    # more or less equivalent to the centroid.
    
    reg.shp <- shp[regional.int, ]
    
    reg.coords <- reg.shp %>%
      st_point_on_surface()
    
    # We reproject each lake to local UTM zone and calculate its area 
    # We use a cluster environment to speed up this process.
    # Adapt to the number of available cores on YOUR machine. The rule of thumb
    # is number of cores minus one, in our case 8 - 1 = 7.
    
    cl <- makeCluster(7) 
    
    # Export R-packages to the cluster, which will be used in the apply loop. Make sure 
    # you installed these packages before.
    
    clusterEvalQ(cl = cl, c(library("sf"),
                            library("stringr"),
                            library("tibble"),
                            library("tidyverse"),
                            library("gfcanalysis")))
    
    # Export the file list of the glacier shapefiles and the extent of the study 
    # regions to the clusters.
    
    clusterExport(cl = cl, list("reg.shp", "min.res", "Author"))
    
    lake.areas.list <- pblapply(1:nrow(reg.shp), cl = cl, function (l) { 
      
      sf_use_s2(FALSE)
      
      lake <- reg.shp[l, ]
      
      # We identify the local UTM zone of the lake using the function utm_zone()
      # This function still works although rgeos has been retired.
      
      xy <- lake %>% 
        st_centroid() %>% 
        st_coordinates()
      
      loc.utm <- utm_zone(x = xy[1], y = xy[2], proj4string = T)
      
      # We reproject the lake to the local UTM zone and calculate the lake area.
      
      lake.utm <- lake %>%
        st_transform(parse_number(loc.utm))
      
      lake.tib <- tibble(AreaUTM = st_area(lake.utm))
      
      # We draw buffers of half and one pixel uncertainties.
      
      half.pix.inner <- st_buffer(lake.utm, dist = -0.5* min.res, 
                                  singleSide = T) %>% st_area()
      
      one.pix.inner  <- st_buffer(lake.utm,  dist = -1* min.res, 
                                  singleSide = T) %>% st_area()
      
      half.pix.outer <- st_buffer(lake.utm, dist = 0.5* min.res, 
                                  singleSide = T) %>% st_area()
      
      one.pix.outer  <- st_buffer(lake.utm, dist = min.res, 
                                  singleSide = T) %>% st_area()
      
      # Add this information to the tibble generated above.
      
      lake.tib <- lake.tib %>%
        mutate(min_05_px = half.pix.inner,
               min_1_px  = one.pix.inner,
               max_05_px = half.pix.outer,
               max_1_px = one.pix.outer, 
               Lat = xy[, "Y"],
               Lon = xy[, "X"],
               Study = Author)
      
      return(lake.tib)
      
    }) 
    
    # We close the cluster environment and free up some memory.
    
    stopCluster(cl)
    gc()
    
    # We bind all individual tibbles into one long tibble.
    
    reg.shp <- reg.shp %>% 
      mutate(do.call(rbind, lake.areas.list))
    
    # Finally, we write that information into a separate directory that
    # we created before we entered the loop.
    
    out.dir <- paste0(out.folder.list[j], "/", Author, "/")
    
    # We write that vector file to disk in the same format as we had read it.
    
    if(!dir.exists(out.dir)) dir.create(out.dir)
    
    st_write(reg.shp,
             paste0(out.dir, "UTMArea_", basename(tab$File[i])),
             delete_layer = T,
             delete_dsn = T)
    
  }
  
  sf_use_s2(TRUE)
  
  message(paste0("Region: ", tab$Region[i], 
                 ", Author: ", tab$Authors[i], 
                 ", File: " , tab$File[i], 
                 " processed (", i, " out of ", nrow(tab), ")"))
  
}


################################################################################
### Calculate the rate of change in total lake area between two time steps #####

# We navigate to the folder where we had stored the lake inventories that we
# had split to the extent of the study region above.

setwd("D:/data/BoxUP/Published_lake_databases/Lakes_split")

# We list all directories with lake inventories.

shp.dirs <- list.dirs(full.names = T, recursive = F) 

# We iterate over these directories, again in parallel, to obtain the rate of 
# change in lake area.

cl <- makeCluster(6) 

# Export R-packages to the cluster, which will be used in the apply loop. 
# Make sure you installed these packages before.

clusterEvalQ(cl = cl, c(library("sf"),
                        library("stringr"),
                        library("tibble"),
                        library("tidyverse")))

# Export the file list of the directory list that contains the lake shapefiles 
# from different regions to the clusters.

clusterExport(cl = cl, list("shp.dirs", "tab"))

lakes.per.region <- pblapply(shp.dirs, cl = cl, function(reg) {
  
  # List the shapefiles in the directory.
  
  reg.shp <- list.files(pattern      = ".shp$", 
                        path         = reg, 
                        full.names   = T,
                        include.dirs = T, 
                        recursive    = T)
  
  if (length(reg.shp) == 0) return()
  
  reg.shp <- reg.shp[grep("Wangchuk_andBolch_2020", reg.shp, invert = T)]
  
  # We loop over every shapefile in that folder and simulate the size of 
  # each glacier lake based on the uncertainty from the sensor resolution.
  
  la <- lapply(reg.shp, function (x) {
    
    # We read the shapefile to memory.
    
    shp <- st_read(x) %>%
      st_drop_geometry()
    
    # The column names were abbreviated in some cases.
    # Those need to be renamed.
    
    if (any(c("mn_05_p", "mn_1_px" ,"mx_05_p", "mx_1_px") %in% colnames(shp))) { 
      
      shp <- shp %>% 
        rename(min_05_px  = mn_05_p,  
               min_1_px = mn_1_px,  
               max_05_px = mx_05_p,   
               max_1_px = mx_1_px) 
      
    }
    
    # We join the info from the glacial lake overview table with the 
    # attribute table of the shapefile. 
    
    shp <- shp %>%
      transmute(AreaUTM, 
                MMU, 
                Year,
                min_05_px, 
                min_1_px, 
                max_05_px, 
                max_1_px) %>%
      mutate(Study = basename(dirname(x))) %>%
      left_join( x = ., 
                 y = tab %>% dplyr::select(Authors,
                                           Year_of_inventory,
                                           Year_Min,
                                           Year_Max), 
                 by = c("Study" = "Authors",
                        "Year" = "Year_of_inventory"))
    
    # For each lake, we randomly draw 1000 lake areas assuming a mapping error 
    # between -0.5 and 0.5 pixels of the original sensor resolution.
    
    size.sim.05 <- list()
    
    for(i in 1:nrow(shp)) {
      
      size.sim.05[[i]] <-  runif(1000, shp$min_05_px[i], shp$max_05_px[i])
      
    }
    
    # We sum all lake areas to obtain the total area within this region.
    
    size.sim.05.vec <- do.call(rbind, size.sim.05) %>% colSums()
    
    # Same as above: For each lake, randomly draw 1000 lake areas assuming a 
    # mapping error between -1 and 1 pixels of the original sensor resolution.
    
    size.sim.1 <- list()
    
    for(i in 1:nrow(shp)) {
      
      size.sim.1[[i]] <-  runif(1000, shp$min_1_px[i], shp$max_1_px[i])
      
    }
    
    # Sum all lake areas to obtain the total area within this region.
    
    size.sim.1.vec <- do.call(rbind, size.sim.1) %>% colSums()
    
    # Return this output.
    
    return(tibble(area_05_px = size.sim.05.vec,
                  area_1_px  = size.sim.1.vec,
                  Study = unique(shp$Study),
                  Year_Min = unique(shp$Year_Min),
                  Year_Max = unique(shp$Year_Max),
                  MMU = unique(shp$MMU),
                  Year = unique(shp$Year)))
    
  })
  
  # Eventually, we obtain 1000 simulations of total lake area for each inventory.
  # We combine this these simulations to one long tibble.
  
  la <- do.call(rbind, la) %>% 
    as_tibble() 
  
  # For every inventory with two or more time steps, we would like to know the
  # rate of change.
  
  Study.list <- list()
  
  for (i in unique(la$Study)) {
    
    # Sort the inventories by study, and then by year, starting from the
    # oldest inventory to the youngest.
    
    study <- la[la$Study == i, ] %>% 
      ungroup() %>%
      arrange(Year)
    
    uni.years <- unique(study$Year)
    
    if (length(uni.years) < 2) next
    
    # Generate a matrix (converted to a tibble) that stores the information
    # on lake area changes between two successive time steps.
    
    change.tib <- matrix(NA, 
                         nrow = length(uni.years)-1,
                         ncol = 11) %>%
      as_tibble(.name_repair = "unique")
    
    colnames(change.tib) <- c("central_year_of_change", 
                              "lowest_year",
                              "year_low_central",
                              "year_high_central",
                              "highest_year", 
                              "perc_ch_median",
                              "perc_ch_min_05_px",
                              "perc_ch_min_1_px",
                              "perc_ch_max_05_px",
                              "perc_ch_max_1_px", 
                              "study")
    
    # From a given inventory, we always look at the next available and 
    # calculate the difference in lake area, divided by the length of the time
    # window. This yields the percent (or relative) change in lake area.
    
    for (j in 1: (length(uni.years)-1)) {
      
      year.one <- study %>% filter(Year == uni.years[j])
      year.two <- study %>% filter(Year == uni.years[j+1])
      
      change.tib$central_year_of_change[j] <- (uni.years[j] + uni.years[j+1]) / 2
      change.tib$lowest_year[j]            <- unique(year.one$Year_Min)
      change.tib$year_low_central[j]       <- unique(year.one$Year)
      change.tib$year_high_central[j]      <- unique(year.two$Year)
      change.tib$highest_year[j]           <- unique(year.two$Year_Max)
      
      ch.05 <- (((year.one$area_05_px - year.two$area_05_px) / 
                   year.one$area_05_px) * -1) / (uni.years[j+1] - uni.years[j]) 
      
      ch.1 <-  (((year.one$area_1_px - year.two$area_1_px) /
                   year.one$area_1_px) * -1) / (uni.years[j+1] - uni.years[j]) 
      
      change.tib$perc_ch_median[j] <- c(ch.05, ch.1) %>% 
        median() * 100 %>% 
        round(digits = 1)
      
      pc.05 <- quantile(ch.05, c(0.0275, 0.975)) * 100  %>% 
        round(digits = 1)
      pc.1 <- quantile(ch.1, c(0.0275, 0.975)) * 100  %>% 
        round(digits = 1)
      
      change.tib$perc_ch_min_05_px[j] <- pc.05[1]
      change.tib$perc_ch_max_05_px[j] <- pc.05[2]
      
      change.tib$perc_ch_min_1_px[j] <- pc.1[1]
      change.tib$perc_ch_max_1_px[j] <- pc.1[2]
      
      change.tib$study <- study$Study[j]
      
    }
    
    Study.list[[i]] <- change.tib
    
  }
  
  o <- do.call(rbind, Study.list) %>%
    mutate(region = reg,
           study = basename(study))
  
  return(o)
  
})

# We close the cluster environment and save this object to disk.

stopCluster(cl)
gc()

saveRDS(lakes.per.region, "lakes_per_region.rds")

# We remove inventories that had been compiled before our study period. 
# We also remove rates of change that are unduly high. This can be a consequence
# of only a handful of lakes mapped in a given inventory.

growth <- do.call(rbind, lakes.per.region) %>%
  filter(year_low_central >= 1985,
         perc_ch_min_1_px > -50,
         perc_ch_max_1_px < 25) %>%
  mutate(region = str_replace_all(region, pattern = "./", replacement = "")) 

# We count the number of inventories in each region. This information
# will be displayed in a figure later.

nobs <- growth %>%
  group_by(region) %>%
  summarise(n = n(),
            nobs = paste0("no = ", n())) %>%
  add_row(region = "Global", 
          n = nrow(growth),
          nobs = paste0("no = ", nrow(growth)))  %>%
  mutate(region = str_replace_all(region, "_", " "),
         region_f = factor(region, 
                           levels = c("Global",
                                      "Alaska Range",
                                      "W Chugach Mountains",
                                      "Saint Elias Mountains",
                                      "Coast Ranges",
                                      "C and N Andes",
                                      "Patagonia",
                                      "Iceland",
                                      "Scandinavia",
                                      "Alps",
                                      "Hissar Alay and Tien Shan",
                                      "HK Pamir Karakoram",
                                      "Himalayas",
                                      "Tibet and Hengduan Shan")))

# We do the same statistic for the number of studies used in every region.

nstud <- growth %>%
  group_by(region) %>%
  summarise(nstud = paste0("ns = ",length(unique(study)))) %>%
  add_row(region  = "Global",
          nstud   = paste0("ns = ", length(unique(growth$study))))  %>%
  mutate(region = str_replace_all(region, "_", " "),
         region_f = factor(region, 
                           levels = c("Global",
                                      "Alaska Range",
                                      "W Chugach Mountains",
                                      "Saint Elias Mountains",
                                      "Coast Ranges",
                                      "C and N Andes",
                                      "Patagonia",
                                      "Iceland",
                                      "Scandinavia",
                                      "Alps",
                                      "Hissar Alay and Tien Shan",
                                      "HK Pamir Karakoram",
                                      "Himalayas",
                                      "Tibet and Hengduan Shan")))

# We now would like to know the average rate of change per region. We do this by
# using a hierarchical linear regression model. One model has variable intercepts, 
# but fixed slopes. That is, every region has a different base rate of change, but the 
# rate of change is not allowed to vary with time. The other model has both 
# varying intercepts and varying slopes, and is therefore to indicate a deceleration
# or acceleration of lake growth.

# We scale the data to mean zero and unit standard deviation.

dat <- growth %>%
  transmute(central_year_of_change, perc_ch_median, region) %>%
  mutate(x_scale = scale(central_year_of_change)[ ,1],
         y_scale = scale(perc_ch_median)[ ,1])

# Obtain the original mean and standard deviation of predictor and response.

sd_y   <- sd(dat$perc_ch_median)
mean_y <- mean(dat$perc_ch_median) 
sd_x   <- sd(dat$central_year_of_change)
mean_x <- mean(dat$central_year_of_change)

# We specify weakly informed, Student t-distributed priors for the model.

bprior <- prior(student_t(3, 0, 2), class = "b") +
  prior(student_t(3, 0, 2), class = "Intercept") +
  prior(student_t(3, 0.5, 2), class = "sd") 

# Varying intercepts, varying slopes: We model the trend in percent lake area 
# change for all regions in a given time interval.

mod <- brm(bf(y_scale ~ x_scale + (x_scale | region)),
           family = student(),
           data = dat,
           prior = bprior,
           iter    = 4000, 
           warmup  = 1000,
           chains  = 4,
           cores   = 4, 
           control = list(adapt_delta = 0.95,
                          max_treedepth = 15),
           backend = "cmdstanr",
           threads = threading(3))

# Assess model fit

mod
plot(mod)
pp_check(mod)

# Looks good!

# Now with varying intercepts, but fixed slopes: We model only the average rate 
# in lake area change.

bprior <- prior(student_t(3, 0, 2), class = "Intercept") +
  prior(student_t(3, 0.5, 2), class = "sd") 

mod.intercept <- brm(bf(y_scale ~ 1 + (1 | region)),
                     family = student(),
                     data = dat,
                     # prior = bprior,
                     iter    = 4000, 
                     warmup  = 1000,
                     chains  = 4,
                     cores   = 4, 
                     control = list(adapt_delta = 0.95,
                                    max_treedepth = 15),
                     backend = "cmdstanr",
                     threads = threading(3))

# Assess model fit

mod.intercept
plot(mod.intercept)
pp_check(mod.intercept)

# Looks fine!

# Define the time period, for which new predictions will be made.
# ... on original scale

seq.x.orig <- seq(from = 1985,
                  to   = 2023,
                  length.out = 132)

# ... and on scaled range

seq.x.scaled <- (seq.x.orig - mean_x) / sd_x

# Obtain the standardized posterior predictive distribution for new observations
# and convert the predictions to original scale.
# Create an a tibble that we will use to specify the predictions.

pred.tibble <- tibble(x_scale = rep(seq.x.scaled, times = length(unique(dat$region))),
                      central_year_of_change  = rep(seq.x.orig, times = length(unique(dat$region))),
                      region  = rep(unique(dat$region), each  = length(seq.x.scaled)))

# Make 500 predictions for each observation.
# Convert predictions to original scale.

# On regional level...

post_epred <- add_epred_draws(mod,
                              ndraws = 500,
                              re_formula = NULL,
                              value = "y_scale",
                              newdata    = pred.tibble) %>%
  ungroup() %>%
  mutate(central_year_of_change = (x_scale * sd_x) + mean_x,
         perc_ch_median = (y_scale * sd_y) + mean_y)

# And on global level.

post_epred_global <- add_epred_draws(mod,
                                     ndraws  = 500,
                                     re_formula = y_scale ~ x_scale,
                                     value   = "y_scale",
                                     newdata = tibble(x_scale = seq.x.scaled)) %>%
  ungroup() %>%
  mutate(central_year_of_change = (x_scale * sd_x) + mean_x,
         perc_ch_median = (y_scale * sd_y) + mean_y,
         region = "Global")

# We combine regional and global predictions.

post_epred_all <- rbind(post_epred,
                        post_epred_global)

# We extract posterior trends in lake growth for each region. 
# We retransform the data to the original scale

slopes.cond <- tidy_draws(mod) %>% 
  pivot_longer(starts_with(paste0("r_", "region")),
               names_to = "param",
               values_to = "post") %>% 
  filter(str_detect(param, 'x_scale')) %>%
  transmute(b_x_scale, param, post) %>%
  mutate(param = str_replace_all(param, "[.]", " "),
         region = str_extract_all(param,"(?<=\\[).+(?=,)") %>% unlist(),
         post_orig = b_x_scale + post,
         post_orig = post_orig  * sd_y / sd_x) 

# We extract posterior trends for the population-level (global) model. 

slopes.all <- mod %>%
  gather_draws(b_x_scale) %>% 
  mutate(region = "Global", param = "x_scale") %>% 
  mutate(post_orig = .value  * sd_y / sd_x) %>%
  ungroup() %>%
  dplyr::select(!starts_with("."))

# We combine slopes from pooled estimate and from group levels.

slopes <- bind_rows(slopes.cond ,
                    slopes.all) %>%
  mutate(first_year = 1985,
         last_year  = 2023,
         response   = "perc_ch_median") 

# Plot the regional trends in lake area change.

mod_param.slope <- slopes %>%
  ggplot(aes(x = post_orig,
             y = region)) +
  stat_pointinterval( .width = c(0.95), 
                      point_size = 1,
                      size = 0.8) + 
  theme_bw() +
  geom_vline(xintercept = 0) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "bottom") +
  guides(color = guide_legend(reverse = TRUE)) 

# We extract the intercepts from the model without slopes, i.e. the average 
# rates of lake growth

# ... for each region

intercepts.cond <- tidy_draws(mod.intercept) %>% 
  pivot_longer(starts_with(paste0("r_", "region")),
               names_to = "param",
               values_to = "post") %>% 
  filter(str_detect(param, 'Intercept')) %>%
  transmute(b_Intercept, param, post) %>%
  mutate(param = str_replace_all(param, "[.]", " "),
         region = str_extract_all(param,"(?<=\\[).+(?=,)") %>% unlist(),
         post_orig = b_Intercept + post,
         post_orig = (post_orig  * sd_y) + mean_y) 

# ... and globally.

intercepts.global <- tidy_draws(mod.intercept) %>% 
  transmute(b_Intercept) %>%
  rename(post = b_Intercept) %>%
  mutate(b_Intercept = NA,
         param = "r_region[Global,Intercept]",
         region = "Global",
         post_orig = (post  * sd_y) + mean_y) 

# We combine the regional and global lake growth rates to one tibble.

intercepts.all <- rbind(intercepts.cond,
                        intercepts.global )

# We plot the posterior distribution of the annual lake growth.

annual.lake.growth <- intercepts.all  %>%
  mutate(region = str_replace_all(region, "_", " ")) %>%
  left_join(x = ., y = nobs, by = "region") %>%
  mutate(reg_obs = paste0(region, " (", n, ")")) %>%
  mutate(region2 = factor(reg_obs, 
                          levels = c(unique(grep("Global", reg_obs, value = T)),
                                     unique(grep("Alaska Range", reg_obs, value = T)),
                                     unique(grep("W Chugach Mountains", reg_obs, value = T)),
                                     unique(grep("Saint Elias Mountains", reg_obs, value = T)),
                                     unique(grep("Coast Ranges", reg_obs, value = T)),
                                     unique(grep("C and N Andes", reg_obs, value = T)),
                                     unique(grep("Patagonia", reg_obs, value = T)),
                                     unique(grep("Iceland", reg_obs, value = T)),
                                     unique(grep("Scandinavia", reg_obs, value = T)),
                                     unique(grep("Alps", reg_obs, value = T)),
                                     unique(grep("Hissar Alay and Tien Shan", reg_obs, value = T)),
                                     unique(grep("HK Pamir Karakoram", reg_obs, value = T)),
                                     unique(grep("Himalayas", reg_obs, value = T)),
                                     unique(grep("Tibet and Hengduan Shan", reg_obs, value = T))))) %>%
  mutate(region2 = fct_rev(region2)) %>%
  ggplot(aes(x = post_orig,
             y = region2)) +
  stat_halfeye(.width = 0.95,
               slab_size = 1,
               interval_size = 2,
               interval_color = "black",
               slab_fill = "lightblue") + 
  theme_bw() +
  xlim(c(-0.5, 3)) +
  labs(y = "Region",
       x = "Annual change in total lake area [%]\nsince the late 1980s") +
  theme( axis.text = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title = element_text(size = 7),
         strip.text = element_text(size = 7),
         legend.position = "none") 

# We extract statistics on the posterior rate (median and 95% HDI) in lake growth. 
# This information will be displayed in a figure later.

change.stats <- intercepts.all %>% group_by(region) %>% 
  summarise(median = median(post_orig), 
            q025 = quantile(post_orig, 0.025), 
            q975 = quantile(post_orig, 0.975),
            low = median - q025,
            high = q975 - median,
            text_median = paste0(round(median, digits = 1), "%"),
            text_hl =     paste0("(",round(q025, digits = 1), "-",
                                 round(q975, digits = 1), "%)" )) %>%
  mutate(region = str_replace_all(region, "_", " "),
         region_f = factor(region, 
                           levels = c("Global",
                                      "Alaska Range",
                                      "W Chugach Mountains",
                                      "Saint Elias Mountains",
                                      "Coast Ranges",
                                      "C and N Andes",
                                      "Patagonia",
                                      "Iceland",
                                      "Scandinavia",
                                      "Alps",
                                      "Hissar Alay and Tien Shan",
                                      "HK Pamir Karakoram",
                                      "Himalayas",
                                      "Tibet and Hengduan Shan")))

# We assess, if the global and regional rates in lake area change are credibly 
# larger than zero.

intercepts.all %>%
  ggplot(aes(x = post_orig,
             y = region)) +
  stat_pointinterval( .width = c(0.95), 
                      point_size = 1,
                      size = 0.8) + 
  theme_bw() +
  geom_vline(xintercept = 0) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "bottom") +
  guides(color = guide_legend(reverse=TRUE)) 

# We now would like to plot the trends in lake growth and the average rate of 
# lake area change.

growth.global <- growth %>% mutate(region = "Global") 

# Remove underscores from the region names in the tibble.

growth.all <- rbind(growth, growth.global) %>%
  mutate(region = str_replace_all(region, "_", " "),
         region_f = factor(region, 
                           levels = c("Global",
                                      "Alaska Range",
                                      "W Chugach Mountains",
                                      "Saint Elias Mountains",
                                      "Coast Ranges",
                                      "C and N Andes",
                                      "Patagonia",
                                      "Iceland",
                                      "Scandinavia",
                                      "Alps",
                                      "Hissar Alay and Tien Shan",
                                      "HK Pamir Karakoram",
                                      "Himalayas",
                                      "Tibet and Hengduan Shan")))

post_epred_all <- post_epred_all %>%
  mutate(region = str_replace_all(region, "_", " "),
         region_f = factor(region, 
                           levels = c("Global",
                                      "Alaska Range",
                                      "W Chugach Mountains",
                                      "Saint Elias Mountains",
                                      "Coast Ranges",
                                      "C and N Andes",
                                      "Patagonia",
                                      "Iceland",
                                      "Scandinavia",
                                      "Alps",
                                      "Hissar Alay and Tien Shan",
                                      "HK Pamir Karakoram",
                                      "Himalayas",
                                      "Tibet and Hengduan Shan")))

# The trend (acceleration/ decelartion) is in green. Individual
# rates of change, including their uncertainties, are shows as black cross hairs.

lake.growth.plot <- ggplot(data = growth.all, 
                           mapping = aes(x = central_year_of_change, 
                                         y = perc_ch_median)) +
  stat_lineribbon(data = post_epred_all,
                  aes(x = central_year_of_change,
                      y = perc_ch_median),
                  fill = "darkolivegreen3",
                  linewidth = 0.5,
                  color = "darkolivegreen4",
                  .width = c(0.95),
                  alpha = 0.6) +
  geom_linerange(aes(y = perc_ch_median,
                     xmin  = lowest_year,
                     xmax  = highest_year),
                 linewidth = 0.3) +
  geom_linerange(aes(x = central_year_of_change,
                     ymin  = perc_ch_min_1_px,
                     ymax  = perc_ch_max_1_px),
                 linewidth  = 0.3) +
  geom_linerange(data = change.stats,
                 aes(x = 2025,
                     ymin  = q025,
                     ymax  = q975,
                     color = "darkolivegreen4"),
                 inherit.aes = F,
                 linewidth  = 1) +
  geom_point(data = change.stats,
             aes(x = 2025,
                 y = median, 
                 color = "darkolivegreen4"), 
             inherit.aes = F,
             shape = 16, 
             size = 0.9) +
  geom_point(shape = 16, 
             size = 0.9) +
  scale_color_identity() +
  geom_hline(yintercept = 0,
             color = "grey50") +
  facet_wrap(~region_f) +
  xlim(c(1985, 2025)) +
  ylim(c(-25, 25)) + 
  theme(legend.position = "none") +
  xlab("Period covered in lake inventory") +
  ylab("Annual change in total lake area [%]")  +
  geom_text(data  = nstud,
            aes(x = 1985, y = -17, label = nstud),
            size = 2.3,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left") +
  geom_text(data  = nobs,
            aes(x = 1985, y = -23, label = nobs),
            size = 2.3,
            colour = "gray10",
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "left") +
  geom_text(aes(x = 2025,
                y = -17,
                label = text_median),
            data = change.stats,
            color = "darkolivegreen4",
            size = 2.3,
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "right") +
  geom_text(aes(x = 2025,
                y = -23,
                label = text_hl),
            data = change.stats,
            color = "darkolivegreen4",
            size = 2.3,
            inherit.aes = FALSE,
            parse = FALSE,
            hjust = "right") +
  theme_bw() +
  theme(axis.text = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.text = element_text(size = 7),
        strip.background = element_blank())

# We Write this plot to disk (Extended Data Figure .

ggsave(
  filename = "lake_area_change.pdf",
  plot = lake.growth.plot ,
  width = 180,
  height = 140,
  units = "mm"
)