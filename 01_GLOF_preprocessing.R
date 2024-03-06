################################################################################
#######         Preprocess tabular data on GLOF occurrences and           ######
#######                 extract descriptive statistics                    ######
#######                                                                   ######
#######                            by Georg Veh                           ######
#######                           March 28, 2023                          ######
#######                       revised Oct 26, 2023                        ######
#######                checked and comments added March 04, 2024          ######
################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

# Important: we run brms models using the cmdstan backend. This needs to be
# installed separately using the instructions given here: 
# https://mc-stan.org/cmdstanr/ 

require(tidyverse)
require(tidybayes)
require(modelr)
require(scales)
require(readODS)
require(sf)
require(lubridate)
require(brms)
require(gridExtra)
require(ggpubr)
require(ggExtra)
require(see)

# Set YOUR working directory folder where to find all files, necessary to run 
# this script. Change the location appropriately.

setwd("D:/data/BoxUP/Work 2022/GLOFsize")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

# Tell R the name of the open-office spreadsheet that contains reported GLOFs per 
# region in separate sheets.

glof.file <- "glofdatabase_2023_12_06.ods"

# Get names of the sheets in the Open Office document. 
# Exclude 'Global', 'Other', and 'Greenland'.

sheetnames <- list_ods_sheets(glof.file)

sheetnames <- sheetnames[!(sheetnames %in% "Other")]
sheetnames <- sheetnames[!(sheetnames %in% "Greenland")]

# We create an empty list which will later have one entry per mountain region.
# Every list element will be a table that holds information about reported GLOFs.

region.list <- list()

# Load the table of regional GLOF reports into memory. 
# We iterate over the names of the spreadsheet.

for(r in sheetnames) {
  
  data <- read_ods(glof.file, sheet = r) %>%
    as_tibble(.name_repair = "unique")
  
  data <- data[-c(1:2), 1:57] 
  data$region <- r 
  
  y <- as.numeric(str_sub(data$Date, 1, 4))
  
  # Extract the years of reported GLOF occurrences. 
  
  for (i in 1 : length(y)) {
    
    # Some years are NA because these GLOFs have no fixed date of occurrence, 
    # but a range of possible dates. 
    # If there is NA in the 'Date' column, we first check whether there is a  
    # given range of dates. If so, we then randomly sample for the range of 
    # plausible years. 
    
    if (is.na(y[i])) {
      
      min.date <- as.numeric(str_sub(data$Date_Min[i], 1, 4))
      max.date <- as.numeric(str_sub(data$Date_Max[i], 1, 4))
      
      if((!is.na(min.date)) & (!is.na(max.date))) {
        
        obs.period <- min.date:max.date 
        
        if (length(obs.period) == 1 ) {
          
          random.year <- obs.period
          
        } else {   random.year <- sample(obs.period, size = 1)}
        
        y[i] <- random.year
        
      } 
    } 
  }
  
  # Append the year column to the data.frame. We also make sure that
  # every column has the format we need to process the data.
  
  data <- data %>% 
    mutate(rounded_year = y) %>%
    mutate(
      Longitude = as.numeric(Longitude),
      Latitude  = as.numeric(Latitude),
      Mean_Lake_Volume_VL = as.numeric(Mean_Lake_Volume_VL),
      Min_VL    = as.numeric(Min_VL), 
      Max_VL    = as.numeric(Max_VL), 
      Mean_Flood_Volume_V0 = as.numeric(Mean_Flood_Volume_V0), 
      Min_V0    = as.numeric(Min_V0), 
      Max_V0    = as.numeric(Max_V0),
      Peak_discharge_Qp = as.numeric(Peak_discharge_Qp), 
      Min_Qp    = as.numeric(Min_Qp), 
      Max_Qp    = as.numeric(Max_Qp),
      Reference = as.character(Reference),
      D_buildings = as.character(D_buildings),
      reported_impacts = as.character(reported_impacts),
      economic_losses	 = as.character(economic_losses),
      D_buildings	     = as.character(D_buildings),
      D_bridges	       = as.character(D_bridges),
      D_roads_paths	   = as.character(D_roads_paths),
      D_railroads	  	 = as.character(D_railroads),
      D_utilities	     = as.character(D_utilities),
      D_flood_protection = as.character(D_flood_protection),
      D_environmental	 = as.character(D_environmental),
      resettlement	   = as.character(resettlement),
      reported_fatalities = as.character(reported_fatalities),
      Image_date_after    = as.character(Image_date_after),
      Image_date_before   = as.character(Image_date_before),
      Lake_area_before    = as.numeric(Lake_area_before), 
      Certainty_level_before = as.numeric(Certainty_level_before),
      Lake_area_after = as.numeric(Lake_area_after), 
      Certainty_level_after = as.numeric(Certainty_level_after)
    )
  
  region.list[[r]] <- data
  
}

# We concatenate all tibbles of regional GLOF occurrences to one long tibble.

all.glofs <- bind_rows(region.list) 

# We may write this table to disk to simply read the preprocessed table later.

saveRDS(all.glofs, "all_glofs_tibble.RDS")
# all.glofs <- readRDS("all_glofs_tibble.RDS")

all.glofs <- all.glofs %>%
  mutate(la_diff = Lake_area_before - Lake_area_after)

################################################################################

# Our focus is on temporal trends in lake areas before the outburst. Put 
# differently, we would like to know whether lakes that had burst out have
# become larger, as lakes have grown in past decades.

# The table needs some preprocessing.
# We trim the table to the study period 1990 - 2023. We remove all lake types 
# that are not moraine, ice, or bedrock. We calculate the difference in 
# lake area in square kilometers.
# Finally, we allocate space for three empty columns that specify the region and
# the dam type.

la.sf <- all.glofs %>%
  mutate(Glacier = replace_na(Glacier, "unknown"),
         Lake    = replace_na(Lake,    "unknown")) %>%
  filter(rounded_year >= 1990,
         rounded_year <= 2023,
         Glacier != "Gl. Perito Moreno",
         !is.na(Latitude),
         Lake_type != "water pocket",
         Lake_type != "unknown",
         Lake_type != "other – colluvial material" ,
         Lake_type != "subglacial",
         Lake_type != "snow",
         Lake_type != "bedrock/ice",
         !str_detect(Lake_type, "volc")) %>%
  mutate(la_diff = (Lake_area_before - Lake_area_after)/ 10^6) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  mutate(RegO1 = NA,
         RegO2 = NA,
         Lake_type_simple = NA)

# Lake_type simple distinguishes between the two main dam types in our study 
# regions: moraine and bedrock versus ice-dammed and supraglacial lakes.

la.sf$Lake_type_simple[str_detect(la.sf$Lake_type, pattern = "moraine|bedrock")] <- "moraine_bedrock"
la.sf$Lake_type_simple[is.na(la.sf$Lake_type_simple)] <- "glacier_supraglacial"

# We add the region codes from the Randolph Glacier Inventory, V6.0 to all 
# reported GLOFs. Navigate to the GLIMS repository (https://www.glims.org/RGI/rgi60_dl.html) 
# and download the complete RGI inventory, which includes both glacier and
# regional outlines. 

rgi.O1 <- st_read("rgi06/00_rgi60_regions/00_rgi60_O1Regions.shp") %>%
  st_make_valid()

# The coordinate reference system of the RGI is WGS84, which has the EPSG code 4326.

st_crs(rgi.O1) <- 4326

# We add the large scale region (O1) according to RGI V6.0 to each GLOF.
# We intersect each reported GLOF (a point feature) with the RGI regions (a 
# polygon feature) and add the 'FULL_NAME' of the RGI region to the reported
# GLOF.

rgi.O1 <- st_make_valid(rgi.O1)

sp.mat.O1 <- st_intersects(rgi.O1, la.sf, sparse = F)

la.sf <- lapply(1:nrow(sp.mat.O1), FUN = function (x) {
  
  out <- la.sf[sp.mat.O1[x, ], ]
  
  out$RegO1 <- rgi.O1$FULL_NAME[x]
  
  if(nrow(out) == 0) {return(NULL)} else { return(out) }
  
})

la.sf <- la.sf[!sapply(la.sf, is.null)]
la.sf <- do.call(rbind, la.sf)

# We do exactly the same for the fine scale region (O2) according to RGI V6.0.

rgi.O2 <- st_read("rgi06/00_rgi60_regions/00_rgi60_O2Regions.shp") %>%
  st_make_valid()

st_crs(rgi.O2) <- 4326

rgi.O2 <- st_make_valid(rgi.O2)

sp.mat.O2 <- st_intersects( rgi.O2, la.sf, sparse = F)

la.sf <- lapply(1:nrow(sp.mat.O2), FUN = function (x) {
  
  out <- la.sf[sp.mat.O2[x, ], ]
  
  out$RegO2 <- rgi.O2$FULL_NAME[x]
  
  if(nrow(out) == 0) {return(NULL)} else { return(out)}
  
})

la.sf <- la.sf[!sapply(la.sf, is.null)]
la.sf <- do.call(rbind, la.sf)

# Some regions have only few lakes and reported outbursts. We group them
# to 13 larger regions.

la.sf <- la.sf %>%
  filter(RegO2 != "Svalbard") %>%
  mutate(Continent = NA,
         pch = NA) %>%
  mutate(RegO2_adj = RegO2) %>%
  mutate(RegO2_adj = str_replace_all(RegO2_adj, "\\(|\\)", "")) %>%
  mutate(
    RegO2_adj = str_replace_all(RegO2_adj, "C Andes",      "C_and_N_Andes"),
    RegO2_adj = str_replace_all(RegO2_adj, "Low-latitude Andes", "C_and_N_Andes"),
    RegO2_adj = str_replace_all(RegO2_adj, "C Himalaya",   "Himalayas"),
    RegO2_adj = str_replace_all(RegO2_adj, "E Himalaya",   "Himalayas"),
    RegO2_adj = str_replace_all(RegO2_adj, "W Himalaya",   "Himalayas"),
    RegO2_adj = str_replace_all(RegO2_adj, "C Himalaya",   "Himalayas"),
    RegO2_adj = str_replace_all(RegO2_adj, "Karakoram",    "HK_Pamir_Karakoram"),
    RegO2_adj = str_replace_all(RegO2_adj, "Pamir Safed Khirs/W Tarim", "HK_Pamir_Karakoram"),
    RegO2_adj = str_replace_all(RegO2_adj, "Hindu Kush",   "HK_Pamir_Karakoram"),
    RegO2_adj = str_replace_all(RegO2_adj, "Hengduan Shan", "Tibet_and_Hengduan_Shan"),
    RegO2_adj = str_replace_all(RegO2_adj, "Inner Tibet",   "Tibet_and_Hengduan_Shan"),
    RegO2_adj = str_replace_all(RegO2_adj, "S and E Tibet", "Tibet_and_Hengduan_Shan"),
    RegO2_adj = str_replace_all(RegO2_adj, "SE Scandinavia", "Scandinavia"),
    RegO2_adj = str_replace_all(RegO2_adj, "SW Scandinavia", "Scandinavia"),
    RegO2_adj = str_replace_all(RegO2_adj, "N Scandinavia", "Scandinavia"),
    RegO2_adj = str_replace_all(RegO2_adj, "Hissar Alay",   "Hissar_Alay_and_Tien_Shan"),
    RegO2_adj = str_replace_all(RegO2_adj, "W Tien Shan",   "Hissar_Alay_and_Tien_Shan"),
    RegO2_adj = str_replace_all(RegO2_adj, "E Tien Shan Dzhungaria", "Hissar_Alay_and_Tien_Shan"),
    RegO2_adj = str_replace_all(RegO2_adj, "N Coast Ranges", "Coast_Ranges"),
    RegO2_adj = str_replace_all(RegO2_adj, "S Coast Ranges", "Coast_Ranges"),
    RegO2_adj = str_replace_all(RegO2_adj, "W Chugach Mtns Talkeetna", "W_Chugach_Mountains"),
    RegO2_adj = str_replace_all(RegO2_adj, "Alaska Ra Wrangell/Kilbuck", "Alaska_Range"),
    RegO2_adj = str_replace_all(RegO2_adj, "St Elias Mtns", "Saint_Elias_Mountains")) %>%
  mutate(RegO2_adj_f = factor(str_replace_all(RegO2_adj, "_", " "),
                              levels = c("Alaska Range", "W Chugach Mountains", "Saint Elias Mountains", 
                                         "Coast Ranges", "Iceland", "Scandinavia", "Alps", 
                                         "Hissar Alay and Tien Shan", "HK Pamir Karakoram", 
                                         "Tibet and Hengduan Shan", "Himalayas", "C and N Andes", 
                                         "Patagonia")))

# We also add the continent to the mountain regions.

la.sf$Continent[grep(pattern = "Alaska_Range|W_Chugach_Mountains|Saint_Elias_Mountains|Coast_Ranges", 
                     x = la.sf$RegO2_adj)] <- "NW_North_America"
la.sf$pch[grep(pattern = "NW_North_America", x = la.sf$Continent)] <- 15

la.sf$Continent[grep(pattern = "Iceland|Scandinavia|Alps", 
                     x = la.sf$RegO2_adj)] <- "Europe"
la.sf$pch[grep(pattern = "Europe", x = la.sf$Continent)] <- 16

la.sf$Continent[grep(pattern = "Hissar_Alay_and_Tien_Shan|HK_Pamir_Karakoram|Tibet_and_Hengduan_Shan|Himalayas", 
                     x = la.sf$RegO2_adj)] <- "High_Mountain_Asia"
la.sf$pch[grep(pattern = "High_Mountain_Asia", x = la.sf$Continent)] <- 17

la.sf$Continent[grep(pattern = "C_and_N_Andes|Patagonia", 
                     x = la.sf$RegO2_adj)] <- "South_America"
la.sf$pch[grep(pattern = "South_America", x = la.sf$Continent)] <- 18

# Write the preprocessed GLOF data table to disk.

saveRDS(la.sf, "la_sf.RDS")

# Add 'global' as another region to the tibble

la.sf.global <- la.sf %>%
  mutate(Continent = "Global", 
         RegO2_adj_f = factor("Global"))

# We attach the global tibble to the regional one. 

la.sf.all <- rbind(la.sf, la.sf.global)

# We remove all cases that have no lake area before the GLOF. 

la.sf.global.regional.nona <- la.sf.all %>%
  filter(!is.na(Lake_area_before)) %>%
  st_drop_geometry() %>%
  mutate(Lake_type_simple = replace(Lake_type_simple, 
                                    Lake_type_simple == "glacier_supraglacial", 
                                    "Glacier & supraglacial"),
         Lake_type_simple = replace(Lake_type_simple, 
                                    Lake_type_simple == "moraine_bedrock", 
                                    "Moraine & bedrock")) 

# We write a geopackage to disk that contains the largest reported GLOF per site.
# This file is part of figure 1 (i.e. the open circles).

# la.sf.all %>%
#   filter(!is.na(Lake_area_before)) %>%
#   filter(RegO2_adj_f != "Global") %>%
#   mutate(Lake_type_simple = replace(Lake_type_simple,
#                                     Lake_type_simple == "glacier_supraglacial",
#                                     "Glacier & supraglacial"),
#          Lake_type_simple = replace(Lake_type_simple,
#                                     Lake_type_simple == "moraine_bedrock",
#                                     "Moraine & bedrock")) %>%
#   group_by(RGI_Glacier_Id, Lake, Lake_type_simple) %>%
#   summarise(n = n(),
#             max_a = max(Lake_area_before)) %>%
#   st_write(dsn = "Max_area_and_number_per_lake.gpkg", delete_dsn = T)

# We load a shapefile that contains the outlines of the 13 (partly dissolved) 
# RGI regions. 

dissolved.rgi <- st_read("rgi06/rgiO2_dissolved_outlines.shp") %>%
  rename(RegO2_adj_f = FULL_NAME) %>%
  mutate(RegO2_adj_f = str_replace_all(RegO2_adj_f, "_", " "))

# We now count the number of GLOFs per region.
# First, across all lake types... 

counts.per.region.all <- la.sf.all %>%
  filter(!is.na(Lake_area_before)) %>%
  filter(RegO2_adj_f != "Global") %>%
  group_by(RegO2_adj_f) %>%
  summarise(n_global = n())

# then, only for glacier-dammed and supraglacial lakes...

counts.per.region.gs <- la.sf.global.regional.nona %>% 
  filter(RegO2_adj_f != "Global",
         Lake_type_simple == "Glacier & supraglacial") %>%
  group_by(RegO2_adj_f) %>%
  summarise(n_glac_sup = n())

# and finally, only for moraine- and bedrock-dammed lakes.

counts.per.region.mb <- la.sf.global.regional.nona %>% 
  filter(RegO2_adj_f != "Global",
         Lake_type_simple == "Moraine & bedrock") %>%
  group_by(RegO2_adj_f) %>%
  summarise(n_mor_bed = n())

# We extract the centroid of each region (i.e. generate a point feature) 
# and add these regional counts as an attribute to each centroid.

counts.global <- la.sf.all %>%
  filter(!is.na(Lake_area_before)) %>%
  filter(RegO2_adj_f == "Global") %>%
  summarise(n_global = n()) %>%
  st_centroid() %>%
  mutate(RegO2_adj_f = "Global",
         Region = "Global",
         n_glac_sup = sum(counts.per.region.gs$n_glac_sup),
         n_mor_bed = sum(counts.per.region.mb$n_mor_bed)) 

# We write a geopackage to disk that contains the regional number of GLOF 
# occurrences. This file is part of figure 1 (i.e. the pie charts).

# dissolved.rgi %>%
#   left_join(y = counts.per.region.all %>% st_drop_geometry()) %>%
#   left_join(y = counts.per.region.gs %>% st_drop_geometry()) %>%
#   left_join(y = counts.per.region.mb %>% st_drop_geometry()) %>%
#   replace_na(list(n_global = 0, n_glac_sup = 0, n_mor_bed = 0)) %>%
#   st_centroid() %>%
#   rbind(counts.global) %>%
#   st_write("Counts_per_region.gpkg",
#            delete_dsn = T)

################################################################################
### Descriptive statistics #####################################################

# We now generate a number of descriptive statistics that will be mentioned
# in the associate manuscript.

# Total number of observations, including lakes without mapped areas

nrow(la.sf)

# Number of lakes with lake area before the GLOF.

la.sf.nona.before <- la.sf %>%
  filter(!is.na(Lake_area_before)) %>%
  st_drop_geometry() 

nrow(la.sf.nona.before)

# Number of unverified GLOFs

nrow(la.sf) - nrow(la.sf.nona.before) 

# Percentage of verified lakes

nrow(la.sf.nona.before) * 100 /nrow(la.sf) 

# Percentage of lakes with area after

la.sf.nona.after <- la.sf.nona.before %>%
  filter(!is.na(Lake_area_after))  

nrow(la.sf.nona.after) * 100 / nrow(la.sf.nona.before)

# Number of observations per region and lake type

nobs <- la.sf.global.regional.nona %>%
  group_by(RegO2_adj_f, Lake_type_simple) %>%
  summarise(n_orig = n()) %>%
  ungroup()

nobs %>% View()

# Number of occurrences per lake type

table(la.sf.nona.before$Lake_type)

# All occurrences per region 

nobs %>% 
  group_by(RegO2_adj_f) %>% 
  summarise(tot = sum(n_orig)) %>% 
  View()

# Share of moraine-dam failures in global GLOF occurrence

n.moraine <- nobs %>% filter(RegO2_adj_f == "Global", 
                             Lake_type_simple == "Moraine & bedrock")

n.moraine$n_orig / nrow(la.sf.nona.before)

# Largest lake area before GLOF

la.sf.nona.before %>%
  slice_max(Lake_area_before) %>% View()

# Largest difference in lake area

la.sf.nona.after %>%
  slice_max(la_diff) %>% View()

# Largest moraine-dam failure

la.sf.nona.before %>% 
  filter(Lake_type == "moraine") %>%
  slice_max(Lake_area_before) %>%
  View()

# Largest moraine-dam failure

la.sf.nona.before %>% 
  filter(Lake_type == "bedrock") %>%
  slice_max(Lake_area_before) %>%
  View()

# Largest difference in lake area for moraine-dammed lakes:

la.sf.nona.after %>%
  filter(Lake_type_simple == "moraine_bedrock") %>%
  slice_max(la_diff) %>% 
  View()

# Lake areas before: ice vs. moraine-dammed lakes:

la.sf.nona.before %>%
  group_by(Lake_type_simple) %>%
  summarise(med = median(Lake_area_before)/10^6,
            q25 = quantile(Lake_area_before, 0.25)/ 10^6,
            q75 = quantile(Lake_area_before, 0.75)/ 10^6) %>%
  mutate(qplus = q75 - med,
         qminus = med - q25)

# Partial vs. complete drainage

la.sf.nona.after %>%
  group_by(Lake_type_simple) %>%
  summarise(n = n(),
            n0 = sum(Lake_area_after == 0)) %>%
  mutate(ratio = (n0/n)*100 )


# Number of GLOFs first reported by us

all.refs <- table(la.sf.nona.before$Reference) 

# ... in NW North America

la.sf.nona.before %>%
  filter(Continent == "NW_North_America") %>%
  group_by(Reference) %>%
  summarise(nref = n()) %>%
  arrange(desc(nref)) %>%
  View()

# ... in Patagonia

la.sf.nona.before %>%
  filter(RegO2_adj == "Patagonia") %>%
  mutate(Glacier_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  group_by(Glacier_lake, Reference) %>%
  summarise(nref = n()) %>%
  arrange(desc(nref)) %>%
  View()

# Length of the GLOF detection window

observation.window <- la.sf.nona.after %>%
  filter(!is.na(Image_date_before),
         !is.na(Image_date_after)) %>%
  mutate(Image_date_before = ymd(Image_date_before),
         Image_date_after = ymd(Image_date_after)) %>%
  mutate(observation_window =  as.numeric(Image_date_after - Image_date_before) ) 

observation.window %>%   
  summarise(median = quantile(observation_window, 0.5),
            q25 =   quantile(observation_window, 0.25),
            q75 =   quantile(observation_window, 0.75))

# Observation window longer than half a year

length(observation.window$observation_window[observation.window$observation_window > 182.5])/
  length(observation.window$observation_window)

observation.window %>%
  filter(nchar(Date) == 10) %>%
  mutate(Date = ymd(Date)) %>%
  mutate(Period_before = Date - Image_date_before,
         Period_after  = Image_date_after - Date) %>%
  summarise(median_before = quantile(Period_before, 0.5),
            q025_before =   quantile(Period_before, 0.25),
            q975_before =   quantile(Period_before, 0.75),
            median_after =  quantile(Period_after, 0.5),
            q025_after =    quantile(Period_after, 0.25),
            q975_after =    quantile(Period_after, 0.75))

# Lake outbursts between 1990 and 2015 to compare with Carrivick and Tweed (2016)

la.sf.nona.before %>%
  filter(rounded_year <= 2015) %>%
  nrow()

# Number of satellite images from different sensors

table(c(la.sf.nona.before$Satellite_before, la.sf.nona.before$Satellite_after))

# Certainty during mapping

table(la.sf.nona.before$Certainty_level_before)
table(la.sf.nona.after$Certainty_level_after)

# Lakes with repeated outbursts

la.sf.global.regional.nona %>%
  filter(region != "Global") %>%
  mutate(Glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  group_by(Glacier_and_lake) %>%
  summarise(nobs = n()) %>%
  arrange(desc(nobs)) %>%
  View()

# Number of outbursts per year

glofs.per.lake <- la.sf.nona.before %>%
  filter(rounded_year >= 1990, 
         rounded_year <= 2023) %>%
  mutate(glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  group_by(Continent, RegO2_adj_f, Lake_type_simple, glacier_and_lake) %>%
  summarise(n_glofs = n()) %>%
  arrange(desc(n_glofs)) 

# Ice-dammed lakes with more than 2 reported outbursts since 1990.

glofs.per.lake %>% 
  filter(Lake_type_simple == "glacier_supraglacial",
         n_glofs >= 2)

# Ice-dammed lakes with more than 10 reported outbursts since 1990.

glofs.per.lake %>% 
  filter(Lake_type_simple == "glacier_supraglacial",
         n_glofs >= 10)

# Sources of ice-dam failures in Patagonia.

glofs.per.lake %>% 
  filter(Lake_type_simple == "glacier_supraglacial",
         n_glofs >= 1, 
         RegO2_adj_f == "Patagonia")

# Sources of moraine- and bedrock-dam failures in High Mountain Asia.

glofs.per.lake %>% 
  filter(Lake_type_simple == "moraine_bedrock",
         n_glofs >= 1, 
         Continent == "High_Mountain_Asia")

# Sources of ice-dam failures in NW North America.

glofs.per.lake %>% 
  filter(Continent == "NW_North_America",
         Lake_type_simple == "glacier_supraglacial",
         RegO2_adj_f != "Coast Ranges")

# Moraine- and bedrock-dammed lakes that burst at least twice.

glofs.per.lake %>% 
  filter(n_glofs > 1,
         Lake_type_simple == "moraine_bedrock")


################################################################################
### Compare GLOF area with flood volume (V0) and peak discharge (Qp) ###########

# From all reported GLOFs, we extract only those that have both a reported 
# flood volume (V0) and a mapped lake area before the GLOF.
# We log10-transform the data, and then standardise lake areas and V0 to
# zero mean and unit standard deviation.

V0 <- all.glofs %>%
  filter(!is.na(Lake_area_before),
         !is.na(Mean_Flood_Volume_V0),
         Mean_Flood_Volume_V0 > 10^-3,
         la_diff > 0) %>%
  mutate(la_log = log10(Lake_area_before),
         la_diff_log = log10(la_diff),
         V0_log = log10(Mean_Flood_Volume_V0),
         la_log_scale = scale_this(la_log),
         la_diff_log_scale = scale_this(la_diff_log),
         V0_log_scale = scale_this(V0_log))


# From all reported GLOFs, we extract only those that have both a reported 
# peak discharge (Qp) and a mapped lake area before the GLOF.
# We log10-transform the data, and then standardise lake areas and Qp to
# zero mean and unit standard deviation.

Qp <- all.glofs %>%
  filter(!is.na(Lake_area_before),
         !is.na(Peak_discharge_Qp),
         la_diff > 0) %>%
  mutate(la_log = log10(Lake_area_before),
         la_diff_log = log10(la_diff),
         Qp_log = log10(Peak_discharge_Qp),
         la_log_scale = scale_this(la_log),
         la_diff_log_scale = scale_this(la_diff_log),
         Qp_log_scale = scale_this(Qp_log))

# We fit robust regression models of lake area vs. flood volume or peak discharge.

# We set weakly informed priors on the intercept and slope

bprior <- prior(student_t(3, 0, 2.5), class = "b") +
  prior(student_t(3, 0, 2.5), class = "Intercept")

# We fit the V0 model.

V0.mod <- brm(V0_log_scale ~ la_log_scale,
              data = V0,
              prior = bprior, 
              warmup  = 1000,
              iter = 4000,
              chains  = 4,
              cores   = 4, 
              control = list(adapt_delta = 0.95,
                             max_treedepth = 15),
              backend = "cmdstanr",
              threads = threading(3))

# To obtain the posterior predictive distribution, we define a range of lake
# areas, still as standardised data.

scaled.area <- (seq_range(c(3, 9), n = 200) - mean(V0$la_log)) / sd(V0$la_log)

# For each datum, we draw 2000 values from the posterior predictive distribution.
# Then, we re-transform the data to the original value range.

preds.V0 <-  add_predicted_draws(
  object = V0.mod, 
  newdata = data.frame(la_log_scale = scaled.area),
  value = "V0_log_scale", 
  ndraws = 2000,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(Lake_area_before =     10^((la_log_scale * sd(V0$la_log)) + mean(V0$la_log)),
         Mean_Flood_Volume_V0 = 10^((V0_log_scale * sd(V0$V0_log)) + mean(V0$V0_log)))

# Finally, we plot the data and the posterior predictive distribution.

plot.V0.la <- preds.V0 %>%
  ggplot(aes(x = Lake_area_before, y = Mean_Flood_Volume_V0)) +
  scale_fill_manual(name = "Posterior rate", values = "cyan") +
  stat_lineribbon(.width = 0.95,
                  point_interval = mean_qi) +
  geom_point2(data = V0,
              aes(x = Lake_area_before, y = Mean_Flood_Volume_V0),
              alpha = 0.6) +
  theme_bw() +
  labs(y = "Reported GLOF volume [Mm³]",
       x = "Lake area before GLOF [m²]") +
  scale_y_continuous(limits = c(10^-4, 10^6.5),
                     trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

plot.V0.la <- plot.V0.la +
  annotate(geom = "text",
           x = 10^3,
           y = 10^6,
           label = paste0("n = ", nrow(V0)),
           hjust = 0,
           color = "grey20",
           size = 0.36 * 7) 

# We fit exactly the the same model to Qp versus lake area.

Qp.mod <- brm(Qp_log_scale ~ la_log_scale,
              data = Qp,
              prior = bprior, 
              warmup  = 1000,
              iter = 4000,
              chains  = 4,
              cores   = 4, 
              control = list(adapt_delta = 0.95,
                             max_treedepth = 15),
              backend = "cmdstanr",
              threads = threading(3))

# Posterior predictive distribution of the Qp model.

scaled.area <- (seq_range(c(3, 9), n = 200) - mean(Qp$la_log)) / sd(Qp$la_log)

preds.Qp <-  add_predicted_draws(
  object = Qp.mod, 
  newdata = data.frame(la_log_scale = scaled.area),
  value = "Qp_log_scale", 
  ndraws = 4000,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(Lake_area_before =     10^((la_log_scale * sd(Qp$la_log)) + mean(Qp$la_log)),
         Peak_discharge_Qp = 10^((Qp_log_scale * sd(Qp$Qp_log)) + mean(Qp$Qp_log)))

# We plot the posterior predictive distribution of Qp for a given lake area.

plot.Qp.la <- preds.Qp %>%
  ggplot(aes(x = Lake_area_before, y = Peak_discharge_Qp)) +
  scale_fill_manual(name = "Posterior rate", values = "lightgreen") +
  stat_lineribbon(.width = 0.95,
                  point_interval = mean_qi) +
  geom_point2(data = Qp,
              aes(x = Lake_area_before, y = Peak_discharge_Qp),
              alpha = 0.6) +
  theme_bw() +
  labs(y = "Reported GLOF peak discharge [m³ s-1]",
       x = "Lake area before GLOF [m²]") +
  scale_y_continuous(limits = c(10^-1, 10^7),
                     trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

plot.Qp.la  <- plot.Qp.la  +
  annotate(geom = "text",
           x = 10^3,
           y = 10^6,
           label = paste0("n = ", nrow(Qp)),
           hjust = 0,
           color = "grey20",
           size = 0.36 * 7) 

# We fit exactly the same models for the loss in lake area due to the GLOF.
# No further comments added.

V0.mod.la.diff <- brm(V0_log_scale ~ la_diff_log_scale,
                      data = V0,
                      prior = bprior, 
                      warmup  = 1000,
                      iter = 4000,
                      chains  = 4,
                      cores   = 4, 
                      control = list(adapt_delta = 0.95,
                                     max_treedepth = 15),
                      backend = "cmdstanr",
                      threads = threading(3))

# Posterior predictive distribution of the V0 model versus loss in lake area.

scaled.area <- (seq_range(c(3, 9), n = 200) - mean(V0$la_log)) / sd(V0$la_log)

preds.V0 <-  add_predicted_draws(
  object = V0.mod.la.diff, 
  newdata = data.frame(la_diff_log_scale = scaled.area),
  value = "V0_log_scale", 
  ndraws = 4000,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(la_diff =     10^((la_diff_log_scale * sd(V0$la_diff_log)) + mean(V0$la_diff_log)),
         Mean_Flood_Volume_V0 = 10^((V0_log_scale * sd(V0$V0_log)) + mean(V0$V0_log)))

plot.V0.la.diff <- preds.V0 %>%
  ggplot(aes(x = la_diff, y = Mean_Flood_Volume_V0)) +
  scale_fill_manual(name = "Posterior rate", values = "cyan") +
  stat_lineribbon(.width = 0.95,
                  point_interval = mean_qi) +
  geom_point2(data = V0,
              aes(x = la_diff, y = Mean_Flood_Volume_V0),
              alpha = 0.6) +
  theme_bw() +
  labs(y = "Reported GLOF volume [Mm³]",
       x = "Difference in lake area [m²]") +
  scale_y_continuous(limits = c(10^-4, 10^6.5),
                     trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

plot.V0.la.diff <- plot.V0.la.diff +
  annotate(geom = "text",
           x = 10^3,
           y = 10^6,
           label = paste0("n = ", nrow(V0)),
           hjust = 0,
           color = "grey20",
           size = 0.36 * 7) 

# Qp model

Qp.mod.diff <- brm(Qp_log_scale ~ la_diff_log_scale,
                   data = Qp,
                   prior = bprior, 
                   warmup  = 1000,
                   iter = 4000,
                   chains  = 4,
                   cores   = 4, 
                   control = list(adapt_delta = 0.95,
                                  max_treedepth = 15),
                   backend = "cmdstanr",
                   threads = threading(3))


# Posterior predictive distribution of the Qp model versus loss in lake area.

scaled.area <- (seq_range(c(3, 9), n = 200) - mean(Qp$la_diff_log)) / sd(Qp$la_diff_log)

preds.Qp <-  add_predicted_draws(
  object = Qp.mod.diff, 
  newdata = data.frame(la_diff_log_scale = scaled.area),
  value = "Qp_log_scale", 
  ndraws = 4000,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(Lake_area_before = 10^((la_diff_log_scale * sd(Qp$la_diff_log)) + mean(Qp$la_diff_log)),
         Peak_discharge_Qp = 10^((Qp_log_scale * sd(Qp$Qp_log)) + mean(Qp$Qp_log)))

plot.Qp.la.diff <- preds.Qp %>%
  ggplot(aes(x = Lake_area_before, y = Peak_discharge_Qp)) +
  scale_fill_manual(name = "Posterior rate", values = "lightgreen") +
  stat_lineribbon(.width = 0.95,
                  point_interval = mean_qi) +
  geom_point2(data = Qp,
              aes(x = Lake_area_before, y = Peak_discharge_Qp),
              alpha = 0.6) +
  theme_bw() +
  labs(y = "Reported GLOF peak discharge [m³ s-1]",
       x = "Difference in lake area [m²]") +
  scale_y_continuous(limits = c(10^-1, 10^7),
                     trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

plot.Qp.la.diff <- plot.Qp.la.diff +
  annotate(geom = "text",
           x = 10^3,
           y = 10^6,
           label = paste0("n = ", nrow(Qp)),
           hjust = 0,
           color = "grey20",
           size = 0.36 * 7) 

# Combine all four figures (V0 vs. lake area before the GLOF, 
#                           V0 vs. loss in lake area,
#                           Qp vs. lake area before the GLOF,
#                           Qp vs. loss in lake area).

Qp.V0.plot <- ggpubr::ggarrange(plot.V0.la, plot.V0.la.diff, 
                                plot.Qp.la, plot.Qp.la.diff, 
                                ncol = 2,
                                nrow = 2,
                                labels = c("A", "B", "C", "D"),
                                align = "hv",
                                font.label = list(size = 8,
                                                  color = "black",
                                                  face = "plain"))

# This will be Extended Data Figure 1.

ggsave("Qp_V0_plot.pdf",
       Qp.V0.plot,
       width = 130,
       height = 130,
       units = "mm")