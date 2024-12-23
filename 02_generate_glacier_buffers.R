################################################################################
#######         Draw 5-km buffers around glaciers in the RGI 6.0          ######
#######                                                                   ######
#######             checked and comments added March 04, 2024             ######
#######                   checked again Dec 23, 2024                      ######
################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

require(sf)
require(tidyverse)
library(doParallel)

# Set YOUR working directory folder where to find all files, necessary to run 
# this script. Change the location appropriately. To run this script on your own
# machine, you need to download the entire Randolph Glacier Inventory V6.0 
# (https://www.glims.org/RGI/rgi60_dl.html) and set the unzipped directory as
# your working directory.

setwd("D:/data/BoxUP/Work 2022/GLOFsize/rgi06")

# We read the RGI O2 regions to memory. We merge the original regions to larger 
# regions given that some regions have only few glaciers, glacier lakes, and 
# reported GLOFs. To do so, we need the RGI60 in a folder (here called rgi06) on
# our machine. 

rgi <- st_read("./00_rgi60_regions/00_rgi60_O2Regions.shp") %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  mutate(FULL_NAME = str_replace_all(FULL_NAME, "\\(|\\)", "")) %>%
  mutate(
    FULL_NAME = str_replace_all(FULL_NAME, "C Andes", "C_and_N_Andes"),
    FULL_NAME = str_replace_all(FULL_NAME, "Low-latitude Andes", "C_and_N_Andes"),
    FULL_NAME = str_replace_all(FULL_NAME, "C Himalaya", "Himalayas"),
    FULL_NAME = str_replace_all(FULL_NAME, "E Himalaya", "Himalayas"),
    FULL_NAME = str_replace_all(FULL_NAME, "W Himalaya", "Himalayas"),
    FULL_NAME = str_replace_all(FULL_NAME, "C Himalaya", "Himalayas"),
    FULL_NAME = str_replace_all(FULL_NAME, "Karakoram", "HK_Pamir_Karakoram"),
    FULL_NAME = str_replace_all(FULL_NAME, "Pamir Safed Khirs/W Tarim", "HK_Pamir_Karakoram"),
    FULL_NAME = str_replace_all(FULL_NAME, "Hindu Kush", "HK_Pamir_Karakoram"),
    FULL_NAME = str_replace_all(FULL_NAME, "Hengduan Shan", "Tibet_and_Hengduan_Shan"),
    FULL_NAME = str_replace_all(FULL_NAME, "Inner Tibet", "Tibet_and_Hengduan_Shan"),
    FULL_NAME = str_replace_all(FULL_NAME, "S and E Tibet", "Tibet_and_Hengduan_Shan"),
    FULL_NAME = str_replace_all(FULL_NAME, "SE Scandinavia", "Scandinavia"),
    FULL_NAME = str_replace_all(FULL_NAME, "SW Scandinavia", "Scandinavia"),
    FULL_NAME = str_replace_all(FULL_NAME, "N Scandinavia", "Scandinavia"),
    FULL_NAME = str_replace_all(FULL_NAME, "Hissar Alay", "Hissar_Alay_and_Tien_Shan"),
    FULL_NAME = str_replace_all(FULL_NAME, "W Tien Shan", "Hissar_Alay_and_Tien_Shan"),
    FULL_NAME = str_replace_all(FULL_NAME, "E Tien Shan Dzhungaria", "Hissar_Alay_and_Tien_Shan"),
    FULL_NAME = str_replace_all(FULL_NAME, "N Coast Ranges", "Coast_Ranges"),
    FULL_NAME = str_replace_all(FULL_NAME, "S Coast Ranges", "Coast_Ranges"),
    FULL_NAME = str_replace_all(FULL_NAME, "W Chugach Mtns Talkeetna", "W_Chugach_Mountains"),
    FULL_NAME = str_replace_all(FULL_NAME, "Alaska Ra Wrangell/Kilbuck", "Alaska_Range"),
    FULL_NAME = str_replace_all(FULL_NAME, "St Elias Mtns", "Saint_Elias_Mountains")) 

# We only keep the 13 study regions, and remove all others.

filt <- rgi$FULL_NAME %in%
  c("Alps", "Iceland", "C_and_N_Andes", "Himalayas", 
    "HK_Pamir_Karakoram", "Tibet_and_Hengduan_Shan", "Scandinavia",
    "Hissar_Alay_and_Tien_Shan", "Coast_Ranges", "W_Chugach_Mountains", 
    "Alaska_Range", "Saint_Elias_Mountains", "Patagonia") 

rgi <- rgi[filt, ]

# We dissolve regions (i.e. remove boundaries for regions that have the same 
# value in the attribute 'FULL_NAME').

rgi.dissolve <- group_by(rgi, FULL_NAME) %>% 
  summarise()

# We write these regions to disk.

st_write(rgi.dissolve, "rgiO2_dissolved_outlines.shp")

# We list all folders that contain polygons of regional glacier inventories.

glacier.folders <- list.dirs(recursive = F)

# We create a folder in which we will drop dissolved buffers around all glaciers 
# per region.

dir.create("glacier_buffers")

# We iterate over the regions in parallel, i.e. every region is processed on
# one core on the computer. This skript was run on a machine with 8 cores. If
# you have more or fewer cores, please adjust to the number of cores minus one.

ncores <- 8 - 1

cl <- makePSOCKcluster(ncores)
registerDoParallel(cl)

buffer.list <- foreach(i = 1:nrow(rgi), 
                       .packages = c('sf', 'stringr', 'dplyr')) %dopar% {
                         
                         # We obtain the regional RGI code for O1 and O2 regions.
                         
                         O1 <- str_sub(rgi$RGI_CODE[i], 1, 2)
                         O2 <- str_sub(rgi$RGI_CODE[i], 4, 5)
                         
                         # We select the folder that contains the data 
                         # for the corresponding O1 region.
                         
                         rgi.glacier.folder <- glacier.folders[str_detect(glacier.folders, O1)]
                         
                         # We read the glacier shapefile into memory.
                         
                         rgi.file <- list.files(path = rgi.glacier.folder, 
                                                pattern = ".shp$", 
                                                full.names = T) %>% 
                           st_read() %>%
                           st_make_valid()
                         
                         # We ensure that the glaciers intersect with the O1 region.
                         
                         filt <- st_intersects(rgi[i, ], rgi.file, sparse = F)
                         
                         sub.glaciers <- rgi.file[filt[1, ], ]
                         
                         # We draw a buffer around the glaciers in that region.
                         
                         sub.glaciers.buffer <- st_buffer(sub.glaciers, dist = 5000) %>% 
                           st_make_valid() %>%
                           st_union() %>%
                           st_make_valid() %>%
                           st_as_sf() %>%
                           mutate(FULL_NAME = rgi$FULL_NAME[i],
                                  RGI_CODE = rgi$RGI_CODE[i])
                         
                         # Finally, we write the buffer to disk.
                         
                         st_write(sub.glaciers.buffer, paste0("./rgi06/glacier_buffers/", 
                                                              rgi$RGI_CODE[i], 
                                                              "_buffer.shp"))
                         
                         # We also return output from the loop. Every regional
                         # glacier buffer will be a separate list item. 
                         
                         return(sub.glaciers.buffer)
                         
                       }


# We close all connections to the cluster and free some memory.

stopCluster(cl)
gc()

# We combine the list output from the foreach-loop to one large polygon that
# contains the glacier buffers for each region.

buffer.shp <- do.call(rbind, buffer.list)

# We dissolve (sf function 'dissolve') the buffered glaciers the attribute
# 'FULL_NAME'.

buffer.dissolve <- group_by(buffer.shp, FULL_NAME) %>% 
  summarise()

# We finally write the glacier buffers to disk.

st_write(buffer.dissolve, "dissolved_buffer.shp")

saveRDS(buffer.dissolve, 
        "dissolved_buffer.RDS")

# Some of the glacier buffers overlap. In QGIS, we clipped all overlapping buffers
# using the outlines of the RGI O2 regions, and only retained one of the 
# overlapping parts. This file is called 
# "glacier_buffers_split_by_O2_no_fid_correct_FULLNAME_2.gpkg".
