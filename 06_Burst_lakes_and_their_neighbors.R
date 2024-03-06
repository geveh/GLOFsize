################################################################################
#######   Calculate distances between burst lakes and their neighbours    ######
#######                                                                   ######
#######                            by Georg Veh                           ######
#######                checked and comments added March 04, 2024          ######
################################################################################

# DISCLAIMER: This script is meant to make our findings reproducible. However,
# we do not redistribute the original lake inventories, which might be subject
# to differing licenses. Therefore, please contact the authors of the underlying
# studies to run this script.

# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

require(tidyverse)
require(sf)
require(units)
require(foreach)
require(parallel)
require(doParallel)

# The goal of this script is to find, for every burst lake, the next larger  
# intact lake, or in other words, the nearest lake that has not burst out yet.

# Set working directory

setwd("D:/data/BoxUP/Work 2022/GLOFsize")

# We load glacier buffers split by O2 regions.

buf.by.reg <- st_read("rgi06/glacier_buffers_split_by_O2_no_fid_correct_FULLNAME_2.gpkg") %>%
  st_make_valid() 

buf.by.reg.13reg <- buf.by.reg %>% 
  group_by(FULL_NAME) %>% 
  summarise() %>% 
  st_make_valid()

# We load all verified GLOF sites.

all.glofs.global <- readRDS("reported_GLOFs_with_geometry.rds")

# We analyse the core regions in the global GLOF inventory.
# We list all inventories from all regions.

all.shp <- list.files(path = "D:/data/BoxUP/Published_lake_databases/Lakes_split",
                      pattern = ".gpkg$|.shp$", 
                      full.names = T, 
                      include.dirs = T,
                      recursive = T) 

# For each lake in each glacier lake inventory, we want to know where the lake is 
# and what size it had. This will allow us to calculate distances between lakes
# and compare their sizes.

dfs <- lapply(all.shp, function (x) {
  
  s <- st_read(x) 
  tib <- tibble(Lat = s$Lat,
                Lon = s$Lon,
                Area = s$AreaUTM,
                Region = basename(dirname(dirname(x))),
                Year = s$Year,
                Study = s$Study)
  
})

# Bind the regional lake inventories and convert them to a point vector file 
# according to the centroidal coordinate.

dfs.bind <- do.call(rbind, dfs) %>%
  st_as_sf(coords = c("Lon", "Lat"), 
           crs = 4326)

# We remove all lakes that could have been counted twice, for whatever reason.

dfs.distinct <- distinct(dfs.bind) 

# How many lakes were mapped?

nrow(dfs.distinct)

# If a lake had repeated outbursts, we select the largest outburst (i.e. area 
# before the GLOF) per site.

glof.sources <- all.glofs.global %>%
  group_by(geometry) %>%
  slice_max(Lake_area_before) %>%
  ungroup()

# We now calculate the distance between every burst lake and the nearest larger 
# intact lake. We iterate over every burst lake in a parallel loop. 
# Set the number of clusters to the number of cores on your machine, but make
# sure to leave at least one core free for other tasks.

cl <- makePSOCKcluster(7)
registerDoParallel(cl)

lake.distances <- foreach(i = 1:nrow(glof.sources), 
                          .packages=c('sf', 'dplyr', 'units'),
                          .combine = 'rbind') %dopar% {
                            
                            # We begin by drawing a buffer of 100 km around the  
                            # GLOF source location to find potential region(s)
                            # to look for unburst lakes in the neighborhood of 
                            # the burst lake.
                            
                            glof.buf.100 <- st_buffer(glof.sources[i, ], 100000)
                            
                            # We select only those regions in this buffer.
                            
                            reg.int.tf <- st_intersects(glof.buf.100, 
                                                        buf.by.reg.13reg, 
                                                        sparse = F)
                            
                            overlapping.regions <- buf.by.reg.13reg$FULL_NAME[reg.int.tf[1, ]]
                            
                            # We select only inventories in these regions that 
                            # mapped lakes as of 2015 or later.
                            
                            dfs.distinct.sub <- dfs.distinct[dfs.distinct$Region %in% overlapping.regions, ] %>%
                              filter(Year >= 2015)
                            
                            # We remove all lakes that burst out in that 100 km 
                            # buffer.
                            
                            dist.to.all.glofs <- st_distance(glof.sources[i, ], glof.sources) %>% 
                              drop_units() 
                            
                            glof.sources.in.reg <- glof.sources[dist.to.all.glofs[1, ] < 100000, ]
                            
                            # We loop over these burst lakes, draw a buffer of 
                            # 1000 m, and remove all lakes within that search
                            # radius from the total lake population.
                            
                            for(g in 1:nrow(glof.sources.in.reg)) {
                              
                              glof.sources.in.reg[g, ]
                              
                              dist.to.g <- st_distance(glof.sources.in.reg[g, ], dfs.distinct.sub)
                              
                              dfs.distinct.sub <- dfs.distinct.sub[drop_units(dist.to.g[1, ]) > 1000, ]
                              
                            }
                            
                            # Now we calculate the distance between that burst 
                            # lake and the entire lake population in its surrounding.
                            
                            dists <- st_distance(glof.sources[i, ], dfs.distinct.sub)
                            dfs.distinct.sub$dists <- drop_units(dists[1, ])
                            
                            # We return the lake that is larger and has the  
                            # closest distance to a burst lake.
                            
                            return(dfs.distinct.sub %>% 
                                     filter(Area > glof.sources$Lake_area_before[i],
                                            dists > 1000) %>%
                                     arrange(dists) %>%
                                     slice_head() %>%
                                     mutate(Lake_type_simple = glof.sources$Lake_type_simple[i]))
                            
                          }

# We stop the cluster.

stopCluster(cl)

# We convert the distances from meters to kilometers.

ld <- lake.distances %>% 
  mutate(value = dists/1000)

# We add the global distances (i.e. all distances from all regions).

ld <- bind_rows(ld,
                ld %>% mutate(Region = "Global"))

# We show the number of lakes that had at least one outburst per region.

nobs.ld <- ld %>%
  st_drop_geometry() %>%
  group_by(Region) %>%
  summarise(nobs = n()) %>%
  st_drop_geometry() %>%
  mutate(Region = str_replace_all(Region, "_", " "))

# We plot the distances to the next intact glacier lake as histogram and
# add the individual data as jittered points on top. 

distances.plot <-  ld %>%
  st_drop_geometry() %>%
  mutate(Lake_type_simple = case_when(Lake_type_simple == "glacier_supraglacial" ~ "Glacier & supraglacial",
                                      Lake_type_simple == "moraine_bedrock" ~ "Moraine & bedrock"),
         Region = str_replace_all(Region, "_", " ")) %>%
  left_join(x = ., y = nobs.ld, by = "Region") %>%
  mutate(reg_obs = paste0(Region, " (", nobs, ")")) %>%
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
  ggplot(aes(x = value, y = region2)) +
  geom_boxplot(fill = "azure2", 
               outlier.shape = NA, 
               width = 0.5,
               linewidth = 0.2) +
  geom_jitter2(aes(color = Lake_type_simple), 
               size = 0.5, height = 0.15,
               alpha = 0.5) +
  scale_colour_manual(name = "Dam type",
                      values = c("Glacier & supraglacial" = "navy", 
                                 "Moraine & bedrock" = "#ee7600"))+ 
  theme_bw() +
  scale_x_continuous(limits = c(1, 1000),
                     trans  = log10_trans(),
                     breaks = c(1, 10, 100, 1000),
                     labels = label_number(drop0trailing = TRUE)) +
  theme( axis.text    = element_text(size = 7),
         axis.text.x  = element_text(size = 7),
         axis.title   = element_text(size = 7),
         legend.title = element_text(size = 7),
         legend.text  = element_text(size = 7),
         legend.position = "bottom") +
  xlab("Distance of burst lakes to\nthe next larger intact lake [km]") +
  ylab("") + 
  guides(colour = guide_legend(ncol = 2))

# We calculate the median distance to the next larger lake.

ld %>% 
  st_drop_geometry() %>% 
  group_by(Region) %>% 
  summarise(m_dist = median(value))

# We save the plot to disk (Figure 4B).

ggsave(filename = "distances_plot.pdf",
       distances.plot,
       width = 180,
       height = 110,
       units = "mm")
