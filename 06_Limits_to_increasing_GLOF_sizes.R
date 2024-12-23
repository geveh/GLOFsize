################################################################################
#######                       Extract all lakes >1km² and                 ######
#######                  assess limits to increasing GLOF sizes           ######
#######                            by Georg Veh                           ######
#######                                                                   ######
#######                checked and comments added March 06, 2024          ######
################################################################################

# DISCLAIMER: This script is meant to make our findings reproducible. However,
# we do not redistribute the original lake inventories, which might be subject
# to differing licenses. Therefore, please contact the authors of the underlying
# studies to run this script.

# Important: we run brms models using the cmdstan backend. This needs to be
# installed separately using the instructions given here: 
# https://mc-stan.org/cmdstanr/ 

# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

require(tidyverse)
require(sf)
require(brms)
require(terra)

# The goal of this script is to extract all lakes that are currently (2015 or 
# later) larger than 1 km², using existing glacier lake inventories.
# All lakes will then be manually mapped in QGIS given their geometric properties.

setwd("D:/data/BoxUP/Published_lake_databases")

# We labelled every lake in QGIS using high-resolution aerial images from Bing 
# Basemaps, satellite images in Google Earth, and Planet images to...
# - add a unique ID for overlapping lakes, if the same lake was mapped in more
#   than one inventory;
# - define the dam type: I = glacier and supraglacial, M = moraine & bedrock,
#                        R = Reservoir for hydropower or freshwater;
# - check whether the lake is in contact with the glacier (Y) or not (N). Glacier
#   contact refers not to the dam itself, but whether any other margin of the
#   glacier has a direct contact to an adjacent glacier;
# - check whether the lake has an established outflow (Y), i.e. a continuous, ideally
#   perennial stream of water leaving the lake, or not (N);
# - whether the outflow is managed (Y) or not (N) by a dyke, artificial dam or
#   tunnel;
# - whether the lake had a reported GLOF before our study period (Y) or not (N).

# Read the labelled inventory of large lakes back to memory.

gt1.labelled <- st_read("Zhang_lakes_2020_gt1km_join_region.gpkg") %>%
  mutate(Uni_ID = 1: nrow(.))

rgi.regions <- st_read("D:/data/BoxUP/Work 2022/GLOFsize/rgi06/dissolved_buffer.shp")

# Besides these static predictors, we also analyse the individual changes in 
# lake area.

sf_use_s2(FALSE)

# We list all folders that contain regional glacier lake inventories.

region.folders <- list.dirs("./Zhang_split", recursive = F)

# We now read all shapefiles from these folders to disk. The key criterion
# is to find out whether the region has inventories mapped between the late
# 1980s and 2015, because we would like to estimate the trend in lake growth 
# in the past decades. Importantly, we do not assume that a lake needs to be
# larger than 1 km², because it could have substantially grown in our study 
# period and just exceeded the 1 km² threshold in the past few years.

all.shps.1990 <- lapply(region.folders, function (x) {
  
  # We list the shapefiles in a given region.
  
  all.files <- list.files(path = x,   
                          pattern = "_1990.shp$", 
                          full.names = T, 
                          include.dirs = T,
                          recursive = T) 
  
  all.files <- lapply(all.files, st_read) %>% 
    do.call(rbind, .)
  
  all.files.regions <- st_join(
    all.files,
    rgi.regions,
    join = st_intersects,
    suffix = c("", ""),
    left = TRUE,
    largest = TRUE
  )
  
  return(all.files.regions)
  
})

# We combine all regional lake inventories to one large global lake inventory.

reg.shps <- do.call(rbind, all.shps.1990)

saveRDS(reg.shps, "reg_shps.RDS")
# reg.shps <- readRDS("reg_shps.RDS")

# For all lakes > 1 km² as of 2020, we now extract the overlapping lakes from the
# preceding years.

lake.ts.list <- list()

for(i in 1:nrow(gt1.labelled)) {
  
  # We first select the corresponding region of that lake.
  
  all.lakes.reg <- reg.shps %>% 
    filter(FULL_NAME == gt1.labelled$FULL_NAME[i])
  
  # Then, we search for lakes that intersect with the lake extent from 2020.
  
  int <- st_intersects(gt1.labelled[i, ], all.lakes.reg, sparse = F)
  
  # If a lake consisted of several parts in a given year, we dissolve the
  # boundaries to one lake extent.
  
  if(!(all(int[1, ] == F))) {
  
  overlapping.lakes  <- all.lakes.reg[int[1, ], ] %>% 
    st_make_valid() %>%
    summarize(Area = sum(Area),
              Year = mean(Year))
  
  # Now combine the entire time series of lake sizes into one simple feature.
  
  lake.ts.all <- rbind(overlapping.lakes,
                       gt1.labelled[i, ] %>% 
                         dplyr::select(c(Area, Year)) %>% 
                         st_set_geometry("geometry")) 
  
  lake.ts.clean <- lake.ts.all %>%
    mutate(Uni_ID = i) %>%
    ungroup() %>%
    st_drop_geometry() 
  
    } else {
      
      lake.ts.clean <- gt1.labelled[i, ] %>% 
        dplyr::select(c(Area, Year)) %>%
        st_drop_geometry() %>%
        add_row(Area = 0, Year = 1990) %>%
        mutate(Uni_ID = i)
       
    }
  
  if (i %in% seq(1, nrow(gt1.labelled), by = 25)) {message(paste0(i, " out of ", nrow(gt1.labelled)))}
  
  lake.ts.list[[i]] <- lake.ts.clean
  
}


lake.ts.list <- do.call(rbind, lake.ts.list) 

unique.ids <- sort(unique(lake.ts.list$Uni_ID))

lake.trends.list <- list()

# We loop over all unique IDs and estimate whether there is a credible trend
# in lake area. We add a label for growth (G), unchanged (U), or decreasing (D).

for (u in unique.ids) {
  
  uni.id <- u
  
  uni.lake <- lake.ts.list %>%
    filter(Uni_ID == uni.id)
  
  start.val <- uni.lake$Area[uni.lake$Year < 2015]
  end.val   <- uni.lake$Area[uni.lake$Year > 2015]
  
  perc.ch <- (end.val - start.val)/ start.val *100
  
    if (perc.ch < -2.5) {trend <- "D"
    } else if ((perc.ch > -2.5) & (perc.ch < 2.5)) {
       trend <- "U"} else if (perc.ch > 2.5) { 
         trend <- "G"} else if (is.infinite(perc.ch)) {
           trend <- "G"
         }

  
  lake.trends.list[[u]] <- tibble(Uni_ID = u,
                                  trend = trend)
  
  }

# We combine all trends in lake area to one long tibble.

lake.trends <- do.call(rbind, lake.trends.list) 

# Select lakes without reported outbursts ('unburst') between 1990 and 2023.

data <- gt1.labelled  %>% 
  left_join(., lake.trends, by = "Uni_ID") %>%
  distinct(geom, .keep_all = T)

rep.glof <- readRDS("D:/data/BoxUP/Work 2022/GLOFsize/la_sf.RDS")

rep.glof.buf <- rep.glof %>% 
  vect() %>%
  buffer(width = 500) %>% 
  st_as_sf()

tfvec <- st_intersects(data, rep.glof.buf, sparse = F) %>%
  apply(MARGIN = 1, function (x) !any(x == T))

data <- data[tfvec, ]

# We add the information on change in lake area to the labelled shapefile.
# We rearrange the data such that every region shows the share of values in each 
# of the six key criteria.

data.rearranged <- data %>% 
  st_drop_geometry() %>%
  group_by(FULL_NAME, Uni_ID) %>%
  reframe(value_lake_type = unique(Lake_type),
          value_Glacier_contact = unique(Glacier_contact),
          value_GLOF_bef90 = unique(GLOF_bef90), 
          value_Managed = unique(Managed),
          value_Estab_outflow = unique(Estab_outflow),
          value_trend = unique(trend)) %>% 
  ungroup() %>% 
  pivot_longer(cols = c(value_lake_type,
                        value_Glacier_contact,
                        value_GLOF_bef90,
                        value_Managed,
                        value_Estab_outflow, 
                        value_trend)) %>%
  mutate(FULL_NAME = as_factor(FULL_NAME),
         value  = as_factor(value),
         name   = as_factor(name)) %>%
  group_by(FULL_NAME, value, name, .drop = F) %>%
  summarise(n_type = n()) %>% 
  ungroup() %>%
  group_by(FULL_NAME, name) %>%
  arrange(desc(n_type)) %>%
  mutate(prop = n_type / sum(n_type) * 100) 

# We count the number of observations in each region.

tot.obs <- data.rearranged %>%
  ungroup() %>%
  filter(name == "value_Managed") %>%
  group_by(FULL_NAME) %>%
  summarise(tot_obs = sum(n_type))

#  We save this information to disk.

saveRDS(tot.obs, "all_observed_lakes_per_region.RDS")

data.rearranged <- left_join(data.rearranged, tot.obs, by = "FULL_NAME") %>%
  mutate(FULL_NAME = str_replace_all(FULL_NAME, "_", " "))

saveRDS(data, "glof_predictors.RDS")

# data <- readRDS("glof_predictors.RDS")
# tot.obs <- readRDS("all_observed_lakes_per_region.RDS")

# Every value has a different color.

cols <- c("Y" = "#fde725", 
          "G" = "blue",
          "I" = "navy",
          "R" = "#21918c",
          "U" = "grey",
          "N" = "#440154", 
          "M" = "#ee7600", 
          "D" = "red")

# We finally plot all the controls on GLOF size in one barplot. Every predictor
# should be one column. The values in the columns are scaled between 0 and 100,
# as they show the share of each group in a given predictor.

glof.predictors.plot <- data.rearranged %>% 
  ungroup() %>%
  mutate(value = factor(value, c("Y" , "G" , "I" , "R" , "U" , "N" , "M" , "D" )),
         name  = factor(name, c("value_lake_type", 
                                "value_trend",
                                "value_Glacier_contact",
                                "value_Estab_outflow", 
                                "value_Managed", 
                                "value_GLOF_bef90")),
         FULL_NAME_n = paste0(FULL_NAME, " (", tot_obs, ")")) %>%
  filter(name != "value_GLOF_bef90") %>%
  ggplot(aes(x = "", 
             y = prop, 
             fill = value,
             width = tot_obs)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() + 
  theme(legend.position = "none") +
  coord_flip() +
  facet_grid(vars(fct_reorder(FULL_NAME_n, tot_obs, .desc = T)), 
             vars(name), 
             labeller = labeller(name = c(value_lake_type = "Lake type",
                                          value_Glacier_contact = "Contact to\nglacier",
                                          value_GLOF_bef90 = "GLOF\nbefore 1990",
                                          value_Managed = "Managed\noutlet",
                                          value_Estab_outflow = "Established\noutflow",
                                          value_trend = "Lake area\nchange")),
             space = "free", 
             scales = "free",
             switch = "y") +
  theme_bw() +
  geom_hline(yintercept = 50, linetype = 1, 
             color = "grey50", linewidth = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = cols, na.value = "grey90") +
  theme( strip.background = element_blank(),
         axis.text    = element_text(size = 7),
         axis.text.x  = element_text(size = 7),
         axis.title   = element_text(size = 7),
         strip.text.y = element_text(size = 7, hjust = 1),
         strip.text.y.left = element_text(angle = 0),
         strip.text.x = element_text(size = 7),
         plot.title   = element_text(size = 7, hjust = 0.5, angle = 0),
         legend.title = element_text(size = 7),
         legend.text  = element_text(size = 7),
         legend.position = "bottom",
         panel.spacing = unit(0.2, "lines")) +
  xlab("") +
  ylab("") + 
  guides(colour = guide_legend(ncol = 2))

# We save this plot (Figure 5).

ggsave(
  filename = "glof_predictors_plot.pdf",
  plot = glof.predictors.plot ,
  width = 180,
  height = 150,
  units = "mm"
)


################################################################################
### Statistics on controls in GLOF size ########################################

data.rearranged %>% 
  filter(name == "value_lake_type") %>% 
  group_by(FULL_NAME, value) %>% 
  summarise(s = sum(n_type)) %>% 
  filter(s != 0) %>% View()

# Number of ice-dammed lakes gt 1 km².

data %>% 
  st_drop_geometry() %>%
  group_by(Uni_ID) %>%
  reframe(Dam = unique(Lake_type)) %>%
  group_by(Dam) %>%
  summarise(n())

# Number of reservoirs per region.

data %>% 
  st_drop_geometry() %>%
  group_by(FULL_NAME, Uni_ID) %>%
  reframe(Dam = unique(Lake_type))  %>%
  group_by(FULL_NAME, Dam) %>%
  summarise(n()) %>% View()

# Lakes detached from parent glacier.

data %>% 
  st_drop_geometry() %>%
  group_by(Uni_ID) %>%
  reframe(detached = unique(Glacier_contact))  %>%
  group_by(detached) %>%
  summarise(n())

# Lakes with managed outlets.

data %>% 
  st_drop_geometry() %>%
  group_by(FULL_NAME, Uni_ID) %>%
  reframe(managed_outflow = unique(Managed))  %>%
  group_by(FULL_NAME, managed_outflow) %>%
  summarise(n()) %>% View()

# Lakes with historical outbursts before 1990.

data %>% 
  st_drop_geometry() %>%
  group_by(Uni_ID) %>%
  reframe(hist_glof = unique(GLOF_bef90))  %>%
  group_by(hist_glof) %>%
  summarise(n())

# Number of lakes with estab. outflows per region

estab.outflow.per.region <- data %>% 
  st_drop_geometry() %>% 
  group_by(FULL_NAME) %>% 
  filter(Estab_outflow == "Y") %>% 
  summarise(n_estab = n())

n.per.region <- data %>% 
  st_drop_geometry() %>% 
  group_by(FULL_NAME) %>% 
  summarise(n_lakes = n())

left_join(estab.outflow.per.region, n.per.region, by = "FULL_NAME") %>%
  mutate(ratio = (n_estab/ n_lakes)*100)

reservoirs.per.region <- data %>% 
  st_drop_geometry() %>% 
  group_by(FULL_NAME) %>% 
  filter(Lake_type == "R") %>% 
  summarise(n_estab = n())
