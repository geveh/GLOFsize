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

# The goal of this script is to extract all lakes that are currently (2015 or 
# later) larger than 1 km², using existing glacier lake inventories.
# All lakes will then be manually mapped in QGIS given their geometric properties.

setwd("D:/data/BoxUP/Published_lake_databases")

# We analyse core regions in the global GLOF inventory.
# First, we list all inventories from all regions.

all.shp <- list.files(path = "D:/data/BoxUP/Published_lake_databases/Lakes_split",
                      pattern = ".gpkg$|.shp$", 
                      full.names = T, 
                      include.dirs = T,
                      recursive = T) 

# We read the vector files and keep only lakes that are today (2015 or later)
# larger than 1 km².

lakes.gt.1km <- lapply(all.shp, function (x) {
  
  s <- st_read(x) 
  
  uni.year <- unique(s$Year)
  uni.study <- unique(s$Study)
  
  if(uni.year < 2015) return()
  
  if(all(s$AreaUTM < 10^6)) return()
  
  s %>% 
    filter(AreaUTM > 10^6) %>%
    transmute(AreaUTM) %>%
    mutate(Region = basename(dirname(dirname(x))),
           Year = uni.year,
           Study = uni.study) %>%
    st_set_geometry("geometry")
  
})

# We remove all NULL entries (i.e. where no lake was mapped after 2015, or 
# all were smaller than 1 km²).

lakes.gt.1km.shp <- lakes.gt.1km[sapply(lakes.gt.1km,  function (x) !is.null(x))]

lakes.gt.1km.shp <- do.call(rbind, lakes.gt.1km.shp)

# We write this shapefile to disk.

st_write(lakes.gt.1km.shp, "lakes_gt_1km.gpkg", delete_dsn = T)

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

gt1.labelled <- st_read("lakes_gt_1km.gpkg")

# Besides these static predictors, we also analyse the individual changes in 
# lake area.

sf_use_s2(FALSE)

# We dissolve the boundary for all lakes that have the same ID.

gt1.uni <- gt1.labelled %>% 
  st_make_valid() %>%
  group_by(Uni_ID, Region) %>% 
  summarize() %>% 
  ungroup() %>% 
  distinct(geom, .keep_all = T)

# We list all folders that contain regional glacier lake inventories.

region.folders <- list.dirs("./Lakes_split", recursive = F)

# We now read all shapefiles from these folders to disk. The key criterion
# is to find out whether the region has inventories mapped between the late
# 1980s and 2015, because we would like to estimate the trend in lake growth 
# in the past decades. Importantly, we do not assume that a lake needs to be
# larger than 1 km², because it could have substantially grown in our study 
# period and just exceeded the 1 km² threshold in the past few years.

all.shpsbef2015 <- lapply(region.folders, function (x) {
  
  # We list the shapefiles in a given region.
  
  all.files <- list.files(path = x,   
                          pattern = ".gpkg$|.shp$", 
                          full.names = T, 
                          include.dirs = T,
                          recursive = T)
  
  # For every shapefile, we check when the lakes have been mapped.
  
  shapes <- lapply(all.files, function(y) {
    
    s <- st_read(y) 
    
    uni.year <- unique(s$Year)
    uni.study <- unique(s$Study)
    
    # We break the loop, if the shapefile was not created between 1980 and 2015.
    
    if(!((uni.year < 2015) & (uni.year > 1980))) return()
    
    s %>% 
      transmute(AreaUTM) %>%
      mutate(Region = basename(x),
             Year = uni.year,
             Study = uni.study) %>%
      st_set_geometry("geometry") }
    
  )
  
  # We remove empty entries from the list of shapefiles. Those occur, when
  # the if-condition above was FALSE. 
  
  shapes <- shapes[sapply(shapes, function(z) !is.null(z))]
  
  # We combine all shapefiles in that region.
  
  shapes.bef2015 <- do.call(rbind, shapes)
  
  return(shapes.bef2015)
  
})

# We combine all regional lake inventories to one large global lake inventory.

reg.shps <- do.call(rbind, all.shpsbef2015)

# For all lakes > 1 km² as of 2015, we now extract the overlapping lakes from the
# preceding years.

lake.ts.list <- list()

for(i in 1:nrow(gt1.uni)) {
  
  # We first select the corresponding region of that lake.
  
  all.lakes.reg <- reg.shps %>% 
    filter(Region == gt1.uni$Region[i])
  
  # Then, we search for lakes that intersect with the lake extent from 2015 or
  # later.
  
  int <- st_intersects(gt1.uni[i, ], all.lakes.reg, sparse = F)
  
  # If a lake consisted of several parts in a given year, we dissolve the
  # boundaries to one lake extent.
  
  overlapping.lakes  <- all.lakes.reg[int[1, ], ] %>% 
    st_make_valid() %>%
    group_by(Region, Year, Study) %>% 
    summarize(AreaUTM = sum(AreaUTM))
  
  # We also extract the original lake outlines as of 2015 or later (i.e.
  # before we dissolved them). 
  
  uni.lake.lab <- gt1.labelled %>%
    filter(Uni_ID == gt1.uni$Uni_ID[i]) %>%
    dplyr::select(c(AreaUTM, Region, Year, Study)) %>%
    st_set_geometry("geometry")
  
  # Now combine the entire time series of lake sizes into one simple feature.
  
  lake.ts.all <- rbind(overlapping.lakes,
                       uni.lake.lab) 
  
  # We do some rough filtering of outliers due to potential errors in the database.
  # All lakes that are two standard deviations smaller and larger than the
  # median lake area are removed from the time series.
  
  med.area <- median(lake.ts.all$AreaUTM)
  sd.area  <- sd(lake.ts.all$AreaUTM)
  
  lake.ts.clean <- lake.ts.all %>%
    filter(AreaUTM < (med.area + 2*sd.area),
           AreaUTM > (med.area - 2*sd.area)) %>%
    mutate(Uni_ID = gt1.uni$Uni_ID[i]) %>%
    ungroup() %>%
    st_drop_geometry()
  
  if (i %in% seq(1, 700, by = 25)) {message(paste0(i, " out of ", nrow(gt1.uni)))}
  
  lake.ts.list[[i]] <- lake.ts.clean
  
}

# We investigate whether the lakes had credible changes in lake area in our
# study period. We only focus on lakes that have more than 3 observations in 
# our study period.

all.lakes.ts <- do.call(rbind, lake.ts.list) %>%
  ungroup() %>%
  mutate(AreaUTM = log10(AreaUTM/ 10^6)) %>%
  filter(AreaUTM <= 250) %>%
  group_by(Uni_ID) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(x_scale = scale(Year),
         y_scale = scale(AreaUTM))

lake.trends.list <- list()

unique.ids <- sort(unique(all.lakes.ts$Uni_ID))

# We loop over all unique IDs and estimate whether there is a credible trend
# in lake area. We add a label for growth (G), unchanged (U), or decreasing (D).

for (u in unique.ids) {
  
  uni.id <- u
  
  uni.lake <- all.lakes.ts %>%
    filter(Uni_ID == uni.id)
  
  x <- uni.lake$Year %>% as_vector() %>% unname()
  y <- uni.lake$AreaUTM %>% as_vector() %>% unname()
  
  if (length(unique(x)) == 1) next
  if (length(unique(y)) == 1) next  
  
  # We scale predictor and response variables to a mean of zero and unit
  # standard deviation.
  
  x_scale <- scale(x)[ ,1]
  y_scale <- scale(y)[ ,1]
  
  # We obtain the mean and standard deviation on original scale. We need these
  # values to convert the data back to the original scale later.
  
  sd_y     <- sd(y)
  mean_y   <- mean(y) 
  sd_x   <- sd(x)
  mean_x <- mean(x)
  
  # We create a tibble that contains the ID, the lake area and the year as
  # standardised data.
  
  dat.new <- tibble(Uni_ID = factor(uni.lake$Uni_ID), 
                    x_scale = x_scale,
                    y_scale = y_scale)
  
  # We specify priors very weakly informed priors.
  
  bprior <- prior(student_t(3, 0, 2.5), class = "Intercept") +
    prior(student_t(3, 0, 5), class = "b") +
    prior(student_t(3, 0, 2.5), class = "sigma")
  
  # We model the trend GLOF size with time using a Student's t-distributed
  # likelihood.
  
  # The if-condition just assesses, if we're in the first iteration of the loop.
  
  if (u == min(unique.ids)) {
    
    mod1 <- brm(y_scale ~ x_scale, 
                family = student(),
                data = dat.new ,
                prior = bprior,
                cores  = 3,
                chains = 3,
                warmup = 1000,
                iter   = 4000,
                control = list(adapt_delta = 0.92,
                               max_treedepth = 15) ,
                backend = "cmdstanr",
                threads = threading(3))
    
    post.sum <- posterior_summary(mod1) 
    
    # Eventually, our only goal is to know whether the posterior slope (parameter
    # "b_x_scale") is credibly larger or smaller than zero, indicating a change
    # in lake area.
    # To this end, we add "G" for growth, "D" for decrease, and "U" for unchanged. 
    
    q_low  <- post.sum["b_x_scale", "Q2.5"]
    q_high <- post.sum["b_x_scale", "Q97.5"]
    
    if (q_high < 0) {trend <- "D"} else if (q_low > 0) {trend <- "G"} else { trend <- "U"}
    
  } else {
    
    # It's not necessary to recompile the model for every lake, because the 
    # model structure is always the same. Therefore, we just "update" the model
    # with new data.
    
    mod2 <- update(mod1, newdata = dat.new)
    
    post.sum <- posterior_summary(mod2) 
    
    q_low  <- post.sum["b_x_scale", "Q2.5"]
    q_high <- post.sum["b_x_scale", "Q97.5"]
    
    if (q_high < 0) {trend <- "D"} else if (q_low > 0) {trend <- "G"} else { trend <- "U"}
    
  }
  
  lake.trends.list[[u]] <- tibble(Uni_ID = u,
                                  trend = trend)
  
}

# We combine all trends in lake area to one long tibble.

lake.trends <- do.call(rbind, lake.trends.list) 

# We add the information on change in lake area to the labelled shapefile.
# We rearrange the data such that every region shows the share of values in each 
# of the six key criteria.

data <- gt1.labelled  %>% 
  left_join(., lake.trends, by = "Uni_ID") %>%
  distinct(geom, .keep_all = T) %>% 
  st_drop_geometry() %>%
  group_by(Region, Uni_ID) %>%
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
  mutate(Region = as_factor(Region),
         value  = as_factor(value),
         name   = as_factor(name)) %>%
  group_by(Region, value, name, .drop = F) %>%
  summarise(n_type = n()) %>% 
  ungroup() %>%
  group_by(Region, name) %>%
  arrange(desc(n_type)) %>%
  mutate(prop = n_type / sum(n_type) * 100) 

# We count the number of observations in each region.

tot.obs <- data %>%
  ungroup() %>%
  filter(name == "value_Managed") %>%
  group_by(Region) %>%
  summarise(tot_obs = sum(n_type))

#  We save this information to disk.

saveRDS(tot.obs, "all_observed_lakes_per_region.RDS")

data <- left_join(data, tot.obs, by = "Region") %>%
  mutate(Region = str_replace_all(Region, "_", " "))

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

glof.predictors.plot <- data %>% 
  ungroup() %>%
  mutate(value = factor(value, c("Y" , "G" , "I" , "R" , "U" , "N" , "M" , "D" )),
         name  = factor(name, c("value_lake_type", 
                                "value_trend",
                                "value_Glacier_contact",
                                "value_Estab_outflow", 
                                "value_Managed", 
                                "value_GLOF_bef90")),
         Region_n = paste0(Region, " (", tot_obs, ")")) %>%
  filter(name != "value_GLOF_bef90") %>%
  ggplot(aes(x = "", 
             y = prop, 
             fill = value,
             width = tot_obs)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  theme_bw() + 
  theme(legend.position = "none") +
  coord_flip() +
  facet_grid(vars(fct_reorder(Region_n, tot_obs, .desc = T)), 
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

data %>% filter(name == "value_lake_type") %>% 
  group_by(Region, value) %>% 
  summarise(s = sum(n_type)) %>% 
  filter(s != 0) %>% View()

################################################################################
### Statistics on controls in GLOF size ########################################

# Number of ice-dammed lakes gt 1 km².

gt1.labelled %>% 
  st_drop_geometry() %>%
  group_by(Uni_ID) %>%
  reframe(Dam = unique(Lake_type)) %>%
  group_by(Dam) %>%
  summarise(n())

# Number of reservoirs per region.

gt1.labelled %>% 
  st_drop_geometry() %>%
  group_by(Region, Uni_ID) %>%
  reframe(Dam = unique(Lake_type))  %>%
  group_by(Region, Dam) %>%
  summarise(n()) %>% View()

# Lakes detached from parent glacier.

gt1.labelled %>% 
  st_drop_geometry() %>%
  group_by(Uni_ID) %>%
  reframe(detached = unique(Glacier_contact))  %>%
  group_by(detached) %>%
  summarise(n())

# Lakes with managed outlets.

gt1.labelled %>% 
  st_drop_geometry() %>%
  group_by(Region, Uni_ID) %>%
  reframe(managed_outflow = unique(Managed))  %>%
  group_by(Region, managed_outflow) %>%
  summarise(n()) %>% View()

# Lakes with historical outbursts before 1990.

gt1.labelled %>% 
  st_drop_geometry() %>%
  group_by(Uni_ID) %>%
  reframe(hist_glof = unique(GLOF_bef90))  %>%
  group_by(hist_glof) %>%
  summarise(n())

# Number of lakes with positive trends.

data %>% 
  ungroup() %>%
  mutate(name = as.character(name)) %>% 
  filter(name == "value_trend",
         (value == "G") | (value == "U") | (value == "D") ) %>% 
  ungroup() %>% 
  group_by(value) %>% 
  summarise(sum(n_type))