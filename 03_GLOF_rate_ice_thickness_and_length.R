################################################################################
#######    Regional GLOF rate and local glacier thickness and length      ######
#######                compared to the size of burst lakes                ######
#######                                                                   ######
#######                            by Georg Veh                           ######
#######             checked and comments added March 04, 2024             ######
################################################################################

# Load the following packages, or use install.packages("nameofpackage"), if some 
# of them are not pre-installed. In some cases you need to restart your R session.

require(tidyverse)
require(scales)
require(sf)
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

# Load the preprocessed table of reported GLOFs.

la.sf <- readRDS("la_sf.RDS")

################################################################################

# We generate a table of reported GLOFs that identifies each lake with its
# parent glacier. We replace special cases such as blank spaces, dots or hyphens
# in this new column, and in the RGI column, because those cause trouble in 
# the brms model later.

rep.glof <- la.sf %>%
  st_drop_geometry() %>%
  mutate(glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  filter(!is.na(Lake_area_before)) %>%
  group_by(Lake_type_simple, glacier_and_lake) %>%
  ungroup() %>%
  mutate(area_scale = scale_this(Lake_area_before),
         year_scale  = scale_this(rounded_year)) %>%
  mutate(glacier_and_lake = str_replace_all(glacier_and_lake, pattern = " ", replacement = "_"),
         glacier_and_lake = str_replace_all(glacier_and_lake, pattern = "[[.]]", replacement = "_"),
         glacier_and_lake = str_replace_all(glacier_and_lake, pattern = "-", replacement = "_"),
         RegO2_adj_f      = str_replace_all(as.character(RegO2_adj_f), " ", "_")) 

# Save this table to disk. We need it later.

saveRDS(rep.glof, "reported_GLOFs.rds")

# For cases where we need the geometry (i.e. the location of the GLOF site),
# we also export an RDS file that includes the coordinate of the GLOF.

la.sf %>%
  mutate(glacier_and_lake = paste0(RGI_Glacier_Id, "_", Lake)) %>%
  filter(!is.na(Lake_area_before)) %>%
  group_by(Lake_type_simple, glacier_and_lake) %>%
  ungroup() %>%
  mutate(area_scale = scale_this(Lake_area_before),
         year_scale  = scale_this(rounded_year)) %>%
  mutate(glacier_and_lake = str_replace_all(glacier_and_lake, pattern = " ", replacement = "_"),
         glacier_and_lake = str_replace_all(glacier_and_lake, pattern = "[[.]]", replacement = "_"),
         glacier_and_lake = str_replace_all(glacier_and_lake, pattern = "-", replacement = "_"),
         RegO2_adj_f      = str_replace_all(as.character(RegO2_adj_f), " ", "_")) %>%
  saveRDS("reported_GLOFs_with_geometry.rds")


################################################################################
### GLOF size versus GLOF reporting rate #######################################

# We calculate the regional reporting rate for both dam types, 
# as well as the regional GLOF size (25th, 50th, and 75th percentile of the lake
# area before the GLOF).

reg.lake.sizes <- rep.glof %>%
  group_by(RegO2_adj_f, Lake_type_simple) %>%
  filter(n() >= 5) %>%
  summarise(freq = n()/ length(1990:2023),
            med = quantile(Lake_area_before, 0.5),
            q25 = quantile(Lake_area_before, 0.25),
            q75 = quantile(Lake_area_before, 0.75)) %>%
  arrange(desc(freq))

# What is the average regional GLOF reporting rate, distinguished by dam type?

reg.lake.sizes %>% 
  ungroup() %>%
  group_by(Lake_type_simple) %>%
  summarise(med_freq = median(freq))

# We plot the rate of reported GLOFs versus the lake area before the GLOF,
# distinguished by region and dam type.

mag.freq <- reg.lake.sizes %>%
  mutate(Lake_type_simple = str_replace_all(Lake_type_simple, 
                                            "glacier_supraglacial", 
                                            "Glacier & supraglacial"),
         Lake_type_simple = str_replace_all(Lake_type_simple, 
                                            "moraine_bedrock", 
                                            "Moraine & bedrock")) %>%
  ggplot(mapping = aes(x = freq, 
                       y = med,
                       color = Lake_type_simple)) +
  geom_point2(alpha = 0.7) +
  geom_linerange(aes(x = freq, 
                     ymin = q25,
                     ymax = q75,
                     alpha = 0.7)) +
  scale_color_manual(name = "Dam type", 
                     values = c("navy",  "#ee7600")) +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  xlab("Annual rate of reported GLOFs per region") +
  ylab("Lake area before the outburst [m²] by region") + 
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.title= element_text(size = 7, face = "bold"), 
         legend.text = element_text(size = 7),
         legend.position = c(0.65, 0.2), 
         aspect.ratio = 1) 

# We add a marginal boxplot on top to show differences in reporting rates
# between dam types.

mag.freq.mar <- ggMarginal(mag.freq, 
                           type = "boxplot",
                           margins = "x",
                           groupColour = T,
                           size = 7,
                           lwd = 0.5)


################################################################################
### GLOF size versus glacier thickness #########################################

# If a lake had several outbursts, we select only the largest in our study period.

max.lake <- rep.glof %>%
  group_by(RGI_Glacier_Id, Lake, Lake_type_simple) %>%
  summarise(max_area = max(Lake_area_before))

# Daniel Farinotti provided consensus estimates of glacier thickness and volume
# in table format. The original publication is: Farinotti et al. (2019):
# "A consensus estimate for the ice thickness distribution of all glaciers on 
# Earth." Nature Geoscience 12, 3, 168-173.
# The table tells us the mean thickness ("Hcomp") of each glacier glacier 
# according to its RGIId.

thickness.data <- read_table(file = "RGI-wide_composites_stats_GV.txt")

# Read the 5 km buffer around glaciers into memory.

buf.o2 <- st_read("rgi06/glacier_buffers_split_by_O2_no_fid_correct_FULLNAME_2.gpkg")

# List all vector outlines of glaciers in this folder.

f.list <- list.files(pattern = "_rgi60_", recursive = T, full.names = T) %>%
  as_tibble() %>%
  filter(grepl(".shp", value)) %>%
  filter(!grepl("00", value))

# We iterate over all glacier shapefiles to obtain the glacier name and
# its region. We then add the corresponding ice thickness to each glacier.
# The output will be a table of all glaciers in each region that had an estimate
# of ice thickness.

reg.invs <- lapply(f.list$value, function (x) {
  
  reg.inv <- st_read(x)
  
  O1 <- str_pad(reg.inv$O1Region, width = 2, side = "left", pad = "0")
  O2 <- str_pad(reg.inv$O2Region, width = 2, side = "left", pad = "0")
  
  O1_O2 <- paste0(O1, "-", O2)
  
  sub.rgi <- reg.inv[O1_O2 %in% buf.o2$RGI_CODE, ]
  sub.rgi$RGI_CODE <- O1_O2[O1_O2 %in% buf.o2$RGI_CODE]
  sub.rgi$RGI_CODE <- O1_O2[O1_O2 %in% buf.o2$RGI_CODE]
  
  sub.rgi <- left_join(x = sub.rgi, 
                       y = buf.o2 %>% 
                         st_drop_geometry() %>% 
                         dplyr::select(RGI_CODE, FULL_NAME), 
                       by = "RGI_CODE")
  
  sub.rgi <- left_join(x = sub.rgi, 
                       y = thickness.data %>% 
                         dplyr::select(RGIId, Hcomp),
                       by = "RGIId") %>%
    st_drop_geometry()
  
})

# We combine all tables to one large table.

reg.invs.bind <- do.call(rbind, reg.invs) %>% 
  as_tibble()

# We save this result to disk.

saveRDS(reg.invs.bind, "reg_invs_bind.rds")
# reg.invs.bind <- readRDS("reg_invs_bind.rds")

# We add the estimated mean ice thickness to all glaciers that generated GLOFs
# between 1990 and 2023.

mag.vs.thick <- left_join(max.lake, 
                          reg.invs.bind %>%
                            dplyr::select(RGIId, Hcomp), 
                          by = c("RGI_Glacier_Id" = "RGIId")) %>%
  drop_na() %>% 
  filter(Hcomp > 0)

# We plot the mean ice thickness against the size of the GLOF.

mag.vs.thick.plot <- mag.vs.thick %>% 
  ggplot(mapping = aes(x = Hcomp, 
                       y = max_area,
                       color = Lake_type_simple)) +
  geom_point2(alpha = 0.6) +
  scale_color_manual(values = c("navy", "#ee7600")) +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(limits = c(10^0.95, 10^3),
                     labels = label_number(drop0trailing = TRUE),
                     trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x)) +
  theme_bw() +
  xlab("Average thickness of parent glacier [m]") +
  ylab("Largest lake area before the outburst [m²] per glacier") + 
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none", 
         aspect.ratio = 1) 

# We add a marginal boxplot to show the differences in ice thickness between
# glacier-dammed and moraine- and bedrock-dammed lakes.

mag.vs.thick.plot.mar <- ggMarginal(mag.vs.thick.plot,
                                    margins = "x", 
                                    size = 7, 
                                    type = "boxplot", 
                                    groupColour = T,
                                    outlier.shape = NA, 
                                    lwd = 0.5)

# We calculate the ice thickness of all glaciers that generated GLOFs from 
# glacier-dammed and supraglacial lakes.

ice.thick <- mag.vs.thick %>%
  filter(Lake_type_simple == "glacier_supraglacial") %>% 
  ungroup() %>% 
  summarise(med = median(Hcomp), 
            q25 = quantile(Hcomp, 0.25),
            q75 = quantile(Hcomp, 0.75))

# We calculate the ice thickness of all glaciers that generated GLOFs from 
# moraine- and bedrock-dammed lakes.

mb.thick <- mag.vs.thick %>%
  filter(Lake_type_simple == "moraine_bedrock") %>% 
  ungroup() %>% 
  summarise(med = median(Hcomp), 
            q25 = quantile(Hcomp, 0.25),
            q75 = quantile(Hcomp, 0.75))

# What is the average ice thickness of all glaciers in our study regions?

reg.invs.bind  %>%
  filter(!is.na(Hcomp)) %>% 
  filter(Hcomp > 0)  %>% 
  summarise(med = median(Hcomp), 
            q25 = quantile(Hcomp, 0.25),
            q75 = quantile(Hcomp, 0.75))

# How large (in terms of the percentile) were glaciers that generated ice-dam
# failures compared to all other glaciers?

v <- reg.invs.bind  %>%
  filter(!is.na(Hcomp)) %>% 
  filter(Hcomp > 0) %>% pull(Hcomp)  

ecdf(v)(ice.thick$med)

################################################################################
### GLOF size versus glacier length ############################################

# Finally, we compare the length of glaciers that burst out against the size
# of GLOFs. As above, we distinguish between glacier-dammed and moraine- and 
# bedrock-dammed lakes.

# We join the length (Lmax) of the glaciers in the RGI inventory to all lakes
# that have burst out in our study period.

mag.vs.length <- left_join(max.lake, 
                           reg.invs.bind %>%
                             dplyr::select(RGIId, Lmax), 
                           by = c("RGI_Glacier_Id" = "RGIId")) %>%
  drop_na() %>% 
  filter(Lmax > 0)

# We generate a plot of GLOF size versus glacier length.

mag.vs.length.plot <- mag.vs.length %>%
  mutate(Lmax = Lmax/ 1000) %>%
  ggplot(mapping = aes(x = Lmax, 
                       y = max_area,
                       color = Lake_type_simple)) +
  geom_point2(alpha = 0.6) +
  scale_color_manual(values = c("navy", "#ee7600"))  +
  scale_y_continuous(# limits = c(10^3, 10^9),
    trans  = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))) +
  scale_x_continuous(limits = c(0.1, 300),
                     trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = label_number(drop0trailing = TRUE)) +
  theme_bw() +
  xlab("Length of parent glacier [km]") +
  ylab("Lake area before the outburst [m²]") + 
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none", aspect.ratio = 1) 

# Add marginal boxplots on both axis to show the differences in glacier length
# and GLOF size between glacier-dammed and moraine- and bedrock-dammed lakes.

mag.vs.length.plot.mar <-  ggMarginal(mag.vs.length.plot,
                                      margins = "both", 
                                      size = 7, 
                                      type = "boxplot", 
                                      groupColour = T,
                                      outlier.shape = NA, 
                                      lwd = 0.5)

# We calculate the length of all glaciers that generated GLOFs from 
# glacier-dammed and supraglacial lakes.

ice.length <- mag.vs.length %>%
  filter(Lake_type_simple == "glacier_supraglacial") %>% 
  ungroup() %>% 
  summarise(med = median(Lmax), 
            q25 = quantile(Lmax, 0.25),
            q75 = quantile(Lmax, 0.75))

# We calculate the length of all glaciers that generated GLOFs from 
# moraine- and bedrock-dammed lakes.

mb.length <- mag.vs.length %>%
  filter(Lake_type_simple == "moraine_bedrock") %>% 
  ungroup() %>% 
  summarise(med = median(Lmax), 
            q25 = quantile(Lmax, 0.25),
            q75 = quantile(Lmax, 0.75))

# What is the length of all glaciers in our regions?

reg.invs.bind  %>%
  filter(!is.na(Lmax)) %>% 
  filter(Lmax > 0)  %>% 
  summarise(med = median(Lmax), 
            q25 = quantile(Lmax, 0.25),
            q75 = quantile(Lmax, 0.75))

v <- reg.invs.bind  %>%
  filter(!is.na(Lmax)) %>% 
  filter(Lmax > 0) %>% pull(Lmax)  

ecdf(v)(ice.length$med)

# Combine the three plots (GLOF size versus regional GLOF reporting rate, 
#                          GLOF size versus ice thickness, and
#                          GLOF size versus glacier length).

arr.plot <- ggarrange(mag.freq.mar, 
                      mag.vs.thick.plot.mar, 
                      mag.vs.length.plot.mar, 
                      ncol = 3)

# Save the plot to disk.

ggsave("Lake_area_vs_all.pdf",
       arr.plot,
       width = 180,
       units = "mm")