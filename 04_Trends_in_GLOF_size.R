################################################################################
#######        Estimate trends in GLOF size between 1990 and 2023         ######
#######                                                                   ######
#######                            by Georg Veh                           ######
#######                checked and comments added March 04, 2024          ######
#######  checked again, drainage ratio and comments added, Dec 23, 2024   ######
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
require(brms)
require(ggpubr)
require(see)

# Set YOUR working directory folder where to find all files, necessary to run 
# this script. Change the location appropriately.

setwd("D:/data/BoxUP/Work 2022/GLOFsize/")

# Useful functions

scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Read the R file containing all reported and verified GLOFs between 1990 and
# 2023. This is the output of the first R script.

rep.glof <- readRDS("reported_GLOFs.rds")

################################################################################
### Trends in size of burst glacier-dammed and supraglacial lakes ##############

# We select all ice-dammed lakes that burst out at least 5 times between 1990
# and 2023.

gs <- rep.glof %>% 
  filter(Lake_type_simple == "glacier_supraglacial") %>%
  group_by(glacier_and_lake) %>%
  filter(n() >= 5)   %>%
  ungroup() %>%
  mutate(area_scale = scale_this(log10(Lake_area_before)),
         year_scale  = scale_this(rounded_year)) %>%
  mutate(RegO2_adj_f2 = str_replace_all(RegO2_adj_f, "_", " ")) %>%
  mutate(RegO2_adj_f_n = factor(RegO2_adj_f2, 
                                levels = c(unique(grep("Alaska Range", RegO2_adj_f2, value = T)),
                                           unique(grep("W Chugach Mountains", RegO2_adj_f2, value = T)),
                                           unique(grep("Saint Elias Mountains", RegO2_adj_f2, value = T)),
                                           unique(grep("Coast Ranges", RegO2_adj_f2, value = T)),
                                           unique(grep("C and N Andes", RegO2_adj_f2, value = T)),
                                           unique(grep("Patagonia", RegO2_adj_f2, value = T)),
                                           unique(grep("Iceland", RegO2_adj_f2, value = T)),
                                           unique(grep("Scandinavia", RegO2_adj_f2, value = T)),
                                           unique(grep("Alps", RegO2_adj_f2, value = T)),
                                           unique(grep("Hissar Alay and Tien Shan", RegO2_adj_f2, value = T)),
                                           unique(grep("HK Pamir Karakoram", RegO2_adj_f2, value = T)),
                                           unique(grep("Himalayas", RegO2_adj_f2, value = T)),
                                           unique(grep("Tibet and Hengduan Shan", RegO2_adj_f2, value = T)))))

# We set weakly informed priors on the model intercept, slopes, group-level 
# standard and residual standard deviations.

bprior <- prior(student_t(3, 0, 2), class = "b") +
  prior(student_t(3, 0, 2), class = "Intercept") +
  prior(student_t(3, 0.5, 2), class = "sd") 

# Fit the model of GLOF size versus time for ice-dammed lakes.

mod.gs <- brm(area_scale ~ year_scale + (year_scale | RegO2_adj_f/glacier_and_lake),
              data = gs,
              family = student(), 
              warmup  = 1000,
              iter = 4000,
              prior = bprior,
              chains  = 4,
              cores   = 4, 
              control = list(adapt_delta = 0.95,
                             max_treedepth = 15),
              backend = "cmdstanr",
              threads = threading(3))

# Show the summary of the model. 

mod.gs
plot(mod.gs)

pp_check(mod.gs)

# Looks very promising!

################################################################################

# After fitting model, we summarise the model on different hierarchical levels:
# - the grand mean (population level);
# - the group level (regions);
# - the local level (individual ice-dammed lakes).

# Grand mean: glacier and supraglacial lakes
# Define a range of time steps for which we would like to obtain draws from the
# expected mean of the posterior distribution. For each time step, we obtain
# 1500 draws, and convert the standardised data back to the original scale.

gs.scaled.years <- (seq_range(c(1990, 2023), n = 100) - mean(gs$rounded_year)) / sd(gs$rounded_year)

pred.grid.gs.grand <- add_epred_draws(
  object = mod.gs, 
  newdata = data.frame(year_scale = gs.scaled.years),
  value = "area_scale", 
  ndraws = 1500,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(rounded_year = (year_scale * sd(gs$rounded_year)) + mean(gs$rounded_year),
         Lake_area_before = 10^((area_scale * sd(log10(gs$Lake_area_before))) + mean(log10(gs$Lake_area_before))),
         Lake_type_simple = "glacier_supraglacial") 

# Plot the trend in the grand mean of ice-dammed lakes.

plot.trend.grand.mean.gs <- pred.grid.gs.grand %>%
  ggplot(aes(x = rounded_year, y = Lake_area_before)) +
  geom_jitter2(data = gs,
               aes(x = rounded_year, y = Lake_area_before)) +
  scale_fill_manual(name = "Posterior rate", values = "navy") +
  stat_lineribbon(aes( y = Lake_area_before), .width = 0.95,
                  point_interval = mean_qi) +
  theme_bw() +
  labs(x = "Year",
       y = "Lake area before GLOF [m²]") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7))  +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

# Obtain the marginal posterior distribution of GLOF size from ice-dammed lakes 
# in 1990 and 2023.

pred.grid.gs.grand %>%
  filter((rounded_year == 1990) | (rounded_year == 2023)) %>%
  group_by(rounded_year) %>%
  mutate(Lake_area_before = Lake_area_before/ 10^6) %>%
  summarise(med = median(Lake_area_before),
            q025 = quantile(Lake_area_before, 0.025),
            q975 = quantile(Lake_area_before, 0.975),
            qminus = med - q025,
            qplus = q975 - med)

################################################################################

# Regional effects: glacier and supraglacial lakes

# We first extract posterior draws on the population level. 

fixef_samples <- fixef(mod.gs, pars = "year_scale", summary = F) %>% 
  as_tibble()

# Then, we extract random effects for the posterior samples for 
# the slope parameter ("year_scale") and the region. The random effects measure
# the deviation from the population level.

ranef_samples_regions <- ranef(mod.gs, 
                               groups = "RegO2_adj_f", 
                               pars = "year_scale",  
                               summary = F) 

ranef_samples_regions <- ranef_samples_regions$RegO2_adj_f[, , "year_scale"]

# Convert random effects to a tibble.

ranef_samples_regions <- ranef_samples_regions  %>% 
  as_tibble() 

# As the random effects measure the deviations from the grand mean, we need to 
# add them to the fixed (population-level) effects.

fixef_and_ranef <- ranef_samples_regions + fixef_samples$year_scale 

fixef_and_ranef <- fixef_and_ranef %>% 
  as_tibble() 

# We summarise the regional trend of ice-dam failures into credibly decreasing, 
# increasing, or unchanged.

region.summary.gs <- pivot_longer(fixef_and_ranef, 
                                  cols = colnames(fixef_and_ranef),
                                  names_to = "RegO2_adj_f") %>%
  group_by(RegO2_adj_f) %>% 
  summarise(quantslow = quantile(value, 0.025), 
            median    = median(value),
            quantsup  = quantile(value, 0.975)) %>%
  mutate(ribcol = case_when(quantslow < 0 & quantsup < 0 ~ "red",
                            quantslow > 0 & quantsup > 0 ~ "blue",
                            .default = "grey"))

# We would like to generate draws from the expected posterior distribution.
# Therefore we generate a sequence of time steps for each region.

pred.grid.regions.gs <-  lapply(seq_range(gs$year_scale, n = 51), 
                                function (x) { (unique(gs[,  "RegO2_adj_f"])) %>%
                                    mutate(year_scale = x) })

pred.grid.regions.gs <- do.call(rbind, pred.grid.regions.gs) %>% 
  ungroup()

# We then draw 500 times from the expected posterior predictive distribution for
# each time step and region. We then convert the standardised data back to the
# original scale.

pred.grid.regions.gs.o <- add_epred_draws(
  object = mod.gs, 
  newdata = pred.grid.regions.gs,
  value = "area_scale", 
  ndraws = 500,
  re_formula =  ~ year_scale + (year_scale | RegO2_adj_f)) %>%
  mutate(rounded_year = (year_scale * sd(gs$rounded_year)) + mean(gs$rounded_year),
         Lake_area_before = 10^((area_scale * sd(log10(gs$Lake_area_before))) + mean(log10(gs$Lake_area_before))),
         Lake_type_simple = "glacier_supraglacial") 

# We add information on the region and the continent to this tibble, and remove
# underscores from names to nicely show them in plot.

pred.grid.regions.gs.o2 <- pred.grid.regions.gs.o %>%
  left_join(x = .,
            y = gs %>% dplyr::select(RegO2_adj_f, Continent),
            by = "RegO2_adj_f",
            relationship = "many-to-many") %>%
  left_join(x = ., y = region.summary.gs, by = "RegO2_adj_f") %>%
  mutate(RegO2_adj_f = str_replace_all(RegO2_adj_f, "_", " ")) %>%
  mutate(RegO2_adj_f_n = factor(RegO2_adj_f, 
                                levels = c(unique(grep("Alaska Range", RegO2_adj_f, value = T)),
                                           unique(grep("W Chugach Mountains", RegO2_adj_f, value = T)),
                                           unique(grep("Saint Elias Mountains", RegO2_adj_f, value = T)),
                                           unique(grep("Coast Ranges", RegO2_adj_f, value = T)),
                                           unique(grep("C and N Andes", RegO2_adj_f, value = T)),
                                           unique(grep("Patagonia", RegO2_adj_f, value = T)),
                                           unique(grep("Iceland", RegO2_adj_f, value = T)),
                                           unique(grep("Scandinavia", RegO2_adj_f, value = T)),
                                           unique(grep("Alps", RegO2_adj_f, value = T)),
                                           unique(grep("Hissar Alay and Tien Shan", RegO2_adj_f, value = T)),
                                           unique(grep("HK Pamir Karakoram", RegO2_adj_f, value = T)),
                                           unique(grep("Himalayas", RegO2_adj_f, value = T)),
                                           unique(grep("Tibet and Hengduan Shan", RegO2_adj_f, value = T)))))

# We plot the posterior trend in GLOF size of ice-dammed lakes with the
# original data on top.

plot.trend.year.regions.gs <- pred.grid.regions.gs.o2 %>%
  ggplot(aes(x = rounded_year, y = Lake_area_before)) +
  geom_jitter2(data = gs,
               aes(x = rounded_year, y = Lake_area_before),
               alpha = 0.3,
               color = "navy") +
  stat_lineribbon(aes(x = rounded_year, 
                      y = Lake_area_before,
                      fill = ribcol), 
                  .width = 0.95,
                  alpha = 0.6,
                  point_interval = mean_qi) +
  scale_fill_identity() +
  facet_wrap(~RegO2_adj_f_n, ncol = 5) +
  theme_bw() +
  labs(x = "Year",
       y = "Lake area before GLOF [m²]") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         strip.background = element_blank(),
         legend.position = "none")  + 
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

# We write this plot to disk.

ggsave("regional_trends_ice_dammed_lakes.pdf",
       plot = plot.trend.year.regions.gs,
       width = 180,
       height = 90,
       units = "mm")

################################################################################

# Local effects: glacier and supraglacial lakes

# Finally, we also extract the local effects, that is the trend in GLOF size
# for all ice-dammed lakes that burst out at least 5 times. This means,
# we need to go one hierarchical level deeper.

ranef_samples_gl <- ranef(mod.gs, 
                          groups = "RegO2_adj_f:glacier_and_lake", 
                          pars = "year_scale",  
                          summary = F) 

ranef_samples_gl <- ranef_samples_gl$`RegO2_adj_f:glacier_and_lake`[ , , "year_scale"] %>%
  as_tibble()

# As above, we need to add the local effect to the regional effect and the 
# global effects.

for(x in colnames(fixef_and_ranef)) {
  
  reg_ef <- fixef_and_ranef[ , x] %>% pull()
  
  idx <- grep(x, colnames(ranef_samples_gl))
  
  ranef_samples_gl[ , idx] <- (ranef_samples_gl[ , idx] + reg_ef) %>% 
    as_tibble()
  
}

# Summarise whether the size of ice-dammed lakes before the GLOFs has increased,
# decreased, or remained unchanged during our study period.
# We also color code these trends to intuitively visualize them in figures.
# We use the 95% HDI as usual, and the the 65% interval, which is equivalent 
# to a variance of one standard deviation in a normal distribution. 

local.summary <- pivot_longer(ranef_samples_gl, 
                              cols = colnames(ranef_samples_gl), 
                              names_to = "glacier_and_lake") %>%
  group_by(glacier_and_lake) %>% 
  summarise(quantslow = quantile(value, 0.025) * sd(log10(gs$Lake_area_before)) / sd(gs$rounded_year) *10, 
            quantsup  = quantile(value, 0.975) * sd(log10(gs$Lake_area_before)) / sd(gs$rounded_year)*10,
            median    = median(value)  * sd(log10(gs$Lake_area_before)) / sd(gs$rounded_year)*10,
            q17       = quantile(value, 0.17)  * sd(log10(gs$Lake_area_before)) / sd(gs$rounded_year)*10,
            q83       = quantile(value, 0.83) * sd(log10(gs$Lake_area_before)) / sd(gs$rounded_year)*10) %>% 
  mutate(ribcol = case_when(quantslow < 0 & quantsup < 0 ~ "red",
                            quantslow > 0 & quantsup > 0 ~ "blue",
                            .default = "grey")) %>%
  mutate(Trend = case_when(ribcol == "red"  ~ "Decrease",
                           ribcol == "blue" ~ "Increase",
                           ribcol == "grey" ~ "No change"))

# Again, remove underscores from the names.

for(x in unique(gs$RegO2_adj_f)) {
  
  local.summary <- local.summary %>%
    mutate(glacier_and_lake = str_replace_all(glacier_and_lake,  paste0(x, "_"), "")) 
  
}

# We compare these trends to the thinning of ice-dammed lakes. To do so,
# we need to bring the glacier name to conventional RGI ID. 

local.summary <- local.summary %>%
  mutate(rgiid = str_replace_all(glacier_and_lake, "_", "."), 
         rgiid = str_sub(rgiid, 1, 14)) 

str_sub(local.summary$rgiid, 6, 6) <- "-"

# We then load the rates of glacier elevation change from the study: 
# Hugonnet et al. (2021): "Accelerated global glacier mass loss in the early 
# twenty-first century." Nature 592, 7856, 726-731. Data are available here: 
# https://doi.org/10.6096/13 (Go to "DOWNLOAD" -> "TIMESERIES", and download
# all Per glacier time series for the regions that cover our study regions.)

# First, list all files from Hugonnet.

dhdt <- list.files(path = "Hugonnet_2021/", 
                   pattern = "rgi60_pergla_rates.csv$", 
                   full.names = T,
                   recursive = T) 

# For each glacier, we extract only the rate of elevation change between 
# 2000 and 2020.

dhdt.all <- lapply(dhdt, function (x) {
  read.table(x, sep = ",", header = T) %>% 
    as_tibble() %>%
    filter(period == "2000-01-01_2020-01-01") %>%
    select(rgiid, dhdt, err_dhdt)})

# Bind the regional tables of elevation changes to one long table.

dhdt.bind <- do.call(rbind, dhdt.all)

# Add the elevation change statistics to the glaciers that had outbursts 
# in our study period.

summary.with.dhdt <- left_join(local.summary, dhdt.bind, by = "rgiid") %>%
  mutate(dhdt_low = dhdt - err_dhdt,
         dhdt_up = dhdt + err_dhdt)

# Plot the trend in glacier elevation change versus the trend in GLOF size. 

# First, as a scatter plot.

dhdt.vs.la.gs <- ggplot(summary.with.dhdt, 
                        aes(x = median, y = dhdt, color = ribcol)) +
  geom_point(size = 0.8) +
  geom_linerange(aes(
    ymin = dhdt_low,
    ymax = dhdt_up),
    alpha = 0.7,
    linewidth = 0.5)+
  geom_linerange(aes(xmin = q17,  
                     xmax = q83),
                 alpha = 0.7,
                 linewidth = 0.5)+
  scale_color_identity() +
  theme_bw() + 
  ylab("Annual glacier elevation change (2000-2019) [m yr-1]") +
  xlab("Decadal change in lake area in\norders of magnitude [log10(m2) dec-1]")  +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7)) +
  ylim(c(-4.1,1.1))

# Then, aggregated as boxplots into decreasing, unchanged, and increasing trends
# in GLOF size.

dhdt.vs.la.gs.boxplot <- summary.with.dhdt %>% 
  mutate(Trend_fct = factor(Trend, levels = c("Decrease", "No change", "Increase"))) %>%
  ggplot(aes(x = Trend_fct, y = dhdt, fill = ribcol)) +
  geom_boxplot() +
  scale_fill_identity() +
  theme_bw() + 
  ylab("") +
  xlab("Posterior trend in lake area")  +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7)) +
  ylim(c(-4.1,1.1))

# Combine these two plots to one.

dhdt.total <- ggarrange(dhdt.vs.la.gs, dhdt.vs.la.gs.boxplot,
                        ncol = 2,
                        labels = c("a", "b"),
                        align = "hv",
                        font.label = list(size = 8,
                                          color = "black",
                                          face = "plain"),
                        widths = c(1, 1))

# Write this figure to disk (Extended Data Figure 4).

ggsave(filename = "dhdt_vs_glacier_dammed_lake_area.pdf",
       plot     = dhdt.total,
       width = 160,
       height = 80,
       units = "mm")

################################################################################

# We also want to draw the trend in pre-GLOF lake area for each lake that had 
# repeated outbursts. 

# We derive a range of time steps for each lake, for which we want to draw
# samples from the expected posterior distribution.

conds <- gs %>% 
  group_by(RegO2_adj_f, glacier_and_lake) %>% 
  summarise(min_r = min(year_scale), 
            max_r = max(year_scale)) %>%
  ungroup()

# We obtain the posterior distribution of GLOF size for each lake.

pred.grid.local.gs <-  lapply(seq_range(gs$year_scale, n = 51), 
                              function (x) { (unique(gs[, c("glacier_and_lake", "RegO2_adj_f")])) %>% 
                                  mutate(year_scale = x) })

pred.grid.local.gs  <- do.call(rbind, pred.grid.local.gs ) %>% 
  ungroup()

# As above, gather the predictions, and convert them to the original scale.

pred.grid.local.gs.o <- add_epred_draws(
  object = mod.gs, 
  newdata = pred.grid.local.gs,
  value = "area_scale", 
  ndraws = 500,
  re_formula = area_scale ~ year_scale + (year_scale | RegO2_adj_f/glacier_and_lake)) %>%
  mutate(rounded_year = (year_scale * sd(gs$rounded_year)) + mean(gs$rounded_year),
         Lake_area_before = 10^((area_scale * sd(log10(gs$Lake_area_before))) + mean(log10(gs$Lake_area_before))),
         Lake_type_simple = "glacier_supraglacial") 

preds.sub.gs <- list()

for (i in 1:nrow(conds)) {
  
  preds.sub.gs[[i]] <- filter(pred.grid.local.gs.o, 
                              RegO2_adj_f == conds$RegO2_adj_f[i],
                              glacier_and_lake == conds$glacier_and_lake[i],
                              year_scale >= conds$min_r[i],  
                              year_scale <= conds$max_r[i])
}

preds.sub.gs <- bind_rows(preds.sub.gs) %>% 
  ungroup()

# We color the local trends in ice-dam failures by red (decreasing),
# grey (unchanged), and blue (increasing).

preds.sub.gs  <- left_join(x = preds.sub.gs,
                           y = local.summary,
                           by = "glacier_and_lake") %>%
  mutate(ribcol_fct = factor(ribcol))

# Plot the posterior rates in GLOF size for each ice-dammed lake with repeated
# failures.

plot.trend.local.gs <- preds.sub.gs  %>%
  mutate(glacier_and_lake = str_replace_all(glacier_and_lake, "_", " "),
         glacier_and_lake = str_replace_all(glacier_and_lake, "RGI60 ", "")) %>%
  ggplot(aes(x = rounded_year, y = Lake_area_before)) +
  stat_lineribbon(mapping = aes(x = rounded_year, 
                                y = Lake_area_before, 
                                fill = ribcol),
                  .width = 0.95,
                  point_interval = mean_qi) +
  scale_fill_identity() +
  geom_point2(data = gs  %>%
                mutate(glacier_and_lake = str_replace_all(glacier_and_lake, "_", " "),
                       glacier_and_lake = str_replace_all(glacier_and_lake, "RGI60 ", "")),
              mapping = aes(x = rounded_year, y = Lake_area_before),
              alpha = 0.6) +
  facet_wrap(~glacier_and_lake, scales = "free", ncol = 5) +
  theme_bw() +
  labs(x = "Year",
       y = "Lake area before GLOF [m²]") +
  theme_bw() +
  theme( axis.text   = element_text(size = 6),
         axis.text.x = element_text(size = 6),
         axis.title  = element_text(size = 6),
         strip.text  = element_text(size = 6),
         strip.background = element_blank(),
         legend.position = "none") +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

# Write this plot to disk (Extended Data Figure 3).

# ggsave("local_trends_idl.pdf",
#        plot.trend.local.gs, 
#        width = 180, 
#        height = 509, 
#        units = "mm")

################################################################################
### Trends in size of burst moraine- and bedrock-dammed lakes ##################

# We do basically the same analysis for moraine- and bedrock-dammed lakes.
# The only difference is that there is no additional level for individual lakes
# given that moraine- and bedrock-dammed lakes rarely burst twice. Put differently,
# this model is simpler, and has fewer data than that of ice-dammed lakes.

mb <- rep.glof %>% 
  filter(Lake_type_simple == "moraine_bedrock") %>%
  group_by(RegO2_adj_f) %>%
  filter(n() >= 5) %>% 
  filter(Lake != "Lago Greve") %>%
  mutate(area_scale = scale_this(log10(Lake_area_before)),
         year_scale  = scale_this(rounded_year)) %>%
  mutate(RegO2_adj_f2 = str_replace_all(RegO2_adj_f, "_", " ")) %>%
  mutate(RegO2_adj_f_n = factor(RegO2_adj_f2, 
                                levels = c(unique(grep("Alaska Range", RegO2_adj_f2, value = T)),
                                           unique(grep("W Chugach Mountains", RegO2_adj_f2, value = T)),
                                           unique(grep("Saint Elias Mountains", RegO2_adj_f2, value = T)),
                                           unique(grep("Coast Ranges", RegO2_adj_f2, value = T)),
                                           unique(grep("C and N Andes", RegO2_adj_f2, value = T)),
                                           unique(grep("Patagonia", RegO2_adj_f2, value = T)),
                                           unique(grep("Iceland", RegO2_adj_f2, value = T)),
                                           unique(grep("Scandinavia", RegO2_adj_f2, value = T)),
                                           unique(grep("Alps", RegO2_adj_f2, value = T)),
                                           unique(grep("Hissar Alay and Tien Shan", RegO2_adj_f2, value = T)),
                                           unique(grep("HK Pamir Karakoram", RegO2_adj_f2, value = T)),
                                           unique(grep("Himalayas", RegO2_adj_f2, value = T)),
                                           unique(grep("Tibet and Hengduan Shan", RegO2_adj_f2, value = T))))) %>%
  ungroup()

# Set weakly informed priors on standardised data.

bprior <- prior(student_t(3, 0, 2), class = "b") +
  prior(student_t(3, 0, 2), class = "Intercept") +
  prior(student_t(3, 0.5, 2), class = "sd") 

mod.mb <- brm(area_scale ~ year_scale + (year_scale | RegO2_adj_f),
              data = mb,
              family = student(), 
              warmup  = 1000,
              iter = 4000,
              prior = bprior,
              chains  = 4,
              cores   = 4, 
              control = list(adapt_delta = 0.97,
                             max_treedepth = 15),
              backend = "cmdstanr",
              threads = threading(3))

# Assess model performance

mod.mb

plot(mod.mb)
pp_check(mod.mb)

# Looks fine!

###############################################################################

# Grand mean: moraine- and bedrock-dammed lakes

# We define a sequence of time steps, for which we would like to have predictions.

mb.scaled.years <- (seq_range(c(1990, 2023), n = 100) - mean(mb$rounded_year)) / sd(mb$rounded_year)

# After drawing samples, we convert the data to the original scale.

pred.grid.grand.mb <- add_epred_draws(
  object = mod.mb, 
  newdata = data.frame(year_scale = mb.scaled.years),
  value = "area_scale", 
  ndraws = 2000,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(rounded_year = (year_scale * sd(mb$rounded_year)) + mean(mb$rounded_year),
         Lake_area_before = 10^((area_scale * sd(log10(mb$Lake_area_before))) + mean(log10(mb$Lake_area_before))),
         Lake_type_simple = "moraine_bedrock") 
  
# We plot the trend of the grand mean, i.e. the population level.

plot.trend.grand.mean.mb <- pred.grid.grand.mb %>%
  ggplot(aes(x = rounded_year, y = Lake_area_before)) +
  geom_point(data = mb,
             aes(x = rounded_year, y = Lake_area_before),
             shape = 16) +
  scale_fill_manual(name = "Posterior trend", values = "#52c8c8c8") +
  stat_lineribbon(aes( y = Lake_area_before), .width = 0.95,
                  point_interval = mean_qi) +
  theme_bw() +
  labs(x = "Year",
       y = "Lake area before GLOF [m²]") +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7))

###############################################################################

# Regional effects: moraine- and bedrock-dammed lakes

# We extract fixed effects posterior samples

fixef_samples <- fixef(mod.mb, pars = "year_scale", summary = F) %>% 
  as_tibble()

# Extract random effects posterior samples for the slope (parameter "year_scale")
# and the given region.

ranef_samples_regions <- ranef(mod.mb, 
                               groups = "RegO2_adj_f", 
                               pars = "year_scale",  
                               summary = F) 

ranef_samples_regions <- ranef_samples_regions$RegO2_adj_f[, , "year_scale"]

ranef_samples_regions <- ranef_samples_regions  %>% 
  as_tibble() 

fixef_and_ranef <- ranef_samples_regions + fixef_samples$year_scale 

fixef_and_ranef <- fixef_and_ranef %>% 
  as_tibble() 

# We assess which regions had positive, negative, or unchanged trends in 
# the size of moraine- and bedrock-dam failures.

pivot_longer(fixef_and_ranef, cols = colnames(fixef_and_ranef )) %>%
  group_by(name) %>% 
  summarise(quantslow = quantile(value, 0.025), 
            quantsup = quantile(value, 0.975)) %>% 
  View()

# We add a color code to better visualise these trends.

region.summary.mb <- pivot_longer(fixef_and_ranef, 
                                  cols = colnames(fixef_and_ranef),
                                  names_to = "RegO2_adj_f") %>%
  group_by(RegO2_adj_f) %>% 
  summarise(quantslow = quantile(value, 0.025), 
            median    = median(value),
            quantsup  = quantile(value, 0.975)) %>%
  mutate(ribcol = case_when(quantslow < 0 & quantsup < 0 ~ "red",
                            quantslow > 0 & quantsup > 0 ~ "blue",
                            .default = "grey"))

# These lines are equivalent to that of the regional trends in ice-dam 
# failures and therefore left uncommented.

pred.grid.regions.mb <-  lapply(mb.scaled.years, 
                                function (x) { (unique(mb[,  "RegO2_adj_f"])) %>%
                                    mutate(year_scale = x) })

pred.grid.regions.mb <- do.call(rbind, pred.grid.regions.mb) %>% 
  ungroup()

pred.grid.regions.mb.o <- add_epred_draws(
  object = mod.mb, 
  newdata = pred.grid.regions.mb,
  value = "area_scale", 
  ndraws = 500,
  re_formula =  ~  year_scale + (year_scale | RegO2_adj_f)) %>%
  mutate(rounded_year = (year_scale * sd(mb$rounded_year)) + mean(mb$rounded_year),
         Lake_area_before = 10^((area_scale * sd(log10(mb$Lake_area_before))) + mean(log10(mb$Lake_area_before))),
         Lake_type_simple = "moraine_bedrock") 

pred.grid.regions.mb.o2 <- pred.grid.regions.mb.o %>%
  left_join(x = .,
            y = mb %>% dplyr::select(RegO2_adj_f, Continent),
            by = "RegO2_adj_f",
            relationship = "many-to-many") %>%
  left_join(x = ., y = region.summary.mb, by = "RegO2_adj_f") %>%
  mutate(RegO2_adj_f = str_replace_all(RegO2_adj_f, "_", " ")) %>%
  mutate(RegO2_adj_f_n = factor(RegO2_adj_f, 
                                levels = c(unique(grep("Alaska Range", RegO2_adj_f, value = T)),
                                           unique(grep("W Chugach Mountains", RegO2_adj_f, value = T)),
                                           unique(grep("Saint Elias Mountains", RegO2_adj_f, value = T)),
                                           unique(grep("Coast Ranges", RegO2_adj_f, value = T)),
                                           unique(grep("C and N Andes", RegO2_adj_f, value = T)),
                                           unique(grep("Patagonia", RegO2_adj_f, value = T)),
                                           unique(grep("Iceland", RegO2_adj_f, value = T)),
                                           unique(grep("Scandinavia", RegO2_adj_f, value = T)),
                                           unique(grep("Alps", RegO2_adj_f, value = T)),
                                           unique(grep("Hissar Alay and Tien Shan", RegO2_adj_f, value = T)),
                                           unique(grep("HK Pamir Karakoram", RegO2_adj_f, value = T)),
                                           unique(grep("Himalayas", RegO2_adj_f, value = T)),
                                           unique(grep("Tibet and Hengduan Shan", RegO2_adj_f, value = T)))))

# Plot the regional trends in moraine- and bedrock-dam failures.

plot.trend.year.regions.mb <- pred.grid.regions.mb.o2 %>%
  ungroup() %>%
  ggplot(aes(x = rounded_year, y = Lake_area_before)) +
  geom_jitter2(data = mb,
               aes(x = rounded_year, y = Lake_area_before),
               alpha = 0.3,
               color = "#ee7600") +
  stat_lineribbon(aes(x = rounded_year, 
                      y = Lake_area_before,
                      fill = ribcol), 
                  .width = 0.95,
                  alpha = 0.6,
                  point_interval = mean_qi) +
  scale_fill_identity() +
  facet_wrap(~RegO2_adj_f_n, ncol = 5) +
  theme_bw() +
  labs(x = "Year",
       y = "Lake area before GLOF [m²]") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         strip.background = element_blank(),
         legend.position = "none")  + 
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)))

# We combine the regional trends of ice-dammed and moraine- and bedrock-
# dammed lakes.

regional.trends.gs.mb <- ggarrange(plot.trend.year.regions.gs,
                                   plot.trend.year.regions.mb,
                                   ncol = 1,
                                   nrow = 2,
                                   labels = c("a", "b"),
                                   align = "hv",
                                   font.label = list(size = 8,
                                                     color = "black",
                                                     face = "plain"))

# We write this plot to disk. 

ggsave(filename = "regional_trends_gs_mb.pdf",
       regional.trends.gs.mb,
       width = 180,
       height = 140,
       units = "mm")

################################################################################
### Global trend of ice-dam and moraine- and bedrock-dam failures ##############

# We combine the grand mean draws for both lake types.

pred.grids.all.grand <- rbind(pred.grid.gs.grand, 
                              pred.grid.grand.mb) %>%
  mutate(RegO2_adj_f = "All regions")

# And we plot the two trends in two different colors.

p <- ggplot(data = rep.glof,
            aes(x = rounded_year, y = Lake_area_before, color = Lake_type_simple)) +
  geom_jitter2(alpha = 0.3) +
  stat_lineribbon(data = pred.grids.all.grand,
                  aes( y = Lake_area_before,
                       fill = Lake_type_simple), 
                  .width = 0.95,
                  alpha = 0.5,
                  point_interval = mean_qi) +
  theme_bw() +
  scale_fill_discrete(type = c("navy", "#ee7600")) +
  scale_color_discrete(type = c("navy", "#ee7600")) +
  labs(x = "Year",
       y = "Lake area before GLOF [m²]") +
  scale_y_continuous(trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")

# We add a marginal histogram to emphasize the difference in lake size before
# the outburst for ice-dammed lakes in contrast to moraine- and bedrock-dammed 
# lakes.

p.mar <- ggMarginal(p, 
                    groupColour = TRUE, 
                    groupFill = TRUE, 
                    type = "boxplot", 
                    margins = "y",
                    outlier.shape = NA)

# We save this plot to disk (Figure 3A).

ggsave("Global_ice_moraine.pdf",
       p.mar, 
       width = 70,
       height = 80,
       units = "mm")

# We extract the global and regional sample sizes of GLOFs per dam type.
# This information will enter the next figure.

samp.size <- rep.glof %>% 
  group_by(RegO2_adj_f, Lake_type_simple) %>% 
  summarise(n = n()) %>% 
  spread(Lake_type_simple, n, fill = 0) %>%
  rename(n_gs = "glacier_supraglacial",
         n_mb = "moraine_bedrock") %>%
  ungroup() %>%
  add_row(RegO2_adj_f = "All regions",
          n_gs = nrow(rep.glof %>% dplyr::filter(Lake_type_simple == "glacier_supraglacial")),
          n_mb = nrow(rep.glof %>% dplyr::filter(Lake_type_simple == "moraine_bedrock")))

# We plot the marginal posterior distribution of GLOF size for every region 
# in 1990 and 2023. We distinguish between the two major dam types.
# We add the sample size to the name of the region, and include the acronym of
# that region.

regional.trends <- rbind(pred.grid.regions.gs.o %>% ungroup(),
                         pred.grid.regions.mb.o %>% ungroup(),
                         pred.grids.all.grand) %>%
  filter((rounded_year == 1990) | (rounded_year == 2023)) %>%
  left_join(x = ., y = samp.size, by = "RegO2_adj_f") %>%
  mutate(RegO2_adj_f = str_replace_all(RegO2_adj_f, "_", " ")) %>%
  mutate(RegO2_adj_f = str_replace_all(RegO2_adj_f, "Alaska Range", "Alaska Range (ALR)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "W Chugach Mountains", "W Chugach Mountains (WCM)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Saint Elias Mountains", "Saint Elias Mountains (SEM)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Coast Ranges", "Coast Ranges (COR)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "C and N Andes", "C and N Andes (CNA)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Patagonia", "Patagonia (PAT)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Iceland", "Iceland (ICL)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Scandinavia", "Scandinavia (SCA)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Alps", "Alps (ALP)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Hissar Alay and Tien Shan", "Hissar Alay and Tien Shan (HATS)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "HK Pamir Karakoram", "HK Pamir Karakoram (HKPK)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Tibet and Hengduan Shan", "Tibet and Hengduan Shan (THS)"),
         RegO2_adj_f = str_replace_all(RegO2_adj_f, "Himalayas", "Himalayas (HIM)")) %>%
  mutate(RegO2_adj_f = paste0(RegO2_adj_f, " (", n_mb, " | ", n_gs, ")")) %>%
  mutate(RegO2_adj_f_n = factor(RegO2_adj_f, 
                                levels = c(unique(grep("All regions", RegO2_adj_f, value = T)),
                                           unique(grep("Alaska Range", RegO2_adj_f, value = T)),
                                           unique(grep("W Chugach Mountains", RegO2_adj_f, value = T)),
                                           unique(grep("Saint Elias Mountains", RegO2_adj_f, value = T)),
                                           unique(grep("Coast Ranges", RegO2_adj_f, value = T)),
                                           unique(grep("C and N Andes", RegO2_adj_f, value = T)),
                                           unique(grep("Patagonia", RegO2_adj_f, value = T)),
                                           unique(grep("Iceland", RegO2_adj_f, value = T)),
                                           unique(grep("Scandinavia", RegO2_adj_f, value = T)),
                                           unique(grep("Alps", RegO2_adj_f, value = T)),
                                           unique(grep("Hissar Alay and Tien Shan", RegO2_adj_f, value = T)),
                                           unique(grep("HK Pamir Karakoram", RegO2_adj_f, value = T)),
                                           unique(grep("Tibet and Hengduan Shan", RegO2_adj_f, value = T)),
                                           unique(grep("Himalayas", RegO2_adj_f, value = T))))) %>%
  ggplot(aes(y = as.character(rounded_year),
             x = Lake_area_before,
             group = Lake_type_simple)) +
  facet_grid(RegO2_adj_f_n~.,  switch = "y") +
  stat_slab( mapping = aes(fill = Lake_type_simple)) + 
  scale_x_continuous(limits =  c(10^3.5, 10^7),
                     trans  = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x))) +
  scale_fill_discrete(type = c("navy", "#ee7600")) +
  theme_bw() +
  theme( strip.background = element_blank(),
         axis.text    = element_text(size = 7),
         axis.text.x  = element_text(size = 7),
         axis.title   = element_text(size = 7),
         strip.text.y = element_text(size = 7, 
                                     hjust = 1),
         legend.title = element_text(size = 7),
         legend.text  = element_text(size = 7),
         legend.position = "bottom",
         strip.text.y.left = element_text(angle = 0),
         panel.spacing = unit(0, "lines"),
         panel.border = element_blank()) +
  guides(colour = guide_legend(ncol = 2)) +
  xlab("") +
  ylab("")

# We save this plot to disk (Figure 3C).

ggsave("Regional_trends.pdf", 
       regional.trends, 
       width = 95, 
       height = 170, 
       units = "mm")

################################################################################
### Drainage ratio #############################################################

# How many moraine- and bedrock-dammed lakes remain a hazard potential because 
# they have not (yet) drained completely?

(rep.glof %>% 
  filter(!is.na(Lake_area_after),
         Lake_area_after > 0, 
         Lake_type_simple == "moraine_bedrock") %>%
   nrow()) / 
  (rep.glof %>% 
     filter(!is.na(Lake_area_after), 
            Lake_type_simple == "moraine_bedrock") %>% 
     nrow()) 


rep.glof %>%
  mutate(drainage_ratio = ((Lake_area_before - Lake_area_after)/Lake_area_before)*100) %>%
  group_by(Lake_type_simple) %>% 
  filter(glacier_and_lake != "NA_unknown",
         drainage_ratio >= 0) %>%
  summarise(med_drain = median(drainage_ratio, na.rm = T),
            p25 = quantile(drainage_ratio, 0.25, na.rm = T),
            p75 = quantile(drainage_ratio, 0.75, na.rm = T))


# We select all moraine- and bedrock-dammed lakes and calculate their drainage 
# ratio, that is the percentage decrease in lake area (lake area before 
# minus lake area after) due to the outburst.
# To calculate a temporal trend in drainage ratio with time, we scale the drainage
# ratio and the year of the GLOF.

mb.drainage.ratio <- rep.glof %>%
  mutate(drainage_ratio = ((Lake_area_before - Lake_area_after)/Lake_area_before)*100) %>%
  filter(glacier_and_lake != "NA_unknown",
         drainage_ratio >= 0,
         Lake_type_simple == "moraine_bedrock",
         Lake_area_before < 2 * 10^7) %>%
  mutate(drainage_ratio_scale = scale_this(drainage_ratio),
         Lake_area_before_scale = scale_this(Lake_area_before),
         la_diff_scale = scale_this(la_diff),
         year_scale = scale_this(rounded_year),
         area_scale = scale_this(Lake_area_before))

bprior <- prior(student_t(3, 0, 2), class = "b") +
  prior(student_t(3, 0, 2), class = "Intercept")

################

# We first derive a model that calculates the temporal trend in the difference in lake area
# before and after GLOF (in km²). 

mod.la.diff <- brm(la_diff_scale ~ year_scale ,
              data = mb.drainage.ratio,
              family = student(), 
              warmup  = 1000,
              iter = 4000,
              prior = bprior,
              chains  = 4,
              cores   = 4, 
              control = list(adapt_delta = 0.95,
                             max_treedepth = 15),
              backend = "cmdstanr",
              threads = threading(3))

# Show the summary of the model. 

mod.la.diff 
plot(mod.la.diff )

pp_check(mod.la.diff)

pred.grid.mb.diff.la.grand <- add_epred_draws(
  object = mod.la.diff , 
  newdata = data.frame(year_scale = seq_range(mb.drainage.ratio$year_scale, n = 100)),
  value = "la_diff_scale", 
  ndraws = 1500,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(rounded_year = (year_scale * sd(mb.drainage.ratio$rounded_year)) + mean(mb.drainage.ratio$rounded_year),
         la_diff = (la_diff_scale * sd(mb.drainage.ratio$la_diff)) + mean(mb.drainage.ratio$la_diff)) 

# Plot the trend in the decrease in lake area with time.

plot.trend.mb.diff.la.grand <- pred.grid.mb.diff.la.grand  %>%
  ggplot(aes(x = rounded_year, y = la_diff)) +
  scale_fill_manual(name = "Posterior rate", values = "#ee7600") +
  stat_lineribbon(aes( y = la_diff), .width = 0.95,
                  point_interval = mean_qi) +
  geom_point(data = mb.drainage.ratio %>% filter(la_diff < 0.5),
             aes(x = rounded_year, y = la_diff)) +
  theme_bw() +
  labs(x = "Year",
       y = "Difference in lake area before and after the GLOF [km²]") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")  

#################

# We then fit a model that calculates the temporal trend in the drainage ratio. 

mod.drainage.ratio <- brm(drainage_ratio_scale ~ area_scale ,
                   data = mb.drainage.ratio,
                   family = student(), 
                   warmup  = 1000,
                   iter = 4000,
                   prior = bprior,
                   chains  = 4,
                   cores   = 4, 
                   control = list(adapt_delta = 0.95,
                                  max_treedepth = 15),
                   backend = "cmdstanr",
                   threads = threading(3))

# Show the summary of the model. 

mod.drainage.ratio 
plot(mod.drainage.ratio )

pp_check(mod.drainage.ratio)

pred.grid.mb.drainage.ratio.grand <- add_epred_draws(
  object = mod.drainage.ratio , 
  newdata = data.frame(area_scale = seq_range(mb.drainage.ratio$area_scale, n = 100)),
  value = "drainage_ratio_scale", 
  ndraws = 1500,
  re_formula = NA) %>%
  ungroup() %>%
  mutate(Lake_area_before = (area_scale * sd(mb.drainage.ratio$Lake_area_before)) + mean(mb.drainage.ratio$Lake_area_before),
         drainage_ratio = (drainage_ratio_scale * sd(mb.drainage.ratio$drainage_ratio)) + mean(mb.drainage.ratio$drainage_ratio)) 

# Plot the trend in the drainage ratio.

plot.trend.mb.drainage.ratio.grand <- pred.grid.mb.drainage.ratio.grand  %>%
  ggplot(aes(x = Lake_area_before, y = drainage_ratio)) +
  scale_fill_manual(name = "Posterior rate", values = "#ee7600") +
  stat_lineribbon(aes( y = drainage_ratio), .width = 0.95,
                  point_interval = mean_qi) +
  geom_point(data = mb.drainage.ratio,
             aes(x = Lake_area_before, y = drainage_ratio)) +
  theme_bw() +
  labs(x = "Lake area before the GLOF",
       y = "Drainage ratio [% drained lake area]") +
  theme( axis.text   = element_text(size = 7),
         axis.text.x = element_text(size = 7),
         axis.title  = element_text(size = 7),
         strip.text  = element_text(size = 7),
         legend.position = "none")  + 
  xlim(c(3000, 1000000))

# Combine both plots, i.e. the decrease in lake area and the drainage ratio.

drainage.ratio.plot <- ggarrange(plot.trend.mb.diff.la.grand, 
          plot.trend.mb.drainage.ratio.grand ,
          ncol = 2,
          labels = c("a", "b"),
          align = "hv",
          font.label = list(size = 10,
                            color = "black",
                            face = "plain"),
          widths = c(1, 1))

ggsave(filename = "drainage_ratio_plot.pdf",
       plot     = drainage.ratio.plot,
       width = 160,
       height = 80,
       units = "mm")
