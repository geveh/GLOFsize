# GLOFsize
# Code for *"Increasingly smaller outbursts despite globally growing glacier lakes"*

## Overview

**This repository contains seven scripts to estimate trends in the volume (*V*<sub>0</sub>), peak discharge (*Q*<sub>p</sub>), timing (day of year *doy*), and source elevation (*Z*) of ice-dam failures on regional and local (i.e. lake-) level. In addition, we investigate the consequences of melting glacier dams on the magnitudes of GLOFs.**

- [01_GLOF_preprocessing.R](#01_glof_preprocessingr)
- [02_generate_glacier_buffers.R](#02_generate_glacier_buffersr)
- [03_GLOF_rate_ice_thickness_and_length.R](#03_glof_rate_ice_thickness_and_lengthr)
- [04_Trends_in_GLOF_size.R](#04_trends_in_glof_sizer)
- [05_Rates_of_lake_growth.R](#05_rates_of_lake_growthr)
- [06_Burst_lakes_and_their_neighbors.R](#06_burst_lakes_and_their_neighborsr)
- [07_Limits_to_increasing_GLOF_sizes.r](#07_limits_to_increasing_glof_sizesr)

The codes are written in the statistical programming language **R** (https://www.r-project.org/), Version 4.2.2, and called within
the Graphical User Interface **RStudio** (https://posit.co/downloads/) under a Microsoft Windows 10 operating system. 
Please install both **R and RStudio** on your machine to successfully run the codes and produce figures and R data objects.

The R codes depend on a number of packages, listed at the beginning of all scripts. Please install those packages before running the scripts. 
The comments within the scripts provide further details on model dependencies and usage of functions. 

Each script will call one or more input data object(s), which are available via ***Zenodo***.  
We also use data on glacier outlines and previously published lake outlines. Please download the data from the web sources provided in the scripts.  
Please put all input files into the same folder, and change the working directory (which is set at beginning of each the script) according to your folder structure. The scripts can be executed one after the other, with the user generating output that is used as input for the next script.
The scripts (and parts thereof) can also be run independent of each other using the input files (in most cases *.RDS* files) from Zenodo.
Each script will produce output in form of a figure (displayed in the associate manuscript and Extended Data figures) or R-objects.

## Scripts

### 01_GLOF_preprocessing.R

**Script to preprocess a raw OpenOffice table of reported glacier lake outburst floods.**

*Mandatory input data*: 
- "glofdatabase_2023_12_06.ods" (table with all reported GLOFs. Compiliation as of December 06, 2023)
- Randolph Glacier Inventory V6.0 (https://www.glims.org/RGI/rgi60_dl.html)


*Main outputs*: 
- "all_glofs_tibble.RDS" (R-object of all reported GLOFs in the global GLOF database; data are not trimmed to the period 1990-2023)
- "la_sf.RDS" (All GLOFs with mapped lake areas before the outburst in the period 1990-2023)
- "Qp_V0_plot.pdf" (Four-panel scatterplot of peak discharge and flood volumes versus lake area before the GLOF and GLOF-related losses in lake area) 

---

### 02_generate_glacier_buffers.R

**Script to generate a 5-km buffer around glaciers in the RGI V6.0.**

*Mandatory input data*: 
- The complete RGI V6.0 in ESRI shapefile format.

*Main outputs*: 
- "rgiO2_dissolved_outlines.shp" (Merged outlines of the RGI O2 according to the 13 study regions) 
- "*RegionXYZ*_buffer.shp") (Buffers around glaciers for the 13 study regions in ESRI shapefile format)
- "glacier_buffers_split_by_O2_no_fid_correct_FULLNAME_2.gpkg" (A slightly manually corrected version of the regional glacier buffers with overlapping buffers removed)

---

### 03_GLOF_rate_ice_thickness_and_length.R

**Script to calculate the regional GLOF rate, local glacier thicknesses and glacier length. Those diagnostics are compared to the size of burst lakes.**

*Mandatory input data*: 
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)

*Output*: 
- "doy_trends_per_region.RDS" (R-object with regression models of *doy* versus time for all dated GLOFs in the six regions)
- "doy_change.pdf" (Plot of the temporal trends in *doy* for each region, including the posterior differences in *doy* between 2021 and 1900)
- "doy_trends_per_glacier.RDS"  (R-object with regression models of *doy* versus time for lakes with repeat GLOFs)
- "doy_local.pdf" (Plot of local changes in *doy* versus time)
- "post_trend_doy_per_lake.pdf"  (Plot of local  posterior differences in *doy* for each lake)

---

### 04_glacier_volumes_and_ice_loss.R

**Script to obtain the total volumes of glaciers and their volume loss between 2000 and 2019 in 100-m elevation bins.**

*Mandatory input data (Data sources from external repositories are provided in the script)*: 
- Folder "Region_extents" (Contains the ESRI shapefile *Extent_pol.shp* to display the extent of the study regions)
- Glacier outlines from the Randolph Glacier Inventory (RGI)
- Glacier surface DEMs from Farinotti et al. (2019)
- Glacier volume DEMs from Farinotti et al. (2019)
- Glacier elevation change data from Hugonnet et al. (2021)

*Output*: 
- "Regional_glacier_and_melt_volumes.rds" (R-object containing the total volume of glacier volume and volume change between 2000 and 2019 in 100-m elevation bins)

---

### 05_trends_in_Z.R

**Script to estimate regional trends in the source elevation (*Z*) of ice-dammed failures.**

*Mandatory input data*: 
- Digital Elevation models from ALOS World 3D - 30m (AW3D30)
- Files from the "GDL_database" (We created a merged lake inventory in ESRI shapefile format from regional lake databases. This lake database is available upon request)   
- "Regional_glacier_and_melt_volumes.rds" (R-object containing the total volume of glacier volume and volume change between 2000 and 2019 in 100-m elevation bins)

*Major outputs*: 
- "gdl_database_centroid.RDS" (R-object of glacier lake centroids in the six study regions)
- "glofs_ice_with_z.RDS" (R-object of first reported GLOF from a given lake and its elevation)
- "Z_trends_per_region.RDS" (R-object with a hierarchical regression models of *Z* versus time for dated GLOFs in the six regions between 1900 and 2021)
- "elev_trend.pdf" (Plot of the change in GLOF source elevation for six regions between 1900 and 2021, including the posterior regression slope)
- "Lake_GLOF_elevation.pdf" / "Lake_GLOF_elevation.png" (Plot of the elevation distribution of historic burst ice-dammed lakes and present-day ice-dammed lakes for six regions between 1900 and 2021)

---

### 06_magnitudes_vs_elev_change.R

**Script to estimate local trends of  V<sub>0</sub> and  Q<sub>p</sub> with elevation change of the glacier dam.**

*Mandatory input data*: 

- Folder "dh_pergla_cut" (Tables of cumulative elevation change (in m) for glaciers with repeat GLOFs between 2000 and 2019)
- "all_glofs_tibble.RDS" (R-object with a preprocessed table of all reported GLOFs)
- "all_glofs_V0_tibble.RDS" (Table of lakes with repeat GLOFs and reported V<sub>0</sub>)
- "all_glofs_qp_tibble.RDS" (Table of lakes with repeat GLOFs and reported Q<sub>p</sub>)
- Folder "Region_extents" (Contains the ESRI shapefile *Extent_pol.shp* to display the extent of the study regions)

*Output*: 

- "local_Qp_vs_dhdt_model.RDS" (R-Object containing a hierarchical model of local changes in Q<sub>p</sub> versus glacier elevation change)
- "local_V0_vs_dhdt_model.RDS" (R-Object containing a hierarchical model of local changes in V<sub>0</sub> versus glacier elevation change)
- "map_and_trends.pdf" (Map of lakes with repeat GLOFs between 2000 and 2019; local trends of V<sub>0</sub> and Q<sub>p</sub> with cumulative changes in glacier dam elevation)
- "dam_thinning_rats.shp" (ESRI shapefile showing mean annual elevation change of glacier dams with repeat outbursts between 2000 and 2019)
- "elev_change_per_glacier.pdf" / "elev_change_per_glacier.png" (Plot of cumulative elevation change for each glacier that produced repeated GLOFs between 2000 and 2019)

---


## Input data

Please visit the repository on Zenodo to obtain the input files.


## References

Georg Veh1, Björn G. Wang, Anika Zirzow, Christoph Schmidt, Natalie Lützow, Frederic Steppat, Guoqing Zhang, Kristin Vogel, Marten Geertsema, Oliver Korup, and John J. Clague: *Increasingly smaller outbursts despite globally growing glacier lakes*. Submitted


## See also

[The Glacier Lake Outburst Flood Database V4.0](http://glofs.geoecology.uni-potsdam.de)

## Contact

**Georg Veh**  
Postdoctoral researcher in the research group on natural hazards  
Institute of Environmental Sciences and Geography  
University of Potsdam  
georg.veh@uni-potsdam.de  
https://www.uni-potsdam.de/de/umwelt/forschung/ag-naturgefahren.html
