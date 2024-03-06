# GLOFsize
# Code for *"Increasingly smaller outbursts despite globally growing glacier lakes"*

## Overview

**This repository contains seven scripts to estimate trends in GLOF size between 1990 and 2023. In addition, we investigate the global and regional change in glacier lake area and focus on controls that limit increases in GLOF size**

- [01_GLOF_preprocessing.R](#01_glof_preprocessingr)
- [02_generate_glacier_buffers.R](#02_generate_glacier_buffersr)
- [03_GLOF_rate_ice_thickness_and_length.R](#03_glof_rate_ice_thickness_and_lengthr)
- [04_Trends_in_GLOF_size.R](#04_trends_in_glof_sizer)
- [05_Rates_of_lake_growth.R](#05_rates_of_lake_growthr)
- [06_Burst_lakes_and_their_neighbors.R](#06_burst_lakes_and_their_neighborsr)
- [07_Limits_to_increasing_GLOF_sizes.R](#07_limits_to_increasing_glof_sizesr)

The codes are written in the statistical programming language **R** (https://www.r-project.org/), Version 4.2.2, and called within
the Graphical User Interface **RStudio** (https://posit.co/downloads/) under a Microsoft Windows 10 operating system. 
Please install both **R and RStudio** on your machine to successfully run the codes and produce figures and R data objects.

The R codes depend on a number of packages, listed at the beginning of all scripts. Please install those packages before running the scripts. 
The comments within the scripts provide further details on model dependencies and usage of functions. 

Each script will call several input data objects, which are available via ***Zenodo***.  
We also use data on glacier outlines, glacier thickness, glacier elevation changes, and previously published lake outlines. Please download the data from the web sources provided in the scripts.  
Please put all input files into the same folder, and change the working directory (which is set at beginning of each the script) according to your folder structure. The scripts can be executed one after the other, with the user generating output that is used as input for the next script.
The scripts (and parts thereof) can also be run independent of each other using the input files (in most cases *.RDS* files) from Zenodo.
Each script will produce output in form of a figure (displayed in the associate manuscript and Extended Data figures) or R-objects.

## Scripts

### 01_GLOF_preprocessing.R

**Script to preprocess a raw OpenOffice table of reported glacier lake outburst floods.**

*Mandatory input data*: 
- "glofdatabase_2023_12_06.ods" (table with all reported GLOFs. Compiliation as of December 06, 2023)
- Randolph Glacier Inventory V6.0 (https://www.glims.org/RGI/rgi60_dl.html)

*Outputs*: 
- "all_glofs_tibble.RDS" (R-object of all reported GLOFs in the global GLOF database; data are not trimmed to the period 1990-2023)
- "la_sf.RDS" (R-object containing all GLOFs with mapped lake areas before the outburst in the period 1990-2023)
- "Qp_V0_plot.pdf" (Four-panel scatterplot of peak discharge and flood volumes versus lake area before the GLOF and GLOF-related losses in lake area) 

---

### 02_generate_glacier_buffers.R

**Script to generate a 5-km buffer around glaciers in the RGI V6.0.**

*Mandatory input data*: 
- The complete RGI V6.0 in ESRI shapefile format.

*Outputs*: 
- "rgiO2_dissolved_outlines.shp" (Merged outlines of the RGI O2 according to the 13 study regions in ESRI shapefile format) 
- "*RegionXYZ*_buffer.shp" (Separate buffers around glaciers for the 13 study regions in ESRI shapefile format)
- "dissolved_buffer.shp" (All buffers around glaciers as one ESRI shapefile)
- "glacier_buffers_split_by_O2_no_fid_correct_FULLNAME_2.gpkg" (A slightly manually corrected version of the regional glacier buffers with overlapping buffers removed in geopackage format)

---

### 03_GLOF_rate_ice_thickness_and_length.R

**Script to calculate the regional GLOF rate, local glacier thicknesses and glacier length. Those diagnostics are compared to the size of burst lakes.**

*Mandatory input data*: 
- "la_sf.RDS" (R-object containing all GLOFs with mapped lake areas before the outburst in the period 1990-2023)
- "RGI-wide_composites_stats_GV.txt" (A text table storing the mean thickness ("Hcomp") of each glacier glacier according to its RGIId. Data provided by Daniel Farinotti.)
- "glacier_buffers_split_by_O2_no_fid_correct_FULLNAME_2.gpkg" (A slightly manually corrected version of the regional glacier buffers with overlapping buffers removed in geopackage format)

*Outputs*: 
- "reported_GLOFs.rds" (R-object containing a table of reported GLOFs in the period 1990-2023 with machine readable names of glaciers and lakes)
- "reported_GLOFs_with_geometry.rds" (R-object containing a *simple features* (sf) object of reported GLOFs in the period 1990-2023 with machine readable names of glaciers and lakes)
- "reg_invs_bind.rds" (A table of all glaciers in each region that had an estimate of ice thickness)
- "Lake_area_vs_all.pdf" (A three panel plot of GLOF rate, local glacier thickness, and glacier length versus GLOF size)
  
---

### 04_Trends_in_GLOF_size.R

**Script to estimate trends in GLOF size between 1990 and 2023.**

*Mandatory input data*: 
- "reported_GLOFs.rds" (R-object containing a table of reported GLOFs in the period 1990-2023 with machine readable names of glaciers and lakes)
- "*RegionXX*_rgi60_pergla_rates.csv" (Regional rate of glacier elevation change between 2000 and 2020, data are from Hugonnet et al., 2021, https://doi.org/10.6096/13)

*Outputs*: 
- "regional_trends_ice_dammed_lakes.pdf" (Figure showing the change in GLOF size from ice-dammed dammed lakes between 1990 and 2023 for each region)
- "dhdt_vs_glacier_dammed_lake_area.pdf" (A two-panel figure showing the relationship between the trend in GLOF size and glacier elevation change)
- "local_trends_idl.pdf" (A multi-panel figure showing the trend in GLOF size for all glacier-dammed lakes between 1990 and 2023)
- "regional_trends_gs_mb.pdf" (A two-panel figure, each consisting of regional panels showing the trend in ice-dam and moraine-and bedrock-dam failures)
- "Global_ice_moraine.pdf" (Figure showing the grand mean of trends in ice-dammed and supraglacial lakes, and moraine- and bedrock-dammed lakes between 1990 and 2023)
- "Regional_trends.pdf" (Figure showing the marginal posterior distribution of GLOF size for every region in 1990 and 2023, distinguished by the two major dam types)

---

### 05_Rates_of_lake_growth.R

**Script to calculate regional rates of glacier lake growth.**

*Mandatory input data*: 
- "dissolved_buffer.shp" (All buffers around glaciers as one ESRI shapefile)
- Outlines of glacier lakes from previous lake inventories (We do not share these data because they might be subject to different licenses. Please contact the authors)
- "Glacier_lakes_global.ods" (OpenOffice spreadsheet of previously prublished glacier lake inventories, including reference to the underlying study, year, and satellite image used to map glacier lakes)

*Outputs*: 
- "UTMArea_XXX" (Glacier lake outlines split to the extent of the 13 study region. Individual lake areas are calculated in the local UTM projection)
- "lakes_per_region.rds" (R-object of rates of lake growth in each region)
- "lake_area_change.pdf" (Figure showing the rate of lake area change in each region)

---

### 06_Burst_lakes_and_their_neighbors.R

**Script to calculate distances between burst lakes and their neighbours.**

*Mandatory input data*: 
- "glacier_buffers_split_by_O2_no_fid_correct_FULLNAME_2.gpkg" (A slightly manually corrected version of the regional glacier buffers with overlapping buffers removed in geopackage format)
- "reported_GLOFs_with_geometry.rds" (R-object containing a *simple features* (sf) object of reported GLOFs in the period 1990-2023 with machine readable names of glaciers and lakes)
- Outlines of glacier lakes from previous lake inventories (We do not share these data because they might be subject to different licenses. Please contact the authors)

*Output*: 

- "distances_plot.pdf" (Figure showing the distances to the next intact glacier lake per region) 

---

### 07_Limits_to_increasing_GLOF_sizes.R

**Script to extract all lakes >1km² from previous lake inventories and assess limits to increasing GLOF sizes.**

*Mandatory input data*: 
- Outlines of glacier lakes from previous lake inventories (We do not share these data because they might be subject to different licenses. Please contact the authors)
- "lakes_gt_1km.gpkg" (All lakes >1km² as of 2015 or later, labelled with descriptors of their geometries)

*Output*: 

- "lakes_gt_1km.gpkg" (All lakes >1km² as of 2015 or later, need to be manually labelled with their geometric properties)
- "glof_predictors_plot.pdf" (Figure showing the geometric characteristics of lakes >1 km²) 

---

## Input data

Please visit the repository on Zenodo to obtain the input files.


## References

Georg Veh, Björn G. Wang, Anika Zirzow, Christoph Schmidt, Natalie Lützow, Frederic Steppat, Guoqing Zhang, Kristin Vogel, Marten Geertsema, Oliver Korup, and John J. Clague: *Increasingly smaller outbursts despite globally growing glacier lakes*. Submitted


## See also

[The Glacier Lake Outburst Flood Database V4.0](http://glofs.geoecology.uni-potsdam.de)

## Contact

**Georg Veh**  
Postdoctoral researcher in the research group on natural hazards  
Institute of Environmental Sciences and Geography  
University of Potsdam  
georg.veh@uni-potsdam.de  
https://www.uni-potsdam.de/de/umwelt/forschung/ag-naturgefahren.html
