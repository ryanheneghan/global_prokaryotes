**Scripts and data to generate statistical models and figures for "The global distribution and climate resilience of marine heterotrophic prokaryotes".**

*Note: the data folder only contains data required to generate figures and statistical models. This data, alongside a 185MB netcdf of prokaryote abundance, cell carbon and total biomass (as well as environmental variables from World Ocean Atlas 2018 and satellite used as inputs into prokaryote statistical models) is available from the Zenodo repository linked to this paper ().*

*Also, the words "bacteria" and "prokaryote" (and abbreviations of these words) are used interchangeably throughout these scripts and data files.*


**/data summary:**\
raw_observations/prokaryote_abundance.csv : all observations of prokaryotic abundance used in study, with environmental variables appended (as described in Methods of paper)\
raw_observations/prokaryote_cellcarbon.csv : all observations of prokaryotic cell carbon used in study, from Malaspina 2010 circumnavigation. Provided by co-authors Josep Gasol and Xose Anxelu G. Moran\
raw_observations/prokaryote_gge.csv : all observations of prokaryotic gross growth efficiency, provided by co-author Gerhard Herndl and Carol Robinson\
raw_observations/prokaryote_sgr.csv : all observations of prokaryotic specific-growth rates, provided by co-author Josep Gasol

bacteria_epipelagic.csv : mean (across top 200m) temperature, nitrate, aou, prokaryotic biomass, cell carbon, respiration and abundance\
bacteria_mesopelagic.csv : mean (across 200-1000m depth) temperature, nitrate, aou, prokaryotic biomass, cell carbon and abundance\
bacteria_bathypelagic.csv : mean (across >1000m depth) temperature, nitrate, aou, prokaryotic biomass, cell carbon and abundance\
bathy_data.csv : mean bathymetry of 1 degree global ocean grid cells (from GEBCO)\
epipelagic_all_heterotroph_biomass.csv : total epipelagic microzoo, macrozoo, mesozoo, rhizaria, fish and prokaryote biomass by grid cell, and fraction of total heterotroph biomass that prokaryotes. \
malaspina_locations.csv : lat-lon locations of malaspina circumnavigation sample locations\
node_purity.csv : output of random forest model for prokaryotic abundance\
surface_chl.csv : lat-lon surface chlorophyll, averaged over 2002-2016, from MODIS-AQUA.\
year_data_processed_epipelagic.csv : mean global environmental variables and total heterotroph biomass and respiration from 1980-2100 under SSP3-7.0, across four different climate models

**/figures summary:**\
Folder where figures are saved.

**/models summary:**\
gam_abundance.RData : gam of prokaryotic abundance\
gam_cellcarbon.RData : gam of prokaryotic cell carbon\
gam_chl_sgr.RData : gam of prokaryotic sgr (in surface waters with chlorophyll measurements only)\
gam_sst_sgr.RData : gam of prokaryotic sgr (in all waters with temperature measurement, not just surface waters)\
glmm_cellcarbon.RData : glmm of prokaryotic cell carbon\
lm_abundance.RData : parametric model of prokaryotic abundance\
lmer_chl_sgr.RData : lmer of prokaryotic sgr (in surface waters with chlorophyll measurements only)\
lmer_sst_sgr.RData : lmer of prokaryotic sgr (in all waters with temperature measurement, not just surface waters)

**/scripts summary:**\
figurexx.R - generates figure xx\
plot_map_func.R - script to generate global maps, with Robinson projection (shapefiles in ./land_shapefiles)\
prokaryote_abundance_stat_model.R - script to build and evaluate prokaryotic abundance statistical models\
prokaryote_carbon_demand_stat_model.R - script to build and evaluate prokaryotic specific growth rate statistical models, and calculate median prokaryotic growth efficiency\
prokaryote_cell_carbon_stat_model.R - script to build and evaluate prokaryotic cell carbon statistical models

**/source_data summary:**\
Folder where source_data for figures are saved. This is different to files in ./data: it's the data in the exact format it appears in the figures (as required by Nature Communications).
