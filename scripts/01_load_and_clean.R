######################################################
# AMT Spatial Conservation Prioritisation for Mammals
#
# Part 1: Load and clean data for PGLS models
#
# Code developed by Emmeline Norris, 01/06/2024
# 
######################################################
install.packages('cowplot')
library(sf)
library(terra)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(ape)
library(usdm)
library(ggthemes)
library(ggspatial)
library(cowplot)

# Disable scientific notation
options(scipen = 999)

#### AUSTRALIAN MONSOONAL TROPICS AND STATE BOUNDARIES ####

# Generate map of Australia (to plot in ggplot)
aust <- ne_countries(scale = "large", returnclass = "sf", country = 'australia') %>%
  st_transform(crs = st_crs(3577))
# Check the CRS proj4string
st_crs(aust)$proj4string
# "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
ggplot(data = aust) + 
  geom_sf()

# Define AMT boundaries by IBRA Bioregions
# Load IBRA bioregions
ibra <- read_sf("data/input-data/ibra7_regions/ibra7_regions.shp")

# Create a dataframe of the 89 IBRA region names, with empty columns for max and min latitude
ibra_names <- data.frame(name = ibra$REG_NAME_7, north = NA, south = NA)

# Extract max and min latitudes for all the IBRA bioregions
for(i in 1:length(ibra_names$name)){
  ibra_names$north[i] <- ext(ibra[ibra$REG_NAME_7==ibra_names$name[i],])[4]
  ibra_names$south[i] <- ext(ibra[ibra$REG_NAME_7==ibra_names$name[i],])[3]
}

# Filter bioregions to those between 17°-20°S latitude (AMT)
# Exclude Wet Tropics which is biogeographically unique compared to AMT
# Exclude Indian Tropical Islands which are geographically isolated from AMT
ibra_amt <- ibra_names %>%
  filter(ibra_names$north > (-17) & ibra_names$south > (-21.72)) %>%
  filter(!(name %in% c("Wet Tropics", "Indian Tropical Islands")))

# Perform a join and create a new column indicating if the bioregion is in the AMT
ibra <- ibra %>%
  mutate(is_amt = if_else(REG_NAME_7 %in% ibra_amt$name, "Yes", "No"))

# Plot the IBRA bioregions with AMT bioregions coloured green
ggplot() +
  geom_sf(data = ibra, aes(fill = is_amt, color = "gray40")) +
  scale_fill_manual(values = c("Yes" = "green", "No" = "grey")) +
  theme_classic()     
# This now corresponds to the boundaries of the AMT as defined in Brundrett (2017).
# DOI: https://doi.org/10.1007/978-3-319-56363-3_17

# Dissolve the AMT IBRA polygons into a single polygon representing the boundary of the AMT
# Buffer IBRA bioregions by 0.1 m to make geometry valid when doing union, then transform to AEA
amt <- ibra %>%
  filter(is_amt == "Yes") %>%
  st_buffer(0.1) %>%
  st_union() %>%
  st_make_valid() %>%
  st_transform(crs = st_crs(3577)) %>%
  st_intersection(aust) %>%
  st_as_sf()

ibra_amt <- ibra %>%
  filter(is_amt == "Yes")

# Crop bounding box to remove Badu, Moa, Boigu and Saibai islands
bbox <- st_bbox(c(xmin = -1315000, ymin = -2360000, xmax = 1619000, ymax = -1150000), crs = st_crs(amt))
amt <- st_crop(amt, bbox)
# Check plot
ggplot(data = amt) + 
  geom_sf()
# Write to shapefile
st_write(amt, "data/output-data/shp/01_amt.shp", append = F)

# Extract the Wet Tropics bioregion from IBRA regions and transform to AEA
wt <- ibra %>%
  filter(REG_NAME_7 == "Wet Tropics") %>%
  st_transform(crs = st_crs(3577)) %>%
  st_intersection(aust) %>%
  st_as_sf()
# Check plot
ggplot(data = wt) + 
  geom_sf()
# Write to shapefile
st_write(wt, "data/output-data/shp/01_wet-tropics.shp", append = F)

#### PRODUCE PLOT SHOWING IBRA AMT DELINEATION in WGS PROJECTION ####

# Plot AMT bioregions
ggplot() +
  geom_sf(data = aust, fill = 'gray95') +
  geom_sf(data = amt, fill = '#c4d5c2') +
  geom_sf(data = wt, fill = 'gray20') +
  geom_sf(data = ibra, fill = 'transparent', color = 'gray50') +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(crs = 4326, xlim = c(110, 160), ylim = c(-44, -11)) +
  theme_minimal() +
  theme(axis.title = element_text(size=12),
        axis.text = element_text(size = 10))

ggsave("01_amt_ibra_bioregions.png", units="cm", width=20, height=16, dpi=600, path = "results/fig/amt", bg  = 'white')


# Now plot AMT with important places:

# Create a data frame for cities
cities <- data.frame(
  name = c("DARWIN", "CAIRNS"),
  lat = c(-12.4634, -16.9186),
  lon = c(130.8456, 145.7781))

# Create a data frame for cities
city_names <- data.frame(
  name = c("DARWIN", "CAIRNS"),
  lat = c(-12.2, -16.5),
  lon = c(129.7, 146.6))

# Create a data frame for towns
towns <- data.frame(
  name = c("Broome", "Katherine"),
  lat = c(-17.9610, -14.4641),
  lon = c(122.2355, 132.2638))
# Create a data frame for towns
town_names <- data.frame(
  name = c("Broome", "Katherine"),
  lat = c(-17.7, -14.3),
  lon = c(121.4, 133.3))

# Create a data frame for regions
regions <- data.frame(
  name = c("Cape\nYork\nPeninsula", "Kimberley", 
           "Top End"),
  lat = c(-16, -16, -12.8),
  lon = c(142.9,127, 133))

# Create a data frame for gaps
gaps <- data.frame(
  name = c("Bonaparte\nGap", "Gulf of\nCarpentaria"),
  lat = c(-13.3, -15.0000),
  lon = c(128.3, 139.0000))

# Extract centroids for bioregions
ibra_amt <- ibra_amt %>%
  mutate(centroid = st_centroid(geometry))

# Crop bounding box to remove Badu, Moa, Boigu and Saibai islands
aust_bbox <- st_bbox(c(xmin = 119, ymin = -23, xmax = 148, ymax = -8), crs = st_crs(4326))

aust_bbox_bound <- st_as_sfc(st_bbox(aust_bbox))

aust_crop <- aust %>%
  st_transform(crs = st_crs(4326)) %>%
  st_crop(aust_bbox)

aust_name <- data.frame(name = c("AUSTRALIA"),
                        lat = c(-20), lon = c(133))

# Create sf objects for cities, towns, regions, and gaps
cities_sf <- st_as_sf(cities, coords = c("lon", "lat"), crs = 4326)
city_names_sf <- st_as_sf(city_names, coords = c("lon", "lat"), crs = 4326)
towns_sf <- st_as_sf(towns, coords = c("lon", "lat"), crs = 4326)
town_names_sf <- st_as_sf(town_names, coords = c("lon", "lat"), crs = 4326)
regions_sf <- st_as_sf(regions, coords = c("lon", "lat"), crs = 4326)
gaps_sf <- st_as_sf(gaps, coords = c("lon", "lat"), crs = 4326)
aust_name_sf <- st_as_sf(aust_name, coords = c("lon", "lat"), crs = 4326)

# Create the map
main_map <- ggplot() +
  geom_sf(data = aust_bbox_bound, fill = "#b4e7ff", color = "#b4e7ff") +
  geom_sf(data = aust_crop, fill = "#f8ebdd", color = "#948d84") +
  geom_sf_text(data = aust_name_sf, aes(label = name), size = 5, fontface = "bold", color = 'gray70') +
  geom_sf(data = wt, fill = '#3e4a3d') +
  geom_sf(data = ibra_amt, fill = '#c4d5c2', color = '#8da78a', alpha = 0.5) +
  geom_sf_text(data = ibra_amt, aes(label = REG_CODE_7), size = 3, color = "#8c662a", alpha = 0.7, fontface = "bold") +
  geom_sf(data = cities_sf, color = "red", size = 3) +
  geom_sf_text(data = city_names_sf, aes(label = name), size = 3, fontface = "bold") +
  geom_sf(data = towns_sf, color = "black", size = 2) +
  geom_sf_text(data = town_names_sf, aes(label = name), size = 2.5) +
  geom_sf_text(data = regions_sf, aes(label = name), size = 4, fontface = "bold", color = "#3e4a3d") +
  geom_sf_text(data = gaps_sf, aes(label = name), size = 3, fontface = "italic", color = 'gray40') +
  xlab("Longitude") + ylab("Latitude") +
  theme_classic() +
  annotation_scale(location = "bl", width_hint = 0.5) +
  coord_sf(expand = FALSE, crs = st_crs(4326), xlim = c(119, 148), ylim = c(-8, -22)) +  
  theme(axis.ticks.length=unit(-0.25, "cm"),
        axis.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "none",
        #axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size = 0.5),
        axis.line = element_line(size = 0.5),
        #axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
plot(main_map)

ggsave("01_amt_POI_20240903.png", units="cm", width=25, height=10, dpi=300, path = "results/fig/amt", bg  = 'white')


# Create the inset map showing AMT in Australia
inset_map <- ggplot() +
  geom_sf(data = aust, fill = "#f8ebdd", color = "#948d84") +
  geom_sf(data = aust_bbox_bound, 
          fill = "darkgreen",
          color = 'transparent',
          alpha = 0.2) +
  coord_sf(crs = st_crs(4326), 
           xlim = c(112, 155), 
           ylim = c(-6, -45),
           expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = "black"))
plot(inset_map)


ggdraw(main_map) +
  draw_plot(
    {inset_map},
    # The distance along a (0,1) x-axis to draw the left edge of the plot
    x = 0.01, 
    # The distance along a (0,1) y-axis to draw the bottom edge of the plot
    y = 0.58,
    # The width and height of the plot expressed as proportion of the entire ggdraw object
    width = 0.28, 
    height = 0.28)

ggsave("01_amt_mapinset_20240903.png", units="cm", width=25, height=14, dpi=600, path = "results/fig/amt", bg  = 'white')



#### RASTER GRID FOR AUSTRALIA AND AMT ####

# Generate raster grids for Australia and the AMT
# Define the extent of the raster grid as the extent of the polygon of Australia
r_ext <- ext(amt) 

# Create a raster grid with a 1 km x 1 km spatial resolution
r <- rast(amt, ext = r_ext, resolution = 1000)
res(r) # Check resolution
crs(r) # Check CRS is Albers Equal Area

# Rasterize polygon of AMT to raster grid and crop extent
r_amt <- rasterize(amt, r)
plot(r_amt)

# Export the raster as a TIFF file
writeRaster(r_amt, filename="data/output-data/tif/01_r_amt.tif", overwrite=TRUE)


#### IUCN SPECIES DISTRIBUTIONS AND THREAT STATUS ####

# Load IUCN polygons shapefile for all Australian mammals
# IUCN data were filtered by:
# Taxonomy: Mammalia
# Land Region: Australia
# Habitats: all terrestrial habitats
# Country Legends: Extant (resident), Extant & Origin Uncertain (Pteropus macrotis)
# Total = 283 species
iucn <- read_sf("data/input-data/iucn/Polys23/data_0.shp") %>%
  st_transform(crs = st_crs(3577))

# Load IUCN species summary csv (includes taxonomic details and threat category)
iucn_summary = read_csv("data/input-data/iucn/Summary/simple_summary.csv")

# Change SCI_NAMES column name to match scientificName in summary file
colnames(iucn)[3] <- 'scientificName'
head(iucn)

# Intersect with AMT polygon to find species that inhabit the AMT
# Note: this crops distribution areas to AMT boundary
iucn_amt <- st_intersection(iucn, amt)

# Create new sf data frame for AMT species that includes all range polygons, including outside AMT
iucn_amt_full <- iucn %>%
  filter(scientificName %in% iucn_amt$scientificName) 

# Merge summary data with IUCN species distribution spatial dataset
iucn_amt_full <- inner_join(iucn_summary, iucn_amt_full, by = 'scientificName')

# Create ordinal 'threat status' variable
# Convert redlistCategory combined with populationTrend to a numeric ordinal  variable signifying the relative threat status of each species
# Remove species listed under the 'Data Deficient' Redlist category
iucn_amt_full <- iucn_amt_full %>%
  filter(!(redlistCategory == 'Data Deficient')) %>%
  mutate(ordinalThreat = case_when(
    redlistCategory == 'Least Concern' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 1,
    redlistCategory == 'Least Concern' & populationTrend == 'Decreasing' ~ 2,
    redlistCategory == 'Near Threatened' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 3,
    redlistCategory == 'Near Threatened' & populationTrend == 'Decreasing' ~ 4,
    redlistCategory == 'Vulnerable' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 5,
    redlistCategory == 'Vulnerable' & populationTrend == 'Decreasing' ~ 6,
    redlistCategory == 'Endangered' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 7,
    redlistCategory == 'Endangered' & populationTrend == 'Decreasing' ~ 8,
    redlistCategory == 'Critically Endangered' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 9,
    redlistCategory == 'Critically Endangered' & populationTrend == 'Decreasing' ~ 10,
    TRUE ~ NA_integer_
  )) 

amt_species <- unique(iucn_amt_full$scientificName) # 176 species

# Remove species listed under Redlist criterion B (geographic range size) to reduce circularity in extinction risk models (14 species)
iucn_mod <- iucn_amt_full %>%
  filter(is.na(redlistCriteria) | !str_detect(redlistCriteria, "B")) %>%
  st_as_sf()
# Check number of species remaining
model_species_v1 <- unique(iucn_mod$scientificName) # 162 species 

# Export data frame of remaining 162 AMT species to use when downloading phylogenies from VertLife
amt_species_df <- as.data.frame(unique(iucn_mod$scientificName))
colnames(amt_species_df)[1] <- "scientificName"
# Add column of species' ordinal threat status
amt_species_df <- amt_species_df %>%
  mutate(ordinalThreat = iucn_mod$ordinalThreat[match(scientificName, iucn_mod$scientificName)])

# Create a data frame of unioned species distribution polygons
## Initialise an empty list to store the unioned polygons
unioned_polys <- list()
for (species in model_species_v1) {
  # Select polygons for the current species
  species_polys <- iucn_mod[iucn_mod$scientificName == species, ]
  # Intersect with aust polygon to remove parts of distribution outside Australia
  species_polys <- st_intersection(species_polys, aust)
  # Make geometries valid
  species_polys <- st_make_valid(species_polys)
  # Perform the union operation
  curr_union <- st_union(species_polys)
  # Store the result in the list
  unioned_polys[[species]] <- st_sf(scientificName = species, geometry = st_sfc(curr_union))
}

# Combine the list into a single sf data frame
iucn_union <- do.call(rbind, unioned_polys)

# Merge data frame containing species' ordinal threat status with the unioned species distribution polygons
iucn_union <- iucn_union %>%
  left_join(amt_species_df, by = "scientificName")

# Map average and net ordinal threat status across the AMT
# Initialize an empty list to store the raster layers
threat_rasters <- list()
for (i in 1:length(iucn_union$scientificName)) {
  curr_threat <- iucn_union$ordinalThreat[i]
  curr_vect <- vect(iucn_union[i, "geometry"]) # Convert the current geometry to a terra vector object
  curr_rast <- rasterize(curr_vect, r_amt, field = curr_threat, background = NA) # Rasterize the vector object using the template raster and threat vals
  threat_rasters[[i]] <- curr_rast # Store the raster in the list
}
# Combine the list of rasters into a single raster stack
threat_stack <- rast(threat_rasters)

# Calculate the net current threat status, excluding NA values
sum_curr_threat <- sum(threat_stack, na.rm = TRUE) %>%
  mask(r_amt)
plot(sum_curr_threat)

# Calculate the mean current threat status, excluding NA values
avg_curr_threat <- mean(threat_stack, na.rm = TRUE) %>%
  mask(r_amt)
plot(avg_curr_threat)


#### MAMMAL PHYLOGENY ####

# Import subset of 100 trees from a Bayesian analysis pruned to Australian 
# terrestrial mammals (271 sp); Downloaded from http://vertlife.org/data/mammals/
phy <- read.nexus("data/input-data/mammal_phy/Pruned/100_terrestrial_mammal_tree_matching_IUCN.nex")

# Extract names of species contained in the mammal phylogeny
phy_names <- phy[[1]]$tip.label # 271 species
phy_names <- gsub("_"," ", phy_names) # Remove underscores in species names

# Find species names that are in both the phylogeny and the subset of 162 AMT species to be included in model
amt_phy_names <- intersect(model_species_v1, phy_names) # 146 species (16 species removed)

# List the species excluded from the list of 162 AMT mammal species
setdiff(model_species_v1, amt_phy_names)
# List the species excluded from the mammal phylogeny
setdiff(phy_names, amt_phy_names)

# Filter unioned IUCN polygons to include only species that are in the mammal phylogeny
iucn_union_mod <- iucn_union %>%
  filter(scientificName %in% amt_phy_names) %>%
  arrange(scientificName)  # order rows alphabetically by scientificName

# Find species remaining 
model_species_v2 <- unique(iucn_union_mod$scientificName) # 146 species

####'COMBINE' MAMMAL TRAIT DATABASE ####

# Load 'COMBINE' mammal trait database
combine <- read.csv("data/input-data/combine/COMBINE_archives/trait_data_imputed.csv")
colnames(combine)[5] <- "scientificName"

combine <- combine %>%
  # Filter to 148 AMT mammals that will be included in the model
  filter(scientificName %in% model_species_v2) %>% 
  # Select columns of interest 
  # scientificName, adult_mass_g, age_first_reproduction_d, litter_size_n, litters_per_year_n, generation_length_d
  dplyr::select(5,7,15,18,19,24) %>% 
  # Remove rows where trait data are missing (none in this case)
  drop_na() 

# Log transform and standardise all trait variables from the COMBINE database
st_combine <- combine %>%
  mutate(
    adult_mass_g_tr = log(adult_mass_g),
    adult_mass_g_tr_st = as.vector(scale(adult_mass_g_tr)),
    age_first_reproduction_d_tr = log(age_first_reproduction_d),
    age_first_reproduction_d_tr_st = as.vector(scale(age_first_reproduction_d_tr)),
    litter_size_n_tr = log(litter_size_n),
    litter_size_n_tr_st = as.vector(scale(litter_size_n_tr)),
    litters_per_year_n_tr = log(litters_per_year_n),
    litters_per_year_n_tr_st = as.vector(scale(litters_per_year_n_tr)),
    generation_length_d_tr = log(generation_length_d),
    generation_length_d_tr_st = as.vector(scale(generation_length_d_tr))
  )

#### RED FOX OCCURRENCE DATA AND PREDICTED DISTRIBUTION ####

# Load red fox occurrence data from Atlas of Living Australia (ALA)
fox_records <- read.csv("data/input-data/ala/V_vulpes/V_vulpes_records-2023-05-21.csv", 
  stringsAsFactors=FALSE) %>%
  dplyr::select(103,105,108) # Select columns for latitude, longitude, & coordinate uncertainty

# Remove points without coordinate data:
fox_records <- fox_records %>%
  filter(!(is.na(decimalLatitude & decimalLongitude))) %>% # Remove NA coordinates
  filter(decimalLatitude>(-40) & decimalLatitude<(-15)) %>% # Exclude Tasmanian records
  filter(coordinateUncertaintyInMeters <= 2000) %>% # Remove records with coordinate uncertainty > 2km
  dplyr::select(!3) # Remove coordinate uncertainty column

# Convert the occurrence record data frame to an sf points object
fox_points <- fox_records %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% # ALA records are in WGS84
  st_transform(crs = st_crs(3577)) %>% 
  distinct() # Remove duplicate points

# Generate minimum convex hull polygon for fox distribution
fox_convex_hull <- st_convex_hull(st_union(fox_points)) %>%
  st_intersection(aust) 
fox_convex_hull <- st_as_sf(as.data.frame(fox_convex_hull))

# Plot the convex hull polygon and species points
ggplot() +
  geom_sf(data = aust, color = "black", alpha = 0.5) +
  geom_sf(data = amt, color = "blue", fill = 'transparent', alpha = 0.5) +
  geom_sf(data = fox_convex_hull, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_sf(data = fox_points, color = "red", size = 1) +
  theme_bw()

# Calculate proportion of proportion of AMT species' ranges (Australian distribution) that overlaps with fox distribution
outfile_fox <- data.frame("scientificName" = model_species_v2,"foxAreaProp" = NA)
for (i in 1:length(model_species_v2)) {
  curr_shape <- iucn_union_mod[iucn_union_mod$scientificName == model_species_v2[i],]
  if (nrow(curr_shape) > 0) {
    fox_area <- st_intersection(fox_convex_hull, curr_shape)
    if (nrow(fox_area) > 0) {
      fox_area <- st_area(fox_area) / st_area(curr_shape)
    } else {
      fox_area <- 0
    }
  } else {
    fox_area <- 0
  }
  outfile_fox$foxAreaProp[i] <- fox_area
}
outfile_fox$foxAreaProp <- round(outfile_fox$foxAreaProp, 4)

# Logit transform and standardise the proportion of species' range overlap with foxes
st_fox <- outfile_fox %>%
  mutate(
    foxAreaProp_tr = ifelse(foxAreaProp == 0, log((foxAreaProp + 0.0001) / (1 - (foxAreaProp + 0.0001))),
                            ifelse(foxAreaProp == 1, log((foxAreaProp - 0.0001) / (1 - (foxAreaProp - 0.0001))),
                                   log(foxAreaProp / (1 - foxAreaProp)))),
    foxAreaProp_tr_st = as.vector(scale(foxAreaProp_tr))
  )

#### CANE TOAD OCCURRENCE DATA AND PREDICTED DISTRIBUTION ####

# Load cane toad occurrence data from Atlas of Living Australia (ALA)
toad_records <- read.csv("data/input-data/ala/R_marina/R_marina_records-2023-05-21.csv", 
                        stringsAsFactors=FALSE) %>%
  dplyr::select(103,105,108) # Select columns for latitude, longitude, & coordinate uncertainty

# Remove points without coordinate data:
toad_records <- toad_records %>%
  filter(!(is.na(decimalLatitude & decimalLongitude))) %>% # Remove NA coordinates
  filter(decimalLatitude>(-34.5) & decimalLatitude<(-15)) %>% # Exclude single record in South Australia
  filter(decimalLongitude>(120) & decimalLongitude<(160)) %>% # Exclude two records in far west of Western Australia, and one point in ocean
  filter(coordinateUncertaintyInMeters <= 2000) %>% # Remove records with coordinate uncertainty > 2km
  dplyr::select(!3) # Remove coordinate uncertainty column

# Convert the occurrence record data frame to an sf points object
toad_points <- toad_records %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% # ALA records are in WGS84
  st_transform(crs = st_crs(3577)) %>% 
  st_intersection(aust) %>%
  distinct() # Remove duplicate points

# Generate minimum convex hull polygon for cane toad distribution
toad_convex_hull <- st_convex_hull(st_union(toad_points)) %>%
  st_intersection(aust) 
toad_convex_hull <- st_as_sf(as.data.frame(toad_convex_hull))

# Plot the convex hull polygon and species points
ggplot() +
  geom_sf(data = aust, color = "black", alpha = 0.5) +
  geom_sf(data = toad_convex_hull, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_sf(data = toad_points, color = "red", size = 1) +
  theme_bw()

# Calculate proportion of proportion of AMT species' ranges (Australian distribution) that overlaps with cane toad convex hull
outfile_toad <- data.frame("scientificName" = model_species_v2,"toadAreaProp" = NA)
for (i in 1:length(model_species_v2)) {
  curr_shape <- iucn_union_mod[iucn_union_mod$scientificName == model_species_v2[i],]
  if (nrow(curr_shape) > 0) {
    toad_area <- st_intersection(toad_convex_hull, curr_shape)
    if (nrow(toad_area) > 0) {
      toad_area <- st_area(toad_area) / st_area(curr_shape)
    } else {
      toad_area <- 0
    }
  } else {
    toad_area <- 0
  }
  outfile_toad$toadAreaProp[i] <- toad_area
}
outfile_toad$toadAreaProp <- round(outfile_toad$toadAreaProp, 4)

# Logit transform and standardise the proportion of species' range overlap with cane toads
st_toad <- outfile_toad %>%
  mutate(
    toadAreaProp_tr = ifelse(toadAreaProp == 0, log((toadAreaProp + 0.0001) / (1 - (toadAreaProp + 0.0001))),
                             ifelse(toadAreaProp == 1, log((toadAreaProp - 0.0001) / (1 - (toadAreaProp - 0.0001))),
                                    log(toadAreaProp / (1 - toadAreaProp)))),
    toadAreaProp_tr_st = as.vector(scale(toadAreaProp_tr))
  )

#### CALCULATE PROPORTION RANGE OVERLAP WITH AMT ####

# Calculate total species distribution area 
outfile_overlap <- data.frame(scientificName = model_species_v2, area_total = NA, area_amt = NA, percent_overlap = NA)
for (i in 1:length(model_species_v2)) {
  curr_shape <- iucn_union_mod[i, "geometry"]
  curr_area <- st_area(curr_shape)
  outfile_overlap$area_total[i] <- curr_area
}

# Calculate area of species distribution that intersects with the AMT
for (i in 1:length(model_species_v2)) {
  curr_shape <- iucn_union_mod[i, "geometry"]
  curr_intersect <- st_intersection(curr_shape, amt)
  curr_area <- st_area(curr_intersect)
  outfile_overlap$area_amt[i] <- curr_area
}

# Calculate the proportion of each species' distributions that overlap with the AMT
outfile_overlap <- outfile_overlap %>%
  mutate(percent_overlap = (area_amt/area_total)*100)

# Remove species with < 2% range overlap with AMT 5
filtered_overlap <- outfile_overlap %>%
  filter(percent_overlap >= 2)

# Log transform and standardise geographic range area
st_overlap <- filtered_overlap %>%
  mutate(percentOverlap_st = as.vector(scale(percent_overlap)))


# Find species remaining 
model_species_v3 <- unique(filtered_overlap$scientificName) # 128 species

# Filter unioned IUCN polygons to include 128 mammal species with >2% AMT overlap
iucn_union_mod <- iucn_union_mod %>%
  filter(scientificName %in% model_species_v3) %>%
  arrange(scientificName)

# Write to shapefile
st_write(iucn_union_mod, "data/output-data/shp/01_iucn_union_mod.shp", append = F)

#### NAFI FIRE FREQUENCY DATA 2000 - 2022 ####

# Load geoTIFF of fire frequency 2000 - 2022 and project to raster grid of AMT
firefreq_00_22 <- rast("data/input-data/nafi/2000-2022/Fire_Freq_00-22/FF_2000_22_gda94.tif")
firefreq_00_22 <- project(firefreq_00_22, r_amt) %>%
                    mask(r_amt)
plot(firefreq_00_22)

# Calculate mean annual fire frequency for 1 x 1 km grid cells in species' distributions
outfile_firefreq <- data.frame("scientificName" = model_species_v3,"mean_firefreq"=NA)
for (i in 1:length(model_species_v3)) {
  # Extract the current polygon
  curr_vect <- vect(iucn_union_mod[i, "geometry"]) 
  # Mask the fire frequency raster by the current polygon
  curr_rast <- mask(firefreq_00_22, curr_vect)
  # Calculate the mean of the raster values within the polygon
  mean_firefreq <- global(curr_rast, fun = mean, na.rm = TRUE)$mean
  # Append the results to the output data frame
  outfile_firefreq$mean_firefreq[i] <- mean_firefreq
}

# Standardise total fire frequency
st_firefreq <- outfile_firefreq %>%
  mutate(
    mean_firefreq_st = as.vector(scale(mean_firefreq))
  )


# Load geoTIFF of late (post-July) fire frequency 2000 - 2022 and project to raster grid of AMT
firefreqLate_00_22 <- rast("data/input-data/nafi/2000-2022/Fire_Freq_00-22_post_July/FFL_2000_22_gda94.tif")
firefreqLate_00_22 <- project(firefreqLate_00_22, r_amt) %>%
                        mask(r_amt)
plot(firefreqLate_00_22)

# Calculate mean annual frequency of late dry season fires for 1 x 1 km grid cells in species' distributions
outfile_latefire <- data.frame("scientificName" = iucn_union_mod$scientificName,"mean_latefire"=NA)
for (i in 1:nrow(iucn_union_mod)) {
  curr_vect <- vect(iucn_union_mod[i, "geometry"]) 
  curr_rast <- mask(firefreqLate_00_22, curr_vect)
  mean_latefire <- global(curr_rast, fun = mean, na.rm = TRUE)$mean
  outfile_latefire$mean_latefire[i] <- mean_latefire
}

# Standardise late-dry season fire frequency
st_latefire <- outfile_latefire %>%
  mutate(
    mean_latefire_st = as.vector(scale(mean_latefire))
  )

#### WORLDCLIM BIOCLIMATIC VARIABLES ####

## BIO8: MEAN TEMPERATURE OF THE WETTEST QUARTER ##
bio8 <- rast("data/input-data/worldclim/biovar/wc2.1_2.5m_bio_8.tif")
bio8 <- project(bio8, r_amt) %>%
  mask(r_amt)

# Calculate the mean temperature of wettest quarter in species' distribution areas
outfile_bio8 <- data.frame("scientificName" = iucn_union_mod$scientificName, "area_km2" = NA, "mean_bio8" = NA)
for(i in 1:nrow(iucn_union_mod)){
  curr_shape <- iucn_union_mod[i, "geometry"]
  outfile_bio8$area_km2[i] <- as.numeric(st_area(curr_shape)) / 1e6
  curr_rast <- rasterize(curr_shape, r_amt)
  curr_rast <- mask(bio8, curr_rast)
  values <- as.vector(curr_rast[])
  outfile_bio8$mean_bio8[i] <- mean(values, na.rm = TRUE)
}

# Standardise temperature of the wettest quarter
st_bio8 <- outfile_bio8 %>%
  mutate(
    mean_bio8_st = as.vector(scale(mean_bio8))
  )

## BIO17: MEAN PRECIPITATION OF THE DRIEST QUARTER ##
bio17 <- rast("data/input-data/worldclim/biovar/wc2.1_2.5m_bio_17.tif")
bio17 <- project(bio17, r_amt) %>%
  mask(r_amt)
plot(bio17)

# Calculate the mean precipitation of the driest quarter in species' distribution areas
outfile_bio17 <- data.frame("scientificName" = iucn_union_mod$scientificName, "area_km2" = NA, "mean_bio17" = NA)
for(i in 1:nrow(iucn_union_mod)){
  curr_shape <- iucn_union_mod[i, "geometry"]
  outfile_bio17$area_km2[i] <- as.numeric(st_area(curr_shape)) / 1e6
  curr_rast <- rasterize(curr_shape, r_amt)
  curr_rast <- mask(bio17, curr_rast)
  values <- as.vector(curr_rast[])
  outfile_bio17$mean_bio17[i] <- mean(values, na.rm = TRUE)
}

# Log transform and standardise precipitation of the driest quarter
st_bio17 <- outfile_bio17 %>%
  mutate(
    mean_bio17_tr = log(mean_bio17),
    mean_bio17_tr_st = as.vector(scale(mean_bio17_tr))
  )

#### HUMAN INFLUENCE INDEX ####

# Load Human Influence Index data for oceania
hii <- rast("data/input-data/hii/hii-oceania-geo-grid/hii_ocean_grid/hii_oceania/w001001.adf")
hii <- project(hii, r_amt) %>%
  mask(r_amt)
plot(hii)

# Calculate median HII for species range polygons
outfile_hii <- data.frame("scientificName" = iucn_union_mod$scientificName, "mean_hii" = NA)
for(i in 1:nrow(iucn_union_mod)){
  curr_shape <- iucn_union_mod[i, "geometry"]
  curr_rast <- rasterize(curr_shape, r_amt)
  curr_rast <- mask(hii, curr_rast)
  values <- as.vector(curr_rast[])
  outfile_hii$mean_hii[i] <- mean(values, na.rm = TRUE)
}

# Square-root transform and standardise human influence index
st_hii <- outfile_hii %>%
  mutate(
    mean_hii_tr = sqrt(mean_hii),
    mean_hii_tr_st = as.vector(scale(mean_hii_tr))
  )

#### AMT MAMMALS GEOGRAPHIC RANGE AREA ####

# Calculate geographic range area for native terrestrial mammals in the AMT
outfile_range <- data.frame("scientificName" = iucn_union_mod$scientificName,"rangeArea_km2" = NA)
for (i in 1:nrow(iucn_union_mod)) {
  curr_shape <- iucn_union_mod[i, "geometry"]
  outfile_range$rangeArea_km2[i] <- as.numeric(st_area(curr_shape)) / 1e6
}

# Log transform and standardise geographic range area
st_range <- outfile_range %>%
  mutate(
    rangeArea_km2_tr = log(rangeArea_km2),
    rangeArea_km2_tr_st = as.vector(scale(rangeArea_km2_tr)))

#### EXTRACT RANGE POLYGON CENTROIDS ####

# Extract centroid coordinates for species range polygons
outfile_centroid <- data.frame("scientificName" = iucn_union_mod$scientificName, centroid_lat = NA, centroid_long = NA)
for(i in 1:nrow(iucn_union_mod)){
  curr_shape <- iucn_union_mod[i, "geometry"]
  centroid <- st_centroid(curr_shape)
  outfile_centroid$centroid_long[i] <- st_coordinates(centroid)[1]
  outfile_centroid$centroid_lat[i] <- st_coordinates(centroid)[2]
}
centroids <- st_as_sf(outfile_centroid, coords = c("centroid_long", "centroid_lat")) %>%
  st_set_crs(crs(amt))

# Check that centroids look correct when plotted
ggplot() +
  geom_sf(data = aust, color = "black", alpha = 0.5) +
  geom_sf(data = centroids, color = "blue") +
  theme_bw()

#### CREATE A PREDICTOR DATA FRAME FOR AMT MAMMALS ####

# Filter data frame of species names and ordinal threat status to include 143 mammals in the phylogeny
amt_species_df <- amt_species_df %>%
  filter(scientificName %in% model_species_v3) %>%
  arrange(scientificName)
# Write to csv
write_csv(amt_species_df, "data/output-data/tbl/01_amt_species_df.csv", append = F)

# List of data frames to merge
data_frames <- list(
  amt_species_df,
  outfile_centroid,
  st_overlap,
  st_combine,
  st_range,
  st_fox,
  st_toad,
  st_firefreq,
  st_latefire,
  st_bio8,
  st_bio17,
  st_hii)

# Merge all data frames
amt_dat <- reduce(data_frames, full_join, by = "scientificName") %>%
  filter(scientificName %in% model_species_v3)

# Filter to keep species name, centroid coords, threat status, and transformed/standardised predictors
amt_dat_st <- amt_dat[, c(1:8,15,17,19,21,23,26,29,32,34,36,39,43,46)]
# Remove rows with incomplete data
amt_dat_st <- amt_dat_st[complete.cases(amt_dat_st),] # No rows removed
# Substitute space with underscore in species scientificName to match phylogeny
amt_dat_st$scientificName <- gsub(" ", "_", amt_dat_st$scientificName)

# Prune phylogenetic trees to AMT mammals included in final data frame 'amt_dat'
amt_phy <-lapply(phy, keep.tip, tip = amt_dat_st$scientificName)
# Save the list of 100 phylogenetic trees pruned to AMT mammals
saveRDS(amt_phy, file = "data/output-data/phy/01_amt_phy.rds")

# Test for correlation between predictor variables for 138 AMT mammals (cut-off = 0.7)
corr_mat <- cor(amt_dat_st[,unlist(lapply(amt_dat_st, is.numeric))], use = "complete.obs")

# Run collinearity diagnostics using Variance Inflation Factor (VIF) cut-off of 5
v_dat <- amt_dat_st[, c(8:21)]
vif <- vifstep(v_dat)
vif
# Bio17 variable (precipitation of driest quarter) found to have collinearity problem
# Generation length has high VIF (5.76)

# Remove precipitation of wettest quarter and generation length from data frame of variables for extinction risk model
amt_dat_st <- amt_dat_st[ , -c(13,20)]

# Run collinearity diagnostics again
v_dat <- amt_dat_st[, c(8:19)]
vif <- vifstep(v_dat)
vif # No variables have a collinearity problem and none have VIF > 5

# Write final data frame containing response and predictor variables to .csv file for use in 02_pgls_model.R
write_csv(amt_dat_st, "data/output-data/tbl/01_amt_dat_st.csv", append = F)

