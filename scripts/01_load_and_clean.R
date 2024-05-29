######################################################
# AMT Spatial Conservation Prioritisation for Mammals
# Part 1: Load and clean data
#
# Code developed by Emmeline Norris, 24/05/2024
# 
######################################################

library(sf)
library(terra)
library(tidyverse)
library(rnaturalearth)
library(ggplot2)
library(ape)

# Disable scientific notation
options(scipen = 999)

#### STUDY AREA: AUSTRALIAN MONSOONAL TROPICS AND STATE BOUNDARIES ####

# Generate map of Australia (to plot in ggplot)
aust <- ne_countries(scale = "large", returnclass = "sf", country = 'australia') %>%
  st_transform(crs = st_crs(3577))
# Check the CRS proj4string
st_crs(aust)$proj4string
# "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"
ggplot(data = aust) + 
  geom_sf()

# Load state & territory borders
states <- read_sf("data/input-data/states/STE_2021_AUST_GDA2020.shp") %>%
  st_transform(crs = st_crs(3577))
st_crs(states)$proj4string # "+proj=longlat +ellps=GRS80 +no_defs"
ggplot(data = states) + 
  geom_sf()

# Extract polygons for QLD, NT and WA (AMT states)
wa <- states %>%
  filter(STE_NAME21 == "Western Australia")
nt <- states %>%
  filter(STE_NAME21 == "Northern Territory")
qld <- states %>%
  filter(STE_NAME21 == "Queensland")

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
  geom_sf(data = ibra, aes(fill = is_amt)) +
  scale_fill_manual(values = c("Yes" = "green", "No" = "grey")) +
  theme_classic()     
# This now corresponds to the boundaries of the AMT as defined in 
# https://link.springer.com/chapter/10.1007/978-3-319-56363-3_17

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

# Check plot
ggplot(data = amt) + 
  geom_sf()
# Write to shapefile
st_write(amt, "data/output-data/shp/amt.shp", append = F)


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
st_write(wt, "data/output-data/shp/wet-tropics-bioregion.shp", append = F)


# Generate raster grids for Australia and the AMT
# Define the extent of the raster grid as the extent of the polygon of Australia
r_ext <- ext(aust) 

# Create a raster grid with a 1 km x 1 km spatial resolution
r <- rast(aust, ext = r_ext, resolution = 1000)
r_aust <- rasterize(aust, r)
res(r) # 10km x 10km spatial resolution
plot(r_aust)

# Rasterize grid to polygon of AMT and crop extent
r_amt <- rasterize(amt, r) %>%
  crop(amt)
plot(r_amt)


#### IUCN SPECIES DISTRIBUTION MAPS ####

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

# Change SCI_NAMES column name to match scientificNames in summary file
colnames(iucn)[3] <- 'scientificName'

# Intersect with AMT polygon to find species that inhabit the AMT
iucn_amt <- st_intersection(iucn, amt)

# Create new sf df for AMT species that includes all range polygons, including outside AMT
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
    redlistCategory == 'Least Concern' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 0,
    redlistCategory == 'Least Concern' & populationTrend == 'Decreasing' ~ 1,
    redlistCategory == 'Near Threatened' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 2,
    redlistCategory == 'Near Threatened' & populationTrend == 'Decreasing' ~ 3,
    redlistCategory == 'Vulnerable' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 4,
    redlistCategory == 'Vulnerable' & populationTrend == 'Decreasing' ~ 5,
    redlistCategory == 'Endangered' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 6,
    redlistCategory == 'Endangered' & populationTrend == 'Decreasing' ~ 7,
    redlistCategory == 'Critically Endangered' & (is.na(populationTrend) | populationTrend %in% c('Stable', 'Increasing', 'Unknown')) ~ 8,
    redlistCategory == 'Critically Endangered' & populationTrend == 'Decreasing' ~ 9,
    redlistCategory == 'Extinct' ~ 10,
    TRUE ~ NA_integer_
  )) 

amt_species <- unique(iucn_amt_full$scientificName) # 178 species

# Remove species listed under Redlist criterion B (geographic range size) to reduce circularity in extinction risk models (14 species)
iucn_mod <- iucn_amt_full %>%
  filter(is.na(redlistCriteria) | !str_detect(redlistCriteria, "B")) # 164 species remaining
# Write to shapefile
st_write(iucn_mod, "data/output-data/shp/iucn_model_polys.shp", append = F)

# Export data frame of remaining 164 AMT species to use when downloading phylogenies from VertLife
amt_species_df <- as.data.frame(unique(iucn_mod$scientificName))
colnames(amt_species_df)[1] <- "scientificName"
# Add column of species' ordinal threat status
amt_species_df <- amt_species_df %>%
  mutate(ordinalThreat = iucn_mod$ordinalThreat[match(scientificName, iucn_mod$scientificName)])
# Write to csv
write_csv(amt_species_df, "data/output-data/tbl/amt_model_species.csv", append = F)

# Coerce iucn_mod to sf data frame
iucn_mod <- st_as_sf(as.data.frame(iucn_mod))
model_species <- unique(iucn_mod$scientificName) # 164 species

# Create a data frame of unioned species distribution polygons
# Initialize an empty list to store the unioned polygons
unified_polys <- list()
for (species in model_species) {
  # Select polygons for the current species
  species_polys <- iucn_mod[iucn_mod$scientificName == species, ]
  # Intersect with aust polygon to remove parts of distribution outside Australia
  species_polys <- st_intersection(species_polys, aust)
  # Make geometries valid
  species_polys <- st_make_valid(species_polys)
  # Perform the union operation
  curr_union <- st_union(species_polys)
  # Store the result in the list
  unified_polys[[species]] <- st_sf(scientificName = species, geometry = st_sfc(curr_union))
}

# Combine the list into a single sf data frame
iucn_union <- do.call(rbind, unified_polys)

# Merge the amt_species_df containing species' ordinal threat status with the unified IUCN species distribution polygons
iucn_union <- iucn_union %>%
  left_join(amt_species_df, by = "scientificName")

# Write to shapefile
st_write(iucn_union, "data/output-data/shp/iucn_model_union.shp", append = F)

# Map average and net threat status across the AMT
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

# Calculate the mean and sum of threat_stack, excluding NA values
sum_curr_threat <- sum(threat_stack, na.rm = TRUE) %>%
  mask(r_amt)
avg_curr_threat <- mean(threat_stack, na.rm = TRUE) %>%
  mask(r_amt)
plot(sum_curr_threat)
plot(avg_curr_threat)


#### MAMMAL PHYLOGENY & 'COMBINE' MAMMAL TRAIT DATABASE ####

# Import subset of 100 trees from a Bayesian analysis pruned to Australian 
# terrestrial mammals (271 sp); Downloaded from http://vertlife.org/data/mammals/
phy <- read.nexus("data/input-data/mammal_phy/Pruned/100_terrestrial_mammal_tree_matching_IUCN.nex")

# Extract names of species contained in the mammal phylogeny
phy_names <- phy[[1]]$tip.label # 271 species
phy_names <- gsub("_"," ", phy_names) # Remove underscores in species names

# Find species names that are in the phylogeny and subset of AMT species to be included in model
amt_phy_names <- intersect(model_species, phy_names) # 148 species

# List the species excluded from the original list of 176 AMT mammal species
setdiff(model_species, amt_phy_names)
# List the species excluded from the mammal phylogeny
setdiff(phy_names, amt_phy_names)

# Filter unified IUCN polygons df to include only species that are in the mammal phylogeny
iucn_union_mod <- iucn_union %>%
  filter(scientificName %in% amt_phy_names) %>%
  arrange(scientificName) %>% # order rows alphabetically by scientificName
  st_transform(crs = st_crs(3577))

# Find species remaining (148 species)
model_phy_species <- unique(iucn_union_mod$scientificName)

# Load 'COMBINE' mammal trait database
combine <- read.csv("data/input-data/combine/COMBINE_archives/trait_data_imputed.csv")
colnames(combine)[5] <- "scientificName"

combine <- combine %>%
  # Filter to 148 AMT mammals that will be included in the model
  filter(scientificName %in% model_phy_species) %>% 
  # Select columns of interest (scientificName, adult_mass_g, age_first_reproduction_d, litter_size_n)
  select(5,7,15,18) %>% 
  # Remove rows where trait data are missing (none in this case)
  drop_na() 

#### RED FOX OCCURRENCE DATA AND PREDICTED DISTRIBUTION ####

# Load red fox occurrence data from Atlas of Living Australia (ALA)
fox_records <- read.csv("data/input-data/ala/V_vulpes/V_vulpes_records-2023-05-21.csv", 
  stringsAsFactors=FALSE) %>%
  select(103,105,108) # Select columns for latitude, longitude, & coordinate uncertainty

# Remove points without coordinate data:
fox_records <- fox_records %>%
  filter(!(is.na(decimalLatitude & decimalLongitude))) %>% # Remove NA coordinates
  filter(decimalLatitude>(-40) & decimalLatitude<(-15)) %>% # Exclude Tasmanian records
  filter(coordinateUncertaintyInMeters <= 2000) %>% # Remove records with coordinate uncertainty > 2km
  select(!3) # Remove coordinate uncertainty column

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
  geom_sf(data = fox_convex_hull, fill = "lightblue", color = "black", alpha = 0.5) +
  geom_sf(data = fox_points, color = "red", size = 1) +
  geom_sf(data = iucn_union_mod[1,], fill = "green", color = "black") +
  geom_sf(data = amt, color = "black", alpha = 0.5) +
  theme_bw()

# Calculate proportion of proportion of AMT species' ranges (Australian distribution) that overlaps with fox distribution
outfile_fox <- data.frame("scientificName" = model_phy_species,"foxAreaProp" = NA)
for (i in 1:length(model_phy_species)) {
  curr_shape <- iucn_union_mod[iucn_union_mod$scientificName == model_phy_species[i],]
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

#### CANE TOAD OCCURRENCE DATA AND PREDICTED DISTRIBUTION ####

# Load cane toad occurrence data from Atlas of Living Australia (ALA)
toad_records <- read.csv("data/input-data/ala/R_marina/R_marina_records-2023-05-21.csv", 
                        stringsAsFactors=FALSE) %>%
  select(103,105,108) # Select columns for latitude, longitude, & coordinate uncertainty

# Remove points without coordinate data:
toad_records <- toad_records %>%
  filter(!(is.na(decimalLatitude & decimalLongitude))) %>% # Remove NA coordinates
  filter(decimalLatitude>(-34.5) & decimalLatitude<(-15)) %>% # Exclude single record in South Australia
  filter(decimalLongitude>(120) & decimalLongitude<(160)) %>% # Exclude two records in far west of Western Australia, and one point in ocean
  filter(coordinateUncertaintyInMeters <= 2000) %>% # Remove records with coordinate uncertainty > 2km
  select(!3) # Remove coordinate uncertainty column

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
  geom_sf(data = iucn_union_mod[99,], fill = "green", color = "black") +
  geom_sf(data = amt, color = "black", alpha = 0.5) +
  theme_bw()

# Calculate proportion of proportion of AMT species' ranges (Australian distribution) that overlaps with cane toad convex hull
outfile_toad <- data.frame("scientificName" = model_phy_species,"toadAreaProp" = NA)
for (i in 1:length(model_phy_species)) {
  curr_shape <- iucn_union_mod[iucn_union_mod$scientificName == model_phy_species[i],]
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

