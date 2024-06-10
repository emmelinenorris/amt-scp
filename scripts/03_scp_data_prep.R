######################################################
# AMT Spatial Conservation Prioritisation for Mammals
#
# Part 3: SCP Planning Unit data preparation
#
# Code developed by Emmeline Norris, 01/06/2024
# 
######################################################

library(sf)
library(terra)
library(tidyverse)
library(ggplot2)

options(scipen = 999)

#### GENERATE HEXAGONAL TESSELLATION OF PLANNING UNITS ####

# Load AMT polygon ceated in 01_load_and_clean
amt <- st_read("data/output-data/shp/01_amt.shp")

# Load amt raster from 01_load_and_clean
r_amt <- rast("data/output-data/tif/01_r_amt.tif")

# Create hexagonal planning unit grid for the AMT
# Extract coordinates of AMT polygon bounding box
ext <- st_as_sfc(st_bbox(amt))
# Sample 6000 hexagons across the extent of the AMT 
centre_pts <- st_sample(ext, size = 6000, type = "hexagonal") %>%
  st_transform(crs = st_crs(3577)) %>%
  st_as_sf() 

# Extract coordinates from the sampled hexagons
coords <- st_coordinates(centre_pts)
# Calculate number of unique x and y coordinates
num_x <- length(unique(coords[, "X"]))
num_y <- length(unique(coords[, "Y"]))

# Create a hexagonal tesselation of planning units covering the AMT bounding box
hex_grid <- st_make_grid(amt,
                         n = c(num_x, num_y), # Specify number of hex units in x and y direction
                         what = 'polygons',
                         square = FALSE, # make them hexagons
                         flat_topped = TRUE) %>% # flat-topped, not pointy
            st_as_sf()

# Create a column of sequential id numbers
hex_grid$id <- seq.int(nrow(hex_grid))

# Intersect hex_grid with AMT boundary polygon
amt_grid <- st_intersection(hex_grid, amt) %>%
  mutate(area = st_area(.)) %>% # Calculate area of each hex unit
  dplyr::select(!2) %>% # remove Fid column
  st_collection_extract("POLYGON") # Extract polygon features from GeometryCollection
# 9909 hex units; 166 km^2 area for entire units

# Check size of hex units is spatially constant
# Only hex units on AMT boundary should have a smaller area
ggplot() + 
  geom_sf(data = amt_grid,
          aes(fill = units::drop_units(area)))

# Write to shapefile
st_write(amt_grid, "data/output-data/shp/03_amt_grid.shp", append = F)

# Add columns to the hex data frame for later designating the status of PUs (i.e. locked in or not) and PU cost
amt_grid$capad <- NA
amt_grid$status <- NA # whether PU is already in the reserve system
amt_grid$locked_in <- NA # TRUE/FALSE whether PU is already in the reserve system
amt_grid$locked_out <- NA # TRUE/FALSE whether PU is to be excluded from the reserve system
amt_grid$cost <- NA

#### CAPAD PROTECTED AREAS #####

# Load CAPAD protected area shapefile
capad <- read_sf("data/input-data/capad/CAPAD2022_terrestrial.shp") %>%
  filter(st_is_valid(.)) %>% # Remove invalid polygons
  st_transform(crs = st_crs(3577))

# Intersect CAPAD protected area polygons with the AMT
capad_amt <- st_intersection(capad, amt) %>%
  st_collection_extract("POLYGON")

# Plot by protected area type
ggplot() +
  geom_sf(data = amt) +
  geom_sf(data = capad_amt, aes(fill = TYPE))

# Calculate the proportion of the AMT designated as protected (as per CAPAD)
capad_area <- (st_area(st_union(capad_amt))/st_area(amt))*100 # 22.6% protected

# Generate a protected area raster that retains categories from original shapefile
pa_type <- unique(capad_amt$TYPE)
pa_df <- data.frame(id = 1:length(pa_type), type = pa_type)
capad_amt$id <- pa_df$id[match(capad_amt$TYPE, pa_df$type)]

# Rasterize maximum PA value to AMT raster grid 
# Note: only options for 'fun' are min, max, mean, or sum. 
# Using mean and sum produces meaningless results, so I have opted to use max consistently throughout 
pa_vect <- vect(capad_amt)
pa_ras <- rasterize(pa_vect, r_amt, field = "id", fun = 'max') 
levels(pa_ras) <- pa_df
plot(pa_ras)

# Assign broad protected area categories
pa_df$category <- ifelse(pa_df$type == "Indigenous Protected Area",
                         "Indigenous Protected Area",
                         "Protected Area")

#### PASTORAL AND OTHER LAND TENURE ####

# Import DAWE geotiff of land use (c. 2015-16)
tenure <- rast("data/input-data/dawe/AUSTEN_250m_2015_16_alb.tif")
# Project raster values to amt raster using 'nearest neighbour' method
tenure <- project(tenure, r_amt, method = 'near') %>%
  mask(r_amt)
plot(tenure)

# Generate a data frame with the "Tenure Level 4" land use types and corresponding id codes
tenure_df <- data.frame(
  id = c(0, 1001, 1002, 2111, 2121, 2131, 2141, 2142, 2151, 2161, 2211, 2212, 2221, 2231, 2232, 2301),
  type = c(
    "No data/unresolved", "Freehold", "Freehold - Indigenous", "Freeholding lease", 
    "Pastoral perpetual lease", "Other perpetual lease", "Pastoral term lease", 
    "Pastoral term lease - Indigenous", "Other term lease", "Other lease",
    "Nature conservation reserve", "Nature conservation reserve - Indigenous",
    "Multiple-use public forest", "Other Crown purposes", 
    "Other Crown purposes - Indigenous", "Other Crown land"))
levels(tenure) <- tenure_df
plot(tenure)

# Assign land type categories
tenure_df$category <- ifelse(
  grepl("pastoral", tenure_df$type, ignore.case = TRUE),
  "Pastoral term or perpetual lease",
  ifelse(
    tenure_df$type == "Freehold - Indigenous",
    "Freehold - Indigenous",
    "Other freehold, term, perpetual lease or Crown purposes"
  )
)

#### NATIVE TITLE LAND ####

# Load shapefile of Native Title Determinations
ntd <- read_sf("data/input-data/nntt/NTD_Register_Nat.shp") %>%
  filter(DETOUTCOME %in% c("Native title exists in parts of the determination area",
                           "Native title exists in the entire determination area")) %>%
  st_transform(crs = st_crs(3577)) 

# Intersect NNTT Native Title Determination polygons with the AMT
ntd_amt <- st_intersection(ntd, amt) %>%
  st_collection_extract("POLYGON")

# Plot by Determination Outcome
ggplot() +
  geom_sf(data = amt) +
  geom_sf(data = ntd_amt, aes(fill = DETOUTCOME))

# Generate protected area raster that retains NT determination categories from original shapefile
ntd_type <- unique(ntd_amt$DETOUTCOME)
ntd_df <- data.frame(id = c(300, 400), type = ntd_type) # Give id numbers different to PA and pastoral ids
ntd_amt$id <- ntd_df$id[match(ntd_amt$DETOUTCOME, ntd_df$type)] 

# Rasterize using modal function to find the most frequent value (mode) for each 10 x 10 km grid cell
ntd_vect <- vect(ntd_amt)
ntd_ras <- rasterize(ntd_amt, r_amt, field = "id", fun = 'max')
levels(ntd_ras) <- ntd_df
plot(ntd_ras)

# Assign land type category
ntd_df$category <- 'Native Title land'

#### INDIGENOUS LAND USE AGREEMENTS ####

# Load shapefile of Indigenous Land Use Agreements
ilua <- read_sf("data/input-data/nntt/ILUA_Registered_Notified_Nat.shp") %>%
  st_transform(crs = st_crs(3577)) 

# Intersect ILUA polygons with the AMT
ilua_amt <- st_intersection(ilua, amt) %>%
  st_collection_extract("POLYGON")

# Plot by Subject Matter
ggplot() +
  geom_sf(data = amt) +
  geom_sf(data = ilua_amt, aes(fill = SUBJMATTER))

# Generate protected area raster with categories from original shapefile
ilua_type <- unique(ilua_amt$SUBJMATTER)
ilua_df <- data.frame(id = c(101:269), type = ilua_type) # Give id numbers different to PA and pastoral ids
ilua_amt$id <- ilua_df$id[match(ilua_amt$SUBJMATTER, ilua_df$type)] 

# Rasterize using modal function to find the most frequent value (mode) for each 10 x 10 km grid cell
ilua_vect <- vect(ilua_amt)
ilua_ras <- rasterize(ilua_vect, r_amt, field = "id", fun = 'max')
levels(ilua_ras) <- ilua_df
plot(ilua_ras)

# Assign broad land type categories
ilua_df$category <- ifelse(grepl("pastoral", ilua_df$type, ignore.case = TRUE),
                           "Indigenous Land Use Agreement - Pastoral",
                           "Indigenous Land Use Agreement")

#### DESIGNATE PROTECTED STATUS OF PUs ####

# Calculate the modal CAPAD protected area type for raster values within each planning unit polygon
amt_grid$capad <- terra::extract(pa_ras, amt_grid, method = "bilinear", fun = modal, na.rm = TRUE)[[2]]
# Check output looks correct
ggplot() + 
  geom_sf(data = amt_grid,
          aes(fill = capad))

# Create a binary variable indicating if a PU is primarily protected (1) or not (0)
amt_grid <- amt_grid %>%
  mutate(status = ifelse(!is.na(capad) & capad > 0, 1, 0))
# Check output looks correct
ggplot() + 
  geom_sf(data = amt_grid,
          aes(fill = status))

# Generate logical (T/F) variable to lock in protected areas since they are already in the reserve system
amt_grid$locked_in <- as.logical(amt_grid$status)

#### MAP ALL LAND TENURE (PROTECTED, NATIVE TITLE, PASTORAL, and OTHER) IN PUs ####

# Create a list of all raster layers relating to land tenure
landcat_rasters <- list(pa_ras, ntd_ras, ilua_ras, tenure)

# Combine the list of rasters into a single raster stack
landcat_stack <- rast(landcat_rasters)

# Extract cell values from each raster layer in landcat_stack for each planning unit polygon
pu_landcat <- terra::extract(landcat_stack, amt_grid, method = 'bilinear', fun = modal, na.rm = T) 
colnames(pu_landcat) <- c("id", "capad", "ntd", "ilua", "tenure") # Rename the columns

# Add a new column to calculate dominant tenure
pu_landcat$main_id <- NA

# Iterate through each row according to the following rules to populate the 'main_cat' column
  # Protected area status takes precedence over all other land categories
  # Native Title determination takes precedence over ILUA and other tenure status
  # ILUA takes precedence over other tenure status
  # If the PU is not protected, Native Title, or under an ILUA, then attribute the other 'tenure' value
for (i in 1:nrow(pu_landcat)) {
  if (!is.na(pu_landcat$capad[i]) && pu_landcat$capad[i] > 0) {
    pu_landcat$main_id[i] <- pu_landcat$capad[i]
  } else if (!is.na(pu_landcat$ntd[i]) && pu_landcat$ntd[i] > 0) {
    pu_landcat$main_id[i] <- pu_landcat$ntd[i]
  } else if (!is.na(pu_landcat$ilua[i]) && pu_landcat$ilua[i] > 0) {
    pu_landcat$main_id[i] <- pu_landcat$ilua[i]
  } else {
    pu_landcat$main_id[i] <- pu_landcat$tenure[i]
  }
}

# Merge all id and land type data frames for different types of tenure
landcat_df <- rbind(pa_df, ntd_df, ilua_df, tenure_df)

# Add new columns to calculate the predominant land tenure type and broad category  of cells within each PU
pu_landcat$main_type <- NA
pu_landcat$category <- NA

# Iterate through each row to assign land tenure type to corresponding id number
for (i in 1:nrow(pu_landcat)) {
  matching_row <- landcat_df[landcat_df$id == pu_landcat$main_id[i],]
  if (nrow(matching_row) > 0) {
    pu_landcat$main_type[i] <- matching_row$type
    pu_landcat$category[i] <- matching_row$category
  }
}

# Check if PUs are classified as pastoral term or perpetual lease in DAWE data
pastoral <- pu_landcat$tenure %in% c(2121, 2141, 2142)

# Update main_type and category based on condition if it is both Native Title land and under pastoral use
pu_landcat <- pu_landcat %>%
  mutate(
    main_type = case_when(
      (main_id == 300 & pastoral) | (main_id == 400 & pastoral) ~ "Native Title land - Pastoral use",
      TRUE ~ main_type),
    category = case_when(
      (main_id == 300 & pastoral) | (main_id == 400 & pastoral) ~ "Native Title land - Pastoral use",
      category == "Indigenous Land Use Agreement" & pastoral ~ "Indigenous Land Use Agreement - Pastoral use",
      TRUE ~ category)
    )

# Bind columns 6 to 8 of pu_landcat to amt_grid
amt_grid_v2 <- amt_grid %>%
  bind_cols(pu_landcat %>% dplyr::select(6:8)) %>%
  filter(if_any(8:10, ~ !is.na(.))) # Remove 20 rows with no land classification

# Check it looks correct
ggplot() + 
  geom_sf(data = amt_grid_v2,
          aes(fill = category)) 

# Convert 'land_type' column to a factor with levels in the order specified below
amt_grid_v2$category <- factor(amt_grid_v2$category, levels = c("Protected Area", 
                                                                      "Indigenous Protected Area", 
                                                                      "Freehold - Indigenous", 
                                                                      "Native Title land", 
                                                                      "Indigenous Land Use Agreement", 
                                                                      "Native Title land - Pastoral use",
                                                                      "Indigenous Land Use Agreement - Pastoral use",
                                                                      "Pastoral term or perpetual lease", 
                                                                      "Other freehold, term, perpetual lease or Crown purposes"))

# Plot by land category
ggplot() +
  geom_sf(data = amt_grid_v2, aes(fill = category), color = NA) +
  scale_fill_manual(values = c("#003c30",
                               "#01665e", 
                               "#35978f", 
                               "#80cdc1", 
                               "#c7eae5", 
                               "#f6e8c3", 
                               "#dfc27d", 
                               "#bf812d", 
                               "#8c510a")) +
  geom_sf(data = amt, fill = "transparent", color = "gray20") +
  theme_minimal()

# Create a new binary column indicating whether the PU can be voluntarily declared as an IPA (i.e., is native title land) and is not already protected
amt_grid_v2 <- amt_grid_v2 %>%
  mutate(declare_ipa = if_else(main_id == 300 | main_id == 400, 1, 0))

#### MAP MEDIAN EMPLOYEE INCOME FOR LGAs to PU LAYER ####

# Read, transform, and crop LGA shapefiles for QLD, NT, and WA
# Queensland
lga_qld <- read_sf("data/input-data/lgas/qld/qld_lga.shp") %>%
  st_transform(crs = st_crs(3577)) %>%
  st_intersection(amt)
# Northern Territory
lga_nt <- read_sf("data/input-data/lgas/nt/nt_lga.shp") %>%
  st_transform(crs = st_crs(3577)) %>%
  st_buffer(1500) %>% # buffer by 1.5km to remove Victoria River which is not featured on 'amt' polygon and results in some PUs with NA values
  st_intersection(amt)
# Western Australia
lga_wa <- read_sf("data/input-data/lgas/wa/wa_lga.shp") %>%
  st_transform(crs = st_crs(3577)) %>%
  st_intersection(amt)

# Merge the LGA shapefiles for QLD, NT, and WA
lga_amt <- rbind(lga_qld, lga_nt, lga_wa) 
# Plot LGAs
ggplot(data = lga_amt) + 
  geom_sf() +
  geom_sf(
    data = amt,
    fill = "transparent",
    color = "red",
    size = 0.5)

# Read in Median Employee Income (MEI) 2011-2018 data (ABS, 2018)
med_income <- read.csv("data/input-data/lgas/LGA_med_income.csv")

# Merge income data with LGA polygons
lga_amt <- lga_amt %>%
  left_join(med_income[,2:3], by = 'ABB_NAME')

# Plot MED_INCOME 
ggplot(lga_amt) +
  geom_sf(aes(fill = MED_INCOME), color = 'transparent') +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(fill = "Median Income", title = "Median Income by LGA")

# Rasterize & specify 'mean' function to calculate average median employee income value 
lga_amt_vect <- vect(lga_amt)
lga_amt_r <- rasterize(lga_amt_vect, r_amt, field = "MED_INCOME", fun = 'mean')
plot(lga_amt_r)

# Extract median employee income value for each PU
# Specify 'bilinear' method and 'modal' function to interpolate modal value from the four nearest cells
pu_medincome <- terra::extract(lga_amt_r, amt_grid_v2, method = 'bilinear', fun = modal, na.rm = T)

# Merge with existing PU data in new data frame
amt_grid_v3 <- amt_grid_v2 %>%
  bind_cols(pu_medincome %>% dplyr::select(2))

# Plot median employee income by LGA
ggplot() +
  geom_sf(data = amt) +
  geom_sf(data = amt_grid_v3, aes(fill = MED_INCOME)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(fill = "Median Employee Income", title = "Median Employee Income by LGA")

#### GENERATE COST PROXY LAYER ####

# Load Human Influence Index data for Oceania
hii <- rast("data/input-data/hii/hii-oceania-geo-grid/hii_ocean_grid/hii_oceania/w001001.adf")
hii <- project(hii, r_amt) %>%
  mask(r_amt)
plot(hii)

# Extract maximum Human Influence Index (HII) values for each PU
# Broadacre land values shown to increase in proximity to infrastructure such as roads, railways, & urban centres which increase profitability
# (Chancellor et al 2019 - Measuring Australian broadacre farmland value)
### Number of buildings on land positively associated with land value
### Transport cost to nearest city negtively associated with land value
pu_hii <- terra::extract(hii, amt_grid_v3, method = 'bilinear', fun = max, na.rm = T)

# Merge with existing PU data in new data frame
amt_grid_v4 <- amt_grid_v3 %>%
  bind_cols(pu_hii %>% dplyr::select(2))

# Plot HII for PUs
ggplot() +
  geom_sf(data = amt) +
  geom_sf(data = amt_grid_v4, aes(fill = hii_oceania)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(fill = "Human Influence Index (HII)")

# Create function to normalise data with NA handling
normalise <- function(x) {
  if(all(is.na(x))) {
    return(rep(NA, length(x))) # Return NA if all values are NA
  }
  range_x <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if(range_x == 0) {
    return(rep(NA, length(x))) 
  }
  return((x - min(x, na.rm = TRUE)) / range_x)
}

# Normalise both HII and MEI to a 0-1 range using min-max scaling
amt_grid_v4$normal_hii <- normalise(amt_grid_v4$hii_oceania)
amt_grid_v4$normal_mei <- normalise(amt_grid_v4$MED_INCOME)

# Compute cost proxy values for each PU giving weightings to HII and MEI
# α and β are weights you can adjust to determine the relative influence of HII and MEI on the cost (these should sum to 1)
alpha <- 0.6
beta <- 0.4

amt_grid_v4$cost <- alpha*amt_grid_v4$normal_hii + beta*amt_grid_v4$normal_mei

# Boxplot of range of cost values for each land category
ggplot(amt_grid_v4, aes(x = category, y = cost, fill = category)) +
  geom_boxplot() +
  labs(x = "Land Category", y = "Cost Proxy Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

# Adjust cost proxy to decrease the cost proxy for PUs on native title land or inalienable Indigenous freehold. This factor should be less than 1 to reduce the cost
# Multiple native title/Indigenous freehold land by 0.3 to remove overlap in interquartile range of cost proxy values between land that is and is not native title
ntd_adjust <- 0.3
amt_grid_v4$cost_adjusted <- ifelse(amt_grid_v4$declare_ipa == 1 | amt_grid_v4$category == "Freehold - Indigenous",
                                     amt_grid_v4$cost*ntd_adjust, 
                                     amt_grid_v4$cost)
hist(amt_grid_v4$cost_adjusted)

# Square root transformation to remove skew
amt_grid_v4$cost_adjusted <- sqrt(amt_grid_v4$cost_adjusted)
hist(amt_grid_v4$cost_adjusted)

# Boxplot of adjusted proxy cost values for land where declare_ipa == 1 and declare_ipa == 0
# No overlap in IQR between native title land and other tenure
ggplot(amt_grid_v4, aes(x = category, y = cost_adjusted, fill = category)) +
  geom_boxplot() +
  labs(x = "Tenure", y = "Cost Proxy Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

# Plot cost proxy values for PUs 
ggplot() +
  geom_sf(data = amt_grid_v4, aes(fill = cost_adjusted), color = NA) +
  scale_fill_viridis_c() +
  geom_sf(data = capad_amt, fill = "transparent", color = "red") + # to plot PA boundaries
  geom_sf(data = amt, fill = "transparent", color = "gray20") +
  theme_minimal() +
  labs(fill = "Cost proxy values")

#### CREATE SPEC_DAT DATA FRAME ####

## SPEC_DAT = data frame of conservation feature identifiers (i.e., the species to be protected)
## id
## name
## representation target

# Load shapefile of IUCN distribution polygons for AMT mammals with threat status data generated in 02_pgls_model
# Note: iucn_threat contains geometry for entire species distribution
iucn_threat <-  read_sf("data/output-data/shp/02_iucn_threat.shp")
colnames(iucn_threat)[1:4] <- c("scientificName", "ordinalThreat", "mean_fitted", "latent_risk")
# Create vector of AMT species names
amt_names <- iucn_threat$scientificName

# Create raster stack of 141 species' distribution polygons including unique identifier and species scientific name
# Initialize an empty list to store the raster layers
spec_rasters <- list()
for (i in 1:length(iucn_threat$scientificName)) {
  curr_vect <- vect(iucn_threat[i, "geometry"]) # Convert the current geometry to a terra vector object
  curr_rast <- rasterize(curr_vect, r_amt, background = NA) # Rasterize the vector object using the template raster and threat vals
  # Associate species name with corresponding raster 
  names(curr_rast) <- paste(iucn_threat$scientificName[i])
  spec_rasters[[i]] <- curr_rast # Store the raster in the list
}
# Combine the list of rasters into a single raster stack
spec_stack <- rast(spec_rasters)
# Export the raster stack as a TIFF file
writeRaster(spec_stack, filename="data/output-data/tif/03_species_stack.tif", overwrite=TRUE)

# Create spec_dat table with feature identifiers, names, and representation targets (prop) for each
spec_dat = data.frame(id = seq_len(nlyr(spec_stack)),
                      scientificName = amt_names,
                      prop = NA) 

# Crop iucn_threat sf data frame to AMT boundary
iucn_amt <- st_intersection(iucn_threat, amt) %>%
  dplyr::select(!5)

# Merge with iucn_amt sf data frame
spec_dat <- merge(spec_dat, iucn_amt, by = 'scientificName')

# Calculate total species distribution area 
spec_dat$area_total <- NA
for (i in 1:length(iucn_threat$scientificName)) {
  curr_shape <- iucn_threat[i, "geometry"]
  curr_area <- st_area(curr_shape)
  spec_dat$area_total[i] <- curr_area
}

# Calculate area of species distribution that intersects with the AMT
spec_dat$area_amt <- NA
for (i in 1:length(iucn_amt$scientificName)) {
  curr_shape <- iucn_amt[i, "geometry"]
  curr_area <- st_area(curr_shape)
  spec_dat$area_amt[i] <- curr_area
}

# Calculate the proportion of each species' distributions that overlap with the AMT
spec_dat <- spec_dat %>%
  mutate(prop = pmin(area_amt/area_total, 1))

# Write to csv
write_csv(spec_dat, "data/output-data/tbl/03_spec_dat.csv", append = F)

#### CALCULATE SPECIES REPRESENTATION TARGETS ####

# TARGET 1: Calculate the CURRENT RISK representation target
# Exclude species with ordinal threat status < 4 (i.e., exclude Near Threatened or Least Concern species)
spec_dat <- spec_dat %>%
  mutate(threatenedSpecies = if_else(ordinalThreat > 3, ordinalThreat, as.numeric(NA)))
# Fifteen species in this category

# Use min-max scaling to normalise current risk scores to between 0 - 1
spec_dat$threatenedSpecies_scaled <- normalise(spec_dat$threatenedSpecies)

# Assign weights to current risk and the proportion of species distribution in the AMT to calculate representation target
# Note: the representation target is the proportion of the species distribution that the prioritizr algorithm aims to secure
alpha <- 0.7 # higher weight for threatened species
beta <- 0.3 # include lower weight for proportion of distribution in AMT to prioritise protection of threatened species that are both endemic to the AMT and have small distributions
spec_dat$target_1 <- (alpha*spec_dat$threatenedSpecies_scaled) + (beta*spec_dat$prop)
# Convert NA values to zeroes
spec_dat <- spec_dat %>%
  mutate(
    target_1 = if_else(target_1 > 0.9, 0.9, target_1),
    target_1 = replace_na(target_1, 0))


# TARGET 2: Calculate the POSITIVE LATENT RISK protection target for species
# Calculate the 90th percentile of latent risk values
stats::quantile(spec_dat$latent_risk, probs = 0.9, na.rm = T) # 1.40

# Create new column with only those species in the top 90th percentile of positive latent risk values 
spec_dat <- spec_dat %>%
  mutate(positiveLatentRisk = if_else(latent_risk > 1.403645, latent_risk, as.numeric(NA)))
# Fifteen species in this category, again

# Use min-max scaling to normalise latent risk scores to between 0 - 1
spec_dat$latentRisk_scaled <- normalise(spec_dat$positiveLatentRisk)

# Assign weights to positive latent risk and proportion of species distribution in the AMT to calculate representation target
# (the representation target is the proportion of the species distribution area that the prioritizr algorithm aims to secure)
alpha <- 0.7 # higher weight for species with a higher positive latent risk value
beta <- 0.3 # include lower weight for proportion of distribution in AMT to prioritise protection of threatened species that are both endemic to the AMT and have small distributions
spec_dat$target_2 <- (alpha*spec_dat$latentRisk_scaled) + (beta*spec_dat$prop) 
# Convert NA values to zeroes and make any values > 0.9 equal 0.9 to ensure algorithm can run
spec_dat <- spec_dat %>%
  mutate(
    target_2 = if_else(target_2 > 0.9, 0.9, target_2),
    target_2 = replace_na(target_2, 0))

#### CALCULATE SPECIES RICHNESS TARGET ####

# Map average and net species richness across the AMT (only including final subset of 141 AMT mammals)
# Initialize an empty list to store the raster layers
richness_rasters <- list()
for (i in 1:length(iucn_threat$scientificName)) {
  curr_vect <- vect(iucn_threat[i, "geometry"]) 
  curr_rast <- rasterize(curr_vect, r_amt) 
  # Associate species name with corresponding raster 
  names(curr_rast) <- paste(iucn_threat$scientificName[i])
  values(curr_rast)[which(values(curr_rast) > 0)] <- 1
  values(curr_rast)[which(is.na(values(curr_rast)))] <- 0
  richness_rasters[[i]] <- curr_rast 
}
# Combine the list of rasters into a single raster stack
richness_stack <- rast(richness_rasters)

# Calculate the net current threat status, excluding NA values
sum_richness <- sum(richness_stack, na.rm = TRUE) %>%
  mask(r_amt)
plot(sum_richness)

# Extract species richness for PUs
pu_richness <- extract(sum_richness, amt_grid_v4, method = 'bilinear', fun = modal, na.rm = T)

# Merge with existing PU data in new data frame
amt_grid_v5 <- cbind(amt_grid_v4, pu_richness)
colnames(amt_grid_v5)[18] <- 'spec_richness'

# Write final PU grid to shapefile
st_write(amt_grid_v5, "data/output-data/shp/03_amt_pu_grid.shp", append = F)

# Plot species richness
ggplot() +
  geom_sf(data = amt_grid_v5, aes(fill = spec_richness)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(fill = "Species richness")


## TARGET 3: prioritise PUs in the top 90th percentile of species richness (41-62 species)
# Calculate quantiles for species richness values across PUs 
quantile(pu_richness[,2], probs = 0.9, na.rm = T)
# 90th percentile = 41 species (i.e., only 10% of PUs have greater than 41 species)

# Define new breaks for species richness
richness_breaks <- c(0, 40, 62)

# Adjusted representation targets based on the two groups
# I.e., if PU has < 36 species, do not prioritise. If > 41 species, must be included in optimal reserve network
richness_targets <- c(0, 0.8)

# Initialize an empty list to store the individual raster layers for the two groups
richness_layers <- list()
# Loop through the breaks to create new raster layers
for (i in 1:(length(richness_breaks) - 1)) {
  lower_bound <- richness_breaks[i]
  upper_bound <- richness_breaks[i + 1]
  # Create a new raster layer for each range, setting cells within the range to 1 and others to 0
  curr_layer <- app(sum_richness, fun = function(x) {
    ifelse(x >= lower_bound & x < upper_bound, 1, 0)
  })
  # Add the current raster layer to the list
  richness_layers[[i]] <- curr_layer
}

# Stack the raster layers into a single raster stack
richness_lyr_stack <- rast(richness_layers)
# Visualise areas in each richness category
plot(richness_lyr_stack[[1]]) # low
plot(richness_lyr_stack[[2]]) # high

# Adjust the layer names to reflect the new groups
rich_group <- c("rich_0_40", "rich_41_62")

# Assign the new names to the layers in the stack
names(richness_lyr_stack) <- rich_group
# Export the raster as a TIFF file
writeRaster(richness_lyr_stack, filename="data/output-data/tif/03_richness_stack.tif", overwrite=TRUE)

# Create a data frame with layer names and the corresponding representation targets
richness_df <- data.frame(rich_group, target_3 = richness_targets)
# Write to csv
write_csv(richness_df, "data/output-data/tbl/03_richness_df.csv", append = F)





