######################################################
# AMT Spatial Conservation Prioritisation for Mammals
#
# Part 2: Phylogenetic Least-Squares Models
#
# Code developed by Emmeline Norris, 01/06/2024
# 
######################################################
install.packages('patchwork')
library(ape)
library(caper)
library(MuMIn)
library(lmtest)
library(rr2)
library(purrr)
library(sf)
library(terra)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Disable scientific notation
options(scipen = 999)

#### CREATE COMPARATIVE DATA OBJECT FOR PGLS MODELS ####

# Load data frame of response and predictor variables generated in 01_load_and_clean.R
amt_dat_st <- read.csv("data/output-data/tbl/01_amt_dat_st.csv")
# Load list containing 100 phylogenetic trees pruned to AMT mammals generated in 01_load_and_clean.R
amt_phy <- readRDS("data/output-data/phy/01_amt_phy.rds")

# Create a "comparative.data" object that puts the phylogenetic tree and predictor df together in one object and matches the species names
# It also removes any species with missing values for any of the variables so that the df is 100% complete
i <- 1
compdat <- comparative.data(phy = amt_phy[[i]], data = amt_dat_st, names.col="scientificName", na.omit = T) # na.omit = F prevents it from deleting species w missing val
# The number of tip labels in the tree and the number of observations in the data frame should be the same
str(compdat)

#### BUILD PHYLOGENETIC LEAST-SQUARES MODELS ####

# Create getMatrices function generates the matrices needed by the fitspphylo function
# Input is a comparative.data object which has a phylogeny ($phy) and a dataset ($data) with identical species sets
getMatrices <- function(compdat, spmat=NULL, phylomat=NULL){
  if(is.null(spmat))spmat = as.matrix(dist(compdat$data[,c(1:2)], diag=T, upper=F))
  if(is.null(phylomat))phylomat = vcv.phylo(compdat$phy, corr=TRUE)
  return(list(spmat = spmat, phylomat = phylomat))
}

# The fitspphylo function was written by Xia Hua to fit gls models with both a phylo and spatial matrix 
# Input arguments are a formula (see below), data frame, spatial distance matrix, phylo distance matrix, p-matrix (see below)
fitspphylo <- function (formula, data, spmatrix, phylomatrix, p) {
  cal <- function (p, formula, data, spmatrix, phylomatrix) {
    spmatrix <- spmatrix/max(spmatrix)
    spmatrix <- exp(-(spmatrix/p[2])^2)
    mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
    res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
    if (inherits(res,"try-error")) {
      out <- -10000
    } else {
      out <- res$logLik
      if (res$logLik>0) {out <- -10000}
    }
    -out
  }
  cal2 <- function (p,formula,data,spmatrix,phylomatrix) {
    spmatrix <- spmatrix/max(spmatrix)
    spmatrix <- exp(-(spmatrix/p[2])^2)
    mat <- as.matrix((p[3]*(1-p[1])*spmatrix+(1-p[1])*(1-p[3])*phylomatrix+p[1]*diag(dim(phylomatrix)[1])))
    res <- try(gls(model=formula,data=data,correlation=corSymm(mat[lower.tri(mat)],fixed=T),method="ML"),silent=T)
    res
  }
  library(nloptr)
  library(nlme)
  p.res <- sbplx(p,cal,formula=formula,data=data,spmatrix=spmatrix,phylomatrix=phylomatrix,lower=c(0.01,0,0),upper=c(1,1,1))
  lm.res <- cal2(p=p.res$par,formula=formula,data=data,spmatrix=spmatrix,phylomatrix=phylomatrix)
  list(p.res=p.res,lm.res=lm.res)
}

# Generate spatial distance matrix among species centroids (call amt_dat_st columns containing centroid lat & lon)
spmat <- as.matrix(dist(amt_dat_st[,c(3,4)], diag=T, upper=F))

# Generate phylogenetic distance matrix from compdat object
phylomat <- getMatrices(compdat)$phylomat

# Starting values for p 
# The 3 values are the overall contribution of autocorrelation, decay rate of spatial correlation, % spatial vs phylo
p <- c(0.5, 0.5, 0.5)


#### EXECUTE FULL AND REDUCED PGLS MODELS ####

## FULL MODEL ##

# Set up a formula object to specify a global (12 predictor variables) 
pred_full <- c("adult_mass_g_tr_st",
                     "I(adult_mass_g_tr_st^2)",
                     "litter_size_n_tr_st",
                     "litters_per_year_n_tr_st",
                     "age_first_reproduction_d_tr_st",
                     "foxAreaProp_tr_st",
                     "toadAreaProp_tr_st",
                     "mean_bio8_st",
                     "mean_firefreq_st",
                     "mean_latefire_st",
                     "mean_hii_tr_st",
                     "rangeArea_km2_tr_st")

full_model <- paste("ordinalThreat ~", paste(pred_full, collapse="+"))
fmla <- as.formula(full_model)
full_model <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
full_model$p.res # the maximum likelihood estimates of the 3 parameters in p
summary(full_model$lm.res)  # the model summary with coefficient estimates, etc.

# Model diagnostics
qqnorm(full_model$lm.res) # Residuals appear normally distributed in centre with slightly heavier tails than expected for normal dist
plot(full_model$lm.res) # Plot indicates heteroscedasticity
AICc(full_model$lm.res) # 486.8173


## REDUCED MODELS ##

# 1
# Remove predictors based on highest p-values to see if it affects model fit
# Remove Bio8 precipitation variable (p = 0.8692)
pred_red_1 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "litter_size_n_tr_st",
                "litters_per_year_n_tr_st",
                "age_first_reproduction_d_tr_st",
                "foxAreaProp_tr_st",
                "toadAreaProp_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st",
                "mean_hii_tr_st",
                "rangeArea_km2_tr_st")

red_model_1 <- paste("ordinalThreat ~", paste(pred_red_1, collapse="+"))
fmla <- as.formula(red_model_1)
red_model_1 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_1$p.res # the maximum likelihood estimates of the 3 parameters in p
summary(red_model_1$lm.res)

# Model diagnostics
qqnorm(red_model_1$lm.res)
plot(red_model_1$lm.res)
AICc(red_model_1$lm.res) # 484.3238


# 2
# Remove litter size variable (P = 0.5464)
pred_red_2 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "litters_per_year_n_tr_st",
                "age_first_reproduction_d_tr_st",
                "foxAreaProp_tr_st",
                "toadAreaProp_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st",
                "mean_hii_tr_st",
                "rangeArea_km2_tr_st")

red_model_2 <- paste("ordinalThreat ~", paste(pred_red_2, collapse="+"))
fmla <- as.formula(red_model_2)
red_model_2 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_2$p.res 
summary(red_model_2$lm.res)

# Model diagnostics
qqnorm(red_model_2$lm.res)
plot(red_model_2$lm.res)
AICc(red_model_2$lm.res) # 482.2471


# 3
# Remove HII variable (p = 0.6008)
pred_red_3 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "litters_per_year_n_tr_st",
                "age_first_reproduction_d_tr_st",
                "foxAreaProp_tr_st",
                "toadAreaProp_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st",
                "rangeArea_km2_tr_st")
red_model_3 <- paste("ordinalThreat ~", paste(pred_red_3, collapse="+"))
fmla <- as.formula(red_model_3)
red_model_3 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_3$p.res # the maximum likelihood estimates of the 3 parameters in p
summary(red_model_3$lm.res)

# Model diagnostics
qqnorm(red_model_3$lm.res)
plot(red_model_3$lm.res)
AICc(red_model_3$lm.res) # 480.1106


# 4
# Remove litters per year variable (p = 0.1276)
pred_red_4 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "age_first_reproduction_d_tr_st",
                "foxAreaProp_tr_st",
                "toadAreaProp_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st",
                "rangeArea_km2_tr_st")

red_model_4 <- paste("ordinalThreat ~", paste(pred_red_4, collapse="+"))
fmla <- as.formula(red_model_4)
red_model_4 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_4$p.res # the maximum likelihood estimates of the 3 parameters in p
summary(red_model_4$lm.res)

# Model diagnostics
qqnorm(red_model_4$lm.res)
plot(red_model_4$lm.res)
AICc(red_model_4$lm.res) # 480.2436 (slight increase)


# 5
# Remove toad area variable (p = 0.1431)
pred_red_5 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "age_first_reproduction_d_tr_st",
                "foxAreaProp_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st",
                "rangeArea_km2_tr_st")

red_model_5 <- paste("ordinalThreat ~", paste(pred_red_5, collapse="+"))
fmla <- as.formula(red_model_5)
red_model_5 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_5$p.res # the maximum likelihood estimates of the 3 parameters in p
summary(red_model_5$lm.res)

# Model diagnostics
qqnorm(red_model_5$lm.res)
plot(red_model_5$lm.res)
AICc(red_model_5$lm.res) # 480.2049 (still higher than model 3)

# 6
# Remove fox area variable (p = 0.1279)
pred_red_6 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "age_first_reproduction_d_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st",
                "rangeArea_km2_tr_st")

red_model_6 <- paste("ordinalThreat ~", paste(pred_red_6, collapse="+"))
fmla <- as.formula(red_model_6)
red_model_6 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_6$p.res # the maximum likelihood estimates of the 3 parameters in p
summary(red_model_6$lm.res)

# Model diagnostics
qqnorm(red_model_6$lm.res)
plot(red_model_6$lm.res)
AICc(red_model_6$lm.res) # 480.3724 (still higher than model 3)

# 7
# Remove range area variable (p = 0.1016)
pred_red_7 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "age_first_reproduction_d_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st")

red_model_7 <- paste("ordinalThreat ~", paste(pred_red_7, collapse="+"))
fmla <- as.formula(red_model_7)
red_model_7 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_7$p.res # the maximum likelihood estimates of the 3 parameters in p
summary(red_model_7$lm.res)

# Model diagnostics
qqnorm(red_model_7$lm.res)
plot(red_model_7$lm.res)
AICc(red_model_7$lm.res) # 480.9421 (still higher than model 3)


## Reduced model #3 has lowest AICc value by a very small amount
# Use model 3 

# Construct a null model for comparison
fmla_null <- as.formula('ordinalThreat ~ 1') # intercept-only model without predictors
model_null <- fitspphylo(formula=fmla_null, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)

# Likelihood ratio test for full and reduced models
lmtest::lrtest(full_model$lm.res, red_model_3$lm.res)
# P-value of 0.8651 indicates that the additional vars in full model do not improve model fit

# Likelihood ratio test for reduced models #5 and null model
lmtest::lrtest(red_model_3$lm.res, model_null$lm.res)
# P-value of < 0.0005 indicates the reduced model is significantly better at explaining variability in data compared to null model

# Calculate R^2
r2 <- R2_lik(mod = red_model_3$lm.res, mod.r = model_null$lm.res)
r2*100 # r^2 = 32.80%

# Calculate delta AICc for full and minimum adequate model 
d_AICc <- AICc(full_model$lm.res) - AICc(red_model_3$lm.res)
d_AICc # 6.71

#### PREDICT LATENT EXTINCTION RISK FOR AMT MAMMALS ####

# Load data frame containing AMT species names and ordinal threat status, generated in 01_load_and_clean
amt_species_df <- read.csv("data/output-data/tbl/01_amt_species_df.csv")
# Replace space between genus and species with underscore
amt_species_df$scientificName <- gsub(" ", "_", amt_species_df$scientificName)

# Run best fitting reduced model with 100 phylogenies
for (i in 1:100){ 
  compdat <- comparative.data(phy=amt_phy[[1]], data= amt_dat_st, names.col="scientificName")
  red_model_3 <- paste("ordinalThreat ~", paste(pred_red_3, collapse="+"))
  fmla <- as.formula(red_model_3)
  red_model_3 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
  saveRDS(red_model_3$lm.res, file=paste0("data/output-data/mod/02_model", i, ".rds"))
}

# Create a list of file paths to load the 100 models
model_list <- list.files(path = "data/output-data/mod", pattern = "*.rds", full.names = TRUE )
all_models <- purrr::map(model_list, readRDS)

# Extract 'fitted' element from all models in list
all_fitted <- lapply(all_models, fitted)
# Create new data frame to output fitted extinction risk values, ordinal threat, latent threat 
all_threat <- amt_species_df 

# Iterate through each object in the all_fitted list
for (i in 1:length(all_fitted)) {
  # Create a new data frame for each object
  model_df <- data.frame(scientificName = names(all_fitted[[i]]), fitted_risk = unlist(all_fitted[[i]]))
  # Generate the data frame name dynamically
  model_name <- paste0("model_", i)
  # Assign the data frame to the generated name
  assign(model_name, model_df)
}

# Iterate through each model data frame and merge the 'fitted_risk' column
for (i in 1:100) {
  model_name <- paste0("model_", i)
  col_name <- paste0("fitted_", i)
  all_threat <- merge(all_threat, get(model_name), by = "scientificName", all = TRUE)
  colnames(all_threat)[colnames(all_threat) == "fitted_risk"] <- col_name
}

# Calculate the mean fitted risk for each row
all_threat$mean_fitted <- rowMeans(all_threat[, 3:102])
# Calculate the latent risk for each row
all_threat$latent_risk <- all_threat$mean_fitted - all_threat$ordinalThreat
# Replace underscores with spaces between genus and species
all_threat$scientificName <- gsub("_", " ", all_threat$scientificName)
# Filter df to include only important columns
all_threat <- all_threat %>%
  dplyr::select(1,2,103,104)

#### PRODUCE TABLE OF ALL THREAT VALUES AND RANGE SIZE ####

# Load IUCN species summary csv (includes taxonomic details and threat category)
iucn_summary = read_csv("data/input-data/iucn/Summary/simple_summary.csv")

# Merge all_threat df with Red List category and population trend columns of summary table
all_threat_df <- all_threat %>%
  left_join(iucn_summary %>% 
              dplyr::select(scientificName, redlistCategory, populationTrend), by = "scientificName")

#Reorder columns
all_threat_df <- all_threat_df %>%
  dplyr::select(1, 5, 6, 2, 3, 4, everything())

# Write to csv
write_csv(all_threat_df, "data/output-data/tbl/02_all_threat_df.csv", append = F)

# Create new df also containing untransformed predictor variables for each species range
species_info <- all_threat_df %>%
  left_join(amt_dat %>% 
              dplyr::select(scientificName, 24,7,9:13,27,30,33,35,38,41,44), by = "scientificName")

# Write to csv
write_csv(species_info, "data/output-data/tbl/02_species_covariates.csv", append = F)

#### CREATE MAPS OF SUM AND MEAN LATENT RISK FROM MODELS ####

# Load shapefile of IUCN distribution polygons for mammals in the AMT generated in 01_load_and_clean
iucn_union_mod <-  read_sf("data/output-data/shp/01_iucn_union_mod.shp")
colnames(iucn_union_mod)[1] <- "scientificName"

# Merge summary data with IUCN species distribution spatial dataset
iucn_threat <- inner_join(iucn_union_mod, all_threat, by = 'scientificName') %>%
  dplyr::select(!2)
# Write to shapefile for use in 03_scp_data_prep.R
st_write(iucn_threat, "data/output-data/shp/02_iucn_threat.shp", append = F)

# Load amt raster from 01_load_and_clean
r_amt <- rast("data/output-data/tif/01_r_amt.tif")

# Initialize an empty list to store the raster layers
latent_risk_rasters <- list()
for (i in 1:length(iucn_threat$scientificName)) {
  curr_latent_risk <- iucn_threat$latent_risk[i]
  curr_vect <- vect(iucn_threat[i, "geometry"]) # Convert the current geometry to a terra vector object
  curr_rast <- rasterize(curr_vect, r_amt, field = curr_latent_risk, background = NA) # Rasterize the vector object using the template raster and threat vals
  latent_risk_rasters[[i]] <- curr_rast # Store the raster in the list
}
# Combine the list of rasters into a single raster stack
latent_stack <- rast(latent_risk_rasters)

# Calculate the net latent extinction risk for the AMT, excluding NA values
sum_latent_risk <- sum(latent_stack, na.rm = TRUE) %>%
  mask(r_amt)

# Reproject to WGS84 (EPSG:4326) geographic coordinate system for plotting
sum_latent_risk <- project(sum_latent_risk, "EPSG:4326")

# Convert the raster data to a data frame
sum_latent_rast_df <- as.data.frame(sum_latent_risk, xy = TRUE)

# Plot the raster data using ggplot2
l1 <- ggplot(sum_latent_rast_df, aes(x = x, y = y, fill = sum)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D",
                       breaks = seq(-10, 20, by = 10), 
                       labels = seq(-10, 20, by = 10)) +
  labs(fill = "Sum") +
  theme_void() +
  coord_fixed(ratio = 1 / cos(mean(sum_latent_rast_df$y) * pi / 180)) +  # Adjust for latitude
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank())

# Calculate the mean latent extinction risk for the AMT, excluding NA values
avg_latent_risk <- mean(latent_stack, na.rm = TRUE) %>%
  mask(r_amt)

# Reproject to WGS84 (EPSG:4326) geographic coordinate system for plotting
avg_latent_risk <- project(avg_latent_risk, "EPSG:4326")

# Convert the raster data to a data frame
avg_latent_rast_df <- as.data.frame(avg_latent_risk, xy = TRUE)

# Plot the raster data using ggplot2
l2 <- ggplot(avg_latent_rast_df, aes(x = x, y = y, fill = mean)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D",
                       breaks = seq(-4, 2, by = 1), 
                       labels = seq(-4, 2, by = 1)) +
  labs(fill = "Mean") +
  theme_void() +
  coord_fixed(ratio = 1 / cos(mean(avg_latent_rast_df$y) * pi / 180)) +  # Adjust for latitude
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank()) 


#### CREATE MAPS OF SUM AND MEAN POSITIVE LATENT RISK ####

# Extract species with positive latent risk values
pos_ler <- iucn_threat[which(iucn_threat$latent_risk > 0),]
pos_ler_rasters <- list()

for (i in 1:length(pos_ler$scientificName)) {
  curr_latent_risk <- pos_ler$latent_risk[i]
  curr_vect <- vect(pos_ler[i, "geometry"]) # Convert the current geometry to a terra vector object
  curr_rast <- rasterize(curr_vect, r_amt, field = curr_latent_risk, background = NA) # Rasterize the vector object using the template raster and threat vals
  pos_ler_rasters[[i]] <- curr_rast # Store the raster in the list
}

# Combine the list of rasters into a single raster stack
pos_latent_stack <- rast(pos_ler_rasters)

# Calculate the net latent extinction risk for the AMT, excluding NA values
sum_pos_latent_risk <- sum(pos_latent_stack, na.rm = TRUE) %>%
  mask(r_amt)

# Reproject to WGS84 (EPSG:4326) geographic coordinate system for plotting
sum_pos_latent_risk <- project(sum_pos_latent_risk, "EPSG:4326")

# Convert the raster data to a data frame
sum_pos_latent_rast_df <- as.data.frame(sum_pos_latent_risk, xy = TRUE)

# Plot the raster data using ggplot2
pl1 <- ggplot(sum_pos_latent_rast_df, aes(x = x, y = y, fill = sum)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D",
                       breaks = seq(0, 40, by = 10), 
                       labels = seq(0, 40, by = 10)) +
  labs(fill = "Sum") +
  theme_void() +
  coord_fixed(ratio = 1 / cos(mean(sum_pos_latent_rast_df$y) * pi / 180)) +  # Adjust for latitude
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank()) +
  ggtitle('(a) Total positive latent risk')


plot(pl1)

# Calculate the mean latent extinction risk for the AMT, excluding NA values
avg_pos_latent_risk <- mean(pos_latent_stack, na.rm = TRUE) %>%
  mask(r_amt)

# Reproject to WGS84 (EPSG:4326) geographic coordinate system for plotting
avg_pos_latent_risk <- project(avg_pos_latent_risk, "EPSG:4326")

# Convert the raster data to a data frame
avg_pos_latent_rast_df <- as.data.frame(avg_pos_latent_risk, xy = TRUE)

# Plot the raster data using ggplot2
pl2 <- ggplot(avg_pos_latent_rast_df, aes(x = x, y = y, fill = mean)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D",
                       breaks = seq(0, 2, by = 0.5), 
                       labels = seq(0, 2, by = 0.5)) +
  labs(fill = "Mean") +
  theme_void() +
  coord_fixed(ratio = 1 / cos(mean(avg_pos_latent_rast_df$y) * pi / 180)) +  # Adjust for latitude
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank())  + 
  ggtitle('(b) Mean positive latent risk')

plot(pl2)


#### CREATE MAPS OF SUM AND MEAN CURRENT RISK ####

# Initialize an empty list to store the raster layers
curr_risk_rasters <- list()
for (i in 1:length(iucn_threat$scientificName)) {
  curr_risk <- iucn_threat$ordinalThreat[i]
  curr_vect <- vect(iucn_threat[i, "geometry"]) # Convert the current geometry to a terra vector object
  curr_rast <- rasterize(curr_vect, r_amt, field = curr_risk, background = NA) # Rasterize the vector object using the template raster and threat vals
  curr_risk_rasters[[i]] <- curr_rast # Store the raster in the list
}
# Combine the list of rasters into a single raster stack
curr_risk_stack <- rast(curr_risk_rasters)

# Calculate the net latent extinction risk for the AMT, excluding NA values
sum_curr_risk <- sum(curr_risk_stack, na.rm = TRUE) %>%
  mask(r_amt)

# Reproject to WGS84 (EPSG:4326) geographic coordinate system for plotting
sum_curr_risk <- project(sum_curr_risk, "EPSG:4326")

# Convert the raster data to a data frame
sum_curr_risk_df <- as.data.frame(sum_curr_risk, xy = TRUE)

# Plot the raster data using ggplot2
c1 <- ggplot(sum_curr_risk_df, aes(x = x, y = y, fill = sum)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D",
                       breaks = seq(0, 120, by = 20), 
                       labels = seq(0, 120, by = 20)) +
  labs(fill = "Sum") +
  theme_void() +
  coord_fixed(ratio = 1 / cos(mean(sum_curr_risk_df$y) * pi / 180)) +  # Adjust for latitude
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank()) + 
  ggtitle('(c) Total current risk')

plot(c1)

# Calculate the mean latent extinction risk for the AMT, excluding NA values
avg_curr_risk <- mean(curr_risk_stack, na.rm = TRUE) %>%
  mask(r_amt)

# Reproject to WGS84 (EPSG:4326) geographic coordinate system for plotting
avg_curr_risk <- project(avg_curr_risk, "EPSG:4326")

# Convert the raster data to a data frame
avg_curr_risk_df <- as.data.frame(avg_curr_risk, xy = TRUE)

# Plot the raster data using ggplot2
c2 <- ggplot(avg_curr_risk_df, aes(x = x, y = y, fill = mean)) +
  geom_raster() +
  scale_fill_viridis_c(option = "D", 
                       breaks = seq(0, 6, by = 1), 
                       labels = seq(0, 6, by = 1)) +
  labs(fill = "Mean") +
  theme_void() +
  coord_fixed(ratio = 1 / cos(mean(avg_curr_risk_df$y) * pi / 180)) + 
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank()) + 
  ggtitle('(d) Mean current risk')

plot(c2)

# Combine the four maps into a 2x2 grid
combined_plot <- (pl1 | pl2) / (c1 | c2)

# Print the combined plot
print(combined_plot)

ggsave("02_combined_risk_plots_PLER.png", units="cm", width=28, height=15, dpi=600, path = "results/fig/pgls", bg  = 'white')

