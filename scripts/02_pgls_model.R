######################################################
# AMT Spatial Conservation Prioritisation for Mammals
#
# Part 2: Phylogenetic Least-Squares Models
#
# Code developed by Emmeline Norris, 01/06/2024
# 
######################################################

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
AICc(full_model$lm.res) # 524.8646


## REDUCED MODELS ##

# 1
# Remove predictors based on highest p-values to see if it affects model fit
# Remove temperature of wettest quarter (bio8) variable (p = 0.5158)
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
AICc(red_model_1$lm.res) # 522.8644


# 2
# Remove human influence index (hii) variable (p = 0.5329)
pred_red_2 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "litter_size_n_tr_st",
                "litters_per_year_n_tr_st",
                "age_first_reproduction_d_tr_st",
                "foxAreaProp_tr_st",
                "toadAreaProp_tr_st",
                "mean_firefreq_st",
                "mean_latefire_st",
                "rangeArea_km2_tr_st")

red_model_2 <- paste("ordinalThreat ~", paste(pred_red_2, collapse="+"))
fmla <- as.formula(red_model_2)
red_model_2 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
red_model_2$p.res 
summary(red_model_2$lm.res)

# Model diagnostics
qqnorm(red_model_2$lm.res)
plot(red_model_2$lm.res)
AICc(red_model_2$lm.res) # 520.8624


# 3
# Remove litter size variable (p = 0.4814)
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
AICc(red_model_3$lm.res) # 519.0112


# 4
# Remove cane toad range overlap variable (p = 0.2302)
pred_red_4 <- c("adult_mass_g_tr_st",
                "I(adult_mass_g_tr_st^2)",
                "litters_per_year_n_tr_st",
                "age_first_reproduction_d_tr_st",
                "foxAreaProp_tr_st",
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
AICc(red_model_4$lm.res) # 518.2123


# 5
# Remove litters per year variable (p = 0.1509)
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
AICc(red_model_5$lm.res) # 518.1067


# 6
# Remove fox range overlap variable (p = 0.1028)
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
AICc(red_model_6$lm.res) # 518.6559 - AICc increases slightly after removal of fox range overlap var


# 7
# Remove geographic range area variable (p = 0.0954)
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
AICc(red_model_7$lm.res) # 519.3437


## Reduced model #5 has lowest AICc value (518.1067) by a very small amount

# Construct a null model for comparison
fmla_null <- as.formula('ordinalThreat ~ 1') # intercept-only model without predictors
model_null <- fitspphylo(formula=fmla_null, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)

# Likelihood ratio test for full and reduced models
lmtest::lrtest(full_model$lm.res, red_model_5$lm.res)
# P-value of 0.3918 indicates that the additional vars in full model do not improve model fit

# Likelihood ratio test for reduced models #5 and null model
lmtest::lrtest(red_model_5$lm.res, model_null$lm.res)
# P-value of < 0.05 indicates the reduced model is significantly better at explaining variability in data compared to null model

# Calculate R^2
r2 <- R2_lik(mod = red_model_5$lm.res, mod.r = model_null$lm.res)
r2*100 # r^2 = 30.02903

# Calculate delta AICc for full and minimum adequate model 
d_AICc <- AICc(full_model$lm.res) - AICc(red_model_5$lm.res)
d_AICc # 6.757928

#### PREDICT LATENT EXTINCTION RISK FOR AMT MAMMALS ####

# Load data frame containing AMT species names and ordinal threat status, generated in 01_load_and_clean
amt_species_df <- read.csv("data/output-data/tbl/01_amt_species_df.csv")
# Replace space between genus and species with underscore
amt_species_df$scientificName <- gsub(" ", "_", amt_species_df$scientificName)

# Run best fitting reduced model with 100 phylogenies
for (i in 1:100){ 
  compdat <- comparative.data(phy=amt_phy[[1]], data= amt_dat_st, names.col="scientificName")
  red_model_5 <- paste("ordinalThreat ~", paste(pred_red_5, collapse="+"))
  fmla <- as.formula(red_model_5)
  red_model_5 <- fitspphylo(formula=fmla, data=compdat$data, spmatrix=spmat, phylomatrix=phylomat, p=p)
  saveRDS(red_model_5$lm.res, file=paste0("data/output-data/mod/02_model", i, ".rds"))
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


#### CREATE MAP OF MEAN LATENT RISK FROM MODELS ####

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
plot(sum_latent_risk)
# Calculate the mean latent extinction risk for the AMT, excluding NA values
avg_latent_risk <- mean(latent_stack, na.rm = TRUE) %>%
  mask(r_amt)
plot(avg_latent_risk)




