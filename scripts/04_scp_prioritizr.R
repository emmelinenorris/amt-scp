###############################################################
# AMT Spatial Conservation Prioritisation for Mammals
#
# Part 4: Spatial Conservation Prioritisation using Prioritizr
#
# Code developed by Emmeline Norris, 08/06/2024
# 
###############################################################

# Install latest version of prioritizr package
install.packages("prioritizr", repos = "https://cran.rstudio.com/")
# Install Gurobi optimiser (v 11.0-1)
install.packages("c:/gurobi1101/win64/R/gurobi_11.0-1.zip", repos = NULL)
# Install slam package (required to run gurobi)
install.packages("slam", repos = "https://cloud.r-project.org")
# Install HiGHS solver
install.packages("highs", repos = "https://cran.rstudio.com/")

library(ggplot2)
library(prioritizr)
library(gurobi)
library(highs)
library(sf)
library(terra)
library(tidyverse)
library(viridis)

# Disable scientific notation
options(scipen = 999)


#### LOAD DATA REQUIRED FOR PRIORITIZR AND PLOTS ####

# Load spec_dat table of conservation features
spec_dat <- read_csv("data/output-data/tbl/03_spec_dat.csv")

# Export the raster stack of species distributions (i.e. the 'conservation features' as a TIFF file)
species <- rast("data/output-data/tif/03_species_stack.tif")

# Load shapefile of hexagonal planning unit tesselation for the AMT
pu_dat <- st_read("data/output-data/shp/03_amt_pu_grid.shp")
colnames(pu_dat)[6:23] <- c("pastoral","status","locked_in", "locked_out_ob2",
                            "locked_out_ob3", "cost","cost_adjusted", "cost_equal", "main_id", 
                            "main_type", "category", "declare_ipa", "median_income",
                            "hii_oceania", "normal_hii", "normal_mei",
                            "id_spec_richness", "spec_richness")

# For Objective 1 (add to existing PA network), lock out freehold land etc.
pu_dat <- pu_dat %>%
  mutate(locked_out_ob1 = if_else(category == "Other freehold, term, perpetual lease or Crown purposes",TRUE, FALSE))

# Filter planning units that are locked out in Objective 1
locked_out_ob1 <- pu_dat %>%
  filter(category == "Other freehold, term, perpetual lease or Crown purposes")

# For Objective 2 (add to existing PA network from only non-Indigenous land), lock out native title & Indigenous freehold
pu_dat <- pu_dat %>%
  mutate(locked_out_ob2 = if_else(category == "Other freehold, term, perpetual lease or Crown purposes",
                                  TRUE, as.logical(declare_ipa)))

# Filter planning units that are locked out in Objective 1
locked_out_ob2 <- pu_dat %>%
  filter(category == "Other freehold, term, perpetual lease or Crown purposes" | declare_ipa == 1)

# For Objective 3 (add to existing PA network from only Indigenous land), lock out pastoral land
pu_dat <- pu_dat %>%
  mutate(locked_out_ob3 = if_else(category == "Other freehold, term, perpetual lease or Crown purposes",
                                  TRUE, as.logical(pastoral)))

# Filter planning units that are locked out in Objective 1
locked_out_ob3 <- pu_dat %>%
  filter(category == "Other freehold, term, perpetual lease or Crown purposes" | pastoral == 1)

# Convert locked_in column back to logical variable
pu_dat$locked_in <- as.logical(pu_dat$locked_in)


# Calculate number of PUs in the PAN
pa_pus <- length(pu_dat$status[which(pu_dat$status==1)])  # 3,469 PUs in already locked into existing protected area network

# Calculate the total cost of a solution.
pa_cost <- sum(pu_dat$cost[which(pu_dat$status == 1)], na.rm = T) # cost of area in PAN = 733.3

# Load shapefile of AMT boundary
amt <- st_read("data/output-data/shp/01_amt.shp")

# Load shapefile of Wet Tropics boundary
wt <- st_read("data/output-data/shp/01_wet-tropics.shp")

#### OBJECTIVE 1 (add to existing PA network); CRITERIA 1 (currently threatened species) ####

cp1.1 = problem(pu_dat, species, cost_column = "cost") %>%
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_1) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_locked_out_constraints("locked_out_ob1") %>% # lock out freehold land etc.
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s1.1 = solve(cp1.1)
s1.1 <- st_as_sf(s1.1) # Convert to sf 
s1.1x <- s1.1[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp1.1, s1.1x)[[2]] - pa_pus # TOTAL PUs = 319

# Calculate total cost of solution
eval_cost_summary(cp1.1, s1.1x)[[2]] - pa_cost # TOTAL COST = 52.13425

# Calculate how well feature representation targets are met by a solution and 
# the proportion of species' distributions covered
tcs1.1 <- eval_target_coverage_summary(cp1.1, s1.1x) # 100% of targets for threatened species met or exceeded

# Calculate total area of PA network (including existing PAs)
sum(st_area(s1.1[which(s1.1x$solution_1 == 1),])) # Approx. 500,093 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp1.1, s1.1x) # Boundary length = 29,577 km

# Calculate the number of selected PUs in each land use category
s1.1_df <- s1.1 %>%
  filter(solution_1 == 1) %>%
  group_by(category) %>%
  summarise(Count = n()) %>%
  mutate(status = case_when(
    category %in% c("Protected Area", "Indigenous Protected Area") ~ "Protected",
    TRUE ~ "Not Protected"))

# Calculate the proportion of PUs in each land use category (excluding protected areas, which are 'locked in')
s1.1_df <- s1.1_df %>%
  mutate(prop_selected = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order below
s1.1_df$category <- factor(s1.1_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Native Title land - Pastoral use",
                                                        "Indigenous Land Use Agreement",
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease"))

# Plot optimal solution
ggplot(data = s1.1_df) +
  geom_sf(mapping = aes(fill = status), color = "transparent") +
  scale_fill_manual(values = c("Not Protected" = "#af8dc3", "Protected" = "#7fbf7b"),
                    labels = c("Selected; Not Protected",
                               "Protected Area")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "black", color = "black", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))


# Plot optimal solution
ggplot(data = s1.1_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_viridis(discrete = T, alpha = 0.85, option = "D", 
                    labels = c("Protected Area", 
                               "Indigenous Protected Area", 
                               "Freehold - Indigenous", 
                               "Native Title land",
                               "Native Title land - Pastoral use",
                               "Indigenous Land Use Agreement",
                               "Indigenous Land Use Agreement - Pastoral use",
                               "Pastoral term or perpetual lease")) +
  geom_sf(data = locked_out_ob1, fill = "gray70", color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray35", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp1_1.png", units="cm", width=30, height=15, dpi=300, path = "results/fig/scp", bg  = 'white')


#### OBJECTIVE 1 (add to existing PA network); CRITERIA 2 (positive latent risk) ####

cp1.2 = problem(pu_dat, species, cost_column = "cost") %>%
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_2) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_locked_out_constraints("locked_out_ob1") %>% # lock out freehold land etc.
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s1.2 = solve(cp1.2)
s1.2 <- st_as_sf(s1.2) # Convert to sf 
s1.2x <- s1.2[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp1.2, s1.2x)[[2]] - pa_pus # TOTAL PUs = 872

# Calculate total cost of solution
eval_cost_summary(cp1.2, s1.2x)[[2]] - pa_cost # TOTAL COST = 142.0343

# Calculate how well feature representation targets are met by a solution and the proportion of species' distributions covered
tcs1.2 <- eval_target_coverage_summary(cp1.2, s1.2x)
# 100% of targets for high positive latent risk species met or exceeded

# Calculate total area of PA network (including existing PAs)
sum(st_area(s1.2[which(s1.2x$solution_1 == 1),])) # Approx. 580,436.8 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp1.2, s1.2x) # Boundary length = 31,830 km

# Calculate the number of selected PUs in each land use category
s1.2_df <- s1.2 %>%
  filter(solution_1 == 1) %>%
  group_by(category) %>%
  summarise(Count = n()) %>%
  mutate(status = case_when(
    category %in% c("Protected Area", "Indigenous Protected Area") ~ "Protected",
    TRUE ~ "Not Protected"))

# Calculate the proportion of PUs in each land use category (excluding protected areas, which are 'locked in')
s1.2_df <- s1.2_df %>%
  mutate(prop_selected = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order
s1.2_df$category <- factor(s1.2_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Native Title land - Pastoral use",
                                                        "Indigenous Land Use Agreement",
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease"))

# Plot optimal solution
ggplot(data = s1.2_df) +
  geom_sf(mapping = aes(fill = status), color = "transparent") +
  scale_fill_manual(values = c("Not Protected" = "#af8dc3", "Protected" = "#7fbf7b"),
                    labels = c("Selected; Not Protected",
                               "Protected Area")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "black", color = "black", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# Plot optimal solution
ggplot(data = s1.2_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_viridis(discrete = T, alpha = 0.85, option = "D", 
                     labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land",
                                "Native Title land - Pastoral use",
                                "Indigenous Land Use Agreement",
                                "Indigenous Land Use Agreement - Pastoral use",
                                "Pastoral term or perpetual lease")) +
  geom_sf(data = locked_out_ob1, fill = "gray70", color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray35", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp1_2.png", units="cm", width=30, height=15, dpi=300, path = "results/fig/scp", bg  = 'white')


#### OBJECTIVE 1 (add to existing PA network); CRITERIA 3 (species richness) ####

cp1.3 = problem(pu_dat, species, cost_column = "cost") %>%
  add_max_features_objective() %>% # sets same minimum objective
  add_relative_targets(0.1) %>% # specifies relative targets species richness classes
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_locked_out_constraints("locked_out_ob1") %>% # lock out freehold land etc.
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s1.3 = solve(cp1.3)
s1.3 <- st_as_sf(s1.3) # Convert to sf 
s1.3x <- s1.3[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp1.3, s1.3x)[[2]] - pa_pus # TOTAL PUs = 322 

# Calculate total cost of solution
eval_cost_summary(cp1.3, s1.3x)[[2]] - pa_cost # TOTAL COST = 77.99684

# Calculate how well feature representation targets are met by a solution
tcs1.3 <- eval_target_coverage_summary(cp1.3, s1.3x)
# Species richness target exceeded

# Calculate total area of PA network (including existing PAs)
sum(st_area(s1.3[which(s1.3x$solution_1 == 1),]))/1000000 # Approx. 501,056.3 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp1.3, s1.3x) # Boundary length = 28,366.12 km

# Calculate the number of selected PUs in each land use category
s1.3_df <- s1.3 %>%
  filter(solution_1 == 1) %>%
  group_by(category) %>%
  summarise(Count = n()) %>%
  mutate(status = case_when(
    category %in% c("Protected Area", "Indigenous Protected Area") ~ "Protected",
    TRUE ~ "Not Protected"))

# Calculate the proportion of PUs in each land use category (excluding protected areas, which are 'locked in')
s1.3_df <- s1.3_df %>%
  mutate(prop_selected = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order
s1.3_df$category <- factor(s1.3_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Indigenous Land Use Agreement", 
                                                        "Native Title land - Pastoral use",
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease", 
                                                        "Other freehold, term, perpetual lease or Crown purposes"))

# Plot optimal solution
ggplot(data = s1.3_df) +
  geom_sf(mapping = aes(fill = status), color = "transparent") +
  scale_fill_manual(values = c("Not Protected" = "#af8dc3", "Protected" = "#7fbf7b"),
                    labels = c("Selected; Not Protected",
                               "Protected Area")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "black", color = "black", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# Plot optimal solution
ggplot(data = s1.3_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(values = c("#35978f",  # Protected Area
                               "#2d6969",  # Indigenous Protected Area
                               "#45581c",  # Freehold - Indigenous
                               "#72922f",  # Native Title land
                               "#afd57e",  # Indigenous Land Use Agreement
                               "#976222",  # Native Title land - Pastoral use
                               "#d09855",  # Indigenous Land Use Agreement - Pastoral use
                               "#e4c49e",  # Pastoral term or perpetual lease
                               "#e15fab"), # Other freehold, term, perpetual lease or Crown purposes
                    labels = c("Protected Area", 
                               "Indigenous Protected Area", 
                               "Freehold - Indigenous", 
                               "Native Title land",
                               "Indigenous Land Use Agreement",
                               "Native Title land - Pastoral use",
                               "Indigenous Land Use Agreement - Pastoral use",
                               "Pastoral term or perpetual lease", 
                               "Other freehold, term, perpetual lease or Crown purposes")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "black", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp1_3.tiff", units="cm", width=30, height=15, dpi=300, compression = 'lzw', path = "results/fig/scp", bg  = 'white')


#### OBJECTIVE 2 (growing IPA network through voluntary declaration of Native Title land); CRITERIA 1 (threatened species) ####

## CRITERIA 1: PRIORITISE CURRENTLY THREATENED SPECIES

cp2.1 = problem(pu_dat, species, cost_column = "cost_adjusted") %>% # Use cost adjusted to reduce cost of declaring Native Title/Indigenous Freehold land
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_1) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s2.1 = solve(cp2.1)
s2.1 <- st_as_sf(s2.1) # Convert to sf 
s2.1x <- s2.1[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp2.1, s2.1x)[[2]] - pa_pus # TOTAL PUs = 313

# Calculate total cost of solution
eval_cost_summary(cp2.1, s2.1x)[[2]] - pa_cost # TOTAL COST = 28.31674

# Calculate how well feature representation targets are met by a solution and the proportion of species' distributions covered
tcs2.1 <- eval_target_coverage_summary(cp2.1, s2.1x)
# 100% of targets for threatened species met or exceeded

# Calculate total area of PA network (including existing PAs)
sum(st_area(s2.1[which(s2.1x$solution_1 == 1),]))/1000000 # Approx. 499,479 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp2.1, s2.1x) # Boundary length = 30,249 km

# Calculate the number of selected PUs in each land use category
s2.1_df <- s2.1 %>%
  filter(solution_1 == 1) %>%
  group_by(category) %>%
  summarise(Count = n()) %>%
  mutate(status = case_when(
    category %in% c("Protected Area", "Indigenous Protected Area") ~ "Protected",
    TRUE ~ "Not Protected"))

# Calculate the proportion of PUs in each land use category (excluding protected areas, which are 'locked in')
s2.1_df <- s2.1_df %>%
  mutate(prop_selected = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order
s2.1_df$category <- factor(s2.1_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Indigenous Land Use Agreement", 
                                                        "Native Title land - Pastoral use",
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease", 
                                                        "Other freehold, term, perpetual lease or Crown purposes"))

# Plot optimal solution
ggplot(data = s2.1_df) +
  geom_sf(mapping = aes(fill = status), color = "transparent") +
  scale_fill_manual(values = c("Not Protected" = "#af8dc3", "Protected" = "#7fbf7b"),
                    labels = c("Selected; Not Protected",
                               "Protected Area")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "black", color = "black", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# Plot optimal solution
ggplot(data = s2.1_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(values = c("#35978f",  # Protected Area
                               "#2d6969",  # Indigenous Protected Area
                               "#45581c",  # Freehold - Indigenous
                               "#72922f",  # Native Title land
                               "#afd57e",  # Indigenous Land Use Agreement
                               "#976222",  # Native Title land - Pastoral use
                               "#d09855",  # Indigenous Land Use Agreement - Pastoral use
                               "#e4c49e",  # Pastoral term or perpetual lease
                               "#e15fab"), # Other freehold, term, perpetual lease or Crown purposes
                    labels = c("Protected Area", 
                               "Indigenous Protected Area", 
                               "Freehold - Indigenous", 
                               "Native Title land",
                               "Indigenous Land Use Agreement",
                               "Native Title land - Pastoral use",
                               "Indigenous Land Use Agreement - Pastoral use",
                               "Pastoral term or perpetual lease", 
                               "Other freehold, term, perpetual lease or Crown purposes")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "black", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp2_1.tiff", units="cm", width=30, height=15, dpi=300, compression = 'lzw', path = "results/fig/scp", bg  = 'white')


#### OBJECTIVE 2 (growing IPA network through voluntary declaration of Native Title land); CRITERIA 2 (positive latent risk) ####

## CRITERIA 2: PRIORITISE HIGH POSITIVE LATENT RISK SPECIES

cp2.2 = problem(pu_dat, species, cost_column = "cost_adjusted") %>%
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_2) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s2.2 = solve(cp2.2)
s2.2 <- st_as_sf(s2.2) # Convert to sf 
s2.2x <- s2.2[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp2.2, s2.2x)[[2]] - pa_pus # TOTAL PUs = 850

# Calculate total cost of solution
eval_cost_summary(cp2.2, s2.2x)[[2]] - pa_cost # TOTAL COST = 83.24743

# Calculate how well feature representation targets are met by a solution and the proportion of species' distributions covered
tcs2.2 <- eval_target_coverage_summary(cp2.2, s2.2x)
# 100% of targets for threatened species met or exceeded

# Calculate total area of PA network (including existing PAs)
sum(st_area(s2.2[which(s2.2x$solution_1 == 1),]))/1000000 # Approx. 578,661.7 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp2.2, s2.2x) # Boundary length = 30,248.46 km

# Calculate the number of selected PUs in each land use category
s2.2_df <- s2.2 %>%
  filter(solution_1 == 1) %>%
  group_by(category) %>%
  summarise(Count = n()) %>%
  mutate(status = case_when(
    category %in% c("Protected Area", "Indigenous Protected Area") ~ "Protected",
    TRUE ~ "Not Protected"))

# Calculate the proportion of PUs in each land use category (excluding protected areas, which are 'locked in')
s2.2_df <- s2.2_df %>%
  mutate(prop_selected = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order
s2.2_df$category <- factor(s2.2_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Indigenous Land Use Agreement", 
                                                        "Native Title land - Pastoral use",
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease", 
                                                        "Other freehold, term, perpetual lease or Crown purposes"))

# Plot optimal solution
ggplot(data = s2.2_df) +
  geom_sf(mapping = aes(fill = status), color = "transparent") +
  scale_fill_manual(values = c("Not Protected" = "#af8dc3", "Protected" = "#7fbf7b"),
                    labels = c("Selected; Not Protected",
                               "Protected Area")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "black", color = "black", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# Plot optimal solution
ggplot(data = s2.2_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(values = c("#35978f",  # Protected Area
                               "#2d6969",  # Indigenous Protected Area
                               "#45581c",  # Freehold - Indigenous
                               "#72922f",  # Native Title land
                               "#afd57e",  # Indigenous Land Use Agreement
                               "#976222",  # Native Title land - Pastoral use
                               "#d09855",  # Indigenous Land Use Agreement - Pastoral use
                               "#e4c49e",  # Pastoral term or perpetual lease
                               "#e15fab"), # Other freehold, term, perpetual lease or Crown purposes
                    labels = c("Protected Area", 
                               "Indigenous Protected Area", 
                               "Freehold - Indigenous", 
                               "Native Title land",
                               "Indigenous Land Use Agreement",
                               "Native Title land - Pastoral use",
                               "Indigenous Land Use Agreement - Pastoral use",
                               "Pastoral term or perpetual lease", 
                               "Other freehold, term, perpetual lease or Crown purposes")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "black", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp2_2.tiff", units="cm", width=30, height=15, dpi=300, compression = 'lzw', path = "results/fig/scp", bg  = 'white')


#### OBJECTIVE 2 (growing IPA network through voluntary declaration of Native Title land); CRITERIA 3 (species richness) ####

cp2.3 = problem(pu_dat, richness_stack, cost_column = "cost_adjusted") %>%
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(richness_df$target_3) %>% # specifies relative targets species richness classes
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s2.3 = solve(cp2.3)
s2.3 <- st_as_sf(s2.3) # Convert to sf 
s2.3x <- s2.3[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp2.3, s2.3x)[[2]] - pa_pus # TOTAL PUs = 332

# Calculate total cost of solution
eval_cost_summary(cp2.3, s2.3x)[[2]] - pa_cost # TOTAL COST = 47.23799

# Calculate how well feature representation targets are met by a solution
tcs2.3 <- eval_target_coverage_summary(cp2.3, s2.3x)
# Species richness target exceeded

# Calculate total area of PA network (including existing PAs)
sum(st_area(s2.3[which(s2.3x$solution_1 == 1),]))/1000000 # Approx. 501,929.5 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp2.3, s2.3x) # Boundary length = 28,383.649 km

# Calculate the number of selected PUs in each land use category
s2.3_df <- s2.3 %>%
  filter(solution_1 == 1) %>%
  group_by(category) %>%
  summarise(Count = n()) %>%
  mutate(status = case_when(
    category %in% c("Protected Area", "Indigenous Protected Area") ~ "Protected",
    TRUE ~ "Not Protected"))

# Calculate the proportion of PUs in each land use category (excluding protected areas, which are 'locked in')
s2.3_df <- s2.3_df %>%
  mutate(prop_selected = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order
s2.3_df$category <- factor(s2.3_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Indigenous Land Use Agreement", 
                                                        "Native Title land - Pastoral use",
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease", 
                                                        "Other freehold, term, perpetual lease or Crown purposes"))

# Plot optimal solution
ggplot(data = s2.3_df) +
  geom_sf(mapping = aes(fill = status), color = "transparent") +
  scale_fill_manual(values = c("Not Protected" = "#af8dc3", "Protected" = "#7fbf7b"),
                    labels = c("Selected; Not Protected",
                               "Protected Area")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "black", color = "black", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

# Plot optimal solution
ggplot(data = s2.3_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(values = c("#35978f",  # Protected Area
                               "#2d6969",  # Indigenous Protected Area
                               "#45581c",  # Freehold - Indigenous
                               "#72922f",  # Native Title land
                               "#afd57e",  # Indigenous Land Use Agreement
                               "#976222",  # Native Title land - Pastoral use
                               "#d09855",  # Indigenous Land Use Agreement - Pastoral use
                               "#e4c49e",  # Pastoral term or perpetual lease
                               "#e15fab"), # Other freehold, term, perpetual lease or Crown purposes
                    labels = c("Protected Area", 
                               "Indigenous Protected Area", 
                               "Freehold - Indigenous", 
                               "Native Title land",
                               "Indigenous Land Use Agreement",
                               "Native Title land - Pastoral use",
                               "Indigenous Land Use Agreement - Pastoral use",
                               "Pastoral term or perpetual lease", 
                               "Other freehold, term, perpetual lease or Crown purposes")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "black", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp2_3.tiff", units="cm", width=30, height=15, dpi=300, compression = 'lzw', path = "results/fig/scp", bg  = 'white')


#### CALCULATE PLANNING UNIT SELECTION FREQUENCY ACROSS ALL SCENARIOS ####

# Create a table containing the PU selections (binary) for each SCP scenario
sel_freq <- tibble(
  s1.1 = s1.1$solution_1,
  s1.2 = s1.2$solution_1,
  s1.3 = s1.3$solution_1,
  s2.1 = s2.1$solution_1,
  s2.2 = s2.2$solution_1,
  s2.3 = s2.3$solution_1) %>%
  mutate(sum_values = rowSums(across(everything()))) %>%
  bind_cols(pu_dat) %>%
  filter(status != 1) %>%
  st_as_sf()

# Create sf data frame containing only PUs in PA network
pa_network <- pu_dat %>%
  filter(status == 1) %>%
  st_as_sf()

# Plot selection frequencies
ggplot(data = sel_freq) +
  geom_sf(mapping = aes(fill = factor(sum_values)), color = "transparent") +
  scale_fill_viridis_d(option = "viridis", 
                       direction = 1, 
                       name = "Selection frequency", 
                       breaks = 0:5, 
                       labels = as.character(0:5),
                       begin = 0.15,
                       end = 0.95) +
  geom_sf(data = pa_network, fill = 'grey80', color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "black", color = "black", size = 0.4) +
  theme_minimal()

ggsave("04_sel_freq.tiff", units="cm", width=30, height=15, dpi=300, compression = 'lzw', path = "results/fig/scp", bg  = 'white')



#### CALCULATE PROPORTION OF SPECIES AMT DISTRIBUTION COVERED BY AMT PROTECTED AREA NETWORK ####

# Load shapefile of IUCN distribution polygons for AMT mammals with threat status data generated in 02_pgls_model
# Note: iucn_threat contains geometry for entire species distribution
iucn_threat <-  read_sf("data/output-data/shp/02_iucn_threat.shp")
colnames(iucn_threat)[1:4] <- c("scientificName", "ordinalThreat", "mean_fitted", "latent_risk")
iucn_amt <- iucn_threat %>% 
  st_intersection(amt)

# Load CAPAD protected areas cropped to AMT boundary
capad_amt <- read_sf("data/output-data/shp/03_capad_amt.shp")

ggplot() +
  geom_sf(data = amt, color = "black", alpha = 0.5) +
  geom_sf(data = capad_amt, color = "black", fill = 'lightgreen', alpha = 0.5) +
  geom_sf(data = iucn_amt[7,], fill = "red", color = "black", alpha = 0.5)

# Calculate area of species' AMT distribution that is covered by existing AMT PA network
spec_dat$area_pa <- NA
spec_dat$prop_protected <- NA
for (i in 1:length(iucn_amt$scientificName)) {
  curr_shape <- iucn_amt[i, "geometry"]
  curr_intersect <- st_intersection(curr_shape, capad_amt)
  curr_area <- st_area(curr_intersect)
  spec_dat$area_pa[i] <- curr_area
}

# Calculate the proportion of each species' AMT distributions that covered by AMT PA network
spec_dat <- spec_dat %>%
  mutate(prop_protected = pmin(area_pa/area_amt, 1))



