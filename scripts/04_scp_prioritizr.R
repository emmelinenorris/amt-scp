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
library(patchwork)
library(cowplot)

# Disable scientific notation
options(scipen = 999)

#### LOAD DATA REQUIRED FOR PRIORITIZR AND PLOTS ####

# Load AMT and WT polygons created in 01_load_and_clean
amt <- st_read("data/output-data/shp/01_amt.shp")
wt <- st_read("data/output-data/shp/01_wet-tropics.shp")

# Load amt raster from 01_load_and_clean
r_amt <- rast("data/output-data/tif/01_r_amt.tif")

# Load spec_dat table of conservation features
spec_dat <- read_csv("data/output-data/tbl/03_spec_dat.csv")

# Export the raster stack of species distributions (i.e. the 'conservation features' as a TIFF file)
species <- rast("data/output-data/tif/03_species_stack.tif")

# Load shapefile of IUCN distribution polygons for mammals in the AMT generated in 01_load_and_clean
iucn_union_mod <-  read_sf("data/output-data/shp/01_iucn_union_mod.shp")
colnames(iucn_union_mod)[1] <- "scientificName"

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

# For Objective 2 (add to existing PA network from only Indigenous land), lock out pastoral land
pu_dat <- pu_dat %>%
  mutate(locked_out_ob2 = if_else(category == "Other freehold, term, perpetual lease or Crown purposes" |
                                    category == "Indigenous Land Use Agreement" |
                                    category == "Indigenous Land Use Agreement - Pastoral use" |
                                    category == "Pastoral term or perpetual lease",
                                  TRUE, FALSE))

# Filter planning units that are locked out in Objective 2
locked_out_ob2 <- pu_dat %>%
  filter(category == "Other freehold, term, perpetual lease or Crown purposes" |
           category == "Indigenous Land Use Agreement" |
           category == "Indigenous Land Use Agreement - Pastoral use" |
           category == "Pastoral term or perpetual lease")

# Convert locked_in column back to logical variable
pu_dat$locked_in <- as.logical(pu_dat$locked_in)

# Filter planning units that are in the national reserve system (non-IPA)
protected_area <- pu_dat %>%
  filter(category == "Protected Area" | category == "Indigenous Protected Area")

# Filter planning units that are in the IPA network
ipa <- pu_dat %>%
  filter(category == "Indigenous Protected Area")

# Calculate number of PUs in the PAN
pa_pus <- length(pu_dat$status[which(pu_dat$status==1)])  # 3,469 PUs in already locked into existing protected area network

# Calculate the total cost of the existing PAN
pa_cost <- sum(pu_dat$cost[which(pu_dat$status == 1)], na.rm = T) # cost of area in PAN = 733.3

# Calculate the total area of the existing PAN
pa_area <- pu_dat %>%
  filter(status == 1) %>%
  st_area() %>%
  sum()

#### OBJECTIVE 1 (add to existing PA network); CRITERIA 1 (currently threatened species) ####

cp1.1 = problem(pu_dat, species, cost_column = "cost") %>%
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_1) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  #add_locked_out_constraints("locked_out_ob1") %>% # lock out freehold land etc.
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  #add_top_portfolio(number_solutions = 3) %>% # generate portfolio of top 10 optimal solutions
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s1.1 = solve(cp1.1)
s1.1 <- st_as_sf(s1.1) # Convert to sf 
s1.1x <- s1.1[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp1.1, s1.1x)[[2]] - pa_pus # TOTAL PUs = 228

# Calculate total cost of solution
eval_cost_summary(cp1.1, s1.1x)[[2]] - pa_cost # TOTAL COST = 19.88

# Calculate how well feature representation targets are met by a solution and 
# the proportion of species' distributions covered
tcs1.1 <- eval_target_coverage_summary(cp1.1, s1.1x) # 100% of targets for threatened species met or exceeded
# Calculate overall species representation as the mean proportion of a species' range covered by the network
mean(tcs1.1$relative_held) # 49.82%

# Calculate total area of PA network (including existing PAs)
sum(st_area(s1.1[which(s1.1x$solution_1 == 1),]))/1000000 # Approx. 487228 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp1.1, s1.1x) # Boundary length = 29541.7 km

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
  mutate(prop_selected_1.1 = ifelse(
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
s1.1_plot <-
ggplot(data = s1.1_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_viridis(discrete = T, alpha = 0.85, begin = 0.1, option = "D", 
                    labels = c("Protected Area", 
                               "Indigenous Protected Area", 
                               "Freehold - Indigenous", 
                               "Native Title land", 
                               "Native Title land - Pastoral use",
                               "Indigenous Land Use Agreement",
                               "Pastoral term or perpetual lease")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "none",
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
  #add_locked_out_constraints("locked_out_ob1") %>% # lock out freehold land etc.
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s1.2 = solve(cp1.2)
s1.2 <- st_as_sf(s1.2) # Convert to sf 
s1.2x <- s1.2[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp1.2, s1.2x)[[2]] - pa_pus # TOTAL PUs = 805

# Calculate total cost of solution
eval_cost_summary(cp1.2, s1.2x)[[2]] - pa_cost # TOTAL COST = 104.90

# Calculate how well feature representation targets are met by a solution and the proportion of species' distributions covered
tcs1.2 <- eval_target_coverage_summary(cp1.2, s1.2x)
# Calculate overall species representation as the mean proportion of a species' range covered by the network
mean(tcs1.2$relative_held) # 53.65%

# Calculate total area of PA network (including existing PAs)
sum(st_area(s1.2[which(s1.2x$solution_1 == 1),]))/1000000 # Approx. 571393 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp1.2, s1.2x[,1]) # Boundary length = 33,758 km

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
  mutate(prop_selected_1.2 = ifelse(
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
                                                        "Pastoral term or perpetual lease",
                                                        "Other freehold, term, perpetual lease or Crown purposes"))

# Plot optimal solution
ggplot(data = s1.2_df) +
  geom_sf(mapping = aes(fill = status), color = "transparent") +
  scale_fill_manual(values = c("Not Protected" = "#af8dc3", "Protected" = "#7fbf7b"),
                    labels = c("Selected; Not Protected",
                               "Protected Area")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 0.8) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 11),
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
s1.2_plot <-
ggplot(data = s1.2_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_viridis(discrete = T, alpha = 0.85, begin = 0.1, option = "D", 
                     labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land",
                                "Native Title land - Pastoral use",
                                "Indigenous Land Use Agreement",
                                "Indigenous Land Use Agreement - Pastoral use",
                                "Pastoral term or perpetual lease",
                                "Other freehold, term, perpetual lease or Crown purposes")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        #legend.position = "none",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp1_2.png", units="cm", width=30, height=15, dpi=300, path = "results/fig/scp", bg  = 'white')



#### OBJECTIVE 2 (expanding PA network from only Indigenous land); CRITERIA 1 (threatened species) ####

cp2.1 = problem(pu_dat, species, cost_column = "cost") %>% # Use cost adjusted to reduce cost of declaring Native Title/Indigenous Freehold land
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_1) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_locked_out_constraints("locked_out_ob2") %>%
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0 will result in the solver stopping only when it has found the optimal solution

# Generate the optimal solution
s2.1 = solve(cp2.1)
s2.1 <- st_as_sf(s2.1) 
s2.1x <- s2.1[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp2.1, s2.1x)[[2]] - pa_pus # TOTAL PUs = 228

# Calculate total cost of solution
eval_cost_summary(cp2.1, s2.1x)[[2]] - pa_cost # TOTAL COST = 20.32826

# Calculate how well feature representation targets are met by a solution and the proportion of species' distributions covered
tcs2.1 <- eval_target_coverage_summary(cp2.1, s2.1x[,1])
# Calculate overall species representation as the mean proportion of a species' range covered by the network
mean(tcs2.1$relative_held) # 49.80%

# Calculate total area of PA network (including existing PAs)
sum(st_area(s2.1[which(s2.1x$solution_1 == 1),]))/1000000 # Approx. 487232.2 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp2.1, s2.1x[,1]) # Boundary length = 28,891 km

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
  mutate(prop_selected_2.1 = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order
s2.1_df$category <- factor(s2.1_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Native Title land - Pastoral use"))

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
s2.1_plot <-
ggplot(data = s2.1_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_viridis(discrete = T, alpha = 0.85, begin = 0.1, option = "D", 
                     labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land", 
                                "Native Title land - Pastoral use")) +
  geom_sf(data = locked_out_ob2, fill = "gray70", color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

ggsave("04_sp2_1.png", units="cm", width=30, height=15, dpi=300, path = "results/fig/scp", bg  = 'white')


#### OBJECTIVE 2 (expanding PA network from only Indigenous land); CRITERIA 2 (positive latent risk) ####

cp2.2 = problem(pu_dat, species, cost_column = "cost") %>%
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_2) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  add_locked_out_constraints("locked_out_ob2") %>%
  add_binary_decisions() %>% 
  add_gurobi_solver(gap = 0) 

# Generate the optimal solution
s2.2 = solve(cp2.2)
s2.2 <- st_as_sf(s2.2) 
s2.2x <- s2.2[, "solution_1", drop = FALSE]

# Calculate number of selected planning units additional to existing PA network
eval_n_summary(cp2.2, s2.2x)[[2]] - pa_pus # TOTAL PUs = 853

# Calculate total cost of solution
eval_cost_summary(cp2.2, s2.2x)[[2]] - pa_cost # TOTAL COST = 128

# Calculate how well feature representation targets are met by a solution and the proportion of species' distributions covered
tcs2.2 <- eval_target_coverage_summary(cp2.2, s2.2x[,1])
# Calculate overall species representation as the mean proportion of a species' range covered by the network
mean(tcs2.2$relative_held) # 55.26%

# Calculate total area of PA network (including existing PAs)
sum(st_area(s2.2[which(s2.2x$solution_1 == 1),]))/1000000 # Approx. 577524 km^2

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp2.2, s2.2x[,1]) # Boundary length = 33,402 km

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
  mutate(prop_selected_2.2 = ifelse(
    category %in% c("Protected Area", "Indigenous Protected Area"), 
    NA, 
    Count / sum(Count[!(category %in% c("Protected Area", "Indigenous Protected Area"))])))

# Convert 'land_type' column to a factor with levels in the specific order
s2.2_df$category <- factor(s2.2_df$category, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Native Title land - Pastoral use"))

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
s2.2_plot 
  
ggplot(data = s2.2_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_viridis(discrete = T, alpha = 0.85, begin = 0.1, end = 0.75, option = "D", 
                     labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land", 
                                "Native Title land - Pastoral use")) +
  geom_sf(data = locked_out_ob2, fill = "gray70", color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        #legend.position = "none",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

#### CREATE PATCHWORK OF PLOTS SHOWING OPTIMAL SOLUTIONS ####

# Plot optimal solution
s1.1_plot <-
  ggplot(data = s1.1_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land", 
                                "Native Title land - Pastoral use",
                                "Indigenous Land Use Agreement",
                                "Pastoral term or perpetual lease"),
                     values = c("#433e85",
                                "#38588c",
                                "#25858e",
                                "#2ab07f",
                                "#52c569",
                                "#86d549",
                                "#fde725")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

s1.2_plot <-
  ggplot(data = s1.2_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land",
                                "Native Title land - Pastoral use",
                                "Indigenous Land Use Agreement",
                                "Indigenous Land Use Agreement - Pastoral use",
                                "Pastoral term or perpetual lease",
                                "Other freehold, term, perpetual lease or Crown purposes"),
                    values = c("#433e85",
                               "#38588c",
                               "#25858e",
                               "#2ab07f",
                               "#52c569",
                               "#86d549",
                               "#c2df23",
                               "#fde725",
                               "#fb9f3a")) +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

# Plot optimal solution
s2.1_plot <-
  ggplot(data = s2.1_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land", 
                                "Native Title land - Pastoral use"),
                    values = c("#433e85",
                               "#38588c",
                               "#25858e",
                               "#2ab07f",
                               "#52c569")) +
  geom_sf(data = locked_out_ob2, fill = "gray70", color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

# Plot optimal solution
s2.2_plot <-
ggplot(data = s2.2_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(labels = c("Protected Area", 
                                "Indigenous Protected Area", 
                                "Freehold - Indigenous", 
                                "Native Title land", 
                                "Native Title land - Pastoral use"),
                     values = c("#433e85",
                                "#38588c",
                                "#25858e",
                                "#2ab07f",
                                "#52c569")) +
  geom_sf(data = locked_out_ob2, fill = "gray70", color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-30))

legend_plot <-
  ggplot(data = s1.2_df) +
  geom_sf(mapping = aes(fill = factor(category)), color = "transparent") +
  scale_fill_manual(labels = c("Protected area", 
                               "Indigenous protected area", 
                               "Indigenous freehold", 
                               "Native title land",
                               "Native title land - pastoral use",
                               "Indigenous land use agreement",
                               "Indigenous land use agreement - pastoral use",
                               "Pastoral term or perpetual lease",
                               "Other freehold, term, perpetual lease or Crown purposes",
                               ""),
                    values = c("#433e85",
                               "#38588c",
                               "#25858e",
                               "#2ab07f",
                               "#52c569",
                               "#86d549",
                               "#c2df23",
                               "#fde725",
                               "#fb9f3a")) +
  theme_classic() +
  theme(legend.text = element_text(size = 11),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0))
# Extract the legend using cowplot
legend <- cowplot::get_legend(legend_plot)
# Plot the legend
grid::grid.newpage()
grid::grid.draw(legend)

ggsave("04_solutions_legend.png", legend, units="cm", width=13, height=8, dpi=300, path = "results/fig/scp", bg  = 'white')

solutions_plot <- s1.1_plot + s1.2_plot + s2.1_plot + s2.2_plot +
  plot_annotation(tag_levels = 'a') 

# Print the combined plot
print(solutions_plot)

ggsave("04_optimal_solution_plots.png", units="cm", width=30, height=15, dpi=300, path = "results/fig/scp", bg  = 'white')

#### CALCULATE PROPORTION OF SPECIES' AMT DISTRIBUTION COVERED BY AMT PROTECTED AREA NETWORK ####

# Load shapefile of IUCN distribution polygons for AMT mammals with threat status data generated in 02_pgls_model
# Note: iucn_threat contains geometry for entire species distribution
iucn_threat <-  read_sf("data/output-data/shp/02_iucn_threat.shp")
colnames(iucn_threat)[1:4] <- c("scientificName", "ordinalThreat", "mean_fitted", "latent_risk")
iucn_amt <- iucn_threat %>% 
  st_intersection(amt)

# Load CAPAD protected areas cropped to AMT boundary
capad_amt <- read_sf("data/output-data/shp/03_capad_amt.shp")
capad_union <- st_union(capad_amt)

ggplot() +
  geom_sf(data = amt, color = "black", alpha = 0.5) +
  geom_sf(data = capad_amt, color = "black", fill = 'lightgreen', alpha = 0.5) +
  geom_sf(data = iucn_amt[128,], fill = "red", color = "black", alpha = 0.5)

# Calculate area of species' AMT distribution that is covered by existing AMT PA network
spec_dat$area_pa <- NA
spec_dat$prop_protected <- NA
for (i in 1:length(iucn_amt$scientificName)) {
  curr_shape <- iucn_amt[i, "geometry"]
  curr_intersect <- st_intersection(curr_shape, capad_union)
  curr_area <- st_area(curr_intersect)
  spec_dat$area_pa[i] <- curr_area
}

# Calculate the proportion of each species' AMT distributions that covered by AMT PA network
spec_dat <- spec_dat %>%
  mutate(prop_protected = area_pa/area_amt)

mean(spec_dat$prop_protected)*100


#### PLOT THE SELECTION FREQUENCY OF PUs ACROSS ALL SCENARIOS ####

# Assuming geometry is the same across all data frames
geometry <- s1.1x$geometry

# Combine and rename columns without geometry
sel_freq <- bind_cols(
  data.frame(s1.1 = s1.1$solution_1),
  data.frame(s1.2 = s1.2$solution_1),
  data.frame(s2.1 = s2.1$solution_1),
  data.frame(s2.2 = s2.2$solution_1)
)

# Add the geometry column and row sums
sel_freq <- sel_freq %>%
  mutate(geometry = geometry) %>%
  mutate(row_sum = rowSums(across(c(s1.1, s1.2, s2.1, s2.2)))) %>%
  filter(!is.na(row_sum)) %>%
  st_as_sf() 


ggplot(data = sel_freq) +
  geom_sf(mapping = aes(fill = as.factor(row_sum)), color = "transparent") +
  scale_fill_manual(values = c("0" = "white", "1" = "#25858e", "2" = "#18be8f", "3" = "#86d549", "4" = "#e1e920")) +
  geom_sf(data = protected_area, fill = "gray70", color = "transparent") +
  geom_sf(data = amt, fill = "transparent", color = "black", size = 1) +
  geom_sf(data = wt, fill = "gray10", color = "black", size = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 10),
        legend.key = element_rect(colour = "gray40"),
        legend.title = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,-10,-20))

ggsave("04_selection_frequency.png", units="cm", width=18, height=8, dpi=300, path = "results/fig/scp", bg  = 'white')

#### PLOT PROPORTION OF LAND TENURE IN PUS ####

s1.1_tenure <- s1.1_df %>%
  bind_rows(data.frame(category = "Indigenous Land Use Agreement - Pastoral use", prop_selected_1.1 = NA)) %>%
  bind_rows(data.frame(category = "Other freehold, term, perpetual lease or Crown purposes", prop_selected_1.1 = NA)) %>%
  as.data.frame() %>%
  dplyr::select(1,5) 

s1.2_tenure <- s1.2_df%>%
  as.data.frame() %>%
  dplyr::select(1,5) 
  
s2.1_tenure <- s2.1_df %>%
  bind_rows(data.frame(category = "Indigenous Land Use Agreement", prop_selected_2.1 = NA)) %>%
  bind_rows(data.frame(category = "Indigenous Land Use Agreement - Pastoral use", prop_selected_2.1 = NA)) %>%
  bind_rows(data.frame(category = "Other freehold, term, perpetual lease or Crown purposes", prop_selected_2.1 = NA)) %>%
  bind_rows(data.frame(category = "Pastoral term or perpetual lease", prop_selected_2.1 = NA)) %>%
  as.data.frame() %>%
  dplyr::select(1,5) 
  
s2.2_tenure <- s2.2_df %>%
  bind_rows(data.frame(category = "Indigenous Land Use Agreement", prop_selected_2.2 = NA)) %>%
  bind_rows(data.frame(category = "Indigenous Land Use Agreement - Pastoral use", prop_selected_2.2 = NA)) %>%
  bind_rows(data.frame(category = "Other freehold, term, perpetual lease or Crown purposes", prop_selected_2.2 = NA)) %>%
  bind_rows(data.frame(category = "Pastoral term or perpetual lease", prop_selected_2.2 = NA)) %>%
  as.data.frame() %>%
  dplyr::select(1,5) 

combined_tenure <- s1.1_tenure %>%
  left_join(s1.2_tenure, by = "category") %>%
  left_join(s2.1_tenure, by = "category") %>%
  left_join(s2.2_tenure, by = "category") 


# Reshape data into long format for plotting
combined_long <- combined_tenure %>%
  pivot_longer(cols = starts_with("prop_selected"),
               names_to = "dataset",
               values_to = "prop_selected")

# Create objective column
combined_long <- combined_long %>%
  mutate(objective = case_when(
    grepl("1.1|1.2", dataset) ~ "Objective 1",
    grepl("2.1|2.2", dataset) ~ "Objective 2",
    TRUE ~ NA_character_
  ))

# Create criteria column
combined_long <- combined_long %>%
  mutate(criteria = case_when(
    grepl("1.1|2.1", dataset) ~ "Current risk",
    grepl("1.2|2.2", dataset) ~ "Latent risk",
    TRUE ~ NA_character_
  )) %>%
  mutate(prop_selected = replace(prop_selected, is.na(prop_selected), 0)) %>%
  filter(!(category %in% c("Protected Area", "Indigenous Protected Area")))

# Convert 'land_type' column to a factor with levels in the specific order
combined_long$category <- factor(combined_long$category, levels = c("Freehold - Indigenous", 
                                                        "Native Title land", 
                                                        "Native Title land - Pastoral use",
                                                        "Indigenous Land Use Agreement",
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease",
                                                        "Other freehold, term, perpetual lease or Crown purposes"))

# Plot using ggplot2
ggplot(combined_long, aes(x = criteria, y = prop_selected, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(. ~ objective, switch="both") + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(name = "Land tenure",
                    labels = c("Indigenous freehold", 
                               "Native title land", 
                               "Native title land - pastoral use",
                               "Indigenous land use agreement",
                               "Indigenous land use agreement - pastoral use",
                               "Pastoral term or perpetual lease",
                               "Other freehold, term, perpetual lease or Crown purposes"),
                    values = c("#25858e",
                               "#2ab07f",
                               "#52c569",
                               "#86d549",
                               "#c2df23",
                               "#fde725",
                               "#fb9f3a")) +
  labs(y = "Proportion of selected PUs") +
  theme_classic() +
  theme(legend.text = element_text(size = 14),
        legend.key = element_rect(colour = "transparent"),
        legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14),
        legend.margin=margin(0,0,0,0))

ggsave("04_proportion_tenure_boxplot.png", units="cm", width=30, height=15, dpi=300, path = "results/fig/scp", bg  = 'white')
