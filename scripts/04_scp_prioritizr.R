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

library(prioritizr)
library(gurobi)
library(highs)

# Load spec_dat table of conservation features
spec_dat <- read_csv("data/output-data/tbl/03_spec_dat.csv")

# Export the raster stack of species distributions (i.e. the 'conservation features' as a TIFF file)
species <- rast("data/output-data/tif/03_species_stack.tif")

# Load shapefile of hexagonal planning unit tesselation for the AMT
pu_dat <- st_read("data/output-data/shp/03_amt_pu_grid.shp")
colnames(pu_dat)[5:18] <- c("locked_in", "locked_out", "cost", "main_id", 
                            "main_type", "category", "declare_ipa", "median_income",
                            "hii_oceania", "normal_hii", "normal_mei", "cost_adjusted",
                            "id_spec_richness", "spec_richness")
# Convert locked_in column back to logical variable
pu_dat$locked_in <- as.logical(pu_dat$locked_in)

# Load raster stack of species richness (0-40 and 41-62)
richness_stack <- rast("data/output-data/tif/03_richness_stack.tif")

# Load data frame containing species richness layer attributes
richness_df <- read_csv("data/output-data/tbl/03_richness_df.csv")

#### OBJECTIVE 1: EXPAND EXISTING PROTECTED AREA NETWORK; NO RESTRICTIONS ####

cp1.1 = problem(pu_dat, species, cost_column = "cost") %>%
  add_min_set_objective() %>% # sets same minimum objective
  add_relative_targets(spec_dat$target_1) %>% # specifies relative targets for each feature
  add_locked_in_constraints("locked_in") %>% # lock in current protected area network (PAN)
  #add_locked_out_constraints("locked_out") %>% # nothing locked out 
  add_binary_decisions() %>% # solution produced involves either a 'selected' planning unit or 'unselected' planning unit
  add_gurobi_solver(gap = 0)  # value of 0.01 will result in the solver stopping when it has found a solution that is optimal to 1%
  #add_boundary_penalties(penalty=0.01) # this time use a more reasonable boundary penalty

# Generate the optimal solution
s1.1 = solve(cp1.1)
s1.1 <- st_as_sf(s1.1) # Convert to sf 
s1.1x <- s1.1[, "solution_1", drop = FALSE]

# PU SUMMARY
# Calculate number of PUs in the PAN
pa_pus <- length(pu_dat$status[which(pu_dat$status==1)]) # 3,410 PUs in PAN
# Calculate number of selected planning units additional to PAN
eval_n_summary(cp1.1, s1.1x)[[2]] - pa_pus # TOTAL PUs = 310 vs 5904 with boundary penalty (bp)

# COST SUMMARY
# Calculate the total cost of a solution.
pa_cost <- sum(hex_grid_sf_v4$cost[which(hex_grid_sf_v4$status == 1)], na.rm = T) # cost of area in PAN = 576.77
# Calculate total cost of solution
eval_cost_summary(cp1.1, s1.1x)[[2]] - pa_cost # TOTAL COST = 172.23 vs 1221.969 with bp

# TARGET COVERAGE SUMMARY
# Calculate how well feature representation targets are met by a solution and proportion of species' distribution covered
tcs1.1 <- eval_target_coverage_summary(cp1.1, s1.1x)
# 100% of priority species have target met by solution

# AREA SUMMARY
# Calculate total area of reserve network
sum(st_area(s1.1[which(s1.1x$solution_1 == 1),]))

# Calculate the exposed boundary length (perimeter) associated with a solution
eval_boundary_summary(cp1.1, s1.1x) # Boundary length = 32,657,849 vs 16,375,270 with bp

# Calculate the number and proportion of selected PUs in each land use category
s1.1_df <- s1.1 %>%
  filter(solution_1 == 1) %>%
  #filter(!(CATEGORY == 'Protected Area')) %>%
  #filter(!CATEGORY == 'Indigenous Protected Area') %>%
  group_by(CATEGORY) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

# Convert 'land_type' column to a factor with levels in the specific order
s1.1_df$CATEGORY <- factor(s1.1_df$CATEGORY, levels = c("Protected Area", 
                                                        "Indigenous Protected Area", 
                                                        "Freehold - Indigenous", 
                                                        "Native Title Land", 
                                                        "Native Title Land - Pastoral use", 
                                                        "Indigenous Land Use Agreement", 
                                                        "Indigenous Land Use Agreement - Pastoral use",
                                                        "Pastoral term or perpetual lease", 
                                                        "Other freehold, term, perpetual lease or Crown purposes"))

# Plot optimal solution
ggplot(data = s1.1_df) +
  geom_sf(mapping = aes(fill = factor(CATEGORY)), color = "transparent") +
  scale_fill_manual(values = c("#003c30",
                               "#01665e", 
                               "#35978f", 
                               "#80cdc1", 
                               "#c7eae5", 
                               "#f6e8c3", 
                               "#dfc27d", 
                               "#bf812d", 
                               "#8c510a",
                               "#b2182b"),
                    labels = c("Protected Area", 
                               "Indigenous Protected Area", 
                               "Freehold - Indigenous", 
                               "Native Title Land", 
                               "Native Title Land - Pastoral use", 
                               "Indigenous Land Use Agreement", 
                               "Indigenous Land Use Agreement - Pastoral use",
                               "Pastoral term or perpetual lease", 
                               "Other freehold, term, perpetual lease or Crown purposes")) +
  geom_sf(data = mines, fill = "#b2182b", color = "transparent", size = 0.4) +
  geom_sf(data = amt, fill = "transparent", color = "gray20", size = 0.4) +
  geom_sf(data = wt, fill = "#6d6657", color = "#6d6657", size = 0.4) +
  theme_minimal() +
  theme(legend.text = element_text(size = 12),
        legend.key = element_rect(colour = "black"),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

#ggsave("sp1.1.tiff", units="cm", width=20, height=15, dpi=300, compression = 'lzw',
#       path = "C:/Users/u7364298/OneDrive - Australian National University/Figures/SCP Maps/Current Risk")





