# ==============================================================================
# Check species growth
# Author: Raquel
# Description: compare observed vs predicted growth
# ==============================================================================

library(tidyverse)  # For data manipulation and ggplot2
library(patchwork)  # For combining plots

# Define species list for consistent processing
species_info <- tribble(
  ~file_name,             ~species_name,          ~display_name,
  "Cod-size-at-age",      "atlantic_cod",         "Atlantic cod",
  "turbot-size-at-age",   "turbot",               "Turbot",
  "redfish-size-at-age",  "redfish",              "Redfish",
  "plaice-size-at-age",   "american_plaice",      "American plaice",
  "witch-size-at-age",    "witch_flounder",       "Witch flounder",
  "yellowtail-size-at-age", "yellowtail_flounder", "Yellowtail flounder",
  "capelin-size-at-age",  "capelin",              "Capelin"
)

# Function to read and preprocess data for a single species
read_species_data <- function(file_name, species_name) {
  data_path <- file.path("data/age-length-key_2005", paste0(file_name, ".csv"))
  
  df <- read.csv(data_path, header = TRUE) %>%
    na.omit() %>%
    mutate(species = species_name)
  
  # Special handling for capelin
  if (species_name == "capelin") {
    df <- df %>%
      mutate(length = length / 10) %>%
      filter(age < 7)
  }
  
  return(df)
}

# Read all species data
species_data_list <- map2(species_info$file_name, species_info$species_name, read_species_data)
names(species_data_list) <- species_info$species_name

# Combine all species data
all_species_data <- bind_rows(species_data_list)

# Save combined data
write_csv(all_species_data, file = "data/species_age_weight.csv")

## ---------------------------------------------------------------
## Load model parameters and calculate size from length
## ---------------------------------------------------------------

# Load model parameters 
#params1 <- readRDS("rds/params_adjusted_RL080.rds")
#params1 <- readRDS("rds/GB_model_optimized_Rmax_erepro_Z0.rds")
params1<- params_optim #params_adjusted #

params1<- readRDS("rds/params_adjusted_growth_FINAL.rds")

# Function to calculate size from length using model parameters
calculate_size <- function(df, species_name, params) {
  a <- params@species_params$a[params@species_params$species == species_name]
  b <- params@species_params$b[params@species_params$species == species_name]
  
  df %>% mutate(size = a * length^b)
}

# Apply size calculation to each species dataset
species_data_with_size <- map2(
  species_data_list,
  species_info$species_name,
  ~calculate_size(.x, .y, params1)
)

## ---------------------------------------------------------------
## Generate growth curves and create plots
## ---------------------------------------------------------------

# Create a function to get modeled growth curves and generate plot for each species
create_species_plot <- function(species_name, display_name, empirical_data) {
  # Get modeled growth data
  modeled_data <- as.data.frame(getGrowthCurves(params1, species = species_name))
  
  # Convert to long format for plotting
  modeled_long <- modeled_data %>%
    gather(age, size, factor_key = TRUE) %>%
    mutate(age = as.numeric(age))
  
  # Create plot
  ggplot() +
    geom_line(data = modeled_long, aes(x = age, y = size, group = 1), 
              color = "#D55E00", size = 1) +
    geom_point(data = empirical_data, aes(x = age, y = size), 
               color = "#0072B2", shape = 1, alpha = 0.6) +
    labs(title = display_name, x = "Age", y = "Size") +
    xlim(0, 20) +
    theme_minimal()
}

# Plot growth curves for all species using the model's built-in function
plotGrowthCurves(params1, species_panel = TRUE, max_age = 20)

# Create individual plots for each species
species_plots <- pmap(
  list(
    species_info$species_name, 
    species_info$display_name,
    species_data_with_size
  ),
  create_species_plot
)
names(species_plots) <- species_info$species_name

# Combine all plots using patchwork
combined_growth_plot <- wrap_plots(species_plots) & theme_publication()

# Save the final plot
ggsave(
  filename = "figures/FINAL/Growth_curves_params_FINAL.png",
  plot = combined_growth_plot,
  width = 22, height = 18, units = "cm"
)

