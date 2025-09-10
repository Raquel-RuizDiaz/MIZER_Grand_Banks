# ==============================================================================
# Time Series Validation Plots: Biomass and Catch
# Author: Raquel Ruiz Diaz
# Description: Create time series plots comparing observed vs predicted biomass and catch
# ==============================================================================

# Load required packages ----
library(tidyverse)
library(reshape2)
library(ggplot2)
library(patchwork)

# Configuration ----
CONFIG <- list(
  # File paths
  BIO_OBS_FILE = "data/biomass_historical.csv",
  CATCH_OBS_FILE = "data/Catches_historical_edit3.csv",
  EFFORT_FILE = "data/Fmort_ratio_historical_2025.csv",
  
  # Output settings
  OUTPUT_DIR = "figures/FINAL/",
  PLOT_WIDTH = 26,
  PLOT_HEIGHT = 20,
  PLOT_UNITS = "cm",
  PLOT_DPI = 300,
  
  # Plot layout
  FACET_NCOL = 3
)

# Species names
SPECIES_NAMES <- c(
  "american_plaice", "atlantic_cod", "redfish", "turbot",
  "witch_flounder", "yellowtail_flounder", "capelin",
  "sand_lance", "thorny_skate"
)

# Create output directory
dir.create(CONFIG$OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Custom theme ----
theme_timeseries <- function() {
  theme_bw() +
    theme(
      strip.background = element_rect(fill = "white"),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 9),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

# ==============================================================================
# DATA PREPARATION FUNCTIONS
# ==============================================================================

#' Prepare model simulation data
#' 
#' @param params_file Path to mizer parameters file
#' @param effort_file Path to effort data file
#' @return Mizer simulation object
prepare_simulation <- function(params_file = "rds/params_adjusted_growth_FINAL.rds",
                               effort_file = CONFIG$EFFORT_FILE) {
  
  cat("Loading model parameters...\n")
  params_multi <- readRDS(params_file)
  
  # Set up initial effort
  effort_initial <- params_multi@initial_effort
  names(effort_initial) <- SPECIES_NAMES
  
  # Run initial equilibrium simulation
  cat("Running equilibrium simulation...\n")
  sim_optim <- project(params_multi, effort = effort_initial, t_max = 500)
  
  # Read historical effort data
  cat("Loading historical effort data...\n")
  effort_historical <- read.csv(effort_file, sep = ";", row.names = 1)
  effort_historical <- as.matrix(effort_historical)
  
  # Run simulation with time-varying effort
  cat("Running historical simulation...\n")
  sim_historical <- project(
    params_multi,
    effort = effort_historical,
    initial_n = sim_optim@n[100,,],
    initial_npp = sim_optim@n_pp[100,]
  )
  
  return(sim_historical)
}

#' Convert model results to long format dataframe
#' 
#' @param model_data Matrix with years as rows, species as columns
#' @param data_type Character string describing the data type
#' @return Long format dataframe
convert_to_long_format <- function(model_data, data_type = "predicted") {
  
  as.data.frame(model_data) %>%
    tibble::rownames_to_column(var = "Year") %>%
    mutate(Year = as.numeric(Year)) %>%
    reshape2::melt(
      id.vars = "Year",
      variable.name = "Species",
      value.name = "value"
    ) %>%
    mutate(source = data_type)
}

#' Load and process observed data
#' 
#' @param file_path Path to observed data file
#' @param year_col Name of year column in the data
#' @return Long format dataframe with observed data
load_observed_data <- function(file_path, year_col = "survey.year") {
  
  if (!file.exists(file_path)) {
    stop(paste("Data file not found:", file_path))
  }
  
  cat("Loading observed data from:", file_path, "\n")
  
  obs_data <- read.csv(file_path, sep = ";", stringsAsFactors = FALSE)
  
  # Handle different year column names
  if (year_col %in% colnames(obs_data)) {
    obs_data <- obs_data %>% rename(Year = !!year_col)
  } else if ("year" %in% colnames(obs_data)) {
    obs_data <- obs_data %>% rename(Year = year)
  }
  
  # Convert to long format
  obs_data %>%
    mutate(Year = as.numeric(Year)) %>%
    reshape2::melt(
      id.vars = "Year",
      variable.name = "Species",
      value.name = "value"
    ) %>%
    mutate(source = "observed") %>%
    filter(!is.na(value))
}

# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================

#' Create time series comparison plot
#' 
#' @param combined_data Combined observed and predicted data
#' @param title Plot title
#' @param y_label Y-axis label
#' @return ggplot object
create_timeseries_plot <- function(combined_data, title, y_label) {
  
  # Filter species with observed data
  species_with_obs <- combined_data %>%
    filter(source == "observed") %>%
    group_by(Species) %>%
    summarise(has_data = sum(!is.na(value)) > 0, .groups = "drop") %>%
    filter(has_data) %>%
    pull(Species)
  
  plot_data <- combined_data %>% 
    filter(Species %in% species_with_obs) %>%
    mutate(Species = as.character(Species))
  
  # Create plot
  ggplot(plot_data, aes(x = Year, y = value, color = source, linetype = source)) +
    geom_line(size = 1) +
    geom_point(
      data = filter(plot_data, source == "observed"), 
      size = 1.5,
      alpha = 0.8
    ) +
    facet_wrap(~Species, scales = "free_y", ncol = CONFIG$FACET_NCOL) +
    scale_color_manual(
      values = c("observed" = "black", "predicted" = "red"),
      labels = c("observed" = "Observed", "predicted" = "Predicted")
    ) +
    scale_linetype_manual(
      values = c("observed" = "solid", "predicted" = "dashed"),
      labels = c("observed" = "Observed", "predicted" = "Predicted")
    ) +
    labs(
      title = title,
      x = "Year",
      y = y_label,
      color = "Data Source",
      linetype = "Data Source"
    ) +
    theme_timeseries()
}

# ==============================================================================
# MAIN ANALYSIS FUNCTION
# ==============================================================================

#' Create biomass and catch time series validation plots
#' 
#' @param sim_results Mizer simulation results
#' @param bio_obs_file Path to observed biomass data
#' @param catch_obs_file Path to observed catch data
#' @param save_plots Logical, whether to save plots
#' @return List with biomass and catch plots
create_validation_timeseries <- function(sim_results,
                                         bio_obs_file = CONFIG$BIO_OBS_FILE,
                                         catch_obs_file = CONFIG$CATCH_OBS_FILE,
                                         save_plots = TRUE) {
  
  cat("=== CREATING TIME SERIES VALIDATION PLOTS ===\n\n")
  
  # Extract model predictions
  cat("Extracting model predictions...\n")
  predicted_biomass <- getBiomass(sim_results)
  predicted_catch <- getYield(sim_results)
  
  # Convert model data to long format
  pred_bio_long <- convert_to_long_format(predicted_biomass, "predicted")
  pred_catch_long <- convert_to_long_format(predicted_catch, "predicted")
  
  # Load observed data
  cat("Loading observed data...\n")
  obs_bio_long <- load_observed_data(bio_obs_file, "survey.year")
  obs_catch_long <- load_observed_data(catch_obs_file, "year")
  
  # Combine observed and predicted data
  biomass_combined <- rbind(obs_bio_long, pred_bio_long)
  catch_combined <- rbind(obs_catch_long, pred_catch_long)
  
  # Create plots
  cat("Creating biomass time series plot...\n")
  biomass_plot <- create_timeseries_plot(
    biomass_combined,
    "Biomass Over Time: Observed vs Predicted",
    "Biomass (kg)"
  )
  
  cat("Creating catch time series plot...\n")
  catch_plot <- create_timeseries_plot(
    catch_combined,
    "Catches Over Time: Observed vs Predicted",
    "Catch (kg)"
  )
  
  # Save plots if requested
  if (save_plots) {
    cat("Saving plots...\n")
    
    # Save biomass plot
    bio_filename <- file.path(CONFIG$OUTPUT_DIR, "biomass_time_series.png")
    ggsave(
      filename = bio_filename,
      plot = biomass_plot,
      width = CONFIG$PLOT_WIDTH,
      height = CONFIG$PLOT_HEIGHT,
      units = CONFIG$PLOT_UNITS,
      dpi = CONFIG$PLOT_DPI
    )
    
    # Save catch plot
    catch_filename <- file.path(CONFIG$OUTPUT_DIR, "catch_time_series.png")
    ggsave(
      filename = catch_filename,
      plot = catch_plot,
      width = CONFIG$PLOT_WIDTH,
      height = CONFIG$PLOT_HEIGHT,
      units = CONFIG$PLOT_UNITS,
      dpi = CONFIG$PLOT_DPI
    )
    
    # Save PDF versions for publication
    ggsave(
      filename = file.path(CONFIG$OUTPUT_DIR, "biomass_time_series.pdf"),
      plot = biomass_plot,
      width = CONFIG$PLOT_WIDTH,
      height = CONFIG$PLOT_HEIGHT,
      units = CONFIG$PLOT_UNITS
    )
    
    ggsave(
      filename = file.path(CONFIG$OUTPUT_DIR, "catch_time_series.pdf"),
      plot = catch_plot,
      width = CONFIG$PLOT_WIDTH,
      height = CONFIG$PLOT_HEIGHT,
      units = CONFIG$PLOT_UNITS
    )
    
    cat("Plots saved to:", CONFIG$OUTPUT_DIR, "\n")
    cat("  - biomass_time_series.png/pdf\n")
    cat("  - catch_time_series.png/pdf\n")
  }
  
  # Print summary
  cat("\n=== SUMMARY ===\n")
  cat("Biomass data - Species with observations:", 
      length(unique(biomass_combined$Species[biomass_combined$source == "observed"])), "\n")
  cat("Catch data - Species with observations:", 
      length(unique(catch_combined$Species[catch_combined$source == "observed"])), "\n")
  cat("Time series validation plots created successfully!\n")
  
  return(list(
    biomass_time = biomass_plot,
    catch_time = catch_plot,
    data = list(
      biomass_combined = biomass_combined,
      catch_combined = catch_combined
    )
  ))
}

# ==============================================================================
# EXAMPLE USAGE
# ==============================================================================

# # Run the complete analysis
sim_results <- prepare_simulation()
validation_plots <- create_validation_timeseries(sim_results)

# # View the plots
validation_plots$biomass_time
validation_plots$catch_time
# 
# # Access the underlying data if needed
biomass_data <- validation_plots$data$biomass_combined
catch_data <- validation_plots$data$catch_combined