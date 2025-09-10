# ==============================================================================
# SCENARIO SIMULATION - COD, CAPELIN, AND SAND LANCE
# Author: Raquel
# Description: 300-year simulations with 10-year averaging for publication figures
# ==============================================================================

# Load required libraries
library(mizer)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(tibble)
library(purrr)
library(stringr)
library(mizerExperimental)

# =============================================================================
# SETUP AND CONFIGURATION
# =============================================================================

# Species configuration
species_names <- c(
  "american_plaice", "atlantic_cod", "redfish", "turbot",
  "witch_flounder", "yellowtail_flounder", "capelin",
  "sand_lance", "thorny_skate"
)

# Key species indices
COD_IDX <- 2
CAPELIN_IDX <- 7
SANDLANCE_IDX <- 8

# Effort levels for scenarios
effort_levels <- c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0)

# Load calibrated model parameters
params_multi <- readRDS("rds/params_adjusted_growth_FINAL.rds")
params_baseline <- params_multi
effort_baseline <- params_multi@initial_effort
names(effort_baseline) <- species_names

# =============================================================================
# BASELINE SIMULATION (300 years)
# =============================================================================

cat("Running 300-year baseline simulation...\n")
sim_baseline <- project(
  params_baseline,
  effort = effort_baseline,
  dt = 0.1,
  t_max = 300
)

# Extract equilibrium conditions from end of baseline
initial_n <- sim_baseline@n[nrow(sim_baseline@n),,]
initial_n_pp <- sim_baseline@n_pp[nrow(sim_baseline@n_pp),]
initial_effort <- sim_baseline@effort[dim(sim_baseline@effort)[1],]

# =============================================================================
# SCENARIO CREATION FUNCTIONS
# =============================================================================

# Standard scenarios for cod and capelin (multiplier-based)
create_standard_scenarios <- function(baseline_effort, target_species_idx,
                                      effort_multipliers = effort_levels) {
  scenarios <- list()
  scenario_names <- paste0(names(baseline_effort)[target_species_idx],
                           "_F", effort_multipliers)
  
  for (i in seq_along(effort_multipliers)) {
    new_effort <- baseline_effort
    new_effort[target_species_idx] <- baseline_effort[target_species_idx] * effort_multipliers[i]
    scenarios[[i]] <- new_effort
  }
  
  names(scenarios) <- scenario_names
  return(scenarios)
}

# Specialized scenarios for sand lance (specific effort levels)
create_sandlance_scenarios <- function(baseline_effort,
                                       specific_efforts = c(0.14, 0.18, 0.28),
                                       effort_labels = c("1.5", "2.0", "3.0")) {
  scenarios <- list()
  scenario_names <- paste0("sand_lance_F", effort_labels)
  
  for (i in seq_along(specific_efforts)) {
    new_effort <- baseline_effort
    new_effort[SANDLANCE_IDX] <- specific_efforts[i]
    scenarios[[i]] <- new_effort
  }
  
  names(scenarios) <- scenario_names
  return(scenarios)
}

# =============================================================================
# CREATE ALL SCENARIOS
# =============================================================================

# Generate scenarios for each target species
cod_scenarios <- create_standard_scenarios(initial_effort, COD_IDX)
capelin_scenarios <- create_standard_scenarios(initial_effort, CAPELIN_IDX)
sandlance_scenarios <- create_sandlance_scenarios(initial_effort)

# Combine all scenarios
all_scenarios <- c(cod_scenarios, capelin_scenarios, sandlance_scenarios)

# =============================================================================
# SIMULATION EXECUTION
# =============================================================================

run_scenario_simulations <- function(scenarios, params, initial_conditions,
                                     t_max = 300, dt = 0.1) {
  results <- list()
  
  for (i in seq_along(scenarios)) {
    cat("Running scenario", i, "of", length(scenarios), ":", names(scenarios)[i], "\n")
    
    results[[i]] <- project(
      params,
      effort = scenarios[[i]],
      dt = dt,
      t_max = t_max,
      initial_n = initial_conditions$n,
      initial_n_pp = initial_conditions$n_pp
    )
  }
  
  names(results) <- names(scenarios)
  return(results)
}

# Run all simulations
initial_conditions <- list(n = initial_n, n_pp = initial_n_pp)

# Running scenario simulations
sim_results <- run_scenario_simulations(all_scenarios, params_baseline, initial_conditions)

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

# Calculate comprehensive changes (biomass, feeding, predation mortality)
calculate_comprehensive_changes <- function(scenario_sim, baseline_sim,
                                            averaging_period = c(291, 300)) {
  
  # 1. Biomass changes
  scenario_biomass <- as.data.frame(getBiomass(scenario_sim))
  baseline_biomass <- as.data.frame(getBiomass(baseline_sim))
  
  scenario_bio_mean <- colMeans(scenario_biomass[averaging_period[1]:averaging_period[2], ])
  baseline_bio_mean <- colMeans(baseline_biomass[averaging_period[1]:averaging_period[2], ])
  
  biomass_change <- (scenario_bio_mean - baseline_bio_mean) / baseline_bio_mean * 100
  
  # 2. Feeding level changes
  scenario_fl <- getFeedingLevel(scenario_sim)
  baseline_fl <- getFeedingLevel(baseline_sim)
  
  # Handle array dimensions consistently
  if (is.array(scenario_fl) && length(dim(scenario_fl)) == 3) {
    scenario_fl_species <- apply(scenario_fl, c(1, 2), mean, na.rm = TRUE)
    baseline_fl_species <- apply(baseline_fl, c(1, 2), mean, na.rm = TRUE)
    
    scenario_fl_mean <- colMeans(scenario_fl_species[averaging_period[1]:averaging_period[2], , drop = FALSE])
    baseline_fl_mean <- colMeans(baseline_fl_species[averaging_period[1]:averaging_period[2], , drop = FALSE])
  } else {
    scenario_fl_mean <- colMeans(scenario_fl[averaging_period[1]:averaging_period[2], , drop = FALSE])
    baseline_fl_mean <- colMeans(baseline_fl[averaging_period[1]:averaging_period[2], , drop = FALSE])
  }
  
  feeding_change <- (scenario_fl_mean - baseline_fl_mean) / baseline_fl_mean * 100
  
  # 3. Predation mortality changes
  scenario_mort <- getPredMort(scenario_sim)
  baseline_mort <- getPredMort(baseline_sim)
  
  # Handle array dimensions consistently
  if (is.array(scenario_mort) && length(dim(scenario_mort)) == 3) {
    scenario_mort_species <- apply(scenario_mort, c(1, 2), mean, na.rm = TRUE)
    baseline_mort_species <- apply(baseline_mort, c(1, 2), mean, na.rm = TRUE)
    
    scenario_mort_mean <- colMeans(scenario_mort_species[averaging_period[1]:averaging_period[2], , drop = FALSE])
    baseline_mort_mean <- colMeans(baseline_mort_species[averaging_period[1]:averaging_period[2], , drop = FALSE])
  } else {
    scenario_mort_mean <- colMeans(scenario_mort[averaging_period[1]:averaging_period[2], , drop = FALSE])
    baseline_mort_mean <- colMeans(baseline_mort[averaging_period[1]:averaging_period[2], , drop = FALSE])
  }
  
  predmort_change <- (scenario_mort_mean - baseline_mort_mean) / baseline_mort_mean * 100
  
  # Combine results
  results <- data.frame(
    species = names(biomass_change),
    biomass_change = biomass_change,
    feeding_change = feeding_change,
    predmort_change = predmort_change,
    stringsAsFactors = FALSE
  )
  
  return(results)
}

# MISSING FUNCTION: System-level metrics calculation
# System-level metrics calculation - BIOMASS ONLY
calculate_system_metrics <- function(scenario_sim, baseline_sim, averaging_period = c(291, 300)) {
  
  # Large fish proportion
  scenario_lf <- getProportionOfLargeFish(scenario_sim, threshold_w = 500)
  baseline_lf <- getProportionOfLargeFish(baseline_sim, threshold_w = 500)
  
  if (is.matrix(scenario_lf)) {
    scenario_lf_mean <- colMeans(scenario_lf[averaging_period[1]:averaging_period[2], ])
    baseline_lf_mean <- colMeans(baseline_lf[averaging_period[1]:averaging_period[2], , drop = FALSE])
  } else {
    scenario_lf_mean <- mean(scenario_lf[averaging_period[1]:averaging_period[2]])
    baseline_lf_mean <- mean(baseline_lf[averaging_period[1]:averaging_period[2]])
  }
  
  lf_change <- (scenario_lf_mean - baseline_lf_mean) / baseline_lf_mean * 100
  
  # Total system biomass
  scenario_total_bio <- apply(getBiomass(scenario_sim), 1, sum)
  baseline_total_bio <- apply(getBiomass(baseline_sim), 1, sum)
  
  scenario_bio_mean <- mean(scenario_total_bio[averaging_period[1]:averaging_period[2]])
  baseline_bio_mean <- mean(baseline_total_bio[averaging_period[1]:averaging_period[2]])
  
  total_bio_change <- (scenario_bio_mean - baseline_bio_mean) / baseline_bio_mean * 100
  
  return(data.frame(
    metric = c("large_fish_proportion", "total_biomass"),
    change = c(lf_change, total_bio_change),
    stringsAsFactors = FALSE
  ))
}

# Add scenario metadata
add_scenario_metadata <- function(data) {
  data %>%
    mutate(
      target_species = case_when(
        str_detect(scenario, "^atlantic_cod") ~ "atlantic_cod",
        str_detect(scenario, "^capelin") ~ "capelin",
        str_detect(scenario, "^sand_lance") ~ "sand_lance",
        TRUE ~ "unknown"
      ),
      effort_multiplier = as.numeric(str_extract(scenario, "(?<=_F)[0-9.]+"))
    )
}

# =============================================================================
# COMPILE RESULTS
# =============================================================================

cat("\nCalculating comprehensive changes...\n")
comprehensive_results <- map_dfr(sim_results,
                                 ~calculate_comprehensive_changes(.x, sim_baseline),
                                 .id = "scenario")

# Add metadata
comprehensive_results <- add_scenario_metadata(comprehensive_results)

# Add baseline sand lance case (F1.0 with zero change)
baseline_sandlance <- data.frame(
  scenario = "sand_lance_F1.0",
  species = species_names,
  biomass_change = 0,
  feeding_change = 0,
  predmort_change = 0,
  target_species = "sand_lance",
  effort_multiplier = 1.0,
  stringsAsFactors = FALSE
)

comprehensive_results <- rbind(comprehensive_results, baseline_sandlance)

# Calculate system metrics
cat("Calculating system metrics...\n")
system_results <- map_dfr(sim_results, ~calculate_system_metrics(.x, sim_baseline),
                          .id = "scenario")

# Add metadata to system results
system_results <- add_scenario_metadata(system_results)

# Add baseline sand lance case to system results
baseline_sandlance_system <- data.frame(
  scenario = "sand_lance_F1.0",
  metric = c("large_fish_proportion", "total_biomass", "spectrum_slope"),
  change = c(0, 0, 0),
  target_species = "sand_lance",
  effort_multiplier = 1.0,
  stringsAsFactors = FALSE
)

system_results <- rbind(system_results, baseline_sandlance_system)

# Save results
#write.csv(comprehensive_results, "results/comprehensive_ecosystem_results.csv", row.names = FALSE)

# =============================================================================
# VISUALIZATION
# =============================================================================

# Updated species color scheme based on functional groups
create_functional_species_colors <- function() {
  colors <- c(
    # Top predator - Red
    "atlantic_cod" = "#D32F2F",
    
    # Forage fish - Orange tones
    "capelin" = "#F57C00",
    "sand_lance" = "#FFB74D",
    
    # Flatfish - Green tones
    "american_plaice"      = "#1B5E20", # Deep forest green
    "yellowtail_flounder"  = "#388E3C", # Medium green
    "witch_flounder"       = "#66BB6A", # Light green
    "turbot"               = "#A5D6A7", # Pale/minty green
    
    
    # Other predators - Blue/Purple tones
    "redfish" = "#5E35B1",
    "thorny_skate" = "chocolate4"
  )
  return(colors)
}


# Updated target vs non-target interaction plot
plot_target_vs_nontarget_updated <- function(biomass_data, target_filter = NULL, max_effort = 3.0) {
  
  if (!is.null(target_filter)) {
    biomass_data <- biomass_data %>% filter(target_species %in% target_filter)
  }
  
  # Get target species changes
  target_changes <- biomass_data %>%
    filter(species == target_species, effort_multiplier <= max_effort) %>%
    select(scenario, target_species, target_change = biomass_change, effort_multiplier)
  
  # Get non-target species changes and join with target changes
  nontarget_changes <- biomass_data %>%
    filter(species != target_species, effort_multiplier <= max_effort) %>%
    left_join(target_changes, by = c("scenario", "target_species", "effort_multiplier")) %>%
    mutate(
      target_label = case_when(
        target_species == "atlantic_cod" ~ "Atlantic Cod",
        target_species == "capelin" ~ "Capelin",
        target_species == "sand_lance" ~ "Sand Lance",
        TRUE ~ str_to_title(str_replace_all(target_species, "_", " "))
      ),
      species_label = str_to_title(str_replace_all(species, "_", " ")),
      fishing_effort_label = case_when(
        effort_multiplier == 0.3 ~ "0.3x",
        effort_multiplier == 0.5 ~ "0.5x",
        effort_multiplier == 0.7 ~ "0.7x",
        effort_multiplier == 1.0 ~ "1.0x",
        effort_multiplier == 1.5 ~ "1.5x",
        effort_multiplier == 2.0 ~ "2.0x",
        effort_multiplier == 3.0 ~ "3.0x",
        TRUE ~ paste0(effort_multiplier, "x")
      )
    )
  
  # Get functional colors
  functional_colors <- create_functional_species_colors()
  
  # Create text labels data frame with proper y-positioning
  text_labels <- target_changes %>%
    mutate(
      target_label = case_when(
        target_species == "atlantic_cod" ~ "Atlantic Cod",
        target_species == "capelin" ~ "Capelin",
        target_species == "sand_lance" ~ "Sand Lance",
        TRUE ~ str_to_title(str_replace_all(target_species, "_", " "))
      ),
      effort_label = paste0(effort_multiplier, "x")
    )
  
  # Get y-axis range for each panel to position text
  y_range <- nontarget_changes %>%
    group_by(target_label) %>%
    summarise(
      y_max = max(biomass_change, na.rm = TRUE),
      y_min = min(biomass_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_text = y_max + (y_max - y_min) * 0.05) # Position text 5% above max
  
  # Join text positioning with text labels
  text_labels <- text_labels %>%
    left_join(y_range, by = "target_label")
  
  # Create the plot
  p <- ggplot(nontarget_changes, aes(x = target_change, y = biomass_change)) +
    # Add vertical lines for each fishing scenario
    #geom_vline(data = target_changes, aes(xintercept = target_change),
    #           color = "gray60", linetype = "solid", alpha = 0.7) +
    
    # Add species response lines (increased line width)
    geom_line(aes(group = interaction(species_label, target_label), color = species),
              size = 1.2, alpha = 0.8) +
    
    # Add open squares for data points
    geom_point(shape = 22, size = 3, stroke = 0.8, fill = "white", color = "black") +
    
    # Add text labels for fishing scenarios with proper positioning
    geom_text(data = text_labels,
              aes(x = target_change, y = y_text, label = effort_label),
              size = 3, color = "gray40", vjust = 0) +
    
    facet_wrap(~ target_label, scales = "free") +
    scale_x_continuous(limits = c(-20, 10), breaks = seq(-20, 10, by = 5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    
    # Updated color scheme
    scale_color_manual(values = functional_colors, name = "Species",
                       labels = function(x) str_to_title(str_replace_all(x, "_", " "))) +
    
    labs(
      x = expression(Delta * " Target Species Biomass (%)"),
      y = expression(Delta * " Non-target Species Biomass (%)")
    ) +
    
    # Reference lines at zero
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    
    theme_bw() +
    theme(
      # Remove minor grid lines and light gray lines
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      
      strip.text = element_text(size = 16, face = "bold"),
      legend.position = "right",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13)
    )
  
  return(p)
}

# Generate the updated plot
species_colors_functional <- create_functional_species_colors()
fig_interactions_updated <- plot_target_vs_nontarget_updated(comprehensive_results, max_effort = 3.0)
print(fig_interactions_updated)

ggsave("figures/target_vs_nontarget_interactions_FINAL_050925.png",
       plot = fig_interactions_updated,
       device = "png",
       width = 16,
       height = 6,
       units = "in",
       dpi = 600)

# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

cat("\nAnalysis Complete!\n")
cat("Results saved to: results/comprehensive_ecosystem_results.csv\n")
cat("Figure saved to: figures/target_vs_nontarget_interactions.png\n")

# Display summary statistics
cat("\nSummary Statistics:\n")
summary_stats <- comprehensive_results %>%
  group_by(target_species) %>%
  summarise(
    scenarios = n_distinct(scenario),
    species_analyzed = n_distinct(species),
    avg_biomass_change = mean(abs(biomass_change), na.rm = TRUE),
    max_biomass_change = max(abs(biomass_change), na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_stats)

# =============================================================================
# POPULATION-LEVEL CHANGES
# =============================================================================

# Modified Function: Comprehensive ecosystem response plot with biomass change labels
# Now includes sand lance scenarios
plot_ecosystem_responses_simple <- function(data, target_filter = NULL, max_effort = 3.0,
                                            show_title = TRUE, show_legend = TRUE) {
  
  if (!is.null(target_filter)) {
    data <- data %>% filter(target_species %in% target_filter)
  }
  
  plot_data <- data %>%
    filter(effort_multiplier <= max_effort) %>%
    filter(species != target_species) %>%  # Remove target species
    select(scenario, target_species, species, effort_multiplier,
           biomass_change, predmort_change) %>%
    pivot_longer(cols = c(biomass_change, predmort_change),
                 names_to = "metric", values_to = "change") %>%
    mutate(
      target_label = case_when(
        target_species == "atlantic_cod" ~ "Atlantic cod",
        target_species == "capelin" ~ "Capelin",
        target_species == "sand_lance" ~ "Sand lance",
        TRUE ~ str_to_title(str_replace_all(target_species, "_", " "))
      ),
      species_label = str_to_title(str_replace_all(species, "_", " ")),
      metric_label = case_when(
        metric == "biomass_change" ~ "Biomass",
        metric == "predmort_change" ~ "Predation Mortality",
        TRUE ~ metric
      ),
      # Updated labeling system to handle all three target species
      biomass_change_label = case_when(
        # For cod and capelin (traditional fishing scenarios)
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 0.3 ~ "Recovery- high",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 0.5 ~ "Recovery- medium",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 0.7 ~ "Recovery- low",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 1.0 ~ "Baseline",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 1.5 ~ "Depletion- low",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 2.0 ~ "Depletion- medium",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 3.0 ~ "Depletion- high",
        
        # For sand lance (specific effort scenarios)
        target_species == "sand_lance" & effort_multiplier == 1.0 ~ "Baseline",
        target_species == "sand_lance" & effort_multiplier == 1.5 ~ "Depletion- low",
        target_species == "sand_lance" & effort_multiplier == 2.0 ~ "Depletion- medium",
        target_species == "sand_lance" & effort_multiplier == 3.0 ~ "Depletion- high",
        
        TRUE ~ paste0(effort_multiplier, "x")
      ),
      # Create ordered factor levels for consistent legend
      biomass_change_label = factor(biomass_change_label, levels = c(
        "Recovery- high", "Recovery- medium", "Recovery- low", "Baseline",
        "Depletion- low", "Depletion- medium", "Depletion- high"
      ))
    )
  
  p <- ggplot(plot_data, aes(x = change, y = reorder(species_label, change),
                             shape = biomass_change_label, fill = biomass_change_label)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    facet_grid(target_label ~ metric_label, scales = "free_x") +
    scale_shape_manual(values = rep(21, 10), name = "Biomass\nScenario") +  # Increased to 10 for more scenarios
    scale_fill_manual(
      values = c(
        "Recovery- high" = "#2166AC",     # Dark blue
        "Recovery- medium" = "#5AADDD",   # Medium blue  
        "Recovery- low" = "#92C5DE",      # Light blue
        "Baseline" = "#F7F7F7",          # White/light gray
        "Depletion- low" = "#FDBF6F",    # Light orange
        "Depletion- medium" = "#F46D43",  # Medium red
        "Depletion- high" = "#A50026"     # Dark red
      ),
      name = "Biomass\nScenario",
      drop = FALSE  # Keep all levels even if not present
    ) +
    labs(x = "Relative Change (%)", y = "Species") +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold", size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      plot.title = element_text(size = 18, face = "bold")
    )
  
  if (show_title) {
    p <- p + labs(title = "Ecosystem Responses to Target Species Fishing")
  }
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + theme(legend.position = "bottom")
  }
  
  return(p)
}

# Generate individual plots for each target species
fig1_cod <- plot_ecosystem_responses_simple(
  comprehensive_results,
  target_filter = "atlantic_cod",
  max_effort = 3.0,
  show_title = FALSE,
  show_legend = FALSE
)

fig2_capelin <- plot_ecosystem_responses_simple(
  comprehensive_results,
  target_filter = "capelin",
  max_effort = 3.0,
  show_title = FALSE,
  show_legend = TRUE
)

fig3_sandlance <- plot_ecosystem_responses_simple(
  comprehensive_results,
  target_filter = "sand_lance",
  max_effort = 3.0,
  show_title = FALSE,
  show_legend = FALSE
)

# Create comprehensive summary figure with all three target species
summary_figure_all <- fig1_cod / fig3_sandlance / fig2_capelin

# Alternative: Create a single plot with all target species
fig_all_targets <- plot_ecosystem_responses_simple(
  comprehensive_results,
  target_filter = c("atlantic_cod", "capelin", "sand_lance"),
  max_effort = 3.0,
  show_title = FALSE,
  show_legend = TRUE
)

# Display the plots
print(summary_figure_all)
print(fig_all_targets)

# Save the plots
ggsave("figures/ecosystem_responses_all_targets_stacked.png",
       plot = summary_figure_all,
       device = "png",
       width = 10,
       height = 12,
       units = "in",
       dpi = 600)

ggsave("figures/ecosystem_responses_all_targets_combined_180525.png",
       plot = fig_all_targets,
       device = "png",
       width = 12,
       height = 10,
       units = "in",
       dpi = 600)

# =============================================================================
# COMMUNITY-LEVEL CHANGES
# =============================================================================

# Modified Function: System-level metrics plot - TOTAL BIOMASS ONLY with target species on Y-axis
plot_system_metrics <- function(data, target_filter = NULL, max_effort = 3.0) {
  
  if (!is.null(target_filter)) {
    data <- data %>% filter(target_species %in% target_filter)
  }
  
  plot_data <- data %>%
    filter(effort_multiplier <= max_effort) %>%
    filter(metric == "total_biomass") %>%  # ONLY total biomass
    mutate(
      target_label = case_when(
        target_species == "atlantic_cod" ~ "Atlantic Cod",
        target_species == "capelin" ~ "Capelin",
        target_species == "sand_lance" ~ "Sand Lance",
        TRUE ~ str_to_title(str_replace_all(target_species, "_", " "))
      ),
      # Create biomass change labels for target species
      biomass_change_label = case_when(
        # For cod and capelin (traditional fishing scenarios)
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 0.3 ~ "Recovery- high",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 0.5 ~ "Recovery- medium",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 0.7 ~ "Recovery- low",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 1.0 ~ "Baseline",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 1.5 ~ "Depletion- low",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 2.0 ~ "Depletion- medium",
        target_species %in% c("atlantic_cod", "capelin") & effort_multiplier == 3.0 ~ "Depletion- high",
        
        # For sand lance (specific effort scenarios)
        target_species == "sand_lance" & effort_multiplier == 1.0 ~ "Baseline",
        target_species == "sand_lance" & effort_multiplier == 1.5 ~ "Depletion- low",
        target_species == "sand_lance" & effort_multiplier == 2.0 ~ "Depletion- medium",
        target_species == "sand_lance" & effort_multiplier == 3.0 ~ "Depletion- high",
        
        TRUE ~ paste0(effort_multiplier, "x")
      ),
      biomass_change_label = factor(biomass_change_label, levels = c(
        "Recovery- high", "Recovery- medium", "Recovery- low", "Baseline",
        "Depletion- low", "Depletion- medium", "Depletion- high"
      ))
    )
  
  ggplot(plot_data, aes(x = change, y = target_label,
                        shape = biomass_change_label, fill = biomass_change_label)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
    scale_shape_manual(values = rep(21, 7), name = "Fishing\nScenario") +
    scale_fill_manual(
      values = c(
        "Recovery- high" = "#2166AC",     # Dark blue
        "Recovery- medium" = "#5AADDD",   # Medium blue  
        "Recovery- low" = "#92C5DE",      # Light blue
        "Baseline" = "#F7F7F7",          # White/light gray
        "Depletion- low" = "#FDBF6F",    # Light orange
        "Depletion- medium" = "#F46D43",  # Medium red
        "Depletion- high" = "#A50026"     # Dark red
      ),
      name = "Fishing\nScenario",
      drop = FALSE
    ) +
    labs(x = "Relative Change in Community Biomass (%)",
         y = "Target Species") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 16),
      axis.text.y = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      panel.grid.minor = element_blank()
    )
}

# Figure 4: System-level responses (TOTAL BIOMASS ONLY)
fig4_system <- plot_system_metrics(
  system_results,
  max_effort = 3.0
)

# Save community-level figure
ggsave("figures/community_biomass090925.png",
       plot = fig4_system,
       device = "png",
       width = 8,
       height = 5,
       units = "in",
       dpi = 600)

