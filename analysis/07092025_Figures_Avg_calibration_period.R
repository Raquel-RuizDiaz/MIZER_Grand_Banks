# ==============================================================================
# MIZER Grand Banks- 10 year average Validation and Analysis
# Author: Raquel Ruiz Diaz
# ==============================================================================

# Load required packages ----
library(mizer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(scales)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggrepel)

# Set up directories ----
fig_dir <- "figures/FINAL/"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Custom theme for publication-quality plots ----
theme_manuscript <- function() {
  theme_bw() +
    theme(
      text = element_text(family = "Arial", size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "top",
      legend.box = "horizontal",
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
}

theme_publication <- function() {
  theme_bw() +
    theme(
      text = element_text(family = "Arial", size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      legend.box.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
      legend.box.margin = margin(6, 6, 6, 6),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(size = 11, face = "bold")
    )
}

# Load and prepare data ----
params_adjusted <- readRDS("rds/params_adjusted_growth_FINAL.rds")
effort <- params_adjusted@initial_effort

# Run simulation ----
sim_optim <- project(params_adjusted, effort = effort, t_max = 100)

# Display basic simulation plot
plot(sim_optim)

# ==============================================================================
# BIOMASS VALIDATION ANALYSIS
# ==============================================================================

# Prepare validation data ----
create_validation_data <- function(sim_data) {
  # Get observed vs model data
  data_sim <- plotBiomassObservedVsModel(sim_data, return_data = TRUE)
  
  # Calculate validation metrics
  R <- cor(data_sim$observed, data_sim$model)
  pbias <- sum(data_sim$observed - data_sim$model) / sum(data_sim$observed)
  
  # Add relative differences and outlier identification
  data_sim$diff <- (data_sim$observed - data_sim$model) / data_sim$observed
  data_sim$outlier <- abs(data_sim$diff) > 1.5
  data_sim$label <- ifelse(data_sim$outlier, as.character(data_sim$species), "")
  
  # Order species by observed biomass
  species_order <- data_sim %>%
    arrange(desc(observed)) %>%
    pull(species)
  data_sim$species <- factor(data_sim$species, levels = species_order)
  
  return(list(data = data_sim, R = R, pbias = pbias))
}

# Create validation plots ----
create_biomass_barplot <- function(data_sim) {
  # Color palette - colorblind friendly
  model_colors <- c("observed" = "#0072B2", "model" = "#D55E00")
  
  # Transform to long format
  long_data <- data_sim %>%
    pivot_longer(cols = c(model, observed), names_to = "type", values_to = "biomass") %>%
    select(species, biomass, type) %>%
    mutate(species = factor(species, levels = levels(data_sim$species)))
  
  # Create bar plot
  ggplot(long_data, aes(x = species, y = biomass, fill = type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = model_colors,
                      labels = c("Model Predicted", "Observed")) +
    scale_y_continuous(labels = comma_format()) +
    labs(x = "",
         y = expression("Biomass [kg]"),
         fill = "") +
    theme_manuscript() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

create_scatter_plot <- function(data_sim, R) {
  ggplot(data_sim, aes(x = log10(observed), y = log10(model))) +
    geom_abline(slope = 1, intercept = 0, color = "darkgray", linetype = "dashed", size = 1) +
    geom_point(shape = 21, size = 3, fill = "#0072B2", color = "black", alpha = 0.8) +
    labs(x = expression(log[10]~"Observed Biomass [kg]"),
         y = expression(log[10]~"Predicted Biomass [kg]")) +
    annotate("text", 
             x = min(log10(data_sim$observed)) + 0.1*(max(log10(data_sim$observed))-min(log10(data_sim$observed))),
             y = max(log10(data_sim$model)) - 0.1*(max(log10(data_sim$model))-min(log10(data_sim$model))),
             label = paste("R =", round(R, 2)),
             hjust = 0, size = 5, fontface = "bold") +
    geom_text_repel(aes(label = label),
                    box.padding = 0.5,
                    max.overlaps = 10,
                    segment.color = "gray50") +
    theme_manuscript()
}

create_residual_plot <- function(data_sim, pbias) {
  ggplot(data_sim, aes(x = species, y = diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray", size = 1) +
    geom_point(shape = 21, size = 3, fill = "#D55E00", color = "black") +
    geom_segment(aes(xend = species, yend = 0), alpha = 0.5, color = "gray50") +
    scale_y_continuous(limits = c(-2.5, 2.5), breaks = seq(-2.5, 2.5, 0.5)) +
    labs(x = "", y = "Relative Difference (Obs - Pred)/Obs") +
    annotate("text", 
             x = length(unique(data_sim$species)) - 2,
             y = 2.25,
             label = paste("Prop. Bias =", round(pbias, 4)),
             hjust = 0, size = 5, fontface = "bold") +
    theme_manuscript() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.y = element_text(face = "plain"))
}

# Generate validation analysis ----
validation_results <- create_validation_data(sim_optim)
data_sim <- validation_results$data
R <- validation_results$R
pbias <- validation_results$pbias

# Display validation plots
plotBiomassObservedVsModel(sim_optim)

# Create individual plots
bar_plot <- create_biomass_barplot(data_sim)
scatter_plot <- create_scatter_plot(data_sim, R)
residual_plot <- create_residual_plot(data_sim, pbias)

# Combine plots
combined_plot <- (scatter_plot + residual_plot) +
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = 'A')

full_plot <- (bar_plot / (scatter_plot + residual_plot)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Time Average Model Validation",
    tag_levels = 'A'
  ) &
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))

# ==============================================================================
# ECOLOGICAL COMPONENT ANALYSIS
# ==============================================================================

# Function to create and save ecological plots
create_ecological_plots <- function(params, sim_data, save_dir) {
  
  # Create individual plots
  plots <- list(
    energy_budget = plotEnergyBudget(params) + theme_publication(),
    death = plotDeath(params, proportion = FALSE) + theme_publication(),
    diet = plotDiet(params) + theme_publication(),
    pred_mort = plotPredMort(params) + theme_publication(),
    resource_pred = plotResourcePred(params, proportion = TRUE) + theme_publication(),
    feeding_level = plotFeedingLevel(params) + theme_publication(),
    growth_curves = plotGrowthCurves(params, species_panel = TRUE)
  )
  
  # Save plots with appropriate dimensions
  plot_specs <- list(
    energy_budget = list(width = 30, height = 18),
    death = list(width = 30, height = 18),
    diet = list(width = 30, height = 18),
    pred_mort = list(width = 20, height = 12),
    resource_pred = list(width = 20, height = 12),
    feeding_level = list(width = 20, height = 12)
  )
  
  # Save individual plots
  for (plot_name in names(plots)) {
    if (plot_name %in% names(plot_specs)) {
      specs <- plot_specs[[plot_name]]
      ggsave(
        filename = paste0(save_dir, "plot", stringr::str_to_title(plot_name), ".png"),
        plot = plots[[plot_name]],
        width = specs$width,
        height = specs$height,
        units = "cm",
        dpi = 300
      )
    }
  }
  
  # Create combined panel
  diet_death_pred_panel <- plots$diet + plots$death + plots$resource_pred
  
  return(list(plots = plots, combined_panel = diet_death_pred_panel))
}

# Generate ecological plots
ecological_results <- create_ecological_plots(params_adjusted, sim_optim, fig_dir)
ecological_results$plots$pred_mort #other plots in here

# ==============================================================================
# SAVE ALL FIGURES
# ==============================================================================

# Save biomass validation figures
ggsave(paste0(fig_dir, "Biomass_barplot_FINAL.png"), bar_plot,
       width = 34, height = 18, units = "cm", dpi = 300)

ggsave(paste0(fig_dir, "Biomass_calibration_FINAL.png"), combined_plot,
       width = 32, height = 14, units = "cm", dpi = 300)

ggsave(paste0(fig_dir, "AvgBiomass_calibration.png"), full_plot,
       width = 40, height = 30, units = "cm", dpi = 300)

# Save combined ecological panel
ggsave(paste0(fig_dir, "diet_death_pred.png"), ecological_results$combined_panel,
       width = 22, height = 18, units = "cm", dpi = 300)
