# ==============================================================================
# Fishing Sensitivity Analysis for MIZER
# Author: Raquel
# Description: Analysis of Beverton-Holt recruitment and yield curves
# ==============================================================================

library(mizer)
library(ggplot2)
library(patchwork)
library(dplyr)

#' Plot Beverton-Holt stock-recruitment relationship
plotBevertonHolt <- function(params, species) {
  # Extract parameters for selected species
  select <- species_params(params)$species == species
  erepro <- species_params(params)$erepro[select]
  w0 <- params@w[params@w_min_idx[select]]
  
  # Calculate steady state values
  E_R_ss <- getRDI(params)[select] / erepro * 2 * w0
  R_dd_ss <- getRDD(params)[select]
  R_max <- species_params(params)$R_max[select]
  
  # Create sequence for plotting
  E_R <- seq(0, 2 * E_R_ss, length.out = 50)
  R_di <- erepro * E_R / 2 / w0  # Density independent recruitment
  R_dd <- R_di / (1 + R_di / R_max)  # Density dependent recruitment
  
  # Prepare data for plotting
  df <- reshape2::melt(data.frame(E_R, R_dd, R_di, R_max = rep(R_max, length(E_R))), 
                       id.vars = "E_R")
  
  # Create plot
  ggplot(df) +
    geom_line(aes(x = E_R, y = value, linetype = variable)) +
    geom_point(aes(x = E_R_ss, y = R_dd_ss), size = 3, color = "red") +
    ylim(NA, 1.1 * R_max) +
    ylab("Reproduction rate [eggs/year]") +
    xlab("Energy invested [g/year]") +
    ggtitle(species) +
    theme_bw()
}

#' Create Beverton-Holt plots for all species
create_BH_plots <- function(params) {
  # Get list of species
  species_list <- species_params(params)$species
  
  # Create individual plots
  plots <- list()
  for (sp in species_list) {
    plots[[sp]] <- plotBevertonHolt(params, sp)
  }
  
  # Combine plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = 4) & 
    theme_bw() & 
    theme(plot.title = element_text(size = 10))
  
  return(combined_plot)
}

#' Create yield vs F plots for all species
create_yield_plots <- function(params, 
                               F_range_default = seq(0.001, 1, 0.01),
                               F_range_special = seq(0.001, 0.5, 0.001)) {
  
  # Get list of species
  species_list <- species_params(params)$species
  
  # Species that need lower F range
  #special_F_species <- c("atlantic_cod", "turbot", "thorny_skate")
  special_F_species <- c("thorny_skate")
  
  # Create individual plots
  plots <- list()
  for (sp in species_list) {
    # Select appropriate F range based on species
    F_range <- if (sp %in% special_F_species) F_range_special else F_range_default
    
    # Create plot with error handling
    tryCatch({
      plots[[sp]] <- plotYieldVsF(params, species = sp, F_range = F_range) +
        ggtitle(sp) +
        theme_bw()
    }, error = function(e) {
      message(paste("Error plotting yield curve for", sp, ":", e$message))
    })
  }
  
  # Remove NULL elements (from errors)
  plots <- plots[!sapply(plots, is.null)]
  
  # Combine plots using patchwork
  combined_plot <- wrap_plots(plots, ncol = 4) &
    theme_bw() &
    theme(plot.title = element_text(size = 10))
  
  return(combined_plot)
}

#' Run complete fishing sensitivity analysis
run_fishing_sensitivity <- function(params, output_dir = "Figures/") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Display current reproduction ratio
  cat("Reproduction ratio (RDI/RDD):\n")
  print(getRDI(params)/getRDD(params))
  
  cat("Reproduction level:\n")
  print(getReproductionLevel(params))
  
  cat("erepro and R_max values:\n")
  print(species_params(params) %>% select(species, erepro, R_max))
  
  # Generate Beverton-Holt plots
  BH_plot <- create_BH_plots(params)
  
  # Save Beverton-Holt plot
  ggsave(
    file.path(output_dir, "BevertonHolt_species_panel.png"),
    BH_plot,
    width = 40, 
    height = 30, 
    units = "cm"
  )
  
  # Generate yield vs F plots
  yield_plot <- create_yield_plots(params)
  
  # Save yield curve plot
  ggsave(
    file.path(output_dir, "YieldCurve_species_panel.png"),
    yield_plot,
    width = 40, 
    height = 30, 
    units = "cm"
  )
  
  return(list(
    BH_plot = BH_plot,
    yield_plot = yield_plot
  ))
}

# Example usage:
params1<- readRDS("rds/params_adjusted_growth_FINAL.rds")
getReproductionLevel(params1)

results <- run_fishing_sensitivity(params1, output_dir = "figures/FINAL/Fishing_sensitivity")
results$BH_plot
results$yield_plot
