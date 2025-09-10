#======== A. LOAD LIBRARIES =========================================================
# A. LOAD LIBRARIES ================================================================
rm(list = ls())

# Load libraries, install if needed
library(ggplot2)
library(devtools)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tidyverse)
library(viridis)
library(patchwork)

#lets' change the plotting colours
library(viridis)
params_uncalibrated@linecolour[1:11] <-plasma(11)
params_uncalibrated@linecolour["Resource"] <-"seagreen3"

bio_obs_file = "data/biomass_historical.csv"
effort <- read.csv("data/Fmort_ratio_historical_2025.csv", sep = ";", row.names = 1)

#======== B. PLOT F AND (RELATIVE) ABUNDANCE ========================================
# Read stock assessment and survey index data. See Appendix S1.
Biomass <- read.csv("data/biomass_historical.csv", sep = ";")
head(Biomass)
Biomass_long<- gather(Biomass, species, biomass, american_plaice:thorny_skate)
Biomass_long<- rename(Biomass_long, year = survey.year)

#Fmort <- read.csv("data/Catches_historical.csv", sep = ";")
Fmort <- read.csv("data/Fmort_ratio_historical_2025.csv", sep=";")
head(Fmort)
Fmort_long<- gather(Fmort, species, Fmort, american_plaice:thorny_skate)
Fmort_long<- rename(Fmort_long, year = survey.year)
#Fmort_long<- rename(Fmort_long, year = Year)

dat<- left_join(Biomass_long, Fmort_long)

# Create dataframes with mean Biomass and F for calibration period plotting
min_cal_yr <- 2000
max_cal_yr <- 2010

ref_time <- data.frame(Year = c(min_cal_yr, max_cal_yr),
                       B  = c(0, max(dat$biomass, na.rm = T)),
                       Fm    = c(0, max(dat$Fmort, na.rm = T)))

mean_B_F <- dat %>%
  filter(year >= min_cal_yr & year <= max_cal_yr) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(mean_B = mean(biomass),
                mean_F   = mean(Fmort))

# Define palette
create_species_colors <- function() {
  colors <- c(
    "american_plaice" = "#8C5400",
    "atlantic_cod" = "#5F53E6",
    "redfish" = "#B5BC30",
    "turbot" = "#EB62EB",
    "witch_flounder" = "#002F00",
    "yellowtail_flounder" = "#3F006E",
    "capelin" = "#84B48C",
    "sand_lance" = "#EF2B62",
    "thorny_skate" = "#007E7A"
  )
  return(colors)
}

col <- create_species_colors()

# Plot Biomass and F, overlay mean for each period
p1 <- ggplot(dat, aes(year, biomass, color = species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE,
            aes(xmin = min(Year),
                xmax = max(Year),
                ymin = 0,
                ymax = max(B)),
            fill  = "black", alpha = 0.1) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(data=mean_B_F, aes(x=year, y=mean_B), size = 1.7, shape = 15) +
  theme(aspect.ratio = 1) +
  ylab("Biomass (kg)") +
  scale_color_manual(values = col) +
  #scale_color_viridis(discrete = TRUE, option = "cividis") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(1995, 2020, 5),
                     limits = c(1995, 2020)) +
  guides(color = FALSE) +
  theme_classic() + 
  theme(
    #legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )


p2 <- ggplot(dat, aes(year, Fmort, color = species)) +
  geom_rect(data = ref_time, inherit.aes = FALSE,
            aes(xmin = min(Year),
                xmax = max(Year),
                ymin = 0,
                ymax = max(Fm)),
            fill  = "black", alpha = 0.1) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(data=mean_B_F, aes(x=year, y=mean_F), size = 1.7, shape = 15) +
  theme(aspect.ratio = 1) +
  ylab("F") +
  scale_color_manual(values = col) +
  #scale_color_viridis(discrete = TRUE, option = "cividis") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(1995, 2020, 5),
                     limits = c(1995, 2020)) +
  theme_classic() + 
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )


pWord<- p1 + p2 + plot_annotation(tag_levels = "A")

ggsave("figures/Biomass_Fmort_timeSeries_180725.png",
       plot = pWord,
       device = "png",
       width = 12,
       height = 6,
       units = "in",
       dpi = 600)


# Get mean SSB and F for Appendix
mean_B_F %>%
  group_by(species) %>%
  summarize(mean(mean_B),
            mean(mean_F))
