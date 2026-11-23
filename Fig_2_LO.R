# Figure 2: DOC and greenhouse gas maps
# Setup and libraries ----
par(mfrow = c(1, 1))

library(tidyverse)
library(sf)
library(cowplot)
library(grid)
library(rnaturalearth)
library(scales)
library(terra)

# Load data ----
# Run this script from the folder containing the input data and raster.
data <- read.csv(
  "Metadata_Malaysia_LO.csv",
  sep = ",",
  header = TRUE,
  fill = TRUE
)

# Keep study rivers

sim <- data %>%
  filter(River == "Simunjan")

seb <- data %>%
  filter(River == "Sebuyau")

mal <- data %>%
  filter(River == "Maludam")

podat <- rbind.data.frame(
  mal,
  seb,
  sim
)

# Keep freshwater samples

podat <- podat %>%
  filter(Sal <= 0.5)

# Extract station number

podat <- podat %>%
  mutate(
    Station = as.numeric(
      stringr::str_extract(
        FT_Name,
        "(?<=st)\\d+"
      )
    )
  )

# Split by season

podat_sw <- podat[
  podat$Season == "Pre-monsoon",
]

podat_ne <- podat[
  podat$Season == "Post-monsoon",
]

# Prepare forest raster ----
# Same raster and crop as Figure 3.

forest <- rast(
  "REGIONBorneo_ForestArea_2.tif"
)

# Original study extent
# Keep this unchanged.

study_extent <- vect(
  matrix(
    c(
      110.73, 1.20,
      111.10, 1.20,
      111.10, 1.70,
      110.73, 1.70,
      110.73, 1.20
    ),
    ncol = 2,
    byrow = TRUE
  ),
  type = "polygons",
  crs = "EPSG:4326"
)

# Transform extent to raster CRS

study_extent_utm <- project(
  study_extent,
  crs(forest)
)

# Crop raster

forest_crop <- crop(
  forest,
  study_extent_utm
)

# Retain forest classes 1, 2 and 3

forest_binary <- forest_crop

forest_binary[
  !forest_binary %in% c(1, 2, 3)
] <- NA

forest_binary[
  !is.na(forest_binary)
] <- 1

# Convert raster to data frame

forest_df <- as.data.frame(
  forest_binary,
  xy = TRUE,
  na.rm = TRUE
)

names(forest_df)[3] <- "forest"

# Calculate seasonal site means ----
prepare_sites <- function(dat) {
  
  pdx <- dat[
    dat$Sal <= 0.5,
  ]
  
  
  # Create less lengthy DOC variable
  if (
    !"DOC" %in% names(pdx) &&
    "DOC.conc..µM." %in% names(pdx)
  ) {
    
    pdx$DOC <- pdx$DOC.conc..µM.
  }
  
  
  # Average replicates by station
  pdx_sites <- pdx %>%
    
    group_by(
      River,
      Season,
      Station
    ) %>%
    
    summarise(
      
      Longitude = mean(
        Longitude,
        na.rm = TRUE
      ),
      
      Latitude = mean(
        Latitude,
        na.rm = TRUE
      ),
      
      DOC = mean(
        DOC,
        na.rm = TRUE
      ),
      
      across(
        where(is.numeric),
        ~ mean(
          .x,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
    )
  
  
  # Convenient CH4 variable
  pdx_sites$CH4 <-
    pdx_sites$CH4..nM.
  
  
  return(
    pdx_sites
  )
}

# Calculate seasonal station means

pdx_sites_sw <- prepare_sites(
  podat_sw
)

pdx_sites_ne <- prepare_sites(
  podat_ne
)

# Offset Simunjan point for display ----
# Only the plotted coordinates of Simunjan Station 3 (NE monsoon) are shifted because points are on top of eachother
# Original coordinates are retained for reference.
add_plot_coords <- function(dat) {
  
  dat %>%
    mutate(
      lon_offset = case_when(
        River == "Simunjan" &
          Season == "Post-monsoon" &
          Station == 3 ~ -0.003,
        TRUE ~ 0
      ),
      
      lat_offset = case_when(
        River == "Simunjan" &
          Season == "Post-monsoon" &
          Station == 3 ~ 0.003,
        TRUE ~ 0
      ),
      
      Longitude_plot = Longitude + lon_offset,
      Latitude_plot = Latitude + lat_offset
    )
}

# Apply plotting offsets

pdx_sites_sw <- add_plot_coords(pdx_sites_sw)
pdx_sites_ne <- add_plot_coords(pdx_sites_ne)

# Prepare variable-specific data ----
# CH4

ch4_sw <- pdx_sites_sw[
  !is.na(pdx_sites_sw$CH4),
]

ch4_ne <- pdx_sites_ne[
  !is.na(pdx_sites_ne$CH4),
]

# N2O

n2o_sw <- pdx_sites_sw[
  !is.na(pdx_sites_sw$N2O..nM.),
]

n2o_ne <- pdx_sites_ne[
  !is.na(pdx_sites_ne$N2O..nM.),
]

# CO2

co2_sw <- pdx_sites_sw %>%
  filter(
    !is.na(CO2..ppm.),
    !is.nan(CO2..ppm.)
  )

# Prepare map projection ----
sf_use_s2(FALSE)

# Geographic bounds for country polygon

lat <- podat$Latitude
long <- podat$Longitude

world <- ne_countries(
  scale = "large",
  continent = "asia",
  returnclass = "sf"
)

asia_cropped <- st_crop(
  world,
  xmin = min(long, na.rm = TRUE) - 0.1,
  xmax = max(long, na.rm = TRUE) + 0.1,
  ymin = min(lat, na.rm = TRUE) - 0.1,
  ymax = max(lat, na.rm = TRUE) + 0.1
)

# Same projected CRS as Figure 3

map_crs <- st_crs(
  crs(forest_binary)
)

asia_utm <- st_transform(
  asia_cropped,
  map_crs
)

# Convert station plotting coordinates to map CRS

add_utm_coords <- function(dat) {
  
  dat_sf <- st_as_sf(
    dat,
    coords = c(
      "Longitude_plot",
      "Latitude_plot"
    ),
    crs = 4326,
    remove = FALSE
  )
  
  
  dat_utm <- st_transform(
    dat_sf,
    map_crs
  )
  
  
  xy <- st_coordinates(
    dat_utm
  )
  
  
  dat$X_utm <- xy[, 1]
  dat$Y_utm <- xy[, 2]
  
  
  dat
}

# Transform seasonal site data

pdx_sites_sw <- add_utm_coords(
  pdx_sites_sw
)

pdx_sites_ne <- add_utm_coords(
  pdx_sites_ne
)

ch4_sw <- add_utm_coords(
  ch4_sw
)

ch4_ne <- add_utm_coords(
  ch4_ne
)

n2o_sw <- add_utm_coords(
  n2o_sw
)

n2o_ne <- add_utm_coords(
  n2o_ne
)

co2_sw <- add_utm_coords(
  co2_sw
)

# SAME MAP EXTENT AS FINAL FIGURE 3

forest_ext <- ext(
  forest_binary
)

# Final Fig. 3 setting:
# +700 m to EAST

map_xlim <- c(
  xmin(forest_ext),
  xmax(forest_ext) + 700
)

# Final Fig. 3 setting:
# -3000 m to SOUTH

map_ylim <- c(
  ymin(forest_ext) - 6000,
  ymax(forest_ext)
)

# Figure formatting ----
axis_text_size <- 16
axis_title_size <- 18
title_size <- 20
legend_text_size <- 13
legend_title_size <- 15

point_size <- 9
label_size <- 8

# Main map theme

big_theme <- theme(
  
  axis.text = element_text(
    size = axis_text_size,
    color = "black"
  ),
  
  axis.title = element_text(
    size = axis_title_size,
    face = "bold"
  ),
  
  plot.title = element_text(
    size = title_size,
    face = "bold",
    hjust = 0.5
  ),
  
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  
  panel.background = element_rect(
    fill = "white",
    color = NA
  ),
  
  panel.border = element_rect(
    color = "black",
    fill = NA,
    linewidth = 1.2
  ),
  
  text = element_text(
    size = axis_text_size
  ),
  
  plot.margin = margin(
    2, 2, 2, 2
  )
)

# Legend theme

legend_theme <- theme(
  
  legend.position = "right",
  
  legend.key.height = unit(
    1.4,
    "cm"
  ),
  
  legend.key.width = unit(
    0.45,
    "cm"
  ),
  
  legend.text = element_text(
    size = legend_text_size
  ),
  
  legend.title = element_text(
    size = legend_title_size,
    face = "bold",
    margin = margin(b = 7)
  )
)

# Map plotting function ----
make_map_plot <- function(
    data,
    fill_var,
    fill_limits = NULL,
    fill_name,
    plot_title = NULL,
    panel_label,
    show_x = FALSE,
    show_y = TRUE,
    add_legend = FALSE,
    labeler = label_number(
      accuracy = 0.01
    ),
    breaks = waiver()
) {
  
  p <- ggplot() +
    
    # Land
    geom_sf(
      data = asia_utm,
      fill = "grey90",
      color = NA
    ) +
    
    # Forest
    geom_raster(
      data = forest_df,
      aes(
        x = x,
        y = y
      ),
      fill = "#C7D2C0",
      inherit.aes = FALSE
    ) +
    
    # Sampling sites
    geom_point(
      data = data,
      aes(
        x = X_utm,
        y = Y_utm,
        fill = .data[[fill_var]]
      ),
      pch = 21,
      color = "black",
      size = point_size
    ) +
    
    # Colour scale
    scale_fill_continuous(
      type = "viridis",
      limits = fill_limits,
      name = fill_name,
      labels = labeler,
      breaks = breaks,
      oob = scales::squish
    ) +
    
    # SAME projected map system as Figure 3
    coord_sf(
      crs = map_crs,
      datum = st_crs(4326),
      xlim = map_xlim,
      ylim = map_ylim,
      expand = FALSE,
      ndiscr = 50
    ) +
    
    big_theme +
    
    # Panel label
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = panel_label,
      hjust = -0.1,
      vjust = 1.3,
      size = label_size,
      fontface = "bold"
    ) +
    
    # Always reserve identical axis/title space
    labs(
      x = "Longitude (°E)",
      y = "Latitude (°N)",
      title =
        if (is.null(plot_title)) {
          " "
        } else {
          plot_title
        }
    )
  
  
  # Hide x axis but preserve its space
  if (!show_x) {
    
    p <- p +
      theme(
        
        axis.title.x =
          element_text(
            color = "transparent"
          ),
        
        axis.text.x =
          element_text(
            color = "transparent"
          ),
        
        axis.ticks.x =
          element_line(
            color = "transparent"
          )
      )
  }
  
  
  # Hide y title but preserve its space
  if (!show_y) {
    
    p <- p +
      theme(
        
        axis.title.y =
          element_text(
            color = "transparent"
          )
      )
  }
  
  
  # Legend
  if (add_legend) {
    
    p <- p +
      legend_theme
    
  } else {
    
    p <- p +
      theme(
        legend.position = "none"
      )
  }
  
  
  p
}

# Export sampling coordinates for ArcGIS ----
export_coords <- function(dat) {
  dat %>%
    select(
      River,
      Season,
      Station,
      Longitude,
      Latitude,
      Longitude_plot,
      Latitude_plot
    ) %>%
    distinct() %>%
    arrange(River, Station)
}

# SW monsoon (September 2017)
coords_SW <- export_coords(pdx_sites_sw)

write.csv(
  coords_SW,
  "Sampling_coordinates_SW_monsoon.csv",
  row.names = FALSE
)

# NE monsoon (March 2017)
coords_NE <- export_coords(pdx_sites_ne)

write.csv(
  coords_NE,
  "Sampling_coordinates_NE_monsoon.csv",
  row.names = FALSE
)

# Create map panels ----
# A: DOC, SW monsoon

p1 <- make_map_plot(
  pdx_sites_sw,
  "DOC",
  c(1000, 4387),
  expression(bold(DOC ~ " (" * mu * "M)")),
  plot_title = "SW monsoon",
  panel_label = "A",
  show_x = FALSE,
  show_y = TRUE,
  add_legend = FALSE,
  breaks = c(1000, 2000, 3000, 4000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

# B: DOC, Monsoon

p2 <- make_map_plot(
  pdx_sites_ne,
  "DOC",
  c(1000, 4387),
  expression(bold(DOC ~ " (" * mu * "M)")),
  plot_title = "NE monsoon",
  panel_label = "B",
  show_x = FALSE,
  show_y = FALSE,
  add_legend = FALSE,
  breaks = c(1000, 2000, 3000, 4000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

# C: CH4, SW monsoon

p3 <- make_map_plot(
  ch4_sw,
  "CH4",
  c(400, 2000),
  expression(bold(CH[4] * " (nM)")),
  plot_title = NULL,
  panel_label = "C",
  show_x = FALSE,
  show_y = TRUE,
  add_legend = FALSE,
  breaks = c(400, 800, 1200, 1600, 2000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

# D: CH4, Monsoon

p4 <- make_map_plot(
  ch4_ne,
  "CH4",
  c(400, 2000),
  expression(bold(CH[4] * " (nM)")),
  plot_title = NULL,
  panel_label = "D",
  show_x = FALSE,
  show_y = FALSE,
  add_legend = FALSE,
  breaks = c(400, 800, 1200, 1600, 2000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

# E: N2O, SW monsoon

p5 <- make_map_plot(
  n2o_sw,
  "N2O..nM.",
  c(2.5, 10),
  expression(bold(N[2] * O * " (nM)")),
  plot_title = NULL,
  panel_label = "E",
  show_x = FALSE,
  show_y = TRUE,
  add_legend = FALSE,
  breaks = c(4, 6, 8, 10),
  labeler = label_number(
    accuracy = 1
  )
)

# F: N2O, Monsoon

p6 <- make_map_plot(
  n2o_ne,
  "N2O..nM.",
  c(2.5, 10),
  expression(bold(N[2] * O * " (nM)")),
  plot_title = NULL,
  panel_label = "F",
  show_x = TRUE,
  show_y = FALSE,
  add_legend = FALSE,
  breaks = c(4, 6, 8, 10),
  labeler = label_number(
    accuracy = 1
  )
)

# G: CO2, SW monsoon

p7 <- make_map_plot(
  co2_sw,
  "CO2..ppm.",
  c(3000, 8000),
  expression(bold(CO[2] * " (ppm)")),
  plot_title = NULL,
  panel_label = "G",
  show_x = TRUE,
  show_y = TRUE,
  add_legend = FALSE,
  breaks = c(3000, 4000, 5000, 6000, 7000, 8000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

# Create legends ----
doc_legend_plot <- make_map_plot(
  pdx_sites_ne,
  "DOC",
  c(1000, 4387),
  expression(bold(DOC ~ " (" * mu * "M)")),
  panel_label = "",
  add_legend = TRUE,
  breaks = c(1000, 2000, 3000, 4000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

ch4_legend_plot <- make_map_plot(
  ch4_ne,
  "CH4",
  c(400, 2000),
  expression(bold(CH[4] * " (nM)")),
  panel_label = "",
  add_legend = TRUE,
  breaks = c(400, 800, 1200, 1600, 2000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

n2o_legend_plot <- make_map_plot(
  n2o_ne,
  "N2O..nM.",
  c(2.5, 10),
  expression(bold(N[2] * O * " (nM)")),
  panel_label = "",
  add_legend = TRUE,
  breaks = c(4, 6, 8, 10),
  labeler = label_number(
    accuracy = 1
  )
)

co2_legend_plot <- make_map_plot(
  co2_sw,
  "CO2..ppm.",
  c(3000, 8000),
  expression(bold(CO[2] * " (ppm)")),
  panel_label = "",
  add_legend = TRUE,
  breaks = c(3000, 4000, 5000, 6000, 7000, 8000),
  labeler = label_number(
    accuracy = 1,
    big.mark = ""
  )
)

doc_legend <- cowplot::get_legend(
  doc_legend_plot
)

ch4_legend <- cowplot::get_legend(
  ch4_legend_plot
)

n2o_legend <- cowplot::get_legend(
  n2o_legend_plot
)

co2_legend <- cowplot::get_legend(
  co2_legend_plot
)

# Format longitude labels ----
# Same approach as final Figure 3:
#
# - remove °E from tick labels
# - keep °E only in axis title
# - hide every second longitude label

format_longitude_labels <- function(p) {
  
  g <- ggplotGrob(p)
  
  axis_index <- grep(
    "axis-b",
    g$layout$name
  )
  
  
  if (length(axis_index) > 0) {
    
    axis_grob <-
      g$grobs[[axis_index[1]]]
    
    
    format_text <- function(x) {
      
      if (inherits(x, "text")) {
        
        labs <- x$label
        
        is_lon <- grepl(
          "degree.*E",
          as.character(labs)
        )
        
        
        if (any(is_lon)) {
          
          lon_positions <- which(
            is_lon
          )
          
          
          new_labels <- vapply(
            labs[lon_positions],
            function(z) {
              
              txt <- paste(
                deparse(z),
                collapse = ""
              )
              
              sub(
                '^"([^"]+)".*$',
                "\\1",
                txt
              )
            },
            character(1)
          )
          
          
          labs[lon_positions] <-
            new_labels
          
          
          if (
            length(lon_positions) > 3
          ) {
            
            remove_positions <-
              lon_positions[
                seq(
                  2,
                  length(lon_positions),
                  by = 2
                )
              ]
            
            labs[remove_positions] <-
              ""
          }
          
          
          x$label <- labs
        }
      }
      
      
      if (!is.null(x$children)) {
        
        x$children <- lapply(
          x$children,
          format_text
        )
      }
      
      
      if (!is.null(x$grobs)) {
        
        x$grobs <- lapply(
          x$grobs,
          format_text
        )
      }
      
      
      x
    }
    
    
    g$grobs[[axis_index[1]]] <-
      format_text(
        axis_grob
      )
  }
  
  
  g
}

# Apply to panels that display longitude

p6_fixed <- format_longitude_labels(
  p6
)

p7_fixed <- format_longitude_labels(
  p7
)

# Arrange figure rows ----
#
# Same architecture as Figure 3:
#
#     map | map | legend

# Row 1: DOC

row1 <- plot_grid(
  p1,
  p2,
  doc_legend,
  nrow = 1,
  rel_widths = c(
    1,
    1,
    0.27
  ),
  align = "h",
  axis = "tb"
)

# Row 2: CH4

row2 <- plot_grid(
  p3,
  p4,
  ch4_legend,
  nrow = 1,
  rel_widths = c(
    1,
    1,
    0.27
  ),
  align = "h",
  axis = "tb"
)

# Row 3: N2O

row3 <- plot_grid(
  p5,
  p6_fixed,
  n2o_legend,
  nrow = 1,
  rel_widths = c(
    1,
    1,
    0.27
  ),
  align = "h",
  axis = "tb"
)

# Row 4: CO2

# Row 4: CO2
# No NE monsoon CO2 data, so place legend directly beside panel G.

row4 <- plot_grid(
  p7_fixed,
  co2_legend,
  NULL,
  nrow = 1,
  rel_widths = c(
    1,
    0.27,
    1
  ),
  align = "h",
  axis = "tb"
)

# Assemble figure ----
panel_grid <- plot_grid(
  row1,
  row2,
  row3,
  row4,
  ncol = 1,
  rel_heights = c(
    1,
    1,
    1,
    1
  )
)

# Export Figure 2 ----
# Same final PDF width as Figure 3: 11 inches.

pdf(
  "Fig2.pdf",
  width = 11,
  height = 22
)

print(
  panel_grid
)

dev.off()

# Display in RStudio

panel_grid
