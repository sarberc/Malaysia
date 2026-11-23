# Figure 3: DOM molecular composition maps of CRAM, NOSC and AImod

# Setup and libraries ----
par(mfrow = c(1, 1))

library(tidyverse)
library(vegan)
library(sf)
library(cowplot)
library(grid)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(terra)

# Load data ----
# Run from the directory containing the input CSV files and forest raster.

data <- read.csv(
  "Metadata_Malaysia_LO.csv",
  sep = ",",
  header = TRUE,
  fill = TRUE
)

ct <- read.csv(
  "DOM_Malaysia_LO.csv",
  sep = ",",
  header = TRUE,
  fill = TRUE
)

ch <- ct[, 1:84]
dom <- ct[, 85:ncol(ct)]

# Load and prepare forest raster ----
forest <- rast("REGIONBorneo_ForestArea_2.tif")

# Original Figure 3 raster extent
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

# Transform extent to forest raster CRS
study_extent_utm <- project(
  study_extent,
  crs(forest)
)

# Crop forest raster
forest_crop <- crop(
  forest,
  study_extent_utm
)

# Combine forest classes
# Classes 1, 2, and 3 are treated as forest.

forest_binary <- forest_crop

forest_binary[
  !forest_binary %in% c(1, 2, 3)
] <- NA

forest_binary[
  !is.na(forest_binary)
] <- 1

# Convert forest raster to data frame
forest_df <- as.data.frame(
  forest_binary,
  xy = TRUE,
  na.rm = TRUE
)

names(forest_df)[3] <- "forest"

# Calculate dom indices ----
# Relative abundance normalisation
domn <- decostand(
  dom,
  method = "total",
  MARGIN = 2
)

podom <- as.data.frame(
  t(domn)
)

# Remove Rajang
keep <- data$River != "Rajang"

podom <- podom[keep, ]

dom <- as.data.frame(
  t(podom)
)

data <- data[keep, ]
# Heteroatom-containing compounds
phet.N <- dom %>%
  filter(ch$N > 0)

phet.S <- dom %>%
  filter(ch$S > 0)

count.N <- phet.N
count.N[count.N > 0] <- 1

count.S <- phet.S
count.S[count.S > 0] <- 1

data$S_ab <- colSums(phet.S)
data$N_ab <- colSums(phet.N)

data$N_num <- colSums(count.N)
data$S_num <- colSums(count.S)

het.dom <- rbind.data.frame(
  phet.N,
  phet.S
)

hetd <- as.data.frame(
  t(het.dom)
)

podat <- data

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

# General DOM info
podat$rich <- specnumber(
  podom
)

podat$diversity <- diversity(
  podom,
  index = "shannon"
)

podat$Het_Perc <- rowSums(
  hetd
) * 100

podat$weighted.OC <- colSums(
  ch$O.C * dom
)

podat$weighted.HC <- colSums(
  ch$H.C * dom
)

podat$highly_unsaturated <- colSums(
  ch$Highly.unsaturated * dom
)

podat$N_containing <- data$N_ab * 100

podat$AImod <- colSums(
  ch$AI.mod * dom
)

# Calculate cram ----
dbe <- 1 + (
  2 * ch$C -
    ch$H +
    ch$N +
    ch$P
) / 2

cram.dom <- as.data.frame(
  t(podom)
) %>%
  filter(
    dbe / ch$C >= 0.3,
    dbe / ch$C <= 0.68,
    dbe / ch$H >= 0.2,
    dbe / ch$H <= 0.95,
    dbe / ch$O >= 0.77,
    dbe / ch$O <= 1.75
  )

podat$CRAM <- colSums(
  cram.dom
)

# Calculate nosc ----
nosc <- -1 * (
  (
    4 * ch$C +
      ch$H -
      3 * ch$N -
      2 * ch$O +
      5 * ch$P -
      2 * ch$S
  ) /
    ch$C
) / 4

nosc.dom <- data.frame(
  nosc * t(podom)
)

podat$NOSC <- colSums(
  nosc.dom
)



# Prepare data for mapping ----
# Remove samples without coordinates
podat <- podat[
  !is.na(podat$Latitude),
]

podat$DOC <- podat$DOC.conc..µM.

lat <- podat$Latitude
long <- podat$Longitude

# Split by sampling campaign
# Pre-monsoon = SW (September); Post-monsoon = NE (March)
podat_sw <- podat[
  podat$Season == "Pre-monsoon",
]

podat_ne <- podat[
  podat$Season == "Post-monsoon",
]

# Average replicates: SW monsoon
pdx_sites_sw <- podat_sw %>%
  dplyr::group_by(
    River,
    Season,
    Station,
    Longitude,
    Latitude
  ) %>%
  dplyr::summarise(
    CRAM = mean(CRAM * 100, na.rm = TRUE),
    NOSC = mean(NOSC, na.rm = TRUE),
    AImod = mean(AImod, na.rm = TRUE),
    .groups = "drop"
  )

# Average replicates: NE monsoon
pdx_sites_ne <- podat_ne %>%
  dplyr::group_by(
    River,
    Season,
    Station,
    Longitude,
    Latitude
  ) %>%
  dplyr::summarise(
    CRAM = mean(CRAM * 100, na.rm = TRUE),
    NOSC = mean(NOSC, na.rm = TRUE),
    AImod = mean(AImod, na.rm = TRUE),
    .groups = "drop"
  )

# Offset Simunjan points for display ----
# Same station-specific offsets as Figure 2.
# Only Simunjan Station 3 during the NE monsoon is shifted.
# Original geographic coordinates remain unchanged.

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

# Prepare map projection and extent ----
sf_use_s2(FALSE)

# Base country polygon
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

# Use native CRS of forest raster
map_crs <- st_crs(
  crs(forest_binary)
)

asia_utm <- st_transform(
  asia_cropped,
  map_crs
)

# Convert plotting coordinates to map CRS
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

pdx_sites_sw <- add_utm_coords(
  pdx_sites_sw
)

pdx_sites_ne <- add_utm_coords(
  pdx_sites_ne
)

# Map extent
#
# Important:
# +3000 m EAST gives Maludam breathing room.
# -3000 m SOUTH gives space below the 1.2°N tick.

forest_ext <- ext(
  forest_binary
)

map_xlim <- c(
  xmin(forest_ext),
  xmax(forest_ext) + 700
)

map_ylim <- c(
  ymin(forest_ext) - 6000,
  ymax(forest_ext)
)

# Figure formatting ----
axis_text_size <- 16
axis_title_size <- 18
title_size <- 20

legend_text_size <- 15
legend_title_size <- 18

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

# Map-plotting function ----
make_dom_plot <- function(
    data,
    fill_var,
    fill_limits = NULL,
    fill_name,
    plot_title = NULL,
    panel_label,
    show_x = FALSE,
    show_y = TRUE,
    add_legend = FALSE,
    labeler = label_number(accuracy = 0.01),
    breaks = waiver()
) {
  
  p <- ggplot() +
    
    geom_sf(
      data = asia_utm,
      fill = "grey90",
      color = NA
    ) +
    
    geom_raster(
      data = forest_df,
      aes(
        x = x,
        y = y
      ),
      fill = "#C7D2C0",
      inherit.aes = FALSE
    ) +
    
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
    
    scale_fill_continuous(
      type = "viridis",
      limits = fill_limits,
      name = fill_name,
      labels = labeler,
      breaks = breaks,
      oob = scales::squish
    ) +
    
    coord_sf(
      crs = map_crs,
      datum = st_crs(4326),
      xlim = map_xlim,
      ylim = map_ylim,
      expand = FALSE,
      ndiscr = 50
    ) +
    
    big_theme +
    
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
    
    labs(
      x = "Longitude (°E)",
      y = "Latitude (°N)",
      title = if (is.null(plot_title)) " " else plot_title
    )
  
  
  # Keep x-axis space but hide labels on rows 1 and 2
  if (!show_x) {
    
    p <- p +
      theme(
        axis.title.x =
          element_text(color = "transparent"),
        
        axis.text.x =
          element_text(color = "transparent"),
        
        axis.ticks.x =
          element_line(color = "transparent")
      )
  }
  
  
  # Keep y-axis title space on right-hand panels
  if (!show_y) {
    
    p <- p +
      theme(
        axis.title.y =
          element_text(color = "transparent")
      )
  }
  
  
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

# Create six map panels ----
# A: CRAM, SW monsoon
p1 <- make_dom_plot(
  pdx_sites_sw,
  "CRAM",
  c(53, 70),
  "CRAM",
  plot_title = "SW monsoon",
  panel_label = "A",
  show_x = FALSE,
  show_y = TRUE,
  add_legend = FALSE,
  labeler = label_number(
    accuracy = 1,
    suffix = "%"
  ),
  breaks = c(
    54, 58, 62, 66, 70
  )
)

# B: CRAM, NE monsoon
p2 <- make_dom_plot(
  pdx_sites_ne,
  "CRAM",
  c(53, 70),
  "CRAM",
  plot_title = "NE monsoon",
  panel_label = "B",
  show_x = FALSE,
  show_y = FALSE,
  add_legend = FALSE,
  labeler = label_number(
    accuracy = 1,
    suffix = "%"
  ),
  breaks = c(
    54, 58, 62, 66, 70
  )
)

# C: NOSC, SW monsoon
p3 <- make_dom_plot(
  pdx_sites_sw,
  "NOSC",
  c(-0.97, -0.90),
  "NOSC",
  plot_title = NULL,
  panel_label = "C",
  show_x = FALSE,
  show_y = TRUE,
  add_legend = FALSE,
  labeler = label_number(
    accuracy = 0.01
  ),
  breaks = c(
    -0.96,
    -0.94,
    -0.92,
    -0.90
  )
)

# D: NOSC, NE monsoon
p4 <- make_dom_plot(
  pdx_sites_ne,
  "NOSC",
  c(-0.97, -0.90),
  "NOSC",
  plot_title = NULL,
  panel_label = "D",
  show_x = FALSE,
  show_y = FALSE,
  add_legend = FALSE,
  labeler = label_number(
    accuracy = 0.01
  ),
  breaks = c(
    -0.96,
    -0.94,
    -0.92,
    -0.90
  )
)

# E: AImod, SW monsoon
p5 <- make_dom_plot(
  pdx_sites_sw,
  "AImod",
  c(0.34, 0.4214),
  bquote(bold(AI[mod])),
  plot_title = NULL,
  panel_label = "E",
  show_x = TRUE,
  show_y = TRUE,
  add_legend = FALSE,
  labeler = label_number(
    accuracy = 0.01
  ),
  breaks = c(
    0.34,
    0.37,
    0.40,
    0.42
  )
)

# F: AImod, NE monsoon
p6 <- make_dom_plot(
  pdx_sites_ne,
  "AImod",
  c(0.34, 0.4214),
  bquote(bold(AI[mod])),
  plot_title = NULL,
  panel_label = "F",
  show_x = TRUE,
  show_y = FALSE,
  add_legend = FALSE,
  labeler = label_number(
    accuracy = 0.01
  ),
  breaks = c(
    0.34,
    0.37,
    0.40,
    0.42
  )
)

# Create separate legends ----
cram_legend_plot <- make_dom_plot(
  pdx_sites_ne,
  "CRAM",
  c(53, 70),
  "CRAM",
  panel_label = "",
  add_legend = TRUE,
  labeler = label_number(
    accuracy = 1,
    suffix = "%"
  ),
  breaks = c(
    54, 58, 62, 66, 70
  )
)

nosc_legend_plot <- make_dom_plot(
  pdx_sites_ne,
  "NOSC",
  c(-0.97, -0.90),
  "NOSC",
  panel_label = "",
  add_legend = TRUE,
  labeler = label_number(
    accuracy = 0.01
  ),
  breaks = c(
    -0.96,
    -0.94,
    -0.92,
    -0.90
  )
)

aimod_legend_plot <- make_dom_plot(
  pdx_sites_ne,
  "AImod",
  c(0.34, 0.4214),
  bquote(bold(AI[mod])),
  panel_label = "",
  add_legend = TRUE,
  labeler = label_number(
    accuracy = 0.01
  ),
  breaks = c(
    0.34,
    0.37,
    0.40,
    0.42
  )
)

# Extract legends
cram_legend <- cowplot::get_legend(
  cram_legend_plot
)

nosc_legend <- cowplot::get_legend(
  nosc_legend_plot
)

aimod_legend <- cowplot::get_legend(
  aimod_legend_plot
)

# Format longitude labels ----
# Removes °E from bottom tick labels because the units are already shown
# in the x-axis title.
#
# Also hides every second longitude tick label to prevent crowding.

format_longitude_labels <- function(p) {
  
  g <- ggplotGrob(p)
  
  axis_index <- grep(
    "axis-b",
    g$layout$name
  )
  
  
  if (length(axis_index) > 0) {
    
    axis_grob <- g$grobs[[axis_index[1]]]
    
    
    format_text <- function(x) {
      
      if (inherits(x, "text")) {
        
        labs <- x$label
        
        # coord_sf stores longitude labels as plotmath expressions
        is_lon <- grepl(
          "degree.*E",
          as.character(labs)
        )
        
        
        if (any(is_lon)) {
          
          lon_positions <- which(
            is_lon
          )
          
          
          # Extract only the numeric part
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
          
          
          # Replace expression labels with plain numeric labels
          labs[lon_positions] <- new_labels
          
          
          # Hide every second longitude label
          if (length(lon_positions) > 3) {
            
            remove_positions <- lon_positions[
              seq(
                2,
                length(lon_positions),
                by = 2
              )
            ]
            
            labs[remove_positions] <- ""
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

# Apply formatting only to bottom row

p5_fixed <- format_longitude_labels(
  p5
)

p6_fixed <- format_longitude_labels(
  p6
)

# Arrange figure rows ----
# Layout:
#
# map | map | legend
#
# Keeping the legends in a third column ensures both map panels
# have identical widths.

# Row 1: CRAM
row1 <- plot_grid(
  p1,
  p2,
  cram_legend,
  nrow = 1,
  rel_widths = c(
    1,
    1,
    0.20
  ),
  align = "h",
  axis = "tb"
)

# Row 2: NOSC
row2 <- plot_grid(
  p3,
  p4,
  nosc_legend,
  nrow = 1,
  rel_widths = c(
    1,
    1,
    0.20
  ),
  align = "h",
  axis = "tb"
)

# Row 3: AImod
row3 <- plot_grid(
  p5_fixed,
  p6_fixed,
  aimod_legend,
  nrow = 1,
  rel_widths = c(
    1,
    1,
    0.20
  ),
  align = "h",
  axis = "tb"
)

# Assemble final figure ----
panel_grid <- plot_grid(
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(
    1,
    1,
    1
  )
)

# Export figure 3 ----
pdf(
  "Fig_3.pdf",
  width = 11,
  height = 17.5
)

print(
  panel_grid
)

dev.off()

# Display in RStudio
panel_grid

