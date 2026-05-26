# Figure 2 - Map plots for DOC, CH4, N2O, CO2)
# Inter-Monsoon vs Monsoon


rm(list = ls(all = TRUE))
par(mfrow = c(1, 1))

# ---- Libraries ----
library(ggplot2)
library(dplyr)
library(sf)
library(cowplot)
library(grid)
library(rnaturalearth)
library(rnaturalearthdata)
library(osmdata)
library(scales)

# 1. ---- Load data ----
#setwd()

#data <- read.csv("Metadata_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)
# Convert oxygen consumption rate from µM min^-1 to µM h^-1

sim <- data %>% filter(River == "Simunjan")
seb <- data %>% filter(River == "Sebuyau")
mal <- data %>% filter(River == "Maludam")

podat <- rbind.data.frame(mal, seb, sim)
podat = podat %>% filter(Sal<=0.5)


podat <- podat %>%
  mutate(
    Station = as.numeric(
      stringr::str_extract(FT_Name, "(?<=st)\\d+")
    )
  )

# ---- Split by season ----
podat1 <- podat[podat$Season == "Pre-monsoon", ]   # Inter-Monsoon / September
podat2 <- podat[podat$Season == "Post-monsoon", ]  # Monsoon / March

# ---- Map boundary setup ----
lat <- podat$Latitude
long <- podat$Longitude

world <- ne_countries(scale = "large", continent = "asia", returnclass = "sf")

asia_cropped <- st_crop(
  world,
  xmin = min(long, na.rm = TRUE) - 0.1,
  xmax = max(long, na.rm = TRUE) + 0.1,
  ymin = min(lat, na.rm = TRUE) - 0.1,
  ymax = max(lat, na.rm = TRUE) + 0.1
)

my_points.df <- data.frame(lon = long, lat = lat)
my_points.sf <- st_as_sf(my_points.df, coords = c("lon", "lat"), crs = 4326)

my_bbox <- c(
  xmin = min(my_points.df$lon, na.rm = TRUE),
  xmax = max(my_points.df$lon, na.rm = TRUE),
  ymin = min(my_points.df$lat, na.rm = TRUE),
  ymax = max(my_points.df$lat, na.rm = TRUE)
)

my_bbox.m <- matrix(
  c(
    my_bbox["xmin"], my_bbox["xmin"], my_bbox["xmax"], my_bbox["xmax"], my_bbox["xmin"],
    my_bbox["ymax"], my_bbox["ymin"], my_bbox["ymin"], my_bbox["ymax"], my_bbox["ymax"]
  ),
  ncol = 2
)

my_bbox.sf <- st_geometry(st_polygon(x = list(my_bbox.m)))
st_crs(my_bbox.sf) <- 4326

my_bbox_buff_25000.sf <- my_bbox.sf %>%
  st_transform(crs = 32632) %>%
  st_buffer(dist = 25000) %>%
  st_transform(crs = 4326)

# ---- OSM rivers ----
osm_rivers.sf <- opq(bbox = st_bbox(my_bbox_buff_25000.sf)) %>%
  add_osm_feature(key = "waterway", value = "river") %>%
  osmdata_sf()

osm_rivers.sf2 <- osm_rivers.sf$osm_lines

riv <- osm_rivers.sf2[osm_rivers.sf2$osm_id == 720489447, ] #manually IDing rivers in the region

#2. --- Function to prepare seasonal site means ---

prepare_sites <- function(dat) {
  pdx <- dat[dat$Sal < 0.5, ]
  if (!"DOC" %in% names(pdx) && "DOC.conc..µM." %in% names(pdx)) {
    pdx$DOC <- pdx$DOC.conc..µM.
  }
  
  pdx_sites <- pdx %>%
    dplyr::group_by(Longitude, Latitude) %>%
    dplyr::summarise(
      DOC = mean(DOC, na.rm = TRUE),
      across(
        -c(DOC),
        ~ if(is.numeric(.x)) mean(.x, na.rm = TRUE) else first(.x)
      ),
      .groups = "drop"
    )
  
  pdx_sites$CH4 <- pdx_sites$CH4..nM.
  
  return(pdx_sites)
}


# ---- Inter-monsoon and monsoon site means ----
pdx_sites_pre  <- prepare_sites(podat1)
pdx_sites_post <- prepare_sites(podat2)

# ---- Offset closely spaced Simunjan points for plotting only ----
# Offsets are for visual separation only; original coordinates are unchanged.
sim_offsets <- tibble::tribble(
  ~River,      ~Season,        ~Station, ~lon_offset, ~lat_offset,
  "Simunjan", "Pre-monsoon",   1,         0.000,       0.000,
  "Simunjan", "Pre-monsoon",   2,        -0.004,      -0.020,
  "Simunjan", "Pre-monsoon",   3,         0.004,       0.070,
  "Simunjan", "Pre-monsoon",   4,         0.004,       0.006,
  "Simunjan", "Post-monsoon",  1,         0.000,      -0.006,
  "Simunjan", "Post-monsoon",  2,        -0.004,      -0.020,
  "Simunjan", "Post-monsoon",  3,         0.004,       0.020,
  "Simunjan", "Post-monsoon",  4,         0.004,       0.006,
  "Simunjan", "Post-monsoon",  5,         0.006,       0.006
)

add_plot_coords <- function(dat) {
  dat %>%
    left_join(sim_offsets, by = c("River", "Season", "Station")) %>%
    mutate(
      lon_offset = tidyr::replace_na(lon_offset, 0),
      lat_offset = tidyr::replace_na(lat_offset, 0),
      Longitude_plot = Longitude + lon_offset,
      Latitude_plot  = Latitude + lat_offset
    )
}

pdx_sites_pre  <- add_plot_coords(pdx_sites_pre)
pdx_sites_post <- add_plot_coords(pdx_sites_post)

# ---- Remove NA values only for plotting specific variables ----
pdx2_pre  <- pdx_sites_pre[!is.na(pdx_sites_pre$CH4..nM.), ]
pdx2_post <- pdx_sites_post[!is.na(pdx_sites_post$CH4..nM.), ]

pdx3_pre  <- pdx_sites_pre[!is.na(pdx_sites_pre$N2O..nM.), ]
pdx3_post <- pdx_sites_post[!is.na(pdx_sites_post$N2O..nM.), ]


#3. Figure setup and plotting function --------------------------------------


axis_text_size <- 16
axis_title_size <- 18
title_size <- 20
legend_text_size <- 15
legend_title_size <- 18
point_size <- 9
label_size <- 8

big_theme <- theme(
  axis.text = element_text(size = axis_text_size, color = "black"),
  axis.title = element_text(size = axis_title_size, face = "bold"),
  plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = "black"),
  text = element_text(size = axis_text_size)
)

legend_theme <- theme(
  legend.position = "right",
  legend.key.size = unit(1.2, "lines"),
  legend.text = element_text(size = legend_text_size),
  legend.title = element_text(size = legend_title_size, face = "bold")
)

make_map_plot <- function(data, fill_var, fill_limits = NULL, fill_name,
                          plot_title = NULL, panel_label,
                          show_x = FALSE, show_y = TRUE,
                          add_legend = FALSE) {
  
  p <- ggplot() +
    geom_sf(data = riv, color = "skyblue", lwd = 0.5, alpha = 1) +
    geom_sf(data = asia_cropped, fill = "gray46") +
    geom_sf(data = osm_rivers.sf2, colour = "skyblue", lwd = 0.5, alpha = 1) +
    geom_point(
      data = data,
      aes_string(x = "Longitude_plot", y = "Latitude_plot", fill = fill_var),
      pch = 21,
      color = "black",
      size = point_size
    )+
    scale_fill_continuous(
      type = "viridis",
      limits = fill_limits,
      name = fill_name,
      oob = scales::squish
    ) +
    guides(
      fill = guide_colorbar(
        title.theme = element_text(size = legend_title_size, face = "bold")
      )
    ) +
    scale_x_continuous(
      limits = c(110.73, 111.1),
      breaks = seq(110.75, 111.1, by = 0.15),
      labels = function(x) paste0(x, "\u00B0", "E")
    ) +
    scale_y_continuous(
      limits = c(1.2, 1.7),
      labels = function(x) paste0(x, "\u00B0", "N")
    ) +
    labs(
      x = if (show_x) "Longitude (°E)" else NULL,
      y = if (show_y) "Latitude (°N)" else NULL,
      title = plot_title
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
    )
  
  if (add_legend) {
    p <- p + legend_theme
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}


# 4. Create plots ---------------------------------------------------------

p1 <- make_map_plot(
  pdx_sites_pre, "DOC", c(1000, 4387),
  expression(bold(DOC ~ " (" * mu * "M)")),
  "Intermonsoon", "A",
  show_x = FALSE, show_y = TRUE
)

p2 <- make_map_plot(
  pdx_sites_post, "DOC", c(1000, 4387),
  expression(bold(DOC ~ " (" * mu * "M)")),
  "Monsoon", "B",
  show_x = FALSE, show_y = FALSE, add_legend = TRUE
)

p3 <- make_map_plot(
  pdx2_pre, "CH4", c(400, 2000),
  expression(bold(CH[4] * " (nM)")),
  NULL, "C",
  show_x = FALSE, show_y = TRUE
)

p4 <- make_map_plot(
  pdx2_post, "CH4", c(400, 2000),
  expression(bold(CH[4] * " (nM)")),
  NULL, "D",
  show_x = FALSE, show_y = FALSE, add_legend = TRUE
)

p5 <- make_map_plot(
  pdx3_pre, "N2O..nM.", c(2.5, 10),
  expression(bold(N[2] * O * " (nM)")),
  NULL, "E",
  show_x = FALSE, show_y = TRUE
)

p6 <- make_map_plot(
  pdx3_post, "N2O..nM.", c(2.5, 10),
  expression(bold(N[2] * O * " (nM)")),
  NULL, "F",
  show_x = TRUE, show_y = FALSE, add_legend = TRUE
)

p7 <- make_map_plot(
  pdx_sites_pre %>% filter(!is.nan(CO2..ppm.)),
  "CO2..ppm.",
  c(3000, 8000),
  expression(bold(CO[2]*" (ppm)")),
  NULL,
  "G",
  show_x = TRUE,
  show_y = TRUE,
  add_legend = TRUE
)
blank <- ggplot() + theme_void()

##NO CO2 for monsoon

# 5. Arrange figure grid and export ---------------------------------------

left_col <- plot_grid(
  p1, p3, p5, p7,
  ncol = 1,
  align = "v"
)

right_col <- plot_grid(
  p2, p4, p6, blank,
  ncol = 1,
  align = "v"
)

left_col <- left_col + theme(plot.margin = margin(0, -8, 0, 0))
right_col <- right_col + theme(plot.margin = margin(0, 0, 0, -8))

panel_grid <- plot_grid(
  left_col,
  right_col,
  ncol = 2,
  align = "hv",
  rel_widths = c(1, 1)
)

pdf("Fig2.pdf", width = 14, height = 26)
print(panel_grid)
dev.off()
