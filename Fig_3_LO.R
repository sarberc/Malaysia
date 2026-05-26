# Figure 3: DOM molecular composition maps
# CRAM, MLB, and AImod
# Inter-monsoon vs Monsoon

rm(list = ls(all = TRUE))
par(mfrow = c(1, 1))

# ---- Libraries ----
library(tidyverse)
library(vegan)
library(sf)
library(cowplot)
library(grid)
library(rnaturalearth)
library(rnaturalearthdata)
library(osmdata)
library(scales)

#1. ---- Load data ----
#setwd()

#data <- read.csv("Metadata_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)
#ct <- read.csv("DOM_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)

ch = ct[,1:84]
dom = ct[,85:ncol(ct)]

# 2. DOM indices and compound classes -------------------------------------
# Relative abundance normalisation
domn <- decostand(dom, method = "total", MARGIN = 2)
podom <- as.data.frame(t(domn))

keep <- data$River != "Rajang"

podom <- podom[keep, ]
dom <- as.data.frame(t(podom))
data <- data[keep, ]

data = data %>% filter(River != "Rajang")


phet.N <- dom %>% filter(ch$N > 0)
phet.S <- dom %>% filter(ch$S > 0)
count.N = phet.N
count.N[count.N>0]=1
count.S = phet.S
count.S[count.S>0]=1


data$S_ab <- colSums(phet.S)
data$N_ab <- colSums(phet.N)
data$N_num=colSums(count.N)
data$S_num=colSums(count.S)


het.dom <- rbind.data.frame(phet.N, phet.S)
hetd <- as.data.frame(t(het.dom))

podat <- data

podat <- podat %>%
  mutate(
    Station = as.numeric(stringr::str_extract(FT_Name, "(?<=st)\\d+"))
  )

podat$rich <- specnumber(podom)
podat$diversity <- diversity(podom, index = "shannon")
podat$Het_Perc <- rowSums(hetd) * 100
podat$weighted.OC <- colSums(ch$O.C * dom)
podat$weighted.HC <- colSums(ch$H.C * dom)
podat$highly_unsaturated <- colSums(ch$Highly.unsaturated * dom)
podat$N_containing <- data$N_ab*100
podat$AImod = colSums(ch$AI.mod*dom)

# CRAM
dbe <- 1 + (2 * ch$C - ch$H + ch$N + ch$P) / 2

cram.dom <- as.data.frame(t(podom)) %>%
  filter(
    dbe / ch$C >= 0.3,
    dbe / ch$C <= 0.68,
    dbe / ch$H >= 0.2,
    dbe / ch$H <= 0.95,
    dbe / ch$O >= 0.77,
    dbe / ch$O <= 1.75
  )

podat$CRAM <- colSums(cram.dom)

# Molecular lability boundary
mlb <- dom %>% filter(ch$H.C >= 1.5)
podat$MLB <- colSums(mlb) * 100

# NOSC
nosc <- -1 * (((4 * ch$C + ch$H - 3 * ch$N - 2 * ch$O + 5 * ch$P - 2 * ch$S) / ch$C) / 4)
nosc.dom <- data.frame(nosc * t(podom))
podat$NOSC <- colSums(nosc.dom)

#Compound classes
lignin <- as.data.frame(t(podom)) %>%
  filter(ch$O.C >= 0.29, ch$O.C <= 0.65, ch$H.C >= 0.7, ch$H.C < 1.5)

podat$lignin_like <- colSums(lignin)

# Keep only main study rivers for tables/statistics
study_rivers <- c("Maludam", "Sebuyau", "Simunjan")

podat_main <- podat 

#3. Prepare data for mapping ------------------------------------------------


podat <- podat[!is.na(podat$Latitude), ]
podat$DOC <- podat$DOC.conc..µM.

lat <- podat$Latitude
long <- podat$Longitude

podat_pre <- podat[podat$Season == "Pre-monsoon", ]   # Inter-monsoon
podat_post <- podat[podat$Season == "Post-monsoon", ]  # Monsoon

pdx_sites_pre <- podat_pre %>%
  dplyr::group_by(River, Season, Station, Longitude, Latitude) %>%
  dplyr::summarise(
    CRAM = mean(CRAM * 100, na.rm = TRUE),
    MLB = mean(MLB, na.rm = TRUE),
    AImod = mean(AImod, na.rm = TRUE),
    .groups = "drop"
  )

pdx_sites_post <- podat_post %>%
  dplyr::group_by(River, Season, Station, Longitude, Latitude) %>%
  dplyr::summarise(
    CRAM = mean(CRAM * 100, na.rm = TRUE),
    MLB = mean(MLB, na.rm = TRUE),
    AImod = mean(AImod, na.rm = TRUE),
    .groups = "drop"
  )

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


# ---- Map boundary setup ----

sf_use_s2(FALSE)

world <- ne_countries(scale = "large", continent = "asia", returnclass = "sf")

asia_cropped <- st_crop(
  world,
  xmin = min(long, na.rm = TRUE) - 0.1,
  xmax = max(long, na.rm = TRUE) + 0.1,
  ymin = min(lat, na.rm = TRUE) - 0.1,
  ymax = max(lat, na.rm = TRUE) + 0.1
)

my_points.df <- data.frame(lon = long, lat = lat)

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

osm_rivers.sf <- opq(bbox = st_bbox(my_bbox_buff_25000.sf)) %>%
  add_osm_feature(key = "waterway", value = "river") %>%
  osmdata_sf()

osm_rivers.sf2 <- osm_rivers.sf$osm_lines
riv <- osm_rivers.sf2[osm_rivers.sf2$osm_id == 720489447, ]

# ---- Plot settings ----

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

make_dom_plot <- function(data, fill_var, fill_limits = NULL, fill_name,
                          plot_title = NULL, panel_label,
                          show_x = FALSE, show_y = TRUE,
                          add_legend = FALSE,
                          labeler = label_number(accuracy = 0.01),
                          breaks = waiver()) {
  
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
    ) +
    scale_fill_continuous(
      type = "viridis",
      limits = fill_limits,
      name = fill_name,
      labels = labeler,
      breaks = breaks,
      oob = scales::squish
    ) +
    scale_x_continuous(
      limits = c(110.73, 111.1),
      breaks = seq(110.75, 111.1, by = 0.15),
      labels = function(x) x
    ) +
    scale_y_continuous(
      limits = c(1.2, 1.7),
      labels = function(x) x
    ) +
    big_theme +
    theme(plot.margin = margin(2, 2, 2, 2)) +
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
      x = if (show_x) "Longitude (°E)" else NULL,
      y = if (show_y) "Latitude (°N)" else NULL,
      title = plot_title
    )
  
  if (!show_x) p <- p + theme(axis.title.x = element_blank())
  if (!show_y) p <- p + theme(axis.title.y = element_blank())
  
  if (add_legend) {
    p <- p + legend_theme
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# 4. Make plots ------------------------------------------------------------


p1 <- make_dom_plot(
  pdx_sites_pre, "CRAM", c(53, 70), "CRAM (%)",
  plot_title = "Intermonsoon",
  panel_label = "A",
  show_x = FALSE,
  show_y = TRUE,
  labeler = label_number(accuracy = 1, suffix = "%"),
  breaks = c(54, 58, 62, 66, 70)
)

p2 <- make_dom_plot(
  pdx_sites_post, "CRAM", c(53, 70), "CRAM (%)",
  plot_title = "Monsoon",
  panel_label = "B",
  show_x = FALSE,
  show_y = FALSE,
  add_legend = TRUE,
  labeler = label_number(accuracy = 1, suffix = "%"),
  breaks = c(54, 58, 62, 66, 70)
)

p3 <- make_dom_plot(
  pdx_sites_pre, "MLB", c(0.20, 0.5), "MLB (%)",
  plot_title = NULL,
  panel_label = "C",
  show_x = FALSE,
  show_y = TRUE,
  labeler = label_number(accuracy = 0.01, suffix = "%"),
  breaks = c(0.20, 0.30, 0.40, 0.50)
)

p4 <- make_dom_plot(
  pdx_sites_post, "MLB", c(0.20, 0.5), "MLB (%)",
  plot_title = NULL,
  panel_label = "D",
  show_x = FALSE,
  show_y = FALSE,
  add_legend = TRUE,
  labeler = label_number(accuracy = 0.01, suffix = "%"),
  breaks = c(0.20, 0.30, 0.40, 0.50)
)

p5 <- make_dom_plot(
  pdx_sites_pre, "AImod", c(0.34, 0.4214), bquote(bold(AI[mod])),
  plot_title = NULL,
  panel_label = "E",
  show_x = TRUE,
  show_y = TRUE,
  labeler = label_number(accuracy = 0.01),
  breaks = c(0.34, 0.37, 0.40, 0.42)
)

p6 <- make_dom_plot(
  pdx_sites_post, "AImod", c(0.34, 0.4214), bquote(bold(AI[mod])),
  plot_title = NULL,
  panel_label = "F",
  show_x = TRUE,
  show_y = FALSE,
  add_legend = TRUE,
  labeler = label_number(accuracy = 0.01),
  breaks = c(0.34, 0.37, 0.40, 0.42)
)


# 5. Arrange figure grid and export ---------------------------------------


blank <- ggplot() + theme_void()

left_col <- plot_grid(
  p1, p3, p5, blank,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1, 1, 1)
)

right_col <- plot_grid(
  p2, p4, p6, blank,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1, 1, 1)
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

pdf("Fig_3.pdf", width = 14, height = 26)
print(panel_grid)
dev.off()
