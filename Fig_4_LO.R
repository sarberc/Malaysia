rm(list = ls(all = TRUE))
par(mfrow = c(1, 1))

# Libraries ------------------------------------------------------------

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
# Load data ------------------------------------------------------------

#setwd()

#data <- read.csv("Metadata_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)
#ct <- read.csv("DOM_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)

# ---- Chemistry and unnormalised DOM matrix ----------------------------

chem <- ct[, 1:84]
dom_raw <- ct[, 85:ncol(ct)]

# samples x formulas
dom_mat <- as.data.frame(t(dom_raw))
dom_mat[] <- lapply(dom_mat, as.numeric)

# add metadata columns
dom_meta <- bind_cols(
  data %>% select(Sample_Name, River, Season),
  dom_mat
)

# mean unnormalised peaks by River x Season
dom_mean <- dom_meta %>%
  dplyr::group_by(River, Season) %>%
  dplyr::summarise(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# ---- Helper function: mean raw spectra, MDL filter, normalise, top 75% ----
make_river_vk <- function(dom_mean, chem, river, season) {
  
  row_dat <- dom_mean %>%
    filter(River == river, Season == season)
  
  if (nrow(row_dat) == 0) {
    warning(paste("No data for", river, season))
    return(NULL)
  }
  
  mean_intensity <- row_dat %>%
    select(-River, -Season) %>%
    unlist(use.names = FALSE)
  
  # MDL after averaging
  mean_intensity[mean_intensity < chem$MDL] <- 0
  
  out <- chem %>%
    transmute(
      MZ = mz,
      Formula = formula,
      OC = O.C,
      HC = H.C,
      C = C,
      H = H,
      O = O,
      N = N,
      S = S,
      Intensity_raw = mean_intensity
    ) %>%
    filter(Intensity_raw > 0) %>%
    mutate(
      Intensity = Intensity_raw / sum(Intensity_raw, na.rm = TRUE) * 100,
      compound_group = case_when(
        N > 0 & S == 0 ~ "CHON",
        S > 0 & N == 0 ~ "CHOS",
        N > 0 & S > 0  ~ "CHONS",
        TRUE           ~ "CHO"
      )
    ) 
  
  return(out)
}

# Plot Van Krevelen ----------------------------------
# Plot Van Krevelen ----------------------------------
plot_vk_ns <- function(df, panel_label = NULL, show_x = TRUE, show_y = TRUE) {
  
  ggplot(df, aes(x = OC, y = HC)) +
    
    geom_point(data = df %>% filter(compound_group == "CHO"),
               shape = 21, fill = "gray75", colour = "black",
               size = 2, alpha = 0.35) +
    
    geom_point(data = df %>% filter(compound_group == "CHON"),
               shape = 21, fill = "gray25", colour = "black",
               size = 4, stroke = 0.15, alpha = 0.7) +
    
    geom_point(data = df %>% filter(compound_group == "CHOS"),
               shape = 21, fill = "purple1", colour = "black",
               size = 6, stroke = 0.15, alpha = 1) +
    
    coord_cartesian(xlim = c(0.2, 0.9), ylim = c(0.4, 1.65)) +
    
    annotate("text",
             x = 0.22, y = 1.62,
             label = panel_label,
             hjust = 0,
             vjust = 1,
             size = 10,
             fontface = "bold") +
    
    labs(
      x = ifelse(show_x, "O/C", ""),
      y = ifelse(show_y, "H/C", "")
    ) +
    
    theme_bw(base_size = 24) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "gray90"),
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 24),
      legend.position = "none",
      aspect.ratio = 0.78
    )
}

# Plots ----------------------------
ma_pre <- make_river_vk(dom_mean, chem, "Maludam", "Pre-monsoon")
ma_post <- make_river_vk(dom_mean, chem, "Maludam",  "Post-monsoon")

se_pre  <- make_river_vk(dom_mean, chem, "Sebuyau",  "Pre-monsoon")
se_post <- make_river_vk(dom_mean, chem, "Sebuyau",  "Post-monsoon")

si_pre  <- make_river_vk(dom_mean, chem, "Simunjan", "Pre-monsoon")
si_post <- make_river_vk(dom_mean, chem, "Simunjan", "Post-monsoon")

# Arrange figure --------------------------------------------------------

x1 <- plot_vk_ns(ma_pre,  "A", show_x = FALSE, show_y = TRUE)
x2 <- plot_vk_ns(ma_post, "B", show_x = FALSE, show_y = FALSE)

x3 <- plot_vk_ns(se_pre,  "C", show_x = FALSE, show_y = TRUE)
x4 <- plot_vk_ns(se_post, "D", show_x = FALSE, show_y = FALSE)

x5 <- plot_vk_ns(si_pre,  "E", show_x = TRUE, show_y = TRUE)
x6 <- plot_vk_ns(si_post, "F", show_x = TRUE, show_y = FALSE)

# ---- Column titles ----

column_titles <- ggarrange(
  ggplot() +
    annotate("text", x = 0.58, y = 0.5,
             label = "Inter-monsoon",
             size = 11) +
    theme_void(),
  
  ggplot() +
    annotate("text", x = 0.42, y = 0.5,
             label = "Monsoon",
             size = 11) +
    theme_void(),
  
  ncol = 2
)

# ---- River label rows ----

river_label_plot <- function(label) {
  ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.35,
      label = label,
      fontface = "italic",
      colour = "gray40",
      size = 11
    ) +
    coord_cartesian(clip = "off") +
    theme_void()
}

maludam_label <- river_label_plot("Maludam")
seb_label     <- river_label_plot("Sebuyau")
sim_label     <- river_label_plot("Simunjan")

# ---- Plot rows ----

row1 <- ggarrange(x1, x2, ncol = 2)
row2 <- ggarrange(x3, x4, ncol = 2)
row3 <- ggarrange(x5, x6, ncol = 2)

# ---- Manual legend ----

legend_plot <- ggplot(
  data.frame(
    compound_group = c("CHO", "CHON", "CHOS"),
    x = 1:3,
    y = 1
  ),
  aes(x = x, y = y, fill = compound_group)
) +
  geom_point(shape = 21, size = 6, colour = "black") +
  scale_fill_manual(
    values = c(
      "CHO" = "gray75",
      "CHON" = "gray25",
      "CHOS" = "purple1"
    ),
    name = NULL
  ) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 24),
    legend.key.size = unit(1.2, "cm")
  ) +
  coord_cartesian(clip = "off") +
  ylim(10, 11)

# ---- Final figure ----

fig_vk <- ggarrange(
  column_titles,
  maludam_label,
  row1,
  seb_label,
  row2,
  sim_label,
  row3,
  legend_plot,
  ncol = 1,
  heights = c(0.08, 0.08, 1, 0.08, 1, 0.08, 1, 0.16)
)

print(fig_vk)

pdf("Fig_4y.pdf", width = 14, height = 18)
print(fig_vk)
dev.off()
