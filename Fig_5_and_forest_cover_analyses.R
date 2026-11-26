# Forest cover, DOM composition and greenhouse gases
# Figure 5 and Table S2

library(tidyverse)
library(vegan)

# Load data ----

data <- read.csv("Metadata_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)
ct <- read.csv("DOM_Malaysia_LO_new.csv", sep = ",", header = TRUE, fill = TRUE)
forest <- read.csv(
  "sampling_locations_forest_change_1973_2015.csv"
)

data_landuse <- data %>%
  left_join(
    forest %>%
      dplyr::select(
        River, Latitude, Longitude,
        Station_ID, WS_ID,
        Forest_remaining, Forest_lost
      ),
    by = c("River", "Latitude", "Longitude")
  )

ch <- ct[, 1:84]
dom <- ct[, 85:ncol(ct)]

stopifnot(ncol(dom) == nrow(data))
colnames(dom) <- data$FT_Name
rownames(dom) <- ch$mz

# Check that the forest cover join did not duplicate samples
stopifnot(nrow(data_landuse) == nrow(data))

# DOM molecular properties ----
# Relative abundance normalisation
domn <- decostand(dom, method = "total", MARGIN = 2)
podom <- as.data.frame(t(domn))

keep <- data$River != "Rajang"

podom <- podom[keep, ]
dom <- as.data.frame(t(podom))
data <- data_landuse[keep, , drop = FALSE]

podat <- data

podat$weighted.OC <- colSums(ch$O.C * dom)
podat$weighted.HC <- colSums(ch$H.C * dom)
podat$weighted.NC <- colSums((ch$N / ch$C) * dom)
podat$weighted.SC <- colSums((ch$S / ch$C) * dom)

podat$AImod <- colSums(ch$AI.mod * dom)

# CRAM
dbe <- 1 + (2 * ch$C - ch$H + ch$N + ch$P) / 2

cram_keep <- dbe / ch$C >= 0.3 & dbe / ch$C <= 0.68 &
  dbe / ch$H >= 0.2 & dbe / ch$H <= 0.95 &
  dbe / ch$O >= 0.77 & dbe / ch$O <= 1.75
cram_keep[is.na(cram_keep)] <- FALSE
cram.dom <- dom[cram_keep, , drop = FALSE]

podat$CRAM <- colSums(cram.dom)

# NOSC
nosc <- -1 * (((4 * ch$C + ch$H - 3 * ch$N - 2 * ch$O + 5 * ch$P - 2 * ch$S) / ch$C) / 4)
nosc.dom <- data.frame(nosc * t(podom))
podat$NOSC <- colSums(nosc.dom)

ideg <- read.csv("ideg.csv")


podat <- podat %>%
  left_join(ideg, by = "FT_Name")

stopifnot(nrow(podat) == 28)
stopifnot(!anyNA(podat$Ideg))

# t-Peak relative intensity (%) from most-likely-match formula assignments
# dom has already been normalised to total peak intensity for each sample.
tpeaks <- read.csv("t_peaks.csv")
tpeak_rows_likely <- which(ch$formula_isotopefree %in% tpeaks$formula)

podat$tpeaks_pct <- 100 * colSums(
  dom[tpeak_rows_likely, , drop = FALSE],
  na.rm = TRUE
)

# Heteroatom classes from Figure 4 ----

heteroatom_data <- readRDS(
  "Fig_4_LO.R"
)

podat_corr <- podat %>%
  left_join(
    heteroatom_data,
    by = "Sample_Name",
    relationship = "one-to-one"
  )

# Table S2: forest cover vs DOM properties ----

dom_vars_corr <- c(
  "weighted.OC",
  "weighted.HC",
  "weighted.NC",
  "weighted.SC",
  "AImod",
  "CRAM",
  "NOSC",
  "RA_CHON",
  "RA_CHOS",
  "RA_CHONS",
  "tpeaks_pct",
  "Ideg"
)

forest_dom_corr <- purrr::map_dfr(dom_vars_corr, function(v) {
  
  d <- podat_corr %>%
    dplyr::select(Forest_remaining, dplyr::all_of(v)) %>%
    tidyr::drop_na()
  
  test <- cor.test(
    d$Forest_remaining,
    d[[v]],
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    DOM_property = v,
    n = nrow(d),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

print(forest_dom_corr, n = Inf)
write.csv(forest_dom_corr, "Table_S2_forest_DOM_correlations_with_tpeaks_Ideg.csv", row.names = FALSE)

# Forest cover vs greenhouse gases ----

ghg_vars <- c("CO2..ppm.", "CH4..nM.", "N2O..nM.")

forest_ghg_corr <- purrr::map_dfr(ghg_vars, function(v) {
  
  d <- podat_corr %>%
    dplyr::select(Forest_remaining, dplyr::all_of(v)) %>%
    tidyr::drop_na()
  
  test <- cor.test(
    d$Forest_remaining,
    d[[v]],
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    GHG = v,
    n = nrow(d),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_value)

print(forest_ghg_corr, n = Inf)

write.csv(forest_ghg_corr, "Table_S2_forest_GHG_correlations.csv", row.names = FALSE)

# Figure 5: Forest cover vs DOM molecular characteristics

fig_vars <- c(
  "CRAM",
  "NOSC",
  "weighted.OC",
  "Ideg"
)

fig_labels <- c(
  "A.  CRAM (%)",
  "B.  NOSC",
  "C.  Weighted O/C",
  "D.  I\u1d05\u1d07\u0262"
)

# Prepare plotting data
fig_data <- podat_corr %>%
  dplyr::select(
    River, Season, Forest_remaining,
    dplyr::all_of(fig_vars)
  ) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(fig_vars),
    names_to = "Property",
    values_to = "Value"
  ) %>%
  mutate(
    Value = if_else(Property == "CRAM", Value * 100, Value),
    Property = factor(
      Property,
      levels = fig_vars,
      labels = fig_labels
    ),
    Season = factor(
      Season,
      levels = c("Pre-monsoon", "Post-monsoon"),
      labels = c("SW monsoon", "NE monsoon")
    )
  )

fig_stats <- forest_dom_corr %>%
  dplyr::filter(DOM_property %in% fig_vars)
print(fig_stats)

# Plot
p_forest <- ggplot(
  fig_data,
  aes(
    Forest_remaining, Value,
    colour = River,
    shape = Season
  )
) +
  geom_point(size = 2.7, alpha = 0.85) +

  facet_wrap(
    ~ Property,
    scales = "free_y",
    ncol = 2
  ) +
  
  scale_colour_manual(values = c(
    "Maludam" = "#2166AC",
    "Sebuyau" = "#1B9E77",
    "Simunjan" = "#D95F02"
  )) +
  
  scale_shape_manual(values = c(
    "SW monsoon" = 16,
    "NE monsoon" = 17
  )) +
  
  labs(
    x = "Remaining forest cover (%)",
    y = NULL,
    colour = "River",
    shape = "Season"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      hjust = 0,
      size = 12
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    ),
    legend.position = "bottom"
  )

print(p_forest)

ggsave(
  "Figure5_forest_cover_DOM.png",
  p_forest,
  width = 8,
  height = 6,
  dpi = 300
)
