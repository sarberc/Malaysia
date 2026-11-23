# Figure 4: Heteroatom-containing molecular formula classes ----
# Relative molecular abundance of CHON, CHOS and CHONS formulas

# 1. Libraries ----
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(grid)

# 2. Load data ----
data <- read.csv("Metadata_Malaysia_LO.csv")
ct <- read.csv("DOM_Malaysia_LO_new.csv")

# Molecular formula information and sample intensities
ch <- ct[, 1:84]
dom <- ct[, 85:ncol(ct)]
colnames(dom) <- data$FT_Name

# 3. Calculate relative molecular abundance ----
# Normalize intensities within each sample
dom_rel <- vegan::decostand(dom, method = "total", MARGIN = 2)

# Retain the three study rivers
keep <- data$River %in% c("Maludam", "Sebuyau", "Simunjan")
data_fig <- data[keep, ]
dom_rel <- dom_rel[, keep, drop = FALSE]

# Mutually exclusive formula classes
is_CHON  <- ch$N > 0 & ch$S == 0
is_CHOS  <- ch$N == 0 & ch$S > 0
is_CHONS <- ch$N > 0 & ch$S > 0

# Sum normalized intensities and express as percentages
data_fig$CHON_percent <- colSums(dom_rel[is_CHON, , drop = FALSE]) * 100
data_fig$CHOS_percent <- colSums(dom_rel[is_CHOS, , drop = FALSE]) * 100
data_fig$CHONS_percent <- colSums(dom_rel[is_CHONS, , drop = FALSE]) * 100

# 4. Prepare plotting data ----
fig4_dat <- data_fig %>%
  select(FT_Name, River, Season,
         CHON_percent, CHOS_percent, CHONS_percent) %>%
  pivot_longer(
    cols = c(CHON_percent, CHOS_percent, CHONS_percent),
    names_to = "Class",
    values_to = "Relative_abundance"
  ) %>%
  mutate(
    Class = recode(
      Class,
      CHON_percent = "CHON",
      CHOS_percent = "CHOS",
      CHONS_percent = "CHONS"
    ),
    River = factor(River, levels = c("Maludam", "Sebuyau", "Simunjan")),
    Season = factor(
      Season,
      levels = c("Pre-monsoon", "Post-monsoon"),
      labels = c("SW monsoon", "NE monsoon")
    ),
    Class = factor(Class, levels = c("CHON", "CHOS", "CHONS"))
  )

panel_labels <- data.frame(
  Class = factor(c("CHON", "CHOS", "CHONS"),
                 levels = c("CHON", "CHOS", "CHONS")),
  label = c("A", "B", "C")
)

# 5. Plot Figure 4 ----
fig4 <- ggplot(
  fig4_dat,
  aes(x = River, y = Relative_abundance, fill = Season)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.75,
    colour = "black",
    linewidth = 0.5
  ) +
  geom_point(
    aes(group = Season),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    shape = 21,
    size = 3,
    stroke = 0.4,
    colour = "black",
    alpha = 0.9
  ) +
  facet_wrap(~ Class, nrow = 1, scales = "free_y") +
  geom_text(
    data = panel_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.4,
    vjust = 1.2,
    size = 5.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "SW monsoon" = "#3F7C78",
      "NE monsoon" = "#E07A4F"
    )
  ) +
  labs(
    x = NULL,
    y = "Relative molecular abundance (%)",
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    axis.text = element_text(colour = "black", size = 13),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.y = element_text(size = 14, margin = margin(r = 8)),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 13),
    panel.spacing = unit(1.2, "lines")
  )

fig4
