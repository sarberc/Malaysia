
# Tables, statistics and data analysis for Malaysia DOM paper ------------------------------------

rm(list = ls(all = TRUE))
par(mfrow = c(1, 1))

# ---- Libraries ----
library(tidyverse)
library(vegan)
library(geosphere)

#1. ---- Load data ----
setwd(...)

data <- #read.csv("Metadata_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)
ct <- #read.csv("DOM_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)

ch = ct[,1:84]
dom = ct[,85:ncol(ct)]

colnames(dom) = data$FT_Name
rownames(dom) = ch$mz

# Convert oxygen consumption rate from µM min^-1 to µM h^-1
data$K <- data$K * 60 


# 2. DOM indices and compound classes -------------------------------------
# Relative abundance normalisation
domn <- decostand(dom, method = "total", MARGIN = 2)
podom <- as.data.frame(t(domn))

keep <- data$River != "Rajang"

podom <- podom[keep, ]
dom <- as.data.frame(t(podom))
data <- data[keep, ]


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


phet.het <- dom %>% filter(ch$N > 0 & ch$S > 0)

podat <- data

podat$rich <- specnumber(podom)
podat$diversity <- diversity(podom, index = "shannon")
podat$weighted.OC <- colSums(ch$O.C * dom)
podat$weighted.HC <- colSums(ch$H.C * dom)
podat$highly_unsaturated <- colSums(ch$Highly.unsaturated * dom)
podat$AImod <- colSums(ch$AI.mod * dom)
podat$N_containing <- data$N_ab*100
podat$Het_Perc <- colSums(phet.het) * 100

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

#3. Setting up tables -------------------------------------------------------
make_summary_table <- function(dat, vars) {
  dat %>%
    dplyr::group_by(Season, River) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(vars),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n = dplyr::n(),
      .groups = "drop"
    )
}

make_range_table <- function(dat, vars) {
  dat %>%
    dplyr::group_by(Season, River) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(vars),
        list(
          min = ~ min(.x, na.rm = TRUE),
          max = ~ max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      n = n(),
      .groups = "drop"
    )
}


#4. Table 1 -----------------------------------------------------------------


table1_vars <- c(
  "Temp",
  "DOC",
 "N2O..nM.",
  "CH4..nM.",
  "CO2..ppm.",
   "K"
)

table1_vars <- table1_vars[table1_vars %in% names(podat_main)]

table1_summary <- make_summary_table(podat_main, table1_vars)

table1_pre <- table1_summary %>%
  filter(Season == "Pre-monsoon")

table1_post <- table1_summary %>%
  filter(Season == "Post-monsoon")

#5.  Table 2 -----------------------------------------------------------------


table2_vars <- c(
  "rich",
  "diversity",
  "weighted.OC",
  "weighted.HC",
  "MLB",
  "CRAM",
  "highly_unsaturated",
  "N_ab",
  "S_ab",
  "N_num",
  "S_num",
  "AImod"
 )

table2_vars <- table2_vars[table2_vars %in% names(podat_main)]

table2_summary <- make_summary_table(podat_main, table2_vars)

table2_pre <- table2_summary %>%
  filter(Season == "Pre-monsoon")

table2_post <- table2_summary %>%
  filter(Season == "Post-monsoon")

#6. Values for text : ranges ------------------------------------------------


range_vars <- c(
  "Temp",
  "DOC",
  "Oxygen",
  "N2O..nM.",
  "CH4..nM.",
  "CO2..ppm.",
  "K",
  "MLB",
  "CRAM",
  "weighted.OC",
  "weighted.HC",
  "AImod"
)

range_vars <- range_vars[range_vars %in% names(podat_main)]

range_table <- make_range_table(podat_main, range_vars)


# 7. Statistical tests ----------------------------------------------------

run_anova_tukey <- function(dat, variable) {
  
  dat2 <- dat %>%
    select(River, all_of(variable)) %>%
    filter(!is.na(.data[[variable]]), !is.na(River))
  
  if (n_distinct(dat2$River) < 2) {
    return(list(
      variable = variable,
      message = paste("Skipped:", variable, "- fewer than 2 rivers with data")
    ))
  }
  
  model <- aov(as.formula(paste(variable, "~ River")), data = dat2)
  
  list(
    variable = variable,
    anova = summary(model),
    tukey = TukeyHSD(model)
  )
}

stats_vars <- c(
  "Temp",
  "K",
  "DOC",
  "CO2..ppm.",
  "N2O..nM.",
  "CH4..nM.",
  "rich",
  "diversity",
  "weighted.OC",
  "weighted.HC",
  "MLB",
  "CRAM",
  "AImod"
)

stats_vars <- stats_vars[stats_vars %in% names(podat_main)]

# Overall river differences using all freshwater samples
stats_overall <- lapply(stats_vars, function(v) run_anova_tukey(podat_main, v))
names(stats_overall) <- stats_vars

# River differences within each season
stats_by_season <- list()

for (season_name in unique(podat_main$Season)) {
  dat_season <- podat_main %>% filter(Season == season_name)
  stats_by_season[[season_name]] <- lapply(stats_vars, function(v) run_anova_tukey(dat_season, v))
  names(stats_by_season[[season_name]]) <- stats_vars
}

# Print selected tests used in manuscript
stats_overall$K
stats_overall$CRAM
stats_overall$weighted.OC

if ("Post-monsoon" %in% names(stats_by_season)) {
  stats_by_season[["Post-monsoon"]]$CRAM
  stats_by_season[["Post-monsoon"]]$weighted.OC
}

if ("Pre-monsoon" %in% names(stats_by_season)) {
  stats_by_season[["Pre-monsoon"]]$CRAM
  stats_by_season[["Pre-monsoon"]]$weighted.OC
}

# 8. Statistical results table ---------------------------------------------

run_stats_table <- function(dat, variable, season_label = "Overall") {
  
  dat2 <- dat %>%
    dplyr::select(River, dplyr::all_of(variable)) %>%
    dplyr::filter(!is.na(.data[[variable]]), !is.na(River))
  
  if (dplyr::n_distinct(dat2$River) < 2) {
    return(NULL)
  }
  
  model <- aov(as.formula(paste(variable, "~ River")), data = dat2)
  
  anova_p <- summary(model)[[1]][["Pr(>F)"]][1]
  
  tukey <- TukeyHSD(model)$River %>%
    as.data.frame() %>%
    tibble::rownames_to_column("comparison") %>%
    dplyr::mutate(
      Season = season_label,
      Variable = variable,
      ANOVA_p = anova_p,
      Tukey_p_adj = `p adj`
    ) %>%
    dplyr::select(
      Season,
      Variable,
      ANOVA_p,
      comparison,
      diff,
      lwr,
      upr,
      Tukey_p_adj
    )
  
  return(tukey)
}

stats_table_overall <- purrr::map_dfr(
  stats_vars,
  ~ run_stats_table(podat_main, .x, season_label = "Overall")
)

stats_table_by_season <- purrr::map_dfr(
  unique(podat_main$Season),
  function(s) {
    dat_s <- podat_main %>% filter(Season == s)
    
    purrr::map_dfr(
      stats_vars,
      ~ run_stats_table(dat_s, .x, season_label = s)
    )
  }
)

stats_table_all <- stats_table_by_season
stats_table_all <- dplyr::bind_rows(
  stats_table_overall,
  stats_table_by_season
) %>%
  mutate(
    ANOVA_p_round = signif(ANOVA_p, 3),
    Tukey_p_adj_round = signif(Tukey_p_adj, 3),
    significant_ANOVA = ANOVA_p < 0.05,
    significant_Tukey = Tukey_p_adj < 0.05
  )



# Linear regressions and supp figs ----------------------------------------

lm_MLB_NOSC <- lm(MLB ~ NOSC, data = podat_main)
summary(lm_MLB_NOSC)

lm_CRAM_lignin <- lm(CRAM ~ lignin_like, data = podat_main)
summary(lm_CRAM_lignin)

p_MLB_NOSC <- ggplot(podat_main, aes(x = NOSC, y = MLB, fill = River)) +
  geom_point(
    pch = 21,
    colour = "black",
    size = 4,
    alpha = 0.85
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = TRUE,
    colour = "black",
    fill = "grey80"
  ) +
  theme_bw() +
  labs(
    x = "NOSC",
    y = "MLB (%)",
    fill = "River"
  ) +
  theme(
    text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

pdf("Fig_NOSC", width = 6, height = 5)
print(p_MLB_NOSC)
dev.off()


# ---- CRAM vs lignin-like molecular formulas ----
p_CRAM_lignin <- ggplot(podat_main, aes(x = lignin_like, y = CRAM, fill = River)) +
  geom_point(
    pch = 21,
    colour = "black",
    size = 4,
    alpha = 0.85
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = TRUE,
    colour = "black",
    fill = "grey80"
  ) +
  theme_bw() +
  labs(
    x = "Lignin-like molecular formulas",
    y = "CRAM",
    fill = "River"
  ) +
  theme(
    text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

pdf("Fig_CRAM_lignin.pdf", width = 6, height = 5)
print(p_CRAM_lignin)
dev.off()