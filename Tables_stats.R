
# Tables, statistics and data analysis for Malaysia DOM paper ------------------------------------


# Load data and libraries -------------------------------------------------
rm(list = ls(all = TRUE))
par(mfrow = c(1, 1))

# ---- Libraries ----
library(tidyverse)
library(vegan)
library(geosphere)

# ---- Load data ----
#setwd of data folder

ct <- #read.csv of DOM crosstab file
data <- #read.csv of metadata file
  
# Convert oxygen consumption rate from µM min^-1 to µM h^-1
data$K <- data$K * 60 

# 1. Process DOM data -----------------------------------------------------
ch <- ct[, 1:83] #chemical attributes on crosstab
dom <- ct[, 84:178]#Peak intensities

dom[is.na(dom)] <- 0

#Average duplicate samples
dom3 <- dom %>%
  group_by(data$Sample_Name) %>%
  summarise_all(mean)

data2 <- data[!duplicated(data$Sample_Name), ]

# Apply method detection limit to averaged samples
mdl <- ct$MDL
dom <- dom3[, 2:ncol(dom3)]
dom[dom < mdl] <- 0

# Relative abundance normalisation
domn <- decostand(t(dom), method = "total", MARGIN = 2)
dom <- as.data.frame(domn)
colnames(dom) <- data2$Sample_Name

# Keep metadata aligned with DOM table
dom <- dom[, 1:54]
data2 <- data2[1:54, ]


# 2. DOM indices and compound classes -------------------------------------
podom <- as.data.frame(t(dom))

phet.N <- dom %>% filter(ct$N > 0)
phet.S <- dom %>% filter(ct$S > 0)

data2$S_ab <- colSums(phet.S)
data2$N_ab <- colSums(phet.N)

het.dom <- rbind.data.frame(phet.N, phet.S)
hetd <- as.data.frame(t(het.dom))

podat <- data2

podat$rich <- specnumber(podom)
podat$diversity <- diversity(podom, index = "shannon")
podat$Het_Perc <- rowSums(hetd) * 100
podat$weighted.OC <- colSums(ct$O.C * dom)
podat$weighted.HC <- colSums(ct$H.C * dom)

# Use AImod if present; otherwise fall back to AI
if ("AI.mod" %in% names(ct)) {
  podat$AImod <- colSums(ct$AI.mod * dom)
} else {
  podat$AImod <- colSums(ct$AI * dom)
}

podat$highly_unsaturated <- colSums(ct$Highly.unsaturated * dom)
podat$N_containing <- colSums(ct$Unsaturated.with.N * dom)
podat$saturated <- colSums(ct$Saturated * dom)

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

podat_main <- podat %>%
  filter(River %in% study_rivers, Sal <= 0.5)


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

write.csv(table1_summary, "table1_R.csv")

#5.  Table 2 -----------------------------------------------------------------


table2_vars <- c(
  "rich",
  "diversity",
  "weighted.OC",
  "weighted.HC",
  "MLB",
  "CRAM",
  "N_containing",
  "S_containing",
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

stats_table_all <- stats_table_all %>%
  mutate(
    ANOVA_p_round = signif(ANOVA_p, 3),
    Tukey_p_adj_round = signif(Tukey_p_adj, 3),
    significant = Tukey_p_adj < 0.05
  )

write.csv(stats_table_all, "river_statistics_p_values.csv", row.names = FALSE)


# 9. Linear regressions described in text ---------------------------------
lm_MLB_NOSC <- lm(MLB ~ NOSC, data = podat_main)
summary(lm_MLB_NOSC)

lm_CRAM_lignin <- lm(CRAM ~ lignin_like, data = podat_main)
summary(lm_CRAM_lignin)

# Optional export of regression summaries as text
sink("Regression_summaries.txt")
cat("MLB ~ NOSC\n")
print(summary(lm_MLB_NOSC))
cat("\n\nCRAM ~ lignin_like\n")
print(summary(lm_CRAM_lignin))
sink()

