
# Tables, statistics and ranges listed in Malaysian rivers DOM paper ------------------------------------
##note my labeling is now referred to differently in manuscript : pre-monsoon is actually SW monsoon, post-monsoon is NE monsoon.

par(mfrow = c(1, 1))

# ---- Libraries ----
library(tidyverse)
library(vegan)
library(geosphere)

library(conflicted)

conflicts_prefer(
  dplyr::select,
  dplyr::filter,
  dplyr::lag
)

#1. ---- Load data ----
# Run from the directory containing the input CSV files.

data <- read.csv("Metadata_Malaysia_LO.csv", sep = ",", header = TRUE, fill = TRUE)
ct <- read.csv("DOM_Malaysia_LO_new.csv", sep = ",", header = TRUE, fill = TRUE)

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
podat$weighted.NC <- colSums((ch$N / ch$C) * dom)
podat$weighted.SC <- colSums((ch$S / ch$C) * dom)
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

# NOSC
nosc <- -1 * (((4 * ch$C + ch$H - 3 * ch$N - 2 * ch$O + 5 * ch$P - 2 * ch$S) / ch$C) / 4)
nosc.dom <- data.frame(nosc * t(podom))
podat$NOSC <- colSums(nosc.dom)

# Study dataset (Rajang excluded above)
podat_main <- podat

ideg <- read.csv("Ideg.csv")

stopifnot(!anyDuplicated(ideg$FT_Name))
stopifnot(all(podat_main$FT_Name %in% ideg_corrected$FT_Name))

podat_main <- podat_main %>%
  select(-any_of("Ideg")) %>%
  left_join(
    ideg_corrected %>% select(FT_Name, Ideg),
    by = "FT_Name"
  )

stopifnot(all(is.finite(podat_main$Ideg)))
# Terrestrial-associated peaks (t-Peaks)
# Relative intensity as % of total assigned FT-ICR-MS peak intensity
tpeaks <- read.csv("t_peaks.csv")

tpeak_rows <- which(
  ch$formula_isotopefree %in% tpeaks$formula
)

podat_main$tpeaks_pct <- 100 * colSums(
  dom[tpeak_rows, , drop = FALSE],
  na.rm = TRUE
)

# 3. Setting up tables -------------------------------------------------------
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


# 4. Environmental summary table (this is now table 2) -----------------------------------------------------------------


table1_vars <- c(
  "Temp",
  "pH",
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

# 5. Table of DOM molecular characteristics(this is now table 3) -----------------------------------------------------------------


table2_vars <- c(
  "rich",
  "diversity",
  "weighted.OC",
  "weighted.HC",
  "weighted.NC",
  "weighted.SC",
  "NOSC",
  "CRAM",
  "highly_unsaturated",
  "N_ab",
  "S_ab",
  "N_num",
  "S_num",
  "AImod",
  "tpeaks_pct"
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
  "pH",
  "DOC",
  "Oxygen",
  "N2O..nM.",
  "CH4..nM.",
  "CO2..ppm.",
  "K",
  "CRAM",
  "NOSC",
  "Ideg",
  "weighted.OC",
  "weighted.HC",
  "weighted.NC",
  "weighted.SC",
  "AImod"
)

range_vars <- range_vars[range_vars %in% names(podat_main)]
range_table <- make_range_table(podat_main, range_vars)
print(range_table, width = Inf)

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
stats_vars <- unique(c(table1_vars, table2_vars))
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



# 9. Statistical comparisons among rivers ----------------------------------

run_river_stats <- function(dat, variable) {
  
  dat2 <- dat %>%
    select(
      River,
      Season,
      all_of(variable)
    ) %>%
    filter(
      !is.na(.data[[variable]]),
      !is.na(River),
      !is.na(Season)
    ) %>%
    mutate(
      River = droplevels(factor(River)),
      Season = droplevels(factor(Season))
    )
  
  n_rivers  <- nlevels(dat2$River)
  n_seasons <- nlevels(dat2$Season)
  
  # Need at least two rivers to make a river comparison
  if (n_rivers < 2) {
    return(
      tibble(
        Variable = variable,
        N = nrow(dat2),
        River_F = NA_real_,
        River_p = NA_real_,
        Season_F = NA_real_,
        Season_p = NA_real_,
        Maludam_vs_Sebuyau = NA_real_,
        Maludam_vs_Simunjan = NA_real_,
        Sebuyau_vs_Simunjan = NA_real_
      )
    )
  }
  
  # If both seasons are represented:
  # two-way ANOVA accounting for season
  if (n_seasons >= 2) {
    test_used <- "Two-way ANOVA"
    model <- aov(
      as.formula(paste(variable, "~ River + Season")),
      data = dat2
    )
    
    anova_out <- summary(model)[[1]]
    
    river_F  <- anova_out["River", "F value"]
    river_p  <- anova_out["River", "Pr(>F)"]
    season_F <- anova_out["Season", "F value"]
    season_p <- anova_out["Season", "Pr(>F)"]
    
  } else {
    
    # Only one season available:
    # compare rivers using one-way ANOVA
    test_used <- "One-way ANOVA"
    model <- aov(
      as.formula(paste(variable, "~ River")),
      data = dat2
    )
    
    anova_out <- summary(model)[[1]]
    
    river_F <- anova_out["River", "F value"]
    river_p <- anova_out["River", "Pr(>F)"]
    
    season_F <- NA_real_
    season_p <- NA_real_
  }
  
  
  # Tukey comparisons among rivers
  tukey <- TukeyHSD(model, "River")$River %>%
    as.data.frame() %>%
    tibble::rownames_to_column("comparison")
  
  
  get_tukey_p <- function(comp) {
    
    value <- tukey %>%
      filter(comparison == comp) %>%
      pull(`p adj`)
    
    if (length(value) == 0) NA_real_ else value
  }
  
  
  tibble(
    Variable = variable,
    N = nrow(dat2),
    Test = test_used,
    River_F = river_F,
    River_p = river_p,
    Season_F = season_F,
    Season_p = season_p,
    Maludam_vs_Sebuyau =
      get_tukey_p("Sebuyau-Maludam"),
    Maludam_vs_Simunjan =
      get_tukey_p("Simunjan-Maludam"),
    Sebuyau_vs_Simunjan =
      get_tukey_p("Simunjan-Sebuyau")
  )
}

# ---- Run statistics for all summary table variables ----------------

# Combine variables from the summary tables
supp_stats_vars <- unique(c(table1_vars, table2_vars, "Ideg"))

# Keep only variables present in dataset
supp_stats_vars <- supp_stats_vars[
  supp_stats_vars %in% names(podat_main)
]

# Run statistics
table_S_stats <- purrr::map_dfr(
  supp_stats_vars,
  ~ run_river_stats(podat_main, .x)
)

write.csv(table_S_stats, "table_S1.csv")

# ---- Format supplementary statistics table -------------------------------

format_p <- function(p) {
  
  if (is.na(p)) {
    return(NA_character_)
  }
  
  if (p < 0.001) {
    return("<0.001")
  }
  
  sprintf("%.3f", p)
}


table_S_formatted <- table_S_stats %>%
  mutate(
    River_F = round(River_F, 2),
    River_p = vapply(River_p, format_p, character(1)),
    
    Season_F = round(Season_F, 2),
    Season_p = vapply(Season_p, format_p, character(1)),
    
    Maludam_vs_Sebuyau =
      vapply(Maludam_vs_Sebuyau, format_p, character(1)),
    
    Maludam_vs_Simunjan =
      vapply(Maludam_vs_Simunjan, format_p, character(1)),
    
    Sebuyau_vs_Simunjan =
      vapply(Sebuyau_vs_Simunjan, format_p, character(1))
  )
