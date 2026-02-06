##------------------------------------------------------###
##---------------------Fig.7-----------------------------####
##------------------------------------------------------###

## 0) Make script runnable anywhere (local / GitHub Actions)
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg) else NA_character_
script_dir <- if (!is.na(script_path)) normalizePath(dirname(script_path)) else getwd()
setwd(script_dir)

dir.create("res", showWarnings = FALSE, recursive = TRUE)

message("Working dir: ", getwd())
message("Files in example/: ", paste(list.files("."), collapse = ", "))

## 1) Package loading (robust)
load_or_stop <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "), "\n",
      "Install them (locally) or add them to GitHub Actions dependencies."
    )
  }
  invisible(TRUE)
}

## NOTE: bayesplot/rstanarm are not actually used in the code you pasted,
## so we load them conditionally to avoid CI failure.
core_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "boot", "emmeans", "ggrepel",
  "cowplot", "gridExtra", "knitr", "stringr"
)

load_or_stop(core_pkgs)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(boot)
  library(emmeans)
  library(ggrepel)
  library(cowplot)
  library(gridExtra)
  library(knitr)
  library(stringr)
})

if (requireNamespace("bayesplot", quietly = TRUE)) {
  suppressPackageStartupMessages(library(bayesplot))
} else {
  message("Note: 'bayesplot' not installed; skipping (not required for this script).")
}
if (requireNamespace("rstanarm", quietly = TRUE)) {
  suppressPackageStartupMessages(library(rstanarm))
} else {
  message("Note: 'rstanarm' not installed; skipping (not required for this script).")
}

theme_set(theme_bw(base_size = 10) + theme(legend.position = "right"))

## 2) Source helper functions (must exist in same folder)
##    - forest.plot.R
##    - my_ggbiplot.R
if (!file.exists("forest.plot.R")) stop("forest.plot.R not found in example/")
if (!file.exists("my_ggbiplot.R"))  stop("my_ggbiplot.R not found in example/")

source("forest.plot.R")
source("my_ggbiplot.R")

## 3) Input data (CSV must exist in same folder)
csv_file <- "1st invivo study young vs aged all data tables.csv"
if (!file.exists(csv_file)) stop("Data file not found: ", csv_file)

dat0 <- read.csv(csv_file, check.names = FALSE)

## convert the Mass of dead mouse to NA
dat0$MassDay6 <- suppressWarnings(as.numeric(dat0$MassDay6))

dat0 <- dat0 %>%
  mutate(
    wtLossThroughD4 = (MassDay4 - MassDay0) / MassDay0,
    wtLossThroughD5 = (MassDay5 - MassDay0) / MassDay0,
    wtLossThroughD6 = (MassDay6 - MassDay0) / MassDay0
  )

## Palettes
trts <- c("Saline", "CDN nanovax", "CpG nanovax", "MPLA nanovax", "R848 nanovax")
trt.shapePal <- c(1, 15, 17, 4, 18)
trt.colorPal <- c(
  rgb(220,105,125,maxColorValue = 255),
  rgb(47, 37, 217,maxColorValue = 255),
  rgb(112,112,112,maxColorValue = 255),
  rgb(230,159,  0,maxColorValue = 255),
  rgb( 86,180,233,maxColorValue = 255)
)

trt.Pal <- data.frame(
  trts = trts,
  trt.colorPal = trt.colorPal,
  trt.shapePal = trt.shapePal,
  stringsAsFactors = FALSE
)

arrow.length <- 1.2
label.adjust <- 1.5

##------------------------------------------------------###
##---------------------Fig.7(A,B,C,D): PCA --------------###
##------------------------------------------------------###

## Prepare response variables
data <- dat0 %>%
  mutate(
    sum.score = (Bronchial.bronchiolar.inflammation.score +
                   Eptihelial.Change.score +
                   Interstitial.inflammation.score +
                   Hemorrhage.score +
                   Edema.score)
  )

data$Treatment <- factor(data$Treatment, levels = trts)

dat.all <- data %>%
  dplyr::select(
    Age, ID, Treatment, sum.score, ViralLoadmL,
    Abs.Titre.background.0.13, Abs.Titre.Background.0.07,
    PenHDay2, PenHDay3, PenHDay4
  )

headers <- c("Age", "ID", "Treatment")

dat.all <- dat.all %>%
  mutate(
    logViral = log10(ViralLoadmL + 1),
    logAbs   = 0.5 * (log10(Abs.Titre.background.0.13) + log10(Abs.Titre.Background.0.07)),
    PenH234  = (PenHDay2 + PenHDay3 + PenHDay4) / 3,
    Histopath = sum.score
  ) %>%
  dplyr::select(
    Age, ID, Treatment, logViral, logAbs, PenH234, Histopath
  )

scale_inputs <- TRUE

pca_plots <- list()

for (age in c("Aged", "Young")) {

  temp <- dat.all %>% dplyr::filter(Age == age)

  ## drop columns that are all NA & then drop rows with any NA
  keep_cols <- vapply(temp, function(x) sum(!is.na(x)) > 0, logical(1))
  temp <- as.data.frame(temp[, keep_cols, drop = FALSE]) %>% tidyr::drop_na()

  if (nrow(temp) < 3) {
    warning("Too few rows after NA drop for age = ", age, "; skipping PCA.")
    next
  }

  res.pca <- prcomp(
    temp %>% dplyr::select(-all_of(headers)),
    center = TRUE,
    scale. = scale_inputs
  )

  ## Flip axes for consistent orientation (optional)
  if (mean(res.pca$rotation[,1]) < 0) {
    res.pca$rotation[,1] <- -res.pca$rotation[,1]
    res.pca$x[,1] <- -res.pca$x[,1]
  }
  if (mean(res.pca$rotation[,2]) < 0) {
    res.pca$rotation[,2] <- -res.pca$rotation[,2]
    res.pca$x[,2] <- -res.pca$x[,2]
  }

  temp$PC1 <- res.pca$x[,1]
  temp$PC2 <- res.pca$x[,2]

  out <- temp %>%
    dplyr::group_by(Treatment) %>%
    dplyr::summarize(
      meanPC1 = mean(PC1), sdPC1 = sd(PC1),
      meanPC2 = mean(PC2), sdPC2 = sd(PC2),
      n = dplyr::n(),
      .groups = "drop"
    )

  res0 <- out %>% dplyr::filter(Treatment == trts[1])
  out2 <- out %>%
    mutate(
      disPC1 = meanPC1 - res0$meanPC1[1],
      disPC2 = meanPC2 - res0$meanPC2[1],
      disEucli = sqrt(disPC1^2 + disPC2^2)
    ) %>%
    arrange(disPC1)

  ## Save tables (FIXED parentheses + safer filenames)
  write.csv(
    out2 %>% dplyr::select(Treatment, meanPC1, sdPC1, n, disPC1),
    file = file.path("res", paste0("PCA_combine_", age, ".csv")),
    row.names = FALSE
  )

  write.csv(
    temp,
    file = file.path("res", paste0("PCA_combine_detail_", age, ".csv")),
    row.names = FALSE
  )

  cat("\n=== ", age, " ===\n")
  print(knitr::kable(out2 %>% dplyr::select(Treatment, meanPC1, sdPC1, n, disPC1)))

  cc <- summary(res.pca)
  var_pc1 <- round(cc$importance[2,1] * 100, 1)
  var_pc2 <- round(cc$importance[2,2] * 100, 1)

  temp.Pal <- trt.Pal

  fig1 <- my_ggbiplot(
    res.pca,
    ellipse = TRUE,
    groups = factor(temp$Treatment, levels = trts),
    labels.size = 4,
    varname.size = 4,
    varname.abbrev = FALSE,
    alpha = 0,
    arrow.length = arrow.length,
    label.adjust = label.adjust
  ) +
    geom_point(
      aes(
        shape = factor(temp$Treatment, levels = trts),
        colour = factor(temp$Treatment, levels = trts)
      )
    ) +
    scale_shape_manual(values = as.numeric(temp.Pal$trt.shapePal)) +
    scale_color_manual(values = temp.Pal$trt.colorPal) +
    theme_bw(base_size = 15) +
    theme(legend.position = "right") +
    labs(
      color = "Adjuvants",
      shape = "Adjuvants",
      x = paste0("PC1 (", var_pc1, "%)"),
      y = paste0("PC2 (", var_pc2, "%)")
    ) +
    coord_equal(ratio = 0.7) +
    ggtitle(age)

  ## Save PCA plot
  ggsave(
    filename = file.path("res", paste0("Fig7_PCA_", age, ".png")),
    plot = fig1,
    width = 7, height = 6, dpi = 300
  )

  pca_plots[[age]] <- fig1
}

##------------------------------------------------------###
##---------------------Fig.7(F): weighted objective -----###
##------------------------------------------------------###

## In your pasted code, you read:
## data = read.csv("res/PP(CA_combine_detail_AgedYoung.csv")
## That file is not created in the loop above.
## We'll rebuild "Aged+Young combined detail" by binding two detail files.

detail_files <- file.path("res", c("PCA_combine_detail_Aged.csv", "PCA_combine_detail_Young.csv"))
if (!all(file.exists(detail_files))) {
  stop("Missing PCA detail file(s). Expected:\n", paste(detail_files, collapse = "\n"))
}

data_detail <- bind_rows(
  read.csv(detail_files[1], check.names = FALSE),
  read.csv(detail_files[2], check.names = FALSE)
)

## Standardize (scale) selected columns
data_detail$logAbs    <- as.numeric(scale(data_detail$logAbs))
data_detail$logViral  <- as.numeric(scale(data_detail$logViral))
data_detail$PenH234   <- as.numeric(scale(data_detail$PenH234))
data_detail$Histopath <- as.numeric(scale(data_detail$Histopath))

## Bring MassDay0 from original dat0 (match by row order is risky).
## Better: merge by ID if ID is unique.
if ("ID" %in% names(dat0) && "ID" %in% names(data_detail)) {
  mass_map <- dat0 %>% dplyr::select(ID, MassDay0)
  data_detail <- data_detail %>% left_join(mass_map, by = "ID")
} else {
  ## fallback
  data_detail$MassDay0 <- dat0$MassDay0[seq_len(nrow(data_detail))]
}

## Harmonize Treatment labels if needed
data_detail$Treatment <- as.character(data_detail$Treatment)
data_detail$Treatment[data_detail$Treatment == "Control (no treatment)"] <- "Saline"
data_detail$Treatment <- factor(data_detail$Treatment, levels = trts)

## Objective functions
obj1 <- function(logViral, logAbs, PenH234, Histopath) {
  (logViral + PenH234 + Histopath) / 3 - logAbs
}

datanew <- data_detail %>%
  mutate(
    Y1 = obj1(logViral, logAbs, PenH234, Histopath)
  )

## Linear model + pairwise contrasts
o.un <- lm(Y1 ~ Age + Treatment + Age:Treatment + MassDay0, data = datanew)

## NOTE: you used lsmeans() in your script, which is now superseded.
## emmeans package supports lsmeans() alias, but we'll use emmeans directly (more stable).
emm_overall <- emmeans(o.un, ~ Treatment)
contr_overall <- pairs(emm_overall)

## Overall forest plot
p1 <- forest.plot(contr_overall, -2, 8, nudge_x = 2, title = "Aged+Young")

## By Age forest plots
emm_by_age <- emmeans(o.un, ~ Treatment | Age)
contr_by_age <- pairs(emm_by_age)

## Split into Aged / Young
contr_df <- as.data.frame(contr_by_age)
contr_aged  <- contr_by_age[contr_df$Age == "Aged"]
contr_young <- contr_by_age[contr_df$Age == "Young"]

p2 <- forest.plot(contr_aged,  -2, 8, nudge_x = 2, title = "Aged")
p3 <- forest.plot(contr_young, -2, 8, nudge_x = 2, title = "Young")

## Save forest plots
ggsave(file.path("res", "Fig7_forest_overall.png"), plot = p1, width = 7, height = 5, dpi = 300)
ggsave(file.path("res", "Fig7_forest_Aged.png"),   plot = p2, width = 7, height = 5, dpi = 300)
ggsave(file.path("res", "Fig7_forest_Young.png"),  plot = p3, width = 7, height = 5, dpi = 300)

## Save contrasts tables
write.csv(as.data.frame(contr_overall), file.path("res", "Fig7_contrasts_overall.csv"), row.names = FALSE)
write.csv(as.data.frame(contr_by_age),  file.path("res", "Fig7_contrasts_by_age.csv"),  row.names = FALSE)

##------------------------------------------------------###
##---------------------Fig.7(E): point + CI plot --------###
##------------------------------------------------------###

res_means <- emmeans(o.un, ~ Treatment | Age)
plot_res <- as.data.frame(res_means)

pE <- ggplot(plot_res, aes(color = Treatment, group = interaction(Age, Treatment))) +
  geom_point(aes(Treatment, datanew$Y1[match(plot_res$Treatment, datanew$Treatment)]), alpha = 0.0) + # keep structure; no raw overlay
  geom_point(aes(Treatment, emmean), size = 2) +
  geom_errorbar(aes(x = Treatment, ymin = lower.CL, ymax = upper.CL), width = 0.4) +
  theme_bw(base_size = 12) +
  ylab("Y = 1/3(logViral + PenH234 + Histopath) - logAbs") +
  scale_color_manual(name = "", labels = trts, values = trt.colorPal) +
  facet_wrap(~Age) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))

ggsave(file.path("res", "Fig7_E_objective_plot.png"), plot = pE, width = 9, height = 4, dpi = 300)

message("Done. Outputs saved to: ", normalizePath("res"))
