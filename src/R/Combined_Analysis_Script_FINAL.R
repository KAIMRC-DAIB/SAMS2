################################################################################
# COMBINED ANALYSIS SCRIPT — FINAL 26/06
################################################################################
library(vegan)
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(stringr)
library(tibble)
library(patchwork)
library(writexl)
library(emmeans)
library(pheatmap)
library(ppcor)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(randomForest)
library(conflicted)
library(cowplot)
conflicted::conflicts_prefer(
  dplyr::filter, dplyr::select, dplyr::mutate, dplyr::arrange,
  dplyr::summarise, dplyr::rename, dplyr::slice_head, dplyr::pull,
  pheatmap::pheatmap,
  ggplot2::margin
)
################################################################################
# 1. SETUP — load data, recode Group, define confounders
################################################################################
## ---- 1.1 File paths -- EDIT THESE ----
mp_file   <- "~/Desktop/metaphlan_species_only_prop.55plus.xlsx"
meta_file <- "~/Desktop/Metadata_55plus _179_sample.xlsx"
path_file <- "~/Desktop/pathabundance_relab_named_clean.xlsx"
## ---- 1.2 Colors ----
group_colors <- c("Healthy" = "#1A9850", "Resilient" = "#2166AC", "Unhealthy" = "#B2182B")
palette_groups <- c(Healthy = "#1A9850", Resilient = "#2C7BB6", Unhealthy = "#D7191C")
## ---- 1.3 Taxonomic side: MetaPhlAn abundance ----
metaphlan_raw <- readxl::read_excel(mp_file)
tax_col <- names(metaphlan_raw)[1]
metaphlan_df <- as.data.frame(metaphlan_raw)
rownames(metaphlan_df) <- metaphlan_df[[tax_col]]
metaphlan_df[[tax_col]] <- NULL
abund_mat <- as.matrix(metaphlan_df)
storage.mode(abund_mat) <- "numeric"
if (max(abund_mat, na.rm = TRUE) > 1) abund_mat <- abund_mat / 100
otu_table <- abund_mat[grepl("s__", rownames(abund_mat)), , drop = FALSE]
otu_table <- t(otu_table)
## ---- 1.4 Metadata ----
meta <- readxl::read_excel(meta_file)
meta <- meta %>%
  mutate(unified_sample_ID = stringr::str_trim(as.character(unified_sample_ID))) %>%
  rename_with(make.names)
meta$SampleID <- as.character(meta$unified_sample_ID)
if (!"Group" %in% colnames(meta)) stop("Column 'Group' was not found in metadata.")
meta$Group_code <- meta$Group
meta$Group <- factor(
  as.numeric(stringr::str_trim(as.character(meta$Group_code))),
  levels = c(1, 2, 3),
  labels = c("Healthy", "Resilient", "Unhealthy")
)
if (any(is.na(meta$Group))) warning("Some Group values became NA. Check Group coding in metadata.")
cat("\n[Setup] Group counts (meta):\n")
print(table(meta$Group, useNA = "ifany"))
## ---- 1.5 Confounders: taxonomic side ----
confounders_taxa <- c("BMI", "Socialized", "Smoking", "wakefulness", "Sex")
missing_confounders_taxa <- setdiff(confounders_taxa, colnames(meta))
if (length(missing_confounders_taxa) > 0) {
  stop(paste("Missing taxonomic confounders:", paste(missing_confounders_taxa, collapse = ", ")))
}
numeric_vars_taxa <- c("BMI")
categorical_vars_taxa <- setdiff(confounders_taxa, numeric_vars_taxa)
meta2 <- meta %>%
  dplyr::select(SampleID, Group, dplyr::all_of(confounders_taxa)) %>%
  dplyr::mutate(BMI = as.numeric(BMI))
for (v in categorical_vars_taxa) meta2[[v]] <- factor(meta2[[v]])
## ---- 1.6 Match samples: taxonomic side ----
common_ids_taxa <- intersect(colnames(abund_mat), meta2$SampleID)
if (length(common_ids_taxa) < 5) stop("Too few matching sample IDs between MetaPhlAn and metadata.")
abund_mat <- abund_mat[, common_ids_taxa, drop = FALSE]
meta2 <- meta2 %>% dplyr::filter(SampleID %in% common_ids_taxa)
meta2 <- meta2[match(common_ids_taxa, meta2$SampleID), ]
cat("\n[Setup] Taxonomic side -- matched samples:", length(common_ids_taxa), "\n")
## ---- 1.7 Functional side: pathway abundance (reuses meta) ----
path_raw <- readxl::read_excel(path_file)
colnames(path_raw) <- trimws(colnames(path_raw))
path_anno <- path_raw %>% dplyr::select(pathway_id, feature, metacyc_link)
# Authoritative list of functional-feature (pathway) column names.
# Defined here so every downstream section (Fig 1 sec 2.7, etc.) can use it.
pathways <- path_anno$pathway_id
path_mat <- path_raw %>% dplyr::select(-feature, -pathway_id, -metacyc_link)
rownames(path_mat) <- path_raw$pathway_id
path_t <- as.data.frame(t(path_mat))
path_t$SampleID <- rownames(path_t)
data_path_raw <- merge(meta, path_t, by = "SampleID")
# Group recode -- the ONLY place this should happen. Uses Group_code,
# not Group, as the source column. If Group ever shows as all-NA after
# this point, re-run THIS block, not a copy of it elsewhere.
data_path_raw$Group <- factor(
  as.numeric(stringr::str_trim(as.character(data_path_raw$Group_code))),
  levels = c(1, 2, 3),
  labels = c("Healthy", "Resilient", "Unhealthy")
)
data_path_raw$Sex <- factor(data_path_raw$Sex, levels = c(1, 2), labels = c("Female", "Male"))
data_path_raw$Smoking     <- factor(data_path_raw$Smoking)
data_path_raw$wakefulness <- factor(data_path_raw$wakefulness)
data_path_raw$Socialized  <- factor(data_path_raw$Socialized)
data_path_raw$BMI         <- as.numeric(data_path_raw$BMI)
# Sanity check: how many pathway columns survived the merge intact?
pathways_present <- intersect(pathways, colnames(data_path_raw))
cat("\n[Setup] pathways defined:", length(pathways),
    "| present in data_path_raw:", length(pathways_present), "\n")
if (length(pathways_present) < length(pathways)) {
  warning("[Setup] ", length(pathways) - length(pathways_present),
          " pathway name(s) did not survive the merge into data_path_raw. ",
          "Downstream functional analyses use only the present ones.")
}
cat("\n[Setup] Functional side -- Group counts:\n")
print(table(data_path_raw$Group, useNA = "ifany"))
if (any(is.na(data_path_raw$Group))) {
  stop(
    "data_path_raw$Group is NA after recoding. Check that Group_code ",
    "contains clean 1/2/3 values before this point in the script."
  )
}
## ---- 1.8 Confounders: functional side ----
confounders_func <- c("Sex", "BMI", "Smoking", "wakefulness", "Socialized")
missing_confounders_func <- setdiff(confounders_func, colnames(data_path_raw))
if (length(missing_confounders_func) > 0) {
  stop(paste("Missing functional confounders:", paste(missing_confounders_func, collapse = ", ")))
}
## ---- 1.9 data_taxa_raw: wide taxa table for pathway-taxa correlations ----
taxa_t <- as.data.frame(t(abund_mat))
taxa_t$SampleID <- rownames(taxa_t)
data_taxa_raw <- merge(meta, taxa_t, by = "SampleID")
data_taxa_raw$Group <- factor(
  as.numeric(stringr::str_trim(as.character(data_taxa_raw$Group_code))),
  levels = c(1, 2, 3),
  labels = c("Healthy", "Resilient", "Unhealthy")
)
data_taxa_raw$Sex <- factor(data_taxa_raw$Sex, levels = c(1, 2), labels = c("Female", "Male"))
data_taxa_raw$Smoking     <- factor(data_taxa_raw$Smoking)
data_taxa_raw$wakefulness <- factor(data_taxa_raw$wakefulness)
data_taxa_raw$Socialized  <- factor(data_taxa_raw$Socialized)
data_taxa_raw$BMI         <- as.numeric(data_taxa_raw$BMI)
cat("\n[Setup] data_taxa_raw created -- samples:", nrow(data_taxa_raw),
    " columns:", ncol(data_taxa_raw), "\n")
## ---- 1.10 Taxonomy helper functions ----
extract_species <- function(x) {
  sapply(strsplit(x, "\\|"), function(z) {
    sp <- z[grepl("^s__", z)]
    if (length(sp) == 0) NA_character_ else sp[1]
  })
}
clean_species_name <- function(x) { x <- sub(".*s__", "", x); gsub("_", " ", x) }
extract_phylum <- function(x) {
  sapply(strsplit(x, "\\|"), function(z) {
    ph <- z[grepl("^p__", z)]
    if (length(ph) == 0) NA_character_ else ph[1]
  })
}
clean_phylum_name <- function(x) { x <- sub(".*p__", "", x); gsub("_", " ", x) }
## ---- 1.11 Spearman correlation CI helper (Fisher z-transformation) ----
# Standard approach for an approximate 95% CI on a Spearman rho: Fisher
# z-transform rho, build a normal-theory CI on the z scale using
# n - 3 as the effective sample size, then back-transform. Used for
# both the pairwise partial Spearman analysis (Section 4) and the
# pathway-taxa Spearman correlations (Section 6).
spearman_ci <- function(rho, n, conf_level = 0.95) {
  if (is.na(rho) || is.na(n) || n < 4) return(c(lower = NA_real_, upper = NA_real_))
  z <- atanh(rho)  # Fisher z-transform: 0.5 * log((1+rho)/(1-rho))
  se_z <- 1 / sqrt(n - 3)
  z_crit <- qnorm(1 - (1 - conf_level) / 2)
  lower <- tanh(z - z_crit * se_z)
  upper <- tanh(z + z_crit * se_z)
  c(lower = lower, upper = upper)
}
cat("\n[Setup] Done.\n")
################################################################################
# 2. FIGURE 1 — Alpha diversity + taxonomic beta diversity (PCoA)
#
# FONT CONTROL: All text uses `fig_font` (set in 2.0). Use the SAME value
#               across Figure 1, Figure 2, and every figure for one typeface.
################################################################################
## ---- 2.0 Global font setting ----
# One place to control the typeface for the whole figure.
#   ""          -> system default sans (renders italic reliably; no showtext)
#   "MyFont"    -> a custom family you register via showtext below
#
# Keep this value identical to Figure 2 so all figures match.
fig_font <- ""
# Optional custom-font registration — ONLY runs if fig_font is not "".
# If you switch to a named font, register it here WITH an italic variant so
# figures that use italic labels (e.g. Figure 2 taxa names) stay italic.
# Registering only a regular variant is exactly what silently kills italics.
if (nzchar(fig_font)) {
  library(showtext)
  library(sysfonts)
  font_add(
    fig_font,
    regular = "YourFont-Regular.ttf",   # <- edit path
    italic  = "YourFont-Italic.ttf"     # <- edit path (keeps italics working)
  )
  showtext_auto()
  showtext_opts(dpi = 300)   # CRITICAL — must match ggsave dpi
}
## ---- 2.0a Journal-grade text theme (all black, consistent sizes) ----
theme_journal <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      text             = element_text(family = fig_font, color = "black"),
      axis.text        = element_text(family = fig_font, color = "black",
                                      size   = base_size),
      axis.title       = element_text(family = fig_font, color = "black",
                                      size   = base_size + 2, face = "bold"),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      axis.ticks       = element_line(color = "black", linewidth = 0.5),
      strip.text       = element_text(family = fig_font, color = "black",
                                      size   = base_size + 2, face = "bold"),
      strip.background = element_rect(fill = "grey92", color = NA),
      legend.text      = element_text(family = fig_font, color = "black",
                                      size   = base_size),
      legend.title     = element_text(family = fig_font, color = "black",
                                      size   = base_size + 2, face = "bold"),
      legend.key.size  = unit(0.6, "cm"),
      panel.border     = element_rect(color = "black", fill = NA,
                                      linewidth = 0.6),
      plot.title       = element_text(family = fig_font, color = "black",
                                      size   = base_size + 4, face = "bold",
                                      hjust  = 0),
      plot.caption     = element_text(family = fig_font, color = "black",
                                      size   = base_size - 2, hjust = 0),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
## ---- 2.1 Alpha diversity indices ----
shannon  <- vegan::diversity(otu_table, index = "shannon")
simpson  <- vegan::diversity(otu_table, index = "simpson")
richness <- vegan::specnumber(otu_table)
evenness <- ifelse(richness > 1, shannon / log(richness), NA)
alpha_df <- data.frame(
  SampleID = rownames(otu_table),
  Shannon  = shannon,
  Simpson  = simpson,
  Richness = richness,
  Evenness = evenness
)
alpha_df <- alpha_df %>%
  left_join(meta[, c("SampleID", "Group", confounders_taxa)], by = "SampleID")
alpha_df$Group <- factor(alpha_df$Group,
                         levels = c("Healthy", "Resilient", "Unhealthy"))
for (v in categorical_vars_taxa) alpha_df[[v]] <- factor(alpha_df[[v]])
for (v in numeric_vars_taxa)     alpha_df[[v]] <- as.numeric(alpha_df[[v]])
## ---- 2.2 Confounder-adjusted alpha models (emmeans pairwise, BH-corrected) ----
metrics                 <- c("Shannon", "Simpson", "Richness", "Evenness")
confounder_formula_taxa <- paste(confounders_taxa, collapse = " + ")
adj_list <- list()
for (metric in metrics) {
  vars_needed <- c(metric, "Group", confounders_taxa)
  dat <- alpha_df[, vars_needed, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat <- droplevels(dat)
  if (nrow(dat) < 5 || nlevels(dat$Group) < 2) next

  factor_vars <- c("Group", categorical_vars_taxa)
  bad_factors <- factor_vars[
    sapply(factor_vars, function(v) nlevels(dat[[v]]) < 2)
  ]
  if (length(bad_factors) > 0) next

  fit <- lm(
    as.formula(paste(metric, "~ Group +", confounder_formula_taxa)),
    data = dat
  )
  em  <- emmeans(fit, specs = "Group", nuisance = confounders_taxa)

  # pairs() with adjust = "BH" returns BH-corrected p-values (= q-values).
  # infer = c(TRUE, TRUE) adds 95% CI columns (lower.CL / upper.CL).
  tmp <- as.data.frame(pairs(em, adjust = "BH", infer = c(TRUE, TRUE)))
  tmp <- tmp %>% dplyr::rename(q_value = p.value)
  tmp$Metric <- metric
  adj_list[[metric]] <- tmp
}
adj_results <- bind_rows(adj_list)
## ---- 2.3 Significant results for bracket annotations ----
if (nrow(adj_results) > 0) {
  significant_results <- adj_results %>%
    tidyr::separate(contrast,
                    into = c("group1", "group2"),
                    sep  = " - ") %>%
    filter(q_value <= 0.05) %>%
    mutate(
      # Bracket label shows BH-corrected q-value explicitly
      p.label = ifelse(
        q_value < 0.001,
        "q < 0.001",
        paste0("q = ", sprintf("%.3f", q_value))
      ),
      sig = case_when(
        q_value < 0.001 ~ "***",
        q_value < 0.01  ~ "**",
        q_value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )
} else {
  significant_results <- data.frame()
}
metric_max <- alpha_df %>%
  summarise(
    Shannon  = max(Shannon,  na.rm = TRUE),
    Simpson  = max(Simpson,  na.rm = TRUE),
    Richness = max(Richness, na.rm = TRUE),
    Evenness = max(Evenness, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "max_value")
if (nrow(significant_results) > 0) {
  significant_results <- significant_results %>%
    left_join(metric_max, by = "Metric")
}
## ---- 2.4 Alpha plot function ----
plot_alpha <- function(data, metric, y_label, sig_df,
                       show_legend = FALSE) {
  df_sig <- sig_df %>% filter(Metric == metric)

  if (nrow(df_sig) > 0) {
    df_sig <- df_sig %>%
      arrange(desc(q_value)) %>%
      mutate(
        x1        = as.numeric(factor(group1, levels = levels(data$Group))),
        x2        = as.numeric(factor(group2, levels = levels(data$Group))),
        row_id    = row_number(),
        y_bracket = max_value * (1.05 + 0.11 * (row_id - 1)),
        y_pvalue  = y_bracket + 0.03 * max_value,
        y_star    = y_bracket - 0.02 * max_value
      )
  }

  data$metric_label <- metric

  p <- ggplot(data, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(
      alpha         = 0.7,
      color         = "black",
      linewidth     = 0.4,
      outlier.shape = NA
    ) +
    geom_point(
      aes(color = Group),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
      size     = 1.8,
      alpha    = 0.75
    )

  if (nrow(df_sig) > 0) {
    p <- p +
      geom_segment(
        data = df_sig,
        aes(x = x1, xend = x2,
            y = y_bracket, yend = y_bracket),
        inherit.aes = FALSE,
        linewidth   = 0.4,
        color       = "black"
      ) +
      geom_segment(
        data = df_sig,
        aes(x = x1, xend = x1,
            y = y_bracket, yend = y_bracket - 0.015 * y_bracket),
        inherit.aes = FALSE,
        linewidth   = 0.4,
        color       = "black"
      ) +
      geom_segment(
        data = df_sig,
        aes(x = x2, xend = x2,
            y = y_bracket, yend = y_bracket - 0.015 * y_bracket),
        inherit.aes = FALSE,
        linewidth   = 0.4,
        color       = "black"
      ) +
      # BH-corrected q-value label on bracket
      geom_text(
        data        = df_sig,
        aes(x = (x1 + x2) / 2, y = y_pvalue, label = p.label),
        inherit.aes = FALSE,
        size        = 5,
        family      = fig_font,
        color       = "black"
      ) +
      geom_text(
        data        = df_sig,
        aes(x = (x1 + x2) / 2, y = y_star, label = sig),
        inherit.aes = FALSE,
        size        = 6,
        family      = fig_font,
        color       = "black"
      )
  }

  p +
    facet_grid(~ metric_label) +
    scale_fill_manual(values  = group_colors) +
    scale_color_manual(values = group_colors) +
    labs(y = y_label, x = NULL) +
    theme_journal(base_size = 14) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.title.x    = element_blank(),
      axis.title.y    = if (nzchar(y_label)) {
        element_text(face = "bold", family = fig_font, color = "black", size = 15)
      } else {
        element_blank()
      },
      legend.position = ifelse(show_legend, "right", "none")
    )
}
## ---- 2.4a Main Figure 1 Panel A: Shannon + Richness only ----
p_shannon  <- plot_alpha(
  data        = alpha_df,
  metric      = "Shannon",
  y_label     = "Diversity value",
  sig_df      = significant_results,
  show_legend = FALSE
)
p_richness <- plot_alpha(
  data        = alpha_df,
  metric      = "Richness",
  y_label     = "",          # no y-axis label — shared label on Shannon panel
  sig_df      = significant_results,
  show_legend = TRUE
)
figureA_main <- (p_shannon | p_richness)
## ---- 2.4b Supplementary: Simpson + Evenness (same style) ----
p_simpson  <- plot_alpha(
  data        = alpha_df,
  metric      = "Simpson",
  y_label     = "Diversity value",
  sig_df      = significant_results,
  show_legend = FALSE
)
p_evenness <- plot_alpha(
  data        = alpha_df,
  metric      = "Evenness",
  y_label     = "",          # no y-axis label — shared label on Simpson panel
  sig_df      = significant_results,
  show_legend = TRUE
)
figureA_supp <- (p_simpson | p_evenness) +
  plot_annotation(
    title   = "Supplementary Figure S1. Alpha diversity — Simpson index and Pielou's evenness",
    caption = paste0(
      "Confounder-adjusted linear models with BMI, socialization, smoking, ",
      "wakefulness, and sex as covariates.\n",
      "Pairwise contrasts via emmeans with Benjamini-Hochberg (BH) correction. ",
      "Brackets show BH-corrected q-values. ",
      "Only significant contrasts (q < 0.05) labelled."
    ),
    theme = theme(
      text         = element_text(family = fig_font, color = "black"),
      plot.title   = element_text(family = fig_font, color = "black",
                                  face = "bold", size = 14),
      plot.caption = element_text(family = fig_font, color = "black",
                                  size = 10, hjust = 0)
    )
  )
print(figureA_supp)
## ---- 2.5 Taxonomic beta diversity + PCoA ----
dist_bc <- vegan::vegdist(otu_table, method = "bray")
pcoa    <- cmdscale(dist_bc, eig = TRUE, k = 2)
df_pcoa <- data.frame(
  SampleID = rownames(otu_table),
  PC1      = pcoa$points[, 1],
  PC2      = pcoa$points[, 2],
  Group    = alpha_df$Group
)
var_PC1 <- round(100 * pcoa$eig[1] / sum(pcoa$eig), 1)
var_PC2 <- round(100 * pcoa$eig[2] / sum(pcoa$eig), 1)
meta_beta <- meta %>%
  filter(SampleID %in% rownames(otu_table)) %>%
  mutate(Group = alpha_df$Group[match(SampleID, alpha_df$SampleID)])
meta_beta <- meta_beta[match(rownames(otu_table), meta_beta$SampleID), ]
for (v in categorical_vars_taxa) meta_beta[[v]] <- factor(meta_beta[[v]])
for (v in numeric_vars_taxa)     meta_beta[[v]] <- as.numeric(meta_beta[[v]])
meta_beta$Group <- factor(meta_beta$Group,
                          levels = c("Healthy", "Resilient", "Unhealthy"))
if (!identical(meta_beta$SampleID, rownames(otu_table)))
  stop("meta_beta sample order does not match otu_table order.")
## ---- 2.6 PERMDISP dispersion check ----
disp              <- vegan::betadisper(dist_bc, meta_beta$Group)
permdisp_test     <- anova(disp)
permdisp_pairwise <- TukeyHSD(disp)
disp_df <- data.frame(
  SampleID             = meta_beta$SampleID,
  Group                = meta_beta$Group,
  Distance_to_centroid = disp$distances
)
distance_kw <- kruskal.test(Distance_to_centroid ~ Group, data = disp_df)
cat("\n[Fig1B] Taxonomic PERMDISP Kruskal-Wallis: chi-sq =",
    round(distance_kw$statistic, 3),
    ", df =", distance_kw$parameter,
    ", p =",  signif(distance_kw$p.value, 3), "\n")
## ---- 2.7 Both PERMANOVAs + pooled FDR q-value (computed once here) ----
# 2.7a Taxonomic PERMANOVA
global_permanova_adj <- vegan::adonis2(
  dist_bc ~ Group + BMI + Socialized + Smoking + wakefulness + Sex,
  data = meta_beta, permutations = 999
)
p_val_adj <- global_permanova_adj$`Pr(>F)`[1]
r2_adj    <- global_permanova_adj$R2[1]
# 2.7b Functional pathway matrix + metadata
# intersect() guards against any pathway name that did not survive the
# merge into data_path_raw -- all_of() would otherwise error on a missing col.
pathway_matrix_func <- data_path_raw %>%
  dplyr::select(all_of(intersect(pathways, colnames(data_path_raw)))) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
  as.data.frame()
rownames(pathway_matrix_func) <- data_path_raw$SampleID
meta_beta_func <- data_path_raw %>%
  dplyr::select(SampleID, Group, Sex, BMI, Smoking, wakefulness, Socialized) %>%
  dplyr::mutate(
    Group       = droplevels(factor(Group)),
    Sex         = droplevels(factor(Sex)),
    Smoking     = droplevels(factor(Smoking)),
    wakefulness = droplevels(factor(wakefulness)),
    Socialized  = droplevels(factor(Socialized)),
    BMI         = as.numeric(BMI)
  )
# 2.7c Functional Bray-Curtis distance
dist_bc_func <- vegan::vegdist(pathway_matrix_func, method = "bray")
# 2.7d Functional PERMANOVA
global_permanova_func <- vegan::adonis2(
  dist_bc_func ~ Group + Sex + BMI + Smoking + wakefulness + Socialized,
  data = meta_beta_func, permutations = 999
)
p_val_func <- global_permanova_func$`Pr(>F)`[1]
r2_func    <- global_permanova_func$R2[1]
# 2.7e Pool both p-values and apply FDR correction once
global_permanova_pvals <- c(
  Taxonomic_beta  = p_val_adj,
  Functional_beta = p_val_func
)
global_permanova_qvals <- p.adjust(global_permanova_pvals, method = "fdr")
cat("\n================ POOLED GLOBAL PERMANOVA Q-VALUES ================\n")
print(data.frame(
  Test    = names(global_permanova_pvals),
  p_value = signif(global_permanova_pvals, 3),
  q_value = signif(global_permanova_qvals, 3)
))
q_val_adj_pooled <- global_permanova_qvals["Taxonomic_beta"]
q_val_func       <- global_permanova_qvals["Functional_beta"]
label_text_taxa <- paste0(
  "q = ", signif(q_val_adj_pooled, 2),
  "\nR² = ", signif(r2_adj, 2)
)
label_text_func <- paste0(
  "q = ", signif(q_val_func, 2),
  "\nR² = ", signif(r2_func, 2)
)
## ---- 2.8 Centroids + Panel B ----
centroids <- df_pcoa %>%
  group_by(Group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
# build_figureB reused for Figure 8 functional PCoA in Section 7
build_figureB <- function(df_points, df_centroids,
                          x_lab, y_lab, label, colors) {
  ggplot(df_points, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2.0, alpha = 0.85) +
    stat_ellipse(type = "norm", level = 0.68, linewidth = 0.6) +
    # Arrows from origin to centroid only — no group name text on plot
    geom_segment(
      data      = df_centroids,
      aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Group),
      arrow     = arrow(length = unit(0.015, "npc"), type = "closed"),
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    # Pooled FDR q-value + R² annotation top-left
    annotate(
      "text",
      x      = min(df_points$PC1, na.rm = TRUE),
      y      = max(df_points$PC2, na.rm = TRUE),
      label  = label,
      hjust  = 0, vjust = 1,
      size   = 5,
      family = fig_font,
      color  = "black"
    ) +
    scale_color_manual(values = colors) +
    labs(x = x_lab, y = y_lab) +
    theme_journal(base_size = 14) +
    theme(legend.position = "right")
}
figureB <- build_figureB(
  df_points    = df_pcoa,
  df_centroids = centroids,
  x_lab        = paste0("PC1 (", var_PC1, "%)"),
  y_lab        = paste0("PC2 (", var_PC2, "%)"),
  label        = label_text_taxa,
  colors       = group_colors
)
## ---- 2.9 Assemble Figure 1 (manual tags for consistent size) ----
figureA_final <- ggdraw(figureA_main) +
  draw_label("A", x = 0.01, y = 0.99,
             hjust      = 0,
             vjust      = 1,
             fontfamily = fig_font,
             fontface   = "bold",
             color      = "black",
             size       = 20)
figureB_final <- ggdraw(figureB) +
  draw_label("B", x = 0.01, y = 0.99,
             hjust      = 0,
             vjust      = 1,
             fontfamily = fig_font,
             fontface   = "bold",
             color      = "black",
             size       = 20)
Figure1 <- plot_grid(
  figureA_final,
  figureB_final,
  ncol        = 1,
  rel_heights = c(1, 1.3)
)
print(Figure1)
## ---- 2.10 Export Figure 1 ----
# ggsave(
#   filename    = "~/Desktop/Figure1.tiff",
#   plot        = Figure1,
#   device      = "tiff",
#   dpi         = 300,
#   width       = 260,          # was 220 mm
#   height      = 300,          # was 260 mm
#   units       = "mm",
#   compression = "lzw"
# )
#
# ## ---- 2.11 Export Supplementary Figure S1 ----
# ggsave(
#   filename    = "~/Desktop/FigureS1_Simpson_Evenness.tiff",
#   plot        = figureA_supp,
#   device      = "tiff",
#   dpi         = 300,
#   width       = 260,          # was 220 mm
#   height      = 150,          # was 130 mm
#   units       = "mm",
#   compression = "lzw"
# )
#
# cat("\n[Section 2 complete]\n")
# cat("  Figure1  saved -> ~/Desktop/Figure1.tiff\n")
# cat("  FigureS1 saved -> ~/Desktop/FigureS1_Simpson_Evenness.tiff\n")
# cat("  label_text_taxa and label_text_func ready for Sections 7 and 8\n")
## ---- 2.12 Pairwise adjusted PERMANOVA + FDR (taxonomic) ----
pairwise_permanova_adj <- function(distmat, metadata, group_var, confounders) {
  groups   <- levels(metadata[[group_var]])
  dist_mat <- as.matrix(distmat)
  results  <- data.frame()
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      g1 <- groups[i]; g2 <- groups[j]
      subset_meta <- metadata %>%
        filter(.data[[group_var]] %in% c(g1, g2))
      subset_meta[[group_var]] <- droplevels(subset_meta[[group_var]])
      subset_ids  <- intersect(subset_meta$SampleID, rownames(dist_mat))
      subset_meta <- subset_meta[match(subset_ids, subset_meta$SampleID), ]
      sub_dist    <- as.dist(dist_mat[subset_ids, subset_ids])
      fml <- as.formula(
        paste("sub_dist ~", group_var, "+",
              paste(confounders, collapse = " + "))
      )
      ad <- vegan::adonis2(fml, data = subset_meta, permutations = 999)
      results <- rbind(results, data.frame(
        Group1  = g1,       Group2  = g2,
        F_model = ad$F[1],  R2      = ad$R2[1],
        p_value = ad$`Pr(>F)`[1]
      ))
    }
  }
  results$Q_value <- p.adjust(results$p_value, method = "fdr")
  results
}
pairwise_results_adj <- pairwise_permanova_adj(
  dist_bc, meta_beta, "Group", confounders_taxa
)
#
################################################################################
#3. FIGURE 2
################################################################################
## ---- 3.1 Required libraries ----
library(tidyverse)
library(cowplot)
library(grid)
## ---- 3.2 Long abundance table ----
abund_long <- as.data.frame(abund_mat) %>%
  rownames_to_column("Taxon") %>%
  pivot_longer(
    cols      = -Taxon,
    names_to  = "SampleID",
    values_to = "Abundance"
  ) %>%
  left_join(meta2, by = "SampleID") %>%
  filter(!is.na(Group))
abund_long$Group <- factor(
  abund_long$Group,
  levels = c("Healthy", "Resilient", "Unhealthy")
)
## ---- 3.3 Safe adjusted categorical model function ----
run_group_model <- function(dat, taxon_name_col = NULL, confounders) {

  keep_cols <- intersect(
    c("Abundance", "Group", confounders),
    colnames(dat)
  )

  dat <- dat[, keep_cols, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat <- droplevels(dat)

  if (!"Group" %in% colnames(dat)) {
    return(data.frame(
      p_group                 = NA_real_,
      fc_unhealthy_vs_healthy = NA_real_
    ))
  }

  dat$Group <- factor(
    dat$Group,
    levels = c("Healthy", "Resilient", "Unhealthy")
  )
  dat$Group <- droplevels(dat$Group)

  valid_confounders <- c()

  for (v in confounders) {

    if (!v %in% colnames(dat)) next

    x <- dat[[v]]

    if (is.factor(x) || is.character(x) || is.logical(x)) {

      x <- droplevels(factor(x))

      if (nlevels(x) >= 2) {
        dat[[v]]          <- x
        valid_confounders <- c(valid_confounders, v)
      }

    } else if (is.numeric(x)) {

      if (length(unique(x[is.finite(x)])) >= 2) {
        valid_confounders <- c(valid_confounders, v)
      }
    }
  }

  if (nrow(dat) < 5 || nlevels(dat$Group) < 2) {
    return(data.frame(
      p_group                 = NA_real_,
      fc_unhealthy_vs_healthy = NA_real_
    ))
  }

  rhs <- paste(c("Group", valid_confounders), collapse = " + ")

  fit <- tryCatch(
    lm(as.formula(paste("Abundance ~", rhs)), data = dat),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(data.frame(
      p_group                 = NA_real_,
      fc_unhealthy_vs_healthy = NA_real_
    ))
  }

  a <- anova(fit)

  m_healthy   <- mean(dat$Abundance[dat$Group == "Healthy"],   na.rm = TRUE)
  m_unhealthy <- mean(dat$Abundance[dat$Group == "Unhealthy"], na.rm = TRUE)

  data.frame(
    p_group = if ("Group" %in% rownames(a)) {
      a["Group", "Pr(>F)"]
    } else {
      NA_real_
    },
    fc_unhealthy_vs_healthy = ifelse(
      is.na(m_healthy) || m_healthy == 0,
      NA_real_,
      m_unhealthy / m_healthy
    )
  )
}
## ---- 3.4 Taxa name cleaning ----
# Returns plain-text display names (no plotmath, no HTML).
# Italic is applied via theme(strip.text = element_text(face = "italic"))
# and theme(legend.text = element_text(face = "italic")).
taxa_font_size <- 10
taxa_geom_size <- taxa_font_size * 0.35
clean_taxon_display <- function(x) {

  x <- as.character(x)

  x <- gsub("</?i>", "", x)        # strip any stray HTML italic tags
  x <- gsub(".*\\|s__", "", x)     # remove clade prefixes up to s__
  x <- gsub(".*\\|p__", "", x)
  x <- gsub("^s__", "", x)
  x <- gsub("^p__", "", x)
  x <- gsub("^g__", "", x)
  x <- gsub("_", " ", x)           # underscores -> spaces
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)

  x[x == ""] <- NA_character_

  x
}
## ---- 3.5 Theme ----
# strip.text uses face = "italic" — this is the reliable italic method.
# legend.text also uses face = "italic" for taxa legends (panels A & C).
# Group labels stay plain (overridden per-panel where needed).
theme_fig2 <- function(base_size = 16) {

  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = fig_font, color = "black"),
      axis.text = element_text(
        family = fig_font, color = "black", size = base_size
      ),
      axis.title = element_text(
        family = fig_font, color = "black", size = base_size, face = "bold"
      ),
      axis.line  = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      panel.border      = element_blank(),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA),
      strip.background  = element_rect(fill = "grey92", color = NA),
      # ── face = "italic" WITHOUT forcing a family ──
      # A family whose italic variant isn't registered silently kills the
      # italic. Omitting family here lets the default font's italic render.
      strip.text = element_text(
        face  = "italic",
        color = "black",
        size  = taxa_font_size
      ),
      legend.text = element_text(
        face  = "italic",   # italic for taxa legend labels
        color = "black",
        size  = taxa_font_size
      ),
      legend.title    = element_blank(),
      legend.key.size = unit(0.32, "cm"),
      legend.spacing.y = unit(0.08, "cm"),
      legend.background = element_rect(fill = "white", color = NA, linewidth = 0),
      plot.caption = element_text(
        family = fig_font, color = "black", size = base_size - 2, hjust = 0
      )
    )
}
# Italic guide for taxa legends (panels A & C) — no family forced
italic_guide <- function() {
  guide_legend(
    ncol      = 1,
    keyheight = unit(0.32, "cm"),
    keywidth  = unit(0.32, "cm"),
    label.theme = element_text(
      face = "italic", color = "black", size = taxa_font_size
    )
  )
}
# Plain guide for group legends (panels B & D)
plain_guide <- function() {
  guide_legend(
    ncol      = 1,
    keyheight = unit(0.4, "cm"),
    keywidth  = unit(0.4, "cm"),
    label.theme = element_text(
      family = fig_font, face = "plain", color = "black", size = 13
    )
  )
}
################################################################################
# 3.6 PHYLUM-LEVEL DATA + MODELS
################################################################################
phylum_long <- abund_long %>%
  mutate(
    Phylum_raw = extract_phylum(Taxon),
    Phylum     = clean_taxon_display(clean_phylum_name(Phylum_raw))
  ) %>%
  filter(!is.na(Phylum), Phylum != "")
phylum_pct <- phylum_long %>%
  group_by(SampleID, Group, Phylum) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(
    total_abund = sum(Abundance, na.rm = TRUE),
    Abundance = ifelse(total_abund > 0, Abundance / total_abund * 100, NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::select(-total_abund) %>%
  left_join(meta2, by = c("SampleID", "Group"))
top_phyla <- phylum_pct %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  pull(Phylum)
phylum_top <- phylum_pct %>%
  filter(Phylum %in% top_phyla) %>%
  mutate(Group = factor(Group, levels = c("Healthy", "Resilient", "Unhealthy")))
lm_results_p <- phylum_top %>%
  group_by(Phylum) %>%
  group_modify(~ run_group_model(.x, "Phylum", confounders_taxa)) %>%
  ungroup() %>%
  mutate(
    qval  = p.adjust(p_group, method = "fdr"),
    label = paste0(
      "\nFC = ", round(fc_unhealthy_vs_healthy, 2),
      "\nFDR q = ", formatC(qval, digits = 3)
    )
  )
sig_phyla_filtered <- lm_results_p %>%
  left_join(
    phylum_top %>%
      group_by(Phylum) %>%
      summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop"),
    by = "Phylum"
  ) %>%
  filter(qval < 0.05, mean_abund >= 5)
phylum_sig <- phylum_top %>%
  filter(Phylum %in% sig_phyla_filtered$Phylum)
## Order significant phyla by mean abundance in Unhealthy
phylum_unhealthy_order <- phylum_sig %>%
  filter(Group == "Unhealthy") %>%
  group_by(Phylum) %>%
  summarise(mean_unhealthy = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_unhealthy)) %>%
  pull(Phylum)
phylum_ymax <- phylum_sig %>%
  group_by(Phylum) %>%
  summarise(
    ymax = ifelse(all(is.na(Abundance)), NA_real_, max(Abundance, na.rm = TRUE)),
    .groups = "drop"
  )
sig_phyla_ann <- sig_phyla_filtered %>%
  left_join(phylum_ymax, by = "Phylum")
# Plain factor levels — italic comes from theme, not plotmath
phylum_sig <- phylum_sig %>%
  mutate(
    Phylum_plot = factor(Phylum, levels = phylum_unhealthy_order)
  )
sig_phyla_ann <- sig_phyla_ann %>%
  mutate(
    Phylum_plot = factor(Phylum, levels = phylum_unhealthy_order)
  )
################################################################################
# 3.7 PANEL A — PHYLUM STACKED BAR
################################################################################
df_bar_p <- phylum_top %>%
  group_by(Group, Phylum) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE) * 100) %>%
  ungroup()
phylum_order <- df_bar_p %>%
  group_by(Phylum) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  pull(Phylum)
# Plain factor levels for Panel A legend
df_bar_p <- df_bar_p %>%
  mutate(
    Phylum_plot = factor(Phylum, levels = rev(phylum_order))
  )
soft_pastel_phyla <- c(
  "#6A5B80", "#606060", "#E3D6E9", "#F2D3B3", "#FFF1B8",
  "#E6D97A", "#F29472", "#C6A3B3", "#9AB3C5", "#B7C9A8",
  "#FFB347", "#C77FA5", "#B3A3CB", "#51687D", "#7A3E63",
  "#B5EAD7", "#6D8A55", "#3E503C", "#3E5C76", "#D76A49"
)
phylum_colors <- setNames(
  rep(soft_pastel_phyla, length.out = length(phylum_order)),
  rev(phylum_order)   # plain names now (no plotmath)
)
pA <- ggplot(
  df_bar_p,
  aes(x = Group, y = Abundance, fill = Phylum_plot)
) +
  geom_bar(stat = "identity", width = 0.65, color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = phylum_colors,
    guide  = italic_guide()   # italic applied via label.theme in guide
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_fig2(base_size = 16) +
  theme(
    axis.text.x = element_text(
      face = "bold", size = 14, family = fig_font,
      color = "black", angle = 0, hjust = 1
    ),
    legend.position   = "right",
    legend.background = element_rect(fill = "white", color = NA, linewidth = 0),
    legend.margin     = margin(0, 0, 0, 4),
    plot.margin       = margin(5, 5, 5, 5),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )
################################################################################
# 3.8 PANEL B — FDR-SIGNIFICANT PHYLA
################################################################################
if (nrow(phylum_sig) == 0) {

  pB <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label  = "No phyla passed FDR q < 0.05\nand mean abundance ≥ 5%",
      size   = taxa_geom_size, family = fig_font, color = "black"
    ) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

} else {

  pB <- ggplot(phylum_sig, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(
      alpha = 0.6, outlier.shape = NA, color = "black",
      linewidth = 0.4, width = 0.6
    ) +
    geom_jitter(aes(color = Group), width = 0.2, size = 1.8, alpha = 0.7) +
    # ── No labeller needed; strip.text face = "italic" handles italics ──
    facet_wrap(~ Phylum_plot, scales = "free_y", ncol = 1) +
    geom_text(
      data = sig_phyla_ann,
      aes(
        x     = 1,
        y     = Inf,
        label = paste0(
          "q = ", formatC(qval, digits = 3),
          "   FC = ", round(fc_unhealthy_vs_healthy, 2)
        )
      ),
      inherit.aes = FALSE,
      hjust  = 0, vjust = 1.4,
      size   = taxa_geom_size, family = fig_font, color = "black"
    ) +
    scale_fill_manual(values  = group_colors, guide = plain_guide()) +
    scale_color_manual(values = group_colors, guide = "none") +
    labs(x = NULL, y = "Relative abundance (%)") +
    theme_fig2(base_size = 16) +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      # Override legend.text to plain for group legend in this panel
      legend.text      = element_text(
        family = fig_font, face = "plain", color = "black", size = 13
      ),
      legend.position  = "none",
      strip.background = element_rect(fill = "grey92", color = NA),
      panel.spacing    = unit(0.8, "lines"),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
################################################################################
# 3.9 SPECIES-LEVEL DATA + MODELS
################################################################################
species_long <- abund_long %>%
  mutate(
    Species_raw = extract_species(Taxon),
    Species     = clean_taxon_display(Species_raw)
  ) %>%
  filter(!is.na(Species), Species != "")
species_pct <- species_long %>%
  group_by(SampleID, Group, Species) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(SampleID) %>%
  mutate(
    total_abund = sum(Abundance, na.rm = TRUE),
    Abundance = ifelse(total_abund > 0, Abundance / total_abund * 100, NA_real_)
  ) %>%
  ungroup() %>%
  dplyr::select(-total_abund) %>%
  left_join(meta2, by = c("SampleID", "Group"))
top30_species <- species_pct %>%
  group_by(Species) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 30) %>%
  pull(Species)
species_top30 <- species_pct %>%
  filter(Species %in% top30_species) %>%
  mutate(Group = factor(Group, levels = c("Healthy", "Resilient", "Unhealthy")))
lm_results_s <- species_top30 %>%
  group_by(Species) %>%
  group_modify(~ run_group_model(.x, "Species", confounders_taxa)) %>%
  ungroup() %>%
  mutate(
    qval  = p.adjust(p_group, method = "fdr"),
    label = paste0(
      "\nFC = ", round(fc_unhealthy_vs_healthy, 2),
      "\nFDR q = ", formatC(qval, digits = 3)
    )
  )
sig_species_omnibus <- lm_results_s %>% filter(qval < 0.05)
species_sig <- species_top30 %>%
  filter(Species %in% sig_species_omnibus$Species)
## Order significant species by mean abundance in Unhealthy
species_unhealthy_order <- species_sig %>%
  filter(Group == "Unhealthy") %>%
  group_by(Species) %>%
  summarise(mean_unhealthy = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_unhealthy)) %>%
  pull(Species)
species_ymax <- species_sig %>%
  group_by(Species) %>%
  summarise(
    ymax = ifelse(all(is.na(Abundance)), NA_real_, max(Abundance, na.rm = TRUE)),
    .groups = "drop"
  )
sig_species_ann <- sig_species_omnibus %>%
  left_join(species_ymax, by = "Species")
# Plain factor levels — italic from theme
species_sig <- species_sig %>%
  mutate(
    Species_plot = factor(Species, levels = species_unhealthy_order)
  )
sig_species_ann <- sig_species_ann %>%
  mutate(
    Species_plot = factor(Species, levels = species_unhealthy_order)
  )
################################################################################
# 3.10 PANEL C — TOP 30 SPECIES STACKED BAR
################################################################################
df_bar_s <- species_top30 %>%
  group_by(Group, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Abundance = Abundance / sum(Abundance, na.rm = TRUE) * 100) %>%
  ungroup()
species_order <- df_bar_s %>%
  group_by(Species) %>%
  summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_abund) %>%
  pull(Species)
# Plain factor levels for Panel C legend
df_bar_s <- df_bar_s %>%
  mutate(
    Species_plot = factor(Species, levels = species_order)
  )
soft_pastel_species <- c(
  "#606060", "#B5EAD7", "#E3D6E9", "#F2D3B3", "#FFF1B8",
  "#E6D97A", "#F29472", "#C6A3B3", "#9AB3C5", "#B7C9A8",
  "#FFB5A7", "#B28DFF", "#D5AAFF", "#AFCBFF", "#7FB285",
  "#D4A5A5", "#AAC0AA", "#97C1A9", "#FFD6A5", "#BCBABA",
  "#FFB347", "#C77FA5", "#B3A3CB", "#51687D", "#6A5B80",
  "#7A3E63", "#D76A49", "#6D8A55", "#3E503C", "#3E5C76"
)
species_colors <- setNames(
  rep(soft_pastel_species, length.out = length(species_order)),
  species_order   # plain names now
)
pC <- ggplot(
  df_bar_s,
  aes(x = Group, y = Abundance, fill = Species_plot)
) +
  geom_bar(stat = "identity", width = 0.65, color = "white", linewidth = 0.2) +
  scale_fill_manual(
    values = species_colors,
    guide  = italic_guide()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_fig2(base_size = 16) +
  theme(
    axis.text.x = element_text(
      face = "bold", size = 14, family = fig_font,
      color = "black", angle = 0, hjust = 1
    ),
    legend.position   = "right",
    legend.background = element_rect(fill = "white", color = NA, linewidth = 0),
    legend.margin     = margin(0, 0, 0, 4),
    plot.margin       = margin(5, 5, 5, 5),
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA)
  )
################################################################################
# 3.11 PANEL D — FDR-SIGNIFICANT SPECIES
################################################################################
if (nrow(species_sig) == 0) {

  pD <- ggplot() +
    annotate(
      "text", x = 0, y = 0,
      label  = "No top-30 species passed FDR q < 0.05",
      size   = taxa_geom_size, family = fig_font, color = "black"
    ) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

} else {

  pD <- ggplot(species_sig, aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(
      alpha = 0.6, outlier.shape = NA, color = "black",
      linewidth = 0.4, width = 0.6
    ) +
    geom_jitter(aes(color = Group), width = 0.2, size = 1.8, alpha = 0.7) +
    # ── No labeller needed; strip.text face = "italic" handles italics ──
    facet_wrap(~ Species_plot, scales = "free_y", ncol = 3) +
    geom_text(
      data = sig_species_ann,
      aes(
        x     = 1,
        y     = Inf,
        label = paste0(
          "q = ", formatC(qval, digits = 3),
          "   FC = ", round(fc_unhealthy_vs_healthy, 2)
        )
      ),
      inherit.aes = FALSE,
      hjust  = 0, vjust = 1.4,
      size   = taxa_geom_size, family = fig_font, color = "black"
    ) +
    scale_fill_manual(values  = group_colors, guide = plain_guide()) +
    scale_color_manual(values = group_colors, guide = "none") +
    labs(x = NULL, y = "Relative abundance (%)") +
    theme_fig2(base_size = 16) +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      # Override legend.text to plain for group legend
      legend.text      = element_text(
        family = fig_font, face = "plain", color = "black", size = 13
      ),
      legend.position  = "right",
      strip.background = element_rect(fill = "grey92", color = NA),
      panel.spacing    = unit(0.8, "lines"),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
################################################################################
# 3.12 ASSEMBLE FIGURE 2
################################################################################
legend_A <- cowplot::get_legend(pA)
legend_C <- cowplot::get_legend(pC)
# Standalone group legend (plain text, not italic)
group_legend_df <- data.frame(
  Group = factor(
    c("Healthy", "Resilient", "Unhealthy"),
    levels = c("Healthy", "Resilient", "Unhealthy")
  ),
  x = 1, y = 1
)
group_legend_plot <- ggplot(
  group_legend_df,
  aes(x = x, y = y, fill = Group)
) +
  geom_point(shape = 21, size = 4, color = "black", stroke = 0.3) +
  scale_fill_manual(values = group_colors, guide = plain_guide()) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.text       = element_text(
      family = fig_font, face = "plain", color = "black", size = 10
    ),
    legend.title      = element_blank(),
    legend.background = element_rect(fill = "white", color = NA)
  )
legend_BD <- cowplot::get_legend(group_legend_plot)
pA_noleg <- pA + theme(legend.position = "none")
pC_noleg <- pC + theme(legend.position = "none")
pB_noleg <- pB + theme(legend.position = "none")
pD_noleg <- pD + theme(legend.position = "none")
pA_draw <- ggdraw(pA_noleg) +
  draw_label(
    "A", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
    fontfamily = fig_font, fontface = "bold", color = "black", size = 20
  )
pB_draw <- ggdraw(pB_noleg) +
  draw_label(
    "B", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
    fontfamily = fig_font, fontface = "bold", color = "black", size = 20
  )
pC_draw <- ggdraw(pC_noleg) +
  draw_label(
    "C", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
    fontfamily = fig_font, fontface = "bold", color = "black", size = 20
  )
pD_draw <- ggdraw(pD_noleg) +
  draw_label(
    "D", x = 0.01, y = 0.99, hjust = 0, vjust = 1,
    fontfamily = fig_font, fontface = "bold", color = "black", size = 20
  )
# Left column: A over C
left_plots <- plot_grid(
  pA_draw, pC_draw,
  ncol = 1, nrow = 2,
  align = "v", axis = "lr",
  rel_heights = c(1, 1)
)
# Taxa legends column: A legend over C legend
left_legends <- plot_grid(
  legend_A, legend_C,
  ncol = 1, nrow = 2,
  rel_heights = c(1, 1)
)
left_legends_white <- ggdraw(left_legends) +
  theme(plot.background = element_rect(fill = "white", color = NA))
# Right column: B above D
right_plots <- plot_grid(
  pB_draw, pD_draw,
  ncol = 1, nrow = 2,
  align = "v", axis = "lr",
  rel_heights = c(1, 2.5)
)
# Final layout: | bar plots | taxa legends | B above D | group legend |
Figure2 <- plot_grid(
  left_plots,
  left_legends_white,
  right_plots,
  legend_BD,
  ncol       = 4,
  nrow       = 1,
  rel_widths = c(1.2, 0.8, 1.6, 0.25)
)
print(Figure2)
################################################################################
# 3.13 EXPORT FIGURE 2
################################################################################
#
# ggsave(
#   filename    = "~/Desktop/Figure2.tiff",
#   plot        = Figure2,
#   device      = "tiff",
#   dpi         = 300,
#   width       = 560,
#   height      = 300,
#   units       = "mm",
#   compression = "lzw",
#   bg          = "white"
# )
#
# cat("\n[Section 3 complete]\n")
# cat("  Figure2 saved -> ~/Desktop/Figure2.tiff\n")
## ---- 3.12 Heatmap of omnibus FDR-significant species ----
significant_species_from_fig2 <- sig_species_omnibus %>%
  pull(Species) %>%
  unique()
# Strip HTML tags for downstream use
significant_species_from_fig2 <- gsub("<.*?>", "",
                                      significant_species_from_fig2)
significant_species_from_fig2 <- significant_species_from_fig2[
  !is.na(significant_species_from_fig2)
]
cat("\n[Fig2] Number of omnibus-significant species:",
    length(significant_species_from_fig2), "\n")
if (length(significant_species_from_fig2) == 0) {
  message("No significant species found at qval < 0.05 -- skipping heatmap.")
} else {
  species_heat_df <- species_pct %>%
    filter(Species %in% significant_species_from_fig2) %>%
    mutate(Group = factor(Group,
                          levels = c("Healthy", "Resilient", "Unhealthy")))

  species_group_mean <- species_heat_df %>%
    group_by(Species, Group) %>%
    summarise(mean_abundance = mean(Abundance, na.rm = TRUE),
              .groups = "drop")

  heat_mat <- species_group_mean %>%
    pivot_wider(names_from  = Group,
                values_from = mean_abundance,
                values_fill = 0) %>%
    column_to_rownames("Species") %>%
    as.matrix()
  heat_mat <- heat_mat[, c("Healthy", "Resilient", "Unhealthy"),
                       drop = FALSE]

  heat_mat_z <- t(scale(t(heat_mat)))
  heat_mat_z[is.na(heat_mat_z)] <- 0
  colnames(heat_mat_z) <- c("Healthy", "Resilient", "Unhealthy")

  heat_colors <- colorRampPalette(
    c("#313695", "white", "#D73027")
  )(100)

  pheatmap(
    heat_mat_z,
    color          = heat_colors,
    breaks         = seq(-2, 2, length.out = 101),
    cluster_rows   = TRUE,
    cluster_cols   = FALSE,
    labels_row     = rownames(heat_mat_z),
    fontsize_row   = taxa_font_size,
    fontsize_col   = taxa_font_size,
    fontsize       = taxa_font_size,
    border_color   = NA,
    main           = "FDR-significant species by group (omnibus model)",
    angle_col      = 90,
    treeheight_row = 20
  )
}
################################################################################
# 3B. FIGURE 3 — Core gut microbiome (shared and unique core taxa)
library(VennDiagram)
library(gridExtra)
library(grid)
library(cowplot)
# Belt-and-suspenders: re-assert ggplot2 scales over any masks.
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(
    ggplot2::scale_fill_manual,
    ggplot2::scale_fill_identity,
    ggplot2::scale_color_manual,
    ggplot2::scale_colour_manual,
    .quiet = TRUE
  )
}
## ---- 3B.0 Global font setting ----
fig_font <- ""
if (nzchar(fig_font)) {
  library(showtext)
  library(sysfonts)
  font_add(
    fig_font,
    regular = "YourFont-Regular.ttf",
    italic  = "YourFont-Italic.ttf"
  )
  showtext_auto()
  showtext_opts(dpi = 300)
}
# Font for base-grid graphics (VennDiagram, pheatmap), which reject "".
grid_font <- if (nzchar(fig_font)) fig_font else "sans"
## ---- 3B.1 Parameters ----
min_abund_cutoff <- 0.001   # 0.1% relative abundance
min_prevalence   <- 0.5     # present in >=50% of samples within a group
# Match Figure 2's taxa label size so species labels are consistent across figures
taxa_font_size <- 9         # slightly smaller than Figure 2's 10
# Figure 3 export dimensions (mm) — edit these to resize the whole figure
fig3_width_mm  <- 250
fig3_height_mm <- 95
# Panel label (A/B) point size. Figures 1 & 2 use 20, but those figures are
# much larger; on this shorter figure 20pt looks oversized, so scale it down.
panel_label_size <- 11
## ---- 3B.2 Build sample x species table directly from abund_mat/meta2 ----
species_mat_raw <- abund_mat[grepl("s__", rownames(abund_mat)), , drop = FALSE]
species_names_clean <- clean_species_name(extract_species(rownames(species_mat_raw)))
species_names_clean <- ifelse(
  is.na(species_names_clean) | species_names_clean == "" |
    grepl("uncult|unclass|metagenome|unknown|bacterium|strain|GGB|UBA|CAG|SGB|MAG|sp\\.|genome",
          species_names_clean, ignore.case = TRUE),
  "Unknown",
  species_names_clean
)
species_mat_df <- as.data.frame(t(species_mat_raw), check.names = FALSE)
colnames(species_mat_df) <- species_names_clean
species_collapsed <- as.data.frame(
  do.call(cbind, lapply(split.default(species_mat_df, colnames(species_mat_df)), function(x) {
    x <- as.data.frame(x)
    if (ncol(x) == 1) x[[1]] else rowMeans(x, na.rm = TRUE)
  })),
  check.names = FALSE
)
species_clean <- cbind(
  SampleID = rownames(species_collapsed),
  species_collapsed
)
species_clean <- merge(meta2[, c("SampleID", "Group")], species_clean, by = "SampleID")
species_clean$Group <- factor(species_clean$Group, levels = c("Healthy", "Resilient", "Unhealthy"))
meta_cols_core <- c("SampleID", "Group")
taxa_cols_core <- setdiff(colnames(species_clean), meta_cols_core)
kept_taxa <- taxa_cols_core[
  sapply(species_clean[, taxa_cols_core, drop = FALSE],
         function(x) any(x >= min_abund_cutoff, na.rm = TRUE))
]
species_clean <- species_clean[, c(meta_cols_core, kept_taxa), drop = FALSE]
cat("\n[Fig3] Samples:", nrow(species_clean), " | Species columns kept:", length(kept_taxa), "\n")
## ---- 3B.3 Identify core species per group ----
get_core_taxa <- function(df) {
  taxa_cols <- setdiff(names(df), meta_cols_core)
  taxa_df <- df[, taxa_cols, drop = FALSE]
  present <- taxa_df >= min_abund_cutoff
  prevalence <- colMeans(present, na.rm = TRUE)
  names(prevalence[prevalence >= min_prevalence])
}
core_species <- lapply(split(species_clean, species_clean$Group), get_core_taxa)
shared_species <- Reduce(intersect, core_species)
shared_tbl <- species_clean[, c(meta_cols_core, intersect(shared_species, colnames(species_clean))), drop = FALSE]
cat("[Fig3] Core species per group:", paste(names(core_species), sapply(core_species, length), sep = "=", collapse = ", "), "\n")
cat("[Fig3] Shared core species across all 3 groups:", length(shared_species), "\n")
## ---- 3B.4 Shared core taxa prevalence table (used by SUPPLEMENTARY heatmap) ----
if (length(shared_species) > 0) {
  heatmap_tbl_core <- shared_tbl %>%
    pivot_longer(cols = -all_of(meta_cols_core), names_to = "Taxon", values_to = "Abundance") %>%
    group_by(Taxon, Group) %>%
    summarise(
      prevalence = 100 * mean(Abundance >= min_abund_cutoff, na.rm = TRUE),
      avg_abundance = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  heatmap_tbl_core <- data.frame()
}
################################################################################
# PANEL A : Venn diagram of core-species overlap  -> captured as a grob
#           Cleaner, smaller styling: thin outlines, soft fills, counts only,
#           no per-region clutter. Title handled in assembly, not inside Venn.
################################################################################
group_levels_core <- names(core_species)
venn_colors <- group_colors[group_levels_core]
if (any(is.na(venn_colors))) {
  warning("Some groups missing from group_colors; using default ggplot2 hues as fallback.")
  fallback_n <- sum(is.na(venn_colors))
  venn_colors[is.na(venn_colors)] <- scales::hue_pal()(fallback_n)
}
# Silence VennDiagram's futile.logger chatter + stray .log file
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
venn.plot <- venn.diagram(
  x = core_species,
  filename        = NULL,
  # ── geometry: a touch of margin so circles don't touch the panel edge ──
  margin          = 0.12,
  # ── fills: soft, semi-transparent; thin coloured outlines ──
  fill            = adjustcolor(venn_colors, alpha.f = 0.35),
  col             = venn_colors,
  lwd             = 1.4,
  alpha           = 0.5,
  # ── region counts: the only numbers shown, kept modest ──
  cex             = 0.8,
  fontface        = "plain",
  fontfamily      = grid_font,
  # ── category (group) labels: smaller, coloured to match, sit just outside ──
  cat.cex         = 0.75,
  cat.fontface    = "bold",
  cat.fontfamily  = grid_font,
  cat.col         = venn_colors,
  cat.default.pos = "outer",
  cat.dist        = c(0.06, 0.06, 0.045),
  # no `main` here -> the "A)" label is drawn in assembly so it aligns with B
  main            = NULL
)
# Wrap the gList as a single gTree, then a cowplot panel.
# draw_grob with explicit x/y/width/height keeps the Venn from ballooning:
# it's placed in a centred sub-region rather than filling the whole panel.
venn_gtree <- grid::gTree(children = venn.plot)
panelA <- cowplot::ggdraw() +
  cowplot::draw_grob(venn_gtree, x = 0.04, y = 0.02,
                     width = 0.92, height = 0.92)
################################################################################
# PANEL B : Unique core taxa per group (prevalence bar plots)
################################################################################
## ---- 3B.5 Unique core taxa per group ----
unique_taxa_list <- list()
for (g in names(core_species)) {
  others <- setdiff(names(core_species), g)
  uniq <- setdiff(core_species[[g]], unlist(core_species[others]))
  if (length(uniq) > 0) {
    unique_taxa_list[[g]] <- data.frame(Group = rep(g, length(uniq)), Taxon = uniq, stringsAsFactors = FALSE)
  }
}
if (length(unique_taxa_list) > 0) {
  unique_taxa_df <- do.call(rbind, unique_taxa_list)
} else {
  unique_taxa_df <- data.frame(Group = character(), Taxon = character())
  message("[Fig3] No unique taxa found in any group.")
}
cat("[Fig3] Unique core taxa per group:",
    if (nrow(unique_taxa_df) > 0) paste(table(unique_taxa_df$Group), collapse = ", ") else "none", "\n")
## ---- 3B.6 Panel B plot ----
if (nrow(unique_taxa_df) > 0) {

  unique_prev_df <- species_clean %>%
    pivot_longer(cols = -all_of(meta_cols_core), names_to = "Taxon", values_to = "Abundance") %>%
    semi_join(unique_taxa_df, by = c("Group", "Taxon")) %>%
    group_by(Group, Taxon) %>%
    summarise(
      Prevalence = 100 * mean(Abundance >= min_abund_cutoff, na.rm = TRUE),
      MeanAbundance = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(unique_taxa_df, by = c("Group", "Taxon"))

  unique_prev_df <- unique_prev_df %>%
    arrange(Group, Prevalence, MeanAbundance) %>%
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)))

  # Canonicalise Group (trim + refactor), then attach the literal hex colour
  # as its OWN column. scale_fill_identity() below uses these hex strings
  # directly -- no name lookup, nothing for a masked scale to break.
  unique_prev_df$Group <- factor(
    trimws(as.character(unique_prev_df$Group)),
    levels = c("Healthy", "Resilient", "Unhealthy")
  )
  unique_prev_df$BarColor <- group_colors[as.character(unique_prev_df$Group)]

  # Diagnostics
  cat("\n[Fig3 PanelB] rows:", nrow(unique_prev_df),
      "| groups:", paste(levels(droplevels(unique_prev_df$Group)), collapse = ", "), "\n")
  if (anyNA(unique_prev_df$BarColor)) {
    warning("[Fig3] Some BarColor values are NA -- a Group had no match in group_colors.")
  }

  p_unique <- ggplot(unique_prev_df,
                     aes(x = Taxon, y = Prevalence, fill = BarColor)) +
    geom_col(color = "black", linewidth = 0.25, width = 0.75) +
    coord_flip() +
    facet_wrap(~ Group, scales = "free_y") +
    # fill values ARE colours -> identity scale, mask-proof
    ggplot2::scale_fill_identity() +
    labs(x = NULL, y = "Prevalence (%)") +   # no title; corner label added in assembly
    # ── SAME BASE AS FIGURE 2: theme_classic (Fig 2 uses theme_classic via
    #    theme_fig2), NOT theme_minimal. This is what makes the taxa font
    #    actually match rather than just sharing a size number. ──
    theme_classic(base_size = 16) +
    theme(
      text = element_text(family = fig_font, color = "black"),
      # species names: identical to Figure 2 taxa style
      # (italic, black, taxa_font_size, NO forced family)
      axis.text.y  = element_text(face = "italic", color = "black",
                                  size = taxa_font_size),
      axis.text.x  = element_text(family = fig_font, color = "black",
                                  size = taxa_font_size),
      axis.title   = element_text(family = fig_font, color = "black",
                                  size = taxa_font_size, face = "bold"),
      axis.line    = element_line(color = "black", linewidth = 0.5),
      axis.ticks   = element_line(color = "black", linewidth = 0.4),
      # group strip labels: PLAIN (these are group names, not Latin taxa).
      # Only the species names on the y-axis are italic.
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text   = element_text(face = "plain", color = "black",
                                  size = taxa_font_size),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      legend.position  = "none"
    )

  print(p_unique)

} else {
  message("[Fig3] No unique taxa to plot in Panel B.")
  p_unique <- cowplot::ggdraw() +
    cowplot::draw_label("No unique core taxa to plot",
                        fontfamily = fig_font, size = 12)
}
################################################################################
# 3B.7 ASSEMBLE FIGURE 3  (Panel A left | Panel B right)
################################################################################
# A/B corner labels via plot_grid's built-in `labels` (robust; avoids the
# nested-ggdraw issue that was dropping the labels).
Figure3 <- cowplot::plot_grid(
  panelA,
  p_unique,
  ncol        = 2,
  rel_widths  = c(0.8, 2),       # Venn compact; bar panel wider -> crisp bars
  labels      = c("A", "B"),
  label_size  = panel_label_size,
  label_fontfamily = fig_font,
  label_fontface   = "bold",
  label_x     = 0.01,
  label_y     = 0.99,
  hjust       = 0,
  vjust       = 1
)
print(Figure3)
## ---- 3B.8 Export combined Figure 3 ----
# Use ragg::agg_tiff for crisp, properly anti-aliased text (italic species
# names render sharp instead of blurry). Falls back to base tiff if ragg
# isn't installed -- install with install.packages("ragg") for best quality.
# if (requireNamespace("ragg", quietly = TRUE)) {
#   ggsave(
#     filename = "~/Desktop/Figure3.tiff",
#     plot     = Figure3,
#     device   = ragg::agg_tiff,
#     dpi      = 300,
#     width    = fig3_width_mm, height = fig3_height_mm, units = "mm",
#     bg       = "white"
#   )
#   cat("[Fig3] Saved with ragg (sharp text) -> ~/Desktop/Figure3.tiff\n")
# } else {
#   message("[Fig3] Package 'ragg' not found -- using base tiff (text may look ",
#           "softer). Run install.packages('ragg') for sharper labels.")
#   ggsave(
#     filename    = "~/Desktop/Figure3.tiff",
#     plot        = Figure3,
#     device      = "tiff",
#     dpi         = 300,
#     width       = fig3_width_mm, height = fig3_height_mm, units = "mm",
#     compression = "lzw",
#     bg          = "white"
#   )
#   cat("[Fig3] Saved with base tiff -> ~/Desktop/Figure3.tiff\n")
# }
################################################################################
# SUPPLEMENTARY — Shared core taxa prevalence HEATMAP (separate export)
################################################################################
## ---- 3B.9 Build + render + export supplementary heatmap ----
if (nrow(heatmap_tbl_core) > 0 && any(is.finite(heatmap_tbl_core$prevalence))) {

  heatmap_df_core <- as.data.frame(heatmap_tbl_core[, c("Taxon", "Group", "prevalence")])

  heatmap_mat_core <- reshape(heatmap_df_core, idvar = "Taxon", timevar = "Group", direction = "wide")
  colnames(heatmap_mat_core) <- gsub("^prevalence\\.", "", colnames(heatmap_mat_core))
  rownames(heatmap_mat_core) <- heatmap_mat_core$Taxon
  heatmap_mat_core$Taxon <- NULL
  heatmap_mat_core <- as.matrix(heatmap_mat_core)
  heatmap_mat_core[is.na(heatmap_mat_core)] <- 0

  row_order_core <- rownames(heatmap_mat_core)[order(rowMeans(heatmap_mat_core, na.rm = TRUE), decreasing = TRUE)]
  heatmap_mat_core <- heatmap_mat_core[row_order_core, , drop = FALSE]

  labels_row_italic <- parse(text = paste0("italic('", rownames(heatmap_mat_core), "')"))

  supp_heatmap <- pheatmap(
    heatmap_mat_core,
    color        = colorRampPalette(c("white", "#B7C9A8", "#2166AC", "#D76A49", "#762A83"))(100),
    cluster_cols = FALSE, cluster_rows = FALSE, show_rownames = TRUE,
    labels_row   = labels_row_italic, fontsize_row = 9,
    fontfamily   = grid_font,
    cellwidth    = 25, cellheight = 12, legend = TRUE,
    legend_breaks = seq(50, 100, 10),
    main         = "Supplementary Figure S2. Shared Core Taxa (Prevalence %)",
    border_color = "gray80",
    silent       = FALSE
  )

  tiff(
    filename    = "~/Desktop/FigureS_CoreSharedHeatmap.tiff",
    width       = 180, height = 220, units = "mm",
    res         = 300, compression = "lzw"
  )
  grid::grid.newpage()
  grid::grid.draw(supp_heatmap$gtable)
  dev.off()

  cat("[Fig3] Supplementary heatmap saved -> ~/Desktop/FigureS_CoreSharedHeatmap.tiff\n")

} else {
  message("[Fig3] No shared core taxa under current thresholds; skipping supplementary heatmap.")
}
cat("\n[Section 3B complete]\n")
cat("  Figure 3 = Panel A (Venn) + Panel B (unique core taxa) in one figure.\n")
cat("  Shared core taxa heatmap exported as a SUPPLEMENTARY figure.\n")
################################################################################
# 5. FIGURE 4 — Differentially abundant microbial pathways
#
#   Split by direction of change in the Unhealthy group:
#     Panel A = pathways INCREASED in Unhealthy
#     Panel B = pathways DECREASED in Unhealthy
#   Boxplots of relative abundance (%) for the top differentially abundant
#   pathways (pairwise emmeans, BH-corrected). Black bars = covariate-
#   adjusted estimated marginal means. q-values/stars on significant pairs.
#
#   Style matches Figures 1-3: uses fig_font and shared sizes.
################################################################################
## ---- 5.0 Figure 4 style + parameters ----
# fig_font is defined in Section 2.0; reuse it. Sizes mirror the other figures.
pseudo                     <- 1e-6
pathway_prevalence_cutoff  <- 0.20
top_n_pathways             <- 12
fig4_base_size   <- 12          # base theme size
fig4_strip_size  <- 8           # facet (pathway) label size — smaller so full names fit
fig4_axis_size   <- 11          # axis text/title size
fig4_stat_size   <- 2.7         # q-value/stars bracket text (geom size)
panel_label_size <- 16          # A/B corner labels (point size; tuned to fig size)
# How many characters before a pathway name wraps to the next line.
# Larger value = fewer line breaks; smaller = more wrapping to fit width.
fig4_label_wrap  <- 30
# Figure 4 export dimensions (mm)
fig4_width_mm  <- 260
fig4_height_mm <- 300
## ---- 5.0a Shared theme (fig_font, journal-style, like other figures) ----
theme_fig4 <- function(base_size = fig4_base_size) {
  theme_classic(base_size = base_size) +
    theme(
      text             = element_text(family = fig_font, color = "black"),
      plot.title       = element_text(family = fig_font, face = "bold", hjust = 0.5),
      axis.title       = element_text(family = fig_font, face = "bold", color = "black",
                                      size = fig4_axis_size),
      axis.text        = element_text(family = fig_font, color = "black",
                                      size = fig4_axis_size - 1),
      axis.line        = element_line(color = "black", linewidth = 0.5),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      legend.title     = element_text(family = fig_font, face = "bold"),
      legend.text      = element_text(family = fig_font, color = "black"),
      legend.key.size  = unit(0.30, "cm"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.35),
      strip.text       = element_text(family = fig_font, face = "bold",
                                      size = fig4_strip_size, lineheight = 0.95,
                                      margin = ggplot2::margin(4, 4, 4, 4)),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.35),
      axis.text.x      = element_text(family = fig_font, angle = 0, hjust = 1, vjust = 1),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
wrap_three_line_label <- function(x, width = 24, max_lines = 4) {
  vapply(x, function(lbl) {
    if (is.na(lbl) || !nzchar(lbl)) return("")
    wrapped <- strwrap(lbl, width = width)
    if (length(wrapped) == 0) return("")
    if (length(wrapped) <= max_lines) return(paste(wrapped, collapse = "\n"))
    paste(c(wrapped[seq_len(max_lines - 1)],
            paste(wrapped[max_lines:length(wrapped)], collapse = " ")),
          collapse = "\n")
  }, character(1))
}
make_pathway_label_lookup <- function(pathway_ids, path_anno, width = 24) {
  pathway_labels <- path_anno %>%
    dplyr::select(pathway_id, feature) %>%
    dplyr::filter(pathway_id %in% pathway_ids) %>%
    dplyr::mutate(feature_clean = wrap_three_line_label(feature, width = width))
  label_lookup <- pathway_labels$feature_clean
  names(label_lookup) <- pathway_labels$pathway_id
  out <- label_lookup[pathway_ids]
  missing_label <- is.na(out)
  out[missing_label] <- wrap_three_line_label(pathway_ids[missing_label], width = width)
  out
}
## ---- 5.1 Define pathway columns ----
metadata_vars <- c("SampleID", "unified_sample_ID", "Group", "Group_code",
                   "Sample.Name", "Age", "Sex", "BMI", "Smoking",
                   "wakefulness", "Socialized")
pathways_all <- setdiff(colnames(data_path_raw), metadata_vars)
## ---- 5.2 Prevalence filter ----
path_prev <- sapply(data_path_raw[pathways_all],
                    function(x) mean(as.numeric(x) > 0, na.rm = TRUE))
pathways <- names(path_prev[path_prev >= pathway_prevalence_cutoff])
if (any(is.na(pathways)))
  stop("NA found in 'pathways' -- check metadata_vars covers all non-pathway columns.")
cat("\n[Pathways] Kept after prevalence filter:", length(pathways),
    "of", length(pathways_all), "\n")
## ---- 5.3 Prepare data for model and plot ----
data_path_model <- data_path_raw
data_path_plot  <- data_path_raw
data_path_model[pathways] <- log10(data_path_model[pathways] + pseudo)
data_path_plot[pathways]  <- data_path_plot[pathways] * 100
## ---- 5.4 Adjusted LM + emmeans pairwise comparisons ----
pairwise_results <- list()
emmeans_results  <- list()
for (p in pathways) {
  vars_needed <- c(p, "Group", confounders_func)
  dat <- data_path_model[, vars_needed, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat <- droplevels(dat)
  if (nrow(dat) < 5 || nlevels(dat$Group) < 2) next

  factor_vars <- c("Group", "Sex", "Smoking", "wakefulness", "Socialized")
  bad_factors <- factor_vars[sapply(factor_vars,
                                    function(v) v %in% colnames(dat) && nlevels(dat[[v]]) < 2)]
  if (length(bad_factors) > 0) next

  formula_p <- as.formula(paste0("`", p, "` ~ Group + ",
                                 paste(confounders_func, collapse = " + ")))
  fit <- tryCatch(lm(formula_p, data = dat), error = function(e) NULL)
  if (is.null(fit)) next

  emm <- tryCatch(emmeans(fit, specs = "Group"), error = function(e) NULL)
  if (is.null(emm)) next

  pw_raw <- as.data.frame(pairs(emm, adjust = "none", infer = c(TRUE, TRUE))) %>%
    rename(raw_p = p.value, raw_lower.CL = lower.CL, raw_upper.CL = upper.CL)
  pw_bh  <- as.data.frame(pairs(emm, adjust = "BH", infer = c(TRUE, TRUE))) %>%
    rename(pairwise_q = p.value, q_lower.CL = lower.CL, q_upper.CL = upper.CL)
  pw <- pw_raw %>%
    dplyr::select(contrast, estimate, SE, df, t.ratio, raw_p, raw_lower.CL, raw_upper.CL) %>%
    left_join(pw_bh %>% dplyr::select(contrast, pairwise_q, q_lower.CL, q_upper.CL),
              by = "contrast") %>%
    mutate(pathway = p)

  emm_df <- as.data.frame(emm) %>% mutate(pathway = p)

  pairwise_results[[p]] <- pw
  emmeans_results[[p]]  <- emm_df
}
pairwise_results <- bind_rows(pairwise_results)
emmeans_results  <- bind_rows(emmeans_results)
cat("[Pathways] Pathways with model results:",
    length(unique(pairwise_results$pathway)), "\n")
## ---- 5.5 FDR correction ----
pairwise_results <- pairwise_results %>%
  mutate(all_tests_q = p.adjust(raw_p, method = "BH"),
         stars = case_when(pairwise_q < 0.001 ~ "***",
                           pairwise_q < 0.01  ~ "**",
                           pairwise_q < 0.05  ~ "*",
                           TRUE               ~ "ns"))
cat("[Pathways] Significant pairwise comparisons (q < 0.05):",
    sum(pairwise_results$pairwise_q < 0.05, na.rm = TRUE),
    "of", nrow(pairwise_results), "\n")
## ---- 5.6 Select top pathways by best pairwise q-value ----
top_pathway_table <- pairwise_results %>%
  group_by(pathway) %>% arrange(pairwise_q, .by_group = TRUE) %>% slice(1) %>% ungroup() %>%
  arrange(pairwise_q) %>% left_join(path_anno, by = c("pathway" = "pathway_id"))
top_pathways <- top_pathway_table %>% slice_head(n = top_n_pathways) %>% pull(pathway)
pathway_label_lookup <- make_pathway_label_lookup(top_pathways, path_anno, width = fig4_label_wrap)
## ---- 5.6a Direction of change in Unhealthy -> the A/B split key ----
pathway_direction <- emmeans_results %>%
  dplyr::filter(pathway %in% top_pathways, Group %in% c("Healthy", "Unhealthy")) %>%
  dplyr::select(pathway, Group, emmean) %>%
  tidyr::pivot_wider(names_from = Group, values_from = emmean) %>%
  dplyr::mutate(
    unhealthy_direction = dplyr::case_when(
      Unhealthy > Healthy ~ "Increased in Unhealthy",
      Unhealthy < Healthy ~ "Decreased in Unhealthy",
      TRUE                ~ "No change"
    )
  ) %>%
  dplyr::left_join(top_pathway_table %>% dplyr::select(pathway, pairwise_q), by = "pathway") %>%
  dplyr::arrange(pairwise_q)
pathways_increased <- pathway_direction %>%
  dplyr::filter(unhealthy_direction == "Increased in Unhealthy") %>% pull(pathway)
pathways_decreased <- pathway_direction %>%
  dplyr::filter(unhealthy_direction != "Increased in Unhealthy") %>% pull(pathway)
cat("[Fig4] Increased in Unhealthy:", length(pathways_increased),
    "| Decreased/other:", length(pathways_decreased), "\n")
## ---- 5.7 Builder: one direction panel (boxplots + emmeans bars + stats) ----
build_pathway_panel <- function(panel_pathways, ncol_facets) {

  if (length(panel_pathways) == 0) {
    return(cowplot::ggdraw() +
             cowplot::draw_label("No pathways in this direction",
                                 fontfamily = fig_font, size = fig4_axis_size))
  }

  lbl_levels <- pathway_label_lookup[panel_pathways]

  plot_df <- data_path_plot %>%
    dplyr::select(SampleID, Group, all_of(panel_pathways)) %>%
    pivot_longer(cols = all_of(panel_pathways),
                 names_to = "pathway", values_to = "abundance") %>%
    mutate(pathway_label = factor(pathway_label_lookup[pathway], levels = lbl_levels))

  emm_plot <- emmeans_results %>%
    filter(pathway %in% panel_pathways) %>%
    mutate(pathway_label = factor(pathway_label_lookup[pathway], levels = lbl_levels),
           emmean_plot = (10^emmean - pseudo) * 100,
           x_num = case_when(Group == "Healthy"   ~ 1,
                             Group == "Resilient" ~ 2,
                             Group == "Unhealthy" ~ 3))

  stat_df <- pairwise_results %>%
    filter(pathway %in% panel_pathways) %>%
    mutate(
      pathway_label = factor(pathway_label_lookup[pathway], levels = lbl_levels),
      group1 = str_trim(str_split_fixed(contrast, "-", 2)[, 1]),
      group2 = str_trim(str_split_fixed(contrast, "-", 2)[, 2]),
      label = dplyr::case_when(
        is.na(pairwise_q) ~ "ns",
        pairwise_q < 0.05 ~ paste0(stars, " q=", signif(pairwise_q, 3)),
        TRUE ~ "ns"
      )
    )

  y_pos <- plot_df %>% group_by(pathway) %>%
    summarise(ymax = max(abundance, na.rm = TRUE),
              ymin = min(abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(yrange = ifelse(ymax - ymin == 0, ymax + 0.1, ymax - ymin))

  stat_df <- stat_df %>% left_join(y_pos, by = "pathway") %>%
    group_by(pathway) %>% arrange(pairwise_q, .by_group = TRUE) %>%
    mutate(y.position = ymax + row_number() * 0.12 * yrange) %>% ungroup()

  ggplot(plot_df, aes(x = Group, y = abundance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75, width = 0.65) +
    geom_jitter(aes(color = Group), width = 0.15, alpha = 0.65, size = 1.6) +
    geom_segment(data = emm_plot,
                 aes(x = x_num - 0.18, xend = x_num + 0.18,
                     y = emmean_plot, yend = emmean_plot),
                 inherit.aes = FALSE, color = "black", linewidth = 1) +
    stat_pvalue_manual(stat_df, label = "label",
                       xmin = "group1", xmax = "group2", y.position = "y.position",
                       tip.length = 0.01, bracket.size = 0.4, size = fig4_stat_size,
                       family = fig_font) +
    facet_wrap(~ pathway_label, scales = "free_y", ncol = ncol_facets) +
    # Extra headroom at the top so the highest q-value bracket clears the
    # facet strip title instead of overlapping it.
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
    scale_fill_manual(values  = palette_groups) +
    scale_color_manual(values = palette_groups) +
    theme_fig4() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "right",
          plot.margin = ggplot2::margin(8, 10, 8, 6)) +
    labs(fill = "Group", color = "Group")
}
## ---- 5.8 Build Panel A (increased) and Panel B (decreased) ----
# Choose facet columns by count so neither panel looks stretched.
ncolA <- if (length(pathways_increased) >= 6) 3 else 2
ncolB <- if (length(pathways_decreased) >= 4) 2 else max(1, length(pathways_decreased))
panelA_fig4 <- build_pathway_panel(pathways_increased, ncolA)
panelB_fig4 <- build_pathway_panel(pathways_decreased, ncolB)
# Shared y-axis title (rotated), one per panel via cowplot
add_y_title <- function(p) {
  cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label("Relative abundance (%)",
                                            angle = 90, fontfamily = fig_font,
                                            fontface = "bold", size = fig4_axis_size),
    p, ncol = 2, rel_widths = c(0.045, 1)
  )
}
panelA_y <- add_y_title(panelA_fig4)
panelB_y <- add_y_title(panelB_fig4)
## ---- 5.9 Assemble Figure 4 (A increased on top, B decreased below) ----
# Heights proportional to facet rows so each box keeps a similar size.
rowsA <- ceiling(length(pathways_increased) / ncolA)
rowsB <- ceiling(max(1, length(pathways_decreased)) / ncolB)
Figure4 <- cowplot::plot_grid(
  panelA_y,
  panelB_y,
  ncol        = 1,
  rel_heights = c(rowsA, rowsB),
  labels      = c("A", "B"),
  label_size  = panel_label_size,
  label_fontfamily = fig_font,
  label_fontface   = "bold",
  label_x     = 0.005,
  label_y     = 0.995,
  hjust       = 0,
  vjust       = 1
)
print(Figure4)
## ---- 5.10 Export Figure 4 ----
# if (requireNamespace("ragg", quietly = TRUE)) {
#   ggsave(
#     filename = "~/Desktop/Figure4.tiff",
#     plot     = Figure4,
#     device   = ragg::agg_tiff,
#     dpi      = 300,
#     width    = fig4_width_mm, height = fig4_height_mm, units = "mm",
#     bg       = "white"
#   )
#   cat("[Fig4] Saved with ragg (sharp text) -> ~/Desktop/Figure4.tiff\n")
# } else {
#   message("[Fig4] Package 'ragg' not found -- using base tiff. ",
#           "Run install.packages('ragg') for sharper labels.")
#   ggsave(
#     filename    = "~/Desktop/Figure4.tiff",
#     plot        = Figure4,
#     device      = "tiff",
#     dpi         = 300,
#     width       = fig4_width_mm, height = fig4_height_mm, units = "mm",
#     compression = "lzw",
#     bg          = "white"
#   )
#   cat("[Fig4] Saved with base tiff -> ~/Desktop/Figure4.tiff\n")
# }
#
# cat("\n[Section 5 / Figure 4 complete]\n")
# cat("  Panel A:", length(pathways_increased), "pathways increased in Unhealthy\n")
# cat("  Panel B:", length(pathways_decreased), "pathways decreased/other in Unhealthy\n")
################################################################################
################################################################################
# 6. PATHWAY-TAXA CORRELATIONS — Scenario A (omnibus, per-cohort combined
#    panel) + Scenario B (pairwise, pooled heatmap)
################################################################################
clean_taxa_names_addon <- function(x) { x <- gsub(".*\\|s__", "", x); x <- gsub("^s__", "", x); gsub("_", " ", x) }
wrap_three_line_label_addon <- function(x, width = 22, max_lines = 3) {
  vapply(x, function(lbl) {
    lbl <- stringr::str_replace_all(as.character(lbl), "_", " ")
    lbl <- stringr::str_squish(lbl)
    wrapped <- strwrap(lbl, width = width)
    if (length(wrapped) == 0) return("")
    if (length(wrapped) <= max_lines) return(paste(wrapped, collapse = "\n"))
    paste(c(wrapped[seq_len(max_lines - 1)], paste(wrapped[max_lines:length(wrapped)], collapse = " ")), collapse = "\n")
  }, character(1))
}
make_pathway_label_lookup_addon <- function(pathway_ids, path_anno, width = 22) {
  pathway_labels <- path_anno %>% dplyr::select(pathway_id, feature) %>% dplyr::filter(pathway_id %in% pathway_ids) %>%
    dplyr::mutate(feature_clean = wrap_three_line_label_addon(feature, width = width))
  label_lookup <- pathway_labels$feature_clean
  names(label_lookup) <- pathway_labels$pathway_id
  out <- label_lookup[pathway_ids]
  missing_label <- is.na(out)
  out[missing_label] <- wrap_three_line_label_addon(pathway_ids[missing_label], width = width)
  out
}
get_residuals_tbl_addon <- function(df, feature, confounders, min_n = 10) {
  vars_needed <- intersect(c("SampleID", feature, confounders), colnames(df))
  dat <- df[, vars_needed, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat <- droplevels(dat)
  if (nrow(dat) < min_n) return(tibble(SampleID = character(), resid = numeric()))

  usable_covars <- confounders[sapply(confounders, function(v) {
    if (!v %in% colnames(dat)) return(FALSE)
    x <- dat[[v]]
    if (is.numeric(x)) length(unique(na.omit(x))) > 1 else nlevels(droplevels(as.factor(x))) > 1
  })]

  if (length(usable_covars) == 0) {
    return(tibble(SampleID = dat$SampleID, resid = dat[[feature]] - mean(dat[[feature]], na.rm = TRUE)))
  }

  form <- as.formula(paste0("`", feature, "` ~ ", paste(usable_covars, collapse = " + ")))
  fit <- lm(form, data = dat)
  tibble(SampleID = dat$SampleID, resid = residuals(fit))
}
match_selected_species_to_taxa_columns <- function(selected_species, data_taxa_raw) {
  selected_species <- unique(na.omit(selected_species))
  taxa_all <- setdiff(colnames(data_taxa_raw), metadata_vars)
  taxa_map <- tibble(taxa = taxa_all, taxa_clean = clean_taxa_names_addon(taxa_all))
  taxa_map %>% filter(taxa_clean %in% selected_species) %>% distinct(taxa, taxa_clean)
}
run_pathway_taxa_correlations <- function(scenario_name, selected_taxa_table, pathways_to_use = top_pathways, min_n = 10) {
  if (nrow(selected_taxa_table) == 0) { warning("No taxa matched for scenario: ", scenario_name); return(tibble()) }
  top_taxa_local <- selected_taxa_table$taxa

  data_taxa_model_local <- data_taxa_raw
  data_taxa_model_local[top_taxa_local] <- log10(data_taxa_model_local[top_taxa_local] + pseudo)

  cor_list <- list()
  for (p in pathways_to_use) {
    path_resid <- get_residuals_tbl_addon(data_path_model, p, confounders_func, min_n) %>% dplyr::rename(pathway_resid = resid)
    for (t in top_taxa_local) {
      taxa_resid <- get_residuals_tbl_addon(data_taxa_model_local, t, confounders_taxa, min_n) %>% dplyr::rename(taxa_resid = resid)
      tmp <- dplyr::inner_join(path_resid, taxa_resid, by = "SampleID")
      if (nrow(tmp) < min_n) next
      if (length(unique(tmp$pathway_resid)) < 3) next
      if (length(unique(tmp$taxa_resid)) < 3) next
      ct <- tryCatch(cor.test(tmp$pathway_resid, tmp$taxa_resid, method = "spearman", exact = FALSE), error = function(e) NULL)
      if (is.null(ct)) next
      taxa_clean_now <- selected_taxa_table$taxa_clean[match(t, selected_taxa_table$taxa)]
      ci <- spearman_ci(unname(ct$estimate), nrow(tmp))
      cor_list[[length(cor_list) + 1]] <- tibble(scenario = scenario_name, pathway = p, taxa = t, taxa_clean = taxa_clean_now,
                                                 rho = unname(ct$estimate), p_value = ct$p.value, n = nrow(tmp),
                                                 rho_lower.CL = ci["lower"], rho_upper.CL = ci["upper"])
    }
  }

  out <- dplyr::bind_rows(cor_list)
  if (nrow(out) == 0) { warning("No valid correlations for scenario: ", scenario_name); return(tibble()) }

  out %>% mutate(q_value = p.adjust(p_value, method = "BH"),
                 stars = case_when(q_value < 0.001 ~ "***", q_value < 0.01 ~ "**", q_value < 0.05 ~ "*", TRUE ~ "")) %>%
    arrange(q_value, p_value)
}
make_correlation_heatmap <- function(cor_df, figure_title, file_prefix, width = 14, height = 8) {
  if (nrow(cor_df) == 0) { message("Skipping heatmap for ", file_prefix, " -- correlation table is empty."); return(NULL) }

  cor_mat <- cor_df %>% dplyr::select(taxa_clean, pathway, rho) %>% group_by(taxa_clean, pathway) %>%
    summarise(rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pathway, values_from = rho) %>% tibble::column_to_rownames("taxa_clean") %>% as.matrix()

  star_mat <- cor_df %>% dplyr::select(taxa_clean, pathway, stars) %>% group_by(taxa_clean, pathway) %>%
    summarise(stars = dplyr::case_when(any(stars == "***", na.rm = TRUE) ~ "***", any(stars == "**", na.rm = TRUE) ~ "**",
                                       any(stars == "*", na.rm = TRUE) ~ "*", TRUE ~ ""), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = pathway, values_from = stars, values_fill = "") %>% tibble::column_to_rownames("taxa_clean") %>% as.matrix()

  star_mat <- star_mat[rownames(cor_mat), colnames(cor_mat), drop = FALSE]

  pw_label_lookup <- make_pathway_label_lookup_addon(colnames(cor_mat), path_anno, width = 22)
  colnames(cor_mat) <- pw_label_lookup[colnames(cor_mat)]
  colnames(star_mat) <- pw_label_lookup[colnames(star_mat)]

  col_fun <- circlize::colorRamp2(c(-0.8, 0, 0.8), c("#2C7BB6", "white", "#D7191C"))

  ht <- ComplexHeatmap::Heatmap(cor_mat, name = "rho", col = col_fun, cluster_rows = TRUE, cluster_columns = TRUE,
                                show_row_names = TRUE, show_column_names = TRUE, row_names_side = "right",
                                row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
                                column_names_gp = grid::gpar(fontsize = 7, lineheight = 0.9),
                                column_names_max_height = grid::unit(7, "cm"), column_names_rot = 90,
                                column_title = figure_title, column_title_gp = grid::gpar(fontsize = 13, fontface = "bold"),
                                heatmap_legend_param = list(title = "Spearman\nrho", at = c(-0.5, 0, 0.5)),
                                cell_fun = function(j, i, x, y, width, height, fill) {
                                  if (!is.na(star_mat[i, j]) && star_mat[i, j] != "") {
                                    grid::grid.text(star_mat[i, j], x, y, gp = grid::gpar(fontsize = 8, fontface = "bold"))
                                  }
                                })
  draw(ht, heatmap_legend_side = "right")
  invisible(ht)
}
################################################################################
# 6.0b  SCENARIO A ADD-ON — split by cohort, transposed, one combined panel
#
# Instead of one pooled Scenario A heatmap across all samples, this
# recomputes the omnibus pathway-taxa correlations SEPARATELY within
# each Group (Healthy / Resilient / Unhealthy), then draws all three as
# one combined ComplexHeatmap panel:
#   - Rows    = pathways, fixed `top_pathways` order, labelled ONCE
#               (on the last / right-most panel only).
#   - Columns = species, fixed alphabetical order, REPEATED under every
#               one of the three cohort sub-panels.
# Each cohort runs its own confounder-adjusted residual model and its
# own BH correction, so q-values/stars are NOT comparable across the
# three sub-panels (each is corrected only within its own cohort).
# Splitting by cohort also shrinks n a lot (e.g. Healthy/Unhealthy are
# much smaller than Resilient) -- pathway-taxa pairs that fail the
# min_n = 10 complete-case threshold in a given cohort are shown as a
# blank/grey cell for that cohort rather than being dropped, so all
# three panels keep the same row/column layout.
################################################################################
run_pathway_taxa_correlations_by_group <- function(scenario_name,
                                                   selected_taxa_table,
                                                   pathways_to_use = top_pathways,
                                                   group_value,
                                                   group_col = "Group",
                                                   min_n = 10) {
  if (nrow(selected_taxa_table) == 0) {
    warning("No taxa matched for scenario: ", scenario_name, " / group: ", group_value)
    return(tibble())
  }
  top_taxa_local <- selected_taxa_table$taxa
  
  if (!group_col %in% colnames(data_path_model)) {
    stop("Column '", group_col, "' not found in data_path_model.")
  }
  if (!group_col %in% colnames(data_taxa_raw)) {
    stop("Column '", group_col, "' not found in data_taxa_raw.")
  }
  
  # Subset BOTH the pathway table and the taxa table to this cohort only,
  # then log-transform the taxa columns exactly as the pooled version does.
  data_path_model_grp <- data_path_model %>% dplyr::filter(.data[[group_col]] == group_value)
  data_taxa_model_grp <- data_taxa_raw    %>% dplyr::filter(.data[[group_col]] == group_value)
  data_taxa_model_grp[top_taxa_local] <- log10(data_taxa_model_grp[top_taxa_local] + pseudo)
  
  n_avail <- nrow(data_path_model_grp)
  cat("  [", scenario_name, "|", group_value, "] samples available before residual/complete-case filtering:", n_avail, "\n")
  
  cor_list <- list()
  for (p in pathways_to_use) {
    path_resid <- get_residuals_tbl_addon(data_path_model_grp, p, confounders_func, min_n) %>%
      dplyr::rename(pathway_resid = resid)
    for (t in top_taxa_local) {
      taxa_resid <- get_residuals_tbl_addon(data_taxa_model_grp, t, confounders_taxa, min_n) %>%
        dplyr::rename(taxa_resid = resid)
      tmp <- dplyr::inner_join(path_resid, taxa_resid, by = "SampleID")
      if (nrow(tmp) < min_n) next
      if (length(unique(tmp$pathway_resid)) < 3) next
      if (length(unique(tmp$taxa_resid)) < 3) next
      ct <- tryCatch(cor.test(tmp$pathway_resid, tmp$taxa_resid, method = "spearman", exact = FALSE), error = function(e) NULL)
      if (is.null(ct)) next
      taxa_clean_now <- selected_taxa_table$taxa_clean[match(t, selected_taxa_table$taxa)]
      ci <- spearman_ci(unname(ct$estimate), nrow(tmp))
      cor_list[[length(cor_list) + 1]] <- tibble(
        scenario     = paste0(scenario_name, "_", group_value),
        pathway      = p,
        taxa         = t,
        taxa_clean   = taxa_clean_now,
        rho          = unname(ct$estimate),
        p_value      = ct$p.value,
        n            = nrow(tmp),
        rho_lower.CL = ci["lower"],
        rho_upper.CL = ci["upper"]
      )
    }
  }
  
  out <- dplyr::bind_rows(cor_list)
  if (nrow(out) == 0) {
    warning("No valid correlations for scenario: ", scenario_name, " / group: ", group_value)
    return(tibble())
  }
  
  out %>%
    mutate(q_value = p.adjust(p_value, method = "BH"),
           stars = case_when(q_value < 0.001 ~ "***", q_value < 0.01 ~ "**", q_value < 0.05 ~ "*", TRUE ~ "")) %>%
    arrange(q_value, p_value)
}

make_correlation_heatmap_combined_alpha <- function(cor_df_list, cohort_labels,
                                                    all_species, top_pathways,
                                                    main_title = NULL) {
  col_fun <- circlize::colorRamp2(c(-0.8, 0, 0.8), c("#2C7BB6", "white", "#D7191C"))
  ht_list <- NULL
  
  for (i in seq_along(cor_df_list)) {
    g <- cohort_labels[i]
    cor_df <- cor_df_list[[i]]
    
    if (is.null(cor_df) || nrow(cor_df) == 0) {
      cor_mat <- matrix(NA_real_, nrow = length(all_species), ncol = length(top_pathways),
                        dimnames = list(all_species, top_pathways))
      star_mat <- matrix("", nrow = length(all_species), ncol = length(top_pathways),
                         dimnames = list(all_species, top_pathways))
    } else {
      cor_mat <- cor_df %>%
        dplyr::select(taxa_clean, pathway, rho) %>%
        group_by(taxa_clean, pathway) %>%
        summarise(rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = pathway, values_from = rho) %>%
        tibble::column_to_rownames("taxa_clean") %>%
        as.matrix()
      
      star_mat <- cor_df %>%
        dplyr::select(taxa_clean, pathway, stars) %>%
        group_by(taxa_clean, pathway) %>%
        summarise(stars = dplyr::case_when(any(stars == "***", na.rm = TRUE) ~ "***",
                                           any(stars == "**",  na.rm = TRUE) ~ "**",
                                           any(stars == "*",   na.rm = TRUE) ~ "*",
                                           TRUE ~ ""), .groups = "drop") %>%
        tidyr::pivot_wider(names_from = pathway, values_from = stars, values_fill = "") %>%
        tibble::column_to_rownames("taxa_clean") %>%
        as.matrix()
      star_mat <- star_mat[rownames(cor_mat), colnames(cor_mat), drop = FALSE]
      
      missing_species <- setdiff(all_species, rownames(cor_mat))
      if (length(missing_species) > 0) {
        pad <- matrix(NA_real_, nrow = length(missing_species), ncol = ncol(cor_mat),
                      dimnames = list(missing_species, colnames(cor_mat)))
        cor_mat <- rbind(cor_mat, pad)
        pad_star <- matrix("", nrow = length(missing_species), ncol = ncol(star_mat),
                           dimnames = list(missing_species, colnames(star_mat)))
        star_mat <- rbind(star_mat, pad_star)
      }
      missing_pathways <- setdiff(top_pathways, colnames(cor_mat))
      if (length(missing_pathways) > 0) {
        pad <- matrix(NA_real_, nrow = nrow(cor_mat), ncol = length(missing_pathways),
                      dimnames = list(rownames(cor_mat), missing_pathways))
        cor_mat <- cbind(cor_mat, pad)
        pad_star <- matrix("", nrow = nrow(star_mat), ncol = length(missing_pathways),
                           dimnames = list(rownames(star_mat), missing_pathways))
        star_mat <- cbind(star_mat, pad_star)
      }
      cor_mat  <- cor_mat[all_species, top_pathways, drop = FALSE]
      star_mat <- star_mat[all_species, top_pathways, drop = FALSE]
    }
    
    # ---- TRANSPOSE: pathways -> rows (vertical axis), species -> columns
    # (horizontal axis). Up to this point cor_mat/star_mat were built as
    # species (rows) x pathway (columns), reindexed to the fixed
    # all_species / top_pathways order -- t() swaps that to
    # pathway (rows) x species (columns) while keeping that same fixed,
    # shared ordering on both sides.
    cor_mat  <- t(cor_mat)
    star_mat <- t(star_mat)
    
    pw_label_lookup <- make_pathway_label_lookup_addon(rownames(cor_mat), path_anno, width = 22)
    rownames(cor_mat)  <- pw_label_lookup[rownames(cor_mat)]
    rownames(star_mat) <- pw_label_lookup[rownames(star_mat)]
    
    is_last_panel  <- (i == length(cor_df_list))
    is_first_panel <- (i == 1)
    
    # capture this panel's own star_mat so cell_fun doesn't all point at
    # the same (last-iteration) object -- classic R closure-in-loop trap
    star_mat_i <- star_mat
    
    ht <- ComplexHeatmap::Heatmap(
      cor_mat, name = "rho", col = col_fun,
      # Row (pathway) order must be FIXED and identical across all three
      # sub-panels, because pathway labels only print once (on the last
      # panel, below, at the far right of the whole combined figure). If
      # rows were independently clustered per panel, panel 2/3 could show
      # a different top-to-bottom pathway order than the panel carrying
      # the labels -- silently mislabeling cells. Using the fixed
      # `top_pathways` order (already applied before transposing) keeps
      # every panel aligned.
      cluster_rows = FALSE,
      show_row_names = is_last_panel,      # pathway names print once only
      row_names_side = "right",
      row_names_gp = grid::gpar(fontsize = 7, lineheight = 0.9),
      # Species stay in the fixed alphabetical order set by all_species
      # (same reasoning as before: this axis repeats on every panel, so
      # keeping it alphabetized/fixed -- rather than re-clustered per
      # panel -- is what makes the three panels visually comparable).
      cluster_columns = FALSE,
      show_column_names = TRUE,            # species names repeat on every sub-panel
      column_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
      column_names_rot = 90,
      column_title = g, column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
      na_col = "grey90",
      show_heatmap_legend = is_first_panel,
      heatmap_legend_param = list(title = "Spearman\nrho", at = c(-0.5, 0, 0.5)),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (!is.na(star_mat_i[i, j]) && star_mat_i[i, j] != "") {
          grid::grid.text(star_mat_i[i, j], x, y, gp = grid::gpar(fontsize = 8, fontface = "bold"))
        }
      }
    )
    
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }
  
  ComplexHeatmap::draw(
    ht_list,
    ht_gap = grid::unit(3, "mm"),
    heatmap_legend_side = "right",
    column_title = main_title,
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  invisible(ht_list)
}
## ---- 6.1 Scenario A: omnibus Group-effect species, split by cohort ----
taxa_table_omnibus <- match_selected_species_to_taxa_columns(significant_species_from_fig2, data_taxa_raw)
cohorts_for_split <- c("Healthy", "Resilient", "Unhealthy")
cor_omnibus_by_group <- setNames(
  lapply(cohorts_for_split, function(g) {
    run_pathway_taxa_correlations_by_group(
      scenario_name        = "omnibus_group_effect",
      selected_taxa_table  = taxa_table_omnibus,
      pathways_to_use      = top_pathways,
      group_value          = g
    )
  }),
  cohorts_for_split
)
for (g in cohorts_for_split) {
  cat("\n[Pathway-taxa corr | Scenario A |", g, "] correlation rows:", nrow(cor_omnibus_by_group[[g]]), "\n")
}
# cor_omnibus kept for backward compatibility with anything downstream
# (e.g. Section 9's Excel export) that still expects the old pooled table.
cor_omnibus <- dplyr::bind_rows(cor_omnibus_by_group)
## Full, fixed, alphabetical species list shared by all three cohort panels.
all_species_alpha <- sort(unique(taxa_table_omnibus$taxa_clean))
make_correlation_heatmap_combined_alpha(
  cor_df_list   = cor_omnibus_by_group[cohorts_for_split],
  cohort_labels = cohorts_for_split,
  all_species   = all_species_alpha,
  top_pathways  = top_pathways,
  main_title    = "Pathway-species correlations (omnibus Group-significant species) by cohort"
)
################################################################################
# 7. FIGURE 8 — Functional beta diversity (pooled q-value with Fig 1B)
################################################################################
## ---- 7.1 Functional pathway matrix ----
pathway_matrix_func <- data_path_raw %>% dplyr::select(all_of(pathways)) %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric))
pathway_matrix_func <- as.data.frame(pathway_matrix_func)
rownames(pathway_matrix_func) <- data_path_raw$SampleID
## ---- 7.2 Metadata for functional beta diversity ----
meta_beta_func <- data_path_raw %>% dplyr::select(SampleID, Group, Sex, BMI, Smoking, wakefulness, Socialized)
meta_beta_func$Group       <- droplevels(factor(meta_beta_func$Group))
meta_beta_func$Sex         <- droplevels(factor(meta_beta_func$Sex))
meta_beta_func$Smoking     <- droplevels(factor(meta_beta_func$Smoking))
meta_beta_func$wakefulness <- droplevels(factor(meta_beta_func$wakefulness))
meta_beta_func$Socialized  <- droplevels(factor(meta_beta_func$Socialized))
meta_beta_func$BMI         <- as.numeric(meta_beta_func$BMI)
## ---- 7.3 Bray-Curtis distance + PCoA ----
dist_bc_func <- vegan::vegdist(pathway_matrix_func, method = "bray")
pcoa_func <- cmdscale(dist_bc_func, eig = TRUE, k = 2)
var_PC1_func <- round(100 * pcoa_func$eig[1] / sum(pcoa_func$eig[pcoa_func$eig > 0]), 1)
var_PC2_func <- round(100 * pcoa_func$eig[2] / sum(pcoa_func$eig[pcoa_func$eig > 0]), 1)
df_pcoa_func <- data.frame(SampleID = data_path_raw$SampleID, PC1 = pcoa_func$points[, 1], PC2 = pcoa_func$points[, 2], Group = meta_beta_func$Group)
## ---- 7.4 PERMDISP dispersion check ----
disp_func <- vegan::betadisper(dist_bc_func, meta_beta_func$Group)
permdisp_anova_func <- anova(disp_func)
permdisp_perm_func  <- vegan::permutest(disp_func, permutations = 999)
permdisp_pairwise_func <- TukeyHSD(disp_func)
disp_df_func <- data.frame(SampleID = meta_beta_func$SampleID, Group = meta_beta_func$Group, Distance_to_centroid = disp_func$distances)
distance_kw_func <- kruskal.test(Distance_to_centroid ~ Group, data = disp_df_func)
cat("\n[Fig8] Functional PERMDISP Kruskal-Wallis: chi-sq =", round(distance_kw_func$statistic, 3),
    ", df =", distance_kw_func$parameter, ", p =", signif(distance_kw_func$p.value, 3), "\n")
cat("(Report as one-line Methods/Results statement; figure not required.)\n")
## ---- 7.5 Adjusted global PERMANOVA + pooled q-value with Fig 1B ----
global_permanova_func <- vegan::adonis2(dist_bc_func ~ Group + Sex + BMI + Smoking + wakefulness + Socialized,
                                         data = meta_beta_func, permutations = 999)
p_val_func <- global_permanova_func$`Pr(>F)`[1]
r2_func    <- global_permanova_func$R2[1]
if (!exists("p_val_adj")) stop("p_val_adj not found. Section 2 (Figure 1B) must run before this section.")
global_permanova_pvals <- c(Taxonomic_beta = p_val_adj, Functional_beta = p_val_func)
global_permanova_qvals <- p.adjust(global_permanova_pvals, method = "fdr")
cat("\n================ POOLED GLOBAL PERMANOVA Q-VALUES ================\n")
print(data.frame(Test = names(global_permanova_pvals), p_value = signif(global_permanova_pvals, 3), q_value = signif(global_permanova_qvals, 3)))
q_val_func <- global_permanova_qvals["Functional_beta"]
label_text_func <- paste0("q = ", signif(q_val_func, 2), "\nR² = ", signif(r2_func, 2))
# Re-build Figure 1B's label_text now that the pooled q-value exists,
# so both beta-diversity panels report matching, genuinely pooled q.
q_val_adj_pooled <- global_permanova_qvals["Taxonomic_beta"]
label_text <- paste0("q = ", signif(q_val_adj_pooled, 2), "\nR² = ", signif(r2_adj, 2))
figureB <- ggplot(df_pcoa, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.7) +
  geom_segment(
    data = centroids,
    aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Group),
    arrow = arrow(length = unit(0.02, "npc")),
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = min(df_pcoa$PC1, na.rm = TRUE),
    y = max(df_pcoa$PC2, na.rm = TRUE),
    label = label_text,
    hjust = 0,
    vjust = 1,
    size = 4,
    color = "black"
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "B",
    x = paste0("PC1 (", var_PC1, "%)"),
    y = paste0("PC2 (", var_PC2, "%)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "right",
    legend.title = element_blank()
  )
print(figureB)
## ---- 7.6 Centroids + Figure 8 PCoA plot ----
centroids_func <- df_pcoa_func %>% dplyr::group_by(Group) %>%
  dplyr::summarise(PC1 = mean(PC1, na.rm = TRUE), PC2 = mean(PC2, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(Group_short = dplyr::case_when(Group == "Healthy" ~ "H", Group == "Resilient" ~ "R", Group == "Unhealthy" ~ "U", TRUE ~ as.character(Group)))
figure8_beta_function <- ggplot(df_pcoa_func, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.8, alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.7) +
  geom_segment(data = centroids_func, aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Group),
               arrow = arrow(length = unit(0.02, "npc")), linewidth = 0.9, show.legend = FALSE) +
  geom_text(data = centroids_func, aes(x = PC1, y = PC2, label = Group_short), color = "black", fontface = "bold", size = 4, vjust = -0.8, show.legend = FALSE) +
  annotate("text", x = min(df_pcoa_func$PC1, na.rm = TRUE), y = max(df_pcoa_func$PC2, na.rm = TRUE), label = label_text_func, hjust = 0, vjust = 1, size = 4) +
  scale_color_manual(values = palette_groups) +
  labs(title = "Functional beta diversity", x = paste0("PCoA1 (", var_PC1_func, "%)"), y = paste0("PCoA2 (", var_PC2_func, "%)")) +
  theme_minimal(base_size = 13) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), legend.position = "right", legend.title = element_blank(), plot.title = element_text(face = "bold"))
print(figure8_beta_function)
## ---- 7.7 Pairwise adjusted PERMANOVA (functional) ----
pairwise_permanova_func_adj <- function(distmat, metadata, group_var, confounders) {
  groups <- unique(metadata[[group_var]])
  dist_mat <- as.matrix(distmat)
  results <- data.frame()
  for (i in 1:(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      g1 <- groups[i]; g2 <- groups[j]
      subset_meta <- metadata %>% dplyr::filter(!!rlang::sym(group_var) %in% c(g1, g2))
      subset_meta[[group_var]] <- droplevels(factor(subset_meta[[group_var]]))
      subset_ids <- intersect(subset_meta$SampleID, rownames(dist_mat))
      subset_meta <- subset_meta %>% dplyr::filter(SampleID %in% subset_ids)
      subset_meta <- subset_meta[match(subset_ids, subset_meta$SampleID), ]
      sub_dist <- as.dist(dist_mat[subset_ids, subset_ids])
      usable_confounders <- confounders[sapply(confounders, function(v) {
        if (!v %in% colnames(subset_meta)) return(FALSE)
        x <- subset_meta[[v]]
        if (is.numeric(x)) length(unique(na.omit(x))) > 1 else nlevels(droplevels(factor(x))) > 1
      })]
      rhs <- paste(c(group_var, usable_confounders), collapse = " + ")
      fml <- as.formula(paste("sub_dist ~", rhs))
      ad <- vegan::adonis2(fml, data = subset_meta, permutations = 999)
      results <- rbind(results, data.frame(Group1 = as.character(g1), Group2 = as.character(g2), F_model = ad$F[1], R2 = ad$R2[1], p_value = ad$`Pr(>F)`[1]))
    }
  }
  results$Q_value <- p.adjust(results$p_value, method = "fdr")
  results
}
pairwise_permanova_func <- pairwise_permanova_func_adj(dist_bc_func, meta_beta_func, "Group", confounders_func)
# Note: dispersion figure intentionally skipped per author decision --
# the Kruskal-Wallis test in Section 7.4 (disp_df_func, distance_kw_func)
# already covers the PERMDISP diagnostic; a one-line Methods/Results
# statement is sufficient without a dedicated figure.
################################################################################
# 8. GROUP-BASED ADJUSTED ABUNDANCE + RANDOM FOREST (combined figure, A/B/C)
#    Method: Panel A re-displays Figure 4's adjusted means (no new
#    model). Panels B/C are multivariate classification importance --
#    complementary, corroborating evidence, not a hypothesis test (see
#    Methods Roadmap at top of script). Recommended placement:
#    SUPPLEMENTARY -- Panel A duplicates Figure 4 information and
#    Panel C trains on a reduced pairwise sample (Healthy n=20,
#    Unhealthy n=19 vs n=85 for the 3-class model).
################################################################################
set.seed(42)  # for reproducible RF results
## ---- 8A. Group-based adjusted abundance dot plot ----
cat("emmeans_results columns:\n")
print(colnames(emmeans_results))
ci_cols_present <- c("lower.CL", "upper.CL") %in% colnames(emmeans_results)
if (all(ci_cols_present)) {
  adj_abund_df <- emmeans_results %>%
    dplyr::filter(pathway %in% top_pathways) %>%
    dplyr::mutate(pathway_label = pathway_label_lookup[pathway],
                  emmean_pct = (10^emmean - pseudo) * 100,
                  lower_pct  = (10^lower.CL - pseudo) * 100,
                  upper_pct  = (10^upper.CL - pseudo) * 100)
} else {
  warning("Could not find lower.CL/upper.CL in emmeans_results -- dot plot will show points only, no error bars.")
  adj_abund_df <- emmeans_results %>%
    dplyr::filter(pathway %in% top_pathways) %>%
    dplyr::mutate(pathway_label = pathway_label_lookup[pathway],
                  emmean_pct = (10^emmean - pseudo) * 100,
                  lower_pct = NA_real_, upper_pct = NA_real_)
}
pathway_order_8a <- adj_abund_df %>%
  group_by(pathway_label) %>% summarise(spread = max(emmean_pct) - min(emmean_pct), .groups = "drop") %>%
  arrange(desc(spread)) %>% pull(pathway_label)
adj_abund_df$pathway_label <- factor(adj_abund_df$pathway_label, levels = rev(pathway_order_8a))
adj_abund_df$Group <- factor(adj_abund_df$Group, levels = c("Healthy", "Resilient", "Unhealthy"))
fig8A_adjusted_abundance <- ggplot(adj_abund_df, aes(x = emmean_pct, y = pathway_label, color = Group)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  { if (all(ci_cols_present))
      geom_errorbarh(aes(xmin = lower_pct, xmax = upper_pct), position = position_dodge(width = 0.6), height = 0.2)
  } +
  scale_color_manual(values = palette_groups) +
  labs(title = "A", x = "Confounder-adjusted estimated abundance (%)", y = NULL, color = "Group") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.y = element_text(size = 8))
## ---- 8B. Random forest variable importance ----
get_residuals_for_rf <- function(df, feature, confounders, min_n = 10) {
  vars_needed <- intersect(c("SampleID", "Group", feature, confounders), colnames(df))
  dat <- df[, vars_needed, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat <- droplevels(dat)
  if (nrow(dat) < min_n) return(NULL)
  usable_covars <- confounders[sapply(confounders, function(v) {
    if (!v %in% colnames(dat)) return(FALSE)
    x <- dat[[v]]
    if (is.numeric(x)) length(unique(na.omit(x))) > 1 else nlevels(droplevels(as.factor(x))) > 1
  })]
  if (length(usable_covars) == 0) {
    resid_vals <- dat[[feature]] - mean(dat[[feature]], na.rm = TRUE)
  } else {
    form <- as.formula(paste0("`", feature, "` ~ ", paste(usable_covars, collapse = " + ")))
    fit <- lm(form, data = dat)
    resid_vals <- residuals(fit)
  }
  tibble(SampleID = dat$SampleID, !!feature := resid_vals)
}
rf_resid_list <- lapply(top_pathways, function(p) get_residuals_for_rf(data_path_model, p, confounders_func))
names(rf_resid_list) <- top_pathways
rf_resid_list <- rf_resid_list[!sapply(rf_resid_list, is.null)]
rf_features_df <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "SampleID"), rf_resid_list)
rf_data <- rf_features_df %>%
  dplyr::inner_join(data_path_model %>% dplyr::select(SampleID, Group), by = "SampleID") %>%
  dplyr::filter(!is.na(Group)) %>% droplevels()
cat("\n[RF] Samples available:", nrow(rf_data), " | Features:", length(top_pathways), "\n")
run_rf_importance <- function(data, group_levels_keep, scenario_label, top_pathways, pathway_label_lookup) {
  dat <- data %>% dplyr::filter(Group %in% group_levels_keep) %>% droplevels()
  if (nrow(dat) < 10 || nlevels(dat$Group) < 2) { warning("Not enough data for RF scenario: ", scenario_label); return(NULL) }
  feature_cols <- intersect(top_pathways, colnames(dat))
  rf_fit <- randomForest::randomForest(x = dat[, feature_cols, drop = FALSE], y = dat$Group, ntree = 1000, importance = TRUE)
  imp <- randomForest::importance(rf_fit, type = 1)
  imp_df <- tibble(pathway = rownames(imp), MeanDecreaseAccuracy = imp[, "MeanDecreaseAccuracy"]) %>%
    mutate(pathway_label = pathway_label_lookup[pathway], scenario = scenario_label) %>%
    arrange(desc(MeanDecreaseAccuracy))
  list(model = rf_fit, importance = imp_df)
}
rf_3class <- run_rf_importance(rf_data, c("Healthy", "Resilient", "Unhealthy"), "3-class (H/R/U)", top_pathways, pathway_label_lookup)
if (!is.null(rf_3class)) {
  print(rf_3class$model)
  cat("\n[RF 3-class] OOB error rate:", round(rf_3class$model$err.rate[nrow(rf_3class$model$err.rate), "OOB"] * 100, 1), "%\n")
}
# CAVEAT: pairwise model trains on a reduced sample (Healthy + Unhealthy
# only, n approx 39 vs n approx 85 for the 3-class model). Treat pairwise
# importance rankings as more exploratory; mention sample size if reported.
rf_pairwise <- run_rf_importance(rf_data, c("Healthy", "Unhealthy"), "Pairwise (Healthy vs Unhealthy)", top_pathways, pathway_label_lookup)
if (!is.null(rf_pairwise)) {
  print(rf_pairwise$model)
  cat("\n[RF pairwise] OOB error rate:", round(rf_pairwise$model$err.rate[nrow(rf_pairwise$model$err.rate), "OOB"] * 100, 1), "%\n")
  cat("[RF pairwise] n =", nrow(rf_data %>% dplyr::filter(Group %in% c("Healthy", "Unhealthy"))), "(treat as exploratory)\n")
}
rf_importance_combined <- dplyr::bind_rows(
  if (!is.null(rf_3class)) rf_3class$importance else NULL,
  if (!is.null(rf_pairwise)) rf_pairwise$importance else NULL
)
plot_rf_dotplot <- function(imp_df, title_label) {
  imp_df <- imp_df %>% arrange(MeanDecreaseAccuracy) %>% mutate(pathway_label = factor(pathway_label, levels = pathway_label))
  ggplot(imp_df, aes(x = MeanDecreaseAccuracy, y = pathway_label)) +
    geom_segment(aes(x = 0, xend = MeanDecreaseAccuracy, y = pathway_label, yend = pathway_label), color = "grey70", linewidth = 0.5) +
    geom_point(size = 3, color = "#D7191C") +
    labs(title = title_label, x = "Mean decrease in accuracy", y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.y = element_text(size = 8))
}
if (!is.null(rf_3class))   fig8B_rf_3class   <- plot_rf_dotplot(rf_3class$importance,   "B")
if (!is.null(rf_pairwise)) fig8B_rf_pairwise <- plot_rf_dotplot(rf_pairwise$importance, "C")
## ---- 8C. Combined A/B/C figure ----
if (!is.null(rf_3class) && !is.null(rf_pairwise)) {
  fig8_combined <- (fig8A_adjusted_abundance | (fig8B_rf_3class / fig8B_rf_pairwise)) +
    plot_annotation(tag_levels = "A")
  print(fig8_combined)
}
cat("\n[Section 8 complete] Objects: fig8A_adjusted_abundance, rf_3class, ",
    "rf_pairwise, rf_importance_combined, fig8_combined\n")
################################################################################
# 9. SUPPLEMENTARY TABLES — all export/write_xlsx calls collected here
################################################################################
# Uncomment to write files. Each list element becomes one Excel sheet.
#
# NOTE: your pasted script had a second copy of this block further down
# with these write_xlsx() calls NOT commented out (i.e. set to run live).
# I've kept them commented here, matching the "Uncomment to write files"
# instruction and the safer default -- uncomment deliberately if you want
# these to actually write to ~/Desktop.
#
# write_xlsx(
#   list(
#     "S1_Alpha_Data"                = alpha_df,
#     "S2_Alpha_Adjusted_Models"     = adj_results,
#     "S3_PERMDISP_Taxonomic"        = as.data.frame(permdisp_test),
#     "S4_PERMANOVA_Global_Taxa"     = as.data.frame(global_permanova_adj),
#     "S5_PERMANOVA_Pairwise_Taxa"   = pairwise_results_adj,
#     "S6_PCoA_Coordinates_Taxa"     = df_pcoa
#   ),
#   path = "~/Desktop/Supplementary_Figure1_Alpha_Beta.xlsx"
# )
#
# write_xlsx(
#   list(
#     "S7_Phylum_All_Results"      = lm_results_p,
#     "S8_Phylum_Significant"      = sig_phyla_filtered,
#     "S9_Species_Top30_Results"   = lm_results_s,
#     "S10_Species_Sig_Omnibus"    = sig_species_omnibus,
#     "S11_Phylum_Bar_Data"        = df_bar_p,
#     "S12_Species_Bar_Data"       = df_bar_s
#   ),
#   path = "~/Desktop/Supplementary_Figure2_Phyla_Species.xlsx"
# )
#
# write_xlsx(
#   list(
#     "S13_PartialSpearman_Species_All" = partial_spearman_species,
#     "S14_PartialSpearman_Species_Sig" = sig_species_pairwise,
#     "S15_PartialSpearman_Phylum_All"  = partial_spearman_phylum,
#     "S16_PartialSpearman_Phylum_Sig"  = sig_phylum_pairwise
#   ),
#   path = "~/Desktop/Supplementary_PartialSpearman.xlsx"
# )
#
# write_xlsx(
#   list(
#     "S17_Pathway_All_Pairwise"   = pairwise_results,
#     "S18_Pathway_Top_Table"      = top_pathway_table,
#     "S19_Pathway_Adjusted_Means" = emmeans_results
#   ),
#   path = "~/Desktop/Supplementary_Functional_Pathways.xlsx"
# )
#
# write_xlsx(
#   list(
#     "S20_PathwayTaxa_Corr_Omnibus"  = cor_omnibus,
#     "S21_PathwayTaxa_Corr_Pairwise" = cor_pairwise
#   ),
#   path = "~/Desktop/Supplementary_Pathway_Taxa_Correlations.xlsx"
# )
#
# write_xlsx(
#   list(
#     "S22_PERMDISP_Functional"      = as.data.frame(permdisp_anova_func),
#     "S23_PERMANOVA_Global_Func"    = as.data.frame(global_permanova_func),
#     "S24_PERMANOVA_Pairwise_Func"  = pairwise_permanova_func,
#     "S25_PCoA_Coordinates_Func"    = df_pcoa_func,
#     "S26_Distance_to_Centroid"     = disp_df_func
#   ),
#   path = "~/Desktop/Supplementary_Figure8_Functional_Beta.xlsx"
# )
#
# write_xlsx(
#   list(
#     "S27_Adjusted_Abundance_Dotplot" = adj_abund_df,
#     "S28_RF_Importance_3class"       = if (!is.null(rf_3class)) rf_3class$importance else data.frame(),
#     "S29_RF_Importance_Pairwise"     = if (!is.null(rf_pairwise)) rf_pairwise$importance else data.frame()
#   ),
#   path = "~/Desktop/Supplementary_RandomForest_AdjustedAbundance.xlsx"
# )
################################################################################
# END OF MAIN SCRIPT (Sections 1-9)
################################################################################

################################################################################
################################################################################
# APPENDED SCRIPT: FIGURE 5 (Option A / Option B) — self-contained add-on
#
# This was a separate, self-contained script in your paste (its own
# library() calls, its own "_fig5"-suffixed helper functions so nothing
# collides with Sections 1-9 above). I've kept it as-is, appended at the
# end, since nothing above depends on it and it guards its own
# prerequisites with exists()/stop() checks.
#
# IMPORTANT: this block calls stop() if `partial_spearman_species` or
# `sig_species_pairwise` don't exist yet (see Section 6.1 "Safety
# checks" below) -- same as Section 6.2 above (Scenario B) already
# required them. Neither your pasted script nor this reconstruction
# defines those two objects; they come from a "Section 4" partial-
# Spearman analysis that wasn't included in what you pasted. Run
# whatever script builds partial_spearman_species / sig_species_pairwise
# (and partial_spearman_phylum / sig_phylum_pairwise, referenced in
# Section 9's commented S13-S16 export) before this block, or before
# Section 6.2 above -- both need it.
################################################################################
suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(cowplot)
  library(grid)
  library(ComplexHeatmap)
  library(circlize)
})
################################################################################
# 0. SETTINGS
################################################################################
fig5_out_dir <- "~/Desktop"
fig5_optionA_file <- file.path(fig5_out_dir, "Figure5_OptionA.tiff")
fig5_optionB_file <- file.path(fig5_out_dir, "Figure5_OptionB.tiff")
fig5_width_mm  <- 280
fig5_height_mm <- 360
fig5_dpi       <- 300
if (!exists("fig_font")) fig_font <- ""
if (!exists("palette_groups")) {
  palette_groups <- c(
    Healthy   = "#1A9850",
    Resilient = "#2C7BB6",
    Unhealthy = "#D7191C"
  )
}
if (!exists("pseudo")) pseudo <- 1e-6
if (!exists("metadata_vars")) {
  metadata_vars <- c(
    "SampleID", "unified_sample_ID", "Group", "Group_code",
    "Sample.Name", "Age", "Sex", "BMI", "Smoking",
    "wakefulness", "Socialized"
  )
}
if (!exists("confounders_func")) {
  confounders_func <- c("Sex", "BMI", "Smoking", "wakefulness", "Socialized")
}
if (!exists("confounders_taxa")) {
  confounders_taxa <- c("BMI", "Socialized", "Smoking", "wakefulness", "Sex")
}
################################################################################
# 1. SAFETY CHECKS
################################################################################
needed_fig5 <- c(
  "data_path_raw",
  "path_anno",
  "species_pct",
  "data_taxa_raw"
)
missing_fig5 <- needed_fig5[!vapply(needed_fig5, exists, logical(1))]
if (length(missing_fig5) > 0) {
  stop(
    "Missing object(s): ",
    paste(missing_fig5, collapse = ", "),
    "\nRun the setup, Figure 2, Figure 4, and partial Spearman sections first."
  )
}
if (!"pathway_id" %in% colnames(path_anno)) {
  stop("path_anno must contain column: pathway_id")
}
if (!"feature" %in% colnames(path_anno)) {
  stop("path_anno must contain column: feature")
}
if (!all(c("SampleID", "Group") %in% colnames(data_path_raw))) {
  stop("data_path_raw must contain SampleID and Group columns.")
}
if (!all(c("SampleID", "Group", "Species", "Abundance") %in% colnames(species_pct))) {
  stop("species_pct must contain SampleID, Group, Species, and Abundance columns.")
}
################################################################################
# 2. HELPER FUNCTIONS
################################################################################
theme_fig5 <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      text             = element_text(family = fig_font, color = "black"),
      axis.text        = element_text(family = fig_font, color = "black"),
      axis.title       = element_text(family = fig_font, face = "bold", color = "black"),
      axis.line        = element_line(color = "black", linewidth = 0.5),
      axis.ticks       = element_line(color = "black", linewidth = 0.4),
      legend.text      = element_text(family = fig_font, color = "black"),
      legend.title     = element_blank(),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.45),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
open_tiff_fig5 <- function(filename, width_mm, height_mm, dpi = 300) {
  width_in  <- width_mm / 25.4
  height_in <- height_mm / 25.4

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(
      filename    = filename,
      width       = width_in,
      height      = height_in,
      units       = "in",
      res         = dpi,
      compression = "lzw",
      background  = "white"
    )
  } else {
    grDevices::tiff(
      filename    = filename,
      width       = width_in,
      height      = height_in,
      units       = "in",
      res         = dpi,
      compression = "lzw",
      bg          = "white"
    )
  }
}
wrap_label_fig5 <- function(x, width = 24, max_lines = 4) {
  vapply(x, function(lbl) {
    lbl <- as.character(lbl)
    lbl <- stringr::str_replace_all(lbl, "_", " ")
    lbl <- stringr::str_squish(lbl)

    wrapped <- strwrap(lbl, width = width)

    if (length(wrapped) == 0) return("")

    if (length(wrapped) <= max_lines) {
      return(paste(wrapped, collapse = "\n"))
    }

    paste(
      c(
        wrapped[seq_len(max_lines - 1)],
        paste(wrapped[max_lines:length(wrapped)], collapse = " ")
      ),
      collapse = "\n"
    )
  }, character(1))
}
make_pathway_label_lookup_fig5 <- function(pathway_ids, path_anno, width = 24) {
  pathway_ids <- unique(as.character(pathway_ids))

  pathway_labels <- path_anno %>%
    dplyr::select(pathway_id, feature) %>%
    dplyr::filter(pathway_id %in% pathway_ids) %>%
    dplyr::mutate(feature_clean = wrap_label_fig5(feature, width = width))

  lookup <- pathway_labels$feature_clean
  names(lookup) <- pathway_labels$pathway_id

  out <- lookup[pathway_ids]
  missing_label <- is.na(out)

  out[missing_label] <- wrap_label_fig5(pathway_ids[missing_label], width = width)

  out
}
clean_species_from_taxa_col_fig5 <- function(x) {
  x <- as.character(x)
  x <- gsub(".*\\|s__", "", x)
  x <- gsub("^s__", "", x)
  x <- gsub("_", " ", x)
  x <- stringr::str_squish(x)
  x
}
stars_from_q_fig5 <- function(q) {
  dplyr::case_when(
    is.na(q)   ~ "",
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    TRUE      ~ ""
  )
}
zscore_rows_fig5 <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  mat_z
}
get_residuals_tbl_fig5 <- function(df, feature, confounders, min_n = 10) {
  vars_needed <- intersect(c("SampleID", feature, confounders), colnames(df))

  dat <- df[, vars_needed, drop = FALSE]
  dat <- dat[complete.cases(dat), , drop = FALSE]
  dat <- droplevels(dat)

  if (nrow(dat) < min_n) {
    return(tibble::tibble(SampleID = character(), resid = numeric()))
  }

  usable_covars <- confounders[vapply(confounders, function(v) {
    if (!v %in% colnames(dat)) return(FALSE)

    x <- dat[[v]]

    if (is.numeric(x)) {
      length(unique(stats::na.omit(x))) > 1
    } else {
      nlevels(droplevels(as.factor(x))) > 1
    }
  }, logical(1))]

  if (length(usable_covars) == 0) {
    return(
      tibble::tibble(
        SampleID = dat$SampleID,
        resid    = dat[[feature]] - mean(dat[[feature]], na.rm = TRUE)
      )
    )
  }

  form <- stats::as.formula(
    paste0("`", feature, "` ~ ", paste(usable_covars, collapse = " + "))
  )

  fit <- stats::lm(form, data = dat)

  tibble::tibble(
    SampleID = dat$SampleID,
    resid    = stats::residuals(fit)
  )
}
match_species_to_taxa_columns_fig5 <- function(selected_species, data_taxa_raw) {
  selected_species <- unique(stats::na.omit(as.character(selected_species)))

  taxa_all <- setdiff(colnames(data_taxa_raw), metadata_vars)

  taxa_map <- tibble::tibble(
    taxa       = taxa_all,
    taxa_clean = clean_species_from_taxa_col_fig5(taxa_all)
  )

  taxa_map %>%
    dplyr::filter(taxa_clean %in% selected_species) %>%
    dplyr::distinct(taxa, taxa_clean)
}
if (!exists("spearman_ci")) {
  spearman_ci <- function(rho, n, conf_level = 0.95) {
    if (is.na(rho) || is.na(n) || n < 4) {
      return(c(lower = NA_real_, upper = NA_real_))
    }

    z <- atanh(rho)
    se_z <- 1 / sqrt(n - 3)
    z_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

    c(
      lower = tanh(z - z_crit * se_z),
      upper = tanh(z + z_crit * se_z)
    )
  }
}
capture_heatmap_grob_fig5 <- function(ht) {
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      padding = grid::unit(c(2, 2, 2, 2), "mm")
    )
  )
}
################################################################################
# 3. PREPARE PATHWAY OBJECTS
################################################################################
pathways_fig5 <- intersect(path_anno$pathway_id, colnames(data_path_raw))
if (length(pathways_fig5) == 0) {
  stop("No pathway columns from path_anno$pathway_id were found in data_path_raw.")
}
# Use top_pathways from your Figure 4 if it exists.
# Otherwise, choose the first 12 pathways from pathway annotation.
if (!exists("top_pathways") || length(top_pathways) == 0) {
  message("top_pathways not found. Using first 12 pathway columns as fallback.")
  top_pathways <- head(pathways_fig5, 12)
}
top_pathways <- intersect(top_pathways, pathways_fig5)
if (length(top_pathways) == 0) {
  stop("top_pathways did not match pathway columns in data_path_raw.")
}
pathway_label_lookup_fig5 <- make_pathway_label_lookup_fig5(
  pathway_ids = top_pathways,
  path_anno   = path_anno,
  width       = 24
)
# Make model-scale pathway table if your Figure 4 object does not already exist.
if (!exists("data_path_model")) {
  data_path_model <- data_path_raw
  data_path_model[pathways_fig5] <- lapply(
    data_path_model[pathways_fig5],
    function(x) log10(as.numeric(x) + pseudo)
  )
}
################################################################################
# 4. PANEL A FOR BOTH OPTIONS — FUNCTIONAL BETA DIVERSITY
################################################################################
pathway_matrix_func_fig5 <- data_path_raw %>%
  dplyr::select(dplyr::all_of(pathways_fig5)) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
  as.data.frame()
rownames(pathway_matrix_func_fig5) <- data_path_raw$SampleID
meta_beta_func_fig5 <- data_path_raw %>%
  dplyr::select(SampleID, Group, dplyr::all_of(intersect(confounders_func, colnames(data_path_raw)))) %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Healthy", "Resilient", "Unhealthy")),
    BMI   = if ("BMI" %in% colnames(.)) as.numeric(BMI) else BMI
  )
factor_covars_fig5 <- intersect(
  c("Sex", "Smoking", "wakefulness", "Socialized"),
  colnames(meta_beta_func_fig5)
)
for (v in factor_covars_fig5) {
  meta_beta_func_fig5[[v]] <- droplevels(factor(meta_beta_func_fig5[[v]]))
}
complete_rows_fig5 <- complete.cases(meta_beta_func_fig5) &
  complete.cases(pathway_matrix_func_fig5)
meta_beta_func_fig5 <- meta_beta_func_fig5[complete_rows_fig5, , drop = FALSE]
pathway_matrix_func_fig5 <- pathway_matrix_func_fig5[complete_rows_fig5, , drop = FALSE]
rownames(pathway_matrix_func_fig5) <- meta_beta_func_fig5$SampleID
dist_bc_func_fig5 <- vegan::vegdist(pathway_matrix_func_fig5, method = "bray")
pcoa_func_fig5 <- stats::cmdscale(dist_bc_func_fig5, eig = TRUE, k = 2)
eig_pos_fig5 <- pcoa_func_fig5$eig[pcoa_func_fig5$eig > 0]
var_PC1_func_fig5 <- round(100 * pcoa_func_fig5$eig[1] / sum(eig_pos_fig5), 1)
var_PC2_func_fig5 <- round(100 * pcoa_func_fig5$eig[2] / sum(eig_pos_fig5), 1)
df_pcoa_func_fig5 <- data.frame(
  SampleID = meta_beta_func_fig5$SampleID,
  PC1      = pcoa_func_fig5$points[, 1],
  PC2      = pcoa_func_fig5$points[, 2],
  Group    = meta_beta_func_fig5$Group
)
# PERMDISP
disp_func_fig5 <- vegan::betadisper(dist_bc_func_fig5, meta_beta_func_fig5$Group)
permdisp_anova_func_fig5 <- anova(disp_func_fig5)
permdisp_perm_func_fig5  <- vegan::permutest(disp_func_fig5, permutations = 999)
cat("\n[Figure 5 Panel A] Functional PERMDISP ANOVA:\n")
print(permdisp_anova_func_fig5)
cat("\n[Figure 5 Panel A] Functional PERMDISP permutation test:\n")
print(permdisp_perm_func_fig5)
# PERMANOVA
usable_func_covars_fig5 <- confounders_func[vapply(confounders_func, function(v) {
  if (!v %in% colnames(meta_beta_func_fig5)) return(FALSE)

  x <- meta_beta_func_fig5[[v]]

  if (is.numeric(x)) {
    length(unique(stats::na.omit(x))) > 1
  } else {
    nlevels(droplevels(factor(x))) > 1
  }
}, logical(1))]
permanova_rhs_fig5 <- paste(c("Group", usable_func_covars_fig5), collapse = " + ")
global_permanova_func_fig5 <- vegan::adonis2(
  stats::as.formula(paste("dist_bc_func_fig5 ~", permanova_rhs_fig5)),
  data = meta_beta_func_fig5,
  permutations = 999
)
p_val_func_fig5 <- global_permanova_func_fig5$`Pr(>F)`[1]
r2_func_fig5    <- global_permanova_func_fig5$R2[1]
# If your previous script already pooled taxonomic + functional PERMANOVA q-values,
# use that q-value. Otherwise use the functional PERMANOVA p-value as a single-test q.
if (exists("global_permanova_qvals") &&
    "Functional_beta" %in% names(global_permanova_qvals)) {
  q_val_func_fig5 <- global_permanova_qvals["Functional_beta"]
} else {
  q_val_func_fig5 <- stats::p.adjust(p_val_func_fig5, method = "BH")
}
label_text_func_fig5 <- paste0(
  "q = ", signif(q_val_func_fig5, 2),
  "\nR² = ", signif(r2_func_fig5, 2)
)
centroids_func_fig5 <- df_pcoa_func_fig5 %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  )
fig5_panelA <- ggplot(df_pcoa_func_fig5, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.2, alpha = 0.88) +
  stat_ellipse(type = "norm", level = 0.68, linewidth = 0.65) +
  geom_segment(
    data = centroids_func_fig5,
    aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Group),
    arrow = arrow(length = unit(0.016, "npc"), type = "closed"),
    linewidth = 0.75,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = min(df_pcoa_func_fig5$PC1, na.rm = TRUE),
    y = max(df_pcoa_func_fig5$PC2, na.rm = TRUE),
    label = label_text_func_fig5,
    hjust = 0,
    vjust = 1,
    size = 4.2,
    family = fig_font,
    color = "black"
  ) +
  scale_color_manual(values = palette_groups) +
  labs(
    x = paste0("PCoA1 (", var_PC1_func_fig5, "%)"),
    y = paste0("PCoA2 (", var_PC2_func_fig5, "%)")
  ) +
  theme_fig5(base_size = 12) +
  theme(
    legend.position = "right",
    plot.margin = margin(6, 6, 6, 6)
  )
fig5_panelA_grob <- ggplotGrob(fig5_panelA)
################################################################################
# 5. FIGURE 5 OPTION A — PANEL B
#    Group-level abundance heatmap of omnibus FDR-significant species
################################################################################
get_omnibus_species_fig5 <- function() {
  if (exists("sig_species_omnibus")) {
    if ("qval" %in% colnames(sig_species_omnibus)) {
      out <- sig_species_omnibus %>%
        dplyr::filter(!is.na(qval), qval < 0.05) %>%
        dplyr::pull(Species) %>%
        unique()
      return(out)
    }
  }

  if (exists("significant_species_from_fig2")) {
    return(unique(stats::na.omit(significant_species_from_fig2)))
  }

  stop(
    "No omnibus significant species object found. ",
    "Need sig_species_omnibus or significant_species_from_fig2 from Figure 2."
  )
}
omnibus_species_fig5 <- get_omnibus_species_fig5()
omnibus_species_fig5 <- unique(stats::na.omit(omnibus_species_fig5))
if (length(omnibus_species_fig5) == 0) {
  stop("No FDR-significant omnibus species found for Figure 5 Option A.")
}
species_heat_omnibus_fig5 <- species_pct %>%
  dplyr::filter(Species %in% omnibus_species_fig5) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Healthy", "Resilient", "Unhealthy")))
heat_mat_omnibus_fig5 <- species_heat_omnibus_fig5 %>%
  dplyr::group_by(Species, Group) %>%
  dplyr::summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = mean_abundance,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Species") %>%
  as.matrix()
heat_mat_omnibus_fig5 <- heat_mat_omnibus_fig5[
  ,
  intersect(c("Healthy", "Resilient", "Unhealthy"), colnames(heat_mat_omnibus_fig5)),
  drop = FALSE
]
heat_mat_omnibus_z_fig5 <- zscore_rows_fig5(heat_mat_omnibus_fig5)
col_fun_z_fig5 <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#2C7BB6", "white", "#D7191C")
)
fig5A_panelB_ht <- ComplexHeatmap::Heatmap(
  heat_mat_omnibus_z_fig5,
  name = "Row\nz-score",
  col = col_fun_z_fig5,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
  column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Row\nz-score",
    at = c(-2, 0, 2)
  )
)
fig5A_panelB_grob <- capture_heatmap_grob_fig5(fig5A_panelB_ht)
################################################################################
# 6. FIGURE 5 OPTION B — PANEL B
#    Partial Spearman associations with pairwise health-status comparisons
################################################################################
if (!exists("partial_spearman_species")) {
  stop("partial_spearman_species not found. Run your pairwise partial Spearman section first.")
}
if (!exists("sig_species_pairwise")) {
  stop("sig_species_pairwise not found. Run your pairwise partial Spearman section first.")
}
cat("\n[Fig5B Panel B] partial_spearman_species columns:\n")
print(colnames(partial_spearman_species))
cat("\n[Fig5B Panel B] sig_species_pairwise columns:\n")
print(colnames(sig_species_pairwise))
sig_species_pairwise_fig5 <- sig_species_pairwise %>%
  dplyr::filter(!is.na(q_value), q_value < 0.05)
if (nrow(sig_species_pairwise_fig5) == 0) {
  stop("No FDR-significant pairwise species found for Figure 5 Option B.")
}
pairwise_species_keep_fig5 <- unique(sig_species_pairwise_fig5$Species)
partial_spearman_species_fig5 <- partial_spearman_species
################################################################################
# 6.1 Detect comparison / contrast column safely
################################################################################
possible_contrast_cols <- c(
  "contrast", "Contrast",
  "comparison", "Comparison",
  "pairwise", "Pairwise",
  "pairwise_comparison", "Pairwise_comparison",
  "Pairwise_Comparison",
  "comparison_name", "Comparison_name",
  "Comparison_Name",
  "group_comparison", "Group_comparison",
  "Group_Comparison"
)
contrast_col_found <- possible_contrast_cols[
  possible_contrast_cols %in% colnames(partial_spearman_species_fig5)
]
if (length(contrast_col_found) >= 1) {

  contrast_col_found <- contrast_col_found[1]

  partial_spearman_species_fig5 <- partial_spearman_species_fig5 %>%
    dplyr::mutate(
      contrast_fig5 = as.character(.data[[contrast_col_found]])
    )

} else if (all(c("group1", "group2") %in% colnames(partial_spearman_species_fig5))) {

  partial_spearman_species_fig5 <- partial_spearman_species_fig5 %>%
    dplyr::mutate(
      contrast_fig5 = paste(as.character(.data$group1), as.character(.data$group2), sep = " vs ")
    )

} else if (all(c("Group1", "Group2") %in% colnames(partial_spearman_species_fig5))) {

  partial_spearman_species_fig5 <- partial_spearman_species_fig5 %>%
    dplyr::mutate(
      contrast_fig5 = paste(as.character(.data$Group1), as.character(.data$Group2), sep = " vs ")
    )

} else {

  stop(
    "\nCould not find a contrast/comparison column.\n\n",
    "Your current partial_spearman_species columns are:\n",
    paste(colnames(partial_spearman_species_fig5), collapse = ", "),
    "\n\nPlease send me this printed column list, or rename your comparison column to 'contrast'."
  )
}
################################################################################
# 6.2 Standardize contrast names
################################################################################
partial_spearman_species_fig5 <- partial_spearman_species_fig5 %>%
  dplyr::mutate(
    contrast_fig5 = stringr::str_squish(as.character(contrast_fig5)),
    contrast_fig5 = stringr::str_replace_all(contrast_fig5, " - ", " vs "),
    contrast_fig5 = stringr::str_replace_all(contrast_fig5, "_vs_", " vs "),
    contrast_fig5 = stringr::str_replace_all(contrast_fig5, "_", " ")
  )
################################################################################
# 6.3 Make clean comparison labels
################################################################################
pairwise_heat_df_fig5 <- partial_spearman_species_fig5 %>%
  dplyr::filter(Species %in% pairwise_species_keep_fig5) %>%
  dplyr::mutate(
    contrast_clean = dplyr::case_when(
      contrast_fig5 %in% c("Healthy vs Resilient", "Resilient vs Healthy") ~ "Resilient\nvs Healthy",
      contrast_fig5 %in% c("Healthy vs Unhealthy", "Unhealthy vs Healthy") ~ "Unhealthy\nvs Healthy",
      contrast_fig5 %in% c("Resilient vs Unhealthy", "Unhealthy vs Resilient") ~ "Unhealthy\nvs Resilient",
      TRUE ~ contrast_fig5
    ),
    contrast_clean = factor(
      contrast_clean,
      levels = c(
        "Resilient\nvs Healthy",
        "Unhealthy\nvs Healthy",
        "Unhealthy\nvs Resilient"
      )
    ),
    stars = stars_from_q_fig5(q_value)
  ) %>%
  dplyr::filter(!is.na(contrast_clean))
if (nrow(pairwise_heat_df_fig5) == 0) {
  stop(
    "pairwise_heat_df_fig5 is empty after cleaning contrast names.\n",
    "Check the printed contrast values using:\n",
    "unique(partial_spearman_species_fig5$contrast_fig5)"
  )
}
################################################################################
# 6.4 Build rho and star matrices
################################################################################
rho_mat_pairwise_fig5 <- pairwise_heat_df_fig5 %>%
  dplyr::select(Species, contrast_clean, rho) %>%
  tidyr::pivot_wider(
    names_from = contrast_clean,
    values_from = rho,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Species") %>%
  as.matrix()
star_mat_pairwise_fig5 <- pairwise_heat_df_fig5 %>%
  dplyr::select(Species, contrast_clean, stars) %>%
  tidyr::pivot_wider(
    names_from = contrast_clean,
    values_from = stars,
    values_fill = ""
  ) %>%
  tibble::column_to_rownames("Species") %>%
  as.matrix()
star_mat_pairwise_fig5 <- star_mat_pairwise_fig5[
  rownames(rho_mat_pairwise_fig5),
  colnames(rho_mat_pairwise_fig5),
  drop = FALSE
]
################################################################################
# 6.5 Heatmap
################################################################################
col_fun_rho_fig5 <- circlize::colorRamp2(
  c(-0.8, 0, 0.8),
  c("#2C7BB6", "white", "#D7191C")
)
fig5B_panelB_ht <- ComplexHeatmap::Heatmap(
  rho_mat_pairwise_fig5,
  name = "Partial\nrho",
  col = col_fun_rho_fig5,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
  column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Partial\nrho",
    at = c(-0.5, 0, 0.5)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(star_mat_pairwise_fig5[i, j]) &&
        star_mat_pairwise_fig5[i, j] != "") {
      grid::grid.text(
        star_mat_pairwise_fig5[i, j],
        x,
        y,
        gp = grid::gpar(fontsize = 9, fontface = "bold")
      )
    }
  }
)
fig5B_panelB_grob <- capture_heatmap_grob_fig5(fig5B_panelB_ht)
################################################################################
# 7. PANEL C FOR BOTH OPTIONS
#    Residualized pathway–species correlations
################################################################################
run_pathway_species_cor_fig5 <- function(selected_species,
                                         scenario_name,
                                         pathways_to_use = top_pathways,
                                         min_n = 10) {
  selected_species <- unique(stats::na.omit(as.character(selected_species)))

  selected_taxa_table <- match_species_to_taxa_columns_fig5(
    selected_species = selected_species,
    data_taxa_raw    = data_taxa_raw
  )

  if (nrow(selected_taxa_table) == 0) {
    warning("No selected species matched data_taxa_raw columns for: ", scenario_name)
    return(tibble::tibble())
  }

  pathways_to_use <- intersect(pathways_to_use, colnames(data_path_model))
  taxa_to_use     <- intersect(selected_taxa_table$taxa, colnames(data_taxa_raw))

  if (length(pathways_to_use) == 0) {
    warning("No selected pathways found in data_path_model for: ", scenario_name)
    return(tibble::tibble())
  }

  if (length(taxa_to_use) == 0) {
    warning("No selected taxa found in data_taxa_raw for: ", scenario_name)
    return(tibble::tibble())
  }

  data_taxa_model_fig5 <- data_taxa_raw

  data_taxa_model_fig5[taxa_to_use] <- lapply(
    data_taxa_model_fig5[taxa_to_use],
    function(x) log10(as.numeric(x) + pseudo)
  )

  cor_list <- list()

  for (p in pathways_to_use) {
    path_resid <- get_residuals_tbl_fig5(
      df          = data_path_model,
      feature     = p,
      confounders = confounders_func,
      min_n       = min_n
    ) %>%
      dplyr::rename(pathway_resid = resid)

    for (t in taxa_to_use) {
      taxa_resid <- get_residuals_tbl_fig5(
        df          = data_taxa_model_fig5,
        feature     = t,
        confounders = confounders_taxa,
        min_n       = min_n
      ) %>%
        dplyr::rename(taxa_resid = resid)

      tmp <- dplyr::inner_join(path_resid, taxa_resid, by = "SampleID")

      if (nrow(tmp) < min_n) next
      if (length(unique(tmp$pathway_resid)) < 3) next
      if (length(unique(tmp$taxa_resid)) < 3) next

      ct <- tryCatch(
        stats::cor.test(
          tmp$pathway_resid,
          tmp$taxa_resid,
          method = "spearman",
          exact  = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(ct)) next

      taxa_clean_now <- selected_taxa_table$taxa_clean[
        match(t, selected_taxa_table$taxa)
      ]

      ci <- spearman_ci(unname(ct$estimate), nrow(tmp))

      cor_list[[length(cor_list) + 1]] <- tibble::tibble(
        scenario     = scenario_name,
        pathway      = p,
        taxa         = t,
        taxa_clean   = taxa_clean_now,
        rho          = unname(ct$estimate),
        p_value      = ct$p.value,
        n            = nrow(tmp),
        rho_lower.CL = ci["lower"],
        rho_upper.CL = ci["upper"]
      )
    }
  }

  out <- dplyr::bind_rows(cor_list)

  if (nrow(out) == 0) {
    warning("No valid pathway–species correlations for: ", scenario_name)
    return(tibble::tibble())
  }

  out %>%
    dplyr::mutate(
      q_value = stats::p.adjust(p_value, method = "BH"),
      stars   = stars_from_q_fig5(q_value)
    ) %>%
    dplyr::arrange(q_value, p_value)
}
make_pathway_species_cor_heatmap_fig5 <- function(cor_df) {
  if (nrow(cor_df) == 0) {
    stop("Correlation table is empty; cannot build Panel C heatmap.")
  }

  cor_mat <- cor_df %>%
    dplyr::select(taxa_clean, pathway, rho) %>%
    dplyr::group_by(taxa_clean, pathway) %>%
    dplyr::summarise(
      rho = mean(rho, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = pathway,
      values_from = rho,
      values_fill = NA_real_
    ) %>%
    tibble::column_to_rownames("taxa_clean") %>%
    as.matrix()

  star_mat <- cor_df %>%
    dplyr::select(taxa_clean, pathway, stars) %>%
    dplyr::group_by(taxa_clean, pathway) %>%
    dplyr::summarise(
      stars = dplyr::case_when(
        any(stars == "***", na.rm = TRUE) ~ "***",
        any(stars == "**",  na.rm = TRUE) ~ "**",
        any(stars == "*",   na.rm = TRUE) ~ "*",
        TRUE ~ ""
      ),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = pathway,
      values_from = stars,
      values_fill = ""
    ) %>%
    tibble::column_to_rownames("taxa_clean") %>%
    as.matrix()

  star_mat <- star_mat[
    rownames(cor_mat),
    colnames(cor_mat),
    drop = FALSE
  ]

  pathway_lookup <- make_pathway_label_lookup_fig5(
    pathway_ids = colnames(cor_mat),
    path_anno   = path_anno,
    width       = 24
  )

  colnames(cor_mat)  <- pathway_lookup[colnames(cor_mat)]
  colnames(star_mat) <- pathway_lookup[colnames(star_mat)]

  col_fun_cor_fig5 <- circlize::colorRamp2(
    c(-0.8, 0, 0.8),
    c("#2C7BB6", "white", "#D7191C")
  )

  ComplexHeatmap::Heatmap(
    cor_mat,
    name = "Spearman\nrho",
    col = col_fun_cor_fig5,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "right",
    row_names_gp = grid::gpar(fontsize = 7.5, fontface = "italic"),
    column_names_gp = grid::gpar(fontsize = 7, lineheight = 0.9),
    column_names_max_height = grid::unit(7, "cm"),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title = "Spearman\nrho",
      at = c(-0.5, 0, 0.5)
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (!is.na(star_mat[i, j]) && star_mat[i, j] != "") {
        grid::grid.text(
          star_mat[i, j],
          x,
          y,
          gp = grid::gpar(fontsize = 8, fontface = "bold")
        )
      }
    }
  )
}
# Option A Panel C data and heatmap
cor_omnibus_fig5 <- run_pathway_species_cor_fig5(
  selected_species = omnibus_species_fig5,
  scenario_name    = "omnibus_group_effect",
  pathways_to_use  = top_pathways,
  min_n            = 10
)
fig5A_panelC_ht <- make_pathway_species_cor_heatmap_fig5(cor_omnibus_fig5)
fig5A_panelC_grob <- capture_heatmap_grob_fig5(fig5A_panelC_ht)
# Option B Panel C data and heatmap
cor_pairwise_fig5 <- run_pathway_species_cor_fig5(
  selected_species = pairwise_species_keep_fig5,
  scenario_name    = "pairwise_health_status",
  pathways_to_use  = top_pathways,
  min_n            = 10
)
fig5B_panelC_ht <- make_pathway_species_cor_heatmap_fig5(cor_pairwise_fig5)
fig5B_panelC_grob <- capture_heatmap_grob_fig5(fig5B_panelC_ht)
################################################################################
# 8. DRAW COMBINED FIGURE 5 OPTIONS
# Layout: A and B side-by-side on top, C below full width
################################################################################
draw_figure5_AB_top_C_bottom <- function(panelA_grob,
                                          panelB_grob,
                                          panelC_grob,
                                          figure_title = NULL) {
  grid::grid.newpage()

  if (is.null(figure_title) || !nzchar(figure_title)) {

    lay <- grid::grid.layout(
      nrow = 2,
      ncol = 2,
      heights = grid::unit(c(0.42, 0.58), "npc"),
      widths  = grid::unit(c(0.50, 0.50), "npc")
    )

    grid::pushViewport(grid::viewport(layout = lay))
    offset <- 0

  } else {

    lay <- grid::grid.layout(
      nrow = 3,
      ncol = 2,
      heights = grid::unit(c(0.055, 0.405, 0.540), "npc"),
      widths  = grid::unit(c(0.50, 0.50), "npc")
    )

    grid::pushViewport(grid::viewport(layout = lay))

    # Figure title
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = 1:2
      )
    )

    grid::grid.text(
      figure_title,
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.55, "npc"),
      gp = grid::gpar(
        fontsize   = 13,
        fontface   = "bold",
        fontfamily = fig_font
      )
    )

    grid::popViewport()
    offset <- 1
  }

  # Panel A: top-left
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 1 + offset,
      layout.pos.col = 1
    )
  )

  grid::grid.draw(panelA_grob)

  grid::grid.text(
    "A",
    x = grid::unit(0.010, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontsize   = 18,
      fontface   = "bold",
      fontfamily = fig_font
    )
  )

  grid::popViewport()

  # Panel B: top-right
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 1 + offset,
      layout.pos.col = 2
    )
  )

  grid::grid.draw(panelB_grob)

  grid::grid.text(
    "B",
    x = grid::unit(0.010, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontsize   = 18,
      fontface   = "bold",
      fontfamily = fig_font
    )
  )

  grid::popViewport()

  # Panel C: bottom, full width
  grid::pushViewport(
    grid::viewport(
      layout.pos.row = 2 + offset,
      layout.pos.col = 1:2
    )
  )

  grid::grid.draw(panelC_grob)

  grid::grid.text(
    "C",
    x = grid::unit(0.010, "npc"),
    y = grid::unit(0.985, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(
      fontsize   = 18,
      fontface   = "bold",
      fontfamily = fig_font
    )
  )

  grid::popViewport()
  grid::popViewport()
}
################################################################################
# 8.1 Preview Figure 5 Option A and Option B in RStudio plot window
################################################################################
draw_figure5_AB_top_C_bottom(
  panelA_grob  = fig5_panelA_grob,
  panelB_grob  = fig5A_panelB_grob,
  panelC_grob  = fig5A_panelC_grob,
  figure_title = "Figure 5 Option A. Functional and taxonomic structure underlying group differences"
)
draw_figure5_AB_top_C_bottom(
  panelA_grob  = fig5_panelA_grob,
  panelB_grob  = fig5B_panelB_grob,
  panelC_grob  = fig5B_panelC_grob,
  figure_title = "Figure 5 Option B. Functional and taxonomic structure underlying pairwise health-status comparisons"
)
################################################################################
# 8.2 Export Figure 5 Option A
################################################################################
open_tiff_fig5(
  filename  = fig5_optionA_file,
  width_mm  = fig5_width_mm,
  height_mm = fig5_height_mm,
  dpi       = fig5_dpi
)
draw_figure5_AB_top_C_bottom(
  panelA_grob  = fig5_panelA_grob,
  panelB_grob  = fig5A_panelB_grob,
  panelC_grob  = fig5A_panelC_grob,
  figure_title = "Figure 5 Option A. Functional and taxonomic structure underlying group differences"
)
grDevices::dev.off()
################################################################################
# 8.3 Export Figure 5 Option B
################################################################################
open_tiff_fig5(
  filename  = fig5_optionB_file,
  width_mm  = fig5_width_mm,
  height_mm = fig5_height_mm,
  dpi       = fig5_dpi
)
draw_figure5_AB_top_C_bottom(
  panelA_grob  = fig5_panelA_grob,
  panelB_grob  = fig5B_panelB_grob,
  panelC_grob  = fig5B_panelC_grob,
  figure_title = "Figure 5 Option B. Functional and taxonomic structure underlying pairwise health-status comparisons"
)
grDevices::dev.off()
################################################################################
# 9. OPTIONAL SUPPLEMENTARY EXCEL TABLES
################################################################################
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      "Fig5_Functional_PERMANOVA" = as.data.frame(global_permanova_func_fig5),
      "Fig5_Functional_PERMDISP_ANOVA" = as.data.frame(permdisp_anova_func_fig5),
      "Fig5A_OmnibusSpecies" = tibble::tibble(Species = omnibus_species_fig5),
      "Fig5A_PathwaySpeciesCorr" = cor_omnibus_fig5,
      "Fig5B_PairwiseSpeciesSig" = sig_species_pairwise_fig5,
      "Fig5B_PathwaySpeciesCorr" = cor_pairwise_fig5
    ),
    path = file.path(fig5_out_dir, "Supplementary_Figure5_OptionA_OptionB.xlsx")
  )
}
cat("\n[Figure 5 complete]\n")
cat("  Option A saved -> ", fig5_optionA_file, "\n", sep = "")
cat("  Option B saved -> ", fig5_optionB_file, "\n", sep = "")
cat("  Supplementary table saved -> ",
    file.path(fig5_out_dir, "Supplementary_Figure5_OptionA_OptionB.xlsx"),
    "\n", sep = "")
################################################################################
# END OF SCRIPT
################################################################################
