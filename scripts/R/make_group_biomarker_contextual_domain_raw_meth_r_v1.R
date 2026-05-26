suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)
set.seed(20260518)

env_path <- function(name, default = NULL, must_work = FALSE) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) {
    return(normalizePath(value, winslash = "/", mustWork = must_work))
  }
  if (!is.null(default)) {
    return(normalizePath(default, winslash = "/", mustWork = must_work))
  }
  stop("Set environment variable ", name, " (see config/paths_template.yml).")
}

biomarker_root <- env_path("EOBC_BIOMARKER_ROOT")
base_dir <- env_path("EOBC_FINAL_ANALYSIS_DIR", file.path(biomarker_root, "final_analysis"))
landscape_dir <- file.path(base_dir, "group_biomarker_landscape_r_v1", "tables")
domain_dir <- file.path(base_dir, "group_marker_domain_specific_r_v1", "tables")
out_dir <- file.path(base_dir, "group_biomarker_contextual_domain_raw_meth_r_v1")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

group_order <- c("G1 | H", "G2 | I", "G3 | L-like", "G4 | L")
modality_order <- c("RNA expression", "TSS methylation")

group_cols <- c(
  "G1 | H" = "#D11224",
  "G2 | I" = "#FF9416",
  "G3 | L-like" = "#8BBE69",
  "G4 | L" = "#1B9278"
)

family_cols <- c(
  "Immune" = "#4EA5F0",
  "Repair" = "#47C56B",
  "Glycolysis / TCA" = "#E6BC18",
  "Fatty acid" = "#F4A259",
  "Kinase signaling" = "#9B6AE8",
  "Hormone signaling" = "#9AAABC"
)

domain_cols <- c(
  "OS" = "#58A9F9",
  "DepMap" = "#FF8A1C",
  "Immune" = "#2EC7D3"
)

os_cols <- c(
  "Protective" = "#58A9F9",
  "Adverse" = "#F4A259",
  "Not significant" = "#AEB8C5"
)

drug_cols <- c(
  "Sensitive" = "#58A9F9",
  "Resistant" = "#FF8A1C",
  "Not significant" = "#AEB8C5"
)

immune_cols <- c(
  "TIL rho positive" = "#17A2A4",
  "TIL rho negative" = "#72A950",
  "TMB rho positive" = "#7B61D1",
  "TMB rho negative" = "#6F8FB7",
  "TIL/TMB rho positive" = "#0E9F6E",
  "TIL/TMB rho negative" = "#4B8B3B",
  "Discordant TIL/TMB rho" = "#B9C4D1",
  "Not significant" = "#AEB8C5"
)

group_shapes <- c(
  "G1 | H" = 21,
  "G2 | I" = 22,
  "G3 | L-like" = 24,
  "G4 | L" = 23
)

theme_eobc <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.28), color = "#111827"),
      plot.subtitle = element_text(color = "#64748B", margin = margin(b = 8)),
      plot.caption = element_text(color = "#64748B", hjust = 0, size = rel(0.72)),
      axis.title = element_text(face = "bold", color = "#111827"),
      axis.text = element_text(color = "#64748B"),
      panel.grid.major = element_line(color = "#DCE6F2", linewidth = 0.42),
      panel.grid.minor = element_line(color = "#ECF2F8", linewidth = 0.22),
      strip.background = element_rect(fill = "#EAF1F8", color = "#D8E3F0"),
      strip.text = element_text(face = "bold", color = "#111827"),
      legend.title = element_text(face = "bold", color = "#111827"),
      legend.text = element_text(color = "#64748B"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(14, 16, 12, 16)
    )
}

save_plot <- function(plot, name, width, height, dpi = 360) {
  ggsave(file.path(plot_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
}

safe_neglog <- function(p, cap = 12) pmin(-log10(pmax(p, 1e-300)), cap)

safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 8 || length(unique(x)) < 3 || length(unique(y)) < 3) {
    return(tibble(rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(rho = unname(ct$estimate), p = ct$p.value)
}

safe_kw <- function(value, group) {
  keep <- is.finite(value) & !is.na(group)
  value <- value[keep]
  group <- droplevels(factor(group[keep]))
  if (length(unique(group)) < 2 || length(value) < 8) return(NA_real_)
  out <- try(kruskal.test(value ~ group)$p.value, silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}

immune_score_status <- function(til_rho, til_fdr, tmb_rho, tmb_fdr, fdr_cutoff = 0.10) {
  til_sig <- is.finite(til_rho) & is.finite(til_fdr) & til_fdr < fdr_cutoff
  tmb_sig <- is.finite(tmb_rho) & is.finite(tmb_fdr) & tmb_fdr < fdr_cutoff
  case_when(
    til_sig & tmb_sig & til_rho >= 0 & tmb_rho >= 0 ~ "TIL/TMB rho positive",
    til_sig & tmb_sig & til_rho < 0 & tmb_rho < 0 ~ "TIL/TMB rho negative",
    til_sig & tmb_sig ~ "Discordant TIL/TMB rho",
    til_sig & til_rho >= 0 ~ "TIL rho positive",
    til_sig & til_rho < 0 ~ "TIL rho negative",
    tmb_sig & tmb_rho >= 0 ~ "TMB rho positive",
    tmb_sig & tmb_rho < 0 ~ "TMB rho negative",
    TRUE ~ "Not significant"
  )
}

read_domain <- function(file) {
  read_csv(file.path(domain_dir, file), show_col_types = FALSE) %>%
    mutate(
      activity_dominant_group = factor(activity_dominant_group, levels = group_order),
      display_layer = if_else(modality == "METH", "TSS methylation", "RNA expression"),
      display_layer = factor(display_layer, levels = modality_order),
      Family6 = factor(Family6, levels = names(family_cols)),
      gene_layer = paste0(gene_label, if_else(modality == "METH", " | TSS methylation", " | RNA")),
      program = paste0(activity_dominant_group, "\n", display_layer)
    )
}

biomarker_long <- read_csv(file.path(landscape_dir, "group_biomarker_long_values.csv"), show_col_types = FALSE) %>%
  mutate(
    group_label = factor(group_label, levels = group_order),
    display_layer = if_else(modality == "METH", "TSS methylation", "RNA expression"),
    display_layer = factor(display_layer, levels = modality_order),
    activity_z = value_z,
    activity_z_clip = pmax(pmin(activity_z, 2.25), -2.25),
    layer_short = if_else(modality == "METH", "M", "R"),
    gene_label = paste0(gene, " [", layer_short, "]"),
    Family6 = factor(Family6, levels = names(family_cols))
  )

immune_metric_path <- file.path(dirname(base_dir), "10_external_validation_TIL_TMB_DepMap_final_v7", "tables", "immune_signature_tcga_til_tmb_merged.csv")
immune_metrics <- read_csv(immune_metric_path, show_col_types = FALSE) %>%
  transmute(
    Sample,
    TIL_score = as.numeric(TIL_score),
    TMB_log1p = as.numeric(tmb_log1p)
  )

tmb_outlier_samples <- immune_metrics %>%
  filter(is.finite(TMB_log1p), TMB_log1p <= 0) %>%
  distinct(Sample) %>%
  arrange(Sample) %>%
  mutate(
    exclusion_scope = "TMB only",
    reason = "TMB_log1p <= 0"
  )
write_csv(tmb_outlier_samples, file.path(table_dir, "tmb_outlier_samples_excluded_tmb_only.csv"))

immune_recalc <- biomarker_long %>%
  select(Sample, gene, modality, activity_z) %>%
  distinct() %>%
  inner_join(immune_metrics, by = "Sample") %>%
  filter(is.finite(activity_z)) %>%
  group_by(gene, modality) %>%
  group_modify(~ {
    til <- safe_cor(.x$activity_z, .x$TIL_score)
    tmb_df <- .x %>% filter(!Sample %in% tmb_outlier_samples$Sample)
    tmb <- safe_cor(tmb_df$activity_z, tmb_df$TMB_log1p)
    tibble(
      rho_TIL_raw = til$rho[1],
      p_TIL_raw = til$p[1],
      rho_TMB_raw = tmb$rho[1],
      p_TMB_raw = tmb$p[1],
      n_TIL = sum(is.finite(.x$activity_z) & is.finite(.x$TIL_score)),
      n_TMB = sum(is.finite(tmb_df$activity_z) & is.finite(tmb_df$TMB_log1p))
    )
  }) %>%
  ungroup() %>%
  mutate(
    til_activity_fdr = p.adjust(p_TIL_raw, method = "BH"),
    tmb_activity_fdr = p.adjust(p_TMB_raw, method = "BH"),
    til_activity_rho = if_else(modality == "METH", -rho_TIL_raw, rho_TIL_raw),
    tmb_activity_rho = if_else(modality == "METH", -rho_TMB_raw, rho_TMB_raw),
    immune_activity_strength = sqrt(til_activity_rho^2 + tmb_activity_rho^2),
    immune_activity_quadrant = case_when(
      til_activity_rho >= 0 & tmb_activity_rho >= 0 ~ "Immune-hot / TMB-high",
      til_activity_rho < 0 & tmb_activity_rho < 0 ~ "Immune-cold / TMB-low",
      til_activity_rho >= 0 & tmb_activity_rho < 0 ~ "TIL-high / TMB-low",
      til_activity_rho < 0 & tmb_activity_rho >= 0 ~ "TMB-shifted",
      TRUE ~ NA_character_
    )
  )
write_csv(immune_recalc, file.path(table_dir, "immune_til_tmb_recalculated_tmb_outliers_removed.csv"))

replace_immune_activity <- function(df) {
  df %>%
    select(-any_of(c(
      "til_activity_rho", "tmb_activity_rho", "til_activity_fdr", "tmb_activity_fdr",
      "immune_activity_strength", "immune_activity_quadrant"
    ))) %>%
    left_join(
      immune_recalc %>%
        select(
          gene, modality, til_activity_rho, tmb_activity_rho, til_activity_fdr, tmb_activity_fdr,
          immune_activity_strength, immune_activity_quadrant, n_TIL, n_TMB
        ),
      by = c("gene", "modality")
    )
}

domain_marker <- read_domain("group_defined_marker_domain_input_tss_low_activity.csv") %>%
  replace_immune_activity() %>%
  mutate(
    os_sig = os_sig %in% TRUE,
    drug_sig = drug_sig %in% TRUE & !str_detect(coalesce(drug_class_clean, ""), regex("Other|exploratory", ignore_case = TRUE)),
    til_sig = is.finite(til_activity_rho) & is.finite(til_activity_fdr) & til_activity_fdr < 0.10,
    tmb_sig = is.finite(tmb_activity_rho) & is.finite(tmb_activity_fdr) & tmb_activity_fdr < 0.10,
    immune_sig = immune_sig %in% TRUE & (til_sig | tmb_sig),
    os_status = case_when(os_sig & os_activity_direction == "Protective" ~ "Protective",
                          os_sig & os_activity_direction == "Adverse" ~ "Adverse",
                          TRUE ~ "Not significant"),
    drug_status = case_when(drug_sig & depmap_activity_direction == "Sensitive" ~ "Sensitive",
                            drug_sig & depmap_activity_direction == "Resistant" ~ "Resistant",
                            TRUE ~ "Not significant"),
    immune_status = if_else(immune_sig, immune_score_status(til_activity_rho, til_activity_fdr, tmb_activity_rho, tmb_activity_fdr), "Not significant")
  )

os_endpoint <- read_domain("OS_group_defined_marker_endpoint_results.csv") %>%
  mutate(
    endpoint_label = factor(endpoint_label, levels = c("Overall OS", "5-year OS", "10-year OS")),
    os_effect = case_when(activity_beta < 0 ~ "Protective", activity_beta > 0 ~ "Adverse", TRUE ~ "Neutral"),
    os_evidence = selected %in% TRUE | p < 0.10
  )

depmap_best <- read_domain("DepMap_group_defined_marker_best_drug_results.csv") %>%
  mutate(
    drug_sig = depmap_fdr < 0.25 & feature_available %in% TRUE &
      !str_detect(coalesce(drug_class_clean, ""), regex("Other|exploratory", ignore_case = TRUE)),
    drug_status = case_when(drug_sig & depmap_activity_direction == "Sensitive" ~ "Sensitive",
                            drug_sig & depmap_activity_direction == "Resistant" ~ "Resistant",
                            TRUE ~ "Not significant"),
    drug_label = if_else(feature_available %in% TRUE, paste0(gene_label, "\n", drug_clean), paste0(gene_label, "\nnot tested"))
  )

immune_best <- read_domain("Immune_group_defined_marker_til_tmb_results.csv") %>%
  replace_immune_activity() %>%
  mutate(
    immune_best_fdr = pmin(til_activity_fdr, tmb_activity_fdr, na.rm = TRUE),
    til_sig = is.finite(til_activity_rho) & is.finite(til_activity_fdr) & til_activity_fdr < 0.10,
    tmb_sig = is.finite(tmb_activity_rho) & is.finite(tmb_activity_fdr) & tmb_activity_fdr < 0.10,
    immune_sig = til_sig | tmb_sig,
    immune_status = if_else(immune_sig, immune_score_status(til_activity_rho, til_activity_fdr, tmb_activity_rho, tmb_activity_fdr), "Not significant")
  )

# Recompute dominance on the displayed feature scale:
# RNA rows are expression; TSS methylation rows use raw beta-value z-scores.
activity_group_means <- biomarker_long %>%
  group_by(modality, display_layer, gene, gene_label, Layer, Family6, target_label, group_label) %>%
  summarise(mean_activity_z = mean(activity_z, na.rm = TRUE), .groups = "drop")

activity_dominance <- activity_group_means %>%
  group_by(modality, display_layer, gene, gene_label, Layer, Family6, target_label) %>%
  summarise(
    activity_dominant_group = as.character(group_label[which.max(mean_activity_z)]),
    activity_suppressed_group = as.character(group_label[which.min(mean_activity_z)]),
    dominant_mean_z = max(mean_activity_z, na.rm = TRUE),
    suppressed_mean_z = min(mean_activity_z, na.rm = TRUE),
    activity_contrast_z = dominant_mean_z - suppressed_mean_z,
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    kw_p = safe_kw(
      biomarker_long$activity_z[biomarker_long$modality == modality & biomarker_long$gene == gene],
      biomarker_long$group_label[biomarker_long$modality == modality & biomarker_long$gene == gene]
    )
  ) %>%
  ungroup() %>%
  group_by(modality) %>%
  mutate(kw_fdr = p.adjust(kw_p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    activity_dominant_group = factor(activity_dominant_group, levels = group_order),
    display_layer = factor(display_layer, levels = modality_order),
    neg_log10_fdr = safe_neglog(kw_fdr, 12),
    state_family_order = paste(display_layer, activity_dominant_group, Family6),
    rank_score = activity_contrast_z * neg_log10_fdr
  ) %>%
  left_join(
    domain_marker %>% select(
      gene, modality, os_sig, drug_sig, immune_sig, os_status, drug_status, immune_status,
      os_best_endpoint, os_activity_direction, drug_clean, drug_class_clean,
      depmap_activity_direction, immune_activity_quadrant,
      til_activity_rho, tmb_activity_rho, til_activity_fdr, tmb_activity_fdr, til_sig, tmb_sig
    ),
    by = c("gene", "modality")
  ) %>%
  mutate(
    os_sig = os_sig %in% TRUE,
    drug_sig = drug_sig %in% TRUE,
    immune_sig = immune_sig %in% TRUE,
    os_status = replace_na(os_status, "Not significant"),
    drug_status = replace_na(drug_status, "Not significant"),
    immune_status = replace_na(immune_status, "Not significant"),
    domain_count = os_sig + drug_sig + immune_sig
  )

write_csv(activity_dominance, file.path(table_dir, "activity_scale_group_dominance_with_domain_evidence.csv"))

rna_meth_alignment <- read_csv(file.path(landscape_dir, "rna_meth_tss_inverse_alignment.csv"), show_col_types = FALSE)
write_csv(rna_meth_alignment, file.path(table_dir, "rna_meth_alignment_reused_for_contextual_figures.csv"))

# ---------------------------------------------------------------------------
# Figure 01A. Publication-style sample heatmap with RNA-METH correlation strip.
# ---------------------------------------------------------------------------
sample_order <- biomarker_long %>%
  distinct(Sample, group_label, Pseudotime, PAM50, Cluster) %>%
  arrange(group_label, Pseudotime, Sample) %>%
  mutate(sample_x = row_number())

row_order <- activity_dominance %>%
  arrange(display_layer, activity_dominant_group, Family6, desc(rank_score), gene) %>%
  mutate(
    feature_id = paste(modality, gene, sep = "___"),
    feature_label = gene
  )

feature_levels <- rev(row_order$feature_id)
feature_labels <- setNames(row_order$feature_label, row_order$feature_id)
n_samples <- nrow(sample_order)

heat_df <- biomarker_long %>%
  inner_join(sample_order %>% select(Sample, sample_x), by = "Sample") %>%
  mutate(feature_id = paste(modality, gene, sep = "___")) %>%
  filter(feature_id %in% feature_levels) %>%
  mutate(feature_id = factor(feature_id, levels = feature_levels))

cor_strip <- row_order %>%
  select(feature_id, gene, modality) %>%
  left_join(rna_meth_alignment %>% select(gene, sample_spearman_rho, group_spearman_rho, inverse_consistency_class), by = "gene") %>%
  pivot_longer(c(sample_spearman_rho, group_spearman_rho), names_to = "cor_type", values_to = "rho") %>%
  mutate(
    feature_id = factor(feature_id, levels = feature_levels),
    cor_type = recode(cor_type, sample_spearman_rho = "sample rho", group_spearman_rho = "group rho"),
    x = n_samples + recode(cor_type, "sample rho" = 4.5, "group rho" = 6.25)
  )

domain_strip <- row_order %>%
  select(feature_id, gene, modality) %>%
  left_join(domain_marker %>% select(gene, modality, os_sig, drug_sig, immune_sig), by = c("gene", "modality")) %>%
  mutate(across(c(os_sig, drug_sig, immune_sig), ~ .x %in% TRUE)) %>%
  pivot_longer(c(os_sig, drug_sig, immune_sig), names_to = "domain", values_to = "present") %>%
  mutate(
    domain = recode(domain, os_sig = "OS", drug_sig = "DepMap", immune_sig = "Immune"),
    feature_id = factor(feature_id, levels = feature_levels),
    x = n_samples + recode(domain, "OS" = 9.10, "DepMap" = 11.05, "Immune" = 13.00),
    domain_fill = if_else(present, domain, "none")
  )

sample_track <- sample_order %>%
  select(sample_x, group_label, Pseudotime) %>%
  mutate(pt_scaled = rescale(Pseudotime, to = c(1.4, 2.15), from = range(Pseudotime, na.rm = TRUE)))

p_track <- ggplot() +
  geom_tile(data = sample_track, aes(sample_x, 1, fill = group_label), height = 0.55, width = 1) +
  geom_line(data = sample_track, aes(sample_x, pt_scaled), color = "#334155", linewidth = 0.28, alpha = 0.75) +
  geom_point(data = sample_track, aes(sample_x, pt_scaled, color = group_label), size = 0.45, alpha = 0.8) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_color_manual(values = group_cols, guide = "none") +
  scale_x_continuous(limits = c(0.5, n_samples + 15.4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.55, 2.35), breaks = c(1, 1.85), labels = c("EOBC", "PT"), expand = c(0, 0)) +
  theme_void(base_family = "Arial") +
  theme(
    axis.text.y = element_text(color = "#64748B", size = 7, face = "bold", hjust = 1),
    plot.margin = margin(0, 14, 0, 14),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_heat <- ggplot() +
  geom_tile(data = heat_df, aes(sample_x, feature_id, fill = activity_z_clip), width = 1, height = 0.92) +
  scale_fill_gradient2(
    low = "#3265A8", mid = "#F8FAFC", high = "#C92A2A",
    midpoint = 0, limits = c(-2.25, 2.25), oob = squish,
    name = "Feature\nz-score"
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = cor_strip, aes(x, feature_id, fill = rho), width = 1.35, height = 0.92, color = "white", linewidth = 0.22) +
  scale_fill_gradient2(
    low = "#3265A8", mid = "#F8FAFC", high = "#C92A2A",
    midpoint = 0, limits = c(-1, 1), oob = squish,
    name = "RNA-METH\nSpearman rho"
  ) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = domain_strip, aes(x, feature_id, fill = domain_fill), width = 1.55, height = 0.80, color = "#CBD5E1", linewidth = 0.22) +
  scale_fill_manual(values = c(domain_cols, "none" = "#F8FAFC"), breaks = names(domain_cols), name = "Domain\nevidence") +
  geom_vline(xintercept = n_samples + 2.45, color = "#111827", linewidth = 0.25) +
  geom_vline(xintercept = n_samples + 7.62, color = "#111827", linewidth = 0.25) +
  facet_grid(display_layer ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(
    limits = c(0.5, n_samples + 15.4),
    breaks = n_samples + c(4.5, 6.25, 9.10, 11.05, 13.00),
    labels = c("sample\nrho", "group\nrho", "OS", "Drug", "Imm"),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = feature_labels) +
  labs(
    title = "A. EOBC biomarker expression and TSS methylation heatmap with RNA-TSS coupling",
    subtitle = "Columns are tumors ordered by EOBC group and pseudotime. RNA rows show expression z-scores; methylation rows show raw TSS beta-value z-scores.",
    x = NULL,
    y = NULL,
    caption = "RNA rows show expression z-scores and TSS methylation rows show raw beta-value z-scores. sample rho is Spearman correlation across matched tumors; group rho is Spearman correlation across the four EOBC group means. Negative RNA-METH rho supports promoter methylation-linked transcriptional repression."
  ) +
  theme_eobc(base_size = 8.2) +
  theme(
    axis.text.x = element_text(size = 7, angle = 0, vjust = 1),
    axis.text.y = element_text(size = 6.6),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    strip.text.y = element_text(angle = 0),
    plot.margin = margin(8, 18, 8, 18)
  )

fig_heat <- p_track / p_heat + plot_layout(heights = c(0.42, 7.2))
save_plot(fig_heat, "Figure_01A_activity_heatmap_with_rna_meth_correlation_R_v1", 13.2, 9.4)

# The faceted version above is useful as a draft, but free-y facets can leave
# excessive whitespace when a full factor order is shared. Rebuild the final
# heatmap as two explicitly stacked panels with aligned x-axes.
make_activity_heat_panel <- function(layer_name, show_x = FALSE, show_legend = TRUE) {
  layer_features <- row_order %>%
    filter(display_layer == layer_name) %>%
    arrange(activity_dominant_group, Family6, desc(rank_score), gene) %>%
    pull(feature_id)
  layer_levels <- rev(layer_features)

  h <- heat_df %>%
    filter(display_layer == layer_name) %>%
    mutate(feature_id = factor(as.character(feature_id), levels = layer_levels))
  cdat <- cor_strip %>%
    filter(as.character(feature_id) %in% layer_features) %>%
    mutate(feature_id = factor(as.character(feature_id), levels = layer_levels))
  ddat <- domain_strip %>%
    filter(as.character(feature_id) %in% layer_features) %>%
    mutate(feature_id = factor(as.character(feature_id), levels = layer_levels))

  ggplot() +
    geom_tile(data = h, aes(sample_x, feature_id, fill = activity_z_clip), width = 1, height = 0.92) +
    scale_fill_gradient2(
      low = "#3265A8", mid = "#F8FAFC", high = "#C92A2A",
      midpoint = 0, limits = c(-2.25, 2.25), oob = squish,
      name = "Feature\nz-score"
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = cdat, aes(x, feature_id, fill = rho), width = 1.35, height = 0.92, color = "white", linewidth = 0.22) +
    scale_fill_gradient2(
      low = "#3265A8", mid = "#F8FAFC", high = "#C92A2A",
      midpoint = 0, limits = c(-1, 1), oob = squish,
      name = "RNA-METH\nSpearman rho"
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = ddat, aes(x, feature_id, fill = domain_fill), width = 1.55, height = 0.80, color = "#CBD5E1", linewidth = 0.22) +
    scale_fill_manual(values = c(domain_cols, "none" = "#F8FAFC"), breaks = names(domain_cols), name = "Domain\nevidence") +
    geom_vline(xintercept = n_samples + 2.45, color = "#111827", linewidth = 0.24) +
    geom_vline(xintercept = n_samples + 7.62, color = "#111827", linewidth = 0.24) +
    scale_x_continuous(
      limits = c(0.5, n_samples + 15.4),
      breaks = if (show_x) n_samples + c(4.5, 6.25, 9.10, 11.05, 13.00) else NULL,
      labels = if (show_x) c("sample\nrho", "group\nrho", "OS", "Drug", "Imm") else NULL,
      expand = c(0, 0)
    ) +
    scale_y_discrete(labels = feature_labels[layer_levels]) +
    labs(x = NULL, y = layer_name) +
    theme_eobc(base_size = 8.0) +
    theme(
      axis.title.y = element_text(face = "bold", size = 9.5, color = "#111827", margin = margin(r = 8)),
      axis.text.y = element_text(size = 6.45),
      axis.text.x = if (show_x) element_text(size = 6.2, angle = 45, hjust = 1, vjust = 1) else element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "horizontal",
      plot.margin = margin(2, 18, 2, 18)
    )
}

p_heat_rna <- make_activity_heat_panel("RNA expression", show_x = FALSE, show_legend = FALSE)
p_heat_tss <- make_activity_heat_panel("TSS methylation", show_x = TRUE, show_legend = TRUE)

fig_heat_refined <- p_track / p_heat_rna / p_heat_tss +
  plot_layout(heights = c(0.42, 3.3, 3.3)) +
  plot_annotation(
    title = "A. EOBC biomarker expression and TSS methylation heatmap with RNA-TSS coupling",
    subtitle = "Columns are tumors ordered by EOBC group and pseudotime. RNA rows show expression z-scores; methylation rows show raw TSS beta-value z-scores.",
    caption = "M denotes raw TSS methylation beta-value z-score. sample rho is Spearman correlation across matched tumors; group rho is Spearman correlation across the four EOBC group means. Negative RNA-METH rho supports promoter methylation-linked transcriptional repression.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, color = "#111827", margin = margin(b = 5)),
      plot.subtitle = element_text(size = 11, color = "#64748B", margin = margin(b = 8)),
      plot.caption = element_text(size = 7.5, color = "#64748B", hjust = 0, margin = margin(t = 6)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 10, 8, 10)
    )
  )
save_plot(fig_heat_refined, "Figure_01A_activity_heatmap_with_rna_meth_correlation_R_v1", 15.6, 8.9)

# ---------------------------------------------------------------------------
# Figure 01B. Full RNA-vs-METH gene-pair correlation matrix.
# ---------------------------------------------------------------------------
rna_wide <- biomarker_long %>%
  filter(modality == "RNA") %>%
  select(Sample, gene, value_z) %>%
  pivot_wider(names_from = gene, values_from = value_z)
meth_wide <- biomarker_long %>%
  filter(modality == "METH") %>%
  select(Sample, gene, value_z) %>%
  pivot_wider(names_from = gene, values_from = value_z)

common_samples <- intersect(rna_wide$Sample, meth_wide$Sample)
rna_mat <- rna_wide %>% filter(Sample %in% common_samples) %>% arrange(Sample)
meth_mat <- meth_wide %>% filter(Sample %in% common_samples) %>% arrange(Sample)
genes_rna <- setdiff(names(rna_mat), "Sample")
genes_meth <- setdiff(names(meth_mat), "Sample")

cor_matrix <- expand_grid(rna_gene = genes_rna, meth_gene = genes_meth) %>%
  rowwise() %>%
  mutate(
    rho = safe_cor(rna_mat[[rna_gene]], meth_mat[[meth_gene]])$rho,
    p = safe_cor(rna_mat[[rna_gene]], meth_mat[[meth_gene]])$p
  ) %>%
  ungroup() %>%
  group_by(rna_gene) %>%
  mutate(row_fdr = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  left_join(activity_dominance %>% filter(modality == "RNA") %>% select(gene, rna_family = Family6, rna_group = activity_dominant_group, rna_rank = rank_score), by = c("rna_gene" = "gene")) %>%
  left_join(activity_dominance %>% filter(modality == "METH") %>% select(gene, meth_family = Family6, meth_group = activity_dominant_group, meth_rank = rank_score), by = c("meth_gene" = "gene"))

gene_order_rna <- activity_dominance %>%
  filter(modality == "RNA") %>%
  arrange(activity_dominant_group, Family6, desc(rank_score), gene) %>%
  pull(gene)
gene_order_meth <- activity_dominance %>%
  filter(modality == "METH") %>%
  arrange(activity_dominant_group, Family6, desc(rank_score), gene) %>%
  pull(gene)

cor_matrix <- cor_matrix %>%
  mutate(
    rna_gene = factor(rna_gene, levels = rev(gene_order_rna)),
    meth_gene = factor(meth_gene, levels = gene_order_meth),
    diagonal = as.character(rna_gene) == as.character(meth_gene)
  )

write_csv(cor_matrix, file.path(table_dir, "rna_by_tss_methylation_gene_pair_correlation_matrix.csv"))

p_cor_matrix <- ggplot(cor_matrix, aes(meth_gene, rna_gene)) +
  geom_tile(aes(fill = rho), color = "white", linewidth = 0.18) +
  geom_tile(data = cor_matrix %>% filter(diagonal), fill = NA, color = "#111827", linewidth = 0.42) +
  geom_text(
    data = cor_matrix %>% filter(diagonal),
    aes(label = sprintf("%.2f", rho)),
    size = 2.15, color = "#111827"
  ) +
  scale_fill_gradient2(
    low = "#3265A8", mid = "#F8FAFC", high = "#C92A2A",
    midpoint = 0, limits = c(-0.8, 0.8), oob = squish,
    name = "Spearman rho\nRNA vs TSS methylation"
  ) +
  labs(
    title = "B. RNA-expression versus TSS-methylation coupling across biomarker genes",
    subtitle = "Rows are RNA features and columns are TSS methylation features. Diagonal cells test the matched gene pair; blue diagonal values support inverse promoter methylation biology.",
    x = "TSS methylation feature",
    y = "RNA feature"
  ) +
  theme_eobc(base_size = 8.6) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 6.4),
    axis.text.y = element_text(size = 6.4),
    legend.position = "right"
  )
save_plot(p_cor_matrix, "Figure_01B_RNA_TSS_methylation_correlation_matrix_R_v1", 8.6, 7.5)

# ---------------------------------------------------------------------------
# Domain overlays on the group-dominance scatter grammar.
# ---------------------------------------------------------------------------
label_top_domain <- function(df, flag_col, n_per_panel = 8) {
  flag <- rlang::ensym(flag_col)
  df %>%
    mutate(.flag = !!flag) %>%
    filter(.flag %in% TRUE | (kw_fdr < 0.05 & activity_contrast_z >= 0.75)) %>%
    group_by(display_layer, activity_dominant_group) %>%
    arrange(desc(.flag), desc(rank_score), .by_group = TRUE) %>%
    slice_head(n = n_per_panel) %>%
    ungroup()
}

base_dominance_plot <- function(df, title, subtitle, stroke_aes, stroke_scale, label_data, label_aes, name) {
  p <- ggplot(df, aes(activity_contrast_z, neg_log10_fdr)) +
    geom_hline(yintercept = -log10(0.10), color = "#FFB26B", linewidth = 0.45, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), color = "#FF6B6B", linewidth = 0.55, linetype = "longdash") +
    geom_point(
      aes(fill = Family6, shape = activity_dominant_group, size = pmax(abs(dominant_mean_z), 0.12), color = {{ stroke_aes }}),
      stroke = 1.1, alpha = 0.92
    ) +
    geom_text_repel(
      data = label_data,
      mapping = label_aes,
      size = 2.55, fontface = "bold",
      min.segment.length = 0,
      segment.color = alpha("#64748B", 0.45),
      box.padding = 0.28,
      point.padding = 0.16,
      max.overlaps = Inf,
      seed = 20260518,
      show.legend = FALSE
    ) +
    facet_wrap(~ display_layer, nrow = 1) +
    scale_fill_manual(values = family_cols, name = "Biomarker family", drop = FALSE) +
    scale_shape_manual(values = group_shapes, name = "Highest EOBC group", drop = FALSE) +
    scale_size_continuous(range = c(2.2, 7.8), name = "|Dominant\nmean z|") +
    stroke_scale +
    coord_cartesian(ylim = c(0, 12.6), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Group-defining contrast (highest mean feature z - lowest mean feature z)",
      y = "-log10(Kruskal-Wallis FDR), capped at 12"
    ) +
    theme_eobc(base_size = 9.6) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  save_plot(p, name, 12.2, 6.8)
  p
}

os_label_data <- activity_dominance %>%
  mutate(label_text = if_else(os_sig, paste0(gene_label, "\n", os_best_endpoint, " / ", os_status), gene_label)) %>%
  label_top_domain(os_sig)

p_os_domain <- base_dominance_plot(
  activity_dominance,
  "C. OS evidence projected onto EOBC group-defining biomarkers",
  "The same group-dominance landscape is reused; point outline shows whether a group marker is linked to protective or adverse OS.",
  os_status,
  scale_color_manual(values = os_cols, name = "OS evidence", drop = FALSE),
  os_label_data,
  aes(label = label_text, color = os_status),
  "Figure_02A_OS_evidence_on_group_dominance_R_v1"
)

depmap_label_data <- activity_dominance %>%
  mutate(label_text = if_else(drug_sig, paste0(gene_label, "\n", drug_clean, " / ", drug_status), gene_label)) %>%
  label_top_domain(drug_sig)

p_depmap_domain <- base_dominance_plot(
  activity_dominance,
  "D. DepMap drug-response evidence projected onto EOBC group-defining biomarkers",
  "Point outline summarizes whether high biomarker activity marks more drug sensitivity or resistance; labels include the best associated drug.",
  drug_status,
  scale_color_manual(values = drug_cols, name = "DepMap evidence", drop = FALSE),
  depmap_label_data,
  aes(label = label_text, color = drug_status),
  "Figure_03A_DepMap_evidence_on_group_dominance_R_v1"
)

immune_label_data <- activity_dominance %>%
  mutate(label_text = if_else(immune_sig, paste0(gene_label, "\n", immune_status), gene_label)) %>%
  label_top_domain(immune_sig)

p_immune_domain <- base_dominance_plot(
  activity_dominance,
  "E. Immune TIL/TMB evidence projected onto EOBC group-defining biomarkers",
  "Only biomarkers with significant continuous TIL-score or TMB-log1p Spearman association are highlighted; labels report the associated score direction.",
  immune_status,
  scale_color_manual(values = immune_cols, name = "Immune evidence", drop = FALSE),
  immune_label_data,
  aes(label = label_text, color = immune_status),
  "Figure_04A_Immune_evidence_on_group_dominance_R_v1"
)

# ---------------------------------------------------------------------------
# Family trajectories with one domain at a time.
# ---------------------------------------------------------------------------
family_activity <- activity_group_means %>%
  group_by(display_layer, Family6, group_label) %>%
  summarise(family_mean_z = mean(mean_activity_z, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    display_layer = factor(display_layer, levels = modality_order),
    Family6 = factor(Family6, levels = names(family_cols)),
    group_label = factor(group_label, levels = group_order)
  )

domain_counts <- domain_marker %>%
  mutate(
    display_layer = factor(display_layer, levels = modality_order),
    Family6 = factor(Family6, levels = names(family_cols)),
    activity_dominant_group = factor(activity_dominant_group, levels = group_order)
  ) %>%
  group_by(display_layer, Family6, activity_dominant_group) %>%
  summarise(
    os_count = sum(os_sig, na.rm = TRUE),
    os_protective = sum(os_sig & os_status == "Protective", na.rm = TRUE),
    os_adverse = sum(os_sig & os_status == "Adverse", na.rm = TRUE),
    drug_count = sum(drug_sig, na.rm = TRUE),
    sensitive_count = sum(drug_sig & drug_status == "Sensitive", na.rm = TRUE),
    resistant_count = sum(drug_sig & drug_status == "Resistant", na.rm = TRUE),
    immune_count = sum(immune_sig, na.rm = TRUE),
    immune_pos_count = sum(immune_sig & str_detect(immune_status, "positive"), na.rm = TRUE),
    immune_neg_count = sum(immune_sig & str_detect(immune_status, "negative"), na.rm = TRUE),
    immune_til_count = sum(immune_sig & str_detect(immune_status, "TIL"), na.rm = TRUE),
    immune_tmb_count = sum(immune_sig & str_detect(immune_status, "TMB"), na.rm = TRUE),
    immune_discordant_count = sum(immune_sig & immune_status == "Discordant TIL/TMB rho", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(group_label = activity_dominant_group)

trajectory_domain_df <- family_activity %>%
  left_join(domain_counts, by = c("display_layer", "Family6", "group_label")) %>%
  mutate(across(ends_with("_count") | c(os_count, drug_count, immune_count, os_protective, os_adverse, sensitive_count, resistant_count, immune_pos_count, immune_neg_count, immune_til_count, immune_tmb_count, immune_discordant_count), ~ replace_na(.x, 0)))

plot_domain_trajectory <- function(df, domain, count_col, balance_col, balance_scale, title, subtitle, name) {
  count_sym <- rlang::ensym(count_col)
  balance_sym <- rlang::ensym(balance_col)
  p <- ggplot(df, aes(group_label, family_mean_z, group = Family6, color = Family6)) +
    geom_hline(yintercept = 0, color = "#94A3B8", linewidth = 0.45, linetype = "longdash") +
    geom_line(linewidth = 1.0, alpha = 0.72) +
    geom_point(aes(size = !!count_sym, fill = !!balance_sym), shape = 21, stroke = 0.75, color = "white") +
    geom_point(data = df %>% filter((!!count_sym) == 0), shape = 21, size = 2.1, fill = "white", color = "#CBD5E1", stroke = 0.55) +
    facet_wrap(~ display_layer, nrow = 1) +
    scale_color_manual(values = family_cols, name = "Biomarker family", drop = FALSE) +
    scale_size_continuous(range = c(2.2, 8.5), breaks = c(0, 1, 3, 6), name = paste0(domain, "\nmarker count")) +
    balance_scale +
    labs(
      title = title,
      subtitle = subtitle,
      x = "EOBC group ordered by group number",
      y = "Family mean feature z-score (RNA expression or raw TSS methylation)"
    ) +
    theme_eobc(base_size = 9.8) +
    theme(legend.position = "bottom", legend.box = "horizontal")
  save_plot(p, name, 11.2, 5.9)
  p
}

trajectory_domain_df <- trajectory_domain_df %>%
  mutate(
    os_balance = case_when(os_protective > os_adverse ~ "Protective-rich",
                           os_adverse > os_protective ~ "Adverse-rich",
                           os_count > 0 ~ "Mixed OS",
                           TRUE ~ "No OS evidence"),
    drug_balance = case_when(sensitive_count > resistant_count ~ "Sensitivity-rich",
                             resistant_count > sensitive_count ~ "Resistance-rich",
                             drug_count > 0 ~ "Mixed drug",
                             TRUE ~ "No drug evidence"),
    immune_balance = case_when(immune_pos_count > immune_neg_count ~ "Positive rho-linked",
                               immune_neg_count > immune_pos_count ~ "Negative rho-linked",
                               immune_discordant_count > 0 ~ "Discordant TIL/TMB rho",
                               immune_count > 0 ~ "Mixed score evidence",
                               TRUE ~ "No immune evidence")
  )

write_csv(trajectory_domain_df, file.path(table_dir, "family_trajectory_with_domain_counts.csv"))

p_traj_os <- plot_domain_trajectory(
  trajectory_domain_df,
  "OS",
  os_count,
  os_balance,
  scale_fill_manual(
    values = c("Protective-rich" = "#58A9F9", "Adverse-rich" = "#F4A259", "Mixed OS" = "#B9C4D1", "No OS evidence" = "white"),
    name = "OS direction", drop = FALSE
  ),
  "F. Family trajectories annotated with OS-linked group markers",
  "Line height is the family-level group program; bubble size counts OS-linked markers whose dominant state is that EOBC group.",
  "Figure_05A_OS_annotated_family_trajectory_R_v1"
)

p_traj_drug <- plot_domain_trajectory(
  trajectory_domain_df,
  "DepMap",
  drug_count,
  drug_balance,
  scale_fill_manual(
    values = c("Sensitivity-rich" = "#58A9F9", "Resistance-rich" = "#FF8A1C", "Mixed drug" = "#B9C4D1", "No drug evidence" = "white"),
    name = "Drug direction", drop = FALSE
  ),
  "G. Family trajectories annotated with DepMap drug-response evidence",
  "Bubble size counts drug-linked markers; fill indicates whether the group program is enriched for sensitivity or resistance evidence.",
  "Figure_05B_DepMap_annotated_family_trajectory_R_v1"
)

p_traj_immune <- plot_domain_trajectory(
  trajectory_domain_df,
  "Immune",
  immune_count,
  immune_balance,
  scale_fill_manual(
    values = c("Positive rho-linked" = "#17A2A4", "Negative rho-linked" = "#72A950", "Discordant TIL/TMB rho" = "#B9C4D1", "Mixed score evidence" = "#D1D9E6", "No immune evidence" = "white"),
    name = "Immune direction", drop = FALSE
  ),
  "H. Family trajectories annotated with TIL/TMB immune evidence",
  "Bubble size counts markers with significant TIL-score or TMB-log1p Spearman association (FDR < 0.10); fill shows signed rho direction.",
  "Figure_05C_Immune_annotated_family_trajectory_R_v1"
)

# ---------------------------------------------------------------------------
# Compact tables for interpretation.
# ---------------------------------------------------------------------------
domain_summary <- domain_marker %>%
  group_by(display_layer, activity_dominant_group, Family6) %>%
  summarise(
    n_markers = n(),
    os_sig = sum(os_sig, na.rm = TRUE),
    drug_sig = sum(drug_sig, na.rm = TRUE),
    immune_sig = sum(immune_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(display_layer, activity_dominant_group, Family6)
write_csv(domain_summary, file.path(table_dir, "program_family_domain_evidence_summary.csv"))

marker_story_routes <- activity_dominance %>%
  arrange(activity_dominant_group, display_layer, desc(domain_count), desc(rank_score)) %>%
  transmute(
    display_layer,
    activity_dominant_group,
    gene,
    gene_label,
    Family6,
    target_label,
    activity_contrast_z,
    kw_fdr,
    os_status,
    os_best_endpoint,
    drug_status,
    drug_clean,
    drug_class_clean,
    immune_status,
    til_activity_rho,
    tmb_activity_rho,
    domain_count
  )
write_csv(marker_story_routes, file.path(table_dir, "marker_story_routes_group_to_domain.csv"))

readme <- c(
  "EOBC group biomarker contextual domain suite v1",
  "",
  "Key design choices:",
  "- RNA rows are plotted as RNA expression z-scores.",
  "- Methylation rows are plotted as raw TSS/promoter beta-value z-scores, so RNA-METH inverse relationships remain visible.",
  "- EOBC group order is fixed as G1 | H, G2 | I, G3 | L-like, G4 | L.",
  "- Domain overlays are intentionally split into OS, DepMap, and Immune figures to avoid overloading a single synthesis plot.",
  "",
  "Main figures:",
  "Figure_01A_activity_heatmap_with_rna_meth_correlation_R_v1: sample heatmap plus RNA-METH rho and domain evidence strips.",
  "Figure_01B_RNA_TSS_methylation_correlation_matrix_R_v1: full RNA-gene by TSS-methylation-gene correlation heatmap.",
  "Figure_02A_OS_evidence_on_group_dominance_R_v1: OS evidence on group marker dominance scatter.",
  "Figure_03A_DepMap_evidence_on_group_dominance_R_v1: drug evidence on group marker dominance scatter.",
  "Figure_04A_Immune_evidence_on_group_dominance_R_v1: immune evidence on group marker dominance scatter.",
  "Figure_05A/B/C: family trajectories annotated one domain at a time."
)
writeLines(readme, file.path(out_dir, "README_contextual_domain_figures_R_v1.md"))

message("Saved contextual domain figures to: ", plot_dir)

# ===========================================================================
# V2 publication refinements
# ---------------------------------------------------------------------------
# The v1 sample-level heatmap is intentionally information-rich, but it is too
# dense for a main-text figure. The v2 figures below keep the visual grammar of
# the group trajectory/dominance plots while making each panel answer one
# manuscript question:
#   1) Which EOBC group does each biomarker represent?
#   2) Does the matched RNA/TSS-methylation pair show inverse coupling?
#   3) Which group-defining markers also carry OS, DepMap, or immune evidence?
# ===========================================================================

short_group <- c(
  "G1 | H" = "G1\nH",
  "G2 | I" = "G2\nI",
  "G3 | L-like" = "G3\nL-like",
  "G4 | L" = "G4\nL"
)

short_endpoint <- c(
  "Overall OS" = "Overall",
  "5-year OS" = "5y",
  "10-year OS" = "10y"
)

clip2 <- function(x, lim = 2.2) pmax(pmin(x, lim), -lim)

domain_plot_context <- activity_dominance %>%
  select(
    gene, modality, gene_label, Family6, target_label, display_layer,
    activity_dominant_group, dominant_mean_z, activity_contrast_z, kw_fdr,
    os_sig, os_status, drug_sig, drug_status, immune_sig, immune_status
  )

os_state_domain_v3 <- os_endpoint %>%
  mutate(selected = selected %in% TRUE) %>%
  group_by(gene, modality) %>%
  arrange(p, desc(selected), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    gene,
    modality,
    os_best_endpoint = as.character(endpoint_label),
    os_endpoint = endpoint,
    os_domain_p = p,
    os_domain_fdr = q,
    domain_x = activity_beta,
    domain_y = safe_neglog(q, 12),
    domain_effect_abs = abs(activity_beta),
    os_selected = selected
  ) %>%
  inner_join(domain_plot_context, by = c("gene", "modality")) %>%
  mutate(
    os_sig = os_sig %in% TRUE | os_selected | (!is.na(os_domain_p) & os_domain_p < 0.05),
    os_status = case_when(
      os_sig & domain_x < 0 ~ "Protective",
      os_sig & domain_x > 0 ~ "Adverse",
      TRUE ~ "Not significant"
    ),
    domain_fdr = os_domain_fdr,
    domain_p = os_domain_p,
    domain_metric = "Best OS endpoint Cox beta"
  )

depmap_state_domain_v3 <- depmap_best %>%
  left_join(
    domain_plot_context %>% select(gene, modality, dominant_mean_z),
    by = c("gene", "modality")
  ) %>%
  mutate(
    drug_sig = drug_sig %in% TRUE,
    domain_x = depmap_activity_rho,
    domain_y = safe_neglog(depmap_fdr, 12),
    domain_fdr = depmap_fdr,
    domain_effect_abs = abs(depmap_activity_rho),
    domain_metric = "Best DepMap drug-response association"
  )

immune_state_domain_v3 <- immune_best %>%
  left_join(
    domain_plot_context %>% select(gene, modality, dominant_mean_z),
    by = c("gene", "modality")
  ) %>%
  mutate(
    immune_best_fdr = pmin(coalesce(til_activity_fdr, 1), coalesce(tmb_activity_fdr, 1)),
    immune_best_metric = if_else(
      coalesce(til_activity_fdr, 1) <= coalesce(tmb_activity_fdr, 1),
      "TIL rho",
      "TMB rho"
    ),
    domain_x = if_else(immune_best_metric == "TIL rho", til_activity_rho, tmb_activity_rho),
    domain_y = safe_neglog(immune_best_fdr, 12),
    domain_fdr = immune_best_fdr,
    domain_effect_abs = abs(domain_x),
    immune_sig = immune_sig %in% TRUE,
    domain_metric = "Best TIL/TMB Spearman association"
  )

write_csv(os_state_domain_v3, file.path(table_dir, "Figure_02A_OS_state_domain_coordinates_R_v3.csv"))
write_csv(depmap_state_domain_v3, file.path(table_dir, "Figure_03A_DepMap_state_domain_coordinates_R_v3.csv"))
write_csv(immune_state_domain_v3, file.path(table_dir, "Figure_04A_Immune_state_domain_coordinates_R_v3.csv"))

make_group_activity_panel <- function(layer_name, show_x = FALSE, show_legend = FALSE) {
  layer_rows <- activity_dominance %>%
    filter(display_layer == layer_name) %>%
    arrange(activity_dominant_group, Family6, desc(rank_score), gene) %>%
    mutate(
      y = rev(row_number()),
      pretty_label = gene
    )

  layer_heat <- activity_group_means %>%
    filter(display_layer == layer_name) %>%
    inner_join(layer_rows %>% select(gene, modality, y), by = c("gene", "modality")) %>%
    mutate(
      x = as.numeric(factor(group_label, levels = group_order)),
      mean_activity_z_clip = clip2(mean_activity_z)
    )

  family_strip <- layer_rows %>%
    transmute(y, Family6)

  cor_strip2 <- layer_rows %>%
    select(y, gene, modality) %>%
    left_join(
      rna_meth_alignment %>% select(gene, sample_spearman_rho, group_spearman_rho),
      by = "gene"
    ) %>%
    pivot_longer(
      c(sample_spearman_rho, group_spearman_rho),
      names_to = "rho_type",
      values_to = "rho"
    ) %>%
    mutate(
      x = if_else(rho_type == "sample_spearman_rho", 5.15, 5.72),
      rho_label = if_else(abs(rho) >= 0.30, sprintf("%.2f", rho), "")
    )

  domain_strip2 <- layer_rows %>%
    select(y, gene, modality) %>%
    left_join(
      domain_marker %>% select(gene, modality, os_sig, drug_sig, immune_sig),
      by = c("gene", "modality")
    ) %>%
    mutate(across(c(os_sig, drug_sig, immune_sig), ~ .x %in% TRUE)) %>%
    pivot_longer(c(os_sig, drug_sig, immune_sig), names_to = "domain", values_to = "present") %>%
    mutate(
      domain = recode(domain, os_sig = "OS", drug_sig = "DepMap", immune_sig = "Immune"),
      x = recode(domain, "OS" = 6.70, "DepMap" = 7.40, "Immune" = 8.10),
      domain_fill = if_else(present, domain, "none")
    )

  ggplot() +
    geom_tile(data = family_strip, aes(0.28, y, fill = Family6), width = 0.22, height = 0.86) +
    scale_fill_manual(values = family_cols, name = "Biomarker family", drop = FALSE) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = layer_heat,
      aes(x, y, fill = mean_activity_z_clip),
      width = 0.96, height = 0.86, color = "white", linewidth = 0.22
    ) +
    scale_fill_gradient2(
      low = "#3265A8", mid = "#F8FAFC", high = "#C92A2A",
      midpoint = 0, limits = c(-2.2, 2.2), oob = squish,
      name = "Mean feature\nz-score"
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = cor_strip2,
      aes(x, y, fill = rho),
      width = 0.48, height = 0.86, color = "white", linewidth = 0.22
    ) +
    geom_text(data = cor_strip2, aes(x, y, label = rho_label), size = 1.85, color = "#111827") +
    scale_fill_gradient2(
      low = "#3265A8", mid = "#F8FAFC", high = "#C92A2A",
      midpoint = 0, limits = c(-0.8, 0.8), oob = squish,
      name = "RNA vs TSS\nSpearman rho"
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = domain_strip2,
      aes(x, y, fill = domain_fill),
      width = 0.52, height = 0.70, color = "#CBD5E1", linewidth = 0.24
    ) +
    scale_fill_manual(
      values = c(domain_cols, "none" = "#F8FAFC"),
      breaks = names(domain_cols),
      name = "Domain\nevidence"
    ) +
    geom_vline(xintercept = 4.55, color = "#111827", linewidth = 0.26) +
    geom_vline(xintercept = 6.18, color = "#111827", linewidth = 0.26) +
    scale_x_continuous(
      limits = c(0.08, 8.42),
      breaks = c(0.28, 1:4, 5.15, 5.72, 6.70, 7.40, 8.10),
      labels = if (show_x) c("Family", short_group, "sample\nrho", "group\nrho", "OS", "Drug", "Imm") else rep("", 10),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = layer_rows$y,
      labels = layer_rows$pretty_label,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(x = NULL, y = layer_name) +
    theme_eobc(base_size = 8.4) +
    theme(
      axis.title.y = element_text(face = "bold", size = 10.2, margin = margin(r = 8)),
      axis.text.y = element_text(size = 6.3),
      axis.text.x = if (show_x) element_text(size = 6.2, angle = 0, vjust = 1) else element_blank(),
      panel.grid = element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      legend.box = "horizontal",
      plot.margin = margin(1, 8, 1, 8)
    )
}

p_group_heat_rna <- make_group_activity_panel("RNA expression", show_x = FALSE, show_legend = FALSE)
p_group_heat_meth <- make_group_activity_panel("TSS methylation", show_x = TRUE, show_legend = TRUE)

fig_group_activity_heat <- p_group_heat_rna / p_group_heat_meth +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "A. EOBC group-defining biomarker expression/methylation map with RNA-TSS methylation coupling",
    subtitle = "Tiles show group mean z-scores for RNA expression and raw TSS methylation. Blue RNA-TSS rho supports promoter methylation-linked repression.",
    caption = "Rows are ordered by the highest EOBC group for each biomarker. Right-side strips show matched RNA/TSS methylation coupling and whether the biomarker has OS, DepMap, or immune evidence.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, color = "#111827"),
      plot.subtitle = element_text(size = 10.4, color = "#64748B", margin = margin(b = 6)),
      plot.caption = element_text(size = 7.4, color = "#64748B", hjust = 0, margin = margin(t = 6)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 10, 8, 10)
    )
  )
save_plot(fig_group_activity_heat, "Figure_01A_EOBC_group_activity_coupling_heatmap_R_v2", 10.8, 10.2)

p_cor_matrix_v2 <- p_cor_matrix +
  labs(
    title = "B. Matched RNA and TSS methylation coupling across biomarker genes",
    subtitle = "Diagonal cells are matched gene pairs. Negative rho values support TSS methylation-linked transcriptional repression."
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = "#111827"),
    plot.subtitle = element_text(size = 8.8, color = "#64748B", margin = margin(b = 8)),
    legend.title = element_text(size = 8.4, face = "bold"),
    legend.text = element_text(size = 7.4)
  )
save_plot(p_cor_matrix_v2, "Figure_01B_RNA_TSS_methylation_correlation_matrix_R_v2", 8.4, 7.2)

domain_sig_label <- function(df, domain = c("OS", "DepMap", "Immune")) {
  domain <- match.arg(domain)
  if (!"domain_x" %in% names(df)) {
    df <- df %>%
      mutate(
        domain_x = activity_contrast_z,
        domain_y = neg_log10_fdr,
        domain_effect_abs = abs(dominant_mean_z)
      )
  }
  if (domain == "OS") {
    df %>%
      filter(os_sig, is.finite(domain_x), is.finite(domain_y)) %>%
      mutate(
        endpoint_short = recode(as.character(os_best_endpoint), !!!short_endpoint, .default = as.character(os_best_endpoint)),
        domain_label = paste0(gene_label, "\n", endpoint_short, " ", os_status),
        domain_status = os_status
      )
  } else if (domain == "DepMap") {
    df %>%
      filter(drug_sig, is.finite(domain_x), is.finite(domain_y)) %>%
      mutate(
        domain_label = paste0(gene_label, "\n", drug_clean, " ", drug_status),
        domain_status = drug_status
      )
  } else {
    df %>%
      filter(immune_sig, is.finite(domain_x), is.finite(domain_y)) %>%
      mutate(
        domain_label = paste0(
          gene_label, "\n", immune_status,
          "\nTIL rho ", sprintf("%+.2f", til_activity_rho),
          " / TMB rho ", sprintf("%+.2f", tmb_activity_rho)
        ),
        domain_status = immune_status
      )
  }
}

domain_axis_limits <- function(plot_df, fdr_thresholds) {
  max_y <- suppressWarnings(max(c(plot_df$domain_y, -log10(fdr_thresholds)), na.rm = TRUE))
  if (!is.finite(max_y)) max_y <- 1
  x_abs <- suppressWarnings(max(abs(plot_df$domain_x), na.rm = TRUE))
  if (!is.finite(x_abs) || x_abs <= 0) x_abs <- 1
  list(
    x = c(-x_abs * 1.14, x_abs * 1.14),
    y = c(0, min(max(max_y + 0.25, 1.5), 12.7))
  )
}

state_domain_scatter <- function(df, domain, title, subtitle, status_values, status_name, file_name) {
  sig_df <- domain_sig_label(df, domain)
  status_col <- if (domain == "OS") {
    "os_status"
  } else if (domain == "DepMap") {
    "drug_status"
  } else {
    "immune_status"
  }

  p <- ggplot(df, aes(activity_contrast_z, neg_log10_fdr)) +
    geom_hline(yintercept = -log10(0.10), color = "#FFB26B", linewidth = 0.35, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), color = "#FF6B6B", linewidth = 0.45, linetype = "longdash") +
    geom_point(
      aes(fill = Family6, shape = activity_dominant_group, size = pmax(abs(dominant_mean_z), 0.10)),
      color = alpha("#94A3B8", 0.52), stroke = 0.72, alpha = 0.36
    ) +
    geom_point(
      data = sig_df,
      aes(fill = Family6, shape = activity_dominant_group, size = pmax(abs(dominant_mean_z), 0.10), color = .data[[status_col]]),
      stroke = 1.35, alpha = 0.98
    ) +
    geom_label_repel(
      data = sig_df,
      aes(label = domain_label, color = .data[[status_col]]),
      size = 2.0,
      fontface = "bold",
      label.size = 0.15,
      label.padding = unit(0.12, "lines"),
      label.r = unit(0.08, "lines"),
      fill = alpha("white", 0.92),
      min.segment.length = 0,
      segment.color = alpha("#64748B", 0.45),
      box.padding = 0.25,
      point.padding = 0.15,
      max.overlaps = Inf,
      seed = 20260518,
      show.legend = FALSE
    ) +
    facet_grid(display_layer ~ activity_dominant_group) +
    scale_fill_manual(values = family_cols, name = "Biomarker family", drop = FALSE) +
    scale_shape_manual(values = group_shapes, name = "Highest EOBC group", drop = FALSE) +
    scale_size_continuous(range = c(1.8, 6.0), name = "|Dominant\nmean z|") +
    scale_color_manual(values = status_values, name = status_name, drop = FALSE) +
    coord_cartesian(ylim = c(0, 12.7), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Group-defining contrast (highest mean feature z - lowest mean feature z)",
      y = "-log10(Kruskal-Wallis FDR), capped at 12"
    ) +
    theme_eobc(base_size = 8.6) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      axis.text.x = element_text(size = 6.6),
      axis.text.y = element_text(size = 6.6),
      strip.text = element_text(size = 8.2, face = "bold"),
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 9.4)
    )
  save_plot(p, file_name, 13.6, 7.6)
  p
}

p_os_state_v2 <- state_domain_scatter(
  activity_dominance,
  "OS",
  "C. OS-linked biomarkers within EOBC group-defining programs",
  "Each panel is one EOBC state and omics layer; colored outlines label group-defining markers with significant OS evidence.",
  os_cols,
  "OS evidence",
  "Figure_02A_OS_state_resolved_group_marker_evidence_R_v2"
)

p_depmap_state_v2 <- state_domain_scatter(
  activity_dominance,
  "DepMap",
  "D. DepMap-linked biomarkers within EOBC group-defining programs",
  "The same group-marker landscape is reused; labels show the best associated drug and whether high biomarker activity marks sensitivity or resistance.",
  drug_cols,
  "DepMap evidence",
  "Figure_03A_DepMap_state_resolved_group_marker_evidence_R_v2"
)

p_immune_state_v2 <- state_domain_scatter(
  activity_dominance,
  "Immune",
  "E. Immune-linked biomarkers within EOBC group-defining programs",
  "Significant TIL/TMB-associated markers are projected onto the EOBC state-specific biomarker dominance landscape.",
  immune_cols,
  "Immune evidence",
  "Figure_04A_Immune_state_resolved_group_marker_evidence_R_v2"
)

v2_readme <- c(
  "Contextual domain figure refinements v2",
  "",
  "Main-text candidate figures:",
  "Figure_01A_EOBC_group_activity_coupling_heatmap_R_v2: compact group-level RNA expression/raw TSS methylation heatmap plus RNA-TSS methylation coupling and domain evidence strips.",
  "Figure_01B_RNA_TSS_methylation_correlation_matrix_R_v2: refined matched RNA/TSS methylation correlation matrix.",
  "Figure_02A_OS_state_resolved_group_marker_evidence_R_v2: OS evidence over state-resolved group-marker dominance.",
  "Figure_03A_DepMap_state_resolved_group_marker_evidence_R_v2: DepMap evidence over state-resolved group-marker dominance.",
  "Figure_04A_Immune_state_resolved_group_marker_evidence_R_v2: immune evidence over state-resolved group-marker dominance.",
  "",
  "Interpretation note: methylation features are TSS/promoter methylation beta values. Group maps show raw methylation z-scores, while RNA-TSS correlation panels retain raw methylation values and expect negative correlations for transcriptional repression."
)
writeLines(v2_readme, file.path(out_dir, "README_contextual_domain_figures_R_v2.md"))

message("Saved v2 contextual domain refinements to: ", plot_dir)

# ------------------------------------------------------------------
# V3 publication polish: shorten titles, fix legend keys, and prevent
# clipping in state-resolved domain panels.
# ------------------------------------------------------------------

fig_group_activity_heat_v3 <- p_group_heat_rna / p_group_heat_meth +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "A. EOBC biomarker expression/methylation map and RNA-TSS coupling",
    subtitle = "Group mean z-scores are shown for RNA expression and raw TSS methylation beta values. Raw RNA-TSS rho at right: negative rho supports promoter methylation-linked repression.",
    caption = "Rows are ordered by the highest EOBC group for each biomarker. sample rho is across matched tumors; group rho is across the four EOBC group means. Domain strips mark OS, DepMap drug-response, or continuous TIL/TMB immune evidence.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, color = "#111827"),
      plot.subtitle = element_text(size = 8.9, color = "#64748B", margin = margin(b = 5)),
      plot.caption = element_text(size = 7.1, color = "#64748B", hjust = 0, margin = margin(t = 5)),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(7, 8, 7, 8)
    )
  )
save_plot(fig_group_activity_heat_v3, "Figure_01A_EOBC_group_activity_coupling_heatmap_R_v3", 12.4, 10.6)

p_cor_matrix_v3 <- p_cor_matrix +
  labs(
    title = "B. RNA versus TSS methylation coupling across biomarker genes",
    subtitle = "Boxed diagonal cells are matched gene pairs; negative Spearman rho is the expected direction for TSS methylation-linked transcriptional repression."
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12.4, color = "#111827"),
    plot.subtitle = element_text(size = 8.3, color = "#64748B", margin = margin(b = 7)),
    legend.title = element_text(size = 8.2, face = "bold"),
    legend.text = element_text(size = 7.2)
  )
save_plot(p_cor_matrix_v3, "Figure_01B_RNA_TSS_methylation_correlation_matrix_R_v3", 8.8, 7.25)

state_domain_scatter_v3 <- function(df, domain, title, subtitle, status_values, status_name, file_name,
                                    x_label, y_label, size_name, fdr_thresholds = c(0.10, 0.05)) {
  plot_df <- df %>%
    filter(is.finite(domain_x), is.finite(domain_y)) %>%
    mutate(
      domain_effect_abs = pmax(domain_effect_abs, 0.03),
      Family6 = factor(Family6, levels = names(family_cols)),
      display_layer = factor(display_layer, levels = modality_order),
      activity_dominant_group = factor(activity_dominant_group, levels = group_order)
    )
  sig_df <- domain_sig_label(plot_df, domain)
  status_col <- if (domain == "OS") {
    "os_status"
  } else if (domain == "DepMap") {
    "drug_status"
  } else {
    "immune_status"
  }
  lims <- domain_axis_limits(plot_df, fdr_thresholds)

  p <- ggplot(plot_df, aes(domain_x, domain_y))
  if (any(abs(fdr_thresholds - 0.25) < 1e-8)) {
    p <- p + geom_hline(yintercept = -log10(0.25), color = "#FFB26B", linewidth = 0.35, linetype = "dotted")
  }
  if (any(abs(fdr_thresholds - 0.10) < 1e-8)) {
    p <- p + geom_hline(yintercept = -log10(0.10), color = "#FFB26B", linewidth = 0.35, linetype = "dotted")
  }
  if (any(abs(fdr_thresholds - 0.05) < 1e-8)) {
    p <- p + geom_hline(yintercept = -log10(0.05), color = "#FF6B6B", linewidth = 0.45, linetype = "longdash")
  }

  p <- p +
    geom_vline(xintercept = 0, color = "#8FA0B4", linewidth = 0.38, linetype = "longdash") +
    geom_point(
      aes(fill = Family6, shape = activity_dominant_group, size = domain_effect_abs),
      color = alpha("#94A3B8", 0.50), stroke = 0.60, alpha = 0.30
    ) +
    geom_point(
      data = sig_df,
      aes(fill = Family6, shape = activity_dominant_group, size = domain_effect_abs, color = .data[[status_col]]),
      stroke = 1.20, alpha = 0.98
    ) +
    geom_label_repel(
      data = sig_df,
      aes(label = domain_label, color = .data[[status_col]]),
      size = 1.85,
      fontface = "bold",
      label.size = 0.13,
      label.padding = unit(0.10, "lines"),
      label.r = unit(0.07, "lines"),
      fill = alpha("white", 0.94),
      min.segment.length = 0,
      segment.color = alpha("#64748B", 0.36),
      box.padding = 0.20,
      point.padding = 0.12,
      max.overlaps = Inf,
      seed = 20260518,
      show.legend = FALSE
    ) +
    facet_grid(display_layer ~ activity_dominant_group) +
    scale_fill_manual(values = family_cols, name = "Biomarker family", drop = FALSE) +
    scale_shape_manual(values = group_shapes, name = "Highest EOBC group", drop = FALSE) +
    scale_size_continuous(range = c(1.5, 5.2), name = size_name) +
    scale_color_manual(values = status_values, name = status_name, drop = FALSE) +
    coord_cartesian(xlim = lims$x, ylim = lims$y, clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(shape = 21, color = "#111827", size = 3.0, alpha = 1, stroke = 0.7)
      ),
      color = guide_legend(
        order = 2,
        override.aes = list(shape = 21, fill = "white", size = 3.0, alpha = 1, stroke = 1.0)
      ),
      size = guide_legend(order = 3),
      shape = guide_legend(order = 4, override.aes = list(fill = "white", color = "#111827", size = 3.0))
    ) +
    theme_eobc(base_size = 8.4) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 2, b = 2),
      legend.box.margin = margin(t = 2, b = 0),
      axis.text.x = element_text(size = 6.4),
      axis.text.y = element_text(size = 6.4),
      strip.text = element_text(size = 8.0, face = "bold"),
      plot.title = element_text(size = 14.8, face = "bold"),
      plot.subtitle = element_text(size = 8.7),
      plot.margin = margin(7, 8, 10, 8)
    )
  save_plot(p, file_name, 14.4, 8.8)
  p
}

p_os_state_v3 <- state_domain_scatter_v3(
  os_state_domain_v3,
  "OS",
  "C. OS evidence within EOBC group-defining biomarker programs",
  "Coordinates show the best OS endpoint Cox beta and FDR; facets retain the EOBC program where each marker is dominant.",
  os_cols,
  "OS evidence",
  "Figure_02A_OS_state_resolved_group_marker_evidence_R_v3",
  "Cox beta for best OS endpoint",
  "-log10(OS Cox FDR), capped at 12",
  "|Cox beta|",
  c(0.10, 0.05)
)

p_depmap_state_v3 <- state_domain_scatter_v3(
  depmap_state_domain_v3,
  "DepMap",
  "D. DepMap drug-response evidence within EOBC group-defining programs",
  "Coordinates show the signed DepMap drug-response association and its FDR; labels retain the best associated drug.",
  drug_cols,
  "DepMap evidence",
  "Figure_03A_DepMap_state_resolved_group_marker_evidence_R_v3",
  "Signed DepMap association with drug response",
  "-log10(DepMap FDR), capped at 12",
  "|DepMap rho|",
  c(0.25, 0.10)
)

p_immune_state_v3 <- state_domain_scatter_v3(
  immune_state_domain_v3,
  "Immune",
  "E. Immune evidence within EOBC group-defining biomarker programs",
  "Coordinates show the strongest TIL/TMB Spearman rho and best FDR; labels report both signed rho values.",
  immune_cols,
  "Immune evidence",
  "Figure_04A_Immune_state_resolved_group_marker_evidence_R_v3",
  "Signed strongest TIL/TMB Spearman rho",
  "-log10(best TIL/TMB FDR), capped at 12",
  "|Best rho|",
  c(0.10, 0.05)
)

v3_readme <- c(
  "Contextual domain figure refinements v3",
  "",
  "V3 is the recommended publication-polish set.",
  "Figure_01A_EOBC_group_activity_coupling_heatmap_R_v3: group-level RNA expression and raw TSS methylation heatmap with RNA-TSS correlation and OS/DepMap/immune evidence strips.",
  "Figure_01B_RNA_TSS_methylation_correlation_matrix_R_v3: matched RNA-vs-TSS methylation correlation matrix.",
  "Figure_02A_OS_state_resolved_group_marker_evidence_R_v3: OS evidence plotted by Cox beta and OS FDR within EOBC state marker programs.",
  "Figure_03A_DepMap_state_resolved_group_marker_evidence_R_v3: DepMap evidence plotted by signed drug-response association and DepMap FDR within EOBC state marker programs.",
  "Figure_04A_Immune_state_resolved_group_marker_evidence_R_v3: immune evidence plotted by strongest TIL/TMB rho and best TIL/TMB FDR within EOBC state marker programs.",
  "",
  "Interpretation note: TSS methylation is shown as raw beta-value z-scores in group maps. Raw RNA-TSS correlations are retained in the correlation panels; negative rho is the biologically expected direction for promoter repression."
)
writeLines(v3_readme, file.path(out_dir, "README_contextual_domain_figures_R_v3.md"))

message("Saved v3 publication-polish contextual domain figures to: ", plot_dir)

# ---------------------------------------------------------------------------
# V2 family trajectory overlays: preserve the line-plot grammar while showing
# where OS, DepMap, and immune evidence lands on each EOBC group program.
# ---------------------------------------------------------------------------
plot_domain_trajectory_v2 <- function(df, domain, count_col, balance_col, balance_values,
                                      title, subtitle, file_name) {
  count_name <- rlang::as_name(rlang::ensym(count_col))
  balance_name <- rlang::as_name(rlang::ensym(balance_col))

  plot_df <- df %>%
    mutate(
      group_index = as.integer(group_label),
      evidence_count = .data[[count_name]],
      evidence_class = .data[[balance_name]]
    )

  label_df <- plot_df %>%
    group_by(display_layer, Family6) %>%
    slice_max(order_by = abs(family_mean_z), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      label_nudge = if_else(family_mean_z >= 0, 0.12, -0.12),
      label_y = family_mean_z + label_nudge
    )

  p <- ggplot(plot_df, aes(group_index, family_mean_z, group = Family6, color = Family6)) +
    geom_hline(yintercept = 0, color = "#94A3B8", linewidth = 0.42, linetype = "longdash") +
    geom_line(linewidth = 1.15, alpha = 0.70, lineend = "round") +
    geom_point(
      shape = 21, size = 2.1, stroke = 0.65,
      fill = "white", color = alpha("#94A3B8", 0.70)
    ) +
    geom_point(
      data = plot_df %>% filter(evidence_count > 0),
      aes(size = evidence_count, fill = evidence_class),
      shape = 21, stroke = 0.90, color = "white", alpha = 0.98
    ) +
    geom_text_repel(
      data = label_df,
      aes(x = group_index, y = label_y, label = Family6, color = Family6),
      inherit.aes = FALSE,
      size = 3.0,
      fontface = "bold",
      min.segment.length = Inf,
      box.padding = 0.18,
      point.padding = 0.08,
      max.overlaps = Inf,
      seed = 20260518,
      show.legend = FALSE
    ) +
    facet_wrap(~ display_layer, nrow = 1) +
    scale_x_continuous(
      breaks = seq_along(group_order),
      labels = group_order,
      limits = c(0.72, 4.28),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    scale_color_manual(values = family_cols, guide = "none", drop = FALSE) +
    scale_fill_manual(values = balance_values, name = paste0(domain, " evidence"), drop = FALSE) +
    scale_size_continuous(
      range = c(2.8, 10.0),
      breaks = c(1, 2, 4, 6),
      name = paste0(domain, "\nmarker count")
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "EOBC group ordered by group number",
      y = "Family mean feature z-score"
    ) +
    guides(
      fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 4.4, color = "white", stroke = 0.8)),
      size = guide_legend(order = 2, override.aes = list(fill = "#E2E8F0", color = "#111827", stroke = 0.8))
    ) +
    theme_eobc(base_size = 9.4) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 2, b = 2),
      legend.box.margin = margin(t = 4, b = 0),
      panel.grid.major.x = element_line(color = "#DCE6F2", linewidth = 0.50),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 24, hjust = 1, size = 8.0, face = "bold"),
      strip.text = element_text(size = 9.0, face = "bold"),
      plot.title = element_text(size = 15.4, face = "bold"),
      plot.subtitle = element_text(size = 9.0),
      plot.margin = margin(10, 14, 12, 14)
    )

  save_plot(p, file_name, 13.8, 6.7)
  p
}

p_traj_os_v2 <- plot_domain_trajectory_v2(
  trajectory_domain_df,
  "OS",
  os_count,
  os_balance,
  c("Protective-rich" = "#58A9F9", "Adverse-rich" = "#F4A259", "Mixed OS" = "#B9C4D1", "No OS evidence" = "white"),
  "F. OS-linked evidence on EOBC biological-family trajectories",
  "Lines show group mean biomarker-family activity; bubbles count group-defining markers with OS evidence, while group separation is supported by Kruskal-Wallis FDR in the state panels.",
  "Figure_05A_OS_annotated_family_trajectory_R_v2"
)

p_traj_depmap_v2 <- plot_domain_trajectory_v2(
  trajectory_domain_df,
  "DepMap",
  drug_count,
  drug_balance,
  c("Sensitivity-rich" = "#58A9F9", "Resistance-rich" = "#FF8A1C", "Mixed drug" = "#B9C4D1", "No drug evidence" = "white"),
  "G. DepMap drug-response evidence on EOBC biological-family trajectories",
  "Lines show group mean biomarker-family activity; bubbles count group-defining markers with DepMap drug sensitivity/resistance evidence (FDR < 0.25).",
  "Figure_05B_DepMap_annotated_family_trajectory_R_v2"
)

p_traj_immune_v2 <- plot_domain_trajectory_v2(
  trajectory_domain_df,
  "Immune",
  immune_count,
  immune_balance,
  c("Positive rho-linked" = "#17A2A4", "Negative rho-linked" = "#72A950", "Discordant TIL/TMB rho" = "#B9C4D1", "Mixed score evidence" = "#D1D9E6", "No immune evidence" = "white"),
  "H. TIL/TMB immune evidence on EOBC biological-family trajectories",
  "Lines show group mean biomarker-family activity; bubbles count markers with significant TIL-score or TMB-log1p Spearman association (FDR < 0.10), with fill showing signed rho direction.",
  "Figure_05C_Immune_annotated_family_trajectory_R_v2"
)

v3_readme_append <- c(
  "",
  "Additional line-plot overlays:",
  "Figure_05A_OS_annotated_family_trajectory_R_v2: family trajectories with OS evidence bubbles; group-defining separation is supported by Kruskal-Wallis FDR in the state panels.",
  "Figure_05B_DepMap_annotated_family_trajectory_R_v2: family trajectories with DepMap drug-response evidence bubbles.",
  "Figure_05C_Immune_annotated_family_trajectory_R_v2: family trajectories with continuous TIL-score/TMB-log1p evidence bubbles; no external hot/cold subtype calls are used."
)
cat(paste0(v3_readme_append, collapse = "\n"), file = file.path(out_dir, "README_contextual_domain_figures_R_v3.md"), append = TRUE)

message("Saved v2 trajectory-domain overlay figures to: ", plot_dir)
