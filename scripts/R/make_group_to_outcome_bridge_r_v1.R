suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

set.seed(20260517)

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
group_dir <- file.path(base_dir, "group_biomarker_landscape_r_v1")
os_dir <- file.path(base_dir, "os_elasticnet_full_candidate_r_v1")
functional_dir <- file.path(base_dir, "full_layer_candidate_validation_r_v1")

out_dir <- file.path(base_dir, "group_to_outcome_bridge_r_v1")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#111827"
muted <- "#64748B"
grid_col <- "#D8E3EF"
protective_col <- "#4EA5F0"
adverse_col <- "#F4A259"

group_order <- c("G1 | H", "G2 | I", "G3 | L-like", "G4 | L")
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
family_shapes <- c(
  "Immune" = 21,
  "Repair" = 22,
  "Glycolysis / TCA" = 24,
  "Fatty acid" = 25,
  "Kinase signaling" = 23,
  "Hormone signaling" = 21
)

theme_eobc <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
      panel.grid.minor = element_line(colour = alpha(grid_col, 0.42), linewidth = 0.20),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = muted),
      strip.text = element_text(face = "bold", colour = ink),
      strip.background = element_rect(fill = "#EEF4FB", colour = "#D6E0EA", linewidth = 0.35),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", colour = ink),
      legend.text = element_text(colour = muted),
      plot.title = element_text(face = "bold", colour = ink, size = rel(1.35)),
      plot.subtitle = element_text(colour = muted, size = rel(0.92)),
      plot.caption = element_text(colour = muted, size = rel(0.72), hjust = 0),
      plot.title.position = "plot"
    )
}

clip_val <- function(x, lim = 1.5) pmax(pmin(x, lim), -lim)
neglog <- function(p) ifelse(is.na(p), NA_real_, -log10(pmax(p, 1e-300)))
safe_num <- function(x) suppressWarnings(as.numeric(x))

group_means <- read_csv(file.path(group_dir, "tables", "group_biomarker_group_means.csv"), show_col_types = FALSE)
dominance_raw <- read_csv(file.path(group_dir, "tables", "group_biomarker_dominant_group_summary.csv"), show_col_types = FALSE)
inverse_alignment <- read_csv(file.path(group_dir, "tables", "rna_meth_tss_inverse_alignment.csv"), show_col_types = FALSE)
os_gene <- read_csv(file.path(os_dir, "tables", "os_gene_univariate_volcano_with_elasticnet_overlay.csv"), show_col_types = FALSE)
os_model <- read_csv(file.path(os_dir, "tables", "os_elasticnet_model_summary.csv"), show_col_types = FALSE)
depmap <- read_csv(file.path(functional_dir, "tables", "depmap_all_candidate_best_hits_by_omics.csv"), show_col_types = FALSE)
immune <- read_csv(file.path(functional_dir, "tables", "immune_all_candidate_til_tmb_correlations_by_omics_wide.csv"), show_col_types = FALSE)

activity_layer_label <- function(modality) {
  ifelse(modality == "RNA", "RNA activity", "TSS methylation-inferred activity")
}

# For TSS methylation, higher methylation is interpreted as lower transcriptional activity.
activity_means <- group_means %>%
  mutate(
    activity_z = if_else(modality == "METH", -mean_z, mean_z),
    activity_z_clip = clip_val(activity_z, 1.5),
    activity_layer = activity_layer_label(modality),
    group_label = factor(group_label, levels = group_order)
  )

activity_dominance <- dominance_raw %>%
  mutate(
    activity_layer = activity_layer_label(modality),
    activity_dominant_group = if_else(modality == "METH", suppressed_group, dominant_group),
    activity_suppressed_group = if_else(modality == "METH", dominant_group, suppressed_group),
    activity_peak_z = if_else(modality == "METH", -suppressed_mean_z, dominant_mean_z),
    activity_low_z = if_else(modality == "METH", -dominant_mean_z, suppressed_mean_z),
    activity_contrast_z = contrast_z,
    group_definition_tier = case_when(
      kw_fdr < 0.05 & activity_contrast_z >= 0.75 ~ "Group-defining",
      kw_fdr < 0.10 & activity_contrast_z >= 0.50 ~ "Supportive",
      TRUE ~ "Exploratory"
    ),
    activity_dominant_group = factor(activity_dominant_group, levels = group_order),
    modality_endpoint = if_else(modality == "RNA", "RNA", "METH")
  ) %>%
  left_join(
    inverse_alignment %>%
      select(
        gene,
        sample_spearman_rho,
        group_spearman_rho,
        rna_matches_meth_low_activity,
        inverse_consistency_class
      ),
    by = "gene"
  )

write_csv(activity_means, file.path(table_dir, "biomarker_group_activity_means_tss_meth_inverted.csv"))
write_csv(activity_dominance, file.path(table_dir, "biomarker_activity_dominance_tss_meth_inverted.csv"))

priority_gene_layers <- activity_dominance %>%
  filter(group_definition_tier %in% c("Group-defining", "Supportive")) %>%
  mutate(
    row_id = paste0(gene_label, " | ", if_else(modality == "RNA", "RNA", "METH activity")),
    group_sort = as.integer(activity_dominant_group),
    tier_sort = if_else(group_definition_tier == "Group-defining", 0, 1)
  ) %>%
  arrange(group_sort, tier_sort, desc(activity_contrast_z), gene)

row_levels <- rev(unique(priority_gene_layers$row_id))

activity_heatmap_df <- activity_means %>%
  mutate(row_id = paste0(gene_label, " | ", if_else(modality == "RNA", "RNA", "METH activity"))) %>%
  semi_join(priority_gene_layers %>% select(gene, modality), by = c("gene", "modality")) %>%
  left_join(
    activity_dominance %>%
      select(gene, modality, activity_dominant_group, activity_contrast_z, kw_fdr, group_definition_tier),
    by = c("gene", "modality")
  ) %>%
  mutate(
    row_id = factor(row_id, levels = row_levels),
    activity_layer = factor(activity_layer, levels = c("RNA activity", "TSS methylation-inferred activity")),
    is_peak = as.character(group_label) == as.character(activity_dominant_group),
    peak_label = case_when(
      is_peak & kw_fdr < 0.05 ~ "**",
      is_peak & kw_fdr < 0.10 ~ "*",
      TRUE ~ ""
    )
  )

p_activity <- ggplot(activity_heatmap_df, aes(group_label, row_id)) +
  geom_tile(aes(fill = activity_z_clip), colour = "white", linewidth = 0.55) +
  geom_point(
    data = activity_heatmap_df %>% filter(is_peak),
    aes(size = activity_contrast_z, colour = activity_dominant_group),
    shape = 21, fill = "white", stroke = 1.05
  ) +
  geom_text(
    data = activity_heatmap_df %>% filter(is_peak, peak_label != ""),
    aes(label = peak_label),
    colour = ink, fontface = "bold", size = 3.2, vjust = 0.55
  ) +
  facet_grid(. ~ activity_layer, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#4C78A8", mid = "#F8FAFC", high = "#D73027",
    midpoint = 0, limits = c(-1.5, 1.5), oob = squish,
    name = "Inferred activity\nz-score"
  ) +
  scale_colour_manual(values = group_cols, name = "Activity-high group", drop = FALSE) +
  scale_size_continuous(range = c(1.8, 5.6), name = "Group contrast") +
  labs(
    title = "A. EOBC group-defining biomarker activity map",
    subtitle = "RNA is shown as expression activity; TSS methylation is inverted so that high activity means low promoter/TSS methylation.",
    x = NULL,
    y = NULL,
    caption = "** FDR<0.05, * FDR<0.10 by Kruskal-Wallis across EOBC groups. Only group-defining/supportive marker-layers are shown."
  ) +
  theme_eobc(9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

ggsave(file.path(plot_dir, "Figure_01_EOBC_group_defining_activity_map_R_v1.png"), p_activity, width = 15.5, height = 11, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_01_EOBC_group_defining_activity_map_R_v1.pdf"), p_activity, width = 15.5, height = 11, device = cairo_pdf, bg = "white")

# OS bridge: beta is also converted to inferred activity direction for methylation.
os_activity <- os_gene %>%
  mutate(
    modality = if_else(modality == "METH", "METH", "RNA"),
    activity_layer = activity_layer_label(modality),
    endpoint_label = factor(endpoint_label, levels = c("Overall OS", "5-year OS", "10-year OS")),
    activity_beta = if_else(modality == "METH", -beta, beta),
    activity_hr = exp(activity_beta),
    os_direction = case_when(
      activity_beta > 0 ~ "Higher inferred activity = adverse",
      activity_beta < 0 ~ "Higher inferred activity = protective",
      TRUE ~ "Neutral"
    ),
    os_neglog10_p = neglog(p),
    os_neglog10_q = neglog(q),
    selected = as.logical(selected)
  ) %>%
  left_join(
    activity_dominance %>%
      select(
        gene, modality, activity_dominant_group,
        activity_contrast_z, kw_fdr, group_definition_tier, sample_spearman_rho, group_spearman_rho,
        inverse_consistency_class
      ),
    by = c("gene", "modality")
  )

write_csv(os_activity, file.path(table_dir, "os_gene_activity_bridge_tss_meth_inverted.csv"))

os_label_df <- os_activity %>%
  filter(
    selected |
      p < 0.05 |
      (group_definition_tier == "Group-defining" & os_neglog10_p >= 1.0)
  ) %>%
  group_by(activity_layer, endpoint_label) %>%
  arrange(desc(selected), desc(os_neglog10_p), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

p_os_volcano <- ggplot(os_activity, aes(activity_beta, os_neglog10_p)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#E8F3FF", alpha = 0.46) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FFF1E8", alpha = 0.46) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#EF6F6C", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#94A3B8", linewidth = 0.45) +
  geom_point(
    aes(fill = activity_dominant_group, shape = Family6, size = stability),
    colour = ink, alpha = 0.78, stroke = 0.38
  ) +
  geom_point(
    data = os_activity %>% filter(selected),
    aes(fill = activity_dominant_group, shape = Family6),
    size = 4.3, colour = "#FFD166", stroke = 1.25, alpha = 1
  ) +
  geom_text_repel(
    data = os_label_df,
    aes(label = gene_label, colour = activity_dominant_group),
    size = 2.55,
    box.padding = 0.28,
    point.padding = 0.20,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  facet_grid(activity_layer ~ endpoint_label, scales = "free_x") +
  scale_fill_manual(values = group_cols, name = "EOBC activity-high group", drop = FALSE, na.value = "#CBD5E1") +
  scale_colour_manual(values = group_cols, drop = FALSE, na.value = muted) +
  scale_shape_manual(values = family_shapes, name = "Biomarker family", drop = FALSE) +
  scale_size_continuous(range = c(1.8, 4.2), limits = c(0, 1), name = "Elastic-net stability") +
  labs(
    title = "B. OS bridge: group-defining biomarkers mapped to survival association",
    subtitle = "Cox beta is interpreted as inferred biomarker activity; methylation beta is sign-inverted because TSS methylation represses expression.",
    x = "Cox beta for inferred activity (negative = protective, positive = adverse)",
    y = "-log10(univariate Cox P)"
  ) +
  theme_eobc(8.5)

ggsave(file.path(plot_dir, "Figure_02A_OS_group_defined_activity_volcano_R_v1.png"), p_os_volcano, width = 15.5, height = 9.8, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_02A_OS_group_defined_activity_volcano_R_v1.pdf"), p_os_volcano, width = 15.5, height = 9.8, device = cairo_pdf, bg = "white")

os_heatmap_df <- os_activity %>%
  semi_join(priority_gene_layers %>% select(gene, modality), by = c("gene", "modality")) %>%
  mutate(
    row_id = paste0(gene_label, " | ", if_else(modality == "RNA", "RNA", "METH activity")),
    row_id = factor(row_id, levels = row_levels),
    endpoint_label = factor(endpoint_label, levels = c("Overall OS", "5-year OS", "10-year OS")),
    activity_beta_clip = clip_val(activity_beta, 1.0),
    selected_label = if_else(selected, "★", ""),
    p_label = case_when(
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

p_os_heatmap <- ggplot(os_heatmap_df, aes(endpoint_label, row_id)) +
  geom_tile(aes(fill = activity_beta_clip), colour = "white", linewidth = 0.55) +
  geom_point(aes(size = os_neglog10_p), shape = 21, fill = "white", colour = ink, alpha = 0.82, stroke = 0.45) +
  geom_text(aes(label = selected_label), colour = "#F59E0B", fontface = "bold", size = 3.6, vjust = 0.55) +
  facet_grid(. ~ activity_layer, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = protective_col, mid = "#F8FAFC", high = adverse_col,
    midpoint = 0, limits = c(-1, 1), oob = squish,
    name = "OS activity beta"
  ) +
  scale_size_continuous(range = c(1.2, 5.0), name = "-log10(P)") +
  labs(
    title = "C. OS evidence matrix for EOBC group-defining biomarkers",
    subtitle = "Blue indicates higher inferred activity is protective; orange indicates higher inferred activity is adverse. Star marks elastic-net selected biomarkers.",
    x = NULL,
    y = NULL
  ) +
  theme_eobc(8.5) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid = element_blank())

ggsave(file.path(plot_dir, "Figure_02B_OS_group_defined_activity_matrix_R_v1.png"), p_os_heatmap, width = 13.5, height = 10.5, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_OS_group_defined_activity_matrix_R_v1.pdf"), p_os_heatmap, width = 13.5, height = 10.5, device = cairo_pdf, bg = "white")

# Functional bridge: DepMap and immune effects are converted to inferred activity for methylation.
depmap_activity <- depmap %>%
  mutate(
    modality = if_else(omics == "METH", "METH", "RNA"),
    activity_layer = activity_layer_label(modality),
    depmap_activity_rho = if_else(modality == "METH", -rho, rho),
    depmap_activity_abs = abs(depmap_activity_rho),
    depmap_activity_direction = case_when(
      is.na(depmap_activity_rho) ~ "Not tested",
      depmap_activity_rho > 0 ~ "Higher inferred activity = resistant",
      depmap_activity_rho < 0 ~ "Higher inferred activity = sensitive",
      TRUE ~ "Neutral"
    ),
    depmap_evidence_tier = case_when(
      is.na(fdr) ~ "Not tested",
      fdr < 0.10 ~ "FDR < 0.10",
      fdr < 0.25 ~ "FDR < 0.25",
      TRUE ~ "Exploratory"
    ),
    depmap_label = if_else(is.na(drug_clean), "", paste0(drug_clean, "\n", drug_target_class))
  ) %>%
  select(
    gene, modality, drug_clean, drug_target_class, depmap_activity_rho,
    depmap_activity_abs, depmap_activity_direction, depmap_evidence_tier, depmap_label,
    depmap_fdr = fdr, depmap_neglog10_fdr = neglog10_fdr, feature_available
  )

immune_activity <- immune %>%
  mutate(
    modality = if_else(omics == "METH", "METH", "RNA"),
    activity_layer = activity_layer_label(modality),
    til_activity_rho = if_else(modality == "METH", -rho_TIL_score, rho_TIL_score),
    tmb_activity_rho = if_else(modality == "METH", -rho_TMB, rho_TMB),
    til_activity_fdr = fdr_TIL_score,
    tmb_activity_fdr = fdr_TMB,
    immune_activity_strength = sqrt(til_activity_rho^2 + tmb_activity_rho^2),
    immune_activity_quadrant = case_when(
      til_activity_rho >= 0 & tmb_activity_rho >= 0 ~ "Immune-hot / TMB-high",
      til_activity_rho < 0 & tmb_activity_rho < 0 ~ "Immune-cold / TMB-low",
      til_activity_rho >= 0 & tmb_activity_rho < 0 ~ "TIL-high / TMB-low",
      til_activity_rho < 0 & tmb_activity_rho >= 0 ~ "TMB-shifted",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    gene, modality, til_activity_rho, tmb_activity_rho, til_activity_fdr, tmb_activity_fdr,
    immune_activity_strength, immune_activity_quadrant
  )

functional_bridge <- activity_dominance %>%
  select(
    gene, modality, gene_label, Family6, target_label, activity_layer,
    activity_dominant_group, activity_contrast_z, kw_fdr, group_definition_tier,
    sample_spearman_rho, group_spearman_rho, inverse_consistency_class
  ) %>%
  left_join(depmap_activity, by = c("gene", "modality")) %>%
  left_join(immune_activity, by = c("gene", "modality")) %>%
  left_join(
    os_activity %>%
      group_by(gene, modality) %>%
      summarise(
        os_min_p = min(p, na.rm = TRUE),
        os_best_endpoint = endpoint_label[which.min(p)][1],
        os_best_activity_beta = activity_beta[which.min(p)][1],
        os_selected_any = any(selected, na.rm = TRUE),
        os_selected_endpoints = paste(endpoint_label[selected], collapse = "; "),
        .groups = "drop"
      ),
    by = c("gene", "modality")
  ) %>%
  mutate(
    os_best_neglog10_p = neglog(os_min_p),
    os_activity_direction = case_when(
      os_best_activity_beta > 0 ~ "Adverse",
      os_best_activity_beta < 0 ~ "Protective",
      TRUE ~ "Neutral"
    )
  )

write_csv(functional_bridge, file.path(table_dir, "group_to_os_drug_immune_activity_bridge.csv"))

functional_long <- functional_bridge %>%
  filter(group_definition_tier %in% c("Group-defining", "Supportive")) %>%
  mutate(
    row_id = paste0(gene_label, " | ", if_else(modality == "RNA", "RNA", "METH activity")),
    row_id = factor(row_id, levels = row_levels),
    activity_dominant_group = factor(activity_dominant_group, levels = group_order)
  ) %>%
  transmute(
    gene, modality, row_id, gene_label, Family6, activity_layer, activity_dominant_group,
    `Drug response` = depmap_activity_rho,
    `TIL association` = til_activity_rho,
    `TMB association` = tmb_activity_rho,
    `Best OS association` = os_best_activity_beta,
    depmap_label,
    drug_clean,
    drug_target_class,
    depmap_evidence_tier,
    os_selected_any,
    os_best_endpoint
  ) %>%
  pivot_longer(
    cols = c(`Drug response`, `TIL association`, `TMB association`, `Best OS association`),
    names_to = "evidence_axis",
    values_to = "signed_effect"
  ) %>%
  mutate(
    evidence_axis = factor(evidence_axis, levels = c("Best OS association", "Drug response", "TIL association", "TMB association")),
    signed_effect_clip = clip_val(signed_effect, 0.85),
    label = case_when(
      evidence_axis == "Drug response" & !is.na(drug_clean) ~ drug_clean,
      evidence_axis == "Best OS association" & os_selected_any ~ paste0("Elastic-net\n", os_best_endpoint),
      evidence_axis %in% c("TIL association", "TMB association") & abs(signed_effect) >= 0.35 ~ sprintf("%.2f", signed_effect),
      TRUE ~ ""
    )
  )

p_func_matrix <- ggplot(functional_long, aes(evidence_axis, row_id)) +
  geom_tile(aes(fill = signed_effect_clip), colour = "white", linewidth = 0.55) +
  geom_text(aes(label = label), size = 2.15, lineheight = 0.82, colour = ink) +
  facet_grid(. ~ activity_layer, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = protective_col, mid = "#F8FAFC", high = adverse_col,
    midpoint = 0, limits = c(-0.85, 0.85), oob = squish,
    name = "Signed activity effect"
  ) +
  labs(
    title = "D. Integrated outcome bridge from EOBC biomarker programs",
    subtitle = "Rows are group-defining marker-layers. Methylation is interpreted as inferred activity (-TSS methylation). Drug: positive=resistance, negative=sensitivity; OS: positive=adverse, negative=protective.",
    x = NULL,
    y = NULL,
    caption = "Drug labels show the best associated DepMap compound for each marker-layer. Numeric immune labels are shown for stronger TIL/TMB associations."
  ) +
  theme_eobc(8.5) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank())

ggsave(file.path(plot_dir, "Figure_03_group_to_OS_drug_immune_bridge_matrix_R_v1.png"), p_func_matrix, width = 14.5, height = 11.2, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_03_group_to_OS_drug_immune_bridge_matrix_R_v1.pdf"), p_func_matrix, width = 14.5, height = 11.2, device = cairo_pdf, bg = "white")

state_summary <- functional_bridge %>%
  filter(group_definition_tier %in% c("Group-defining", "Supportive")) %>%
  group_by(activity_layer, activity_dominant_group) %>%
  summarise(
    n_marker_layers = n(),
    top_families = paste(names(sort(table(Family6), decreasing = TRUE))[1:min(3, n_distinct(Family6))], collapse = ", "),
    os_selected_genes = paste(unique(gene_label[os_selected_any]), collapse = "; "),
    drug_sensitive_genes = paste(unique(gene_label[depmap_activity_rho < -0.35]), collapse = "; "),
    drug_resistant_genes = paste(unique(gene_label[depmap_activity_rho > 0.35]), collapse = "; "),
    immune_hot_genes = paste(unique(gene_label[til_activity_rho > 0.25 & tmb_activity_rho > 0]), collapse = "; "),
    immune_cold_genes = paste(unique(gene_label[til_activity_rho < -0.25 & tmb_activity_rho < 0]), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(activity_layer, activity_dominant_group)

write_csv(state_summary, file.path(table_dir, "state_program_to_outcome_summary.csv"))

summary_md <- file.path(out_dir, "group_to_outcome_bridge_interpretation.md")
writeLines(c(
  "# EOBC group-to-outcome bridge analysis",
  "",
  "## Analysis rule",
  "- RNA features are interpreted as expression activity.",
  "- TSS methylation features are interpreted as inferred transcriptional activity after sign inversion: activity = -TSS methylation.",
  "- Therefore, methylation-high states are not automatically gene-active states; they indicate likely transcriptional suppression at the corresponding TSS feature.",
  "",
  "## Main outputs",
  "- Figure_01_EOBC_group_defining_activity_map_R_v1: group-defining biomarker activity map.",
  "- Figure_02A_OS_group_defined_activity_volcano_R_v1: OS association of group-defining biomarkers.",
  "- Figure_02B_OS_group_defined_activity_matrix_R_v1: compact OS evidence matrix.",
  "- Figure_03_group_to_OS_drug_immune_bridge_matrix_R_v1: integrated OS/drug/immune consequence matrix.",
  "",
  "## State-level summary",
  paste(capture.output(print(state_summary, n = Inf)), collapse = "\n")
), summary_md, useBytes = TRUE)

message("Wrote group-to-outcome bridge analysis to: ", out_dir)
