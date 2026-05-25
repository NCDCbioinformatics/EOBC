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

bridge_dir <- file.path(base_dir, "group_to_outcome_bridge_r_v1", "tables")
validation_dir <- file.path(base_dir, "full_layer_candidate_validation_r_v1", "tables")
out_dir <- file.path(base_dir, "group_marker_domain_specific_r_v1")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#111827"
muted <- "#64748B"
grid_col <- "#D8E3EF"

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
layer_cols <- c("RNA activity" = "#0F766E", "TSS-low activity" = "#C45A11")
os_cols <- c("Protective" = "#4EA5F0", "Adverse" = "#F4A259", "Not significant" = "#D7DEE8")
drug_cols <- c("Sensitive" = "#4EA5F0", "Resistant" = "#F4A259", "Not prioritized" = "#D7DEE8")
immune_cols <- c("TIL-high/TMB-high" = "#40C7D8", "TMB-shifted" = "#8FAAD0", "Immune-cold/TMB-low" = "#8BBE69", "Mixed/other" = "#D7DEE8")

theme_eobc <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
      panel.grid.minor = element_line(colour = alpha(grid_col, 0.35), linewidth = 0.18),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = muted),
      strip.text = element_text(face = "bold", colour = ink),
      strip.background = element_rect(fill = "#EEF4FB", colour = "#D6E0EA", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", colour = ink),
      legend.text = element_text(colour = muted),
      plot.title = element_text(face = "bold", colour = ink, size = rel(1.25)),
      plot.subtitle = element_text(colour = muted, size = rel(0.92)),
      plot.caption = element_text(colour = muted, size = rel(0.72), hjust = 0),
      plot.title.position = "plot"
    )
}

clean_drug_class <- function(x, drug = NULL) {
  x0 <- str_to_lower(coalesce(x, ""))
  d0 <- str_to_lower(coalesce(drug, ""))
  case_when(
    str_detect(x0, "checkpoint|ddr|dna damage|parp|atr|wee|chk|topoisomerase") |
      str_detect(d0, "adavosertib|ceralasertib|alisertib|rucaparib|olaparib") ~ "DDR / checkpoint",
    str_detect(x0, "rtk|mapk|egfr|mek|raf|erk|src|abl") |
      str_detect(d0, "dasatinib|nilotinib|vx-11e|erlotinib") ~ "RTK/MAPK",
    str_detect(x0, "mitotic|cell-cycle|microtubule|aurora|polo|cdk") |
      str_detect(d0, "bi-2536|velban|vinblastine|palbociclib|mk 0457|rebemadlin") ~ "Mitotic / cell-cycle",
    str_detect(x0, "pi3k|akt|mtor") | str_detect(d0, "azd-8055|rapamycin") ~ "PI3K/AKT/mTOR",
    str_detect(x0, "metabolic|wnt") | str_detect(d0, "agi-6780|sb 216763") ~ "Metabolic / WNT",
    str_detect(x0, "cytotoxic|chemotherapy|taxane|docetaxel|topotecan|carmustine|mitoxantrone|elephantin") |
      str_detect(d0, "docetaxel|carmustine|mitoxantrone|elephantin") ~ "Cytotoxic / chemotherapy",
    TRUE ~ "Other / exploratory"
  )
}

short_layer <- function(x) {
  if_else(x == "RNA activity", "RNA activity", "TSS-low activity")
}

format_p <- function(x) {
  if_else(is.na(x), "NA", if_else(x < 1e-4, "<1e-4", sprintf("%.3g", x)))
}

bridge <- read_csv(
  file.path(bridge_dir, "group_to_os_drug_immune_activity_bridge.csv"),
  show_col_types = FALSE
) %>%
  filter(group_definition_tier %in% c("Group-defining", "Supportive")) %>%
  mutate(
    activity_dominant_group = factor(activity_dominant_group, levels = group_order),
    activity_layer = factor(activity_layer, levels = c("RNA activity", "TSS methylation-inferred activity")),
    layer_story = short_layer(as.character(activity_layer)),
    program = paste0(as.character(activity_dominant_group), "\n", layer_story),
    program = factor(
      program,
      levels = c(
        "G1 | H\nTSS-low activity",
        "G2 | I\nRNA activity",
        "G3 | L-like\nRNA activity",
        "G3 | L-like\nTSS-low activity",
        "G4 | L\nRNA activity",
        "G4 | L\nTSS-low activity"
      )
    ),
    gene_layer = paste0(gene_label, " | ", if_else(modality == "RNA", "RNA", "TSS-low")),
    os_sig = coalesce(as.logical(os_selected_any), FALSE) | (!is.na(os_min_p) & os_min_p < 0.05),
    drug_sig = coalesce(as.logical(feature_available), FALSE) & !is.na(depmap_fdr) & depmap_fdr < 0.25,
    immune_sig = (!is.na(til_activity_fdr) & til_activity_fdr < 0.10) | (!is.na(tmb_activity_fdr) & tmb_activity_fdr < 0.10),
    os_direction = case_when(
      !os_sig ~ "Not significant",
      os_best_activity_beta < 0 ~ "Protective",
      os_best_activity_beta > 0 ~ "Adverse",
      TRUE ~ "Not significant"
    ),
    depmap_activity_direction = case_when(
      !coalesce(as.logical(feature_available), FALSE) ~ "Not prioritized",
      depmap_activity_rho < 0 ~ "Sensitive",
      depmap_activity_rho > 0 ~ "Resistant",
      TRUE ~ "Not prioritized"
    ),
    drug_class_clean = clean_drug_class(drug_target_class, drug_clean),
    immune_quadrant_clean = case_when(
      til_activity_rho > 0 & tmb_activity_rho > 0 ~ "TIL-high/TMB-high",
      til_activity_rho < 0 & tmb_activity_rho > 0 ~ "TMB-shifted",
      til_activity_rho < 0 & tmb_activity_rho < 0 ~ "Immune-cold/TMB-low",
      TRUE ~ "Mixed/other"
    )
  )

write_csv(bridge, file.path(table_dir, "group_defined_marker_domain_input_tss_low_activity.csv"))

# -----------------------------------------------------------------------------
# OS domain-specific analysis
# -----------------------------------------------------------------------------
os_long <- read_csv(
  file.path(bridge_dir, "os_gene_activity_bridge_tss_meth_inverted.csv"),
  show_col_types = FALSE
) %>%
  filter(group_definition_tier %in% c("Group-defining", "Supportive")) %>%
  mutate(
    activity_dominant_group = factor(activity_dominant_group, levels = group_order),
    layer_story = short_layer(activity_layer),
    program = paste0(as.character(activity_dominant_group), "\n", layer_story),
    program = factor(program, levels = levels(bridge$program)),
    endpoint_label = case_when(
      endpoint == "overall_os" ~ "Overall OS",
      endpoint == "os_5y" ~ "5-year OS",
      endpoint == "os_10y" ~ "10-year OS",
      TRUE ~ endpoint_label
    ),
    endpoint_label = factor(endpoint_label, levels = c("Overall OS", "5-year OS", "10-year OS")),
    os_direction = case_when(
      p < 0.05 & activity_beta < 0 ~ "Protective",
      p < 0.05 & activity_beta > 0 ~ "Adverse",
      TRUE ~ "Not significant"
    ),
    selected = coalesce(as.logical(selected), FALSE),
    os_label = if_else(
      p < 0.05 | selected,
      paste0(gene_label, "\n", if_else(activity_beta < 0, "protective", "adverse"), " P=", format_p(p)),
      ""
    ),
    gene_layer = paste0(gene_label, " | ", if_else(modality == "RNA", "RNA", "TSS-low"))
  )

write_csv(os_long, file.path(table_dir, "OS_group_defined_marker_endpoint_results.csv"))

os_rows <- os_long %>%
  group_by(program, gene_layer) %>%
  summarise(best_p = min(p, na.rm = TRUE), .groups = "drop") %>%
  arrange(program, best_p, gene_layer) %>%
  mutate(row_id = factor(gene_layer, levels = rev(unique(gene_layer)))) %>%
  select(program, gene_layer, row_id)

os_heat_df <- os_long %>%
  left_join(os_rows, by = c("program", "gene_layer")) %>%
  mutate(
    beta_cap = pmax(pmin(activity_beta, 1.5), -1.5),
    point_label = case_when(
      selected ~ "Elastic-net selected",
      p < 0.05 ~ "Cox P<0.05",
      TRUE ~ "Not prioritized"
    ),
    point_label = factor(point_label, levels = c("Elastic-net selected", "Cox P<0.05", "Not prioritized"))
  )

p_os_heat <- ggplot(os_heat_df, aes(endpoint_label, row_id)) +
  geom_tile(aes(fill = beta_cap), colour = "white", linewidth = 0.55, width = 0.94, height = 0.88) +
  geom_point(
    data = os_heat_df %>% filter(selected | p < 0.05),
    aes(shape = point_label, size = os_neglog10_p),
    fill = "white", colour = ink, stroke = 0.45
  ) +
  facet_grid(program ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#4EA5F0", mid = "white", high = "#F4A259", midpoint = 0,
    limits = c(-1.5, 1.5), oob = squish, name = "Cox beta\n(activity)"
  ) +
  scale_shape_manual(values = c("Elastic-net selected" = 8, "Cox P<0.05" = 21, "Not prioritized" = 21), name = "OS evidence") +
  scale_size_continuous(range = c(2.2, 5.5), name = "-log10(P)") +
  labs(
    title = "A. OS evidence of EOBC group-defining biomarkers",
    subtitle = "Rows are EOBC group marker-layer features. Blue indicates protective activity; orange indicates adverse activity. TSS methylation is sign-inverted as TSS-low inferred activity.",
    x = NULL, y = NULL,
    caption = "OS is analyzed separately from DepMap and immune evidence. Symbols mark endpoint-level Cox P<0.05 or elastic-net selection."
  ) +
  theme_eobc(8.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", colour = ink, size = 9),
    axis.text.y = element_text(colour = ink, size = 7.2),
    strip.text.y = element_text(angle = 0, size = 8.2),
    legend.box = "horizontal"
  )

ggsave(file.path(plot_dir, "Figure_01A_OS_group_marker_endpoint_heatmap_R_v1.png"), p_os_heat, width = 9.8, height = 13.2, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_01A_OS_group_marker_endpoint_heatmap_R_v1.pdf"), p_os_heat, width = 9.8, height = 13.2, device = cairo_pdf, bg = "white")

p_os_volcano <- ggplot(os_long, aes(activity_beta, os_neglog10_p)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#EAF4FE", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FFF1E7", alpha = 0.45) +
  geom_hline(yintercept = -log10(0.05), colour = "#F87171", linewidth = 0.45, linetype = "longdash") +
  geom_vline(xintercept = 0, colour = "#8FA0B4", linewidth = 0.45, linetype = "longdash") +
  geom_point(
    aes(fill = activity_dominant_group, shape = layer_story, size = if_else(selected, 1.3, 1)),
    colour = ink, stroke = 0.32, alpha = 0.88
  ) +
  geom_text_repel(
    data = os_long %>% filter(p < 0.05 | selected),
    aes(label = gene_label, colour = activity_dominant_group),
    size = 2.65, fontface = "bold", max.overlaps = Inf,
    box.padding = 0.25, point.padding = 0.18, min.segment.length = 0,
    segment.alpha = 0.45, show.legend = FALSE
  ) +
  facet_wrap(~ endpoint_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = group_cols, name = "Dominant EOBC group") +
  scale_colour_manual(values = group_cols, guide = "none") +
  scale_shape_manual(values = c("RNA activity" = 21, "TSS-low activity" = 22), name = "Feature layer") +
  scale_size_identity() +
  labs(
    title = "B. OS volcano landscape restricted to EOBC group-defining markers",
    subtitle = "This plot is intentionally OS-only: each point is one marker-layer feature in one endpoint.",
    x = "Cox beta of marker activity score",
    y = "-log10(Cox P)"
  ) +
  theme_eobc(9) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, colour = ink, size = 3.2)),
    shape = guide_legend(override.aes = list(fill = "#FFFFFF", colour = ink, size = 3.2))
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(plot_dir, "Figure_01B_OS_group_marker_volcano_R_v1.png"), p_os_volcano, width = 14.8, height = 5.2, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_01B_OS_group_marker_volcano_R_v1.pdf"), p_os_volcano, width = 14.8, height = 5.2, device = cairo_pdf, bg = "white")

os_summary <- os_long %>%
  group_by(activity_dominant_group, layer_story, endpoint_label) %>%
  summarise(
    n_markers = n(),
    n_protective = sum(p < 0.05 & activity_beta < 0, na.rm = TRUE),
    n_adverse = sum(p < 0.05 & activity_beta > 0, na.rm = TRUE),
    n_selected = sum(selected, na.rm = TRUE),
    top_marker = gene_label[which.min(p)],
    top_beta = activity_beta[which.min(p)],
    top_p = min(p, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(os_summary, file.path(table_dir, "OS_group_program_summary.csv"))

# -----------------------------------------------------------------------------
# DepMap domain-specific analysis
# -----------------------------------------------------------------------------
depmap_df <- bridge %>%
  mutate(
    depmap_y = if_else(!is.na(depmap_fdr), -log10(pmax(depmap_fdr, 1e-300)), NA_real_),
    depmap_priority = case_when(
      !feature_available ~ "Feature unavailable",
      depmap_fdr < 0.10 ~ "FDR < 0.10",
      depmap_fdr < 0.25 ~ "FDR < 0.25",
      TRUE ~ "Exploratory"
    ),
    depmap_priority = factor(depmap_priority, levels = c("FDR < 0.10", "FDR < 0.25", "Exploratory", "Feature unavailable")),
    drug_label = if_else(feature_available, paste0(gene_label, "\n", drug_clean, "\n", drug_class_clean), paste0(gene_label, "\nnot tested")),
    direction_label = depmap_activity_direction
  )
write_csv(depmap_df, file.path(table_dir, "DepMap_group_defined_marker_best_drug_results.csv"))

drug_shape_values <- c(
  "DDR / checkpoint" = 21,
  "RTK/MAPK" = 24,
  "Mitotic / cell-cycle" = 22,
  "PI3K/AKT/mTOR" = 23,
  "Metabolic / WNT" = 25,
  "Cytotoxic / chemotherapy" = 24,
  "Other / exploratory" = 21
)

p_depmap <- ggplot(depmap_df %>% filter(feature_available), aes(depmap_activity_rho, depmap_y)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#EAF4FE", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#FFF1E7", alpha = 0.45) +
  geom_hline(yintercept = -log10(0.25), colour = "#FFB15C", linewidth = 0.40, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.10), colour = "#F87171", linewidth = 0.45, linetype = "longdash") +
  geom_vline(xintercept = 0, colour = "#8FA0B4", linewidth = 0.45, linetype = "longdash") +
  geom_point(
    aes(fill = Family6, shape = drug_class_clean, size = depmap_priority),
    colour = ink, stroke = 0.42, alpha = 0.92
  ) +
  geom_text_repel(
    aes(label = drug_label, colour = Family6),
    size = 2.25, fontface = "bold", lineheight = 0.82, max.overlaps = Inf,
    box.padding = 0.28, point.padding = 0.16, min.segment.length = 0,
    force = 1.4, max.time = 2,
    segment.alpha = 0.35, show.legend = FALSE
  ) +
  facet_wrap(~ program, ncol = 2, scales = "free") +
  scale_fill_manual(values = family_cols, name = "Biomarker family") +
  scale_colour_manual(values = family_cols, guide = "none") +
  scale_shape_manual(values = drug_shape_values, name = "Drug class") +
  scale_size_manual(values = c("FDR < 0.10" = 4.6, "FDR < 0.25" = 3.7, "Exploratory" = 2.5, "Feature unavailable" = 1.6), name = "Evidence tier") +
  labs(
    title = "C. DepMap drug-response evidence of EOBC group-defining biomarkers",
    subtitle = "Each marker-layer is tested independently; x-axis is signed activity association with drug response. Left = high marker activity aligns with sensitivity; right = resistance.",
    x = "Signed DepMap association with drug response",
    y = "-log10(best FDR across AUC / IC50)",
    caption = "TSS methylation features are sign-inverted before plotting. Labels show biomarker, best associated drug, and drug class."
  ) +
  theme_eobc(8.4) +
  coord_cartesian(clip = "off") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8.1),
    legend.position = "bottom",
    plot.margin = margin(10, 30, 10, 30)
  )

ggsave(file.path(plot_dir, "Figure_02A_DepMap_group_marker_landscape_R_v1.png"), p_depmap, width = 16.2, height = 13.2, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_02A_DepMap_group_marker_landscape_R_v1.pdf"), p_depmap, width = 16.2, height = 13.2, device = cairo_pdf, bg = "white")

depmap_summary <- depmap_df %>%
  filter(feature_available) %>%
  mutate(
    evidence_tier = case_when(
      depmap_fdr < 0.10 ~ "FDR < 0.10",
      depmap_fdr < 0.25 ~ "FDR < 0.25",
      TRUE ~ "Exploratory"
    )
  ) %>%
  group_by(activity_dominant_group, layer_story, direction_label, drug_class_clean) %>%
  summarise(
    n = n(),
    n_fdr_025 = sum(depmap_fdr < 0.25, na.rm = TRUE),
    top_marker = gene_label[which.min(depmap_fdr)],
    top_drug = drug_clean[which.min(depmap_fdr)],
    top_fdr = min(depmap_fdr, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(depmap_summary, file.path(table_dir, "DepMap_group_program_summary.csv"))

p_depmap_summary <- depmap_summary %>%
  mutate(
    program = paste0(as.character(activity_dominant_group), "\n", layer_story),
    program = factor(program, levels = levels(bridge$program)),
    direction_label = factor(direction_label, levels = c("Sensitive", "Resistant", "Not prioritized"))
  ) %>%
  filter(direction_label != "Not prioritized") %>%
  ggplot(aes(drug_class_clean, program)) +
  geom_tile(fill = "#F8FAFC", colour = "#E6EEF7", linewidth = 0.45, width = 0.92, height = 0.82) +
  geom_point(
    data = depmap_summary %>%
      mutate(
        program = paste0(as.character(activity_dominant_group), "\n", layer_story),
        program = factor(program, levels = levels(bridge$program)),
        direction_label = factor(direction_label, levels = c("Sensitive", "Resistant", "Not prioritized"))
      ) %>%
      filter(direction_label != "Not prioritized", n_fdr_025 > 0),
    aes(size = n_fdr_025, fill = direction_label),
    shape = 21, colour = ink, stroke = 0.35, alpha = 0.95
  ) +
  geom_text(
    data = depmap_summary %>%
      mutate(
        program = paste0(as.character(activity_dominant_group), "\n", layer_story),
        program = factor(program, levels = levels(bridge$program)),
        direction_label = factor(direction_label, levels = c("Sensitive", "Resistant", "Not prioritized"))
      ) %>%
      filter(direction_label != "Not prioritized", n_fdr_025 > 0),
    aes(label = as.character(n_fdr_025)),
    colour = "white", fontface = "bold", size = 3
  ) +
  scale_fill_manual(values = drug_cols, name = "Direction") +
  scale_size_continuous(range = c(4, 10), breaks = c(1, 3, 6), name = "FDR<0.25\nmarker count") +
  labs(
    title = "D. Drug-class burden among group-defining biomarkers",
    subtitle = "Counts summarize FDR-prioritized DepMap associations by EOBC program and drug class.",
    x = NULL, y = NULL
  ) +
  theme_eobc(8.8) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(file.path(plot_dir, "Figure_02B_DepMap_group_marker_drug_class_summary_R_v1.png"), p_depmap_summary, width = 13.2, height = 5.8, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_group_marker_drug_class_summary_R_v1.pdf"), p_depmap_summary, width = 13.2, height = 5.8, device = cairo_pdf, bg = "white")

# -----------------------------------------------------------------------------
# Immune domain-specific analysis
# -----------------------------------------------------------------------------
immune_df <- bridge %>%
  mutate(
    immune_best_fdr = pmin(coalesce(til_activity_fdr, 1), coalesce(tmb_activity_fdr, 1), na.rm = TRUE),
    immune_priority = case_when(
      immune_best_fdr < 0.05 ~ "FDR < 0.05",
      immune_best_fdr < 0.10 ~ "FDR < 0.10",
      TRUE ~ "Exploratory"
    ),
    immune_priority = factor(immune_priority, levels = c("FDR < 0.05", "FDR < 0.10", "Exploratory")),
    immune_label = paste0(
      gene_label, "\n",
      "TIL ", sprintf("%.2f", til_activity_rho), " / TMB ", sprintf("%.2f", tmb_activity_rho)
    )
  )
write_csv(immune_df, file.path(table_dir, "Immune_group_defined_marker_til_tmb_results.csv"))

p_immune <- ggplot(immune_df, aes(til_activity_rho, tmb_activity_rho)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "#EAF4FE", alpha = 0.65) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#E8FAF7", alpha = 0.75) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#F0FAED", alpha = 0.75) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "#F7F9FC", alpha = 0.75) +
  geom_hline(yintercept = 0, colour = "#8FA0B4", linewidth = 0.45, linetype = "longdash") +
  geom_vline(xintercept = 0, colour = "#8FA0B4", linewidth = 0.45, linetype = "longdash") +
  geom_point(
    aes(fill = Family6, shape = layer_story, size = immune_priority),
    colour = ink, stroke = 0.42, alpha = 0.92
  ) +
  geom_text_repel(
    aes(label = immune_label, colour = Family6),
    size = 2.25, fontface = "bold", lineheight = 0.82, max.overlaps = Inf,
    box.padding = 0.28, point.padding = 0.16, min.segment.length = 0,
    force = 1.4, max.time = 2,
    segment.alpha = 0.35, show.legend = FALSE
  ) +
  facet_wrap(~ program, ncol = 2, scales = "free") +
  scale_fill_manual(values = family_cols, name = "Biomarker family") +
  scale_colour_manual(values = family_cols, guide = "none") +
  scale_shape_manual(values = c("RNA activity" = 21, "TSS-low activity" = 22), name = "Feature layer") +
  scale_size_manual(values = c("FDR < 0.05" = 4.6, "FDR < 0.10" = 3.7, "Exploratory" = 2.4), name = "Evidence tier") +
  labs(
    title = "E. Immune TIL/TMB evidence of EOBC group-defining biomarkers",
    subtitle = "Each facet is an EOBC state program. Axes show Spearman correlations with TIL and TMB; methylation is plotted as TSS-low inferred activity.",
    x = "Spearman rho with TIL score",
    y = "Spearman rho with TMB"
  ) +
  theme_eobc(8.4) +
  coord_cartesian(clip = "off") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 8.1),
    legend.position = "bottom",
    plot.margin = margin(10, 30, 10, 30)
  )

ggsave(file.path(plot_dir, "Figure_03A_Immune_group_marker_TIL_TMB_landscape_R_v1.png"), p_immune, width = 16.2, height = 13.2, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_03A_Immune_group_marker_TIL_TMB_landscape_R_v1.pdf"), p_immune, width = 16.2, height = 13.2, device = cairo_pdf, bg = "white")

immune_summary <- immune_df %>%
  group_by(activity_dominant_group, layer_story, immune_quadrant_clean) %>%
  summarise(
    n = n(),
    n_fdr_010 = sum(immune_best_fdr < 0.10, na.rm = TRUE),
    top_marker = gene_label[which.min(immune_best_fdr)],
    top_til = til_activity_rho[which.min(immune_best_fdr)],
    top_tmb = tmb_activity_rho[which.min(immune_best_fdr)],
    top_fdr = min(immune_best_fdr, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(immune_summary, file.path(table_dir, "Immune_group_program_summary.csv"))

p_immune_summary <- immune_summary %>%
  mutate(
    program = paste0(as.character(activity_dominant_group), "\n", layer_story),
    program = factor(program, levels = levels(bridge$program)),
    immune_quadrant_clean = factor(immune_quadrant_clean, levels = names(immune_cols))
  ) %>%
  ggplot(aes(immune_quadrant_clean, program)) +
  geom_tile(fill = "#F8FAFC", colour = "#E6EEF7", linewidth = 0.45, width = 0.92, height = 0.82) +
  geom_point(
    data = immune_summary %>%
      mutate(
        program = paste0(as.character(activity_dominant_group), "\n", layer_story),
        program = factor(program, levels = levels(bridge$program)),
        immune_quadrant_clean = factor(immune_quadrant_clean, levels = names(immune_cols))
      ) %>%
      filter(n_fdr_010 > 0),
    aes(size = n_fdr_010, fill = immune_quadrant_clean),
    shape = 21, colour = ink, stroke = 0.35, alpha = 0.95
  ) +
  geom_text(
    data = immune_summary %>%
      mutate(
        program = paste0(as.character(activity_dominant_group), "\n", layer_story),
        program = factor(program, levels = levels(bridge$program)),
        immune_quadrant_clean = factor(immune_quadrant_clean, levels = names(immune_cols))
      ) %>%
      filter(n_fdr_010 > 0),
    aes(label = as.character(n_fdr_010)),
    colour = "white", fontface = "bold", size = 3
  ) +
  scale_fill_manual(values = immune_cols, name = "Immune context") +
  scale_size_continuous(range = c(4, 10), breaks = c(1, 3, 6), name = "FDR<0.10\nmarker count") +
  labs(
    title = "F. Immune-context burden among group-defining biomarkers",
    subtitle = "Counts summarize FDR-prioritized TIL/TMB associations by EOBC program.",
    x = NULL, y = NULL
  ) +
  theme_eobc(8.8) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )

ggsave(file.path(plot_dir, "Figure_03B_Immune_group_marker_context_summary_R_v1.png"), p_immune_summary, width = 12.8, height = 5.8, dpi = 330, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_group_marker_context_summary_R_v1.pdf"), p_immune_summary, width = 12.8, height = 5.8, device = cairo_pdf, bg = "white")

readme <- c(
  "# EOBC group-marker domain-specific analyses v1",
  "",
  "This folder intentionally separates the downstream analyses instead of making one final integrated conclusion figure.",
  "",
  "Core principle:",
  "- Group-defining/supportive biomarker-layer features are fixed first.",
  "- OS, DepMap drug response, and immune/TIL-TMB evidence are analyzed in separate domain-specific figures.",
  "- TSS/promoter methylation is sign-inverted and interpreted as TSS-low inferred activity, because high TSS methylation generally suppresses RNA expression.",
  "",
  "Figures:",
  "- Figure_01A/01B: OS-only endpoint heatmap and volcano landscape.",
  "- Figure_02A/02B: DepMap-only drug-response landscape and drug-class burden.",
  "- Figure_03A/03B: Immune-only TIL/TMB landscape and immune-context burden.",
  "",
  "Tables:",
  "- OS_group_defined_marker_endpoint_results.csv",
  "- DepMap_group_defined_marker_best_drug_results.csv",
  "- Immune_group_defined_marker_til_tmb_results.csv",
  "- Domain-specific summary tables."
)
writeLines(readme, file.path(out_dir, "README_group_marker_domain_specific_R_v1.md"))

message("Saved group-marker domain-specific outputs to: ", out_dir)
