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

set.seed(20260515)

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

base_dir <- env_path("EOBC_BIOMARKER_ROOT")
analysis_dir <- file.path(base_dir, "final_analysis", "group_biomarker_landscape_r_v1")
table_dir <- file.path(analysis_dir, "tables")
plot_dir <- file.path(analysis_dir, "plots")
final_fig_dir <- file.path(base_dir, "final_analysis", "final_fig")
final_revision_dir <- file.path(final_fig_dir, "final_paper_figures_20260522_v2")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_revision_dir, recursive = TRUE, showWarnings = FALSE)

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

theme_eobc <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
      panel.grid.minor = element_line(colour = alpha(grid_col, 0.45), linewidth = 0.20),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = muted),
      strip.text = element_text(face = "bold", colour = ink),
      strip.background = element_rect(fill = "#EEF4FB", colour = "#D6E0EA", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", colour = ink),
      legend.text = element_text(colour = muted),
      plot.title = element_text(face = "bold", colour = ink, size = rel(1.45)),
      plot.subtitle = element_text(colour = muted, size = rel(0.95)),
      plot.caption = element_text(colour = muted, size = rel(0.72), hjust = 0),
      plot.title.position = "plot"
    )
}

save_plot_all <- function(plot, basename, width, height, dpi = 360) {
  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  for (dest_dir in unique(c(plot_dir, final_fig_dir, final_revision_dir))) {
    ggsave(file.path(dest_dir, paste0(basename, ".png")),
           plot, width = width, height = height, dpi = dpi, bg = "white")
    ggsave(file.path(dest_dir, paste0(basename, ".pdf")),
           plot, width = width, height = height, bg = "white", device = pdf_device)
  }
}

biomarker_long <- read_csv(file.path(table_dir, "group_biomarker_long_values.csv"),
                           show_col_types = FALSE) %>%
  mutate(
    group_label = factor(group_label, levels = group_order),
    group_index = as.numeric(group_label),
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")),
    layer_short = if_else(Layer == "Methylation", "M", "R"),
    gene_short = paste0(gene, " [", layer_short, "]")
  )

group_means <- read_csv(file.path(table_dir, "group_biomarker_group_means.csv"),
                        show_col_types = FALSE) %>%
  mutate(
    group_label = factor(group_label, levels = group_order),
    group_index = as.numeric(group_label),
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")),
    layer_short = if_else(Layer == "Methylation", "M", "R"),
    gene_short = paste0(gene, " [", layer_short, "]")
  )

dominant <- read_csv(file.path(table_dir, "group_biomarker_dominant_group_summary.csv"),
                     show_col_types = FALSE) %>%
  mutate(
    dominant_group = factor(dominant_group, levels = group_order),
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")),
    layer_short = if_else(Layer == "Methylation", "M", "R"),
    gene_short = paste0(gene, " [", layer_short, "]"),
    neg_log10_fdr = -log10(pmax(kw_fdr, 1e-12)),
    evidence_tier = case_when(
      kw_fdr < 0.05 & contrast_z >= 1.0 ~ "Core",
      kw_fdr < 0.10 & contrast_z >= 0.75 ~ "Supportive",
      kw_fdr < 0.10 ~ "Weak",
      TRUE ~ "Exploratory"
    ),
    rank_score = contrast_z * pmin(neg_log10_fdr, 12)
  )

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

trend_stats <- biomarker_long %>%
  group_by(modality, modality_label, gene, gene_short, Layer, Family6, target_label) %>%
  summarise(
    n = n(),
    pseudotime_rho = safe_cor(value_z, Pseudotime)$rho,
    pseudotime_p = safe_cor(value_z, Pseudotime)$p,
    group_index_rho = safe_cor(value_z, group_index)$rho,
    group_index_p = safe_cor(value_z, group_index)$p,
    .groups = "drop"
  ) %>%
  group_by(modality) %>%
  mutate(
    pseudotime_fdr = p.adjust(pseudotime_p, method = "BH"),
    group_index_fdr = p.adjust(group_index_p, method = "BH")
  ) %>%
  ungroup() %>%
  left_join(
    dominant %>%
      select(modality, gene, dominant_group, contrast_z, kw_fdr, evidence_tier, rank_score),
    by = c("modality", "gene")
  )

write_csv(trend_stats, file.path(table_dir, "group_biomarker_ordered_trend_scores.csv"))

## Biomarker-only PCA: are the four EOBC states visible from the selected markers alone?
pca_one_modality <- function(df, mod_label) {
  wide <- df %>%
    filter(modality_label == mod_label) %>%
    select(Sample, group_label, gene, value_z) %>%
    distinct() %>%
    pivot_wider(names_from = gene, values_from = value_z) %>%
    arrange(group_label, Sample)
  meta <- wide %>% select(Sample, group_label)
  mat <- wide %>% select(-Sample, -group_label) %>% as.matrix()
  mat[!is.finite(mat)] <- 0
  pca <- prcomp(mat, center = FALSE, scale. = FALSE)
  var_exp <- (pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100
  scores <- as_tibble(pca$x[, 1:2]) %>%
    bind_cols(meta) %>%
    mutate(
      modality_label = mod_label,
      pc1_lab = sprintf("PC1 (%.1f%%)", var_exp[1]),
      pc2_lab = sprintf("PC2 (%.1f%%)", var_exp[2])
    )
  attr(scores, "var_exp") <- var_exp
  scores
}

pca_scores <- bind_rows(
  pca_one_modality(biomarker_long, "RNA feature"),
  pca_one_modality(biomarker_long, "Methylation feature")
) %>%
  mutate(
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")),
    group_label = factor(group_label, levels = group_order)
  )

write_csv(pca_scores, file.path(table_dir, "group_biomarker_pca_scores.csv"))

fig_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = group_label, fill = group_label)) +
  stat_ellipse(type = "norm", geom = "polygon", alpha = 0.08, colour = NA, level = 0.68) +
  stat_ellipse(type = "norm", linewidth = 0.55, alpha = 0.75, level = 0.68) +
  geom_point(size = 2.25, alpha = 0.85, stroke = 0) +
  facet_wrap(~ modality_label, ncol = 2, scales = "free") +
  scale_colour_manual(values = group_cols, name = "EOBC group") +
  scale_fill_manual(values = group_cols, guide = "none") +
  labs(
    title = "F. Biomarker-only sample map recapitulates EOBC group structure",
    subtitle = "PCA was calculated independently for RNA and methylation features using only the selected biomarker panel.",
    x = "PC1",
    y = "PC2"
  ) +
  theme_eobc(base_size = 9.4) +
  theme(legend.position = "bottom")

save_plot_all(fig_pca, "Figure_06_EOBC_biomarker_only_PCA_R_v1", 10.8, 5.3)

## State-specific biomarker trajectories.
trajectory_markers <- dominant %>%
  filter(kw_fdr < 0.10, contrast_z >= 0.75) %>%
  group_by(modality_label, dominant_group) %>%
  arrange(desc(rank_score), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  select(modality, modality_label, gene, dominant_group, rank_score, kw_fdr, neg_log10_fdr, contrast_z)

traj_df <- group_means %>%
  inner_join(trajectory_markers,
             by = c("modality", "modality_label", "gene")) %>%
  mutate(
    dominant_group = factor(dominant_group, levels = group_order),
    facet_lab = paste(modality_label, dominant_group, sep = "\n"),
    facet_lab = factor(
      facet_lab,
      levels = as.vector(t(outer(c("RNA feature", "Methylation feature"), group_order, paste, sep = "\n")))
    ),
    is_peak = group_label == dominant_group
  )

traj_label_df <- traj_df %>%
  filter(is_peak) %>%
  group_by(facet_lab) %>%
  arrange(desc(rank_score), .by_group = TRUE) %>%
  slice_head(n = 6) %>%
  ungroup()

fig_traj <- ggplot(traj_df, aes(x = group_label, y = mean_z, group = gene, colour = Family6)) +
  geom_hline(yintercept = 0, colour = "#94A3B8", linetype = "dashed", linewidth = 0.38) +
  geom_line(linewidth = 0.85, alpha = 0.70) +
  geom_point(aes(size = if_else(is_peak, neg_log10_fdr, 1.2)),
             shape = 21, fill = "white", stroke = 0.65, alpha = 0.95) +
  geom_text_repel(
    data = traj_label_df,
    aes(label = gene_short),
    size = 2.35,
    fontface = "bold",
    box.padding = 0.18,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    segment.size = 0.16,
    max.overlaps = Inf,
    seed = 20260515,
    show.legend = FALSE
  ) +
  facet_wrap(~ facet_lab, ncol = 4) +
  scale_colour_manual(values = family_cols, name = "Biomarker family") +
  scale_size_continuous(range = c(1.4, 5.5), breaks = c(1, 3, 6, 9, 12),
                        name = "Peak -log10\nKW FDR") +
  coord_cartesian(ylim = c(-1.35, 1.35), clip = "off") +
  labs(
    title = "G. Group-wise biomarker trajectories reveal EOBC state-specific programs",
    subtitle = "Each line follows group mean z-score across G1-G4; RNA uses expression z-score and methylation uses raw TSS beta-value z-score without direction flipping.",
    x = "EOBC group ordered by group number",
    y = "Group mean z-score",
    caption = "Markers shown here pass Kruskal-Wallis FDR < 0.10 and dominant-vs-suppressed contrast >= 0.75; point size at the peak encodes -log10(KW FDR)."
  ) +
  theme_eobc(base_size = 8.0) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    strip.text = element_text(size = 8.4),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

save_plot_all(fig_traj, "Figure_07_EOBC_group_biomarker_trajectory_R_v1", 14.6, 8.7)

## Family-level state trajectories.
family_program_scores <- read_csv(file.path(table_dir, "group_biomarker_family_program_scores.csv"),
                                  show_col_types = FALSE) %>%
  mutate(
    group_label = factor(group_label, levels = group_order),
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature"))
  )

family_label_df <- family_program_scores %>%
  group_by(modality_label, Family6) %>%
  filter(abs(family_mean_z) == max(abs(family_mean_z), na.rm = TRUE)) %>%
  slice_head(n = 1) %>%
  ungroup()

fig_family_traj <- ggplot(family_program_scores,
                          aes(x = group_label, y = family_mean_z, group = Family6, colour = Family6)) +
  geom_hline(yintercept = 0, colour = "#94A3B8", linetype = "dashed", linewidth = 0.40) +
  geom_line(linewidth = 1.05, alpha = 0.82) +
  geom_point(aes(size = n_significant_markers), shape = 21, fill = "white", stroke = 0.75) +
  geom_text_repel(
    data = family_label_df,
    aes(label = Family6),
    size = 2.55,
    fontface = "bold",
    box.padding = 0.18,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    segment.size = 0.16,
    max.overlaps = Inf,
    seed = 20260515,
    show.legend = FALSE
  ) +
  facet_wrap(~ modality_label, ncol = 2) +
  scale_colour_manual(values = family_cols, name = "Biomarker family") +
  scale_size_continuous(range = c(1.6, 6.2), breaks = c(0, 2, 4, 6),
                        name = "Dominant\nmarker count") +
  labs(
    title = "H. Biological-family trajectories across the four EOBC states",
    subtitle = "Line height is the average group mean z-score within each biomarker family; methylation uses raw TSS beta-value z-score.",
    x = "EOBC group ordered by group number",
    y = "Family mean z-score",
    caption = "Point size counts state-dominant markers passing Kruskal-Wallis FDR < 0.10 and contrast >= 0.75 in that family/state."
  ) +
  theme_eobc(base_size = 9.2) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    legend.box = "horizontal"
  )

save_plot_all(fig_family_traj, "Figure_08_EOBC_family_trajectory_R_v1", 11.8, 6.1)

## State-specific dominance landscape, easier to read than the all-in-one landscape.
dom_state_df <- dominant %>%
  mutate(
    neg_log10_fdr_cap = pmin(neg_log10_fdr, 12),
    label = if_else(kw_fdr < 0.10 & contrast_z >= 0.75, gene_short, NA_character_)
  )

fig_dom_state <- ggplot(dom_state_df,
                        aes(x = contrast_z, y = neg_log10_fdr_cap, fill = Family6)) +
  geom_hline(yintercept = -log10(0.10), colour = "#FDBA74", linetype = "dotted", linewidth = 0.52) +
  geom_hline(yintercept = -log10(0.05), colour = "#F87171", linetype = "longdash", linewidth = 0.52) +
  geom_point(
    aes(size = abs(dominant_mean_z)),
    shape = 21,
    colour = ink,
    stroke = 0.55,
    alpha = 0.93
  ) +
  geom_text_repel(
    aes(label = label, colour = Family6),
    size = 2.15,
    fontface = "bold",
    box.padding = 0.16,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    segment.size = 0.15,
    force = 1.5,
    max.overlaps = Inf,
    seed = 20260515,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  facet_grid(modality_label ~ dominant_group) +
  scale_fill_manual(values = family_cols, name = "Biomarker family") +
  scale_colour_manual(values = family_cols, guide = "none") +
  scale_size_continuous(range = c(1.8, 5.5), name = "|Dominant\nmean z|") +
  scale_y_continuous(limits = c(0, 12.6), breaks = c(0, 3, 6, 9, 12),
                     expand = expansion(mult = c(0.03, 0.08))) +
  coord_cartesian(clip = "off") +
  labs(
    title = "I. State-resolved dominance landscape of EOBC biomarkers",
    subtitle = "Splitting the dominance plot by highest EOBC group makes each state-defining marker program explicit.",
    x = "Dominant-vs-suppressed group contrast",
    y = "-log10(Kruskal-Wallis FDR), capped at 12"
  ) +
  theme_eobc(base_size = 7.8) +
  theme(
    strip.text = element_text(size = 7.8),
    legend.box = "horizontal",
    plot.margin = margin(10, 16, 10, 10)
  )

save_plot_all(fig_dom_state, "Figure_09_EOBC_state_resolved_dominance_landscape_R_v1", 14.5, 7.4)

message("Trajectory-focused group biomarker figures complete.")
message("Plots written to: ", plot_dir)
message("Final figure copies written to: ", final_fig_dir)
message("Tables written to: ", table_dir)
