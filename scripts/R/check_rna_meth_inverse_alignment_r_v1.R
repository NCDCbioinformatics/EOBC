suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

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
final_analysis_dir <- env_path("EOBC_FINAL_ANALYSIS_DIR", file.path(biomarker_root, "final_analysis"))
out_dir <- file.path(final_analysis_dir, "group_biomarker_landscape_r_v1")
table_dir <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#111827"
muted <- "#64748B"
grid_col <- "#D8E3EF"
family_cols <- c(
  "Immune" = "#4EA5F0",
  "Repair" = "#47C56B",
  "Glycolysis / TCA" = "#E6BC18",
  "Fatty acid" = "#F4A259",
  "Kinase signaling" = "#9B6AE8",
  "Hormone signaling" = "#9AAABC"
)
group_order <- c("G1 | H", "G2 | I", "G3 | L-like", "G4 | L")

theme_eobc <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = ink),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.35),
      panel.grid.minor = element_line(colour = scales::alpha(grid_col, 0.45), linewidth = 0.20),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = muted),
      strip.text = element_text(face = "bold", colour = ink),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", colour = ink),
      legend.text = element_text(colour = muted),
      plot.title = element_text(face = "bold", colour = ink, size = rel(1.35)),
      plot.subtitle = element_text(colour = muted, size = rel(0.9)),
      plot.title.position = "plot"
    )
}

long_values <- read_csv(file.path(table_dir, "group_biomarker_long_values.csv"), show_col_types = FALSE)
group_means <- read_csv(file.path(table_dir, "group_biomarker_group_means.csv"), show_col_types = FALSE)

sample_inverse <- long_values %>%
  select(Sample, gene, gene_label, Family6, modality, value_z) %>%
  filter(modality %in% c("RNA", "METH")) %>%
  pivot_wider(names_from = modality, values_from = value_z) %>%
  group_by(gene, gene_label, Family6) %>%
  summarise(
    n_overlap = sum(is.finite(RNA) & is.finite(METH)),
    sample_spearman_rho = if_else(
      n_overlap >= 10,
      suppressWarnings(cor(RNA, METH, method = "spearman", use = "complete.obs")),
      NA_real_
    ),
    sample_pearson_r = if_else(
      n_overlap >= 10,
      suppressWarnings(cor(RNA, METH, method = "pearson", use = "complete.obs")),
      NA_real_
    ),
    .groups = "drop"
  )

profile_wide <- group_means %>%
  select(gene, gene_label, Family6, modality, group_label, mean_z) %>%
  filter(modality %in% c("RNA", "METH")) %>%
  pivot_wider(names_from = modality, values_from = mean_z)

group_inverse <- profile_wide %>%
  group_by(gene, gene_label, Family6) %>%
  summarise(
    n_groups = sum(is.finite(RNA) & is.finite(METH)),
    group_spearman_rho = if_else(
      n_groups >= 4,
      suppressWarnings(cor(RNA, METH, method = "spearman", use = "complete.obs")),
      NA_real_
    ),
    group_pearson_r = if_else(
      n_groups >= 4,
      suppressWarnings(cor(RNA, METH, method = "pearson", use = "complete.obs")),
      NA_real_
    ),
    .groups = "drop"
  )

peak_alignment <- profile_wide %>%
  mutate(group_label = factor(group_label, levels = group_order)) %>%
  group_by(gene, gene_label, Family6) %>%
  summarise(
    rna_peak_group = as.character(group_label[which.max(RNA)]),
    meth_high_group_raw = as.character(group_label[which.max(METH)]),
    meth_low_group_inferred_activity_peak = as.character(group_label[which.min(METH)]),
    rna_peak_z = max(RNA, na.rm = TRUE),
    meth_high_z = max(METH, na.rm = TRUE),
    meth_low_z = min(METH, na.rm = TRUE),
    rna_matches_meth_low_activity = rna_peak_group == meth_low_group_inferred_activity_peak,
    rna_matches_raw_meth_high = rna_peak_group == meth_high_group_raw,
    .groups = "drop"
  )

alignment <- sample_inverse %>%
  left_join(group_inverse, by = c("gene", "gene_label", "Family6")) %>%
  left_join(peak_alignment, by = c("gene", "gene_label", "Family6")) %>%
  mutate(
    sample_inverse_flag = sample_spearman_rho < 0,
    group_inverse_flag = group_spearman_rho < 0,
    inverse_consistency_class = case_when(
      sample_inverse_flag & group_inverse_flag & rna_matches_meth_low_activity ~ "Strong inverse / activity-consistent",
      group_inverse_flag & rna_matches_meth_low_activity ~ "Group inverse / activity-consistent",
      sample_inverse_flag | group_inverse_flag ~ "Partial inverse",
      TRUE ~ "Concordant or context-specific"
    )
  ) %>%
  arrange(desc(group_inverse_flag), desc(sample_inverse_flag), group_spearman_rho)

write_csv(alignment, file.path(table_dir, "rna_meth_tss_inverse_alignment.csv"))

summary_tbl <- alignment %>%
  summarise(
    n_genes = n(),
    n_sample_inverse = sum(sample_inverse_flag, na.rm = TRUE),
    pct_sample_inverse = n_sample_inverse / n_genes,
    n_group_inverse = sum(group_inverse_flag, na.rm = TRUE),
    pct_group_inverse = n_group_inverse / n_genes,
    n_activity_peak_match = sum(rna_matches_meth_low_activity, na.rm = TRUE),
    pct_activity_peak_match = n_activity_peak_match / n_genes,
    median_sample_rho = median(sample_spearman_rho, na.rm = TRUE),
    median_group_rho = median(group_spearman_rho, na.rm = TRUE)
  )
write_csv(summary_tbl, file.path(table_dir, "rna_meth_tss_inverse_alignment_summary.csv"))

label_df <- alignment %>%
  filter(
    abs(group_spearman_rho) >= 0.8 |
      abs(sample_spearman_rho) >= 0.35 |
      inverse_consistency_class == "Strong inverse / activity-consistent"
  )

p1 <- ggplot(alignment, aes(sample_spearman_rho, group_spearman_rho)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "#E8F3FF", alpha = 0.55) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "#94A3B8", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#94A3B8", linewidth = 0.45) +
  geom_point(aes(fill = Family6, size = abs(group_pearson_r)), shape = 21, colour = ink, alpha = 0.88, stroke = 0.45) +
  geom_text_repel(
    data = label_df,
    aes(label = gene_label, colour = Family6),
    size = 3.0,
    box.padding = 0.35,
    point.padding = 0.22,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = family_cols, drop = FALSE) +
  scale_colour_manual(values = family_cols, drop = FALSE) +
  scale_size_continuous(range = c(2.5, 8), name = "|Group-level Pearson r|") +
  labs(
    title = "RNA-METH inverse sanity check for TSS methylation biomarkers",
    subtitle = "Negative correlation supports the expected TSS methylation-to-expression repression model.",
    x = "Sample-level Spearman rho: RNA z vs TSS methylation z",
    y = "Group-profile Spearman rho: RNA group mean z vs methylation group mean z",
    fill = "Biomarker family"
  ) +
  coord_cartesian(xlim = c(-0.8, 0.8), ylim = c(-1.05, 1.05), clip = "off") +
  theme_eobc(10)

peak_counts <- peak_alignment %>%
  mutate(
    rna_peak_group = factor(rna_peak_group, levels = group_order),
    meth_low_group_inferred_activity_peak = factor(meth_low_group_inferred_activity_peak, levels = group_order)
  ) %>%
  count(rna_peak_group, meth_low_group_inferred_activity_peak)

p2 <- ggplot(peak_counts, aes(rna_peak_group, meth_low_group_inferred_activity_peak)) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 1.1) +
  geom_text(aes(label = n), colour = ink, fontface = "bold", size = 4) +
  scale_fill_gradient(low = "#F8FAFC", high = "#4EA5F0", name = "Gene count") +
  labs(
    title = "RNA peak versus methylation-inferred activity peak",
    subtitle = "For TSS methylation, the biologically active methylation peak is the lowest methylation group.",
    x = "RNA-high group",
    y = "Methylation-low group (-METH activity peak)"
  ) +
  coord_fixed() +
  theme_eobc(10) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

combined <- p1 / p2 + plot_layout(heights = c(2.1, 1))

ggsave(file.path(plot_dir, "Figure_10_RNA_METH_TSS_inverse_alignment_R_v1.png"), combined, width = 10.5, height = 10.5, dpi = 320, bg = "white")
ggsave(file.path(plot_dir, "Figure_10_RNA_METH_TSS_inverse_alignment_R_v1.pdf"), combined, width = 10.5, height = 10.5, device = cairo_pdf, bg = "white")

message("Wrote RNA-METH inverse alignment checks to: ", out_dir)
