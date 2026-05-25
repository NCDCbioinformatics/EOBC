suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
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
input_dir <- file.path(base_dir, "00_inputs_detected")
panel_path <- file.path(
  base_dir,
  "09_selected_relaxed_union_family6_final_signedMeth",
  "tables",
  "selected_union_biomarker_genes_family6_v18.csv"
)
clinical_path <- file.path(input_dir, "total_sample_clinical_all.csv")
rna_path <- file.path(input_dir, "TPM_young.csv")
meth_path <- file.path(input_dir, "MET_young_batch_JW.csv")

out_dir <- file.path(base_dir, "final_analysis", "group_biomarker_landscape_r_v1")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

ink <- "#111827"
muted <- "#64748B"
grid_col <- "#D8E3EF"

group_order <- c("H", "I", "L_like", "L")
group_labels <- c(
  "H" = "G1 | H",
  "I" = "G2 | I",
  "L_like" = "G3 | L-like",
  "L" = "G4 | L"
)
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

pam50_cols <- c(
  "Basal" = "#1F78B4",
  "LumA" = "#A6CEE3",
  "LumB" = "#FDBF6F",
  "Her2" = "#E31A1C",
  "Normal" = "#B2DF8A",
  "Unknown" = "#CBD5E1"
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

clip_value <- function(x, lim = 2.5) pmax(pmin(x, lim), -lim)

zscore_by_gene <- function(df) {
  df %>%
    group_by(gene, modality) %>%
    mutate(
      value_z = as.numeric(scale(value)),
      value_z = if_else(is.finite(value_z), value_z, 0),
      value_z_clip = clip_value(value_z, 2.5)
    ) %>%
    ungroup()
}

read_omics_long <- function(path, genes, modality, transform = c("none", "log2")) {
  transform <- match.arg(transform)
  raw <- read_csv(path, show_col_types = FALSE)
  names(raw)[1] <- "gene"
  raw <- raw %>%
    mutate(gene = as.character(gene)) %>%
    filter(gene %in% genes) %>%
    distinct(gene, .keep_all = TRUE)

  raw %>%
    pivot_longer(-gene, names_to = "Sample", values_to = "value") %>%
    mutate(
      Sample = str_replace_all(Sample, '^"|"$', ""),
      value = suppressWarnings(as.numeric(value)),
      value = if (transform == "log2") log2(value + 1) else value,
      modality = modality,
      modality_label = if_else(modality == "RNA", "RNA feature", "Methylation feature")
    ) %>%
    filter(str_detect(Sample, "\\.T$"), is.finite(value))
}

safe_kw <- function(value, group) {
  keep <- is.finite(value) & !is.na(group)
  value <- value[keep]
  group <- droplevels(factor(group[keep]))
  if (length(unique(group)) < 2 || length(value) < 8) return(NA_real_)
  out <- try(kruskal.test(value ~ group)$p.value, silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}

safe_wilcox <- function(value, is_group) {
  keep <- is.finite(value) & !is.na(is_group)
  value <- value[keep]
  is_group <- is_group[keep]
  if (length(unique(is_group)) < 2 || min(table(is_group)) < 3) return(NA_real_)
  out <- try(wilcox.test(value ~ is_group)$p.value, silent = TRUE)
  if (inherits(out, "try-error")) NA_real_ else out
}

panel <- read_csv(panel_path, show_col_types = FALSE) %>%
  mutate(
    gene = as.character(gene),
    layer_short = recode(Layer, "Transcriptome" = "R", "Methylation" = "M", .default = Layer),
    gene_label = paste0(gene, " [", layer_short, "]")
  ) %>%
  distinct(gene, .keep_all = TRUE)
genes <- panel$gene

clinical <- read_csv(clinical_path, show_col_types = FALSE) %>%
  filter(tpye == "Tumor", !is.na(Pseudotime_type), Pseudotime_type %in% group_order) %>%
  mutate(
    Sample = as.character(Row.names),
    group_raw = factor(Pseudotime_type, levels = group_order),
    group_label = factor(group_labels[as.character(group_raw)], levels = group_labels[group_order]),
    Cluster = factor(Cluster),
    PAM50 = coalesce(PAM50, "Unknown"),
    Pseudotime = as.numeric(Pseudotime)
  )

message("Reading RNA and methylation biomarker matrices...")
rna_long <- read_omics_long(rna_path, genes, "RNA", "log2")
meth_long <- read_omics_long(meth_path, genes, "METH", "none")

biomarker_long <- bind_rows(rna_long, meth_long) %>%
  inner_join(
    clinical %>% select(Sample, group_raw, group_label, Cluster, PAM50, Pseudotime, age, Stage_final),
    by = "Sample"
  ) %>%
  left_join(panel %>% select(gene, gene_label, Layer, Family6, target_label), by = "gene") %>%
  zscore_by_gene()

write_csv(biomarker_long, file.path(table_dir, "group_biomarker_long_values.csv"))

group_stats <- biomarker_long %>%
  group_by(modality, modality_label, gene, gene_label, Layer, Family6, target_label, group_label) %>%
  summarise(
    n = n(),
    mean_z = mean(value_z, na.rm = TRUE),
    median_z = median(value_z, na.rm = TRUE),
    sd_z = sd(value_z, na.rm = TRUE),
    .groups = "drop"
  )

kw_stats <- biomarker_long %>%
  group_by(modality, modality_label, gene, gene_label, Layer, Family6, target_label) %>%
  summarise(
    kw_p = safe_kw(value_z, group_label),
    group_mean_range = diff(range(tapply(value_z, group_label, mean, na.rm = TRUE), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  group_by(modality) %>%
  mutate(kw_fdr = p.adjust(kw_p, method = "BH")) %>%
  ungroup()

dominant_stats <- group_stats %>%
  group_by(modality, modality_label, gene, gene_label, Layer, Family6, target_label) %>%
  arrange(desc(mean_z), .by_group = TRUE) %>%
  summarise(
    dominant_group = as.character(first(group_label)),
    dominant_mean_z = first(mean_z),
    suppressed_group = as.character(last(group_label)),
    suppressed_mean_z = last(mean_z),
    contrast_z = dominant_mean_z - suppressed_mean_z,
    .groups = "drop"
  ) %>%
  left_join(kw_stats, by = c("modality", "modality_label", "gene", "gene_label", "Layer", "Family6", "target_label")) %>%
  rowwise() %>%
  mutate(
    dominant_vs_rest_p = safe_wilcox(
      biomarker_long$value_z[biomarker_long$modality == modality & biomarker_long$gene == gene],
      biomarker_long$group_label[biomarker_long$modality == modality & biomarker_long$gene == gene] == dominant_group
    )
  ) %>%
  ungroup() %>%
  group_by(modality) %>%
  mutate(dominant_vs_rest_fdr = p.adjust(dominant_vs_rest_p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    dominant_group = factor(dominant_group, levels = group_labels[group_order]),
    significance = case_when(
      kw_fdr < 0.05 ~ "FDR < 0.05",
      kw_fdr < 0.10 ~ "FDR < 0.10",
      kw_p < 0.05 ~ "P < 0.05",
      TRUE ~ "Exploratory"
    )
  )

write_csv(group_stats, file.path(table_dir, "group_biomarker_group_means.csv"))
write_csv(kw_stats, file.path(table_dir, "group_biomarker_kruskal_tests.csv"))
write_csv(dominant_stats, file.path(table_dir, "group_biomarker_dominant_group_summary.csv"))

gene_order <- dominant_stats %>%
  group_by(gene_label, Family6) %>%
  summarise(max_contrast = max(contrast_z, na.rm = TRUE),
            min_fdr = min(kw_fdr, na.rm = TRUE), .groups = "drop") %>%
  arrange(Family6, min_fdr, desc(max_contrast), gene_label) %>%
  pull(gene_label)

sample_order_tbl <- biomarker_long %>%
  distinct(modality, modality_label, Sample, group_label, Cluster, PAM50, Pseudotime) %>%
  group_by(modality, modality_label, group_label) %>%
  arrange(Pseudotime, Cluster, Sample, .by_group = TRUE) %>%
  mutate(sample_index_group = row_number()) %>%
  ungroup() %>%
  group_by(modality, modality_label) %>%
  arrange(group_label, sample_index_group, .by_group = TRUE) %>%
  mutate(sample_index = row_number()) %>%
  ungroup()

heatmap_df <- biomarker_long %>%
  inner_join(sample_order_tbl %>% select(modality, Sample, sample_index), by = c("modality", "Sample")) %>%
  mutate(
    gene_label = factor(gene_label, levels = rev(gene_order)),
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature"))
  )

sample_group_blocks <- sample_order_tbl %>%
  group_by(modality, modality_label, group_label) %>%
  summarise(xmin = min(sample_index) - 0.5, xmax = max(sample_index) + 0.5,
            xmid = mean(range(sample_index)), n = n(), .groups = "drop") %>%
  mutate(modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")))

annotation_df <- sample_order_tbl %>%
  select(modality, modality_label, sample_index, group_label, PAM50) %>%
  pivot_longer(c(group_label, PAM50), names_to = "annotation", values_to = "annotation_value") %>%
  mutate(
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")),
    annotation = factor(annotation, levels = c("group_label", "PAM50"), labels = c("EOBC group", "PAM50")),
    y = as.numeric(annotation)
  )

ann_cols <- c(group_cols, pam50_cols)

fig_sample_heatmap <- ggplot() +
  geom_rect(
    data = sample_group_blocks,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = group_label),
    alpha = 0.045,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = group_cols, guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = heatmap_df,
    aes(x = sample_index, y = gene_label, fill = value_z_clip),
    width = 1.0,
    height = 0.96
  ) +
  geom_vline(
    data = sample_group_blocks,
    aes(xintercept = xmax),
    colour = "white",
    linewidth = 0.42
  ) +
  geom_text(
    data = sample_group_blocks,
    aes(x = xmid, y = Inf, label = paste0(group_label, "\nn=", n)),
    inherit.aes = FALSE,
    vjust = 1.08,
    size = 2.55,
    fontface = "bold",
    colour = "white",
    lineheight = 0.86
  ) +
  facet_wrap(~ modality_label, ncol = 1, scales = "free_x") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = squish,
    name = "Within-gene\nz-score"
  ) +
  labs(
    title = "A. EOBC group-defining biomarker heatmap across patient samples",
    subtitle = "Samples are ordered by EOBC group number and within-group pseudotime; rows are the 26 candidate biomarkers grouped by biological family.",
    x = "Patient samples ordered within EOBC group",
    y = NULL,
    caption = "RNA values are log2(TPM+1); methylation values are beta values. Both are z-scored within each gene and omics layer."
  ) +
  theme_eobc(base_size = 8.4) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6.2),
    legend.position = "right",
    strip.text = element_text(size = 10.5)
  )

ggsave(file.path(plot_dir, "Figure_01_EOBC_group_biomarker_sample_heatmap_R_v1.png"),
       fig_sample_heatmap, width = 12.8, height = 9.3, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_01_EOBC_group_biomarker_sample_heatmap_R_v1.pdf"),
       fig_sample_heatmap, width = 12.8, height = 9.3, bg = "white", device = cairo_pdf)

mean_heatmap_df <- group_stats %>%
  left_join(kw_stats %>% select(modality, gene, kw_fdr), by = c("modality", "gene")) %>%
  mutate(
    gene_label = factor(gene_label, levels = rev(gene_order)),
    group_label = factor(group_label, levels = group_labels[group_order]),
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature"))
  )

dominant_marker_df <- dominant_stats %>%
  mutate(
    gene_label = factor(gene_label, levels = rev(gene_order)),
    group_label = factor(dominant_group, levels = group_labels[group_order]),
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")),
    mark = case_when(
      kw_fdr < 0.05 ~ "**",
      kw_fdr < 0.10 ~ "*",
      kw_p < 0.05 ~ ".",
      TRUE ~ ""
    )
  )

fig_mean_heatmap <- ggplot(mean_heatmap_df, aes(x = group_label, y = gene_label)) +
  geom_tile(aes(fill = mean_z), colour = "white", linewidth = 0.46) +
  geom_point(
    data = dominant_marker_df,
    aes(x = group_label, y = gene_label, size = contrast_z, colour = dominant_group),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    stroke = 0.85
  ) +
  geom_text(
    data = dominant_marker_df %>% filter(mark != ""),
    aes(x = group_label, y = gene_label, label = mark),
    inherit.aes = FALSE,
    colour = ink,
    fontface = "bold",
    size = 3.3,
    vjust = 0.45
  ) +
  facet_wrap(~ modality_label, ncol = 2) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1.4, 1.4),
    oob = squish,
    name = "Group mean\nz-score"
  ) +
  scale_colour_manual(values = group_cols, name = "Dominant group") +
  scale_size_continuous(range = c(1.8, 5.0), name = "Max-min\ngroup contrast") +
  labs(
    title = "B. EOBC group-level biomarker program map",
    subtitle = "Each tile is the group mean z-score. Open circles mark the group where each biomarker is highest; ** FDR<0.05, * FDR<0.10, . nominal P<0.05 by Kruskal-Wallis.",
    x = NULL,
    y = NULL
  ) +
  theme_eobc(base_size = 8.6) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 6.4),
    strip.text = element_text(size = 10.0),
    legend.box = "horizontal"
  )

ggsave(file.path(plot_dir, "Figure_02_EOBC_group_biomarker_mean_heatmap_R_v1.png"),
       fig_mean_heatmap, width = 12.8, height = 8.2, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02_EOBC_group_biomarker_mean_heatmap_R_v1.pdf"),
       fig_mean_heatmap, width = 12.8, height = 8.2, bg = "white", device = cairo_pdf)

dominance_plot_df <- dominant_stats %>%
  mutate(
    modality_label = factor(modality_label, levels = c("RNA feature", "Methylation feature")),
    neg_log10_fdr_raw = -log10(pmax(kw_fdr, 1e-12)),
    neg_log10_fdr = pmin(neg_log10_fdr_raw, 12)
  ) %>%
  group_by(modality_label) %>%
  mutate(
    rank_for_label = dense_rank(desc(contrast_z + neg_log10_fdr_raw / 10)),
    label = if_else(rank_for_label <= 18 | kw_fdr < 0.001, gene_label, NA_character_)
  ) %>%
  ungroup()

fig_dominance <- ggplot(dominance_plot_df, aes(x = contrast_z, y = neg_log10_fdr)) +
  geom_hline(yintercept = -log10(0.10), colour = "#FDBA74", linetype = "dotted", linewidth = 0.55) +
  geom_hline(yintercept = -log10(0.05), colour = "#F87171", linetype = "longdash", linewidth = 0.55) +
  geom_point(
    aes(fill = Family6, shape = dominant_group, size = abs(dominant_mean_z)),
    colour = "#111827",
    stroke = 0.55,
    alpha = 0.92
  ) +
  geom_text_repel(
    aes(label = label, colour = Family6),
    size = 2.15,
    fontface = "bold",
    box.padding = 0.16,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    segment.size = 0.16,
    force = 1.7,
    max.time = 3,
    max.overlaps = Inf,
    seed = 20260515,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  facet_wrap(~ modality_label, ncol = 2) +
  scale_fill_manual(values = family_cols, drop = FALSE, name = "Biomarker family") +
  scale_colour_manual(values = family_cols, guide = "none") +
  scale_shape_manual(
    values = c("G1 | H" = 21, "G2 | I" = 22, "G3 | L-like" = 24, "G4 | L" = 23),
    name = "Highest group"
  ) +
  scale_size_continuous(range = c(2.3, 7.0), name = "|Dominant\nmean z|") +
  scale_y_continuous(
    limits = c(0, 12.7),
    breaks = c(0, 2.5, 5, 7.5, 10, 12),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "C. Biomarker dominance landscape across EOBC groups",
    subtitle = "Genes far to the right have larger separation between their highest and lowest EOBC groups; y-axis is capped at -log10(FDR)=12 for readability.",
    x = "Dominant-vs-suppressed group contrast (max mean z - min mean z)",
    y = "-log10(Kruskal-Wallis FDR)"
  ) +
  theme_eobc(base_size = 8.8) +
  theme(
    legend.box = "horizontal",
    plot.margin = margin(12, 18, 10, 12)
  )

ggsave(file.path(plot_dir, "Figure_03_EOBC_group_biomarker_dominance_landscape_R_v1.png"),
       fig_dominance, width = 12.4, height = 5.7, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03_EOBC_group_biomarker_dominance_landscape_R_v1.pdf"),
       fig_dominance, width = 12.4, height = 5.7, bg = "white", device = cairo_pdf)

top_summary <- dominant_stats %>%
  arrange(modality, kw_fdr, desc(contrast_z)) %>%
  group_by(modality_label, dominant_group) %>%
  slice_head(n = 5) %>%
  ungroup()
write_csv(top_summary, file.path(table_dir, "group_biomarker_top_markers_by_group.csv"))

analysis_note <- c(
  "EOBC group-centric biomarker analysis",
  "",
  "Grouping:",
  "- Pseudotime_type was mapped to four EOBC states: G1|H, G2|I, G4|L-like, G3|L.",
  "- Samples without Pseudotime_type were excluded from this group-centric analysis.",
  "",
  "Features:",
  "- RNA: log2(TPM+1), z-scored within each biomarker gene.",
  "- Methylation: beta value, z-scored within each biomarker gene.",
  "- All 26 selected candidate biomarkers were tested in both omics layers when available.",
  "",
  "Statistics:",
  "- Kruskal-Wallis test across the four EOBC groups per biomarker and omics layer.",
  "- FDR was adjusted within each omics layer.",
  "- Dominant group is the EOBC group with the highest mean z-score for that biomarker.",
  "",
  "Interpretation:",
  "- These figures should precede OS/drug/immune consequence figures because they establish which biomarkers actually define the EOBC group states."
)
writeLines(analysis_note, file.path(out_dir, "analysis_note.txt"))

message("Done.")
message("Plots written to: ", plot_dir)
message("Tables written to: ", table_dir)
