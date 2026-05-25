#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(survival)
  library(survminer)
})

options(stringsAsFactors = FALSE)
set.seed(20260522)

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
analysis_dir <- file.path(biomarker_root, "final_analysis")
final_fig_dir <- file.path(analysis_dir, "final_fig")
final_revision_dir <- file.path(final_fig_dir, "final_paper_figures_20260522_v2")

out_roots <- unique(c(
  file.path(final_fig_dir, "individual_biomarker_evidence_20260522"),
  file.path(final_revision_dir, "individual_biomarker_evidence_20260522")
))

subdirs <- c("OS_KM", "DepMap_drug_response", "Immune_TIL_TMB", "Immune_TIL", "Immune_TMB", "tables")
invisible(lapply(out_roots, function(root) {
  invisible(lapply(file.path(root, subdirs), dir.create, recursive = TRUE, showWarnings = FALSE))
}))

theme_paper <- function(base_size = 9.0) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", color = "#111827", size = rel(1.12)),
      plot.subtitle = element_text(color = "#4B5563", size = rel(0.86), margin = margin(b = 7)),
      plot.caption = element_text(color = "#64748B", size = rel(0.72), hjust = 0),
      axis.title = element_text(face = "bold", color = "#111827"),
      axis.text = element_text(color = "#475569"),
      legend.title = element_text(face = "bold", color = "#111827"),
      legend.text = element_text(color = "#475569"),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 12, 10, 12)
    )
}

fmt_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

fmt_num <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "")
}

save_plot_all <- function(plot, subdir, stub, width, height, dpi = 450) {
  for (root in out_roots) {
    dest_dir <- file.path(root, subdir)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(dest_dir, paste0(stub, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
    ggsave(file.path(dest_dir, paste0(stub, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
  }
}

save_surv_all <- function(g, subdir, stub, width = 6.2, height = 5.8, res = 450) {
  for (root in out_roots) {
    dest_dir <- file.path(root, subdir)
    dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
    png(file.path(dest_dir, paste0(stub, ".png")), width = width, height = height, units = "in", res = res, bg = "white")
    print(g)
    dev.off()
    cairo_pdf(file.path(dest_dir, paste0(stub, ".pdf")), width = width, height = height, bg = "white")
    print(g)
    dev.off()
  }
}

write_table_all <- function(x, file_name) {
  for (root in out_roots) {
    write_csv(x, file.path(root, "tables", file_name))
  }
}

family_cols <- c(
  "Immune" = "#4EA5F0",
  "Repair" = "#47C56B",
  "Glycolysis / TCA" = "#E6BC18",
  "Fatty acid" = "#F4A259",
  "Kinase signaling" = "#9B6AE8",
  "Hormone signaling" = "#9AAABC"
)

group_cols <- c(
  "G1 | H" = "#D11224",
  "G2 | I" = "#FF9416",
  "G3 | L-like" = "#8BBE69",
  "G4 | L" = "#1B9278"
)

# ---------------------------------------------------------------------------
# OS: individual Kaplan-Meier curves with risk tables.
# ---------------------------------------------------------------------------
km_path <- file.path(analysis_dir, "domain_publication_suite_v2", "tables", "os_single_gene_km_curve_table.csv")
os_sum_path <- file.path(analysis_dir, "domain_publication_suite_v2", "tables", "os_single_gene_survival_summary.csv")

km_df <- read_csv(km_path, show_col_types = FALSE) %>%
  mutate(
    time_months = as.numeric(time_months),
    event = as.integer(event),
    group = factor(group, levels = c("Low", "High"))
  ) %>%
  filter(is.finite(time_months), !is.na(event), !is.na(group))

os_summary <- read_csv(os_sum_path, show_col_types = FALSE) %>%
  mutate(across(c(logrank_p, hr, lower95, upper95, p_value), as.numeric))

endpoint_label <- c(
  "overall_os" = "Overall OS",
  "os_5y" = "5-year OS",
  "os_10y" = "10-year OS"
)

os_pairs <- km_df %>%
  distinct(gene, modality, endpoint) %>%
  arrange(gene, modality, endpoint)

os_manifest <- list()
for (i in seq_len(nrow(os_pairs))) {
  gene_i <- os_pairs$gene[i]
  modality_i <- os_pairs$modality[i]
  endpoint_i <- os_pairs$endpoint[i]
  dat <- km_df %>%
    filter(gene == gene_i, modality == modality_i, endpoint == endpoint_i)
  if (nrow(dat) < 20 || n_distinct(dat$group) < 2) next

  fit <- survfit(Surv(time_months, event) ~ group, data = dat)
  info <- os_summary %>%
    filter(gene == gene_i, modality == modality_i, endpoint == endpoint_i) %>%
    slice(1)
  ep_lab <- endpoint_label[[endpoint_i]]
  if (is.null(ep_lab)) ep_lab <- endpoint_i
  p_label <- paste0(
    "log-rank p = ", fmt_p(info$logrank_p), "\n",
    "HR = ", fmt_num(info$hr), " (95% CI ", fmt_num(info$lower95), "-", fmt_num(info$upper95), ")"
  )
  max_time <- max(dat$time_months, na.rm = TRUE)
  break_by <- ifelse(max_time <= 66, 12, 24)

  g <- ggsurvplot(
    fit,
    data = dat,
    risk.table = TRUE,
    risk.table.height = 0.27,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = TRUE,
    conf.int = FALSE,
    censor = TRUE,
    censor.shape = 124,
    censor.size = 2.5,
    palette = c("Low" = "#2C7FB8", "High" = "#D95F02"),
    break.time.by = break_by,
    pval = p_label,
    pval.coord = c(max_time * 0.05, 0.13),
    pval.size = 3.2,
    risk.table.fontsize = 3.0,
    legend.title = paste0(gene_i, " ", modality_i),
    legend.labs = c("Low", "High"),
    xlab = "Months",
    ylab = "Overall survival probability",
    title = paste0(gene_i, " | ", ep_lab),
    ggtheme = theme_paper(8.2)
  )
  g$plot <- g$plot +
    coord_cartesian(ylim = c(0, 1.02), clip = "off") +
    theme(
      legend.position = c(0.80, 0.18),
      legend.background = element_rect(fill = scales::alpha("white", 0.82), color = "#E5E7EB"),
      plot.title = element_text(size = 10.5, face = "bold"),
      axis.title = element_text(size = 8.6, face = "bold"),
      axis.text = element_text(size = 7.2),
      legend.title = element_text(size = 7.2, face = "bold"),
      legend.text = element_text(size = 7.0)
    )
  g$table <- g$table +
    theme_paper(7.3) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 6.6),
      axis.text.y = element_text(size = 6.6)
    )

  stub <- paste0("Figure_OS_KM_", safe_name(gene_i), "_", safe_name(endpoint_i), "_", safe_name(modality_i))
  save_surv_all(g, "OS_KM", stub, width = 5.9, height = 4.8)
  os_manifest[[length(os_manifest) + 1]] <- tibble(gene = gene_i, modality = modality_i, endpoint = endpoint_i, file_stub = stub)
}
os_manifest <- bind_rows(os_manifest)
write_table_all(os_manifest, "individual_OS_KM_manifest.csv")

# ---------------------------------------------------------------------------
# DepMap: per-gene dose-response curves for prioritized biomarker-drug pairs.
# ---------------------------------------------------------------------------
drug_path <- file.path(biomarker_root, "15_final_manuscript_OS_TILTMB_drug_sankey", "tables", "gene_drug_topk_hits.csv")
drug_df <- read_csv(drug_path, show_col_types = FALSE) %>%
  mutate(
    direction_score = as.numeric(direction_score),
    q_primary = as.numeric(q_primary),
    rank_in_gene = as.numeric(rank_in_gene),
    drug_target_class = na_if(drug_target_class, ""),
    direction_label = if_else(direction_score < 0, "Sensitivity-associated", "Resistance-associated"),
    drug_label = paste0(drug_short, "\n", drug_target_class),
    abs_score = abs(direction_score)
  ) %>%
  filter(
    is.finite(direction_score),
    is.finite(q_primary),
    q_primary <= 0.25,
    !is.na(drug_target_class),
    !str_detect(drug_target_class, regex("Other|Experimental", ignore_case = TRUE))
  ) %>%
  group_by(gene) %>%
  arrange(q_primary, desc(abs_score), .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup()

write_table_all(drug_df, "individual_DepMap_curated_drug_hits_used.csv")

dose_path <- file.path(biomarker_root, "13_final_manuscript_validation_suite_v5_full_rewrite_raw_depmap_SVS_MAF", "tables", "top_gene_dose_response_curve_table.csv")
dose_df <- read_csv(dose_path, show_col_types = FALSE) %>%
  mutate(
    dose_uM = as.numeric(dose_uM),
    log10_uM = as.numeric(log10_uM),
    mean_viability = as.numeric(mean_viability),
    se_viability = as.numeric(se_viability),
    n = as.integer(n),
    group = factor(str_to_lower(group), levels = c("low", "high"), labels = c("Low marker", "High marker")),
    drug_target_class = na_if(drug_target_class, "")
  ) %>%
  filter(
    gene %in% unique(drug_df$gene),
    is.finite(log10_uM),
    is.finite(mean_viability),
    !is.na(group),
    !is.na(drug_target_class),
    !str_detect(drug_target_class, regex("Other|Experimental", ignore_case = TRUE))
  )

thin_dose_curve <- function(df, max_points = 8) {
  df %>%
    group_by(gene, Family6, drug_clean, drug_target_class, direction, group) %>%
    arrange(log10_uM, .by_group = TRUE) %>%
    group_modify(~ {
      x <- .x %>% arrange(log10_uM)
      if (n_distinct(x$log10_uM) <= max_points) return(x)
      x %>%
        mutate(.dose_bin = ntile(row_number(), max_points)) %>%
        group_by(.dose_bin) %>%
        summarise(
          dose_uM = 10^weighted.mean(log10_uM, pmax(n, 1), na.rm = TRUE),
          log10_uM = weighted.mean(log10_uM, pmax(n, 1), na.rm = TRUE),
          mean_viability = weighted.mean(mean_viability, pmax(n, 1), na.rm = TRUE),
          se_viability = sqrt(sum((coalesce(se_viability, 0)^2) * pmax(n, 1), na.rm = TRUE)) /
            sqrt(sum(pmax(n, 1), na.rm = TRUE)),
          n = sum(n, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        select(-.dose_bin)
    }) %>%
    ungroup()
}

curve_pairs <- dose_df %>%
  distinct(gene, Family6, drug_clean, drug_target_class, direction) %>%
  left_join(
    drug_df %>%
      select(gene, drug_clean, direction_score, q_primary, rank_in_gene, direction_label) %>%
      distinct(),
    by = c("gene", "drug_clean")
  ) %>%
  mutate(
    q_sort = if_else(is.finite(q_primary), q_primary, 1),
    rank_sort = if_else(is.finite(rank_in_gene), rank_in_gene, 999),
    direction_label = coalesce(
      direction_label,
      case_when(
        direction == "more_sensitive_when_high" ~ "High marker sensitive",
        direction == "more_resistant_when_high" ~ "High marker resistant",
        TRUE ~ "Direction not assigned"
      )
    )
  ) %>%
  group_by(gene) %>%
  arrange(q_sort, rank_sort, drug_clean, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup()

write_table_all(curve_pairs, "individual_DepMap_dose_response_curve_pairs_used.csv")

for (root in out_roots) {
  depmap_dir <- file.path(root, "DepMap_drug_response")
  old_depmap_files <- list.files(
    depmap_dir,
    pattern = "^Figure_DepMap_drug_response_.*\\.(png|pdf)$",
    full.names = TRUE
  )
  if (length(old_depmap_files) > 0) unlink(old_depmap_files)
}

depmap_manifest <- list()
for (gene_i in sort(unique(curve_pairs$gene))) {
  pair_i <- curve_pairs %>% filter(gene == gene_i) %>% slice_head(n = 1)
  dat <- dose_df %>%
    filter(gene == gene_i, drug_clean == pair_i$drug_clean) %>%
    arrange(group, log10_uM)
  if (nrow(dat) == 0) next
  dat_plot <- thin_dose_curve(dat, max_points = 8)

  direction_tag <- case_when(
    pair_i$direction == "more_sensitive_when_high" ~ "High marker group shows greater drug sensitivity",
    pair_i$direction == "more_resistant_when_high" ~ "High marker group shows greater drug resistance",
    TRUE ~ "High/low marker response curves"
  )
  stat_label <- if (is.finite(pair_i$direction_score)) {
    paste0("rho(AUC)=", fmt_num(pair_i$direction_score, 2), "; q=", fmt_p(pair_i$q_primary))
  } else {
    "AUC association from selected DepMap pair"
  }

  p <- ggplot(dat_plot, aes(log10_uM, mean_viability, color = group, fill = group)) +
    geom_ribbon(
      aes(ymin = pmax(0, mean_viability - coalesce(se_viability, 0)), ymax = pmin(1.05, mean_viability + coalesce(se_viability, 0))),
      alpha = 0.14,
      color = NA
    ) +
    geom_line(linewidth = 1.15) +
    geom_point(size = 2.05, shape = 21, stroke = 0.55, color = "white") +
    scale_color_manual(values = c("Low marker" = "#2563EB", "High marker" = "#D95F02"), name = NULL) +
    scale_fill_manual(values = c("Low marker" = "#2563EB", "High marker" = "#D95F02"), name = NULL) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), labels = number_format(accuracy = 0.01)) +
    scale_x_continuous(breaks = pretty_breaks(n = 6)) +
    labs(
      title = paste0(gene_i, " dose-response curve"),
      subtitle = paste0(pair_i$drug_clean, " | ", pair_i$drug_target_class, "   ", stat_label),
      x = "Drug dose (log10 uM)",
      y = "Mean viability",
      caption = paste0(direction_tag, ". Curves are grouped by marker-high versus marker-low cell-line strata used for AUC-based DepMap prioritization.")
    ) +
    theme_paper(8.6) +
    theme(
      panel.grid.major = element_line(color = "#E2E8F0", linewidth = 0.35),
      panel.grid.minor = element_line(color = "#EEF2F7", linewidth = 0.25),
      legend.position = "bottom",
      legend.justification = "left"
    )

  stub <- paste0("Figure_DepMap_drug_response_", safe_name(gene_i))
  save_plot_all(p, "DepMap_drug_response", stub, width = 5.6, height = 4.2)
  depmap_manifest[[length(depmap_manifest) + 1]] <- tibble(
    gene = gene_i,
    drug_clean = pair_i$drug_clean,
    drug_target_class = pair_i$drug_target_class,
    direction = pair_i$direction,
    rho_auc = pair_i$direction_score,
    q_auc = pair_i$q_primary,
    n_raw_curve_points = nrow(dat),
    n_curve_points = nrow(dat_plot),
    file_stub = stub
  )
}
depmap_manifest <- bind_rows(depmap_manifest)
write_table_all(depmap_manifest, "individual_DepMap_drug_response_manifest.csv")

# ---------------------------------------------------------------------------
# Immune: per-gene RNA expression versus TIL score and TMB score.
# ---------------------------------------------------------------------------
expr_path <- file.path(analysis_dir, "group_biomarker_landscape_r_v1", "tables", "group_biomarker_long_values.csv")
immune_path <- file.path(biomarker_root, "10_external_validation_TIL_TMB_DepMap_final_v7", "tables", "immune_signature_tcga_til_tmb_merged.csv")

expr_df <- read_csv(expr_path, show_col_types = FALSE) %>%
  filter(modality == "RNA") %>%
  transmute(
    Sample,
    gene,
    Family6,
    group_label,
    rna_z = as.numeric(value_z)
  )

immune_metrics <- read_csv(immune_path, show_col_types = FALSE) %>%
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
write_table_all(tmb_outlier_samples, "individual_Immune_TMB_outlier_samples_excluded.csv")

immune_plot_df <- bind_rows(
  expr_df %>%
    inner_join(immune_metrics %>% select(Sample, TIL_score), by = "Sample") %>%
    transmute(
      Sample, gene, Family6, group_label, rna_z,
      metric = "TIL score",
      metric_value = TIL_score
    ),
  expr_df %>%
    inner_join(immune_metrics %>% select(Sample, TMB_log1p), by = "Sample") %>%
    filter(!Sample %in% tmb_outlier_samples$Sample) %>%
    transmute(
      Sample, gene, Family6, group_label, rna_z,
      metric = "TMB log1p",
      metric_value = TMB_log1p
    )
) %>%
  filter(is.finite(rna_z), is.finite(metric_value)) %>%
  mutate(metric = factor(metric, levels = c("TIL score", "TMB log1p")))

cor_summary <- immune_plot_df %>%
  group_by(gene, Family6, metric) %>%
  summarise(
    n = n(),
    rho = suppressWarnings(cor(rna_z, metric_value, method = "spearman")),
    p = suppressWarnings(cor.test(rna_z, metric_value, method = "spearman", exact = FALSE)$p.value),
    .groups = "drop"
  ) %>%
  group_by(metric) %>%
  mutate(q = p.adjust(p, method = "BH")) %>%
  ungroup()

write_table_all(cor_summary, "individual_Immune_expression_TIL_TMB_correlation_summary.csv")

make_immune_panel <- function(dat, stats, gene_i, metric_i, show_legend = TRUE) {
  fam_i <- dat$Family6[1]
  fam_col <- family_cols[[fam_i]]
  if (is.null(fam_col) || is.na(fam_col)) fam_col <- "#64748B"
  stat <- stats %>% filter(gene == gene_i, metric == metric_i) %>% slice(1)
  y_lab <- if (metric_i == "TIL score") "TIL score" else "TMB score (log1p)"
  ggplot(dat, aes(rna_z, metric_value)) +
    geom_point(aes(fill = group_label), shape = 21, size = 2.9, color = "white", stroke = 0.45, alpha = 0.92) +
    geom_smooth(method = "lm", se = TRUE, color = fam_col, fill = scales::alpha(fam_col, 0.16), linewidth = 0.85) +
    scale_fill_manual(values = group_cols, name = "EOBC group", drop = FALSE) +
    labs(
      title = paste0(gene_i, " vs ", y_lab),
      subtitle = paste0("Spearman rho = ", fmt_num(stat$rho), ", p = ", fmt_p(stat$p), ", q = ", fmt_p(stat$q), " (n=", stat$n, ")"),
      x = "RNA expression z-score",
      y = y_lab
    ) +
    theme_paper(8.5) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      panel.grid.major = element_line(color = "#E2E8F0", linewidth = 0.32),
      panel.grid.minor = element_line(color = "#F1F5F9", linewidth = 0.20),
      plot.title = element_text(size = 9.4, face = "bold"),
      plot.subtitle = element_text(size = 7.3)
    )
}

immune_manifest <- list()
for (gene_i in sort(unique(immune_plot_df$gene))) {
  dat_gene <- immune_plot_df %>% filter(gene == gene_i)
  if (nrow(dat_gene) == 0) next

  p_til <- make_immune_panel(dat_gene %>% filter(metric == "TIL score"), cor_summary, gene_i, "TIL score", show_legend = FALSE)
  p_tmb <- make_immune_panel(dat_gene %>% filter(metric == "TMB log1p"), cor_summary, gene_i, "TMB log1p", show_legend = FALSE)
  p_combo <- (p_til | p_tmb) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  stub_combo <- paste0("Figure_Immune_expression_TIL_TMB_", safe_name(gene_i))
  stub_til <- paste0("Figure_Immune_expression_TIL_", safe_name(gene_i))
  stub_tmb <- paste0("Figure_Immune_expression_TMB_", safe_name(gene_i))
  save_plot_all(p_combo, "Immune_TIL_TMB", stub_combo, width = 8.8, height = 3.8)
  save_plot_all(p_til + theme(legend.position = "bottom"), "Immune_TIL", stub_til, width = 4.7, height = 3.8)
  save_plot_all(p_tmb + theme(legend.position = "bottom"), "Immune_TMB", stub_tmb, width = 4.7, height = 3.8)

  immune_manifest[[length(immune_manifest) + 1]] <- tibble(
    gene = gene_i,
    combined_file_stub = stub_combo,
    til_file_stub = stub_til,
    tmb_file_stub = stub_tmb
  )
}
immune_manifest <- bind_rows(immune_manifest)
write_table_all(immune_manifest, "individual_Immune_TIL_TMB_manifest.csv")

manifest <- bind_rows(
  os_manifest %>% transmute(domain = "OS", gene, detail = endpoint, file_stub),
  depmap_manifest %>% transmute(domain = "DepMap", gene, detail = paste0(drug_clean, " dose-response curve"), file_stub),
  immune_manifest %>% transmute(domain = "Immune", gene, detail = "RNA expression vs TIL and TMB", file_stub = combined_file_stub)
)
write_table_all(manifest, "individual_biomarker_evidence_manifest.csv")

message("Saved individual biomarker evidence plots to:")
message("  ", paste(out_roots, collapse = "\n  "))
