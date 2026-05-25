library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(stringr)
library(scales)

set.seed(20260513)

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

biomarker_root_env <- Sys.getenv("EOBC_BIOMARKER_ROOT", unset = "")
root <- env_path(
  "EOBC_PROJECT_ROOT",
  if (nzchar(biomarker_root_env)) dirname(normalizePath(biomarker_root_env, winslash = "/", mustWork = FALSE)) else NULL
)
biomarker_root <- file.path(root, "biomarker")
input_dir <- file.path(biomarker_root, "00_inputs_detected")
manuscript_dir <- file.path(
  biomarker_root,
  "13_final_manuscript_validation_suite_v5_full_rewrite_raw_depmap_SVS_MAF"
)
out_dir <- file.path(
  biomarker_root,
  "final_analysis",
  "full_layer_candidate_validation_r_v1"
)
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

candidate_path <- file.path(manuscript_dir, "tables", "selected_biomarker_panel_used.csv")
depmap_rna_path <- file.path(root, "depmap", "RNA", "dep_RNA.csv")
depmap_meth_path <- file.path(
  root, "depmap", "meth", "Methylation_(1kb_upstream_TSS)_subsetted_NAsdropped.csv"
)
depmap_auc_path <- file.path(
  root, "depmap", "drug", "Drug_sensitivity_AUC_(Sanger_GDSC2)_subsetted_NAsdropped.csv"
)
depmap_ic50_path <- file.path(
  root, "depmap", "drug", "Drug_sensitivity_IC50_(Sanger_GDSC2)_subsetted_NAsdropped.csv"
)
old_depmap_annot_path <- file.path(
  manuscript_dir, "tables", "depmap_gene_drug_association_joint_AUC_IC50_annotated.csv"
)
legacy_immune_path <- file.path(
  biomarker_root,
  "final_analysis",
  "domain_publication_suite_v2",
  "tables",
  "immune_all_biomarker_landscape.csv"
)
tpm_path <- file.path(input_dir, "TPM_young.csv")
patient_meth_path <- file.path(input_dir, "MET_young_batch_JW.csv")
immune_ref_path <- file.path(
  biomarker_root,
  "10_external_validation_TIL_TMB_DepMap",
  "tables",
  "immune_signature_tcga_til_tmb_merged.csv"
)

family_cols <- c(
  "Immune" = "#4EA5F5",
  "Repair" = "#44C06A",
  "Glycolysis / TCA" = "#D8B425",
  "Fatty acid" = "#F0A35E",
  "Kinase signaling" = "#A57AE5",
  "Hormone signaling" = "#90A2B8"
)
layer_cols <- c("RNA" = "#4EA5F5", "METH" = "#F0A35E")
drug_shape_values <- c(
  "DDR / Checkpoint" = 21,
  "RTK/MAPK kinase" = 24,
  "Mitotic / Cell-cycle" = 22,
  "PI3K/AKT/mTOR" = 23,
  "Apoptosis / p53" = 25,
  "Epigenetic / chromatin" = 8,
  "Cytotoxic / chemotherapy" = 3,
  "Metabolic / WNT" = 4,
  "Bioactive / experimental" = 7
)
layer_shape_values <- c("Transcriptome" = 21, "Methylation" = 22)
ink <- "#111827"
muted <- "#64748B"
grid_col <- "#D9E3EF"
protective_bg <- "#EAF4FF"
adverse_bg <- "#FFF3EA"

theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.45),
      panel.grid.minor = element_line(colour = alpha(grid_col, 0.45), linewidth = 0.25),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = muted),
      plot.title = element_text(face = "bold", colour = ink, size = rel(1.18), margin = margin(b = 4)),
      plot.subtitle = element_text(colour = muted, size = rel(0.82), margin = margin(b = 10)),
      strip.text = element_text(face = "bold", colour = ink, size = rel(1.02)),
      legend.title = element_text(face = "bold", colour = ink),
      legend.text = element_text(colour = muted),
      legend.background = element_rect(fill = "white", colour = grid_col),
      legend.key = element_rect(fill = "white", colour = NA),
      plot.margin = margin(14, 18, 12, 18)
    )
}

clean_drug <- function(x) str_trim(str_replace(x, "\\s*\\(GDSC2:[0-9]+\\)$", ""))

classify_drug <- function(x) {
  u <- toupper(x)
  case_when(
    str_detect(u, "CERALASERTIB|ADAVOSERTIB|AZD7762|MK-8776|OLAPARIB|NIRAPARIB|TALAZOPARIB|BERZOSERTIB|KU-55933|KU-57788|VE 821|MIRIN|CRA-032765|PARP|ATR|ATM|WEE|CHK|DNA-PK") ~ "DDR / Checkpoint",
    str_detect(u, "GEFITINIB|ERLOTINIB|LAPATINIB|AFATINIB|OSIMERTINIB|AZD8931|TRAMETINIB|DABRAFENIB|DASATINIB|NILOTINIB|CRIZOTINIB|VX-11E|BVD-523|SCH772984|PD 0325901|BIRB 796|PLX 4720|RTK|MAPK|MEK|BRAF|EGFR|ABL|ALK") ~ "RTK/MAPK kinase",
    str_detect(u, "ALISERTIB|BI-2536|MK 0457|ZM 447439|PALBOCICLIB|DINACICLIB|RIBOCICLIB|RO-3306|284461-73-0|MITOTIC|AURORA|CDK|PLK") ~ "Mitotic / Cell-cycle",
    str_detect(u, "BEZ235|AZD-8055|RAPAMYCIN|PICTILISIB|ALPELISIB|TASELISIB|MK 2206|GSK-2141795|PF 4708671|PI3K|AKT|MTOR|S6K") ~ "PI3K/AKT/mTOR",
    str_detect(u, "NAVITOCLAX|GX15-070|LCL161|REBEMADLIN|APR-246|MDM2|BCL|SMAC|P53") ~ "Apoptosis / p53",
    str_detect(u, "VORINOSTAT|EPZ-004777|PINOMETOSTAT|HDAC|DOT1L|EZH|BET|JQ1") ~ "Epigenetic / chromatin",
    str_detect(u, "CAMPTOTHECIN|CIS-DDP|DOCETAXEL|PACLITAXEL|IRINOTECAN|TOPOTECAN|GEMCITABINE|FLUOROURACIL|CARMUSTINE|VINCRISTINE|VELBAN|MITOXANTRONE|EPIRUBICIN|CYCLOPHOSPHAMIDE|BORTEZOMIB|CYTARABINE|OXALIPLATIN|ELOXATIN|SINULARIN|ELEPHANTIN") ~ "Cytotoxic / chemotherapy",
    str_detect(u, "AGI-6780|SB 216763|XAV-939|APO866|IDH|GSK3|TANKYRASE|WNT|NAMPT") ~ "Metabolic / WNT",
    TRUE ~ "Bioactive / experimental"
  )
}

read_gene_matrix_long <- function(path, genes, omics, value_transform = identity) {
  dat <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  names(dat)[1] <- "gene"
  dat %>%
    filter(.data$gene %in% genes) %>%
    pivot_longer(-gene, names_to = "sample", values_to = "feature_value") %>%
    mutate(
      feature_value = value_transform(as.numeric(feature_value)),
      omics = omics
    )
}

read_depmap_drug_long <- function(path, metric) {
  dat <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  meta_cols <- c("depmap_id", "cell_line_display_name", "lineage_1", "lineage_2",
                 "lineage_3", "lineage_6", "lineage_4")
  dat %>%
    select(any_of(meta_cols), everything()) %>%
    pivot_longer(
      cols = -any_of(meta_cols),
      names_to = "drug_raw",
      values_to = "response_value"
    ) %>%
    mutate(
      response_value = as.numeric(response_value),
      metric = metric,
      drug_clean = clean_drug(drug_raw)
    )
}

safe_spearman <- function(x, y, min_n = 8) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(n = n, rho = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(n = n, rho = unname(test$estimate), p = test$p.value)
}

candidate <- read_csv(candidate_path, show_col_types = FALSE) %>%
  transmute(
    gene = as.character(gene),
    original_layer = Layer,
    Family6,
    target_label,
    gene_label = paste0(gene, " [", if_else(Layer == "Methylation", "M", "R"), "]")
  ) %>%
  distinct(gene, .keep_all = TRUE)
genes <- candidate$gene

message("Building DepMap RNA features for all candidate genes...")
dep_rna <- read_csv(depmap_rna_path, show_col_types = FALSE, progress = FALSE)
names(dep_rna)[1] <- "gene"
depmap_rna_long <- dep_rna %>%
  filter(gene %in% genes) %>%
  pivot_longer(-gene, names_to = "cell_line_display_name", values_to = "feature_value") %>%
  mutate(feature_value = as.numeric(feature_value), omics = "RNA") %>%
  select(omics, gene, cell_line_display_name, feature_value)

message("Building DepMap TSS-methylation gene-level features for all candidate genes...")
dep_meth <- read_csv(depmap_meth_path, show_col_types = FALSE, progress = FALSE)
dep_meth_features <- lapply(genes, function(g) {
  cols <- names(dep_meth)[str_detect(names(dep_meth), paste0("^", fixed(g), "_"))]
  if (length(cols) == 0) {
    return(tibble(
      omics = "METH", gene = g, cell_line_display_name = character(),
      feature_value = numeric(), n_depmap_meth_probes = integer()
    ))
  }
  tibble(
    omics = "METH",
    gene = g,
    cell_line_display_name = dep_meth$cell_line_display_name,
    feature_value = rowMeans(as.data.frame(dep_meth[, cols]), na.rm = TRUE),
    n_depmap_meth_probes = length(cols)
  )
}) %>% bind_rows()

depmap_features <- bind_rows(
  depmap_rna_long %>% mutate(n_depmap_meth_probes = NA_integer_),
  dep_meth_features
) %>%
  left_join(candidate, by = "gene")

message("Running all-candidate DepMap drug-response correlations by RNA and methylation layer...")
drug_long <- bind_rows(
  read_depmap_drug_long(depmap_auc_path, "AUC"),
  read_depmap_drug_long(depmap_ic50_path, "IC50")
) %>%
  select(metric, depmap_id, cell_line_display_name, drug_raw, drug_clean, response_value)

depmap_assoc <- depmap_features %>%
  inner_join(drug_long, by = "cell_line_display_name", relationship = "many-to-many") %>%
  group_by(omics, gene, Family6, original_layer, target_label, gene_label,
           metric, drug_raw, drug_clean) %>%
  group_modify(~ safe_spearman(.x$feature_value, .x$response_value, min_n = 8)) %>%
  ungroup() %>%
  group_by(omics, metric) %>%
  mutate(fdr = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    direction = case_when(
      is.na(rho) ~ "not_tested",
      rho < 0 ~ "High biomarker = sensitive",
      rho > 0 ~ "High biomarker = resistant",
      TRUE ~ "neutral"
    ),
    neglog10_fdr = -log10(pmax(fdr, 1e-300))
  )

old_drug_map <- if (file.exists(old_depmap_annot_path)) {
  read_csv(old_depmap_annot_path, show_col_types = FALSE) %>%
    select(drug_raw, drug_clean, drug_target_class) %>%
    distinct() %>%
    filter(!is.na(drug_raw), !is.na(drug_target_class)) %>%
    group_by(drug_raw, drug_clean) %>%
    summarise(drug_target_class = first(drug_target_class), .groups = "drop")
} else {
  tibble(drug_raw = character(), drug_clean = character(), drug_target_class = character())
}

candidate_layer_grid <- tidyr::crossing(
  gene = genes,
  omics = c("RNA", "METH")
) %>%
  left_join(candidate, by = "gene")

depmap_best <- depmap_assoc %>%
  filter(is.finite(fdr), is.finite(rho)) %>%
  arrange(omics, gene, fdr, p, desc(abs(rho)), desc(n)) %>%
  group_by(omics, gene) %>%
  slice(1) %>%
  ungroup() %>%
  left_join(old_drug_map, by = c("drug_raw", "drug_clean")) %>%
  mutate(
    original_drug_target_class = drug_target_class,
    curated_drug_class = classify_drug(drug_clean),
    drug_target_class = case_when(
      is.na(drug_clean) | drug_clean == "" ~ NA_character_,
      !is.na(curated_drug_class) & curated_drug_class != "" ~ curated_drug_class,
      !is.na(original_drug_target_class) & original_drug_target_class != "" ~ original_drug_target_class,
      TRUE ~ "Bioactive / experimental"
    )
  )

depmap_available <- depmap_features %>%
  group_by(omics, gene) %>%
  summarise(
    n_feature_values = sum(is.finite(feature_value)),
    n_depmap_meth_probes = suppressWarnings(max(n_depmap_meth_probes, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(n_depmap_meth_probes = if_else(is.infinite(n_depmap_meth_probes), NA_integer_, as.integer(n_depmap_meth_probes)))

depmap_best_full <- candidate_layer_grid %>%
  left_join(depmap_available, by = c("omics", "gene")) %>%
  left_join(
    depmap_best %>%
      select(omics, gene, metric, drug_raw, drug_clean, drug_target_class, n, rho, p, fdr,
             direction, neglog10_fdr),
    by = c("omics", "gene")
  ) %>%
  mutate(
    feature_available = !is.na(n_feature_values) & n_feature_values >= 8,
    plot_y = if_else(is.finite(neglog10_fdr), neglog10_fdr, NA_real_),
    effect_abs = if_else(is.finite(rho), abs(rho), NA_real_),
    layer_label = factor(
      recode(omics, RNA = "RNA feature", METH = "Methylation feature"),
      levels = c("RNA feature", "Methylation feature")
    )
  )

write_csv(depmap_assoc, file.path(table_dir, "depmap_all_candidate_all_drug_correlations_by_omics.csv"))
write_csv(depmap_best_full, file.path(table_dir, "depmap_all_candidate_best_hits_by_omics.csv"))

message("Running all-candidate immune TIL/TMB correlations by patient RNA and methylation layer...")
patient_rna_long <- read_gene_matrix_long(tpm_path, genes, "RNA", value_transform = function(x) log2(x + 1)) %>%
  rename(cell_or_sample = sample)
patient_meth_long <- read_gene_matrix_long(patient_meth_path, genes, "METH", value_transform = identity) %>%
  rename(cell_or_sample = sample)

immune_ref <- read_csv(immune_ref_path, show_col_types = FALSE) %>%
  transmute(
    cell_or_sample = Sample,
    TIL_score = as.numeric(TIL_score),
    TMB = as.numeric(tmb_count)
  )

patient_features <- bind_rows(patient_rna_long, patient_meth_long) %>%
  select(omics, gene, cell_or_sample, feature_value) %>%
  left_join(candidate, by = "gene")

immune_long <- patient_features %>%
  inner_join(immune_ref, by = "cell_or_sample") %>%
  pivot_longer(c(TIL_score, TMB), names_to = "immune_metric", values_to = "immune_value") %>%
  group_by(omics, gene, Family6, original_layer, target_label, gene_label, immune_metric) %>%
  group_modify(~ safe_spearman(.x$feature_value, .x$immune_value, min_n = 10)) %>%
  ungroup() %>%
  group_by(omics, immune_metric) %>%
  mutate(fdr = p.adjust(p, method = "BH")) %>%
  ungroup()

immune_wide <- immune_long %>%
  select(omics, gene, immune_metric, n, rho, p, fdr) %>%
  pivot_wider(
    names_from = immune_metric,
    values_from = c(n, rho, p, fdr),
    names_glue = "{.value}_{immune_metric}"
  ) %>%
  left_join(candidate_layer_grid, by = c("omics", "gene")) %>%
  mutate(
    rho_TIL_score = as.numeric(rho_TIL_score),
    rho_TMB = as.numeric(rho_TMB),
    immune_strength = sqrt(coalesce(rho_TIL_score, 0)^2 + coalesce(rho_TMB, 0)^2),
    immune_best_fdr = pmin(coalesce(fdr_TIL_score, 1), coalesce(fdr_TMB, 1)),
    immune_tier = case_when(
      immune_best_fdr < 0.10 ~ "FDR < 0.10",
      immune_best_fdr < 0.25 ~ "FDR < 0.25",
      TRUE ~ "Exploratory"
    ),
    immune_quadrant = case_when(
      rho_TIL_score >= 0 & rho_TMB >= 0 ~ "Immune-hot / TMB-high",
      rho_TIL_score < 0 & rho_TMB < 0 ~ "Immune-cold / TMB-low",
      rho_TIL_score >= 0 & rho_TMB < 0 ~ "TIL-high / TMB-low",
      rho_TIL_score < 0 & rho_TMB >= 0 ~ "TMB-shifted",
      TRUE ~ "Not tested"
    ),
    layer_label = factor(
      recode(omics, RNA = "RNA feature", METH = "Methylation feature"),
      levels = c("RNA feature", "Methylation feature")
    )
  )

write_csv(immune_long, file.path(table_dir, "immune_all_candidate_til_tmb_correlations_by_omics_long.csv"))
write_csv(immune_wide, file.path(table_dir, "immune_all_candidate_til_tmb_correlations_by_omics_wide.csv"))

dep_xmax <- max(1, ceiling(max(abs(depmap_best_full$rho), na.rm = TRUE) * 10) / 10)
dep_ymax <- max(1.75, ceiling(max(depmap_best_full$plot_y, na.rm = TRUE) * 10) / 10 + 0.15)

dep_plot_data <- depmap_best_full %>%
  mutate(
    shape_layer = if_else(omics == "RNA", "RNA", "METH"),
    missing_label = if_else(!feature_available, gene, NA_character_),
    family_short = recode(
      Family6,
      "Immune" = "Imm",
      "Repair" = "Rep",
      "Glycolysis / TCA" = "Gly/TCA",
      "Fatty acid" = "FA",
      "Kinase signaling" = "Kin",
      "Hormone signaling" = "Horm",
      .default = Family6
    ),
    drug_label = str_trunc(coalesce(drug_clean, "not tested"), width = 18),
    drug_class_plot = case_when(
      str_detect(coalesce(drug_target_class, ""), "DDR|Checkpoint") ~ "DDR / Checkpoint",
      str_detect(coalesce(drug_target_class, ""), "RTK|MAPK|kinase") ~ "RTK/MAPK kinase",
      str_detect(coalesce(drug_target_class, ""), "Mitotic|Cell-cycle") ~ "Mitotic / Cell-cycle",
      str_detect(coalesce(drug_target_class, ""), "PI3K|AKT|mTOR") ~ "PI3K/AKT/mTOR",
      str_detect(coalesce(drug_target_class, ""), "Apoptosis|p53") ~ "Apoptosis / p53",
      str_detect(coalesce(drug_target_class, ""), "Epigenetic|chromatin") ~ "Epigenetic / chromatin",
      str_detect(coalesce(drug_target_class, ""), "Cytotoxic|chemotherapy") ~ "Cytotoxic / chemotherapy",
      str_detect(coalesce(drug_target_class, ""), "Metabolic|WNT") ~ "Metabolic / WNT",
      TRUE ~ "Bioactive / experimental"
    ),
    drug_class_plot = factor(drug_class_plot, levels = names(drug_shape_values)),
    hit_tier = case_when(
      coalesce(fdr, 1) < 0.10 ~ "FDR < 0.10",
      coalesce(fdr, 1) < 0.25 ~ "FDR < 0.25",
      TRUE ~ "Exploratory"
    )
  ) %>%
  group_by(layer_label) %>%
  arrange(desc(coalesce(plot_y, -Inf)), desc(coalesce(effect_abs, -Inf)), .by_group = TRUE) %>%
  mutate(label_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    label = if_else(
      feature_available & is.finite(rho),
      paste0(gene, " [", family_short, "]\n", drug_label),
      NA_character_
    )
  )
dep_label_data <- dep_plot_data %>% filter(!is.na(label))

dep_landscape <- ggplot(dep_plot_data, aes(x = rho, y = plot_y)) +
  annotate("rect", xmin = -dep_xmax, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = protective_bg, alpha = 0.72) +
  annotate("rect", xmin = 0, xmax = dep_xmax, ymin = -Inf, ymax = Inf,
           fill = adverse_bg, alpha = 0.72) +
  geom_vline(xintercept = 0, colour = "#8FA3B8", linetype = "longdash", linewidth = 0.55) +
  geom_hline(yintercept = -log10(0.10), colour = "#F87171", linetype = "longdash", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.25), colour = "#FDBA74", linetype = "dotted", linewidth = 0.55) +
  geom_point(
    aes(fill = Family6, colour = Family6, shape = drug_class_plot, alpha = hit_tier),
    stroke = 0.50, size = 3.05, na.rm = TRUE
  ) +
  geom_text_repel(
    data = dep_label_data,
    aes(label = label, colour = Family6),
    size = 1.72,
    fontface = "bold",
    lineheight = 0.82,
    box.padding = 0.14,
    point.padding = 0.09,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    segment.size = 0.16,
    max.overlaps = Inf,
    force = 2.25,
    max.time = 4,
    seed = 20260513,
    na.rm = TRUE
  ) +
  annotate("text", x = -dep_xmax * 0.92, y = dep_ymax * 0.93, label = "More sensitive",
           colour = layer_cols[["RNA"]], hjust = 0, size = 2.95, fontface = "bold") +
  annotate("text", x = dep_xmax * 0.92, y = dep_ymax * 0.93, label = "More resistant",
           colour = layer_cols[["METH"]], hjust = 1, size = 2.95, fontface = "bold") +
  facet_wrap(~ layer_label, ncol = 2) +
  scale_fill_manual(values = family_cols, drop = FALSE) +
  scale_colour_manual(values = family_cols, guide = "none", drop = FALSE) +
  scale_shape_manual(values = drug_shape_values, drop = FALSE, name = "Drug class") +
  scale_alpha_manual(values = c("FDR < 0.10" = 0.95, "FDR < 0.25" = 0.8, "Exploratory" = 0.48),
                     guide = "none") +
  coord_cartesian(xlim = c(-dep_xmax, dep_xmax), ylim = c(-0.03, dep_ymax), clip = "off") +
  labs(
    title = "DepMap drug-response landscape by omics layer",
    subtitle = "All 26 candidate genes were re-screened independently as RNA-expression and TSS-methylation features; labels show biomarker class and the best associated drug.",
    x = "Signed DepMap association with drug response",
    y = expression(-log[10]("best FDR across AUC / IC50")),
    fill = "Biomarker family"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 3.2, alpha = 1), nrow = 1, order = 1),
    shape = guide_legend(override.aes = list(colour = "#334155", fill = "grey88", size = 3.2, alpha = 1), nrow = 2, order = 2)
  ) +
  theme_pub(9.6) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key.size = unit(0.38, "cm"),
    legend.text = element_text(size = 7.6, colour = muted),
    legend.title = element_text(size = 8.2, face = "bold", colour = ink),
    panel.spacing = unit(1.1, "lines"),
    plot.title = element_text(size = 13.2, face = "bold", colour = ink),
    plot.subtitle = element_text(size = 8.4, colour = muted),
    plot.title.position = "plot"
  )

dep_missing <- dep_plot_data %>%
  filter(!feature_available) %>%
  mutate(chip_label = paste0(gene, " [", omics, "]"))

dep_caption_text <- if (nrow(dep_missing) > 0) {
  str_wrap(
    paste0(
      "DepMap feature unavailable and therefore not tested: ",
      paste(dep_missing$chip_label, collapse = ", "),
      "."
    ),
    width = 145
  )
} else {
  NULL
}
dep_final <- dep_landscape +
  plot_annotation(caption = dep_caption_text) &
  theme(plot.caption = element_text(size = 7.2, colour = muted, hjust = 0, margin = margin(t = 4)))

ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v1.png"),
       dep_final, width = 16.5, height = 9.2, dpi = 320, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v1.pdf"),
       dep_final, width = 16.5, height = 9.2, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v2.png"),
       dep_final, width = 13.8, height = 6.6, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v2.pdf"),
       dep_final, width = 13.8, height = 6.6, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v3.png"),
       dep_final, width = 13.6, height = 6.2, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v3.pdf"),
       dep_final, width = 13.6, height = 6.2, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v4_clean.png"),
       dep_final, width = 12.2, height = 5.7, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v4_clean.pdf"),
       dep_final, width = 12.2, height = 5.7, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v5_drug_labels.png"),
       dep_final, width = 13.4, height = 5.9, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v5_drug_labels.pdf"),
       dep_final, width = 13.4, height = 5.9, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v6_drug_class_shape.png"),
       dep_final, width = 13.8, height = 6.15, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v6_drug_class_shape.pdf"),
       dep_final, width = 13.8, height = 6.15, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v7_refined_drug_classes.png"),
       dep_final, width = 14.2, height = 6.35, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v7_refined_drug_classes.pdf"),
       dep_final, width = 14.2, height = 6.35, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v8_all_labels.png"),
       dep_final, width = 14.8, height = 6.75, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02B_DepMap_all_candidates_RNA_METH_rescreen_R_v8_all_labels.pdf"),
       dep_final, width = 14.8, height = 6.75, bg = "white", device = cairo_pdf)

immune_plot_data <- immune_wide %>%
  mutate(
    shape_layer = if_else(omics == "RNA", "RNA", "METH"),
    label = gene
  ) %>%
  group_by(layer_label) %>%
  arrange(desc(immune_strength), .by_group = TRUE) %>%
  mutate(label_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    label = gene
  )
immune_label_data <- immune_plot_data %>% filter(!is.na(label))

immune_ranges <- immune_plot_data %>%
  group_by(layer_label) %>%
  summarise(
    x_min_raw = min(rho_TIL_score, na.rm = TRUE),
    x_max_raw = max(rho_TIL_score, na.rm = TRUE),
    y_min_raw = min(rho_TMB, na.rm = TRUE),
    y_max_raw = max(rho_TMB, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_pad = pmax((x_max_raw - x_min_raw) * 0.10, 0.045),
    y_pad = pmax((y_max_raw - y_min_raw) * 0.14, 0.045),
    x_min = pmin(x_min_raw, 0) - x_pad,
    x_max = pmax(x_max_raw, 0) + x_pad,
    y_min = pmin(y_min_raw, 0) - y_pad,
    y_max = pmax(y_max_raw, 0) + y_pad
  )

immune_hot_rect <- immune_ranges %>% transmute(layer_label, xmin = 0, xmax = x_max, ymin = 0, ymax = y_max)
immune_cold_rect <- immune_ranges %>% transmute(layer_label, xmin = x_min, xmax = 0, ymin = y_min, ymax = 0)
immune_tmb_rect <- immune_ranges %>% transmute(layer_label, xmin = x_min, xmax = 0, ymin = 0, ymax = y_max)
immune_low_rect <- immune_ranges %>% transmute(layer_label, xmin = 0, xmax = x_max, ymin = y_min, ymax = 0)

immune_quad_labels <- bind_rows(
  immune_ranges %>%
    transmute(layer_label, label = "Immune-hot / TMB-high",
              x = x_max - x_pad * 0.5, y = y_max - y_pad * 0.65,
              hjust = 1, colour = "#47C7D9"),
  immune_ranges %>%
    transmute(layer_label, label = "Immune-cold / TMB-low",
              x = x_min + x_pad * 0.5, y = y_min + y_pad * 0.65,
              hjust = 0, colour = "#8BC56C"),
  immune_ranges %>%
    transmute(layer_label, label = "TMB-shifted",
              x = x_min + x_pad * 0.5, y = y_max - y_pad * 0.65,
              hjust = 0, colour = "#8FA9C9")
)

immune_landscape <- ggplot(immune_plot_data, aes(x = rho_TIL_score, y = rho_TMB)) +
  geom_rect(data = immune_hot_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#E7FAF7", alpha = 0.82) +
  geom_rect(data = immune_cold_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#EFF8EE", alpha = 0.82) +
  geom_rect(data = immune_tmb_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#EEF5FF", alpha = 0.72) +
  geom_rect(data = immune_low_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "#F8FAFC", alpha = 0.74) +
  geom_vline(xintercept = 0, colour = "#8FA3B8", linetype = "longdash", linewidth = 0.55) +
  geom_hline(yintercept = 0, colour = "#8FA3B8", linetype = "longdash", linewidth = 0.55) +
  geom_point(
    aes(fill = Family6, shape = shape_layer, alpha = immune_tier),
    size = 3.15,
    colour = "black", stroke = 0.45
  ) +
  geom_text_repel(
    data = immune_label_data,
    aes(label = label, colour = Family6),
    size = 1.78,
    fontface = "bold",
    box.padding = 0.14,
    point.padding = 0.09,
    min.segment.length = 0,
    segment.colour = "#94A3B8",
    segment.size = 0.16,
    max.overlaps = Inf,
    force = 2.15,
    max.time = 4,
    seed = 20260513
  ) +
  geom_text(
    data = immune_quad_labels,
    aes(x = x, y = y, label = label, hjust = hjust, colour = I(colour)),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 2.9
  ) +
  facet_wrap(~ layer_label, ncol = 2, scales = "free") +
  scale_fill_manual(values = family_cols, drop = FALSE) +
  scale_colour_manual(values = family_cols, guide = "none", drop = FALSE) +
  scale_shape_manual(values = c("RNA" = 21, "METH" = 22), guide = "none") +
  scale_alpha_manual(values = c("FDR < 0.10" = 0.96, "FDR < 0.25" = 0.78, "Exploratory" = 0.48),
                     guide = "none") +
  coord_cartesian(clip = "off") +
  labs(
    title = "TIL/TMB immune landscape by omics layer",
    subtitle = "All 26 candidate genes were re-screened independently as patient RNA and methylation features in TCGA-overlapping samples.",
    x = "Spearman rho with TIL score",
    y = "Spearman rho with TMB",
    fill = "Biomarker family"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size = 3.4, alpha = 1), nrow = 1)
  ) +
  theme_pub(9.6) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key.size = unit(0.38, "cm"),
    legend.text = element_text(size = 7.6, colour = muted),
    legend.title = element_text(size = 8.2, face = "bold", colour = ink),
    panel.spacing = unit(1.1, "lines"),
    plot.title = element_text(size = 13.2, face = "bold", colour = ink),
    plot.subtitle = element_text(size = 8.4, colour = muted),
    plot.title.position = "plot"
  )

ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v1.png"),
       immune_landscape, width = 16.5, height = 9.1, dpi = 320, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v1.pdf"),
       immune_landscape, width = 16.5, height = 9.1, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v2.png"),
       immune_landscape, width = 13.8, height = 6.4, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v2.pdf"),
       immune_landscape, width = 13.8, height = 6.4, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v3.png"),
       immune_landscape, width = 13.6, height = 6.1, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v3.pdf"),
       immune_landscape, width = 13.6, height = 6.1, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v4_clean.png"),
       immune_landscape, width = 12.2, height = 5.6, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v4_clean.pdf"),
       immune_landscape, width = 12.2, height = 5.6, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v5_corr_axis.png"),
       immune_landscape, width = 12.2, height = 5.6, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v5_corr_axis.pdf"),
       immune_landscape, width = 12.2, height = 5.6, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v7_free_axis.png"),
       immune_landscape, width = 12.2, height = 5.6, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v7_free_axis.pdf"),
       immune_landscape, width = 12.2, height = 5.6, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v8_all_labels.png"),
       immune_landscape, width = 12.8, height = 5.9, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03B_Immune_all_candidates_RNA_METH_rescreen_R_v8_all_labels.pdf"),
       immune_landscape, width = 12.8, height = 5.9, bg = "white", device = cairo_pdf)

if (file.exists(legacy_immune_path)) {
  legacy_immune <- read_csv(legacy_immune_path, show_col_types = FALSE) %>%
    mutate(
      immune_available = as.logical(immune_available),
      Layer = factor(Layer, levels = c("Transcriptome", "Methylation")),
      Family6 = factor(Family6, levels = names(family_cols)),
      short_label = paste0(gene, " [", if_else(as.character(Layer) == "Methylation", "M", "R"), "]"),
      immune_signif_tier = factor(immune_signif_tier, levels = c("FDR < 0.10", "FDR < 0.25", "Exploratory"))
    )

  legacy_available <- legacy_immune %>% filter(immune_available)
  legacy_missing <- legacy_immune %>% filter(!immune_available)

  x_rng <- range(legacy_available$rho_til, na.rm = TRUE)
  y_rng <- range(legacy_available$rho_tmb, na.rm = TRUE)
  x_pad <- 0.06
  y_pad <- 0.06

  legacy_immune_plot <- ggplot(legacy_available, aes(rho_til, rho_tmb)) +
    annotate("rect", xmin = x_rng[1] - x_pad, xmax = 0, ymin = 0, ymax = y_rng[2] + y_pad,
             fill = "#EEF5FF", alpha = 0.72) +
    annotate("rect", xmin = 0, xmax = x_rng[2] + x_pad, ymin = 0, ymax = y_rng[2] + y_pad,
             fill = "#E7FAF7", alpha = 0.82) +
    annotate("rect", xmin = x_rng[1] - x_pad, xmax = 0, ymin = y_rng[1] - y_pad, ymax = 0,
             fill = "#EFF8EE", alpha = 0.82) +
    annotate("rect", xmin = 0, xmax = x_rng[2] + x_pad, ymin = y_rng[1] - y_pad, ymax = 0,
             fill = "#F8FAFC", alpha = 0.74) +
    geom_vline(xintercept = 0, colour = "#8FA3B8", linetype = "longdash", linewidth = 0.55) +
    geom_hline(yintercept = 0, colour = "#8FA3B8", linetype = "longdash", linewidth = 0.55) +
    geom_point(
      aes(fill = Family6, shape = Layer, size = immune_size_norm, alpha = immune_signif_tier),
      colour = "black", stroke = 0.65
    ) +
    geom_label_repel(
      aes(label = short_label, colour = Family6),
      fill = alpha("white", 0.92),
      label.size = 0.14,
      size = 2.45,
      fontface = "bold",
      box.padding = 0.24,
      point.padding = 0.16,
      segment.size = 0,
      max.overlaps = Inf,
      seed = 20260513,
      show.legend = FALSE
    ) +
    annotate("text", x = x_rng[1] - 0.02, y = y_rng[2] + 0.035,
             label = "TMB-shifted", hjust = 0, colour = "#90A7C8", size = 3.1, fontface = "bold") +
    annotate("text", x = 0.12, y = y_rng[2] + 0.035,
             label = "Immune-hot / TMB-high", hjust = 0, colour = "#47C7D9", size = 3.1, fontface = "bold") +
    annotate("text", x = -0.015, y = y_rng[1] + 0.035,
             label = "Immune-cold / TMB-low", hjust = 1, colour = "#8BC56C", size = 2.9, fontface = "bold") +
    scale_fill_manual(values = family_cols, drop = TRUE) +
    scale_colour_manual(values = family_cols, guide = "none", drop = TRUE) +
    scale_shape_manual(values = layer_shape_values, drop = FALSE) +
    scale_alpha_manual(values = c("FDR < 0.10" = 0.96, "FDR < 0.25" = 0.78, "Exploratory" = 0.50),
                       drop = FALSE) +
    scale_size_area(max_size = 8.8, limits = c(0, 1), breaks = c(0.5, 0.75, 1.0)) +
    coord_cartesian(
      xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
      ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
      clip = "off"
    ) +
    labs(
      title = "Immune-evidence TIL/TMB landscape",
      subtitle = paste0(
        "Layer-specific immune-correlation evidence restored from the prior analysis (",
        nrow(legacy_available),
        " biomarkers with usable TIL/TMB evidence; ",
        nrow(legacy_missing),
        " retained as not-carried candidates)."
      ),
      x = "Spearman rho with TIL score",
      y = "Spearman rho with TMB",
      fill = "Biomarker family",
      shape = "Evidence layer",
      alpha = "Evidence tier",
      size = "Immune strength"
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = 21, size = 3.4, alpha = 1), order = 1),
      shape = guide_legend(override.aes = list(size = 3.4, fill = "grey90", alpha = 1), order = 2),
      alpha = guide_legend(override.aes = list(shape = 21, fill = "grey90", size = 3.4), order = 3),
      size = guide_legend(override.aes = list(shape = 21, fill = "grey90", alpha = 1), order = 4)
    ) +
    theme_pub(9.8) +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.spacing.y = unit(0.12, "cm"),
      legend.key.size = unit(0.42, "cm"),
      legend.text = element_text(size = 7.8, colour = muted),
      legend.title = element_text(size = 8.5, face = "bold", colour = ink),
      plot.title = element_text(size = 13.5, face = "bold", colour = ink),
      plot.subtitle = element_text(size = 8.2, colour = muted),
      plot.title.position = "plot"
    )

  ggsave(file.path(plot_dir, "Figure_03B_Immune_evidence_restored_landscape_R_v6.png"),
         legacy_immune_plot, width = 10.8, height = 6.1, dpi = 360, bg = "white")
  ggsave(file.path(plot_dir, "Figure_03B_Immune_evidence_restored_landscape_R_v6.pdf"),
         legacy_immune_plot, width = 10.8, height = 6.1, bg = "white", device = cairo_pdf)

  immune_compare <- legacy_available %>%
    mutate(omics = if_else(as.character(Layer) == "Methylation", "METH", "RNA")) %>%
    select(omics, gene, Layer, Family6, prior_rho_til = rho_til, prior_rho_tmb = rho_tmb,
           prior_best_fdr = immune_best_fdr) %>%
    left_join(
      immune_wide %>% select(omics, gene, rescreen_rho_til = rho_TIL_score,
                             rescreen_rho_tmb = rho_TMB, rescreen_best_fdr = immune_best_fdr),
      by = c("gene", "omics")
    )
  write_csv(immune_compare, file.path(table_dir, "immune_prior_evidence_vs_all_candidate_rescreen.csv"))
}

dep_summary <- depmap_best_full %>%
  group_by(omics) %>%
  summarise(
    n_candidates = n(),
    n_tested = sum(feature_available, na.rm = TRUE),
    n_fdr_10 = sum(fdr < 0.10, na.rm = TRUE),
    n_fdr_25 = sum(fdr < 0.25, na.rm = TRUE),
    strongest_gene = gene[which.max(if_else(is.finite(plot_y), plot_y, -Inf))],
    strongest_drug = drug_clean[which.max(if_else(is.finite(plot_y), plot_y, -Inf))],
    .groups = "drop"
  )
immune_summary <- immune_wide %>%
  group_by(omics) %>%
  summarise(
    n_candidates = n(),
    n_til_fdr_10 = sum(fdr_TIL_score < 0.10, na.rm = TRUE),
    n_tmb_fdr_10 = sum(fdr_TMB < 0.10, na.rm = TRUE),
    n_any_fdr_10 = sum(immune_best_fdr < 0.10, na.rm = TRUE),
    strongest_gene = gene[which.max(immune_strength)],
    .groups = "drop"
  )
write_csv(dep_summary, file.path(table_dir, "depmap_layer_rescreen_summary.csv"))
write_csv(immune_summary, file.path(table_dir, "immune_layer_rescreen_summary.csv"))

message("Done.")
message("Plots written to: ", plot_dir)
message("Tables written to: ", table_dir)
