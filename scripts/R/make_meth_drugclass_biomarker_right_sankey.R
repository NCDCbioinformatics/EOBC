suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggalluvial)
  library(scales)
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
base_dir <- env_path("EOBC_FINAL_ANALYSIS_DIR", file.path(biomarker_root, "final_analysis"))
context_dir <- file.path(base_dir, "group_biomarker_contextual_domain_raw_meth_r_v1", "tables")
validation_dir <- file.path(base_dir, "full_layer_candidate_validation_r_v1", "tables")
out_dir <- file.path(base_dir, "eobc_integrated_story_sankey_r_v5")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
final_fig_dir <- file.path(base_dir, "final_fig")
final_revision_dir <- file.path(final_fig_dir, "final_paper_figures_20260522_v2")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_revision_dir, recursive = TRUE, showWarnings = FALSE)

dominance_path <- file.path(context_dir, "activity_scale_group_dominance_with_domain_evidence.csv")
alignment_path <- file.path(context_dir, "rna_meth_alignment_reused_for_contextual_figures.csv")
depmap_path <- file.path(validation_dir, "depmap_all_candidate_all_drug_correlations_by_omics.csv")
expr_immune_path <- file.path(base_dir, "group_biomarker_landscape_r_v1", "tables", "group_biomarker_long_values.csv")
immune_score_path <- file.path(
  dirname(base_dir),
  "10_external_validation_TIL_TMB_DepMap_final_v7",
  "tables",
  "immune_signature_tcga_til_tmb_merged.csv"
)
gene_rf_path <- file.path(
  dirname(base_dir),
  "05_model_output_archive",
  "biomarker_gene_moment_RF_SEM_manual_groups_v2",
  "tables",
  "gene_driver_summary_all_targets.csv"
)

group_order <- c("G1 | H", "G2 | I", "G3 | L-like", "G4 | L")
group_display <- c(
  "G1 | H" = "G1 | H\nHypermethylated\nEOBC state",
  "G2 | I" = "G2 | I\nIntermediate\nEOBC state",
  "G3 | L-like" = "G3 | L-like\nLuminal-like\nEOBC state",
  "G4 | L" = "G4 | L\nLuminal\nEOBC state"
)

group_node_cols <- c(
  "G1 | H" = "#D81B60",
  "G2 | I" = "#F59E0B",
  "G3 | L-like" = "#00A6A6",
  "G4 | L" = "#7CB342",
  "G1 | H\nHypermethylated\nEOBC state" = "#D81B60",
  "G2 | I\nIntermediate\nEOBC state" = "#F59E0B",
  "G3 | L-like\nLuminal-like\nEOBC state" = "#00A6A6",
  "G4 | L\nLuminal\nEOBC state" = "#7CB342"
)

family_cols <- c(
  "Immune" = "#2B8CBE",
  "Repair" = "#2CA25F",
  "Glycolysis / TCA" = "#D4A017",
  "Fatty acid" = "#E67E22",
  "Kinase signaling" = "#7B61D1",
  "Hormone signaling" = "#7F8FA6"
)

evidence_cols <- c(
  "OS | Protective" = "#58A9F9",
  "OS | Adverse" = "#F4A259",
  "Drug | DDR/checkpoint\nSensitive" = "#2563EB",
  "Drug | DDR/checkpoint\nResistant" = "#D95F02",
  "Drug | RTK/MAPK\nSensitive" = "#2563EB",
  "Drug | RTK/MAPK\nResistant" = "#D95F02",
  "Drug | Mitotic/cell-cycle\nSensitive" = "#2563EB",
  "Drug | Mitotic/cell-cycle\nResistant" = "#D95F02",
  "Drug | Cytotoxic/DNA damage\nSensitive" = "#2563EB",
  "Drug | Cytotoxic/DNA damage\nResistant" = "#D95F02",
  "Drug | PI3K/mTOR\nSensitive" = "#2563EB",
  "Drug | PI3K/mTOR\nResistant" = "#D95F02",
  "Drug | p53/MDM2\nSensitive" = "#2563EB",
  "Drug | p53/MDM2\nResistant" = "#D95F02",
  "Immune | TIL/TMB positive" = "#0E9F6E",
  "Immune | TIL/TMB weak" = "#B9C4D1",
  "Immune | TIL/TMB negative" = "#4B8B3B",
  "Immune | TIL rho positive" = "#17A2A4",
  "Immune | TIL rho negative" = "#72A950",
  "Immune | TMB rho positive" = "#7B61D1",
  "Immune | TMB rho negative" = "#6F8FB7",
  "Immune | TIL/TMB rho positive" = "#0E9F6E",
  "Immune | TIL/TMB rho negative" = "#4B8B3B",
  "Immune | Discordant TIL/TMB rho" = "#B9C4D1"
)

biomarker_set_cols <- c(
  "Multi-domain" = "#6D28D9",
  "OS + Immune" = "#DB2777",
  "OS + DepMap" = "#F97316",
  "DepMap + Immune" = "#059669",
  "OS-only" = "#EF4444",
  "DepMap-only" = "#2563EB",
  "Immune-only" = "#06B6D4"
)

rf_dom_cols <- c(
  "RF M-dominant" = "#0F766E",
  "RF R-dominant" = "#2563EB",
  "RF M/R-balanced" = "#7C3AED",
  "RF n/a" = "#94A3B8"
)

meth_rna_cols <- c(
  "Meth-RNA inverse" = "#3B82F6",
  "Meth-RNA positive" = "#F97316",
  "Meth-RNA weak" = "#CBD5E1"
)

pick_col <- function(palette, key, default = "#94A3B8") {
  key <- as.character(key)
  out <- unname(palette[key])
  if (length(out) == 0 || is.na(out)) default else out
}

mix_cols <- function(cols, weights) {
  keep <- !is.na(cols) & nzchar(cols) & is.finite(weights) & weights > 0
  cols <- cols[keep]
  weights <- weights[keep]
  if (length(cols) == 0) return("#94A3B8")
  weights <- weights / sum(weights)
  rgb_mat <- grDevices::col2rgb(cols)
  rgb_mix <- as.numeric(rgb_mat %*% weights)
  grDevices::rgb(rgb_mix[1], rgb_mix[2], rgb_mix[3], maxColorValue = 255)
}

boost_col <- function(col, sat = 1.22, val = 1.03) {
  hsv_mat <- grDevices::rgb2hsv(grDevices::col2rgb(col))
  grDevices::hsv(
    h = hsv_mat["h", ],
    s = pmin(1, hsv_mat["s", ] * sat),
    v = pmin(1, hsv_mat["v", ] * val)
  )
}

evidence_color_for_node <- function(x) {
  x <- as.character(x)
  case_when(
    str_detect(x, "^OS \\| Protective") ~ "#58A9F9",
    str_detect(x, "^OS \\| Adverse") ~ "#F4A259",
    str_detect(x, "^Drug \\| .*Sensitive") ~ "#2563EB",
    str_detect(x, "^Drug \\| .*Resistant") ~ "#D95F02",
    str_detect(x, "^Immune \\| TIL/TMB positive") ~ "#0E9F6E",
    str_detect(x, "^Immune \\| TIL/TMB weak") ~ "#B9C4D1",
    str_detect(x, "^Immune \\| TIL/TMB negative") ~ "#4B8B3B",
    str_detect(x, "^Immune \\| TIL rho positive") ~ "#17A2A4",
    str_detect(x, "^Immune \\| TIL rho negative") ~ "#72A950",
    str_detect(x, "^Immune \\| TMB rho positive") ~ "#7B61D1",
    str_detect(x, "^Immune \\| TMB rho negative") ~ "#6F8FB7",
    str_detect(x, "^Immune \\| TIL/TMB rho positive") ~ "#0E9F6E",
    str_detect(x, "^Immune \\| TIL/TMB rho negative") ~ "#4B8B3B",
    str_detect(x, "^Immune \\| Discordant") ~ "#B9C4D1",
    TRUE ~ "#94A3B8"
  )
}

blend_route_color <- function(group, family, evidence, biomarker_set, rf_dom) {
  mapply(
    function(g, f, e, b, r) {
      boost_col(mix_cols(
        c(
          pick_col(group_node_cols, g),
          pick_col(family_cols, f),
          evidence_color_for_node(e),
          pick_col(biomarker_set_cols, b),
          pick_col(rf_dom_cols, r)
        ),
        c(0.20, 0.30, 0.28, 0.12, 0.10)
      ), sat = 1.28, val = 1.04)
    },
    as.character(group),
    as.character(family),
    as.character(evidence),
    as.character(biomarker_set),
    as.character(rf_dom),
    USE.NAMES = FALSE
  )
}

route_stage_colors <- function(group, family, evidence, biomarker_set, rf_dom, meth_rna_class) {
  cols <- mapply(
    function(g, f, e, b, r, m) {
      group_col <- pick_col(group_node_cols, g)
      family_col <- pick_col(family_cols, f)
      evidence_col <- evidence_color_for_node(e)
      biomarker_col <- boost_col(mix_cols(
        c(pick_col(meth_rna_cols, m), pick_col(rf_dom_cols, r), pick_col(biomarker_set_cols, b)),
        c(0.44, 0.34, 0.22)
      ), sat = 1.12, val = 1.02)

      c(
        stage1_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col),
          c(0.88, 0.08, 0.04)
        ), sat = 1.08, val = 1.02),
        stage2_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col, evidence_col),
          c(0.44, 0.46, 0.08, 0.02)
        ), sat = 1.10, val = 1.02),
        stage3_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col, evidence_col),
          c(0.18, 0.28, 0.42, 0.12)
        ), sat = 1.13, val = 1.02),
        stage4_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col, evidence_col),
          c(0.08, 0.16, 0.26, 0.50)
        ), sat = 1.15, val = 1.02)
      )
    },
    as.character(group),
    as.character(family),
    as.character(evidence),
    as.character(biomarker_set),
    as.character(rf_dom),
    as.character(meth_rna_class),
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )

  tibble::as_tibble(t(cols), .name_repair = "minimal") %>%
    mutate(across(everything(), as.character))
}

truthy <- function(x) {
  if (is.logical(x)) {
    return(replace_na(x, FALSE))
  }
  str_to_upper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
}

clean_label <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

pretty_pathway_label <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x == "Base Excision\nRepair" ~ "DNA-repair axis\nBase-excision repair",
    x == "Cell Cycle\nCheckpoint" ~ "DNA-repair axis\nCell-cycle checkpoint",
    x == "Checkpoint\nInhibition" ~ "DNA-damage axis\nCheckpoint inhibition",
    x == "Chromatin\nRemodeling" ~ "Epigenetic-repair axis\nChromatin remodeling",
    x == "DNA Damage\nResponse" ~ "DNA-repair axis\nDNA-damage response",
    x == "Homologous\nRecombination" ~ "DNA-repair axis\nHomologous recombination",
    x == "Mismatch Repair" ~ "DNA-repair axis\nMismatch repair",
    x == "Fatty Acid\nMetabolism" ~ "Metabolic axis\nFatty-acid metabolism",
    x == "Glycolysis" ~ "Metabolic axis\nGlycolysis",
    x == "TCA Cycle" ~ "Metabolic axis\nTCA cycle",
    x == "mTOR Signaling" ~ "Oncogenic signaling\nmTOR signaling",
    x == "PI3K-AKT Pathway" ~ "Oncogenic signaling\nPI3K-AKT pathway",
    x == "RAS-MAPK Pathway" ~ "Oncogenic signaling\nRAS-MAPK pathway",
    x == "WNT Signaling" ~ "Developmental signaling\nWNT signaling",
    x == "Hormone\nSignaling" ~ "Hormone axis\nHormone signaling",
    str_detect(x, "^Immune \\| ") ~ str_replace(x, "^Immune \\| ", "Immune axis\n"),
    TRUE ~ str_replace_all(x, " \\| ", "\n")
  )
}

pretty_evidence_label <- function(x) {
  out <- as.character(x)
  out[out == "OS | Protective"] <- "Clinical outcome\nFavorable OS"
  out[out == "OS | Adverse"] <- "Clinical outcome\nAdverse OS"
  immune_idx <- str_detect(out, "^Immune \\| ")
  out[immune_idx] <- out[immune_idx] %>%
    str_replace("^Immune \\| ", "TIL/TMB immune axis\n") %>%
    str_replace("TIL/TMB positive", "Concordant positive") %>%
    str_replace("TIL/TMB weak", "Weak / discordant") %>%
    str_replace("TIL/TMB negative", "Concordant negative") %>%
    str_replace("TIL rho positive", "TIL rho positive") %>%
    str_replace("TIL rho negative", "TIL rho negative") %>%
    str_replace("TMB rho positive", "TMB rho positive") %>%
    str_replace("TMB rho negative", "TMB rho negative") %>%
    str_replace("TIL/TMB rho positive", "TIL/TMB rho positive") %>%
    str_replace("TIL/TMB rho negative", "TIL/TMB rho negative") %>%
    str_replace("Discordant TIL/TMB rho", "Discordant TIL/TMB rho")
  drug_idx <- str_detect(out, "^Drug \\| ")
  out[drug_idx] <- out[drug_idx] %>%
    str_replace("^Drug \\| ", "DepMap therapeutic evidence\n") %>%
    str_replace("DDR/checkpoint", "DDR / checkpoint") %>%
    str_replace("RTK/MAPK", "RTK / MAPK") %>%
    str_replace("Mitotic/cell-cycle", "Mitotic / cell-cycle") %>%
    str_replace("Cytotoxic/DNA damage", "Cytotoxic / DNA-damage") %>%
    str_replace("PI3K/mTOR", "PI3K / mTOR") %>%
    str_replace("p53/MDM2", "p53 / MDM2") %>%
    str_replace("Experimental/natural", "Experimental / natural product") %>%
    str_replace("Other/exploratory", "Exploratory / other") %>%
    str_replace("\nSensitive$", "\nSensitivity-associated") %>%
    str_replace("\nResistant$", "\nResistance-associated")
  out
}

pretty_biomarker_label <- function(x) {
  as.character(x) %>%
    str_replace_all("RNA inv strong", "Meth-RNA inverse") %>%
    str_replace_all("RNA inv\\b", "Meth-RNA inverse") %>%
    str_replace_all("RNA ctx\\b", "Meth-RNA positive") %>%
    str_replace_all("RNA weak\\b", "Meth-RNA weak")
}

immune_subclass_from_gene <- function(gene, target_label) {
  target_label <- clean_label(coalesce(target_label, ""))
  case_when(
    str_detect(target_label, regex("Adaptive", ignore_case = TRUE)) ~ "Immune | Adaptive Immunity",
    str_detect(target_label, regex("Innate", ignore_case = TRUE)) ~ "Immune | Innate Immunity",
    str_detect(target_label, regex("Antigen", ignore_case = TRUE)) ~ "Immune | Antigen Processing",
    str_detect(target_label, regex("B.?Cell|B cell", ignore_case = TRUE)) ~ "Immune | B Cell Markers",
    str_detect(target_label, regex("NK", ignore_case = TRUE)) ~ "Immune | NK Cell Markers",
    str_detect(target_label, regex("Macrophage", ignore_case = TRUE)) ~ "Immune | Macrophage Activation",
    str_detect(target_label, regex("Treg|Regulatory", ignore_case = TRUE)) ~ "Immune | Treg Markers",
    str_detect(target_label, regex("Cytokine", ignore_case = TRUE)) ~ "Immune | Cytokine Signaling",
    str_detect(target_label, regex("Apoptosis", ignore_case = TRUE)) ~ "Immune | Apoptosis Regulators",
    gene %in% c("CLEC7A", "MMP10", "TRAF3IP3") ~ "Immune | Innate Immunity",
    gene %in% c("MYO1G", "LRRC10B") ~ "Immune | Treg Markers",
    gene %in% c("RASAL3", "CPLX1", "KLHDC7B", "MLPH", "C1QTNF6") ~ "Immune | Cytokine Signaling",
    gene %in% c("TFF1", "BCL2A1") ~ "Immune | Macrophage Activation",
    TRUE ~ "Immune | Immune program"
  )
}

pathway_subclass <- function(family, gene, target_label) {
  target_label <- clean_label(coalesce(target_label, ""))
  out <- case_when(
    family == "Immune" ~ immune_subclass_from_gene(gene, target_label),
    family == "Repair" & str_detect(target_label, regex("Base Excision", ignore_case = TRUE)) ~ "Base Excision\nRepair",
    family == "Repair" & str_detect(target_label, regex("Cell Cycle", ignore_case = TRUE)) ~ "Cell Cycle\nCheckpoint",
    family == "Repair" & str_detect(target_label, regex("Checkpoint", ignore_case = TRUE)) ~ "Checkpoint\nInhibition",
    family == "Repair" & str_detect(target_label, regex("Chromatin", ignore_case = TRUE)) ~ "Chromatin\nRemodeling",
    family == "Repair" & str_detect(target_label, regex("DNA Damage", ignore_case = TRUE)) ~ "DNA Damage\nResponse",
    family == "Repair" & str_detect(target_label, regex("Homologous", ignore_case = TRUE)) ~ "Homologous\nRecombination",
    family == "Repair" ~ paste0("Repair | ", target_label),
    family == "Glycolysis / TCA" & str_detect(target_label, regex("TCA", ignore_case = TRUE)) ~ "TCA Cycle",
    family == "Glycolysis / TCA" & str_detect(target_label, regex("Glycol", ignore_case = TRUE)) ~ "Glycolysis",
    family == "Glycolysis / TCA" ~ "Glycolysis / TCA",
    family == "Fatty acid" ~ "Fatty Acid\nMetabolism",
    family == "Kinase signaling" & str_detect(target_label, regex("mTOR", ignore_case = TRUE)) ~ "mTOR Signaling",
    family == "Kinase signaling" & str_detect(target_label, regex("PI3K|AKT", ignore_case = TRUE)) ~ "PI3K-AKT Pathway",
    family == "Kinase signaling" & str_detect(target_label, regex("RAS|MAPK", ignore_case = TRUE)) ~ "RAS-MAPK Pathway",
    family == "Kinase signaling" ~ clean_label(target_label),
    family == "Hormone signaling" & str_detect(target_label, regex("WNT", ignore_case = TRUE)) ~ "WNT Signaling",
    family == "Hormone signaling" ~ "Hormone\nSignaling",
    TRUE ~ clean_label(target_label)
  )
  out %>%
    str_replace_all("Immune \\| ", "Immune | ") %>%
    str_replace_all(" / ", "/") %>%
    str_squish()
}

classify_drug <- function(drug) {
  d <- str_to_upper(coalesce(drug, ""))
  case_when(
    str_detect(d, "ADAVOSERTIB|CERALASERTIB|ATR|CHK|WEE1|CHECKPOINT") ~ "DDR/checkpoint",
    str_detect(d, "DASATINIB|ERLOTINIB|NILOTINIB|VX-11E|EGFR|ERK|MEK|RAF|MAPK|RTK") ~ "RTK/MAPK",
    str_detect(d, "ALISERTIB|AURORA|BI-2536|PLK|PALBOCICLIB|CDK|MK 0457|MK-0457|VELBAN|VINBLASTINE|MITOTIC") ~ "Mitotic/cell-cycle",
    str_detect(d, "DOCETAXEL|CARMUSTINE|MITOXANTRONE|TOPOTECAN|CAMPTOTHECIN|TAXANE|CHEMO") ~ "Cytotoxic/DNA damage",
    str_detect(d, "RAPAMYCIN|AZD-8055|MTOR|PI3K|AKT") ~ "PI3K/mTOR",
    str_detect(d, "REBEMADLIN|MDM2|P53") ~ "p53/MDM2",
    str_detect(d, "SINULARIN|ELEPHANTIN|5XE|284461|BDP-|AK175558|CRA-032765|AGI-6780") ~ "Experimental/natural",
    TRUE ~ "Other/exploratory"
  )
}

evidence_from_immune <- function(quadrant, status) {
  q <- str_to_lower(coalesce(status, quadrant, ""))
  case_when(
    str_detect(q, "til/tmb rho positive") ~ "Immune | TIL/TMB rho positive",
    str_detect(q, "til/tmb rho negative") ~ "Immune | TIL/TMB rho negative",
    str_detect(q, "discordant") ~ "Immune | Discordant TIL/TMB rho",
    str_detect(q, "til rho positive") ~ "Immune | TIL rho positive",
    str_detect(q, "til rho negative") ~ "Immune | TIL rho negative",
    str_detect(q, "tmb rho positive") ~ "Immune | TMB rho positive",
    str_detect(q, "tmb rho negative") ~ "Immune | TMB rho negative",
    TRUE ~ NA_character_
  )
}

evidence_from_til_tmb_pair <- function(til_rho, tmb_rho, weak_cutoff = 0.20) {
  case_when(
    is.finite(til_rho) & is.finite(tmb_rho) &
      til_rho > 0 & tmb_rho > 0 &
      pmin(abs(til_rho), abs(tmb_rho)) >= weak_cutoff ~ "Immune | TIL/TMB positive",
    is.finite(til_rho) & is.finite(tmb_rho) &
      til_rho < 0 & tmb_rho < 0 &
      pmin(abs(til_rho), abs(tmb_rho)) >= weak_cutoff ~ "Immune | TIL/TMB negative",
    is.finite(til_rho) | is.finite(tmb_rho) ~ "Immune | TIL/TMB weak",
    TRUE ~ NA_character_
  )
}

evidence_from_recomputed_til_tmb <- function(til_rho, tmb_rho, til_q, tmb_q, weak_cutoff = 0.15, q_cutoff = 0.25) {
  best_q <- pmin(coalesce(til_q, Inf), coalesce(tmb_q, Inf))
  case_when(
    is.finite(til_rho) & is.finite(tmb_rho) &
      til_rho > 0 & tmb_rho > 0 &
      pmin(abs(til_rho), abs(tmb_rho)) >= weak_cutoff &
      best_q <= q_cutoff ~ "Immune | TIL/TMB positive",
    is.finite(til_rho) & is.finite(tmb_rho) &
      til_rho < 0 & tmb_rho < 0 &
      pmin(abs(til_rho), abs(tmb_rho)) >= weak_cutoff &
      best_q <= q_cutoff ~ "Immune | TIL/TMB negative",
    best_q <= q_cutoff ~ "Immune | TIL/TMB weak",
    TRUE ~ NA_character_
  )
}

spearman_summary <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 4) {
    return(tibble(n = n, rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(n = n, rho = unname(ct$estimate), p = ct$p.value)
}

evidence_from_os <- function(status, direction) {
  s <- str_to_lower(coalesce(status, direction, ""))
  case_when(
    str_detect(s, "adverse|risk|poor") ~ "OS | Adverse",
    str_detect(s, "protect|favorable|good") ~ "OS | Protective",
    TRUE ~ "OS | Protective"
  )
}

dominance <- read_csv(dominance_path, show_col_types = FALSE)
alignment <- read_csv(alignment_path, show_col_types = FALSE)
depmap_all <- read_csv(depmap_path, show_col_types = FALSE)

gene_rf_summary <- if (file.exists(gene_rf_path)) {
  read_csv(gene_rf_path, show_col_types = FALSE) %>%
    transmute(
      gene_clean = gene,
      layer_short = case_when(
        str_detect(str_to_lower(Layer), "methyl") ~ "METH",
        str_detect(str_to_lower(Layer), "transcript") ~ "RNA",
        TRUE ~ NA_character_
      ),
      rf_importance = as.numeric(rf_importance)
    ) %>%
    filter(!is.na(layer_short), !is.na(gene_clean)) %>%
    group_by(gene_clean, layer_short) %>%
    summarise(rf_score = max(rf_importance, na.rm = TRUE), .groups = "drop") %>%
    mutate(rf_score = if_else(is.finite(rf_score), rf_score, 0)) %>%
    pivot_wider(names_from = layer_short, values_from = rf_score, names_prefix = "rf_")
} else {
  tibble(gene_clean = character(), rf_METH = double(), rf_RNA = double())
}
if (!"rf_METH" %in% names(gene_rf_summary)) gene_rf_summary$rf_METH <- numeric(nrow(gene_rf_summary))
if (!"rf_RNA" %in% names(gene_rf_summary)) gene_rf_summary$rf_RNA <- numeric(nrow(gene_rf_summary))
gene_rf_summary <- gene_rf_summary %>%
  mutate(
    rf_METH = coalesce(rf_METH, 0),
    rf_RNA = coalesce(rf_RNA, 0),
    rf_total = rf_METH + rf_RNA,
    rf_dom = case_when(
      rf_total <= 0 ~ "RF n/a",
      rf_METH >= rf_RNA * 1.25 ~ "RF M-dominant",
      rf_RNA >= rf_METH * 1.25 ~ "RF R-dominant",
      TRUE ~ "RF M/R-balanced"
    ),
    rf_layer_label = if_else(
      rf_total <= 0,
      "RF n/a",
      paste0(
        case_when(
          rf_dom == "RF M-dominant" ~ "METH-dom",
          rf_dom == "RF R-dominant" ~ "RNA-dom",
          rf_dom == "RF M/R-balanced" ~ "METH/RNA-bal",
          TRUE ~ "RF n/a"
        ),
        " | RF M", sprintf("%.1f", rf_METH),
        "/R", sprintf("%.1f", rf_RNA)
      )
    )
  )
write_csv(gene_rf_summary, file.path(table_dir, "meth_group_pathway_evidence_biomarker_right_sankey_rf_layer_summary.csv"))

immune_recomputed <- if (file.exists(expr_immune_path) && file.exists(immune_score_path)) {
  expr_immune <- read_csv(expr_immune_path, show_col_types = FALSE) %>%
    filter(modality == "RNA") %>%
    transmute(gene_clean = gene, Sample, rna_z = as.numeric(value_z))
  immune_score <- read_csv(immune_score_path, show_col_types = FALSE) %>%
    transmute(
      Sample,
      TIL_score = as.numeric(TIL_score),
      tmb_log1p = as.numeric(tmb_log1p)
    )
  tmb_outlier_samples <- immune_score %>%
    filter(!is.na(tmb_log1p), tmb_log1p <= 0) %>%
    distinct(Sample) %>%
    pull(Sample)
  calc_metric <- function(metric_name, metric_col, excluded_samples = character()) {
    immune_metric <- immune_score %>%
      filter(!Sample %in% excluded_samples) %>%
      transmute(Sample, metric_value = .data[[metric_col]])
    expr_immune %>%
      inner_join(immune_metric, by = "Sample") %>%
      group_by(gene_clean) %>%
      group_modify(~ spearman_summary(.x$rna_z, .x$metric_value)) %>%
      ungroup() %>%
      mutate(metric = metric_name)
  }
  bind_rows(
    calc_metric("TIL", "TIL_score"),
    calc_metric("TMB", "tmb_log1p", tmb_outlier_samples)
  ) %>%
    group_by(metric) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup()
} else {
  tibble(gene_clean = character(), metric = character(), n = integer(), rho = double(), p = double(), q = double())
}
write_csv(immune_recomputed, file.path(table_dir, "meth_group_pathway_evidence_biomarker_right_sankey_immune_til_tmb_recomputed.csv"))

immune_recomputed_wide <- immune_recomputed %>%
  filter(metric %in% c("TIL", "TMB")) %>%
  select(gene_clean, metric, n, rho, p, q) %>%
  pivot_wider(
    names_from = metric,
    values_from = c(n, rho, p, q),
    names_glue = "{metric}_{.value}"
  )

alignment_prepped <- alignment %>%
  transmute(
    gene_clean = gene,
    sample_spearman_rho,
    group_spearman_rho,
    inverse_consistency_class = coalesce(inverse_consistency_class, ""),
    rna_matches_meth_low_activity = truthy(rna_matches_meth_low_activity),
    rna_matches_raw_meth_high = truthy(rna_matches_raw_meth_high),
    sample_inverse_flag = truthy(sample_inverse_flag),
    group_inverse_flag = truthy(group_inverse_flag),
    rna_support_class = case_when(
      str_detect(inverse_consistency_class, regex("Strong", ignore_case = TRUE)) ~ "Meth-RNA inverse",
      rna_matches_meth_low_activity | sample_inverse_flag | group_inverse_flag ~ "Meth-RNA inverse",
      rna_matches_raw_meth_high ~ "Meth-RNA positive",
      TRUE ~ "Meth-RNA weak"
    )
  )

meth_markers <- dominance %>%
  filter(display_layer == "TSS methylation") %>%
  mutate(
    gene_clean = gene,
    EOBC_group = factor(activity_dominant_group, levels = group_order),
    Family6 = if_else(Family6 %in% names(family_cols), Family6, "Immune"),
    pathway_program = pathway_subclass(Family6, gene_clean, target_label),
    marker_label = gene_clean,
    group_weight = pmax(0.18, pmin(activity_contrast_z, 2.2)),
    neg_log10_fdr = pmin(coalesce(neg_log10_fdr, 0), 12),
    os_sig = truthy(os_sig),
    drug_sig = truthy(drug_sig),
    immune_sig = truthy(immune_sig)
  ) %>%
  left_join(alignment_prepped, by = "gene_clean") %>%
  left_join(
    gene_rf_summary %>%
      select(gene_clean, rf_METH, rf_RNA, rf_dom, rf_layer_label),
    by = "gene_clean"
  ) %>%
  mutate(
    rna_support_class = coalesce(rna_support_class, "Meth-RNA weak"),
    meth_rna_class = rna_support_class,
    rf_METH = coalesce(rf_METH, 0),
    rf_RNA = coalesce(rf_RNA, 0),
    rf_dom = coalesce(rf_dom, "RF n/a"),
    rf_layer_label = coalesce(rf_layer_label, "RF n/a"),
    biomarker_node = paste0(marker_label, "\n", rna_support_class, "\n", rf_layer_label)
  )

depmap_routes <- depmap_all %>%
  filter(str_to_upper(omics) == "METH") %>%
  filter(gene %in% meth_markers$gene_clean) %>%
  mutate(
    gene_clean = gene,
    drug_class = classify_drug(drug_clean),
    drug_direction = case_when(
      str_detect(str_to_lower(direction), "sensitive") ~ "Sensitive",
      str_detect(str_to_lower(direction), "resistant") ~ "Resistant",
      rho < 0 ~ "Sensitive",
      TRUE ~ "Resistant"
    ),
    evidence_node = paste0("Drug | ", drug_class, "\n", drug_direction),
    drug_label = paste0(drug_clean, " (", drug_direction, ")"),
    fdr = as.numeric(fdr),
    neglog10_fdr = -log10(pmax(fdr, 1e-300))
  ) %>%
  filter(!is.na(fdr), fdr <= 0.25, !drug_class %in% c("Other/exploratory", "Experimental/natural")) %>%
  group_by(gene_clean, evidence_node, drug_class, drug_direction) %>%
  summarise(
    n_drugs = n(),
    best_fdr = min(fdr, na.rm = TRUE),
    best_neglog10_fdr = max(neglog10_fdr, na.rm = TRUE),
    top_drugs = paste(head(unique(drug_label[order(fdr)]), 4), collapse = "; "),
    .groups = "drop"
  ) %>%
  left_join(
    meth_markers %>%
      select(
        gene_clean, EOBC_group, Family6, pathway_program, biomarker_node,
        meth_rna_class,
        rf_METH, rf_RNA, rf_dom, rf_layer_label,
        activity_contrast_z, group_weight, neg_log10_fdr
      ),
    by = "gene_clean"
  ) %>%
  mutate(
    route_type = "DepMap",
    route_weight = group_weight * pmin(1.35, 0.60 + 0.13 * n_drugs + 0.06 * pmin(best_neglog10_fdr, 3))
  )

os_routes <- meth_markers %>%
  filter(os_sig) %>%
  transmute(
    gene_clean,
    EOBC_group,
    Family6,
    pathway_program,
    evidence_node = evidence_from_os(os_status, os_activity_direction),
    biomarker_node,
    meth_rna_class,
    rf_METH,
    rf_RNA,
    rf_dom,
    rf_layer_label,
    route_type = "OS",
    route_weight = group_weight * 1.15,
    activity_contrast_z,
    neg_log10_fdr,
    n_drugs = NA_integer_,
    best_fdr = NA_real_,
    best_neglog10_fdr = NA_real_,
    top_drugs = NA_character_
  )

os_km_filter_path <- file.path(
  final_revision_dir,
  "individual_biomarker_evidence_20260522",
  "tables",
  "supplement_OS_KM_sankey_aligned_source_table.csv"
)
if (file.exists(os_km_filter_path)) {
  os_keep <- read_csv(os_km_filter_path, show_col_types = FALSE) %>%
    mutate(individual_logrank_p = as.numeric(individual_logrank_p)) %>%
    filter(is.finite(individual_logrank_p), individual_logrank_p < 0.05) %>%
    pull(gene) %>%
    unique()
  os_routes <- os_routes %>% filter(gene_clean %in% os_keep)
}

immune_routes <- meth_markers %>%
  filter(immune_sig) %>%
  left_join(immune_recomputed_wide, by = "gene_clean") %>%
  mutate(
    recomputed_immune_evidence = evidence_from_recomputed_til_tmb(TIL_rho, TMB_rho, TIL_q, TMB_q)
  ) %>%
  transmute(
    gene_clean,
    EOBC_group,
    Family6,
    pathway_program,
    evidence_node = recomputed_immune_evidence,
    biomarker_node,
    meth_rna_class,
    rf_METH,
    rf_RNA,
    rf_dom,
    rf_layer_label,
    TIL_n,
    TMB_n,
    TIL_rho,
    TMB_rho,
    TIL_q,
    TMB_q,
    route_type = "Immune",
    route_weight = group_weight * 1.00,
    activity_contrast_z,
    neg_log10_fdr,
    n_drugs = NA_integer_,
    best_fdr = NA_real_,
    best_neglog10_fdr = NA_real_,
    top_drugs = NA_character_
  ) %>%
  filter(!is.na(evidence_node))

depmap_routes_for_bind <- depmap_routes %>%
  transmute(
    gene_clean,
    EOBC_group,
    Family6,
    pathway_program,
    evidence_node,
    biomarker_node,
    meth_rna_class,
    rf_METH,
    rf_RNA,
    rf_dom,
    rf_layer_label,
    route_type,
    route_weight,
    activity_contrast_z,
    neg_log10_fdr,
    n_drugs,
    best_fdr,
    best_neglog10_fdr,
    top_drugs
  )

routes_full <- bind_rows(os_routes, immune_routes, depmap_routes_for_bind) %>%
  filter(!is.na(EOBC_group)) %>%
  mutate(
    EOBC_group = factor(as.character(EOBC_group), levels = group_order),
    Family6 = factor(Family6, levels = names(family_cols)),
    rf_dom = factor(rf_dom, levels = names(rf_dom_cols)),
    route_weight = pmax(route_weight, 0.08)
  )

biomarker_sets <- routes_full %>%
  group_by(gene_clean) %>%
  summarise(
    has_os = any(route_type == "OS"),
    has_depmap = any(route_type == "DepMap"),
    has_immune = any(route_type == "Immune"),
    n_domains = sum(c(has_os, has_depmap, has_immune)),
    biomarker_set = case_when(
      n_domains >= 3 ~ "Multi-domain",
      has_os & has_immune ~ "OS + Immune",
      has_os & has_depmap ~ "OS + DepMap",
      has_depmap & has_immune ~ "DepMap + Immune",
      has_os ~ "OS-only",
      has_depmap ~ "DepMap-only",
      has_immune ~ "Immune-only",
      TRUE ~ "Unassigned"
    ),
    .groups = "drop"
  )

routes_full <- routes_full %>%
  left_join(biomarker_sets %>% select(gene_clean, biomarker_set), by = "gene_clean") %>%
  mutate(biomarker_set = factor(biomarker_set, levels = names(biomarker_set_cols)))

priority_genes <- routes_full %>%
  group_by(gene_clean) %>%
  summarise(
    n_functional_domains = n_distinct(route_type),
    max_contrast = max(activity_contrast_z, na.rm = TRUE),
    has_os = any(route_type == "OS"),
    has_depmap = any(route_type == "DepMap"),
    has_immune = any(route_type == "Immune"),
    .groups = "drop"
  ) %>%
  filter(has_os | n_functional_domains >= 2 | max_contrast >= 1.35) %>%
  pull(gene_clean)

routes_priority <- routes_full %>%
  filter(gene_clean %in% priority_genes)

biomarker_summary <- routes_full %>%
  group_by(gene_clean, biomarker_node, EOBC_group, Family6, pathway_program, biomarker_set, meth_rna_class, rf_METH, rf_RNA, rf_dom, rf_layer_label) %>%
  summarise(
    n_routes = n(),
    has_os = any(route_type == "OS"),
    has_depmap = any(route_type == "DepMap"),
    has_immune = any(route_type == "Immune"),
    evidence_domains = paste(sort(unique(route_type)), collapse = "; "),
    depmap_evidence = paste(na.omit(unique(top_drugs)), collapse = " | "),
    max_activity_contrast_z = max(activity_contrast_z, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(EOBC_group, desc(has_os + has_depmap + has_immune), desc(max_activity_contrast_z))

write_csv(routes_full, file.path(table_dir, "meth_group_pathway_evidence_biomarker_right_sankey_routes_full.csv"))
write_csv(routes_priority, file.path(table_dir, "meth_group_pathway_evidence_biomarker_right_sankey_routes_priority.csv"))
write_csv(biomarker_summary, file.path(table_dir, "meth_group_pathway_evidence_biomarker_right_sankey_biomarker_summary.csv"))
write_csv(depmap_routes, file.path(table_dir, "meth_group_pathway_evidence_biomarker_right_sankey_depmap_drugclass_routes.csv"))

plot_sankey <- function(df, title, subtitle, width, height, outfile_stub, label_size = 2.25) {
  df_plot <- df %>%
    mutate(
      `EOBC methylation state` = unname(group_display[as.character(EOBC_group)]),
      `EOBC methylation state` = if_else(
        is.na(`EOBC methylation state`),
        as.character(EOBC_group),
        `EOBC methylation state`
      ),
      `Biological program` = pathway_program,
      `Functional evidence` = evidence_node,
      `Convergent biomarker` = biomarker_node
    )

  state_levels <- unname(group_display[group_order])
  program_levels <- df_plot %>%
    distinct(`Biological program`, Family6) %>%
    arrange(Family6, `Biological program`) %>%
    pull(`Biological program`)
  evidence_levels <- df_plot %>%
    distinct(`Functional evidence`) %>%
    mutate(
      evidence_order = case_when(
        str_detect(`Functional evidence`, "^Drug \\| ") ~ 1L,
        str_detect(`Functional evidence`, "^Immune \\| ") ~ 2L,
        str_detect(`Functional evidence`, "^OS \\| ") ~ 3L,
        TRUE ~ 4L
      )
    ) %>%
    arrange(evidence_order, `Functional evidence`) %>%
    pull(`Functional evidence`)
  terminal_levels <- df_plot %>%
    distinct(`Convergent biomarker`, biomarker_set, Family6, gene_clean) %>%
    arrange(biomarker_set, Family6, gene_clean) %>%
    pull(`Convergent biomarker`)

  df_plot <- df_plot %>%
    mutate(
      route_id = paste0("route_", row_number()),
      `EOBC methylation state` = factor(`EOBC methylation state`, levels = state_levels),
      `Biological program` = factor(`Biological program`, levels = program_levels),
      `Functional evidence` = factor(`Functional evidence`, levels = evidence_levels),
      `Convergent biomarker` = factor(`Convergent biomarker`, levels = terminal_levels),
      route_color = blend_route_color(EOBC_group, Family6, `Functional evidence`, biomarker_set, rf_dom),
      route_color_id = route_color
    )
  route_stage_cols <- route_stage_colors(
    df_plot$EOBC_group,
    df_plot$Family6,
    df_plot$`Functional evidence`,
    df_plot$biomarker_set,
    df_plot$rf_dom,
    df_plot$meth_rna_class
  )
  df_plot <- bind_cols(df_plot, route_stage_cols)

  axis_names <- c("EOBC methylation state", "Biological program", "Convergent biomarker", "Functional evidence")
  axis_fill_names <- c("stage1_color", "stage2_color", "stage3_color", "stage4_color")

  all_strata <- unique(c(
    as.character(df_plot$`EOBC methylation state`),
    as.character(df_plot$`Biological program`),
    as.character(df_plot$`Convergent biomarker`),
    as.character(df_plot$`Functional evidence`)
  ))
  stratum_cols <- setNames(rep("#F8FAFC", length(all_strata)), all_strata)

  group_cols_for_plot <- group_node_cols[names(group_node_cols) %in% names(stratum_cols)]
  stratum_cols[names(group_cols_for_plot)] <- scales::alpha(group_cols_for_plot, 0.90)

  program_map <- df_plot %>%
    distinct(`Biological program`, Family6) %>%
    filter(!is.na(Family6), `Biological program` %in% names(stratum_cols))
  stratum_cols[as.character(program_map$`Biological program`)] <- scales::alpha(family_cols[as.character(program_map$Family6)], 0.84)

  biomarker_map <- df_plot %>%
    group_by(`Convergent biomarker`) %>%
    arrange(desc(route_weight), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    select(`Convergent biomarker`, biomarker_set, Family6, EOBC_group, rf_dom, meth_rna_class) %>%
    filter(!is.na(biomarker_set), `Convergent biomarker` %in% names(stratum_cols))
  stratum_cols[as.character(biomarker_map$`Convergent biomarker`)] <- mapply(
    function(b, f, g, r, m) {
      boost_col(mix_cols(
        c(pick_col(meth_rna_cols, m), pick_col(rf_dom_cols, r), pick_col(biomarker_set_cols, b), pick_col(family_cols, f), pick_col(group_node_cols, g)),
        c(0.38, 0.24, 0.18, 0.13, 0.07)
      ), sat = 1.24, val = 1.03)
    },
    as.character(biomarker_map$biomarker_set),
    as.character(biomarker_map$Family6),
    as.character(biomarker_map$EOBC_group),
    as.character(biomarker_map$rf_dom),
    as.character(biomarker_map$meth_rna_class),
    USE.NAMES = FALSE
  )

  evidence_names <- unique(as.character(df_plot$`Functional evidence`))
  stratum_cols[evidence_names] <- evidence_color_for_node(evidence_names)

  df_lodes <- bind_rows(lapply(seq_along(axis_names), function(i) {
    df_plot %>%
      transmute(
        x = factor(axis_names[[i]], levels = axis_names),
        stratum = factor(as.character(.data[[axis_names[[i]]]]), levels = all_strata),
        alluvium = route_id,
        route_weight = route_weight,
        biomarker_set = biomarker_set,
        lode_fill = .data[[axis_fill_names[[i]]]]
      )
  }))

  lode_color_values <- setNames(unique(df_lodes$lode_fill), unique(df_lodes$lode_fill))
  route_color_values <- setNames(unique(df_plot$route_color_id), unique(df_plot$route_color_id))
  group_fill_cols <- group_node_cols[group_order]
  fill_cols <- c(
    lode_color_values,
    route_color_values,
    biomarker_set_cols,
    family_cols,
    group_fill_cols,
    stratum_cols[setdiff(
      names(stratum_cols),
      c(names(biomarker_set_cols), names(family_cols), names(group_fill_cols), names(route_color_values))
    )]
  )
  column_bands <- tibble(
    xmin = c(0.82, 1.82, 2.82, 3.82),
    xmax = c(1.18, 2.18, 3.18, 4.18),
    fill = c("#FFF7ED", "#F0F9FF", "#F8FAFC", "#F0FDFA")
  )

  p <- ggplot(
    df_lodes,
    aes(
      x = x,
      stratum = stratum,
      alluvium = alluvium,
      y = route_weight
    )
  ) +
    geom_rect(
      data = column_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = column_bands$fill,
      alpha = 0.34,
      color = NA
    ) +
    geom_flow(fill = "#0F172A", alpha = 0.055, width = 0.128, knot.pos = 0.50, color = NA, show.legend = FALSE) +
    geom_flow(
      aes(fill = lode_fill),
      aes.flow = "forward",
      alpha = 0.52,
      width = 0.102,
      knot.pos = 0.47,
      color = scales::alpha("white", 0.30),
      linewidth = 0.09,
      show.legend = FALSE
    ) +
    geom_flow(
      aes(fill = lode_fill),
      aes.flow = "backward",
      alpha = 0.48,
      width = 0.092,
      knot.pos = 0.53,
      color = NA,
      show.legend = FALSE
    ) +
    geom_flow(
      aes(fill = biomarker_set),
      alpha = 0.08,
      width = 0.026,
      knot.pos = 0.48,
      color = NA,
      show.legend = TRUE,
      aes.flow = "forward"
    ) +
    geom_stratum(
      fill = "#0F172A",
      alpha = 0.13,
      width = 0.145,
      color = NA,
      show.legend = FALSE
    ) +
    geom_stratum(
      aes(fill = after_stat(stratum)),
      width = 0.112,
      color = "#0F172A",
      linewidth = 0.62,
      show.legend = FALSE
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = label_size,
      lineheight = 0.86,
      color = "#111827",
      fontface = "bold"
    ) +
    scale_x_discrete(
        limits = axis_names,
        labels = c("EOBC state\n(raw TSS-METH)", "Biological\nprogram", "Terminal\nbiomarker", "Evidence\ndomain"),
        expand = c(0.038, 0.020)
      ) +
    scale_fill_manual(values = fill_cols, breaks = names(biomarker_set_cols), drop = FALSE, na.value = "#E5E7EB") +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Evidence-flow weight",
      fill = "Biomarker evidence set",
      caption = paste(
        "Ribbon hues are route-specific blends that transition across columns: EOBC-state at entry, EOBC+program in the middle, biomarker Meth-RNA/RF layer near terminal genes, and evidence-domain at the endpoint.",
        "Thin inner strands retain terminal evidence-set class. Labels show Meth-RNA class plus RF methylation/RNA scores.",
        "OS routes require KM log-rank p < 0.05. DepMap uses curated drug classes only. TMB-log1p zero outliers are excluded only for TMB correlations.",
        sep = "\n"
      )
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold", color = "#111827", margin = margin(t = 12)),
      axis.title.y = element_text(size = 13, face = "bold", color = "#111827", margin = margin(r = 10)),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 25.5, face = "bold", color = "#0F172A", margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12.8, color = "#64748B", margin = margin(b = 18)),
      plot.caption = element_text(size = 9.3, color = "#64748B", hjust = 0, margin = margin(t = 16)),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold", color = "#111827"),
      legend.text = element_text(size = 10.5, color = "#64748B"),
      plot.margin = margin(20, 30, 18, 30)
    ) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 0.90)))

  pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
  for (dest_dir in unique(c(plot_dir, final_fig_dir, final_revision_dir))) {
    ggsave(file.path(dest_dir, paste0(outfile_stub, ".png")), p, width = width, height = height, dpi = 320, bg = "white")
    ggsave(file.path(dest_dir, paste0(outfile_stub, ".pdf")), p, width = width, height = height, device = pdf_device, bg = "white")
  }
  p
}

plot_sankey(
  routes_full,
  title = "EOBC raw TSS-methylation biomarker-set evidence Sankey",
  subtitle = "Only harmonized functional evidence routes are shown: KM-significant OS, curated DepMap therapeutic classes, and continuous TIL/TMB associations.",
  width = 26.2,
  height = 15.2,
  outfile_stub = "Figure_10C_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_full_R_v12",
  label_size = 1.55
)

plot_sankey(
  routes_priority,
  title = "Priority EOBC raw TSS-methylation biomarker-set evidence Sankey",
  subtitle = "Focused view after harmonizing OS KM evidence and removing weak drug classes or unsupported immune categories.",
  width = 24,
  height = 12.8,
  outfile_stub = "Figure_10D_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_priority_R_v12",
  label_size = 1.75
)

message("Saved Sankey plots to: ", plot_dir)
message("Saved Sankey final copies to: ", final_fig_dir)
message("Saved Sankey revision copies to: ", final_revision_dir)
message("Saved Sankey tables to: ", table_dir)
