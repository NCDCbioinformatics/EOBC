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

group_short <- c(
  "G1 | H" = "G1",
  "G2 | I" = "G2",
  "G3 | L-like" = "G3",
  "G4 | L" = "G4"
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
  "Kinase signaling" = "#9B6AE8",
  "Developmental signaling" = "#B8A6C9"
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
  "Immune | posi" = "#17A2A4",
  "Immune | nega" = "#72A950"
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
    str_detect(x, "^DepMap \\| Sensitive") ~ "#2563EB",
    str_detect(x, "^DepMap \\| Resistant") ~ "#D95F02",
    str_detect(x, "^Drug \\| .*Sensitive") ~ "#2563EB",
    str_detect(x, "^Drug \\| .*Resistant") ~ "#D95F02",
    str_detect(x, "^Immune \\| posi") ~ "#17A2A4",
    str_detect(x, "^Immune \\| nega") ~ "#72A950",
    str_detect(x, regex("^Immune \\| .*positive", ignore_case = TRUE)) ~ "#17A2A4",
    str_detect(x, regex("^Immune \\| .*negative", ignore_case = TRUE)) ~ "#72A950",
    str_detect(x, regex("^Immune \\| .*weak", ignore_case = TRUE)) ~ "#94A3B8",
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
        c(0.64, 0.15, 0.12, 0.05, 0.04)
      ), sat = 1.18, val = 1.03)
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
        c(0.48, 0.32, 0.20)
      ), sat = 1.06, val = 1.02)

      c(
        stage1_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col, evidence_col),
          c(0.92, 0.05, 0.02, 0.01)
        ), sat = 1.10, val = 1.02),
        stage2_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col, evidence_col),
          c(0.74, 0.19, 0.05, 0.02)
        ), sat = 1.10, val = 1.02),
        stage3_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col, evidence_col),
          c(0.66, 0.09, 0.19, 0.06)
        ), sat = 1.10, val = 1.02),
        stage4_color = boost_col(mix_cols(
          c(group_col, family_col, biomarker_col, evidence_col),
          c(0.62, 0.06, 0.08, 0.24)
        ), sat = 1.12, val = 1.02)
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

first_nonempty <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) NA_character_ else x[[1]]
}

format_p_compact <- function(x) {
  case_when(
    !is.finite(x) ~ "p NA",
    x < 0.001 ~ "p <0.001",
    TRUE ~ paste0("p ", formatC(x, format = "f", digits = 3))
  )
}

depmap_class_label <- function(x) {
  x <- coalesce(as.character(x), "Drug class n/a")
  x %>%
    str_replace_all("Mitotic/Cell-cycle", "Mitotic/cell-cycle") %>%
    str_replace_all("DDR/Checkpoint", "DDR/checkpoint") %>%
    str_replace_all("Other/Experimental", "Other/experimental")
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
    x == "Developmental\nSignaling" ~ "Developmental signaling\nDevelopmental signaling",
    x == "Hormone\nSignaling" ~ "Developmental signaling\nDevelopmental signaling",
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
    str_replace("^Immune \\| ", "TIL/TMB immune axis\n")
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
    str_replace_all("\\s*\\[[MR]\\]", "") %>%
    str_replace_all("RNA inv strong", "Meth-RNA inverse") %>%
    str_replace_all("RNA inv\\b", "Meth-RNA inverse") %>%
    str_replace_all("RNA ctx\\b", "Meth-RNA positive") %>%
    str_replace_all("RNA weak\\b", "Meth-RNA weak")
}

compact_biomarker_label <- function(x) {
  labels <- pretty_biomarker_label(x)
  vapply(
    strsplit(labels, "\n", fixed = TRUE),
    function(parts) paste(head(parts, 2), collapse = "\n"),
    character(1)
  )
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
    family == "Developmental signaling" & str_detect(target_label, regex("WNT", ignore_case = TRUE)) ~ "WNT Signaling",
    family == "Developmental signaling" ~ "Developmental\nSignaling",
    family == "Hormone signaling" & str_detect(target_label, regex("WNT", ignore_case = TRUE)) ~ "WNT Signaling",
    family == "Hormone signaling" ~ "Developmental\nSignaling",
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
    str_detect(q, "positive") ~ "Immune | posi",
    str_detect(q, "negative") ~ "Immune | nega",
    TRUE ~ NA_character_
  )
}

evidence_from_til_tmb_pair <- function(til_rho, tmb_rho, weak_cutoff = 0.20) {
  case_when(
    is.finite(til_rho) & is.finite(tmb_rho) &
      til_rho > 0 & tmb_rho > 0 &
      pmin(abs(til_rho), abs(tmb_rho)) >= weak_cutoff ~ "Immune | posi",
    is.finite(til_rho) & is.finite(tmb_rho) &
      til_rho < 0 & tmb_rho < 0 &
      pmin(abs(til_rho), abs(tmb_rho)) >= weak_cutoff ~ "Immune | nega",
    TRUE ~ NA_character_
  )
}

evidence_from_recomputed_til_tmb <- function(til_rho, tmb_rho, til_q, tmb_q, weak_cutoff = 0.15, q_cutoff = 0.25) {
  til_sig <- is.finite(til_rho) & is.finite(til_q) & til_q <= q_cutoff
  tmb_sig <- is.finite(tmb_rho) & is.finite(tmb_q) & tmb_q <= q_cutoff
  til_abs <- if_else(is.finite(til_rho), abs(til_rho), -Inf)
  tmb_abs <- if_else(is.finite(tmb_rho), abs(tmb_rho), -Inf)
  best_rho <- case_when(
    til_sig & tmb_sig & til_abs >= tmb_abs ~ til_rho,
    til_sig & tmb_sig ~ tmb_rho,
    til_sig ~ til_rho,
    tmb_sig ~ tmb_rho,
    TRUE ~ NA_real_
  )
  case_when(
    (til_sig | tmb_sig) & is.finite(best_rho) & best_rho >= 0 ~ "Immune | posi",
    (til_sig | tmb_sig) & is.finite(best_rho) & best_rho < 0 ~ "Immune | nega",
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

feature_group_summary <- dominance %>%
  mutate(
    gene_clean = gene,
    layer_short = case_when(
      display_layer == "TSS methylation" ~ "METH",
      display_layer == "RNA expression" ~ "RNA",
      TRUE ~ NA_character_
    ),
    feature_group = unname(group_short[as.character(activity_dominant_group)]),
    feature_fdr = suppressWarnings(as.numeric(kw_fdr)),
    feature_sig = is.finite(feature_fdr) & feature_fdr <= 0.05
  ) %>%
  filter(layer_short %in% c("METH", "RNA")) %>%
  group_by(gene_clean, layer_short) %>%
  arrange(feature_fdr, .by_group = TRUE) %>%
  summarise(
    feature_group = first(feature_group),
    feature_fdr = first(feature_fdr),
    feature_sig = any(feature_sig, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = layer_short,
    values_from = c(feature_group, feature_fdr, feature_sig),
    names_glue = "{layer_short}_{.value}"
  ) %>%
  mutate(
    METH_feature_sig = coalesce(METH_feature_sig, FALSE),
    RNA_feature_sig = coalesce(RNA_feature_sig, FALSE),
    meth_feature_label = if_else(
      METH_feature_sig & !is.na(METH_feature_group),
      paste0("METH ", METH_feature_group),
      "METH n.s."
    ),
    rna_feature_label = if_else(
      RNA_feature_sig & !is.na(RNA_feature_group),
      paste0("RNA ", RNA_feature_group),
      "RNA n.s."
    ),
    feature_layer_label = paste(meth_feature_label, rna_feature_label, sep = " | ")
  )
write_csv(feature_group_summary, file.path(table_dir, "meth_group_pathway_evidence_biomarker_right_sankey_feature_group_summary.csv"))

feature_markers_raw <- dominance %>%
  filter(display_layer %in% c("TSS methylation", "RNA expression")) %>%
  mutate(
    gene_clean = gene,
    EOBC_group = factor(activity_dominant_group, levels = group_order),
    source_layer = case_when(
      display_layer == "TSS methylation" ~ "METH",
      display_layer == "RNA expression" ~ "RNA",
      TRUE ~ NA_character_
    ),
    Family6 = recode(as.character(Family6), "Hormone signaling" = "Developmental signaling"),
    Family6 = if_else(Family6 %in% names(family_cols), Family6, "Immune"),
    pathway_program = pathway_subclass(Family6, gene_clean, target_label),
    marker_label = str_remove(gene_clean, "\\s*\\[[MR]\\]$"),
    kw_fdr = suppressWarnings(as.numeric(kw_fdr)),
    feature_sig = is.finite(kw_fdr) & kw_fdr <= 0.05,
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
  left_join(
    feature_group_summary %>%
      select(gene_clean, feature_layer_label, METH_feature_group, RNA_feature_group, METH_feature_sig, RNA_feature_sig),
    by = "gene_clean"
  ) %>%
  mutate(
    rna_support_class = coalesce(rna_support_class, "Meth-RNA weak"),
    meth_rna_class = rna_support_class,
    rf_METH = coalesce(rf_METH, 0),
    rf_RNA = coalesce(rf_RNA, 0),
    rf_dom = coalesce(rf_dom, "RF n/a"),
    rf_layer_label = coalesce(rf_layer_label, "RF n/a"),
    feature_layer_label = coalesce(feature_layer_label, "METH n.s. | RNA n.s."),
    biomarker_node = paste0(marker_label, "\n", feature_layer_label, "\n", rna_support_class, "\n", rf_layer_label)
  )

feature_marker_layers <- feature_markers_raw %>%
  filter(feature_sig | os_sig | drug_sig | immune_sig) %>%
  distinct(gene_clean, EOBC_group, source_layer)

feature_markers <- feature_markers_raw %>%
  filter(feature_sig | os_sig | drug_sig | immune_sig) %>%
  arrange(gene_clean, EOBC_group, desc(group_weight), desc(neg_log10_fdr)) %>%
  group_by(gene_clean, EOBC_group) %>%
  summarise(
    source_layer_label = paste(sort(unique(source_layer)), collapse = "+"),
    Family6 = first(Family6),
    pathway_program = first(pathway_program),
    biomarker_node = first(biomarker_node),
    feature_layer_label = first(feature_layer_label),
    meth_rna_class = first(meth_rna_class),
    rf_METH = first(rf_METH),
    rf_RNA = first(rf_RNA),
    rf_dom = first(rf_dom),
    rf_layer_label = first(rf_layer_label),
    group_weight = max(group_weight, na.rm = TRUE),
    activity_contrast_z = max(activity_contrast_z, na.rm = TRUE),
    neg_log10_fdr = max(neg_log10_fdr, na.rm = TRUE),
    os_sig = any(os_sig, na.rm = TRUE),
    drug_sig = any(drug_sig, na.rm = TRUE),
    immune_sig = any(immune_sig, na.rm = TRUE),
    os_status = first_nonempty(os_status[os_sig]),
    os_activity_direction = first_nonempty(os_activity_direction[os_sig]),
    immune_activity_quadrant = first_nonempty(immune_activity_quadrant[immune_sig]),
    immune_status = first_nonempty(immune_status[immune_sig]),
    .groups = "drop"
  ) %>%
  mutate(
    group_weight = if_else(is.finite(group_weight), group_weight, 0.18),
    activity_contrast_z = if_else(is.finite(activity_contrast_z), activity_contrast_z, 0),
    neg_log10_fdr = if_else(is.finite(neg_log10_fdr), neg_log10_fdr, 0)
  )
write_csv(feature_markers, file.path(table_dir, "meth_rna_union_group_marker_features_for_sankey.csv"))

depmap_routes <- depmap_all %>%
  mutate(
    gene_clean = gene,
    source_layer = str_to_upper(omics),
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
  filter(gene_clean %in% feature_markers$gene_clean) %>%
  filter(!is.na(fdr), fdr <= 0.25, !drug_class %in% c("Other/exploratory", "Experimental/natural")) %>%
  group_by(gene_clean, evidence_node, drug_class, drug_direction) %>%
  summarise(
    n_drugs = n(),
    evidence_layers = paste(sort(unique(source_layer)), collapse = "+"),
    best_fdr = min(fdr, na.rm = TRUE),
    best_neglog10_fdr = max(neglog10_fdr, na.rm = TRUE),
    top_drugs = paste(head(unique(drug_label[order(fdr)]), 4), collapse = "; "),
    .groups = "drop"
  ) %>%
  left_join(
    feature_markers %>%
      select(
        gene_clean, EOBC_group, Family6, pathway_program, biomarker_node,
        source_layer_label,
        feature_layer_label,
        meth_rna_class,
        rf_METH, rf_RNA, rf_dom, rf_layer_label,
        activity_contrast_z, group_weight, neg_log10_fdr
      ),
    by = "gene_clean",
    relationship = "many-to-many"
  ) %>%
  mutate(
    route_type = "DepMap",
    route_weight = group_weight * pmin(1.35, 0.60 + 0.13 * n_drugs + 0.06 * pmin(best_neglog10_fdr, 3))
  )

os_gene_evidence <- feature_markers %>%
  filter(os_sig) %>%
  group_by(gene_clean) %>%
  summarise(
    evidence_node = evidence_from_os(first_nonempty(os_status), first_nonempty(os_activity_direction)),
    .groups = "drop"
  )

os_routes <- feature_markers %>%
  inner_join(os_gene_evidence, by = "gene_clean") %>%
  transmute(
    gene_clean,
    EOBC_group,
    Family6,
    pathway_program,
    evidence_node,
    biomarker_node,
    source_layer_label,
    feature_layer_label,
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

immune_routes <- feature_markers %>%
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
    source_layer_label,
    feature_layer_label,
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
    source_layer_label,
    feature_layer_label,
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
  group_by(gene_clean, biomarker_node, source_layer_label, feature_layer_label, EOBC_group, Family6, pathway_program, biomarker_set, meth_rna_class, rf_METH, rf_RNA, rf_dom, rf_layer_label) %>%
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

domain_coord_dir <- file.path(base_dir, "group_biomarker_contextual_domain_raw_meth_r_v1", "tables")
os_story_summary <- read_csv(file.path(domain_coord_dir, "OS_significant_biomarker_summary_no_group.csv"), show_col_types = FALSE)
depmap_story_summary <- read_csv(file.path(domain_coord_dir, "DepMap_significant_biomarker_summary_no_group.csv"), show_col_types = FALSE)
immune_story_summary <- read_csv(file.path(domain_coord_dir, "Immune_significant_biomarker_summary_no_group.csv"), show_col_types = FALSE)

story_source_routes <- bind_rows(
  os_story_summary %>%
    filter(domain_status %in% c("Protective", "Adverse")) %>%
    mutate(
      Family6_source = recode(as.character(Family6), "Hormone signaling" = "Developmental signaling"),
      Family6_source = if_else(Family6_source %in% names(family_cols), Family6_source, "Immune")
    ) %>%
    transmute(
      gene_clean = gene,
      EOBC_group = activity_dominant_group,
      source_layer = modality,
      Family6_source,
      pathway_program_source = pathway_subclass(Family6_source, gene_clean, target_label),
      route_type = "OS",
      evidence_summary = paste0("OS | ", domain_status),
      evidence_display = evidence_summary,
      route_weight = pmin(0.86, 0.32 + 0.065 * pmin(domain_y, 5)),
      evidence_details = paste(os_best_endpoint, domain_status, format_p_compact(domain_p)),
      max_activity_contrast_z = activity_contrast_z
    ),
  depmap_story_summary %>%
    filter(domain_status %in% c("Sensitive", "Resistant")) %>%
    mutate(
      Family6_source = recode(as.character(Family6), "Hormone signaling" = "Developmental signaling"),
      Family6_source = if_else(Family6_source %in% names(family_cols), Family6_source, "Immune")
    ) %>%
    transmute(
      gene_clean = gene,
      EOBC_group = activity_dominant_group,
      source_layer = modality,
      Family6_source,
      pathway_program_source = pathway_subclass(Family6_source, gene_clean, target_label),
      route_type = "DepMap",
      evidence_summary = paste0("DepMap | ", domain_status),
      evidence_display = paste0("DepMap | ", depmap_class_label(drug_class_clean), " ", domain_status),
      route_weight = pmin(0.82, 0.30 + 0.060 * pmin(domain_y, 5)),
      evidence_details = paste(drug_clean, domain_status, format_p_compact(domain_p)),
      max_activity_contrast_z = activity_contrast_z
    ),
  immune_story_summary %>%
    filter(domain_status %in% c("posi", "nega")) %>%
    mutate(
      Family6_source = recode(as.character(Family6), "Hormone signaling" = "Developmental signaling"),
      Family6_source = if_else(Family6_source %in% names(family_cols), Family6_source, "Immune")
    ) %>%
    transmute(
      gene_clean = gene,
      EOBC_group = activity_dominant_group,
      source_layer = modality,
      Family6_source,
      pathway_program_source = pathway_subclass(Family6_source, gene_clean, target_label),
      route_type = "Immune",
      evidence_summary = case_when(
        domain_status == "posi" ~ "Immune | TIL/TMB positive",
        domain_status == "nega" ~ "Immune | TIL/TMB negative",
        TRUE ~ "Immune | TIL/TMB weak"
      ),
      evidence_display = evidence_summary,
      route_weight = pmin(0.82, 0.30 + 0.060 * pmin(domain_y, 5)),
      evidence_details = paste0(
        "TIL rho ", sprintf("%.2f", til_activity_rho),
        " / TMB rho ", sprintf("%.2f", tmb_activity_rho),
        " / ", format_p_compact(domain_p)
      ),
      max_activity_contrast_z = activity_contrast_z
    )
) %>%
  filter(EOBC_group %in% group_order, !is.na(evidence_summary))

routes_story <- story_source_routes %>%
  left_join(
    feature_markers %>%
      select(
        gene_clean, EOBC_group, Family6, pathway_program, biomarker_node,
        source_layer_label, feature_layer_label, meth_rna_class,
        rf_METH, rf_RNA, rf_dom, rf_layer_label
      ),
    by = c("gene_clean", "EOBC_group")
  ) %>%
  left_join(
    feature_group_summary %>%
      distinct(gene_clean, .keep_all = TRUE) %>%
      select(gene_clean, fallback_feature_layer_label = feature_layer_label),
    by = "gene_clean"
  ) %>%
  left_join(
    gene_rf_summary %>%
      distinct(gene_clean, .keep_all = TRUE) %>%
      select(
        gene_clean,
        fallback_rf_METH = rf_METH,
        fallback_rf_RNA = rf_RNA,
        fallback_rf_dom = rf_dom,
        fallback_rf_layer_label = rf_layer_label
      ),
    by = "gene_clean"
  ) %>%
  left_join(
    alignment_prepped %>%
      distinct(gene_clean, .keep_all = TRUE) %>%
      select(gene_clean, fallback_meth_rna_class = rna_support_class),
    by = "gene_clean"
  ) %>%
  mutate(
    Family6 = coalesce(as.character(Family6), Family6_source),
    pathway_program = coalesce(pathway_program, pathway_program_source),
    source_layer_label = coalesce(source_layer_label, source_layer),
    feature_layer_label = coalesce(
      feature_layer_label,
      fallback_feature_layer_label,
      paste0(source_layer, " ", unname(group_short[EOBC_group]))
    ),
    meth_rna_class = coalesce(meth_rna_class, fallback_meth_rna_class, "Meth-RNA weak"),
    rf_METH = coalesce(rf_METH, fallback_rf_METH, 0),
    rf_RNA = coalesce(rf_RNA, fallback_rf_RNA, 0),
    rf_dom = coalesce(rf_dom, fallback_rf_dom, "RF n/a"),
    rf_layer_label = coalesce(rf_layer_label, fallback_rf_layer_label, "RF n/a"),
    source_group_label = paste0(source_layer, " ", unname(group_short[EOBC_group])),
    story_feature_label = if_else(
      str_detect(feature_layer_label, fixed(source_group_label)),
      feature_layer_label,
      paste0(source_group_label, " evidence")
    ),
    biomarker_node = paste0(
      gene_clean,
      "\n", story_feature_label,
      "\n", feature_layer_label,
      "\n", meth_rna_class,
      "\n", rf_layer_label
    )
  ) %>%
  group_by(
    gene_clean, EOBC_group, Family6, pathway_program, biomarker_node,
    source_layer_label, feature_layer_label, meth_rna_class,
    rf_METH, rf_RNA, rf_dom, rf_layer_label, route_type, evidence_summary, evidence_display
  ) %>%
  summarise(
    route_weight = max(route_weight, na.rm = TRUE),
    n_underlying_routes = n(),
    evidence_details = paste(unique(na.omit(evidence_details)), collapse = " | "),
    max_activity_contrast_z = max(max_activity_contrast_z, na.rm = TRUE),
    .groups = "drop"
  )

story_biomarker_sets <- routes_story %>%
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

routes_story <- routes_story %>%
  left_join(story_biomarker_sets %>% select(gene_clean, biomarker_set), by = "gene_clean") %>%
  mutate(
    route_weight = if_else(is.finite(route_weight), route_weight, 0.32),
    max_activity_contrast_z = if_else(is.finite(max_activity_contrast_z), max_activity_contrast_z, 0),
    EOBC_group = factor(as.character(EOBC_group), levels = group_order),
    Family6 = factor(Family6, levels = names(family_cols)),
    biomarker_set = factor(biomarker_set, levels = names(biomarker_set_cols)),
    rf_dom = factor(rf_dom, levels = names(rf_dom_cols))
  )

priority_story_genes <- story_biomarker_sets %>%
  filter(has_os | n_domains >= 3) %>%
  pull(gene_clean)

routes_story_priority <- routes_story %>%
  filter(gene_clean %in% priority_story_genes)

write_csv(routes_story, file.path(table_dir, "meth_rna_group_biomarker_domain_sankey_story_routes_full.csv"))
write_csv(routes_story_priority, file.path(table_dir, "meth_rna_group_biomarker_domain_sankey_story_routes_priority.csv"))

plot_sankey <- function(df, title, subtitle, width, height, outfile_stub, label_size = 2.05) {
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
      `Convergent biomarker` = pretty_biomarker_label(biomarker_node)
    )

  state_levels <- rev(unname(group_display[group_order]))
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

  axis_names <- c("EOBC methylation state", "Biological program", "Functional evidence", "Convergent biomarker")
  axis_fill_names <- c("stage1_color", "stage2_color", "stage4_color", "stage3_color")

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
        c(0.18, 0.14, 0.08, 0.10, 0.50)
      ), sat = 1.14, val = 1.02)
    },
    as.character(biomarker_map$biomarker_set),
    as.character(biomarker_map$Family6),
    as.character(biomarker_map$EOBC_group),
    as.character(biomarker_map$rf_dom),
    as.character(biomarker_map$meth_rna_class),
    USE.NAMES = FALSE
  )

  evidence_names <- unique(as.character(df_plot$`Functional evidence`))
  evidence_map <- df_plot %>%
    group_by(`Functional evidence`) %>%
    summarise(
      group_mix = mix_cols(vapply(EOBC_group, function(g) pick_col(group_node_cols, g), character(1)), route_weight),
      evidence_tone = evidence_color_for_node(first(as.character(`Functional evidence`))),
      .groups = "drop"
    ) %>%
    mutate(
      evidence_fill = mapply(
        function(g_col, e_col) {
          boost_col(mix_cols(c(g_col, e_col), c(0.72, 0.28)), sat = 1.10, val = 1.02)
        },
        group_mix,
        evidence_tone,
        USE.NAMES = FALSE
      )
    )
  stratum_cols[as.character(evidence_map$`Functional evidence`)] <- evidence_map$evidence_fill

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
    fill = c("#FFF7ED", "#F0F9FF", "#F0FDFA", "#F8FAFC")
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
      alpha = 0.66,
      width = 0.102,
      knot.pos = 0.47,
      color = scales::alpha("white", 0.36),
      linewidth = 0.09,
      show.legend = FALSE
    ) +
    geom_flow(
      aes(fill = lode_fill),
      aes.flow = "backward",
      alpha = 0.38,
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
      lineheight = 0.82,
      color = "#111827",
      fontface = "bold"
    ) +
    scale_x_discrete(
        limits = axis_names,
        labels = c("EOBC state\n(METH/RNA markers)", "Biological\nprogram", "Evidence\ndomain", "Terminal\nbiomarker"),
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
        "Ribbon hues keep EOBC-state identity across METH/RNA union marker routes, then receive subtle biological-program, biomarker Meth-RNA/RF, and evidence-domain tints.",
        "Terminal biomarker labels show significant METH/RNA EOBC feature groups, Meth-RNA class, and RF methylation/RNA scores.",
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

plot_story_sankey <- function(df, title, subtitle, width, height, outfile_stub, label_size = 2.05) {
  if (!"evidence_display" %in% names(df)) {
    df <- df %>% mutate(evidence_display = evidence_summary)
  }

  df_plot <- df %>%
    mutate(
      `EOBC state` = unname(group_display[as.character(EOBC_group)]),
      `EOBC state` = if_else(is.na(`EOBC state`), as.character(EOBC_group), `EOBC state`),
      `Terminal biomarker` = compact_biomarker_label(biomarker_node),
      evidence_display = coalesce(evidence_display, evidence_summary),
      evidence_label_short = evidence_display %>%
        str_remove("^OS \\|\\s*") %>%
        str_remove("^DepMap \\|\\s*") %>%
        str_remove("^Immune \\|\\s*"),
      `Evidence domain` = paste0(unname(group_short[as.character(EOBC_group)]), " | ", route_type, "\n", evidence_label_short)
    )

  state_levels <- unname(group_display[group_order])
  biomarker_levels <- df_plot %>%
    distinct(`Terminal biomarker`, gene_clean, EOBC_group, Family6, biomarker_set) %>%
    arrange(EOBC_group, biomarker_set, Family6, gene_clean) %>%
    pull(`Terminal biomarker`) %>%
    unique()
  evidence_levels <- df_plot %>%
    distinct(`Evidence domain`, EOBC_group, route_type, evidence_summary, evidence_display) %>%
    mutate(
      group_rank = match(as.character(EOBC_group), group_order),
      route_rank = case_when(
        route_type == "OS" ~ 1L,
        route_type == "Immune" ~ 2L,
        route_type == "DepMap" ~ 3L,
        TRUE ~ 9L
      ),
      direction_rank = case_when(
        str_detect(evidence_summary, "Protective|positive|Sensitive") ~ 1L,
        str_detect(evidence_summary, "Adverse|negative|Resistant") ~ 2L,
        TRUE ~ 3L
      )
    ) %>%
    arrange(group_rank, route_rank, direction_rank, evidence_display, `Evidence domain`) %>%
    pull(`Evidence domain`)

  df_plot <- df_plot %>%
    mutate(
      route_id = paste0("story_route_", row_number()),
      `EOBC state` = factor(`EOBC state`, levels = state_levels),
      `Terminal biomarker` = factor(`Terminal biomarker`, levels = biomarker_levels),
      `Evidence domain` = factor(`Evidence domain`, levels = evidence_levels)
    )

  stage_cols <- mapply(
    function(g, f, e, b, r, m) {
      group_col <- pick_col(group_node_cols, g)
      family_col <- pick_col(family_cols, f)
      evidence_col <- evidence_color_for_node(e)
      biomarker_col <- boost_col(mix_cols(
        c(group_col, family_col, pick_col(meth_rna_cols, m), pick_col(rf_dom_cols, r), pick_col(biomarker_set_cols, b)),
        c(0.72, 0.12, 0.06, 0.05, 0.05)
      ), sat = 1.08, val = 1.02)

      c(
        stage1_color = boost_col(mix_cols(c(group_col, family_col), c(0.96, 0.04)), sat = 1.10, val = 1.02),
        stage2_color = boost_col(mix_cols(c(group_col, biomarker_col, family_col), c(0.60, 0.30, 0.10)), sat = 1.08, val = 1.02),
        stage3_color = boost_col(mix_cols(c(group_col, evidence_col, family_col, biomarker_col), c(0.62, 0.22, 0.08, 0.08)), sat = 1.10, val = 1.02)
      )
    },
    as.character(df_plot$EOBC_group),
    as.character(df_plot$Family6),
    as.character(df_plot$evidence_summary),
    as.character(df_plot$biomarker_set),
    as.character(df_plot$rf_dom),
    as.character(df_plot$meth_rna_class),
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  )
  df_plot <- bind_cols(df_plot, tibble::as_tibble(t(stage_cols), .name_repair = "minimal"))

  axis_names <- c("EOBC state", "Terminal biomarker", "Evidence domain")
  axis_fill_names <- c("stage1_color", "stage2_color", "stage3_color")
  all_strata <- unique(c(
    as.character(df_plot$`EOBC state`),
    as.character(df_plot$`Terminal biomarker`),
    as.character(df_plot$`Evidence domain`)
  ))

  stratum_cols <- setNames(rep("#F8FAFC", length(all_strata)), all_strata)
  group_cols_for_plot <- group_node_cols[names(group_node_cols) %in% names(stratum_cols)]
  stratum_cols[names(group_cols_for_plot)] <- scales::alpha(group_cols_for_plot, 0.92)

  biomarker_map <- df_plot %>%
    group_by(`Terminal biomarker`) %>%
    summarise(
      group_mix = mix_cols(vapply(EOBC_group, function(g) pick_col(group_node_cols, g), character(1)), route_weight),
      family_mix = mix_cols(vapply(Family6, function(f) pick_col(family_cols, f), character(1)), route_weight),
      set_mix = mix_cols(vapply(biomarker_set, function(b) pick_col(biomarker_set_cols, b), character(1)), route_weight),
      meth_mix = mix_cols(vapply(meth_rna_class, function(m) pick_col(meth_rna_cols, m), character(1)), route_weight),
      .groups = "drop"
    ) %>%
    mutate(
      biomarker_fill = mapply(
        function(g_col, f_col, s_col, m_col) {
          boost_col(mix_cols(c(g_col, f_col, s_col, m_col), c(0.74, 0.12, 0.06, 0.08)), sat = 1.08, val = 1.02)
        },
        group_mix, family_mix, set_mix, meth_mix,
        USE.NAMES = FALSE
      )
    )
  stratum_cols[as.character(biomarker_map$`Terminal biomarker`)] <- biomarker_map$biomarker_fill

  evidence_map <- df_plot %>%
    group_by(`Evidence domain`) %>%
    summarise(
      group_mix = mix_cols(vapply(EOBC_group, function(g) pick_col(group_node_cols, g), character(1)), route_weight),
      evidence_tone = evidence_color_for_node(first(as.character(evidence_summary))),
      .groups = "drop"
    ) %>%
    mutate(
      evidence_fill = mapply(
        function(g_col, e_col) boost_col(mix_cols(c(g_col, e_col), c(0.72, 0.28)), sat = 1.08, val = 1.02),
        group_mix, evidence_tone,
        USE.NAMES = FALSE
      )
    )
  stratum_cols[as.character(evidence_map$`Evidence domain`)] <- evidence_map$evidence_fill

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
  group_fill_cols <- group_node_cols[group_order]
  fill_cols <- c(
    lode_color_values,
    biomarker_set_cols,
    group_fill_cols,
    stratum_cols[setdiff(names(stratum_cols), c(names(biomarker_set_cols), names(group_fill_cols)))]
  )
  column_bands <- tibble(
    xmin = c(0.82, 1.82, 2.82),
    xmax = c(1.18, 2.18, 3.18),
    fill = c("#FFF7ED", "#F0F9FF", "#F0FDFA")
  )

  p <- ggplot(
    df_lodes,
    aes(x = x, stratum = stratum, alluvium = alluvium, y = route_weight)
  ) +
    geom_rect(
      data = column_bands,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = column_bands$fill,
      alpha = 0.38,
      color = NA
    ) +
    geom_flow(fill = "#0F172A", alpha = 0.052, width = 0.14, knot.pos = 0.50, color = NA, show.legend = FALSE) +
    geom_flow(
      aes(fill = lode_fill),
      aes.flow = "forward",
      alpha = 0.76,
      width = 0.108,
      knot.pos = 0.48,
      color = scales::alpha("white", 0.38),
      linewidth = 0.10,
      show.legend = FALSE
    ) +
    geom_flow(
      aes(fill = biomarker_set),
      aes.flow = "forward",
      alpha = 0.045,
      width = 0.028,
      knot.pos = 0.49,
      color = NA,
      show.legend = TRUE
    ) +
    geom_stratum(fill = "#0F172A", alpha = 0.12, width = 0.16, color = NA, show.legend = FALSE) +
    geom_stratum(
      aes(fill = after_stat(stratum)),
      width = 0.122,
      color = "#0F172A",
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = label_size,
      lineheight = 0.82,
      color = "#111827",
      fontface = "bold"
    ) +
    scale_x_discrete(
      limits = axis_names,
      labels = c("EOBC state\n(METH/RNA markers)", "Evidence-linked\nbiomarker", "Functional\nevidence"),
      expand = c(0.06, 0.028)
    ) +
    scale_fill_manual(values = fill_cols, breaks = names(biomarker_set_cols), drop = FALSE, na.value = "#E5E7EB") +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = NULL,
      fill = "Biomarker evidence set",
      caption = paste(
        "Condensed Sankey of Figure 1C-E and Figure S2A-C evidence: only final significant OS, DepMap, and immune-linked marker routes are shown.",
        "Functional evidence nodes are split by EOBC state; DepMap endpoints are further split by drug class and sensitivity/resistance.",
        "Ribbons keep EOBC-state color as the main hue, then receive subtle biomarker/family/Meth-RNA/RF and functional-evidence tinting.",
        "OS uses nominal Cox p < 0.05; DepMap uses AUC/IC50 dose-response curve p < 0.05; immune routes require concordant non-zero TIL/TMB directions after excluding six zero-TMB samples for TMB correlations.",
        sep = "\n"
      )
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 15, face = "bold", color = "#111827", margin = margin(t = 12)),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 25.5, face = "bold", color = "#0F172A", margin = margin(b = 5)),
      plot.subtitle = element_text(size = 13.0, color = "#64748B", margin = margin(b = 18)),
      plot.caption = element_text(size = 9.5, color = "#64748B", hjust = 0, margin = margin(t = 16)),
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

copy_plot_alias <- function(source_stub, alias_stub) {
  for (dest_dir in unique(c(plot_dir, final_fig_dir, final_revision_dir))) {
    for (ext in c("png", "pdf")) {
      src <- file.path(dest_dir, paste0(source_stub, ".", ext))
      dst <- file.path(dest_dir, paste0(alias_stub, ".", ext))
      if (file.exists(src)) {
        file.copy(src, dst, overwrite = TRUE)
      }
    }
  }
}

plot_story_sankey(
  routes_story,
  title = "EOBC biomarker functional-evidence Sankey",
  subtitle = "Condensed view of Figure 1C-E and Figure S2A-C: final significant evidence is grouped by EOBC state, with DepMap split by drug class.",
  width = 20.5,
  height = 12.0,
  outfile_stub = "Figure_10C_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_full_R_v14",
  label_size = 1.60
)

plot_story_sankey(
  routes_story_priority,
  title = "Priority EOBC biomarker functional-evidence Sankey",
  subtitle = "Focused marker view retaining OS-linked and multi-domain biomarkers; endpoint evidence remains grouped by EOBC state.",
  width = 18.5,
  height = 9.5,
  outfile_stub = "Figure_10D_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_priority_R_v14",
  label_size = 1.78
)

copy_plot_alias(
  source_stub = "Figure_10C_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_full_R_v14",
  alias_stub = "Figure_10C_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_full_R_v12"
)

copy_plot_alias(
  source_stub = "Figure_10D_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_priority_R_v14",
  alias_stub = "Figure_10D_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_priority_R_v12"
)

copy_plot_alias(
  source_stub = "Figure_10D_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_priority_R_v13",
  alias_stub = "Figure_10D_EOBC_METH_group_pathway_evidence_biomarker_right_sankey_priority_R_v8"
)

message("Saved Sankey plots to: ", plot_dir)
message("Saved Sankey final copies to: ", final_fig_dir)
message("Saved Sankey revision copies to: ", final_revision_dir)
message("Saved Sankey tables to: ", table_dir)
