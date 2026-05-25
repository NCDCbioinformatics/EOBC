suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(survival)
  library(glmnet)
  library(scales)
})

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

base_dir <- env_path("EOBC_BIOMARKER_ROOT")
input_dir <- file.path(base_dir, "00_inputs_detected")
panel_path <- file.path(
  base_dir,
  "09_selected_relaxed_union_family6_final_signedMeth",
  "tables",
  "selected_union_biomarker_genes_family6_v18.csv"
)
rna_path <- file.path(input_dir, "TPM_young.csv")
meth_path <- file.path(input_dir, "MET_young_batch_JW.csv")
clinical_path <- file.path(input_dir, "total_sample_clinical_all.csv")

out_dir <- file.path(base_dir, "final_analysis", "os_elasticnet_full_candidate_r_v1")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

days_per_month <- 365.25 / 12
endpoints <- tibble(
  endpoint = c("overall_os", "os_5y", "os_10y"),
  endpoint_label = c("Overall OS", "5-year OS", "10-year OS"),
  horizon_days = c(Inf, 365.25 * 5, 365.25 * 10),
  x_month_max = c(300, 60, 120)
)

ink <- "#111827"
muted <- "#64748B"
grid_col <- "#D8E3EF"
axis_col <- "#94A3B8"
protective_blue <- "#55A9F7"
adverse_orange <- "#F4A261"
rna_green <- "#0B6E4F"
meth_orange <- "#BF5B17"
low_grey <- "#8FA3B5"

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
      panel.grid.major = element_line(colour = grid_col, linewidth = 0.38),
      panel.grid.minor = element_line(colour = alpha(grid_col, 0.45), linewidth = 0.22),
      axis.title = element_text(face = "bold", colour = ink),
      axis.text = element_text(colour = muted),
      strip.text = element_text(face = "bold", colour = ink, size = rel(1.03)),
      strip.background = element_rect(fill = "#EEF4FB", colour = "#D6E0EA", linewidth = 0.35),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", colour = ink),
      legend.text = element_text(colour = muted),
      plot.title = element_text(face = "bold", size = rel(1.45), colour = ink),
      plot.subtitle = element_text(colour = muted, size = rel(0.95)),
      plot.caption = element_text(colour = muted, size = rel(0.72), hjust = 0),
      plot.title.position = "plot"
    )
}

safe_stage_number <- function(x) {
  x <- as.character(x)
  out <- str_extract(x, "\\d+")
  suppressWarnings(as.numeric(out))
}

zscore_matrix <- function(mat) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  center <- colMeans(mat, na.rm = TRUE)
  scale <- apply(mat, 2, sd, na.rm = TRUE)
  scale[is.na(scale) | scale == 0] <- 1
  z <- sweep(sweep(mat, 2, center, "-"), 2, scale, "/")
  z[!is.finite(z)] <- 0
  z
}

read_omics_matrix <- function(path, genes, transform = c("none", "log2")) {
  transform <- match.arg(transform)
  raw <- read_csv(path, show_col_types = FALSE)
  names(raw)[1] <- "gene"
  raw <- raw %>%
    mutate(gene = as.character(gene)) %>%
    filter(gene %in% genes) %>%
    distinct(gene, .keep_all = TRUE)

  missing_genes <- setdiff(genes, raw$gene)
  mat_df <- raw %>% arrange(match(gene, genes))
  row_genes <- mat_df$gene
  mat_df <- mat_df %>% select(-gene)
  names(mat_df) <- str_replace_all(names(mat_df), '^"|"$', "")
  mat <- t(as.matrix(suppressWarnings(as.data.frame(lapply(mat_df, as.numeric)))))
  colnames(mat) <- row_genes
  rownames(mat) <- names(mat_df)
  mat <- mat[str_detect(rownames(mat), "\\.T$"), , drop = FALSE]

  if (transform == "log2") {
    mat <- log2(mat + 1)
  }

  list(
    matrix = zscore_matrix(mat),
    present_genes = intersect(genes, row_genes),
    missing_genes = missing_genes
  )
}

make_endpoint <- function(clinical, endpoint_name) {
  ep <- endpoints %>% filter(endpoint == endpoint_name)
  horizon <- ep$horizon_days[[1]]
  clinical %>%
    mutate(
      time_days = if (is.infinite(horizon)) overall_days else pmin(overall_days, horizon),
      event = if_else(death_event == 1 & overall_days <= horizon, 1, 0)
    ) %>%
    filter(is.finite(time_days), time_days > 0, !is.na(event))
}

balanced_foldid <- function(event, nfolds = 5) {
  n <- length(event)
  nfolds <- max(3, min(nfolds, sum(event == 1), sum(event == 0), n))
  foldid <- integer(n)
  for (cls in c(0, 1)) {
    idx <- which(event == cls)
    idx <- sample(idx)
    foldid[idx] <- rep(seq_len(nfolds), length.out = length(idx))
  }
  foldid
}

fit_cv_cox <- function(x, time, event, alpha = 0.5) {
  y <- Surv(time, event)
  nfolds <- max(3, min(5, sum(event == 1), sum(event == 0)))
  foldid <- balanced_foldid(event, nfolds = nfolds)
  cv.glmnet(
    x = x,
    y = y,
    family = "cox",
    alpha = alpha,
    standardize = FALSE,
    type.measure = "deviance",
    nfolds = nfolds,
    foldid = foldid
  )
}

extract_model_coef <- function(cvfit) {
  coefs_1se <- as.matrix(coef(cvfit, s = "lambda.1se"))[, 1]
  lambda_used <- "lambda.1se"
  if (sum(abs(coefs_1se) > 1e-9) == 0) {
    coefs_1se <- as.matrix(coef(cvfit, s = "lambda.min"))[, 1]
    lambda_used <- "lambda.min"
  }
  list(coef = coefs_1se, lambda_used = lambda_used)
}

fallback_univariate <- function(x, time, event) {
  y <- Surv(time, event)
  out <- lapply(seq_len(ncol(x)), function(j) {
    fit <- try(coxph(y ~ x[, j]), silent = TRUE)
    if (inherits(fit, "try-error")) {
      return(tibble(gene = colnames(x)[j], beta = NA_real_, p = NA_real_))
    }
    s <- summary(fit)
    tibble(
      gene = colnames(x)[j],
      beta = unname(coef(fit)[1]),
      p = s$coefficients[1, "Pr(>|z|)"]
    )
  }) %>% bind_rows()
  best <- out %>% arrange(p) %>% slice(1)
  coefs <- rep(0, ncol(x))
  names(coefs) <- colnames(x)
  coefs[best$gene] <- best$beta
  coefs
}

stability_selection <- function(x, time, event, b = 100, alpha = 0.5) {
  genes <- colnames(x)
  counts <- setNames(rep(0, length(genes)), genes)
  coef_sum <- setNames(rep(0, length(genes)), genes)
  valid <- 0
  n <- nrow(x)

  for (i in seq_len(b)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    if (sum(event[idx] == 1) < 3 || length(unique(event[idx])) < 2) next
    cvb <- try(fit_cv_cox(x[idx, , drop = FALSE], time[idx], event[idx], alpha = alpha), silent = TRUE)
    if (inherits(cvb, "try-error")) next
    eb <- extract_model_coef(cvb)
    sel <- abs(eb$coef) > 1e-9
    counts[sel] <- counts[sel] + 1
    coef_sum <- coef_sum + eb$coef
    valid <- valid + 1
  }

  tibble(
    gene = genes,
    bootstrap_valid = valid,
    stability = if (valid > 0) counts / valid else 0,
    bootstrap_mean_coef = if (valid > 0) coef_sum / valid else 0
  )
}

cox_metrics <- function(dat, adjusted = FALSE) {
  dat <- dat %>% mutate(risk_std = as.numeric(scale(score)))
  if (adjusted) {
    model_dat <- dat %>%
      select(time_days, event, risk_std, age_num, stage_num) %>%
      filter(complete.cases(.))
    if (nrow(model_dat) < 20 || sum(model_dat$event) < 8) {
      return(tibble(hr = NA_real_, p = NA_real_, cindex = NA_real_, n_used = nrow(model_dat)))
    }
    fit <- try(coxph(Surv(time_days, event) ~ risk_std + age_num + stage_num, data = model_dat), silent = TRUE)
  } else {
    model_dat <- dat %>% select(time_days, event, risk_std)
    fit <- try(coxph(Surv(time_days, event) ~ risk_std, data = model_dat), silent = TRUE)
  }
  if (inherits(fit, "try-error")) {
    return(tibble(hr = NA_real_, p = NA_real_, cindex = NA_real_, n_used = nrow(model_dat)))
  }
  s <- summary(fit)
  tibble(
    hr = unname(s$coefficients["risk_std", "exp(coef)"]),
    p = unname(s$coefficients["risk_std", "Pr(>|z|)"]),
    cindex = unname(s$concordance[1]),
    n_used = nrow(model_dat)
  )
}

logrank_p <- function(dat) {
  fit <- survdiff(Surv(time_days, event) ~ risk_group, data = dat)
  pchisq(fit$chisq, length(fit$n) - 1, lower.tail = FALSE)
}

make_km_curve <- function(dat) {
  sf <- survfit(Surv(time_months, event) ~ risk_group, data = dat)
  ss <- summary(sf)
  curve <- tibble(
    time_months = ss$time,
    survival = ss$surv,
    n_risk = ss$n.risk,
    n_event = ss$n.event,
    n_censor = ss$n.censor,
    risk_group = str_remove(as.character(ss$strata), "^risk_group=")
  )
  baseline <- dat %>%
    distinct(risk_group) %>%
    transmute(time_months = 0, survival = 1, n_risk = NA_integer_,
              n_event = 0, n_censor = 0, risk_group)
  bind_rows(baseline, curve) %>% arrange(risk_group, time_months)
}

run_screen <- function(modality, modality_label, xmat, clinical, endpoint_name, panel, b = 100) {
  endpoint_df <- make_endpoint(clinical, endpoint_name)
  samples <- intersect(endpoint_df$Sample, rownames(xmat))
  endpoint_df <- endpoint_df %>% filter(Sample %in% samples) %>% arrange(match(Sample, samples))
  x <- xmat[endpoint_df$Sample, , drop = FALSE]
  keep <- apply(x, 2, sd, na.rm = TRUE) > 0
  x <- x[, keep, drop = FALSE]
  time <- endpoint_df$time_days
  event <- endpoint_df$event

  cvfit <- fit_cv_cox(x, time, event)
  extracted <- extract_model_coef(cvfit)
  coefs <- extracted$coef
  fallback_used <- FALSE
  if (sum(abs(coefs) > 1e-9) == 0) {
    coefs <- fallback_univariate(x, time, event)
    extracted$lambda_used <- "univariate_fallback"
    fallback_used <- TRUE
  }

  score <- as.numeric(x %*% coefs)
  if (sd(score, na.rm = TRUE) == 0 || all(!is.finite(score))) {
    coefs <- fallback_univariate(x, time, event)
    score <- as.numeric(x %*% coefs)
    extracted$lambda_used <- "univariate_fallback"
    fallback_used <- TRUE
  }

  score <- as.numeric(scale(score))
  if (median(score, na.rm = TRUE) > 0) {
    # Keep the Cox score direction: larger score should represent higher predicted hazard.
    score <- score
  }

  dat <- endpoint_df %>%
    mutate(
      score = score,
      risk_group = if_else(score >= median(score, na.rm = TRUE), "High-risk score", "Low-risk score"),
      time_months = time_days / days_per_month
    )

  unadj <- cox_metrics(dat, adjusted = FALSE)
  adj <- cox_metrics(dat, adjusted = TRUE)
  lr_p <- logrank_p(dat)
  stab <- stability_selection(x, time, event, b = b)

  ep_label <- endpoints %>% filter(endpoint == endpoint_name) %>% pull(endpoint_label)
  screen_id <- paste(modality, endpoint_name, sep = "_")
  screen_label <- paste(modality_label, ep_label, sep = " | ")

  coef_tbl <- tibble(
    modality = modality,
    modality_label = modality_label,
    endpoint = endpoint_name,
    endpoint_label = ep_label,
    screen_id = screen_id,
    screen_label = screen_label,
    gene = names(coefs),
    coefficient = as.numeric(coefs),
    selected = abs(coefs) > 1e-9
  ) %>%
    left_join(stab, by = "gene") %>%
    left_join(panel %>% select(gene, Layer, Family6, target_label), by = "gene")

  selected_genes <- coef_tbl %>% filter(selected) %>% arrange(desc(abs(coefficient))) %>% pull(gene)
  if (length(selected_genes) == 0) selected_genes <- "(none)"

  model_summary <- tibble(
    modality = modality,
    modality_label = modality_label,
    endpoint = endpoint_name,
    endpoint_label = ep_label,
    screen_id = screen_id,
    screen_label = screen_label,
    n_samples = nrow(dat),
    events = sum(dat$event),
    n_genes_available = ncol(x),
    n_selected = sum(coef_tbl$selected),
    selected_genes = paste(selected_genes, collapse = " + "),
    lambda_used = extracted$lambda_used,
    fallback_used = fallback_used,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    cv_min_deviance = min(cvfit$cvm, na.rm = TRUE),
    hr_per_sd = unadj$hr,
    cox_p = unadj$p,
    cindex = unadj$cindex,
    adjusted_hr_per_sd = adj$hr,
    adjusted_cox_p = adj$p,
    adjusted_cindex = adj$cindex,
    logrank_p = lr_p,
    high_n = sum(dat$risk_group == "High-risk score"),
    low_n = sum(dat$risk_group == "Low-risk score")
  )

  km_curve <- make_km_curve(dat) %>%
    mutate(
      modality = modality,
      modality_label = modality_label,
      endpoint = endpoint_name,
      endpoint_label = ep_label,
      screen_id = screen_id,
      screen_label = screen_label
    )

  score_tbl <- dat %>%
    select(Sample, time_days, time_months, event, age_num, stage_num, score, risk_group) %>%
    mutate(
      modality = modality,
      modality_label = modality_label,
      endpoint = endpoint_name,
      endpoint_label = ep_label,
      screen_id = screen_id,
      screen_label = screen_label
    )

  list(
    model_summary = model_summary,
    coef_tbl = coef_tbl,
    km_curve = km_curve,
    score_tbl = score_tbl
  )
}

panel <- read_csv(panel_path, show_col_types = FALSE) %>%
  mutate(gene = as.character(gene)) %>%
  distinct(gene, .keep_all = TRUE)
genes <- panel$gene

cached_tables <- c(
  file.path(table_dir, "os_elasticnet_model_summary.csv"),
  file.path(table_dir, "os_elasticnet_coefficients_and_stability.csv"),
  file.path(table_dir, "os_elasticnet_patient_scores.csv"),
  file.path(table_dir, "os_elasticnet_km_curve_long.csv")
)

if (all(file.exists(cached_tables)) && Sys.getenv("EOBC_OS_FORCE_REFIT") != "1") {
  message("Using cached elastic-net OS tables. Set EOBC_OS_FORCE_REFIT=1 to refit models.")
  model_summary <- read_csv(cached_tables[[1]], show_col_types = FALSE)
  coef_tbl <- read_csv(cached_tables[[2]], show_col_types = FALSE)
  score_tbl <- read_csv(cached_tables[[3]], show_col_types = FALSE)
  km_curve <- read_csv(cached_tables[[4]], show_col_types = FALSE)
} else {
  clinical <- read_csv(clinical_path, show_col_types = FALSE) %>%
    filter(tpye == "Tumor") %>%
    mutate(
      Sample = as.character(Row.names),
      overall_days = as.numeric(overall_survival),
      death_event = case_when(alive == "No" ~ 1, alive == "Yes" ~ 0, TRUE ~ NA_real_),
      age_num = as.numeric(age),
      stage_num = safe_stage_number(Stage_final)
    ) %>%
    filter(is.finite(overall_days), !is.na(death_event), overall_days > 0)

  message("Reading RNA and methylation matrices...")
  rna <- read_omics_matrix(rna_path, genes, transform = "log2")
  meth <- read_omics_matrix(meth_path, genes, transform = "none")

  screen_grid <- tidyr::crossing(
    modality = c("RNA", "METH"),
    endpoint = endpoints$endpoint
  )

  message("Running elastic-net Cox screens with bootstrap stability selection...")
  results <- vector("list", nrow(screen_grid))
  for (i in seq_len(nrow(screen_grid))) {
    modality <- screen_grid$modality[[i]]
    endpoint_name <- screen_grid$endpoint[[i]]
    modality_label <- if_else(modality == "RNA", "RNA", "Methylation")
    xmat <- if (modality == "RNA") rna$matrix else meth$matrix
    message("  ", modality_label, " | ", endpoint_name)
    results[[i]] <- run_screen(
      modality = modality,
      modality_label = modality_label,
      xmat = xmat,
      clinical = clinical,
      endpoint_name = endpoint_name,
      panel = panel,
      b = 100
    )
  }

  model_summary <- bind_rows(lapply(results, `[[`, "model_summary"))
  coef_tbl <- bind_rows(lapply(results, `[[`, "coef_tbl"))
  km_curve <- bind_rows(lapply(results, `[[`, "km_curve"))
  score_tbl <- bind_rows(lapply(results, `[[`, "score_tbl"))

  write_csv(model_summary, cached_tables[[1]])
  write_csv(coef_tbl, cached_tables[[2]])
  write_csv(score_tbl, cached_tables[[3]])
  write_csv(km_curve, cached_tables[[4]])
}

screen_levels <- model_summary %>%
  mutate(screen_label = factor(screen_label, levels = c(
    "RNA | Overall OS", "RNA | 5-year OS", "RNA | 10-year OS",
    "Methylation | Overall OS", "Methylation | 5-year OS", "Methylation | 10-year OS"
  ))) %>%
  arrange(screen_label) %>%
  pull(screen_label)

coef_plot <- coef_tbl %>%
  mutate(
    screen_label = factor(screen_label, levels = screen_levels),
    gene_label = paste0(gene, " [", recode(Layer, "Transcriptome" = "R", "Methylation" = "M", .default = Layer), "]"),
    coef_plot = if_else(selected, coefficient, 0),
    support_plot = pmax(stability, if_else(selected, 0.06, 0)),
    selected_label = factor(if_else(selected, "Selected", "Not selected"), levels = c("Not selected", "Selected"))
  )

gene_order <- coef_plot %>%
  group_by(gene_label, Family6) %>%
  summarise(max_support = max(support_plot, na.rm = TRUE),
            max_abs_coef = max(abs(coefficient), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(max_support), desc(max_abs_coef), Family6, gene_label) %>%
  pull(gene_label)

coef_plot <- coef_plot %>%
  mutate(gene_label = factor(gene_label, levels = rev(gene_order)))

fig_stability <- ggplot(coef_plot, aes(x = screen_label, y = gene_label)) +
  geom_tile(aes(fill = coefficient), colour = "white", linewidth = 0.35, alpha = 0.94) +
  geom_point(
    aes(size = support_plot, shape = selected_label),
    colour = "#111827",
    fill = "white",
    stroke = 0.35,
    alpha = 0.92
  ) +
  geom_point(
    data = coef_plot %>% filter(selected),
    aes(colour = Family6),
    size = 1.25,
    shape = 8,
    stroke = 0.55,
    show.legend = FALSE
  ) +
  scale_fill_gradient2(
    low = protective_blue,
    mid = "white",
    high = adverse_orange,
    midpoint = 0,
    name = "Cox coefficient\n(high feature value)"
  ) +
  scale_colour_manual(values = family_cols, drop = FALSE) +
  scale_size_continuous(
    range = c(0.25, 5.4),
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    name = "Bootstrap\nselection rate"
  ) +
  scale_shape_manual(values = c("Not selected" = 21, "Selected" = 24), name = "Final model") +
  labs(
    title = "A. Full-candidate OS Cox model: selected biomarkers and bootstrap stability",
    subtitle = "Elastic-net Cox uses all 26 biomarkers simultaneously, allowing the model to choose the effective signature size without a fixed 3-gene cap.",
    x = NULL,
    y = NULL,
    caption = "Star marks the final elastic-net model coefficient; point size shows bootstrap selection rate across 100 resamples."
  ) +
  theme_eobc(base_size = 8.6) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 6.5),
    panel.grid = element_blank(),
    legend.box = "horizontal"
  )

ggsave(file.path(plot_dir, "Figure_01_OS_elasticnet_stability_matrix_R_v1.png"),
       fig_stability, width = 10.8, height = 8.8, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_01_OS_elasticnet_stability_matrix_R_v1.pdf"),
       fig_stability, width = 10.8, height = 8.8, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_01_OS_elasticnet_stability_matrix_R_v2.png"),
       fig_stability, width = 10.8, height = 8.8, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_01_OS_elasticnet_stability_matrix_R_v2.pdf"),
       fig_stability, width = 10.8, height = 8.8, bg = "white", device = cairo_pdf)

summary_plot <- model_summary %>%
  mutate(
    screen_label = factor(screen_label, levels = screen_levels),
    modality_colour = if_else(modality == "RNA", rna_green, meth_orange),
    complexity_label = paste0(n_selected, " genes"),
    perf_label = paste0("HR=", number(hr_per_sd, accuracy = 0.01),
                        "\nC=", number(cindex, accuracy = 0.01),
                        "\nP=", pvalue(cox_p, accuracy = 0.001))
  )

fig_complexity <- ggplot(summary_plot, aes(x = screen_label, y = n_selected, fill = modality)) +
  geom_col(width = 0.64, colour = "#263241", linewidth = 0.38) +
  geom_text(aes(label = complexity_label), vjust = -0.45, size = 3.1, fontface = "bold", colour = ink) +
  geom_text(aes(y = pmax(n_selected * 0.50, 0.55), label = str_wrap(selected_genes, width = 22)),
            size = 2.55, colour = "white", fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(
    values = c("RNA" = rna_green, "METH" = meth_orange),
    breaks = c("RNA", "METH"),
    labels = c("RNA", "Methylation")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)), breaks = pretty_breaks()) +
  labs(
    title = "B. Model-selected OS signature size by screen",
    subtitle = "The selected size is data-driven by penalized Cox regularization, not pre-fixed to one, two, or three genes.",
    x = NULL,
    y = "Number of selected biomarkers",
    fill = "Omics layer"
  ) +
  theme_eobc(base_size = 8.8) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"))

ggsave(file.path(plot_dir, "Figure_02_OS_elasticnet_model_complexity_R_v1.png"),
       fig_complexity, width = 10.3, height = 4.8, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02_OS_elasticnet_model_complexity_R_v1.pdf"),
       fig_complexity, width = 10.3, height = 4.8, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_02_OS_elasticnet_model_complexity_R_v2.png"),
       fig_complexity, width = 10.3, height = 4.8, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_02_OS_elasticnet_model_complexity_R_v2.pdf"),
       fig_complexity, width = 10.3, height = 4.8, bg = "white", device = cairo_pdf)

km_plot_data <- km_curve %>%
  left_join(model_summary %>% select(screen_id, selected_genes, hr_per_sd, cox_p, logrank_p, high_n, low_n, x_month_max = endpoint), by = "screen_id") %>%
  left_join(endpoints %>% select(endpoint, x_month_max), by = c("endpoint" = "endpoint")) %>%
  mutate(
    screen_label = factor(screen_label, levels = screen_levels),
    line_group = case_when(
      risk_group == "Low-risk score" ~ "Low-risk score",
      modality == "RNA" ~ "High-risk score (RNA model)",
      TRUE ~ "High-risk score (methylation model)"
    )
  )

km_ann <- model_summary %>%
  left_join(endpoints %>% select(endpoint, x_month_max), by = "endpoint") %>%
  mutate(
    screen_label = factor(screen_label, levels = screen_levels),
    x = x_month_max * 0.03,
    y = 0.605,
    label = paste0(
      str_wrap(selected_genes, width = 25),
      "\nHR=", number(hr_per_sd, accuracy = 0.01),
      " | Cox P=", pvalue(cox_p, accuracy = 0.001),
      "\nLog-rank P=", pvalue(logrank_p, accuracy = 0.001),
      " | High/Low n=", high_n, "/", low_n
    )
  )

fig_km <- ggplot(km_plot_data, aes(x = time_months, y = survival, colour = line_group)) +
  geom_step(linewidth = 0.82) +
  geom_point(
    data = km_plot_data %>% filter(n_censor > 0),
    shape = 3,
    size = 1.15,
    stroke = 0.35,
    show.legend = FALSE
  ) +
  geom_label(
    data = km_ann,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0,
    size = 2.15,
    lineheight = 0.86,
    colour = ink,
    fill = alpha("white", 0.88),
    linewidth = 0.22,
    label.padding = unit(0.11, "lines")
  ) +
  facet_wrap(~ screen_label, ncol = 3, scales = "free_x") +
  scale_colour_manual(
    values = c(
      "High-risk score (RNA model)" = rna_green,
      "High-risk score (methylation model)" = meth_orange,
      "Low-risk score" = low_grey
    ),
    breaks = c("High-risk score (RNA model)", "High-risk score (methylation model)", "Low-risk score"),
    name = NULL
  ) +
  coord_cartesian(ylim = c(0.56, 1.02), clip = "off") +
  labs(
    title = "C. Kaplan-Meier validation of full-candidate elastic-net OS scores",
    subtitle = "High-risk and low-risk groups are defined by the median model score within each omics layer and endpoint.",
    x = "Months",
    y = "Survival probability"
  ) +
  theme_eobc(base_size = 8.6) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    strip.text = element_text(size = 9.4)
  )

ggsave(file.path(plot_dir, "Figure_03_OS_elasticnet_KM_atlas_R_v1.png"),
       fig_km, width = 13.8, height = 8.2, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03_OS_elasticnet_KM_atlas_R_v1.pdf"),
       fig_km, width = 13.8, height = 8.2, bg = "white", device = cairo_pdf)
ggsave(file.path(plot_dir, "Figure_03_OS_elasticnet_KM_atlas_R_v2.png"),
       fig_km, width = 13.8, height = 8.2, dpi = 360, bg = "white")
ggsave(file.path(plot_dir, "Figure_03_OS_elasticnet_KM_atlas_R_v2.pdf"),
       fig_km, width = 13.8, height = 8.2, bg = "white", device = cairo_pdf)

analysis_note <- c(
  "EOBC OS elastic-net full-candidate analysis",
  "",
  "Rationale:",
  "- Previous exhaustive Cox screen intentionally evaluated only subset sizes 1-3, as recorded in os_phase1_subset_screen_v1/logs/run_info.txt.",
  "- This analysis removes the arbitrary 3-gene cap by entering all 26 candidate biomarkers simultaneously into an elastic-net Cox model.",
  "- The selected signature size is therefore model-driven. Bootstrap stability selection estimates whether each biomarker is repeatedly selected.",
  "- With only 10-25 OS events per endpoint, this regularized model is more defensible than unpenalized exhaustive testing of all 67,108,863 possible subsets.",
  "",
  "Outputs:",
  "- Figure_01_OS_elasticnet_stability_matrix_R_v1: selected genes, coefficients, and bootstrap stability.",
  "- Figure_02_OS_elasticnet_model_complexity_R_v1: selected model size per screen.",
  "- Figure_03_OS_elasticnet_KM_atlas_R_v1: KM validation of each model score.",
  "- Tables contain coefficients/stability, patient scores, KM curves, and model-level metrics."
)
writeLines(analysis_note, file.path(out_dir, "analysis_note.txt"))

message("Done.")
message("Plots written to: ", plot_dir)
message("Tables written to: ", table_dir)
