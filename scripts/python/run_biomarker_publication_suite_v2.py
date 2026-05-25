from __future__ import annotations

import importlib.util
import math
import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test


SCRIPT_DIR = Path(__file__).resolve().parent
OS_CONSENSUS_SCRIPT = SCRIPT_DIR / "run_os_consensus_analysis_v3.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import module at {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


oscons = load_module(OS_CONSENSUS_SCRIPT, "os_consensus_v3")
phase1_script = SCRIPT_DIR / "run_os_phase1_subset_screen.py"
if not phase1_script.exists():
    phase1_script = oscons.PHASE1_ROOT.parent / "run_os_phase1_subset_screen.py"
phase1 = load_module(phase1_script, "os_phase1_base")

BIOMARKER_ROOT = Path(os.environ.get("EOBC_BIOMARKER_ROOT", oscons.BIOMARKER_ROOT)).expanduser()
FINAL_ANALYSIS_ROOT = Path(
    os.environ.get("EOBC_FINAL_ANALYSIS_DIR", oscons.FINAL_ANALYSIS_ROOT)
).expanduser()
VALIDATION_TABLES = Path(
    os.environ.get(
        "EOBC_VALIDATION_TABLES_DIR",
        BIOMARKER_ROOT
        / "13_final_manuscript_validation_suite_v5_full_rewrite_raw_depmap_SVS_MAF"
        / "tables",
    )
).expanduser()
OUTPUT_ROOT = Path(
    os.environ.get("EOBC_DOMAIN_PUBLICATION_OUT", FINAL_ANALYSIS_ROOT / "domain_publication_suite_v2")
).expanduser()
PLOTS_DIR = OUTPUT_ROOT / "plots"
TABLES_DIR = OUTPUT_ROOT / "tables"

FAMILY_COLORS = {
    "Immune": "#2a9d8f",
    "Repair": "#d62828",
    "Glycolysis / TCA": "#f77f00",
    "Hormone signaling": "#577590",
    "Fatty acid": "#f4a261",
    "Kinase signaling": "#7b2cbf",
}
FAMILY_ABBR = {
    "Immune": "Imm",
    "Repair": "Rep",
    "Glycolysis / TCA": "Gly/TCA",
    "Hormone signaling": "Horm",
    "Fatty acid": "FA",
    "Kinase signaling": "Kin",
}
LAYER_EDGE = {
    "Transcriptome": "#0b6e4f",
    "Methylation": "#8d5524",
}
LAYER_MARKERS = {"Transcriptome": "o", "Methylation": "s"}
MODALITY_COLORS = {"RNA": "#0b6e4f", "METH": "#bf5b17"}
PROTECTIVE_COLOR = "#69b3ff"
ADVERSE_COLOR = "#f4a261"
SENSITIVE_COLOR = "#3a86ff"
RESISTANT_COLOR = "#ff7f11"
KM_LOW_COLOR = "#94a3b8"
TITLE_BAR = "#334155"
GRID_COLOR = "#dbe4ee"
TEXT_COLOR = "#111827"
MUTED_TEXT = "#6b7280"

SCREEN_ORDER = [
    ("RNA", "overall_os"),
    ("RNA", "os_5y"),
    ("RNA", "os_10y"),
    ("METH", "overall_os"),
    ("METH", "os_5y"),
    ("METH", "os_10y"),
]
SCREEN_LABELS = {
    ("RNA", "overall_os"): "RNA\nOverall",
    ("RNA", "os_5y"): "RNA\n5y",
    ("RNA", "os_10y"): "RNA\n10y",
    ("METH", "overall_os"): "METH\nOverall",
    ("METH", "os_5y"): "METH\n5y",
    ("METH", "os_10y"): "METH\n10y",
}
ENDPOINT_DISPLAY = {
    "overall_os": "Overall OS",
    "os_5y": "5-year OS",
    "os_10y": "10-year OS",
}


def ensure_dirs() -> None:
    for path in [OUTPUT_ROOT, PLOTS_DIR, TABLES_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_style() -> None:
    plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans", "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["axes.facecolor"] = "white"
    plt.rcParams["figure.facecolor"] = "white"
    plt.rcParams["savefig.facecolor"] = "white"


def save_figure_bundle(fig: plt.Figure, output_path: Path, dpi: int = 320) -> None:
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")


def lighten(color: str, factor: float = 0.55) -> Tuple[float, float, float]:
    rgb = np.array(to_rgb(color))
    out = rgb + (1.0 - rgb) * factor
    return tuple(np.clip(out, 0.0, 1.0))


def darken(color: str, factor: float = 0.78) -> Tuple[float, float, float]:
    rgb = np.array(to_rgb(color))
    return tuple(np.clip(rgb * factor, 0.0, 1.0))


def layer_tag(layer: str) -> str:
    return "R" if layer == "Transcriptome" else "M"


def gene_label_with_meta(row: pd.Series) -> str:
    return f"{row['gene']} [{layer_tag(str(row['Layer']))}]"


def format_p(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3f}"


def header_bar(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.0,
        1.03,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        color="white",
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": TITLE_BAR,
            "edgecolor": TITLE_BAR,
        },
    )


def draw_chip_strip(ax: plt.Axes, df: pd.DataFrame, title: str) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(
        0.01,
        0.92,
        title,
        ha="left",
        va="top",
        fontsize=10,
        color=MUTED_TEXT,
        fontweight="bold",
    )
    if df.empty:
        ax.text(0.01, 0.48, "All biomarkers had usable validation data.", fontsize=9.3, color=MUTED_TEXT)
        return

    gap = 0.012
    chips: List[Tuple[str, str, float]] = []
    rows: List[List[Tuple[str, str, float]]] = [[]]
    current_x = 0.01

    for _, row in df.iterrows():
        family = str(row["Family6"])
        label = f"{row['gene']} [{FAMILY_ABBR.get(family, family)}]"
        width = min(0.22, 0.045 + 0.0115 * len(label))
        chips.append((label, family, width))
        if current_x + width > 0.985 and rows[-1]:
            rows.append([])
            current_x = 0.01
        rows[-1].append((label, family, width))
        current_x += width + gap

    n_rows = max(1, len(rows))
    top_y = 0.70 if n_rows <= 2 else 0.76
    bottom_y = 0.08
    row_step = (top_y - bottom_y) / max(n_rows - 1, 1) if n_rows > 1 else 0.0
    h = 0.20 if n_rows <= 2 else min(0.18, row_step * 0.72 if row_step else 0.18)
    font_size = 8.4 if n_rows <= 2 else 8.0

    for row_idx, row_items in enumerate(rows):
        x = 0.01
        y = top_y - row_idx * row_step
        for label, family, width in row_items:
            face = lighten(FAMILY_COLORS.get(family, "#94a3b8"), 0.82)
            edge = FAMILY_COLORS.get(family, "#94a3b8")
            ax.add_patch(
                FancyBboxPatch(
                    (x, y),
                    width,
                    h,
                    boxstyle="round,pad=0.018,rounding_size=0.04",
                    linewidth=1.2,
                    facecolor=face,
                    edgecolor=edge,
                )
            )
            ax.text(
                x + width / 2,
                y + h / 2,
                label,
                ha="center",
                va="center",
                fontsize=font_size,
                color=TEXT_COLOR,
            )
            x += width + gap


def annotate_scatter_points(
    ax: plt.Axes,
    df: pd.DataFrame,
    xcol: str,
    ycol: str,
    labelcol: str,
    sublabelcol: str | None = None,
    sort_cols: List[str] | None = None,
    offset_map: Dict[str, Tuple[float, float]] | None = None,
    collision_x_frac: float = 0.11,
    collision_y_frac: float = 0.06,
    label_fontsize: float = 9.2,
    sublabel_fontsize: float = 7.9,
) -> None:
    if df.empty:
        return
    use = df.copy()
    if sort_cols:
        use = use.sort_values(sort_cols, ascending=[False] * len(sort_cols))
    else:
        use = use.sort_values(ycol, ascending=False)

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_range = x_max - x_min
    y_range = y_max - y_min
    placed: List[Tuple[float, float]] = []

    for _, row in use.iterrows():
        x = float(row[xcol])
        y = float(row[ycol])
        if pd.isna(x) or pd.isna(y):
            continue
        row_key = str(row.get("gene", row[labelcol]))
        if offset_map and row_key in offset_map:
            dx_frac, dy_frac = offset_map[row_key]
            tx = x + dx_frac * x_range
            ty = y + dy_frac * y_range
        else:
            dx = 0.038 * x_range if x >= 0 else -0.038 * x_range
            tx = x + dx
            ty = y + 0.024 * y_range
        for step in range(18):
            conflict = any(
                abs(tx - px) < collision_x_frac * x_range and abs(ty - py) < collision_y_frac * y_range
                for px, py in placed
            )
            if not conflict:
                break
            ty += (0.028 * y_range) * (1 if step % 2 == 0 else -1)
        placed.append((tx, ty))

        leader = darken(FAMILY_COLORS.get(str(row["Family6"]), "#64748b"), 0.9)
        ax.plot([x, tx], [y, ty], color=leader, lw=1.0, alpha=0.85, zorder=4)
        ha = "left" if tx >= x else "right"
        ax.text(
            tx,
            ty,
            str(row[labelcol]),
            ha=ha,
            va="bottom",
            fontsize=label_fontsize,
            color=TEXT_COLOR,
            zorder=5,
        )
        if sublabelcol is not None:
            ax.text(
                tx,
                ty - 0.027 * y_range,
                str(row[sublabelcol]),
                ha=ha,
                va="top",
                fontsize=sublabel_fontsize,
                color=MUTED_TEXT,
                zorder=5,
            )


def ensure_os_consensus_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    metrics_path = oscons.TABLES_DIR / "os_gene_screen_metrics.csv"
    summary_path = oscons.TABLES_DIR / "os_gene_consensus_summary.csv"
    if metrics_path.exists() and summary_path.exists():
        metrics = pd.read_csv(metrics_path)
        summary = pd.read_csv(summary_path)
        meta = oscons.load_meta()
        return meta, metrics, summary

    oscons.ensure_dirs()
    meta = oscons.load_meta()
    volcano, top_hits, km_tables = oscons.load_screen_tables()
    metrics = oscons.build_screen_gene_metrics(meta, volcano, top_hits, km_tables)
    summary = oscons.build_gene_summary(metrics)
    return meta, metrics, summary


def load_core_gene_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta, metrics, summary = ensure_os_consensus_tables()
    gene_table = meta.merge(
        summary,
        on=["gene", "gene_label", "Layer", "Family6", "target_label"],
        how="left",
    )
    gene_table["row_label"] = gene_table.apply(gene_label_with_meta, axis=1)
    return gene_table, metrics


def build_depmap_all_biomarker_table(gene_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    depmap = pd.read_csv(VALIDATION_TABLES / "depmap_gene_drug_association_joint_AUC_IC50_annotated.csv")
    depmap = depmap[depmap["gene"].isin(gene_table["gene"])].copy()
    depmap["best_fdr"] = depmap[["auc_fdr", "ic50_fdr"]].min(axis=1, skipna=True)
    depmap["best_p"] = depmap[["p_auc", "p_ic50"]].min(axis=1, skipna=True)
    depmap["direction_score"] = pd.to_numeric(depmap["direction_score"], errors="coerce")
    depmap["combined_score"] = pd.to_numeric(depmap["combined_score"], errors="coerce")
    depmap["best_rank_value"] = depmap["combined_score"].fillna(-np.inf)

    best = (
        depmap.sort_values(
            ["gene", "best_rank_value", "best_fdr", "best_p"],
            ascending=[True, False, True, True],
        )
        .groupby("gene", as_index=False)
        .head(1)
        .copy()
    )

    best = best.rename(
        columns={
            "drug_clean": "best_drug",
            "drug_target_class": "best_target_class",
            "direction_score": "depmap_signed_score",
            "combined_score": "depmap_combined_score",
            "consistent_direction": "depmap_direction_text",
        }
    )
    best["depmap_neglog10"] = -np.log10(best["best_fdr"].fillna(best["best_p"]).clip(lower=1e-300))
    best["depmap_direction"] = np.where(
        best["depmap_signed_score"] < 0,
        "More sensitive when high",
        "More resistant when high",
    )
    best["depmap_signif_tier"] = np.select(
        [
            best["best_fdr"] < 0.10,
            best["best_fdr"] < 0.25,
        ],
        ["FDR < 0.10", "FDR < 0.25"],
        default="Exploratory",
    )

    merged = gene_table.merge(
        best[
            [
                "gene",
                "best_drug",
                "best_target_class",
                "depmap_signed_score",
                "depmap_combined_score",
                "best_fdr",
                "best_p",
                "depmap_neglog10",
                "depmap_direction",
                "depmap_signif_tier",
            ]
        ],
        on="gene",
        how="left",
    )
    merged["depmap_available"] = merged["best_drug"].notna()

    score = merged["depmap_combined_score"].fillna(0.0).abs()
    denom = float(score.max()) if float(score.max()) > 0 else 1.0
    merged["depmap_size_norm"] = score / denom
    merged["depmap_label"] = merged["gene"]
    merged["depmap_sublabel"] = merged["best_drug"].fillna("Not carried into tumor-intrinsic DepMap screen")
    merged.to_csv(TABLES_DIR / "depmap_all_biomarker_best_hits.csv", index=False)

    missing = merged[~merged["depmap_available"]].copy()
    missing.to_csv(TABLES_DIR / "depmap_biomarkers_without_intrinsic_model.csv", index=False)
    return merged, depmap


def build_immune_all_biomarker_table(gene_table: pd.DataFrame) -> pd.DataFrame:
    immune = pd.read_csv(VALIDATION_TABLES / "immune_gene_correlations_with_TIL_and_TMB.csv")
    immune = immune[immune["gene"].isin(gene_table["gene"])].copy()
    immune["immune_best_fdr"] = immune[["til_fdr", "tmb_fdr"]].min(axis=1, skipna=True)
    immune["immune_strength"] = np.sqrt(immune["rho_til"].pow(2) + immune["rho_tmb"].pow(2))
    immune["immune_signif_tier"] = np.select(
        [
            immune["immune_best_fdr"] < 0.10,
            immune["immune_best_fdr"] < 0.25,
        ],
        ["FDR < 0.10", "FDR < 0.25"],
        default="Exploratory",
    )
    immune["immune_quadrant"] = np.select(
        [
            (immune["rho_til"] > 0) & (immune["rho_tmb"] > 0),
            (immune["rho_til"] < 0) & (immune["rho_tmb"] < 0),
            (immune["rho_til"] < 0) & (immune["rho_tmb"] > 0),
            (immune["rho_til"] > 0) & (immune["rho_tmb"] < 0),
        ],
        [
            "Immune-hot / TMB-high",
            "Immune-cold / TMB-low",
            "TMB-shifted",
            "TIL-high / TMB-low",
        ],
        default="Mixed / borderline",
    )

    merged = gene_table.merge(
        immune[
            [
                "gene",
                "rho_til",
                "rho_tmb",
                "til_fdr",
                "tmb_fdr",
                "immune_best_fdr",
                "immune_strength",
                "immune_signif_tier",
                "immune_quadrant",
            ]
        ],
        on="gene",
        how="left",
    )
    merged["immune_available"] = merged["rho_til"].notna() & merged["rho_tmb"].notna()
    denom = float(merged["immune_strength"].fillna(0).max())
    merged["immune_size_norm"] = merged["immune_strength"].fillna(0) / (denom if denom > 0 else 1.0)
    merged["immune_label"] = merged["gene"]
    merged.to_csv(TABLES_DIR / "immune_all_biomarker_landscape.csv", index=False)
    return merged


def compute_single_gene_os_table(gene_table: pd.DataFrame) -> pd.DataFrame:
    panel = phase1.load_panel(phase1.SELECTED_PANEL)
    genes = panel["gene"].astype(str).tolist()
    clinical = phase1.load_clinical(phase1.CLINICAL_PATH)
    rna = phase1.zscore_df(phase1.load_omics_matrix(phase1.RNA_PATH, genes))
    meth = phase1.zscore_df(phase1.load_omics_matrix(phase1.METH_PATH, genes))
    layer_map = panel.set_index("gene")["Layer"].to_dict()
    rows: List[Dict[str, object]] = []

    for gene in genes:
        modality = "RNA" if str(layer_map[gene]) == "Transcriptome" else "METH"
        zmat = rna if modality == "RNA" else meth
        if gene not in zmat.columns:
            continue
        modality_clinical = clinical[clinical["Sample"].isin(set(zmat.index))].copy()
        for endpoint in phase1.HORIZONS:
            endpoint_df = phase1.build_endpoint_df(modality_clinical, endpoint)
            endpoint_df = endpoint_df[endpoint_df["Sample"].isin(zmat.index)].copy()
            endpoint_df["score"] = endpoint_df["Sample"].map(zmat[gene])
            endpoint_df["time_months"] = endpoint_df["time_days"] / phase1.DAYS_PER_MONTH
            median_score = endpoint_df["score"].median()
            endpoint_df["group"] = np.where(endpoint_df["score"] >= median_score, "High", "Low")

            uni = phase1.fit_univariate_score_model(endpoint_df[["time_days", "event", "score"]].dropna())
            adj = phase1.fit_adjusted_score_model(
                endpoint_df[["time_days", "event", "score", "age_num", "stage_num"]].dropna()
            )
            high = endpoint_df[endpoint_df["group"] == "High"].copy()
            low = endpoint_df[endpoint_df["group"] == "Low"].copy()
            logrank_p = logrank_test(
                high["time_months"],
                low["time_months"],
                event_observed_A=high["event"],
                event_observed_B=low["event"],
            ).p_value

            rows.append(
                {
                    "gene": gene,
                    "modality": modality,
                    "endpoint": endpoint,
                    "endpoint_label": ENDPOINT_DISPLAY[endpoint],
                    "n_high": int(high.shape[0]),
                    "n_low": int(low.shape[0]),
                    "events_high": int(high["event"].sum()),
                    "events_low": int(low["event"].sum()),
                    "score_median": float(median_score),
                    "logrank_p": float(logrank_p),
                    **uni,
                    **adj,
                }
            )

    stats = pd.DataFrame(rows)
    stats["rank_metric"] = stats["p_value"].fillna(1.0)
    best = (
        stats.sort_values(["gene", "rank_metric", "logrank_p"], ascending=[True, True, True])
        .groupby("gene", as_index=False)
        .head(1)
        .copy()
    )
    best = best.merge(
        gene_table[
            [
                "gene",
                "gene_label",
                "Layer",
                "Family6",
                "target_label",
                "consensus_strength",
                "top1_hits",
                "consensus_class",
                "best_screen",
            ]
        ],
        on="gene",
        how="left",
    )
    best["display_rank_score"] = (
        -np.log10(best["p_value"].clip(lower=1e-300))
        + 0.15 * best["top1_hits"].fillna(0)
        + 0.08 * best["consensus_strength"].fillna(0)
    )
    best["km_priority"] = best["p_value"].fillna(1.0)
    best.to_csv(TABLES_DIR / "os_single_gene_survival_summary.csv", index=False)
    return best


def plot_os_gene_consensus_matrix_v2(gene_table: pd.DataFrame, metrics: pd.DataFrame) -> None:
    gene_order = gene_table["gene"].tolist()
    n_rows = len(gene_order)
    fig, ax = plt.subplots(figsize=(12.4, max(8.3, 0.33 * n_rows)))
    ax.set_xlim(-1.25, len(SCREEN_ORDER) - 0.35)
    ax.set_ylim(-0.7, n_rows - 0.3)
    ax.axvspan(-0.5, 2.5, color="#eef6ff", alpha=0.95, zorder=0)
    ax.axvspan(2.5, 5.5, color="#fff2e8", alpha=0.95, zorder=0)

    for y in np.arange(-0.5, n_rows, 1):
        ax.axhline(y, color="#edf2f7", lw=0.8, zorder=0)
    for x in np.arange(-0.5, len(SCREEN_ORDER), 1):
        ax.axvline(x, color="#d7e3f0", lw=0.8, zorder=0)

    max_abs = float(metrics["screen_contrib"].abs().max()) if not metrics.empty else 1.0
    contrib_cmap = LinearSegmentedColormap.from_list(
        "os_v2_diverging",
        ["#69b3ff", "#f8fafc", "#f4a261"],
        N=256,
    )
    norm = Normalize(vmin=-max_abs, vmax=max_abs)
    metric_map = metrics.set_index(["gene", "screen"])
    y_positions = {gene: n_rows - 1 - idx for idx, gene in enumerate(gene_order)}

    for _, row in gene_table.iterrows():
        gene = str(row["gene"])
        y = y_positions[gene]
        family = str(row["Family6"])
        layer = str(row["Layer"])
        ax.add_patch(
            Rectangle(
                (-1.18, y - 0.38),
                0.18,
                0.76,
                facecolor=FAMILY_COLORS.get(family, "#94a3b8"),
                edgecolor="none",
                zorder=1,
            )
        )
        ax.text(
            -0.92,
            y,
            f"{row['gene_label']} [{FAMILY_ABBR.get(family, family)}]",
            ha="right",
            va="center",
            fontsize=9.4,
            color=TEXT_COLOR,
        )
        for x_idx, (modality, endpoint) in enumerate(SCREEN_ORDER):
            screen = f"{modality}_{endpoint}"
            if (gene, screen) not in metric_map.index:
                continue
            cell = metric_map.loc[(gene, screen)]
            support = float(cell["support_top30"])
            contrib = float(cell["screen_contrib"])
            size = 26 + 650 * support
            face = contrib_cmap(norm(contrib))
            ax.scatter(
                x_idx,
                y,
                s=size,
                c=[face],
                marker=LAYER_MARKERS.get(layer, "o"),
                edgecolors=FAMILY_COLORS.get(family, "#64748b"),
                linewidths=1.35,
                zorder=3,
            )
            if bool(cell["top1_hit"]):
                ax.scatter(
                    x_idx,
                    y,
                    s=max(48, size * 0.18),
                    marker="*",
                    c="#ffd166",
                    edgecolors="white",
                    linewidths=0.35,
                    zorder=4,
                )

    ax.text(1.0, n_rows - 0.05, "RNA screens", ha="center", va="bottom", fontsize=11, fontweight="bold", color=MODALITY_COLORS["RNA"])
    ax.text(4.0, n_rows - 0.05, "Methylation screens", ha="center", va="bottom", fontsize=11, fontweight="bold", color=MODALITY_COLORS["METH"])
    ax.set_xticks(range(len(SCREEN_ORDER)))
    ax.set_xticklabels([SCREEN_LABELS[key] for key in SCREEN_ORDER], fontsize=10)
    ax.set_yticks([])
    ax.tick_params(length=0)
    header_bar(ax, "A. Gene-level OS consensus extracted from exhaustive subset screening")

    for spine in ax.spines.values():
        spine.set_visible(False)

    size_vals = [0.10, 0.30, 0.60]
    size_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor="#94a3b8",
            markeredgecolor="#94a3b8",
            markersize=math.sqrt(26 + 650 * val) / 1.8,
            label=f"{int(val * 100)}% of top 30",
        )
        for val in size_vals
    ]
    layer_handles = [
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="white", markeredgecolor="#4b5563", markersize=8, label="Transcriptome"),
        Line2D([0], [0], marker="s", linestyle="", color="w", markerfacecolor="white", markeredgecolor="#4b5563", markersize=8, label="Methylation"),
        Line2D([0], [0], marker="*", linestyle="", color="#ffd166", markeredgecolor="#ffffff", markersize=10, label="Present in screen-leading signature"),
    ]
    family_handles = [Patch(facecolor=FAMILY_COLORS[f], edgecolor=FAMILY_COLORS[f], label=f) for f in FAMILY_COLORS]
    legend_style = {
        "frameon": True,
        "fancybox": True,
        "framealpha": 0.96,
        "facecolor": "white",
        "edgecolor": "#e2e8f0",
    }
    fig.legend(
        handles=size_handles,
        fontsize=8.4,
        loc="lower left",
        bbox_to_anchor=(0.07, 0.065),
        ncol=3,
        title="Top-subset support",
        title_fontsize=8.6,
        **legend_style,
    )
    fig.legend(
        handles=family_handles,
        fontsize=8.4,
        loc="lower center",
        bbox_to_anchor=(0.50, 0.065),
        ncol=3,
        title="Biological class",
        title_fontsize=8.6,
        **legend_style,
    )
    fig.legend(
        handles=layer_handles,
        fontsize=8.3,
        loc="lower right",
        bbox_to_anchor=(0.92, 0.058),
        ncol=1,
        title="Layer / top hit",
        title_fontsize=8.6,
        **legend_style,
    )

    sm = plt.cm.ScalarMappable(norm=norm, cmap=contrib_cmap)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.026, pad=0.02)
    cbar.ax.set_ylabel("Signed OS consensus in top subsets  (protective  ←   →  adverse)", fontsize=8.7)
    cbar.ax.tick_params(labelsize=8)

    ax.text(
        0.0,
        -0.22,
        "Cell fill encodes screen-specific consensus direction; biological-class palette is carried by row strips and marker outlines.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.6,
        color=MUTED_TEXT,
    )
    fig.subplots_adjust(left=0.20, right=0.91, top=0.96, bottom=0.17)
    save_figure_bundle(fig, PLOTS_DIR / "Figure_01A_OS_gene_consensus_matrix_v2.png")
    plt.close(fig)


def plot_os_temporal_scatter_v2(gene_table: pd.DataFrame) -> None:
    df = gene_table.copy()
    fig, ax = plt.subplots(figsize=(10.4, 8.4))
    header_bar(ax, "B. Temporal OS behavior of each biomarker")

    x = df["early_score"].astype(float)
    y = df["late_score"].astype(float)
    x_pad = max(0.08, float(np.nanmax(np.abs(x))) * 0.12)
    y_pad = max(0.08, float(np.nanmax(np.abs(y))) * 0.12)
    x_lim = (float(np.nanmin(x)) - x_pad, float(np.nanmax(x)) + x_pad)
    y_lim = (float(np.nanmin(y)) - y_pad, float(np.nanmax(y)) + y_pad)

    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    ax.axvspan(x_lim[0], 0, color="#eaf4ff", alpha=0.75, zorder=0)
    ax.axvspan(0, x_lim[1], color="#fff0e4", alpha=0.75, zorder=0)
    ax.axhspan(y_lim[0], 0, color="#f5fbff", alpha=0.45, zorder=0)
    ax.axhspan(0, y_lim[1], color="#fff7ef", alpha=0.45, zorder=0)
    ax.axvline(0, color="#6b7280", lw=1.1, ls="--")
    ax.axhline(0, color="#6b7280", lw=1.1, ls="--")

    size_raw = df["consensus_strength"].fillna(0.0).astype(float)
    denom = float(size_raw.max()) if float(size_raw.max()) > 0 else 1.0
    df["bubble_size"] = 120 + 1600 * (size_raw / denom)

    for _, row in df.iterrows():
        family = str(row["Family6"])
        layer = str(row["Layer"])
        face = FAMILY_COLORS.get(family, "#94a3b8")
        ax.scatter(
            float(row["early_score"]),
            float(row["late_score"]),
            s=float(row["bubble_size"]),
            marker=LAYER_MARKERS.get(layer, "o"),
            c=[face],
            edgecolors=LAYER_EDGE.get(layer, "#334155"),
            linewidths=1.4,
            alpha=0.92,
            zorder=3,
        )
        if float(row.get("top1_hits", 0) or 0) > 0:
            ax.scatter(
                float(row["early_score"]),
                float(row["late_score"]),
                s=float(row["bubble_size"]) * 1.18,
                facecolors="none",
                edgecolors="#1f2937",
                linewidths=1.2,
                zorder=4,
            )

    label_df = df.sort_values(
        ["top1_hits", "consensus_strength", "support_sum"],
        ascending=[False, False, False],
    ).head(10)
    label_df = label_df.assign(label_text=label_df["row_label"])
    temporal_offsets = {
        "PRR15": (-0.08, -0.03),
        "BCL2A1": (-0.09, -0.02),
        "SPDEF": (-0.06, -0.06),
        "ZMYND10": (-0.08, -0.01),
        "MARVELD2": (-0.07, 0.02),
        "SKA3": (-0.03, 0.025),
        "AURKB": (0.05, 0.015),
        "RASAL3": (-0.06, -0.02),
        "TRAF3IP3": (0.08, -0.04),
        "C1QTNF6": (0.03, 0.03),
        "MYO1G": (0.05, 0.04),
        "MLPH": (0.05, 0.03),
        "KLHDC7B": (0.03, 0.03),
        "TFF1": (0.05, -0.03),
    }
    annotate_scatter_points(
        ax,
        label_df,
        "early_score",
        "late_score",
        "label_text",
        sublabelcol="consensus_class",
        sort_cols=["consensus_strength"],
        offset_map=temporal_offsets,
        collision_x_frac=0.09,
        collision_y_frac=0.05,
        label_fontsize=8.9,
        sublabel_fontsize=7.5,
    )

    ax.set_xlabel("Short-term OS consensus (5-year subset screens)", fontsize=12, fontweight="bold")
    ax.set_ylabel("Long-term OS consensus (10-year subset screens)", fontsize=12, fontweight="bold")
    ax.grid(color=GRID_COLOR, lw=0.9)
    ax.set_axisbelow(True)
    ax.text(0.02, 0.96, "5y/10y protective", transform=ax.transAxes, ha="left", va="top", fontsize=10.2, color="#3b82f6")
    ax.text(0.98, 0.96, "Persistent adverse", transform=ax.transAxes, ha="right", va="top", fontsize=10.2, color="#d97706")
    ax.text(0.02, 0.05, "5y protective / 10y adverse", transform=ax.transAxes, ha="left", va="bottom", fontsize=9.5, color="#b45309")
    ax.text(0.98, 0.05, "5y adverse / 10y protective", transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5, color="#2563eb")

    family_handles = [Patch(facecolor=FAMILY_COLORS[f], edgecolor=FAMILY_COLORS[f], label=f) for f in FAMILY_COLORS]
    layer_handles = [
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#9ca3af", markeredgecolor=LAYER_EDGE["Transcriptome"], markersize=8, label="Transcriptome"),
        Line2D([0], [0], marker="s", linestyle="", color="w", markerfacecolor="#9ca3af", markeredgecolor=LAYER_EDGE["Methylation"], markersize=8, label="Methylation"),
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="none", markeredgecolor="#1f2937", markersize=10, label="Appears in a screen-leading signature"),
    ]
    legend_style = {
        "frameon": True,
        "fancybox": True,
        "framealpha": 0.96,
        "facecolor": "white",
        "edgecolor": "#e2e8f0",
    }
    fig.legend(
        handles=family_handles,
        fontsize=8.5,
        loc="lower left",
        bbox_to_anchor=(0.08, 0.045),
        ncol=3,
        title="Biological class",
        title_fontsize=8.7,
        **legend_style,
    )
    fig.legend(
        handles=layer_handles,
        fontsize=8.4,
        loc="lower right",
        bbox_to_anchor=(0.94, 0.038),
        ncol=1,
        title="Layer / top hit ring",
        title_fontsize=8.7,
        **legend_style,
    )
    fig.text(
        0.015,
        0.015,
        "Bubble size scales with total consensus strength across all six subset screens.",
        ha="left",
        va="top",
        fontsize=8.7,
        color=MUTED_TEXT,
    )
    fig.subplots_adjust(left=0.11, right=0.98, top=0.95, bottom=0.19)
    save_figure_bundle(fig, PLOTS_DIR / "Figure_01B_OS_temporal_scatter_v2.png")
    plt.close(fig)


def plot_os_single_gene_km_atlas(best_os: pd.DataFrame, gene_table: pd.DataFrame) -> None:
    selected = (
        best_os.sort_values(["p_value", "logrank_p", "consensus_strength"], ascending=[True, True, False])
        .head(8)
        .copy()
    )
    clinical = phase1.load_clinical(phase1.CLINICAL_PATH)
    genes = gene_table["gene"].tolist()
    rna = phase1.zscore_df(phase1.load_omics_matrix(phase1.RNA_PATH, genes))
    meth = phase1.zscore_df(phase1.load_omics_matrix(phase1.METH_PATH, genes))

    fig, axes = plt.subplots(2, 4, figsize=(16.8, 8.8))
    axes = axes.ravel()
    curve_rows: List[pd.DataFrame] = []

    for ax, (_, row) in zip(axes, selected.iterrows()):
        gene = str(row["gene"])
        modality = str(row["modality"])
        endpoint = str(row["endpoint"])
        family = str(row["Family6"])
        layer = str(row["Layer"])
        zmat = rna if modality == "RNA" else meth
        modality_clinical = clinical[clinical["Sample"].isin(set(zmat.index))].copy()
        endpoint_df = phase1.build_endpoint_df(modality_clinical, endpoint)
        endpoint_df = endpoint_df[endpoint_df["Sample"].isin(zmat.index)].copy()
        endpoint_df["score"] = endpoint_df["Sample"].map(zmat[gene])
        endpoint_df["time_months"] = endpoint_df["time_days"] / phase1.DAYS_PER_MONTH
        median_score = endpoint_df["score"].median()
        endpoint_df["group"] = np.where(endpoint_df["score"] >= median_score, "High", "Low")
        endpoint_df["gene"] = gene
        endpoint_df["modality"] = modality
        endpoint_df["endpoint"] = endpoint
        curve_rows.append(endpoint_df[["gene", "modality", "endpoint", "time_months", "event", "group", "score"]].copy())

        high = endpoint_df[endpoint_df["group"] == "High"].copy()
        low = endpoint_df[endpoint_df["group"] == "Low"].copy()
        high_color = FAMILY_COLORS.get(family, "#475569")
        border = FAMILY_COLORS.get(family, "#475569")

        kmf = KaplanMeierFitter()
        kmf.fit(high["time_months"], high["event"], label="High biomarker")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=high_color, linewidth=2.4)
        kmf.fit(low["time_months"], low["event"], label="Low biomarker")
        kmf.plot_survival_function(ax=ax, ci_show=False, color=KM_LOW_COLOR, linewidth=2.2)
        if ax.get_legend() is not None:
            ax.get_legend().remove()

        ax.set_title(
            f"{gene} [{layer_tag(layer)}] | {ENDPOINT_DISPLAY[endpoint]}",
            fontsize=11.2,
            fontweight="bold",
            pad=10,
            color=TEXT_COLOR,
        )
        ax.text(
            0.02,
            0.97,
            f"{family} | {row['target_label']}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.4,
            color=border,
            fontweight="bold",
        )
        ax.text(
            0.02,
            0.83,
            f"HR={row['hr']:.2f} | Cox P={format_p(float(row['p_value']))}\n"
            f"Adj.P={format_p(float(row['p_value_adj']))} | Log-rank P={format_p(float(row['logrank_p']))}\n"
            f"High/Low n={int(row['n_high'])}/{int(row['n_low'])}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=7.9,
            color=TEXT_COLOR,
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": lighten(border, 0.88),
                "edgecolor": lighten(border, 0.45),
                "linewidth": 0.8,
            },
        )
        ax.set_xlim(0, phase1.HORIZONS[endpoint] / phase1.DAYS_PER_MONTH if phase1.HORIZONS[endpoint] else 300)
        ax.set_ylim(0.55 if endpoint != "overall_os" else 0.50, 1.02)
        ax.grid(color=GRID_COLOR, lw=0.8)
        ax.set_xlabel("Months", fontsize=9.5, fontweight="bold")
        ax.set_ylabel("Survival probability", fontsize=9.5, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_color("#cbd5e1")

    for ax in axes[len(selected) :]:
        ax.axis("off")

    handles = [
        Line2D([0], [0], color="#111827", lw=2.4, label="High biomarker (panel color)"),
        Line2D([0], [0], color=KM_LOW_COLOR, lw=2.2, label="Low biomarker"),
    ]
    fig.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=2, fontsize=9.2)
    fig.suptitle(
        "C. Leading single-gene Kaplan–Meier profiles among OS-prioritized biomarkers",
        fontsize=18,
        fontweight="bold",
        y=0.992,
    )
    fig.text(
        0.5,
        0.972,
        "Genes are ranked by their best raw single-gene Cox signal on the modality in which each biomarker was selected.",
        ha="center",
        va="top",
        fontsize=9.4,
        color=MUTED_TEXT,
    )
    fig.suptitle(
        "C. Leading single-gene Kaplan-Meier profiles among OS-prioritized biomarkers",
        fontsize=18,
        fontweight="bold",
        y=0.992,
    )
    fig.subplots_adjust(top=0.84, bottom=0.08, left=0.07, right=0.985, wspace=0.24, hspace=0.33)
    save_figure_bundle(fig, PLOTS_DIR / "Figure_01C_OS_single_gene_km_atlas_v1.png")
    plt.close(fig)

    if curve_rows:
        pd.concat(curve_rows, ignore_index=True).to_csv(TABLES_DIR / "os_single_gene_km_curve_table.csv", index=False)


def plot_depmap_landscape_all(depmap_best: pd.DataFrame) -> None:
    available = depmap_best[depmap_best["depmap_available"]].copy()
    missing = depmap_best[~depmap_best["depmap_available"]].copy()
    fig = plt.figure(figsize=(11.8, 9.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[5.6, 1.45], hspace=0.18)
    ax = fig.add_subplot(gs[0, 0])
    strip = fig.add_subplot(gs[1, 0])
    header_bar(ax, "A. Biomarker-centered DepMap response landscape")

    if not available.empty:
        x_min = float(available["depmap_signed_score"].min()) - 0.45
        x_max = float(available["depmap_signed_score"].max()) + 0.35
        y_max = float(available["depmap_neglog10"].max()) + 0.22
    else:
        x_min, x_max, y_max = -1.0, 1.0, 1.0
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.02, y_max)
    ax.axvspan(x_min, 0, color="#eaf2ff", alpha=0.9, zorder=0)
    ax.axvspan(0, x_max, color="#fff0e3", alpha=0.9, zorder=0)
    ax.axvline(0, color="#64748b", lw=1.2, ls="--")
    ax.grid(color=GRID_COLOR, lw=0.9)
    ax.set_axisbelow(True)
    ax.text(0.015, 0.965, "More sensitive", transform=ax.transAxes, ha="left", va="top", fontsize=11, color=SENSITIVE_COLOR)
    ax.text(0.985, 0.965, "More resistant", transform=ax.transAxes, ha="right", va="top", fontsize=11, color=RESISTANT_COLOR)

    for _, row in available.iterrows():
        family = str(row["Family6"])
        layer = str(row["Layer"])
        tier = str(row["depmap_signif_tier"])
        face = FAMILY_COLORS.get(family, "#94a3b8")
        size = 140 + 560 * float(row["depmap_size_norm"])
        if tier == "FDR < 0.10":
            alpha = 0.95
            facecolor = face
            linewidth = 1.8
        elif tier == "FDR < 0.25":
            alpha = 0.78
            facecolor = lighten(face, 0.20)
            linewidth = 1.7
        else:
            alpha = 0.48
            facecolor = lighten(face, 0.60)
            linewidth = 1.4
        ax.scatter(
            float(row["depmap_signed_score"]),
            float(row["depmap_neglog10"]),
            s=size,
            c=[facecolor],
            marker=LAYER_MARKERS.get(layer, "o"),
            edgecolors=darken(face, 0.85),
            linewidths=linewidth,
            alpha=alpha,
            zorder=3,
        )

    label_df = available.copy()
    label_df["label_text"] = label_df["gene"]
    label_df["sub_text"] = label_df["best_drug"].fillna("") + " | " + label_df["best_target_class"].fillna("NA")
    annotate_scatter_points(
        ax,
        label_df,
        "depmap_signed_score",
        "depmap_neglog10",
        "label_text",
        sublabelcol="sub_text",
        sort_cols=["depmap_neglog10", "depmap_combined_score"],
    )

    ax.set_xlabel("Signed DepMap association", fontsize=12, fontweight="bold")
    ax.set_ylabel("-log10(best FDR across AUC / IC50)", fontsize=12, fontweight="bold")

    family_handles = [Patch(facecolor=FAMILY_COLORS[f], edgecolor=FAMILY_COLORS[f], label=f) for f in FAMILY_COLORS]
    layer_handles = [
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#94a3b8", markeredgecolor="#334155", markersize=8, label="Transcriptome"),
        Line2D([0], [0], marker="s", linestyle="", color="w", markerfacecolor="#94a3b8", markeredgecolor="#334155", markersize=8, label="Methylation"),
    ]
    tier_handles = [
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#64748b", markeredgecolor="#334155", alpha=0.95, markersize=8, label="FDR < 0.10"),
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#cbd5e1", markeredgecolor="#334155", alpha=0.78, markersize=8, label="FDR < 0.25"),
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#f8fafc", markeredgecolor="#94a3b8", alpha=0.48, markersize=8, label="Exploratory"),
    ]
    legend_style = {
        "frameon": True,
        "fancybox": True,
        "framealpha": 0.96,
        "facecolor": "white",
        "edgecolor": "#e2e8f0",
    }
    leg1 = ax.legend(handles=family_handles, fontsize=8.4, loc="lower left", bbox_to_anchor=(0.015, 0.02), ncol=2, title="Biomarker class", title_fontsize=8.6, **legend_style)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=layer_handles, fontsize=8.4, loc="lower center", bbox_to_anchor=(0.55, 0.02), ncol=2, title="Layer", title_fontsize=8.6, **legend_style)
    ax.add_artist(leg2)
    ax.legend(handles=tier_handles, fontsize=8.3, loc="lower right", bbox_to_anchor=(0.985, 0.02), ncol=1, title="Evidence tier", title_fontsize=8.6, **legend_style)

    strip_df = missing.sort_values(["Family6", "gene"]).copy()
    draw_chip_strip(strip, strip_df, "Biomarkers not carried into the tumor-intrinsic DepMap screen")
    fig.subplots_adjust(top=0.96, bottom=0.08, left=0.08, right=0.98)
    save_figure_bundle(fig, PLOTS_DIR / "Figure_02A_drug_depmap_landscape_all_biomarkers_v2.png")
    plt.close(fig)


def plot_depmap_dose_response_atlas(depmap_best: pd.DataFrame) -> None:
    dose = pd.read_csv(VALIDATION_TABLES / "top_gene_dose_response_curve_table.csv")
    available = depmap_best[depmap_best["depmap_available"]].copy()
    available = available.sort_values(["best_fdr", "depmap_combined_score"], ascending=[True, False]).reset_index(drop=True)

    n = available.shape[0]
    ncols = 5
    nrows = math.ceil(n / ncols) if n else 1
    fig, axes = plt.subplots(nrows, ncols, figsize=(18.0, 3.35 * nrows), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, row) in zip(axes, available.iterrows()):
        gene = str(row["gene"])
        drug = str(row["best_drug"])
        family = str(row["Family6"])
        border = FAMILY_COLORS.get(family, "#475569")
        direction = str(row["depmap_direction"])
        sub = dose[(dose["gene"] == gene) & (dose["drug_clean"] == drug)].copy()
        agg = (
            sub.groupby(["group", "log10_uM"], as_index=False)
            .agg(
                mean_viability=("mean_viability", "mean"),
                se_viability=("se_viability", "mean"),
                n=("n", "max"),
            )
            .sort_values(["group", "log10_uM"])
        )

        for group_name, color in [("high", border), ("low", KM_LOW_COLOR)]:
            part = agg[agg["group"] == group_name].copy()
            if part.empty:
                continue
            x = part["log10_uM"].to_numpy(dtype=float)
            y = part["mean_viability"].to_numpy(dtype=float)
            se = part["se_viability"].fillna(0).to_numpy(dtype=float)
            ax.plot(x, y, color=color, linewidth=2.2, label=group_name.title())
            ax.fill_between(x, y - se, y + se, color=color, alpha=0.15)

        ax.set_title(f"{gene} [{layer_tag(str(row['Layer']))}]", fontsize=11, fontweight="bold", color=TEXT_COLOR, pad=6)
        ax.text(0.02, 0.93, drug, transform=ax.transAxes, ha="left", va="top", fontsize=8.7, color=TEXT_COLOR, fontweight="bold")
        ax.text(0.02, 0.83, str(row["best_target_class"]), transform=ax.transAxes, ha="left", va="top", fontsize=7.9, color=MUTED_TEXT)
        badge_color = SENSITIVE_COLOR if direction == "More sensitive when high" else RESISTANT_COLOR
        badge_text = "High biomarker = sensitive" if direction == "More sensitive when high" else "High biomarker = resistant"
        ax.text(
            0.98,
            0.95,
            badge_text,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.7,
            color="white",
            bbox={
                "boxstyle": "round,pad=0.22",
                "facecolor": badge_color,
                "edgecolor": badge_color,
            },
        )
        ax.text(
            0.98,
            0.08,
            f"best FDR={format_p(float(row['best_fdr']))}",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.6,
            color=MUTED_TEXT,
        )
        ax.grid(color=GRID_COLOR, lw=0.75)
        ax.set_axisbelow(True)
        ax.set_ylim(0.18, 1.04)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_color(lighten(border, 0.45))
            ax.spines[spine].set_linewidth(1.1)

    for ax in axes[n:]:
        ax.axis("off")

    for ax in axes[-ncols:]:
        ax.set_xlabel("Drug dose (log10 uM)", fontsize=9.5, fontweight="bold")
    for ax in axes[::ncols]:
        ax.set_ylabel("Mean viability", fontsize=9.5, fontweight="bold")

    handles = [
        Line2D([0], [0], color="#111827", lw=2.2, label="High biomarker group (panel color)"),
        Line2D([0], [0], color=KM_LOW_COLOR, lw=2.2, label="Low biomarker group"),
    ]
    fig.legend(handles=handles, frameon=False, loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=2, fontsize=9.2)
    fig.suptitle(
        "B. Dose-response curves for biomarker–drug pairs prioritized by DepMap",
        fontsize=18,
        fontweight="bold",
        y=0.992,
    )
    fig.text(
        0.5,
        0.972,
        "All biomarkers with tumor-intrinsic DepMap hits are shown; each panel uses the strongest drug association for that gene.",
        ha="center",
        va="top",
        fontsize=9.4,
        color=MUTED_TEXT,
    )
    fig.suptitle(
        "B. Dose-response curves for biomarker-drug pairs prioritized by DepMap",
        fontsize=18,
        fontweight="bold",
        y=0.992,
    )
    fig.subplots_adjust(top=0.84, bottom=0.08, left=0.06, right=0.99, wspace=0.24, hspace=0.36)
    save_figure_bundle(fig, PLOTS_DIR / "Figure_02B_drug_dose_response_atlas_v1.png")
    plt.close(fig)


def plot_immune_landscape_all(immune_df: pd.DataFrame) -> None:
    available = immune_df[immune_df["immune_available"]].copy()
    missing = immune_df[~immune_df["immune_available"]].copy()
    fig = plt.figure(figsize=(11.8, 11.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[5.0, 2.35], hspace=0.14)
    ax = fig.add_subplot(gs[0, 0])
    strip = fig.add_subplot(gs[1, 0])
    header_bar(ax, "A. Biomarker-centered TIL/TMB immune landscape")

    x_min = float(available["rho_til"].min()) - 0.08 if not available.empty else -0.6
    x_max = float(available["rho_til"].max()) + 0.08 if not available.empty else 0.6
    y_min = float(available["rho_tmb"].min()) - 0.08 if not available.empty else -0.6
    y_max = float(available["rho_tmb"].max()) + 0.08 if not available.empty else 0.6
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    ax.add_patch(Rectangle((x_min, 0), 0 - x_min, y_max, facecolor="#eef3ff", edgecolor="none", alpha=0.78, zorder=0))
    ax.add_patch(Rectangle((0, 0), x_max, y_max, facecolor="#ecfdf5", edgecolor="none", alpha=0.78, zorder=0))
    ax.add_patch(Rectangle((x_min, y_min), 0 - x_min, 0 - y_min, facecolor="#f0fdf4", edgecolor="none", alpha=0.78, zorder=0))
    ax.add_patch(Rectangle((0, y_min), x_max, 0 - y_min, facecolor="#f8fafc", edgecolor="none", alpha=0.78, zorder=0))
    ax.axvline(0, color="#64748b", lw=1.2, ls="--")
    ax.axhline(0, color="#64748b", lw=1.2, ls="--")
    ax.grid(color=GRID_COLOR, lw=0.9)
    ax.set_axisbelow(True)

    for _, row in available.iterrows():
        family = str(row["Family6"])
        layer = str(row["Layer"])
        tier = str(row["immune_signif_tier"])
        face = FAMILY_COLORS.get(family, "#94a3b8")
        size = 140 + 700 * float(row["immune_size_norm"])
        if tier == "FDR < 0.10":
            facecolor = face
            alpha = 0.95
            linewidth = 1.8
        elif tier == "FDR < 0.25":
            facecolor = lighten(face, 0.22)
            alpha = 0.82
            linewidth = 1.7
        else:
            facecolor = lighten(face, 0.58)
            alpha = 0.58
            linewidth = 1.4
        ax.scatter(
            float(row["rho_til"]),
            float(row["rho_tmb"]),
            s=size,
            c=[facecolor],
            marker=LAYER_MARKERS.get(layer, "o"),
            edgecolors=darken(face, 0.88),
            linewidths=linewidth,
            alpha=alpha,
            zorder=3,
        )

    label_df = available.copy()
    label_df["label_text"] = label_df["gene"]
    label_df["sub_text"] = label_df["immune_quadrant"]
    immune_offsets = {
        "AURKB": (0.015, 0.03),
        "BCL2A1": (0.03, 0.025),
        "CLEC7A": (0.02, 0.02),
        "KLHDC7B": (0.018, 0.018),
        "MYO1G": (0.012, 0.02),
        "SLC44A4": (0.012, 0.014),
        "LRRC10B": (0.02, 0.02),
        "STMN3": (0.016, 0.012),
        "ZMYND10": (0.016, 0.016),
        "LRRC56": (0.016, 0.012),
        "TFF1": (0.018, 0.012),
        "MLPH": (0.018, 0.014),
        "MMP10": (0.014, 0.016),
        "CPLX1": (0.016, 0.018),
    }
    annotate_scatter_points(
        ax,
        label_df,
        "rho_til",
        "rho_tmb",
        "label_text",
        sublabelcol="sub_text",
        sort_cols=["immune_best_fdr", "immune_strength"],
        offset_map=immune_offsets,
        collision_x_frac=0.08,
        collision_y_frac=0.05,
        label_fontsize=8.8,
        sublabel_fontsize=7.4,
    )

    ax.set_xlabel("Spearman rho with TIL score", fontsize=12, fontweight="bold")
    ax.set_ylabel("Spearman rho with TMB", fontsize=12, fontweight="bold")
    ax.text(0.02, 0.96, "TMB-shifted", transform=ax.transAxes, ha="left", va="top", fontsize=11, color="#64748b")
    ax.text(0.98, 0.96, "Immune-hot / TMB-high", transform=ax.transAxes, ha="right", va="top", fontsize=11, color="#0f766e")
    ax.text(0.02, 0.05, "Immune-cold / TIL-TMB low", transform=ax.transAxes, ha="left", va="bottom", fontsize=10, color="#4d7c0f")
    ax.text(0.98, 0.05, "TIL-high / TMB-low", transform=ax.transAxes, ha="right", va="bottom", fontsize=10, color="#64748b")

    family_handles = [Patch(facecolor=FAMILY_COLORS[f], edgecolor=FAMILY_COLORS[f], label=f) for f in FAMILY_COLORS]
    layer_handles = [
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#94a3b8", markeredgecolor="#334155", markersize=8, label="Transcriptome"),
        Line2D([0], [0], marker="s", linestyle="", color="w", markerfacecolor="#94a3b8", markeredgecolor="#334155", markersize=8, label="Methylation"),
    ]
    tier_handles = [
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#64748b", markeredgecolor="#334155", alpha=0.95, markersize=8, label="FDR < 0.10"),
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#cbd5e1", markeredgecolor="#334155", alpha=0.82, markersize=8, label="FDR < 0.25"),
        Line2D([0], [0], marker="o", linestyle="", color="w", markerfacecolor="#f8fafc", markeredgecolor="#94a3b8", alpha=0.58, markersize=8, label="Exploratory"),
    ]
    legend_style = {
        "frameon": True,
        "fancybox": True,
        "framealpha": 0.96,
        "facecolor": "white",
        "edgecolor": "#e2e8f0",
    }
    leg1 = ax.legend(handles=family_handles, fontsize=8.4, loc="lower left", bbox_to_anchor=(0.015, 0.02), ncol=2, title="Biomarker class", title_fontsize=8.6, **legend_style)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=layer_handles, fontsize=8.4, loc="lower center", bbox_to_anchor=(0.56, 0.02), ncol=2, title="Layer", title_fontsize=8.6, **legend_style)
    ax.add_artist(leg2)
    ax.legend(handles=tier_handles, fontsize=8.3, loc="lower right", bbox_to_anchor=(0.985, 0.02), ncol=1, title="Evidence tier", title_fontsize=8.6, **legend_style)

    draw_chip_strip(strip, missing.sort_values(["Family6", "gene"]), "Biomarkers not carried into the TIL/TMB immune-correlation screen")
    fig.subplots_adjust(top=0.96, bottom=0.08, left=0.08, right=0.98)
    save_figure_bundle(fig, PLOTS_DIR / "Figure_03A_immune_tiltmb_landscape_all_biomarkers_v2.png")
    plt.close(fig)


def write_summary(
    gene_table: pd.DataFrame,
    depmap_best: pd.DataFrame,
    immune_df: pd.DataFrame,
    single_gene_os: pd.DataFrame,
) -> None:
    best_os = single_gene_os.sort_values(["p_value", "logrank_p"]).head(8)
    depmap_avail = depmap_best[depmap_best["depmap_available"]]
    immune_avail = immune_df[immune_df["immune_available"]]

    lines = [
        "# Biomarker Publication Suite v2",
        "",
        "## Figure set",
        "",
        "- `Figure_01A_OS_gene_consensus_matrix_v2`: gene-level OS consensus matrix using exhaustive subset screening across 6 OS screens.",
        "- `Figure_01B_OS_temporal_scatter_v2`: short-term versus long-term OS consensus behavior of each biomarker.",
        "- `Figure_01C_OS_single_gene_km_atlas_v1`: leading single-gene KM plots ranked by best raw single-gene Cox signal.",
        "- `Figure_02A_drug_depmap_landscape_all_biomarkers_v2`: all biomarkers in the tumor-intrinsic DepMap landscape, including a strip for biomarkers not carried into the screen.",
        "- `Figure_02B_drug_dose_response_atlas_v1`: dose-response curves for all biomarker genes with DepMap drug hits.",
        "- `Figure_03A_immune_tiltmb_landscape_all_biomarkers_v2`: biomarker-centered TIL/TMB landscape plus a strip for genes not carried into the immune screen.",
        "",
        "## Quick summary",
        "",
        f"- Total biomarker genes: {gene_table.shape[0]}",
        f"- Tumor-intrinsic DepMap hits available: {depmap_avail.shape[0]} genes",
        f"- TIL/TMB immune-correlation hits available: {immune_avail.shape[0]} genes",
        "",
        "## Leading single-gene OS profiles",
        "",
    ]
    for _, row in best_os.iterrows():
        lines.append(
            f"- `{row['gene']}` ({row['endpoint_label']}, {row['modality']}): HR={row['hr']:.2f}, Cox P={format_p(float(row['p_value']))}, "
            f"Adj.P={format_p(float(row['p_value_adj']))}, Log-rank P={format_p(float(row['logrank_p']))}"
        )

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- DepMap plot now retains all biomarker genes; genes without tumor-intrinsic validation are shown in the lower strip instead of being silently dropped.",
            "- Immune plot uses the same all-biomarker logic for TIL/TMB correlation availability.",
            "- Single-gene OS KM curves are ranked by raw Cox signal because age/stage-adjusted single-gene significance was limited in this cohort.",
        ]
    )
    (OUTPUT_ROOT / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ensure_dirs()
    setup_style()

    gene_table, metrics = load_core_gene_table()
    depmap_best, _ = build_depmap_all_biomarker_table(gene_table)
    immune_df = build_immune_all_biomarker_table(gene_table)
    single_gene_os = compute_single_gene_os_table(gene_table)

    plot_os_gene_consensus_matrix_v2(gene_table, metrics)
    plot_os_temporal_scatter_v2(gene_table)
    plot_os_single_gene_km_atlas(single_gene_os, gene_table)
    plot_depmap_landscape_all(depmap_best)
    plot_depmap_dose_response_atlas(depmap_best)
    plot_immune_landscape_all(immune_df)

    write_summary(gene_table, depmap_best, immune_df, single_gene_os)


if __name__ == "__main__":
    main()
