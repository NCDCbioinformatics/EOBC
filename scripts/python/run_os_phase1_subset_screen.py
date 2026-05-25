from __future__ import annotations

import argparse
import itertools
import math
import os
from pathlib import Path
from textwrap import fill
from typing import Dict, Iterable, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test


def env_path(name: str, default: Path | None = None) -> Path:
    value = os.environ.get(name)
    if value:
        return Path(value).expanduser()
    if default is not None:
        return Path(default).expanduser()
    raise RuntimeError(f"Set {name}; see config/paths_template.yml.")


BASE_DIR = env_path("EOBC_BIOMARKER_ROOT")
INPUT_DIR = BASE_DIR / "00_inputs_detected"
SELECTED_PANEL = (
    BASE_DIR
    / "09_selected_relaxed_union_family6_final_signedMeth"
    / "tables"
    / "selected_union_biomarker_genes_family6_v18.csv"
)
OUTPUT_ROOT = BASE_DIR / "final_analysis" / "os_phase1_subset_screen_v1"

RNA_PATH = INPUT_DIR / "TPM_young.csv"
METH_PATH = INPUT_DIR / "MET_young_batch_JW.csv"
CLINICAL_PATH = INPUT_DIR / "total_sample_clinical_all.csv"

HORIZONS = {
    "overall_os": None,
    "os_5y": 365.25 * 5,
    "os_10y": 365.25 * 10,
}

MODALITY_COLORS = {
    "RNA": "#0b6e4f",
    "METH": "#bf5b17",
}

ENDPOINT_COLORS = {
    "overall_os": "#264653",
    "os_5y": "#2a9d8f",
    "os_10y": "#e76f51",
}

SIZE_COLORS = {
    1: "#2a9d8f",
    2: "#e9c46a",
    3: "#f4a261",
    4: "#e76f51",
}

DAYS_PER_MONTH = 365.25 / 12.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Phase 1 EOBC prognostic subset screen with RNA/METH biomarker panels."
    )
    parser.add_argument(
        "--max-subset-size",
        type=int,
        default=3,
        help="Maximum gene-set size for exhaustive subset screening.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top subsets to retain per modality-endpoint condition.",
    )
    parser.add_argument(
        "--top-frequency-n",
        type=int,
        default=30,
        help="Number of top subsets used for gene-frequency heatmaps.",
    )
    return parser.parse_args()


def setup_plot_style() -> None:
    plt.rcParams["font.family"] = ["Malgun Gothic", "DejaVu Sans", "sans-serif"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["figure.facecolor"] = "#fcfbf7"
    plt.rcParams["axes.facecolor"] = "#fcfbf7"
    plt.rcParams["savefig.facecolor"] = "#fcfbf7"
    plt.rcParams["axes.edgecolor"] = "#6c757d"
    plt.rcParams["axes.titleweight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    sns.set_theme(
        style="whitegrid",
        rc={
            "axes.spines.right": False,
            "axes.spines.top": False,
            "grid.alpha": 0.15,
            "grid.color": "#5c677d",
        },
    )


def ensure_output_dirs(root: Path) -> Dict[str, Path]:
    tables = root / "tables"
    plots = root / "plots"
    km_dir = plots / "km_per_condition"
    logs = root / "logs"
    for p in (root, tables, plots, km_dir, logs):
        p.mkdir(parents=True, exist_ok=True)
    return {"root": root, "tables": tables, "plots": plots, "km": km_dir, "logs": logs}


def save_figure_bundle(fig: plt.Figure, output_path: Path, dpi: int = 320) -> None:
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")


def safe_stage_number(stage_value: str) -> float:
    if pd.isna(stage_value):
        return np.nan
    digits = "".join(ch for ch in str(stage_value) if ch.isdigit())
    return float(digits) if digits else np.nan


def clean_sample_columns(columns: Iterable[str]) -> List[str]:
    return [str(col).strip('"') for col in columns]


def load_panel(panel_path: Path) -> pd.DataFrame:
    panel = pd.read_csv(panel_path)
    panel["gene"] = panel["gene"].astype(str)
    return panel


def load_clinical(clinical_path: Path) -> pd.DataFrame:
    clinical = pd.read_csv(clinical_path)
    clinical = clinical[clinical["tpye"] == "Tumor"].copy()
    clinical["Sample"] = clinical["Row.names"].astype(str)
    clinical["overall_days"] = pd.to_numeric(clinical["overall_survival"], errors="coerce")
    clinical["death_event"] = clinical["alive"].map({"No": 1, "Yes": 0})
    clinical["recurrence_event"] = clinical["recurrence"].map({"Yes": 1, "No": 0})
    clinical["stage_num"] = clinical["Stage_final"].apply(safe_stage_number)
    clinical["age_num"] = pd.to_numeric(clinical["age"], errors="coerce")
    clinical["cohort_group"] = clinical["race"].fillna("Unknown").replace(
        {
            "Korean": "Korean",
            "White": "White",
            "Blackorafricanamerican": "Black",
        }
    )
    clinical = clinical.dropna(subset=["overall_days", "death_event"]).copy()
    return clinical


def load_omics_matrix(matrix_path: Path, selected_genes: List[str]) -> pd.DataFrame:
    raw = pd.read_csv(matrix_path, index_col=0)
    raw.index = raw.index.astype(str)
    raw.columns = clean_sample_columns(raw.columns)
    gene_order = [g for g in selected_genes if g in raw.index]
    matrix = raw.loc[gene_order].T
    matrix.index.name = "Sample"
    matrix = matrix.apply(pd.to_numeric, errors="coerce")
    tumor_cols = [c for c in matrix.index if c.endswith(".T")]
    matrix = matrix.loc[tumor_cols]
    return matrix


def zscore_df(df: pd.DataFrame) -> pd.DataFrame:
    centered = df - df.mean(axis=0)
    scale = df.std(axis=0, ddof=0).replace(0, np.nan)
    z = centered.divide(scale, axis=1)
    return z.fillna(0.0)


def build_endpoint_df(clinical: pd.DataFrame, endpoint_name: str) -> pd.DataFrame:
    endpoint = clinical[["Sample", "overall_days", "death_event", "age_num", "stage_num", "cohort_group"]].copy()
    horizon = HORIZONS[endpoint_name]
    if horizon is None:
        endpoint["time_days"] = endpoint["overall_days"]
        endpoint["event"] = endpoint["death_event"].astype(int)
    else:
        endpoint["time_days"] = endpoint["overall_days"].clip(upper=horizon)
        endpoint["event"] = ((endpoint["death_event"] == 1) & (endpoint["overall_days"] <= horizon)).astype(int)
    endpoint["endpoint"] = endpoint_name
    return endpoint


def fit_univariate_score_model(data: pd.DataFrame) -> Dict[str, float]:
    result = {
        "n": int(data.shape[0]),
        "events": int(data["event"].sum()),
        "beta": np.nan,
        "hr": np.nan,
        "lower95": np.nan,
        "upper95": np.nan,
        "p_value": np.nan,
        "c_index": np.nan,
    }
    if data.shape[0] < 20 or data["event"].sum() < 3:
        return result
    if np.isclose(data["score"].std(ddof=0), 0):
        return result
    cph = CoxPHFitter(penalizer=0.03)
    try:
        cph.fit(data[["time_days", "event", "score"]], duration_col="time_days", event_col="event")
        row = cph.summary.loc["score"]
        result.update(
            {
                "beta": float(row["coef"]),
                "hr": float(row["exp(coef)"]),
                "lower95": float(row["exp(coef) lower 95%"]),
                "upper95": float(row["exp(coef) upper 95%"]),
                "p_value": float(row["p"]),
                "c_index": float(cph.concordance_index_),
            }
        )
    except Exception:
        pass
    return result


def fit_adjusted_score_model(data: pd.DataFrame) -> Dict[str, float]:
    result = {
        "n_adj": int(data.shape[0]),
        "events_adj": int(data["event"].sum()),
        "beta_adj": np.nan,
        "hr_adj": np.nan,
        "lower95_adj": np.nan,
        "upper95_adj": np.nan,
        "p_value_adj": np.nan,
    }
    use = data.dropna(subset=["age_num", "stage_num"]).copy()
    if use.shape[0] < 30 or use["event"].sum() < 6:
        return result
    if np.isclose(use["score"].std(ddof=0), 0):
        return result
    cph = CoxPHFitter(penalizer=0.08)
    try:
        cph.fit(
            use[["time_days", "event", "score", "age_num", "stage_num"]],
            duration_col="time_days",
            event_col="event",
        )
        row = cph.summary.loc["score"]
        result.update(
            {
                "n_adj": int(use.shape[0]),
                "events_adj": int(use["event"].sum()),
                "beta_adj": float(row["coef"]),
                "hr_adj": float(row["exp(coef)"]),
                "lower95_adj": float(row["exp(coef) lower 95%"]),
                "upper95_adj": float(row["exp(coef) upper 95%"]),
                "p_value_adj": float(row["p"]),
            }
        )
    except Exception:
        pass
    return result


def derive_gene_signs(z_matrix: pd.DataFrame, endpoint_df: pd.DataFrame) -> Dict[str, float]:
    signs: Dict[str, float] = {}
    merged = endpoint_df.set_index("Sample").join(z_matrix, how="inner")
    for gene in z_matrix.columns:
        data = merged[["time_days", "event", gene]].dropna().rename(columns={gene: "score"})
        model = fit_univariate_score_model(data)
        beta = model["beta"]
        signs[gene] = 1.0 if np.isnan(beta) or np.isclose(beta, 0) else float(np.sign(beta))
    return signs


def screen_subsets_for_condition(
    z_matrix: pd.DataFrame,
    endpoint_df: pd.DataFrame,
    modality: str,
    endpoint_name: str,
    max_subset_size: int,
) -> pd.DataFrame:
    endpoint_join = endpoint_df.set_index("Sample").join(z_matrix, how="inner")
    gene_signs = derive_gene_signs(z_matrix.loc[endpoint_join.index], endpoint_df.loc[endpoint_df["Sample"].isin(endpoint_join.index)])
    genes = list(z_matrix.columns)
    x = endpoint_join[genes].to_numpy(dtype=float)
    time_days = endpoint_join["time_days"].to_numpy(dtype=float)
    event = endpoint_join["event"].to_numpy(dtype=int)
    sample_ids = endpoint_join.index.to_numpy()
    sign_array = np.array([gene_signs[g] for g in genes], dtype=float)

    rows: List[Dict[str, object]] = []
    for subset_size in range(1, max_subset_size + 1):
        combos = list(itertools.combinations(range(len(genes)), subset_size))
        print(f"[{modality} | {endpoint_name}] subset size {subset_size}: {len(combos)} combinations")
        for i, idx_tuple in enumerate(combos, start=1):
            idx = np.array(idx_tuple, dtype=int)
            score = (x[:, idx] * sign_array[idx]).mean(axis=1)
            data = pd.DataFrame(
                {
                    "Sample": sample_ids,
                    "time_days": time_days,
                    "event": event,
                    "score": score,
                    "age_num": endpoint_join["age_num"].to_numpy(dtype=float),
                    "stage_num": endpoint_join["stage_num"].to_numpy(dtype=float),
                }
            )
            fit = fit_univariate_score_model(data[["time_days", "event", "score"]])
            if np.isnan(fit["p_value"]):
                continue
            rows.append(
                {
                    "modality": modality,
                    "endpoint": endpoint_name,
                    "subset_size": subset_size,
                    "subset_label": " + ".join(genes[j] for j in idx),
                    "genes": "|".join(genes[j] for j in idx),
                    "n": fit["n"],
                    "events": fit["events"],
                    "beta": fit["beta"],
                    "hr": fit["hr"],
                    "lower95": fit["lower95"],
                    "upper95": fit["upper95"],
                    "p_value": fit["p_value"],
                    "neglog10_p": -math.log10(max(fit["p_value"], 1e-300)),
                    "c_index": fit["c_index"],
                    "gene_signs": "|".join(str(int(gene_signs[genes[j]])) for j in idx),
                }
            )
            if i % 500 == 0:
                print(f"  processed {i}/{len(combos)} combinations")
    result = pd.DataFrame(rows)
    if not result.empty:
        result = result.sort_values(
            by=["p_value", "c_index", "subset_size"],
            ascending=[True, False, True],
        ).reset_index(drop=True)
        result["rank"] = np.arange(1, len(result) + 1)
    return result


def add_adjusted_metrics(
    top_results: pd.DataFrame,
    z_matrix: pd.DataFrame,
    endpoint_df: pd.DataFrame,
) -> pd.DataFrame:
    endpoint_join = endpoint_df.set_index("Sample").join(z_matrix, how="inner")
    rows: List[Dict[str, float]] = []
    for _, row in top_results.iterrows():
        genes = row["genes"].split("|")
        signs = np.array([float(s) for s in row["gene_signs"].split("|")], dtype=float)
        score = (endpoint_join[genes].to_numpy(dtype=float) * signs).mean(axis=1)
        data = endpoint_join[["time_days", "event", "age_num", "stage_num"]].copy()
        data["score"] = score
        adj = fit_adjusted_score_model(data.reset_index(drop=True))
        rows.append(adj)
    adj_df = pd.DataFrame(rows)
    return pd.concat([top_results.reset_index(drop=True), adj_df], axis=1)


def compute_km_stats(
    z_matrix: pd.DataFrame,
    endpoint_df: pd.DataFrame,
    genes: List[str],
    signs: List[float],
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    joined = endpoint_df.set_index("Sample").join(z_matrix[genes], how="inner")
    score = (joined[genes].to_numpy(dtype=float) * np.array(signs, dtype=float)).mean(axis=1)
    km_df = joined[["time_days", "event"]].copy()
    km_df["time_months"] = km_df["time_days"] / DAYS_PER_MONTH
    km_df["score"] = score
    split = float(np.median(score))
    km_df["group"] = np.where(km_df["score"] >= split, "High score", "Low score")

    fit = fit_univariate_score_model(km_df[["time_days", "event", "score"]].copy())
    stats = {
        "split_value": split,
        "logrank_p": np.nan,
        "n_high": int((km_df["group"] == "High score").sum()),
        "n_low": int((km_df["group"] == "Low score").sum()),
        "hr": fit["hr"],
        "p_value": fit["p_value"],
    }
    try:
        high = km_df[km_df["group"] == "High score"]
        low = km_df[km_df["group"] == "Low score"]
        logrank = logrank_test(high["time_days"], low["time_days"], high["event"], low["event"])
        stats["logrank_p"] = float(logrank.p_value)
    except Exception:
        pass
    return km_df.reset_index(), stats


def plot_dataset_overview(summary_df: pd.DataFrame, output_path: Path) -> None:
    plot_df = summary_df.copy()
    plot_df["condition"] = plot_df["modality"] + "\n" + plot_df["endpoint"].str.replace("_", " ", regex=False)

    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(len(plot_df))
    width = 0.36
    ax.bar(
        x - width / 2,
        plot_df["n_samples"],
        width,
        color=[MODALITY_COLORS[m] for m in plot_df["modality"]],
        alpha=0.9,
        label="Samples",
    )
    ax.bar(
        x + width / 2,
        plot_df["events"],
        width,
        color=[ENDPOINT_COLORS[e] for e in plot_df["endpoint"]],
        alpha=0.9,
        label="Events",
    )
    for i, row in plot_df.iterrows():
        ax.text(i - width / 2, row["n_samples"] + 1.2, int(row["n_samples"]), ha="center", va="bottom", fontsize=9)
        ax.text(i + width / 2, row["events"] + 1.2, int(row["events"]), ha="center", va="bottom", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(plot_df["condition"], fontsize=10)
    ax.set_ylabel("Count")
    ax.set_title("EOBC Prognostic Screening Overview", fontsize=16, pad=12)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.set_axisbelow(True)
    fig.tight_layout()
    save_figure_bundle(fig, output_path, dpi=300)
    plt.close(fig)


def plot_top_sets_grid(top_results: pd.DataFrame, output_path: Path) -> None:
    modalities = ["RNA", "METH"]
    endpoints = ["overall_os", "os_5y", "os_10y"]
    fig, axes = plt.subplots(len(endpoints), len(modalities), figsize=(18, 16), sharex=False)

    for r, endpoint in enumerate(endpoints):
        for c, modality in enumerate(modalities):
            ax = axes[r, c]
            panel = top_results[(top_results["modality"] == modality) & (top_results["endpoint"] == endpoint)].copy()
            panel = panel.sort_values("p_value", ascending=False).tail(10)
            if panel.empty:
                ax.axis("off")
                continue
            y = np.arange(panel.shape[0])
            colors = [SIZE_COLORS.get(int(v), "#999999") for v in panel["subset_size"]]
            ax.hlines(y, panel["lower95"], panel["upper95"], color="#8d99ae", linewidth=2, alpha=0.9)
            ax.scatter(panel["hr"], y, s=90, c=colors, edgecolor="#1f2933", linewidth=0.7, zorder=3)
            ax.axvline(1.0, color="#6c757d", linestyle="--", linewidth=1.0)
            labels = [
                label if len(label) <= 36 else label[:33] + "..."
                for label in panel["subset_label"]
            ]
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=9)
            ax.set_xscale("log")
            ax.set_xlabel("Hazard ratio")
            ax.set_title(f"{modality} | {endpoint.replace('_', ' ')}", fontsize=13, pad=10)
            for yy, (_, row) in zip(y, panel.iterrows()):
                ax.text(
                    row["upper95"] * 1.05,
                    yy,
                    f"p={row['p_value']:.2e}\nC={row['c_index']:.3f}",
                    va="center",
                    fontsize=8,
                    color="#2f3e46",
                )
            ax.grid(axis="x", alpha=0.15)
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", label=f"Size {k}", markerfacecolor=v, markeredgecolor="#1f2933", markersize=8)
        for k, v in SIZE_COLORS.items()
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.0))
    fig.suptitle("Top Prognostic Gene Sets by Modality and Endpoint", fontsize=18, y=1.02)
    fig.tight_layout()
    save_figure_bundle(fig, output_path, dpi=320)
    plt.close(fig)


def plot_gene_frequency_heatmap(freq_df: pd.DataFrame, output_path: Path) -> None:
    heat = freq_df.pivot(index="gene", columns="condition", values="frequency").fillna(0.0)
    heat = heat.loc[heat.mean(axis=1).sort_values(ascending=False).index]

    fig, ax = plt.subplots(figsize=(10, max(8, 0.35 * heat.shape[0])))
    sns.heatmap(
        heat,
        cmap=sns.color_palette(["#f3ede2", "#dcefe4", "#5b8a72"], as_cmap=True),
        linewidths=0.4,
        linecolor="#fcfbf7",
        cbar_kws={"label": "Frequency among top subsets"},
        ax=ax,
    )
    ax.set_title("Gene Reuse Frequency in High-Ranking Subsets", fontsize=16, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    save_figure_bundle(fig, output_path, dpi=320)
    plt.close(fig)


def plot_km_panels(km_results: Dict[Tuple[str, str], Dict[str, object]], output_path: Path, km_dir: Path) -> None:
    endpoints = ["overall_os", "os_5y", "os_10y"]
    modalities = ["RNA", "METH"]
    fig, axes = plt.subplots(len(endpoints), len(modalities), figsize=(17, 15), sharey=True)

    for r, endpoint in enumerate(endpoints):
        for c, modality in enumerate(modalities):
            ax = axes[r, c]
            item = km_results.get((modality, endpoint))
            if not item:
                ax.axis("off")
                continue
            km_df = item["km_df"]
            stats = item["stats"]
            top_row = item["top_row"]
            for group_name, color in [("High score", MODALITY_COLORS[modality]), ("Low score", "#7a8c99")]:
                subset = km_df[km_df["group"] == group_name]
                kmf = KaplanMeierFitter()
                kmf.fit(subset["time_months"], subset["event"], label=group_name)
                kmf.plot_survival_function(ax=ax, ci_show=False, color=color, linewidth=2.4)
            ax.set_title(f"{modality} | {endpoint.replace('_', ' ')}", fontsize=13)
            ax.set_xlabel("Months")
            ax.set_ylabel("Survival probability")
            ax.set_ylim(0.0, 1.02)
            subset_label = fill(str(top_row["subset_label"]), width=28)
            subtitle = (
                f"{subset_label}\n"
                f"HR={stats['hr']:.2f} | Cox p={stats['p_value']:.2e} | Log-rank p={stats['logrank_p']:.2e}"
            )
            ax.text(
                0.03,
                0.96,
                subtitle,
                transform=ax.transAxes,
                fontsize=8.6,
                va="top",
                ha="left",
                bbox={"boxstyle": "round,pad=0.28", "facecolor": "#fcfbf7", "edgecolor": "#d9d9d9", "alpha": 0.95},
            )
            ax.legend(frameon=False, loc="lower left", fontsize=8.5)
            single_path = km_dir / f"KM_{modality}_{endpoint}.png"
            single_fig, single_ax = plt.subplots(figsize=(7.8, 5.8))
            for group_name, color in [("High score", MODALITY_COLORS[modality]), ("Low score", "#7a8c99")]:
                subset = km_df[km_df["group"] == group_name]
                kmf = KaplanMeierFitter()
                kmf.fit(subset["time_months"], subset["event"], label=group_name)
                kmf.plot_survival_function(ax=single_ax, ci_show=False, color=color, linewidth=2.5)
            single_ax.set_title(f"{modality} | {endpoint.replace('_', ' ')}", fontsize=13)
            single_ax.set_xlabel("Months")
            single_ax.set_ylabel("Survival probability")
            single_ax.set_ylim(0.0, 1.02)
            single_ax.text(
                0.03,
                0.96,
                subtitle,
                transform=single_ax.transAxes,
                fontsize=9,
                va="top",
                ha="left",
                bbox={"boxstyle": "round,pad=0.30", "facecolor": "#fcfbf7", "edgecolor": "#d9d9d9", "alpha": 0.95},
            )
            single_ax.legend(frameon=False, loc="lower left")
            single_fig.tight_layout()
            save_figure_bundle(single_fig, single_path, dpi=320)
            plt.close(single_fig)

    fig.suptitle("Kaplan-Meier Panels for the Highest-Ranking Gene Sets", fontsize=18, y=1.01)
    fig.tight_layout()
    save_figure_bundle(fig, output_path, dpi=320)
    plt.close(fig)


def write_summary_markdown(
    output_path: Path,
    top_results: pd.DataFrame,
    overview_df: pd.DataFrame,
    max_subset_size: int,
) -> None:
    lines = [
        "# EOBC OS Phase 1 Subset Screen",
        "",
        f"- Maximum subset size screened: `{max_subset_size}`",
        "- Screen design: RNA and methylation analyzed separately.",
        "- Endpoint definitions:",
        "  - `overall_os`: raw overall survival",
        "  - `os_5y`: administratively censored at 5 years",
        "  - `os_10y`: administratively censored at 10 years",
        "- Event mapping:",
        "  - death event: `alive == No`",
        "",
        "## Cohort summary",
        "",
    ]
    for _, row in overview_df.iterrows():
        lines.append(
            f"- `{row['modality']} | {row['endpoint']}`: n={int(row['n_samples'])}, events={int(row['events'])}"
        )
    lines.append("")
    lines.append("## Top set per condition")
    lines.append("")
    for modality in ["RNA", "METH"]:
        for endpoint in ["overall_os", "os_5y", "os_10y"]:
            panel = top_results[(top_results["modality"] == modality) & (top_results["endpoint"] == endpoint)]
            if panel.empty:
                continue
            row = panel.iloc[0]
            lines.extend(
                [
                    f"### {modality} | {endpoint}",
                    "",
                    f"- Subset: `{row['subset_label']}`",
                    f"- Size: `{int(row['subset_size'])}`",
                    f"- HR: `{row['hr']:.3f}` ({row['lower95']:.3f}-{row['upper95']:.3f})",
                    f"- Cox p: `{row['p_value']:.3e}`",
                    f"- C-index: `{row['c_index']:.3f}`",
                    (
                        f"- Adjusted HR (age + stage): `{row['hr_adj']:.3f}` "
                        f"({row['lower95_adj']:.3f}-{row['upper95_adj']:.3f}), p=`{row['p_value_adj']:.3e}`"
                        if not pd.isna(row["hr_adj"])
                        else "- Adjusted model: not stable / not enough complete cases"
                    ),
                    "",
                ]
            )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    setup_plot_style()
    out = ensure_output_dirs(OUTPUT_ROOT)

    panel = load_panel(SELECTED_PANEL)
    genes = panel["gene"].tolist()
    clinical = load_clinical(CLINICAL_PATH)
    rna = load_omics_matrix(RNA_PATH, genes)
    meth = load_omics_matrix(METH_PATH, genes)

    modality_mats = {
        "RNA": zscore_df(rna),
        "METH": zscore_df(meth),
    }

    overview_rows: List[Dict[str, object]] = []
    all_results: List[pd.DataFrame] = []
    top_results_with_adj: List[pd.DataFrame] = []
    km_payload: Dict[Tuple[str, str], Dict[str, object]] = {}
    freq_rows: List[Dict[str, object]] = []
    score_rows: List[pd.DataFrame] = []

    for modality, z_matrix in modality_mats.items():
        available_samples = set(z_matrix.index)
        modality_clinical = clinical[clinical["Sample"].isin(available_samples)].copy()
        for endpoint_name in HORIZONS:
            endpoint_df = build_endpoint_df(modality_clinical, endpoint_name)
            endpoint_df = endpoint_df[endpoint_df["Sample"].isin(z_matrix.index)].copy()
            overview_rows.append(
                {
                    "modality": modality,
                    "endpoint": endpoint_name,
                    "n_samples": int(endpoint_df.shape[0]),
                    "events": int(endpoint_df["event"].sum()),
                }
            )
            screen_path = out["tables"] / f"subset_screen_{modality}_{endpoint_name}.csv"
            if screen_path.exists():
                print(f"Using cached screen for {modality} | {endpoint_name}")
                screen_df = pd.read_csv(screen_path)
            else:
                print(f"Running screen for {modality} | {endpoint_name}: n={endpoint_df.shape[0]}, events={endpoint_df['event'].sum()}")
                screen_df = screen_subsets_for_condition(
                    z_matrix=z_matrix.loc[endpoint_df["Sample"]],
                    endpoint_df=endpoint_df,
                    modality=modality,
                    endpoint_name=endpoint_name,
                    max_subset_size=args.max_subset_size,
                )
                screen_df.to_csv(screen_path, index=False)
            all_results.append(screen_df)

            top_df = screen_df.head(args.top_n).copy()
            top_df = add_adjusted_metrics(top_df, z_matrix=z_matrix, endpoint_df=endpoint_df)
            top_path = out["tables"] / f"top_subsets_{modality}_{endpoint_name}.csv"
            top_df.to_csv(top_path, index=False)
            top_results_with_adj.append(top_df)

            freq_top = screen_df.head(args.top_frequency_n)
            for _, row in freq_top.iterrows():
                for gene in row["genes"].split("|"):
                    freq_rows.append(
                        {
                            "modality": modality,
                            "endpoint": endpoint_name,
                            "condition": f"{modality}\n{endpoint_name.replace('_', ' ')}",
                            "gene": gene,
                        }
                    )

            if not top_df.empty:
                top_row = top_df.iloc[0]
                top_genes = top_row["genes"].split("|")
                top_signs = [float(s) for s in top_row["gene_signs"].split("|")]
                km_df, km_stats = compute_km_stats(
                    z_matrix=z_matrix,
                    endpoint_df=endpoint_df,
                    genes=top_genes,
                    signs=top_signs,
                )
                km_payload[(modality, endpoint_name)] = {
                    "km_df": km_df,
                    "stats": km_stats,
                    "top_row": top_row,
                }
                km_table_path = out["tables"] / f"km_table_{modality}_{endpoint_name}.csv"
                km_df.to_csv(km_table_path, index=False)

                patient_scores = endpoint_df[["Sample", "time_days", "event"]].copy()
                patient_scores["time_months"] = patient_scores["time_days"] / DAYS_PER_MONTH
                patient_scores["score"] = km_df["score"].to_numpy()
                patient_scores["group"] = km_df["group"].to_numpy()
                patient_scores["modality"] = modality
                patient_scores["endpoint"] = endpoint_name
                patient_scores["subset_label"] = top_row["subset_label"]
                score_rows.append(patient_scores)

    overview_df = pd.DataFrame(overview_rows)
    overview_df.to_csv(out["tables"] / "cohort_overview.csv", index=False)

    all_results_df = pd.concat(all_results, ignore_index=True)
    all_results_df.to_csv(out["tables"] / "all_subset_screen_results_combined.csv", index=False)

    top_results_df = pd.concat(top_results_with_adj, ignore_index=True)
    top_results_df.to_csv(out["tables"] / "top_subset_summary_combined.csv", index=False)

    freq_df = (
        pd.DataFrame(freq_rows)
        .groupby(["modality", "endpoint", "condition", "gene"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    freq_df["frequency"] = freq_df["count"] / float(args.top_frequency_n)
    freq_df.to_csv(out["tables"] / "gene_frequency_top_subsets.csv", index=False)

    if score_rows:
        pd.concat(score_rows, ignore_index=True).to_csv(out["tables"] / "patient_scores_top_sets.csv", index=False)

    plot_dataset_overview(overview_df, out["plots"] / "Figure_01_dataset_overview.png")
    plot_top_sets_grid(top_results_df, out["plots"] / "Figure_02_top_sets_grid.png")
    plot_gene_frequency_heatmap(freq_df, out["plots"] / "Figure_03_gene_frequency_heatmap.png")
    plot_km_panels(km_payload, out["plots"] / "Figure_04_km_panels.png", out["km"])

    write_summary_markdown(
        out["root"] / "analysis_summary.md",
        top_results=top_results_df,
        overview_df=overview_df,
        max_subset_size=args.max_subset_size,
    )

    (out["logs"] / "run_info.txt").write_text(
        "\n".join(
            [
                "EOBC OS Phase 1 subset screen completed.",
                f"max_subset_size={args.max_subset_size}",
                f"top_n={args.top_n}",
                f"top_frequency_n={args.top_frequency_n}",
            ]
        ),
        encoding="utf-8",
    )
    print("Analysis complete.")


if __name__ == "__main__":
    main()
