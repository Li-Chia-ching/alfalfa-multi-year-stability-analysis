"""Generate Figures 1–7 from cohort-specific statistical reanalysis outputs.

Input files are expected in the directory supplied with --results-dir.
Figures are written to --output-dir in PDF and PNG formats at 600 dpi.
"""

import argparse
import json
import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

# ---------------- Logging & Optional Dependencies ----------------
logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

try:
    from adjustText import adjust_text

    HAS_ADJUST_TEXT = True
    logger.info("adjustText is available; enhanced label placement enabled.")
except ImportError:
    HAS_ADJUST_TEXT = False
    logger.warning(
        "adjustText is not installed; using simple label offsets.\n"
        "  Install with: pip install adjustText"
    )

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------- Configuration Defaults ----------------
DEFAULT_RESULTS_DIR = Path("04_Results/Revision_Final/Final_Statistical_Reanalysis")
DEFAULT_OUTPUT_DIR = Path("Figures")

# Color Palette (Okabe-Ito colorblind-friendly)
C_BLUE = "#0072B2"
C_ORANGE = "#D55E00"
C_AMBER = "#E69F00"
C_GREEN = "#009E73"
C_PINK = "#CC79A7"
C_SKY = "#56B4E9"
C_YELLOW = "#F0E442"
C_GREY = "#555555"
C_GRAY = "#9AA0A6"

HIGHLIGHT_COLORS: Dict[str, str] = {
    "Sitel": C_BLUE,
    "Xinjiang_Daye": C_ORANGE,
    "ACA542": C_AMBER,
    "L33": C_GREEN,
    "Algonquin": C_PINK,
}

COHORT_COLORS: Dict[str, str] = {
    "2001-established": C_BLUE,
    "2002-established": C_ORANGE,
}
COHORT_LABELS: Dict[str, str] = {
    "2001-established": "2001 cohort",
    "2002-established": "2002 cohort",
}
TRAIT_LABELS: Dict[str, str] = {
    "SummerPH": "Summer plant height (cm)",
    "FW": "Fresh weight (kg 5 m$^{-2}$)",
    "DW": "Dry weight (kg 5 m$^{-2}$)",
}
TRAIT_SHORT: Dict[str, str] = {
    "SummerPH": "Summer height",
    "FW": "Fresh weight",
    "DW": "Dry weight",
}
TRAIT_UNIT: Dict[str, str] = {
    "SummerPH": "cm",
    "FW": "kg 5 m$^{-2}$",
    "DW": "kg 5 m$^{-2}$",
}

# Plot rcParams setup
rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "bold",
    "axes.labelsize": 9,
    "axes.labelweight": "bold",
    "axes.linewidth": 1.0,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "legend.frameon": False,
    "figure.dpi": 100,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ---------------- Helper Functions ----------------
def despine(ax: plt.Axes, keep_left: bool = True, keep_bottom: bool = True) -> None:
    """Remove top and right spines from axis."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if not keep_left:
        ax.spines["left"].set_visible(False)
    if not keep_bottom:
        ax.spines["bottom"].set_visible(False)


def save_fig(fig: plt.Figure, name: str, output_dir: Path) -> None:
    """Save a figure as PDF and PNG format at specified resolution."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        path = output_dir / f"{name}.{ext}"
        fig.savefig(path, format=ext, dpi=600)
    plt.close(fig)
    logger.info("Saved figure: %s.pdf / .png", name)


def load_csv(results_dir: Path, name: str) -> pd.DataFrame:
    """Load a CSV file from the given results directory."""
    file_path = results_dir / name
    if not file_path.is_file():
        raise FileNotFoundError(f"Missing input dataset: {file_path}")
    return pd.read_csv(file_path)


def color_for(genotype: str) -> str:
    """Return highlight color for a given genotype."""
    return HIGHLIGHT_COLORS.get(genotype, C_GRAY)


def panel_letter(i: int) -> str:
    """Return an alphabetical panel label (A, B, C...)."""
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i]


def adjust_labels_with_leader_lines(
    ax: plt.Axes,
    texts: List[Any],
    points: List[Tuple[float, float]],
    force_text: Tuple[float, float] = (0.2, 0.25),
    force_points: Tuple[float, float] = (0.2, 0.25),
    expand: Tuple[float, float] = (1.05, 1.1),
    arrowprops: Optional[Dict[str, Any]] = None,
) -> None:
    """Adjust scatter labels, add leader lines, and expand axes to avoid overlap."""
    if not HAS_ADJUST_TEXT:
        for text, (x, y) in zip(texts, points):
            xlim, ylim = ax.get_xlim(), ax.get_ylim()
            x_range, y_range = xlim[1] - xlim[0], ylim[1] - ylim[0]
            offset_x = 0.05 * x_range if x < xlim[0] + 0.1 * x_range else (
                -0.05 * x_range if x > xlim[1] - 0.1 * x_range else 0.02 * x_range
            )
            offset_y = -0.05 * y_range if y > ylim[1] - 0.1 * y_range else 0.02 * y_range
            text.set_position((x + offset_x, y + offset_y))
    else:
        if arrowprops is None:
            arrowprops = dict(arrowstyle="-", color="black", lw=0.5, alpha=0.6)

        adjust_text(
            texts,
            x=[p[0] for p in points],
            y=[p[1] for p in points],
            ax=ax,
            force_text=force_text,
            force_points=force_points,
            expand=expand,
            arrowprops=arrowprops,
            only_move={"text": "xy"},
            avoid_self=True,
            avoid_points=True,
            avoid_axes=True,
            max_iter=800,
            tol=0.005,
        )

    try:
        fig = ax.get_figure()
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        for text in texts:
            bbox = text.get_window_extent(renderer=renderer)
            bbox_data = bbox.transformed(ax.transData.inverted())
            xmin = min(xmin, bbox_data.xmin)
            xmax = max(xmax, bbox_data.xmax)
            ymin = min(ymin, bbox_data.ymin)
            ymax = max(ymax, bbox_data.ymax)

        x_margin = (xmax - xmin) * 0.08
        y_margin = (ymax - ymin) * 0.08
        ax.set_xlim(xmin - x_margin, xmax + x_margin)
        ax.set_ylim(ymin - y_margin, ymax + y_margin)
    except Exception as exc:
        logger.debug("Axis expansion failed: %s", exc)


# ---------------- Matrix Decompositions (AMMI / GGE) ----------------
def _genotype_year_matrix(results_dir: Path, trait: str, cohort: str = "2001-established") -> pd.DataFrame:
    py = load_csv(results_dir, "plot_year_wide.csv")
    py = py[py["cohort"] == cohort]
    mat = py.groupby(["genotype_ID", "calendar_year"])[trait].mean().unstack().dropna()
    return mat


def _ammi_scores(mat: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray, pd.Index, pd.Index]:
    grand = mat.values.mean()
    row_eff = mat.values.mean(axis=1, keepdims=True) - grand
    col_eff = mat.values.mean(axis=0, keepdims=True) - grand
    resid = mat.values - grand - row_eff - col_eff
    U, s, Vt = np.linalg.svd(resid, full_matrices=False)
    g_scores = U * s
    y_scores = Vt.T
    pc_pct = (s ** 2) / np.sum(resid ** 2) * 100
    return g_scores, y_scores, pc_pct, mat.index, mat.columns


def _gge_scores(mat: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, pd.Index, pd.Index]:
    col_mean = mat.values.mean(axis=0, keepdims=True)
    centered = mat.values - col_mean
    U, s, Vt = np.linalg.svd(centered, full_matrices=False)
    g_scores = U * s
    y_scores = Vt.T
    pc_pct = (s ** 2) / np.sum(centered ** 2) * 100
    g_means = mat.values.mean(axis=1)
    return g_scores, y_scores, pc_pct, g_means, mat.index, mat.columns


# ---------------- Figure Generators ----------------
def fig1_same_stand_age(results_dir: Path, output_dir: Path) -> None:
    """Figure 1: Same-stand-age sensitivity analysis."""
    logger.info("[Fig 1] Generating same-stand-age sensitivity plot...")
    df = load_csv(results_dir, "same_stand_age_sensitivity_final.csv")
    traits = ["FW", "DW", "SummerPH"]
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.2))
    stand_ages = sorted(df["stand_age"].unique())
    x = np.arange(len(stand_ages))
    width = 0.38

    for ax, trait in zip(axes, traits):
        for i, cohort in enumerate(["2001-established", "2002-established"]):
            sub = df[(df["cohort"] == cohort) & (df["stand_age"].isin(stand_ages))].sort_values("stand_age")
            vals = sub[trait].values
            ax.bar(
                x + (i - 0.5) * width,
                vals,
                width,
                color=COHORT_COLORS[cohort],
                label=COHORT_LABELS[cohort],
                edgecolor="black",
                linewidth=0.5,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([f"Stand age {a}" for a in stand_ages])
        ax.set_ylabel(TRAIT_LABELS[trait])
        ax.set_title(f"({'abc'[list(traits).index(trait)]}) {TRAIT_SHORT[trait]}")
        despine(ax)
        ax.tick_params(direction="out", length=3)

    axes[0].legend(loc="upper left", handlelength=1.2, handletextpad=0.5)
    fig.tight_layout()
    save_fig(fig, "Fig1_same_stand_age_sensitivity", output_dir)


def fig2_variance_heritability(results_dir: Path, output_dir: Path) -> None:
    """Figure 2: Variance component partitioning and broad-sense heritability."""
    logger.info("[Fig 2] Generating variance & heritability plot...")
    vc = load_csv(results_dir, "variance_components_final_by_cohort.csv")
    h2 = load_csv(results_dir, "heritability_final_by_cohort.csv")
    traits = ["SummerPH", "FW", "DW"]
    cohorts = ["2001-established", "2002-established"]
    components = [
        ("V_G", "Genotype", C_BLUE),
        ("V_GY", "G x Year", C_ORANGE),
        ("V_plot", "Plot", C_AMBER),
        ("V_e", "Residual", C_GREY),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), gridspec_kw={"width_ratios": [1.5, 1]})
    
    # Panel A: Variance components
    ax = axes[0]
    x = np.arange(len(traits))
    width = 0.38
    for ti, trait in enumerate(traits):
        for ci, cohort in enumerate(cohorts):
            row = vc[(vc["cohort"] == cohort) & (vc["trait"] == trait)]
            if row.empty:
                continue
            bottom = 0.0
            for col, _, color in components:
                val = float(row[col].iloc[0])
                if val > 0:
                    ax.bar(
                        x[ti] + (ci - 0.5) * width,
                        val,
                        width,
                        bottom=bottom,
                        color=color,
                        edgecolor="black",
                        linewidth=0.4,
                    )
                    bottom += val
            if cohort == "2002-established" and trait == "FW":
                ax.scatter(
                    x[ti] + (ci - 0.5) * width,
                    0.05,
                    marker="*",
                    s=60,
                    color=C_PINK,
                    zorder=5,
                    edgecolor="black",
                    linewidth=0.5,
                )
    ax.set_xticks(x)
    ax.set_xticklabels([TRAIT_SHORT[t] for t in traits])
    ax.set_ylabel("Variance component")
    ax.set_title("(A) Variance partitioning by cohort")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, ec="black", lw=0.5) for _, _, c in components]
    handles.append(plt.scatter([], [], marker="*", s=60, color=C_PINK, edgecolor="black", linewidth=0.5))
    ax.legend(
        handles,
        [n for _, n, _ in components] + ["V_G at boundary"],
        loc="upper right",
        ncol=2,
        handlelength=1.2,
        handletextpad=0.4,
        columnspacing=1.0,
    )
    despine(ax)
    ax.tick_params(direction="out", length=3)

    # Panel B: Heritability
    ax = axes[1]
    for ci, cohort in enumerate(cohorts):
        vals = []
        for trait in traits:
            r = h2[(h2["cohort"] == cohort) & (h2["trait"] == trait)]
            v = float(r["H2"].iloc[0]) if not r.empty else np.nan
            vals.append(0.0 if v < 1e-6 else v)
        offset = (ci - 0.5) * width
        ax.bar(
            np.arange(len(traits)) + offset,
            vals,
            width,
            color=COHORT_COLORS[cohort],
            label=COHORT_LABELS[cohort],
            edgecolor="black",
            linewidth=0.5,
        )
        for bi, v in enumerate(vals):
            if cohort == "2002-established" and traits[bi] == "FW":
                ax.text(bi + offset, 0.02, "boundary", ha="center", va="bottom", fontsize=6, color=C_PINK, rotation=90)
    ax.set_xticks(np.arange(len(traits)))
    ax.set_xticklabels([TRAIT_SHORT[t] for t in traits])
    ax.set_ylabel("Broad-sense H$^2$ (genotype-mean)")
    ax.set_ylim(0, 1.0)
    ax.set_title("(B) Heritability by cohort")
    ax.legend(loc="upper right", handlelength=1.2)
    despine(ax)
    ax.tick_params(direction="out", length=3)

    fig.tight_layout()
    save_fig(fig, "Fig2_variance_heritability", output_dir)


def fig3_blup_ranking(results_dir: Path, output_dir: Path) -> None:
    """Figure 3: Cohort-specific BLUP caterpillar ranking plots."""
    logger.info("[Fig 3] Generating BLUP ranking plot...")
    blup = load_csv(results_dir, "BLUP_final_by_cohort.csv")

    meta = {
        "input_file": "BLUP_final_by_cohort.csv",
        "input_columns": list(blup.columns),
        "n_rows": int(len(blup)),
        "counts": {
            f"{c[:4]}_{t}": len(blup[(blup["cohort"] == c) & (blup["trait"] == t)])
            for c in ["2001-established", "2002-established"]
            for t in ["SummerPH", "FW", "DW"]
        },
        "identities": {
            "Xinjiang_Daye_exists": "Xinjiang_Daye" in blup["genotype"].unique(),
            "ACA542_exists": "ACA542" in blup["genotype"].unique(),
        },
        "uncertainty_intervals": "NONE",
    }
    (output_dir / "Fig3_input_metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False))

    fig, axes = plt.subplots(2, 3, figsize=(12, 13.2))
    panel_idx = 0

    for ri, cohort in enumerate(["2001-established", "2002-established"]):
        for ci, trait in enumerate(["SummerPH", "FW", "DW"]):
            ax = axes[ri, ci]
            sub = blup[(blup["cohort"] == cohort) & (blup["trait"] == trait)].sort_values("BLUP", ascending=False).reset_index(drop=True)
            n = len(sub)
            y = np.arange(n)
            blups = sub["BLUP"].values
            labels = sub["genotype"].values
            colors = [color_for(g) for g in labels]

            ax.scatter(blups, y, s=26, c=colors, edgecolors="black", linewidths=0.35, zorder=3)
            bmin, bmax = blups.min(), blups.max()
            span = bmax - bmin

            if trait == "FW" and cohort == "2002-established":
                m = float(np.mean(blups))
                ax.set_xlim(m - 0.35, m + 0.35)
                ax.axvline(m, color=C_GRAY, lw=0.8, ls="--", zorder=1, alpha=0.7)
                ax.text(
                    m + 0.36 - 0.7 * 0.02,
                    n - 0.5,
                    "V$_G$ $\\approx$ 0;\ngenotype\ndiscrimination\nlimited",
                    fontsize=7,
                    color=C_PINK,
                    va="top",
                    ha="right",
                    fontweight="bold",
                )
            else:
                margin = max(span * 0.06, 0.001)
                ax.set_xlim(bmin - margin, bmax + margin)

            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7.5)
            ax.set_ylim(n - 0.5, -0.5)
            ax.set_title(f"({panel_letter(panel_idx)})  {COHORT_LABELS[cohort]} — {TRAIT_SHORT[trait]}", loc="left")
            ax.set_xlabel(f"BLUP ({TRAIT_UNIT[trait]})")
            ax.xaxis.grid(True, which="major", color="#E0E0E0", lw=0.6, zorder=0)
            ax.yaxis.grid(False)
            ax.set_axisbelow(True)
            despine(ax)
            ax.tick_params(axis="y", length=0)
            ax.tick_params(axis="x", direction="out", length=3)
            panel_idx += 1

    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=C_GRAY, markeredgecolor="black", markersize=7, label="Other genotypes")
    ]
    for genotype, color in HIGHLIGHT_COLORS.items():
        handles.append(
            plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markeredgecolor="black", markersize=7, label=genotype)
        )

    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=6,
        bbox_to_anchor=(0.5, 1.005),
        frameon=False,
        handletextpad=0.4,
        columnspacing=1.3,
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.subplots_adjust(hspace=0.32, wspace=0.42)

    output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(output_dir / f"Fig3_BLUP_Caterpillar_v3.{ext}", format=ext, dpi=600)

    try:
        fig.savefig(output_dir / "Fig3_BLUP_Caterpillar_v3.tif", format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    except Exception as exc:
        logger.warning("TIFF export skipped: %s", exc)
    plt.close(fig)


def fig4_ammi_biplot(results_dir: Path, output_dir: Path) -> None:
    """Figure 4: AMMI biplots for temporal stability."""
    logger.info("[Fig 4] Generating AMMI biplot...")
    traits = ["SummerPH", "FW", "DW"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))

    for ax, trait in zip(axes, traits):
        mat = _genotype_year_matrix(results_dir, trait)
        g, y, pct, gnames, ynames = _ammi_scores(mat)

        scale = np.sqrt(np.sum(g[:, 0] ** 2) / np.sum(y[:, 0] ** 2)) if np.sum(y[:, 0] ** 2) > 0 else 1.0
        gx, gy = g[:, 0], g[:, 1] if g.shape[1] > 1 else np.zeros_like(g[:, 0])
        yx, yy = y[:, 0] * scale, (y[:, 1] * scale if y.shape[1] > 1 else np.zeros_like(y[:, 0]))

        ax.scatter(gx, gy, color=C_BLUE, s=18, edgecolor="black", linewidth=0.4, zorder=3)
        ax.scatter(yx, yy, color=C_ORANGE, marker="^", s=40, edgecolor="black", linewidth=0.4, zorder=4)

        texts, points = [], []
        top_idx = np.argsort(-mat.mean(axis=1).values)[:5]
        for i in top_idx:
            t = ax.annotate(gnames[i], (gx[i], gy[i]), fontsize=6, color=C_BLUE, fontweight="bold", xytext=(2, 2), textcoords="offset points")
            texts.append(t)
            points.append((gx[i], gy[i]))

        for j, yn in enumerate(ynames):
            t = ax.annotate(str(yn), (yx[j], yy[j]), fontsize=7, color=C_ORANGE, fontweight="bold", xytext=(3, 3), textcoords="offset points")
            texts.append(t)
            points.append((yx[j], yy[j]))

        ax.margins(0.2)
        adjust_labels_with_leader_lines(ax, texts, points, force_text=(0.15, 0.2), force_points=(0.15, 0.2), expand=(1.05, 1.1))

        ax.axhline(0, color="grey", lw=0.5, ls="--")
        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.set_xlabel(f"IPCA1 ({pct[0]:.1f}%)")
        if g.shape[1] > 1:
            ax.set_ylabel(f"IPCA2 ({pct[1]:.1f}%)")
        ax.set_title(TRAIT_SHORT[trait])
        despine(ax)
        ax.tick_params(direction="out", length=3)

    fig.suptitle("AMMI biplots (2001 cohort; temporal stability at a single site)", y=1.02, fontsize=9, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "Fig4_ammi_biplot", output_dir)


def fig5_gge_biplot(results_dir: Path, output_dir: Path) -> None:
    """Figure 5: GGE biplots for temporal stability."""
    logger.info("[Fig 5] Generating GGE biplot...")
    traits = ["SummerPH", "FW", "DW"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))

    for ax, trait in zip(axes, traits):
        mat = _genotype_year_matrix(results_dir, trait)
        g, y, pct, g_means, gnames, ynames = _gge_scores(mat)

        scale = np.sqrt(np.sum(g[:, 0] ** 2) / np.sum(y[:, 0] ** 2)) if np.sum(y[:, 0] ** 2) > 0 else 1.0
        gx, gy = g[:, 0], g[:, 1] if g.shape[1] > 1 else np.zeros_like(g[:, 0])
        yx, yy = y[:, 0] * scale, (y[:, 1] * scale if y.shape[1] > 1 else np.zeros_like(y[:, 0]))

        g_norm = (g_means - g_means.min()) / (g_means.max() - g_means.min() + 1e-9)
        sizes = 15 + g_norm * 35
        ax.scatter(gx, gy, color=C_BLUE, s=sizes, edgecolor="black", linewidth=0.4, zorder=3, alpha=0.85)
        ax.scatter(yx, yy, color=C_ORANGE, marker="^", s=50, edgecolor="black", linewidth=0.4, zorder=4)

        texts, points = [], []
        top_idx = np.argsort(-g_means)[:5]
        for i in top_idx:
            t = ax.annotate(gnames[i], (gx[i], gy[i]), fontsize=6, color=C_BLUE, fontweight="bold", xytext=(2, 2), textcoords="offset points")
            texts.append(t)
            points.append((gx[i], gy[i]))

        for j, yn in enumerate(ynames):
            t = ax.annotate(str(yn), (yx[j], yy[j]), fontsize=7, color=C_ORANGE, fontweight="bold", xytext=(3, 3), textcoords="offset points")
            texts.append(t)
            points.append((yx[j], yy[j]))

        ax.margins(0.2)
        adjust_labels_with_leader_lines(ax, texts, points, force_text=(0.15, 0.2), force_points=(0.15, 0.2), expand=(1.05, 1.1))

        ax.axhline(0, color="grey", lw=0.5, ls="--")
        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.set_xlabel(f"PC1 ({pct[0]:.1f}%)")
        if g.shape[1] > 1:
            ax.set_ylabel(f"PC2 ({pct[1]:.1f}%)")
        ax.set_title(TRAIT_SHORT[trait])
        despine(ax)
        ax.tick_params(direction="out", length=3)

    fig.suptitle("GGE biplots (2001 cohort; temporal stability at a single site)", y=1.02, fontsize=9, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "Fig5_gge_biplot", output_dir)


def fig6_mean_stability(results_dir: Path, output_dir: Path) -> None:
    """Figure 6: Mean performance vs AMMI stability value (ASV)."""
    logger.info("[Fig 6] Generating mean vs stability plot...")
    traits = ["SummerPH", "FW", "DW"]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))

    for ax, trait in zip(axes, traits):
        mat = _genotype_year_matrix(results_dir, trait)
        g, _, _, gnames, _ = _ammi_scores(mat)
        ss1 = np.sum(g[:, 0] ** 2)
        ss2 = np.sum(g[:, 1] ** 2) if g.shape[1] > 1 else 0
        asv = np.sqrt(ss1 / ss2) * np.abs(g[:, 0]) + np.abs(g[:, 1]) if ss2 > 0 else np.abs(g[:, 0])
        g_means = mat.mean(axis=1).values

        ax.scatter(asv, g_means, color=C_BLUE, s=22, edgecolor="black", linewidth=0.4, zorder=3)

        texts, points = [], []
        top_idx = np.argsort(-g_means)[:5]
        for i in top_idx:
            t = ax.annotate(gnames[i], (asv[i], g_means[i]), fontsize=6, color=C_BLUE, fontweight="bold", xytext=(3, 3), textcoords="offset points")
            texts.append(t)
            points.append((asv[i], g_means[i]))

        ax.margins(0.2)
        adjust_labels_with_leader_lines(ax, texts, points, force_text=(0.15, 0.2), force_points=(0.15, 0.2), expand=(1.05, 1.1))

        ax.set_xlabel("AMMI Stability Value (lower = more stable)")
        ax.set_ylabel(TRAIT_LABELS[trait])
        ax.set_title(TRAIT_SHORT[trait])
        despine(ax)
        ax.tick_params(direction="out", length=3)

    fig.suptitle("Mean performance vs temporal stability (2001 cohort)", y=1.03, fontsize=9, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "Fig6_mean_stability", output_dir)


def fig7_blup_correlation(results_dir: Path, output_dir: Path) -> None:
    """Figure 7: BLUP correlation heatmaps by cohort."""
    logger.info("[Fig 7] Generating BLUP correlation plot...")
    corr = load_csv(results_dir, "BLUP_correlation_by_cohort.csv")
    cohorts = ["2001-established", "2002-established"]
    trait_pairs = [("SummerPH", "FW"), ("SummerPH", "DW"), ("FW", "DW")]
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6))

    for ax, cohort in zip(axes, cohorts):
        mat = np.full((3, 3), np.nan)
        labels = ["SummerPH", "FW", "DW"]
        for i, (t1, t2) in enumerate(trait_pairs):
            r = corr[(corr["cohort"] == cohort) & (corr["trait1"] == t1) & (corr["trait2"] == t2)]
            if not r.empty:
                val = float(r["r"].iloc[0])
                mat[i, 2 - labels.index(t2)] = val
                mat[2 - labels.index(t2), i] = val

        for k in range(3):
            mat[k, k] = 1.0

        im = ax.imshow(mat, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
        ax.set_xticks(range(3))
        ax.set_yticks(range(3))
        ax.set_xticklabels(["Summer\nPH", "FW", "DW"])
        ax.set_yticklabels(["Summer PH", "FW", "DW"])
        ax.set_title(COHORT_LABELS[cohort])

        for ii in range(3):
            for jj in range(3):
                v = mat[ii, jj]
                if not np.isnan(v):
                    color = "white" if abs(v) > 0.6 else "black"
                    ax.text(jj, ii, f"{v:.2f}", ha="center", va="center", color=color, fontsize=8, fontweight="bold")

        ax.tick_params(direction="out", length=3)

    cbar_ax = fig.add_axes([0.93, 0.18, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    cbar.ax.tick_params(direction="out", length=3)
    fig.subplots_adjust(left=0.07, right=0.9, wspace=0.35)
    save_fig(fig, "Fig7_blup_correlation_cohort", output_dir)


# ---------------- CLI Entry Point ----------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate Figures 1-7 from cohort-specific statistical outputs.")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Directory containing the statistical result CSV files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for generated figures and metadata.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    results_dir: Path = args.results_dir
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    required_files = [
        "same_stand_age_sensitivity_final.csv",
        "variance_components_final_by_cohort.csv",
        "heritability_final_by_cohort.csv",
        "BLUP_final_by_cohort.csv",
        "plot_year_wide.csv",
        "BLUP_correlation_by_cohort.csv",
    ]
    missing = [name for name in required_files if not (results_dir / name).is_file()]
    if missing:
        missing_msg = "Required input files are missing:\n" + "\n".join(f"  - {name}" for name in missing)
        logger.error(missing_msg)
        raise FileNotFoundError(missing_msg)

    logger.info("Results directory: %s", results_dir.resolve())
    logger.info("Output directory:  %s", output_dir.resolve())

    figures = [
        ("Fig1", fig1_same_stand_age),
        ("Fig2", fig2_variance_heritability),
        ("Fig3", fig3_blup_ranking),
        ("Fig4", fig4_ammi_biplot),
        ("Fig5", fig5_gge_biplot),
        ("Fig6", fig6_mean_stability),
        ("Fig7", fig7_blup_correlation),
    ]

    for fig_name, fig_func in figures:
        try:
            fig_func(results_dir, output_dir)
        except Exception as exc:
            logger.error("Failed to generate %s: %s", fig_name, exc, exc_info=True)

    logger.info("All processing completed.")


if __name__ == "__main__":
    main()
