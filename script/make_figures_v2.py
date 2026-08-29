"""
Regenerate all 7 manuscript figures from the final cohort-specific results.

Inputs (read from 04_Results/Revision_Final/Final_Statistical_Reanalysis/):
  - BLUP_final_by_cohort.csv         (BLUP values per cohort/trait/genotype)
  - variance_components_final_by_cohort.csv  (V_G, V_GY, V_plot, V_e)
  - heritability_final_by_cohort.csv (H2 per cohort/trait)
  - same_stand_age_sensitivity_final.csv (cohort means by stand age)
  - BLUP_correlation_by_cohort.csv   (Pearson r per cohort)
  - plot_year_wide.csv               (analytical input for AMMI/GGE)

Outputs (PDF + PNG, 600 dpi) to 06_Manuscript/Formal_To_NZCHS/Figures_v2/:
  - Fig1_same_stand_age_sensitivity.{pdf,png}
  - Fig2_variance_heritability.{pdf,png}
  - Fig3_blup_ranking_cohort.{pdf,png}
  - Fig4_ammi_biplot.{pdf,png}
  - Fig5_gge_biplot.{pdf,png}
  - Fig6_mean_stability.{pdf,png}
  - Fig7_blup_correlation_cohort.{pdf,png}

Style: Nature/Springer (no gridlines, bold axes, minimal titles),
Okabe-Ito colorblind-friendly palette, Arial font, despine.
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------- Paths ----------------
# ANONYMIZED PATH
PROJECT_ROOT = Path("./")
RESULTS_DIR = PROJECT_ROOT / "04_Results" / "Revision_Final" / "Final_Statistical_Reanalysis"
FIG_OUT_DIR = PROJECT_ROOT / "06_Manuscript" / "Formal_To_NZCHS" / "Figures_v2"
FIG_OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------- Style ----------------
# Okabe-Ito colorblind-friendly subset
C_BLUE = "#0072B2"
C_ORANGE = "#D55E00"
C_AMBER = "#E69F00"
C_GREEN = "#009E73"
C_PINK = "#CC79A7"
C_SKY = "#56B4E9"
C_YELLOW = "#F0E442"
C_GREY = "#555555"

COHORT_COLORS = {
    "2001-established": C_BLUE,
    "2002-established": C_ORANGE,
}
COHORT_LABELS = {
    "2001-established": "2001 cohort",
    "2002-established": "2002 cohort",
}
TRAIT_LABELS = {
    "SummerPH": "Summer plant height (cm)",
    "FW": "Fresh weight (kg 5 m$^{-2}$)",
    "DW": "Dry weight (kg 5 m$^{-2}$)",
}
TRAIT_SHORT = {
    "SummerPH": "Summer height",
    "FW": "Fresh weight",
    "DW": "Dry weight",
}

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


def despine(ax, keep_left=True, keep_bottom=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if not keep_left:
        ax.spines["left"].set_visible(False)
    if not keep_bottom:
        ax.spines["bottom"].set_visible(False)


def save_fig(fig, name):
    for ext in ("pdf", "png"):
        path = FIG_OUT_DIR / f"{name}.{ext}"
        fig.savefig(path, format=ext, dpi=600)
    plt.close(fig)
    print(f"  saved: {name}.pdf + .png")


def load_csv(name):
    return pd.read_csv(RESULTS_DIR / name)


# ---------------- Fig 1: same-stand-age sensitivity ----------------
def fig1_same_stand_age():
    print("[Fig 1] same-stand-age sensitivity")
    df = load_csv("same_stand_age_sensitivity_final.csv")
    traits = ["FW", "DW", "SummerPH"]
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.2))
    stand_ages = sorted(df["stand_age"].unique())
    x = np.arange(len(stand_ages))
    width = 0.38
    for ax, trait in zip(axes, traits):
        for i, cohort in enumerate(["2001-established", "2002-established"]):
            sub = df[(df["cohort"] == cohort) & (df["stand_age"].isin(stand_ages))].sort_values("stand_age")
            vals = sub[trait].values
            ax.bar(x + (i - 0.5) * width, vals, width,
                   color=COHORT_COLORS[cohort], label=COHORT_LABELS[cohort],
                   edgecolor="black", linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f"Stand age {a}" for a in stand_ages])
        ax.set_ylabel(TRAIT_LABELS[trait])
        ax.set_title(f"({'abc'[list(traits).index(trait)]}) {TRAIT_SHORT[trait]}")
        despine(ax)
        ax.tick_params(direction="out", length=3)
    axes[0].legend(loc="upper left", handlelength=1.2, handletextpad=0.5)
    fig.tight_layout()
    save_fig(fig, "Fig1_same_stand_age_sensitivity")


# ---------------- Fig 2: variance + heritability ----------------
def fig2_variance_heritability():
    print("[Fig 2] variance + heritability")
    vc = load_csv("variance_components_final_by_cohort.csv")
    h2 = load_csv("heritability_final_by_cohort.csv")
    traits = ["SummerPH", "FW", "DW"]
    cohorts = ["2001-established", "2002-established"]
    components = [("V_G", "Genotype", C_BLUE),
                  ("V_GY", "G x Year", C_ORANGE),
                  ("V_plot", "Plot", C_AMBER),
                  ("V_e", "Residual", C_GREY)]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6),
                             gridspec_kw={"width_ratios": [1.5, 1]})
    # Panel A: stacked variance
    ax = axes[0]
    x = np.arange(len(traits))
    n_c = len(cohorts)
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
                    ax.bar(x[ti] + (ci - 0.5) * width, val, width,
                           bottom=bottom, color=color, edgecolor="black",
                           linewidth=0.4)
                    bottom += val
                else:
                    # boundary ~0 marker
                    pass
            # boundary marker for 2002 FW
            if cohort == "2002-established" and trait == "FW":
                ax.scatter(x[ti] + (ci - 0.5) * width, 0.05,
                           marker="*", s=60, color=C_PINK, zorder=5,
                           edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([TRAIT_SHORT[t] for t in traits])
    ax.set_ylabel("Variance component")
    ax.set_title("(A) Variance partitioning by cohort")
    # legend
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, ec="black", lw=0.5)
               for _, _, c in components]
    handles.append(plt.scatter([], [], marker="*", s=60, color=C_PINK,
                               edgecolor="black", linewidth=0.5))
    ax.legend(handles, [n for _, n, _ in components] + ["V_G at boundary"],
              loc="upper right", ncol=2, handlelength=1.2, handletextpad=0.4,
              columnspacing=1.0)
    despine(ax)
    ax.tick_params(direction="out", length=3)

    # Panel B: H2
    ax = axes[1]
    for ci, cohort in enumerate(cohorts):
        vals = []
        for trait in traits:
            r = h2[(h2["cohort"] == cohort) & (h2["trait"] == trait)]
            v = float(r["H2"].iloc[0]) if not r.empty else np.nan
            # boundary: 2002 FW ~ 2.3e-9 -> display as 0
            if v < 1e-6:
                v = 0.0
            vals.append(v)
        offset = (ci - 0.5) * 0.38
        bars = ax.bar(np.arange(len(traits)) + offset, vals, 0.38,
                      color=COHORT_COLORS[cohort], label=COHORT_LABELS[cohort],
                      edgecolor="black", linewidth=0.5)
        for bi, v in enumerate(vals):
            if cohort == "2002-established" and traits[bi] == "FW":
                ax.text(bi + offset, 0.02, "boundary", ha="center", va="bottom",
                        fontsize=6, color=C_PINK, rotation=90)
    ax.set_xticks(np.arange(len(traits)))
    ax.set_xticklabels([TRAIT_SHORT[t] for t in traits])
    ax.set_ylabel("Broad-sense H$^2$ (genotype-mean)")
    ax.set_ylim(0, 1.0)
    ax.set_title("(B) Heritability by cohort")
    ax.legend(loc="upper right", handlelength=1.2)
    despine(ax)
    ax.tick_params(direction="out", length=3)
    fig.tight_layout()
    save_fig(fig, "Fig2_variance_heritability")


# ---------------- Fig 3: BLUP ranking (cohort-specific) ----------------
def fig3_blup_ranking():
    print("[Fig 3] BLUP ranking")
    blup = load_csv("BLUP_final_by_cohort.csv")
    traits = ["FW", "DW", "SummerPH"]
    cohorts = ["2001-established", "2002-established"]
    fig, axes = plt.subplots(2, 3, figsize=(10, 6),
                             sharey="row")
    for ri, cohort in enumerate(cohorts):
        for ti, trait in enumerate(traits):
            ax = axes[ri, ti]
            sub = blup[(blup["cohort"] == cohort) & (blup["trait"] == trait)] \
                .sort_values("BLUP", ascending=False).reset_index(drop=True)
            n = len(sub)
            y = np.arange(n)
            vals = sub["BLUP"].values
            colors = [C_BLUE if i < 5 else C_GREY for i in range(n)]
            # 2002 FW: all shrunk -> use orange pale
            if cohort == "2002-established" and trait == "FW":
                colors = [C_ORANGE] * n
            ax.hbar = ax.barh(y, vals, color=colors, edgecolor="black",
                              linewidth=0.3, height=0.8)
            ax.invert_yaxis()
            # label top 5
            for i in range(min(5, n)):
                ax.text(vals[i] + (vals.max() - vals.min()) * 0.01, i,
                        sub["genotype"].iloc[i], va="center", fontsize=6.5,
                        fontweight="bold" if i < 5 else "normal")
            if cohort == "2002-established" and trait == "FW":
                ax.set_title(f"{TRAIT_SHORT[trait]}\n(BLUP shrunk; uninformative)",
                             fontsize=8)
            else:
                ax.set_title(TRAIT_SHORT[trait], fontsize=9)
            if ti == 0:
                ax.set_ylabel(COHORT_LABELS[cohort])
            despine(ax)
            ax.tick_params(direction="out", length=3)
            ax.set_yticks([])
    fig.tight_layout()
    save_fig(fig, "Fig3_blup_ranking_cohort")


# ---------------- AMMI / GGE score computation ----------------
def _genotype_year_matrix(trait, cohort="2001-established"):
    """Build genotype x year mean matrix for the 2001 cohort (complete 4-year)."""
    py = pd.read_csv(RESULTS_DIR / "plot_year_wide.csv")
    py = py[py["cohort"] == cohort]
    # aggregate to genotype x year mean
    mat = py.groupby(["genotype_ID", "calendar_year"])[trait].mean().unstack()
    # keep genotypes with complete 4 years
    mat = mat.dropna()
    return mat


def _ammi_scores(mat):
    """AMMI: additive model (row+col+grand) then SVD of residuals."""
    grand = mat.values.mean()
    row_eff = mat.values.mean(axis=1, keepdims=True) - grand
    col_eff = mat.values.mean(axis=0, keepdims=True) - grand
    resid = mat.values - grand - row_eff - col_eff
    U, s, Vt = np.linalg.svd(resid, full_matrices=False)
    # genotype scores: U * s ; year scores: Vt.T
    g_scores = U * s  # (n_g, k)
    y_scores = Vt.T  # (n_y, k)
    ss_total = np.sum(resid ** 2)
    pc_ss = s ** 2
    pc_pct = pc_ss / ss_total * 100
    return g_scores, y_scores, pc_pct, mat.index, mat.columns


def _gge_scores(mat):
    """GGE: center by environment (column) means only, then SVD."""
    col_mean = mat.values.mean(axis=0, keepdims=True)
    centered = mat.values - col_mean
    U, s, Vt = np.linalg.svd(centered, full_matrices=False)
    g_scores = U * s
    y_scores = Vt.T
    ss_total = np.sum(centered ** 2)
    pc_ss = s ** 2
    pc_pct = pc_ss / ss_total * 100
    # also keep genotype means for "which-won-where"
    g_means = mat.values.mean(axis=1)
    return g_scores, y_scores, pc_pct, g_means, mat.index, mat.columns


# ---------------- Fig 4: AMMI biplots ----------------
def fig4_ammi_biplot():
    print("[Fig 4] AMMI biplots")
    traits = ["SummerPH", "FW", "DW"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, trait in zip(axes, traits):
        try:
            mat = _genotype_year_matrix(trait)
            g, y, pct, gnames, ynames = _ammi_scores(mat)
        except Exception as e:
            ax.text(0.5, 0.5, f"error: {e}", transform=ax.transAxes)
            continue
        # scale for biplot (symmetric)
        scale = np.sqrt(np.sum(g[:, 0] ** 2) / np.sum(y[:, 0] ** 2)) if np.sum(y[:, 0] ** 2) > 0 else 1.0
        gx = g[:, 0]
        gy = g[:, 1] if g.shape[1] > 1 else np.zeros_like(gx)
        yx = y[:, 0] * scale
        yy = y[:, 1] * scale if y.shape[1] > 1 else np.zeros_like(yx)
        ax.scatter(gx, gy, color=C_BLUE, s=18, edgecolor="black",
                   linewidth=0.4, zorder=3)
        ax.scatter(yx, yy, color=C_ORANGE, marker="^", s=40,
                   edgecolor="black", linewidth=0.4, zorder=4)
        # label top 5 genotypes by mean
        top_idx = np.argsort(-mat.mean(axis=1).values)[:5]
        for i in top_idx:
            ax.annotate(gnames[i], (gx[i], gy[i]), fontsize=6,
                        color=C_BLUE, fontweight="bold",
                        xytext=(2, 2), textcoords="offset points")
        for j, yn in enumerate(ynames):
            ax.annotate(str(yn), (yx[j], yy[j]), fontsize=7,
                        color=C_ORANGE, fontweight="bold",
                        xytext=(3, 3), textcoords="offset points")
        ax.axhline(0, color="grey", lw=0.5, ls="--")
        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.set_xlabel(f"IPCA1 ({pct[0]:.1f}%)")
        if g.shape[1] > 1:
            ax.set_ylabel(f"IPCA2 ({pct[1]:.1f}%)")
        ax.set_title(TRAIT_SHORT[trait])
        despine(ax)
        ax.tick_params(direction="out", length=3)
    fig.suptitle("AMMI biplots (2001 cohort; temporal stability at a single site)",
                 y=1.02, fontsize=9, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "Fig4_ammi_biplot")


# ---------------- Fig 5: GGE biplots ----------------
def fig5_gge_biplot():
    print("[Fig 5] GGE biplots")
    traits = ["SummerPH", "FW", "DW"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    for ax, trait in zip(axes, traits):
        try:
            mat = _genotype_year_matrix(trait)
            g, y, pct, g_means, gnames, ynames = _gge_scores(mat)
        except Exception as e:
            ax.text(0.5, 0.5, f"error: {e}", transform=ax.transAxes)
            continue
        scale = np.sqrt(np.sum(g[:, 0] ** 2) / np.sum(y[:, 0] ** 2)) if np.sum(y[:, 0] ** 2) > 0 else 1.0
        gx = g[:, 0]
        gy = g[:, 1] if g.shape[1] > 1 else np.zeros_like(gx)
        yx = y[:, 0] * scale
        yy = y[:, 1] * scale if y.shape[1] > 1 else np.zeros_like(yx)
        # marker size proportional to genotype mean
        g_norm = (g_means - g_means.min()) / (g_means.max() - g_means.min() + 1e-9)
        sizes = 15 + g_norm * 35
        ax.scatter(gx, gy, color=C_BLUE, s=sizes, edgecolor="black",
                   linewidth=0.4, zorder=3, alpha=0.85)
        ax.scatter(yx, yy, color=C_ORANGE, marker="^", s=50,
                   edgecolor="black", linewidth=0.4, zorder=4)
        top_idx = np.argsort(-g_means)[:5]
        for i in top_idx:
            ax.annotate(gnames[i], (gx[i], gy[i]), fontsize=6,
                        color=C_BLUE, fontweight="bold",
                        xytext=(2, 2), textcoords="offset points")
        for j, yn in enumerate(ynames):
            ax.annotate(str(yn), (yx[j], yy[j]), fontsize=7,
                        color=C_ORANGE, fontweight="bold",
                        xytext=(3, 3), textcoords="offset points")
        ax.axhline(0, color="grey", lw=0.5, ls="--")
        ax.axvline(0, color="grey", lw=0.5, ls="--")
        ax.set_xlabel(f"PC1 ({pct[0]:.1f}%)")
        if g.shape[1] > 1:
            ax.set_ylabel(f"PC2 ({pct[1]:.1f}%)")
        ax.set_title(TRAIT_SHORT[trait])
        despine(ax)
        ax.tick_params(direction="out", length=3)
    fig.suptitle("GGE biplots (2001 cohort; temporal stability at a single site)",
                 y=1.02, fontsize=9, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "Fig5_gge_biplot")


# ---------------- Fig 6: mean vs stability ----------------
def fig6_mean_stability():
    print("[Fig 6] mean vs stability")
    traits = ["SummerPH", "FW", "DW"]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.4))
    for ax, trait in zip(axes, traits):
        try:
            mat = _genotype_year_matrix(trait)
            g, y, pct, gnames, ynames = _ammi_scores(mat)
            # ASV (Purchase) using first 2 IPCA
            ss1 = np.sum(g[:, 0] ** 2)
            ss2 = np.sum(g[:, 1] ** 2) if g.shape[1] > 1 else 0
            if ss2 > 0:
                asv = np.sqrt(ss1 / ss2) * np.abs(g[:, 0]) + np.abs(g[:, 1])
            else:
                asv = np.abs(g[:, 0])
            g_means = mat.mean(axis=1).values
            ax.scatter(asv, g_means, color=C_BLUE, s=22, edgecolor="black",
                       linewidth=0.4, zorder=3)
            top_idx = np.argsort(-g_means)[:5]
            for i in top_idx:
                ax.annotate(gnames[i], (asv[i], g_means[i]), fontsize=6,
                            color=C_BLUE, fontweight="bold",
                            xytext=(3, 3), textcoords="offset points")
            ax.set_xlabel("AMMI Stability Value (lower = more stable)")
            ax.set_ylabel(TRAIT_LABELS[trait])
            ax.set_title(TRAIT_SHORT[trait])
            despine(ax)
            ax.tick_params(direction="out", length=3)
        except Exception as e:
            ax.text(0.5, 0.5, f"error: {e}", transform=ax.transAxes)
    fig.suptitle("Mean performance vs temporal stability (2001 cohort)",
                 y=1.03, fontsize=9, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "Fig6_mean_stability")


# ---------------- Fig 7: BLUP correlation (cohort-specific) ----------------
def fig7_blup_correlation():
    print("[Fig 7] BLUP correlation")
    corr = load_csv("BLUP_correlation_by_cohort.csv")
    cohorts = ["2001-established", "2002-established"]
    trait_pairs = [("SummerPH", "FW"), ("SummerPH", "DW"), ("FW", "DW")]
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6))
    for ax, cohort in zip(axes, cohorts):
        mat = np.full((3, 3), np.nan)
        labels = ["SummerPH", "FW", "DW"]
        for i, (t1, t2) in enumerate(trait_pairs):
            r = corr[(corr["cohort"] == cohort) & (corr["trait1"] == t1) & (corr["trait2"] == t2)]
            if not r.empty:
                mat[i, 2 - labels.index(t2)] = float(r["r"].iloc[0])
                mat[2 - labels.index(t2), i] = float(r["r"].iloc[0])
        # diagonal = 1
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
                    ax.text(jj, ii, f"{v:.2f}", ha="center", va="center",
                            color=color, fontsize=8, fontweight="bold")
        ax.tick_params(direction="out", length=3)
    cbar_ax = fig.add_axes([0.93, 0.18, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax, label="Pearson r")
    cbar.ax.tick_params(direction="out", length=3)
    fig.subplots_adjust(left=0.07, right=0.9, wspace=0.35)
    save_fig(fig, "Fig7_blup_correlation_cohort")


# ---------------- Main ----------------
def main():
    print(f"Output dir: {FIG_OUT_DIR}")
    try:
        fig1_same_stand_age()
    except Exception as e:
        print(f"  Fig1 error: {e}")
    try:
        fig2_variance_heritability()
    except Exception as e:
        print(f"  Fig2 error: {e}")
    try:
        fig3_blup_ranking()
    except Exception as e:
        print(f"  Fig3 error: {e}")
    try:
        fig4_ammi_biplot()
    except Exception as e:
        print(f"  Fig4 error: {e}")
    try:
        fig5_gge_biplot()
    except Exception as e:
        print(f"  Fig5 error: {e}")
    try:
        fig6_mean_stability()
    except Exception as e:
        print(f"  Fig6 error: {e}")
    try:
        fig7_blup_correlation()
    except Exception as e:
        print(f"  Fig7 error: {e}")
    print("Done.")


if __name__ == "__main__":
    main()
