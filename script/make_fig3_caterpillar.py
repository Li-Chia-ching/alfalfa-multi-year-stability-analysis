"""
Fig. 3 — Cohort-specific BLUP caterpillar plot (Cleveland dot plot).

STATUS: Figure-only task. Statistical analysis is FROZEN.
  - Reads frozen BLUP point estimates from BLUP_final_by_cohort.csv.
  - Does NOT recompute genotype identities, BLUPs, model, H2, or variance components.
  - Does NOT refit any model.

UNCERTAINTY INTERVALS:
  The frozen pipeline (01_final_reanalysis.R) extracted ranef(m)$genotype WITHOUT
  condVar=TRUE and did not save model .rds objects. Therefore no model-derived
  SE or 95% CI is available in the frozen output. Per spec ("Never invent
  uncertainty values"), this figure shows BLUP POINT ESTIMATES ONLY, with no
  horizontal uncertainty intervals. The caption states this explicitly.
  Adding 95% CI would require re-running the frozen pipeline with
  condVar=TRUE (estimates unchanged) — pending user approval.

OUTPUT:
  Fig3_BLUP_Caterpillar_v3.png  (600 dpi)
  Fig3_BLUP_Caterpillar_v3.pdf  (vector)
  Fig3_BLUP_Caterpillar_v3.tif  (600 dpi, from same source)
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---------------- Paths ----------------
PROJECT = Path(r"C:\Users\lijia\Documents\R Workplace\Sino-Australian_Alfalfa_Project")
SRC = PROJECT / "04_Results" / "Revision_Final" / "Final_Statistical_Reanalysis" / "BLUP_final_by_cohort.csv"
OUT = PROJECT / "06_Manuscript" / "Formal_To_NZCHS" / "Figures_v2"
OUT.mkdir(parents=True, exist_ok=True)
META_LOG = OUT / "Fig3_input_metadata.json"

# ---------------- Style ----------------
C_GRAY = "#9AA0A6"      # ordinary genotypes
C_BLUE = "#0072B2"      # Sitel
C_ORANGE = "#D55E00"    # Xinjiang_Daye (ACA541)
C_AMBER = "#E69F00"     # ACA542
C_GREEN = "#009E73"     # L33
C_PINK = "#CC79A7"      # Algonquin

HIGHLIGHT = {
    "Sitel":          (C_BLUE, "Sitel"),
    "Xinjiang_Daye":  (C_ORANGE, "Xinjiang_Daye (ACA541)"),
    "ACA542":         (C_AMBER, "ACA542"),
    "L33":            (C_GREEN, "L33"),
    "Algonquin":      (C_PINK, "Algonquin"),
}

TRAITS = ["SummerPH", "FW", "DW"]
TRAIT_TITLE = {
    "SummerPH": "Summer height",
    "FW": "Fresh weight",
    "DW": "Dry weight",
}
TRAIT_UNIT = {
    "SummerPH": "cm",
    "FW": "kg 5 m$^{-2}$",
    "DW": "kg 5 m$^{-2}$",
}
COHORTS = ["2001-established", "2002-established"]
COHORT_LABEL = {"2001-established": "2001 cohort", "2002-established": "2002 cohort"}

rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 11.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.labelweight": "bold",
    "axes.linewidth": 0.9,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 7.5,
    "figure.dpi": 100,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def despine(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)


def label_for(g):
    if g in HIGHLIGHT:
        return HIGHLIGHT[g][1]
    return g


def color_for(g):
    if g in HIGHLIGHT:
        return HIGHLIGHT[g][0]
    return C_GRAY


def panel_letter(i):
    return "ABCDEF"[i]


def make_fig():
    df = pd.read_csv(SRC)
    # metadata log
    meta = {
        "input_file": str(SRC.name),
        "input_columns": list(df.columns),
        "n_rows": int(len(df)),
        "counts": {
            f"{c[:4]}_{t}": int(len(df[(df.cohort == c) & (df.trait == t)]))
            for c in COHORTS for t in TRAITS
        },
        "identities": {
            "Xinjiang_Daye_is_ACA541": "Xinjiang_Daye" in df.genotype.unique(),
            "ACA542_separate": "ACA542" in df.genotype.unique(),
            "L33_in_both": all(
                len(df[(df.genotype == "L33") & (df.cohort == c)]) > 0 for c in COHORTS
            ),
            "Algonquin_in_both": all(
                len(df[(df.genotype == "Algonquin") & (df.cohort == c)]) > 0 for c in COHORTS
            ),
        },
        "uncertainty_intervals": "NONE — frozen pipeline saved BLUP point estimates only; no SE/CI available without refit. Point-only per spec §11/§18.",
    }
    META_LOG.write_text(json.dumps(meta, indent=2, ensure_ascii=False))
    print("metadata log:", META_LOG.name)

    fig, axes = plt.subplots(2, 3, figsize=(12, 13.2))
    panel_idx = 0
    for ri, coh in enumerate(COHORTS):
        for ci, trait in enumerate(TRAITS):
            ax = axes[ri, ci]
            sub = df[(df.cohort == coh) & (df.trait == trait)].copy()
            # sort by BLUP descending -> highest at top
            sub = sub.sort_values("BLUP", ascending=False).reset_index(drop=True)
            n = len(sub)
            y = np.arange(n)  # 0 at top after invert
            blups = sub["BLUP"].values
            labels = [label_for(g) for g in sub["genotype"].values]
            colors = [color_for(g) for g in sub["genotype"].values]

            # plot points (Cleveland dot plot — no stems, no rectangles)
            ax.scatter(blups, y, s=26, c=colors, edgecolors="black",
                       linewidths=0.35, zorder=3)

            # x-axis limits (independent per panel; special handling 2002 FW)
            bmin, bmax = blups.min(), blups.max()
            span = bmax - bmin
            if trait == "FW" and coh == "2002-established":
                # boundary: all BLUPs ~ identical; use a sensible non-degenerate
                # range centered on cohort mean to avoid misleading zoom AND
                # avoid degenerate axis. Differences (10^-9) are NOT magnified.
                m = float(np.mean(blups))
                ax.set_xlim(m - 0.35, m + 0.35)
                # subtle dashed line at cohort mean = shrinkage target
                ax.axvline(m, color=C_GRAY, lw=0.8, ls="--", zorder=1, alpha=0.7)
                ax.text(m + 0.36 - (m + 0.35 - (m - 0.35)) * 0.02, n - 0.5,
                        "V$_G$ $\\approx$ 0;\ngenotype\ndiscrimination\nlimited",
                        fontsize=7, color=C_PINK, va="top", ha="right",
                        fontweight="bold")
            else:
                margin = max(span * 0.06, 0.001)
                ax.set_xlim(bmin - margin, bmax + margin)

            # y-axis: genotype labels
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7.5)
            ax.set_ylim(n - 0.5, -0.5)  # invert: rank 1 at top

            # optional rank annotation for 2001 biomass (FW, DW) only
            if coh == "2001-established" and trait in ("FW", "DW"):
                # rank 1 = Sitel, rank 2 = Xinjiang_Daye (verified top-2)
                for r, g in [(0, "Sitel"), (1, "Xinjiang_Daye")]:
                    if sub["genotype"].iloc[r] == g:
                        ax.text(blups[r], y[r] - 0.28, f"Rank {r+1}",
                                fontsize=6.5, color=HIGHLIGHT[g][0],
                                ha="center", va="top", fontweight="bold")

            # panel letter + title
            ax.set_title(f"({panel_letter(panel_idx)})  {COHORT_LABEL[coh]} — {TRAIT_TITLE[trait]}",
                         loc="left")
            ax.set_xlabel(f"BLUP ({TRAIT_UNIT[trait]})")
            # subtle vertical gridlines only, no horizontal
            ax.xaxis.grid(True, which="major", color="#E0E0E0", lw=0.6, zorder=0)
            ax.yaxis.grid(False)
            ax.set_axisbelow(True)
            despine(ax)
            ax.tick_params(axis="y", length=0)
            ax.tick_params(axis="x", direction="out", length=3)
            panel_idx += 1

    # shared legend for highlights (top, outside panels)
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=c,
                          markeredgecolor="black", markersize=7, label=lbl)
               for _, (c, lbl) in [(0, (C_GRAY, "Other genotypes"))] +
               [(k, v) for k, v in HIGHLIGHT.items()]]
    leg = fig.legend(handles=handles, loc="upper center", ncol=6,
                     bbox_to_anchor=(0.5, 1.005), frameon=False,
                     handletextpad=0.4, columnspacing=1.3, fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.subplots_adjust(hspace=0.32, wspace=0.42)

    # save png, pdf, tif
    for ext in ("png", "pdf"):
        p = OUT / f"Fig3_BLUP_Caterpillar_v3.{ext}"
        fig.savefig(p, format=ext, dpi=600)
        print("  saved:", p.name)
    # TIFF from same source (600 dpi)
    try:
        fig.savefig(OUT / "Fig3_BLUP_Caterpillar_v3.tif", format="tiff",
                    dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        print("  saved: Fig3_BLUP_Caterpillar_v3.tif")
    except Exception as e:
        print("  tiff skipped:", e)
    plt.close(fig)


if __name__ == "__main__":
    make_fig()
    print("Done.")
