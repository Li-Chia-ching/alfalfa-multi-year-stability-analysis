"""
合并后的可视化脚本：整合 make_fig3_caterpillar.py 和 make_figures_v2.py
包含7个Figure的生成功能，并优化Figure 4-6的标签防重叠算法。

输入数据来源：04_Results/Revision_Final/Final_Statistical_Reanalysis/
输出位置：06_Manuscript/Formal_To_NZCHS/Figures_v2/

Style: Nature/Springer (no gridlines, bold axes, minimal titles),
Okabe-Ito colorblind-friendly palette, Arial font, despine.
"""
import os
import sys
import warnings
import json
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

# 尝试导入adjustText用于标签防重叠
try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
    print("✓ adjustText库已安装，将使用标签防重叠算法")
except ImportError:
    HAS_ADJUST_TEXT = False
    print("⚠ adjustText库未安装，将使用简单标签偏移")
    print("  安装命令: pip install adjustText")

warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------- Paths ----------------
PROJECT_ROOT = Path(r"C:\Users\lijia\Documents\R Workplace\Sino-Australian_Alfalfa_Project")
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
C_GRAY = "#9AA0A6"  # 普通基因型颜色

# Figure 3专用颜色映射
HIGHLIGHT_COLORS = {
    "Sitel": C_BLUE,
    "Xinjiang_Daye": C_ORANGE,
    "ACA542": C_AMBER,
    "L33": C_GREEN,
    "Algonquin": C_PINK,
}

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
TRAIT_UNIT = {
    "SummerPH": "cm",
    "FW": "kg 5 m$^{-2}$",
    "DW": "kg 5 m$^{-2}$",
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
    """移除顶部和右侧边框，保持底部和左侧边框"""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if not keep_left:
        ax.spines["left"].set_visible(False)
    if not keep_bottom:
        ax.spines["bottom"].set_visible(False)


def save_fig(fig, name):
    """保存图片为PDF和PNG格式"""
    for ext in ("pdf", "png"):
        path = FIG_OUT_DIR / f"{name}.{ext}"
        fig.savefig(path, format=ext, dpi=600)
    plt.close(fig)
    print(f"  saved: {name}.pdf + .png")


def load_csv(name):
    """加载CSV数据文件"""
    return pd.read_csv(RESULTS_DIR / name)


def label_for(genotype):
    """获取基因型标签（用于Figure 3）"""
    if genotype in HIGHLIGHT_COLORS:
        return genotype
    return genotype


def color_for(genotype):
    """获取基因型颜色（用于Figure 3）"""
    if genotype in HIGHLIGHT_COLORS:
        return HIGHLIGHT_COLORS[genotype]
    return C_GRAY


def panel_letter(i):
    """生成面板字母标签"""
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i]


def adjust_labels_with_leader_lines(ax, texts, points, 
                                   force_text: Tuple[float, float] = (0.2, 0.25),
                                   force_points: Tuple[float, float] = (0.2, 0.25),
                                   expand: Tuple[float, float] = (1.05, 1.1),
                                   arrowprops: Optional[Dict[str, Any]] = None):
    """
    调整标签位置以避免重叠，并添加引线。
    包含终极修复：动态扩展坐标轴，确保任何标签都不会飞出图外。
    """
    if not HAS_ADJUST_TEXT:
        print("  使用简单标签偏移算法")
        for i, (text, (x, y)) in enumerate(zip(texts, points)):
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            x_range = xlim[1] - xlim[0]
            y_range = ylim[1] - ylim[0]
            
            offset_x = 0.02 * x_range
            offset_y = 0.02 * y_range
            
            if x > xlim[1] - 0.1 * x_range:
                offset_x = -0.05 * x_range
            elif x < xlim[0] + 0.1 * x_range:
                offset_x = 0.05 * x_range
                
            if y > ylim[1] - 0.1 * y_range:
                offset_y = -0.05 * y_range
            elif y < ylim[0] + 0.1 * y_range:
                offset_y = 0.05 * y_range
            
            text.set_position((x + offset_x, y + offset_y))
    else:
        print("  使用adjustText标签防重叠算法")
        
        if arrowprops is None:
            arrowprops = dict(
                arrowstyle='-', 
                color='black', 
                lw=0.5,
                alpha=0.6
            )
        
        # 修复1：将 avoid_axes 恢复为 True，让算法尽量在内部消化
        adjust_text(
            texts,
            x=[p[0] for p in points],
            y=[p[1] for p in points],
            ax=ax,
            force_text=force_text,
            force_points=force_points,
            expand=expand,
            arrowprops=arrowprops,
            only_move={'text': 'xy'},
            avoid_self=True,
            avoid_points=True,
            avoid_axes=True, # 已恢复
            max_iter=800,
            tol=0.005
        )

    # 修复2：终极坐标轴自适应（适用于所有图表和算法分支）
    # 在所有标签排布完成后，抓取文本真实的渲染边界，反推需要的坐标轴范围
    try:
        fig = ax.get_figure()
        fig.canvas.draw() # 强制执行一次渲染以获取文本的真实像素尺寸
        renderer = fig.canvas.get_renderer()
        
        # 获取当前的坐标轴限制
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        
        # 遍历所有文本，寻找极值
        for text in texts:
            bbox = text.get_window_extent(renderer=renderer)
            # 将屏幕像素边界转换回数据的XY坐标
            bbox_data = bbox.transformed(ax.transData.inverted())
            xmin = min(xmin, bbox_data.xmin)
            xmax = max(xmax, bbox_data.xmax)
            ymin = min(ymin, bbox_data.ymin)
            ymax = max(ymax, bbox_data.ymax)
            
        # 根据计算出的极值，再外扩 8% 作为安全留白
        x_margin = (xmax - xmin) * 0.08
        y_margin = (ymax - ymin) * 0.08
        
        # 重新设定坐标轴
        ax.set_xlim(xmin - x_margin, xmax + x_margin)
        ax.set_ylim(ymin - y_margin, ymax + y_margin)
        
    except Exception as e:
        print(f"  坐标轴自适应拓展失败: {e}")


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
            if cohort == "2002-established" and trait == "FW":
                ax.scatter(x[ti] + (ci - 0.5) * width, 0.05,
                           marker="*", s=60, color=C_PINK, zorder=5,
                           edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([TRAIT_SHORT[t] for t in traits])
    ax.set_ylabel("Variance component")
    ax.set_title("(A) Variance partitioning by cohort")
    handles = [plt.Rectangle((0, 0), 1, 1, color=c, ec="black", lw=0.5)
               for _, _, c in components]
    handles.append(plt.scatter([], [], marker="*", s=60, color=C_PINK,
                               edgecolor="black", linewidth=0.5))
    ax.legend(handles, [n for _, n, _ in components] + ["V_G at boundary"],
              loc="upper right", ncol=2, handlelength=1.2, handletextpad=0.4,
              columnspacing=1.0)
    despine(ax)
    ax.tick_params(direction="out", length=3)

    ax = axes[1]
    for ci, cohort in enumerate(cohorts):
        vals = []
        for trait in traits:
            r = h2[(h2["cohort"] == cohort) & (h2["trait"] == trait)]
            v = float(r["H2"].iloc[0]) if not r.empty else np.nan
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
    print("[Fig 3] BLUP ranking (caterpillar)")
    blup = load_csv("BLUP_final_by_cohort.csv")
    
    meta = {
        "input_file": "BLUP_final_by_cohort.csv",
        "input_columns": list(blup.columns),
        "n_rows": int(len(blup)),
        "counts": {},
        "identities": {},
        "uncertainty_intervals": "NONE"
    }
    
    for cohort in ["2001-established", "2002-established"]:
        for trait in ["SummerPH", "FW", "DW"]:
            count = len(blup[(blup["cohort"] == cohort) & (blup["trait"] == trait)])
            meta["counts"][f"{cohort[:4]}_{trait}"] = count
    
    meta["identities"]["Xinjiang_Daye_exists"] = "Xinjiang_Daye" in blup["genotype"].unique()
    meta["identities"]["ACA542_exists"] = "ACA542" in blup["genotype"].unique()
    
    meta_log = FIG_OUT_DIR / "Fig3_input_metadata.json"
    meta_log.write_text(json.dumps(meta, indent=2, ensure_ascii=False))
    
    fig, axes = plt.subplots(2, 3, figsize=(12, 13.2))
    panel_idx = 0
    
    for ri, cohort in enumerate(["2001-established", "2002-established"]):
        for ci, trait in enumerate(["SummerPH", "FW", "DW"]):
            ax = axes[ri, ci]
            sub = blup[(blup["cohort"] == cohort) & (blup["trait"] == trait)].copy()
            sub = sub.sort_values("BLUP", ascending=False).reset_index(drop=True)
            n = len(sub)
            y = np.arange(n)
            blups = sub["BLUP"].values
            labels = sub["genotype"].values
            colors = [color_for(g) for g in labels]
            
            ax.scatter(blups, y, s=26, c=colors, edgecolors="black",
                       linewidths=0.35, zorder=3)
            
            bmin, bmax = blups.min(), blups.max()
            span = bmax - bmin
            
            if trait == "FW" and cohort == "2002-established":
                m = float(np.mean(blups))
                ax.set_xlim(m - 0.35, m + 0.35)
                ax.axvline(m, color=C_GRAY, lw=0.8, ls="--", zorder=1, alpha=0.7)
                ax.text(m + 0.36 - (m + 0.35 - (m - 0.35)) * 0.02, n - 0.5,
                        "V$_G$ $\\approx$ 0;\ngenotype\ndiscrimination\nlimited",
                        fontsize=7, color=C_PINK, va="top", ha="right",
                        fontweight="bold")
            else:
                margin = max(span * 0.06, 0.001)
                ax.set_xlim(bmin - margin, bmax + margin)
            
            ax.set_yticks(y)
            ax.set_yticklabels(labels, fontsize=7.5)
            ax.set_ylim(n - 0.5, -0.5)
            ax.set_title(f"({panel_letter(panel_idx)})  {COHORT_LABELS[cohort]} — {TRAIT_SHORT[trait]}",
                         loc="left")
            ax.set_xlabel(f"BLUP ({TRAIT_UNIT[trait]})")
            ax.xaxis.grid(True, which="major", color="#E0E0E0", lw=0.6, zorder=0)
            ax.yaxis.grid(False)
            ax.set_axisbelow(True)
            despine(ax)
            ax.tick_params(axis="y", length=0)
            ax.tick_params(axis="x", direction="out", length=3)
            panel_idx += 1
    
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=C_GRAY,
                          markeredgecolor="black", markersize=7, label="Other genotypes")]
    for genotype, color in HIGHLIGHT_COLORS.items():
        handles.append(plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color,
                                  markeredgecolor="black", markersize=7, label=genotype))
    
    leg = fig.legend(handles=handles, loc="upper center", ncol=6,
                     bbox_to_anchor=(0.5, 1.005), frameon=False,
                     handletextpad=0.4, columnspacing=1.3, fontsize=8.5)
    
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    fig.subplots_adjust(hspace=0.32, wspace=0.42)
    
    for ext in ("png", "pdf"):
        p = FIG_OUT_DIR / f"Fig3_BLUP_Caterpillar_v3.{ext}"
        fig.savefig(p, format=ext, dpi=600)
    
    try:
        tif_path = FIG_OUT_DIR / "Fig3_BLUP_Caterpillar_v3.tif"
        fig.savefig(tif_path, format="tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    except Exception as e:
        print(f"  tiff skipped: {e}")
    plt.close(fig)


# ---------------- AMMI / GGE score computation ----------------
def _genotype_year_matrix(trait, cohort="2001-established"):
    py = pd.read_csv(RESULTS_DIR / "plot_year_wide.csv")
    py = py[py["cohort"] == cohort]
    mat = py.groupby(["genotype_ID", "calendar_year"])[trait].mean().unstack()
    mat = mat.dropna()
    return mat

def _ammi_scores(mat):
    grand = mat.values.mean()
    row_eff = mat.values.mean(axis=1, keepdims=True) - grand
    col_eff = mat.values.mean(axis=0, keepdims=True) - grand
    resid = mat.values - grand - row_eff - col_eff
    U, s, Vt = np.linalg.svd(resid, full_matrices=False)
    g_scores = U * s 
    y_scores = Vt.T 
    ss_total = np.sum(resid ** 2)
    pc_ss = s ** 2
    pc_pct = pc_ss / ss_total * 100
    return g_scores, y_scores, pc_pct, mat.index, mat.columns

def _gge_scores(mat):
    col_mean = mat.values.mean(axis=0, keepdims=True)
    centered = mat.values - col_mean
    U, s, Vt = np.linalg.svd(centered, full_matrices=False)
    g_scores = U * s
    y_scores = Vt.T
    ss_total = np.sum(centered ** 2)
    pc_ss = s ** 2
    pc_pct = pc_ss / ss_total * 100
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
        
        scale = np.sqrt(np.sum(g[:, 0] ** 2) / np.sum(y[:, 0] ** 2)) if np.sum(y[:, 0] ** 2) > 0 else 1.0
        gx = g[:, 0]
        gy = g[:, 1] if g.shape[1] > 1 else np.zeros_like(gx)
        yx = y[:, 0] * scale
        yy = y[:, 1] * scale if y.shape[1] > 1 else np.zeros_like(yx)
        
        ax.scatter(gx, gy, color=C_BLUE, s=18, edgecolor="black", linewidth=0.4, zorder=3)
        ax.scatter(yx, yy, color=C_ORANGE, marker="^", s=40, edgecolor="black", linewidth=0.4, zorder=4)
        
        texts = []
        points = []
        
        top_idx = np.argsort(-mat.mean(axis=1).values)[:5]
        for i in top_idx:
            text = ax.annotate(gnames[i], (gx[i], gy[i]), fontsize=6,
                               color=C_BLUE, fontweight="bold",
                               xytext=(2, 2), textcoords="offset points")
            texts.append(text)
            points.append((gx[i], gy[i]))
        
        for j, yn in enumerate(ynames):
            text = ax.annotate(str(yn), (yx[j], yy[j]), fontsize=7,
                               color=C_ORANGE, fontweight="bold",
                               xytext=(3, 3), textcoords="offset points")
            texts.append(text)
            points.append((yx[j], yy[j]))
        
        # 关键修复点：增加绘图区边界余量，防止文本在调整时被挤出视口外
        ax.margins(0.2)
        
        adjust_labels_with_leader_lines(
            ax, texts, points,
            force_text=(0.15, 0.2),
            force_points=(0.15, 0.2),
            expand=(1.05, 1.1)
        )
        
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
        
        g_norm = (g_means - g_means.min()) / (g_means.max() - g_means.min() + 1e-9)
        sizes = 15 + g_norm * 35
        ax.scatter(gx, gy, color=C_BLUE, s=sizes, edgecolor="black", linewidth=0.4, zorder=3, alpha=0.85)
        ax.scatter(yx, yy, color=C_ORANGE, marker="^", s=50, edgecolor="black", linewidth=0.4, zorder=4)
        
        texts = []
        points = []
        
        top_idx = np.argsort(-g_means)[:5]
        for i in top_idx:
            text = ax.annotate(gnames[i], (gx[i], gy[i]), fontsize=6,
                               color=C_BLUE, fontweight="bold",
                               xytext=(2, 2), textcoords="offset points")
            texts.append(text)
            points.append((gx[i], gy[i]))
        
        for j, yn in enumerate(ynames):
            text = ax.annotate(str(yn), (yx[j], yy[j]), fontsize=7,
                               color=C_ORANGE, fontweight="bold",
                               xytext=(3, 3), textcoords="offset points")
            texts.append(text)
            points.append((yx[j], yy[j]))
        
        # 关键修复点：增加绘图区边界余量
        ax.margins(0.2)
        
        adjust_labels_with_leader_lines(
            ax, texts, points,
            force_text=(0.15, 0.2),
            force_points=(0.15, 0.2),
            expand=(1.05, 1.1)
        )
        
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
            ss1 = np.sum(g[:, 0] ** 2)
            ss2 = np.sum(g[:, 1] ** 2) if g.shape[1] > 1 else 0
            if ss2 > 0:
                asv = np.sqrt(ss1 / ss2) * np.abs(g[:, 0]) + np.abs(g[:, 1])
            else:
                asv = np.abs(g[:, 0])
            g_means = mat.mean(axis=1).values
            
            ax.scatter(asv, g_means, color=C_BLUE, s=22, edgecolor="black", linewidth=0.4, zorder=3)
            
            texts = []
            points = []
            
            top_idx = np.argsort(-g_means)[:5]
            for i in top_idx:
                text = ax.annotate(gnames[i], (asv[i], g_means[i]), fontsize=6,
                                   color=C_BLUE, fontweight="bold",
                                   xytext=(3, 3), textcoords="offset points")
                texts.append(text)
                points.append((asv[i], g_means[i]))
            
            # 关键修复点：增加绘图区边界余量
            ax.margins(0.2)
            
            adjust_labels_with_leader_lines(
                ax, texts, points,
                force_text=(0.15, 0.2),
                force_points=(0.15, 0.2),
                expand=(1.05, 1.1)
            )
            
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


# ---------------- Fig 7: BLUP correlation ----------------
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
    
    if HAS_ADJUST_TEXT:
        print("✓ adjustText库可用，Figure 4-6将使用标签防重叠算法")
    else:
        print("⚠ adjustText库不可用，Figure 4-6将使用简单标签偏移")
        print("  安装命令: pip install adjustText")
    
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
            fig_func()
        except Exception as e:
            print(f"  {fig_name} error: {e}")
    
    print("Done.")

if __name__ == "__main__":
    main()
