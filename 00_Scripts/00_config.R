# 00_config.R
# ---------------------------------------------------------------------------
# Scientific Manuscript Reconstruction -- Sino-Australian Alfalfa Project
# SINGLE SOURCE OF TRUTH: paths, packages, graphical theme.
# All numbered scripts (01_..10_) source this file first.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
  library(showtext)  # 新增：用于接管全局字体渲染
})

# --- Project paths (forward slashes; R on Windows accepts them) -------------
ROOT <- "C:/Users/lijia/Documents/R Workplace/Sino-Australian_Alfalfa_Project"
RAW  <- file.path(ROOT, "01_Raw_Phenotype_Data")
RES  <- file.path(ROOT, "04_Results")
FIG  <- file.path(RES, "Figures")
TAB  <- file.path(RES, "Tables")
STAT <- file.path(RES, "Statistics")
SUPP <- file.path(RES, "Supplementary")
INT  <- file.path(ROOT, "03_Analysis", "interim")

for (d in c(FIG, TAB, STAT, SUPP, INT)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# --- Traits analysed (core trial = 2002-2005) -------------------------------
TRAITS      <- c("SummerPH", "FW", "DW")
TRAIT_UNIT  <- c(SummerPH = "cm", FW = "g plot-1", DW = "g plot-1")
TRAIT_LABEL <- c(SummerPH = "Summer plant height", FW = "Fresh weight", DW = "Dry weight")

# --- Publication-ready figure style (unified across all figures) -----------
# 开启 showtext 自动接管图形设备
showtext_auto()

FIG_W <- 22   # cm (从17扩大至22，给文字和图形留出充足空间)
FIG_H <- 15   # cm (从12扩大至15)
DPI   <- 600

# 显著提升基础字号 (从11提升至15)，确保肉眼可读
FIG_THEME <- theme_classic(base_size = 15) +
  theme(
    text            = element_text(family = "sans", size = 15, colour = "black"),
    axis.title      = element_text(size = 16, face = "bold"), # 坐标轴标题加粗加大
    axis.text       = element_text(size = 14, colour = "black"), # 坐标轴刻度文字
    axis.line       = element_line(linewidth = 0.8, colour = "black"),
    panel.border    = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text      = element_text(size = 14, face = "bold"),
    legend.key      = element_blank(),
    legend.position = "right",
    legend.text     = element_text(size = 13),
    plot.tag        = element_text(face = "bold")
  )

# --- Colour-blind friendly palette (Okabe-Ito) with alpha support ----------
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Apply via PAL list; users can additionally set alpha in geoms (e.g., alpha = 0.7)
PAL <- list(
  scale_colour_manual(values = cb_palette),
  scale_fill_manual(values = cb_palette)
)

# --- Export helper: PDF + SVG + PNG(600 dpi) for every figure --------------
export_fig <- function(g, name, w = FIG_W, h = FIG_H) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), g, width = w, height = h,
         units = "cm", device = "pdf")
  ggsave(file.path(FIG, paste0(name, ".svg")), g, width = w, height = h,
         units = "cm", device = "svg")
  ggsave(file.path(FIG, paste0(name, ".png")), g, width = w, height = h,
         units = "cm", dpi = DPI, device = "png")
  invisible(NULL)
}

# --- Figure + companion CSV export (ensures full reproducibility) ----------
export_fig_with_data <- function(g, plot_data, name, w = FIG_W, h = FIG_H) {
  # 1. Export the figure in three formats
  export_fig(g, name, w, h)
  # 2. Export the underlying data used for plotting
  data_path <- file.path(FIG, paste0(name, "_plot_data.csv"))
  write.csv(plot_data, data_path, row.names = FALSE)
  message(sprintf("Saved: %s.* and %s", name, basename(data_path)))
  invisible(NULL)
}

# Rounding helper for tables
rr <- function(x, n = 3) formatC(round(x, n), format = "f", digits = n)