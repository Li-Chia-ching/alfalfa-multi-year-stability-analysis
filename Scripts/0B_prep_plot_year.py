# -*- coding: utf-8 -*-
"""
build_plot_year_wide.py

Generates plot x year wide-format dataset from canonical cutting-level data.
Includes physical anomaly correction (swapping FW and DW for 2002 if DW > FW).
"""

from pathlib import Path
import pandas as pd

# Define relative paths for data availability and reproducibility
BASE_DIR = Path(__file__).resolve().parent
INPUT_PATH = BASE_DIR / "output" / "canonical_cutting_wide_v2.csv"
OUTPUT_PATH = BASE_DIR / "output" / "plot_year_wide.csv"

# Load canonical cutting-level dataset
df = pd.read_csv(INPUT_PATH)

# Correct physical anomaly for 2002 data where DW > FW (physically impossible)
swap_mask = (
    (df["calendar_year"] == 2002)
    & df["fresh_weight_kg5m2"].notna()
    & df["dry_weight_kg5m2"].notna()
    & (df["dry_weight_kg5m2"] > df["fresh_weight_kg5m2"])
)

# Swap FW and DW values where anomaly is detected
fw_temp = df.loc[swap_mask, "fresh_weight_kg5m2"]
df.loc[swap_mask, "fresh_weight_kg5m2"] = df.loc[swap_mask, "dry_weight_kg5m2"]
df.loc[swap_mask, "dry_weight_kg5m2"] = fw_temp
print(f"Corrected FW/DW swapped rows for 2002: {int(swap_mask.sum())}")

# Aggregate cutting-level metrics to plot x year level using standardized genotype_ID
py = (
    df.groupby(
        [
            "plot",
            "genotype_ID",
            "cohort",
            "establishment_year",
            "calendar_year",
            "stand_age",
        ],
        dropna=False,
    )[["summer_height_cm", "fresh_weight_kg5m2", "dry_weight_kg5m2"]]
    .mean()
    .reset_index()
)

# Rename trait columns for concise presentation
py = py.rename(
    columns={
        "summer_height_cm": "SummerPH",
        "fresh_weight_kg5m2": "FW",
        "dry_weight_kg5m2": "DW",
    }
)

# Export plot x year wide-format dataset
py.to_csv(OUTPUT_PATH, index=False, na_rep="")

# Print summary
print(f"Total plot x year rows: {len(py)}")
print(f"Unique genotypes count: {py['genotype_ID'].nunique()}")
print(f"Calendar years covered: {sorted(py['calendar_year'].unique())}")
print(
    "Non-null trait counts - SummerPH: %d, FW: %d, DW: %d"
    % (py["SummerPH"].notna().sum(), py["FW"].notna().sum(), py["DW"].notna().sum())
)
