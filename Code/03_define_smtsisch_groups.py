# 03_define_smtsisch_groups.py
# Creates LOW vs HIGH SMTSISCH groups using Q1 vs Q3 (drops middle 50%).
# Outputs: meta_groups.csv, group_stats_table.csv, smtsisch_thresholds.txt

import pandas as pd
import numpy as np
from pathlib import Path

IN_META = Path("meta_final.csv")

OUT_META_GROUPS = Path("meta_groups.csv")
OUT_STATS = Path("group_stats_table.csv")
OUT_THRESH = Path("smtsisch_thresholds.txt")
OUT_HIST = Path("smtsisch_histogram.csv")

meta = pd.read_csv(IN_META)

# Required columns
needed = ["SAMPID", "SMTSISCH", "SMRIN", "AGE", "SEX", "DeathType", "SMGEBTCH"]
missing = [c for c in needed if c not in meta.columns]
if missing:
    raise ValueError(f"Missing columns in meta_final.csv: {missing}")

# Ensure numeric
meta["SMTSISCH"] = pd.to_numeric(meta["SMTSISCH"], errors="coerce")
meta["SMRIN"] = pd.to_numeric(meta["SMRIN"], errors="coerce")

# Compute Q1/Q3 on non-missing values
x = meta["SMTSISCH"].dropna()
if x.shape[0] < 20:
    raise ValueError("Too few non-missing SMTSISCH values to form Q1/Q3 groups.")

q1 = float(x.quantile(0.25))
q3 = float(x.quantile(0.75))

# Define groups
meta["Group"] = pd.NA
meta.loc[meta["SMTSISCH"] <= q1, "Group"] = "LOW"
meta.loc[meta["SMTSISCH"] >= q3, "Group"] = "HIGH"

meta_g = meta.dropna(subset=["Group"]).copy()

# Save thresholds used (for report reproducibility)
OUT_THRESH.write_text(
    f"SMTSISCH thresholds (Q1 vs Q3)\n"
    f"n_total={meta.shape[0]}\n"
    f"n_nonmissing_SMTSISCH={x.shape[0]}\n"
    f"Q1 (25%) = {q1}\n"
    f"Q3 (75%) = {q3}\n"
    f"LOW:  SMTSISCH <= Q1\n"
    f"HIGH: SMTSISCH >= Q3\n"
    f"n_LOW={int((meta_g['Group']=='LOW').sum())}\n"
    f"n_HIGH={int((meta_g['Group']=='HIGH').sum())}\n",
    encoding="utf-8"
)

# Simple histogram bins (for sanity check)
bins = np.histogram_bin_edges(x.to_numpy(), bins=20)
hist_counts, hist_bins = np.histogram(x.to_numpy(), bins=bins)
hist_df = pd.DataFrame({
    "bin_left": hist_bins[:-1],
    "bin_right": hist_bins[1:],
    "count": hist_counts
})
hist_df.to_csv(OUT_HIST, index=False)

# Build a clean group stats table
def summarize_group(df: pd.DataFrame, name: str):
    n = df.shape[0]
    out = {
        "Group": name,
        "n": n,
        "SMTSISCH_mean": df["SMTSISCH"].mean(),
        "SMTSISCH_median": df["SMTSISCH"].median(),
        "SMTSISCH_min": df["SMTSISCH"].min(),
        "SMTSISCH_max": df["SMTSISCH"].max(),
        "SMRIN_mean": df["SMRIN"].mean(),
        "SMRIN_median": df["SMRIN"].median(),
    }

    # AGE in GTEx is categorical bins like "50-59"
    if "AGE" in df.columns:
        out["AGE_mode"] = df["AGE"].mode(dropna=True).iloc[0] if df["AGE"].notna().any() else pd.NA

    # SEX coded 1/2 in GTEx subject phenotypes
    if "SEX" in df.columns:
        sex_counts = df["SEX"].value_counts(dropna=False).to_dict()
        out["SEX_male_n"] = int(sex_counts.get(1, 0))
        out["SEX_female_n"] = int(sex_counts.get(2, 0))

    # DeathType (DTHHRDY) 0-4
    if "DeathType" in df.columns:
        d_counts = df["DeathType"].value_counts(dropna=False).to_dict()
        for k in sorted([kk for kk in d_counts.keys() if pd.notna(kk)]):
            out[f"DeathType_{int(k)}_n"] = int(d_counts[k])

    return out

stats = pd.DataFrame([
    summarize_group(meta_g[meta_g["Group"] == "LOW"], "LOW"),
    summarize_group(meta_g[meta_g["Group"] == "HIGH"], "HIGH"),
])

stats.to_csv(OUT_STATS, index=False)

# Save meta_groups.csv (keep everything + Group)
meta_g.to_csv(OUT_META_GROUPS, index=False)

print("Wrote:", OUT_META_GROUPS.resolve())
print("Wrote:", OUT_STATS.resolve())
print("Wrote:", OUT_THRESH.resolve())
print("Wrote:", OUT_HIST.resolve())
print("\nThresholds:")
print(OUT_THRESH.read_text(encoding="utf-8"))
