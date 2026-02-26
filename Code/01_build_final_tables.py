# 01_build_final_tables.py
# Builds:
#   - expr_final.csv (genes x samples)
#   - meta_final.csv (samples x attributes)
# by aligning sample order and merging subject phenotypes (AGE/SEX/DeathType).

import pandas as pd
from pathlib import Path

IN_EXPR = Path("expression_after_outlier_removal.csv")
IN_META = Path("meta_after_outlier_removal.csv")
IN_SUBJ = Path("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

OUT_EXPR = Path("expr_final.csv")
OUT_META = Path("meta_final.csv")
OUT_SUMMARY = Path("step1_build_final_summary.txt")

# Load
expr = pd.read_csv(IN_EXPR, index_col=0)   # genes x samples
meta = pd.read_csv(IN_META)               # samples x attributes

# Extract subject ID from SAMPID (e.g., GTEX-1117F-0626-SM-xxxx -> GTEX-1117F)
def sampid_to_subjid(s: str) -> str:
    parts = str(s).split("-")
    return "-".join(parts[:2]) if len(parts) >= 2 else str(s)

meta["SUBJID"] = meta["SAMPID"].map(sampid_to_subjid)

# Subject phenotypes (keep AGE/SEX/DTHHRDY as DeathType)
subj = pd.read_csv(IN_SUBJ, sep="\t", dtype=str)
subj = subj[["SUBJID", "SEX", "AGE", "DTHHRDY"]].copy()
subj = subj.rename(columns={"DTHHRDY": "DeathType"})

# Merge subject phenotypes onto sample-level meta
meta_m = meta.merge(subj, on="SUBJID", how="left", validate="many_to_one")

# Keep required fields (+ tissue context is helpful for the report)
keep_cols = [
    "SAMPID", "SUBJID",
    "SMTS", "SMTSD",
    "SMTSISCH", "SMRIN",
    "SMNABTCH", "SMGEBTCH",
    "SMNABTCHD", "SMGEBTCHD",
    "SMNABTCHT", "SMGEBTCHT",
    "AGE", "SEX", "DeathType"
]
keep_cols = [c for c in keep_cols if c in meta_m.columns]
meta_final = meta_m[keep_cols].copy()

# Numeric coercions (keep NaN for missing)
for c in ["SMTSISCH", "SMRIN"]:
    if c in meta_final.columns:
        meta_final[c] = pd.to_numeric(meta_final[c], errors="coerce")

if "SEX" in meta_final.columns:
    meta_final["SEX"] = pd.to_numeric(meta_final["SEX"], errors="coerce").astype("Int64")
    meta_final["SEX_label"] = meta_final["SEX"].map({1: "Male", 2: "Female"}).astype("string")

if "DeathType" in meta_final.columns:
    meta_final["DeathType"] = pd.to_numeric(meta_final["DeathType"], errors="coerce").astype("Int64")

# Align meta order to expr columns (and then force expr columns to match)
meta_final = (
    meta_final.set_index("SAMPID")
    .loc[expr.columns]
    .reset_index()
    .rename(columns={"index": "SAMPID"})  # safety
)

expr_final = expr.loc[:, meta_final["SAMPID"].tolist()]
assert list(expr_final.columns) == meta_final["SAMPID"].tolist(), "Alignment failed."

# Write outputs
expr_final.to_csv(OUT_EXPR)
meta_final.to_csv(OUT_META, index=False)

# Summary QC
missing_cols = [c for c in ["AGE", "SEX", "DeathType", "SMTSISCH", "SMRIN"] if c in meta_final.columns]
missing = meta_final[missing_cols].isna().sum()

lines = []
lines.append("STEP 1 — Build expr_final.csv + meta_final.csv")
lines.append(f"Input expr: {IN_EXPR.name}  shape={expr.shape} (genes x samples)")
lines.append(f"Input meta: {IN_META.name}  shape={meta.shape} (samples x attributes)")
lines.append(f"Output expr_final: {OUT_EXPR.name}  shape={expr_final.shape}")
lines.append(f"Output meta_final: {OUT_META.name}  shape={meta_final.shape}")
lines.append("")
lines.append("Missing values in key fields (after merge/alignment):")
for k, v in missing.to_dict().items():
    lines.append(f"  {k}: {int(v)}")
lines.append("")
lines.append("Head(meta_final):")
lines.append(meta_final.head(5).to_string(index=False))

OUT_SUMMARY.write_text("\n".join(lines), encoding="utf-8")
print(f"Wrote: {OUT_EXPR}, {OUT_META}, {OUT_SUMMARY}")
