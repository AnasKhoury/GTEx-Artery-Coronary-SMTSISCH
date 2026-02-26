import pandas as pd
import numpy as np
from sklearn.preprocessing import quantile_transform
from pathlib import Path

# ================================
# 1. Quantile normalization (genes)
# ================================

expr_path = r"C:\Users\97250\Downloads\Lab2_GTex\artery_coronary_filtered_log2.csv"
out_dir = Path(r"C:\Users\97250\Downloads\Lab4_GTex")
out_dir.mkdir(parents=True, exist_ok=True)

# Load filtered + log2 expression matrix from Lab 2
expr = pd.read_csv(expr_path, index_col=0)
print("Expression (log2-filtered) shape:", expr.shape)

# Drop non-expression columns if exist (e.g. 'id')
if "id" in expr.columns:
    expr = expr.drop(columns=["id"])

# Keep only numeric columns
expr = expr.select_dtypes(include=[np.number])
print("Expression (numeric only) shape:", expr.shape)

# Quantile normalization across samples
nq = min(1000, expr.shape[0])  # n_quantiles ≤ number of genes
expr_qn_vals = quantile_transform(
    expr.values,
    axis=0,
    copy=True,
    n_quantiles=nq,
    output_distribution="normal",
    subsample=int(1e9)
)
expr_qn = pd.DataFrame(expr_qn_vals, index=expr.index, columns=expr.columns)

# Save normalized matrix
expr_qn_path = out_dir / "artery_coronary_quantile_normalized.csv"
expr_qn.to_csv(expr_qn_path)
print("Quantile-normalized matrix saved to:", expr_qn_path)
print("Normalized shape:", expr_qn.shape)

# ============================================
# 2–5. Work with meta-data (sample attributes)
# ============================================

meta_path = r"C:\Users\97250\Downloads\GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
meta = pd.read_csv(meta_path, sep="\t")
print("\nMeta-data shape (all tissues):", meta.shape)

# 3. Keep meta-data only for chosen tissue (SMTSD)
tissue_name = "Artery - Coronary"
meta_tissue = meta[meta["SMTSD"] == tissue_name]
print("Meta-data for tissue", tissue_name, ":", meta_tissue.shape[0], "samples")

# 4. Keep only samples that have gene-expression values
expr_samples = set(expr_qn.columns)
meta_tissue_matched = meta_tissue[meta_tissue["SAMPID"].isin(expr_samples)]
print("Samples after matching meta-data with expression:",
      meta_tissue_matched.shape[0])

# (optional) save this matched meta-data
matched_meta_path = out_dir / "meta_artery_coronary_matched.csv"
meta_tissue_matched.to_csv(matched_meta_path, index=False)
print("Matched meta-data saved to:", matched_meta_path)

# 5. Filter samples with RIN > 5.8
meta_final = meta_tissue_matched[meta_tissue_matched["SMRIN"] > 5.8]
print("Samples after RIN > 5.8 filtering:", meta_final.shape[0])

# save final meta-data
final_meta_path = out_dir / "meta_artery_coronary_RIN_gt_5_8.csv"
meta_final.to_csv(final_meta_path, index=False)
print("Final meta-data saved to:", final_meta_path)

# (optional) filter expression matrix to the same samples
final_ids = list(meta_final["SAMPID"])
expr_qn_final = expr_qn[final_ids]
expr_qn_final_path = out_dir / "expr_artery_coronary_QN_RIN_gt_5_8.csv"
expr_qn_final.to_csv(expr_qn_final_path)
print("Final expression matrix saved to:", expr_qn_final_path)
print("Final expression shape:", expr_qn_final.shape)
