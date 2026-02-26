# ===== Lab 2 — GTEx (Artery–Coronary) =====
# Steps: load -> log2(x+1) -> histograms -> filtering -> save
# Edit the 2 paths below to match your computer.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---- EDIT THESE PATHS ----
EXPR_CSV = Path(r"C:\Users\97250\Downloads\gene_expression_data.csv")
OUT_DIR  = Path(r"C:\Users\97250\Downloads\Lab2_GTex")  # results folder
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1) Load expression data (rows=gene IDs, cols=sample IDs)
expression_data = pd.read_csv(EXPR_CSV, index_col=0)
print("Shape before log transform:", expression_data.shape)

# 2) Log2(x+1) transform
expression_log = np.log2(expression_data + 1.0)
print("Shape after log transform:", expression_log.shape)

# 3) Histograms (saved as PNGs for Word)
plt.figure()
plt.hist(expression_data.values.flatten(), bins=150)
plt.xlim(0, 10000)  # zoom in to see the distribution better (teacher comment)
plt.title("Distribution BEFORE log2(x+1) (zoomed to 0–10,000 TPM)")
plt.xlabel("TPM")
plt.ylabel("Frequency")
before_png = OUT_DIR / "hist_before_log.png"
plt.savefig(before_png, dpi=150, bbox_inches="tight")
plt.close()

plt.figure()
plt.hist(expression_log.values.flatten(), bins=100)
plt.title("Distribution AFTER log2(x+1)")
plt.xlabel("log2(TPM + 1)")
plt.ylabel("Frequency")
after_png = OUT_DIR / "hist_after_log.png"
plt.savefig(after_png, dpi=150, bbox_inches="tight")
plt.close()

print(f"Saved histograms to:\n  {before_png}\n  {after_png}")

# 4) Filtering: remove genes with variance == 0 and low mean
#    You can adjust thresholds if your instructor specifies others.
MIN_VAR = 0.0           # keep genes with variance > 0
MIN_MEAN_LOG = 1.0      # keep genes with mean log2(TPM+1) > 1

n_before = expression_log.shape[0]
var_by_gene = expression_log.var(axis=1)
mean_by_gene = expression_log.mean(axis=1)
mask = (var_by_gene > MIN_VAR) & (mean_by_gene > MIN_MEAN_LOG)
filtered = expression_log.loc[mask]

n_after = filtered.shape[0]
n_removed = n_before - n_after

print("\n--- Gene counts ---")
print("Genes before filtering:", n_before)
print("Genes after  filtering:", n_after)
print("Genes removed:", n_removed)

# 5) Save cleaned matrix for next lab
clean_csv = OUT_DIR / "artery_coronary_filtered_log2.csv"
filtered.to_csv(clean_csv)
print(f"\nSaved cleaned matrix:\n  {clean_csv}")
