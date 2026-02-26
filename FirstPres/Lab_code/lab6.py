# ======================================================
# Lab 6 — Outlier Detection + Cleaning + PCA (GTEx)
# Continues from Lab5 outputs (Artery - Coronary)
# ======================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from collections import Counter

# ------------------------------------------------------
# 0. File paths (same structure as Lab4 / Lab5)
# ------------------------------------------------------
expr_path = r"C:\Users\97250\Downloads\Lab4_GTex\expr_artery_coronary_QN_RIN_gt_5_8.csv"
meta_path = r"C:\Users\97250\Downloads\Lab4_GTex\meta_artery_coronary_RIN_gt_5_8.csv"
out_dir = Path(r"C:\Users\97250\Downloads\Lab6_Outliers")
out_dir.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------
# 1. Load expression + meta and align samples
# ------------------------------------------------------
expr = pd.read_csv(expr_path, index_col=0)      # genes × samples
meta = pd.read_csv(meta_path).set_index("SAMPID")

# keep only overlapping samples
common_samples = [s for s in expr.columns if s in meta.index]
expr = expr[common_samples]
meta = meta.loc[common_samples]

print("Expression shape (genes x samples):", expr.shape)
print("Meta shape (samples x fields):", meta.shape)

# transpose for sample-based operations
X = expr.T   # samples × genes

# ------------------------------------------------------
# 2. Outlier Method 1 — PCA-based sample outliers
# ------------------------------------------------------
X_scaled = StandardScaler().fit_transform(X)

pca = PCA(n_components=5, random_state=0)
scores = pca.fit_transform(X_scaled)

pc_df = pd.DataFrame(scores[:, :2], index=X.index, columns=["PC1", "PC2"])

# Z-scores of PC1 and PC2 across samples
pc_z = (pc_df - pc_df.mean()) / pc_df.std(ddof=0)

# Mark outliers as |z| > 3 on PC1 or PC2
mask_pca = (pc_z.abs() > 3).any(axis=1)
outliers_pca = pc_df.index[mask_pca].tolist()

print("\nPCA-based outliers (|z|>3 on PC1/PC2):", outliers_pca)

# Cleaner PCA plot
plt.figure(figsize=(7.5, 6))
is_out = pc_df.index.isin(outliers_pca)

plt.scatter(pc_df.loc[~is_out, "PC1"], pc_df.loc[~is_out, "PC2"], alpha=0.7, label="Non-outliers")
plt.scatter(pc_df.loc[is_out, "PC1"],  pc_df.loc[is_out, "PC2"],  alpha=0.9, label="Outlier")

for s in pc_df.index[is_out]:
    plt.text(pc_df.loc[s, "PC1"], pc_df.loc[s, "PC2"], s, fontsize=7)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA Outlier Detection (PC1/PC2 z-score > 3)")
plt.legend()
plt.tight_layout()
plt.savefig(out_dir / "pca_outliers.png", dpi=200)
plt.close()

# --- Hierarchical clustering outliers (teacher-style: extreme merge distance) ---
dist_corr = pdist(X, metric="correlation")
Z = linkage(dist_corr, method="average")
# ------------------------------------------------------
# FULL dendrogram (all samples) — for completeness
# ------------------------------------------------------
plt.figure(figsize=(18, 7))
dendrogram(
    Z,
    no_labels=True,            # hide labels to keep it readable
    color_threshold=0,
    above_threshold_color="black"
)
plt.title("Hierarchical Clustering of Samples (correlation distance, average linkage) — FULL")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig(out_dir / "dendrogram_full_no_labels.png", dpi=300)
plt.close()


# Readable dendrogram (truncated)
plt.figure(figsize=(12, 5))
dendrogram(Z, truncate_mode="lastp", p=35, show_leaf_counts=True, no_labels=True)
plt.title("Hierarchical Clustering (correlation distance, average linkage) — truncated view")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig(out_dir / "dendrogram_truncated.png", dpi=250)
plt.close()

# Teacher-style outliers: cut very high, keep only the most separated branches
# Try 0.99 (top 1%) first; if it returns 0, use 0.98
cut_h = float(np.quantile(Z[:, 2], 0.99))
clusters = fcluster(Z, t=cut_h, criterion="distance")
cluster_series = pd.Series(clusters, index=X.index)
cluster_sizes = cluster_series.value_counts()

# Outliers = samples in clusters that are very small at this extreme cut
outliers_hc = cluster_series[cluster_series.map(cluster_sizes) <= 2].index.tolist()

print("\nHierarchical outliers (extreme cut, teacher-style):", outliers_hc)
print("Hierarchical cut height used:", cut_h)


# Also save which cluster each sample got (for debugging/report)
hc_assign = pd.DataFrame({
    "SAMPID": cluster_series.index,
    "cluster": cluster_series.values,
    "cluster_size": cluster_series.map(cluster_sizes).values
}).set_index("SAMPID")
hc_assign.to_csv(out_dir / "hierarchical_clusters_assignment.csv")

# ------------------------------------------------------
# 4. Outlier Method 3 — Total expression (per sample, IQR)
# ------------------------------------------------------
# NOTE: You said you want ONLY 2 outlier methods (PCA + Hierarchical Clustering).
# We keep this block in the script for completeness, but we do NOT use it in the final union.
# If you want it fully disabled, set USE_TOTAL_EXPR_METHOD = False.

USE_TOTAL_EXPR_METHOD = False  # <-- set True if you ever want to include the 3rd method again

total_expr = expr.sum(axis=0)  # sum across genes → one value per sample

q1, q3 = total_expr.quantile([0.25, 0.75])
iqr = q3 - q1
lower = q1 - 1.5 * iqr
upper = q3 + 1.5 * iqr

outliers_total = total_expr[(total_expr < lower) | (total_expr > upper)].index.tolist()

print("\nTotal-expression outliers (IQR rule):", outliers_total)

# Boxplot with clear x-axis label (no confusing "1")
plt.figure(figsize=(6.5, 4))
plt.boxplot(total_expr.values, tick_labels=["All samples"])
plt.title("Total Expression per Sample (IQR rule)")
plt.ylabel("Total expression (sum across genes)")
plt.tight_layout()
plt.savefig(out_dir / "total_expression_boxplot.png", dpi=200)
plt.close()

# ------------------------------------------------------
# 5. FINAL OUTLIERS = UNION of selected methods
# ------------------------------------------------------
# IMPORTANT: Only PCA + Hierarchical are included.
# Total-expression method is excluded unless USE_TOTAL_EXPR_METHOD=True.
if USE_TOTAL_EXPR_METHOD:
    final_outliers = sorted(set(outliers_pca) |
                            set(outliers_hc) |
                            set(outliers_total))
    counts = Counter(outliers_pca + outliers_hc + outliers_total)
else:
    final_outliers = sorted(set(outliers_pca) |
                            set(outliers_hc))
    counts = Counter(outliers_pca + outliers_hc)

print("\nCounts per sample (how many methods flagged it):")
for s, c in counts.items():
    print(f"{s}: {c} method(s)")

print("\nFINAL REMOVED OUTLIERS (union of methods):")
for s in final_outliers:
    print(" -", s)

pd.Series(final_outliers, name="outlier_sample").to_csv(
    out_dir / "final_outliers_union.csv", index=False
)

# ------------------------------------------------------
# 6. Remove outliers from expression + meta
# ------------------------------------------------------
expr_clean = expr.drop(columns=final_outliers)
meta_clean = meta.drop(index=final_outliers)

print("\nAfter removing outliers:")
print("Expression shape (genes x samples):", expr_clean.shape)
print("Meta shape (samples x fields):", meta_clean.shape)

expr_clean.to_csv(out_dir / "expression_after_outlier_removal.csv")
meta_clean.to_csv(out_dir / "meta_after_outlier_removal.csv")

# ------------------------------------------------------
# 7. PCA AFTER cleaning (as required by Lab6)
# ------------------------------------------------------
X2 = expr_clean.T
X2_scaled = StandardScaler().fit_transform(X2)

pca2 = PCA(n_components=5, random_state=0)
pc2 = pca2.fit_transform(X2_scaled)

pc2_df = pd.DataFrame(pc2[:, :2], index=X2.index, columns=["PC1", "PC2"])
explained = pca2.explained_variance_ratio_[:2] * 100

plt.figure(figsize=(7.5, 6))
plt.scatter(pc2_df["PC1"], pc2_df["PC2"], alpha=0.75)
plt.xlabel(f"PC1 ({explained[0]:.2f}%)")
plt.ylabel(f"PC2 ({explained[1]:.2f}%)")
plt.title("PCA AFTER outlier removal")
plt.tight_layout()
plt.savefig(out_dir / "pca_after_outlier_removal.png", dpi=200)
plt.close()

# ------------------------------------------------------
# 8. Hierarchical clustering AFTER cleaning (recompute!)
# ------------------------------------------------------
X_clean = expr_clean.T
dist_corr2 = pdist(X_clean, metric="correlation")
Z2 = linkage(dist_corr2, method="average")

plt.figure(figsize=(12, 5))
dendrogram(
    Z2,
    truncate_mode="lastp",
    p=35,
    show_leaf_counts=True,
    no_labels=True
)
plt.title("Hierarchical Clustering AFTER outlier removal — truncated view")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig(out_dir / "dendrogram_after_removal.png", dpi=250)
plt.close()

print("\nAll Lab6 results saved in:", out_dir)
