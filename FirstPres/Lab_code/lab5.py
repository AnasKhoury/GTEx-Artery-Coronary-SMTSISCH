import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA

# =========================
# 0. FILE PATHS – EDIT IF NEEDED
# =========================
expr_path = r"C:\Users\97250\Downloads\Lab4_GTex\expr_artery_coronary_QN_RIN_gt_5_8.csv"
meta_path = r"C:\Users\97250\Downloads\Lab4_GTex\meta_artery_coronary_RIN_gt_5_8.csv"
out_dir = Path(r"C:\Users\97250\Downloads\Lab5_PCA")  # new folder for this lab
out_dir.mkdir(parents=True, exist_ok=True)

# =========================
# 1. LOAD EXPRESSION + META
# =========================
print("Loading expression matrix...")
expr = pd.read_csv(expr_path, index_col=0)   # rows = genes, columns = samples
print("Expression shape (genes x samples):", expr.shape)

print("Loading metadata...")
meta = pd.read_csv(meta_path)
print("Metadata shape:", meta.shape)

# Make sure SAMPID is the index of meta
if "SAMPID" not in meta.columns:
    raise ValueError("SAMPID column not found in metadata!")
meta = meta.set_index("SAMPID")

# Keep only numeric expression columns, drop any 'id' column if exists
if "id" in expr.columns:
    expr = expr.drop(columns=["id"])
expr = expr.select_dtypes(include=[np.number])

# Make sure samples match between expression and metadata
expr_samples = set(expr.columns)
meta_samples = set(meta.index)

common_samples = sorted(expr_samples & meta_samples)
print("Number of common samples:", len(common_samples))

if len(common_samples) == 0:
    raise ValueError("No overlapping samples between expression and meta!")

# Reorder both to the same sample order
expr = expr[common_samples]          # genes x samples
meta = meta.loc[common_samples, :]   # samples x traits

print("Aligned expression shape:", expr.shape)
print("Aligned meta shape:", meta.shape)

# =========================
# 2. PCA (DIMENSION REDUCTION)
# =========================
# PCA works on matrix (samples x features), so transpose:
X = expr.T  # shape: (n_samples, n_genes)
print("Matrix for PCA (samples x genes):", X.shape)

# Compute first 10 PCs (more than enough)
n_components = 10
pca = PCA(n_components=n_components)
X_pca = pca.fit_transform(X)  # shape: (n_samples, n_components)

# Make a DataFrame with PCs
pc_cols = [f"PC{i+1}" for i in range(n_components)]
pcs = pd.DataFrame(X_pca, index=meta.index, columns=pc_cols)
pcs.to_csv(out_dir / "PCA_scores.csv")
print("Saved PCA scores to:", out_dir / "PCA_scores.csv")

# Explained variance
var_ratio = pca.explained_variance_ratio_
print("\nExplained variance ratio:")
for i in range(5):
    print(f"PC{i+1}: {var_ratio[i]*100:.2f}%")

print(f"\nTotal variance explained by first 5 PCs: {var_ratio[:5].sum()*100:.2f}%")

# =========================
# 3. BASIC PC1 vs PC2 SCATTER (NO COLOR)
# =========================
plt.figure(figsize=(7, 6))
plt.scatter(pcs["PC1"], pcs["PC2"], alpha=0.7)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PC1 vs PC2 (all samples)")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(out_dir / "PC1_vs_PC2_basic.png", dpi=150)
plt.close()
print("Saved basic PC1 vs PC2 plot.")

# =========================
# 4. COLOR BY META TRAITS (ONE PLOT PER TRAIT)
# =========================
traits = ["SMGEBTCH", "SMNABTCH", "SMTSISCH", "SMRIN"]

for trait in traits:
    if trait not in meta.columns:
        print(f"[Warning] Trait '{trait}' not found in metadata. Skipping.")
        continue

    data = pd.concat([pcs[["PC1", "PC2"]], meta[trait]], axis=1)

    # Drop samples with missing values for this trait
    data = data.dropna(subset=[trait])
    if data.empty:
        print(f"[Warning] Trait '{trait}' has only NaNs after filtering. Skipping.")
        continue

    # If trait is non-numeric (categorical), encode as integers for coloring
    if not np.issubdtype(data[trait].dtype, np.number):
        data[trait] = data[trait].astype("category").cat.codes

    plt.figure(figsize=(7, 6))
    sc = plt.scatter(
        data["PC1"],
        data["PC2"],
        c=data[trait],
        cmap="viridis",
        alpha=0.8
    )
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(f"PC1 vs PC2 colored by {trait}")
    cbar = plt.colorbar(sc)
    cbar.set_label(trait)
    plt.grid(alpha=0.3)
    plt.tight_layout()

    out_file = out_dir / f"PC1_PC2_colored_by_{trait}.png"
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f"Saved PC1 vs PC2 plot colored by {trait} to:", out_file)

print("\nAll PCA plots and results saved in:", out_dir)
