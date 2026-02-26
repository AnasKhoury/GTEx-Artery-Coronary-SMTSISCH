import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from pathlib import Path

base_dir = Path(r"C:\Users\97250\Downloads\Lab7_confounded_corrected")
meta_path = r"C:\Users\97250\Downloads\Lab4_GTex\meta_artery_coronary_RIN_gt_5_8.csv"

# load matrices
expr_std = pd.read_csv(base_dir / "expr_standardized.csv", index_col=0)
expr_resid = pd.read_csv(base_dir / "expr_residuals_corrected.csv", index_col=0)
meta = pd.read_csv(meta_path)
meta = meta.set_index("SAMPID")

meta_clean = meta.loc[expr_std.columns]

# ---------- helper functions ----------

def plot_pca_continuous(df, pca_obj, color_col, title, out_name):
    plt.figure(figsize=(7,6))
    sc = plt.scatter(df["PC1"], df["PC2"], c=df[color_col],
                     cmap="viridis", edgecolor="k", s=50)
    plt.colorbar(sc, label=color_col)
    var1 = pca_obj.explained_variance_ratio_[0] * 100
    var2 = pca_obj.explained_variance_ratio_[1] * 100
    plt.xlabel(f"PC1 ({var1:.1f}% variance)")
    plt.ylabel(f"PC2 ({var2:.1f}% variance)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(base_dir / out_name, dpi=300)
    plt.close()

def plot_pca_categorical(df, pca_obj, cat_col, title, out_name):
    plt.figure(figsize=(7,6))
    cats = df[cat_col].astype(str).fillna("NA").unique()
    for c in cats:
        mask = df[cat_col].astype(str).fillna("NA") == c
        plt.scatter(df.loc[mask,"PC1"], df.loc[mask,"PC2"],
                    label=c, s=50, edgecolor="k")
    var1 = pca_obj.explained_variance_ratio_[0] * 100
    var2 = pca_obj.explained_variance_ratio_[1] * 100
    plt.xlabel(f"PC1 ({var1:.1f}% variance)")
    plt.ylabel(f"PC2 ({var2:.1f}% variance)")
    plt.title(title)
    plt.legend(title=cat_col, bbox_to_anchor=(1.05,1), loc="upper left")
    plt.tight_layout()
    plt.savefig(base_dir / out_name, dpi=300)
    plt.close()

# ---------- PCA BEFORE (expr_std) ----------

pca_before = PCA(n_components=2)
pcs_before = pca_before.fit_transform(expr_std.T)

pc_before = pd.DataFrame(
    pcs_before,
    index=expr_std.columns,
    columns=["PC1","PC2"]
)
pc_before["SMRIN"] = meta_clean["SMRIN"]
pc_before["SMTSISCH"] = meta_clean["SMTSISCH"]
pc_before["SMGEBTCH"] = meta_clean["SMGEBTCH"]

plot_pca_continuous(pc_before, pca_before, "SMRIN",
                    "PCA Before Confounder Correction (RIN)",
                    "PCA_before_RIN.png")
plot_pca_continuous(pc_before, pca_before, "SMTSISCH",
                    "PCA Before Confounder Correction (Ischemic Time)",
                    "PCA_before_IschemicTime.png")
plot_pca_categorical(pc_before, pca_before, "SMGEBTCH",
                     "PCA Before Confounder Correction (Batch SMGEBTCH)",
                     "PCA_before_Batch.png")

# ---------- PCA AFTER (expr_resid) ----------

pca_after = PCA(n_components=2)
pcs_after = pca_after.fit_transform(expr_resid.T)

pc_after = pd.DataFrame(
    pcs_after,
    index=expr_resid.columns,
    columns=["PC1","PC2"]
)
pc_after["SMRIN"] = meta_clean["SMRIN"]
pc_after["SMTSISCH"] = meta_clean["SMTSISCH"]
pc_after["SMGEBTCH"] = meta_clean["SMGEBTCH"]

plot_pca_continuous(pc_after, pca_after, "SMRIN",
                    "PCA After Confounder Correction (RIN)",
                    "PCA_after_RIN.png")
plot_pca_continuous(pc_after, pca_after, "SMTSISCH",
                    "PCA After Confounder Correction (Ischemic Time)",
                    "PCA_after_IschemicTime.png")
plot_pca_categorical(pc_after, pca_after, "SMGEBTCH",
                     "PCA After Confounder Correction (Batch SMGEBTCH)",
                     "PCA_after_Batch.png")

print("Saved all PCA plots to", base_dir)
