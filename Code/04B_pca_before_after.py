# 04B_pca_before_after.py
# PCA before/after R ComBat, colored by batch (SMGEBTCH)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA

IN_META = Path("meta_final.csv")
IN_EXPR_BEFORE = Path("expr_final.csv")
IN_EXPR_AFTER  = Path("expr_combat.csv")

BATCH_COL = "SMGEBTCH"

OUT_PCA_BEFORE = Path("pca_before.csv")
OUT_PCA_AFTER  = Path("pca_after.csv")
OUT_FIG_BEFORE = Path("pca_before_by_batch.png")
OUT_FIG_AFTER  = Path("pca_after_by_batch.png")

meta = pd.read_csv(IN_META)

# --- Detect the sample-id column (expects GTEX-...-SM-...)
def detect_sample_col(df: pd.DataFrame) -> str:
    candidates = ["SAMPID", "sample_id", "SampleID", "index", "Unnamed: 0"]
    for c in candidates:
        if c in df.columns:
            return c
    # fallback: find a column that looks like GTEx sample IDs
    for c in df.columns:
        s = df[c].astype(str)
        if s.str.contains(r"^GTEX-[A-Z0-9]+-\d{4}-SM-", regex=True, na=False).mean() > 0.5:
            return c
    raise ValueError(f"Could not find sample-id column. Columns are: {list(df.columns)}")

SAMP_COL = detect_sample_col(meta)
print("Using sample-id column:", SAMP_COL)

def load_expr(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)  # genes x samples

def run_pca(expr_gxs: pd.DataFrame, meta_df: pd.DataFrame, out_csv: Path):
    # align meta to expression columns
    meta2 = meta_df.set_index(SAMP_COL).loc[expr_gxs.columns].reset_index()

    # After reset_index(), the sample-id column might be called "index"
    # Force the first column to be named "SAMPID"
    meta2 = meta2.rename(columns={meta2.columns[0]: "SAMPID"})

    X = expr_gxs.T.to_numpy(dtype=float)  # samples x genes

    if not np.isfinite(X).all():
        raise ValueError("Non-finite values found in PCA input. Check expr_combat.csv.")

    pca = PCA(n_components=2, random_state=0)
    Z = pca.fit_transform(X)

    out = pd.DataFrame({
        "SAMPID": meta2["SAMPID"].values,
        "PC1": Z[:, 0],
        "PC2": Z[:, 1],
        "batch": meta2[BATCH_COL].astype(str).values
    })
    out["PC1_var"] = pca.explained_variance_ratio_[0]
    out["PC2_var"] = pca.explained_variance_ratio_[1]
    out.to_csv(out_csv, index=False)
    return out

def plot(df: pd.DataFrame, title: str, out_png: Path):
    fig, ax = plt.subplots(figsize=(8, 6))
    batches = df["batch"].unique()
    batch_to_i = {b: i for i, b in enumerate(batches)}
    colors = df["batch"].map(batch_to_i)

    sca = ax.scatter(df["PC1"], df["PC2"], c=colors, s=20, alpha=0.8)
    ax.set_xlabel(f"PC1 (var={df['PC1_var'].iloc[0]:.3f})")
    ax.set_ylabel(f"PC2 (var={df['PC2_var'].iloc[0]:.3f})")
    ax.set_title(title)
    ax.grid(True, linewidth=0.3, alpha=0.5)

    cbar = plt.colorbar(sca, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Batch index (SMGEBTCH)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

expr_before = load_expr(IN_EXPR_BEFORE)
expr_after  = load_expr(IN_EXPR_AFTER)

# Use only samples with batch present
meta_plot = meta.dropna(subset=[BATCH_COL]).copy()

# Filter expression columns to match meta_plot sample list
expr_before = expr_before.loc[:, meta_plot[SAMP_COL]]
expr_after  = expr_after.loc[:, meta_plot[SAMP_COL]]

p1 = run_pca(expr_before, meta_plot, OUT_PCA_BEFORE)
plot(p1, "PCA BEFORE ComBat (colored by SMGEBTCH)", OUT_FIG_BEFORE)

p2 = run_pca(expr_after, meta_plot, OUT_PCA_AFTER)
plot(p2, "PCA AFTER ComBat (colored by SMGEBTCH)", OUT_FIG_AFTER)

print("Wrote:", OUT_FIG_BEFORE.resolve())
print("Wrote:", OUT_FIG_AFTER.resolve())
print("Wrote:", OUT_FIG_AFTER.resolve())
