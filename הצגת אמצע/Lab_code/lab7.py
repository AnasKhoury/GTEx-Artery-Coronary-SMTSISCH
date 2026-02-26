
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
from scipy.stats import zscore
import matplotlib.pyplot as plt


from pathlib import Path

# =========================
# 0. PATHS 
# =========================

# Expression matrix (after tissue selection and RIN filtering)
expr_path = r"C:\Users\97250\Downloads\Lab4_GTex\expr_artery_coronary_QN_RIN_gt_5_8.csv"

# Metadata for the same samples
meta_path = r"C:\Users\97250\Downloads\Lab4_GTex\meta_artery_coronary_RIN_gt_5_8.csv"

# Output directory for Lab 7 results
out_dir = Path(r"C:\Users\97250\Downloads\Lab7_confounded_corrected")
out_dir.mkdir(parents=True, exist_ok=True)



# =====================================================
# 1. LOAD EXPRESSION + METADATA
# =====================================================

# Expression: rows = genes, columns = samples (SAMPID)
expr = pd.read_csv(expr_path, index_col=0)

# Metadata: must contain SAMPID column (or already indexed by it)
meta = pd.read_csv(meta_path)

if "SAMPID" in meta.columns:
    meta = meta.set_index("SAMPID")

print("Expression shape (genes x samples):", expr.shape)
print("Metadata shape (samples x variables):", meta.shape)


# =====================================================
# 2. ALIGN SAMPLES BETWEEN expr AND meta
# =====================================================

# We keep only samples that appear in BOTH expression and metadata
common_samples = expr.columns.intersection(meta.index)

expr = expr[common_samples]
meta = meta.loc[common_samples]

print("After aligning:", expr.shape, meta.shape)
# Now the order of samples is consistent between expr and meta.


# =====================================================
# 3. GENE FILTERING: LOW EXPRESSION + LOW VARIANCE
# =====================================================

# Remove genes with very low mean expression or almost no variance.
# You can adjust thresholds to match what you did in Lab 2.
min_mean_expr = 0.1   # threshold for mean expression
min_var_expr = 0.1    # threshold for variance

gene_means = expr.mean(axis=1)
gene_vars = expr.var(axis=1)

keep_genes = (gene_means > min_mean_expr) & (gene_vars > min_var_expr)
expr_filt = expr.loc[keep_genes]

print("Genes before filtering:", expr.shape[0])
print("Genes after filtering:", expr_filt.shape[0])


# =====================================================
# 4. REMOVE OUTLIERS (FROM LAB 6)
# =====================================================

# Outliers detected previously (10 by IQR + 1 by PCA)
outlier_samples = [
    "GTEX-13S7M-2226-SM-5S2WB",
    "GTEX-15DDE-1226-SM-6PAN8",
    "GTEX-16YQH-0826-SM-7938J",
    "GTEX-1HSMQ-0626-SM-CL54A",
    "GTEX-1IDJH-1226-SM-E6CJU",
    "GTEX-PWCY-0426-SM-48TCW",
    "GTEX-RWS6-0426-SM-47JXH",
    "GTEX-WHSE-0726-SM-4RTXA",
    "GTEX-XV7Q-0626-SM-4BRV5",
    "GTEX-ZPIC-0326-SM-DO91U",
    "GTEX-14E7W-1026-SM-62LEK"
]

valid_samples = [s for s in expr_filt.columns if s not in outlier_samples]

expr_clean = expr_filt[valid_samples]
meta_clean = meta.loc[valid_samples]

print("Samples before removing outliers:", expr_filt.shape[1])
print("Samples after removing outliers:", expr_clean.shape[1])

# Optional: save intermediate matrix
expr_clean.to_csv(out_dir / "expr_filtered_no_outliers.csv")


# =====================================================
# 5. QUANTILE NORMALIZATION
# =====================================================

def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform quantile normalization across samples (columns).
    All samples will share the same value distribution.
    """
    # sort each column
    sorted_df = pd.DataFrame(
        np.sort(df.values, axis=0),
        index=df.index,
        columns=df.columns
    )
    # mean across columns for each rank
    rank_mean = sorted_df.mean(axis=1)

    # create a copy to fill
    df_qn = df.copy()
    for col in df.columns:
        # order of original values
        order = np.argsort(df[col].values)
        # new column with same shape
        new_col = np.empty_like(df[col].values, dtype=float)
        # assign rank means according to sorted order
        new_col[order] = rank_mean.values
        df_qn[col] = new_col

    return df_qn

expr_qn = quantile_normalize(expr_clean)
print("Quantile normalization done. Shape:", expr_qn.shape)

expr_qn.to_csv(out_dir / "expr_quantile_normalized.csv")


# =====================================================
# 6. STANDARDIZATION (Z-SCORE PER GENE)
# =====================================================

# For each gene, subtract mean and divide by std across samples.
# This puts all genes on the same scale (mean=0, std=1) – important for PCA and regression.
expr_std_values = zscore(expr_qn.values, axis=1, ddof=1)
expr_std = pd.DataFrame(
    expr_std_values,
    index=expr_qn.index,
    columns=expr_qn.columns
)

print("Standardization done. Any NaNs?", np.isnan(expr_std.values).sum())
expr_std.to_csv(out_dir / "expr_standardized.csv")


# =====================================================
# 7. PCA BEFORE CORRECTION (FOR VISUALIZATION)
# =====================================================

pca_before = PCA(n_components=2)
pcs_before = pca_before.fit_transform(expr_std.T)  # samples as rows

pc_before_df = pd.DataFrame(
    pcs_before,
    index=expr_std.columns,
    columns=["PC1", "PC2"]
)

# Attach metadata columns we care about for coloring
pc_before_df["SMRIN"] = meta_clean.loc[pc_before_df.index, "SMRIN"]
pc_before_df["SMTSISCH"] = meta_clean.loc[pc_before_df.index, "SMTSISCH"]
pc_before_df["SMGEBTCH"] = meta_clean.loc[pc_before_df.index, "SMGEBTCH"]
# Add Death Type (Hardy scale) if exists
if "DTHHRDY" in meta_clean.columns:
    pc_before_df["DTHHRDY"] = meta_clean.loc[pc_before_df.index, "DTHHRDY"]



# =====================================================
# 8. BUILD COVARIATE MATRIX FROM METADATA (NUMERIC SAFE)
# =====================================================

# Continuous covariates
numeric_covariates = ["SMRIN", "SMTSISCH"]
# Categorical covariates (batches + death type)
categorical_covariates = ["SMNABTCH", "SMGEBTCH", "DTHHRDY"]


# Keep only covariates that actually exist in metadata
numeric_covariates = [c for c in numeric_covariates if c in meta_clean.columns]
categorical_covariates = [c for c in categorical_covariates if c in meta_clean.columns]

print("Numeric covariates:", numeric_covariates)
print("Categorical covariates:", categorical_covariates)

# Extract these columns
meta_cov = meta_clean[numeric_covariates + categorical_covariates].copy()

# One-hot encode categorical covariates (batch variables)
if len(categorical_covariates) > 0:
    meta_cov = pd.get_dummies(meta_cov, columns=categorical_covariates, drop_first=True)

# Convert all covariate columns to numeric (in case some are still 'object')
meta_cov = meta_cov.apply(pd.to_numeric, errors="coerce")

# Drop columns that are completely NaN (just in case)
meta_cov = meta_cov.dropna(axis=1, how="all")

# Drop samples with any NaNs in selected covariates
mask_complete = meta_cov.notna().all(axis=1)
meta_cov = meta_cov.loc[mask_complete]

# Match expr_std to these samples
expr_std = expr_std.loc[:, meta_cov.index]

print("After aligning with covariates (clean numeric):")
print("expr_std shape:", expr_std.shape)
print("meta_cov shape:", meta_cov.shape)

# Add intercept term
X = meta_cov.copy()
X["intercept"] = 1.0

# Ensure a consistent column order
X = X[list(meta_cov.columns) + ["intercept"]]

# >>> IMPORTANT: force to float matrix <<<
X_mat = X.to_numpy(dtype=float)  # shape: (n_samples, p)

# Precompute (X'X)^(-1) X' to efficiently get beta for each gene
XtX_inv_Xt = np.linalg.pinv(X_mat.T @ X_mat) @ X_mat.T


# =====================================================
# 9. MULTIPLE REGRESSION PER GENE – RESIDUAL MATRIX
# =====================================================

genes = expr_std.index
samples = expr_std.columns
Y = expr_std.values  # shape: (n_genes, n_samples)

residuals = np.empty_like(Y)

for i, gene in enumerate(genes):
    y = Y[i, :]                 # expression of this gene across samples
    beta_hat = XtX_inv_Xt @ y   # regression coefficients
    y_hat = X_mat @ beta_hat    # fitted values
    residuals[i, :] = y - y_hat # residuals

expr_resid = pd.DataFrame(
    residuals,
    index=genes,
    columns=samples
)

print("Residual matrix shape:", expr_resid.shape)

# This is the main output of Lab 7
expr_resid.to_csv(out_dir / "expr_residuals_corrected.csv")

print("✅ Lab 7 pipeline finished. Residual matrix saved.")

# =====================================================
# 10. PCA AFTER CORRECTION
# =====================================================

pca_after = PCA(n_components=2)
pcs_after = pca_after.fit_transform(expr_resid.T)  # samples as rows

pc_after_df = pd.DataFrame(
    pcs_after,
    index=expr_resid.columns,
    columns=["PC1", "PC2"]
)

pc_after_df["SMRIN"] = meta_clean.loc[pc_after_df.index, "SMRIN"]
pc_after_df["SMTSISCH"] = meta_clean.loc[pc_after_df.index, "SMTSISCH"]
pc_after_df["SMGEBTCH"] = meta_clean.loc[pc_after_df.index, "SMGEBTCH"]
# Add Death Type (Hardy scale) if exists
if "DTHHRDY" in meta_clean.columns:
    pc_after_df["DTHHRDY"] = meta_clean.loc[pc_after_df.index, "DTHHRDY"]



# =====================================================
# 11. PCA PLOTS BEFORE VS AFTER CORRECTION
# =====================================================

def plot_pca_continuous(df, pca_obj, color_col, title, output_name=None):
    """
    Scatter plot of PC1 vs PC2 colored by a continuous covariate.
    """
    plt.figure(figsize=(7, 6))
    sc = plt.scatter(df["PC1"], df["PC2"], c=df[color_col],
                     cmap="viridis", edgecolor="k", s=50)
    plt.colorbar(sc, label=color_col)
    var1 = pca_obj.explained_variance_ratio_[0] * 100
    var2 = pca_obj.explained_variance_ratio_[1] * 100
    plt.xlabel(f"PC1 ({var1:.1f}% variance)")
    plt.ylabel(f"PC2 ({var2:.1f}% variance)")
    plt.title(title)
    plt.tight_layout()
    if output_name is not None:
        plt.savefig(out_dir / output_name, dpi=300)
    plt.show()


def plot_pca_categorical(df, pca_obj, cat_col, title, output_name=None):
    """
    Scatter plot of PC1 vs PC2 colored by a categorical covariate.
    """
    plt.figure(figsize=(7, 6))
    categories = df[cat_col].astype(str).fillna("NA").unique()
    for cat in categories:
        mask = df[cat_col].astype(str).fillna("NA") == cat
        plt.scatter(df.loc[mask, "PC1"],
                    df.loc[mask, "PC2"],
                    label=cat,
                    s=50,
                    edgecolor="k")
    var1 = pca_obj.explained_variance_ratio_[0] * 100
    var2 = pca_obj.explained_variance_ratio_[1] * 100
    plt.xlabel(f"PC1 ({var1:.1f}% variance)")
    plt.ylabel(f"PC2 ({var2:.1f}% variance)")
    plt.title(title)
    plt.legend(title=cat_col, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    if output_name is not None:
        plt.savefig(out_dir / output_name, dpi=300)
    plt.show()


# --- BEFORE correction ---

plot_pca_continuous(
    pc_before_df,
    pca_before,
    color_col="SMRIN",
    title="PCA Before Confounder Correction (Colored by RIN)",
    output_name="PCA_before_RIN.png"
)

plot_pca_continuous(
    pc_before_df,
    pca_before,
    color_col="SMTSISCH",
    title="PCA Before Confounder Correction (Colored by Ischemic Time)",
    output_name="PCA_before_IschemicTime.png"
)

plot_pca_categorical(
    pc_before_df,
    pca_before,
    cat_col="SMGEBTCH",
    title="PCA Before Confounder Correction (Colored by Batch SMGEBTCH)",
    output_name="PCA_before_Batch.png"

)
# OPTIONAL: Death type plot (Before)
if "DTHHRDY" in pc_before_df.columns:
    plot_pca_categorical(
        pc_before_df,
        pca_before,
        cat_col="DTHHRDY",
        title="PCA Before Confounder Correction (Colored by Death Type - DTHHRDY)",
        output_name="PCA_before_DeathType.png"
    )


# --- AFTER correction ---

plot_pca_continuous(
    pc_after_df,
    pca_after,
    color_col="SMRIN",
    title="PCA After Confounder Correction (Colored by RIN)",
    output_name="PCA_after_RIN.png"
)

plot_pca_continuous(
    pc_after_df,
    pca_after,
    color_col="SMTSISCH",
    title="PCA After Confounder Correction (Colored by Ischemic Time)",
    output_name="PCA_after_IschemicTime.png"
)

plot_pca_categorical(
    pc_after_df,
    pca_after,
    cat_col="SMGEBTCH",
    title="PCA After Confounder Correction (Colored by Batch SMGEBTCH)",
    output_name="PCA_after_Batch.png"
)
# OPTIONAL: Death type plot (After)
if "DTHHRDY" in pc_after_df.columns:
    plot_pca_categorical(
        pc_after_df,
        pca_after,
        cat_col="DTHHRDY",
        title="PCA After Confounder Correction (Colored by Death Type - DTHHRDY)",
        output_name="PCA_after_DeathType.png"
    )




print("✅ PCA plots saved in:", out_dir)
