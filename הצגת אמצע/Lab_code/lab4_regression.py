import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import linregress

# ----------- paths (EDIT ONLY IF NEEDED) -----------
expr_path = r"C:\Users\97250\Downloads\Lab4_GTex\expr_artery_coronary_QN_RIN_gt_5_8.csv"
meta_path = r"C:\Users\97250\Downloads\Lab4_GTex\meta_artery_coronary_RIN_gt_5_8.csv"

# Subject-level phenotypes (GTEx v8) – contains AGE / SEX(or GENDER) / DTHHRDY
pheno_path = r"C:\Users\97250\Downloads\GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"

out_dir = Path(r"C:\Users\97250\Downloads\Lab4_GTex\Lab4_plots")
out_dir.mkdir(parents=True, exist_ok=True)

# ----------- load data -----------
expr = pd.read_csv(expr_path, index_col=0)     # genes × samples
meta = pd.read_csv(meta_path)                  # samples in rows

print("Expression shape:", expr.shape)
print("Meta shape (before merge):", meta.shape)

# make sure SAMPID exists
if "SAMPID" not in meta.columns:
    raise ValueError("Column 'SAMPID' not found in meta-data file.")

# ==========================================================
# MERGE BLOCK: bring AGE / SEX / DeathType into meta
# ==========================================================
pheno = pd.read_csv(pheno_path, sep="\t", low_memory=False)
pheno.columns = pheno.columns.str.strip()

# Find the sex column name (some GTEx versions use GENDER)
sex_col = None
for cand in ["SEX", "GENDER", "Gender", "sex", "gender"]:
    if cand in pheno.columns:
        sex_col = cand
        break

if sex_col is None:
    raise ValueError(
        "Could not find a sex column in SubjectPhenotypesDS.txt. "
        f"Available columns are: {list(pheno.columns)}"
    )

# Require core columns
required = {"SUBJID", "AGE", "DTHHRDY"}
missing = required - set(pheno.columns)
if missing:
    raise ValueError(f"SubjectPhenotypesDS.txt missing columns: {missing}")

# Build SUBJID from SAMPID: GTEX-1117F-0006-SM-... -> GTEX-1117F
parts = meta["SAMPID"].astype(str).str.split("-", expand=True)
if parts.shape[1] < 2:
    raise ValueError("Unexpected SAMPID format; cannot build SUBJID.")
meta["SUBJID"] = parts[0] + "-" + parts[1]

# Merge subject-level traits into the sample-level meta
meta = meta.merge(
    pheno[["SUBJID", sex_col, "AGE", "DTHHRDY"]].rename(columns={sex_col: "SEX"}),
    on="SUBJID",
    how="left"
)

# Rename DTHHRDY (Hardy scale) to DeathType, as your code expects
meta = meta.rename(columns={"DTHHRDY": "DeathType"})

print("Meta shape (after merge):", meta.shape)
print("Merged traits availability:",
      {k: (k in meta.columns) for k in ["AGE", "SEX", "DeathType"]})

# ==========================================================
# align samples between expression and meta-data
# ==========================================================
sample_ids = [s for s in expr.columns if s in set(meta["SAMPID"])]
expr = expr[sample_ids]
meta = meta[meta["SAMPID"].isin(sample_ids)].set_index("SAMPID").loc[sample_ids]

print("Aligned expression shape:", expr.shape)
print("Aligned meta shape:", meta.shape)

# ========================================
# 1) choose 5 genes with strong corr to SMRIN
# ========================================
trait_for_selection = "SMRIN"
if trait_for_selection not in meta.columns:
    raise ValueError("SMRIN not found in meta-data; needed to select genes.")

x_rin = pd.to_numeric(meta[trait_for_selection], errors="coerce")
valid = x_rin.notna()

corrs = []
for gene in expr.index:
    y = expr.loc[gene, valid]
    if y.var() == 0:
        continue
    r = np.corrcoef(x_rin[valid], y)[0, 1]
    if np.isnan(r):
        continue
    corrs.append((gene, r))

# --- choose TOP 5 POSITIVE correlations with SMRIN ---
positive_corrs = [(gene, r) for gene, r in corrs if r > 0]
positive_corrs = sorted(positive_corrs, key=lambda t: t[1], reverse=True)

top_genes = [g for g, r in positive_corrs[:5]]

print("\nTop 5 genes by POSITIVE correlation with SMRIN:")
for g, r in positive_corrs[:5]:
    print(f"{g}: r = {r:.3f}")

# Safety: if there are fewer than 5 positive genes (rare), fall back gracefully
if len(top_genes) < 5:
    print(f"\n[Warning] Only found {len(top_genes)} genes with positive correlation. "
          f"Filling the rest with strongest |corr| genes.")
    corrs_abs = sorted(corrs, key=lambda t: abs(t[1]), reverse=True)
    for g, r in corrs_abs:
        if g not in top_genes:
            top_genes.append(g)
        if len(top_genes) == 5:
            break

# ========================================
# 2) regression plots for each trait & gene
# ========================================
traits = ["SMGEBTCH", "SMNABTCH", "SMTSISCH", "SMRIN", "AGE", "SEX", "DeathType"]

# Traits that should NOT have a regression line (categorical / binned)
no_reg_line_traits = {"AGE", "SEX", "DeathType"}

for trait in traits:
    if trait not in meta.columns:
        print(f"\n[Warning] Trait '{trait}' not found in meta-data. Skipping.")
        continue

    print(f"\nProcessing trait: {trait}")

    t_series = meta[trait]

    # Build x values correctly depending on the trait type
    if trait == "AGE":
        # GTEx AGE is usually binned like "60-69" => take lower bound
        age_num = t_series.astype(str).str.extract(r"(\d+)")[0]
        x_all = pd.to_numeric(age_num, errors="coerce")
        x_label = "AGE (binned → lower bound)"
    elif pd.api.types.is_numeric_dtype(t_series):
        x_all = pd.to_numeric(t_series, errors="coerce")
        x_label = trait
    else:
        # encode categories as integers for correlation/regression
        codes, uniques = pd.factorize(t_series)
        x_all = pd.Series(codes, index=t_series.index, dtype=float)
        x_label = f"{trait} (encoded categories)"

    for gene in top_genes:
        y_all = expr.loc[gene]

        mask = x_all.notna() & y_all.notna()
        x = x_all[mask]
        y = y_all[mask]

        # need enough points
        if len(x) < 3:
            continue

        # skip if no variance (regression will fail / meaningless)
        if np.var(x) == 0 or np.var(y) == 0:
            continue

        # (Keep r and p in the title like before)
        slope, intercept, r, p, _ = linregress(x, y)

        plt.figure(figsize=(6, 4))

        # Optional: jitter x for categorical/binned traits so points don't sit exactly on top of each other
        if trait in no_reg_line_traits:
            rng = np.random.default_rng(0)
            x_plot = x + rng.normal(0, 0.03, size=len(x))
        else:
            x_plot = x

        plt.scatter(x_plot, y, alpha=0.6)

        # Draw regression line ONLY for continuous traits
        if trait not in no_reg_line_traits:
            line_x = np.linspace(x.min(), x.max(), 100)
            line_y = slope * line_x + intercept
            plt.plot(line_x, line_y, color="red", linewidth=2)

        plt.xlabel(x_label)
        plt.ylabel(f"Expression of {gene}")
        plt.title(f"{gene} vs {trait}\n"
                  f"r = {r:.3f}, p = {p:.2e}")

        plt.tight_layout()
        out_file = out_dir / f"{trait}_{gene}.png"
        plt.savefig(out_file, dpi=150)
        plt.close()

        print(f"  Saved plot: {out_file}")
