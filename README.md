# GTEx Artery–Coronary (SMTSISCH) — Final Project (Bio Data Science)

Anas Khoury 
Bshara Naoum 

**Question:** Does ischemic time (**SMTSISCH**) associate with gene-expression differences in **GTEx Artery–Coronary** tissue?  
**Comparison:** extreme groups  
- **LOW:** SMTSISCH ≤ Q1 → n=53  
- **HIGH:** SMTSISCH ≥ Q3 → n=53  

**Direction convention used throughout**
- **logFC = HIGH − LOW**
- **GSEA:** **NES > 0 = HIGH-enriched**, **NES < 0 = LOW-enriched**

---

## Key Rubric Constraint (Critical)
**Batch is handled ONLY via ComBat (Step 4) and is NOT included as a categorical term in the limma DE model.**  
- ComBat reference: https://rdrr.io/bioc/sva/man/ComBat.html  
- limma reference: https://bioconductor.org/packages/release/bioc/html/limma.html  

---

## Repository Structure

### `FinalProject.docx`
The final report (≤8 pages) submitted to the course website.

### `Code/` (final reproducible pipeline)
Scripts are numbered in the required run order:
1. `01_build_final_tables.py` — build aligned final inputs
2. `02_choose_batch_column.py` — select batch variable (chosen: **SMGEBTCH**)
3. `03_define_smtsisch_groups.py` — define LOW/HIGH extreme groups (Q1/Q3)
4. `04A_combat_sva.R` — ComBat batch correction (sva::ComBat)
5. `04B_pca_before_after.py` — PCA plots before/after ComBat
6. `05_limma_degs.R` — differential expression with limma (NO batch term)
7. `06_gsea_preranked_fgsea.R` — pre-ranked GSEA using fgsea + Hallmark sets

### `FirstStep/` (Step 1–3 outputs: aligned inputs + group definition)
- `expr_final.csv` — aligned expression matrix (genes × samples)
- `meta_final.csv` — aligned metadata
- `batch_diagnostics.csv` — diagnostics for batch candidates
- `chosen_batch_column.txt` — selected batch column (SMGEBTCH)
- `smtsisch_histogram.csv` — SMTSISCH distribution summary
- `meta_groups.csv` — extreme-group metadata (LOW/HIGH)
- `group_stats_table.csv` — thresholds + group sizes
- (optional) `expression_after_outlier_removal.csv`, `meta_after_outlier_removal.csv`

### `Pca/` (Step 4 outputs: ComBat + PCA validation)
- `expr_combat.csv` — ComBat-corrected expression matrix
- `combat_report.txt` — ComBat settings + summary
- `pca_before_by_batch.png` — PCA before correction (colored by SMGEBTCH)
- `pca_after_by_batch.png` — PCA after correction (colored by SMGEBTCH)

### `DifferentialExpression/` (Step 5 outputs: limma)
- `deg_table.csv` — full limma results
- `deg_top50.csv` — top genes list
- `deg_summary.txt` — summary counts (e.g., FDR<0.05)
- `volcano.png` — volcano plot (HIGH vs LOW)
- `top20_genes_bar.png` — Top 20 DE genes (by FDR)

### `GSEA/` (Step 6 outputs: pre-ranked Hallmark GSEA)
- `gsea_ranked_t_HighVsLow.rnk` — ranked list (gene symbol + limma t-stat)
- `gsea_hallmark_full.csv` — all Hallmark results
- `gsea_hallmark_top_up.csv` — top HIGH-enriched (NES>0)
- `gsea_hallmark_top_down.csv` — top LOW-enriched (NES<0)
- `gsea_hallmark_dotplot.png` — dotplot (NES vs -log10(FDR))
- `gsea_enrichment_*.png` — enrichment plots (examples: TNFA/NFKB, OXPHOS, etc.)

### `FirstPres/` (course labs 1–7 — not part of the final pipeline)
This folder contains materials from the **first course presentation / first 7 labs** (early work and course exercises).  
It is **not required** to reproduce the final results, but kept for completeness.

---

## Pipeline Summary (Steps 1–6)

**Step 1 — Build aligned inputs**  
Create aligned expression + metadata tables matched by sample IDs.

**Step 2 — Choose batch variable**  
Evaluate batch candidates and select the batch column used for correction (**SMGEBTCH**).

**Step 3 — Define extreme SMTSISCH groups**  
Select LOW/HIGH groups using quartiles (Q1/Q3) and build `meta_groups.csv`.

**Step 4 — Batch correction (ComBat)**  
Apply ComBat to remove batch effects while preserving biological/technical covariates.  
(ComBat: https://rdrr.io/bioc/sva/man/ComBat.html)

**Step 5 — Differential expression (limma)**  
Run limma on ComBat-corrected expression:  
- Model: `~ Group + AGE + SEX + DeathType + SMRIN` (**NO batch term**)  
- Contrast: **HIGH − LOW**  
- Multiple testing: BH-FDR  
(limma: https://bioconductor.org/packages/release/bioc/html/limma.html)

**Step 6 — Pre-ranked GSEA (Hallmark)**  
Rank genes by limma t-stat (HIGH − LOW) and run fgsea on MSigDB Hallmark sets.  
- fgsea: https://bioconductor.posit.co/packages/release/bioc/manuals/fgsea/man/fgsea.pdf  
- MSigDB Hallmark: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp

---

## How to Reproduce (Windows CMD)

From the repository root:

### 1) Python steps (1–3 + PCA plotting)
```bat
python Code\01_build_final_tables.py
python Code\02_choose_batch_column.py
python Code\03_define_smtsisch_groups.py
python Code\04B_pca_before_after.py

2) R steps (ComBat + limma + GSEA)
Rscript Code\04A_combat_sva.R
Rscript Code\05_limma_degs.R
Rscript Code\06_gsea_preranked_fgsea.R

Requirements (minimal)

Python: pandas, numpy, matplotlib, scikit-learn
R: sva, limma, fgsea, data.table (or tidyverse)

Outputs for Grading (quick list)

Batch correction: Pca/pca_before_by_batch.png, Pca/pca_after_by_batch.png

DE: DifferentialExpression/volcano.png, DifferentialExpression/top20_genes_bar.png

GSEA: GSEA/gsea_hallmark_dotplot.png, plus example enrichment plots in GSEA/

Full tables are included as CSV files in their respective folders.