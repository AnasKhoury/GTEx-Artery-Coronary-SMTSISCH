# GTEx Artery–Coronary Final Project (Data Science Bio)
**Topic:** Effect of ischemic time (**SMTSISCH**) on gene expression in **GTEx Artery–Coronary**  
**Deliverables:** Report (≤8 pages) + GitHub link to this repository

---

## 1) Research Question
**Does ischemic time (SMTSISCH) associate with transcriptomic differences in Artery–Coronary tissue?**  
We compare extreme groups:
- **LOW:** SMTSISCH ≤ Q1 (139) → n=53  
- **HIGH:** SMTSISCH ≥ Q3 (567) → n=53  

Direction convention used throughout the project:
- **logFC = HIGH − LOW**
- In GSEA: **NES > 0 = HIGH-enriched**, **NES < 0 = LOW-enriched**

---

## 2) Key Rubric Constraint (Important)
**Batch is NOT included as a categorical variable in the differential expression (DE) model.**  
Batch is handled **only via ComBat** (Step 4).  
The limma model is: `~ Group + AGE + SEX + DeathType + SMRIN` (NO batch term).

---

## 3) Folder Structure (What is Where)

### `/Report/`
- `FinalProject.docx` (final write-up, ≤8 pages)
- (optional) `FinalProject.pdf`

### `/Code/` (reproducible pipeline scripts)
1. `01_build_final_tables.py`  
2. `02_choose_batch_column.py`  
3. `03_define_smtsisch_groups.py`  
4. `04A_combat_sva.R`  
5. `04B_pca_before_after.py`  
6. `05_limma_degs.R`  
7. `06_gsea_preranked_fgsea.R`

### `/first step/` (cleaned + aligned inputs, grouping)
- `expr_final.csv` — aligned expression matrix (genes × samples)
- `meta_final.csv` — aligned metadata
- `batch_diagnostics.csv` — batch candidate diagnostics
- `chosen_batch_column.txt` — selected batch variable (SMGEBTCH)
- `smtsisch_histogram.csv` — SMTSISCH distribution summary
- `meta_groups.csv` — extreme-group metadata (LOW/HIGH)
- `group_stats_table.csv` — thresholds + group sizes summary
- `Archive_or_unused/` — optional: files produced earlier but not used in final analysis (e.g., outlier-removal versions)

### `/pca/` (ComBat outputs + PCA diagnostics)
- `expr_combat.csv` — ComBat-corrected expression matrix
- `combat_report.txt` — ComBat settings + summary (#samples, #genes, covariates)
- `pca_before_by_batch.png` — PCA before correction (colored by SMGEBTCH)
- `pca_after_by_batch.png` — PCA after correction (colored by SMGEBTCH)

### `/Differential Expression/` (Step 5 outputs)
- `deg_table.csv` — full limma results table
- `deg_top50.csv` — top genes list
- `deg_summary.txt` — summary counts (e.g., FDR<0.05, FDR<0.10)
- `volcano.png` — volcano plot (HIGH vs LOW)
- `top20_genes_bar.png` — Top 20 DE genes (by FDR)

### `/GSEA/` (Step 6 outputs)
- `gsea_ranked_t_HighVsLow.rnk` — ranked list (gene symbol + limma t-stat)
- `gsea_hallmark_full.csv` — all Hallmark results
- `gsea_hallmark_top_up.csv` — top HIGH-enriched pathways (NES>0)
- `gsea_hallmark_top_down.csv` — top LOW-enriched pathways (NES<0)
- `gsea_hallmark_dotplot.png` — dotplot (NES vs -log10(FDR))
- `gsea_enrichment_*.png` — enrichment plots (examples: TNFA/NFKB, OXPHOS, PEROXISOME, MYC_TARGETS)

---

## 4) Pipeline Summary (Steps 1–6)

### Step 1 — Build aligned final inputs
Goal: ensure expression and metadata are aligned by sample IDs.
- Output: `expr_final.csv`, `meta_final.csv`

### Step 2 — Choose batch variable
Goal: pick the best batch variable for correction diagnostics.
- Chosen: **SMGEBTCH**
- Output: `batch_diagnostics.csv`, `chosen_batch_column.txt`

### Step 3 — Define SMTSISCH extreme groups (Q1 vs Q3)
Goal: define LOW and HIGH groups using quartiles (Q1/Q3).
- Output: `smtsisch_thresholds.txt` (or recorded in group stats), `meta_groups.csv`, `group_stats_table.csv`

### Step 4 — ComBat batch correction (sva::ComBat)
Goal: remove batch effects using ComBat while preserving covariates:
- batch: SMGEBTCH
- covariates preserved: SMTSISCH + AGE + SEX + DeathType + SMRIN
- Output: `expr_combat.csv`, `combat_report.txt`, PCA plots before/after

### Step 5 — Differential expression (limma)
Goal: identify DEGs between HIGH vs LOW after ComBat.
- Model: `~ Group + AGE + SEX + DeathType + SMRIN` (NO batch term)
- Contrast: HIGH − LOW
- Multiple testing: BH-FDR
- Output: `deg_table.csv`, plots, summaries

### Step 6 — GSEA (pre-ranked)
Goal: pathway-level interpretation using ranked gene list.
- Ranking: limma **t-statistic** (HIGH − LOW)
- Gene sets: MSigDB **Hallmark (H)**
- Output: `gsea_hallmark_full.csv`, top tables, dotplot, enrichment plots

---

## 5) How to Reproduce (Windows CMD)

### Run scripts in order
Open CMD in the `Project` folder and run:

#### Python steps (1–3, 4B)
```bat
python Code\01_build_final_tables.py
python Code\02_choose_batch_column.py
python Code\03_define_smtsisch_groups.py
python Code\04B_pca_before_after.py