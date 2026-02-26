# 04A_combat_sva.R
# ComBat batch correction using Bioconductor sva::ComBat
# Input:  expr_final.csv (genes x samples), meta_final.csv
# Output: expr_combat.csv (genes x samples) + combat_report.txt

suppressPackageStartupMessages({
  library(data.table)
  library(sva)
})

IN_EXPR <- "expr_final.csv"
IN_META <- "meta_final.csv"
BATCH_COL <- "SMGEBTCH"

OUT_EXPR <- "expr_combat.csv"
OUT_REP  <- "combat_report.txt"

# Read expression (genes x samples)
expr <- fread(IN_EXPR, data.table=FALSE)
gene_id <- expr[[1]]
expr <- as.matrix(expr[, -1, drop=FALSE])
rownames(expr) <- gene_id

# Read metadata
meta <- fread(IN_META, data.table=FALSE)

# Align meta to expr column order
stopifnot(all(colnames(expr) %in% meta$SAMPID))
meta <- meta[match(colnames(expr), meta$SAMPID), ]
stopifnot(identical(colnames(expr), meta$SAMPID))

# Remove samples with missing batch
batch <- meta[[BATCH_COL]]
keep <- !is.na(batch)
expr <- expr[, keep, drop=FALSE]
meta <- meta[keep, , drop=FALSE]
batch <- as.factor(meta[[BATCH_COL]])

# Model matrix (covariates preserved)
# AGE is categorical like "50-59" -> factor
# SEX is 1/2 -> factor
# DeathType is 0-4 -> factor
mod <- model.matrix(~ SMTSISCH + AGE + factor(SEX) + factor(DeathType) + SMRIN, data=meta)

# Run ComBat
combat_expr <- ComBat(dat=expr, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

# Write corrected matrix
out_df <- data.frame(Name=rownames(combat_expr), combat_expr, check.names=FALSE)
fwrite(out_df, OUT_EXPR)

# Report
rep <- c(
  "ComBat (sva::ComBat) batch correction report",
  paste0("Batch column: ", BATCH_COL),
  paste0("Samples used: ", ncol(combat_expr)),
  paste0("Genes used: ", nrow(combat_expr)),
  "",
  "Batch size summary:",
  paste(capture.output(print(summary(table(batch)))), collapse="\n")
)
writeLines(rep, OUT_REP)

cat("Wrote:", OUT_EXPR, "\n")
cat("Wrote:", OUT_REP, "\n")