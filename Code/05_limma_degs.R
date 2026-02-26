# 05_limma_degs.R
# Differential expression with limma on ComBat-corrected expression.
# Model: ~ Group + AGE + SEX + DeathType + SMRIN   (NO batch term)
# Outputs: deg_table.csv, volcano.png, top20_genes_bar.png, deg_summary.txt

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

IN_EXPR <- "expr_combat.csv"      # genes x samples (ComBat output)
IN_META <- "meta_groups.csv"      # only LOW/HIGH samples (Q1 vs Q3)
OUT_DEG <- "deg_table.csv"
OUT_TOP <- "deg_top50.csv"
OUT_SUM <- "deg_summary.txt"
OUT_VOL <- "volcano.png"
OUT_BAR <- "top20_genes_bar.png"

# -----------------------
# Load expression
# -----------------------
expr <- fread(IN_EXPR, data.table=FALSE)
gene_id <- expr[[1]]
expr <- as.matrix(expr[, -1, drop=FALSE])
rownames(expr) <- gene_id

# -----------------------
# Load metadata (groups)
# -----------------------
meta <- fread(IN_META, data.table=FALSE)

# Ensure required cols exist
req <- c("SAMPID","Group","AGE","SEX","DeathType","SMRIN")
miss <- setdiff(req, colnames(meta))
if (length(miss) > 0) stop(paste("Missing columns in meta_groups.csv:", paste(miss, collapse=", ")))

# Align samples
stopifnot(all(meta$SAMPID %in% colnames(expr)))
expr <- expr[, meta$SAMPID, drop=FALSE]
stopifnot(identical(colnames(expr), meta$SAMPID))

# Make factors
meta$Group <- factor(meta$Group, levels=c("LOW","HIGH"))  # HIGH vs LOW
meta$AGE <- factor(meta$AGE)                              # GTEx age bins
meta$SEX <- factor(meta$SEX)                              # 1/2
meta$DeathType <- factor(meta$DeathType)                  # 0-4
meta$SMRIN <- as.numeric(meta$SMRIN)

# -----------------------
# Design matrix
# -----------------------
design <- model.matrix(~ Group + AGE + SEX + DeathType + SMRIN, data=meta)

# Fit limma
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# We want the Group effect: HIGH vs LOW
coef_name <- "GroupHIGH"
if (!(coef_name %in% colnames(fit$coefficients))) {
  stop(paste("Could not find coefficient", coef_name, "in design. Coefs are:", paste(colnames(fit$coefficients), collapse=", ")))
}

deg <- topTable(fit, coef=coef_name, number=Inf, sort.by="P", adjust.method="BH")
deg$Gene <- rownames(deg)

# Save full table
fwrite(deg, OUT_DEG)

# Save top 50
top50 <- head(deg, 50)
fwrite(top50, OUT_TOP)

# Summary text
n_fdr05 <- sum(deg$adj.P.Val < 0.05, na.rm=TRUE)
n_fdr10 <- sum(deg$adj.P.Val < 0.10, na.rm=TRUE)

summ <- c(
  "DEGs (limma) summary",
  paste0("Samples: ", ncol(expr), "  (LOW=", sum(meta$Group=="LOW"), ", HIGH=", sum(meta$Group=="HIGH"), ")"),
  paste0("Genes tested: ", nrow(expr)),
  paste0("Model: ~ Group + AGE + SEX + DeathType + SMRIN  (batch handled via ComBat only)"),
  paste0("Contrast tested: HIGH vs LOW (coef=", coef_name, ")"),
  paste0("Significant DEGs: FDR<0.05 => ", n_fdr05),
  paste0("Significant DEGs: FDR<0.10 => ", n_fdr10)
)
writeLines(summ, OUT_SUM)

# -----------------------
# Volcano plot
# -----------------------
deg$negLog10P <- -log10(deg$P.Value)
deg$Signif <- deg$adj.P.Val < 0.05

p_vol <- ggplot(deg, aes(x=logFC, y=negLog10P)) +
  geom_point(aes(alpha=Signif), size=1.4) +
  scale_alpha_manual(values=c(`FALSE`=0.35, `TRUE`=0.9), guide="none") +
  labs(
    title="Volcano: HIGH vs LOW SMTSISCH (limma)",
    x="logFC (HIGH - LOW)",
    y="-log10(p-value)"
  ) +
  theme_minimal()

ggsave(OUT_VOL, p_vol, width=7, height=5, dpi=200)

# -----------------------
# Top 20 bar plot (by FDR)
# -----------------------
top20 <- head(deg[order(deg$adj.P.Val), ], 20)
top20$Gene <- factor(top20$Gene, levels=rev(top20$Gene))

p_bar <- ggplot(top20, aes(x=Gene, y=logFC)) +
  geom_col() +
  coord_flip() +
  labs(
    title="Top 20 DE genes (by FDR)",
    x="Gene",
    y="logFC (HIGH - LOW)"
  ) +
  theme_minimal()

ggsave(OUT_BAR, p_bar, width=7, height=6, dpi=200)

cat("Wrote:", OUT_DEG, "\n")
cat("Wrote:", OUT_TOP, "\n")
cat("Wrote:", OUT_SUM, "\n")
cat("Wrote:", OUT_VOL, "\n")
cat("Wrote:", OUT_BAR, "\n")