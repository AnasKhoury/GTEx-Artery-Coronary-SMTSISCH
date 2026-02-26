# 06_gsea_preranked_fgsea.R
# Pre-ranked GSEA (Hallmark) using limma t-statistic (HIGH vs LOW)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

IN_DEG   <- "deg_table.csv"
OUT_RNK  <- "gsea_ranked_t_HighVsLow.rnk"
OUT_FULL <- "gsea_hallmark_full.csv"
OUT_UP   <- "gsea_hallmark_top_up.csv"
OUT_DN   <- "gsea_hallmark_top_down.csv"
OUT_DOT  <- "gsea_hallmark_dotplot.png"

# -----------------------
# Load limma table
# -----------------------
deg <- fread(IN_DEG)

stopifnot(all(c("Gene","t","adj.P.Val","logFC") %in% colnames(deg)))

# Strip Ensembl version suffix (ENSG... .#)
deg$ensembl <- sub("\\..*$", "", deg$Gene)

# -----------------------
# Map Ensembl -> HGNC symbol
# Option A (recommended): biomaRt (needs internet)
# -----------------------
use_biomart <- FALSE

if (use_biomart) {
  suppressPackageStartupMessages(library(biomaRt))
  mart <- biomaRt::useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl")

  map <- biomaRt::getBM(
    attributes = c("ensembl_gene_id","hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = unique(deg$ensembl),
    mart       = mart
  )

  deg <- deg %>%
    left_join(map, by=c("ensembl"="ensembl_gene_id")) %>%
    rename(symbol = hgnc_symbol)
} else {
  # Option B: org.Hs.eg.db (works offline if installed)
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })
  deg$symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = unique(deg$ensembl),
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )[deg$ensembl]
}

# Clean symbols
deg <- deg %>%
  mutate(symbol = ifelse(is.na(symbol) | symbol=="", NA, symbol)) %>%
  filter(!is.na(symbol))

# If multiple Ensembl map to same symbol, keep the strongest signal (max |t|)
deg_collapsed <- deg %>%
  dplyr::group_by(symbol) %>%
  dplyr::slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# -----------------------
# Build ranked list (named numeric vector)
# -----------------------
ranks <- deg_collapsed$t
names(ranks) <- deg_collapsed$symbol
ranks <- sort(ranks, decreasing = TRUE)

# Save .rnk (2 columns)
fwrite(
  data.table(gene=names(ranks), score=as.numeric(ranks)),
  OUT_RNK, sep="\t", col.names=FALSE
)

# -----------------------
# Load Hallmark gene sets from MSigDB via msigdbr
# -----------------------
h <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)

pathways <- split(h$gene_symbol, h$gs_name)

# -----------------------
# Run pre-ranked GSEA (fgsea multilevel is recommended)
# -----------------------
set.seed(1)
fg <- fgseaMultilevel(
  pathways = pathways,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500
)

fg <- fg %>%
  as.data.frame() %>%
  arrange(padj) %>%
  mutate(direction = ifelse(NES > 0, "Enriched in HIGH", "Enriched in LOW"))

fwrite(fg, OUT_FULL)

# Top up/down tables
top_up <- fg %>% filter(NES > 0) %>% arrange(padj) %>% head(15)
top_dn <- fg %>% filter(NES < 0) %>% arrange(padj) %>% head(15)
fwrite(top_up, OUT_UP)
fwrite(top_dn, OUT_DN)

# -----------------------
# Dotplot (NES vs -log10(FDR))
# -----------------------
fg2 <- fg %>% mutate(minusLog10FDR = -log10(padj + 1e-300))

p <- ggplot(fg2, aes(x = NES, y = minusLog10FDR)) +
  geom_point(alpha=0.75) +
  labs(
    title = "Hallmark GSEA (pre-ranked by limma t-stat)",
    x = "NES (positive = HIGH-enriched)",
    y = "-log10(FDR)"
  ) +
  theme_minimal()

ggsave(OUT_DOT, p, width=7, height=5, dpi=200)

# -----------------------
# Enrichment plots for top pathways (2 up + 2 down)
# -----------------------
plot_one <- function(gs_name, out_png) {
  g <- plotEnrichment(pathways[[gs_name]], ranks) +
    labs(title = gs_name, subtitle = "Pre-ranked GSEA (fgsea)")
  ggsave(out_png, g, width=7, height=5, dpi=200)
}

top_up_names <- top_up$pathway[1:min(2, nrow(top_up))]
top_dn_names <- top_dn$pathway[1:min(2, nrow(top_dn))]

for (nm in c(top_up_names, top_dn_names)) {
  out <- paste0("gsea_enrichment_", nm, ".png")
  plot_one(nm, out)
}

cat("Wrote:", OUT_RNK, "\n")
cat("Wrote:", OUT_FULL, "\n")
cat("Wrote:", OUT_UP, "\n")
cat("Wrote:", OUT_DN, "\n")
cat("Wrote:", OUT_DOT, "\n")
cat("Wrote: enrichment plots for selected pathways\n")