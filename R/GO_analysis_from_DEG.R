#!/usr/bin/env Rscript

# =========================
# Usage:
#   Rscript GO_analysis_from_DEG.R deg.tsv outprefix up BP
#   Rscript GO_analysis_from_DEG.R deg.tsv outprefix down ALL
#   Rscript GO_analysis_from_DEG.R go_enrich.R deg.tsv outprefix all BP
#   Rscript GO_analysis_from_DEG.R outprefix up BP 0.5 0.1
#
# Args:
#   1) deg file (tsv/csv, must contain: transcript, log2FoldChange, padj)
#   2) outpdf_prefix
#   3) direction: up | down | all
#   4) ont: BP | MF | CC | ALL   (optional, default BP)
#   5) lfc_cutoff (optional, default 1)
#   6) padj_cutoff (optional, default 0.05)
# =========================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Need at least: deg, outpdf_prefix, direction. Example: Rscript GO_analysis_from_DEG.R deg.tsv out up BP")
}

deg <- args[1]
outpdf_prefix <- args[2]
direction <- tolower(args[3])

ont <- ifelse(length(args) >= 4, toupper(args[4]), "BP")
lfc_cutoff <- ifelse(length(args) >= 5, as.numeric(args[5]), 1)
padj_cutoff <- ifelse(length(args) >= 6, as.numeric(args[6]), 0.05)

if (!direction %in% c("up", "down", "all")) stop("direction 必须是 up / down / all")
if (!ont %in% c("BP", "MF", "CC", "ALL")) stop("ont 必须是 BP / MF / CC / ALL")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(enrichplot)
})

# ---- read DEG ----
res <- read.table(deg, header = TRUE, stringsAsFactors = FALSE, sep = "", quote = "", comment.char = "")

needed <- c("transcript", "log2FoldChange", "padj")
miss <- setdiff(needed, colnames(res))
if (length(miss) > 0) stop(paste0("缺少列：", paste(miss, collapse = ", ")))

# ---- pick genes by direction ----
if (direction == "up") {
  genes <- res$transcript[res$log2FoldChange > lfc_cutoff & res$padj < padj_cutoff]
  title_prefix <- "Up-Genes"
} else if (direction == "down") {
  genes <- res$transcript[res$log2FoldChange < lfc_cutoff & res$padj < padj_cutoff]
  title_prefix <- "Down-Genes"
} else {
  genes <- res$transcript[abs(res$log2FoldChange) > lfc_cutoff & res$padj < padj_cutoff]
  title_prefix <- "DE-Genes"
}

genes <- unique(as.character(genes))
genes <- trimws(genes)
genes <- genes[!is.na(genes) & genes != ""]

if (length(genes) == 0) stop("No Genes left. Check cutoff or column names.")

# ---- SYMBOL -> ENTREZ ----
gene.df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

if (is.null(gene.df) || nrow(gene.df) == 0) stop("bitr failed: no matching ENTREZID (check SYMBOL / species)")

entrez <- unique(gene.df$ENTREZID)

# ---- enrichGO runner ----
run_enrich_go <- function(entrez_ids, ont_one) {
  enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = ont_one,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
}

# ---- plot & save (dot + bar) ----
plot_and_save <- function(ego_obj, out_prefix, plot_title, show_n = 20) {
  df <- as.data.frame(ego_obj)
  if (is.null(df) || nrow(df) == 0) {
    message("No significant terms for: ", out_prefix, " (skip plots)")
    return(invisible(NULL))
  }

  n_take <- min(show_n, nrow(df))
  p1 <- dotplot(ego_obj, showCategory = n_take) + 
    ggtitle(plot_title) + 
    scale_y_discrete(labels = rev(df$Description[1:n_take]))
  ggsave(paste0(out_prefix, ".dot.pdf"), p1, width = 10, height = 6)

  p2 <- barplot(ego_obj, showCategory = n_take) + 
    ggtitle(plot_title) + 
    scale_y_discrete(labels = rev(df$Description[1:n_take]))
  ggsave(paste0(out_prefix, ".bar.pdf"), p2, width = 10, height = 6)

  # 同时输出表格
  write.table(df,
              file = paste0(out_prefix, ".enrichGO.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(NULL)
}

# ---- main logic ----
if (ont == "ALL") {
  onts <- c("BP", "MF", "CC")
  for (o in onts) {
    ego_o <- run_enrich_go(entrez, o)
    out_prefix <- paste0(outpdf_prefix, ".", direction, ".ALL.", o)
    plot_title <- paste0("GO ", o, " Enrichment of ", title_prefix)
    plot_and_save(ego_o, out_prefix, plot_title, show_n = 20)
  }
} else {
  # 单一 ontology
  ego <- run_enrich_go(entrez, ont)
  out_prefix <- paste0(outpdf_prefix, ".", direction, ".", ont)
  plot_title <- paste0("GO ", ont, " Enrichment of ", title_prefix)
  plot_and_save(ego, out_prefix, plot_title, show_n = 20)
}

message("Done.")
message("Genes used (SYMBOL): ", length(genes))
message("Mapped ENTREZ used: ", length(entrez))