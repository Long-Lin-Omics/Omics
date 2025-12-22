#!/usr/bin/env Rscript
Sys.setlocale("LC_COLLATE", "C")

# Load libraries
suppressPackageStartupMessages({
  library(argparse)
  library(DESeq2)
  library(tximport)
  library(biomaRt)
  library(pheatmap)
  library(ggplot2)
  library(RColorBrewer)
})

# Create parser
parser <- ArgumentParser(description="Run DESeq2 RNA-seq differential expression analysis")

parser$add_argument("--samples", required=TRUE, help="Samples file (tab-delimited with 'sample', 'path' and 'condition' columns)")
parser$add_argument("--tx2gene", required=TRUE, help="Transcript-to-gene mapping file (2-column TSV)")
parser$add_argument("--genesHeatmap", required=FALSE, help="first column will be used as gene list to plot heatmap.")
parser$add_argument("--count", type="integer", required=FALSE, help="Genes whose reads number among samples above the value will be kept for analysis")
parser$add_argument("--group_size", type="integer", required=FALSE, help="Genes with at least this size of individuals, who have at least --count reads or 10, will be kept for analysis")
parser$add_argument("--ercc_norm", action="store_true", help="Use ERCC spike-ins for normalization")
parser$add_argument("--rev_ref", action="store_true", help="reverse the reference level")
parser$add_argument("--output_prefix", required=TRUE, help="Output prefix for files and plots")

args <- parser$parse_args()

# Load sample table and tx2gene
print('reading sample.list')
samples <- read.table(args$samples, header=TRUE, sep="\t")
samples$condition <- factor(samples$condition)
if (args$rev_ref){
  samples$condition <- factor(samples$condition, levels = rev(levels(samples$condition)))
}
tx2gene <- read.table(args$tx2gene, header=FALSE)
colnames(tx2gene) <- c("TXNAME", "GENEID")

# Load counts via tximport
print('loading counts')
files <- as.character(samples$path)
names(files) <- as.character(samples$sample)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

# DESeqDataSet
dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~condition)

# Optional Filtering
filtered = FALSE
if (!is.null(args$group_size) & !is.null(args$count)){
    keep <- rowSums(counts(dds) >= args$count) >= args$group_size
    print(paste0('Filter genes with --count: ', args$count, ' and --group_size: ', args$group_size))
    filtered = TRUE
}else if (!is.null(args$group_size)){
    keep <- rowSums(counts(dds) >= 10) >= args$group_size
    print(paste0('Filter genes with --count: 10 and --group_size: ', args$group_size))
    filtered = TRUE
}else if(!is.null(args$count)){
    keep <- rowSums(counts(dds)) > args$count
    print(paste0('Filter genes with --count: 10'))
    filtered = TRUE
}else{
    print('No Filter on Genes') 
}


# Optional ERCC normalization
if (args$ercc_norm) {
  ercc_genes <- grep("^ERCC", rownames(dds))
  print(colSums(counts(dds[ercc_genes,])))
  if (filtered){
    raw_gene_n = length(names(dds))
    # always keep all the ercc genes even with filters
    keep <- keep | grepl("^ERCC", rownames(dds))
    dds <- dds[keep,]
    print(paste(sum(keep), 'out of ', raw_gene_n, ' genes were kept!' ))
  }
  print(paste0(length(ercc_genes), ' ERCC genes are used to perform normalization.'))
  dds <- estimateSizeFactors(dds, controlGenes=ercc_genes)
}else{
  ercc_genes <- grep("^ERCC", rownames(dds))
  print('Reads mapped to ERCC genes (not for normalization) :')
  print(colSums(counts(dds[ercc_genes,])))
  if (filtered){
    raw_gene_n = length(names(dds))
    keep <- keep & !grepl("^ERCC", rownames(dds))
    dds <- dds[keep,]
    print(paste(sum(keep), 'out of ', raw_gene_n, ' genes were kept!' ))
  }
}


# Run DESeq2
dds <- DESeq(dds)
res <- results(dds) #  baseMean log2FoldChange     lfcSE        stat      pvalue        padj
resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm") # change log2FoldChange, lfcSE and delete stat

#output
summary_file = paste0(args$output_prefix,'.summary.txt')
deg_tsv = paste0(args$output_prefix,'.deg.txt')
QA_plot = paste0(args$output_prefix,'.QA.png')
heatmap_plot = paste0(args$output_prefix,'.heatmap.png')
pca_plot = paste0(args$output_prefix,'.pca.png')
valcano_plot = paste0(args$output_prefix,'.valcano.png')

# Save DEG results and summary
final_res = cbind(resLFC,res[,c('log2FoldChange','lfcSE','stat')])
colnames(final_res)[6:7] = paste0(colnames(final_res)[6:7],'_no_LFC')

# Get normalized counts
norm_counts <- counts(dds, normalized=TRUE)
if (! is.null(args$genesHeatmap)){
  deg_genes <- read.table(args$genesHeatmap,h=F,stringsAsFactors=F)$V1
  missing = sum(!deg_genes %in% rownames(norm_counts))
  print(paste(missing, "genes not in the matrix"))
  mat <- norm_counts[deg_genes, ]
  mat_z <- t(scale(t(mat)))
  write.table(mat_z,file='nolog.txt')
  geneHeatmp <- paste0(args$output_prefix,'.geenHeatmap.png')
  png(geneHeatmp,width = 1500, height = 1500,res=300)
  pheatmap(mat_z,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize = 10,
         main = paste0("Heatmap of DEGs(n=",length(deg_genes)-missing,") (row z-score)"))
  dev.off()
  log_counts <- log2(norm_counts + 1)
  mat <- log_counts[deg_genes, ]
  mat_z <- t(scale(t(mat)))
  geneHeatmp <- paste0(args$output_prefix,'.geenHeatmap.log.png')
  png(geneHeatmp,width = 1500, height = 1500,res=300)
  pheatmap(mat_z,
         show_rownames = FALSE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize = 10,
         main = paste0("Heatmap of DEGs(n=",length(deg_genes)-missing,") (row log z-score)"))
  dev.off()
  write.table(mat_z,file='log.txt')
}else {

  # Get sample groups
  groupA <- rownames(colData(dds))[dds$condition == levels(dds$condition)[1]]
  groupB <- rownames(colData(dds))[dds$condition == levels(dds$condition)[2]]

  # Calculate mean normalized counts for each group
  meanA <- rowMeans(norm_counts[, groupA, drop=FALSE])
  meanB <- rowMeans(norm_counts[, groupB, drop=FALSE])

  # Add to final results
  final_res <- cbind(final_res, norm_counts, meanA, meanB)

  colnames(final_res)[ (ncol(final_res)-1):ncol(final_res) ] <- paste0("mean_", levels(samples$condition))

  final_res <- cbind(transcript = rownames(final_res), final_res)

  final_res <- final_res[order(final_res$pvalue),]

  write.table(final_res, file=deg_tsv, sep="\t", quote=FALSE, row.names = FALSE)
  
  sink(summary_file)
  print(resultsNames(dds)[2])
  summary(res)
  sink()

  # Generate QA plots: Dispersion, MA before and after LFC 
  png(QA_plot,width = 3000, height = 3000, res=300)
  # Set up a 2x2 plotting area
  par(mfrow = c(2, 2))  # 2 rows, 2 columns
  plotDispEsts(dds)
  plotMA(res, ylim=c(-5,5), main="MA plot (before shrinkage)")
  abline(h=c(-1,1), col="dodgerblue", lwd=2)
  plotMA(resLFC, ylim=c(-5,5), main="MA plot (after shrinkage)")
  abline(h=c(-1,1), col="dodgerblue", lwd=2)
  dev.off()

  # pca and heatmap
  png(heatmap_plot,width = 1500, height = 1500,res=300)
  vsd <- vst(dds)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- vsd$samples
  colnames(sampleDistMatrix) <- vsd$samples
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  # png(pca_plot,width = 1500, height = 1500,res=300)
  a=plotPCA(vsd, intgroup="condition")
  ggsave(pca_plot, plot=a, width=6, height=6, dpi=300)
  # dev.off()

  # valcano plot 
  plotValcano <- function(log2FoldChange, padj, lfc_thres = 1, padj_thres=0.05,main='') {
    tab = data.frame(logFC = log2FoldChange, negLogPval = -log10(padj))#make a data frame with the log2 fold-changes and adjusted p-values
    plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change),
       ylab = expression(-log[10]~pvalue[adj]), main = main, col='grey') 
    signGenesUp = (tab$logFC > lfc_thres   & tab$negLogPval > -log10(padj_thres))
    N_up = sum(signGenesUp,na.rm=T)
    points(tab[signGenesUp, ], pch = 16, cex = 0.6, col = "red")
    signGenesDn = (tab$logFC < -lfc_thres & tab$negLogPval > -log10(padj_thres))
    N_dn = sum(signGenesDn,na.rm=T)
    points(tab[signGenesDn, ], pch = 16, cex = 0.6, col = "blue")
    abline(h = -log10(padj_thres), col = "green3", lty = 2)
    abline(v = c(-lfc_thres, lfc_thres), col = "blue", lty = 2)
    mtext(paste0("FDR = ", padj_thres), side = 4, at = -log10(padj_thres), cex = 0.6, line = 0.5, las = 1)
    mtext(c(paste0("-", lfc_thres, " fold"), paste0("+", lfc_thres, " fold")), side = 3, at = c(-lfc_thres, lfc_thres), cex = 0.6, line = 0.5)
    lims <- par("usr")
    y_len = lims[4] - lims[3]
    x_len = lims[2] - lims[1]
    y4txt = lims[4] - y_len/4
    x4txt = (lims[1] + lims[2]) / 2 + x_len/4
    text(x4txt,y4txt,N_up)
    text(-x4txt,y4txt,N_dn)
  }

  png(valcano_plot,width = 1500, height = 1500,res=300)
  par(mar = c(5, 5, 4, 4))
  plotValcano(final_res$log2FoldChange,final_res$padj,main=sub("^condition_", "",resultsNames(dds)[2]))
  dev.off()


  # Count plots for top 10 DEGs
  # pdf(args$count_plot)
  # par(mfrow=c(2,5))
  # top10 <- head(order(res$pvalue), 10)
  # for (gene in rownames(res)[top10]) {
  #   plotCounts(dds, gene=gene, intgroup="condition", main=gene)
  # }
  # dev.off()
}
