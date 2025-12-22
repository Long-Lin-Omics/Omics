#!/usr/bin/env Rscript

suppressMessages(library(VennDiagram))
suppressMessages(library(dplyr))
suppressMessages(library(grid))

# ==== 1. 获取命令行参数 ====
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 6) {
  stop("Usage: Rscript venn_plot.R file1 file2 output_png set1_name set2_name title_identifier")
}

file1 <- args[1]
file2 <- args[2]
out_png <- args[3]
set1_name <- args[4]
set2_name <- args[5]
title <- args[6]

cat("Input files:\n")
cat("  File1:", file1, "\n")
cat("  File2:", file2, "\n")
cat("Output PNG:", out_png, "\n")
cat("Set names:", set1_name, ",", set2_name, "\n")
cat("Title: Venn Diagram of " , title)

# ==== 2. 自动判断 CSV / TSV ====
read_auto <- function(f) {
  if (grepl("\\.csv$", f, ignore.case = TRUE)) {
    read.csv(f, header = TRUE)
  } else {
    read.table(f, header = TRUE, sep = "\t")
  }
}

df1 <- read_auto(file1)
df2 <- read_auto(file2)

# ==== 3. 可选过滤条件（按需修改） ====
df1_f <- df1 %>% filter(is.na(padj) | padj < 0.05, is.na(log2FoldChange) | abs(log2FoldChange) > 1)
df2_f <- df2 %>% filter(is.na(padj) | padj < 0.05, is.na(log2FoldChange) | abs(log2FoldChange) > 1)

# ==== 4. 第一列作为基因名 ====
genes1 <- df1_f[[1]]
genes2 <- df2_f[[1]]

# ==== 5. 绘制并保存 PNG ====
png(out_png, width = 2000, height = 2000, res = 300)

venn.plot <- venn.diagram(
  x = setNames(
        list(genes1, genes2),
        c(set1_name, set2_name)
      ),
  filename = NULL,
  fill = c("skyblue", "pink"),
  alpha = c(0.5, 0.5),
  cex = 1.5,
  cat.cex = 1.5,
  main = title
)

grid.newpage()
grid.draw(venn.plot)

dev.off()

cat("Venn diagram saved to:", out_png, "\n")