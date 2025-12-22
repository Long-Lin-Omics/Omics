#!/usr/bin/env Rscript

# 从命令行获取参数：bed 输入文件 & png 输出文件
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_bed_len_hist.R input.bed output.png")
}

bed_file <- args[1]
png_file <- args[2]

# 读入 bed 文件（假设至少有 3 列）
bed <- read.table(bed_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# 计算区间长度：第三列 - 第二列
len <- bed[[3]] - bed[[2]]

# 可以去掉非正长度（如果有脏数据）
len <- len[len > 0]

# 画图到 png
png(png_file, width = 800, height = 600)
hist(
  len,
  main = paste0("Interval Length Distribution (n=", length(len), ')'),
  xlab = "Interval length (V3 - V2)",
  ylab = "Frequency",
  breaks = 50
)
dev.off()
