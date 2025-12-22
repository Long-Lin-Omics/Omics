

#!/usr/bin/env Rscript

# 加载包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

# 从命令行读取参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("用法: Rscript plot_IR_polyA.R dorado_polya_per_read read_IR_category output.png")
}

file1 <- args[1]
file2 <- args[2]
outfile <- args[3]

# 读入文件
df2 <- read.table(file1, h=T,stringsAsFactors=F)
df1 <- read.table(file2, h=F,stringsAsFactors=F)
colnames(df1) <- c('readid','IR_category')

# 合并
df <- inner_join(df1, df2, by = "readid")

# 画直方图

df_hist <- df %>%
  group_by(IR_category, pAlength) %>%
  summarise(count = n(), .groups = 'drop')

df_hist <- df_hist[!is.na(df_hist$pAlength),]

# 绘折线图
ggplot(df_hist, aes(x = pAlength, y = count, color = IR_category)) +
  geom_line(linewidth = 1) +
  labs(
    x = "polyA length",
    y = "Count",
    title = "PolyA length histogram (line)"
  ) +
  theme_minimal()
ggsave(outfile, width = 10, height = 6, dpi = 300)