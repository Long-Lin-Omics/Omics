#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggpubr)
  library(stringr)
  library(data.table)
})

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript plot_peak_distances_multi_inlist.R <annotated_peaks.txt> <output.png> <tag for title> <gene_list1.txt> [<gene_list2.txt> ...]")
}

annot_file <- args[1]
output_png <- args[2]
tag <- args[3]
gene_list_files <- args[4:length(args)]

# Read in HOMER annotation file
annot <- fread(annot_file, h=T)

# Check for required columns
required_cols <- c("Gene Name", "Distance to TSS")
if (!all(required_cols %in% colnames(annot))) {
  stop("Annotated peaks file must contain columns: 'Gene Name' and 'Distance to TSS'")
}

# Prepare 'All' group data
annot <- annot %>%
  mutate(Distance = abs(`Distance to TSS`),
            Group = "All",Class='InList')

plot_data <- annot
group_labels = c('All')
# For each gene list file
for (gene_file in gene_list_files) {
  gene_list <- read_lines(gene_file)
  group_label <- str_remove(basename(gene_file), "\\.txt$|\\.bed$|\\.list$")  # e.g. up_DEGs
  group_labels=c(group_labels,group_label)
  # Annotate peaks as 'InList' or 'Other' for this gene list
  subset_data <- annot %>%
    mutate(Class = ifelse(`Gene Name` %in% gene_list, "InList", "Other"),
           Group = group_label) 

  # Combine into plot_data
  plot_data <- bind_rows(plot_data, subset_data)
}

plot_data$Group<- factor(plot_data$Group, levels = group_labels)

inlist_counts <- plot_data %>%
  filter(Class == "InList") %>%
  group_by(Group) %>%
  summarise(n = n())
print(inlist_counts)

# 箱线图 + 标数字
p <- ggplot(plot_data, aes(x = Group, y = log10(Distance), fill = Class)) +
  geom_boxplot() +
  geom_text(data = inlist_counts, aes(x = Group, y = max(log10(plot_data$Distance)) * 1.05, label = paste0("n(InList)=", n)),
            inherit.aes = FALSE, vjust = 0, size = 3.5) +
  labs(
    title = paste0("Distance of Peaks (", tag, ") to Nearest TSS"),
    y = "Distance to TSS (bp, log10)",
    x = "Groups"
  ) 
#   theme_minimal() +
#   theme(
#     text = element_text(size = 14),
#     axis.text.x = element_text(angle = 30, hjust = 1),
#     legend.position = "none"
#   ) +
#   scale_fill_brewer(palette = "Set2")

# Save to file
ggsave(output_png, p, width = 9, height = 5, dpi = 300)

