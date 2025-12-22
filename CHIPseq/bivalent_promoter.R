#!/usr/bin/env Rscript

# Load libraries
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(tools)

# Receive arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 6) {
  stop("Usage: Rscript promoter_overlap_double_plot.R <INO80_file> <K27me3_file> <K4me3_file> <promoter_distance> <output_file> <output_file>")
}

ino80_file <- args[1]
k27me3_file <- args[2]
k4me3_file <- args[3]
promoter_range <- as.numeric(args[4])
output_file <- args[5]
output_file2 <- args[6]

# Read annotated HOMER files
ino80 <- read_tsv(ino80_file, col_types = cols())
k27me3 <- read_tsv(k27me3_file, col_types = cols())
k4me3  <- read_tsv(k4me3_file,  col_types = cols())

# Filter to promoter regions
ino80_prom <- ino80 %>% filter(abs(`Distance to TSS`) <= promoter_range) %>% mutate(Gene = `Gene Name`)
k27_prom   <- k27me3 %>% filter(abs(`Distance to TSS`) <= promoter_range) %>% mutate(Gene = `Gene Name`)
k4_prom    <- k4me3  %>% filter(abs(`Distance to TSS`) <= promoter_range) %>% mutate(Gene = `Gene Name`)

# Check overlaps for INO80 promoter regions
ino80_prom <- ino80_prom %>%
  mutate(K27me3 = Gene %in% k27_prom$Gene,
         K4me3  = Gene %in% k4_prom$Gene) %>%
  mutate(Category = case_when(
    K27me3 & K4me3             ~ "Bivalent",
    !K27me3 & K4me3            ~ "H3K4me3 only",
    K27me3 & !K4me3            ~ "H3K27me3 only",
    TRUE                       ~ "No overlap"
  ))

# Summary counts for INO80 promoters
summary_counts <- ino80_prom %>%
  count(Category)

total_promoters <- nrow(ino80_prom)
bivalent_regions <- summary_counts %>% filter(Category == "Bivalent") %>% pull(n)
if(length(bivalent_regions) == 0) bivalent_regions <- 0  # handle if no bivalent regions

# Identify bivalent promoter genes (K27me3 + K4me3 overlap)
bivalent_promoters <- intersect(k27_prom$Gene, k4_prom$Gene) %>% unique()

# Check how many of them are also bound by INO80
bivalent_df <- tibble(Gene = bivalent_promoters) %>%
  mutate(INO80 = Gene %in% ino80_prom$Gene) %>%
  count(INO80) %>%
  mutate(Category = ifelse(INO80, "Bound by INO80", "Not bound by INO80"))

total_bivalent_promoters <- length(bivalent_promoters)

# Open graphic device based on file extension
file_ext <- file_ext(output_file)
if (file_ext == "pdf") {
  pdf(output_file, width = 8, height = 8)
} else if (file_ext %in% c("png", "tiff", "jpeg", "jpg")) {
  png(output_file, width = 2000, height = 2000, res = 300)
} else {
  stop("Unsupported output file type. Please use .pdf, .png, .jpg, or .tiff")
}

# First plot: INO80 promoter overlaps
print(summary_counts)
p1 <- ggplot(summary_counts, aes(x = paste0("INO80 Promoters(n=",sum(summary_counts$n),")"), y = n, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4) +
  # annotate("text", x = 1, y = total_promoters + 10,
  #          label = paste("Total promoters:", total_promoters), size = 5) +
  # annotate("text", x = 1, y = total_promoters + 50,
  #          label = paste("Bivalent regions:", bivalent_regions), size = 5) +
  scale_fill_manual(values = c("Bivalent" = "#E69F00", "H3K4me3 only" = "#56B4E9",
                               "H3K27me3 only" = "#009E73", "No overlap" = "grey70")) +
  theme_minimal() +
  ylab("Number of Promoter Regions") +
  xlab("") +
  ggtitle(paste("Promoter Overlap within ±", promoter_range, "bp of TSS (INO80 Promoters)"))

print(p1)
dev.off()

file_ext <- file_ext(output_file2)
if (file_ext == "pdf") {
  pdf(output_file2, width = 8, height = 8)
} else if (file_ext %in% c("png", "tiff", "jpeg", "jpg")) {
  png(output_file2, width = 2000, height = 2000, res = 300)
} else {
  stop("Unsupported output file type. Please use .pdf, .png, .jpg, or .tiff")
}
# Second plot: Bivalent promoters overlapped by INO80
print(bivalent_df)
p2 <- ggplot(bivalent_df, aes(x = paste0("Bivalent Promoters(n=",sum(bivalent_df$n),")"), y = n, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            color = "white", size = 5) +
  # annotate("text", x = 1, y = total_bivalent_promoters + 5,
  #          label = paste("Total bivalent promoters:", total_bivalent_promoters), size = 5) +
  scale_fill_manual(values = c("Bound by INO80" = "#D55E00", "Not bound by INO80" = "grey70")) +
  theme_minimal() +
  ylab("Number of Promoter Regions") +
  xlab("") +
  ggtitle(paste("INO80 Binding at Bivalent Promoters within ±", promoter_range, "bp"))

print(p2)

# Close device
dev.off()

cat("Plots saved to", output_file, "\n")
