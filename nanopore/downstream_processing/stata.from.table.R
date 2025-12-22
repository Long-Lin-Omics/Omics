
args <- commandArgs(trailingOnly=TRUE)
input_file <- args[1]

library(data.table)

data <- fread(input_file, header = TRUE)

if (ncol(data) == 5) {
  setnames(data, c("read_id", "transcript", "MAPQ", "barcode", "polya_length"))
} else if (ncol(data) == 6) {
  setnames(data, c("read_id", "transcript", "polya_length", "barcode", "CI", "Sample"))
} else {
  stop("Unexpected number of columns: ", ncol(data))
}

library(rtracklayer)
gtf <- import("/ddn/gs1/project/nextgen/post/hug4/LongLin/rnaseq/reference/Mus_musculus.GRCm38.102.gtf")
gtf <- as.data.frame(gtf)

data$transcript_clean <- sub("\\.\\d+$", "", data$transcript)
gtf$transcript_clean <- sub("\\.\\d+$", "", gtf$transcript_id)
idx <- match(data$transcript_clean, gtf$transcript_clean)
data$gene_or_transcript <- data$transcript
data$gene_or_transcript[!is.na(idx)] <- gtf$gene_name[idx[!is.na(idx)]]

library(dplyr)

data$barcode[is.na(data$barcode)]  = 'Undemuxed'

total_reads = dim(data)[1]
total_mapped = 0
total_mapped_A = 0
total_unmapped = 0
total_unmapped_A = 0
to_print = 'barcode\ttotoal_reads\tmapped_rate\tmapped_reads_count\tmapped_polyA_count\tmapped_polyA_rate\tunmapped_reads_count\tunmmapped_polyA_count\tunmmapped_polyA_rate\n'

for (bc in sort(unique(data$barcode))) {
  barcode_subset <- data %>% filter(barcode == bc)

  if (nrow(barcode_subset) == 0) next

  gene_counts <- barcode_subset %>%
  group_by(gene_or_transcript) %>%
  summarise(
    read_count = n(),  # unique read 数
    reads_count_with_poly_a = sum(polya_length > 0,na.rm=T),
    mean_polya_length = ifelse(reads_count_with_poly_a > 0,
                               mean(polya_length[polya_length > 0], na.rm = TRUE),
                               NA_real_),
    quantile_0 = ifelse(reads_count_with_poly_a > 0,
                        quantile(polya_length[polya_length > 0], 0, na.rm = TRUE),
                        NA_real_),
    quantile_25 = ifelse(reads_count_with_poly_a > 0,
                         quantile(polya_length[polya_length > 0], 0.25, na.rm = TRUE),
                         NA_real_),
    quantile_50 = ifelse(reads_count_with_poly_a > 0,
                         quantile(polya_length[polya_length > 0], 0.5, na.rm = TRUE),
                         NA_real_),
    quantile_75 = ifelse(reads_count_with_poly_a > 0,
                         quantile(polya_length[polya_length > 0], 0.75, na.rm = TRUE),
                         NA_real_),
    quantile_100 = ifelse(reads_count_with_poly_a > 0,
                          quantile(polya_length[polya_length > 0], 1, na.rm = TRUE),
                          NA_real_)
  ) %>%
  arrange(desc(read_count))
  
  mapped_total = sum(gene_counts$read_count[gene_counts$gene_or_transcript!='*'])
  mapped_A_total = sum(gene_counts$reads_count_with_poly_a[gene_counts$gene_or_transcript!='*'])
  unmapped_total = sum(gene_counts$read_count[gene_counts$gene_or_transcript=='*'])
  unmapped_A_total = sum(gene_counts$reads_count_with_poly_a[gene_counts$gene_or_transcript=='*'])
  
  total_mapped = total_mapped + mapped_total
  total_mapped_A = total_mapped_A + mapped_A_total
  total_unmapped = total_unmapped + unmapped_total
  total_unmapped_A = total_unmapped_A + unmapped_A_total

  to_print = paste0(to_print, paste(bc, dim(barcode_subset)[1], mapped_total/dim(barcode_subset)[1], mapped_total, mapped_A_total, mapped_A_total/mapped_total, unmapped_total,  unmapped_A_total, unmapped_A_total/unmapped_total, sep='\t'),'\n')

  gene_count_file <- paste0( bc, "_gene_counts.tsv")
  write.table(gene_counts, file = gene_count_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Saved gene counts for barcode", bc, "to", gene_count_file, "\n")
}


cat('totoal_reads\tmapped_rate\tmapped_reads_count\tmapped_polyA_count\tmapped_polyA_rate\tunmapped_reads_count\tunmmapped_polyA_count\tunmmapped_polyA_rate\n',file='summary.txt')
cat(paste(total_reads, total_mapped/total_reads, total_mapped, total_mapped_A, total_mapped_A/total_mapped, total_unmapped, total_unmapped_A, total_unmapped_A/total_unmapped, sep='\t'),'\n',file='summary.txt',append=T)
cat('\n',file='summary.txt',append=T)
cat(to_print,file='summary.txt',append=T)


library(dplyr)
library(ggplot2)
library(patchwork)

# 准备两个空的 data.frame
polyA_freq_all_allreads <- data.frame()
polyA_freq_all_mapped <- data.frame()

for (bc in sort(unique(data$barcode))) {
  
  # 情况1：所有 reads
  barcode_subset_all <- data %>% filter(barcode == bc, polya_length > 0)
  if (nrow(barcode_subset_all) > 0) {
    freq_table <- barcode_subset_all %>%
      group_by(polya_length) %>%
      summarise(count = n()) %>%
      mutate(barcode = bc) 
    polyA_freq_all_allreads <- rbind(polyA_freq_all_allreads, freq_table)
  }
  
  # 情况2：transcript != "*"
  barcode_subset_mapped <- data %>% filter(barcode == bc, transcript != "*", polya_length > 0)
  if (nrow(barcode_subset_mapped) > 0) {
    freq_table <- barcode_subset_mapped %>%
      group_by(polya_length) %>%
      summarise(count = n()) %>%
      mutate(barcode = bc) 
    polyA_freq_all_mapped <- rbind(polyA_freq_all_mapped, freq_table)
  }
}

median_by_barcode <- polyA_freq_all_allreads %>%
  group_by(barcode) %>%
  summarize(median_polya = median(polya_length, na.rm = TRUE))

cat('\nMedians from All Reads (polya_length > 0):\n',file='summary.txt',append=T)
cat(capture.output(print(median_by_barcode)), file = 'summary.txt', append = TRUE, sep = "\n")

median_by_barcode <- polyA_freq_all_mapped %>%
  group_by(barcode) %>%
  summarize(median_polya = median(polya_length, na.rm = TRUE))

cat('Medians from Mapped Reads (transcript != \'*\'):\n',file='summary.txt',append=T)
cat(capture.output(print(median_by_barcode)), file = 'summary.txt', append = TRUE, sep = "\n")


# 绘图
p1 <- ggplot(polyA_freq_all_allreads, aes(x = polya_length, y = count, color = barcode)) +
  geom_line(size = 0.7) +
  theme_bw() +
  labs(title = "All Reads (polya_length > 0)", x = "PolyA Length", y = "Count")
  # + geom_vline(data = median_by_barcode,
  #            aes(xintercept = median_polya, color = barcode),
  #            linetype = "dashed")



p2 <- ggplot(polyA_freq_all_mapped, aes(x = polya_length, y = count, color = barcode)) +
  geom_line(size = 0.7) +
  theme_bw() +
  labs(title = "Mapped Reads (transcript != '*')", x = "PolyA Length", y = "Count")

# 保存到同一个 PNG
png("polyA_length_two_conditions.png", width = 2500, height = 1500, res = 300)
p2 / p1   # patchwork 用 / 表示上下排列
dev.off()

cat("Plot saved to polyA_length_two_conditions.png\n")




