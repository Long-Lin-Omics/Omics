#!/usr/bin/env Rscript

# Load necessary library
library(rtracklayer)
library(data.table)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if input and output are provided
if (length(args) < 3) {
    stop("Usage: Rscript convert_peaks_to_bigbed.R input_peaks.bed output.bb fa.fai")
}

# Assign input and output file paths
peak_file <- args[1]
bb_file <- args[2]
chrom_sizes_file <- args[3]

if (!file.exists(peak_file)) {
    stop("Input peak file does not exist: ", peak_file)
}

create_empty_bb <- function(output_file, input_file) {
    con <- file(output_file, open = "wb")
    close(con)

    cat(
        "Input peak file is empty; created empty placeholder: ",
        input_file, " → ", output_file, "\n",
        sep = ""
    )

    quit(save = "no", status = 0)
}

if (is.na(file.info(peak_file)$size) ||
    file.info(peak_file)$size == 0) {
    create_empty_bb(bb_file, peak_file)
}

chrom_sizes <- fread(chrom_sizes_file, header = FALSE, col.names = c("chrom", "size","offset","linebases","linewidth"))

# Import peak file (BED/narrowPeak format)
col_names <- c("chrom", "start", "end", "name", "score", "strand",
               "signalValue", "pValue", "qValue", "peak")

peak_data <- fread(peak_file, col.names = col_names, na.strings = ".", fill = TRUE)

if (nrow(peak_data) == 0) {
    create_empty_bb(bb_file, peak_file)
}

gr <- GRanges(seqnames = peak_data$chrom,
              ranges = IRanges(start = peak_data$start, end = peak_data$end),
              strand = "*",
              score = peak_data$score)

chrom_sizes_filtered <- chrom_sizes[chrom_sizes$chrom %in% seqlevels(gr), ]
chrom_sizes_filtered <- chrom_sizes_filtered[match(seqlevels(gr), chrom_sizes_filtered$chrom), ]

seqlengths(gr) <- setNames(chrom_sizes_filtered$size, chrom_sizes_filtered$chrom)


# Export as BigBed
export(gr, bb_file, format = "BigBed")

cat("Conversion complete: ", peak_file, " → ", bb_file, "\n")

