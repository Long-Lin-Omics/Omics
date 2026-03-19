#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    paste0(
      "Usage:\n",
      "  Rscript compare_scatter.R file1.csv file2.csv key1,key2 cmp1,cmp2 [output.png]\n\n",
      "Arguments:\n",
      "  key1,key2  : key columns in file1 and file2 to perform inner join\n",
      "  cmp1,cmp2  : key columns in file1 and file2 to perform plotting\n",
      "  output.png : Output plot file (.pdf or .png recommended), defaulted as scatter_merge.png\n"
    )
  )
}

file1 <- args[1]
file2 <- args[2]
key_cols <- trimws(strsplit(args[3], ",", fixed = TRUE)[[1]])
cmp_cols <- trimws(strsplit(args[4], ",", fixed = TRUE)[[1]])
out_file <- if (length(args) >= 5) args[5] else "scatter_merge.png"

if (length(cmp_cols) != 2) {
  stop("compared columns must be two column names from two files, e.g. column_from_file1,column_from_file2")
}

read_table <- function(path) {
  if (!file.exists(path)) stop(paste("File doesn't exist:", path))
  read.table(path, check.names = FALSE, stringsAsFactors = FALSE,h=T)
}

resolve_merged_col <- function(df, preferred,df_index) {
  candidates <- c(
    preferred,
    ifelse(df_index==1,paste0(preferred, ".x"),paste0(preferred, ".y"))
  )
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) {
    stop(paste(
      "Can't find column:", preferred, "in the merged columns: ",paste(names(df), collapse = ", ")
    ))
  }
  hit[1]
}

df1 <- read_table(file1)
df2 <- read_table(file2)

missing1 <- setdiff(c(key_cols[1], cmp_cols[1]), names(df1))
missing2 <- setdiff(c(key_cols[2], cmp_cols[2]), names(df2))

if (length(missing1) > 0) {
  stop(paste(file1, 'missing column(s):', paste(missing1, collapse = ", ")))
}
if (length(missing2) > 0) {
  stop(paste(file2, "missing column(s):", paste(missing2, collapse = ", ")))
}

n1 <- nrow(df1)
n2 <- nrow(df2)

# inner join
merged <- merge(df1, df2, by.x = key_cols[1], by.y = key_cols[2])

n_intersection <- nrow(merged)
n_unique_keys_intersection <- if (n_intersection > 0) nrow(unique(merged[key_cols[1]])) else 0

x_col <- resolve_merged_col(merged, cmp_cols[1], df_index=1)
y_col <- resolve_merged_col(merged, cmp_cols[2], df_index=2)

x <- suppressWarnings(as.numeric(merged[[x_col]]))
y <- suppressWarnings(as.numeric(merged[[y_col]]))

ok <- complete.cases(x, y)
plot_df <- data.frame(x = x[ok], y = y[ok])

n_plot <- nrow(plot_df)

if (n_plot == 0) {
  stop(" the selected columns for plotting are nulls or mising in both files")
}

# output
if (grepl("\\.pdf$", out_file, ignore.case = TRUE)) {
  pdf(out_file, width = 7, height = 7)
  on.exit(dev.off(), add = TRUE)
} else {
  png(out_file, width = 1200, height = 1200, res = 200)
  on.exit(dev.off(), add = TRUE)
}

xr <- range(plot_df$x, na.rm = TRUE)
yr <- range(plot_df$y, na.rm = TRUE)
lim <- range(c(xr, yr), na.rm = TRUE)

plot(
  plot_df$x, plot_df$y,
  xlab = paste0(basename(file1), " : ", cmp_cols[1]),
  ylab = paste0(basename(file2), " : ", cmp_cols[2]),
  main = "Scatter plot after inner join",
  xlim = lim, ylim = lim,
  asp = 1,
  pch = 16, cex = 0.6,col='steelblue'
)
abline(0, 1, lty = 2, lwd = 2)

legend(
  "topleft",
  bty = "n",
  legend = c(
    paste0("#rows(",file1,"): ", n1),
    paste0("#rows(",file2,"): ", n2),
    paste0("#rows(inner join): ", n_intersection),
    paste0("#rows(unique joined keys): ", n_unique_keys_intersection),
    paste0("#rows(plotting(non NA): ", n_plot)
  )
)

mtext(
  paste0(
    "keys = ", paste(key_cols, collapse = ", "),
    " | compare = ", cmp_cols[1], " vs ", cmp_cols[2]
  ),
  side = 3, line = 0.2, cex = 0.8
)

message("output: ", out_file)
message("merged comp columns: ", x_col, " vs ", y_col)