plot_hist_density <- function(df, col, sample_name = "", title_suffix = "",
                             fill_col = "grey70", line_col = "black") {
  if (is.numeric(col)) {
    colname <- names(df)[col]
  } else if (is.character(col)) {
    colname <- col
  } else {
    stop("col must be column name or index")
  }

  if (!colname %in% names(df)) {
    stop(paste0("Column '", colname, "' not found"))
  }

  x <- df[[col]]
  
  n_reads <- sum(!is.na(x))
  med_val <- median(x, na.rm = TRUE)

  ggplot(df, aes(x = .data[[col]])) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 50,
      fill = fill_col,
      alpha = 0.6
    ) +
    geom_density(color = line_col, linewidth = 1) +
    geom_vline(
      xintercept = med_val,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = med_val,
      y = Inf,
      label = paste0("median = ", round(med_val, 1)),
      vjust = 2,
      hjust = -0.1,
      color = "red",
      size = 3.5
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = paste0("n = ", n_reads),
      hjust = 1.1,
      vjust = 1.5,
      size = 4
    ) +
    labs(
      title = paste(sample_name, title_suffix),
      x = col,
      y = "Density"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank()   
    )
}

