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


plot_hist_density_grouped <- function(df, value_col, group_col,sample_name = "", title_suffix = "",
                             fill_cols = c("steelblue","darkorange"), line_cols = c("darkblue","brown") ){ 
  df2 <- df |>
      dplyr::filter(!is.na(.data[[value_col]]), !is.na(.data[[group_col]]))

  if ( length(unique(df2[[group_col]])) > 2 ) {
    warning('The group column has more than two unique elements!')
  }

  df_median <- df2 %>%
    group_by(.data[[group_col]]) %>%
    summarise(n=n(),median_val = median(.data[[value_col]], na.rm = TRUE))
  
  mm = median(df_median$median_val)
  ggplot(df2, aes(x = .data[[value_col]], )) +
    geom_histogram(aes(y = after_stat(density), fill = .data[[group_col]]), alpha = 0.5, bins = 50, position = "identity") +
    geom_density(aes(color = .data[[group_col]]),linewidth = 1) +   # 轮廓线
    theme_bw() +
    labs(title = paste(sample_name, title_suffix),
         x = paste(df_median[[1]], collapse = " vs "),
         y = "Density") + 
    annotate("text",
           x = Inf, y = Inf,
           label = paste(paste0("N(",df_median[[group_col]],"): "),  df_median[['n']], collapse='\n'),
           hjust = 1.1, vjust = 1.5,
           size = 5) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),   # 去网格
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 14)
    )+
    geom_vline(data = df_median,
             aes(xintercept = median_val, color = .data[[group_col]]),
             linetype = "dashed", linewidth = 1) +
    geom_text(data = df_median,
            aes(x = median_val,
                y = Inf,
                label = round(median_val, 1),
                color = .data[[group_col]],
                hjust = ifelse(median_val < mm, 1.2, -0.2)),
            vjust = 1.5,
            size = 4,
            show.legend = FALSE)+
    scale_fill_manual(values = fill_cols) +
    scale_color_manual(values = line_cols)

}

