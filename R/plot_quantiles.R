plot_quantiles <- function(
    x,
    output_png,
    probs = c(0.25, 0.5, 0.75),
    width = 3400,
    height = 2000
) {
    # check NA
    na_count <- sum(is.na(x))
    if (na_count > 0) {
        warning(
            sprintf(
                "%d NA values found and removed before plotting.",
                na_count
            )
        )
    }

    # remove NA
    x_clean <- x[!is.na(x)]
    if (length(x_clean) == 0) {
        stop("No non-NA values remaining.")
    }

    # quantiles
    qs <- quantile(x_clean, probs = probs)

    png(
        filename = output_png,
        width = width,
        height = height,
        res = 300
    )

    # plot two figures in one png
    par(mfrow = c(1, 2))

    ## Histogram
    hist(
        x_clean,
        breaks = 30,
        col = "lightblue",
        border = "white",
        main = "Histogram with Quantiles",
        xlab = "Value"
    )
    abline(
        v = qs,
        col = c("red", "blue", "red")[seq_along(qs)],
        lwd = 2,
        lty = 2
    )
    legend(
        "topright",
        legend = paste0(names(qs), " Quantile"),
        col = c("red", "blue", "red")[seq_along(qs)],
        lty = 2,
        lwd = 2,
        bty = "n"
    )

    ## Density
    plot(
        density(x_clean),
        main = "Density Plot with Quantiles",
        lwd = 2
    )
    abline(
        v = qs,
        col = c("red", "blue", "red")[seq_along(qs)],
        lwd = 2,
        lty = 2
    )
    legend(
        "topright",
        legend = paste0(names(qs), " Quantile"),
        col = c("red", "blue", "red")[seq_along(qs)],
        lty = 2,
        lwd = 2,
        bty = "n"
    )

    dev.off()

    invisible(
        list(
            n = length(x_clean),
            na_count = na_count,
            quantiles = quantile(x_clean)
        )
    )
}