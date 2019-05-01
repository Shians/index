#' Create main INdEX analysis plots
#'
#' @param index_output the output from running index_analysis
#'
#' @return None
#' @export
#'
#' @examples
plot_index <-
function(index_output) {
    par(mfrow = c(1, 2))
    on.exit(par(mfrow = c(1, 1)))

    bar_col <- RColorBrewer::brewer.pal(9, "Paired")
    bar_col <- bar_col[c(3:4, 5:8, 1:2, 9)]

    # assign colours to points, every point should only belong to one category
    categories <- index_output$category
    categories <- factor(
        categories,
        levels = c("Mixed+-", "Mixed-+", "Intron-", "Intron+", "Exon-", "Exon+", "+", "-", "")
    )
    point_col <- c(bar_col)[categories]

    category_counts <- sapply(
        setdiff(levels(categories), ""),
        function(x) {
            sum(categories == x)
        }
    )

    barplot(
        category_counts,
        col = bar_col,
        las = 2,
        main = "INdEX Categories"
    )

    plot(
        index_output$tops$intron$logFC,
        index_output$tops$exon$logFC,
        pch = 20,
        col = point_col,
        main = "logFC",
        xlab = "Intron logFC",
        ylab = "Exon logFC"
    )
    abline(a = 0, b = 1, col = "#BBBBBB", lwd = 2)
    par(mfrow = c(1, 1))
}
