plot_index <-
function(index_output) {
    par(mfrow = c(1, 2))
    on.exit(par(mfrow = c(1, 1)))

    category <- with(index_output$decide.tests,
         list(
            `Mixed+-` = (Exon == 1) & (Intron == -1),
            `Mixed-+` = (Exon == -1) & (Intron == 1),
            `Intron-` = (Intron == -1) & (Exon == 0),
            `Intron+` = (Intron == 1) & (Exon == 0),
            `Exon-` = (Exon == -1) & (Intron == 0),
            `Exon+` = (Exon == 1) & (Intron == 0),
            `+` = (Exon == 1) & (Intron == 1),
            `-` = (Exon == -1) & (Intron == -1)
        )
    )

    bar_col <- RColorBrewer::brewer.pal(9, "Paired")
    bar_col <- bar_col[c(3:4, 5:8, 1:2, 9)]

    # assign colours to points, every point should only belong to one category
    point_col <- character(length(category[[1]]))
    for (n in names(category)) {
        point_col[category[[n]]] <- n
    }
    point_col <- factor(point_col, levels = c(names(category), ""))
    point_col <- c(bar_col)[point_col]

    barplot(
        sapply(category, sum),
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
