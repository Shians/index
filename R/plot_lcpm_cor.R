plot_lcpm_cor <-
function(index_output) {
    exon_lcpm <- edgeR::cpm(index_output$dges$exon, log = TRUE)
    intron_lcpm <- edgeR::cpm(index_output$dges$intron, log = TRUE)

    par(mfrow = c(1, 2))
    on.exit(par(mfrow = c(1, 1)))

    plot(rowMeans(intron_lcpm), rowMeans(exon_lcpm), xlab = "Intron", ylab = "Exon", col = "#666666")
    abline(a = 0, b = 1, lty = 2, lwd = 1.5)
    boxplot(list(Intron = rowMeans(intron_lcpm), Exon = rowMeans(exon_lcpm)))
}
