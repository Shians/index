filter_dges <-
function(exon, intron, group) {
    get_lib_sizes <- function(exon, intron) {
        list(
            exon = colSums(exon$counts),
            intron = colSums(intron$counts),
            total = colSums(exon$counts + intron$counts)
        )
    }

    # Filter out genes across all count sets
    keep.exon <- edgeR::filterByExpr(exon, group = group)
    keep.intron <- edgeR::filterByExpr(intron, group = group)
    keep.exprs <- keep.exon & keep.intron

    exon <- exon[keep.exprs, ]
    intron <- intron[keep.exprs, ]

    lib.sizes <- get_lib_sizes(exon, intron)
    exon$samples$lib.size <- lib.sizes$total
    intron$samples$lib.size <- lib.sizes$total

    exon <- edgeR::calcNormFactors(exon)
    intron <- edgeR::calcNormFactors(intron)

    list(
        exon = exon,
        intron = intron
    )
}
