plot_voom <-
function(index_output) {
    capitalise <- function(x) {
        substring(x, 1, 1) <- toupper(substring(x, 1, 1))
        x
    }

    par(mfrow = c(1, 2)); on.exit(par(mfrow = c(1, 1)))

    vooms <- index_output$voom

    for (vname in names(vooms)) {
        v <- vooms[[vname]]
        v$voom.xy$main <- capitalise(vname)
        v$voom.xy$pch <- 20
        do.call(plot, v$voom.xy)
        v$voom.line$col <- "red"
        do.call(lines, v$voom.line)
    }
}
