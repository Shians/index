\name{index-package}
\alias{index-package}
\alias{index}
\docType{package}
\title{
\packageTitle{index}
}
\description{
\packageDescription{index}
}
\details{
}
\author{
\packageMaintainer{index}
}
\references{

}

\keyword{ package }
\examples{
library(index)

exon <- readRDS(system.file("extdata/exon_dge.Rds", package = "index"))
intron <- readRDS(system.file("extdata/intron_dge.Rds", package = "index"))
group <- readRDS(system.file("extdata/group.Rds", package = "index"))

x <- index_analysis(exon, intron, group)
plot_voom(x)
plot_lcpm_cor(x)
plot_index(x)
}
