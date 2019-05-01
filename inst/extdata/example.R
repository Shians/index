library(limma)
library(edgeR)
library(Homo.sapiens)

exon <- readRDS("inst/extdata/exon_dge.Rds")
intron <- readRDS("inst/extdata/intron_dge.Rds")
group <- readRDS("inst/extdata/group.Rds")

x <- index_analysis(exon, intron, group)
plot_voom(x)
plot_lcpm_cor(x)
plot_index(x)
