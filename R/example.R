library(limma)
library(edgeR)
library(Homo.sapiens)

exon <- readRDS("data/exon_dge.Rds")
intron <- readRDS("data/intron_dge.Rds")
group <- readRDS("data/group.Rds")

x <- index_analysis(exon, intron, group, design = model.matrix(~0 + group), contrast = c(-1, 1))
plot_voom(x)
plot_lcpm_cor(x)
plot_index(x)
