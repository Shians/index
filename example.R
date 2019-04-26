
# load data ----
library(limma)
library(edgeR)
library(Homo.sapiens)
# Examine Total RNA samples only
keep <-
    c(paste0("R", 1:3, ".000.Total"), paste0("R", 1:3, ".100.Total"))
load("./counts/Gencode/Genebody_genelevel.RData")
x <- DGEList(counts$counts)
x$genes <- counts$annotation
geneid <- strsplit2(rownames(x), split = "\\.")[, 1]
genes <-
    select(Homo.sapiens,
           keys = geneid,
           columns = "SYMBOL",
           keytype = "ENSEMBL")
m <- match(geneid, genes$ENSEMBL)
x$genes$Symbol <- genes$SYMBOL[m]
x <- x[, keep]
genebody <- x

# Exon counts
load("./counts/Gencode/Exon_genelevel.RData")
x <- DGEList(counts$counts)
x$genes <-
    counts$annotation # Length is actually the length of individual exons
nfeatures <- nchar(gsub(";", "", x$genes$Strand))
x$genes <- x$genes[, c("GeneID", "Length")]
x$genes$Length <- as.numeric(as.character(x$genes$Length))
x$genes$NFeatures <- nfeatures
x <- x[, keep]
exon <- x

# Intron counts
load("./counts/Gencode/Delta_genelevel.RData")
x <- DGEList(counts$counts)
x$genes <-
    as.data.frame(cbind(
        GeneID = rownames(x),
        Length = genebody$genes$Length - exon$genes$Length + 1
    ))
x$genes$Length <- as.numeric(as.character(x$genes$Length))
x$genes$NFeatures <- exon$genes$NFeatures - 1
x <- x[, keep]
intron <- x

group <- as.factor(rep(c("000", "100"), each = 3))

x <- index_analysis(exon, intron, group, design = model.matrix(~0 + group), contrast = c(-1, 1))
plot_voom(x)
plot_lcpm_cor(x)
plot_index(x)
