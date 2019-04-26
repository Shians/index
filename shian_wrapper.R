# 19 March 2019 (Last editted 22 March 2019)
# Charity Law
# Compare DE analyses with genebody counts, exon counts and intron counts

# Load packages
library(limma)
library(edgeR)
library(Homo.sapiens)
setwd(
    "/stornext/General/data/user_managed/grpu_law_0/Signalling_in_intronic_region/datasets/Human_celllines_Illumina/Main_hg20/"
)

# Examine Total RNA samples only
keep <-
    c(paste0("R", 1:3, ".000.Total"), paste0("R", 1:3, ".100.Total"))

# Genebody counts
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

# Sample information
group <- as.factor(rep(c("000", "100"), each = 3))
samplenames <- gsub(".Total", "", colnames(x))
nsamples <- length(group)

### FUNCTION START

# Library sizes
lib.genebody <- colSums(genebody$counts)
lib.exon <- colSums(exon$counts)
lib.intron <- colSums(intron$counts)
lib.total <- lib.exon + lib.intron

# Normalised by sequencing depth using total library size
genebody.lcpm <- cpm(genebody, lib.size = lib.total, log = TRUE)
exon.lcpm <- cpm(exon, lib.size = lib.total, log = TRUE)
intron.lcpm <- cpm(intron, lib.size = lib.total, log = TRUE)

# Filter out genes across all count sets
keep.exon <- filterByExpr(exon, group = group)
keep.intron <- filterByExpr(intron, group = group)
table(keep.exon, keep.intron)
keep.exprs <- keep.exon & keep.intron
genebody <- genebody[keep.exprs, ]
exon <- exon[keep.exprs, ]
intron <- intron[keep.exprs, ]
genebody$samples$lib.size <-
    exon$samples$lib.size <- intron$samples$lib.size <- lib.total
genebody.lcpm <- cpm(genebody, log = TRUE)
exon.lcpm <- cpm(exon, log = TRUE)
intron.lcpm <- cpm(intron, log = TRUE)
# Check that uni-exonic genes are removed

# Normalisation (not needed as seen from density plots, but TMM is used anyway)
genebody <- calcNormFactors(genebody)
exon <- calcNormFactors(exon)
intron <- calcNormFactors(intron)

# Log-RPKM values using total library size
genebody.rpkm <- rpkm(genebody, log = F, length = genebody$genes$Length)
exon.rpkm <- rpkm(exon, log = F, length = exon$genes$Length)
intron.rpkm <- rpkm(intron, log = F, length = intron$genes$Length)



# DE model
design <- model.matrix( ~ group)

# limma-voom fit
# Genebody
vg <- voom(genebody, design, plot = TRUE)
fitg <- lmFit(vg, design)
fitg <- eBayes(fitg)
dtg <- decideTests(fitg, p.value = 0.01)[, 2]
tableg <- topTable(fitg,
                   coef = 2,
                   n = Inf,
                   sort.by = "none")
# Exon
ve <- voom(exon, design, plot = TRUE)
fite <- lmFit(ve, design)
fite <- eBayes(fite)
dte <- decideTests(fite, p.value = 0.01)[, 2]
tablee <- topTable(fite,
                   coef = 2,
                   n = Inf,
                   sort.by = "none")
# Intron
vi <- voom(intron, design, plot = TRUE)
fiti <- lmFit(vi, design)
fiti <- eBayes(fiti)
dti <- decideTests(fiti, p.value = 0.01)[, 2]
tablei <- topTable(fiti,
                   coef = 2,
                   n = Inf,
                   sort.by = "none")

# Venn diagram across three comparisons - interesting but difficult to read (in supplementary materials if needed)

dt <- cbind(Genebody = dtg,
            Exon = dte,
            Intron = dti)
geneid <- rownames(dt)

### FUNCTION END

# Genebody results are more similar than intron results when compared to exon results
par(mfrow = c(2, 2))
vennDiagram(dt[, 2:1], main = "DE")
vennDiagram(dt[, 2:3], main = "DE")
library(venneuler)
vc.eg <- vennCounts(dt[, 2:1])
vc.ei <- vennCounts(dt[, 2:3])
plot(venneuler(c(
    "Genebody" = vc.eg[2, "Counts"],
    "Exon" = vc.eg[3, "Counts"],
    "Genebody&Exon" = vc.eg[4, "Counts"]
)),
col = c("limegreen", "dodgerblue"))
plot(venneuler(c(
    "Intron" = vc.ei[2, "Counts"],
    "Exon" = vc.ei[3, "Counts"],
    "Intron&Exon" = vc.ei[4, "Counts"]
)),
col = c("magenta", "dodgerblue"))

# Split results into directions
vennDiagram(
    dt[, 2:1],
    include = c("up", "down"),
    counts.col = c("red", "blue"),
    main = "DE same dir"
)
vennDiagram(
    dt[, 2:3],
    include = c("up", "down"),
    counts.col = c("red", "blue"),
    main = "DE same dir"
)
