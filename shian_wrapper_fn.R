# Based on DE_analyses.R
filter_dges <- function(exon, intron, group) {
    get_lib_sizes <- function(exon, intron) {
        list(
            exon = colSums(exon$counts),
            intron = colSums(intron$counts),
            total = colSums(exon$counts + intron$counts)
        )
    }

    # Filter out genes across all count sets
    keep.exon <- filterByExpr(exon, group = group)
    keep.intron <- filterByExpr(intron, group = group)
    keep.exprs <- keep.exon & keep.intron

    exon <- exon[keep.exprs, ]
    intron <- intron[keep.exprs, ]

    lib.sizes <- get_lib_sizes(exon, intron)
    exon$samples$lib.size <- lib.sizes$total
    intron$samples$lib.size <- lib.sizes$total

    exon <- calcNormFactors(exon)
    intron <- calcNormFactors(intron)

    list(
        exon = exon,
        intron = intron
    )
}

index_analysis <- function(exon, intron, group, design = NULL, contrast = NULL, p.value = 0.01) {
    # Input checks ----
    if (!all.equal(dim(exon), dim(intron))) {
        stop("exon and intron DGEs must have the same dimensions")
    }

    if (!all.equal(dimnames(exon), dimnames(intron))) {
        stop("exon and intron DGEs must have the same dimname()")
    }

    dges <- filter_dges(exon, intron, group)
    exon <- dges$exon
    intron <- dges$intron

    if (is.null(design)) {
        message("creating desing matrix by 'model.matrix(~group)'")
        design <- model.matrix(~group)
    }

    if (nrow(design) != ncol(exon)) {
        stop("'design' must have the same number of rows as columns in DGEList objects")
    }

    if (!is.null(contrast)) {
        if (!is(contrast, "numeric")) {
            stop("'contrast' argument must be numeric values")
        }

        if (length(contrast) != 1 && length(contrast) != ncol(design)) {
            stop("'contrast' argument must be length 1 or equal to the number of columns in design")
        }
    }

    # Helper functions ----
    get_voom <- function(dge) {
        stopifnot(is(dge, "DGEList"))

        voom(dge, design, plot = FALSE, save.plot = TRUE)
    }

    get_fit <- function(v, contrast) {
        stopifnot(is(v, "EList"))

        fit <- lmFit(v, design)

        if (is.null(contrast)) {
            fit <- eBayes(fit)
        } else {
            fit <- contrasts.fit(fit, contrasts = contrast)
            fit <- eBayes(fit)
        }

        return(fit)
    }

    get_dt <- function(fit, p.value) {
        stopifnot(is(p.value, "numeric"))

        decideTests(fit, p.value = p.value)[, ncol(fit)]
    }

    get_top <- function(fit) {
        topTable(fit, coef = ncol(fit), number = Inf, sort.by = "none")
    }

    get_full_dt <- function(fite, fiti, p.value) {
        dte <- get_dt(fite, p.value = p.value)
        dti <- get_dt(fiti, p.value = p.value)

        dt <- cbind(Exon = dte, Intron = dti)
        geneid <- rownames(dt)
        data.frame(dt)
    }

    # Main computations ----
    dges <- list(
        exon = exon,
        intron = intron
    )

    vooms <- lapply(
        dges,
        get_voom
    )

    fits <- lapply(
        vooms,
        get_fit,
        contrast = contrast
    )

    tops <- lapply(
        fits,
        get_top
    )

    list(
        decide.tests = get_full_dt(fits$exon, fits$intron, p.value),
        dges = dges,
        voom = vooms,
        tops = tops
    )
}

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

# run function ----
plot_voom <- function(index_output) {
    par(mfrow = c(1, 2)); on.exit(par(mfrow = c(1, 1)))

    vooms <- index_output$voom

    for (vname in names(vooms)) {
        v <- vooms[[vname]]
        v$voom.xy$main <- vname
        v$voom.xy$pch <- 20
        do.call(plot, v$voom.xy)
        v$voom.line$col <- "red"
        do.call(lines, v$voom.line)
    }
}

plot_lcpm_cor <- function(index_output) {
    exon_lcpm <- cpm(index_output$dges$exon, log = TRUE)
    intron_lcpm <- cpm(index_output$dges$intron, log = TRUE)

    par(mfrow = c(1, 2))
    on.exit(par(mfrow = c(1, 1)))

    plot(rowMeans(intron_lcpm), rowMeans(exon_lcpm))
    boxplot(list(Intron = rowMeans(intron_lcpm), Exon = rowMeans(exon_lcpm)))
}

plot_index <- function(index_output) {
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

x <- index_analysis(exon, intron, group, design = model.matrix(~0 + group), contrast = c(-1, 1))
plot_voom(x)
plot_lcpm_cor(x)
plot_index(x)
