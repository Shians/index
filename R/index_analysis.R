index_analysis <-
function(exon, intron, group, design = NULL, contrast = NULL, p.value = 0.01) {
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

        limma::voom(dge, design, plot = FALSE, save.plot = TRUE)
    }

    get_fit <- function(v, contrast) {
        stopifnot(is(v, "EList"))

        fit <- limma::lmFit(v, design)

        if (!is.null(contrast)) {
            fit <- limma::contrasts.fit(fit, contrasts = contrast)
        }

        fit <- limma::eBayes(fit)
        return(fit)
    }

    get_dt <- function(fit, p.value) {
        stopifnot(is(p.value, "numeric"))

        limma::decideTests(fit, p.value = p.value)[, ncol(fit)]
    }

    get_top <- function(fit) {
        limma::topTable(fit, coef = ncol(fit), number = Inf, sort.by = "none")
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
