#' Perform INtron differences to EXons (INdEX) analysis
#'
#' @param exon the DGEList containing exon information
#' @param intron the DGEList containing intron information
#' @param group the vector of groups
#' @param design the design matrix (if absent, design will be constructed from group argument)
#' @param contrast the contrast vector to test (cannot be a matrix)
#' @param p.value the p-value for decideTests
#'
#' @details There are 9 INdEX categories for genes stored in the \code{category} slot of the output:
#' \enumerate{
#'     \item "+": significant upregulation both exon and intron
#'     \item "-": significant downregulation both exon and intron
#'     \item "Exon+": significant upregulation exons only
#'     \item "Exon-": significant downregulation exons only
#'     \item "Intron-": significant upregulation introns only
#'     \item "Intron-": significant downregulation introns only
#'     \item "Mixed+-": significant upregulation in exons and downregulation in introns
#'     \item "Mixed-+": significant upregulation in introns and downregulation in exons
#'     \item "": No significance in any other category
#' }
#'
#' @return list containing a combined decideTest result, INdEX categories, the input dges, voom objects and top tables for exon and introns.
#' @export
#'
#' @examples
index_analysis <- function(exon, intron, group, design = NULL, contrast = NULL, p.value = 0.05) {
    # Input checks ----
    if (isFALSE(all.equal(dim(exon), dim(intron)))) {
        stop("exon and intron DGEs must have the same dimensions")
    }

    if (isFALSE(all.equal(dimnames(exon), dimnames(intron)))) {
        stop("exon and intron DGEs must have the same dimname()")
    }

    dges <- filter_dges(exon, intron, group)
    exon <- dges$exon
    intron <- dges$intron

    if (is.null(design)) {
        message("creating design matrix by 'model.matrix(~group)'")
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

        dt <- set_colnames(
            cbind(dte, dti),
            c("Exon", "Intron")
        )

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

    decide_tests <- get_full_dt(fits$exon, fits$intron, p.value)

    category <- local({
        category_list <- with(decide_tests,
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

        category <- character(length(category_list[[1]]))
        for (n in names(category_list)) {
            category[category_list[[n]]] <- n
        }
        category
    })

    list(
        decide.tests = decide_tests,
        category = category,
        dges = dges,
        voom = vooms,
        tops = tops
    )
}
