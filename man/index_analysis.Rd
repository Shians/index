\name{index_analysis}
\alias{index_analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
index_analysis(exon, intron, group, design = NULL, contrast = NULL, p.value = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exon}{
%%     ~~Describe \code{exon} here~~
}
  \item{intron}{
%%     ~~Describe \code{intron} here~~
}
  \item{group}{
%%     ~~Describe \code{group} here~~
}
  \item{design}{
%%     ~~Describe \code{design} here~~
}
  \item{contrast}{
%%     ~~Describe \code{contrast} here~~
}
  \item{p.value}{
%%     ~~Describe \code{p.value} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (exon, intron, group, design = NULL, contrast = NULL, 
    p.value = 0.01) 
{
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
        limma::topTable(fit, coef = ncol(fit), number = Inf, 
            sort.by = "none")
    }
    get_full_dt <- function(fite, fiti, p.value) {
        dte <- get_dt(fite, p.value = p.value)
        dti <- get_dt(fiti, p.value = p.value)
        dt <- cbind(Exon = dte, Intron = dti)
        geneid <- rownames(dt)
        data.frame(dt)
    }
    dges <- list(exon = exon, intron = intron)
    vooms <- lapply(dges, get_voom)
    fits <- lapply(vooms, get_fit, contrast = contrast)
    tops <- lapply(fits, get_top)
    list(decide.tests = get_full_dt(fits$exon, fits$intron, p.value), 
        dges = dges, voom = vooms, tops = tops)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
