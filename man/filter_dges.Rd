\name{filter_dges}
\alias{filter_dges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filter_dges(exon, intron, group)
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
function (exon, intron, group) 
{
    get_lib_sizes <- function(exon, intron) {
        list(exon = colSums(exon$counts), intron = colSums(intron$counts), 
            total = colSums(exon$counts + intron$counts))
    }
    keep.exon <- edgeR::filterByExpr(exon, group = group)
    keep.intron <- edgeR::filterByExpr(intron, group = group)
    keep.exprs <- keep.exon & keep.intron
    exon <- exon[keep.exprs, ]
    intron <- intron[keep.exprs, ]
    lib.sizes <- get_lib_sizes(exon, intron)
    exon$samples$lib.size <- lib.sizes$total
    intron$samples$lib.size <- lib.sizes$total
    exon <- edgeR::calcNormFactors(exon)
    intron <- edgeR::calcNormFactors(intron)
    list(exon = exon, intron = intron)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
