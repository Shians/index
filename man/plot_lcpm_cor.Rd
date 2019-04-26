\name{plot_lcpm_cor}
\alias{plot_lcpm_cor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_lcpm_cor(index_output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{index_output}{
%%     ~~Describe \code{index_output} here~~
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
function (index_output) 
{
    exon_lcpm <- edgeR::cpm(index_output$dges$exon, log = TRUE)
    intron_lcpm <- edgeR::cpm(index_output$dges$intron, log = TRUE)
    par(mfrow = c(1, 2))
    on.exit(par(mfrow = c(1, 1)))
    plot(rowMeans(intron_lcpm), rowMeans(exon_lcpm), xlab = "Intron", 
        ylab = "Exon", col = "#666666")
    abline(a = 0, b = 1, lty = 2, lwd = 1.5)
    boxplot(list(Intron = rowMeans(intron_lcpm), Exon = rowMeans(exon_lcpm)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
