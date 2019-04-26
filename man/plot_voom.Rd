\name{plot_voom}
\alias{plot_voom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_voom(index_output)
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
    capitalise <- function(x) {
        substring(x, 1, 1) <- toupper(substring(x, 1, 1))
        x
    }
    par(mfrow = c(1, 2))
    on.exit(par(mfrow = c(1, 1)))
    vooms <- index_output$voom
    for (vname in names(vooms)) {
        v <- vooms[[vname]]
        v$voom.xy$main <- capitalise(vname)
        v$voom.xy$pch <- 20
        do.call(plot, v$voom.xy)
        v$voom.line$col <- "red"
        do.call(lines, v$voom.line)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
