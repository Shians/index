\name{plot_index}
\alias{plot_index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot_index(index_output)
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
    par(mfrow = c(1, 2))
    on.exit(par(mfrow = c(1, 1)))
    category <- with(index_output$decide.tests, list(`Mixed+-` = (Exon == 
        1) & (Intron == -1), `Mixed-+` = (Exon == -1) & (Intron == 
        1), `Intron-` = (Intron == -1) & (Exon == 0), `Intron+` = (Intron == 
        1) & (Exon == 0), `Exon-` = (Exon == -1) & (Intron == 
        0), `Exon+` = (Exon == 1) & (Intron == 0), `+` = (Exon == 
        1) & (Intron == 1), `-` = (Exon == -1) & (Intron == -1)))
    bar_col <- RColorBrewer::brewer.pal(9, "Paired")
    bar_col <- bar_col[c(3:4, 5:8, 1:2, 9)]
    point_col <- character(length(category[[1]]))
    for (n in names(category)) {
        point_col[category[[n]]] <- n
    }
    point_col <- factor(point_col, levels = c(names(category), 
        ""))
    point_col <- c(bar_col)[point_col]
    barplot(sapply(category, sum), col = bar_col, las = 2, main = "INdEX Categories")
    plot(index_output$tops$intron$logFC, index_output$tops$exon$logFC, 
        pch = 20, col = point_col, main = "logFC", xlab = "Intron logFC", 
        ylab = "Exon logFC")
    abline(a = 0, b = 1, col = "#BBBBBB", lwd = 2)
    par(mfrow = c(1, 1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
