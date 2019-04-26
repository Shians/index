\name{intron}
\alias{intron}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("intron")}
\format{
  The format is:
Formal class 'DGEList' [package "edgeR"] with 1 slot
  ..@ .Data:List of 3
  .. ..$ : num [1:58288, 1:6] 74 9 0 533 971 1 226 0 590 0 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:58288] "ENSG00000000003.14" "ENSG00000000005.5" "ENSG00000000419.12" "ENSG00000000457.13" ...
  .. .. .. ..$ : chr [1:6] "R1.000.Total" "R2.000.Total" "R3.000.Total" "R1.100.Total" ...
  .. ..$ :'data.frame':	6 obs. of  3 variables:
  .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1
  .. .. ..$ lib.size    : num [1:6] 9090869 8896734 9069887 8776476 8546322 ...
  .. .. ..$ norm.factors: num [1:6] 1 1 1 1 1 1
  .. ..$ :'data.frame':	58288 obs. of  3 variables:
  .. .. ..$ GeneID   : Factor w/ 58288 levels "ENSG00000000003.14",..: 1 2 3 4 5 6 7 8 9 10 ...
  .. .. ..$ Length   : num [1:58288] 8349 13475 22483 37755 186108 ...
  .. .. ..$ NFeatures: num [1:58288] 9 6 9 14 31 14 20 6 17 9 ...
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(intron)
## maybe str(intron) ; plot(intron) ...
}
\keyword{datasets}
