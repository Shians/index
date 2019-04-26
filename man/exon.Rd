\name{exon}
\alias{exon}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("exon")}
\format{
  The format is:
Formal class 'DGEList' [package "edgeR"] with 1 slot
  ..@ .Data:List of 3
  .. ..$ : int [1:58288, 1:6] 2477 0 1272 360 712 0 437 1491 1115 522 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:58288] "ENSG00000000003.14" "ENSG00000000005.5" "ENSG00000000419.12" "ENSG00000000457.13" ...
  .. .. .. ..$ : chr [1:6] "R1.000.Total" "R2.000.Total" "R3.000.Total" "R1.100.Total" ...
  .. ..$ :'data.frame':	6 obs. of  3 variables:
  .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1
  .. .. ..$ lib.size    : num [1:6] 24811479 24098659 24355219 24263842 21483032 ...
  .. .. ..$ norm.factors: num [1:6] 1 1 1 1 1 1
  .. ..$ :'data.frame':	58288 obs. of  3 variables:
  .. .. ..$ GeneID   : chr [1:58288] "ENSG00000000003.14" "ENSG00000000005.5" "ENSG00000000419.12" "ENSG00000000457.13" ...
  .. .. ..$ Length   : num [1:58288] 4535 1610 1207 6883 5967 ...
  .. .. ..$ NFeatures: int [1:58288] 10 7 10 15 32 15 21 7 18 10 ...
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
data(exon)
## maybe str(exon) ; plot(exon) ...
}
\keyword{datasets}
