\name{MAD-package}
\alias{MAD-package}
\alias{MAD}
\docType{package}
\title{MAD}
\description{
%More about what it does (maybe more than one line)
%~~ DA concise (1-5 lines) description of the package ~~
This package is used to generate a modified augmentated experiemental design for plant breeding. 
}
\details{
\tabular{ll}{
Package: \tab MAD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-06-05\cr
LiMcense: 
%What license is it under?
\tab GPL \cr
}
%~~ An overview of how to use the package, including the most important functions ~~
Curently, this package includes a function design.dma. To use this function, you need to give an accession list, a check list and number of rows of trial field in your experiment. Other parameters could be NULL or user options. 
}
\author{
%Who wrote it#
Tyler Tiede in University of Minnesota wrote design.dma fucntion \cr
Aimin Yan in Cornell University reframed Tyler's code to make it into a R package \cr
Maintainer: Aimin Yan <ay247@cornell.edu>
%~~ The author and/or maintainer of the package ~~
}
\references{
%~~ Literature or other references for background information ~~
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ 
%package 
Modified Augmentated Design
}

\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%\code{\link{agricolae}}
}
\examples{
%~~ simple examples of the most important functions ~~
#test file
#simple example to use design.dma
test.ma.acc.72.pchk.1.schk.3.0<-design.dma(entries=c(seq(1,300,1)),chk.names= c(seq(1,4,1)),nFieldRows=10)
#test.ma.acc.72.pchk.1.schk.3.0[[2]]
}