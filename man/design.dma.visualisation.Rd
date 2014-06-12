\name{design.dma.visualisation}
\alias{design.dma.visualisation}
%- Also NEED an '\alis' for EACH other topic documented here.
\title{
Visualize the field layout of trial generated from a modified augmented experiment design
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function is used to visualize the field layout of a modified augmentated experiment design 
}
\usage{
design.dma<-function(enviro=format(Sys.Date(), "\%x"), entries= NULL, nEntries= NULL, chk.names= NULL, nSecChk= NULL, nFieldRows= NULL, nFieldCols= NULL, nRowsPerBlk=NULL, nColsPerBlk=NULL, nChksPerBlk=2,  plot.start=1001, maxPerChks=0.12, fillWithChk=T, minRowBlkDim=2, minColBlkDim=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{enviro}{
%%     ~~Describe \code{enviro} here~~
Label for experiemnt
}
  \item{entries}{
%%     ~~Describe \code{entries} here~~
Accession list 
}
  \item{nEntries}{
%%     ~~Describe \code{num.entries} here~~
The number of accessions
}
  \item{chk.names}{
%%     ~~Describe \code{chk.names} here~~
Check list
}
  \item{nSecChk}{
%%     ~~Describe \code{num.sec.chk} here~~
The number of secondary checks
}
  \item{nFieldRows}{
%%     ~~Describe \code{num.rows} here~~
Number of rows in the trial field
}
  \item{nFieldCols}{
%%     ~~Describe \code{num.cols} here~~
Number of columns in the trial field
}
  \item{plot.start}{
%%     ~~Describe \code{plot.start} here~~
The start index for plot
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
test.ma.acc.72.pchk.1.schk.3.0<-design.dma(entries=c(seq(1,300,1)),chk.names= c(seq(1,4,1)),nFieldRow=10)
#test.ma.acc.72.pchk.1.schk.3.0[[2]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Modified Augmentated Design}
\keyword{Experiment design }% __ONLY ONE__ keyword per line