\name{ODmatrix3D_summary}
\alias{ODmatrix3D_summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing OD matrix.
}
\description{
Summarizing OD matrix.
}
\usage{
ODmatrix3D_summary(frame, oID, oV, dID, dV, threshold = 500, avaiW)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frame}{
A dataframe of OD matrix.
}
  \item{oID}{
A numeric vector of the IDs of the origin locations
}
  \item{oV}{
A numeric vector of the vertices IDs of the origin locations
}
  \item{dID}{
A numeric vector of the IDs of the destination locations
}
  \item{dV}{
A numeric vector of the vertices IDs of the destination locations
}
  \item{threshold}{
A number indicating the cost threshold for defining the catchment for calculating the availability of destinations.
}
  \item{avaiW}{
A numeric vector defining the weight for each destinatioin for calculating the availability within a threshold travel cost from each origin location. Weight of one for a destination refers to counting a destination as one if it is within the catchment of the origin within a cost budget.
}
}
\details{
This function summarizes the OD matrix produced from the function "ODmatrix3D". This function calculates the minimum cost, mean cost, median cost, and availability for each origin to many destinations. The availability is calculated by the sum of "avaiW" of each destination for each origin within a "threshold" search. avaiW is the weight for each destination. For instance, one means that the destination is counted as 1. The threshold is the cost <= which the destination is counted as available within a certain level of cost.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
