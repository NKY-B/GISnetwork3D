\name{ODmatrix3D}
\alias{ODmatrix3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create OD matrix.
}
\description{
Create OD matrix.
}
\usage{
ODmatrix3D(graph, oID, oV, dID, dV, weight, mode = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
An igraph object created from the function "GISigraph3D"
}
  \item{oID}{
 A numeric vector of the ID of the origin location(s)
}
  \item{oV}{
A numeric vector of the vertex ID of the origin location(s)
}
  \item{dID}{
A numeric vector of the ID of the destination location(s)
}
  \item{dV}{
A numeric vector of the vertex ID of the destination location(s)
}
  \item{weight}{
A numeric vector of weight/ cost for each edge for the least-cost routing.
}
  \item{mode}{
The direction of routing: "all" is undirected, "in" is toward the origin, and "out" is from the origin.
}
}
\details{
This function produce an OD matrix that calculate the cost of the least-cost paths between each origin and destination.
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
