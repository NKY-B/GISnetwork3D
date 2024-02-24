\name{Shortest_3Dpaths_for_pairs}
\alias{Shortest_3Dpaths_for_pairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create one or more paths for each pair of origin and destination.
}
\description{
Create one or more paths for each pair of origin and destination.
}
\usage{
Shortest_3Dpaths_for_pairs(graph, edges, oID, oV, dID, dV, weight, mode = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
An igraph object created from the function "GISigraph3D"
}
  \item{edges}{
A sf line layer with the fields of "from" and "to" to indicate the IDs of vertices for the start and end of each line segment. This object is preferably the "edges" output from the function "build3DNET"
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
  \item{weight}{
	A numeric vector of weight/ cost for each edge for the least-cost routing
}
  \item{mode}{
The direction of routing: "all" is undirected, "in" is toward the origin, and "out" is from the origin
}
}
\details{
This function extracts the edges of shortest path for each OD pair.
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
