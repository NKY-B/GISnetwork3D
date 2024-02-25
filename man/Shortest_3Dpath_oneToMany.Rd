\name{Shortest_3Dpath_oneToMany}
\alias{Shortest_3Dpath_oneToMany}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create paths from one origin to one or more destination(s).
}
\description{
Create paths from one origin to one or more destination(s).
}
\usage{
Shortest_3Dpath_oneToMany(graph, edges, weight, oV, dV, mode = "all")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
An igraph object created from the function "GISigraph3D"
}
  \item{edges}{
A sf line layer with the fields of "from" and "to" to indicate the IDs of vertices for the start and end of each line segment. This object is preferably the "edges" output from the function "build3DNET"
}
  \item{weight}{
A numeric vector of weight/ cost for each edge for the least-cost routing
}
  \item{oV}{
A numeric vector of the vertex ID of the origin location
}
  \item{dV}{
A numeric vector of the vertex/vertices ID(s) of the destination location(s)
}
  \item{mode}{
The direction of routing: "all" is undirected, "in" is toward the origin, and "out" is from the origin
}
}
\details{
This function create separate shortest path for each O and D pair in a one-to-many basis.
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
