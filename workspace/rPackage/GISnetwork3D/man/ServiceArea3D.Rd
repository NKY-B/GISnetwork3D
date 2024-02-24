\name{ServiceArea3D}
\alias{ServiceArea3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create 3D service area.
}
\description{
Create 3D service area.
}
\usage{
ServiceArea3D(graph, edges, vertices, oID, oV, mode = "all", weight, search = 500, costThreshold = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
An igraph object created from the function "GISigraph3D"
}
  \item{edges}{
A sf line layer with the fields of "from" and "to" to indicate the IDs of vertices for the start and end of each line segment. This object is preferably the "edges" output from the function "build3DNET"
}
  \item{vertices}{
A sf point layer with the fields of "vID".
}
  \item{oID}{
A numeric vector of the IDs of the origin locations
}
  \item{oV}{
A numeric vector of the vertices IDs of the origin locations
}
  \item{mode}{
The direction of routing: "all" is undirected, "in" is toward the origin, and "out" is from the origin
}
  \item{weight}{
A numeric vector of weight/ cost for each edge for the least-cost routing
}
  \item{search}{
A number indicating the size of straight-line radius used to first identify possible reachable vertices. This is expressed as 2D distance
}
  \item{costThreshold}{
The catchment size expressed as the maximum 3D network cost.
}
}
\details{
This function identifies vertices reachable within a cost threshold from an origin vertex and also all edges linking to these vertices. Therefore, the output is a list of (1) sf point showing the reachable vertices, and (2) the edges within threshold cost from each origin location.
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
