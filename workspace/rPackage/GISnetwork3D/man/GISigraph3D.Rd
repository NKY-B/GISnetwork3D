\name{GISigraph3D}
\alias{GISigraph3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create igraph from edges and vertices produced from the function "build3DNET".
}
\description{
Create igraph from edges and vertices produced from the function "build3DNET".
}
\usage{
GISigraph3D(edges, vertices, directed = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges}{
Edges produced from the function "build3DNET"
}
  \item{vertices}{
Vertices produced from the function "build3DNET"
}
  \item{directed}{
Boolean: If TRUE, a directed graph is produced.
}
}
\details{
This function convert the "edges" and "vertices" to igraph object for network analysis. The input "edges" and vertices should be the output of "build3DNET". Or a line layer of "edges" and a point layer of "vertices" with the Z dimension. The "edges" must record the vertex ID of the from and to vertices in the first two columns. The "vertices" should have a field "vID" in the first column, which should be a vector of chronological order (1,2,3...n).
}
\value{
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
