\name{build3DNET}
\alias{build3DNET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build network for 3D line vector.
}
\description{
Build network for 3D line vector.
}
\usage{
build3DNET(line, ped = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{line}{
 A sf line feature with the z value.
}
  \item{ped}{
Boolean: Indicate whether the network is built as a pedestrian network. If TRUE, this function will create a network that supports bi-directional routing with slope/directional dependent walking time(s) for each edge.
}
}
\details{
This function produces the edges and vertices of the input line feature.

Each edge is only connected by two vertices as this function first explodes the line feature. The edges record the from and to vertices, of which the index is the vID recorded in the vertices layer. The digitization direction determines the from and to vertices. A unique ID is generated for each edge segment (eID) and each vertex (vID).

If the ped is TRUE, the produced edges will record the slopes and the walking time for both directions. The walking time is based on Tobler's hiking function.

If a column "directCond" is presented in the attribute table of the input line, directed edges will be filtered before producing a bi-directional network when the ped argument is TRUE. Only 1 and 0 are accepted for the directCond. 1 indicates that the edge is directed according to its digitization direction. 0 indicates that it is bi-directed, and bi-directed walking speed and travel time were calculated for these edges.
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
