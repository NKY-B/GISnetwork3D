\name{findVertices3D}
\alias{findVertices3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the nearest network vertex for each input point.
}
\description{
Find the nearest network vertex for each input point.
}
\usage{
findVertices3D(pt, vertices)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pt}{
A sf point feature with Z value.
}
  \item{vertices}{
The vertices of a network produced from the function - "build3DNET"
}
}
\details{
This function finds the nearest vertex from the vertices layer, produced from the function "build3DNET", for each input "pt". The input "pt" must have the z value. The output has the original pt with the vID of the nearest network vertex.
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
