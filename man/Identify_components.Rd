\name{Identify_components}
\alias{Identify_components}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify the components of the network.
}
\description{
Identify the components of the network.
}
\usage{
Identify_components(line)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{line}{
 A sf line feature with the z value.
}
}
\details{
This function identifies the component to which each edge belongs for the input line layer. This function first converts the input line layer into network data using the "build3DNET" function, and it is then converted to igraph using the "GISigraph" function to identify the component of edges. Please be cautious that the output layer is an exploded version of the line, meaning each line segment only has two connected vertices.
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
