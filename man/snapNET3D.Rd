\name{snapNET3D}
\alias{snapNET3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Connect disconnected line within a distance threshod.
}
\description{
Connect disconnected line within a distance threshod.
}
\usage{
snapNET3D(line, threshold = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{line}{
 A sf line feature with the z value.
}
  \item{threshold}{
A distance threshold within which the disconnected vertices are to be linked by a small edge.
}
}
\details{
This function connects the disconnected lines that are separated from each other within a predefined distance in 3D. Please note that the output layer is an exploded version of the line. Each line segment only has two connected vertices.
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
