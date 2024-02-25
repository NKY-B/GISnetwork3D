\name{dist3D}
\alias{dist3D}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the straight-line 3D distance between two points.
}
\description{
Calculate the straight-line 3D distance between two points.
}
\usage{
dist3D(x1, y1, z1, x2, y2, z2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1}{
A numeric vector of X coordinate of the from point.
}
  \item{y1}{
A numeric vector of Y coordinate of the from point.
}
  \item{z1}{
A numeric vector of Z coordinate of the from point.
}
  \item{x2}{
A numeric vector of X coordinate of the to point.
}
  \item{y2}{
A numeric vector of Y coordinate of the to point.
}
  \item{z2}{
A numeric vector of Z coordinate of the to point.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
# Define the from and to coordinates
x1 = 1
y1 = 2
z1 = 3
x2 = 4
y2 = 9
z2 = 12

result = dist3D(x1, y1, z1, x2, y2, z2)
# result = 11.78983

## The function is currently defined as
function (x1, y1, z1, x2, y2, z2)
{
    distance = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
    return(distance)
  }
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
