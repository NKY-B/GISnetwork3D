\name{calculate_slope}
\alias{calculate_slope}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the degree slope.
}
\description{
Calculate the degree slope.
}
\usage{
calculate_slope(x1, y1, z1, x2, y2, z2)
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x1, y1, z1, x2, y2, z2)
{
    run = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    rise = z2 - z1
    degree_slope = atan(rise/run) * (180/pi)
    return(degree_slope)
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
