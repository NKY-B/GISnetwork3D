\name{reverseLine_by_condition}
\alias{reverseLine_by_condition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reverse line by condition.
}
\description{
This function reverses the direction of each line for the sf line features.
}
\usage{
reverseLine_by_condition(line, condition = -999)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{line}{
A sf line feature.
}
  \item{condition}{
A numeric vector (1 or 0) defining which line requires reversion. The input length must be the same as the total number of rows of lines. Value of 1 = requiring reversion. 0 = not requiring reversion. All line segments are reversed if the condition is not specified by default.
}
}
\details{
This function reverses the direction of each line for the sf line features. The output line layer has an additional column called reversed. 0 means no reversion, and 1 means reversed.
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
