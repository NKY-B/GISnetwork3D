\name{tobr_hf_pace}
\alias{tobr_hf_pace}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the walking speed according to the Tobler's hiking function.
}
\description{
Calculate the walking speed according to the Tobler's hiking function.
}
\usage{
tobr_hf_pace(Deg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Deg}{
A numeric vector in degree slope.
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
## The function is currently defined as
function (Deg)
{
    S = tan(Deg * pi/180)
    result = 0.6 * exp(3.5 * abs(S + 0.05))
    return(result)
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
