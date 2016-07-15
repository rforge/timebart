\name{recur.pre.bart}
\title{Data construction for nonparametric survival analysis with BART}
\alias{recur.pre.bart}
\description{
Recurrent event data contained in \eqn{(t_1,\delta_1, ..., t_k,\delta_k, x)} must be translated to data
suitable for the BART model; see \code{recur.bart} for more details.
}
\usage{
recur.pre.bart( times, delta, x.train=NULL, x.test=NULL )
}
\arguments{
   \item{times}{
     Matrix of time to event or right-censoring.\cr
   }

   \item{delta}{
     Matrix of event indicators: 1 is an event while 0 is censored.\cr
   }

   \item{x.train}{
   Explanatory variables for training (in sample) data.\cr
   If provided, must be a matrix
   with (as usual) rows corresponding to observations and columns to variables.\cr
   }

   \item{x.test}{
     Currently, does nothing since there is not obvious
     basis for \eqn{v(t)} and \eqn{N(t-)}.\cr
   }
}

\value{
   \code{recur.pre.bart} returns a list.
Besides the items listed below, the list has
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.
   \item{y.train}{A vector of binary responses.}
   \item{X.train}{A matrix with the rows of the training data.}
   \item{X.test}{A matrix with the zero dimensions.}
}

\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{recur.bart}}
}
\examples{
\dontrun{
require(timebart)


}
}
\keyword{data construction}
