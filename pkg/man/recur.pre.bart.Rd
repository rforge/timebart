\name{recur.pre.bart}
\title{Data construction for nonparametric survival analysis with BART}
\alias{recur.pre.bart}
\description{
Recurrent event data contained in \eqn{(t_1,\delta_1, ..., t_k,\delta_k, x)} must be translated to data
suitable for the BART model; see \code{recur.bart} for more details.
}
\usage{
recur.pre.bart( times, delta, x.train=NULL, probs=c(0.15, 0.15) )
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

   \item{probs}{ \code{X.test} is generated from \code{x.train} in a
     similar fashion to \code{X.train}.  However, observations beyond
     censoring contain \code{NA} for columns two and three containing
     \eqn{v(t)} and \eqn{N(t-)} respectively.  In order to preform
     prediction with \code{X.test}, then these missing values need to be
     handled accordingly.  We provide an automated method for this which
     is returned as \code{X.base}.  It would be difficult to provide an
     efficient method that randomly generates \eqn{v(t)} and \eqn{N(t-)}
     since they are structurally dependent.  Rather, we
     deterministically generate them in such a fashion that they define
     a "middle" baseline, i.e. the median \eqn{N(t-)} as a function of
     \eqn{t} for non-censored subjects (calculated from \code{X.train});
     the median is adjusted upward as the censoring increases with time
     (if it does at all) to a maximum value controlled by the second
     argument contained in \code{probs} since otherwise the "low" risk
     subjects would be over-represented.  Each subject's previous value
     is carried forward unless their computed sojourn time (which is a
     function of the generated \eqn{N(t-)}) is longer than the the
     median \eqn{v(t)} as a function of \eqn{t} (which is also
     calculated from \code{X.train}); the median is adjusted downward as
     the censoring increases with time (if it does at all) to a maximum
     value controlled by the first argument contained in \code{probs}. }

}

\value{
   \code{recur.pre.bart} returns a list.
Besides the items listed below, the list has
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.
   \item{y.train}{A vector of binary responses.}
   \item{X.train}{A matrix with the rows of the training data.}
   \item{X.test}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
   \item{X.base}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
   \item{sojourn}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
   \item{pattern}{Generated from \code{x.train} (see discussion above included in
     the argument \code{probs}).}
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
