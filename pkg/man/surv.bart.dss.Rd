\name{surv.bart.dss}
\title{Decoupling Shrinkage and Selection (DSS) variable selection for nonparametric survival analysis with BART}
\alias{surv.bart.dss}
\alias{mc.surv.bart.dss}
\description{
Here we implement the Decoupling Shrinkage and Selection (DSS) method for variable selection in
nonparametric survival analysis with BART.  The method is
computationally friendly. }

\usage{
    surv.bart.dss(post,
                  mc.cores=2L,   ## used by mc.surv.bart.dss only
                  nice=19L       ## used by mc.surv.bart.dss only
                  )
    mc.surv.bart.dss(post,
                     mc.cores=2L,## used by mc.surv.bart.dss only
                     nice=19L    ## used by mc.surv.bart.dss only
                     )
}

\arguments{

    \item{post}{ The posterior from a previous call to \code{surv.bart} or \code{mc.surv.bart}. }

   \item{mc.cores}{
     Number of cores to employ in parallel: used by \code{mc.surv.bart.dss} only.
   }

   \item{nice}{
     Set the job priority: used by \code{mc.surv.bart.dss} only.  The default
     priority is 19: priorities go from 0 (highest) to 19 (lowest).
   }
}

\value{
   Returns a list of selected variable column numbers and names.
}
\references{

A General Approach to Variable Selection for Bayesian Nonparametric Models
JSM, August 10, 2015. \url{http://rob-mcculloch.org/nonpar-varsel_2015-8-10.pdf}
}
\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{surv.bart}}
}
\examples{
\dontrun{
require(timebart)
}
}
\keyword{nonproportional hazards variable selection}
