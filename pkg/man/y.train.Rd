\name{y.train}
\alias{y.train}
\alias{y.test}
\docType{data}
\title{
A data set used in example of \code{recur.bart}. 
}
\description{
  
Two vectors containing the training and testing set outcomes for a real
data example of recurrent events survival analysis.  There are 488
patients in the cohort: 235 in the training set and 253 in the testing
set.  See the Reference below (and the References therein) for more
detailed information; a brief synopsis follows.

\code{y.train} contains the training set only.  \code{y.test} contains
both the training set and the testing set.  There are 798 unique time
points so there are 235*798=187530 rows of the training set followed by
253*798=201894 rows of the testing set.  For patient's who died prior to
the end of follow-up, no event, i.e., zero, is carried forward.

}
\usage{
data(y.train)
data(y.test)
}
\seealso{x.train}
\references{
Sparapani R., Rein L., Tarima S., Jackson T., Meurer J. (2017)
Nonparametric recurrent events analysis and
Bayesian Additive Regression Trees with an application
to hospital admissions for patients with diabetes.
\emph{Biostatistics} \bold{manuscript submitted}.

%% Xu D, Daniels MJ and Winterstein AG. (2016).
%% Sequential BART for imputation of missing covariates.
%% \emph{Biostatistics} \bold{17}, 589-602.
}
\examples{
data(y.train)
data(y.test)
table(y.train)
table(y.test)
}
\keyword{datasets}
