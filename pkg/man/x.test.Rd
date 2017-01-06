\name{x.test}
\alias{x.test}
\docType{data}
\title{
A data set used in example of \code{recur.bart}. 
}
\description{
  
A matrix containing the testing set for a real data example of
recurrent events survival analysis.  There are 488 patients in the
cohort: 235 in the training set and 253 in the testing set.  See the
Reference below (and the References therein) for more detailed
information; a brief synopsis follows.

\code{x.test} contains both the training set and the testing set.  There
are 798 unique time points so there are 235*798=187530 rows of the
training set followed by 253*798=201894 rows of the testing set.  For
patient's who died prior to the end of follow-up, their external factors
are last value carried forward.  Therefore, we can use \code{x.test}
to estimate the cumulative hazard for all patients for all time points.

}
\usage{data(x.test)}
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
data(x.test)
head(x.test[ , 1:10])
}
\keyword{datasets}
