\name{recur.bart}
\title{Nonparametric survival analysis with BART}
\alias{recur.bart}
\alias{mc.recur.bart}
\description{
Here we have implemented a simple and direct approach to utilize BART in
survival analysis that is very flexible, and is akin to discrete-time
survival analysis.  Following the capabilities of BART, we allow for
maximum flexibility in modeling the dependence of survival times on
covariates.  In particular, we do not impose proportional hazards.

To elaborate, consider data in the usual form:
\eqn{(t_i, \delta_i, {x}_i)}{(t, delta, x)} where \eqn{t_i}{t} is the event time,
\eqn{\delta_i}{delta} is an indicator distinguishing events
(\eqn{\delta=1}{delta=1}) from right-censoring
(\eqn{\delta=0}{delta=0}), \eqn{{x}_i}{x} is a vector of covariates, and
\eqn{i=1, ..., N}{i=1, ..., N (i suppressed for convenience)}
indexes subjects.

We denote the \eqn{K} distinct event/censoring times by
\eqn{0<t_{(1)}<...<t_{(K)}<\infty}{0<t(1)<...< t(K)<infinity} thus
taking \eqn{t_{(j)}}{t(j)} to be the \eqn{j^{th}}{j'th} order
statistic among distinct observation times and, for convenience,
\eqn{t_{(0)}=0}{t(0)=0}. Now consider event indicators \eqn{y_{ij}}{y(j)}
for each subject \eqn{i} at each distinct time \eqn{t_{(j)}}{t(j)}
up to and including the subject's observation time
\eqn{t_i=t_{(n_i)}}{t=t(n)} with
\eqn{n_i=\sum_j I[t_{(j)}\leq t_i]}{n=sum I[t(j)<=t]}.
This means \eqn{y_{ij}=0}{y(j)=0} if \eqn{j<n_i}{j<n} and
\eqn{y_{in_i}=\delta_i}{y(n)=delta}.

We then denote by \eqn{p_{ij}}{p(j)} the probability
of an event at time \eqn{t_{(j)}}{t(j)} conditional on no previous event. We
now write the model for \eqn{y_{ij}}{y(j)} as a nonparametric probit
regression of \eqn{y_{ij}}{y(j)} on the time \eqn{t_{(j)}}{t(j)} and the covariates
\eqn{{x}_i}{x}, and then utilize BART for binary responses.  Specifically,
\eqn{ y_{ij}\ =\ \delta_i I[t_i=t_{(j)}],\ j=1, ..., n_i }{ y(j) = delta I[t=t(j)], j=1, ..., n }; we have
\eqn{p_{ij} = F(\mu_{ij}),\ \mu_{ij} = \mu_0+f(t_{(j)}, {x}_i)}{p(j) =
  F(mu(j)), mu(j) = mu0+f(t(j), x)} where \eqn{F} denotes the standard normal cdf (probit link).
As in the binary
response case, \eqn{f} is the sum of many tree models.  }

\usage{

recur.bart(x.train, y.train=NULL, times=NULL, delta=NULL, x.test=matrix(0.0,0,0),
           k=2.0, power=2.0, base=.95, binaryOffset=NULL,
           ntree=50, ndpost=10000, nskip=250, printevery=100, keepevery=10,
           keeptrainfits=TRUE, usequants=FALSE, numcut=100, printcutoffs=0,
           verbose=TRUE,
           seed=99,    ## mc.recur.bart only
           mc.cores=2, ## mc.recur.bart only
           nice=19L    ## mc.recur.bart only
         )

mc.recur.bart(x.train, y.train=NULL, times=NULL, delta=NULL, x.test=matrix(0.0,0,0),
              k=2.0, power=2.0, base=.95, binaryOffset=NULL,
              ntree=50, ndpost=10000, nskip=250, printevery=100, keepevery=10,
              keeptrainfits=TRUE, usequants=FALSE, numcut=100, printcutoffs=0,
              verbose=TRUE,
              seed=99,    ## mc.recur.bart only
              mc.cores=2, ## mc.recur.bart only
              nice=19L    ## mc.recur.bart only
            )
}

\arguments{

  \item{x.train}{ Explanatory variables for training (in sample)
data.\cr Must be a matrix with (as usual) rows corresponding to
observations and columns to variables.\cr \code{recur.bart} will generate
draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.train (note
that the definition of \code{x.train} is dependent on whether
\code{y.train} has been specified; see below).  }

   \item{y.train}{
   Binary response dependent variable for training (in sample) data.\cr
   If \code{y.train} is \code{NULL}, then \code{y.train} (\code{x.train} and
   \code{x.test}, if specified) are generated
   by a call to \code{recur.pre.bart} (which require that \code{times} and \code{delta} be
   provided: see below); otherwise, \code{y.train} (\code{x.train} and
   \code{x.test}, if specified)
   are utilized as given assuming that the data construction has already been performed.
   }

   \item{times}{
     The time of event or right-censoring.\cr
     If \code{y.train} is \code{NULL}, then \code{times} (and \code{delta})
     must be provided.
   }

   \item{delta}{
     The event indicator: 1 is an event while 0 is censored.\cr
     If \code{y.train} is \code{NULL}, then \code{delta} (and \code{times})
     must be provided.
   }

   \item{x.test}{
   Explanatory variables for test (out of sample) data.\cr
   Must be a matrix and have the same structure as x.train.\cr
   \code{recur.bart} will generate draws of \eqn{f(t, x)} for each \eqn{x} which is a row of x.test.
   }

   \item{k}{
   k is the number of prior standard deviations \eqn{f(t, x)} is away from +/-3.
   The bigger k is, the more conservative the fitting will be.
   }

   \item{power}{
   Power parameter for tree prior.
   }

   \item{base}{
   Base parameter for tree prior.
   }

   \item{binaryOffset}{
     The model is \eqn{P(Y=1 | t, x) = F(f(t, x) + \mu_0)}{P(Y=1 | t, x) = F(f(t, x) + mu0)}
     where \eqn{\mu_0}{mu0} is specified by \code{binaryOffset}.\cr
   The idea is that \eqn{f} is shrunk towards 0, so the offset allows you to shrink towards
   a probability other than .5.\cr
   If \code{binaryOffset=NULL} when \code{times} and \code{delta} were provided, then
   an exponential distribution offset is assumed independent of the covariates, i.e.
   \code{binaryOffset=qnorm(1-exp(-mean.diff*sum(delta)/sum(times)))} where
   \code{mean.diff} is the mean of the differences of the distinct
   ordered adjacent \code{times},
   i.e. \eqn{mean(t_{(1)}-t_{(0)}, ..., t_{(K)}-t_{(K-1)})}{mean(t(1)-t(0), ..., t(K)-t(K-1))}.\cr
   If \code{binaryOffset=NULL} when \code{times} and \code{delta} were not provided, then
   \code{binaryOffset=0}.
   }

   \item{ntree}{
   The number of trees in the sum.
   }

   \item{ndpost}{
   The number of posterior draws after burn in, ndpost/keepevery will actually be returned.
   }

   \item{nskip}{
   Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
   As the MCMC runs, a message is printed every printevery draws.
   }

   \item{keepevery}{
   Every keepevery draw is kept to be returned to the user.\cr
   A \dQuote{draw} will consist of values \eqn{f^*(t, x)}{f*(t, x)}
   at \eqn{x} = rows from the train(optionally) and test data, where \eqn{f^*}{f*} denotes
   the current draw of \eqn{f}.
   }

   \item{keeptrainfits}{
   If true the draws of \eqn{f(t, x)} for \eqn{x} = rows of x.train are returned.
   }

   \item{usequants}{
   Decision rules in the tree are of the form
   \eqn{x \le c}{x <= c}  vs. \eqn{x > c}
   for each variable corresponding to a column of x.train.
   usequants determines how the set of possible c is determined.
   If usequants is true, then the c are
   a subset of the values (xs[i]+xs[i+1])/2 where xs is
   unique sorted values obtained from
   the corresponding column of x.train.
   If usequants is false, the cutoffs are equally spaced across the
   range of values taken on by the corresponding
   column of x.train.
   }

   \item{numcut}{
   The number of possible values of c (see usequants).
   If a single number if given, this is used for all variables.
   Otherwise a vector with length equal to ncol(x.train) is required,
   where the \eqn{i^{th}}{i^th} element gives the number of c used for
   the \eqn{i^{th}}{i^th} variable in x.train.
   If usequants is false, numcut equally spaced cutoffs
   are used covering the range of values in the corresponding
   column of x.train.  If usequants is true, then  min(numcut, the number of unique values in the
   corresponding columns of x.train - 1) c values are used.
   }

   \item{printcutoffs}{
   The number of cutoff rules c to printed to screen before the MCMC is run.
   Give a single integer, the same value will be used for all variables.
   If 0, nothing is printed.
   }

   \item{verbose}{
   Logical, if FALSE supress printing.
   }

   \item{seed}{
     \code{mc.recur.bart} only: seed required for reproducible MCMC.
   }

   \item{mc.cores}{
     \code{mc.recur.bart} only: number of cores to employ in parallel.
   }

   \item{nice}{
     \code{mc.recur.bart} only: set the job priority.  The default
     priority is 19: priorities go from 0 (highest) to 19 (lowest).
   }
}

\value{
   \code{recur.bart} returns a list.
Besides the items listed below,
the list has a \code{binaryOffset} component giving the value used,
 a \code{times} component giving the unique times and \code{K} which is the number of
 unique times.

   \item{yhat.train}{
   A matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row corresponds to a draw \eqn{f^*}{f*} from the posterior of \eqn{f}
   and each column corresponds to a row of x.train.
   The \eqn{(i,j)} value is \eqn{f^*(t, x)}{f*(t, x)} for the \eqn{i^{th}}{i\^th} kept draw of \eqn{f}
   and the \eqn{j^{th}}{j\^th} row of x.train.\cr
   Burn-in is dropped.
   }

   \item{haz.train}{The hazard function, \eqn{h(t|x)}, where x's are the
   rows of the training data.}
   \item{cum.train}{The cumulative hazard function, \eqn{h(t|x)}, where x's are the
   rows of the training data.}
   \item{yhat.test}{Same as yhat.train but now the x's are the rows of the test data.}
   \item{haz.test}{The hazard function, \eqn{h(t|x)}, where x's are the
   rows of the test data.}
   \item{cum.test}{The cumulative hazard function, \eqn{h(t|x)}, where x's are the
   rows of the test data.}
   \item{varcount}{a matrix with (ndpost/keepevery) rows and nrow(x.train) columns.
   Each row is for a draw. For each variable (corresponding to the columns),
   the total count of the number of times
   that variable is used in a tree decision rule (over all trees) is given.}

   Note that yhat.train and yhat.test are
   \eqn{f(t, x)} + \code{binaryOffset}.  If you want draws of the probability
   \eqn{P(Y=1 | t, x)} you need to apply the normal cdf (\code{pnorm})
   to these values.
}
\author{
Rodney Sparapani: \email{rsparapa@mcw.edu}
}
\seealso{
\code{\link{recur.pre.bart}}
}
\examples{
\dontrun{

require(survival)

subset <- -which(bladder1$stop==0)

bladder0 <- bladder1[subset, ]

id <- unique(sort(bladder0$id))

N <- length(id)
K <- max(bladder0$enum)

times <- matrix(0, nrow=N, ncol=K)
dimnames(times)[[1]] <- paste0(id)

delta <- matrix(0, nrow=N, ncol=K)
dimnames(delta)[[1]] <- paste0(id)

X.train <- matrix(0, nrow=N, ncol=2*K+1) ## add time-dependent cols too
dimnames(X.train)[[1]] <- paste0(id)
dimnames(X.train)[[2]] <- c('tx', rep(c('number', 'size'), K))
##dimnames(X.train)[[2]] <- c('Pl=1,B6=2,Th=3', rep(c('number', 'size'), K))

for(i in 1:N) {
    h <- id[i]

    for(j in 1:K) {
        k <- which(bladder0$id==h & bladder0$enum==j)

        if(length(k)==1) {
            times[i, j] <- bladder0$stop[k]
            delta[i, j] <- (bladder0$status[k]==1)*1

            if(j==1) {
                X.train[i, 1] <- as.numeric(bladder0$treatment[k])
                X.train[i, 2] <- bladder0$number[k]
                X.train[i, 3] <- bladder0$size[k]
            }
            else if(delta[i, j]==1) {
                suppressWarnings(X.train[i, 2*j] <- as.numeric(bladder0$rtumor[k]))
                suppressWarnings(X.train[i, 2*j+1] <- as.numeric(bladder0$rsize[k]))
            }
        }
    }
}

pre <- recur.pre.bart(times=times, delta=delta, x.train=X.train)

K <- nrow(pre$X.train)

for(j in 1:K) {
    if(pre$X.train[j, 3]>0) {
        pre$X.train[j, 5] <- pre$X.train[j, 5+pre$X.train[j, 3]*2]
        pre$X.train[j, 6] <- pre$X.train[j, 6+pre$X.train[j, 3]*2]
    }
}

pre$X.train <- pre$X.train[ , 1:6]

pre$X.test <- rbind(pre$X.train, pre$X.train, pre$X.train)

pre$X.test[ , 4] <- c(rep(1, K), rep(2, K), rep(3, K))

post <- recur.bart(y.train=pre$y.train, x.train=pre$X.train, x.test=pre$X.test)

}
}
\keyword{nonparametric recurrent events model}
\keyword{nonproportional hazards}
