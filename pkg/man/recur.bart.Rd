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
           x.test.short=FALSE, keepcall=FALSE,
           k=2.0, power=2.0, base=.95, binaryOffset=NULL,
           ntree=50, ndpost=10000, nskip=250, printevery=100, keepevery=10,
           keeptrainfits=TRUE, usequants=FALSE, numcut=100, printcutoffs=0,
           verbose=TRUE,
           seed=99,    ## mc.recur.bart only
           mc.cores=2, ## mc.recur.bart only
           nice=19L    ## mc.recur.bart only
         )

mc.recur.bart(x.train, y.train=NULL, times=NULL, delta=NULL, x.test=matrix(0.0,0,0),
              x.test.short=FALSE, keepcall=FALSE,
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

   \item{x.test.short}{
     Occasionally, you do not need the entire time grid for
     \code{x.test}.  If so, then for performance reasons, you can set this argument to \code{TRUE}.
   }

   \item{keepcall}{
     Logical; if \code{FALSE}, returned object will have \code{call} set to
          \code{call("NULL")}, otherwise the call used to instantiate BART.
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
\references{
Xu, Daniels and Winterstein.  Sequential BART for imputation of missing
covariates.  Biostatistics 2016 doi: 10.1093/biostatistics/kxw009
with supplement
\url{http://biostatistics.oxfordjournals.org/content/early/2016/03/15/biostatistics.kxw009/suppl/DC1}
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

X.train <- matrix(NA, nrow=N, ncol=2*K+1) ## add time-dependent cols too
dimnames(X.train)[[1]] <- paste0(id)
dimnames(X.train)[[2]] <- c('tx', rep(c('number', 'size'), K))

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
                if(bladder0$rtumor[k]!='.')
                    X.train[i, 2*j] <- as.numeric(bladder0$rtumor[k])
                if(bladder0$rsize[k]!='.')
                    X.train[i, 2*j+1] <- as.numeric(bladder0$rsize[k])
            }
        }
    }
}

pre <- recur.pre.bart(times=times, delta=delta, x.train=X.train)

J <- nrow(pre$X.train)

for(j in 1:J) {
    if(pre$X.train[j, 3]>0) {
        pre$X.train[j, 5] <- pre$X.train[j, 5+pre$X.train[j, 3]*2]
        pre$X.train[j, 6] <- pre$X.train[j, 6+pre$X.train[j, 3]*2]
    }
}

pre$X.train <- pre$X.train[ , 1:6]

J <- nrow(pre$X.base)

for(j in 1:J) {
    if(pre$X.base[j, 3]>0) {
        pre$X.base[j, 5] <- pre$X.base[j, 5+pre$X.base[j, 3]*2]
        pre$X.base[j, 6] <- pre$X.base[j, 6+pre$X.base[j, 3]*2]
    }
}

pre$X.base <- pre$X.base[ , 1:6]

pre$X.base <- rbind(pre$X.base, pre$X.base, pre$X.base)

pre$X.base[ , 4] <- c(rep(1, J), rep(2, J), rep(3, J))

plot(pre$times, pre$pattern, ylim=c(0, max(pre$X.train[ , 3])))
set.seed(10)
u <- runif(N)*0.25
for(i in 1:N) {
    k <- (i-1)*60
    h <- k+1:60
    j <- which(!is.na(pre$X.test[h, 3]))
    points(pre$times[j], pre$X.test[k+j, 3]+u[i], type='l')
    j <- c(max(j), which(is.na(pre$X.test[h, 3])))
    points(pre$times[j], pre$X.base[k+j, 3]+u[i], type='l', col='red')
}

## since BART is essentially a nonparametric function, handling
## missing data can, in certain circumstances, be as easy as creating
## a missing "category" for continuous and ordinal variables
## in the bladder ex, both number and size are recorded as integers
## from 1 to 8 so a reasonable missing category might be 0 or 9
## however, we really want to create a generous amount of
## separation between observed and missing data to give BART a
## reasonable chance of splitting between them; but, perhaps, not too
## generous otherwise you may sacrfice performance within the
## observed data range; the usual caveats apply and it is always
## a good idea to explore more sophisticated methods (see below)    
pre$X.base[which(is.na(pre$X.base[ , 5])), 5] <- -2
pre$X.base[which(is.na(pre$X.base[ , 6])), 6] <- -2
    
## or impute missing data via your favorite method
## such as Sequential BART (see ref above)
## library(MASS)
## library(Rlab)
## library(msm)
## source('~/sequential_bart/load-pbart.R')
## set.seed(21)
## base <- BartMI(bartmpicstem, xx=pre$X.base, yy=NULL, datatype=rep(0, 6),
##                type=0, numskip=20, burn=1000)
## coarsen the imputed data same way as observed example data
## base$imputed5[which(base$imputed5[ , 5]<1), 5] <- 1
## base$imputed5[which(base$imputed5[ , 5]>8), 5] <- 8
## base$imputed5[ , 5] <- round(base$imputed5[ , 5])
## base$imputed5[which(base$imputed5[ , 6]<1), 6] <- 1
## base$imputed5[which(base$imputed5[ , 6]>8), 6] <- 8
## base$imputed5[ , 6] <- round(base$imputed5[ , 6])
## pre$X.base <- base$imputed5
## dimnames(pre$X.base)[[2]] <- c('t', 'v', 'N', 'tx', 'number', 'size')

set.seed(99)
post <- recur.bart(y.train=pre$y.train, x.train=pre$X.train, x.test=pre$X.base)
## depending on your performance, you may want to run in parallel if available
## post <- mc.recur.bart(y.train=pre$y.train, x.train=pre$X.train, x.test=pre$X.base,
##                       mc.cores=16, seed=99)

K <- post$K
M <- 1000
haz.Pl <- matrix(0, nrow=M, ncol=K)
haz.B6 <- matrix(0, nrow=M, ncol=K)
haz.Th <- matrix(0, nrow=M, ncol=K)

for(i in 1:N) {
    j <- (i-1)*K+1:K
    haz.Pl <- haz.Pl+post$haz.test[1:M, j]/N
    haz.B6 <- haz.B6+post$haz.test[1:M, j+N*K]/N
    haz.Th <- haz.Th+post$haz.test[1:M, j+N*K*2]/N
}

haz.Pl.m <- apply(haz.Pl, 2, mean)
haz.B6.m <- apply(haz.B6, 2, mean)
haz.Th.m <- apply(haz.Th, 2, mean)

plot(post$times, haz.Pl.m, type='l',
     main='bladder ex: Placebo (black), B6 (red), Thiotepa (blue)',
     ylab=expression(lambda(t)), xlab='t (months)',
     ylim=c(0, max(c(haz.Pl.m, haz.B6.m, haz.Th.m))))
points(post$times, haz.B6.m, type='l', col='red')
points(post$times, haz.Th.m, type='l', col='blue')

plot(post$times, apply(1*((haz.Th/haz.Pl) < 1), 2, mean),
     type='l', ylim=c(0, 1),
     main='bladder ex: Thiotepa vs. Placebo',
     ylab=expression(P(RR(t)<1)), xlab='t (months)')
     points(c(0, max(post$times)), c(0.5, 0.5), type='l', col='gray')

rr.Th.Pl.m <- apply(haz.Th/haz.Pl, 2, mean)
rr.Th.Pl.l <- apply(haz.Th/haz.Pl, 2, quantile, probs=0.25)
rr.Th.Pl.u <- apply(haz.Th/haz.Pl, 2, quantile, probs=0.75)
plot(post$times, rr.Th.Pl.m, type='l', log='y',
     main='bladder ex: Thiotepa vs. Placebo with 50% CI',
     ylab=expression(RR(t)), xlab='t (months)', ylim=c(0.625, 1.6))
points(post$times, rr.Th.Pl.l, type='l', col='gray')
points(post$times, rr.Th.Pl.u, type='l', col='gray')
points(c(0, max(post$times)), c(1, 1), type='l')

rr.Th.B6.m <- apply(haz.Th/haz.B6, 2, mean)
rr.Th.B6.l <- apply(haz.Th/haz.B6, 2, quantile, probs=0.25)
rr.Th.B6.u <- apply(haz.Th/haz.B6, 2, quantile, probs=0.75)
plot(post$times, rr.Th.B6.m, type='l', log='y',
     main='bladder ex: Thiotepa vs. Vitamin B6 with 50% CI',
     ylab=expression(RR(t)), xlab='t (months)', ylim=c(0.625, 1.6))
points(post$times, rr.Th.B6.l, type='l', col='gray')
points(post$times, rr.Th.B6.u, type='l', col='gray')
points(c(0, max(post$times)), c(1, 1), type='l')

plot(post$times, apply(1*((haz.Th/haz.B6) < 1), 2, mean),
     type='l', ylim=c(0, 1),
     main='bladder ex: Thiotepa vs. Vitamin B6',
     ylab=expression(P(RR(t)<1)), xlab='t (months)')
     points(c(0, max(post$times)), c(0.5, 0.5), type='l', col='gray')

}
}
\keyword{nonparametric recurrent events model}
\keyword{nonproportional hazards}
