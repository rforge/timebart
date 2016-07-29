## 02/20/16

## you call this function before bart()
## this function takes traditional time/delta
## recurrent event variables and regressors (if any)
## and it constructs the corresponding
## X.train, y.train and X.test approriate for use with bart()

recur.pre.bart <- function(
                      times,
                      ## matrix of recur times

                      delta,
                      ## matrix of indicators: 0=no event, 1=event

                      x.train=NULL,
                      ## matrix of covariate regressors
                      ## can be NULL, i.e. KM analog

                      x.test=NULL
                      ## matrix of covariate regressors at X.test settings
                      ## does nothing for now since there is no obvious basis for v(t) and N(t-)
                      ) {
    ## currently does not handle time dependent Xs
    ## can be extended later
    ## most likely via the alternative counting process notation

    N <- nrow(times)
    J <- ncol(times)

    times <- cbind(times, apply(times, 1, max))

    dimnames(times)[[2]][J+1] <- 'stop'

    events <- unique(sort(times))
    ## time grid of events including censoring times

    if(events[1]==0) events <- events[-1]

    K <- length(events)

    y.train <- integer(N) ## y.train is at least N long

    k <- 1

    for(i in 1:N) for(j in 1:K) if(events[j] <= times[i, J+1]) {
        y.train[k] <- 0

        for(h in 1:J) if(y.train[k]==0 & events[j]==times[i, h] & delta[i, h]==1)
            y.train[k] <- 1

        k <- k+1
    }

    m <- length(y.train)

    if(length(x.train)==0) {
        p <- 0
        n <- 1

        X.train <- matrix(nrow=m, ncol=3)

        dimnames(X.train)[[2]] <- c('t', 'v', 'N')
    } else {
        p <- ncol(x.train)

        if(length(x.test)>0) n <- nrow(x.test)

        ##print(list(m, n, p))

        X.train <- matrix(nrow=m, ncol=p+3)

        if(length(dimnames(x.train)[[2]])>0)
            dimnames(X.train)[[2]] <- c('t', 'v', 'N', dimnames(x.train)[[2]])
        else dimnames(X.train)[[2]] <- c('t', 'v', 'N', paste0('x', 1:p))
    }

    k <- 1

    for(i in 1:N) {
        n.t <- 0
        t.0 <- 0

        for(j in 1:K) if(events[j] <= times[i, J+1]) {
            X.train[k, 1:3] <- c(events[j], events[j]-t.0, n.t)

            for(h in 1:J) if(events[j]==times[i, h] & delta[i, h]==1) {
                n.t <- n.t+1
                t.0 <- events[j]
            }

            if(p>0) X.train[k, 4:(3+p)] <- x.train[i, ]

            k <- k+1
        }
    }

## automated X.test creation is not feasible since there is no obvious basis for v(t) and N(t-)
    ## if(p==0 | length(x.test)>0) {
    ##     X.test <- matrix(nrow=K*n, ncol=p+3, dimnames=dimnames(X.train))

    ##     ## we can summarize this, but this is only one of many possibilities
    ##     for(i in 1:n) for(j in 1:K) {
    ##         X.test[(i-1)*K+j, 1:(J+1)] <- c(rep(events[j], J), 0)
    ##         if(p>0) X.test[(i-1)*K+j, (J+2):(J+1+p)] <- x.test[i, ]
    ##     }
    ## }
    ## else X.test <- matrix(0.0, 0L, 0L)

    return(list(y.train=y.train, X.train=X.train, X.test=matrix(0.0, 0L, 0L), times=events, K=K))
                ##X.test=data.matrix(X.test), times=events, K=K))
}
