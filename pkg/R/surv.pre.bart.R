
## you call this function before surv.bart()
## this function takes traditional time/delta
## survival variables and regressors (if any)
## and it constructs the corresponding
## X.train, y.train and X.test approriate for use with bart()

surv.pre.bart <- function(
                      times,
                      ## vector of survival times
                      
                      delta,
                      ## vector of event indicators: 1 event, 0 censoring
    
                      x.train=NULL,
                      ## matrix of covariate regressors
                      ## can be NULL, i.e. KM analog
                      
                      x.test=NULL
                      ## matrix of covariate regressors at X.test settings
                      ) {
    ## currently does not handle time dependent Xs
    ## can be extended later
    ## most likely via the alternative counting process notation

    N <- length(times)

    events <- unique(sort(times))
    ## time grid of events including censoring times

    K <- length(events)

    y.train <- integer(N) ## y.train is at least N long
            
    k <- 1
        
    for(i in 1:N) for(j in 1:K) if(events[j] <= times[i]) {
        y.train[k] <- delta[i]*(times[i] == events[j])

        k <- k+1
    }

    m <- length(y.train)

    if(length(x.train)==0) {
        p <- 0
        n <- 1
        
        X.train <- matrix(nrow=m, ncol=1, dimnames=list(NULL, 't'))
    } else {
        p <- ncol(x.train)
        
        if(length(x.test)>0) n <- nrow(x.test)
        
        X.train <- matrix(nrow=m, ncol=p+1)

        if(length(dimnames(x.train)[[2]])>0)
            dimnames(X.train)[[2]] <- c('t', dimnames(x.train)[[2]])
        else dimnames(X.train)[[2]] <- c('t', paste0('x', 1:p))
    }
    
    k <- 1
    
    for(i in 1:N) for(j in 1:K) if(events[j] <= times[i]) {
        if(p==0) X.train[k, ] <- c(events[j])
        else X.train[k, ] <- c(events[j], x.train[i, ])

        k <- k+1
    }

    if(p==0 | length(x.test)>0) {
        X.test <- matrix(nrow=K*n, ncol=p+1, dimnames=dimnames(X.train))

        for(i in 1:n) for(j in 1:K) {
            if(p==0) X.test[j, ] <- c(events[j])
            else X.test[(i-1)*K+j, ] <- c(events[j], x.test[i, ])
        }
    }
    else X.test <- matrix(nrow=0, ncol=0)*0
    
    return(list(y.train=y.train, X.train=X.train, X.test=X.test, times=events, K=K))
}
