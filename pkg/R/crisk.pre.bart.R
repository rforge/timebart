
## you call this function before crisk.bart()
## this function takes traditional time/delta
## competing risk variables and regressors (if any)
## and it constructs the corresponding
## X.train, y.train, y.train2, and X.test approriate for use with bart()

crisk.pre.bart <- function(
                      times,
                      ## vector of survival times
                      
                      delta,
                      ## vector of event indicators
                      ## 0=censoring, 1=cause of interest, 2=other causes
    
                      x.train=NULL,
                      ## matrix of covariate regressors
                      ## can be NULL, i.e. KM analog
                      
                      x.test=NULL
                      ## matrix of covariate regressors at X.test settings
                      ) {
    ## currently does not handle time dependent Xs
    ## can be extended later
    ## most likely via the alternative counting process notation

    if(!all(unique(sort(delta))==0:2)) 
        stop('delta must be coded as: 0(censored), 1(cause of interest) or 2(other cause)')
        
    pre <- surv.pre.bart(times=times, 1*(delta==1), x.train=x.train, x.test=x.test)

    pre$cond <- which(pre$y.train==0)
    
    pre2 <- surv.pre.bart(times=times, 1*(delta==2))
    
    pre$y.train2 <- pre2$y.train
    pre$binaryOffset2 <- pre2$binaryOffset

    return(pre)
}
