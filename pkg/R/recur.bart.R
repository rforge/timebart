## run BART with recurrent events

recur.bart <- function(
    x.train, y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L),
    x.test.nogrid = FALSE, ## you may not need the whole grid
    keepcall = FALSE, ## the call object can get rather large
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = NULL,
    ntree = 50L,
    ndpost = 10000L, nskip = 250L,
    printevery = 100L,
    keepevery = 10L, keeptrainfits = TRUE,
    usequants = FALSE, numcut = 100L, printcutoffs = 0L,
    verbose = TRUE,
    seed = 99L,    ## only used by mc.recur.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
    )
{
    if(length(y.train)==0) {
        recur <- recur.pre.bart(times, delta, x.train, x.test)

        y.train <- recur$y.train
        x.train <- recur$X.train
        x.test  <- recur$X.test

        times   <- recur$times
        K       <- recur$K

        if(length(binaryOffset)==0) {
            lambda <- sum(delta)/sum(times)
            delta  <- mean(times[2:K]-times[1:(K-1)])
            binaryOffset <- qnorm(1-exp(-lambda*delta))
        }
    }
    else {
        if(length(binaryOffset)==0) binaryOffset <- 0

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    ##cat('dbarts\n')

    post <- bart(x.train=x.train, y.train=y.train, x.test=x.test,
                        sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
                        k=k, keepcall=keepcall,
                        power=power, base=base,
                        binaryOffset=binaryOffset,
                        ntree=ntree,
                        ndpost=ndpost, nskip=nskip,
                        printevery=printevery, keepevery=keepevery, keeptrainfits=keeptrainfits,
                        usequants=usequants, numcut=numcut, printcutoffs=printcutoffs,
                        verbose=verbose)

    ##post$call <- NULL
    post$binaryOffset <- NULL
    post$times <- times
    post$K <- K
    post$x.train <- x.train

    ## if(keepevery>1L) {
    ##     thin <- seq(1, ndpost, keepevery)
    ##     post$yhat.train <- post$yhat.train[thin, ]
    ##     post$varcount <- post$varcount[thin, ]
    ## }

    post$cum.train <- pnorm(post$yhat.train)
    post$haz.train <- post$cum.train

    H <- nrow(x.train)

    for(h in 1:H) {
        j <- which(x.train[h, 1]==times) ## for grid points only

        if(j==1) post$haz.train[ , h] <- post$haz.train[ , h]/times[1]
        else {
            post$haz.train[ , h] <- post$haz.train[ , h]/(times[j]-times[j-1])
            post$cum.train[ , h] <- post$cum.train[ , h-1]+post$cum.train[ , h]
        }
    }

    ## H <- nrow(x.train)
    ## time <- 0

    ## for(h in 1:H) {
    ##     prev <- time
    ##     time <- x.train[h, 1]

    ##     if(time==post$times[1]) post$haz.train[ , h] <- post$haz.train[ , h]/time
    ##     else {
    ##         post$haz.train[ , h] <- post$haz.train[ , h]/(time-prev)
    ##         post$cum.train[ , h] <- post$cum.train[ , h-1]+post$cum.train[ , h]
    ##     }
    ## }

    if(length(x.test)>0) {
        ##if(keepevery>1L) post$yhat.test <- post$yhat.test[thin, ]

        post$x.test <- x.test
        
        post$haz.test <- pnorm(post$yhat.test)

        if(!x.test.nogrid) {
            post$cum.test <- post$haz.test

            H <- nrow(x.test)

            for(h in 1:H) {
                j <- which(x.test[h, 1]==times) ## for grid points only

                if(j==1) post$haz.test[ , h] <- post$haz.test[ , h]/times[1]
                else {
                    post$haz.test[ , h] <- post$haz.test[ , h]/(times[j]-times[j-1])
                    post$cum.test[ , h] <- post$cum.test[ , h-1]+post$cum.test[ , h]
                }
            }
        }

        ## H <- nrow(x.test)
        ## time <- 0

        ## for(h in 1:H) {
        ##     prev <- time
        ##     time <- x.test[h, 1]

        ##     if(time==post$times[1]) post$haz.test[ , h] <- post$haz.test[ , h]/time
        ##     else {
        ##         post$haz.test[ , h] <- post$haz.test[ , h]/(time-prev)
        ##         post$cum.test[ , h] <- post$cum.test[ , h-1]+post$cum.test[ , h]
        ##     }
        ## }

    }

    return(post)
}
