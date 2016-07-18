## run BART and generate survival
## 07/13/16

surv.bart <- function(
    x.train, y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L),
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = NULL,
    ntree = 50L,
    ndpost = 10000L, nskip = 250L,
    printevery = 100L,
    keepevery = 10L, keeptrainfits = TRUE,
    usequants = FALSE, numcut = 100L, printcutoffs = 0L,
    verbose = TRUE,
    id = NULL,     ## only used by surv.bart
    seed = 99L,    ## only used by mc.surv.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
)
{
    if(length(y.train)==0) {
        surv <- surv.pre.bart(times, delta, x.train, x.test)

        y.train <- surv$y.train
        x.train <- surv$X.train
        x.test  <- surv$X.test

        times   <- surv$times
        K       <- surv$K

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

    cat('timebart::surv.bart\n')
    
    post <- dbarts::bart(x.train=x.train, y.train=y.train, x.test=x.test,
                        k=k,
                        power=power, base=base,
                        binaryOffset=binaryOffset,
                        ntree=ntree,
                        ndpost=ndpost, nskip=nskip,
                        printevery=printevery, keepevery=1L, keeptrainfits=keeptrainfits,
                        usequants=usequants, numcut=numcut, printcutoffs=printcutoffs,
                        verbose=verbose)

    post$call <- NULL
    post$binaryOffset <- NULL
    post$id <- id
    post$times <- times
    post$K <- K
    post$x.train <- x.train

    if(keepevery>1L) { ## thinning with dbarts not available
        thin <- seq(1, ndpost, keepevery)
        post$yhat.train <- post$yhat.train[thin, ]
        post$varcount <- post$varcount[thin, ]
    }

    post$surv.train <- 1-pnorm(post$yhat.train)

    H <- nrow(x.train)/K ## the number of different settings

    for(h in 1:H) for(j in 2:K)
                      post$surv.train[ , K*(h-1)+j] <-
                          post$surv.train[ , K*(h-1)+j-1]*post$surv.train[ , K*(h-1)+j]

    post$surv.train.mean <- apply(post$surv.train, 2, mean)
    
    if(length(x.test)>0) {
        post$x.test <- x.test
        H <- nrow(x.test)/K ## the number of different settings

        if(keepevery>1L) post$yhat.test <- post$yhat.test[thin, ]

        post$surv.test <- 1-pnorm(post$yhat.test)

        for(h in 1:H) for(j in 2:K)
                post$surv.test[ , K*(h-1)+j] <-
                    post$surv.test[ , K*(h-1)+j-1]*post$surv.test[ , K*(h-1)+j]

        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    return(post)
}
