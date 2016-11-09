## run BART and generate survival

surv.bart <- function(
    x.train, y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L),
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
    id = NULL,     ## only used by surv.bart
    seed = 99L,    ## only used by mc.surv.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
)
{
    if(length(y.train)==0) {
        pre <- surv.pre.bart(times, delta, x.train, x.test)

        y.train <- pre$y.train
        x.train <- pre$X.train
        x.test  <- pre$X.test

        times   <- pre$times
        K       <- pre$K

        if(length(binaryOffset)==0) binaryOffset <- pre$binaryOffset
    }
    else {
        if(length(binaryOffset)==0) binaryOffset <- 0

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    post <- bart(x.train=x.train, y.train=y.train, x.test=x.test,
                        keepcall=keepcall, k=k,
                        power=power, base=base,
                        binaryOffset=binaryOffset,
                        ntree=ntree,
                        ndpost=ndpost, nskip=nskip,
                        printevery=printevery, keepevery=keepevery, keeptrainfits=keeptrainfits,
                        usequants=usequants, numcut=numcut, printcutoffs=printcutoffs,
                        verbose=verbose)

    post$binaryOffset <- binaryOffset
    post$id <- id
    post$times <- times
    post$K <- K
    post$x.train <- x.train

    ## if(keepevery>1L) { ## manual thinning needed for dbarts < 0.8-6
    ##     thin <- seq(1, ndpost, keepevery)
    ##     if(keeptrainfits) post$yhat.train <- post$yhat.train[thin, ]
    ##     post$varcount <- post$varcount[thin, ]
    ## }

    if(keeptrainfits) {
        post$surv.train <- 1-pnorm(post$yhat.train)

        H <- nrow(x.train)/K ## the number of different settings

        for(h in 1:H) for(j in 2:K) {
                l <- K*(h-1)+j

                post$surv.train[ , l] <- post$surv.train[ , l-1]*post$surv.train[ , l]
                      }

        post$surv.train.mean <- apply(post$surv.train, 2, mean)
    }

    if(length(x.test)>0) {
        post$x.test <- x.test
        H <- nrow(x.test)/K ## the number of different settings

        ##if(keepevery>1L) post$yhat.test <- post$yhat.test[thin, ]

        post$surv.test <- 1-pnorm(post$yhat.test)

        for(h in 1:H) for(j in 2:K) {
                l <- K*(h-1)+j

                post$surv.test[ , l] <- post$surv.test[ , l-1]*post$surv.test[ , l]
                      }

        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    return(post)
}
