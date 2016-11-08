## run BART and generate survival

crisk.bart <- function(
    x.train, y.train1=NULL, y.train2=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L),
    keepcall = FALSE, ## the call object can get rather large
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset1 = NULL,
    binaryOffset2 = NULL,
    ntree = 50L,
    ndpost = 10000L, nskip = 250L,
    printevery = 100L,
    keepevery = 10L, keeptrainfits = TRUE,
    usequants = FALSE, numcut = 100L, printcutoffs = 0L,
    verbose = TRUE,
    id = NULL,     ## only used by crisk.bart
    seed = 99L,    ## only used by mc.crisk.bart
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
)
{
    if(length(y.train1)==0) {
        if(length(binaryOffset1)==0) {
            lambda <- sum(delta==1)/sum(times)
            binaryOffset1 <- qnorm(1-exp(-lambda))
        }

        if(length(binaryOffset2)==0) {
            lambda <- sum(delta==2)/sum(times)
            binaryOffset2 <- qnorm(1-exp(-lambda))
        }

        crisk <- crisk.pre.bart(times, delta, x.train, x.test)

        y.train1 <- crisk$y.train1
        y.train2 <- crisk$y.train2
        x.train <- crisk$X.train
        x.test  <- crisk$X.test

        times   <- crisk$times
        K       <- crisk$K
    }
    else {
        if(length(binaryOffset1)==0) binaryOffset1 <- 0
        if(length(binaryOffset2)==0) binaryOffset2 <- 0

        times <- unique(sort(x.train[ , 1]))
        K     <- length(times)
    }

    post1 <- bart(x.train=x.train, y.train=y.train1, x.test=x.test,
                        keepcall=keepcall, k=k,
                        power=power, base=base,
                        binaryOffset=binaryOffset1,
                        ntree=ntree,
                        ndpost=ndpost, nskip=nskip,
                        printevery=printevery, keepevery=keepevery, keeptrainfits=keeptrainfits,
                        usequants=usequants, numcut=numcut, printcutoffs=printcutoffs,
                        verbose=verbose)

    forced <- which(y.train1==0)
    
    post2 <- bart(x.train=x.train[forced, ], y.train=y.train2[forced, ], x.test=x.test,
                        keepcall=keepcall, k=k,
                        power=power, base=base,
                        binaryOffset=binaryOffset2,
                        ntree=ntree,
                        ndpost=ndpost, nskip=nskip,
                        printevery=printevery, keepevery=keepevery, keeptrainfits=keeptrainfits,
                        usequants=usequants, numcut=numcut, printcutoffs=printcutoffs,
                        verbose=verbose)

    ##post$call <- NULL
    post$binaryOffset <- NULL
    post$id <- id
    post$times <- times
    post$K <- K
    post$x.train <- x.train

    ## if(keepevery>1L) { ## thinning with dbarts not available
    ##     thin <- seq(1, ndpost, keepevery)
    ##     if(keeptrainfits) post$yhat.train <- post$yhat.train[thin, ]
    ##     post$varcount <- post$varcount[thin, ]
    ## }

    if(keeptrainfits) {
        post$surv.train <- 1-pnorm(post$yhat.train)

        H <- nrow(x.train)/K ## the number of different settings

        for(h in 1:H) for(j in 2:K)
                      post$surv.train[ , K*(h-1)+j] <-
                          post$surv.train[ , K*(h-1)+j-1]*post$surv.train[ , K*(h-1)+j]

        post$surv.train.mean <- apply(post$surv.train, 2, mean)
    }

    if(length(x.test)>0) {
        post$x.test <- x.test
        H <- nrow(x.test)/K ## the number of different settings

        ##if(keepevery>1L) post$yhat.test <- post$yhat.test[thin, ]

        post$surv.test <- 1-pnorm(post$yhat.test)

        for(h in 1:H) for(j in 2:K)
                post$surv.test[ , K*(h-1)+j] <-
                    post$surv.test[ , K*(h-1)+j-1]*post$surv.test[ , K*(h-1)+j]

        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    return(post)
}
