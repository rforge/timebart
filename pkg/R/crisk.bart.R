## run BART and generate CIF probability for cause 1

crisk.bart <- function(
    x.train, y.train=NULL, y.train2=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L), cond=NULL,
    keepcall = FALSE, ## the call object can get rather large
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = NULL,
    binaryOffset2 = NULL,
    ntree = 50L,
    ndpost = 10000L, nskip = 250L,
    printevery = 100L,
    keepevery = 10L, keeptrainfits = FALSE,
    usequants = FALSE, numcut = 100L, printcutoffs = 0L,
    verbose = TRUE,
    id = NULL     
)
{
    if(length(y.train)==0) {
        pre <- crisk.pre.bart(times, delta, x.train, x.test)

        y.train <- pre$y.train
        x.train <- pre$X.train
        x.test  <- pre$X.test
        y.train2 <- pre$y.train2

        times   <- pre$times
        K       <- pre$K

        if(length(cond)==0) cond <- pre$cond
        if(length(binaryOffset)==0) binaryOffset <- pre$binaryOffset
        if(length(binaryOffset2)==0) binaryOffset2 <- pre$binaryOffset2
    }
    else {
        if(length(binaryOffset)==0) binaryOffset <- 0
        if(length(binaryOffset2)==0) binaryOffset2 <- 0

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

    post2 <- bart(x.train=x.train[cond, ], y.train=y.train2[cond], x.test=x.test,
                        keepcall=keepcall, k=k,
                        power=power, base=base,
                        binaryOffset=binaryOffset2,
                        ntree=ntree,
                        ndpost=ndpost, nskip=nskip,
                        printevery=printevery, keepevery=keepevery, keeptrainfits=keeptrainfits,
                        usequants=usequants, numcut=numcut, printcutoffs=printcutoffs,
                        verbose=verbose)

    post$binaryOffset <- binaryOffset
    post$binaryOffset2 <- binaryOffset2
    post$id <- id
    post$times <- times
    post$K <- K
    post$x.train <- x.train
    post$cond <- cond
    post$varcount2 <- post2$varcount

    if(length(x.test)>0) {
        post$x.test <- x.test
        H <- nrow(x.test)/K ## the number of different settings

        post$yhat.test2 <- post2$yhat.test

        post$prob.test <- pnorm(post$yhat.test)
        post$surv.test <- (1-post$prob.test)*(1-pnorm(post$yhat.test2))

        for(h in 1:H) for(j in 2:K) {
                l <- K*(h-1)+j
                
                post$prob.test[ , l] <- post$prob.test[ , l-1]+post$surv.test[ , l-1]*post$prob.test[ , l]
                post$surv.test[ , l] <- post$surv.test[ , l-1]*post$surv.test[ , l]
                      }

        post$prob.test.mean <- apply(post$prob.test, 2, mean)
        post$surv.test.mean <- apply(post$surv.test, 2, mean)
    }

    return(post)
}
