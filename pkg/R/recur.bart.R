## run BART with recurrent events
## 7/12/16

recur.bart <- function(
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

    cat('dbarts\n')
    
    post <- dbarts::bart(x.train=x.train, y.train=y.train, x.test=x.test,
                        sigest = NA_real_, sigdf = 3.0, sigquant = 0.90, 
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
    post$times <- times
    post$K <- K
    post$x.train <- x.train

    if(keepevery>1L) {
        thin <- seq(1, ndpost, keepevery)
        post$yhat.train <- post$yhat.train[thin, ]
        post$varcount <- post$varcount[thin, ]
    }

    post$cum.train <- pnorm(post$yhat.train)
    post$haz.train <- post$cum.train

    H <- nrow(x.train)/K ## the number of different settings

    for(h in 1:H) for(i in 1:K) {
                      j <- (h-1)*K+i
                      if(i==1) post$haz.train[ , j] <- post$haz.train[ , j]/post$times[1]
                      else {
                          post$haz.train[ , j] <- post$haz.train[ , j]/(post$times[i]-post$times[i-1])
                          post$cum.train[ , j] <- post$cum.train[ , j-1]+post$cum.train[ , j]
                      }
                  }

    if(length(x.test)>0) {
        if(keepevery>1L) post$yhat.test <- post$yhat.test[thin, ]

        post$x.test <- x.test
        post$cum.test <- pnorm(post$yhat.test)
        post$haz.test <- post$cum.test

        H <- nrow(x.test)/K ## the number of different settings

        for(h in 1:H) for(i in 1:K) {
                          j <- (h-1)*K+i
                          if(i==1) post$haz.test[ , j] <- post$haz.test[ , j]/post$times[1]
                          else {
                              post$haz.test[ , j] <- post$haz.test[ , j]/(post$times[i]-post$times[i-1])
                              post$cum.test[ , j] <- post$cum.test[ , j-1]+post$cum.test[ , j]
                          }
                      }
    }

    return(post)
}
