## run BART with recurrent events in parallel

mc.recur.bart <- function(
    x.train, y.train=NULL, times=NULL, delta=NULL,
    x.test = matrix(0.0, 0L, 0L),
    x.test.short = FALSE, ## you may not need the whole grid
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
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) mc.cores->mc.cores.detected
        ## warning(paste0('The number of cores requested, mc.cores=', mc.cores,
        ##                ',\n exceeds the number of cores detected via detectCores() ',
        ##                'which yields ', mc.cores.detected, ' .'))

    if(length(y.train)==0) {
        recur <- recur.pre.bart(times, delta, x.train, x.test)

        y.train <- recur$y.train
        x.train <- recur$X.train
        x.test  <- recur$X.test

        if(length(binaryOffset)==0) {
            lambda <- sum(delta)/sum(times[ , ncol(times)])
            binaryOffset <- qnorm(1-exp(-lambda))
        }
    }
    else if(length(binaryOffset)==0) binaryOffset <- 0

    mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

    while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery

    ##print(binaryOffset)

    for(i in 1:mc.cores) {
       parallel::mcparallel({psnice(value=nice);
                   recur.bart(x.train=x.train, y.train=y.train,
                              x.test=x.test, x.test.short=x.test.short,
                        k=k, keepcall=keepcall,
                        power=power, base=base,
                        binaryOffset=binaryOffset,
                        ntree=ntree,
                        ndpost=mc.ndpost, nskip=nskip,
                        printevery=printevery, keepevery=keepevery, keeptrainfits=keeptrainfits,
                        usequants=usequants, numcut=numcut, printcutoffs=printcutoffs,
                        verbose=verbose)}, silent=(i!=1))
                                          ## to avoid duplication of output
                                          ## capture stdout from first posterior only
    }

    post.list <- parallel::mccollect()

    post <- post.list[[1]]

    if(mc.cores==1) return(post)
    else {
        for(i in 2:mc.cores) {
            post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)
            post$cum.train <- rbind(post$cum.train, post.list[[i]]$cum.train)
            post$haz.train <- rbind(post$haz.train, post.list[[i]]$haz.train)

            if(length(post$yhat.test)>0) {
                post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)
                post$haz.test <- rbind(post$haz.test, post.list[[i]]$haz.test)
                if(!x.test.short) post$cum.test <- rbind(post$cum.test, post.list[[i]]$cum.test)
            }

            if(length(post$sigma)>0)
                post$sigma <- c(post$sigma, post.list[[i]]$sigma)

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
        }

        post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
        post$cum.train.mean <- apply(post$cum.train, 2, mean)
        post$haz.train.mean <- apply(post$haz.train, 2, mean)

        if(length(post$yhat.test)>0) {
            post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
            post$haz.test.mean <- apply(post$haz.test, 2, mean)
            if(!x.test.short) post$cum.test.mean <- apply(post$cum.test, 2, mean)
        }

        return(post)
    }
}
