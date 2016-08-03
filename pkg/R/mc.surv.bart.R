## run BART and generate survival in parallel
## 07/12/16

mc.surv.bart <- function(
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
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected)
        warning(paste0('The number of cores requested, mc.cores=', mc.cores,
                       ',\n exceeds the number of cores detected via detectCores() ',
                       'which yields ', mc.cores.detected, ' .'))

    if(length(y.train)==0) {
        surv <- surv.pre.bart(times, delta, x.train, x.test)

        y.train <- surv$y.train
        x.train <- surv$X.train
        x.test  <- surv$X.test

        if(length(binaryOffset)==0) {
            lambda <- sum(delta)/sum(times)
            delta  <- mean(surv$times[2:surv$K]-surv$times[1:(surv$K-1)])
            binaryOffset <- qnorm(1-exp(-lambda*delta))
        }
    }
    else if(length(binaryOffset)==0) binaryOffset <- 0

    mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

    while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
              surv.bart(x.train=x.train, y.train=y.train, x.test=x.test,
                        k=k,
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
            if(keeptrainfits) post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)

            if(length(post$surv.train)>0)
                post$surv.train <- rbind(post$surv.train, post.list[[i]]$surv.train)

            if(length(post$yhat.test)>0)
                post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)

            if(length(post$surv.test)>0)
                post$surv.test <- rbind(post$surv.test, post.list[[i]]$surv.test)

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
        }

        if(length(post$yhat.train.mean)>0)
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)

        if(length(post$surv.train.mean)>0)
            post$surv.train.mean <- apply(post$surv.train, 2, mean)

        if(length(post$yhat.test.mean)>0)
            post$yhat.test.mean <- apply(post$yhat.test, 2, mean)

        if(length(post$surv.test.mean)>0)
            post$surv.test.mean <- apply(post$surv.test, 2, mean)

        return(post)
    }
}
