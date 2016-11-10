## run BART and generate CIF probability for cause 1

mc.crisk.bart <- function(
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
    seed = 99L,   
    mc.cores = 2L,
    nice=19L    
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

    mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

    while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery
    
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

    Mx <- 2^31-1
    Nx <- nrow(x.train)
    if(Nx*ndpost>Mx) warning('nrow(x.train)*ndpost>2Gi: due to the 2Gi limit in sendMaster,\n',
                             '(unless this limit was increased) reduce ndpost to ', Mx %/% Nx)
    Nx <- nrow(x.test)
    if(Nx*ndpost>Mx) warning('nrow(x.test)*ndpost>2Gi: due to the 2Gi limit in sendMaster,\n',
                             '(unless this limit was increased) reduce ndpost to ', Mx %/% Nx)
    
        for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
              crisk.bart(x.train=x.train, y.train=y.train, y.train2=y.train2, x.test=x.test,
                        cond=cond, keepcall=keepcall, k=k,
                        power=power, base=base,
                        binaryOffset=binaryOffset,
                        binaryOffset2=binaryOffset2,
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
            if(length(post$yhat.test)>0)
                post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)

            if(length(post$yhat.test2)>0)
                post$yhat.test2 <- rbind(post$yhat.test2, post.list[[i]]$yhat.test2)

            if(length(post$prob.test)>0)
                post$prob.test <- rbind(post$prob.test, post.list[[i]]$prob.test)

            if(length(post$prob.test2)>0)
                post$prob.test2 <- rbind(post$prob.test2, post.list[[i]]$prob.test2)

            if(length(post$prob.test12)>0)
                post$prob.test12 <- rbind(post$prob.test12, post.list[[i]]$prob.test12)

            if(length(post$surv.test)>0)
                post$surv.test <- rbind(post$surv.test, post.list[[i]]$surv.test)

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
            post$varcount2 <- rbind(post$varcount2, post.list[[i]]$varcount2)
        }

        if(length(post$prob.test.mean)>0)
            post$prob.test.mean <- apply(post$prob.test, 2, mean)

        if(length(post$prob.test2.mean)>0)
            post$prob.test2.mean <- apply(post$prob.test2, 2, mean)

        if(length(post$prob.test12.mean)>0)
            post$prob.test12.mean <- apply(post$prob.test12, 2, mean)

        if(length(post$surv.test.mean)>0)
            post$surv.test.mean <- apply(post$surv.test, 2, mean)

        return(post)
    }
}
