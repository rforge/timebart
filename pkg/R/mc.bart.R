## run BART in parallel

mc.bart <- function(
    x.train, y.train=NULL, x.test = matrix(0.0, 0L, 0L),
    keepcall = TRUE, 
    k = 2.0, ## BEWARE: do NOT use k for other purposes below
    power = 2.0, base = 0.95,
    binaryOffset = 0,
    ntree = 200L,
    ndpost = 10000L, nskip = 250L,
    printevery = 100L,
    keepevery = 10L, keeptrainfits = TRUE,
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

    for(i in 1:mc.cores) {
        parallel::mcparallel({psnice(value=nice);
                   bart(x.train=x.train, y.train=y.train, x.test=x.test,
                        keepcall=keepcall, k=k,
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
        keeptestfits <- length(x.test)>0
    
        for(i in 2:mc.cores) {
            if(keeptrainfits) post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)

            if(keeptestfits) post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)

            post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
        }

        if(length(post$yhat.train.mean)>0)
            post$yhat.train.mean <- apply(post$yhat.train, 2, mean)

        if(length(post$yhat.test.mean)>0)
            post$yhat.test.mean <- apply(post$yhat.test, 2, mean)

        return(post)
    }
}
