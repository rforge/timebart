## run BART and generate CIF probability for cause 1

mc.crisk.bart <- function(
    x.train = matrix(0.0, 0L, 0L),
    y.train=NULL, y.train2=NULL, times=NULL, delta=NULL,
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

    H <- 1
    Mx <- 2^31-1
    Nx <- max(nrow(x.train), nrow(x.test))
    
    if(Nx>Mx%/%ndpost) {
        H <- ceiling(ndpost / (Mx %/% Nx))
        ndpost <- ndpost %/% H
        ##nrow*ndpost>2Gi: due to the 2Gi limit in sendMaster
        ##(unless this limit was increased): reducing ndpost 
    }
    
    mc.cores.detected <- detectCores()

    if(mc.cores>mc.cores.detected) {
        message('The number of cores requested, ', mc.cores,
                       ',\n exceeds the number of cores detected via detectCores() ',
                       'reducing to ', mc.cores.detected)
        mc.cores <- mc.cores.detected
    }

    mc.ndpost <- ((ndpost %/% mc.cores) %/% keepevery)*keepevery

    while(mc.ndpost*mc.cores<ndpost) mc.ndpost <- mc.ndpost+keepevery
    
    post.list <- list()
        
    for(h in 1:H) {
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

        post.list[[h]] <- parallel::mccollect()
    }

    if(H==1 & mc.cores==1) return(post.list[[1]][[1]])
    else {
        for(h in 1:H) for(i in mc.cores:1) {
            if(h==1 & i==mc.cores) post <- post.list[[1]][[mc.cores]]
            else {
                if(length(post$yhat.test)>0)
                    post$yhat.test <- rbind(post$yhat.test, post.list[[h]][[i]]$yhat.test)

                if(length(post$yhat.test2)>0)
                    post$yhat.test2 <- rbind(post$yhat.test2, post.list[[h]][[i]]$yhat.test2)

                if(length(post$prob.test)>0)
                    post$prob.test <- rbind(post$prob.test, post.list[[h]][[i]]$prob.test)

                if(length(post$prob.test2)>0)
                    post$prob.test2 <- rbind(post$prob.test2, post.list[[h]][[i]]$prob.test2)

                if(length(post$prob.test12)>0)
                    post$prob.test12 <- rbind(post$prob.test12, post.list[[h]][[i]]$prob.test12)

                if(length(post$surv.test)>0)
                    post$surv.test <- rbind(post$surv.test, post.list[[h]][[i]]$surv.test)

                post$varcount <- rbind(post$varcount, post.list[[h]][[i]]$varcount)
                post$varcount2 <- rbind(post$varcount2, post.list[[h]][[i]]$varcount2)
                      }

            post.list[[h]][[i]] <- NULL
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
