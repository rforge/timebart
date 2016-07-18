mc.surv.bart.gse <- function(
                    x.train, times, delta,
                    P=50L, ## number of permutations 
                    ntree=20L,  
                    C=1,
                    alpha=0.05,
                    k=2.0,
                    power=2.0, base=0.95,
                    binaryOffset=NULL,
                    ndpost=10000L, nskip=50L,
                    printevery=100L, keepevery=1L, keeptrainfits=FALSE,
                    usequants=FALSE, numcut=100L, printcutoffs=0L,
                    verbose=TRUE,
                    seed = 99L,    
                    mc.cores = 2L, 
                    nice=19L       
                    )
{
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
    parallel::mc.reset.stream()

    mc.cores.detected <- parallel::detectCores()

    if(mc.cores>mc.cores.detected)
        warning(paste0('The number of cores requested, mc.cores=', mc.cores,
                       ',\n exceeds the number of cores detected via detectCores() ',
                       'which yields ', mc.cores.detected, ' .'))        

    N <- length(delta)
    K <- ncol(x.train)+1 ## add 1 for time
    
    perm <- matrix(runif(N*P), nrow=N, ncol=P)
    
    post <- list()
    post.list <- list()
    
    j <- 0
    h <- 1
    l <- 0
    
    while(j<=P) {
        for(i in 1:mc.cores) if(j<=P) {
            if(j==0) {
                times. <- times
                delta. <- delta
            }
            else {
                times. <- times[rank(perm[ , j])]
                delta. <- delta[rank(perm[ , j])]
            }
            
            parallel::mcparallel({tools::psnice(value=nice);
                surv.bart(x.train=x.train, times=times., delta=delta.,
                            k=k,
                            power=power, base=base,
                            binaryOffset=binaryOffset,
                            ntree=ntree,
                            ndpost=ndpost, nskip=nskip,
                            printevery=printevery,
                            keepevery=keepevery,
                            keeptrainfits=keeptrainfits,
                            usequants=usequants, numcut=numcut,
                            printcutoffs=printcutoffs,
                            verbose=verbose, id=j)},
                       silent=(j!=0))
            ## to avoid duplication of output
            ## capture stdout from first posterior only

            j <- j+1
        }
        
        post.list <- parallel::mccollect()

        return(post.list)

        for(i in 1:length(post.list)) {
            if(post.list[[i]]$id==0) l <- h
            post[[h]] <- post.list[[i]]$varcount
            h <- h+1
        }

        post.list <- NULL
    }

    varcount <- post[[l]]/apply(post[[l]], 1, sum)

    varcount <- apply(varcount, 2, mean)

    #remove l-th item from list corresponding to unpermuted times/delta
    post[[l]] <- NULL
    
    for(j in 1:P) {
        total <- apply(post[[j]], 1, sum)

        post[[j]] <- post[[j]]/total

        post[[j]] <- apply(post[[j]], 2, mean)

        if(j==1) prob <- matrix(post[[j]], nrow=1, ncol=K)
        else prob <- rbind(prob, post[[j]])
    }

    mu.k <- apply(prob, 2, mean)
    sd.k <- apply(prob, 2, sd)

    cov.prob <- double(K)

    iter <- 0

    while(min(cov.prob)<(1-alpha)) {
        if(iter>0) {
            C <- C*1.01
            cov.prob <- cov.prob*0
        }
        
        for(i in 1:P) for(j in 1:K) {
            cov.prob[j] <- cov.prob[j]+(prob[i, j]<=(mu.k[j]+C*sd.k[j]))/P
        }

        iter <- iter+1
    }

    if(iter==1) 
        warning('Algorithm stopped at iteration 1.  Try again with a smaller C.')
    
    return(list(which=which(varcount>(mu.k+C*sd.k)), prob=varcount,
                C=C, mu.k=mu.k, sd.k=sd.k, iter=iter, perm.prob=post))
}
               
