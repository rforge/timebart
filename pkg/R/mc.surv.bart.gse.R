mc.surv.bart.gse <- function(
                    x.train, y.train=NULL, times=NULL, delta=NULL,
                    varcount, ## from previous bart() fit
                    P=100L,    ## number of permutations 
                    ntree=20L,  
                    C=1,
                    alpha=0.05,
                    k=2.0,
                    power=2.0, base=0.95,
                    binaryOffset=NULL,
                    ndpost=1000L, nskip=250L,
                    printevery=100L, keepevery=1L, keeptrainfits=TRUE,
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

    if(length(y.train)==0) {
        surv <- surv.pre.bart(times, delta, x.train)

        y.train <- surv$y.train
        x.train <- surv$X.train
    }

    N <- length(y.train)
    K <- ncol(x.train)
    P <- max(P, mc.cores)
    
    perm <- matrix(runif(N*P), nrow=N, ncol=P)
    
    post <- list()
    post.list <- list()
    
    j <- 1
    h <- 1
    
    while(j<=P) {
        for(i in 1:mc.cores) if(j<=P) {
            R <- rank(perm[ , j])
            
            y. <- y.train[R]
            
            parallel::mcparallel({tools::psnice(value=nice);
                surv.bart(x.train=x.train, y.train=y.,
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
                            verbose=verbose)},
                       silent=(i!=1))
            ## to avoid duplication of output
            ## capture stdout from first posterior only

            j <- j+1
        }
        
        post.list <- parallel::mccollect()

        for(i in 1:length(post.list)) {
            post[[h]] <- post.list[[i]]
            h <- h+1
        }
    }

    print(paste0('Number of permutations=', length(post)))

    P <- length(post)
    
    for(j in 1:P) {
        total <- apply(post[[j]], 1, sum)

        post[[j]] <- post[[j]]/total

        post[[j]] <- apply(post[[j]], 2, mean)

        if(j==1) prob <- matrix(post[[j]], nrow=1, ncol=K)
        else prob <- rbind(prob, post[[j]])
    }

    mu.k <- apply(prob, 2, mean)
    sd.k <- apply(prob, 2, sd)

    R <- double(K)

    iter <- 0

    while(min(R)<(1-alpha)) {
        if(iter>0) {
            C <- C*1.01
            R <- R*0
        }
        
        for(i in 1:P) for(j in 1:K) {
            R[j] <- R[j]+(prob[i, j]<=(mu.k[j]+C*sd.k[j]))/P
        }

        iter <- iter+1
    }

    if(iter==1) {
        warning('Algorithm stopped at iteration 1.  Try again with a smaller C.')
        return(post.list)
    }

    varcount <- varcount/apply(varcount, 1, sum)

    prob <- apply(varcount, 2, mean)
    
    return(list(which=which(prob>(mu.k+C*sd.k)), prob=prob,
                C=C, mu.k=mu.k, sd.k=sd.k, iter=iter))
}
               
