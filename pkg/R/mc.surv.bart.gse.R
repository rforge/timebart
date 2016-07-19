mc.surv.bart.gse <- function(
                    x.train, times, delta,
                    P=50L, ## number of permutations 
                    R=5L,  ## number of replicates 
                    ntree=20L,  
                    C=1,
                    alpha=0.05,
                    k=2.0,
                    power=2.0, base=0.95,
                    binaryOffset=NULL,
                    ndpost=2000L, nskip=50L,
                    printevery=100L, keepevery=1L, keeptrainfits=FALSE,
                    usequants=FALSE, numcut=100L, printcutoffs=0L,
                    verbose=TRUE,
                    seed = 99L,    
                    mc.cores = 2L, 
                    nice=19L       
                    )
{
    set.seed(seed)

    N <- length(delta)
    K <- ncol(x.train)+1 ## add 1 for time
    
    perm <- matrix(runif(N*P), nrow=N, ncol=P)
    
    prob <- matrix(nrow=P, ncol=K)
    
    h <- 1
    
    for(i in 0:P) {
        if(i==0) {
            times. <- times
            delta. <- delta
        }
        else {
            times. <- times[rank(perm[ , i])]
            delta. <- delta[rank(perm[ , i])]
        }
        
        tmp2 <- matrix(nrow=R, ncol=K)
        
        for(j in 1:R) {
            tmp1 <- mc.surv.bart(x.train=x.train, times=times., delta=delta.,
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
                                 verbose=(i==0 & j==1 & verbose), 
                                 seed=h, mc.cores=mc.cores, nice=nice)$varcount

            tmp2[j, ] <- apply(tmp1, 2, mean)

            h <- h+1
        }

        tmp1 <- apply(tmp2, 2, mean)
        tmp1 <- tmp1/sum(tmp1)
        
        if(i==0) varcount <- tmp1
        else prob[i, ] <- tmp1
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
                C=C, mu.k=mu.k, sd.k=sd.k, iter=iter, perm.prob=prob))
}
               
