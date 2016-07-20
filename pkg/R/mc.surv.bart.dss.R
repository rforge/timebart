## run DSS in parallel
## 07/20/16

mc.surv.bart.dss <- function(
    post,          ## posterior from surv.bart/mc.surv.bart                         
                   ## which includes x.train and yhat.train
    mc.cores = 2L, ## ditto
    nice=19L       ## ditto
)
{
    if(length(post$yhat.train.mean)==0)
        post$yhat.train.mean <- apply(post$yhat.train, 2, mean)

    names. <- dimnames(post$x.train)[[2]]

    K <- length(names.)

    R2 <- double(K)

    h <- 1 ## start with time only
    j <- 1
    
    if(K>1) while(j==1) {
        l <- 0
        
        for(k in 2:K) {
            if(k %in% h) R2[k] <- 0
            else {
                if(l[1]==0) l <- k
                else l <- c(k, l) ## LIFO
                    
                parallel::mcparallel({tools::psnice(value=nice);
                    rpart::rpart(post$yhat.train.mean~post$x.train[ , c(k, h)])})
            }

            if(l[1]==0) j <- 0
            else j <- length(l)
                
            if(j==mc.cores | (j>0 & k==K)) { 
                fit <- parallel::mccollect()

                i <- 1
                    
                for(j in l) { ## LIFO
                    if(sd(predict(fit[[i]]))==0) R2[j] <- 0
                    else R2[j] <- cor(post$yhat.train.mean, predict(fit[[i]]))^2

                    i <- i+1
                }

                fit <- NULL
                l <- 0
            }
        }
        
        k <- which(R2==max(R2))

        j <- length(k)
        
        if(j==1){
            i <- length(h)+1
            h[i] <- k
            print(c(h[i], R2[k]))
            print(names.[h[i]])
        }
        else print(c(k, R2[k]))
    }

    return(list(pick=h, names=names.[h]))
}
