## run DSS
## 07/22/16

surv.bart.dss <- function(
    post,          ## posterior from surv.bart/mc.surv.bart
                   ## which includes x.train and yhat.train
    mc.cores = 2L, ## for mc.surv.bart.dss only
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

    stopifnot(K>1)

    while(j==1) {
        for(k in 2:K) {
            if(k %in% h) R2[k] <- 0
            else {
                fit <- rpart(post$yhat.train.mean~post$x.train[ , c(k, h)])

                if(sd(predict(fit))==0) R2[k] <- 0
                else R2[k] <- cor(post$yhat.train.mean, predict(fit))^2

                fit <- NULL
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
