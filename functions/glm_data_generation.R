pacman::p_load(Matrix)
generate_data_glm <- function(n,p,cov_setting=c("ind","comsym","ar1","group","factor","extreme"),
                              family=gaussian("identity"),
                              ntest=0,a=NULL,ind=NULL,beta=NULL,alpha=NULL,signal_strength=2,avg_exp=0.5,
                              rho=NULL,xmean=rep(0,p)) {
  
  # first select active vars and generate beta
  bSb <- 0
  Ft <- NULL
  cov_setting <- match.arg(cov_setting)
  if (is.null(beta)) {
    if (is.null(a)) {
      if (is.null(ind)) {
        a <- round(4*log(p))
      } else {
        a <- length(ind)
      }
    }
    if (!cov_setting=="extreme") {
      if (is.null(ind)) {
        ind <- sample(1:p,a,replace=FALSE)
      }
      beta <- Matrix(data=c(0),p,1,sparse = TRUE)
      beta[ind] <- sample(c(1,-1),a,replace = TRUE,prob = c(0.6,0.4)) * (4*log(n)/sqrt(n) + abs(rnorm(a)))
    }
  } else {
    ind <- which(beta!=0)
    a <- length(ind)
  }
  
  if (is.null(rho)) {
    if (cov_setting=="ar1") {
      rho <- 0.9
    }
    if (cov_setting=="comsym") {
      rho <- 0.5
    }
  }
  # then generate x with given cov structure and calculate bSb = beta' Sigma beta
  # efficient implementation for high p
  switch(cov_setting,
         ind={
           x <- matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)
           bSb <- sum(beta^2)
         },
         comsym={
           x <- sqrt(rho)*matrix(rep(rnorm((n+ntest),0,1),p),n+ntest,p) + sqrt(1-rho)*matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)
           bSb <- rho*sum(beta)^2 + (1-rho)*sum(beta^2)
         },
         ar1={ # explizit cholesky dec from Ben Bolker 09 Oct 2019 (https://bbolker.github.io/mixedmodels-misc/notes/varmats.html)
           R <- matrix(0,p,p)
           R[1,] <- rho^(0:(p-1))        ## formula for 1st row
           cc <- sqrt(1 - rho^2);        ## scaling factor: c^2 + rho^2 = 1
           R2 <-  cc * R[1,]             ## formula for 2nd row */
           for (j in 2:p) {              ## shift elements in 2nd row for remaining rows
             R[j, j:p] <- R2[1:(p-j+1)] 
           }
           # crossprod(R)
           x <- matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)%*%R
           bSb <- sum((R%*%beta)^2)
         },
         group={
           x <- matrix(c(0),n+ntest,p)
           nb <- p%/%100 ## number of blocks
           
           if (p > 100) {
             R <- matrix(0,100,100)
             R[1,] <- 0.9^(0:(100-1))     
             cc <- sqrt(1 - 0.9^2);        
             R2 <-  cc * R[1,]     
             for (j in 2:100) {  
               R[j, j:100] <- R2[1:(100-j+1)] 
             }
             for (i in 1:(nb-1)) {
               if (i <= (nb-1)/2 ) {
                 x[,(1:100)+(i-1)*100] <- sqrt(0.5)*matrix(rep(rnorm((n+ntest),0,1),100),n+ntest,100) + 
                   sqrt(1-0.5)*matrix(rnorm((n+ntest)*100,0,1),n+ntest,100)
                 bSb <- bSb + 0.5*sum(beta[(1:100)+(i-1)*100])^2 + (1-0.5)*sum(beta[(1:100)+(i-1)*100]^2)
               } else {
                 x[,(1:100)+(i-1)*100] <- matrix(rnorm((n+ntest)*100,0,1),n+ntest,100)%*%R
                 bSb <- bSb + sum((R%*%beta[(1:100)+(i-1)*100])^2)
               }
             }
             x[,-(1:((nb-1)*100))] <- matrix(rnorm((n+ntest)*(p-(nb-1)*100)),n+ntest,p-(nb-1)*100)
             bSb <- bSb + sum(beta[-(1:((nb-1)*100))]^2)
           } else {
             x <- matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)
             bSb <- sum(beta^2)
           }
         },
         factor={
           k <- a
           Ft <- matrix(rnorm(p*k,0,1),k,p)
           x <- matrix(rnorm((n+ntest)*k,0,1),n+ntest,k)%*%Ft + 0.1*matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)
           bSb <- sum((Ft%*%beta)^2) + 0.01*sum(beta^2)
         },
         extreme={
           if (is.null(beta)) {
             ind <- 1:a
             beta <- Matrix(data=c(0),p,1,sparse = TRUE)
             beta[ind] <- 1
           } 
           x <- matrix(c(0),n+ntest,p)
           z <- matrix(rnorm((n+ntest)*p,0,1),n+ntest,p)
           w <- matrix(rnorm((n+ntest)*a,0,1),n+ntest,a)
           x[,1:a] <- (z[,1:a]+w)/sqrt(2)
           x[,-(1:a)] <- (z[,-(1:a)]+rowSums(z[,1:a]))/sqrt(a+1)
           bSb <- sum(beta[ind]^2)
         },
         stop("Not an allowed covariance setting, choose one of ind, comsym, ar1, group, factor, extreme!")
  )
  stopifnot(nrow(x)==(n+ntest),ncol(x)==p)
  x <- scale(x,center = -xmean,scale = FALSE)
  # scale beta such that bSb = signal_strength
  beta <- beta*sqrt(signal_strength)/sqrt(bSb)
  xbeta <- x[,ind,drop=FALSE]%*%beta[ind]
  
  # finally generate responses
  if (is.null(alpha)) {
    # set such that we reach avg expected value
    alpha <- tryCatch( uniroot(function(interc) mean(family$linkinv(interc + xbeta))-avg_exp, c(-max(abs(xbeta))-7,-max(abs(xbeta))+7),extendInt = "yes" )$root,
                       error=function(error_message) {
                         message(paste0("Error in uniroot for alpha: ",error_message))
                         return(0)
                       })
  }

  # y_tmp <- rbinom(n+ntest,1,1/(1+exp(-alpha - xbeta)))
  # cbind(rowMeans(sapply(1:10000,function(i)rbinom(n+ntest,1,1/(1+exp(-alpha - xbeta))) )),
  #       rowMeans(sapply(1:10000,function(i)family$simulate(glm(rep(0,n+ntest)~0,offset = alpha+xbeta,family = family),nsim = 1) )),
  #       1/(1+exp(-alpha-xbeta)))
  
  if ("simulate" %in% names(family)) {
    y <- family$simulate(glm(rep(1,n+ntest)~0,offset = alpha+xbeta,family = family),nsim = 1)
  } else {
    if (family$family=="gaussian") {
      if (family$link =="log") {
        mysd <- 0.1
      } else {
        mysd <- 1
      }
      y <- family$linkinv(alpha+xbeta ) + rnorm(n+ntest,0,mysd)
    } else {
      stop("Family can't be simulated!")
    }
  }
  
  return(list(x=x[1:n,],y=y[1:n],xtest=x[-(1:n),],ytest=y[-(1:n)],
              ind=ind,alpha=alpha,beta=beta,cov_setting=cov_setting,
              Ft=Ft,bSb=bSb))
}

# data <- generate_data_glm(200,2000,family=poisson(log),ntest=100,avg_exp = 4)

