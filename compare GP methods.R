comparemods <- function(z, X, y, Xtest){
  alpha=2.291052
  phi=1.5  #2.369451
  Nclust <- length(unique(z))
  N <- length(z)
  parm_all <- data.frame('t1'=0,'nug1'=0,'t2'=0,'nug2'=0,'t3'=0,'nug3'=0)
  sigM <- SigM <- matrix(c(2,0,0,2), nrow = 2)
  D <- array(0, dim = c(N,N,ncol(X)))
  for (i in 1:ncol(X)) {
    D[,,i] <- as.matrix(dist(X[,i], diag = TRUE, upper = TRUE))
  }
  Dpred <- array(0, dim = c(nrow(Xtest),nrow(Xtest),ncol(Xtest)))
  for (i in 1:ncol(Xtest)) {
    Dpred[,,i] <- as.matrix(dist(Xtest[,i], diag = TRUE, upper = TRUE))
  }
  dd <- array(0, dim = c(nrow(Xtest),N,ncol(X)))
  for (i in 1:ncol(X)) {
    dd[,,i] <- sqrt(laGP::distance(Xtest, X))
  }
  llgp <- function(D, y, theta, nug){
    n <- length(y)
    if (theta >0 & nug >= 0){
      Sig <- covmat(D=D, ls=theta, nug = nug)
      dsig <- as.numeric(determinant(Sig)[1])
      prtheta <- 0
      for (i in 1:length(theta)) {
        prtheta <- prtheta + dinvgamma(theta[i], 1, 1, log = TRUE)
      }
      prnug <- ifelse(nug==0, 0, dinvgamma(nug, 1, 1, log = TRUE))
      out <- (n/2)*log(2*pi) + (n/2)*dsig + 0.5*t(y)%*%solve(Sig)%*%y - prtheta - prnug
    } else {
      out <- 0
    }
    return(-out)
  }
  
  datdf <- data.frame('x'=as.numeric(X), 'y'=y)
  suppressMessages(suppressWarnings({
    precomp <- brm(y~gp(x, scale=FALSE), data = datdf, 
                   prior = c(prior(constant(0), class=Intercept), 
                             prior(inv_gamma(1,1), class=lscale, coef=gpx), 
                             prior(inv_gamma(1,1), class=sdgp)), 
                   chains = 1, iter = 500, silent = 2, refresh = 0)}))
  newmods <- list()
  
  for (i in 1:Nclust) {
    zi <- which(z==i)
    dat_i <- X[z==i,,drop=FALSE]
    y_i <- y[z==i]
    
    myllgp <- function(x){llgp(D=D[zi,zi,,drop=FALSE], y=y[zi], theta=x[1], nug=x[2])}
    draws <- matrix(c(1,0.1), nrow = 1000, ncol = 2, byrow = T)
    for (k in 2:1000) {
      setTimeLimit(elapsed = 0.1, transient = TRUE)
      on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
      
      #tryCatch({slice_elliptical_mv(x = draws[k-1,], log_target=myllgp, 
      #                              mu=c(2,1), Sig = SigM)$x},
      #         error=function(msg){ print(msg) }) 
      drawk <- tryCatch({slice_elliptical_mv(x = draws[k-1,], log_target=myllgp, 
                                             mu=c(2,1), Sig = SigM)$x},
                        error=function(msg){ draws[k-1,] }) 
      draws[k,] <- ifelse(drawk>0, drawk, 0.001)
    }
    setTimeLimit(elapsed = Inf, transient = TRUE)
    sliceparm_i <- colMeans(draws[501:1000,])
    
    newmods[[i]] <- update(precomp, newdata = datdf[zi,])
    brmparm_i <- posterior_summary(newmods[[i]])
    brmparm_i <- brmparm_i[2:3,1]
    
    gp_i <- newGPsep(dat_i, y_i, d=7.84, g=0.081, dK = TRUE)
    mle_i <- tryCatch({
      mleGPsep(gp_i, param = "both", ab = c(1,1,1,1), tmax = c(500,500))
    }, error = function(msg){
      out <- c()
      for (M in 1:(length(ab)/2)) {
        out <- c(out, rgamma(1, 1, 1))
      }
      out
    })
    if (class(mle_i)=="list"){
      mle_i <- as.vector(mle_i$theta)
    }
    
    parm_all[i,] <- c(mle_i, brmparm_i, sliceparm_i)
  }
  
  # Prediction Section
  C <- exp(-rowSums((dd^2)/rep(phi, each=nrow(Xtest)*N), dims = 2)) #changed Dpred to dd
  predsbyc_opt <- matrix(nrow=nrow(C), ncol = Nclust)
  predsbyc_brm <- matrix(nrow=nrow(C), ncol = Nclust)
  predsbyc_slc <- matrix(nrow=nrow(C), ncol = Nclust)
  var_opt <- matrix(nrow = nrow(C), ncol = Nclust)
  var_brm <- matrix(nrow = nrow(C), ncol = Nclust)
  var_slc <- matrix(nrow = nrow(C), ncol = Nclust)
  probsbyc <- matrix(nrow = nrow(C), ncol = Nclust)
  mlepreds <- c()
  brmpreds <- c()
  slcpreds <- c()
  mlebestpreds <- c()
  brmbestpreds <- c()
  slcbestpreds <- c()
  
  for(k in 1:Nclust){
    opt_i <- predGPsep(k-1, Xtest, lite = TRUE)
    predsbyc_opt[,k] <- opt_i$mean
    var_opt[,k] <- opt_i$s2
    brm_i <- fitted(newmods[[k]], data.frame('x'=as.numeric(Xtest)))
    predsbyc_brm[,k] <- brm_i[,1]
    var_brm[,k] <- brm_i[,2]^2
    z_k <- which(z==k)
    pars <- as.numeric(parm_all[k,5:6])
    Kinv <- solve(covmat(D=D[z_k,z_k,,drop=FALSE], 
                         ls=pars[-length(pars)], 
                         nug = pars[length(pars)]))
    tau2 <- drop(t(y[z_k])%*%Kinv%*%y[z_k])/length(z_k)
    K12 <- drop(covmat(D=dd, set1 = 1:nrow(dd), set2 = z_k,
                      ls = pars[-length(pars)]))
    K11 <- covmat(D=Dpred, ls = pars[-length(pars)], nug = pars[length(pars)])
    slc_i <- K12%*%Kinv%*%y[z_k]
    var_slc[,k] <- tau2*diag(K11 - K12%*%Kinv%*%t(K12))
    predsbyc_slc[,k] <- as.numeric(slc_i)
  }
  for(i in 1:nrow(C)){
    lpz <- function(x){#logprob_zi(zi=i,j=x,z=z,C=C,alpha=alpha)
        n <- length(z)
        lp <- log(n-1) + log(sum(C[i, which(z==x)])) -
          log(sum(C[i,])) - log(n-1+alpha)
        return(lp)
      }
    logp <- sapply(1:Nclust, lpz)
    probs <- exp(logp)/sum(exp(logp))
    probsbyc[i,] <- probs
    mlepreds <- c(mlepreds, sum(predsbyc_opt[i,]*probs))
    brmpreds <- c(brmpreds, sum(predsbyc_brm[i,]*probs))
    slcpreds <- c(slcpreds, sum(predsbyc_slc[i,]*probs))
    mlebestpreds <- c(mlebestpreds, predsbyc_opt[i,which.max(probs)])
    brmbestpreds <- c(brmbestpreds, predsbyc_brm[i,which.max(probs)])
    slcbestpreds <- c(slcbestpreds, predsbyc_slc[i,which.max(probs)])
  }
  
  preds_all <- data.frame('mle_pred'=mlepreds, 'STAN_pred'=brmpreds,
                          'slice_pred'=slcpreds, 'mle_best'=mlebestpreds,
                          'STAN_best'=brmbestpreds, 'slice_best'=slcbestpreds)
  return(list('GPparms'=parm_all, 'ests'=preds_all, 'probs'=probsbyc, 'clust.ests'=
                list('mle'=predsbyc_opt, 'stan'=predsbyc_brm, 'slice'=predsbyc_slc),
              'clust.var'=list('mle'=var_opt, 'stan'=var_brm,  'slice'=var_slc)))
}

library(brms)
library(laGP)
library(qslice)
library(invgamma)
cmp1 <- comparemods(z=slc_basic1$z$z200, X=simX, y=simdat$y,
                    Xtest = matrix(seq(0.01,5,by=0.01),nrow = 500))
cmp2 <- comparemods(z=rep(1,100), X=simX, y=simdat$y,
                    Xtest = matrix(seq(0.01,5,by=0.01),nrow = 500))
cmp8 <- comparemods(z=rep(1:10,each=10)[order(order(as.numeric(simX)))], X=simX, y=simdat$y,
                    Xtest = matrix(seq(0.01,5,by=0.01),nrow = 500))

preds <- cmp8$ests
preds$x <- seq(0.01,5,by=0.01)
preds$true <- sapply(seq(0.01,5,by=0.01), simfn)
ggplot(data = preds, aes(x=x)) + geom_line(aes(y=true), col="black") +
  geom_line(aes(y=mle_pred), col="green") + geom_line(aes(y=STAN_pred), col="blue") +
  geom_line(aes(y=slice_pred), col="red") + xlab('X') + ylab('Y') +
  geom_point(data = simdat, aes(x=x,y=y)) +
  ggtitle('GP Comparison with Reduced Phi')
ggplot(data = predc, aes(x=x)) + geom_line(aes(y=true), col="black", size=1.2) +
  geom_line(aes(y=c1), linetype=2, col="red", size=1.2) + 
  geom_line(aes(y=c2), linetype=2, col="orange", size=1.2) +
  geom_line(aes(y=c3), linetype=2, col="yellow", size=1.2) + 
  geom_line(aes(y=c4), linetype=2, col="green", size=1.2) +
  geom_line(aes(y=c5), linetype=2, col="darkgreen", size=1.2) + 
  geom_line(aes(y=c6), linetype=2, col="blue", size=1.2) +
  geom_line(aes(y=c7), linetype=2, col="lightblue", size=1.2) + 
  geom_line(aes(y=c8), linetype=2, col="purple", size=1.2) +
  geom_line(aes(y=c9), linetype=2, col="magenta", size=1.2) + 
  geom_line(aes(y=c10), linetype=2, col="brown", size=1.2) +
  xlab('X') + ylab('Y') + ggtitle('GP Predictions from Ten Experts')
