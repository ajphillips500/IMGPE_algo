# iMGPE Gibbs step functions
###################################

# log posterior for alpha
logpost_alpha <- function(alpha, N, k){
  value <- log(gamma(alpha)) + (k-1.5)*log(alpha) -(1/(2*alpha)) -
    log(gamma(N+alpha))
  return(value)
}

# log posterior for alpha, MCRP
lp_alpha <- function(alpha, z, C){
  n <- length(z)
  sumlp <- 0
  for (i in 1:n) {
    if(i==1){ zprime=c(-1,z[-1]) } else {
      zprime <- c(1:(i-1),-1,z[-c(1:i)])
    }
    if(sum(zprime==z[i])==0){
      sumlp <- sumlp + log(alpha)
    } else {
      numerator <- log(sum(C[z[i],which(zprime==z[i])]))
      denom <- log(sum(C[z[i],-z[i]]))
      sumlp <- sumlp +numerator - denom
    }
  }
  return(sumlp - n*log(n+alpha-1))
}

# log posterior for phi (approximate)
logpost_phi <- function(phi, alpha, z, D, Xq=NULL, lp=TRUE){
  n <- length(z)
  if (is.null(Xq)){
    if (length(dim(D))==2){
      C <- exp(-(D^2)/phi)
    } else {
      C <- exp(-rowSums((D^2)/rep(phi, each=n^2), dims = 2))
    }
  } else {
    # must code different covariance matrix for qualitative variables
    C <- exp(-rowSums((D^2)/rep(phi, each=n^2), dims = 2))
  }
  
  step1 <- unlist(lapply(1:n, logprob_zi, j=NULL, z=z, C=C, alpha=alpha))
  step1[step1== -Inf] <- -1000
  step2 <- max(-1000, log(phi))
  step3 <- max(-1000, log(phi*sqrt(2*pi)))
  step4 <- sum(step1) - 0.5*sum(step2^2) - step3
  #step4 <- sum(step1) + dgamma(phi, 1, 1, log = TRUE)
  if(!lp){ step4 <- exp(step4) }
  return(step4)
}

# CRP functions
#######################################

# log probability of z_i=k | z_-i, etc.
logprob_zi <- function(zi, j=NULL, z, C, alpha){
  if(is.null(j)){ j=z[zi] }
  n <- length(z)
  #lp <- log(n-1) + log(sum(C[zi, which(z==j)[!(which(z==j)==zi)]])) -
  #  log(sum(C[zi,-zi])) - log(n-1+alpha)
  if(sum(z==j)==1){
    lp = log(alpha) - log(n - 1 + alpha)
  } else {
    lp <- log(n-1) + log(sum(C[zi, which(z==j)[!(which(z==j)==zi)]])) -
      log(sum(C[zi,-zi])) - log(n-1+alpha)
  }
  return(lp)
}

# function to construct covariance matrix from Gaussian kernel function
covmat <- function(D, set1=1:dim(D)[1], set2=1:dim(D)[2], ls, nug=0){
  D <- D[set1,set2,,drop=FALSE]
  N1 <- length(set1)
  N2 <- length(set2)
  dim(D) <- c(N1,N2,length(ls))
  C <- exp(-rowSums((D^2)/rep(ls, each=N1*N2), dims = 2)) +
    diag(nug, nrow = N1, ncol = N2)
  return(C)
}

# remove row & col i from precision matrix
inverse.update <- function(x, i) {
  a <- x[-i,-i, drop=FALSE]
  b <- x[-i,i, drop=FALSE]
  c <- x[i,-i, drop=FALSE]
  d <- x[i,i]
  a - b %*% c / d # For production code, should throw an error when d is 0.
}

# Functions for Distance Dependent CRP
##########################################

# Take log without distortion
safelog <- function (x) {
  safelog.f <- function (x)
    if (x == Inf)
      Inf
  else if (x == 0)
    -100
  else
    log(x)
  
  if (length(x) == 1)
    safelog.f(x)
  else
    sapply(x, safelog.f)
}

# Get the log of a sum without distortion
log.sum <- function(v) {
  log.sum.pair <- function(x,y)
  {
    if ((y == -Inf) && (x == -Inf))
    { return(-Inf); }
    if (y < x) return(x+log(1 + exp(y-x)))
    else return(y+log(1 + exp(x-y)));
  }
  if (length(v) == 1)
    return(v)
  r <- v[1];
  for (i in 2:length(v))
    r <- log.sum.pair(r, v[i])
  return(r)
}

# Get all points that point i connects to
connections <- function(i, links)
{
  visited <- c()
  to.visit <- c(i)
  while (length(to.visit) > 0) {
    curr <- to.visit[1]
    visited <- c(visited, curr)
    to.visit <- to.visit[-1]
    pointers <- which(links == curr)
    for (p in pointers) {
      if (!(p %in% visited))
        to.visit <- c(to.visit, p)
    }
  }
  visited
}

# Generate Data from IMGPE Model
##############################################

gendat <- function(seed=12345, X, alp=NULL, phi=NULL, z=NULL, gpparm=NULL){
  require(MASS)
  set.seed(seed)
  # Generate phi and alpha from priors
  if(is.null(alp)){ alpha <- invgamma::rinvgamma(1, 1, 1) }
  if(is.null(phi)){ phi <- rlnorm(1) }
  # Generate z from phi, alpha
  D <- as.matrix(dist(X, diag = TRUE, upper = TRUE))
  C <- exp(-D^2/phi)
  if (is.null(z)){
    z <- rep(1,nrow(X))
    clusters <- c()
    for (i in 2:nrow(X)) {
      n <- i-1
      zi <- z[1:n]
      nclust <- length(unique(zi))
      probs <- c()
      for (j in 1:nclust) {
        newprob <- log(n-1) - log(n-1+alpha) + 
          log(sum(C[i, which(zi==j)])) - log(sum(C[i,zi]))
        probs <- c(probs, newprob)
      }
      probs <- c(probs, log(alpha)-log(n-1+alpha))
      probs <- exp(probs)
      newi <- sample(1:(nclust+1), size = 1, prob = probs)
      z[i] <- newi
    }
  }
  
  # Generate theta from prior for each cluster
  nclust <- length(unique(z))
  theta <- matrix(0, nrow = nclust, ncol = 3)
  if (is.null(gpparm)){
    for (i in 1:nclust) {
      theta[i,] <- c(rgamma(1,1,1), rgamma(1,1,10), 1)
    }
  } else {
    for (i in 1:nclust) {
      theta[i,] <- gpparm
    }
  }
  # Generate y from theta, z
  # Generate y for all X per cluster
  DD <- array(data = D, dim = c(nrow(D), ncol(D), 1))
  #Sig <- matrix(0, nrow = nrow(D), ncol = nrow(D))
  y <- rep(0, nrow(D))
  for(j in 1:nclust) {
    Sigj <- covmat(DD, ls=theta[i,1], nug = theta[i,2])
    Sigj <- Sigj*theta[i,3]
    #Sig[which(z==j), which(z==j)] <- Sigj
    yj <- as.vector(mvrnorm(n=1, mu = rep(0,length(z)), Sigma = Sigj))
    y[which(z==j)] <- yj[which(z==j)]
  }
  #y <- mvtnorm::rmvnorm(1, mean = rep(0,nrow(D)), sigma = Sig)
  
  return(list('X'=X, 'alpha'=alpha, 'phi'=phi, 'z'=z, 'theta'=theta, 'y'=y))
}

# IMGPE MCMC Fitting Functions
##############################################

# IMGPE-EM with GP fitting w/laGP
imgpe.em <- function(X, Xq=NULL, y, parms=NULL, z_init=NULL, draw = 10, Xpred=NULL,
                      maxSize=nrow(X), maxIters=5000){
  require(laGP)
  require(qslice)
  require(numDeriv)
  
  N <- nrow(X)
  d <- ncol(X)
  dq <- ifelse(is.null(Xq), 0, ncol(Xq))
  alpha <- ifelse(!is.null(parms), parms[1], 0.01*N)
  if (!is.null(parms)){
    phi <- parms[-1]
  } else {
    phi <- rep(1, d+dq)
  }
  ab <- c(darg(NULL, X)$ab, garg(NULL, y)$ab)
  if(sum(ab[(length(ab)-1):length(ab)])==0){
    ab[(length(ab)-1):length(ab)] <- c(1, 1)
  }
  z <- z_init
  if(is.null(z)){z <- rep(1,N)}
  n_k <- as.vector(table(z))
  Nclust <- length(n_k)
  D <- array(0, dim = c(N,N,d+dq))
  for (i in 1:d) {
    D[,,i] <- as.matrix(dist(X[,i], diag = T, upper = T))
  }
  if (!is.null(Xq)){
    for (i in 1:dq) {
      Di <- as.matrix(cluster::daisy(Xq[,i,drop=FALSE], metric = "gower"))
      Di[Di > 0] <- 1
      D[,,d+i] <- Di
    }
  }
  if (!is.null(Xpred)){
    Dpred <- array(0, dim = c(nrow(Xpred),N,d+dq))
    for (i in 1:d){
      Dpred[,,i] <- as.matrix(distance(Xpred[,1:d], X))
    }
    if (!is.null(Xq)) {
      for (i in 1:dq){
        Di <- matrix(0, nrow = nrow(Xpred), ncol = N)
        for (A in 1:nrow(Xpred)){
          for (B in 1:N) {
            Di[A,B] <- ifelse(Xpred[A,d+i]==Xq[B,i], 0, 1)
          }
        }
        Dpred[,,d+i] <- Di
      }
    }
  }
  Qlist <- list()
  maxranges <- apply(X, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  
  # prepare results
  res_gatepar <- data.frame(matrix(c(alpha, phi), nrow=1))
  colnames(res_gatepar) <- c("alpha", colnames(X), colnames(Xq))
  res_z <- data.frame(z0=z)
  res_gp <- data.frame()
  res_noise <- data.frame('nug0'=rep(0,N))
  res_spvar <- data.frame('spv0'=rep(0,N))
  if (!is.null(Xpred)){
    res_drawmu <- data.frame('start'=rep(0,nrow(Xpred)))
    res_pred <- data.frame('start'=rep(0,nrow(Xpred)))
  }
  
  for (j in 1:maxIters) {
    for (i in 1:nrow(X)) {
      # Gibbs sampling pass over indicator variables z
      z_i <- z[i] # what is the nth persons table assignment?
      n_k[z_i] <- n_k[z_i] - 1 # remove the nth person from table
      # update the matching precision matrix
      if (j>=2){
        j_index <- which((1:N)[z==z_i]==i)
        Qlist[[z_i]] <- inverse.update(Qlist[[z_i]], j_index)
        if (n_k[z_i] == 0){
          Qlist <- Qlist[-z_i]
          res_gp <- res_gp[-z_i,]
        }
      }
      # if the table becomes empty when the nth person is removed,
      # then that table/cluster is removed.
      if( n_k[z_i] == 0 ) {
        n_k <- n_k[-z_i]    # take out the empty cluster
        z[z>=z_i] <- z[z>=z_i]-1 # move up all other clusters
        Nclust <- Nclust - 1     # decrease total number of clusters by 1
      }
      z[i] <- -1
      
      logp <- rep(NA, Nclust+1)
      C <- exp(-rowSums((D^2)/rep(phi, each=N^2), dims = 2))
      for (k in 1:Nclust) {
        if (length(Qlist) >= k){
          z_k <- which(z==k)
          if (length(z_k) != nrow(Qlist[[k]])){
            print(paste('k:',k))
            print(table(z))
            print(nrow(Qlist[[k]]))
            break
          }
          pars <- as.numeric(res_gp[1,-ncol(res_gp)])
          Qy <- drop(covmat(D=D[,,1:d,drop=FALSE], set1 = i, set2 = z_k,
                            ls = pars))
          #print(paste('zk:',length(z_k),'Qy:',length(Qy),'Qlist:',nrow(Qlist[[k]])))
          mu <- t(Qy)%*%(Qlist[[k]])%*%y[z_k]
          spvar <- drop(t(y[z_k])%*%(Qlist[[k]])%*%y[z_k])/length(z_k)
          noise <- res_gp[k,ncol(res_gp)]
          sig2 <- spvar + noise - t(Qy)%*%(Qlist[[k]])%*%Qy
          sig2 <- ifelse(sig2<=0, 1e-4, sig2)
          lcondprob <- dnorm(y[i], mean = mu, sd = sqrt(sig2), log = TRUE)
        } else {
          nug <- 0
          for (M in 1:(length(ab)/2)) {
            nug <- nug + rgamma(1, ab[2*M-1], ab[2*M])
          }
          lcondprob <- dnorm(y[i], mean = 0, sd = sqrt(nug), log = TRUE)
        }
        logp[k] <- lcondprob + logprob_zi(i, j=k, z=z, C=C, alpha=alpha)
      }
      
      M <- length(ab)
      logp[Nclust+1] <- dnorm(y[i], 0, sd=sqrt(rgamma(1, ab[M-1], ab[M])), log = TRUE) +
        log(alpha) - log(nrow(X)-1+alpha)
      floor <- min(logp[is.finite(logp)],1e-9)
      logp <- ifelse(!is.finite(logp), floor, logp)
      
      # transform unnormalized log probabilities into probabilities
      loc_probs <- exp(logp - max(logp))
      loc_probs[loc_probs==0] <- min(1e-9,min(loc_probs[loc_probs>0]))
      # set probability of full clusters to zero
      loc_probs[c(n_k,0) > maxSize] <- 0
      loc_probs <- loc_probs / sum(loc_probs)
      # draw a sample of which cluster this data point should belong to
      newz <- sample(1:(Nclust+1), 1, replace = TRUE, prob = loc_probs)
      # spawn a new cluster if necessary
      if(newz == Nclust + 1) {
        n_k <- c(n_k, 0)
        Nclust <- Nclust + 1
      }
      z[i] <- newz
      if (length(Qlist) >= newz){
        parms <- drop(as.matrix(res_gp[newz,]))
        cmatnew <- covmat(D=D[,,1:d,drop=FALSE], set1 = which(z==newz), set2 = which(z==newz),
                          ls=parms[-length(parms)], nug=parms[length(parms)])
        Qlist[[newz]] <- solve(cmatnew + diag(1e-9, nrow = nrow(cmatnew),
                                              ncol = nrow(cmatnew)))
      }
      n_k[newz] <- n_k[newz] + 1 # update the cluster n_k
    }
    
    # Fit a GP to each cluster
    invisible(capture.output(g <- deleteGPseps()))
    Qlist <- list()
    res_gp <- data.frame(matrix(ncol = d+1))  
    for (i in 1:Nclust) {
      dat_i <- X[z==i,,drop=FALSE]
      response_i <- y[z==i]
      gp_i <- newGPsep(dat_i, response_i, d=7.84, g=0.081, dK = TRUE)
      mle_i <- tryCatch({
        mleGPsep(gp_i, param = "both", ab = ab, tmax = c(500,500))
      }, error = function(msg){
        out <- c()
        for (M in 1:(length(ab)/2)) {
          out <- c(out, rgamma(1, ab[2*M-1], ab[2*M]))
        }
        out
      })
      if (class(mle_i) == "list"){
        cvm <- covmat(D=D[,,1:d,drop=FALSE], set1 = which(z==i), set2 = which(z==i),
                      ls=mle_i$theta[-length(mle_i$theta)], 
                      nug = mle_i$theta[length(mle_i$theta)])
        Qlist[[i]] <- solve(cvm + diag(1e-9, nrow = sum(z==i), ncol = sum(z==i)))
        res_gp[i,] <- as.vector(mle_i$theta)
      } else {
        # return draws from prior of theta
        Qlist[[i]] <- diag(mle_i[length(mle_i)], nrow = length(response_i))
        res_gp[i,] <- mle_i
      }
    }
    
    # Sample from posterior of y
    if(!is.null(Xpred) & j%%draw == 0){
      C <- exp(-rowSums((Dpred^2)/rep(phi, each=nrow(Xpred)*N), dims = 2))
      predsbyc <- matrix(nrow=nrow(C), ncol = Nclust)
      uncerbyc <- matrix(nrow = nrow(C), ncol = Nclust)
      finalpreds <- c()
      bestpreds <- c()
      for(k in 1:Nclust){
        pred <- predGPsep(k-1, Xpred[,1:d,drop=FALSE], lite = TRUE)
        predsbyc[,k] <- pred$mean
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
        finalpreds <- c(finalpreds, sum(predsbyc[i,]*probs))
        bestpreds <- c(bestpreds, predsbyc[i,which.max(probs)])
      }
      res_drawmu[,j/draw] <- finalpreds
      colnames(res_drawmu)[j/draw] <- paste('draw',j/draw,sep = '')
      res_pred[,j/draw] <- bestpreds
      colnames(res_pred)[j/draw] <- paste('draw',j/draw,sep = '')
    }
    
    # Sample alpha
    lalpha <- function(x){ out <- logpost_alpha(alpha = 2*x, N=N, k=Nclust)
    out <- ifelse(!is.finite(out), -1000, out)
    return(out)}
    pseu <- list(ld = function(x) dgamma(x,1,1,log = TRUE), 
                 q = function(x) qgamma(x,1,1))
    setTimeLimit(elapsed = 5, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
    proposal <- tryCatch({ slice_quantile(alpha, lalpha, pseudo = pseu)$x}, 
                         error=function(msg){ alpha })
    alpha <- ifelse(is.finite(proposal), proposal/2, alpha)
    setTimeLimit(elapsed = Inf, transient = FALSE)
    
    # Sample phi
    prop_cmat <- hessian(logpost_phi, x=phi, alpha=alpha, z=z, D=D, lp=FALSE)
    prop_cmat <- ifelse(is.finite(prop_cmat), prop_cmat, 0)
    sdphi <- (2.38^2/(d+dq))*abs(solve(prop_cmat + diag(1e-9, d+dq, d+dq)))
    for(k in 1:ncol(sdphi)){
      sdphi[k,k] <- min(maxranges[k], sdphi[k,k])
      sdphi[k,k] <- max(1e-5, sdphi[k,k])
    }
    sdphi <- ifelse(is.na(sdphi), 1, sdphi)
    proposal <- rep(-1, d+dq)
    while(!(all(proposal > 0))){proposal <- mvtnorm::rmvnorm(1, phi, as.matrix(sdphi))}
    aprob <- logpost_phi(proposal, alpha=alpha, z = z, D = D) -
      logpost_phi(phi = phi, alpha = alpha, z = z, D = D)
    u <- runif(1)
    aprob <- ifelse(is.nan(aprob), -Inf, ifelse(is.infinite(aprob),-Inf,aprob))
    if(u < exp(aprob)){ phi = proposal }
    
    # Record results
    res_gatepar[j+1,] <- c(alpha, phi)
    res_z[,j+1] <- z
    colnames(res_z)[j+1] <- paste('z',j,sep = '')
    allnugs <- c()
    allvar <- c()
    for(k in 1:N){
      cnum <- z[k]
      allnugs <- c(allnugs, res_gp[cnum, ncol(res_gp)])
      allvar <- c(allvar, drop(t(y[z==cnum])%*%(Qlist[[cnum]])%*%y[z==cnum])/sum(z==cnum))
    }
    res_noise[,j+1] <- allnugs
    res_spvar[,j+1] <- allvar
    colnames(res_noise)[j+1] <- paste('nug',j,sep = '')
    colnames(res_spvar)[j+1] <- paste('spv',j,sep='')
    
    if(j%%10 == 0){
      print(paste("Iter",j,"done."))
    }
  }
  
  # Return the results
  spv <- c()
  for (i in 1:nrow(res_gp)) {
    tau2 <- drop(t(y[z==i])%*%(Qlist[[i]])%*%y[z==i])/sum(z==i)
    spv <- c(spv, tau2)
  }
  res_gp$spvar <- spv
  res_posterior <- 'none'
  if (!is.null(Xpred)){
    #res_posterior <- list('mu'=res_drawmu, 'lower'=res_drawlo, 'upper'=res_drawhi)
    res_posterior <- list('mean'=res_drawmu, 'best_ests'=res_pred)
  }
  return(list('z'=res_z, 'gatepar'=res_gatepar, 'expertprior'=ab, 'spvar'=res_spvar,
              'noise'=res_noise, 'experts'=res_gp, 'draws'=res_posterior))
}

# IMGPE with GP fitting w/STAN
imgpe.stan <- function(X, Xq=NULL, y, parms=NULL, z_init=NULL, draw = 10, Xpred=NULL,
                       maxSize=nrow(X), maxIters=5000){
  require(brms)
  #require(future)
  #require(future.apply)
  require(qslice)
  require(numDeriv)
  
  N <- nrow(X)
  alpha <- ifelse(!is.null(parms), parms[1], 0.01*N)
  if (!is.null(parms)){
    phi <- parms[-1]
  } else {
    phi <- rep(1, ncol(X))
  }
  ab <- c(1,1,1,1)
  z <- z_init
  if(is.null(z)){z <- rep(1,N)}
  n_k <- as.vector(table(z))
  Nclust <- length(n_k)
  D <- array(0, dim = c(N,N,ncol(X)))
  for (i in 1:ncol(X)) {
    D[,,i] <- as.matrix(dist(X[,i], diag = T, upper = T))
  }
  Qlist <- list()
  maxranges <- apply(X, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  
  # prepare results
  res_gatepar <- data.frame(matrix(c(alpha, phi), nrow=1))
  colnames(res_gatepar) <- c("alpha", colnames(X), colnames(Xq))
  res_z <- data.frame(z0=z)
  res_gp <- data.frame()
  res_noise <- data.frame('nug0'=rep(0,N))
  res_spvar <- data.frame('spv0'=rep(0,N))
  if (!is.null(Xpred)){
    res_drawmu <- data.frame('start'=rep(0,nrow(Xpred)))
    res_drawbest <- data.frame('start'=rep(0,nrow(Xpred)))
    #res_drawlo <- data.frame('start'=rep(0,nrow(Xpred)))
    #res_drawhi <- data.frame('start'=rep(0,nrow(Xpred)))
    Dpred <- array(0, dim = c(nrow(Xpred),N,ncol(X)))
    for (i in 1:ncol(X)) {
      Dpred[,,i] <- as.matrix(laGP::distance(Xpred[,i], X[,i]))
    }
  }
  
  # prepare parallel model fitting
  datdf <- data.frame('x'=as.numeric(X), 'y'=y)
  suppressMessages(suppressWarnings({
    precomp <- brm(y~gp(x, scale=FALSE), data = datdf, 
                   prior = c(prior(constant(0), class=Intercept), 
                             prior(inv_gamma(1,1), class=lscale, coef=gpx), 
                             prior(inv_gamma(1,1), class=sdgp)), 
                   chains = 1, iter = 500, silent = 2, refresh = 0)}))
  dflist <- split(datdf, z)
  #plan(multisession(workers = 8))
  fit_fn <- function(df){update(precomp, newdat=df)}
  
  for (j in 1:maxIters) {
    for (i in 1:nrow(X)) {
      # Gibbs sampling pass over indicator variables z
      z_i <- z[i] # what is the nth persons table assignment?
      n_k[z_i] <- n_k[z_i] - 1 # remove the nth person from table
      # update the matching precision matrix
      if (j>=2){
        j_index <- which((1:N)[z==z_i]==i)
        Qlist[[z_i]] <- inverse.update(Qlist[[z_i]], j_index)
        if (n_k[z_i] == 0){
          Qlist <- Qlist[-z_i]
          res_gp <- res_gp[-z_i,]
        }
      }
      # if the table becomes empty when the nth person is removed,
      # then that table/cluster is removed.
      if( n_k[z_i] == 0 ) {
        n_k <- n_k[-z_i]    # take out the empty cluster
        z[z>=z_i] <- z[z>=z_i]-1 # move up all other clusters
        Nclust <- Nclust - 1     # decrease total number of clusters by 1
      }
      z[i] <- -1
      
      logp <- rep(NA, Nclust+1)
      C <- exp(-rowSums((D^2)/rep(phi, each=N^2), dims = 2))
      for (k in 1:Nclust) {
        if (length(Qlist) >= k){
          z_k <- which(z==k)
          if (length(z_k) != nrow(Qlist[[k]])){
            print(paste('k:',k))
            print(table(z))
            print(nrow(Qlist[[k]]))
            break
          }
          pars <- as.numeric(res_gp[1,-ncol(res_gp)])
          Qy <- drop(covmat(D=D, set1 = i, set2 = z_k,
                            ls = pars))
          #print(paste('zk:',length(z_k),'Qy:',length(Qy),'Qlist:',nrow(Qlist[[k]])))
          mu <- t(Qy)%*%(Qlist[[k]])%*%y[z_k]
          spvar <- drop(t(y[z_k])%*%(Qlist[[k]])%*%y[z_k])/length(z_k)
          noise <- res_gp[k,ncol(res_gp)]
          sig2 <- spvar + noise - t(Qy)%*%(Qlist[[k]])%*%Qy
          sig2 <- ifelse(sig2<=0, 1e-4, sig2)
          lcondprob <- dnorm(y[i], mean = mu, sd = sqrt(sig2), log = TRUE)
        } else {
          nug <- 0
          for (M in 1:(length(ab)/2)) {
            nug <- nug + (1/rgamma(1, ab[2*M-1], ab[2*M]))
          }
          lcondprob <- dnorm(y[i], mean = 0, sd = sqrt(nug), log = TRUE)
        }
        logp[k] <- lcondprob + logprob_zi(i, j=k, z=z, C=C, alpha=alpha)
      }
      
      M <- length(ab)
      logp[Nclust+1] <- dnorm(y[i], sd=sqrt(1/rgamma(1, ab[M-1], ab[M])), log = TRUE) +
        log(alpha) - log(nrow(X)-1+alpha)
      floor <- min(logp[is.finite(logp)],1e-9)
      logp <- ifelse(!is.finite(logp), floor, logp)
      
      # transform unnormalized log probabilities into probabilities
      loc_probs <- exp(logp - max(logp))
      loc_probs[loc_probs==0] <- min(1e-9,min(loc_probs[loc_probs>0]))
      # set probability of full clusters to zero
      loc_probs[c(n_k,0) > maxSize] <- 0
      loc_probs <- loc_probs / sum(loc_probs)
      # draw a sample of which cluster this data point should belong to
      newz <- sample(1:(Nclust+1), 1, replace = TRUE, prob = loc_probs)
      # spawn a new cluster if necessary
      if(newz == Nclust + 1) {
        n_k <- c(n_k, 0)
        Nclust <- Nclust + 1
      }
      z[i] <- newz
      if (length(Qlist) >= newz){
        parms <- drop(as.matrix(res_gp[newz,]))
        cmatnew <- covmat(D=D, set1 = which(z==newz), set2 = which(z==newz),
                          ls=parms[-length(parms)], nug=parms[length(parms)])
        Qlist[[newz]] <- solve(cmatnew + diag(1e-9, nrow = nrow(cmatnew),
                                              ncol = nrow(cmatnew)))
      }
      n_k[newz] <- n_k[newz] + 1 # update the cluster n_k
    }
    
    # Fit a GP to each cluster
    Qlist <- list()
    res_gp <- data.frame(matrix(ncol = length(phi)+1))
    newmods <- list()
    #dflist <- split(datdf, z)
    #suppressWarnings({dfsplit <- split(dflist, rep(1:5, each=12))})
    #for (B in 1:length(dflist)) {
    #  suppressMessages(suppressWarnings({
    #newmodset <- future_lapply(dfsplit[[B]], fit_fn, future.seed = TRUE)
    #    newmodset <- fit_fn(dflist[[B]])
    #  }))
    #  newmods <- append(newmods, newmodset)
    #}
    
    for (i in 1:Nclust) {
      zw <- which(z==i)
      newmods[[i]] <- update(precomp, newdata = datdf[zw,])
      pars_all <- posterior_summary(newmods[[i]])
      res_gp[i,] <- pars_all[2:3,1]
      Qlist[[i]] <- solve(covmat(D=D[zw,zw,,drop=FALSE], ls=pars_all[2,1], 
                                 nug = pars_all[3,1]))
    }
    
    # Sample from posterior of y
    if(!is.null(Xpred) & j%%draw == 0){
      C <- exp(-rowSums((Dpred^2)/rep(phi, each=nrow(Xpred)*N), dims = 2))
      predsbyc <- matrix(nrow=nrow(Xpred), ncol = Nclust)
      finalpreds <- c()
      bestpreds <- c()
      for(k in 1:Nclust){
        pred <- fitted(newmods[[k]], data.frame('x'=as.numeric(Xpred)))
        predsbyc[,k] <- pred[,1]
      }
      for(i in 1:nrow(Xpred)){
        lpz <- function(x){#logprob_zi(zi=i,j=x,z=z,C=C,alpha=alpha)
          n <- length(z)
          lp <- log(n-1) + log(sum(C[i, which(z==x)])) -
            log(sum(C[i,])) - log(n-1+alpha)
          return(lp)
        }
        logp <- sapply(1:Nclust, lpz)
        probs <- exp(logp)/sum(exp(logp))
        finalpreds <- c(finalpreds, sum(predsbyc[i,]*probs))
        bestpreds <- c(bestpreds, predsbyc[i,which.max(probs)])
      }
      res_drawmu[,j/draw] <- finalpreds
      colnames(res_drawmu)[j/draw] <- paste('draw',j/draw,sep = '')
      res_drawbest[,j/draw] <- bestpreds
      colnames(res_drawbest)[j/draw] <- paste('draw',j/draw,sep = '')
    }
    
    # Sample alpha
    lalpha <- function(x){ out <- logpost_alpha(alpha = x, N=N, k=Nclust)
    out <- ifelse(!is.finite(out), -1000, out)
    return(out)}
    pseu <- list(ld = function(x) dgamma(x,1,1,log = TRUE), 
                 q = function(x) qgamma(x,1,1))
    setTimeLimit(elapsed = 5, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
    proposal <- tryCatch({ slice_quantile(alpha, lalpha, pseudo = pseu)$x}, 
                         error=function(msg){ alpha })
    alpha <- ifelse(is.finite(proposal), proposal, alpha)
    setTimeLimit(elapsed = Inf, transient = TRUE)
    
    # Sample phi
    prop_cmat <- hessian(logpost_phi, x=phi, alpha=alpha, z=z, D=D, lp=FALSE)
    sdphi <- (2.38^2/length(phi))*abs(solve(prop_cmat + diag(1e-9, length(phi), length(phi))))
    for(k in 1:ncol(sdphi)){
      sdphi[k,k] <- min(maxranges[k], sdphi[k,k])
      sdphi[k,k] <- max(1e-5, sdphi[k,k])
    }
    proposal <- rep(-1, length(phi))
    while(!(all(proposal > 0))){proposal <- mvtnorm::rmvnorm(1, phi, as.matrix(sdphi))}
    aprob <- logpost_phi(proposal, alpha=alpha, z = z, D = D) -
      logpost_phi(phi = phi, alpha = alpha, z = z, D = D)
    u <- runif(1)
    aprob <- ifelse(is.nan(aprob), -Inf, ifelse(is.infinite(aprob),-Inf,aprob))
    if(u < exp(aprob)){ phi = proposal }
    
    # Record results
    res_gatepar[j+1,] <- c(alpha, phi)
    res_z[,j+1] <- z
    colnames(res_z)[j+1] <- paste('z',j,sep = '')
    allnugs <- c()
    allvar <- c()
    for(k in 1:N){
      cnum <- z[k]
      allnugs <- c(allnugs, res_gp[z[k], ncol(res_gp)])
      allvar <- c(allvar, drop(t(y[z==cnum])%*%(Qlist[[cnum]])%*%y[z==cnum])/sum(z==cnum))
    }
    res_noise[,j+1] <- allnugs
    res_spvar[,j+1] <- allvar
    colnames(res_noise)[j+1] <- paste('nug',j,sep = '')
    colnames(res_spvar)[j+1] <- paste('spv',j,sep='')
    
    #if(j%%1 == 0){
    print(paste("Iter",j,"done."))
    #}
  }
  
  # Return the results
  spv <- c()
  for (i in 1:nrow(res_gp)) {
    tau2 <- drop(t(y[z==i])%*%(Qlist[[i]])%*%y[z==i])/sum(z==i)
    spv <- c(spv, tau2)
  }
  res_gp$spvar <- spv
  res_posterior <- 'none'
  if (!is.null(Xpred)){
    #res_posterior <- list('mu'=res_drawmu, 'lower'=res_drawlo, 'upper'=res_drawhi)
    res_posterior <- list('avg'=res_drawmu, 'best'=res_drawbest)
  }
  return(list('z'=res_z, 'gatepar'=res_gatepar, 'expertprior'=ab, 'spvar'=res_spvar,
              'noise'=res_noise, 'experts'=res_gp, 'draws'=res_posterior))
}

# IMGPE with GP fitting w/custom slice sampler
imgpe.slice <- function(X, Xq=NULL, y, parms=NULL, z_init=NULL, draw = 10, Xpred=NULL,
                       maxSize=nrow(X), maxIters=5000){
  require(future)
  require(future.apply)
  require(qslice)
  require(numDeriv)
  require(invgamma)
  #require(lognorm)
  
  N <- nrow(X)
  ndim <- ncol(X)
  alpha <- ifelse(!is.null(parms), parms[1], 0.01*N)
  if (!is.null(parms)){
    phi <- parms[-1]
  } else {
    phi <- rep(1, ncol(X))
  }
  ab <- rep(1, (ndim+1)*2)
  z <- z_init
  if(is.null(z)){z <- rep(1,N)}
  n_k <- as.vector(table(z))
  Nclust <- length(n_k)
  D <- array(0, dim = c(N,N,ncol(X)))
  for (i in 1:ncol(X)) {
    D[,,i] <- as.matrix(dist(X[,i], diag = TRUE, upper = TRUE))
  }
  Qlist <- list()
  maxranges <- apply(X, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  llgp <- function(D, y, theta, spvar=1, nug){
    n <- length(y)
    mu <- mean(y)
    if (theta >0 & nug >= 0){
      Sig <- covmat(D=D, ls=theta, nug = nug)*spvar
      dsig <- as.numeric(determinant(Sig)[1])
      prtheta <- sum(dinvgamma(theta, 1, 1, log = TRUE))
      prnug <- ifelse(nug==0, 0, dinvgamma(nug, 1, 1, log = TRUE))
      out <- (n/2)*log(2*pi) + (1/2)*log(dsig) + 0.5*t(y-mu)%*%solve(Sig)%*%(y-mu) - prtheta - prnug
    } else {
      out <- 0
    }
    return(-out)
  }
  SigM <- matrix(c(2,0,0,2), nrow = 2)
  
  # prepare results
  res_gatepar <- data.frame(matrix(c(alpha, phi), nrow=1))
  colnames(res_gatepar) <- c("alpha", colnames(X), colnames(Xq))
  res_z <- data.frame(z0=z)
  res_gp <- data.frame()
  res_noise <- data.frame('nug0'=rep(0,N))
  #res_spvar <- data.frame('spv0'=rep(0,N))
  if (!is.null(Xpred)){
    Dpred <- array(0, dim = c(nrow(Xpred),nrow(Xpred),ncol(Xpred)))
    for (i in 1:ncol(Xpred)) {
      Dpred[,,i] <- as.matrix(dist(Xpred[,i], diag = TRUE, upper = TRUE))
    }
    dd <- array(0, dim = c(nrow(Xpred),N,ncol(X)))
    for (i in 1:ncol(X)) {
      dd[,,i] <- sqrt(laGP::distance(Xpred, X))
    }
    res_drawmu <- data.frame('start'=rep(0,nrow(Xpred)))
    res_drawbest <- data.frame('start'=rep(0,nrow(Xpred)))
  }
  finaltau <- c()
  
  # prepare parallel model fitting
  plan(multisession(workers = 8))
  
  for (j in 1:maxIters) {
    for (i in 1:nrow(X)) {
      # Gibbs sampling pass over indicator variables z
      z_i <- z[i] # what is the nth persons table assignment?
      n_k[z_i] <- n_k[z_i] - 1 # remove the nth person from table
      # update the matching precision matrix
      if (j>=2){
        j_index <- which((1:N)[z==z_i]==i)
        Qlist[[z_i]] <- inverse.update(Qlist[[z_i]], j_index)
        if (n_k[z_i] == 0){
          Qlist <- Qlist[-z_i]
          res_gp <- res_gp[-z_i,]
        }
      }
      # if the table becomes empty when the nth person is removed,
      # then that table/cluster is removed.
      if( n_k[z_i] == 0 ) {
        n_k <- n_k[-z_i]    # take out the empty cluster
        z[z>=z_i] <- z[z>=z_i]-1 # move up all other clusters
        Nclust <- Nclust - 1     # decrease total number of clusters by 1
      }
      z[i] <- -1
      
      logp <- rep(NA, Nclust+1)
      C <- exp(-rowSums((D^2)/rep(phi, each=N^2), dims = 2))
      for (k in 1:Nclust) {
        if (length(Qlist) >= k){
          z_k <- which(z==k)
          pars <- as.numeric(res_gp[1,-ncol(res_gp)])
          Qy <- drop(covmat(D=D, set1 = i, set2 = z_k,
                            ls = pars))
          #print(paste('zk:',length(z_k),'Qy:',length(Qy),'Qlist:',nrow(Qlist[[k]])))
          mu <- t(Qy)%*%(Qlist[[k]])%*%y[z_k]
          spvar <- drop(t(y[z_k])%*%(Qlist[[k]])%*%y[z_k])/length(z_k)
          noise <- res_gp[k,ncol(res_gp)]
          sig2 <- spvar + noise - t(Qy)%*%(Qlist[[k]])%*%Qy
          sig2 <- ifelse(sig2<=0, 1e-4, sig2)
          lcondprob <- dnorm(y[i], mean = mu, sd = sqrt(sig2), log = TRUE)
        } else {
          nug <- 0
          for (M in 1:(length(ab)/2)) {
            nug <- nug + rinvgamma(1, ab[2*M-1], ab[2*M])
          }
          lcondprob <- dnorm(y[i], mean = 0, sd = sqrt(nug), log = TRUE)
        }
        logp[k] <- lcondprob + logprob_zi(i, j=k, z=z, C=C, alpha=alpha)
      }
      
      M <- length(ab)
      logp[Nclust+1] <- dnorm(y[i], 0, sd=sqrt(rinvgamma(1, ab[M-1], ab[M])), log = TRUE) +
        log(alpha) - log(nrow(X)-1+alpha)
      floor <- min(logp[is.finite(logp)],1e-9)
      logp <- ifelse(!is.finite(logp), floor, logp)
      
      # transform unnormalized log probabilities into probabilities
      loc_probs <- exp(logp - max(logp))
      loc_probs[loc_probs==0] <- min(1e-9,min(loc_probs[loc_probs>0]))
      # set probability of full clusters to zero
      loc_probs[c(n_k,0) > maxSize] <- 0
      loc_probs <- loc_probs / sum(loc_probs)
      # draw a sample of which cluster this data point should belong to
      newz <- sample(1:(Nclust+1), 1, replace = TRUE, prob = loc_probs)
      # spawn a new cluster if necessary
      if(newz == Nclust + 1) {
        n_k <- c(n_k, 0)
        Nclust <- Nclust + 1
      }
      z[i] <- newz
      if (length(Qlist) >= newz){
        parms <- drop(as.matrix(res_gp[newz,]))
        cmatnew <- covmat(D=D, set1 = which(z==newz), set2 = which(z==newz),
                          ls=parms[-length(parms)], nug=parms[length(parms)])
        Qlist[[newz]] <- solve(cmatnew + diag(1e-9, nrow = nrow(cmatnew),
                                              ncol = nrow(cmatnew)))
      }
      n_k[newz] <- n_k[newz] + 1 # update the cluster n_k
    }
    
    # Fit a GP to each cluster
    Qlist <- list()
    res_gp <- data.frame(matrix(ncol = length(phi)+1))
    res_noise[,j+1] <- 0
    
    indexlist <- split(seq_along(z), z)
    fitslice <- function(zw){
      myllgp <- function(x){llgp(D=D[zw,zw,,drop=FALSE], y=y[zw], theta=x[1], nug=x[2])}
      draws <- matrix(c(1,0.1), nrow = 1000, ncol = 2, byrow = T)
      for (k in 2:1000) {
        setTimeLimit(elapsed = 0.1, transient = TRUE)
        on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
        
        tryCatch({slice_elliptical_mv(x = draws[k-1,], log_target=myllgp, 
                                      mu=c(2,1), Sig = SigM)$x},
                 error=function(msg){ print(msg) }) 
        drawk <- tryCatch({slice_elliptical_mv(x = draws[k-1,], log_target=myllgp, 
                                               mu=c(2,1), Sig = SigM)$x},
                          error=function(msg){ draws[k-1,] }) 
        draws[k,] <- ifelse(drawk>0, drawk, 0.001)
      }
      setTimeLimit(elapsed = Inf, transient = TRUE)
      return(colMeans(draws[501:1000,]))
    }
    newGPparms <- future_lapply(indexlist, fitslice, future.seed = TRUE)
    
    for (i in 1:Nclust) {
      zw <- indexlist[[i]]
      res_noise[zw,j+1] <- newGPparms[[i]][ndim+1]
      res_gp[i,] <- newGPparms[[i]]
      Qlist[[i]] <- solve(covmat(D=D[zw,zw,,drop=FALSE], 
                                 ls=newGPparms[[i]][-(ndim+1)], 
                                 nug = newGPparms[[i]][ndim+1]))
    }
    colnames(res_noise)[j+1] <- paste('nug',j,sep = '')
    
    # Sample from posterior of y
    if(!is.null(Xpred) & j%%draw == 0){
      C <- exp(-rowSums((Dpred^2)/rep(phi, each=nrow(Xpred)^2), dims = 2))
      predsbyc <- matrix(nrow=nrow(C), ncol = Nclust)
      finalpreds <- c()
      bestpreds <- c()
      for(k in 1:Nclust){
        z_k <- which(z==k)
        pars <- as.numeric(res_gp[k,])
        Qy <- drop(covmat(D=dd, set1 = 1:nrow(dd), set2 = z_k,
                          ls = pars[-length(pars)]))
        pred <- Qy%*%(Qlist[[k]])%*%y[z_k]
        predsbyc[,k] <- as.numeric(pred)
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
        finalpreds <- c(finalpreds, sum(predsbyc[i,]*probs))
        bestpreds <- c(bestpreds, predsbyc[i,which.max(probs)])
      }
      res_drawmu[,j/draw] <- finalpreds
      colnames(res_drawmu)[j/draw] <- paste('draw',j/draw,sep = '')
      res_drawbest[,j/draw] <- bestpreds
      colnames(res_drawbest)[j/draw] <- paste('draw',j/draw,sep = '')
    }
    
    # Sample alpha
    if(FALSE){
      lalpha <- function(x){ out <- logpost_alpha(alpha = 2*x, N=N, k=Nclust)
      out <- ifelse(!is.finite(out), -1000, out)
      return(out)}
      pseu <- list(ld = function(x) dgamma(x,1,1,log = TRUE), 
                   q = function(x) qgamma(x,1,1))
      setTimeLimit(elapsed = 5, transient = TRUE)
      on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
      proposal <- tryCatch({ slice_quantile(alpha, lalpha, pseudo = pseu)$x}, 
                           error=function(msg){ alpha })
      alpha <- ifelse(is.finite(proposal), proposal, alpha)
      setTimeLimit(elapsed = Inf, transient = TRUE)
    }
    
    # Sample phi
    prop_cmat <- hessian(logpost_phi, x=phi, alpha=alpha, z=z, D=D, lp=FALSE)
    sdphi <- (2.38^2/length(phi))*abs(solve(prop_cmat + diag(1e-9, length(phi), length(phi))))
    for(k in 1:ncol(sdphi)){
      sdphi[k,k] <- min(maxranges[k], sdphi[k,k])
      sdphi[k,k] <- max(1e-5, sdphi[k,k])
    }
    proposal <- rep(-1, length(phi))
    while(!(all(proposal > 0))){proposal <- mvtnorm::rmvnorm(1, phi, as.matrix(sdphi))}
    aprob <- logpost_phi(proposal, alpha=alpha, z = z, D = D) -
      logpost_phi(phi = phi, alpha = alpha, z = z, D = D)
    u <- runif(1)
    aprob <- ifelse(is.nan(aprob), -Inf, ifelse(is.infinite(aprob),-Inf,aprob))
    if(u < exp(aprob)){ phi = proposal }
    
    # Record results
    res_gatepar[j+1,] <- c(alpha, phi)
    res_z[,j+1] <- z
    colnames(res_z)[j+1] <- paste('z',j,sep = '')
    
    if(j%%10 == 0){
      print(paste("Iter",j,"done."))
    }
  }
  
  # Return the results
  spv <- c()
  for (i in 1:nrow(res_gp)) {
    tau2 <- drop(t(y[z==i])%*%(Qlist[[i]])%*%y[z==i])/sum(z==i)
    spv <- c(spv, tau2)
  }
  res_gp$spvar <- spv
  res_posterior <- 'none'
  if (!is.null(Xpred)){
    res_posterior <- list('avg'=res_drawmu, 'best'=res_drawbest)
  }
  return(list('z'=res_z, 'gatepar'=res_gatepar, 'expertprior'=ab, 
              'noise'=res_noise, 'experts'=res_gp, 'draws'=res_posterior))
}

# IMGPE-DD-EM
imgpe.ddem <- function(X, Xq=NULL, y, parms=NULL, z_init=NULL, draw = 10, Xpred=NULL,
                       maxSize=nrow(X), maxIters=5000){
  require(future)
  require(future.apply)
  library(laGP)
  require(qslice)
  require(numDeriv)
  #require(plyr)
  #require(invgamma)
  
  N <- nrow(X)
  ndim <- ncol(X)
  alpha <- ifelse(!is.null(parms), parms[1], 0.01*N)
  if (!is.null(parms)){
    phi <- parms[-1]
  } else {
    phi <- rep(1, ndim)
  }
  ab <- c(1,1,1,1)
  st <- data.frame(idx=1:N, cluster=ifelse(is.null(z_init), 1:N, z_init), customer=1:N)
  n_k <- as.vector(table(st$cluster))
  Nclust <- length(n_k)
  D <- array(0, dim = c(N,N,ncol(X)))
  for (i in 1:ncol(X)) {
    D[,,i] <- as.matrix(dist(X[,i], diag = TRUE, upper = TRUE))
  }
  maxranges <- apply(X, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  llgp <- function(D, y, theta, nug=0.0001){
    n <- length(y)
    if (all(theta >0) & nug >= 0){
      Sig <- covmat(D=D, ls=theta, nug = nug)
      dsig <- as.numeric(determinant(Sig)[1])
      prtheta <- sum(dgamma(theta, 1, 1, log = TRUE))
      prnug <- ifelse(nug==0, 0, dgamma(nug, 1, 1, log = TRUE))
      out <- (n/2)*log(2*pi) + (n/2)*dsig + 0.5*t(y)%*%solve(Sig)%*%y - prtheta - prnug
    } else {
      out <- 0
    }
    return(-out)
  }
  SigM <- matrix(c(2,0,0,2), nrow = 2)
  lhood <- rep(0,N)
  names(lhood) <- as.character(1:N)
  
  # prepare results
  res_gatepar <- data.frame(matrix(c(alpha, phi), nrow=1))
  colnames(res_gatepar) <- c("alpha", colnames(X), colnames(Xq))
  res_z <- data.frame(z0=1:N)
  res_gp <- data.frame(matrix(c(phi,1), nrow = N, ncol = ndim+1, byrow = T))
  res_noise <- data.frame('nug0'=rep(0,N))
  #res_spvar <- data.frame('spv0'=rep(0,N))
  if (!is.null(Xpred)){
    Dpred <- array(0, dim = c(nrow(Xpred),nrow(Xpred),ncol(Xpred)))
    for (i in 1:ncol(Xpred)) {
      Dpred[,,i] <- as.matrix(dist(Xpred[,i], diag = TRUE, upper = TRUE))
    }
    dd <- array(0, dim = c(nrow(Xpred),N,ncol(X)))
    for (i in 1:ncol(X)) {
      dd[,,i] <- sqrt(laGP::distance(Xpred[,i], X[,i]))
    }
    res_drawmu <- data.frame('start'=rep(0,nrow(Xpred)))
  }
  finaltau <- c()
  
  for (j in 1:maxIters) {
    # compute priors for clustering
    logpriorlist <- list()
    candlinkslist <- list()
    logpriorthresh <- -10
    probconn_one <- c()
    for (i in 2:N) {
      logprior_i <- sapply(1:N, function(j) safelog(exp(-sum((X[i,]-X[j,])^2/phi))))
      logprior_i[i] <- log(alpha)
      logprior_i <- logprior_i - log.sum(logprior_i)
      probconn_one <- c(probconn_one, logprior_i[1])
      candlinks_i <- which(logprior_i > logpriorthresh)
      logpriorlist[[i]] <- logprior_i
      candlinkslist[[i]] <- candlinks_i
    }
    C <- exp(-rowSums((D^2)/rep(phi, each=N^2), dims = 2))
    #if(j==5){print(probconn_one)}
    
    # Gibbs sampling pass over indicator variables z
    for (i in 2:N) {
      # remove point i from the data
      # if the table becomes empty when the nth person is removed,
      # then that table/cluster is removed.
      old.cluster <- st$cluster[i]
      old.customer <- st$customer[i]
      conn.i <- connections(i, st$customer)
      if (!(old.customer %in% conn.i)) st$cluster[conn.i] <- i
      st$customer[i] <- i
      
      # update the likelihoods if a table has been split
      if (old.customer != i & !(old.customer %in% conn.i))
      {
        old.idx <- st[which(st$cluster==old.cluster),"idx"]
        lhood[as.character(old.cluster)] <- llgp(D[old.idx,old.idx,,drop=FALSE], y[old.idx], 
                                                 theta = unlist(res_gp[old.cluster,-(ndim+1)]),
                                                 nug = res_gp[old.cluster, ndim+1])
      }
      else
      {
        lhood[as.character(old.cluster)] <- 0
      }
      
      # get log prior
      
      log.prior <- logpriorlist[[i]]
      cand.links <- candlinkslist[[i]]
      
      # compute the likelihood of data point i (and its connectors)
      # with all other tables
      
      cand.clusts <- unique(st$cluster[cand.links])
      new.lhood <- c()
      for (k in 1:length(cand.clusts)) {
        df <- st[st$cluster == cand.clusts[k],]
        pts <- unique(c(df$idx,st[conn.i,"idx"]))
        new.lhood[k] <- llgp(D[pts,pts,,drop=FALSE], y[pts],
                             theta = unlist(res_gp[cand.clusts[k],-(ndim+1)]),
                             nug = res_gp[cand.clusts[k], ndim+1])
      }
      
      names(new.lhood) <- cand.clusts
      
      # set up the old likelihoods
      
      old.lhood <- lhood[as.character(cand.clusts)]
      sum.old.lhood <- sum(old.lhood)
      
      # compute the conditional distribution
      log.prob <-
        log.prior[cand.links] +
        sapply(cand.links,
               function (j) {
                 c.j <- as.character(st$cluster[j])
                 sum.old.lhood - old.lhood[c.j] + new.lhood[c.j] })
      
      prob <- exp(log.prob - log.sum(log.prob))
      if(sum(prob)==0) {
        print('log.prob:')
        print(log.prob)
      }
      if (length(prob)==1)
        st$customer[i] <- cand.links[1]
      else
        st$customer[i] <- sample(cand.links, 1, prob=prob)
      
      # update the clusters
      st$cluster[conn.i] <- st$cluster[st$customer[i]]
      pts_i <- subset(st, cluster == st$cluster[i])$idx
      clus_i <- st$cluster[i]
      res_gp[pts_i,] <- as.list(res_gp[st$customer[i],])
      lhood[as.character(st$cluster[i])] <- llgp(D[pts_i,pts_i,,drop=FALSE], y[pts_i],
                                                 theta = unlist(res_gp[clus_i,-(ndim+1)]),
                                                 nug = res_gp[clus_i, ndim+1])
    }
    #update Nclust, n_k, z
    z <- st$cluster
    n_k <- as.vector(table(z))
    Nclust <- length(n_k)
    
    # Fit a GP to each cluster
    invisible(capture.output(g <- deleteGPseps()))
    res_gp <- data.frame(matrix(ncol = length(phi)+1))
    res_noise[,j+1] <- 0
    
    indexlist <- split(seq_along(z), z)
    clustlist <- unique(st$cluster)
    for (i in 1:Nclust) {
      ci <- which(z == clustlist[i])
      dat_i <- X[ci,,drop=FALSE]
      response_i <- y[ci]
      gp_i <- newGPsep(dat_i, response_i, d=7.84, g=0.081, dK = TRUE)
      mle_i <- tryCatch({
        mleGPsep(gp_i, param = "both", ab = ab, tmax = c(500,500))
      }, error = function(msg){
        out <- c()
        for (M in 1:(length(ab)/2)) {
          out <- c(out, rinvgamma(1, ab[2*M-1], ab[2*M]))
        }
        out
      })
      if (class(mle_i) == "list"){
        res_gp[ci,] <- as.list(mle_i$theta)
      } else {
        # return draws from prior of theta
        res_gp[ci,] <- as.list(mle_i)
      }
     
      res_noise[ci,j+1] <- res_gp[ci[1],ndim+1]
    }
    colnames(res_noise)[j+1] <- paste('nug',j,sep = '')
    
    # Sample from posterior of y
    if(!is.null(Xpred) & j%%draw == 0){
      CC <- exp(-rowSums((dd^2)/rep(phi, each=nrow(Xpred)*N), dims = 2))
      finalpreds <- c()
      predsbyc <- matrix(nrow=nrow(CC), ncol = Nclust)
      for(k in 1:Nclust){
        pred <- predGPsep(k-1, Xpred, lite = TRUE)
        predsbyc[,k] <- pred$mean
      }
      for (i in 1:nrow(CC)) {
        CCi <- CC[i,]
        #CCi[i] <- log(alpha)
        #CCi <- safelog(CCi) - safelog(sum(CCi))
        probs <- as.numeric(tapply(CCi, z, sum))
        probs <- probs/sum(probs) #exp(lprobs)/sum(exp(lprobs))
        finalpreds <- c(finalpreds, sum(predsbyc[i,]*probs))
      }
      res_drawmu[,j/draw] <- finalpreds
      colnames(res_drawmu)[j/draw] <- paste('draw',j/draw,sep = '')
    }
    
    # Sample alpha
    if (TRUE) {
      lalpha <- function(x){ out <- logpost_alpha(alpha = x, N=N, k=Nclust)
      out <- ifelse(!is.finite(out), -1000, out)
      return(out)}
      pseu <- list(ld = function(x) dgamma(x,1,1,log = TRUE), 
                   q = function(x) qgamma(x,1,1))
      setTimeLimit(elapsed = 5, transient = TRUE)
      on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
      proposal <- tryCatch({ slice_quantile(alpha, lalpha, pseudo = pseu)$x}, 
                           error=function(msg){ alpha })
      alpha <- ifelse(is.finite(proposal), proposal, alpha)
      setTimeLimit(elapsed = Inf, transient = TRUE)
    }
    
    # Sample phi
    logpost_phi <- function(phi, D=D, z=z, alpha=alpha, Xq=NULL, lp=TRUE){
      n <- length(z)
      C <- exp(-rowSums((D^2)/rep(phi, each=n^2), dims = 2))
      step1 <- 0
      for (i in 1:nrow(C)) {
        Ci <- C[i,]
        Ci[i] <- log(alpha)
        Ci <- safelog(Ci) - safelog(sum(Ci))
        lprobs <- as.numeric(tapply(Ci, z, sum))
        probs <- exp(lprobs)/sum(exp(lprobs))
        step1 <- step1 + log(probs[which(unique(z)==z[i])])
      }
      step2 <- max(-1000, log(phi))
      step3 <- max(-1000, log(phi*sqrt(2*pi)))
      step4 <- step1 - 0.5*sum(step2^2) - step3
      if(!lp){ step4 <- exp(step4) }
      return(step4)
    }
    
    if (TRUE) {
      prop_cmat <- hessian(logpost_phi, x=phi, D=D, z=z, alpha=alpha, lp=FALSE)
      prop_cmat <- ifelse(!is.finite(prop_cmat), 0, prop_cmat)
      sdphi <- (2.38^2/length(phi))*abs(solve(prop_cmat + diag(1e-9, length(phi), length(phi))))
      for(k in 1:ncol(sdphi)){
        sdphi[k,k] <- min(maxranges[k], sdphi[k,k])
        sdphi[k,k] <- max(1e-5, sdphi[k,k])
      }
      proposal <- rep(-1, length(phi))
      while(!(all(proposal > 0))){proposal <- mvtnorm::rmvnorm(1, phi, as.matrix(sdphi))}
      aprob <- logpost_phi(proposal, alpha=alpha, z = z, D = D) -
        logpost_phi(phi = phi, alpha = alpha, z = z, D = D)
      u <- runif(1)
      aprob <- ifelse(is.nan(aprob), -Inf, ifelse(is.infinite(aprob),-Inf,aprob))
      if(u < exp(aprob)){ phi = proposal }
    }
    
    # Record results
    res_gatepar[j+1,] <- c(alpha, phi)
    res_z[,j+1] <- z
    colnames(res_z)[j+1] <- paste('z',j,sep = '')
    
    if(j%%10 == 0){
      print(paste("Iter",j,"done."))
    }
  }
  
  # Return the results
  spv <- c()
  clusts <- unique(z)
  for (i in clusts) {
    zi <- which(z==i)
    Kinv <- solve(covmat(D[zi,zi,,drop=FALSE], ls=unlist(res_gp[i,-ncol(res_gp)]),
                         nug = res_gp[i,ncol(res_gp)]))
    tau2 <- drop(t(y[z==i])%*%Kinv%*%y[z==i])/sum(z==i)
    spv <- c(spv, tau2)
  }
  res_gp$spvar <- spv
  res_posterior <- 'none'
  if (!is.null(Xpred)){
    res_posterior <- res_drawmu
  }
  return(list('z'=res_z, 'gatepar'=res_gatepar, 'expertprior'=ab, 
              'noise'=res_noise, 'experts'=res_gp, 'draws'=res_posterior))
}

# IMGPE-DD-Slice Sampler
imgpe.ddslc <- function(X, Xq=NULL, y, parms=NULL, z_init=NULL, draw = 10, Xpred=NULL,
                       maxSize=nrow(X), maxIters=5000){
  require(future)
  require(future.apply)
  library(laGP)
  require(qslice)
  require(numDeriv)
  #require(plyr)
  #require(invgamma)
  
  N <- nrow(X)
  ndim <- ncol(X)
  alpha <- ifelse(!is.null(parms), parms[1], 0.01*N)
  if (!is.null(parms)){
    phi <- parms[-1]
  } else {
    phi <- rep(1, ndim)
  }
  ab <- c(1,1,1,1)
  st <- data.frame(idx=1:N, cluster=ifelse(is.null(z_init), 1:N, z_init), customer=1:N)
  n_k <- as.vector(table(st$cluster))
  Nclust <- length(n_k)
  D <- array(0, dim = c(N,N,ncol(X)))
  for (i in 1:ncol(X)) {
    D[,,i] <- as.matrix(dist(X[,i], diag = TRUE, upper = TRUE))
  }
  maxranges <- apply(X, 2, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  llgp <- function(D, y, theta, nug=0){
    n <- length(y)
    if (theta >0 & nug >= 0){
      Sig <- covmat(D=D, ls=theta, nug = nug)
      dsig <- as.numeric(determinant(Sig)[1])
      prtheta <- sum(dgamma(theta, 1, 1, log = TRUE))
      prnug <- ifelse(nug==0, 0, dgamma(nug, 1, 1, log = TRUE))
      out <- (n/2)*log(2*pi) + (n/2)*dsig + 0.5*t(y)%*%solve(Sig)%*%y - prtheta - prnug
    } else {
      out <- 0
    }
    return(-out)
  }
  SigM <- matrix(c(2,0,0,2), nrow = 2)
  lhood <- rep(0,N)
  names(lhood) <- as.character(1:N)
  
  # prepare results
  res_gatepar <- data.frame(matrix(c(alpha, phi), nrow=1))
  colnames(res_gatepar) <- c("alpha", colnames(X), colnames(Xq))
  res_z <- data.frame(z0=1:N)
  res_gp <- data.frame(matrix(c(phi,1), nrow = N, ncol = ndim+1, byrow = T))
  res_noise <- data.frame('nug0'=rep(0,N))
  #res_spvar <- data.frame('spv0'=rep(0,N))
  if (!is.null(Xpred)){
    Dpred <- array(0, dim = c(nrow(Xpred),nrow(Xpred),ncol(Xpred)))
    for (i in 1:ncol(Xpred)) {
      Dpred[,,i] <- as.matrix(dist(Xpred[,i], diag = TRUE, upper = TRUE))
    }
    dd <- array(0, dim = c(nrow(Xpred),N,ncol(X)))
    for (i in 1:ncol(X)) {
      dd[,,i] <- sqrt(laGP::distance(Xpred, X))
    }
    res_drawmu <- data.frame('start'=rep(0,nrow(Xpred)))
  }
  finaltau <- c()
  
  for (j in 1:maxIters) {
    # compute priors for clustering
    logpriorlist <- list()
    candlinkslist <- list()
    logpriorthresh <- -10
    probconn_one <- c()
    for (i in 2:N) {
      logprior_i <- sapply(1:N, function(j) safelog(exp(-sum((X[i,]-X[j,])^2/phi))))
      logprior_i[i] <- log(alpha)
      logprior_i <- logprior_i - log.sum(logprior_i)
      probconn_one <- c(probconn_one, logprior_i[1])
      candlinks_i <- which(logprior_i > logpriorthresh)
      logpriorlist[[i]] <- logprior_i
      candlinkslist[[i]] <- candlinks_i
    }
    C <- exp(-rowSums((D^2)/rep(phi, each=N^2), dims = 2))
    #if(j==5){print(probconn_one)}
    
    # Gibbs sampling pass over indicator variables z
    for (i in 2:N) {
      # remove point i from the data
      # if the table becomes empty when the nth person is removed,
      # then that table/cluster is removed.
      old.cluster <- st$cluster[i]
      old.customer <- st$customer[i]
      conn.i <- connections(i, st$customer)
      if (!(old.customer %in% conn.i)) st$cluster[conn.i] <- i
      st$customer[i] <- i
      
      # update the likelihoods if a table has been split
      if (old.customer != i & !(old.customer %in% conn.i))
      {
        old.idx <- st[which(st$cluster==old.cluster),"idx"]
        lhood[as.character(old.cluster)] <- llgp(D[old.idx,old.idx,,drop=FALSE], y[old.idx], 
                                                 theta = res_gp[old.cluster,-(ndim+1)],
                                                 nug = res_gp[old.cluster, ndim+1])
      }
      else
      {
        lhood[as.character(old.cluster)] <- 0
      }
      
      # get log prior
      
      log.prior <- logpriorlist[[i]]
      cand.links <- candlinkslist[[i]]
      
      # compute the likelihood of data point i (and its connectors)
      # with all other tables
      
      cand.clusts <- unique(st$cluster[cand.links])
      new.lhood <- c()
      for (k in 1:length(cand.clusts)) {
        df <- st[st$cluster == cand.clusts[k],]
        pts <- unique(c(df$idx,st[conn.i,"idx"]))
        new.lhood[k] <- llgp(D[pts,pts,,drop=FALSE], y[pts],
                             theta = res_gp[cand.clusts[k],-(ndim+1)],
                             nug = res_gp[cand.clusts[k], ndim+1])
      }
      
      names(new.lhood) <- cand.clusts
      
      # set up the old likelihoods
      
      old.lhood <- lhood[as.character(cand.clusts)]
      sum.old.lhood <- sum(old.lhood)
      
      # compute the conditional distribution
      log.prob <-
        log.prior[cand.links] +
        sapply(cand.links,
               function (j) {
                 c.j <- as.character(st$cluster[j])
                 sum.old.lhood - old.lhood[c.j] + new.lhood[c.j] })
      
      prob <- exp(log.prob - log.sum(log.prob))
      if(sum(prob)==0) {
        print('log.prob:')
        print(log.prob)
      }
      if (length(prob)==1)
        st$customer[i] <- cand.links[1]
      else
        st$customer[i] <- sample(cand.links, 1, prob=prob)
      
      # update the clusters
      st$cluster[conn.i] <- st$cluster[st$customer[i]]
      pts_i <- subset(st, cluster == st$cluster[i])$idx
      clus_i <- st$cluster[i]
      res_gp[pts_i,] <- as.list(res_gp[st$customer[i],])
      lhood[as.character(st$cluster[i])] <- llgp(D[pts_i,pts_i,,drop=FALSE], y[pts_i],
                                                 theta = res_gp[clus_i,-(ndim+1)],
                                                 nug = res_gp[clus_i, ndim+1])
    }
    #update Nclust, n_k, z
    z <- st$cluster
    n_k <- as.vector(table(z))
    Nclust <- length(n_k)
    
    # Fit a GP to each cluster
    res_noise[,j+1] <- 0
    Qlist <- list()
    indexlist <- split(seq_along(z), z)
    fitslice <- function(zw){
      myllgp <- function(x){llgp(D=D[zw,zw,,drop=FALSE], y=y[zw], theta=x[1], nug=x[2])}
      draws <- matrix(c(1,0.1), nrow = 1000, ncol = 2, byrow = T)
      for (k in 2:1000) {
        setTimeLimit(elapsed = 0.5, transient = TRUE)
        on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
        
        tryCatch({slice_elliptical_mv(x = draws[k-1,], log_target=myllgp, 
                                      mu=c(2,1), Sig = SigM)$x},
                 error=function(msg){ print(msg) }) 
        drawk <- tryCatch({slice_elliptical_mv(x = draws[k-1,], log_target=myllgp, 
                                               mu=c(2,1), Sig = SigM)$x},
                          error=function(msg){ draws[k-1,] }) 
        draws[k,] <- ifelse(drawk>0, drawk, 0.001)
      }
      setTimeLimit(elapsed = Inf, transient = TRUE)
      return(colMeans(draws[501:1000,]))
    }
    newGPparms <- future_lapply(indexlist, fitslice, future.seed = TRUE)
    
    for (i in 1:Nclust) {
      zw <- indexlist[[i]]
      res_noise[zw,j+1] <- newGPparms[[i]][ndim+1]
      res_gp[i,] <- newGPparms[[i]]
      Qlist[[i]] <- solve(covmat(D=D[zw,zw,,drop=FALSE], 
                                 ls=newGPparms[[i]][-(ndim+1)], 
                                 nug = newGPparms[[i]][ndim+1]))
    }
    
    colnames(res_noise)[j+1] <- paste('nug',j,sep = '')
    
    # Sample from posterior of y
    if(!is.null(Xpred) & j%%draw == 0){
      CC <- exp(-rowSums((dd^2)/rep(phi, each=nrow(Xpred)*N), dims = 2))
      finalpreds <- c()
      predsbyc <- matrix(nrow=nrow(CC), ncol = Nclust)
      for(k in 1:Nclust){
        z_k <- which(z==k)
        pars <- as.numeric(res_gp[k,])
        Qy <- drop(covmat(D=dd, set1 = 1:nrow(dd), set2 = z_k,
                          ls = pars[-length(pars)]))
        pred <- Qy%*%(Qlist[[k]])%*%y[z_k]
        predsbyc[,k] <- as.numeric(pred)
      }
      for (i in 1:nrow(CC)) {
        CCi <- CC[i,]
        #CCi[i] <- log(alpha)
        #CCi <- safelog(CCi) - safelog(sum(CCi))
        probs <- as.numeric(tapply(CCi, z, sum))
        probs <- probs/sum(probs) #exp(lprobs)/sum(exp(lprobs))
        finalpreds <- c(finalpreds, sum(predsbyc[i,]*probs))
      }
      res_drawmu[,j/draw] <- finalpreds
      colnames(res_drawmu)[j/draw] <- paste('draw',j/draw,sep = '')
    }
    
    # Sample alpha
    lalpha <- function(x){ out <- logpost_alpha(alpha = x, N=N, k=Nclust)
    out <- ifelse(!is.finite(out), -1000, out)
    return(out)}
    pseu <- list(ld = function(x) dgamma(x,1,1,log = TRUE), 
                 q = function(x) qgamma(x,1,1))
    setTimeLimit(elapsed = 5, transient = TRUE)
    on.exit(setTimeLimit(elapsed = Inf, transient = FALSE))
    proposal <- tryCatch({ slice_quantile(alpha, lalpha, pseudo = pseu)$x}, 
                         error=function(msg){ alpha })
    alpha <- ifelse(is.finite(proposal), proposal, alpha)
    setTimeLimit(elapsed = Inf, transient = TRUE)
    
    # Sample phi
    logpost_phi <- function(phi, D=D, z=z, alpha=alpha, Xq=NULL, lp=TRUE){
      n <- length(z)
      C <- exp(-rowSums((D^2)/rep(phi, each=n^2), dims = 2))
      step1 <- 0
      for (i in 1:nrow(C)) {
        Ci <- C[i,]
        Ci[i] <- log(alpha)
        Ci <- safelog(Ci) - safelog(sum(Ci))
        lprobs <- as.numeric(tapply(Ci, z, sum))
        probs <- exp(lprobs)/sum(exp(lprobs))
        step1 <- step1 + log(probs[which(unique(z)==z[i])])
      }
      step2 <- max(-1000, log(phi))
      step3 <- max(-1000, log(phi*sqrt(2*pi)))
      step4 <- step1 - 0.5*sum(step2^2) - step3
      if(!lp){ step4 <- exp(step4) }
      return(step4)
    }
    
    prop_cmat <- hessian(logpost_phi, x=phi, D=D, z=z, alpha=alpha, lp=FALSE)
    prop_cmat <- ifelse(!is.finite(prop_cmat), 0, prop_cmat)
    sdphi <- (2.38^2/length(phi))*abs(solve(prop_cmat + diag(1e-9, length(phi), length(phi))))
    for(k in 1:ncol(sdphi)){
      sdphi[k,k] <- min(maxranges[k], sdphi[k,k])
      sdphi[k,k] <- max(1e-5, sdphi[k,k])
    }
    proposal <- rep(-1, length(phi))
    while(!(all(proposal > 0))){proposal <- mvtnorm::rmvnorm(1, phi, as.matrix(sdphi))}
    aprob <- logpost_phi(proposal, alpha=alpha, z = z, D = D) -
      logpost_phi(phi = phi, alpha = alpha, z = z, D = D)
    u <- runif(1)
    aprob <- ifelse(is.nan(aprob), -Inf, ifelse(is.infinite(aprob),-Inf,aprob))
    if(u < exp(aprob)){ phi = proposal }
    
    # Record results
    res_gatepar[j+1,] <- c(alpha, phi)
    res_z[,j+1] <- z
    colnames(res_z)[j+1] <- paste('z',j,sep = '')
    
    if(j%%10 == 0){
      print(paste("Iter",j,"done."))
    }
  }
  
  # Return the results
  spv <- c()
  clusts <- unique(z)
  for (i in clusts) {
    zi <- which(z==i)
    Kinv <- solve(covmat(D[zi,zi,,drop=FALSE], ls=res_gp[i,-ncol(res_gp)],
                         nug = res_gp[i,ncol(res_gp)]))
    tau2 <- drop(t(y[z==i])%*%Kinv%*%y[z==i])/sum(z==i)
    spv <- c(spv, tau2)
  }
  res_gp$spvar <- spv
  res_posterior <- 'none'
  if (!is.null(Xpred)){
    res_posterior <- res_drawmu
  }
  return(list('z'=res_z, 'gatepar'=res_gatepar, 'expertprior'=ab, 
              'noise'=res_noise, 'experts'=res_gp, 'draws'=res_posterior))
}

# DDCRP Clustering Algorithm
ddcrp.gibbs <- function(dat, y, alpha, dist.fn, decay.fn, lhood.fn,
                        niter, summary.fn = ncomp.summary,
                        log.prior.thresh=-10,
                        clust.traj=FALSE, cust.traj=FALSE)
{
  require(plyr)
  ### set up summary statistics and trajectories
  
  ndata <- dim(dat)[1]
  msg.inc <- 10^(floor(log10(dim(dat)[1]))-1)
  if (clust.traj)
    clust.traj <- matrix(NA, nrow=niter, ncol=ndata)
  if (cust.traj)
    cust.traj <- matrix(NA, nrow=niter, ncol=ndata)
  score <- numeric(niter)
  map.score <- 0
  
  ### set up initial state, summaries, and cluster likelihoods
  
  msg("setting up the initial state")
  st <- data.frame(idx=1:ndata, cluster=1:ndata, customer=1:ndata)
  lhood <- daply(st, .(cluster), function (df) lhood.fn(dat[df$idx,,drop=FALSE], y[df$idx]))
  
  summary <- summary.fn(dat, 0, st, lhood, alpha)
  
  ### run for niter iterations
  
  for (iter in 1:niter)
  {
    msg(sprintf("iter=%d", iter))
    
    iter.score <- 0
    for (i in seq(2,ndata)) # note: index i = 1 is correct at the outset
    {
      #if ((i %% msg.inc) == 0) msg(sprintf("%04d", i))
      
      ### "remove" the i-th data point from the state
      ### to do this, set its cluster to i, and set its connected data to i
      
      old.cluster <- st$cluster[i]
      old.customer <- st$customer[i]
      conn.i <- connections(i, st$customer)
      st$cluster[conn.i] <- i
      st$customer[i] <- i
      
      ### if this removal splits a table update the likelihoods.
      ### note: splits only happen if c_i^{old} != i
      
      if (old.customer != i)
      {
        ### !!! do we need to use "idx"
        old.idx <- st[which(st$cluster==old.cluster),"idx"]
        lhood[char(old.cluster)] <- lhood.fn(dat[old.idx,,drop=FALSE], y[old.idx])
      }
      else
      {
        lhood[char(old.cluster)] <- 0
      }
      
      ### compute the log prior
      ### (this should be precomputed---see opt.ddcrp.gibbs below)
      
      log.prior <- sapply(1:ndata,
                          function (j) safelog(decay.fn(dist.fn(dat[i,], dat[j,]))))
      log.prior[i] <- log(alpha)
      log.prior <- log.prior - log.sum(log.prior)
      cand.links <- which(log.prior > log.prior.thresh)
      
      ### compute the likelihood of data point i (and its connectors)
      ### with all other tables (!!! do we need to use "idx"?)
      
      cand.clusts <- unique(st$cluster[cand.links])
      
      new.lhood <- daply(subset(st, cluster %in% cand.clusts), .(cluster),
                         function (df)
                           lhood.fn(dat[unique(c(df$idx,st[conn.i,"idx"])),,drop=FALSE],
                                    y[unique(c(df$idx,st[conn.i,"idx"]))]))
      
      if (length(new.lhood)==1) names(new.lhood) <- cand.clusts
      
      ### set up the old likelihoods
      
      old.lhood <- lhood[char(cand.clusts)]
      sum.old.lhood <- sum(old.lhood)
      
      ### compute the conditional distribution
      
      log.prob <-
        log.prior[cand.links] +
        sapply(cand.links,
               function (j) {
                 c.j <- char(st$cluster[j])
                 sum.old.lhood - old.lhood[c.j] + new.lhood[c.j] })
      
      ### sample from the distribution
      
      prob <- exp(log.prob - log.sum(log.prob))
      if (length(prob)==1)
        st$customer[i] <- cand.links[1]
      else
        st$customer[i] <- sample(cand.links, 1, prob=prob)
      
      ### update the score with the prior and update the clusters
      
      iter.score <- iter.score + log.prior[st$customer[i]]
      st$cluster[conn.i] <- st$cluster[st$customer[i]]
      clust.i.idx <- subset(st, cluster == st$cluster[i])$idx
      lhood[char(st$cluster[i])] <- lhood.fn(dat[clust.i.idx,,drop=FALSE], y[clust.i.idx])
    }
    
    ### update the summary
    
    iter.score <- iter.score + sum(lhood)
    score[iter] <- iter.score
    if ((score[iter] > map.score) || (iter==1))
    {
      map.score <- score[iter]
      map.state <- st
    }
    summary <- rbind(summary, summary.fn(dat, iter, st, lhood, alpha))
    if (!is.null(dim(cust.traj))) cust.traj[iter,] <- st$customer
    if (!is.null(dim(clust.traj))) clust.traj[iter,] <- st$cluster
  }
  
  ### return everything
  
  list(summary=summary, cust.traj=cust.traj, clust.traj=clust.traj, score=score,
       map.score = map.score, map.state = map.state)
}

slc.fixa <- imgpe.slice(X=simX, y=simdat$y, draw = 5, Xpred = Xtest, maxIters = 3000)
set.seed(102925)
x1 <- runif(150,0,5)