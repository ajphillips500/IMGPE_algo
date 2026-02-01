# Generate lists of all possible partitions
# N options: 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975
nparts <- data.frame('n'=1:11, 'prt'=c(1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570))
nparts$logprt <- log10(nparts$prt)

partitions <- list('n2'=matrix(c(1,1,1,2),2,2))
for (i in 1:8) {
  n <- i + 2
  pset <- matrix(0, nrow = 1, ncol = n)
  prevset <- partitions[[i]]
  
  for (j in 1:nrow(prevset)) {
    options <- c(unique(prevset[j,]), n)
    set_j <- matrix(rep(prevset[j,], each=length(options)),
                    nrow = length(options), ncol = n-1)
    set_j <- cbind(set_j, options)
    colnames(set_j) <- c()
    pset <- rbind(pset, set_j)
  }
  
  partitions[[i+1]] <- pset[-1,]
  names(partitions)[i+1] <- paste0('n',n)
}

# P(z|alpha, phi)
# returns un-normalized log prob of cluster assn. given alpha and kernel matrix
probza <- function(z, alpha, C){
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

# P(alpha|z)
# returns log prob of alpha, given z. Calculates normalizing constant.
lpalpha <- function(alpha, z, C){
  n <- length(z)
  probs <- apply(partitions[[n-1]], 1, probza, alpha, C)
  probs <- exp(probs)
  probs <- probs/sum(probs) #NC neutralized
  zid <- which(rowSums(partitions[[n-1]] == rep(z, each = nrow(partitions[[n-1]])))==n)
  return(log(probs[zid]) + dgamma(alpha, 1,1, log = TRUE))
}

# Test functions
X <- matrix(runif(5), nrow = 5)
D <- as.matrix(dist(X, diag = TRUE, upper = TRUE))
C <- exp(-D)

# data frame of p(a|z) for each z, n=5
# plot sample of dists for same number of clusters, same num singletons
n5dists <- data.frame('x'=1:500/100)
for (i in 5:6) {
  n7dists$add <- vapply(n7dists$x, lpalpha, FUN.VALUE = c(0), z=zz[i,], C=C2)
  colnames(n7dists)[i+1] <- paste0('z', paste0(zz[i,], collapse = ''))
}

#If nsingleton changes then the shape of the log distribution does not change.
#If only nclust changes the shape does not change.