### This is the marginal likelihood function to be maximized
# it is weighted if weight var is not null
#par - the score vector
#getIndLnL - logical, should the likelihood be calculated only for a subset
#X_ - the matrix of X values, defaults to all X , used for calculating group likelihoods
#i - list of indexes for individual likelihood (used to subset weights and r22)
#repweightvar  - the replicate weight for replicate variance estimation
fn.regression <- function(X_, i=NULL, wv=NULL, rr1, stuDat, nodes) {
  K <- ncol(X_)
  # fix rr1 to just regard this data
  if(!is.null(i)) {
    rr1p <- rr1[,i,drop=FALSE]
  } else{
    rr1p <- rr1
    i <- 1:nrow(X_)
  }
  # fix weights
  if(!is.null(wv)) {
    w <- stuDat[i, wv]
  } else {
    w <- rep(1, nrow(X_))
  }
  dTheta <-nodes[2] - nodes[1]
  function(par) {
    # break up par into beta and residual components
    B <- par[1:K]
    s2 <- par[(K+1)]^2
    # form prediction
    XB <- X_ %*% B
    # residual, per indivudal (outer sapply), per node (inner sapply)
    nodes.minus.XB <- t(t(matrix(nodes,nrow=length(nodes),ncol=nrow(XB))) - as.vector(XB))
    # likelihood of normal distribution for residuals nodes.minus.XB
    rr2 <- rr1p * ((1/(sqrt(2 * pi * s2))) * exp(-((nodes.minus.XB)^2/(2 * s2))))
    # aggregate likelihood, weight, multiply by -2 to make it a deviance
    cs <- pmax(.Machine$double.eps, colSums(rr2))
    return(-2*sum(w * log(cs*dTheta)))
  }
}

# graded response model
grm <- function (theta, d, score, a, D) {
  maxD <- length(d)   
  if(score == 0) {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)])))
  }
  else if(score == maxD) {
    pr <- 1/(1 + exp(-D*a*(theta - d[score])))
  } else {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)]))) - 1/(1 + exp(D*a*(theta - d[score])))
  }
  pr
}

# log of density of binomial where size=1, accepts x in [0,1] for partial credit
ldbinom2 <- function(x,pr) {
  return(x*log(pr) + (1-x)*log(1-pr))
}

# graded partial credit model
gpcm <- function (theta, d, score, a, D) {
  if (score > length(d)) {
    stop (paste0("Score of ", score," higher than maximum (", length(d),")"))
  }
  if(score <= 0) {
    stop (paste0("Score of ", score," lower than minimum (1)"))
  }
  Da <- D * a
  exp(sum(Da * (theta - d[1:score]))) / sum(exp(cumsum(Da * (theta - d))))
}

# helper, get the graded partial credit model likelihood
gpcmLikelihood <- function (theta, d, score, a, D=1.7) {
  exp(sum(D * (theta - d[1:score]))) / sum(exp(cumsum(D * (theta - d))))
}

# helper, get the graded response model likelihood
grmLikelihood <- function (theta, d, score, a, D=1.7) {
  maxD <- length(d)   
  if(score == 0) {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)])))
  }
  else if(score == maxD) {
    pr <- 1/(1 + exp(-D*a*(theta - d[score])))
  } else {
    pr <- 1/(1 + exp(D*a*(theta - d[(score+1)]))) - 1/(1 + exp(D*a*(theta - d[score])))
  }
  pr
}
