print.mml.means <- function(object, ...){
  co <- mml.remap.coef(object$coefficients, object$location, object$scale)
  print(co, ...)
}

mml.remap.coef <- function(coefficients, location, scale, noloc=FALSE) {
  coefficients <- coefficients * scale
  if(!noloc) {
    coefficients[names(coefficients) == "(Intercept)"] <- coefficients[names(coefficients) == "(Intercept)"] + location
  }
  return(coefficients)
}

summary.mml.means <- function(object, gradientHessian=FALSE,
                              varType=c("consistent", "robust", "cluster", "replicate", "Taylor"),
                              clustervar=NULL, jkSumMultiplier=1, # cluster
                              rep_weight=NULL, # replicate
                              stratavar=NULL, psuvar=NULL, singletonFix=c("drop", "use mean", "pool", "with next", "with next clean"),# Taylor
                              ...){
  sumCall <- match.call()
  H_B_prime <- getIHessian.mml.means(object, gradientHessian)
  # check/fix varType argument
  varType <- match.arg(varType)
  singletonFix <- match.arg(singletonFix)
  if(varType=="consistent") {
    VC <- getVarConsistent(object, H_B_prime)
  }
  if(varType=="robust") {
    VC <- getVarRobust(object, H_B_prime)
  }
  if(varType=="cluster") {
    stuDat <- object$stuDat
    if(is.null(clustervar)) {
      stop("You must define a valid clustervar to use cluster variance estimation.")
    }
    if(length(clustervar) != 1) {
      if("ClusterVar__" %in% colnames(stuDat)) {
        stop("Please rename the variable ", dQuote("ClusterVar__"), " on the ", dQuote("stuDat"), " argument.")
      }
      # paste together variables with colons
      # first, remove colons from existing variables, if necessary
      for(i in 1:length(clustervar)) {
        if(inherits(stuDat[clustervar], "character") && sd(nchar(stuDat[clustervar])) > 0) {
          stuDat[clustervar] <- gsub(":", "-", stuDat[clustervar], fixed=TRUE)
        }
        if(inherits(stuDat[clustervar], "factor")) {
          stuDat[clustervar] <- gsub(":", "-", as.character(stuDat[clustervar]), fixed=TRUE)
        }
      }
      stuDat$ClusterVar__ <- apply(stuDat[clustervar], 1, function(x) { paste(x, collapse=":") } )
      clustervar <- "ClusterVar__"
    }
    if(!clustervar %in% colnames(stuDat)) {
      stop(paste0("Could not find clustervar column named ", dQuote(clustervar), " on ", dQuote("stuDat"), " data"))
    }
    if(length(unique(stuDat[,clustervar])) <= 1) {
      stop("There must be more than one cluster for cluster variance estimation.")
    }
    VC <- getVarCluster(object, H_B_prime, clustervar)
  }
  if(varType=="replicate") {
    if(is.null(rep_weight)) {
      stop(paste0("the argument ", dQuote("rep_weight"), " must be defined for varType ", dQuote("replicate"), "."))
    }
    if( any(! rep_weight %in% colnames(stuDat))) {
      stop(paste0("Could not find rep_weight column named ", dQuote(rep_weight), " on ", dQuote("stuDat"), " data."))
    }
    VC <- getVarReplicate(object, H_B_prime, rep_weight, jkSumMultiplier)
  }
  if(varType=="Taylor") {
    stuDat <- object$stuDat
    if(is.null(stratavar) | is.null(psuvar)) {
      stop(paste0("the arguments ", dQuote("stratavar"), " and ", dQuote("psuvar")," must be defined for varType ", dQuote("Taylor"), "."))
    }
    if(!stratavar %in% colnames(stuDat)) {
      stop(paste0("Could not find stratavar column named ", dQuote(stratavar), " on ", dQuote("stuDat"), " data."))
    }
    if(!psuvar %in% colnames(stuDat)) {
      stop(paste0("Could not find stratavar column named ", dQuote(psuvar), " on ", dQuote("stuDat"), " data."))
    }
    VC <- getVarTaylor(object, H_B_prime, stratavar, psuvar, singletonFix)
  }

  se <- sqrt(diag(VC))
  tval <- as.vector(mml.remap.coef(coef(object), object$location, object$scale))/mml.remap.coef(se, object$location, object$scale, noloc=TRUE)
  TAB <- cbind(Estimate = as.vector(mml.remap.coef(coef(object), object$location, object$scale)),
               StdErr = mml.remap.coef(se, object$location, object$scale, noloc=TRUE),
               t.value = tval)
  row.names(TAB) <- names(object$coefficients) 
  return(structure(list("call" = object$call,
                        "summaryCall" = sumCall,
                        "coefficients" = TAB,
                        "converged" = object$Convergence,
                        "LL" = object$LogLik,
                        "convergence" = object$convergence,
                        "iterations" = object$iterations,
                        "VC" = VC,
                        "iHessian" = H_B_prime,
                        "stuDat" = object$stuDat,
                        "X" = object$X,
                        "obs" = object$obs),
                   class="summary.mml.means"))
}

#' @importFrom stats printCoefmat
print.summary.mml.means <- function(x, ...){
  cat(paste0("Call:\n"))
  print(x$call)
  cat(paste0("Summary Call:\n"))
  print(x$summaryCall)
  cat("\n")
  cat("Summary:\n")
  cof <- x$coefficients
  cof1 <- cof[1:(nrow(cof)-1),,drop=FALSE]
  cof2 <- cof[nrow(cof),1:2,drop=FALSE]
  printCoefmat(cof1)
  cat("\n")
  cat("Residual Variance Estimate:\n")
  print(cof2)
  cat("\n")
  cat(paste0("Convergence = ", x$converged, "\n"))
  cat(paste0("Iterations = ", x$iterations, "\n"))
  cat(paste0("LogLike = ", round(x$LL,2), "\n"))
  cat(paste0("observations = ", x$obs, "\n"))
} 

getIHessian.mml.means <- function(object, gradientHessian=FALSE) {
  if(gradientHessian) {
    X <- object$X
    x_ind <- split(X, as.factor(1:nrow(X)), drop=FALSE)
    # dl/dBeta * dl/dBeta ' evaluated for each individual 
    stuDat <- object$stuDat
    stuDat$one <- 1
    weightvar <- object$weightvar
    if(is.null(object$weightvar)) {
      weightvar <- "one"
    }
    vars <- lapply(1:nrow(X), FUN=function(i){stuDat[i,weightvar] * gradgradT(location=coef(object),
                                                        ii=i,
                                                        X_subset=matrix(x_ind[[i]], nrow=1),
                                                        weightvar="one",
                                                        rr1=object$rr1,
                                                        stuDat=stuDat,
                                                        nodes=object$nodes)})
    Wsum <- sum(stuDat[,weightvar])
    H_B_prime <- -1 * Wsum/(Wsum-1) * solve(Reduce("+", vars))
  } else {
    H_B_prime <- solve(-1/2*getHessian(object$lnlf, object$coefficients))
  }
  return(H_B_prime)
}

# estimate covariance matrix (Standard error)
# -1/2 turns a deviance into a likelihood
getVarConsistent <- function(object, H_B_prime) {
  #consistent  estimate
  return(-1 * H_B_prime )
}

getVarRobust <- function(object, H_B_prime) {
  X <- object$X
  x_ind <- split(X, as.factor(1:nrow(X)), drop=FALSE)
  # dl/dBeta * dl/dBeta ' evaluated for each individual 
  object$stuDat$one <- 1
  vars <- lapply(1:nrow(X), FUN=function(i){gradgradT(location=coef(object),
                                                      ii=i,
                                                      X_subset=matrix(x_ind[[i]], nrow=1),
                                                      #weightvar=object$weightvar,
                                                      weightvar="one",
                                                      rr1=object$rr1,
                                                      stuDat=object$stuDat,
                                                      nodes=object$nodes)})
  V <- Reduce("+", vars)
  return(H_B_prime * V * H_B_prime)
}

getVarCluster <- function(object, H_B_prime, clustervar) {
  stuDat <- object$stuDat
  X <- object$X
  # dl/dBeta * dl/dBeta ' evaluated for each group  before being multipled and then summed
  # first get list of each group index 
  clusters <- unique(stuDat[[clustervar]])
  group_index <- lapply(clusters, FUN=function(x) {
    which(stuDat[[clustervar]]==x)
  })
  #this is important to ensure group_index and x_groups are in the same order
  stuDat[[clustervar]] <- factor(stuDat[[clustervar]], levels=clusters) 
  x_groups <- lapply(split(X, stuDat[[clustervar]]), matrix, ncol=ncol(X))
  vars <- lapply(c(1:length(group_index)), FUN=function(group){
    gradgradT(location=coef(object),
              ii=group_index[[group]],
              X_subset = x_groups[[group]],
              weightvar=object$weightvar,
              rr1=object$rr1,
              stuDat=object$stuDat,
              nodes=object$nodes)
  })
  V <-  Reduce("+",vars)
  return(H_B_prime %*% V %*% H_B_prime)
} 

getVarReplicate <- function(object, H_B_prime, rep_weight, jkSumMultiplier=1) {
  X <- object$X
  B0  <- object$coefficients
  B_j <- lapply(rep_weight, FUN=function(x){ #restimate with each weight 
    fn2B <- fn.regression(X_=X, wv=x, rr1=object$rr1)
    return(bobyqa(B0, fn2B)$par)
  })
  rep <- lapply(B_j, function(x){(x-B0) %*% t(x-B0)})
  return(jkSumMultiplier * Reduce("+", rep))
} 

getVarTaylor <- function(object, H_B_prime, stratavar, psuvar, singletonFix=c("drop", "use mean", "pool", "with next", "with next clean")) {
  #find PSU per strata, and warn about dropping strata with one PSU
  stuDat <- object$stuDat
  X <- object$X
  strata <- lapply(sort(unique(stuDat[[stratavar]])), FUN=function(x) {
    list(strat=x,
         psu=sort(unique(stuDat[stuDat[,stratavar]==x, psuvar])))
  })
  n_psu <- lapply(strata, function(x) { length(x$psu)}) 
  if (any(n_psu==1)){
    if(singletonFix == "with next") {
      warning(paste0("Of the ", length(n_psu)," strata, ", sum(n_psu<2) ," strata have only one PSU. All strata with only one PSU will be combined with the next stratum, maintaining their PSU identifier. See the ", dQuote("singletonFix"), " argument for more details and other options."))
      ustrata <- sort(unique(stuDat[[stratavar]]))
      last_occupied <- NA
      for(i in 1:(length(ustrata)-1)) {
        if(n_psu[i] == 1) {
          stuDat[stuDat[,stratavar] == ustrata[i],stratavar] <- ustrata[i+1]
        } else{
          last_occupied <- i
        }
      }
      if(n_psu[length(n_psu)] == 1) {
        stuDat[stuDat[,stratavar] == ustrata[length(n_psu)],stratavar] <- ustrata[last_occupied]
      }
      strata <- strata[n_psu>1] #variance estimation can only happen for Strata with more than one PSU 
    }
    if(singletonFix == "with next clean") {
      warning(paste0("Of the ", length(n_psu)," strata, ", sum(n_psu<2) ," strata have only one PSU. All strata with only one PSU will be combined with the next stratum, after being given a PSU identifier that is not occupid in the new stratum. See the ", dQuote("singletonFix"), " argument for more details and other options."))
      ustrata <- sort(unique(stuDat[[stratavar]]))
      last_occupied <- NA
      for(i in 1:(length(ustrata)-1)) {
        if(n_psu[i] == 1) {
          stuDat[stuDat[,stratavar] == ustrata[i],psuvar] <- stuDat[stuDat[,stratavar] == ustrata[i],psuvar] + max(stuDat[,psuvar])
          stuDat[stuDat[,stratavar] == ustrata[i],stratavar] <- ustrata[i+1]
        } else{
          last_occupied <- i
        }
      }
      if(n_psu[length(n_psu)] == 1) {
        i <- length(n_psu)
        stuDat[stuDat[,stratavar] == ustrata[i],psuvar] <- stuDat[stuDat[,stratavar] == ustrata[i],psuvar] + max(stuDat[,psuvar])
        stuDat[stuDat[,stratavar] == ustrata[i],stratavar] <- ustrata[last_occupied]
      }
      strata <- strata[n_psu>1] #variance estimation can only happen for Strata with more than one PSU 
    }
    if(singletonFix == "pool") {
      warning(paste0("Of the ", length(n_psu)," strata, ", sum(n_psu<2) ," strata have only one PSU. All strata with only one PSU have their value compared to the mean of these strata. See the ", dQuote("singletonFix"), " argument for more details and other options."))
    }
    if(singletonFix == "drop") {
      warning(paste0("Of the ", length(n_psu)," strata, ", sum(n_psu<2) ," strata have only one PSU. All strata with only one PSU are excluded from variance estimation. See the ", dQuote("singletonFix"), " argument for other options."))
      strata <- strata[n_psu>1] #variance estimation can only happen for Strata with more than one PSU 
    }
    if(singletonFix == "use mean") {
      warning(paste0("Of the ", length(n_psu)," strata, ", sum(n_psu<2) ," strata have only one PSU. All strata with only one PSU have their value compared to the mean. See the ", dQuote("singletonFix"), " argument for more details and other options."))
    }
  }
  #split X based on psu for access later
  # loop through strata 
  str <- lapply(strata, FUN=function(st) {
    # number of PSUs in this stratum
    n_a <- length(st$psu)
    # PSU index for units in this stratum and PSU
    group_index <- lapply(st$psu, FUN=function(x) {
      which(stuDat[[psuvar]]==x & stuDat[[stratavar]]==st$strat)
    })
    # extract data for this stratum
    X_strata <- X[stuDat[[stratavar]] %in% st$strat,]
    stuDat_strata <- stuDat[stuDat[[stratavar]] %in% st$strat,]
    #split up data by PSU, lapply after the split, making each one a matrix
    x_groups <- lapply(split(X_strata, stuDat_strata[[psuvar]]), matrix, ncol=ncol(X))
    #vector of scores per psu
    s_p <- lapply(c(1:length(group_index)),FUN=function(k){
      gradInd(location=coef(object),
              ii=group_index[[k]],
              X_subset=x_groups[[k]],
              weightvar=object$weightvar,
              rr1=object$rr1,
              stuDat=object$stuDat,
              nodes=object$nodes)
    })
    st$s_p <- s_p 
    if(n_a > 1) {
      #average score across as psu in strata
      s_a_bar <- Reduce("+",s_p)/n_a
      #(s_p - s_a_bar)*(s_p - s_a_bar)' 
      v_a <- lapply(s_p, FUN=function(s){
        (s - s_a_bar) %*% t(s - s_a_bar)
      })
      st$V_a <- (n_a/(n_a-1))*Reduce("+",v_a)
    }
    return(st)
  }) # end lapply(strata, FUN=function(st) {
  if(singletonFix %in% c("use mean", "pool")) {
    # get all strata's PSUs to find mean
    s_p_all <- list()
    for(stri in str) {
      if( (length(stri$s_p) == 1) | singletonFix == "use mean") {
        s_p_all <- c(s_p_all, stri$s_p)
      }
    }
    browser()
    # this is the overall average, across strata, should be zero
    n_all <- length(s_p_all)
    s_a_barOverall <- Reduce("+", s_p_all)/n_all
  }
  if(singletonFix=="drop") {
    # when we reduce these, we will need a zero matrix
    # so grab a valid V_a and multiply it by 0
    zeroV_a <- NULL
    i <- 1
    while(is.null(zeroV_a)) {
      if(!is.null(str[[i]]$V_a)) {
        zeroV_a <- 0 * str[[i]]$V_a
      }
      i <- i + 1
    }
  }
  # aggregate V_a
  V_a <- lapply(str, FUN=function(st) {
    npsu <- length(st$psu)
    if(npsu > 1) {
      return(st$V_a)
    }
    if(singletonFix %in% c("use mean", "pool")) {
      #(s_p - s_a_barOverall)*(s_p - s_a_barOverall)' 
      s <- st$s_p[[1]] # there is only one element
      v_a <- 2*(s - s_a_barOverall) %*% t(s - s_a_barOverall)
      return(v_a)
    }
    if(singletonFix=="drop") {
      return(zeroV_a)
    }
  })
  #sum variance across over all strata

  V <- Reduce("+",V_a)
  return(H_B_prime %*% V %*% H_B_prime)
}

# from EdSurvey
# author: Paul Bailey
pasteItems <- function(vector, final="and") {
  # no need to do anything if there is one or fewer elements
  if(length(vector) <= 1) {
    return(vector)
  }
  if(length(vector) == 2) {
    return(paste0(vector[1], " ", final, " ", vector[2]))
  }
  v <- vector[-length(vector)]
  f <- vector[length(vector)]
  return(paste0(paste(v, collapse=", "), ", ", final, " ", f))
}

#' @importFrom stats vcov
vcov.mml.means <- function(object, ...){
  vcov(summary(object, ...))
}

vcov.summary.mml.means <- function(object, ...){
  object$VC
}
